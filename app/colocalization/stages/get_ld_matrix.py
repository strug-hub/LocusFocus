import os
from typing import Optional, Tuple

import pandas as pd
import numpy as np
from flask import current_app as app

from app.colocalization.constants import LD_MAT_DIAG_CONSTANT
from app.colocalization.payload import SessionPayload
from app.utils import download_file_with_ext
from app.colocalization.plink import plink_ld_pairwise, plink_ldmat
from app.pipeline import PipelineStage
from app.utils.errors import InvalidUsage, ServerError


class GetLDMatrixStage(PipelineStage):
    """
    Adds an LD matrix to the payload for the given data.
    First checks if the user uploaded their own LD, and if not,
    we create an LD matrix with the user's preferred population.

    If no GWAS data is in the session, then an error is raised.

    Will also determine lead SNP and R2 values for the gwas dataset.

    In the case where an LD matrix is created, the gwas_data is subset
    to only the SNPs in the generated LD matrix.

    Prerequisites:
    - Session is created.
    - `gwas_data` is defined in the session.
    """

    def name(self) -> str:
        return "get-ld-matrix"

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        # Enforce prerequisites
        if payload.gwas_data is None:
            raise ServerError(f"Cannot use GetLDMatrixStage; gwas_data is None")

        # Read from file if it exists. Otherwise, create with PLINK
        ld_matrix, _ = self._read_ld_matrix_file(payload)
        if ld_matrix is None:
            ld_matrix, _ = self._create_ld_matrix(payload)

        payload.ld_matrix = ld_matrix

        return payload

    def _read_ld_matrix_file(
        self, payload: SessionPayload
    ) -> Tuple[Optional[np.matrix], Optional[pd.DataFrame]]:
        """
        Try to read in an LD matrix if one is provided by the user.

        If the user did not provide their own LD matrix, then return None.

        If the user provided an LD matrix, but the format is invalid or otherwise unusable,
        then raise an error.

        Otherwise, return the LD matrix, as well as the BIM file as a DataFrame.
        """
        ld_matrix_filepath = download_file_with_ext(payload.request_files, ["ld"])
        if ld_matrix_filepath is None:
            return None, None

        if payload.gwas_data is None:
            raise ServerError(
                f"Cannot validate user-provided LD matrix; gwas_data is not defined"
            )

        ld_mat = pd.read_csv(
            ld_matrix_filepath, sep="\t", encoding="utf-8", header=None
        )
        # Dimensions check
        if not len(ld_mat.shape) == 2:
            raise InvalidUsage(f"LD matrix input is not a 2D matrix", status_code=410)

        if not ld_mat.shape[0] == ld_mat.shape[1]:
            raise InvalidUsage(f"LD matrix input is not square", status_code=410)

        # Compare LD matrix dimensions to GWAS data
        if not ld_mat.shape[0] == payload.gwas_data_kept.shape[0]:
            if ld_mat.shape[0] == payload.gwas_data.shape[0]:
                # Special case: user uploaded LD matrix matches GWAS data exactly, but not after subsetting
                app.logger.warning(
                    f"LD matrix input has same dimensions as GWAS data, but not after subsetting. LD matrix will be subsetted to match current GWAS data."
                )
                ld_mat = ld_mat.loc[
                    payload.gwas_indices_kept, payload.gwas_indices_kept
                ]
            else:
                raise InvalidUsage(
                    f"GWAS and LD matrix input have different dimensions:\nRaw GWAS Length: {payload.gwas_data.shape[0]}\nGWAS Length after format check: {payload.gwas_data_kept.shape[0]}\nLD matrix shape: {ld_mat.shape}",
                    status_code=410,
                )

        # Ensure custom LD matrix and GWAS files are sorted for accurate matching
        if not payload.gwas_data_kept["POS"].is_monotonic_increasing:
            raise InvalidUsage(
                "GWAS data input is not sorted and may not match with the LD matrix",
                status_code=410,
            )

        payload.r2 = list(ld_mat.iloc[:, payload.get_current_lead_snp_index()])  # type: ignore

        ld_mat = np.matrix(ld_mat)

        # Recreate BIM file from PLINK
        # since the user provided their own LD, we assume they correspond to the provided SNPs
        ld_snps_df = pd.DataFrame(
            {
                "CHROM": payload.gwas_data_kept["CHROM"],
                "CHROM_POS": payload.gwas_data_kept["CHROM_POS"],
                "POS": payload.gwas_data_kept["POS"],
                "ALT": payload.gwas_data_kept["ALT"],
                "REF": payload.gwas_data_kept["REF"],
            }
        )

        ld_snps_df.iloc[:, 0] = x_to_23(list(ld_snps_df.iloc[:, 0]))  # type: ignore
        payload.ld_snps_bim_df = ld_snps_df

        return ld_mat, ld_snps_df

    def _create_ld_matrix(
        self, payload: SessionPayload
    ) -> Tuple[np.matrix, pd.DataFrame]:
        """
        Try to create an LD matrix using PLINK.

        Return the LD matrix, and the BIM file as a DataFrame.
        """
        if payload.gwas_data is None:
            raise ServerError(f"Cannot create LD matrix; gwas_data is not defined")

        chrom, _, _ = payload.get_locus_tuple()

        # Create LD matrix with PLINK
        ld_snps_df, ldmat = plink_ldmat(
            build=payload.get_coordinate(),
            pop=payload.get_ld_population(),
            chrom=chrom,
            snp_positions=list(payload.gwas_data_kept["POS"]),
            outfilename=os.path.join(
                app.config["SESSION_FOLDER"], f"ld-{payload.session_id}"
            ),
            region=payload.get_locus_tuple(),
        )

        # Update lead SNP if needed, and set R2
        temp_ld_mat, new_lead_snp_position = plink_ld_pairwise(
            build=payload.get_coordinate(),
            pop=payload.get_ld_population(),
            chrom=chrom,
            snp_positions=list(payload.gwas_data_kept["POS"]),
            snp_pvalues=list(payload.gwas_data_kept["P"]),
            outfilename=os.path.join(
                app.config["SESSION_FOLDER"], f"ld-{payload.session_id}"
            ),
        )

        # Rename to consistent naming structure (see plink .bim file format)
        ld_snps_df = ld_snps_df.rename(
            columns={0: "CHROM", 1: "CHROM_POS", 3: "POS", 4: "ALT", 5: "REF"}
        ).iloc[
            :, [0, 1, 3, 4, 5]
        ]  # drop third column (all zeroes)

        # Generated LD matrix, need to update GWAS indices
        ld_indices = payload.gwas_data["POS"].isin(ld_snps_df["POS"])
        payload.gwas_indices_kept &= ld_indices

        payload.r2 = list(temp_ld_mat["R2"])
        payload.ld_snps_bim_df = ld_snps_df

        return np.matrix(ldmat), ld_snps_df
