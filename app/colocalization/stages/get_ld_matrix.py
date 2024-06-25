import os
from typing import Optional

import pandas as pd
import numpy as np
from flask import current_app as app

from app.colocalization.payload import SessionPayload
from app.colocalization.utils import download_file
from app.colocalization.plink import plink_ld_pairwise
from app.pipeline import PipelineStage
from app.routes import InvalidUsage


class GetLDMatrixStage(PipelineStage):
    """
    Adds an LD matrix to the payload for the given data.
    First checks if the user uploaded their own LD, and if not,
    we create an LD matrix with the user's preferred population.

    If no GWAS data is in the session, then an error is raised.

    Prerequisites:
    - Session is created.
    - `gwas_data` is defined in the session.
    """

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        # Enforce prerequisites
        if payload.gwas_data is None:
            raise Exception(f"Cannot use GetLDMatrixStage; gwas_data is None")

        # First, read from file if it exists
        ld_matrix = self._read_ld_matrix_file(payload)
        if ld_matrix is not None:
            payload.ld_data = ld_matrix
            return payload

        # Second, create matrix using PLINK if first didnt work
        ld_matrix = self._create_ld_matrix(payload)
        payload.ld_data = ld_matrix

        return payload


    def _read_ld_matrix_file(self, payload: SessionPayload) -> Optional[pd.DataFrame]:
        """
        Try to read in an LD matrix if one is provided by the user.

        If the user did not provide their own LD matrix, then return None.

        If the user provided an LD matrix, but the format is invalid or otherwise unusable,
        then raise an error.

        Otherwise, return the LD matrix as a DataFrame.
        """
        ld_matrix_filepath = download_file(payload.request, ['ld'])
        if ld_matrix_filepath is None:
            return None

        if payload.gwas_data is None:
            raise Exception(f"Cannot validate user-provided LD matrix; gwas_data is not defined")

        ld_mat = pd.read_csv(ld_matrix_filepath, sep="\t", encoding='utf-8', header=None)
        if not (len(ld_mat.shape) == 2 and ld_mat.shape[0] == ld_mat.shape[1] and ld_mat.shape[0] == payload.gwas_data.shape[0]):
            raise InvalidUsage(f"GWAS and LD matrix input have different dimensions:\nGWAS Length: {payload.gwas_data.shape[0]}\nLD matrix shape: {ld_mat.shape}", status_code=410)

        # Ensure custom LD matrix and GWAS files are sorted for accurate matching
        if not np.all(payload.gwas_data["POS"][:-1] <= payload.gwas_data["POS"][1:]):
            raise InvalidUsage('GWAS data input is not sorted and may not match with the LD matrix', status_code=410)

        # We have a user-submitted LD matrix. We need to subset it to the GWAS indices that were kept
        ld_mat = ld_mat.loc[ payload.gwas_indices_kept, payload.gwas_indices_kept ]

        payload.r2 = list(ld_mat.iloc[:, payload.gwas_lead_snp_index]) # type: ignore

        return ld_mat


    def _create_ld_matrix(self, payload: SessionPayload) -> pd.DataFrame:
        """
        Try to create an LD matrix using PLINK.
        """
        if payload.gwas_data is None:
            raise Exception(f"Cannot create LD matrix; gwas_data is not defined")

        chrom, _, _ = payload.get_locus_tuple()

        ld_df, new_lead_snp_position = plink_ld_pairwise(
            build=payload.get_coordinate(),
            pop=payload.get_ld_population(),
            chrom=chrom,
            snp_positions=list(payload.gwas_data["POS"]),
            snp_pvalues=list(payload.gwas_data["P"]),
            outfilename=os.path.join(app.config["SESSION_FOLDER"], f"ld-{payload.session_id}")
        )

        old_lead_snp_position = payload.gwas_data.iloc[payload.gwas_lead_snp_index, :]["POS"]

        if new_lead_snp_position != old_lead_snp_position:
            payload.gwas_lead_snp_index = list(payload.gwas_data["POS"]).index(new_lead_snp_position)

        return ld_df
