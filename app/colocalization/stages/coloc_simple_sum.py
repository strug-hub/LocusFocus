import json
from typing import Optional, Tuple
import numpy as np
import pandas as pd

from app.colocalization.constants import LD_MAT_DIAG_CONSTANT
from app.colocalization.payload import SessionPayload
from app.utils import clean_snps, standardize_snps, write_list, write_matrix
from app.pipeline.pipeline_stage import PipelineStage
from app.scripts import ScriptError, coloc2, simple_sum
from app.utils.errors import InvalidUsage, ServerError
from app.utils.gtex import get_gtex_data


class ColocSimpleSumStage(PipelineStage):
    """
    Stage for preparing and running Simple Sum 2 *and* COLOC2 colocalization.

    Prerequisites:
    - gwas_data is loaded and subsetted
    - ld_matrix is loaded
    - gtex_data is reported
    """

    COLOC2_COLNAMES = [
        "CHR",
        "POS",
        "SNPID",
        "A2",
        "A1",
        "BETA",
        "SE",
        "PVAL",
        "MAF",
        "N",
    ]
    COLOC2_EQTL_COLNAMES = COLOC2_COLNAMES + ["ProbeID"]
    COLOC2_GWAS_COLNAMES = COLOC2_COLNAMES + ["type"]

    def name(self) -> str:
        return "simple-sum"
    
    def description(self) -> str:
        return "Perform Simple Sum colocalization (and COLOC2 if selected)"

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        # Check prerequisites
        errors = []
        if payload.gwas_data is None:
            errors.append("GWAS data is not loaded")
        if payload.ld_matrix is None:
            errors.append("LD Matrix is not loaded")
        if payload.ld_snps_bim_df is None:
            errors.append("LD Matrix BIM file is missing")
        if payload.reported_gtex_data == {} and payload.get_gtex_selection() != (
            [],
            [],
        ):
            # reported_gtex_data isn't actually used, but it's just a proxy for whether GTEx was fetched
            errors.append(
                "GTEx genes/tissues were selected but data has not been reported"
            )

        if len(errors) > 0:
            raise ServerError(
                f"Missing prerequisites for Simple Sum Stage: {', '.join(errors)}"
            )

        p_value_matrix, coloc2eqtl_df = self._build_pvalue_matrix(payload)

        if len(p_value_matrix) == 1:
            # Only GWAS data was provided, no secondary datasets
            raise InvalidUsage("No secondary datasets provided for colocalization. Please select at least one GTEx tissue/gene combination, or upload a secondary dataset.\nIf you would like to run a set-based test for your GWAS data, please use the Set-Based Test form in the navbar instead.")

        self._run_simple_sum(
            p_value_matrix, payload, coloc2eqtl_df,
        )

        return payload

    def _build_pvalue_matrix(
        self, payload: SessionPayload
    ) -> Tuple[np.matrix, Optional[pd.DataFrame]]:
        """
        Create a P-value matrix required for Simple Sum colocalization,
        as well as a COLOC2 dataframe if COLOC2 is enabled.

        The P-value matrix is a rectangular matrix with M rows and N columns, where
        N is the number of SNPs being tested, and
        M is the number of datasets.

        N is also equal to the side-length of the square LD matrix.

        The first row contains the P values of the provided GWAS dataset
        (after subsetting, cleaning, etc.).

        The subsequent rows contain the P values of the selected eQTL datasets from GTEx.
        There are {tissues}*{genes} rows for eQTLs, ordered by tissue then by gene.
        The order is determined by what the user selects in the form. Ordering
        is not important but the grouping is.

        Rows after this are the P values for the
        secondary datasets provided by the user, if any.

        Example:
        - User selects GTEx datasets for
            - tissues "Blood" and "Liver"
            - genes "NUCKS1" and "CDK18"
            - uploads secondary datasets Alpha, Bravo, Charlie
        - We expect 8 rows in the matrix that could be labelled as:
            - GWAS
            - Blood, NUCKS1
            - Blood, CDK18
            - Liver, NUCKS1
            - Liver, CDK18
            - Alpha
            - Bravo
            - Charlie
        """
        assert payload.gwas_data is not None

        p_value_matrix = []
        gtex_version = payload.get_gtex_version()
        regionstr = payload.get_locus()
        coordinate = payload.get_coordinate()

        coloc2eqtl_df = pd.DataFrame({})

        # 1. GWAS Data
        p_value_matrix.append(list(payload.gwas_data_kept["P"]))  # type: ignore

        # 2. GTEx secondary datasets
        std_snp_list = pd.Series(
            clean_snps(list(payload.std_snp_list), regionstr, coordinate)
        )
        ss_std_snp_list = std_snp_list.loc[payload.gwas_indices_kept]

        gtex_tissues, gtex_genes = payload.get_gtex_selection()
        
        if len(gtex_tissues) > 0:
            for tissue in gtex_tissues:
                for agene in gtex_genes:
                    gtex_eqtl_df = get_gtex_data(
                        gtex_version, tissue, agene, ss_std_snp_list
                    )
                    if len(gtex_eqtl_df) > 0:
                        pvalues = list(gtex_eqtl_df["pval"])
                        if payload.coloc2:
                            tempdf = gtex_eqtl_df.rename(
                                columns={
                                    "rs_id": "SNPID",
                                    "pval": "PVAL",
                                    "beta": "BETA",
                                    "se": "SE",
                                    "sample_maf": "MAF",
                                    "chr": "CHR",
                                    "variant_pos": "POS",
                                    "ref": "A2",
                                    "alt": "A1",
                                }
                            )
                            tempdf.dropna(inplace=True)
                            if len(tempdf.index) != 0:
                                numsamples = round(
                                    tempdf["ma_count"].tolist()[0]
                                    / tempdf["MAF"].tolist()[0]
                                )
                                numsampleslist = np.repeat(
                                    numsamples, tempdf.shape[0]
                                ).tolist()
                                tempdf = pd.concat(
                                    [tempdf, pd.Series(numsampleslist, name="N")],
                                    axis=1,
                                )
                                probeid = str(tissue) + ":" + str(agene)
                                probeidlist = np.repeat(
                                    probeid, tempdf.shape[0]
                                ).tolist()
                                tempdf = pd.concat(
                                    [tempdf, pd.Series(probeidlist, name="ProbeID")],
                                    axis=1,
                                )
                                tempdf = tempdf.reindex(
                                    columns=self.COLOC2_EQTL_COLNAMES
                                )
                                coloc2eqtl_df = pd.concat(
                                    [coloc2eqtl_df, tempdf], axis=0
                                )
                    else:
                        pvalues = np.repeat(np.nan, len(ss_std_snp_list))
                    p_value_matrix.append(pvalues)

        # 3. Uploaded secondary datasets
        if (
            payload.secondary_datasets is not None
            and len(payload.secondary_datasets) > 0
        ):
            if payload.coloc2:
                # Saving uploaded secondary datasets for coloc2 run
                for (
                    dataset_title,
                    secondary_dataset,
                ) in payload.secondary_datasets.items():
                    secondary_dataset = pd.DataFrame(secondary_dataset)
                    if secondary_dataset.shape[0] == 0:
                        # print(f'No data for table {table_titles[i]}')
                        pvalues = np.repeat(np.nan, len(ss_std_snp_list))
                        p_value_matrix.append(pvalues)
                        continue
                    try:
                        if not set(self.COLOC2_EQTL_COLNAMES).issubset(
                            secondary_dataset
                        ):
                            raise InvalidUsage(
                                f"You have chosen to run COLOC2. COLOC2 assumes eQTL data as secondary dataset, and you must have all of the following column names: {self.COLOC2_EQTL_COLNAMES}"
                            )
                        secondary_dataset["SNPID"] = clean_snps(
                            secondary_dataset["SNPID"].tolist(), regionstr, coordinate
                        )
                        # secondary_dataset.set_index('SNPID', inplace=True)
                        idx = pd.Index(list(secondary_dataset["SNPID"]))
                        secondary_dataset = (
                            secondary_dataset.loc[~idx.duplicated()]
                            .reset_index()
                            .drop(columns=["index"])
                        )
                        # merge to keep only SNPs already present in the GWAS/primary dataset (SS subset):
                        secondary_data_std_snplist = standardize_snps(
                            secondary_dataset["SNPID"].tolist(), regionstr, coordinate
                        )
                        secondary_dataset = pd.concat(
                            [
                                secondary_dataset,
                                pd.DataFrame(
                                    secondary_data_std_snplist, columns=["SNPID.tmp"]
                                ),
                            ],
                            axis=1,
                        )
                        snp_df = pd.DataFrame(ss_std_snp_list, columns=["SNPID.tmp"])
                        secondary_data = (
                            snp_df.reset_index()
                            .merge(
                                secondary_dataset,
                                on="SNPID.tmp",
                                how="left",
                                sort=False,
                            )
                            .sort_values("index")
                        )
                        pvalues = list(secondary_data["PVAL"])
                        p_value_matrix.append(pvalues)
                        coloc2eqtl_df = pd.concat(
                            [
                                coloc2eqtl_df,
                                secondary_data.reindex(
                                    columns=self.COLOC2_EQTL_COLNAMES
                                ),
                            ],
                            axis=0,
                        )
                    except InvalidUsage as e:
                        e.message = f"[secondary dataset '{dataset_title}'] {e.message}"
                        raise e
            else:
                for (
                    dataset_title,
                    secondary_dataset,
                ) in payload.secondary_datasets.items():
                    secondary_dataset = pd.DataFrame(secondary_dataset)
                    if secondary_dataset.shape[0] == 0:
                        # print(f'No data for table {table_titles[i]}')
                        pvalues = np.repeat(np.nan, len(ss_std_snp_list))
                        p_value_matrix.append(pvalues)
                        continue
                    # remove duplicate SNPs
                    try:
                        secondary_dataset["SNP"] = clean_snps(
                            secondary_dataset["SNP"].tolist(), regionstr, coordinate
                        )
                        idx = pd.Index(list(secondary_dataset["SNP"]))
                        secondary_dataset = (
                            secondary_dataset.loc[~idx.duplicated()]
                            .reset_index()
                            .drop(columns=["index"])
                        )
                        # merge to keep only SNPs already present in the GWAS/primary dataset (SS subset):
                        secondary_data_std_snplist = standardize_snps(
                            secondary_dataset["SNP"].tolist(), regionstr, coordinate
                        )
                        std_snplist_df = pd.DataFrame(
                            secondary_data_std_snplist, columns=["SNP" + ".tmp"]
                        )
                        secondary_dataset = pd.concat(
                            [secondary_dataset, std_snplist_df], axis=1
                        )
                        snp_df = pd.DataFrame(ss_std_snp_list, columns=["SNP" + ".tmp"])
                        secondary_data = (
                            snp_df.reset_index()
                            .merge(
                                secondary_dataset,
                                on="SNP" + ".tmp",
                                how="left",
                                sort=False,
                            )
                            .sort_values("index")
                        )
                        pvalues = list(secondary_data["P"])
                        p_value_matrix.append(pvalues)
                    except InvalidUsage as e:
                        e.message = f"[secondary dataset '{dataset_title}'] {e.message}"
                        raise e

        return np.matrix(p_value_matrix), coloc2eqtl_df

    def _run_simple_sum(
        self,
        p_value_matrix: np.matrix,
        payload: SessionPayload,
        coloc2eqtl_df: Optional[pd.DataFrame] = None,
    ):
        """
        Execute the Simple Sum colocalization R script with the given parameters.
        Also executes COLOC2 if specified by the user.

        Return the DataFrame created by the R script as a result.
        Raise error if the script fails to run.
        """
        assert payload.gwas_data is not None
        assert payload.ld_snps_bim_df is not None
        assert payload.ld_matrix is not None

        if payload.coloc2 and coloc2eqtl_df is None:
            raise InvalidUsage("COLOC2 was selected but no coloc2eqtl_df was provided")

        if len(payload.gwas_data_kept) == 0:
            raise InvalidUsage("No SNPs in Simple Sum subset")

        SS_positions = payload.gwas_data_kept["POS"].to_list()

        # Handle LD matrix
        ld_matrix = payload.ld_matrix
        ld_mat_snps = payload.ld_snps_bim_df["CHROM_POS"].to_list()
        ld_mat_positions = payload.ld_snps_bim_df["POS"].to_list()

        np.fill_diagonal(ld_matrix, np.diag(ld_matrix) + LD_MAT_DIAG_CONSTANT)

        p_matrix_indices = [
            i for i, e in enumerate(SS_positions) if e in ld_mat_positions
        ]
        p_value_matrix = p_value_matrix[:, p_matrix_indices]  # type: ignore

        write_matrix(p_value_matrix, payload.file.p_value_filepath)
        write_matrix(ld_matrix, payload.file.ld_matrix_filepath)
        # Extra files written for LD matrix:
        write_list(ld_mat_snps, payload.file.ld_mat_snps_filepath)
        write_list(ld_mat_positions, payload.file.ld_mat_positions_filepath)

        # Run Simple Sum
        try:
            SSdf = simple_sum(
                payload.file.p_value_filepath,
                payload.file.ld_matrix_filepath,
                payload.file.simple_sum_results_filepath,
                payload.get_p_value_threshold(),
            )
        except ScriptError as e:
            raise InvalidUsage(e.message, status_code=410)

        # Save results
        payload.ss_result_df = SSdf

        # Gather results
        SSPvalues = SSdf["Pss"].tolist()
        num_SNP_used_for_SS = SSdf["n"].tolist()
        comp_used = SSdf["comp_used"].tolist()
        first_stages = SSdf["first_stages"].tolist()
        first_stage_p = SSdf["first_stage_p"].tolist()

        for i in np.arange(len(SSPvalues)):
            if SSPvalues[i] > 0:
                SSPvalues[i] = np.format_float_scientific((-np.log10(SSPvalues[i])), precision=2)  # type: ignore
        SSPvaluesMatGTEx = np.empty(0)
        num_SNP_used_for_SSMat = np.empty(0)
        comp_usedMat = np.empty(0)
        SSPvaluesSecondary = []
        numSNPsSSPSecondary = []
        compUsedSecondary = []

        gtex_tissues, gtex_genes = payload.get_gtex_selection()

        table_titles = []
        if payload.secondary_datasets is not None:
            table_titles = list(payload.secondary_datasets.keys())

        if len(gtex_tissues) > 0:
            SSPvaluesMatGTEx = np.array(
                SSPvalues[0 : (len(gtex_tissues) * len(gtex_genes))]
            ).reshape(len(gtex_tissues), len(gtex_genes))
            num_SNP_used_for_SSMat = np.array(
                num_SNP_used_for_SS[0 : (len(gtex_tissues) * len(gtex_genes))]
            ).reshape(len(gtex_tissues), len(gtex_genes))
            comp_usedMat = np.array(
                comp_used[0 : (len(gtex_tissues) * len(gtex_genes))]
            ).reshape(len(gtex_tissues), len(gtex_genes))
        if len(SSPvalues) > len(gtex_tissues) * len(gtex_genes):
            SSPvaluesSecondary = SSPvalues[
                (len(gtex_tissues) * len(gtex_genes)) : (len(SSPvalues))
            ]
            numSNPsSSPSecondary = num_SNP_used_for_SS[
                (len(gtex_tissues) * len(gtex_genes)) : (len(SSPvalues))
            ]
            compUsedSecondary = comp_used[
                (len(gtex_tissues) * len(gtex_genes)) : (len(SSPvalues))
            ]
        SSPvalues_dict = {
            "Genes": gtex_genes,
            "Tissues": gtex_tissues,
            "Secondary_dataset_titles": table_titles,
            "SSPvalues": SSPvaluesMatGTEx.tolist()  # GTEx pvalues
            # ,'Num_SNPs_Used_for_SS': [int(x) for x in num_SNP_used_for_SS]
            ,
            "Num_SNPs_Used_for_SS": num_SNP_used_for_SSMat.tolist(),
            "Computation_method": comp_usedMat.tolist(),
            "SSPvalues_secondary": SSPvaluesSecondary,
            "Num_SNPs_Used_for_SS_secondary": numSNPsSSPSecondary,
            "Computation_method_secondary": compUsedSecondary,
            "First_stages": first_stages,
            "First_stage_Pvalues": first_stage_p,
        }
        json.dump(SSPvalues_dict, open(payload.file.SSPvalues_filepath, "w"))

        ####################################################################################################

        if payload.coloc2:

            assert coloc2eqtl_df is not None

            # COLOC2
            gwas_to_coloc2_colnames = {
                "CHROM": "CHR",
                "SNP": "SNPID",
                "P": "PVAL",
                "REF": "A2",
                "ALT": "A1",
            }

            coloc2_gwasdf = payload.gwas_data_kept.rename(
                columns=gwas_to_coloc2_colnames
            ).reindex(columns=self.COLOC2_GWAS_COLNAMES)

            # Calculating COLOC2 stats
            coloc2_gwasdf.dropna().to_csv(
                payload.file.coloc2_gwas_filepath,
                index=False,
                encoding="utf-8",
                sep="\t",
            )
            if coloc2_gwasdf.shape[0] == 0 or coloc2eqtl_df.shape[0] == 0:
                raise InvalidUsage(
                    f"Empty datasets for coloc2. Cannot proceed. GWAS numRows: {coloc2_gwasdf.shape[0]}; eQTL numRows: {coloc2eqtl_df.shape[0]}. May be due to inability to match with GTEx variants. Please check position, REF/ALT allele correctness, and or SNP names."
                )
            coloc2eqtl_df.dropna().to_csv(
                payload.file.coloc2_eqtl_filepath,
                index=False,
                encoding="utf-8",
                sep="\t",
            )

            # Run COLOC2
            try:
                coloc2df = coloc2(
                    payload.file.coloc2_gwas_filepath,
                    payload.file.coloc2_eqtl_filepath,
                    payload.file.coloc2_results_filepath,
                )
            except ScriptError as e:
                raise InvalidUsage(e.stdout, status_code=410)

            # save as json:
            coloc2_dict = {
                "ProbeID": coloc2df["ProbeID"].tolist(),
                "PPH4abf": coloc2df["PPH4abf"].tolist(),
            }
        else:
            coloc2_dict = {"ProbeID": [], "PPH4abf": []}

        json.dump(coloc2_dict, open(payload.file.coloc2_filepath, "w"))

        return
