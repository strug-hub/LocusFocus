from typing import List

import numpy as np
import pandas as pd
from flask import current_app

from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.helpers import adjust_snp_column
from app.utils.liftover import run_liftover


class LiftoverGWASFile(PipelineStage):
    """
    Lift over coordinates as needed and update `gwas_indices_kept` to exclude non-lifted-over rows
    """

    def name(self) -> str:
        return "liftover-gwas-file"

    def description(self) -> str:
        return "LiftOver GWAS file (if needed, may be skipped)"

    def invoke(self, payload: SessionPayload) -> object:
        needs_liftover = False

        if payload.get_gtex_version() == "V7":
            liftover_target = "hg19"
        else:
            liftover_target = "hg38"

        needs_liftover = payload.get_coordinate() != liftover_target

        # keep type checker happy
        if payload.gwas_data is not None:

            if needs_liftover:
                current_app.logger.debug("Lifting over GWAS data")

                lifted_over, unlifted_over = run_liftover(
                    payload.gwas_data, liftover_target, chrom_col="CHROM", pos_col="POS"
                )
                current_app.logger.debug(f"Original: {len(payload.gwas_data)} -> Lifted over: {len(lifted_over)}, unlifted over: {len(unlifted_over)}")

                # Handle lead SNP
                lead_snp_is_user_defined = payload.get_lead_snp_name() != ""
                lead_snp_is_rsid = payload.gwas_data.iloc[payload.get_lead_snp_index()]["SNP"].startswith("rs")
                lead_snp_lifted_over = payload.get_lead_snp_index() not in unlifted_over

                if lead_snp_is_user_defined:
                    if not lead_snp_lifted_over:
                        # clear the lead SNP name, warn the user
                        old_lead_snp = payload.get_lead_snp_name()
                        new_lead_snp = lifted_over[lifted_over["P"].argmin()]["SNP"]  # lowest p-value
                        payload.liftover_lead_snp_warning = f"The specified lead SNP '{old_lead_snp}' was not lifted over. The lifted SNP with the lowest p-value is used instead: '{new_lead_snp}'."
                        payload.lead_snp_name = ""

                lifted_over = adjust_snp_column(lifted_over, liftover_target)

                # lifted_over and gwas_data are not the same size so we need to be careful
                changed_indices = payload.gwas_data.index[~payload.gwas_data.index.isin(unlifted_over)]
                payload.gwas_data.loc[changed_indices] = lifted_over.values  # type: ignore

                if lead_snp_is_user_defined and lead_snp_lifted_over and not lead_snp_is_rsid:
                    # update the lead SNP name if it changed due to liftover
                    assert payload.gwas_data_original is not None  # type checker
                    lead_snp_index = payload.gwas_data.index[payload.gwas_data_original["SNP"] == payload.get_lead_snp_name()].tolist()[0]
                    payload.lead_snp_name = str(payload.gwas_data.iloc[lead_snp_index]["SNP"])

                payload.gwas_indices_kept = self.update_indices_kept(
                    payload, unlifted_over
                )

            else:
                current_app.logger.debug("No liftover needed")

        return payload

    def update_indices_kept(
        self,
        payload: SessionPayload,
        unlifted_over: List[int],
    ) -> pd.Series:
        """
        Helper function to update `gwas_indices_kept` to exclude non-lifted-over rows
        """
        assert payload.gwas_data_original is not None

        lifted_over_kept_ = pd.Series(
            np.ones(len(payload.gwas_data_original), dtype=bool)
        )

        lifted_over_kept_[unlifted_over] = False

        return pd.Series(np.logical_and(lifted_over_kept_, payload.gwas_indices_kept))
