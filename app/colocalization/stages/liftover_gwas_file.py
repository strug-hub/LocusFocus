from typing import List

import numpy as np
import pandas as pd
from flask import current_app as app

from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
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
                app.logger.debug("Lifting over GWAS data")

                lifted_over, unlifted_over = run_liftover(
                    payload.gwas_data, liftover_target, chrom_col="CHROM", pos_col="POS"
                )
                app.logger.debug(f"Original: {len(payload.gwas_data)} -> Lifted over: {len(lifted_over)}, unlifted over: {len(unlifted_over)}")

                payload.gwas_data = lifted_over

                payload.unlifted_over_indices = pd.Series(list(unlifted_over))

                payload.gwas_indices_kept = self.update_indices_kept(
                    payload, unlifted_over
                )

            else:
                app.logger.debug("No liftover needed")

        return payload

    def update_indices_kept(
        self,
        payload: SessionPayload,
        unlifted_over: List[int],
    ) -> pd.Series:
        assert payload.gwas_data_original is not None

        lifted_over_kept_ = pd.Series(
            np.ones(len(payload.gwas_data_original), dtype=bool)
        )

        lifted_over_kept_[unlifted_over] = False

        return pd.Series(np.logical_and(lifted_over_kept_, payload.gwas_indices_kept))
