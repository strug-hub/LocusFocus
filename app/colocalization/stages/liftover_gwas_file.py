from typing import List

import numpy as np
import pandas as pd

from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.liftover import run_liftover


class LiftoverGWASFile(PipelineStage):
    """Lift over coordinates as needed and update `gwas_indices_kept` to exclude non-lifted-over rows"""

    def name(self) -> str:
        return "liftover-gwas-file"

    def invoke(self, payload: SessionPayload) -> object:
        # These two will be from form fields, I think
        needs_liftover = False
        liftover_target = "hg19"

        # keep type checker happy
        if payload.gwas_data is not None:

            if needs_liftover:

                lifted_over, unlifted_over = run_liftover(
                    payload.gwas_data, liftover_target
                )

                payload.gwas_data = lifted_over

                payload.unlifted_over_indices = pd.Series(list(unlifted_over))

                payload.gwas_indices_kept = self.update_indices_kept(
                    payload, unlifted_over
                )

        return payload

    def update_indices_kept(
        self,
        payload: SessionPayload,
        unlifted_over: List[int],
    ) -> pd.Series:

        lifted_over_kept_ = pd.Series(
            np.ones(len(payload.gwas_data_original), dtype=bool)
        )

        lifted_over_kept_[unlifted_over] = False

        return pd.Series(np.logical_and(lifted_over_kept_, payload.gwas_indices_kept))
