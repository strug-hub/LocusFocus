import pandas as pd

from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.helpers import adjust_snp_column
from app.utils.liftover import run_liftover
from app.utils.errors import InvalidUsage


class LiftoverSecondaryDatasets(PipelineStage):
    """Lift over coordinates as needed for secondary datasets"""

    def name(self) -> str:
        return "liftover-secondary-datasets"

    def description(self) -> str:
        return "LiftOver secondary datasets (if needed, may be skipped)"

    def invoke(self, payload: SessionPayload) -> object:
        if payload.secondary_datasets is None:
            return payload

        needs_liftover = False

        if payload.get_gtex_version() == "V7":
            raise InvalidUsage(
                "GTEx V7 is no longer available. Please use GTEx V8 or GTEx V10."
            )
        else:
            liftover_target = "hg38"

        needs_liftover = payload.get_secondary_coordinate() != liftover_target

        if needs_liftover:
            payload.secondary_datasets_unlifted_indices = {}

            for table_title, table in payload.secondary_datasets.items():
                dataset = pd.DataFrame(table)

                if dataset.shape[0] == 0:
                    continue

                lifted_over, unlifted_over = run_liftover(
                    dataset, liftover_target, chrom_col="CHROM", pos_col="BP"
                )

                lifted_over = adjust_snp_column(
                    lifted_over,
                    liftover_target,
                    ignore_alleles=True,
                )

                payload.secondary_datasets[table_title] = lifted_over.fillna(
                    -1
                ).to_dict(  # type: ignore
                    orient="records"
                )

                payload.secondary_datasets_unlifted_indices[table_title] = unlifted_over

        return payload
