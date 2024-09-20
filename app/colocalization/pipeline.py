from app.pipeline import Pipeline
from app.pipeline.pipeline_stage import PipelineStage


class ColocalizationPipeline(Pipeline):
    """
    Pipeline class for handling colocalization.
    """

    def __init__(self):
        super().__init__()
        # TODO: add stages

    def invoke_stage(self, stage: PipelineStage, payload: object) -> object:
        return super().invoke_stage(stage, payload)
