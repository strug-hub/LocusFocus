from typing import List # type: ignore
from .pipeline_stage import PipelineStage


class Pipeline():
    """
    Generic pipeline class for creating an ordered series of stages to execute on a given payload.

    Similar implementation as the Chain of Responsibility design pattern.
    """

    def __init__(self):
        self.stages: List[PipelineStage] = []

    def invoke_stage(self, stage: PipelineStage, payload: object) -> object:
        payload = self.pre_stage(stage, payload)
        payload = stage.invoke(payload)
        payload = self.post_stage(stage, payload)
        return payload

    def pre_stage(self, stage: PipelineStage, payload: object):
        """
        Operations performed before the stage's `invoke` method is called.
        Can be extended.
        """
        return payload

    def post_stage(self, stage: PipelineStage, payload: object):
        """
        Operations performed after the stage's `invoke` method is called.
        Can be extended.
        """
        return payload

    def pipe(self, *stages: PipelineStage):
        """Add new stage(s) to the pipeline. 
        Serves as a builder function for the pipeline, and can be chain-called.

        Args:
            stage (PipelineStage): A pipeline stage to add to the pipeline.

        Returns:
            Pipeline: The new pipeline with the new stage added.
        """
        self.stages.extend(stages)
        return self

    def process(self, payload: object) -> object:
        """Execute the pipeline stages in-order on the given payload.

        Args:
            payload (object): A payload to be processed by this pipeline.

        Returns:
            object: A payload that has been processed by all stages in this pipeline.
        """
        for stage in self.stages:
            payload = self.invoke_stage(stage, payload)
        return payload
