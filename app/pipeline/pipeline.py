from typing import List # type: ignore
from .pipeline_stage import PipelineStage
from .payload import Payload


class Pipeline():
    """
    Generic pipeline class for creating an ordered series of stages to execute on a given payload.

    Similar implementation as the Chain of Responsibility design pattern.
    """

    def __init__(self):
        self.stages: List[PipelineStage] = []

    def invoke_stage(self, stage: PipelineStage, payload: Payload) -> Payload:
        # pre-stage steps go here
        new_payload = stage.invoke(payload)
        # post-stage steps go here
        return new_payload

    def pipe(self, stage: PipelineStage):
        """Add a new stage to the pipeline. 
        Serves as a builder function for the pipeline, and can be chain-called.

        Args:
            stage (PipelineStage): A pipeline stage to add to the pipeline.

        Returns:
            Pipeline: The new pipeline with the new stage added.
        """
        self.stages.append(stage)
        return self

    def process(self, payload: Payload) -> Payload:
        """Execute the pipeline stages in-order on the given payload.

        Args:
            payload (object): A payload to be processed by this pipeline.

        Returns:
            object: A payload that has been processed by all stages in this pipeline.
        """
        for stage in self.stages:
            payload = self.invoke_stage(stage, payload)
        return payload
