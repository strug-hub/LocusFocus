from typing import List  # type: ignore

from uuid import uuid4

from .pipeline_stage import PipelineStage


class Pipeline:
    """
    Generic pipeline class for creating an 
    ordered series of stages to execute on a given payload.

    Similar implementation as the Chain of Responsibility design pattern.

    Includes methods for handling pre and post-stage operations.
    """

    def __init__(self):
        self.stages: List[PipelineStage] = []
        self.id = uuid4()

    def invoke_stage(self, stage: PipelineStage, payload: object) -> object:
        payload = self.pre_stage(stage, payload)
        payload = stage.invoke(payload)
        payload = self.post_stage(stage, payload)
        return payload

    def pre_stage(self, stage: PipelineStage, payload: object) -> object:
        """
        Operations performed before the stage's `invoke` method is called.
        Executes before every stage in the pipeline, and can be extended.
        
        Parameters:
            stage (PipelineStage): The stage that's about to be invoked.
            payload (object): The payload to be processed by the stage.
        """
        return payload

    def post_stage(self, stage: PipelineStage, payload: object) -> object:
        """
        Operations performed after the stage's `invoke` method is called.
        Executes immediately after every stage in the pipeline, and can be extended.

        Parameters:
            stage (PipelineStage): The stage that was just invoked.
            payload (object): The payload that was processed by the stage.
        """
        return payload

    def pre_pipeline(self, payload: object) -> object:
        """
        Operations performed before the pipeline's `process` method is called.
        Can be extended.
        """
        return payload

    def post_pipeline(self, payload: object) -> object:
        """ 
        Operations performed after the pipeline's `process` method is called.
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
        self.pre_pipeline(payload)
        for stage in self.stages:
            payload = self.invoke_stage(stage, payload)
        self.post_pipeline(payload)
        return payload
