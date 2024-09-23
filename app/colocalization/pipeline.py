import timeit

from app.pipeline import Pipeline
from app.pipeline.pipeline_stage import PipelineStage
from app.colocalization.stages import *

class ColocalizationPipeline(Pipeline):
    """
    Pipeline class for handling colocalization.
    """

    def __init__(self):
        super().__init__()
        self.pipe(
            CreateSessionStage(),
            CollectUserInputStage(),
            ReadGWASFileStage(enforce_one_chrom=True),
            ReadSecondaryDatasetsStage(),
            ReportGTExDataStage(),
            GetLDMatrixStage(),
            SimpleSumSubsetGWASStage(),
            ColocSimpleSumStage(),
            FinalizeResultsStage()
        )

        self.timers = { f"{stage.name()}": 0.0 for stage in self.stages }

    def pre_stage(self, stage: PipelineStage, payload: object):
        # Timer for each stage
        if stage.name() in self.timers:
            self.timers[stage.name()] = timeit.default_timer()

        return super().pre_stage(stage, payload)
    
    def post_stage(self, stage: PipelineStage, payload: object):
        # Timer for each stage (stop time - start time)
        if stage.name() in self.timers:
            self.timers[stage.name()] = timeit.default_timer() - self.timers[stage.name()]

        return super().post_stage(stage, payload)
