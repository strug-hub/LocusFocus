import timeit

from flask import Request

from app.colocalization.payload import SessionPayload
from app.utils import get_session_filepath
from app.pipeline import Pipeline
from app.pipeline.pipeline_stage import PipelineStage
from app.colocalization.stages import *
from app.utils.errors import LocusFocusError


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
            FinalizeResultsStage(),
        )

        self.timers = {f"{stage.name()}": 0.0 for stage in self.stages}

    def process(self, payload: Request) -> SessionPayload:
        return super().process(payload) # type: ignore

    def pre_stage(self, stage: PipelineStage, payload: object):
        # Timer for each stage (start time)
        if stage.name() in self.timers:
            self.timers[stage.name()] = timeit.default_timer()

        return super().pre_stage(stage, payload)

    def post_stage(self, stage: PipelineStage, payload: object):
        # Timer for each stage (stop time - start time)
        if stage.name() in self.timers:
            self.timers[stage.name()] = (
                timeit.default_timer() - self.timers[stage.name()]
            )

        return super().post_stage(stage, payload)

    def pre_pipeline(self, payload: object):
        # Timer for the entire pipeline
        self.timers["pipeline"] = timeit.default_timer()

        return super().pre_pipeline(payload)

    def post_pipeline(self, payload: SessionPayload):
        # Timer for the entire pipeline (stop time - start time)
        self.timers["pipeline"] = timeit.default_timer() - self.timers["pipeline"]

        # Save timers to file
        timing_file_path = get_session_filepath(f"times-{payload.session_id}.txt")
        with open(timing_file_path, "w") as f:
            f.write("-----------------------------------------------------------\n")
            f.write(" Times Report\n")
            f.write("-----------------------------------------------------------\n")
            for stage_name, timer in self.timers.items():
                f.write(f"'{stage_name}': {timer}\n")
            f.write(f"Total time: {self.timers['pipeline']}\n")

        return super().post_pipeline(payload)

    def invoke_stage(self, stage: PipelineStage, payload: object):
        try:
            return super().invoke_stage(stage, payload)
        except LocusFocusError as e:
            e.message = f"[{stage.name()}] {e.message}"
            raise e
