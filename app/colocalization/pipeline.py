from os import PathLike
import timeit
from typing import List

from flask import current_app as app
from werkzeug.datastructures import ImmutableMultiDict

from app.colocalization.payload import SessionPayload
from app.utils import get_session_filepath
from app.pipeline import Pipeline
from app.pipeline.pipeline_stage import PipelineStage
from app.colocalization import stages
from app.utils.errors import LocusFocusError


class ColocalizationPipeline(Pipeline):
    """
    Pipeline class for handling colocalization.
    """

    def __init__(self, id=None, bound_task=None):
        super().__init__(id=id, bound_task=bound_task)
        self.pipe(
            stages.CreateSessionStage(),
            stages.CollectUserInputStage(),
            stages.ReadGWASFileStage(enforce_one_chrom=False),
            stages.LiftoverGWASFile(),
            stages.ReadSecondaryDatasetsStage(),
            stages.LiftoverSecondaryDatasets(),
            stages.ReportGTExDataStage(),
            stages.SimpleSumSubsetGWASStage(),
            stages.GetLDMatrixStage(),
            stages.ColocSimpleSumStage(),
            stages.FinalizeResultsStage(),
        )
        self.timers = {f"{stage.name()}": 0.0 for stage in self.stages}

    def process(
        self, request_form: ImmutableMultiDict, uploaded_files: List[PathLike]
    ) -> SessionPayload:
        """
        Run the colocalization pipeline with the provided Request form and file upload dicts.

        Args:
            request_form: The form data from the request.
            uploaded_files: The uploaded file data from the request.

        Returns:
            SessionPayload: The final payload object that has been processed by all stages in this pipeline.
        """
        initial_payload = SessionPayload(
            request_form=request_form,
            uploaded_files=uploaded_files,
            session_id=self.id,
        )

        return super().process(initial_payload)  # type: ignore

    def pre_stage(self, stage: PipelineStage, payload: object):
        # Timer for each stage (start time)
        app.logger.debug(f"Starting stage {stage.name()}")
        if stage.name() in self.timers:
            self.timers[stage.name()] = timeit.default_timer()

        return super().pre_stage(stage, payload)

    def post_stage(self, stage: PipelineStage, payload: object):
        # Timer for each stage (stop time - start time)
        app.logger.debug(f"Finished stage {stage.name()}")
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
                percentage = round(timer / self.timers["pipeline"] * 100, 2)

                f.write(f"'{stage_name}': {timer} ({percentage}%)\n")
            f.write(f"Total time: {self.timers['pipeline']}\n")

        return super().post_pipeline(payload)

    def invoke_stage(self, stage: PipelineStage, payload: object) -> object:
        try:
            return super().invoke_stage(stage, payload)
        except LocusFocusError as e:
            # Expected errors
            payload = self.post_pipeline(payload)  # type: ignore
            e.message = f"[{stage.name()}] {e.message}"
            raise e
        except BaseException as e:
            # Unexpected errors
            raise e
