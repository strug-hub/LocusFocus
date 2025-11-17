import json
from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage


class FinalizeResultsStage(PipelineStage):
    """
    Perform final steps for the colocalization process.
    """

    def name(self) -> str:
        return "finalize-results"

    def description(self) -> str:
        return "Finalize results"

    def invoke(self, payload: SessionPayload) -> SessionPayload:
        # Indicate that the request was a success
        payload.success = True

        with open(payload.file.session_filepath, "w", encoding="utf-8") as file:
            json.dump(payload.dump_session_data(), file)

        if payload.report_liftover() != {}:
            with open(payload.file.liftover_filepath, "w", encoding="utf-8") as file:
                json.dump(payload.report_liftover(), file)
        return payload
