import json
from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage


class FinalizeResultsStage(PipelineStage):
    """
    Perform final steps for the colocalization process.
    """

    def name(self) -> str:
        return "finalize-results"

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        # Indicate that the request was a success
        payload.success = True

        json.dump(payload.dump_session_data(), open(payload.file.session_filepath, 'w'))
        return payload
