import json
from app.colocalization.payload import SessionPayload
from app.colocalization.utils import get_session_filepath
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

        # Save data in JSON format for plotting
        sessionfilepath = get_session_filepath(f"form_data-{payload.session_id}.json")

        json.dump(payload.dump_session_data(), open(sessionfilepath, 'w'))
        return payload
