from flask import Request

from app.colocalization.payload import SessionPayload
from app.pipeline import PipelineStage


class CreateSessionStage(PipelineStage):
    """
    Given a Flask request, 
    create a Colocalization payload to use for the rest of the pipeline.
    """

    def invoke(self, request: Request) -> SessionPayload:
        payload = SessionPayload(request)
        return payload