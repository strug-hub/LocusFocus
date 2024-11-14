from datetime import datetime
import json
from flask import Request
from werkzeug.utils import secure_filename

from app.colocalization.payload import SessionPayload
from app.pipeline import PipelineStage
from app.utils.errors import InvalidUsage


class CreateSessionStage(PipelineStage):
    """
    Given a Flask request,
    create a Colocalization payload to use for the rest of the pipeline.
    """

    def name(self) -> str:
        return "create-session"

    def invoke(self, request: Request) -> SessionPayload:
        payload = SessionPayload(request=request)

        self._create_metadata_file(payload)
        self._check_file_upload(payload)

        return payload

    def _create_metadata_file(self, payload: SessionPayload):
        """
        Create JSON dict of session data needed for metadata file.

        The existence of the metadata file is what we use to check whether a session has been started.
        """

        metadata = {}
        metadata.update(
            {
                "datetime": datetime.now().isoformat(),
                "files_uploaded": [
                    file.filename or ""
                    for file in payload.request.files.getlist("files[]")
                ],
                "session_id": str(payload.session_id),
                "type": "default",
            }
        )

        json.dump(metadata, open(payload.file.metadata_filepath, "w"))
        return None

    def _check_file_upload(self, payload: SessionPayload):
        """
        Check if the user has uploaded any files.
        """
        if "files[]" not in payload.request.files:
            raise InvalidUsage(f"No files found in request")
        return None
