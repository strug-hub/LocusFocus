from collections import namedtuple
from datetime import datetime
import json

from app.colocalization.payload import SessionPayload
from app.pipeline import PipelineStage
from app.utils.errors import InvalidUsage


class CreateSessionStage(PipelineStage):
    """
    Given an initial payload containing a Flask request's form data and file data,
    create a Colocalization payload to use for the rest of the pipeline.
    """

    def name(self) -> str:
        return "create-session"

    def invoke(self, payload: SessionPayload) -> SessionPayload:

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
                "files_uploaded": payload.uploaded_files,
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
        if payload.uploaded_files is None or len(payload.uploaded_files) == 0:
            raise InvalidUsage(f"No files found in request")
        return None
