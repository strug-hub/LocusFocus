import os
from uuid import UUID
from flask import Blueprint, jsonify, current_app as app

from app.utils import get_session_filepath


jobs_bp = Blueprint("jobs", __name__)


@jobs_bp.route("/job/status/<job_id>", methods=["GET"])
def get_job_status(job_id):
    """
    Route for getting or updating the status of a job.

    Args:
        job_id (str): The ID of the job to get or update the status of.
    """
    # We currently check job status by checking for created session files.
    # In the future, we should use a database to store job results.
    # https://github.com/strug-hub/LocusFocus/issues/13

    # Check that job_id is a valid UUID
    try:
        if not str(UUID(job_id, version=4)) == job_id:
            return jsonify({"status": "INVALID"}), 400
    except ValueError:
        return jsonify({"status": "INVALID"}), 400

    status = "PENDING"
    if os.path.exists(get_session_filepath(f"metadata-{job_id}.json")):
        status = "RUNNING"

    if os.path.exists(get_session_filepath(f"form_data-{job_id}.json")):
        status = "COMPLETE"

    return jsonify({"status": status})
