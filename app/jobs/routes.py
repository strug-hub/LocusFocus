from flask import Blueprint, jsonify, current_app as app
from celery.result import AsyncResult

from app.utils.errors import InvalidUsage, ServerError


jobs_bp = Blueprint("jobs", __name__)


@jobs_bp.route("/job/status/<job_id>", methods=["GET"])
def get_job_status(job_id):
    """
    Route for getting or updating the status of a job.

    Args:
        job_id (str): The ID of the job to get or update the status of.
    """
    result = AsyncResult(job_id, app=app.extensions["celery"])

    if result.status == "FAILURE":
        error = result.result
        if isinstance(error, InvalidUsage):
            return jsonify({"status": "FAILURE", "message": error.message["message"], "status_code": error.status_code, "payload": error.payload})
        elif isinstance(error, ServerError):
            return jsonify({"status": "FAILURE", "message": error.message["message"], "status_code": error.status_code, "payload": error.payload})
        else:
            # Unexpected error
            app.logger.error("An unexpected error occurred!")
            app.logger.error(error.__repr__(), exc_info=True)
            return jsonify({"status": "FAILURE", "message": "An unexpected error occurred", "status_code": 500})

    return jsonify({"status": result.status})
