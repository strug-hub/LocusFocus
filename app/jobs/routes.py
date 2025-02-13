from flask import Blueprint, jsonify, current_app as app, url_for
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
        
        # message
        try:
            error_message = error.message["message"]
        except TypeError:
            error_message = error.message

        if isinstance(error, InvalidUsage):
            return jsonify(
                {
                    "status": "FAILURE",
                    "error_title": "Invalid Usage Error",
                    "error_message": error_message,
                    "status_code": error.status_code,
                    "payload": error.payload,
                }
            )
        elif isinstance(error, ServerError):
            return jsonify(
                {
                    "status": "FAILURE",
                    "error_title": "Server Error",
                    "error_message": error_message,
                    "status_code": error.status_code,
                    "payload": error.payload,
                }
            )
        else:
            # Unexpected error
            app.logger.error("An unexpected error occurred!")
            return jsonify(
                {
                    "status": "FAILURE",
                    "error_title": "Unexpected Server Error",
                    "error_message": "An unexpected error occurred. Please contact your system administrator.",
                    "status_code": 500,
                }
            )

    if result.status == "SUCCESS":
        return jsonify(
            {
                "status": "SUCCESS",
                "redirect_url": url_for("prev_session_input", old_session_id=result.id),
            }
        )

    # default (in-progress)
    return jsonify({
        "status": result.status,
        "stage_index": result.info.get("stage_index"),
        "stage_count": result.info.get("stage_count"),
    })
