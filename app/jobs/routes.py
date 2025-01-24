from flask import Blueprint, jsonify, current_app as app
from celery.result import AsyncResult


jobs_bp = Blueprint("jobs", __name__)


@jobs_bp.route("/job/status/<job_id>", methods=["GET"])
def get_job_status(job_id):
    """
    Route for getting or updating the status of a job.

    Args:
        job_id (str): The ID of the job to get or update the status of.
    """
    result = AsyncResult(job_id, app=app.extensions["celery"])

    return jsonify({"status": result.status})
