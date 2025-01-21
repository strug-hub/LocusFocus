"""
Celery tasks for LocusFocus.
"""

from typing import Iterable
from uuid import UUID

from celery import shared_task
from celery.result import AsyncResult
from flask import current_app as app


def run_pipeline_async(
    session_id: UUID, pipeline_type: str, payload_params: Iterable
) -> AsyncResult:
    """
    Run a pipeline asynchronously using Celery.
    Return the task ID, which can be used to check the task status.

    pipeline must be one of: "colocalization"
    """
    result = _pipeline_task.apply_async(session_id=session_id, pipeline_type=pipeline_type, payload_params=payload_params, task_id=str(session_id))  # type: ignore
    return result


@shared_task(ignore_result=False)
def _pipeline_task(**kwargs) -> object:
    """
    Celery task for running a pipeline.

    Do not call this function directly. Instead, use `run_pipeline_async`.
    """
    session_id = kwargs.pop("session_id")
    pipeline_type = kwargs.pop("pipeline_type")
    payload_params = kwargs.pop("payload_params")

    pipeline = None
    if pipeline_type == "colocalization":
        from app.colocalization.pipeline import ColocalizationPipeline

        pipeline = ColocalizationPipeline(id=session_id)
    else:
        raise ValueError(f"Invalid pipeline type: '{pipeline_type}'")

    app.logger.debug(f"Running pipeline '{session_id}' with type '{pipeline_type}'")

    result = pipeline.process(*payload_params)

    return result
