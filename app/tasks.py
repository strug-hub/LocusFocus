"""
Celery tasks for LocusFocus.
"""

from typing import Iterable
from celery import shared_task
from celery.result import AsyncResult

from app.pipeline import Pipeline

def run_pipeline_async(pipeline: Pipeline, payload_params: Iterable) -> AsyncResult:
    """
    Run a pipeline asynchronously using Celery.
    Return the task ID, which can be used to check the task status.
    """
    result = _pipeline_task.apply_async(pipeline, payload_params, task_id=pipeline.id)  # type: ignore
    return result


@shared_task(ignore_result=False)
def _pipeline_task(pipeline: Pipeline, payload_params: Iterable) -> object:
    """
    Celery task for running a pipeline.
    
    Do not call this function directly. Instead, use `run_pipeline_async`.
    """

    result = pipeline.process(*payload_params)

    return result
