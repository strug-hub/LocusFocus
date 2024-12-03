"""
Celery tasks for LocusFocus.
"""

from uuid import UUID
from celery import shared_task

from app.pipeline import Pipeline

def run_pipeline_async(pipeline: Pipeline, initial_payload: object) -> UUID:
    """
    Run a pipeline asynchronously using Celery.
    Return the task ID, which can be used to check the task status.
    """
    result = _pipeline_task.apply_async(pipeline, initial_payload, task_id=pipeline.id)  # type: ignore
    return result.id


@shared_task(ignore_result=False)
def _pipeline_task(pipeline: Pipeline, initial_payload: object) -> object:
    """
    Celery task for running a pipeline.
    
    Do not call this function directly. Instead, use `run_pipeline_async`.
    """
    return pipeline.process(initial_payload)
