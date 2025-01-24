"""
Celery tasks for LocusFocus.
"""
from os import PathLike
from typing import List

from celery import shared_task
from celery.result import AsyncResult
from werkzeug.datastructures import ImmutableMultiDict
from flask import current_app as app


def run_pipeline_async(pipeline_type: str, request_form: ImmutableMultiDict, uploaded_files: List[PathLike]) -> AsyncResult:
    """
    Run a pipeline asynchronously using Celery.
    Return the task ID, which can be used to check the task status.

    pipeline must be one of: "colocalization"
    """
    result = _pipeline_task.apply_async((pipeline_type, request_form, uploaded_files))  # type: ignore
    return result


@shared_task(ignore_result=False, bind=True)
def _pipeline_task(self, pipeline_type: str, request_form: ImmutableMultiDict, uploaded_files: List[PathLike]) -> object:
    """
    Celery task for running a pipeline.

    Do not call this function directly. Instead, use `run_pipeline_async`.
    """
    session_id = self.request.id

    pipeline = None
    if pipeline_type == "colocalization":
        from app.colocalization.pipeline import ColocalizationPipeline

        pipeline = ColocalizationPipeline(id=session_id)
    else:
        raise ValueError(f"Invalid pipeline type: '{pipeline_type}'")

    app.logger.debug(f"Running pipeline '{session_id}' with type '{pipeline_type}'")

    result = pipeline.process(request_form, uploaded_files)

    return result
