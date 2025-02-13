"""
Celery tasks for LocusFocus.
"""
from os import PathLike
from typing import List

from celery import shared_task
from celery.result import AsyncResult
from werkzeug.datastructures import ImmutableMultiDict
from flask import current_app as app
from kombu.exceptions import OperationalError


def get_is_celery_running() -> bool:
    """
    Check if Celery is running and available to accept tasks.
    """
    try:
        return app.extensions["celery"].control.ping() is not None
    except IOError:
        return False
    except OperationalError:
        app.logger.error("Could not connect to redis server. Celery is not available.")
        return False
    except Exception:
        app.logger.error("An unexpected error occurred while checking if Celery is running.")
        return False

def run_pipeline_async(pipeline_type: str, request_form: ImmutableMultiDict, uploaded_filepaths: List[PathLike]) -> AsyncResult:
    """
    Run a pipeline asynchronously using Celery.
    Return the task ID, which can be used to check the task status.

    pipeline must be one of: "colocalization"

    Tasks can be looked up using `AsyncResult(task_id, app=app.extensions["celery"])`.
    Read `.state` to check the task status:
    - PENDING: task does not exist (yet)
    - RUNNING: task has been started by the worker
    - FAILURE: task has failed
    - SUCCESS: task has been executed successfully
    """
    result = _pipeline_task.apply_async((pipeline_type, request_form, uploaded_filepaths))  # type: ignore
    return result


@shared_task(ignore_result=False, bind=True)
def _pipeline_task(self, pipeline_type: str, request_form: ImmutableMultiDict, uploaded_filepaths: List[PathLike]) -> object:
    """
    Celery task for running a pipeline.

    Do not call this function directly. Instead, use `run_pipeline_async`.
    """
    session_id = self.request.id

    pipeline = None
    if pipeline_type == "colocalization":
        from app.colocalization.pipeline import ColocalizationPipeline

        pipeline = ColocalizationPipeline(id=session_id, bound_task=self)
    else:
        raise ValueError(f"Invalid pipeline type: '{pipeline_type}'")

    app.logger.debug(f"Running pipeline '{session_id}' with type '{pipeline_type}'")

    result = pipeline.process(request_form, uploaded_filepaths)

    return result.file.get_plot_template_paths(session_id=str(session_id))
