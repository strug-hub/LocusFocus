"""
Utility package for defining abstract classes for Pipelines, PipelineStages,
and other utilities for creating pipelines for LocusFocus.
"""

from app.pipeline.pipeline import Pipeline
from app.pipeline.pipeline_stage import PipelineStage

__all__ = ["Pipeline", "PipelineStage"]
