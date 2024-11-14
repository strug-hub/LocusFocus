"""
Utility package for defining abstract classes for Pipelines, PipelineStages, 
and other utilities for creating pipelines for LocusFocus.
"""

from .pipeline import Pipeline
from .pipeline_stage import PipelineStage

__all__ = ["Pipeline", "PipelineStage"]
