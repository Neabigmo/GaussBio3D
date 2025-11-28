"""
Preset pipelines for common tasks
常见任务的预设流水线
"""

from .protein import flexibility_mgli_pipeline, compute_bfactor_mgli
from .dti import default_dti_mgli_pipeline, compute_pl_complex_mgli

__all__ = [
    "flexibility_mgli_pipeline",
    "compute_bfactor_mgli",
    "default_dti_mgli_pipeline",
    "compute_pl_complex_mgli",
]

