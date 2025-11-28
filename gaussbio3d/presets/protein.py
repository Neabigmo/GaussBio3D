"""
Protein flexibility preset pipeline
蛋白柔性预设流水线
"""

from __future__ import annotations

import numpy as np
from typing import Optional

from ..config import MgliConfig
from ..core.pipeline import MGLIPipeline, PCAProjector
from ..molecules.protein import Protein


def flexibility_mgli_pipeline() -> MGLIPipeline:
    cfg = MgliConfig(
        distance_bins=[5.0 + i for i in range(23)],  # 5–27 Å, step 1 Å
        use_rbf=False,
        signed=False,
        stats=["sum", "median", "std"],
        group_mode_A="element",
        group_mode_B="element",
        max_distance=27.0,
    )
    projector = None
    return MGLIPipeline(config=cfg, projector=projector)


def compute_bfactor_mgli(pdb_path: str, chain_id: Optional[str] = None) -> np.ndarray:
    A = Protein.from_pdb(pdb_path, chain_id=chain_id)
    pipe = flexibility_mgli_pipeline()
    X = pipe.transform([(A, None)])
    return X[0]

