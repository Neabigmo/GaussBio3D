"""
Protein–ligand interaction preset pipeline
蛋白–配体互作预设流水线
"""

from __future__ import annotations

import numpy as np
from typing import Optional

from ..config import MgliConfig
from ..core.pipeline import MGLIPipeline, PCAProjector
from ..molecules.protein import Protein
from ..molecules.ligand import Ligand


def default_dti_mgli_pipeline(n_components: Optional[int] = None) -> MGLIPipeline:
    cfg = MgliConfig(
        distance_bins=[3.0, 5.0, 7.0, 10.0, 15.0],
        use_rbf=False,
        signed=False,
        stats=["sum", "median", "std"],
        group_mode_A="element",
        group_mode_B="element",
        max_distance=15.0,
    )
    projector = PCAProjector(n_components=n_components or 256, whiten=True)
    return MGLIPipeline(config=cfg, projector=projector)


def compute_pl_complex_mgli(protein_pdb: str, ligand_sdf: str, chain_id: Optional[str] = None) -> np.ndarray:
    A = Protein.from_pdb(protein_pdb, chain_id=chain_id)
    B = Ligand.from_sdf(ligand_sdf)
    pipe = default_dti_mgli_pipeline()
    X = pipe.fit_transform([(A, B)])
    return X[0]

