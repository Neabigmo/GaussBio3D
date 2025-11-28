"""
MGLI Pipeline
MGLI流水线

Sklearn-style pipeline to compute multiscale mGLI descriptors with optional projection.
仿sklearn风格的流水线，用于计算多尺度mGLI描述符并可选降维。
"""

from __future__ import annotations

import numpy as np
from typing import Optional, Sequence

from ..config import MgliConfig
from ..core.geometry import Structure
from ..features.descriptor import global_mgli_descriptor


class Projector:
    def fit(self, X: np.ndarray) -> "Projector":
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        return X


class PCAProjector(Projector):
    def __init__(self, n_components: int = 256, whiten: bool = True):
        self.n_components = int(n_components)
        self.whiten = bool(whiten)
        self._model = None

    def fit(self, X: np.ndarray) -> "PCAProjector":
        try:
            from sklearn.decomposition import PCA  # type: ignore
            self._model = PCA(n_components=self.n_components, whiten=self.whiten)
            self._model.fit(X)
        except Exception:
            self._model = None
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self._model is None:
            return X
        return self._model.transform(X)


class MGLIPipeline:
    def __init__(
        self,
        config: Optional[MgliConfig] = None,
        projector: Optional[Projector] = None,
    ):
        self.config = config or MgliConfig()
        self.projector = projector

    def _featurize_one(self, A: Structure, B: Optional[Structure]) -> np.ndarray:
        return global_mgli_descriptor(A, B, self.config)

    def fit(self, pairs: Sequence[tuple[Structure, Optional[Structure]]]) -> "MGLIPipeline":
        X = [self._featurize_one(a, b) for a, b in pairs]
        X = np.vstack(X) if len(X) else np.zeros((0, 0), dtype=float)
        if self.projector is not None and X.size:
            self.projector.fit(X)
        return self

    def transform(self, pairs: Sequence[tuple[Structure, Optional[Structure]]]) -> np.ndarray:
        X = [self._featurize_one(a, b) for a, b in pairs]
        X = np.vstack(X) if len(X) else np.zeros((0, 0), dtype=float)
        if self.projector is not None and X.size:
            X = self.projector.transform(X)
        return X

    def fit_transform(self, pairs: Sequence[tuple[Structure, Optional[Structure]]]) -> np.ndarray:
        self.fit(pairs)
        return self.transform(pairs)


__all__ = [
    "Projector",
    "PCAProjector",
    "MGLIPipeline",
]

