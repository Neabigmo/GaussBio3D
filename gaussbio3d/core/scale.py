"""
Multi-scale weighting schemes for mGLI
多尺度加权方案（mGLI）

Provides binning and RBF-based radial weighting on pairwise GLI matrices.
为成对GLI矩阵提供分箱与RBF径向加权。
"""

from __future__ import annotations

import numpy as np
from typing import Sequence


class ScaleScheme:
    def apply(self, G: np.ndarray, dists: np.ndarray) -> np.ndarray:
        raise NotImplementedError


class BinningScaleScheme(ScaleScheme):
    def __init__(self, edges: Sequence[float]):
        self.edges = np.asarray(edges, dtype=float)
        assert self.edges.ndim == 1 and self.edges.size >= 2

    def apply(self, G: np.ndarray, dists: np.ndarray) -> np.ndarray:
        n, m = G.shape
        K = self.edges.size - 1
        out = np.zeros((K, n, m), dtype=float)
        for k in range(K):
            mask = (dists >= self.edges[k]) & (dists < self.edges[k + 1])
            if np.any(mask):
                out[k][mask] = G[mask]
        return out


class RBFScaleScheme(ScaleScheme):
    def __init__(self, centers: Sequence[float], sigma: float | None = None):
        self.centers = np.asarray(centers, dtype=float)
        assert self.centers.ndim == 1 and self.centers.size >= 1
        if sigma is None:
            if self.centers.size == 1:
                sigma = 1.0
            else:
                gaps = np.diff(np.sort(self.centers))
                sigma = float(np.mean(gaps)) or 1.0
        self.sigma = float(sigma)

    def apply(self, G: np.ndarray, dists: np.ndarray) -> np.ndarray:
        n, m = G.shape
        K = self.centers.size
        out = np.zeros((K, n, m), dtype=float)
        r = dists[None, :, :]
        c = self.centers[:, None, None]
        W = np.exp(-((r - c) ** 2) / (2.0 * self.sigma**2))
        for k in range(K):
            out[k] = G * W[k]
        return out


__all__ = [
    "ScaleScheme",
    "BinningScaleScheme",
    "RBFScaleScheme",
]

