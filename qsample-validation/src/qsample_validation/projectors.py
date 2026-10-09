"""Projector-intertwining checks."""

from __future__ import annotations

import numpy as np

from .szegedy import SpectralWalk


def projector_intertwining_probe_residual(walk: SpectralWalk, seed: int = 0) -> float:
    """Return the residual on one deterministic unit probe vector."""

    generator = np.random.default_rng(seed)
    vector = generator.normal(size=walk.n_states)
    vector = vector + 1.0j * generator.normal(size=walk.n_states)
    vector = vector.astype(np.complex128)
    vector /= np.linalg.norm(vector)

    left = walk.projector_one(walk.apply_t(vector))
    stationary = np.sqrt(walk.chain.probabilities).astype(np.complex128)
    right_logical = stationary * np.vdot(stationary, vector)
    right = walk.apply_t(right_logical)
    return float(np.linalg.norm(left - right))


def projector_intertwining_operator_norm(walk: SpectralWalk) -> float:
    """Return the logical-to-walk operator norm in the spectral basis."""

    selected = np.abs(walk.eigenvalues - 1.0) <= 64.0 * np.finfo(float).eps
    eigenvectors = walk.eigenvectors[:, selected].astype(np.complex128)
    numerical_projector = np.einsum(
        "ik,jk->ij", eigenvectors, eigenvectors.conj(), optimize=False
    )
    stationary = np.sqrt(walk.chain.probabilities).astype(np.complex128)
    stationary /= np.linalg.norm(stationary)
    target_projector = np.outer(stationary, stationary.conj())
    return float(np.linalg.norm(numerical_projector - target_projector, ord=2))
