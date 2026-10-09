"""Exact output-state metrics."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


def normalized_state(state: ComplexArray) -> ComplexArray:
    """Normalize a finite nonzero state without overflow or underflow."""

    array = np.asarray(state, dtype=np.complex128)
    if array.size == 0:
        raise ValueError("state must be nonempty")
    scale = float(np.max(np.abs(array)))
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("state must be finite and nonzero")
    scaled = array / scale
    return scaled / np.linalg.norm(scaled)


def clean_joint_state(probabilities: FloatArray) -> ComplexArray:
    amplitudes = np.sqrt(np.asarray(probabilities, dtype=np.float64))
    n_states = len(amplitudes)
    state = np.zeros((2, n_states, n_states), dtype=np.complex128)
    state[0, :, 0] = amplitudes
    return state


def joint_trace_distance(state: ComplexArray, probabilities: FloatArray) -> float:
    target = clean_joint_state(probabilities)
    normalized = normalized_state(state)
    overlap = np.vdot(target, normalized)
    orthogonal_residual = normalized - overlap * target
    return min(1.0, float(np.linalg.norm(orthogonal_residual)))


def reduced_density_x(state: ComplexArray) -> ComplexArray:
    amplitudes = normalized_state(state).transpose(1, 0, 2)
    matrix = amplitudes.reshape(amplitudes.shape[0], -1)
    return np.einsum("xe,ye->xy", matrix, matrix.conj(), optimize=True)


def reduced_trace_distance(state: ComplexArray, probabilities: FloatArray) -> float:
    amplitudes = np.sqrt(np.asarray(probabilities, dtype=np.float64))
    target = np.outer(amplitudes, amplitudes)
    difference = reduced_density_x(state) - target
    difference = (difference + difference.conj().T) / 2.0
    eigenvalues = np.linalg.eigvalsh(difference)
    return float(0.5 * np.sum(np.abs(eigenvalues)))


def measured_distribution_x(state: ComplexArray) -> FloatArray:
    normalized = normalized_state(state)
    return np.sum(np.abs(normalized) ** 2, axis=(0, 2)).real


def total_variation_distance(state: ComplexArray, probabilities: FloatArray) -> float:
    observed = measured_distribution_x(state)
    return float(0.5 * np.sum(np.abs(observed - probabilities)))


def clean_subspace_leakage(state: ComplexArray) -> float:
    normalized = normalized_state(state)
    if normalized.ndim != 3 or normalized.shape[0] != 2:
        raise ValueError("state must have axes (signal, A, B) with two signal levels")
    signal_one = float(np.sum(np.abs(normalized[1]) ** 2))
    dirty_transition = float(np.sum(np.abs(normalized[0, :, 1:]) ** 2))
    return min(1.0, signal_one + dirty_transition)


def norm_residual(state: ComplexArray) -> float:
    return abs(float(np.vdot(state, state).real) - 1.0)
