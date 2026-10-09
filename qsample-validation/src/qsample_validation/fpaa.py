"""Yoder-Low-Chuang fixed-point amplitude amplification."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .metrics import normalized_state

ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True)
class FPAASequence:
    lower_bound: float
    target_error: float
    length: int
    alpha: NDArray[np.float64]
    beta: NDArray[np.float64]

    @property
    def n_phases(self) -> int:
        return self.length - 1

    @property
    def phase_weight(self) -> float:
        """Return the phase-sensitive compiler weight A_i."""

        return float(
            np.sum(np.abs(np.sin(self.alpha / 2.0)))
            + np.sum(np.abs(np.sin(self.beta / 2.0)))
        )


def yoder_width(length: int, target_error: float) -> float:
    scale = math.acosh(1.0 / target_error) / length
    return math.tanh(scale) ** 2


def minimal_odd_length(lower_bound: float, target_error: float) -> int:
    if not 0.0 < lower_bound <= 1.0:
        raise ValueError("lower_bound must lie in (0, 1]")
    if not 0.0 < target_error < 1.0:
        raise ValueError("target_error must lie in (0, 1)")
    length = 3
    while yoder_width(length, target_error) > lower_bound + 1.0e-15:
        length += 2
    return length


def _arccot(value: float) -> float:
    return math.atan2(1.0, value)


def yoder_sequence(lower_bound: float, target_error: float) -> FPAASequence:
    length = minimal_odd_length(lower_bound, target_error)
    count = (length - 1) // 2
    xi = 1.0 / math.cosh(math.acosh(1.0 / target_error) / length)
    root = math.sqrt(max(0.0, 1.0 - xi * xi))
    alpha = np.array(
        [
            2.0 * _arccot(math.tan(2.0 * math.pi * j / length) * root)
            for j in range(1, count + 1)
        ],
        dtype=np.float64,
    )
    beta = -alpha[::-1].copy()
    return FPAASequence(
        lower_bound=float(lower_bound),
        target_error=float(target_error),
        length=length,
        alpha=alpha,
        beta=beta,
    )


def apply_selective_phase(
    state: ComplexArray, target: ComplexArray, phase: float
) -> ComplexArray:
    return state + (np.exp(1.0j * phase) - 1.0) * target * np.vdot(target, state)


def apply_ideal_fpaa(
    source: ComplexArray,
    target: ComplexArray,
    sequence: FPAASequence,
    input_state: ComplexArray | None = None,
) -> ComplexArray:
    state = np.array(source if input_state is None else input_state, dtype=np.complex128)
    source_state = np.asarray(source, dtype=np.complex128)
    target_state = np.asarray(target, dtype=np.complex128)
    for alpha, beta in zip(sequence.alpha, sequence.beta, strict=True):
        state = apply_selective_phase(state, target_state, float(beta))
        state = apply_selective_phase(state, source_state, float(-alpha))
    return state


def pure_trace_distance(left: ComplexArray, right: ComplexArray) -> float:
    left_normalized = normalized_state(left)
    right_normalized = normalized_state(right)
    overlap = np.vdot(right_normalized, left_normalized)
    orthogonal_residual = left_normalized - overlap * right_normalized
    return min(1.0, float(np.linalg.norm(orthogonal_residual)))
