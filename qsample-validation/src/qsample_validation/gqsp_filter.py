"""Explicit average-power eigenspace filter used by the spectral reference."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True)
class FilterSpec:
    gap: float
    requested_eta: float
    t: int
    n: int
    degree: int
    certified_eta: float


def design_filter(gap: float, eta: float) -> FilterSpec:
    if not 0.0 < gap <= math.pi:
        raise ValueError("gap must lie in (0, pi]")
    if not 0.0 < eta < 1.0:
        raise ValueError("eta must lie in (0, 1)")
    inverse_width = math.e / math.sin(gap / 2.0)
    if not math.isfinite(inverse_width):
        raise ValueError("gap is below the representable filter-design range")
    t = int(math.ceil(inverse_width))
    n = int(math.ceil(-math.log(eta)))
    return FilterSpec(
        gap=float(gap),
        requested_eta=float(eta),
        t=t,
        n=n,
        degree=n * (t - 1),
        certified_eta=math.exp(-n),
    )


def average_polynomial(z: complex | ComplexArray, t: int) -> complex | ComplexArray:
    if t < 1:
        raise ValueError("t must be a positive integer")
    values = np.asarray(z, dtype=np.complex128)
    result = np.empty_like(values)
    unit_circle = np.abs(np.abs(values) - 1.0) <= 32.0 * np.finfo(float).eps
    angles = np.angle(values[unit_circle])
    result[unit_circle] = (
        np.exp(0.5j * (t - 1) * angles)
        * np.sinc(t * angles / (2.0 * np.pi))
        / np.sinc(angles / (2.0 * np.pi))
    )

    zero = values == 0.0
    result[zero] = 1.0 / t
    off_circle = ~unit_circle & ~zero
    if np.any(off_circle):
        logarithms = np.log(values[off_circle])
        result[off_circle] = np.expm1(t * logarithms) / (
            t * np.expm1(logarithms)
        )
    if np.isscalar(z):
        return complex(result.item())
    return result


def evaluate_filter(
    z: complex | ComplexArray, spec: FilterSpec
) -> complex | ComplexArray:
    return average_polynomial(z, spec.t) ** spec.n


def scalar_completion(z: complex, spec: FilterSpec) -> ComplexArray:
    value = complex(evaluate_filter(z, spec))
    complement = math.sqrt(max(0.0, 1.0 - abs(value) ** 2))
    return np.array(
        [[value, -complement], [complement, np.conjugate(value)]],
        dtype=np.complex128,
    )
