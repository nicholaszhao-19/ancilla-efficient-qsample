"""Complete common-interface selective-phase compiler."""

from __future__ import annotations

import cmath
import math

import numpy as np
from numpy.typing import NDArray

from .gqsp_filter import FilterSpec, evaluate_filter
from .szegedy import SpectralWalk

ComplexArray = NDArray[np.complex128]


def _filter_value(z: complex, spec: FilterSpec) -> complex:
    return complex(evaluate_filter(z, spec))


def apply_filter_completion(
    state: ComplexArray,
    walk: SpectralWalk,
    spec: FilterSpec,
    *,
    adjoint: bool = False,
) -> ComplexArray:
    """Apply the documented canonical scalar completion of the GQSP block."""

    signal_zero = np.asarray(state[0], dtype=np.complex128)
    signal_one = np.asarray(state[1], dtype=np.complex128)

    def polynomial(z: complex | ComplexArray) -> complex | ComplexArray:
        return evaluate_filter(z, spec)

    def polynomial_adjoint(z: complex | ComplexArray) -> complex | ComplexArray:
        values = np.asarray(evaluate_filter(z, spec), dtype=np.complex128)
        result = np.conjugate(values)
        if np.isscalar(z):
            return complex(result.item())
        return result

    def complement(z: complex | ComplexArray) -> complex | ComplexArray:
        values = np.asarray(evaluate_filter(z, spec), dtype=np.complex128)
        result = np.sqrt(np.maximum(0.0, 1.0 - np.abs(values) ** 2)).astype(
            np.complex128
        )
        if np.isscalar(z):
            return complex(result.item())
        return result

    if not adjoint:
        output_zero = walk.apply_function(signal_zero, polynomial)
        output_zero -= walk.apply_function(signal_one, complement)
        output_one = walk.apply_function(signal_zero, complement)
        output_one += walk.apply_function(signal_one, polynomial_adjoint)
    else:
        output_zero = walk.apply_function(signal_zero, polynomial_adjoint)
        output_zero += walk.apply_function(signal_one, complement)
        output_one = -walk.apply_function(signal_zero, complement)
        output_one += walk.apply_function(signal_one, polynomial)
    return np.stack((output_zero, output_one))


def apply_compiled_selective_phase(
    state: ComplexArray,
    walk: SpectralWalk,
    spec: FilterSpec,
    phase: float,
) -> ComplexArray:
    lifted = walk.oracle.apply(np.asarray(state, dtype=np.complex128))
    computed = apply_filter_completion(lifted, walk, spec)
    computed[0] *= cmath.exp(1.0j * phase)
    uncomputed = apply_filter_completion(computed, walk, spec, adjoint=True)
    return walk.oracle.apply(uncomputed)


def exact_restricted_error(walk: SpectralWalk, spec: FilterSpec, phase: float) -> float:
    leakages: list[float] = []
    for index, eigenvalue in enumerate(walk.eigenvalues):
        if index == walk.stationary_index:
            continue
        angle = math.acos(float(np.clip(eigenvalue, -1.0, 1.0)))
        leakages.append(abs(_filter_value(cmath.exp(1.0j * angle), spec)))
    return abs(cmath.exp(1.0j * phase) - 1.0) * max(leakages, default=0.0)
