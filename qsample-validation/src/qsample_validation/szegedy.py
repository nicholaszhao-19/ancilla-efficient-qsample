"""Matrix-free Szegedy walk and exact functional calculus."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .markov import MarkovChain
from .oracle import HouseholderOracle

ComplexArray = NDArray[np.complex128]
ComplexFunction = Callable[[complex | ComplexArray], complex | ComplexArray]
ENDPOINT_ATOL = 64.0 * np.finfo(float).eps
MIN_SUPPORTED_SPECTRAL_GAP = float(np.cbrt(np.finfo(float).eps))


@dataclass
class SpectralWalk:
    """Exact walk functions without constructing the dense |Omega|^2 matrix."""

    chain: MarkovChain

    def __post_init__(self) -> None:
        if (
            not np.isfinite(self.chain.spectral_gap)
            or self.chain.spectral_gap < MIN_SUPPORTED_SPECTRAL_GAP
        ):
            raise ValueError(
                "spectral gap is below the stable complex128 range of this "
                "nonorthogonal spectral reference"
            )
        self.oracle = HouseholderOracle.from_transition(self.chain.transition)
        self.sqrt_transition = np.sqrt(self.chain.transition)
        self.eigenvalues = self.chain.eigenvalues
        self.eigenvectors = self.chain.eigenvectors
        self.stationary_index = int(np.argmax(self.eigenvalues))

    @property
    def n_states(self) -> int:
        return self.chain.transition.shape[0]

    def apply_t(self, logical_state: ComplexArray) -> ComplexArray:
        logical = np.asarray(logical_state, dtype=np.complex128)
        return logical[:, None] * self.sqrt_transition

    def apply_t_adjoint(self, walk_state: ComplexArray) -> ComplexArray:
        state = np.asarray(walk_state, dtype=np.complex128)
        return np.einsum("xy,xy->x", self.sqrt_transition, state, optimize=True)

    def apply_walk(self, walk_state: ComplexArray) -> ComplexArray:
        swapped = np.asarray(walk_state, dtype=np.complex128).T.copy()
        lifted = self.oracle.apply(swapped)
        reflected = -lifted
        reflected[:, 0] = lifted[:, 0]
        return self.oracle.apply(reflected)

    def _active_coordinates(
        self, walk_state: ComplexArray
    ) -> tuple[ComplexArray, ComplexArray, ComplexArray, ComplexArray]:
        state = np.asarray(walk_state, dtype=np.complex128)
        overlap_a = np.einsum(
            "xk,x->k",
            self.eigenvectors,
            self.apply_t_adjoint(state),
            optimize=True,
        )
        overlap_b = np.einsum(
            "xk,x->k",
            self.eigenvectors,
            self.apply_t_adjoint(state.T),
            optimize=True,
        )
        coefficients_a = np.zeros(self.n_states, dtype=np.complex128)
        coefficients_b = np.zeros(self.n_states, dtype=np.complex128)
        for index, eigenvalue in enumerate(self.eigenvalues):
            if (
                abs(eigenvalue - 1.0) <= ENDPOINT_ATOL
                or abs(eigenvalue + 1.0) <= ENDPOINT_ATOL
            ):
                coefficients_a[index] = overlap_a[index]
                continue
            denominator = 1.0 - eigenvalue * eigenvalue
            coefficients_a[index] = (
                overlap_a[index] - eigenvalue * overlap_b[index]
            ) / denominator
            coefficients_b[index] = (
                overlap_b[index] - eigenvalue * overlap_a[index]
            ) / denominator
        logical_a = np.einsum(
            "xk,k->x", self.eigenvectors, coefficients_a, optimize=True
        )
        logical_b = np.einsum(
            "xk,k->x", self.eigenvectors, coefficients_b, optimize=True
        )
        active = self.apply_t(logical_a)
        active += self.apply_t(logical_b).T
        return coefficients_a, coefficients_b, active, overlap_a

    def apply_function(
        self,
        walk_state: ComplexArray,
        function: ComplexFunction,
    ) -> ComplexArray:
        """Apply ``function(W)`` exactly from the discriminant spectrum."""

        state = np.asarray(walk_state, dtype=np.complex128)
        coefficients_a, coefficients_b, active, overlap_a = self._active_coordinates(
            state
        )
        transformed_a = np.zeros_like(coefficients_a)
        transformed_b = np.zeros_like(coefficients_b)

        for index, eigenvalue in enumerate(self.eigenvalues):
            if abs(eigenvalue - 1.0) <= ENDPOINT_ATOL:
                transformed_a[index] = complex(function(1.0 + 0.0j)) * overlap_a[index]
                continue
            if abs(eigenvalue + 1.0) <= ENDPOINT_ATOL:
                transformed_a[index] = complex(function(-1.0 + 0.0j)) * overlap_a[index]
                continue

            sine = np.sqrt(max(0.0, 1.0 - eigenvalue * eigenvalue))
            z_plus = complex(eigenvalue, sine)
            z_minus = complex(eigenvalue, -sine)
            f_plus = complex(function(z_plus))
            f_minus = complex(function(z_minus))
            denominator = z_plus - z_minus
            matrix = np.array(
                [[2.0 * eigenvalue, 1.0], [-1.0, 0.0]],
                dtype=np.complex128,
            )
            identity = np.eye(2, dtype=np.complex128)
            matrix_function = (
                f_plus * (matrix - z_minus * identity)
                - f_minus * (matrix - z_plus * identity)
            ) / denominator
            transformed = np.einsum(
                "ij,j->i",
                matrix_function,
                np.array([coefficients_a[index], coefficients_b[index]]),
                optimize=True,
            )
            transformed_a[index], transformed_b[index] = transformed

        transformed_logical_a = np.einsum(
            "xk,k->x", self.eigenvectors, transformed_a, optimize=True
        )
        transformed_logical_b = np.einsum(
            "xk,k->x", self.eigenvectors, transformed_b, optimize=True
        )
        transformed_active = self.apply_t(transformed_logical_a)
        transformed_active += self.apply_t(transformed_logical_b).T

        residual = state - active
        symmetric = (residual + residual.T) / 2.0
        antisymmetric = (residual - residual.T) / 2.0
        transformed_residual = complex(function(-1.0 + 0.0j)) * symmetric
        transformed_residual += complex(function(1.0 + 0.0j)) * antisymmetric
        return transformed_active + transformed_residual

    def projector_one(self, walk_state: ComplexArray) -> ComplexArray:
        def indicator(z: complex | ComplexArray) -> complex | ComplexArray:
            values = np.asarray(z, dtype=np.complex128)
            result = (values == (1.0 + 0.0j)).astype(np.complex128)
            if np.isscalar(z):
                return complex(result.item())
            return result

        return self.apply_function(walk_state, indicator)

    def random_unitarity_residual(self, seed: int = 0) -> float:
        generator = np.random.default_rng(seed)
        state = generator.normal(size=(self.n_states, self.n_states))
        state = state + 1.0j * generator.normal(size=state.shape)
        state = state.astype(np.complex128)
        state /= np.linalg.norm(state)
        forward = self.apply_walk(state)
        backward = self.apply_function(forward, lambda z: 1.0 / z)
        return float(np.linalg.norm(backward - state))
