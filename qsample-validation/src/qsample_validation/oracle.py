"""Deterministic row-wise Householder extension of a transition oracle."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class HouseholderOracle:
    transition: FloatArray
    vectors: FloatArray

    @classmethod
    def from_transition(cls, transition: FloatArray) -> HouseholderOracle:
        targets = np.sqrt(np.asarray(transition, dtype=np.float64))
        origin = np.zeros(transition.shape[0], dtype=np.float64)
        origin[0] = 1.0
        differences = origin[None, :] - targets
        norms = np.linalg.norm(differences, axis=1)
        vectors = np.zeros_like(differences)
        nontrivial = norms > 1.0e-15
        vectors[nontrivial] = differences[nontrivial] / norms[nontrivial, None]
        return cls(transition=np.asarray(transition), vectors=vectors)

    @property
    def n_states(self) -> int:
        return self.transition.shape[0]

    def apply(self, state: ComplexArray) -> ComplexArray:
        """Apply the controlled row Householders on the last two axes (A,B)."""

        array = np.asarray(state, dtype=np.complex128)
        if array.shape[-2:] != (self.n_states, self.n_states):
            raise ValueError("state must end in axes (A,B)")
        overlaps = np.einsum("...xy,xy->...x", array, self.vectors, optimize=True)
        return array - 2.0 * overlaps[..., :, None] * self.vectors

    def clean_column(self) -> ComplexArray:
        state = np.zeros((self.n_states, self.n_states), dtype=np.complex128)
        state[:, 0] = 1.0
        return self.apply(state)
