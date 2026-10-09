"""Finite open-boundary Ising-ladder enumeration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class IsingModel:
    """Disagreement-count Ising model with a fixed configuration ordering."""

    rows: int
    cols: int
    disagreement_cost: float = 1.0
    external_field_h: float = 0.0

    @property
    def n_sites(self) -> int:
        return self.rows * self.cols

    @property
    def n_states(self) -> int:
        return 1 << self.n_sites

    @property
    def edges(self) -> tuple[tuple[int, int], ...]:
        edges: list[tuple[int, int]] = []
        for row in range(self.rows):
            for col in range(self.cols):
                site = row * self.cols + col
                if col + 1 < self.cols:
                    edges.append((site, site + 1))
                if row + 1 < self.rows:
                    edges.append((site, site + self.cols))
        return tuple(edges)

    def configurations(self) -> IntArray:
        """Return row-major lexicographic bit strings mapped to spins -1,+1."""

        state_ids = np.arange(self.n_states, dtype=np.uint64)[:, None]
        shifts = np.arange(self.n_sites - 1, -1, -1, dtype=np.uint64)[None, :]
        bits = ((state_ids >> shifts) & 1).astype(np.int64)
        return 2 * bits - 1

    def index_from_spins(self, spins: IntArray) -> IntArray:
        """Map one or more spin configurations back to lexicographic indices."""

        bits = ((np.asarray(spins, dtype=np.int64) + 1) // 2).astype(np.int64)
        weights = (1 << np.arange(self.n_sites - 1, -1, -1)).astype(np.int64)
        return bits @ weights

    def energies(self) -> FloatArray:
        spins = self.configurations()
        disagreements = np.zeros(self.n_states, dtype=np.float64)
        for left, right in self.edges:
            disagreements += spins[:, left] != spins[:, right]
        field_term = -self.external_field_h * spins.sum(axis=1)
        return self.disagreement_cost * disagreements + field_term

    def flipped_indices(self, site: int, spin: int) -> IntArray:
        """Indices obtained by fixing ``site`` to ``spin`` in every state."""

        configs = self.configurations().copy()
        configs[:, site] = spin
        return self.index_from_spins(configs)

