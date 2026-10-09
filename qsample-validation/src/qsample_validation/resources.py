"""Exact resource accounting for the compiled FPAA unitary."""

from __future__ import annotations

from dataclasses import dataclass

from .fpaa import FPAASequence
from .gqsp_filter import FilterSpec


@dataclass(frozen=True)
class StageResources:
    controlled_walk_calls: int
    wrapper_oracle_calls: int
    transition_oracle_calls: int


def stage_resources(
    sequence: FPAASequence,
    source_filter: FilterSpec,
    target_filter: FilterSpec,
) -> StageResources:
    pairs = (sequence.length - 1) // 2
    controlled_walk_calls = 2 * pairs * (source_filter.degree + target_filter.degree)
    wrapper_oracle_calls = 2 * sequence.n_phases
    transition_oracle_calls = 2 * controlled_walk_calls + wrapper_oracle_calls
    return StageResources(
        controlled_walk_calls=controlled_walk_calls,
        wrapper_oracle_calls=wrapper_oracle_calls,
        transition_oracle_calls=transition_oracle_calls,
    )
