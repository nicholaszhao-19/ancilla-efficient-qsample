"""Sequential approximate compiled-phase FPAA and exact output metrics."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .fpaa import FPAASequence, apply_ideal_fpaa, pure_trace_distance, yoder_sequence
from .gqsp_filter import FilterSpec, design_filter
from .markov import MarkovChain, squared_qsample_overlap
from .metrics import (
    clean_joint_state,
    clean_subspace_leakage,
    joint_trace_distance,
    norm_residual,
    reduced_trace_distance,
    total_variation_distance,
)
from .resources import StageResources, stage_resources
from .selective_phase import (
    apply_compiled_selective_phase,
    exact_restricted_error,
)
from .szegedy import SpectralWalk

ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True)
class StageResult:
    index: int
    beta_source: float
    beta_target: float
    overlap: float
    sequence: FPAASequence
    filter_allocation: Literal["phase_weighted", "uniform"]
    phase_weight: float
    eta: float
    source_filter: FilterSpec
    target_filter: FilterSpec
    resources: StageResources
    local_joint_error: float
    cumulative_joint_error: float
    cumulative_leakage: float
    phase_sensitive_compiler_bound: float


@dataclass(frozen=True)
class EndToEndResult:
    """Compiled production result with a separately labelled ideal control."""

    epsilon: float
    state: ComplexArray
    stages: tuple[StageResult, ...]
    joint_error: float
    reduced_error: float
    total_variation: float
    leakage: float
    norm_residual: float
    ideal_joint_error: float
    compiled_to_ideal_joint_distance: float


def _apply_stage(
    state: ComplexArray,
    source_walk: SpectralWalk,
    target_walk: SpectralWalk,
    source_filter: FilterSpec,
    target_filter: FilterSpec,
    sequence: FPAASequence,
) -> ComplexArray:
    result = np.asarray(state, dtype=np.complex128)
    for alpha, beta in zip(sequence.alpha, sequence.beta, strict=True):
        result = apply_compiled_selective_phase(result, target_walk, target_filter, float(beta))
        result = apply_compiled_selective_phase(result, source_walk, source_filter, float(-alpha))
    return result


def _embed_logical_state(amplitudes: ComplexArray) -> ComplexArray:
    """Embed an arbitrary logical state in the common clean joint space."""

    logical = np.asarray(amplitudes, dtype=np.complex128)
    state = np.zeros((2, logical.size, logical.size), dtype=np.complex128)
    state[0, :, 0] = logical
    return state


def run_end_to_end(
    chains: list[MarkovChain],
    epsilon: float,
    *,
    filter_allocation: Literal["phase_weighted", "uniform"] = "phase_weighted",
) -> EndToEndResult:
    """Run compiled propagation and a non-headline ideal-projector control."""

    if len(chains) < 2:
        raise ValueError("at least two chains are required")
    if not 0.0 < epsilon < 1.0:
        raise ValueError("epsilon must lie in (0, 1)")
    if filter_allocation not in {"phase_weighted", "uniform"}:
        raise ValueError("unknown filter allocation")

    walks = [SpectralWalk(chain) for chain in chains]
    n_stages = len(chains) - 1
    state = clean_joint_state(chains[0].probabilities)
    ideal_state = np.sqrt(chains[0].probabilities).astype(np.complex128)
    stage_results: list[StageResult] = []

    for index in range(n_stages):
        source = chains[index]
        target = chains[index + 1]
        overlap = squared_qsample_overlap(source.probabilities, target.probabilities)
        target_fpaa_error = epsilon / (2.0 * n_stages)
        sequence = yoder_sequence(overlap, target_fpaa_error)
        phase_weight = sequence.phase_weight
        if filter_allocation == "phase_weighted":
            eta = epsilon / (4.0 * n_stages * phase_weight)
        else:
            eta = epsilon / (4.0 * n_stages * (sequence.length - 1))
        source_filter = design_filter(source.walk_phase_gap, eta)
        target_filter = design_filter(target.walk_phase_gap, eta)
        resources = stage_resources(sequence, source_filter, target_filter)
        phase_sensitive_compiler_bound = math.fsum(
            exact_restricted_error(walks[index + 1], target_filter, float(beta))
            + exact_restricted_error(walks[index], source_filter, float(-alpha))
            for alpha, beta in zip(sequence.alpha, sequence.beta, strict=True)
        )

        local_input = clean_joint_state(source.probabilities)
        local_output = _apply_stage(
            local_input,
            walks[index],
            walks[index + 1],
            source_filter,
            target_filter,
            sequence,
        )
        local_error = joint_trace_distance(local_output, target.probabilities)

        state = _apply_stage(
            state,
            walks[index],
            walks[index + 1],
            source_filter,
            target_filter,
            sequence,
        )
        ideal_state = apply_ideal_fpaa(
            np.sqrt(source.probabilities),
            np.sqrt(target.probabilities),
            sequence,
            input_state=ideal_state,
        )
        stage_results.append(
            StageResult(
                index=index,
                beta_source=source.beta,
                beta_target=target.beta,
                overlap=overlap,
                sequence=sequence,
                filter_allocation=filter_allocation,
                phase_weight=phase_weight,
                eta=eta,
                source_filter=source_filter,
                target_filter=target_filter,
                resources=resources,
                local_joint_error=local_error,
                cumulative_joint_error=joint_trace_distance(state, target.probabilities),
                cumulative_leakage=clean_subspace_leakage(state),
                phase_sensitive_compiler_bound=(phase_sensitive_compiler_bound),
            )
        )

    target_probabilities = chains[-1].probabilities
    ideal_joint_state = _embed_logical_state(ideal_state)
    return EndToEndResult(
        epsilon=float(epsilon),
        state=state,
        stages=tuple(stage_results),
        joint_error=joint_trace_distance(state, target_probabilities),
        reduced_error=reduced_trace_distance(state, target_probabilities),
        total_variation=total_variation_distance(state, target_probabilities),
        leakage=clean_subspace_leakage(state),
        norm_residual=norm_residual(state),
        ideal_joint_error=joint_trace_distance(ideal_joint_state, target_probabilities),
        compiled_to_ideal_joint_distance=pure_trace_distance(state, ideal_joint_state),
    )
