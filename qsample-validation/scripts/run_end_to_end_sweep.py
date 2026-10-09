#!/usr/bin/env python
"""Generate compiled spectral-reference propagation data and resources."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from qsample_validation.benchmark import build_from_config, load_config
from qsample_validation.end_to_end import EndToEndResult, run_end_to_end
from qsample_validation.selective_phase import exact_restricted_error
from qsample_validation.szegedy import SpectralWalk
from qsample_validation.validation import (
    validate_chain,
    validate_compiler_error,
    validate_end_to_end,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "ising_2x4.yaml"
END_OUTPUT = ROOT / "outputs" / "data" / "end_to_end_sweep.csv"
STAGE_OUTPUT = ROOT / "outputs" / "data" / "benchmark_stages.csv"
CORRECTED_END_OUTPUT = ROOT / "outputs" / "data" / "corrected_end_to_end_sweep.csv"
CORRECTED_STAGE_OUTPUT = ROOT / "outputs" / "data" / "corrected_benchmark_stages.csv"


def _format_float(value: float) -> str:
    """Keep more than paper precision without exposing unstable last bits."""

    return f"{value:.10g}"


def _write_end_results(results: list[EndToEndResult], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "requested_epsilon",
        "compiled_joint_trace_distance",
        "compiled_reduced_trace_distance",
        "compiled_total_variation_distance",
        "ideal_joint_trace_distance",
        "compiled_to_ideal_joint_trace_distance",
        "compiled_clean_subspace_leakage",
        "compiled_norm_residual",
        "total_controlled_walk_calls",
        "total_wrapper_oracle_calls",
        "total_transition_oracle_calls",
    ]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "requested_epsilon": _format_float(result.epsilon),
                    "compiled_joint_trace_distance": _format_float(result.joint_error),
                    "compiled_reduced_trace_distance": _format_float(result.reduced_error),
                    "compiled_total_variation_distance": _format_float(result.total_variation),
                    "ideal_joint_trace_distance": _format_float(result.ideal_joint_error),
                    "compiled_to_ideal_joint_trace_distance": _format_float(
                        result.compiled_to_ideal_joint_distance
                    ),
                    "compiled_clean_subspace_leakage": _format_float(result.leakage),
                    "compiled_norm_residual": _format_float(result.norm_residual),
                    "total_controlled_walk_calls": sum(
                        stage.resources.controlled_walk_calls for stage in result.stages
                    ),
                    "total_wrapper_oracle_calls": sum(
                        stage.resources.wrapper_oracle_calls for stage in result.stages
                    ),
                    "total_transition_oracle_calls": sum(
                        stage.resources.transition_oracle_calls for stage in result.stages
                    ),
                }
            )


def _write_stage_results(result: EndToEndResult, chains: list, output: Path) -> None:
    fields = [
        "stage",
        "beta_i",
        "beta_i_plus_1",
        "p_i_exact",
        "delta_i",
        "phase_gap_i",
        "delta_i_plus_1",
        "phase_gap_i_plus_1",
        "L_i",
        "epsilon_i_fpaa",
        "phase_weight_A_i",
        "eta_i",
        "t_i",
        "n_i",
        "d_i",
        "t_i_plus_1",
        "n_i_plus_1",
        "d_i_plus_1",
        "compiler_error_i_phi_pi",
        "compiler_error_i_plus_1_phi_pi",
        "compiled_local_joint_error",
        "compiled_cumulative_joint_error",
        "phase_sensitive_compiler_hybrid_bound",
        "controlled_walk_calls",
        "wrapper_oracle_calls",
        "transition_oracle_calls",
        "compiled_cumulative_clean_subspace_leakage",
    ]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for stage in result.stages:
            source = chains[stage.index]
            target = chains[stage.index + 1]
            source_error = exact_restricted_error(
                SpectralWalk(source), stage.source_filter, math.pi
            )
            target_error = exact_restricted_error(
                SpectralWalk(target), stage.target_filter, math.pi
            )
            writer.writerow(
                {
                    "stage": stage.index,
                    "beta_i": _format_float(stage.beta_source),
                    "beta_i_plus_1": _format_float(stage.beta_target),
                    "p_i_exact": _format_float(stage.overlap),
                    "delta_i": _format_float(source.spectral_gap),
                    "phase_gap_i": _format_float(source.walk_phase_gap),
                    "delta_i_plus_1": _format_float(target.spectral_gap),
                    "phase_gap_i_plus_1": _format_float(target.walk_phase_gap),
                    "L_i": stage.sequence.length,
                    "epsilon_i_fpaa": _format_float(stage.sequence.target_error),
                    "phase_weight_A_i": _format_float(stage.phase_weight),
                    "eta_i": _format_float(stage.eta),
                    "t_i": stage.source_filter.t,
                    "n_i": stage.source_filter.n,
                    "d_i": stage.source_filter.degree,
                    "t_i_plus_1": stage.target_filter.t,
                    "n_i_plus_1": stage.target_filter.n,
                    "d_i_plus_1": stage.target_filter.degree,
                    "compiler_error_i_phi_pi": _format_float(source_error),
                    "compiler_error_i_plus_1_phi_pi": _format_float(target_error),
                    "compiled_local_joint_error": _format_float(stage.local_joint_error),
                    "compiled_cumulative_joint_error": (
                        _format_float(stage.cumulative_joint_error)
                    ),
                    "phase_sensitive_compiler_hybrid_bound": _format_float(
                        stage.phase_sensitive_compiler_bound
                    ),
                    "controlled_walk_calls": (stage.resources.controlled_walk_calls),
                    "wrapper_oracle_calls": (stage.resources.wrapper_oracle_calls),
                    "transition_oracle_calls": (stage.resources.transition_oracle_calls),
                    "compiled_cumulative_clean_subspace_leakage": (
                        _format_float(stage.cumulative_leakage)
                    ),
                }
            )


def main() -> None:
    config = load_config(CONFIG)
    _, chains = build_from_config(config)
    atol = float(config["numerics"]["atol"])
    for chain in chains:
        validate_chain(chain, atol)
    epsilons = [float(value) for value in config["numerics"]["sweep_epsilon"]]
    results = [run_end_to_end(chains, epsilon) for epsilon in epsilons]
    for result in results:
        validate_end_to_end(result, atol)
        for stage in result.stages:
            source_error = exact_restricted_error(
                SpectralWalk(chains[stage.index]),
                stage.source_filter,
                math.pi,
            )
            target_error = exact_restricted_error(
                SpectralWalk(chains[stage.index + 1]),
                stage.target_filter,
                math.pi,
            )
            validate_compiler_error(source_error, 2.0 * stage.eta, atol)
            validate_compiler_error(target_error, 2.0 * stage.eta, atol)
    _write_end_results(results, END_OUTPUT)
    reference_epsilon = float(config["numerics"]["benchmark_epsilon"])
    reference = run_end_to_end(chains, reference_epsilon)
    validate_end_to_end(reference, atol)
    _write_stage_results(reference, chains, STAGE_OUTPUT)

    corrected_results = [
        run_end_to_end(chains, epsilon, filter_allocation="uniform")
        for epsilon in epsilons
    ]
    for result in corrected_results:
        validate_end_to_end(result, atol)
    _write_end_results(corrected_results, CORRECTED_END_OUTPUT)
    corrected_reference = run_end_to_end(
        chains,
        reference_epsilon,
        filter_allocation="uniform",
    )
    validate_end_to_end(corrected_reference, atol)
    _write_stage_results(corrected_reference, chains, CORRECTED_STAGE_OUTPUT)
    print(END_OUTPUT)
    print(STAGE_OUTPUT)
    print(CORRECTED_END_OUTPUT)
    print(CORRECTED_STAGE_OUTPUT)


if __name__ == "__main__":
    main()
