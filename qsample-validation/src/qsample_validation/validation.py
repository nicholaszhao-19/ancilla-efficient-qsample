"""Fail-fast acceptance checks for publication output generation."""

from __future__ import annotations

import math

import numpy as np

from .end_to_end import EndToEndResult
from .markov import MarkovChain


def validate_chain(chain: MarkovChain, atol: float) -> None:
    """Raise when a generated chain violates a benchmark invariant."""

    arrays = (
        chain.probabilities,
        chain.transition,
        chain.discriminant,
        chain.eigenvalues,
        chain.eigenvectors,
    )
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise RuntimeError(f"beta={chain.beta}: chain contains non-finite values")
    if float(np.min(chain.probabilities)) < -atol:
        raise RuntimeError(f"beta={chain.beta}: Gibbs probabilities are negative")
    if abs(float(np.sum(chain.probabilities)) - 1.0) > atol:
        raise RuntimeError(f"beta={chain.beta}: Gibbs probabilities are not normalized")
    if float(np.min(chain.transition)) < -atol:
        raise RuntimeError(f"beta={chain.beta}: transition matrix has negative entries")
    if (
        max(
            chain.row_residual,
            chain.stationarity_residual,
            chain.detailed_balance_residual,
        )
        > atol
    ):
        raise RuntimeError(f"beta={chain.beta}: Markov-chain residual exceeds {atol:g}")
    if float(np.max(np.abs(chain.discriminant - chain.discriminant.T))) > atol:
        raise RuntimeError(f"beta={chain.beta}: discriminant is not symmetric")
    if abs(float(chain.eigenvalues[-1]) - 1.0) > atol:
        raise RuntimeError(f"beta={chain.beta}: stationary eigenvalue is not one")
    if not 0.0 < chain.spectral_gap <= 1.0:
        raise RuntimeError(f"beta={chain.beta}: invalid lazy-chain spectral gap")
    if not 0.0 < chain.walk_phase_gap <= math.pi:
        raise RuntimeError(f"beta={chain.beta}: invalid walk phase gap")
    stationary = np.sqrt(chain.probabilities)
    leading = chain.eigenvectors[:, -1]
    alignment_residual = min(
        float(np.linalg.norm(leading - stationary)),
        float(np.linalg.norm(leading + stationary)),
    )
    if alignment_residual > atol:
        raise RuntimeError(f"beta={chain.beta}: stationary eigenspace is not numerically resolved")


def validate_compiler_error(error: float, envelope: float, atol: float) -> None:
    if not math.isfinite(error) or error < 0.0:
        raise RuntimeError("compiler error is not finite and nonnegative")
    if error > envelope + atol:
        raise RuntimeError(f"compiler error {error:.6g} exceeds envelope {envelope:.6g}")


def validate_end_to_end(result: EndToEndResult, atol: float) -> None:
    """Check theorem budgets, metric ordering, and exact resource identities."""

    metrics = (
        result.joint_error,
        result.reduced_error,
        result.total_variation,
        result.leakage,
        result.norm_residual,
        result.ideal_joint_error,
        result.compiled_to_ideal_joint_distance,
    )
    if not all(math.isfinite(value) for value in metrics):
        raise RuntimeError(f"epsilon={result.epsilon}: non-finite output metric")
    if min(metrics) < 0.0:
        raise RuntimeError(f"epsilon={result.epsilon}: negative output metric")
    if (
        max(
            result.joint_error,
            result.reduced_error,
            result.total_variation,
        )
        > result.epsilon + atol
    ):
        raise RuntimeError(f"epsilon={result.epsilon}: final error budget failed")
    if result.reduced_error > result.joint_error + atol:
        raise RuntimeError(f"epsilon={result.epsilon}: trace-distance contractivity failed")
    if result.total_variation > result.reduced_error + atol:
        raise RuntimeError(f"epsilon={result.epsilon}: measurement contractivity failed")
    if result.leakage > 1.0 + atol:
        raise RuntimeError(f"epsilon={result.epsilon}: leakage is not a probability")
    if result.norm_residual > atol:
        raise RuntimeError(f"epsilon={result.epsilon}: state normalization failed")
    if result.ideal_joint_error > result.epsilon / 2.0 + atol:
        raise RuntimeError(f"epsilon={result.epsilon}: ideal FPAA control failed")
    if result.compiled_to_ideal_joint_distance > result.epsilon / 2.0 + atol:
        raise RuntimeError(f"epsilon={result.epsilon}: compiled FPAA control failed")
    if (
        result.joint_error
        > result.ideal_joint_error + result.compiled_to_ideal_joint_distance + atol
    ):
        raise RuntimeError(f"epsilon={result.epsilon}: control triangle failed")

    n_stages = len(result.stages)
    stage_budget = result.epsilon / n_stages
    compiler_stage_budget = result.epsilon / (2.0 * n_stages)
    for stage in result.stages:
        if stage.filter_allocation == "phase_weighted":
            denominator = stage.phase_weight
        elif stage.filter_allocation == "uniform":
            denominator = stage.sequence.length - 1
        else:
            raise RuntimeError(f"stage={stage.index}: unknown filter allocation")
        expected_eta = result.epsilon / (4.0 * n_stages * denominator)
        if not math.isclose(stage.eta, expected_eta, rel_tol=1.0e-12, abs_tol=0.0):
            raise RuntimeError(f"stage={stage.index}: filter allocation drifted")
        if stage.local_joint_error > stage_budget + atol:
            raise RuntimeError(f"stage={stage.index}: local stage budget failed")
        if stage.cumulative_joint_error > (stage.index + 1) * stage_budget + atol:
            raise RuntimeError(f"stage={stage.index}: cumulative stage budget failed")
        if (
            not math.isfinite(stage.phase_sensitive_compiler_bound)
            or stage.phase_sensitive_compiler_bound < 0.0
        ):
            raise RuntimeError(f"stage={stage.index}: invalid phase-sensitive compiler bound")
        if stage.phase_sensitive_compiler_bound > compiler_stage_budget + atol:
            raise RuntimeError(f"stage={stage.index}: phase-sensitive compiler budget failed")
        if stage.phase_sensitive_compiler_bound > 2.0 * stage.phase_weight * stage.eta + atol:
            raise RuntimeError(f"stage={stage.index}: phase-sensitive theorem bound failed")
        resources = stage.resources
        if (
            resources.transition_oracle_calls
            != 2 * resources.controlled_walk_calls + resources.wrapper_oracle_calls
        ):
            raise RuntimeError(f"stage={stage.index}: resource identity failed")
        for spec in (stage.source_filter, stage.target_filter):
            if spec.degree != spec.n * (spec.t - 1):
                raise RuntimeError(f"stage={stage.index}: filter degree identity failed")
