#!/usr/bin/env python
"""Independently audit the saved benchmark and write a machine-readable record."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import platform
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from qsample_validation.benchmark import build_from_config, load_config
from qsample_validation.end_to_end import run_end_to_end
from qsample_validation.validation import validate_end_to_end

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
DATA = ROOT / "outputs" / "data"
OUTPUT = DATA / "publication_audit.json"
BETAS = np.array([0.0, 0.3, 0.6, 0.9, 1.2])
EXPECTED_HISTOGRAM = {0: 2, 2: 14, 3: 24, 4: 48, 5: 80, 6: 48, 7: 24, 8: 14, 10: 2}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _read_csv(path: Path, expected_fields: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        _require(tuple(reader.fieldnames or ()) == expected_fields, f"{path.name}: schema drift")
        rows = list(reader)
    _require(all(None not in row for row in rows), f"{path.name}: malformed row")
    return rows


def _independent_model() -> tuple[np.ndarray, list[tuple[int, int]]]:
    state_ids = np.arange(256, dtype=np.uint64)[:, None]
    shifts = np.arange(7, -1, -1, dtype=np.uint64)[None, :]
    spins = 2 * ((state_ids >> shifts) & 1).astype(np.int64) - 1
    edges = [(row * 4 + col, row * 4 + col + 1) for row in range(2) for col in range(3)]
    edges.extend((col, col + 4) for col in range(4))
    energies = sum((spins[:, left] != spins[:, right]) for left, right in edges)
    return np.asarray(energies, dtype=np.float64), edges


def _independent_chain(
    energies: np.ndarray, beta: float
) -> tuple[np.ndarray, np.ndarray, float, float]:
    weights = np.exp(-beta * energies)
    probabilities = weights / np.sum(weights)
    heat_bath = np.zeros((256, 256), dtype=np.float64)
    for state in range(256):
        for site in range(8):
            mask = 1 << (7 - site)
            minus = state & ~mask
            plus = state | mask
            local = np.exp(
                -beta
                * (
                    np.array([energies[minus], energies[plus]])
                    - min(energies[minus], energies[plus])
                )
            )
            local /= np.sum(local)
            heat_bath[state, minus] += local[0] / 8.0
            heat_bath[state, plus] += local[1] / 8.0
    transition = 0.5 * np.eye(256) + 0.5 * heat_bath
    roots = np.sqrt(probabilities)
    discriminant = roots[:, None] * transition / roots[None, :]
    discriminant = 0.5 * (discriminant + discriminant.T)
    eigenvalues = np.linalg.eigvalsh(discriminant)
    spectral_gap = 1.0 - float(eigenvalues[-2])
    phase_gap = math.acos(float(eigenvalues[-2]))
    return probabilities, eigenvalues, spectral_gap, phase_gap


def _filter_magnitude(phase: float, t: int, n: int) -> float:
    average = abs(np.sinc(t * phase / (2.0 * np.pi)) / np.sinc(phase / (2.0 * np.pi)))
    return float(average**n)


def _independent_fpaa(
    overlap: float, target_error: float
) -> tuple[int, np.ndarray, np.ndarray, float]:
    length = 3
    while math.tanh(math.acosh(1.0 / target_error) / length) ** 2 > overlap + 1.0e-15:
        length += 2
    count = (length - 1) // 2
    xi = 1.0 / math.cosh(math.acosh(1.0 / target_error) / length)
    root = math.sqrt(max(0.0, 1.0 - xi * xi))
    alpha = np.array(
        [
            2.0
            * math.atan2(
                1.0,
                math.tan(2.0 * math.pi * index / length) * root,
            )
            for index in range(1, count + 1)
        ]
    )
    beta = -alpha[::-1]
    phase_weight = float(np.sum(np.abs(np.sin(alpha / 2.0))) + np.sum(np.abs(np.sin(beta / 2.0))))
    return length, alpha, beta, phase_weight


def _apply_ideal_fpaa_control(
    state: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
    alpha: np.ndarray,
    beta: np.ndarray,
) -> np.ndarray:
    """Apply the ideal FPAA unitary without calling production routines."""

    result = np.asarray(state, dtype=np.complex128).copy()
    source_state = np.asarray(source, dtype=np.complex128)
    target_state = np.asarray(target, dtype=np.complex128)
    for source_phase, target_phase in zip(alpha, beta, strict=True):
        result += (np.exp(1.0j * target_phase) - 1.0) * target_state * np.vdot(target_state, result)
        result += (
            (np.exp(-1.0j * source_phase) - 1.0) * source_state * np.vdot(source_state, result)
        )
    return result


def _pure_trace_distance(left: np.ndarray, right: np.ndarray) -> float:
    """Evaluate pure-state distance by a stable projection residual."""

    left_state = np.asarray(left, dtype=np.complex128).copy()
    right_state = np.asarray(right, dtype=np.complex128).copy()
    left_state /= np.linalg.norm(left_state)
    right_state /= np.linalg.norm(right_state)
    overlap = np.vdot(right_state, left_state)
    return min(1.0, float(np.linalg.norm(left_state - overlap * right_state)))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _stable_json(value: Any) -> Any:
    """Remove numerically meaningless process-level last-bit variation."""

    if isinstance(value, float):
        return float(f"{value:.12g}")
    if isinstance(value, dict):
        return {key: _stable_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_stable_json(item) for item in value]
    return value


def _artifact_manifest() -> dict[str, str]:
    paths = [
        ROOT / "configs" / "ising_2x4.yaml",
        ROOT / "pyproject.toml",
        ROOT / "requirements-repro.txt",
        *sorted(ROOT.glob("*.md")),
        *sorted((ROOT / "src" / "qsample_validation").glob("*.py")),
        *sorted((ROOT / "scripts").glob("*.py")),
        *sorted((ROOT / "tests").glob("*.py")),
        *sorted(DATA.glob("*.csv")),
        *sorted((ROOT / "outputs" / "figures").glob("*")),
        *sorted((ROOT / "paper").glob("*")),
        REPOSITORY / "SAMPLE_Final.tex",
        REPOSITORY / "SAMPLE_Final.pdf",
    ]
    _require(all(path.is_file() for path in paths), "publication artifact is missing")
    return {str(path.relative_to(REPOSITORY)): _sha256(path) for path in paths if path != OUTPUT}


def main() -> None:
    chain_fields = (
        "beta",
        "spectral_gap",
        "walk_phase_gap",
        "overlap_to_next",
        "row_residual",
        "stationarity_residual",
        "detailed_balance_residual",
        "projector_intertwining_operator_norm",
    )
    compiler_fields = (
        "beta",
        "walk_phase_gap",
        "requested_eta",
        "certified_eta",
        "t",
        "n",
        "degree",
        "scaled_degree",
        "phase",
        "restricted_operator_error",
        "theorem_envelope",
    )
    end_fields = (
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
    )
    stage_fields = (
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
    )
    chain_rows = _read_csv(DATA / "benchmark_chains.csv", chain_fields)
    compiler_rows = _read_csv(DATA / "compiler_sweep.csv", compiler_fields)
    end_rows = _read_csv(DATA / "end_to_end_sweep.csv", end_fields)
    stage_rows = _read_csv(DATA / "benchmark_stages.csv", stage_fields)
    _require(
        (len(chain_rows), len(compiler_rows), len(end_rows), len(stage_rows)) == (5, 65, 6, 4),
        "unexpected CSV row count",
    )

    energies, edges = _independent_model()
    histogram = {int(key): int(value) for key, value in Counter(energies).items()}
    _require(len(edges) == 10, "independent ladder does not have ten edges")
    _require(histogram == EXPECTED_HISTOGRAM, "energy histogram mismatch")

    independent = [_independent_chain(energies, beta) for beta in BETAS]
    probabilities = [item[0] for item in independent]
    gaps = np.array([item[2] for item in independent])
    phase_gaps = np.array([item[3] for item in independent])
    overlaps = np.array(
        [
            float(np.dot(np.sqrt(probabilities[index]), np.sqrt(probabilities[index + 1])) ** 2)
            for index in range(4)
        ]
    )
    saved_gaps = np.array([float(row["spectral_gap"]) for row in chain_rows])
    saved_phases = np.array([float(row["walk_phase_gap"]) for row in chain_rows])
    saved_overlaps = np.array([float(row["overlap_to_next"]) for row in chain_rows[:-1]])
    gap_difference = float(np.max(np.abs(gaps - saved_gaps)))
    phase_difference = float(np.max(np.abs(phase_gaps - saved_phases)))
    overlap_difference = float(np.max(np.abs(overlaps - saved_overlaps)))
    _require(
        max(gap_difference, phase_difference, overlap_difference) < 2.0e-13,
        "independent chain audit failed",
    )

    independent_eigenvalues = {
        float(beta): independent[index][1] for index, beta in enumerate(BETAS)
    }
    max_compiler_ratio = 0.0
    max_compiler_difference = 0.0
    for row in compiler_rows:
        beta = float(row["beta"])
        eta = float(row["requested_eta"])
        gap = float(row["walk_phase_gap"])
        t = int(row["t"])
        n = int(row["n"])
        degree = int(row["degree"])
        _require(t == math.ceil(math.e / math.sin(gap / 2.0)), "filter t drift")
        _require(n == math.ceil(-math.log(eta)), "filter n drift")
        _require(degree == n * (t - 1), "filter degree drift")
        phases = np.arccos(np.clip(independent_eigenvalues[beta][:-1], -1.0, 1.0))
        independent_error = 2.0 * max(_filter_magnitude(float(value), t, n) for value in phases)
        saved_error = float(row["restricted_operator_error"])
        envelope = float(row["theorem_envelope"])
        max_compiler_difference = max(max_compiler_difference, abs(independent_error - saved_error))
        max_compiler_ratio = max(max_compiler_ratio, saved_error / envelope)
        _require(saved_error <= envelope + 1.0e-14, "compiler envelope failed")
    _require(max_compiler_difference < 2.0e-13, "independent compiler audit failed")

    reference_epsilon = 1.0e-3
    expected_transition_calls: list[int] = []
    independent_phase_bounds: list[float] = []
    max_phase_bound_difference = 0.0
    ideal_control_state = np.sqrt(probabilities[0]).astype(np.complex128)
    for index, row in enumerate(stage_rows):
        length, alpha, beta, phase_weight = _independent_fpaa(
            overlaps[index], reference_epsilon / 8.0
        )
        eta = reference_epsilon / (16.0 * phase_weight)
        _require(int(row["L_i"]) == length, "FPAA length drift")
        _require(
            abs(float(row["phase_weight_A_i"]) - phase_weight) < 1.0e-9,
            "phase weight drift",
        )
        _require(abs(float(row["eta_i"]) - eta) < 1.0e-13, "eta allocation drift")
        source_degree = int(row["d_i"])
        target_degree = int(row["d_i_plus_1"])
        controlled = 2 * ((length - 1) // 2) * (source_degree + target_degree)
        wrapper = 2 * (length - 1)
        transition_calls = 2 * controlled + wrapper
        _require(int(row["controlled_walk_calls"]) == controlled, "controlled-walk count drift")
        _require(int(row["wrapper_oracle_calls"]) == wrapper, "wrapper count drift")
        _require(int(row["transition_oracle_calls"]) == transition_calls, "query count drift")
        expected_transition_calls.append(transition_calls)

        source_beta = float(row["beta_i"])
        target_beta = float(row["beta_i_plus_1"])
        source_phases = np.arccos(np.clip(independent_eigenvalues[source_beta][:-1], -1.0, 1.0))
        target_phases = np.arccos(np.clip(independent_eigenvalues[target_beta][:-1], -1.0, 1.0))
        source_pi_error = 2.0 * max(
            _filter_magnitude(float(value), int(row["t_i"]), int(row["n_i"]))
            for value in source_phases
        )
        target_pi_error = 2.0 * max(
            _filter_magnitude(
                float(value),
                int(row["t_i_plus_1"]),
                int(row["n_i_plus_1"]),
            )
            for value in target_phases
        )
        independent_phase_bound = source_pi_error * math.fsum(
            abs(math.sin(float(value) / 2.0)) for value in alpha
        ) + target_pi_error * math.fsum(abs(math.sin(float(value) / 2.0)) for value in beta)
        saved_phase_bound = float(row["phase_sensitive_compiler_hybrid_bound"])
        max_phase_bound_difference = max(
            max_phase_bound_difference,
            abs(independent_phase_bound - saved_phase_bound),
        )
        _require(
            saved_phase_bound <= reference_epsilon / 8.0 + 1.0e-14,
            "phase-sensitive compiler budget failed",
        )
        independent_phase_bounds.append(independent_phase_bound)
        ideal_control_state = _apply_ideal_fpaa_control(
            ideal_control_state,
            np.sqrt(probabilities[index]),
            np.sqrt(probabilities[index + 1]),
            alpha,
            beta,
        )
    _require(
        max_phase_bound_difference < 2.0e-13,
        "independent phase-sensitive compiler audit failed",
    )

    for row in end_rows:
        epsilon = float(row["requested_epsilon"])
        joint = float(row["compiled_joint_trace_distance"])
        reduced = float(row["compiled_reduced_trace_distance"])
        total_variation = float(row["compiled_total_variation_distance"])
        ideal_joint = float(row["ideal_joint_trace_distance"])
        compiled_to_ideal = float(row["compiled_to_ideal_joint_trace_distance"])
        _require(total_variation <= reduced + 1.0e-10, "TVD contractivity failed")
        _require(reduced <= joint + 1.0e-10, "trace-distance contractivity failed")
        _require(
            max(joint, reduced, total_variation) <= epsilon + 1.0e-10,
            "endpoint budget failed",
        )
        _require(
            ideal_joint <= epsilon / 2.0 + 1.0e-10,
            "ideal FPAA control budget failed",
        )
        _require(
            compiled_to_ideal <= epsilon / 2.0 + 1.0e-10,
            "compiled FPAA control budget failed",
        )
        _require(
            joint <= ideal_joint + compiled_to_ideal + 1.0e-10,
            "ideal/compiled control triangle failed",
        )
        _require(
            int(row["total_transition_oracle_calls"])
            == 2 * int(row["total_controlled_walk_calls"]) + int(row["total_wrapper_oracle_calls"]),
            "endpoint resource identity failed",
        )

    config = load_config(ROOT / "configs" / "ising_2x4.yaml")
    _, chains = build_from_config(config)
    fresh = run_end_to_end(chains, reference_epsilon)
    validate_end_to_end(fresh, float(config["numerics"]["atol"]))
    saved_reference = next(
        row for row in end_rows if float(row["requested_epsilon"]) == reference_epsilon
    )
    metric_pairs = (
        (
            fresh.joint_error,
            float(saved_reference["compiled_joint_trace_distance"]),
        ),
        (
            fresh.reduced_error,
            float(saved_reference["compiled_reduced_trace_distance"]),
        ),
        (
            fresh.total_variation,
            float(saved_reference["compiled_total_variation_distance"]),
        ),
        (
            fresh.ideal_joint_error,
            float(saved_reference["ideal_joint_trace_distance"]),
        ),
        (
            fresh.compiled_to_ideal_joint_distance,
            float(saved_reference["compiled_to_ideal_joint_trace_distance"]),
        ),
        (
            fresh.leakage,
            float(saved_reference["compiled_clean_subspace_leakage"]),
        ),
        (
            fresh.norm_residual,
            float(saved_reference["compiled_norm_residual"]),
        ),
    )
    fresh_difference = max(abs(left - right) for left, right in metric_pairs)
    _require(fresh_difference < 2.0e-12, "fresh end-to-end replay differs from CSV")

    independent_ideal_joint_error = _pure_trace_distance(
        ideal_control_state, np.sqrt(probabilities[-1])
    )
    ideal_joint_state = np.zeros_like(fresh.state)
    ideal_joint_state[0, :, 0] = ideal_control_state
    independent_compiled_to_ideal = _pure_trace_distance(fresh.state, ideal_joint_state)
    ideal_control_difference = max(
        abs(independent_ideal_joint_error - fresh.ideal_joint_error),
        abs(independent_compiled_to_ideal - fresh.compiled_to_ideal_joint_distance),
    )
    _require(
        ideal_control_difference < 2.0e-12,
        "independent ideal-projector control failed",
    )
    resource_totals = {
        "controlled_walk_calls": sum(
            stage.resources.controlled_walk_calls for stage in fresh.stages
        ),
        "wrapper_oracle_calls": sum(stage.resources.wrapper_oracle_calls for stage in fresh.stages),
        "transition_oracle_calls": sum(
            stage.resources.transition_oracle_calls for stage in fresh.stages
        ),
    }
    _require(
        resource_totals
        == {
            "controlled_walk_calls": 10516,
            "wrapper_oracle_calls": 32,
            "transition_oracle_calls": 21064,
        },
        "headline resource totals drifted",
    )

    record: dict[str, Any] = {
        "status": "PASS",
        "scope": "deterministic complex128 spectral reference; not gate-level synthesis",
        "independent_checks": {
            "energy_histogram": histogram,
            "max_spectral_gap_difference": gap_difference,
            "max_phase_gap_difference": phase_difference,
            "max_overlap_difference": overlap_difference,
            "max_compiler_error_difference": max_compiler_difference,
            "max_compiler_error_to_envelope_ratio": max_compiler_ratio,
            "max_phase_sensitive_hybrid_bound_difference": (max_phase_bound_difference),
            "stage_phase_sensitive_compiler_hybrid_bounds": (independent_phase_bounds),
            "stage_transition_oracle_calls": expected_transition_calls,
        },
        "fresh_reference": {
            "epsilon": reference_epsilon,
            "compiled_joint_trace_distance": fresh.joint_error,
            "compiled_reduced_trace_distance": fresh.reduced_error,
            "compiled_total_variation_distance": fresh.total_variation,
            "ideal_joint_trace_distance": fresh.ideal_joint_error,
            "compiled_to_ideal_joint_trace_distance": (fresh.compiled_to_ideal_joint_distance),
            "compiled_clean_subspace_leakage": fresh.leakage,
            "compiled_norm_residual": fresh.norm_residual,
            **resource_totals,
            "max_difference_from_saved_csv": fresh_difference,
        },
        "ideal_projector_control": {
            "purpose": (
                "diagnostic only; all paper metrics use approximate compiled selective phases"
            ),
            "independent_ideal_joint_trace_distance": (independent_ideal_joint_error),
            "independent_compiled_to_ideal_joint_trace_distance": (independent_compiled_to_ideal),
            "max_difference_from_production_control": (ideal_control_difference),
            "total_phase_sensitive_compiler_hybrid_bound": math.fsum(independent_phase_bounds),
            "allocated_total_compiler_budget": reference_epsilon / 2.0,
        },
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": importlib.metadata.version("numpy"),
            "matplotlib": importlib.metadata.version("matplotlib"),
            "pyyaml": importlib.metadata.version("PyYAML"),
            "pytest": importlib.metadata.version("pytest"),
            "ruff": importlib.metadata.version("ruff"),
            "blas": "Apple Accelerate (recorded by numpy.show_config)",
        },
        "artifact_sha256": _artifact_manifest(),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    stable_record = _stable_json(record)
    OUTPUT.write_text(
        json.dumps(stable_record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"PASS: {OUTPUT}")


if __name__ == "__main__":
    main()
