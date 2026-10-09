"""Configuration and construction helpers for the fixed paper benchmark."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml

from .ising import IsingModel
from .markov import MarkovChain, build_chain


def validate_config(config: dict[str, Any]) -> None:
    """Reject configuration drift that the frozen paper benchmark cannot honor."""

    expected_values: tuple[tuple[tuple[str, ...], Any], ...] = (
        (("model", "rows"), 2),
        (("model", "cols"), 4),
        (("model", "boundary"), "open"),
        (("model", "energy_convention"), "disagreement_count"),
        (("model", "disagreement_cost"), 1.0),
        (("model", "conventional_ising_equivalent_J"), 0.5),
        (("model", "external_field_h"), 0.0),
        (("model", "spin_values"), [-1, 1]),
        (
            ("model", "state_ordering"),
            "row_major_lexicographic_bits_0_then_1",
        ),
        (("schedule", "beta"), [0.0, 0.3, 0.6, 0.9, 1.2]),
        (("dynamics", "type"), "random_scan_single_site_heat_bath"),
        (("dynamics", "resampling"), "exact_conditional_gibbs"),
        (("dynamics", "laziness"), 0.5),
        (("filter", "family"), "claudon_average_power"),
        (("filter", "phase"), "pi"),
        (("filter", "completion"), "canonical_spectral"),
        (("numerics", "dtype"), "complex128"),
    )
    errors: list[str] = []
    for path, expected in expected_values:
        value: Any = config
        try:
            for key in path:
                value = value[key]
        except (KeyError, TypeError):
            errors.append(f"{'.'.join(path)} is missing")
            continue
        if value != expected:
            errors.append(f"{'.'.join(path)} must be {expected!r}, got {value!r}")

    try:
        numerics = config["numerics"]
        benchmark_epsilon = float(numerics["benchmark_epsilon"])
        sweep_epsilon = [float(value) for value in numerics["sweep_epsilon"]]
        atol = float(numerics["atol"])
        rtol = float(numerics["rtol"])
        if not 0.0 < benchmark_epsilon < 1.0:
            errors.append("numerics.benchmark_epsilon must lie in (0, 1)")
        if not sweep_epsilon or any(
            not math.isfinite(value) or not 0.0 < value < 1.0
            for value in sweep_epsilon
        ):
            errors.append("numerics.sweep_epsilon must contain finite values in (0, 1)")
        if benchmark_epsilon not in sweep_epsilon:
            errors.append("numerics.sweep_epsilon must include benchmark_epsilon")
        if not math.isfinite(atol) or atol <= 0.0:
            errors.append("numerics.atol must be finite and positive")
        if not math.isfinite(rtol) or rtol <= 0.0:
            errors.append("numerics.rtol must be finite and positive")
    except (KeyError, TypeError, ValueError):
        errors.append("numerics tolerances and epsilon values must be present and numeric")

    if errors:
        raise ValueError("invalid frozen benchmark configuration: " + "; ".join(errors))


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError("benchmark configuration must be a mapping")
    validate_config(config)
    return config


def build_from_config(config: dict[str, Any]) -> tuple[IsingModel, list[MarkovChain]]:
    validate_config(config)
    model_config = config["model"]
    model = IsingModel(
        rows=int(model_config["rows"]),
        cols=int(model_config["cols"]),
        disagreement_cost=float(model_config["disagreement_cost"]),
        external_field_h=float(model_config["external_field_h"]),
    )
    laziness = float(config["dynamics"]["laziness"])
    chains = [
        build_chain(model, float(beta), laziness=laziness)
        for beta in config["schedule"]["beta"]
    ]
    return model, chains
