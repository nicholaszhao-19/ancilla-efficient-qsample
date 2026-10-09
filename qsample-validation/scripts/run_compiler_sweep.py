#!/usr/bin/env python
"""Generate the exact restricted selective-phase compiler sweep."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from qsample_validation.benchmark import build_from_config, load_config
from qsample_validation.gqsp_filter import design_filter
from qsample_validation.selective_phase import exact_restricted_error
from qsample_validation.szegedy import SpectralWalk
from qsample_validation.validation import validate_chain, validate_compiler_error

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "ising_2x4.yaml"
OUTPUT = ROOT / "outputs" / "data" / "compiler_sweep.csv"


def main() -> None:
    config = load_config(CONFIG)
    _, chains = build_from_config(config)
    atol = float(config["numerics"]["atol"])
    requested_eta = np.logspace(-1, -7, 13)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
    ]
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for chain in chains:
            validate_chain(chain, atol)
            walk = SpectralWalk(chain)
            for eta in requested_eta:
                spec = design_filter(chain.walk_phase_gap, float(eta))
                error = exact_restricted_error(walk, spec, math.pi)
                validate_compiler_error(error, 2.0 * float(eta), atol)
                writer.writerow(
                    {
                        "beta": f"{chain.beta:.17g}",
                        "walk_phase_gap": f"{chain.walk_phase_gap:.17g}",
                        "requested_eta": f"{eta:.17g}",
                        "certified_eta": f"{spec.certified_eta:.17g}",
                        "t": spec.t,
                        "n": spec.n,
                        "degree": spec.degree,
                        "scaled_degree": (
                            f"{spec.degree * chain.walk_phase_gap:.17g}"
                        ),
                        "phase": f"{math.pi:.17g}",
                        "restricted_operator_error": f"{error:.17g}",
                        "theorem_envelope": f"{2.0 * eta:.17g}",
                    }
                )
    print(OUTPUT)


if __name__ == "__main__":
    main()
