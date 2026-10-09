#!/usr/bin/env python
"""Generate exact chain- and overlap-level benchmark data."""

from __future__ import annotations

import csv
from pathlib import Path

from qsample_validation.benchmark import build_from_config, load_config
from qsample_validation.markov import squared_qsample_overlap
from qsample_validation.projectors import projector_intertwining_operator_norm
from qsample_validation.szegedy import SpectralWalk
from qsample_validation.validation import validate_chain

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "ising_2x4.yaml"
OUTPUT = ROOT / "outputs" / "data" / "benchmark_chains.csv"


def main() -> None:
    config = load_config(CONFIG)
    _, chains = build_from_config(config)
    atol = float(config["numerics"]["atol"])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "beta",
        "spectral_gap",
        "walk_phase_gap",
        "overlap_to_next",
        "row_residual",
        "stationarity_residual",
        "detailed_balance_residual",
        "projector_intertwining_operator_norm",
    ]
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index, chain in enumerate(chains):
            validate_chain(chain, atol)
            projector_residual = projector_intertwining_operator_norm(
                SpectralWalk(chain)
            )
            if projector_residual > atol:
                raise RuntimeError(
                    f"beta={chain.beta}: projector residual exceeds {atol:g}"
                )
            overlap = ""
            if index + 1 < len(chains):
                overlap = squared_qsample_overlap(
                    chain.probabilities, chains[index + 1].probabilities
                )
            writer.writerow(
                {
                    "beta": f"{chain.beta:.17g}",
                    "spectral_gap": f"{chain.spectral_gap:.17g}",
                    "walk_phase_gap": f"{chain.walk_phase_gap:.17g}",
                    "overlap_to_next": (
                        f"{overlap:.17g}" if isinstance(overlap, float) else overlap
                    ),
                    "row_residual": f"{chain.row_residual:.17g}",
                    "stationarity_residual": f"{chain.stationarity_residual:.17g}",
                    "detailed_balance_residual": (
                        f"{chain.detailed_balance_residual:.17g}"
                    ),
                    "projector_intertwining_operator_norm": (
                        f"{projector_residual:.17g}"
                    ),
                }
            )
    print(OUTPUT)


if __name__ == "__main__":
    main()
