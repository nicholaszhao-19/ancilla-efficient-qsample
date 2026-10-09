#!/usr/bin/env python3
"""Render room-readable validation panels for the QCE26 Beamer deck.

The data and semantics are identical to Fig. 3 in Corrected_Proceedings.pdf.
Only typography, line weight, and the one-panel-per-slide layout change.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "outputs" / "data"
FIGURES = ROOT / "outputs" / "figures"

BLUE = "#1769AA"
TEAL = "#00897B"
ORANGE = "#D97706"
INK = "#102A43"
GRID = "#D9E2EC"


def read_numeric_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def configure() -> None:
    source_sans_dir = Path("/usr/local/texlive/2025/texmf-dist/fonts/opentype/adobe/sourcesanspro")
    if source_sans_dir.is_dir():
        for font_path in source_sans_dir.glob("SourceSansPro-*.otf"):
            font_manager.fontManager.addfont(font_path)
        font_family = "Source Sans Pro"
    else:
        font_family = "DejaVu Sans"
    plt.rcParams.update(
        {
            "font.family": font_family,
            "mathtext.fontset": "stixsans",
            "font.size": 18,
            "axes.labelsize": 20,
            "axes.titlesize": 20,
            "legend.fontsize": 15.5,
            "xtick.labelsize": 17.5,
            "ytick.labelsize": 17.5,
            "axes.linewidth": 1.15,
            "lines.linewidth": 2.6,
            "lines.markersize": 8.5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def finish(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        FIGURES / f"{stem}.pdf",
        bbox_inches="tight",
        pad_inches=0.04,
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(
        FIGURES / f"{stem}.png",
        dpi=240,
        bbox_inches="tight",
        pad_inches=0.04,
    )
    plt.close(fig)


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, which="major", color=GRID, linewidth=0.8, alpha=0.9)
    ax.grid(False, which="minor")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(colors=INK, length=5)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def make_end_to_end() -> None:
    rows = read_numeric_csv(DATA / "corrected_end_to_end_sweep.csv")
    rows.sort(key=lambda row: row["requested_epsilon"])
    epsilon = [row["requested_epsilon"] for row in rows]

    fig, ax = plt.subplots(figsize=(6.6, 3.6), constrained_layout=True)
    ax.loglog(
        epsilon,
        [row["compiled_joint_trace_distance"] for row in rows],
        marker="o",
        markerfacecolor="white",
        markeredgewidth=2.0,
        color=BLUE,
        label=r"Trace distance $D_{\rm tr}$",
        zorder=4,
    )
    ax.loglog(
        epsilon,
        [row["compiled_total_variation_distance"] for row in rows],
        marker="^",
        linestyle="-.",
        color=TEAL,
        label=r"Total variation $d_{\rm TV}$",
        zorder=3,
    )
    ax.loglog(
        epsilon,
        epsilon,
        color=INK,
        linestyle=":",
        linewidth=2.3,
        label=r"Guarantee line $y=\varepsilon$",
        zorder=2,
    )
    ax.set_xlabel(r"Target accumulated error $\varepsilon$")
    ax.set_ylabel("Observed final error")
    ax.legend(
        loc="upper left",
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.94,
        handlelength=2.5,
    )
    style_axis(ax)
    finish(fig, "slide_end_to_end_corrected")


def make_compiler() -> None:
    rows = read_numeric_csv(DATA / "compiler_sweep.csv")
    fig, ax = plt.subplots(figsize=(6.6, 3.6), constrained_layout=True)

    selected_beta = [0.0, 0.6, 1.2]
    colors = [BLUE, ORANGE, TEAL]
    markers = ["o", "s", "^"]
    for beta, color, marker in zip(selected_beta, colors, markers, strict=True):
        series = [row for row in rows if abs(row["beta"] - beta) < 1.0e-12]
        series.sort(key=lambda row: row["requested_eta"])
        ax.loglog(
            [row["requested_eta"] for row in series],
            [row["restricted_operator_error"] for row in series],
            marker=marker,
            markerfacecolor="white" if marker != "^" else color,
            markeredgewidth=1.8,
            color=color,
            label=(
                rf"$\beta={beta:g}$; "
                rf"$\Delta_{{\mathcal{{W}}}}={series[0]['walk_phase_gap']:.3f}$"
            ),
        )

    envelope = [row for row in rows if abs(row["beta"]) < 1.0e-12]
    envelope.sort(key=lambda row: row["requested_eta"])
    ax.loglog(
        [row["requested_eta"] for row in envelope],
        [row["theorem_envelope"] for row in envelope],
        color=INK,
        linestyle="--",
        linewidth=2.3,
        label=r"Analytic envelope $2\varepsilon^{\mathcal{W}}$",
    )
    ax.set_xlabel(r"Filter tolerance $\varepsilon^{\mathcal{W}}$")
    ax.set_ylabel("Restricted operator error")
    ax.invert_xaxis()
    ax.legend(
        loc="lower left",
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.94,
        handlelength=2.3,
    )
    style_axis(ax)
    finish(fig, "slide_compiler_corrected")


def main() -> None:
    configure()
    make_end_to_end()
    make_compiler()
    print(FIGURES / "slide_end_to_end_corrected.pdf")
    print(FIGURES / "slide_compiler_corrected.pdf")


if __name__ == "__main__":
    main()
