#!/usr/bin/env python
"""Render the proceedings figure and generated LaTeX benchmark table."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "outputs" / "data"
FIGURES = ROOT / "outputs" / "figures"
PAPER = ROOT / "paper"


def _read_numeric_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def _scientific(value: float, digits: int = 2) -> str:
    if value == 0.0:
        return "$0$"
    exponent = int(np.floor(np.log10(abs(value))))
    coefficient = value / (10.0**exponent)
    return rf"${coefficient:.{digits}f}\!\times\!10^{{{exponent}}}$"


def _resolved_nonnegative(value: float) -> str:
    if value <= 10.0 * np.finfo(float).eps:
        return r"$<2.3\!\times\!10^{-15}$"
    return _scientific(value)


def make_validation_figure(
    *,
    end_to_end_name: str = "end_to_end_sweep.csv",
    output_stem: str = "fig3_validation",
    corrected_proceedings_notation: bool = False,
) -> None:
    compiler = _read_numeric_csv(DATA / "compiler_sweep.csv")
    end_to_end = _read_numeric_csv(DATA / end_to_end_name)
    FIGURES.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.5,
            "legend.fontsize": 6.8,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.2,
            "lines.markersize": 3.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    colors = ["#0072B2", "#D55E00", "#009E73"]
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 2.55))

    selected_beta = [0.0, 0.6, 1.2]
    markers = ["o", "s", "^"]
    for beta, color, marker in zip(selected_beta, colors, markers, strict=True):
        rows = [row for row in compiler if abs(row["beta"] - beta) < 1.0e-12]
        rows.sort(key=lambda row: row["requested_eta"])
        axes[0].loglog(
            [row["requested_eta"] for row in rows],
            [row["restricted_operator_error"] for row in rows],
            marker=marker,
            color=color,
            label=(
                rf"$\beta={beta:g}$, "
                rf"$\Delta_{{\mathcal{{W}}}}={rows[0]['walk_phase_gap']:.3f}$"
                if corrected_proceedings_notation
                else rf"$\beta={beta:g}$, $\Delta={rows[0]['walk_phase_gap']:.3f}$"
            ),
        )
    envelope_rows = [row for row in compiler if abs(row["beta"] - selected_beta[0]) < 1.0e-12]
    envelope_rows.sort(key=lambda row: row["requested_eta"])
    envelope_label = (
        r"$2\varepsilon^{\mathcal{W}}$"
        if corrected_proceedings_notation
        else r"$2\eta$"
    )
    axes[0].loglog(
        [row["requested_eta"] for row in envelope_rows],
        [row["theorem_envelope"] for row in envelope_rows],
        color="#333333",
        linestyle="--",
        label=envelope_label,
    )
    if corrected_proceedings_notation:
        axes[0].set_xlabel(r"compiler error $\varepsilon^{\mathcal{W}}$")
        axes[0].set_ylabel(
            r"$\|\widetilde{\mathcal{S}}_{\pi}(\Pi_0)\mathcal{J}"
            r"-\mathcal{J}\mathcal{S}_\pi(\Pi_{\pi_i})\|$"
        )
    else:
        axes[0].set_xlabel(r"requested filter tolerance $\eta$")
        axes[0].set_ylabel(
            r"$\|\widetilde{\mathcal{S}}_{i,\pi}^{(\eta)}\mathcal{J}"
            r"-\mathcal{J}\mathcal{S}_\pi(\Pi_{\pi_i})\|$"
        )
    axes[0].set_title(
        "(a) Compiler spectral norm"
        if corrected_proceedings_notation
        else "(a) Compiler spectral error"
    )
    axes[0].grid(True, which="major", linewidth=0.4, alpha=0.28)
    axes[0].legend(loc="lower left", frameon=False)
    axes[0].invert_xaxis()

    end_to_end.sort(key=lambda row: row["requested_epsilon"])
    epsilon = [row["requested_epsilon"] for row in end_to_end]
    axes[1].loglog(
        epsilon,
        [row["compiled_joint_trace_distance"] for row in end_to_end],
        marker="o",
        markerfacecolor="none",
        markeredgewidth=1.1,
        color=colors[0],
        label=(
            r"$D_{\rm tr}$"
            if corrected_proceedings_notation
            else r"$D_{\rm tr}$, compiled joint"
        ),
        zorder=4,
    )
    if not corrected_proceedings_notation:
        axes[1].loglog(
            epsilon,
            [row["compiled_reduced_trace_distance"] for row in end_to_end],
            marker="s",
            color=colors[1],
            linestyle="--",
            label=r"$D_{\rm tr}$, compiled reduced $A$",
            zorder=3,
        )
    axes[1].loglog(
        epsilon,
        [row["compiled_total_variation_distance"] for row in end_to_end],
        marker="^",
        color=colors[2],
        linestyle="-.",
        label=(
            r"$d_{\rm TV}$"
            if corrected_proceedings_notation
            else r"$d_{\rm TV}$, compiled"
        ),
    )
    axes[1].loglog(
        epsilon,
        epsilon,
        color="#333333",
        linestyle=":",
        label=r"$y=\varepsilon$",
    )
    if corrected_proceedings_notation:
        axes[1].set_xlabel(r"error accumulation $\varepsilon$")
    else:
        axes[1].set_xlabel(r"requested final error $\varepsilon$")
    axes[1].set_ylabel("observed final error")
    axes[1].set_title(
        "(b) Total error vs. accumulated error"
        if corrected_proceedings_notation
        else "(b) Approximate compiled FPAA"
    )
    axes[1].grid(True, which="major", linewidth=0.4, alpha=0.28)
    axes[1].legend(loc="upper left", frameon=False)

    fig.subplots_adjust(left=0.085, right=0.995, bottom=0.20, top=0.90, wspace=0.31)
    fig.savefig(
        FIGURES / f"{output_stem}.pdf",
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(FIGURES / f"{output_stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_benchmark_table() -> None:
    rows = _read_numeric_csv(DATA / "benchmark_stages.csv")
    PAPER.mkdir(parents=True, exist_ok=True)
    lines = [
        r"\begin{table*}[!t]",
        r"\caption{Fully enumerated open $2\times4$ disagreement-count Ising "
        r"benchmark. "
        r"Panel (a) gives the stage inputs. Panel (b) fixes "
        r"$\varepsilon=10^{-3}$, $\varepsilon_i^{\rm FP}=1.25\times10^{-4}$, "
        r"and $A_i=2.28249$, and reports full-state propagation under the "
        r"approximate compiled phases $\widetilde{\mathcal S}$. The ideal "
        r"projectors are used only for sequence design and an audit control. "
        r"$N_{\rm cW}$ counts controlled "
        r"walk calls, $N_{\rm wrap}$ counts outer $O_{\mathcal P_j}^{(\dagger)}$ calls, "
        r"and $Q_i=2N_{\rm cW}+N_{\rm wrap}$. "
        r"$\widetilde D_i^{\rm local}$ starts "
        r"stage $i$ from a fresh $\mathcal J\lvert\pi_i\rangle$; "
        r"$\widetilde D_i^{\rm cum}$ "
        r"compares the sequential prefix through stage $i$ with "
        r"$\mathcal J\lvert\pi_{i+1}\rangle$.}",
        r"\label{tab:exact-benchmark}",
        r"\centering",
        r"\footnotesize",
        r"\textbf{(a) Stage inputs.}\quad",
        r"\setlength{\tabcolsep}{4.3pt}",
        r"\begin{tabular}{@{}cccccc@{}}",
        r"\toprule",
        r"$i$ & $\beta_i\!\to\!\beta_{i+1}$ & $p_i=q_i$ & "
        r"$(\delta_i,\delta_{i+1})$ & $(\Delta_i,\Delta_{i+1})$ & $L_i$\\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            rf"{int(row['stage'])} & "
            rf"${row['beta_i']:.1f}\!\to\!{row['beta_i_plus_1']:.1f}$ & "
            rf"${row['p_i_exact']:.6f}$ & "
            rf"$({row['delta_i']:.6f},{row['delta_i_plus_1']:.6f})$ & "
            rf"$({row['phase_gap_i']:.6f},{row['phase_gap_i_plus_1']:.6f})$ & "
            rf"${int(row['L_i'])}$\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\\[2pt]",
            r"\textbf{(b) Compiled transport.}\quad",
            r"\setlength{\tabcolsep}{3.4pt}",
            r"\begin{tabular}{@{}ccccccccc@{}}",
            r"\toprule",
            r"$i$ & $(d_{i,i},d_{i,i+1})$ & $\eta_i$ & "
            r"$\widetilde D_i^{\rm local}$ & $\widetilde D_i^{\rm cum}$ & "
            r"$N_{\rm cW}$ & $N_{\rm wrap}$ & $Q_i$ & leakage\\",
            r"\midrule",
        ]
    )
    for row in rows:
        lines.append(
            rf"{int(row['stage'])} & "
            rf"$({int(row['d_i'])},{int(row['d_i_plus_1'])})$ & "
            rf"{_scientific(row['eta_i'])} & "
            rf"{_scientific(row['compiled_local_joint_error'])} & "
            rf"{_scientific(row['compiled_cumulative_joint_error'])} & "
            rf"${int(row['controlled_walk_calls'])}$ & "
            rf"${int(row['wrapper_oracle_calls'])}$ & "
            rf"${int(row['transition_oracle_calls'])}$ & "
            rf"{_resolved_nonnegative(row['compiled_cumulative_clean_subspace_leakage'])}\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table*}",
            "",
        ]
    )
    (PAPER / "table2_benchmark.tex").write_text("\n".join(lines), encoding="utf-8")


def make_captions() -> None:
    text = r"""# Paper captions

## Figure 1

Common-interface architecture of the end-to-end QSample preparation. The
prescribed schedule advances through stationary states with certified adjacent
squared overlap \(q_i\geq p_i\). At stage \(i\), the source- and target-chain
filters for \(j\in\{i,i+1\}\) are synthesized using the supplied certified
phase-gap bounds. The resulting compiled phases form the FPAA stage
\(\widetilde U_i^{\rm FP}\).

## Figure 2

Complete common-interface selective-phase compiler. The transition oracle
lifts the clean embedding into the walk space, the one-signal-qubit GQSP
filter computes the zero-eigenphase flag, \(P_\\phi\) applies the requested
phase, and both computations are unperformed. On the clean column the
restricted error is at most
\(|e^{i\\phi}-1|\\eta\). Open control circles denote control on the signal
state \(\ket0_s\).

## Figure 3

One transport stage in the fixed common embedding. In chronological order,
the Yoder-Low-Chuang FPAA unitary alternates a target-chain phase
\(\\widetilde{\\mathcal S}_{i+1,\\beta_j}\) with a source-chain phase
\(\\widetilde{\\mathcal S}_{i,-\\alpha_j}\). Each compiler begins and ends in the same
\(\mathcal J\)-encoded space.

## Figure 4

Deterministic double-precision validation on the open
\(2\\times4\) disagreement-count Ising ladder. (a) The restricted operator
error of the complete
common-interface compiler at \(\\phi=\\pi\) lies below the \(2\\eta\) theorem
bound for three representative walk gaps. (b) Full-state spectral-reference
propagation through all four approximate compiled-selective-phase FPAA stages;
joint trace distance, the reduced-\(A\) trace distance, and measurement TV
distance all remain below the requested final error. Ideal projectors are used
only for sequence design and a separately labelled control. The filter
polynomial is evaluated on the complete walk spectrum with the canonical
completion specified in the reproducibility package; no gate-level
angle-synthesis claim is made.
"""
    (PAPER / "captions.md").write_text(text, encoding="utf-8")


def main() -> None:
    make_validation_figure()
    make_validation_figure(
        end_to_end_name="corrected_end_to_end_sweep.csv",
        output_stem="fig3_validation_corrected",
        corrected_proceedings_notation=True,
    )
    make_benchmark_table()
    make_captions()
    print(FIGURES / "fig3_validation.pdf")
    print(FIGURES / "fig3_validation.png")
    print(FIGURES / "fig3_validation_corrected.pdf")
    print(FIGURES / "fig3_validation_corrected.png")
    print(PAPER / "table2_benchmark.tex")


if __name__ == "__main__":
    main()
