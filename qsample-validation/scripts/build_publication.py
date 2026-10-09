#!/usr/bin/env python
"""Run every acceptance check and rebuild every proceedings artifact."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
PAPER = ROOT / "paper"
MAIN_TEX = REPOSITORY / "SAMPLE_Final.tex"
MAIN_PDF = REPOSITORY / "SAMPLE_Final.pdf"


def _run(command: list[str], *, cwd: Path, environment: dict[str, str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, env=environment, check=True)


def _compile_latex(
    source: Path,
    destination: Path,
    *,
    cwd: Path,
    environment: dict[str, str],
    expected_pages: int,
) -> None:
    with tempfile.TemporaryDirectory(prefix="qsample-latex-") as temporary:
        output_directory = Path(temporary)
        _run(
            [
                "latexmk",
                "-pdf",
                "-g",
                "-silent",
                "-interaction=nonstopmode",
                "-halt-on-error",
                f"-outdir={output_directory}",
                source.name,
            ],
            cwd=cwd,
            environment=environment,
        )
        built = output_directory / source.with_suffix(".pdf").name
        if not built.is_file():
            raise RuntimeError(f"latexmk did not create {built.name}")
        log = built.with_suffix(".log").read_text(encoding="utf-8", errors="replace")
        warning_patterns = (
            r"LaTeX Warning:",
            r"Package \S+ Warning:",
            r"Overfull \\[hv]box",
            r"pdfTeX warning",
        )
        if any(re.search(pattern, log) for pattern in warning_patterns):
            raise RuntimeError(f"{source.name}: final LaTeX log contains a warning")
        pdf_information = subprocess.run(
            ["pdfinfo", str(built)],
            cwd=cwd,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        page_match = re.search(r"^Pages:\s+(\d+)$", pdf_information, re.MULTILINE)
        if page_match is None or int(page_match.group(1)) != expected_pages:
            raise RuntimeError(
                f"{source.name}: expected {expected_pages} PDF page(s)"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(built, destination)


def main() -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    environment["SOURCE_DATE_EPOCH"] = "946684800"
    environment["FORCE_SOURCE_DATE"] = "1"
    environment["TZ"] = "UTC"
    environment["PIP_NO_CACHE_DIR"] = "1"
    existing_texinputs = environment.get("TEXINPUTS", "")
    environment["TEXINPUTS"] = (
        str(REPOSITORY) + os.pathsep + existing_texinputs
    )

    with tempfile.TemporaryDirectory(prefix="qsample-mpl-") as matplotlib_cache:
        environment["MPLCONFIGDIR"] = matplotlib_cache
        _run(
            [sys.executable, "-W", "error", "-m", "pytest", "-q"],
            cwd=ROOT,
            environment=environment,
        )
        _run(
            [sys.executable, "-m", "ruff", "check", "."],
            cwd=ROOT,
            environment=environment,
        )
        _run(
            [sys.executable, "-m", "pip", "check"],
            cwd=ROOT,
            environment=environment,
        )
        for script in (
            "run_benchmark.py",
            "run_compiler_sweep.py",
            "run_end_to_end_sweep.py",
            "make_paper_outputs.py",
        ):
            _run(
                [sys.executable, str(ROOT / "scripts" / script)],
                cwd=ROOT,
                environment=environment,
            )

    for name in ("fig0_overview", "fig1_compiler", "fig2_architecture"):
        _compile_latex(
            PAPER / f"{name}.tex",
            PAPER / f"{name}.pdf",
            cwd=PAPER,
            environment=environment,
            expected_pages=1,
        )
    _compile_latex(
        MAIN_TEX,
        MAIN_PDF,
        cwd=REPOSITORY,
        environment=environment,
        expected_pages=9,
    )
    _run(
        [sys.executable, str(ROOT / "scripts" / "run_publication_audit.py")],
        cwd=ROOT,
        environment=environment,
    )
    print(f"PASS: {MAIN_PDF}")


if __name__ == "__main__":
    main()
