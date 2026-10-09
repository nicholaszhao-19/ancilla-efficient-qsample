# Reproducibility record

## GitHub release environment (2026-10-09)

- macOS 15.7.8, Apple arm64;
- CPython 3.12.14;
- NumPy 1.26.4 linked to OpenBLAS (wheel build: 0.3.23.dev);
- Matplotlib 3.10.3, PyYAML 6.0.3;
- pytest 8.3.5, Ruff 0.15.14;
- setuptools 80.9.0 and wheel 0.45.1 for isolated package builds;
- latexmk 4.86a and pdfTeX 1.40.27 from TeX Live 2025;
- Poppler `pdfinfo` 26.05.0.

`requirements-repro.txt` pins the Python packages used for the recorded run.
The operating system, BLAS implementation, Python minor version, and TeX
distribution remain part of the numerical and typesetting environment and
must be recorded separately; the package does not claim bitwise identity
across different eigensolvers or TeX installations.

## Complete rebuild

From `qsample-validation`, run:

```bash
python scripts/build_publication.py
```

The command promotes warnings to errors, runs all tests and lint checks,
checks installed-package consistency, regenerates every CSV, table, caption,
PNG, and vector PDF, recompiles both schematic figures and the proceedings
paper, rejects LaTeX warnings and page-count drift, and finally executes an
independent numerical audit. It fixes PDF timestamps and trailer identifiers
and uses a temporary Matplotlib cache.

The final audit record is `outputs/data/publication_audit.json`. It contains:

- an independently reconstructed energy histogram, transition spectra,
  adjacent overlaps, compiler errors, every actual-phase hybrid bound, FPAA
  allocations, and resource counts;
- a fresh full-state replay at the headline \(\varepsilon=10^{-3}\);
- an independently implemented ideal-projector control and its full-state
  distance from the approximate compiled trajectory;
- exact environment versions;
- SHA-256 hashes for source, configuration, data, figures, table, and paper.

Floating-point comparisons use stated tolerances. Hashes establish the exact
artifact set for one recorded build; they are not substitutes for the
tolerance-based numerical checks.

For a clean package-only check, build a wheel in a fresh environment with the
two versions declared under `[build-system]` in `pyproject.toml`, then run
`python -m pip wheel . --no-deps`. The audited wheel is a pure-Python
`py3-none-any` artifact.
