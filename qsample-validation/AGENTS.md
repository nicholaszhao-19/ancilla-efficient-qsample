# Repository rules

- Treat `SPEC.md` as the mathematical source of truth.
- Do not change the energy, state-ordering, oracle-extension, FPAA, or filter
  conventions without updating `SPEC.md`, `ASSUMPTIONS.md`, and the tests.
- Run `python scripts/build_publication.py` for a publication release; focused
  changes must at least pass `PYTHONWARNINGS=error python -m pytest` and
  `ruff check .`.
- Never describe the current backend as a gate-level GQSP simulation. It is a
  deterministic double-precision spectral reference using the intended
  polynomial and a documented canonical unitary completion.
- Save every plotted dataset as CSV and every paper figure as vector PDF.
- Record benchmark commands and residuals in `REPORT.md` and progress in
  `STATUS.md`.
