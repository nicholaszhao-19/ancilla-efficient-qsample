# QSample validation

This package reproduces the exact \(2\times4\) Ising benchmark and validates
the corrected common-interface selective-phase compiler and sequential FPAA
theorem used in the accompanying IEEE manuscript.

The backend performs deterministic full-state propagation in `complex128` and
evaluates the explicit GQSP filter on the complete walk spectrum. It is a
**double-precision spectral reference**, not a gate-level simulation: a
documented canonical unitary completion fixes off-clean signal dynamics, while
GQSP rotation angles are not synthesized.

All headline, figure, and table errors come from approximate compiled
selective-phase FPAA. Ideal projectors are retained only as a labelled
numerical control.

See `SPEC.md` for definitions, `IMPLEMENT.md` for commands, `REPORT.md`
for final residuals, and `REPRODUCIBILITY.md` for exact provenance.

Quick start:

```bash
python -m pip install -e ".[test]" -c requirements-repro.txt
python scripts/build_publication.py
```

The build command fails on test warnings, lint violations, broken package
requirements, numerical acceptance failures, LaTeX warnings, or page-count
drift. It writes an independent audit and SHA-256 manifest to
`outputs/data/publication_audit.json`.
