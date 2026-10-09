# Implementation runbook

Install the package and validation tools from this directory:

```bash
python -m pip install -e ".[test]" -c requirements-repro.txt
```

Run the complete acceptance and artifact build:

```bash
python scripts/build_publication.py
```

The build uses temporary writable Matplotlib and LaTeX directories,
deterministic PDF timestamps, and regenerates all three schematic PDFs, the
validation figure, and the final nine-page proceedings PDF. It requires
`latexmk`, `pdflatex`, and Poppler `pdfinfo` on `PATH`.

For focused development, the numerical generators remain individually
callable:

```bash
PYTHONPATH=src python scripts/run_benchmark.py
PYTHONPATH=src python scripts/run_compiler_sweep.py
PYTHONPATH=src python scripts/run_end_to_end_sweep.py
MPLCONFIGDIR=/tmp/qsample-mpl PYTHONPATH=src python scripts/make_paper_outputs.py
PYTHONPATH=src python scripts/run_publication_audit.py
```

The scripts produce:

- `outputs/data/benchmark_chains.csv`: chain invariants and residuals;
- `outputs/data/compiler_sweep.csv`: restricted compiler errors;
- `outputs/data/end_to_end_sweep.csv`: explicitly labelled compiled
  joint/reduced/TVD metrics plus the ideal-projector diagnostic control;
- `outputs/data/benchmark_stages.csv`: compiled stage errors, phase-sensitive
  compiler hybrid bounds, allocations, and exact counts;
- `outputs/data/publication_audit.json`: independent checks, environment, and
  artifact hashes;
- `outputs/figures/fig3_validation.{pdf,png}`;
- `paper/table2_benchmark.tex` and `paper/captions.md`.

Stop and document a blocker if stochasticity, stationarity, detailed balance,
unitarity, FPAA ordering, the compiler envelope, or the end-to-end theorem
budget fails its stated tolerance.
