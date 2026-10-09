# Validation report

## Outcome

The fully enumerated open \(2\times4\) disagreement-count Ising benchmark
passes the complete acceptance and independent audit suites. The
implementation propagates the full signal--system--transition state through
the common-interface compiler and all four approximate compiled-phase FPAA
stages using deterministic complex128 spectral functional calculus.

At requested final error \(\varepsilon=10^{-3}\):

- joint trace distance: \(1.0816783801\times10^{-4}\);
- reduced-\(A\) trace distance: \(1.0816783746\times10^{-4}\);
- measurement total-variation distance: \(1.8887983925\times10^{-5}\);
- clean-subspace leakage: \(1.1766847306\times10^{-16}\);
- state-norm residual: \(2.2870594307\times10^{-14}\);
- total controlled-walk calls: \(10516\);
- total outer wrapper calls: \(32\);
- total transition-oracle calls: \(21064\).

The ideal-projector control has joint distance
\(1.0816783746\times10^{-4}\), while the compiled-to-ideal full-state distance
is \(1.0847510\times10^{-8}\). The control is diagnostic only: every headline,
figure, and table error above comes from the compiled trajectory. The exact
phase-sensitive compiler hybrid bounds by stage are
\(1.661,2.472,2.339,1.140\times10^{-7}\), each far below its
\(1.25\times10^{-4}\) allocation.

The stable projection-residual evaluation gives joint distance larger than
reduced distance by \(5.44\times10^{-13}\), as required by contractivity.
The paper rounds both to \(1.082\times10^{-4}\). Leakage is computed by
directly summing positive off-clean probabilities; its exact value is below
the paper's \(2.3\times10^{-15}\) display threshold.

## Chain and walk checks

Across all five temperatures, the maximum residuals are:

| Check | Maximum residual |
|---|---:|
| transition row sum | \(4.44\times10^{-16}\) |
| stationarity | \(1.39\times10^{-17}\) |
| detailed balance | \(1.08\times10^{-19}\) |
| stationary-projector intertwining, operator norm | \(1.70\times10^{-14}\) |

The projector quantity is the full logical-space operator norm, not a seeded
probe-vector residual. The exact stage inputs, filter degrees, allocations,
and resource counts are saved in `outputs/data/benchmark_stages.csv`.

## Independent reconstruction

`scripts/run_publication_audit.py` reconstructs the model without using the
production Ising or Markov-chain builders. It independently checks:

- the energy histogram
  \(\{0:2,2:14,3:24,4:48,5:80,6:48,7:24,8:14,10:2\}\);
- every \(256\times256\) transition spectrum and adjacent Gibbs overlap;
- all 65 compiler-sweep rows and their theorem envelopes;
- every FPAA length, phase weight, filter allocation, degree, and query count;
- every actual Yoder phase's spectrum-restricted hybrid bound;
- an independently implemented ideal-projector trajectory and its distance
  from the compiled full state;
- CSV schemas, row counts, contractivity, and final budgets;
- a fresh full-state replay at \(\varepsilon=10^{-3}\).

Maximum independent differences are \(9.99\times10^{-16}\) for spectral gaps,
\(1.07\times10^{-14}\) for phase gaps, zero at stored precision for overlaps,
and \(1.37\times10^{-14}\) for compiler errors. The phase-sensitive hybrid
bounds agree within \(4.95\times10^{-17}\); the independent ideal controls
agree within \(2.47\times10^{-17}\). The largest saved compiler error is only
\(0.105\) of its analytic envelope. The machine-readable result and artifact
SHA-256 manifest are in
`outputs/data/publication_audit.json`.

## Acceptance suite

- `python -W error -m pytest -q`: 18 passed.
- `python -m ruff check .`: passed.
- `python -m pip check`: no broken requirements.
- package build backend: setuptools 80.9.0 and wheel 0.45.1 (pinned).
- independent publication audit: passed.
- every compiler sweep point satisfies its \(2\eta\) envelope.
- the production-path regression enforces all 32 approximate compiled calls
  in target-\(\beta_j\), source-\(-\alpha_j\) order.
- every phase-sensitive stage hybrid bound satisfies the compiler allocation.
- every joint, reduced, and TVD endpoint satisfies its requested
  \(\varepsilon\).
- the complete nine-page paper builds without LaTeX warnings or page-count
  drift.
- PDF preflight finds all fonts embedded, no encryption, no active link
  annotations, no bookmark outline, and title/author metadata matching the
  visible plain-text title and byline.
- A notation-only camera-ready pass introduced the general density-operator
  trace distance before its pure-state specialization, moved the clean
  embedding ahead of Figures 1--2, corrected the introductory restricted-map
  compiler bound and GQSP complement identity, and removed the numerical
  \(\ket{\psi_i}\)/\(\mathcal L_i\) collisions. Saved numerical CSV values
  were not regenerated or changed.
- The filter-degree asymptotic is uniform on the stated tolerance domain;
  Fejér--Riesz factorization explicitly supplies the complementary
  polynomial, and the asserted \(A_i\geq1\) allocation condition is proved.
- Figure 1 is a vector schedule-to-stage overview adapted from the protected
  IEEE QALG baseline without modifying that baseline. Its notation is
  reconciled with the proceedings definitions: exact squared overlap \(q_i\),
  certified bound \(p_i\), source/target walks and gaps, compiled phases, and
  the common \(J\)-embedding. Figure 2 layout audit: both circuit panels use
  \(0.84\) of the same
  minipage linewidth, have a \(7\) mm separation, and carry unambiguous
  left-side labels `(a)` and `(b)` matched by the IEEE caption. The redundant
  in-figure headings and displayed equations were removed. Figure 3's
  duplicate bottom explanatory sentence was also removed from the artwork
  while its explanation remains in the caption. Both rebuilt vector figures
  and their camera-ready pages were visually checked for alignment, clipping,
  and legibility.

An additional adversarial audit covered 48 random reversible chains,
degenerate spectra, exact \(\lambda=-1\), 360 compiler/operator comparisons,
24 dense completion comparisons, extreme representable filter tolerances,
state scalings from \(10^{-200}\) to \(10^{200}\), and high-temperature
failure guards. No publication-blocking discrepancy remained.

## Numerical scope

The top-left block is the explicit average-power polynomial used in the paper.
The full-state reference fixes the pointwise unitary completion in `SPEC.md`;
it does not synthesize the angle-level complementary GQSP polynomial. The
proved compiler envelope is completion-independent.

The results use floating-point diagonalization, not interval arithmetic, and
exclude finite-precision rotation-synthesis error. Chains below the documented
complex128 spectral-gap threshold, underflowed Gibbs distributions, and
unresolved stationary eigenspaces fail explicitly. The benchmark's minimum
gap \(4.36\times10^{-3}\) is far inside the accepted range.

Numerical acceptance is tolerance-based. Different BLAS/eigensolver builds
can change insignificant last bits; saved endpoint floats are rounded to 10
significant digits, well beyond paper precision, to avoid serializing
unstable noise.

## GitHub release environment (2026-10-09)

- macOS 15.7.8, arm64
- Python 3.12.14
- NumPy 1.26.4 with OpenBLAS (wheel build: 0.3.23.dev)
- Matplotlib 3.10.3
- PyYAML 6.0.3
- pytest 8.3.5
- Ruff 0.15.14
- setuptools 80.9.0 and wheel 0.45.1 (clean build backend)
- latexmk 4.86a
- pdfTeX 1.40.27 (TeX Live 2025)
- Poppler `pdfinfo` 26.05.0

## GitHub release verification (2026-10-09)

`python scripts/build_publication.py` completed successfully in the release
environment above: 18 warnings-as-errors tests, Ruff, package consistency,
all numerical generators, three schematic PDFs, the nine-page manuscript and
the independent numerical audit. Plotting-script import order and one long
legend string were reformatted without changing the numerical algorithms.
The saved headline at requested epsilon 0.001 is joint trace distance
0.000108167838, reduced trace distance 0.0001081678375, TVD
0.00001888798392, leakage 1.176685122e-16, and norm residual
2.420286194e-14, with 21,064 transition-oracle calls.
