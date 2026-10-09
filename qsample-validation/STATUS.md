# Status

- Complete: mathematical specification, hardened implementation, independent
  reconstruction, generated artifacts, and manuscript integration.
- Benchmark normalization is explicit: disagreement-count cost 1 equals
  conventional Ising \(J=1/2\) up to a constant.
- Acceptance suite: 18 tests pass with warnings treated as errors; Ruff,
  package consistency, independent audit, and the warning-free nine-page
  paper build pass.
- Camera-ready preflight uses the same plain-text title and author in the
  typeset byline and PDF metadata while emitting no active link annotations
  or bookmark outline.
- The notation audit now defines trace distance first for density operators,
  places the clean embedding before both compiler figures, makes the GQSP
  complement identity explicit at \(\mathcal W\), and separates the numerical
  normalized-state and leakage symbols from \(\ket{\psi}\) and the FPAA
  length \(L_i\). The numerical CSV artifacts are unchanged.
- The manuscript uses a filter-degree bound uniform over
  \(\varepsilon^{\mathcal W}\in(0,1)\), explicitly supplies the
  Fejér--Riesz complementary polynomial, and proves \(A_i\geq1\).
- Figure 1 is a vector end-to-end overview adapted from the protected IEEE
  QALG baseline with proceedings-consistent overlap, walk, compiler, and
  common-embedding notation; the baseline itself is unchanged. Figure 2's two
  circuits use the same \(0.84\)-linewidth display width, a
  \(7\) mm inter-panel separation, and left-side `(a)`/`(b)` labels matched by
  the IEEE caption. Its redundant headings and equations are removed. Figure
  3's duplicate bottom sentence is removed from the artwork while its caption
  retains the explanation. The vector PDFs and final page renders are
  verified.
- At \(\varepsilon=10^{-3}\): joint/reduced trace distance
  \(1.082\times10^{-4}\), TVD \(1.889\times10^{-5}\), clean leakage
  \(1.177\times10^{-16}\), and norm residual \(2.29\times10^{-14}\).
- These are approximate compiled-phase FPAA outputs. The ideal-projector
  control is \(1.0816783746\times10^{-4}\), and the compiled-to-ideal
  full-state distance is \(1.0847510\times10^{-8}\).
- Exact phase-sensitive compiler hybrid bounds are
  \(1.66,2.47,2.34,1.14\times10^{-7}\), each below its
  \(1.25\times10^{-4}\) stage allocation; their
  \(7.6122\times10^{-7}\) total is now reported in the manuscript.
- Maximum row-sum, stationarity, detailed-balance, and operator-norm
  projector-intertwining residuals are \(4.44\times10^{-16}\),
  \(1.39\times10^{-17}\), \(1.08\times10^{-19}\), and
  \(1.70\times10^{-14}\), respectively.
- Numerical scope: deterministic complex128 spectral reference with a
  canonical completion; no angle-level, gate-level, or fault-tolerant
  synthesis claim.
- Release command: `python scripts/build_publication.py`.
