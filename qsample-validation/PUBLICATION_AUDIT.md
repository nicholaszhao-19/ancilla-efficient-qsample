# Objective publication audit

## Verdict

**Ready for publication within the stated scope: a deterministic
double-precision spectral-reference validation.** No discrepancy was found in
the benchmark physics, Markov chains, spectra, overlaps, FPAA construction,
filter degrees, phase-sensitive compiler bounds, compiled propagation,
final-state metrics, or resource counts.

This verdict does not extend to angle-level or gate-level GQSP synthesis,
physical runtime, fault-tolerant cost, or asymptotic empirical scaling.

## Material issues found and resolved

1. Small joint errors and leakage were evaluated by cancellation-prone
   subtraction. Stable orthogonal-residual and positive-sum formulas now
   resolve the headline values and restore the expected contractivity order.
2. A seeded random-vector projector probe was described as a residual
   maximum. It is now replaced everywhere by the full operator norm.
3. Output scripts wrote artifacts without enforcing all claimed invariants.
   They now fail before publication output is accepted when residuals,
   compiler envelopes, endpoint budgets, metric ordering, allocations, or
   resource identities fail.
4. Unsupported configuration fields were silently ignored. The complete
   frozen benchmark schema is now validated, and configuration drift is
   rejected.
5. Near-one filter evaluation, extreme \(\eta\), unresolved spectral gaps,
   Gibbs underflow, and arbitrarily scaled states exposed numerical edge
   failures. Stable evaluation and explicit complex128-domain guards now
   prevent silent misreporting.
6. Reproducibility lacked one complete build and machine-verifiable evidence.
   The release command now rebuilds every artifact, fixes PDF timestamps,
   enforces an nine-page warning-free paper, and writes independent numerical
   checks plus an artifact hash manifest. A clean PEP 517 wheel build and
   installation also pass with the pinned build backend.
7. Regression coverage omitted the headline \(\varepsilon=10^{-3}\), stable
   small-error metrics, global completion unitarity, dense walk equivalence,
   configuration drift, and exact resource totals. These cases are now
   covered; all 18 tests pass with warnings promoted to errors.
8. The rounded target error alone could not prove that the production path
   used approximate compiled phases rather than ideal projectors. The
   regression now enforces all 32 target/source compiled calls and their
   phases, the CSV schemas explicitly label compiled metrics, and the
   independent audit records an ideal-projector control plus every
   phase-sensitive hybrid bound.

The fixes change only insignificant numerical digits and the rigor of the
reported projector quantity. The paper-level headline remains
\(1.082\times10^{-4}\) joint/reduced trace distance,
\(1.889\times10^{-5}\) TVD, leakage below \(2.3\times10^{-15}\), and
\(21064\) transition-oracle calls.

At \(\varepsilon=10^{-3}\), the ideal-projector control is
\(1.0816783746\times10^{-4}\), the compiled-to-ideal state distance is
\(1.0847510\times10^{-8}\), and the total exact phase-sensitive hybrid bound is
\(7.6122\times10^{-7}\) against a \(5\times10^{-4}\) allocation. These values
support, but do not replace, the compiled headline results.

## Remaining bounded caveats

- The experiment is one fully enumerated \(2\times4\) instance. It validates
  the implementation and stated example, not empirical scaling.
- Exact enumerated overlaps and phase gaps instantiate the paper's
  \(p_i=q_i\) reference case; a user-supplied promise-bound API is not tested
  here.
- Resource results are exact symbolic query counts for the specified
  primitives, not compiled gate counts or wall-clock measurements.
- Floating-point diagonalization is not interval-certified. Cross-platform
  equality is tolerance-based, although the audited build records exact
  environment and artifact hashes.
- The GitHub release includes paper citation metadata and preserves the
  bundled typesetting licences. No additional reuse licence for the research
  code is granted by this release, and no separate software archival DOI has
  been minted.
