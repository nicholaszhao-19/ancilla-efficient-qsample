# Assumptions and resolved conventions

1. **Hamiltonian normalization.** The corrected manuscript explicitly fixes
   disagreement-count energy. The attached draft configuration's
   `coupling_J: 1.0` would mean a different temperature scale if interpreted
   as conventional \(-J\sum s_us_v\). We retain the manuscript convention and
   record its conventional equivalent \(J=1/2\).
2. **State ordering.** The old formula-only figure left ordering unstated.
   Row-major lexicographic bits are now fixed. Permuting this order leaves all
   reported scalar results invariant.
3. **Oracle extension.** The manuscript specified only the clean column.
   Full-state reference propagation uses the row-wise Householder extension in
   `SPEC.md`.
4. **FPAA convention.** Source and target signs/order were checked against the
   primary Yoder-Low-Chuang definition and by direct two-dimensional tests.
   Production propagation uses the approximate compiled phases; the ideal
   two-dimensional routine is only a sequence/control oracle.
5. **GQSP backend.** The exact intended polynomial is evaluated on the walk
   spectrum. A canonical scalar unitary completion fixes leaked-signal
   dynamics. Rotation-angle synthesis is not claimed; the paper and plots call
   this a deterministic double-precision spectral reference.
6. **Query primitive.** One transition-oracle call is one controlled or
   uncontrolled use of \(O_{\mathcal P}\) or \(O_{\mathcal P}^\dagger\). Tables separately report
   controlled-walk calls and outer wrapper calls.
7. **Numerical domain.** The backend is a complex128 reference for the frozen
   five-chain schedule, not a general arbitrary-gap solver. It rejects
   spectral gaps below \(\epsilon_{\rm mach}^{1/3}\), Gibbs probabilities that
   underflow to zero, and unresolved stationary eigenspaces rather than
   silently extrapolating beyond the validated regime.
