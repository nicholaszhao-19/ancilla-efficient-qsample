# Paper captions

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
