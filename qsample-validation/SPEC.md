# Mathematical specification

## Benchmark

The benchmark is the convention already stated in the corrected manuscript:
an open \(2\times4\) ladder with ten edges and

\[
H(x)=\sum_{(u,v)\in E}\mathbf 1[x_u\ne x_v].
\]

The coefficient of each disagreement is one. In the conventional
\(-J\sum_{(u,v)}s_us_v\) notation this is \(J=1/2\), up to the additive
constant \(|E|/2\). The configuration file therefore does **not** claim
conventional Ising coupling \(J=1\).

Configurations are lexicographically ordered row-major bit strings, with bit
0 mapped to spin \(-1\) and bit 1 to spin \(+1\). This choice cannot alter
overlaps or spectra, but it fixes every saved vector and matrix entry.

At inverse temperature \(\beta\), random-scan heat-bath dynamics chooses one
site uniformly and exactly resamples it from the conditional Gibbs
distribution. The implemented chain is

\[
\mathcal P_\beta\text{ has transition matrix }\;
P_\beta=\tfrac12(I+G_\beta).
\]

## Walk and common interface

\[
T|x\rangle=|x\rangle\sum_y\sqrt{P(x,y)}|y\rangle,\qquad
W=(2TT^\dagger-I)\,\mathrm{SWAP}.
\]

The clean common embedding is

\[
J|v\rangle=|0\rangle_a|v\rangle_A|0\rangle_B.
\]

For each row \(x\), the transition oracle is the deterministic Householder
extension mapping \(|0\rangle_B\) to
\(\sum_y\sqrt{P(x,y)}|y\rangle\). It is real, Hermitian, and unitary.

## GQSP filter and selective phase

The exact benchmark phase gap is
\(\Delta_{\mathcal W_\beta}=\arccos\lambda_2(P_\beta)\). The symbol
\(\gamma_i\) is reserved for a supplied certified lower bound satisfying
\(0<\gamma_i\leq\Delta_{\mathcal W_i}\). For certified phase-gap lower bound
\(\gamma\) and requested filter tolerance \(\eta\),
use the explicit Claudon average-power polynomial

\[
\Upsilon_{t,n}(z)=
\left(\frac1t\sum_{k=0}^{t-1}z^k\right)^n,\quad
t=\left\lceil\frac{e}{\sin(\gamma/2)}\right\rceil,\quad
n=\left\lceil\ln(1/\eta)\right\rceil .
\]

Its degree is \(d=n(t-1)\), \(\Upsilon(1)=1\), and for
\(\gamma\le|\theta|\le\pi\),

\[
|\Upsilon(e^{i\theta})|
\le [t\sin(\gamma/2)]^{-n}\le e^{-n}\le\eta.
\]

The spectral reference fixes the scalar unitary completion

\[
C(z)=
\begin{pmatrix}
\Upsilon(z)&-\sqrt{1-|\Upsilon(z)|^2}\\
\sqrt{1-|\Upsilon(z)|^2}&\overline{\Upsilon(z)}
\end{pmatrix}.
\]

This completion makes the double-precision spectral-reference propagation
deterministic and unitary, but no claim is made that it is the complementary
polynomial returned by angle-level GQSP synthesis. The top-left polynomial
block is the paper's intended GQSP block, and the theorem's error envelope is
independent of the chosen completion.

The physical selective phase is

\[
\widetilde{\mathcal S}_{i,\phi}^{(\eta)}
=(I_a\otimes O_{\mathcal P_i}^\dagger)
C_i^\dagger(P_\phi\otimes I)C_i
(I_a\otimes O_{\mathcal P_i}).
\]

## FPAA

Use the Yoder-Low-Chuang convention translated to
\(\mathcal S_\phi(\Pi)=I+(e^{i\phi}-1)\Pi\). Chronological action at a stage is

\[
\mathcal S_{\beta_1}(\Pi_t),\ \mathcal S_{-\alpha_1}(\Pi_s),\ldots,
\mathcal S_{\beta_r}(\Pi_t),\ \mathcal S_{-\alpha_r}(\Pi_s).
\]

No other phase ordering is permitted without re-running the two-dimensional
acceptance tests.

For \(\ell\) stages and requested final error \(\varepsilon\), the implementation
uses the theorem allocation

\[
\varepsilon_i=\varepsilon/\ell,\qquad
\varepsilon_i^{\mathrm{FP}}=\varepsilon_i/2,\qquad
A_i=\sum_j\left(
\left|\sin\frac{\alpha_{i,j}}2\right|
+\left|\sin\frac{\beta_{i,j}}2\right|\right),\qquad
\eta_i=\frac{\varepsilon_i}{4A_i}.
\]

Consequently, the phase-weighted compiler increment
\(2A_i\eta_i=\varepsilon_i/2\), rather than the looser
\((L_i-1)\)-weighted allocation.

The production benchmark propagates every compiled FPAA unitary with the
approximate operators \(\widetilde{\mathcal S}_{i,\phi}^{(\eta_i)}\). The ideal
two-dimensional FPAA implementation is used only to design/check the Yoder
phase sequence and as a separately labelled diagnostic control; it never
supplies a reported table or figure error.

## Output metrics

The report computes the following metrics from the approximate compiled
full-state trajectory:

- joint pure-state trace distance to \(J|\pi_\ell\rangle\);
- exact reduced-state trace distance
  \(\frac12\|\rho_A-|\pi_\ell\rangle\langle\pi_\ell|\|_1\);
- total variation distance after measuring \(A\);
- clean-subspace leakage, evaluated without subtractive cancellation as
  \[
  \sum_{a,x,y:\,(a,y)\ne(0,0)}
  |\langle a,x,y|\Psi_\ell\rangle|^2.
  \]

All spectral propagation uses binary64/complex128 arithmetic. The
audit additionally records the ideal-projector endpoint and the
compiled-to-ideal pure-state distance, without treating either as a primary
simulation result. The
nonorthogonal active-space coordinates are accepted only when the enumerated
spectral gap is at least
\(\epsilon_{\rm mach}^{1/3}\); the benchmark minimum
\(4.36\times10^{-3}\) is more than two orders of magnitude larger. This guard
prevents unsupported nearly degenerate chains from returning silent
low-accuracy results.
