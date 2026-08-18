# Assumption Verification Checklist

This checklist enumerates the theoretical assumptions required for the main result (Theorem 1) and the supporting lemmas in the paper. Before running the numerical experiments or applying the theoretical result, verify that each condition below is satisfied for the chosen specification of the model primitives.
$\sigma$
---

## Notation

| Symbol | Description | Reference |
|--------|-------------|-----------|
| \(b: \mathbb{R}^d \to \mathbb{R}^d\) | drift coefficient | Equation (1), Section 1 |
| \(\sigma: \mathbb{R}^d \to \mathbb{R}^{d \times m}\) | diffusion coefficient | Equation (1), Section 1 |
| \(K(u) = u^{H-1/2}/\Gamma(H+1/2)\) | fractional Volterra kernel | Equation (1), Section 1 |
| \(H\) | Hurst parameter | Equation (1), Section 1 |
| \(X_0\) | initial condition | Equation (1), Section 1 |
| \(\theta_t\) | sampling intensity process | Equation (3), Section 1 |
| \(\tau_k^n\) | irregular sampling times | Equation (3), Section 1 |
| \(\eta_n(t)\) | sampling map: \(\eta_n(t) = \tau_k^n\) for \(\tau_k^n \le t < \tau_{k+1}^n\) | Section 2 |
| \(W\) | \(m\)-dimensional Brownian motion driving the SVE | Equation (1), Section 1 |
| \(B\) | \(m^2\)-dimensional Brownian motion on an extended space, independent of \(W\) | Theorem 1 |
| \(X^n_t\) | Euler scheme solution | Equation (2), Section 1 |
| \(U^n_t = n^H (X_t - X^n_t)\) | normalized error process | Theorem 1 |
| \(U_t\) | limiting error process | Theorem 1 |
| \(V^{n,k,j}_t\) | auxiliary martingale process | Section 4 |
| \(\Delta^n_t\) | remainder term in the decomposition of \(U^n\) | Section 4 |
| \(\Phi^n_t\), \(\Phi_t\) | auxiliary limiting terms in Lemma rev3 | Section 4 |

---

## Assumption Summary

| # | Object | Required Condition | Status | Notes |
|---|--------|---------------------|--------|-------|
| A1 | \(b\) | continuous on \(\mathbb{R}^d\) | ☐ | |
| A2 | \(b\) | bounded first derivatives: \(\sup_{x \in \mathbb{R}^d} \| \nabla b(x) \| < \infty\) | ☐ | |
| A3 | \(b\) | uniformly continuous derivatives | ☐ | |
| A4 | \(\sigma\) | continuous on \(\mathbb{R}^d\) | ☐ | |
| A5 | \(\sigma\) | bounded first derivatives: \(\sup_{x \in \mathbb{R}^d} \| \nabla \sigma(x) \| < \infty\) | ☐ | |
| A6 | \(\sigma\) | uniformly continuous derivatives | ☐ | |
| A7 | \(H\) | \(H \in (0, 1/2]\) | ☐ | See equation (1); the case \(H = 1/2\) reduces to classical SDE |
| A8 | \(\theta_t\) | strictly positive: \(\theta_t > 0\) a.s. for all \(t\) | ☐ | |
| A9 | \(\theta_t\) | adapted to the filtration \(\{\mathcal{F}_t\}\) | ☐ | |
| A10 | \(\theta_t\) | càdlàg in \(t\) a.s. | ☐ | |
| A11 | \(1/\theta_t\) | bounded: there exists \(K_0 < \infty\) s.t. \(|1/\theta_t| \le K_0\) a.s. for all \(t \in [0,T]\) | ☐ | Used in Lemmas 4.2–4.7; can be relaxed via localization (see Remark after Theorem 1) |
| A12 | \(1/\theta_t\) | a.s. Riemann integrable on \([0,T]\) | ☐ | Assumption in Theorem 1 |
| A13 | \(X_0\) | \(\mathbb{E}[|X_0|^p] < \infty\) for all \(p \ge 1\) (or as needed for moment estimates) | ☐ | Used in Lemma 2.2 (moment bounds) |
| A14 | Sampling grid | \(\tau_0^n = 0\), \(\tau_{k+1}^n = \tau_k^n + 1/(n\theta(\tau_k^n))\) | ☐ | Defined in equation (3) |
| A15 | \(\eta_n(t)\) | \(\eta_n(t) = \tau_k^n\) for \(\tau_k^n \le t < \tau_{k+1}^n\) | ☐ | Sampling map used in the Euler scheme (Section 2) |

---

## Additional Conditions for Specific Lemmas

Some lemmas require stronger or additional assumptions beyond the main theorem. These are listed below for completeness.

| # | Lemma / Result | Additional Condition | Notes |
|---|----------------|----------------------|-------|
| L1 | Lemma 2.2 (moment bounds for \(X\) and \(X^n\)) | \(b(0), \sigma(0)\) finite; \(p \ge 1\) | Follows from A1–A6 and A13 |
| L2 | Lemma 2.3 (Hölder continuity of \(X\) and \(X^n\)) | \(p > H^{-1}\) | Required for Kolmogorov continuity criterion |
| L3 | Lemma 2.4 (convergence of the irregular grid integral) | \(H\) continuous adapted; \(H^n \to H\) in \(L^2\) | New lemma replacing Fukasawa–Ugai Lemma C.2 |
| L4 | Lemma 4.1 (tightness of \(U^n\)) | All assumptions A1–A12 | See Section 4 proof |
| L5 | Lemma 4.7 (pathwise uniqueness for the limit equation) | Lipschitz continuity of \(b, \sigma\) (follows from A2, A5) | The \(\theta\)-term cancels in the difference argument |

---

## Verification Status Summary

| Category | Status |
|----------|--------|
| Coefficient assumptions (A1–A6) | ☐ All verified |
| Hurst parameter (A7) | ☐ Verified |
| Sampling intensity (A8–A12) | ☐ All verified |
| Initial condition (A13) | ☐ Verified |
| Grid definition (A14–A15) | ☐ Verified |

---

## Instructions for Use

1. For each assumption A1–A15, verify analytically or numerically that the condition holds for your specific choice of \(b\), \(\sigma\), \(\theta\), and \(H\).
2. Mark the corresponding checkbox as ☐ Pass or ☐ Fail.
3. If any assumption fails, the theoretical result (Theorem 1) may not apply; consider modifying the specification or using a localization argument as mentioned in the remark following Theorem 1.
4. For the additional lemma-specific conditions, confirm that the required bounds hold in the verification of each lemma.

---

*Last updated: August 2026*
