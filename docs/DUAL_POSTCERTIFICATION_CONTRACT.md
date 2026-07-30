# Custom RG/RDM dual replay contract

## Scope

This contract rigorously replays the dual arithmetic of an assembled Float64 compressed-LTI or RDM8 SDP. A replay becomes a physical ground-state certificate only when the same artifact also proves that its exact dyadic coarse-graining coefficients satisfy the recursive physical-witness identities.

The retained implementation rounds the coarse-grainer tensor `B` once to a 13-bit dyadic grid and derives every map and quadratic SDP coefficient from that common tensor. Independently rounded maps are not accepted into the certificate chain.

## Exact dyadic physical assembly

For the accepted D = 4 tensor, every `B` entry is an integer divided by \(2^{13}\). The largest absolute numerator is 6,547. Integer audits bound the `W2` contraction numerator by 171,452,836, the `MR`/`ML` quadratic numerator by 85,726,418, and every quadratic `W2` coefficient numerator by 3,707,608,477,307,858, below the exact `Float64` integer limit \(2^{53}\). The recorded `Float64` contractions therefore equal their exact rational contractions, and `W2`, `MR`, `ML`, the linked RDM maps, and the recursive witness identities all come from one exact tensor.

The replay transcript records this audit and accepts a physical claim only when `dyadic_compatibility_audit.status = "PROVEN"` and `physical_relaxation_compatibility_proven = true`. The proof applies to both the compressed and RDM8 formulations.

## Dual reconstruction

Write the primal SDP as

```math
\min_X \sum_g \langle C_g,X_g\rangle\quad\text{subject to}\quad A(X)=b,\qquad X_g\succeq0.
```

Every scalar equality multiplier returned by Mosek is converted exactly to a rational number. The code then reconstructs

```math
S_g=C_g-A_g^\ast y,\qquad d=b^\mathsf{T}y
```

directly from the JuMP objective and equality coefficients. This avoids assuming that the solver-reported dual residual itself is a rigorous error bound.

## Trace-bounded eigenvalue correction

For each logical PSD variable group, let `T_g` be a proven trace upper bound and let `ell_g` be a rigorous lower enclosure of the smallest eigenvalue over all U(1) blocks in that group. QMBCertify's Arblib helper computes `ell_g` from the exact rational slack matrices. Every primal-feasible point then satisfies

```math
\langle S_g,X_g\rangle\ge \min(0,\ell_g)\operatorname{tr}(X_g)\ge \min(0,\ell_g)T_g.
```

Consequently,

```math
E_{\mathrm{cert}}=d+\sum_g\min(0,\ell_g)T_g
```

is a rigorous lower bound on the SDP optimum and therefore on the physical energy covered by the relaxation contract.

The trace bounds are:

- `tr(rho3) = 1`;
- `tr(rho4) = 1` when the RDM8 lift is enabled;
- `tr(omega4) <= ||W2||_2^2`;
- `tr(omega(m+1)) <= min(||MR||_2^2, ||ML||_2^2) tr(omega(m))`.

The implementation upper-bounds each squared spectral norm by the exact rational induced-norm product `||A||_1 ||A||_infinity`. This is conservative and does not depend on floating-point singular values.

## Published fields

A successful replay records `dual_certificate.status = "CERTIFIED"`, exact rational strings for the reconstructed dual constant, every group correction, and the final endpoint, plus downward-safe floating-point summaries. These fields may be used for a physical lower-bound or finite-gap claim only when `physical_relaxation_compatibility_proven = true`.

Each certificate also records a hashed binary replay transcript containing the assembled primal objective and equalities, the PSD block inventory, all equality multipliers, and the coarse-grainer. Validation independently rebuilds the SDP from the stored coarse-grainer, requires exact equality with the transcript's primal data, reconstructs every slack, reruns the rigorous eigenvalue enclosures, recomputes `min(0, ell_g) T_g`, and compares every exact summary field. A JSON summary without a valid replay transcript is not accepted as a strict certificate.

`residual_adjusted_lower_bound` and `published_lower_bound` are legacy numerical-diagnostic field names retained for historical schema compatibility. They are not rigorous certificates and must not be used for strict claims.

For a comparison of two solves, the raw objective separation must exceed the sum of their independent residual scales. The difference between two certified lower bounds is also recorded, but it is a comparison of certified endpoints rather than a proof that the two exact SDP optima are separated by that amount.

## Validated points

- H7 D = 4, `n_super = 6`: the exact-dyadic RDM8 endpoint exceeds the exact-dyadic compressed endpoint by `2.127e-5`, and physical compatibility is proven.
- H7 finite periodic N = 20, D = 4, `n_super = 10`: the exact-dyadic RDM8 endpoint is `-0.44528862186638035`; together with the exact Rayleigh upper certificate it gives a rigorous finite gap of `6.9295372641e-5`. The six-decimal downward-published endpoint gives the more conservative gap `6.9673506261e-5`.

The reference H7 runs used names `20260729-h7-dyadic-d4-n6-k13` and `20260729-h7-finite-n20-dyadic-k13`. Generated artifacts are not versioned; their runnable entrypoints are `scripts/heisenberg1d/rg_rdm8_local.jl` and `scripts/heisenberg1d/rg_rdm8_finite_n20.jl`.
