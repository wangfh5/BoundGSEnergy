# Full-model SOHS post-certification contract

## Scope

`h11_full_sohs_certificate.jl` certifies the complete finite-chain structured-NPA model assembled by the pinned QMBCertify source, including its base Gram, PSD state-optimality, local-RDM positivity, and linear state-optimality terms. This route does not use RG coarse-graining. It is an independent parent-method certificate and an N = 20 accuracy gate, not the N = 200 challenge result.

For the periodic spin-1/2 Heisenberg chain, QMBCertify assembles the coefficient identity

\[
H_N/N-\ell I=\sum_g v_g^\dagger Q_gv_g+\mathcal P+\mathcal R+\mathcal L,
\]

where every \(Q_g\) is positive semidefinite, \(\mathcal P\) and \(\mathcal R\) are the PSD state-optimality and local-RDM dual terms, and \(\mathcal L\) is a linear state-optimality term whose expectation vanishes in a ground state. Taking the expectation in a translation-twirled ground state therefore gives \(e_{\mathrm{PBC}}(N)\ge\ell\).

## Exact replay

The run-time instrumentation gives the pinned public entrypoint a private name, adds a model-only return before optimization, returns the JuMP model after solving, and multiplies the coefficient-matching equations by a recorded positive power of two. The reproduction entrypoint fixes and records the random seed used by QMBCertify's algebraic independence filter. Both `bound_gsp.jl` and the complete QMBCertify Julia source tree are pinned by hard-coded SHA-256 values. The scale, complete variable vector, objective, coefficient equations, PSD block inventory, and canonical Pauli words are stored in a hashed replay transcript.

Validation does not trust that transcript as the model definition. It reruns the pinned package with the recorded full configuration and seed in model-only mode, then requires exact equality of the source hashes, objective expression and maximization sense, complete coefficient system, variable count and indexing, PSD cone inventory, and canonical Pauli-word support before inserting the stored solver values. The model is rejected unless JuMP contains exactly one vector of coefficient-matching equalities and PSD matrix cones; an unsupported constraint type cannot be silently omitted. The status gate separately requires both `OPTIMAL` and `FEASIBLE_POINT`.

Every stored `Float64` coefficient and solver value is interpreted as its exact dyadic rational. Approximate Fourier coefficients remain valid because an arbitrary exact linear combination of Pauli words is still an admissible SOHS polynomial; their positivity does not require the floating-point transformation to be exactly unitary. Scaling all coefficient equations by \(2^{14}=16384\) leaves the feasible set unchanged, and division by the same power of two is exact in binary floating point before rational conversion.

For each exact rational PSD block \(Q_g\), Arblib encloses its minimum eigenvalue from below. If the enclosure is \(\lambda_g<0\), replay replaces \(Q_g\) by \(Q_g-\lambda_gI\). It then recomputes every coefficient equation exactly. If the residual coefficient of canonical Pauli word \(P_w\) is \(r_w\), then \(|\langle P_w\rangle|\le1\), so

\[
\left|\left\langle\sum_wr_wP_w\right\rangle\right|\le\sum_w|r_w|.
\]

The reported endpoint

\[
\ell_{\mathrm{cert}}=\ell_{\mathrm{exact}}-\sum_w|r_w|
\]

is therefore a rigorous lower bound. After rebuilding and matching the physical model, the validator reloads the hashed transcript, reruns all eigenvalue enclosures and exact coefficient arithmetic, and rejects stale summaries, altered values, or a transcript whose model structure was changed even when its file digest was recomputed.

## N = 20 target-gate result

The accepted reference run used the name `20260729-repro-julia115-full-sohs-n20`. Generated artifacts are not versioned. The certificate generator requires and records Julia `1.11.5`. MOSEK terminates `OPTIMAL` with primal status `FEASIBLE_POINT` and recorded objective sense `MAX_SENSE` for order 4, `extra = 1`, RDM9, PSD state optimality level 3, and linear state optimality enabled. The model has 59,413 scalar variables, 3,432 coefficient equations, 49 PSD blocks, and largest PSD dimension 126. All exact rational PSD blocks already have nonnegative rigorous eigenvalue enclosures, so no diagonal repair is needed.

The numeric objective is `-0.4452254675515185`. Exact replay subtracts an identity-residual \(\ell_1\) correction of `2.9462514233e-6` and certifies

\[
e_{\mathrm{PBC}}(20)\ge -0.44522841380294176.
\]

The independently replayed exact dyadic Rayleigh quotient gives

\[
e_{\mathrm{PBC}}(20)\le -0.4452193264937391.
\]

Thus the rigorous finite-energy interval has width `9.0873092027e-6`, passing the N = 20 `1e-5` gate. The downward-safe published lower bound is `-0.445228414`. This result validates the structured-NPA route locally; the challenge target remains N = 200.
