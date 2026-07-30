# Structured-NPA + RG fusion contract

## Scope

`structured_rg_local.jl` implements a periodic `N = 20`, `D = 4` validation model containing both the accepted full structured-NPA constraints and the exact-compatible RG hierarchy. Its purpose is to establish a physical, exactly replayable intersection of the two methods before any larger-system run. The finite witness requires `2 * n_super <= N`; the implementation enforces this condition in argument parsing, model construction, result validation, and transcript rebuilding.

## Unified primal model

The structured-NPA configuration is order 4, `extra = 1`, RDM9, PSD state optimality level 3, and linear state optimality. Dualizing its SOHS maximization gives a moment minimization with canonical Pauli moments \(m_w\). QMBCertify scales its coefficient equations by \(2^{14}\), and the implementation records the exact relation between its internal variables \(y_w\) and physical moments:

\[
m_w=-2^{14}y_w.
\]

The accepted link reconstructs the eight-spin density matrix for four two-spin super-sites affinely from the same canonical moments:

\[
\rho_4^{\mathrm{structured}}=2^{-8}\sum_{P\in\{I,X,Y,Z\}^{\otimes8}}m_P P.
\]

An independent U(1)-blocked variable \(\rho_4\succeq0\) is constrained entrywise to equal this structured density matrix. The coefficient audit enumerates every Pauli contribution, proves that each addend and final affine coefficient lies on the declared dyadic grid, and bounds every accumulated integer numerator below \(2^{53}\). The stored Float64 equations therefore equal their exact-rational construction rather than an independently rounded approximation.

The charge-preserving coarse-grainer is quantized to a 13-bit dyadic grid. The first RG variable \(\omega_4\succeq0\) is connected to the independent physical RDM by

\[
\omega_4=(I_4\otimes W_2\otimes I_4)\rho_4(I_4\otimes W_2\otimes I_4)^\mathsf T.
\]

For depths \(k\ge5\), U(1)-blocked variables \(\omega_k\succeq0\) satisfy both recursive marginal equations derived from \(M_R\) and \(M_L\). The existing dyadic audit proves that `W2`, `MR`, `ML`, their Kronecker maps, and the symmetric-variable quadratic coefficients are assembled exactly in Float64. Every physical periodic-chain state therefore supplies feasible structured moments, the physical \(\rho_4\), and all recursive RG witnesses, so the fused minimization remains a lower-bound relaxation of the physical problem.

The accepted unnormalized numerical formulation multiplies the local-RDM, \(\omega_4\)-link, and recursive equality families by the exact factors `(1, 16, 1)`. The normalized formulation uses `(1, 32, 8)`: at the N=20 control point these factors exactly compensate the \((\alpha_4,\alpha_5)=(2,8)\) variable-column scaling, while the fixed factors remain bounded at larger depth. These nonzero power-of-two factors leave the feasible set unchanged and are included in the transcript, deterministic model rebuild, and tamper tests.

For size-stable RDM8 runs, the implementation can additionally use the exact dyadic variable substitution \(\widetilde{\omega}_k=\omega_k/\alpha_k\), where \(\alpha_k\) is the smallest power of two greater than or equal to the analytic trace bound \(T_k\). The first link becomes \(\widetilde{\omega}_4=C\rho_4C^\mathsf{T}/\alpha_4\), and each recursion carries the exact coefficient \(\alpha_{k-1}/\alpha_k\). Every \(\alpha_k\) is positive, so this is a bijective reparameterization of the same PSD feasible set rather than a relaxation change. Exact replay uses the corresponding trace bound \(T_k/\alpha_k\), and transcript version 5 binds the normalization mode, every dyadic scale, and the normalized equality scales. The unnormalized accepted N=20 certificate remains transcript version 4.

## Exact certificate

The fused moment model is dualized back to a SOHS-like maximization. Its transcript records the pinned Julia and QMBCertify versions, objective and `MAX_SENSE`, canonical moment equations, all SOHS PSD blocks, the exact local-RDM coefficient audit, the dyadic coarse-grainer and charge labels, and every affine RG PSD slack.

Validation rebuilds the complete model without solving and requires exact equality of every recorded structure field. Every floating-point solution value is interpreted as an exact dyadic rational. SOHS blocks are repaired with rigorous Arblib minimum-eigenvalue enclosures. The remaining Pauli-word identity residual is subtracted in \(\ell_1\) norm because every Pauli word has operator norm one. For each RG primal group with trace upper bound \(T\) and rigorous dual-slack eigenvalue lower bound \(\lambda\), replay adds the conservative correction

\[
T\min(0,\lambda).
\]

The trace correction covers `hybrid_rho4`, `hybrid_omega4`, and every deeper `hybrid_omega` group. The finite upper endpoint is an independently replayed exact-rational Rayleigh quotient. A result is accepted only when both solves are `OPTIMAL`, both primal statuses are `FEASIBLE_POINT`, the transcript matches the rebuilt model, both coefficient audits pass, and the exact lower endpoint does not cross the exact upper endpoint.

## N = 20 validation

The accepted reference run is generated locally and ignored by Git.

| Route | Scalar variables | Raw lower bound | Strict lower bound | Strict gap |
| --- | ---: | ---: | ---: | ---: |
| Matched structured NPA | 59,413 | `-0.4452253797415892` | `-0.4452272566443788` | `7.930150639712252e-6` |
| Structured NPA + `D = 4`, `n_super = 5` RG/RDM8 | 76,501 | `-0.4452249879529044` | `-0.44522454084723984` | `5.214353500759827e-6` |

Both solves terminate `OPTIMAL` with `FEASIBLE_POINT` primal status. The fused raw objective improves by `3.9178868482814266e-7`, while exact mixed replay improves the matched strict endpoint by `2.715797138952425e-6`. The fused gap is below `1e-5` and is also substantially tighter than the independent exact-compatible RG/RDM8 gap. The complete optional validator independently rebuilds the 76,501-variable model, replays both certificates, and rejects altered coefficients or equality scales even when the transcript hash is recomputed.

An exploratory fixed-depth `N = 30`, `n_super = 5` screen remains `OPTIMAL` and has a positive raw gain of `2.1437405645086116e-7`, but its strict gain is `-1.3332406440951713e-5` because the accumulated moment residual grows. This is not a formal finite-size certificate. The next trend gate is therefore `N = 30`, `n_super = 8`, which tests increasing RG coverage at fixed system size before any `N = 50` or `N = 200` run.
