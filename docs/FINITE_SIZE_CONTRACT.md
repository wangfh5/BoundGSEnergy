# Finite-size contract for the 1D Heisenberg chain

## Target problem

For an even periodic chain of \(N\) spin-\(1/2\) sites,

\[
H_N=\sum_{i=1}^{N}\mathbf S_i\cdot\mathbf S_{i+1},\qquad \mathbf S_{N+1}\equiv\mathbf S_1,
\]

the target is the energy density

\[
e_{\mathrm{PBC}}(N)=\frac{1}{N}\min_{\rho_N\succeq0,\ \operatorname{tr}\rho_N=1}\operatorname{tr}(H_N\rho_N).
\]

Translation twirling preserves the objective, so the minimization may be restricted to translation-invariant states. The boundary bond \((N,1)\) is therefore included explicitly in \(H_N\), while the translation-reduced objective uses one representative of each bond orbit.

## Compressed-LTI lower bound

Group physical sites \((2j-1,2j)\) into \(M=N/2\) super-sites. The local super-site objective used by `t2t3_compressed_lti.jl` is

\[
h_{\mathrm{super}}=\frac12(h_{12}\otimes I+I\otimes h_{12})+h_{23},
\]

and the script minimizes \(\operatorname{tr}(h_{\mathrm{super}}\rho^{(2)})/2\). Pair-translation invariance makes this exactly the physical energy per site: the first two terms average the intra-block bond and the last term is the inter-block bond, including the periodic boundary bond after translation twirling.

Every translation-invariant periodic state on \(M\) super-sites supplies compatible marginals through depth \(n_{\mathrm{super}}\le M\). Applying the completely positive MPS coarse-graining maps to those marginals produces a feasible point of the compressed relaxation. Consequently,

\[
e_{\mathrm{compressed}}(n_{\mathrm{super}},D)\le e_{\mathrm{LTI}}(n_{\mathrm{super}})\le e_{\mathrm{PBC}}(N)
\]

for every even \(N\ge 2n_{\mathrm{super}}\). The second inequality is the finite-ring relation stated in Kull et al., PRX 14, 021008, Sec. V B: the LTI problem is obtained from any periodic ring of at least the same support by dropping constraints.

No equality between hierarchy support and system size is assumed. A compressed point may support several finite-size claims, but its numerical value remains the same relaxation result.

## U(1)-blocked restriction

The H3 implementation further restricts every local matrix to blocks of fixed total \(2S^z\). This does not remove all physical witnesses needed for the lower-bound argument. Global U(1) twirling, \(\mathcal T(\rho)=\frac{1}{2\pi}\int_0^{2\pi}e^{i\theta S^z_{\mathrm{tot}}}\rho e^{-i\theta S^z_{\mathrm{tot}}}\,d\theta\), preserves positivity, trace, translation invariance, and the Heisenberg objective because \(H_N\) commutes with \(S^z_{\mathrm{tot}}\). The SU(2)-symmetric coarse-grainer is U(1)-covariant; the implementation verifies zero forbidden tensor entries and exact charge preservation for `W2`, `MR`, `ML`, both Kronecker link maps, and `HLOC`. Applying those maps to the twirled marginals therefore lands exactly in the declared charge blocks. Every periodic physical state can thus be replaced by an objective-equivalent witness feasible for the U(1)-blocked formulation, so the blocked optimum remains a lower bound on \(e_{\mathrm{PBC}}(N)\).

## Required metadata

Every finite-size interpretation records:

- `system_scope = "finite_periodic_chain"`
- `system_size = N`
- `support_sites = 2*n_super`
- `n_super`
- `coarse_grainer_fingerprint`
- `status`
- `numerical_bound_schema_version`
- achieved primal/dual residuals and objectives
- raw objectives and residual-adjusted numerical diagnostics
- strict dual-certificate fields when the result is used for a rigorous lower-bound or finite-gap claim

A point can be attached to a finite system only when `system_size` is even and `support_sites <= system_size`. Historical infinite-chain run records retain their original metadata; a finite-size result cites them as provenance or reruns the same parameter point in a fresh result folder.

## Acceptance checks

For the first rung, \(N=20\):

1. Exact diagonalization evaluates the full 20-bond periodic Hamiltonian in the \(S^z_{\mathrm{tot}}=0\) sector.
2. QMBCertify supplies an uncompressed structured-NPA lower bound for the same periodic chain.
3. The compressed ladder uses \(n_{\mathrm{super}}\le10\), terminates `OPTIMAL`, is non-decreasing with hierarchy depth at fixed \(D\), and remains below the exact finite-chain energy.
4. Numeric compressed runs record achieved primal and dual residuals and both solver objectives. Residual-adjusted values are diagnostics only. A strict custom-SDP lower bound must pass [DUAL_POSTCERTIFICATION_CONTRACT.md](DUAL_POSTCERTIFICATION_CONTRACT.md); the QMBCertify comparison uses its rational post-certification step.

## Validated finite ladder

The stored XDiag Lanczos vector has exact dyadic Rayleigh quotient `-0.4452193264937391`, which is a rigorous variational upper bound on the finite-chain ground energy density. H7 replaces the independently rounded recursive maps with a common 13-bit dyadic D = 4 tensor whose entire coefficient assembly is exact in `Float64`. At `n_super = 10`, the compressed endpoint is `-0.4453263421613551` and the RDM8 endpoint is `-0.44528862186638035`, giving a rigorous RDM8 finite gap of `6.9295372641e-5`. The separately reported six-decimal published bound gives the conservative gap `6.9673506261e-5`.

The independent H11 structured-NPA model uses order 4, one extra basis range, RDM9, PSD state optimality level 3, and linear state optimality. Its full-model SOHS replay certifies `-0.44522836660051746`, giving the rigorous interval width `9.0401067783e-6`. This passes the N = 20 `1e-5` validation gate; its proof and transcript contract are recorded in [FULL_SOHS_CERTIFICATE.md](FULL_SOHS_CERTIFICATE.md).

For N = 50 and 100, finite-periodic Bethe equations provide the references. The implementation reproduces the N = 20 XDiag energy to `1.21e-13` before using the same equations at larger N. The OPTIMAL D = 4 endpoint objectives are monotone at fixed coarse-grainer, but their displayed residual-adjusted values remain numerical diagnostics until the dual post-step is run at those sizes.
