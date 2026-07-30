# Research plan — next gates

Updated 2026-07-30. The completed evidence and current numerical values are maintained in [PROGRESS.md](PROGRESS.md); [PLAN.md](PLAN.md) remains the ratified challenge scope.

## Immediate objective

Measure whether the strict fusion gain survives increasing system size when the RG support grows with the system. Keep the accepted order 4, `extra = 1`, RDM9, PSD state optimality level 3, linear state optimality, 13-bit `D = 4` map, RDM8 link, and exact equality scales `(1, 16, 1)` fixed. Run `N = 30`, `n_super = 8` before any `N = 50` or `N = 200` calculation.

## Priority 1 — two-axis size trend

The accepted `N = 20`, `n_super = 5` fused certificate is `-0.44522454084723984`, improves the matched structured-NPA certificate by `2.715797138952425e-6`, and has strict gap `5.214353500759827e-6`. Complete model rebuild, exact mixed replay, and coefficient/scale tamper rejection pass.

The exploratory `N = 30`, `n_super = 5` point has positive raw gain `2.1437405645086116e-7` but negative strict gain because its moment residual grows. This rules out the assumption that a fixed shallow RG cone automatically becomes more valuable at larger `N`. Test two axes:

- Fixed depth: `(N, n_super) = (20, 5), (30, 5), (40, 5)`.
- Approximately fixed physical coverage: `(30, 8), (40, 10)`.

Record raw and strict improvements, identity/Hamiltonian/other moment residuals, RG slack correction, scalar variables, PSD blocks, wall time, and scheduler-measured maximum RSS. Accept the large-size-benefit hypothesis only if the fixed-coverage strict gain remains positive and grows beyond the combined numerical uncertainty.

## Priority 2 — unified `N = 50` gate

After the `N = 30`, `n_super = 8` gate passes, compare `N = 50`, `n_super = 5` with `n_super = 13`. Require `OPTIMAL`, positive strict improvement under same-process exact replay, and residual scaling consistent with a future formal certificate. Do not infer a certificate from the existing `SLOW_PROGRESS` structured-NPA screen or from raw monotonicity alone.

Continue using only charge-preserving, exactly representable coarse-grainers so every recursive map and SDP coefficient retains the physical witness identities.

## Priority 3 — higher-`D` and high-performance compute

The retained `D = 5`, `N = 100`, `n_super = 40` screen is numerical and high-memory. Continue to `N = 200` only after:

- a smaller strict certificate validates the same map family;
- the predicted certified gap can meet `1e-5`;
- a high-memory job has explicit resource and checkpoint settings;
- the N = 200 post-certificate path is implemented, so the result cannot be mistaken for a residual-adjusted diagnostic.

Current variable-count extrapolation is approximately `V_structured(N) = 10283 + 2456.5N`. The fused model adds 17,088 variables at `n_super = 5` and about 1,184 per additional depth. This gives about 150,196 variables for `N = 50, n_super = 5`, 296,701 for `N = 100, n_super = 25`, and 571,951 for `N = 200, n_super = 50`; conservative memory requests are 64 GB, 128 GB, and 256 GB respectively.

## Stop conditions

- Stop the proportional-depth ladder if `N = 30`, `n_super = 8` has nonpositive strict improvement or loses `OPTIMAL` termination.
- Stop a map search if it improves only the compressed objective and not the exact-compatible RDM8 certificate.
- Stop local RG scaling before memory exhaustion; the retained N = 100 observation already establishes the need for a high-memory environment.
- Report a numerical screen as such whenever exact post-certification is absent.
