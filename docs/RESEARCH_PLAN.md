# Research plan — next gates

Updated 2026-07-29. The completed evidence and current numerical values are maintained in [PROGRESS.md](PROGRESS.md); [PLAN.md](PLAN.md) remains the ratified challenge scope.

## Immediate objective

Produce an `OPTIMAL`, exactly post-certified lower bound for periodic `N = 50` whose rigorous gap is close enough to `1e-5` to justify scaling. Do not launch `N = 200` merely because a numerical entrypoint exists.

## Priority 1 — structured-NPA conditioning and constraint selection

The accepted `N = 20` configuration is order 4, `extra = 1`, RDM9, PSD state optimality level 3, and linear state optimality. Its direct `N = 50` screen remains locally tractable but returns `SLOW_PROGRESS` and has a raw reference gap near `4e-5`.

Run cheap `N = 20` screens before every `N = 50` attempt. Accept a candidate only if it remains `OPTIMAL`, improves by more than the combined numerical uncertainty, and does not create a disproportionate variable or cone increase. Exact SOHS replay is required only after the numerical gate passes.

## Priority 2 — exactly representable RG-map improvement

Search only within charge-preserving, exactly representable tensors so that every derived map and SDP coefficient retains the physical witness identities. Score candidates first at local `D = 4`, `n_super = 6`, then at finite `N = 20`, `n_super = 10`. Require improvement in the RDM8 certificate itself, not only in the compressed baseline.

Do not accept a floating-point map whose dual arithmetic replays but whose exact recursive compatibility is unproven.

## Priority 3 — higher-`D` and high-performance compute

The retained `D = 5`, `N = 100`, `n_super = 40` screen is numerical and high-memory. Continue to `N = 200` only after:

- a smaller strict certificate validates the same map family;
- the predicted certified gap can meet `1e-5`;
- a high-memory job has explicit resource and checkpoint settings;
- the N = 200 post-certificate path is implemented, so the result cannot be mistaken for a residual-adjusted diagnostic.

## Stop conditions

- Stop a constraint family if its `N = 20` improvement is below the combined residual scale or it loses `OPTIMAL` termination without a material accuracy gain.
- Stop a map search if it improves only the compressed objective and not the exact-compatible RDM8 certificate.
- Stop local RG scaling before memory exhaustion; the retained N = 100 observation already establishes the need for a high-memory environment.
- Report a numerical screen as such whenever exact post-certification is absent.
