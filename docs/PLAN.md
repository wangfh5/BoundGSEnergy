# PLAN — BoundGSEnergy / challenge #49

Ratified 2026-07-27. The concise stage result is in the [Challenge Report](challenge-report/report.html), current evidence is maintained in [PROGRESS.md](PROGRESS.md), and the next decision gates are maintained in [RESEARCH_PLAN.md](RESEARCH_PLAN.md).

## 1. What we compute

Implement the NPA + coarse-graining combination: block neighbouring spins into super-sites following PRX 14, 021008, Sec. III-D-2, then add structured-SDP constraints and strict post-certification. Start with certified lower bounds for the periodic 1D spin-1/2 Heisenberg chain.

## 2. Sizes

The ladder is `N = 20 → 50 → 100 → 200`, where `N = 200` is the first challenge target. Each rung must validate accuracy, certification, and cost before the next is attempted.

The finite-size contract proves that an exactly compatible compressed point at depth `n_super` lower-bounds every even periodic chain with `N >= 2n_super`. `support_sites = 2n_super` and the claimed `system_size = N` remain separate metadata fields, and accuracy is always measured against the reference for the claimed finite system.

## 3. Success criteria

- **Floor:** reproduce known `N = 20–50` results and demonstrate a measurable accuracy or cost benefit from coarse-graining.
- **Target:** certify the periodic `N = 200` energy density with gap at most `1e-5`.
- **Stretch:** start the 1D `J1–J2` chain at `N = 100` with gap at most `1e-3`.

## 4. Stop conditions

- Stop a method branch if coarse-graining shows no accuracy or cost advantage at `N = 50`.
- Stop local scaling before memory exhaustion and move to a high-memory environment only after a smaller certified gate predicts that the target is reachable.
- Never promote a residual-adjusted numerical diagnostic to a certified lower bound without exact post-certification.
