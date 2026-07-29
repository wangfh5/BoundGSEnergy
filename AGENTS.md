# BoundGSEnergy — agent guide

Research repo for challenge QuantumBFS/quantum.harness#49: certified
ground-state energy bounds for spin-1/2 Heisenberg systems via
NPA + RG coarse-graining. Registered as PR #219.

## Layout

- `scripts/` — runnable Julia (see README "Reproduce"). Project env is the
  repo-root `Project.toml`/`Manifest.toml`: run as `julia --project=. ...`
- `results/<run>/` — generated run artifacts; ignored by Git, and new runs get a fresh dated folder rather than overwriting an old one
- `knowledge/` — survey library; `ref.bib` is the bibliography source of
  truth; `.raw/` and `.figures/` are gitignored local caches
- `docs/` — PLAN.md (ratified plan), discussion logs

## Conventions

- Every claimed bound must come from a clean OPTIMAL solver status; quote
  digits at the residual level, never beyond
- The sandbox ceiling: a certified lower bound must sit at or below its
  reference value — a violation means an index/model bug, always
- Commit the code, pinned environment, tests, documentation, and reproduction commands; do not commit generated results or coarse-grainer caches
- Run the 1D certificate pipeline with exactly Julia 1.11.5
- Skills are Ion remote dependencies (`ion add` to sync); do not vendor
  copies of skills into this repo
