# BoundGSEnergy

Certified ground-state energy bounds for quantum spin systems via the
**NPA SDP hierarchy + renormalization-group coarse-graining**.

Team **BoundGSEnergy** (Fo-Hong Wang) — entry for challenge
[**QuantumBFS/quantum.harness#49**](https://github.com/QuantumBFS/quantum.harness/issues/49),
registered as [PR #219](https://github.com/QuantumBFS/quantum.harness/pull/219)
(Harnessing Quantum 2026).

## What is here

- `scripts/` — all runnable code:
  - `t1_vumps.jl` — VUMPS uniform-MPS coarse-grainer for the 1D Heisenberg
    chain (MPSKit), D = 2, 3.
  - `t2t3_compressed_lti.jl` — the Kull–Schuch compressed-LTI relaxation
    (PRX 14, 021008, Sec. II, Eq. 14) in the super-spin formulation,
    JuMP + Mosek. Reproduces the Fig. 2b behavior.
  - `heisenberg2d_sdp_table8.jl`, `heisenberg2d_sdp_l4_sweep.jl` — the
    precursor reproduction of arXiv:2604.01555 Table 8 (structured NPA,
    QMBCertify.jl): certified lower bounds for the 2D square-lattice
    Heisenberg model, L = 4, 6, 8.
- `results/` — run artifacts (JSON numbers + figures):
  - `20260727-kullschuch-mve/` — compressed-LTI reproduction of
    Kull–Schuch Fig. 2b: bounds below the Bethe e₀ at all n, monotone,
    saturation plateaus 4.69e-3 (D=2) and 2.01e-3 (D=3).
  - `20260727-183425-wang2026-table8-heisenberg2d/` — Table 8 run
    (`run.json` + self-contained `report.html`).
- `knowledge/polynomial-optimization/` — a 34-reference survey library
  (rendered full text + INDEX.md + NOTES.md field map) built with the
  harness's `/survey` pipeline. `knowledge/ref.bib` is the bibliography
  source of truth.
- `docs/PLAN.md` — the ratified plan: what we compute, sizes, success
  criteria, kill criteria.
- `Ion.toml` / `Ion.lock` — the research skills this repo's agents use,
  pinned as remote dependencies on `QuantumBFS/quantum.harness` and
  `QuantumBFS/sci-brain` (install with `ion add`).

## Reproduce

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # pinned Manifest
julia --project=. scripts/t1_vumps.jl                 # coarse-grainer MPS
julia --project=. scripts/t2t3_compressed_lti.jl      # compressed-LTI ladder
```

Mosek needs an academic license at `~/mosek/mosek.lic` (free:
https://www.mosek.com/products/academic-licenses/).

## Status / roadmap

- [x] Reproduce arXiv:2604.01555 Table 8 with public QMBCertify (finding:
      shipped code reproduces the PRX-2024 "old" column exactly, not the
      2026 "new" column)
- [x] MVE: Kull–Schuch compressed-LTI relaxation, Heisenberg chain,
      Fig. 2b behavior reproduced (D = 2, 3)
- [ ] #49 rung 1: 1D Heisenberg, 200 spins, 1e-5 accuracy (needs D ≳ 7
      and/or coarse-grainer optimization)
- [ ] The combination proper: NPA moment/symmetry constraints inside the
      compressed relaxation
- [ ] Rungs 2–4: J1–J2 chains, 2D lattices
