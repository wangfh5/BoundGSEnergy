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
    chain (MPSKit), parameterized over bond dimensions. The no-argument run is
    the bounded D = 2, 3 baseline; larger D sweeps require explicit `--Ds=...`.
  - `t2t3_compressed_lti.jl` — the Kull–Schuch compressed-LTI relaxation
    (PRX 14, 021008, Sec. II, Eq. 14) in the super-spin formulation,
    JuMP + Mosek, with dated result folders, coarse-grainer fingerprints, and
    atomic incremental JSON writes.
  - `heisenberg2d_sdp_table8.jl`, `heisenberg2d_sdp_l4_sweep.jl` — the
    precursor reproduction of arXiv:2604.01555 Table 8 (structured NPA,
    QMBCertify.jl): certified lower bounds for the 2D square-lattice
    Heisenberg model, L = 4, 6, 8.
- `results/` — run artifacts (JSON numbers + figures):
  - `20260727-kullschuch-mve/` — compressed-LTI reproduction of
    Kull–Schuch Fig. 2b: bounds below the Bethe e₀ at all n, monotone,
    saturation plateaus 4.69e-3 (D=2) and 2.01e-3 (D=3).
  - `20260728-kullschuch-d-sweep/` — higher-D follow-up: VUMPS converged at
    D = 4, 5, 6, 7; compressed-LTI remained OPTIMAL for D=4 at 8 and 12
    physical spins and for D=5 at 8 physical spins, but runtime/RAM rose
    sharply before the paper's n_eff regime.
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
julia --project=. scripts/t1_vumps.jl                 # bounded D = 2, 3 baseline
julia --project=. scripts/t2t3_compressed_lti.jl      # bounded D = 2, 3 ladder
```

Higher-D probes are explicit:

```bash
julia --project=. scripts/t1_vumps.jl --Ds=4,5,6,7
julia --project=. scripts/t2t3_compressed_lti.jl --Ds=4,5 --ns=4,6 --run-name=20260728-kullschuch-d-sweep
```

Mosek needs an academic license at `~/mosek/mosek.lic` (free:
https://www.mosek.com/products/academic-licenses/).

## Status / roadmap

- [x] Reproduce arXiv:2604.01555 Table 8 with public QMBCertify (finding:
      shipped code reproduces the PRX-2024 "old" column exactly, not the
      2026 "new" column)
- [x] MVE: Kull–Schuch compressed-LTI relaxation, Heisenberg chain,
      Fig. 2b behavior reproduced (D = 2, 3)
- [x] VUMPS coarse-grainers extended to D = 4, 5, 6, 7 (all converged with
      2-site Float64 uMPS; see `scripts/vumps_summary.json`)
- [ ] #49 rung 1: 1D Heisenberg, 200 spins, 1e-5 accuracy (needs D ≳ 7
      and/or coarse-grainer optimization; current local compressed-LTI cost
      wall appears already near D=4, n=6 and D=5, n=4)
- [ ] The combination proper: NPA moment/symmetry constraints inside the
      compressed relaxation
- [ ] Rungs 2–4: J1–J2 chains, 2D lattices
