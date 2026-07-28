# 1D Heisenberg chain

This folder contains the current research path for the antiferromagnetic spin-1/2 Heisenberg chain,

```math
H = \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+1},
```

with `J = 1` and exact thermodynamic energy density `e0 = 1/4 - ln(2)`.

## Scripts

- `t1_vumps.jl` builds converged two-site `Float64` VUMPS coarse-grainers. Reusable tensors are written to `artifacts/coarse_grainers/heisenberg1d/`.
- `t2t3_compressed_lti.jl` solves the super-spin compressed-LTI relaxation and writes incremental run artifacts under `results/<run-name>/`.
- `test_metadata.jl` checks legacy-result migration and infinite-system scope metadata without running an SDP.

Run these scripts from the repository root:

```bash
julia --project=. scripts/heisenberg1d/t1_vumps.jl
julia --project=. scripts/heisenberg1d/t2t3_compressed_lti.jl
julia --project=. scripts/heisenberg1d/test_metadata.jl
```

Only a clean `OPTIMAL` solver status counts as a certified lower bound. Every certified energy must remain at or below `e0`, and a fixed-`D` ladder must be non-decreasing with hierarchy depth.
