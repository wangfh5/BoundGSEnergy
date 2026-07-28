# 2D Heisenberg model

This folder contains the structured-NPA baseline for the spin-1/2 square-lattice Heisenberg model with periodic boundaries.

## Scripts

- `heisenberg2d_sdp_table8.jl` reproduces the public QMBCertify Table 8 baseline for `L = 4, 6, 8`.
- `heisenberg2d_sdp_l4_sweep.jl` probes which QMBCertify strengthenings move the `L = 4` bound.

Run scripts from the repository root. The Table 8 reproduction accepts an optional result directory:

```bash
julia --project=. scripts/heisenberg2d/heisenberg2d_sdp_table8.jl
julia --project=. scripts/heisenberg2d/heisenberg2d_sdp_l4_sweep.jl
```

These scripts remain a baseline for the challenge. A certificate-preserving 2D coarse-grained relaxation is a later research stage and must not be confused with truncating an effective Hamiltonian.
