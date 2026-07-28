# 2026-07-28 higher-D follow-up

This run continues the 2026-07-27 compressed-LTI MVE by extending the uniform-MPS coarse-grainer bond dimension beyond `D = 3`.

## VUMPS coarse-grainers

All new coarse-grainers converged with the existing 2-site `Float64` setup (`tol = 1e-12`, `maxiter = 400`):

- `D = 4`: `E/N = -0.4410581979`, variational gap `2.09e-3`, `||B|| = 9.3e-13`
- `D = 5`: `E/N = -0.4420478402`, variational gap `1.10e-3`, `||B|| = 8.2e-13`
- `D = 6`: `E/N = -0.4424950823`, variational gap `6.52e-4`, `||B|| = 6.8e-13`
- `D = 7`: `E/N = -0.4426482472`, variational gap `4.99e-4`, `||B|| = 7.5e-13`

The cached tensors live in `scripts/vumps_mps.jls`, with a machine-readable summary in `scripts/vumps_summary.json`.

## Compressed-LTI cost probe

Certified points obtained in this folder:

- `D = 4`, `n = 4` super-spins (`8` physical spins): `E = -0.4563866762`, gap `1.324e-2`, `OPTIMAL`, `16.2 s`
- `D = 4`, `n = 6` super-spins (`12` physical spins): `E = -0.4489520041`, gap `5.805e-3`, `OPTIMAL`, `211.9 s`
- `D = 5`, `n = 4` super-spins (`8` physical spins): `E = -0.4563866761`, gap `1.324e-2`, `OPTIMAL`, `162.8 s`

All certified points remain below the Bethe ansatz reference `e0 = 1/4 - ln 2`.

## Interpretation

The bottleneck in this repo is no longer obtaining converged higher-`D` uMPS coarse-grainers. The local blocker is the compressed SDP cost: moving from `D = 4, n = 4` to `D = 4, n = 6` already raises wall time by about `13x`, and `D = 5, n = 4` is already slower than `D = 4, n = 6`. On this machine, that makes the paper's `n_eff(7) ~ 60` regime inaccessible without a better relaxation.

The most plausible next lever remains the one highlighted by Kull et al. in PRX 14, 021008 Sec. IV and Appendix C: optimize the coarse-graining maps and/or formulate symmetry-aware compressed variables, rather than only increasing `D` inside the current super-spin SDP.
