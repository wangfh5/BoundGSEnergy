# 1D Heisenberg scripts

These scripts study the periodic antiferromagnetic spin-1/2 Heisenberg chain

```math
H = \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+1}
```

with `J = 1`. Run every entrypoint from the repository root with the pinned runtime as `julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/<script>.jl`.

## Certified routes

- `t2t3_compressed_lti.jl` contains the compressed-LTI SDP, exact U(1) block construction, RDM lifts, numerical-bound schema, and exact-rational dual post-certification utilities.
- `rdm8_common.jl` contains only shared constants and plotting for the exact-compatible RG/RDM8 runs.
- `rg_rdm8_local.jl` builds one 13-bit dyadic `D = 4` tensor, derives every recursive map and SDP coefficient from it, and certifies the local `n_super = 6` compressed/RDM8 comparison.
- `rg_rdm8_finite_n20.jl` applies the same exact-compatible map at finite periodic `N = 20`, `n_super = 10`, validates its H0/H1/H3 provenance chain, and compares both lower certificates with an exact-rational Rayleigh upper bound.
- `h10_qmbcertify_n20.jl` is the reusable QMBCertify configuration screen used by the full-model route.
- `h11_full_sohs_certificate.jl` solves the selected structured-NPA model and validates its exact SOHS replay against a deterministically rebuilt, source-pinned model.
- `finite_upper_certificate.jl` constructs and validates the exact-rational finite-energy upper certificate used by both routes.

## Baselines and numerical ladder

- `t1_vumps.jl` deterministically builds two-site VUMPS coarse-grainers with a recorded per-dimension seed and writes them under `artifacts/coarse_grainers/heisenberg1d/`.
- `h0_finite_n20.jl` checks the finite-periodic contract against exact diagonalization and a structured-NPA baseline.
- `symmetry_probe.jl` compares dense and exact U(1)-blocked `D = 4` formulations.
- `h3_finite_n20_d4.jl`, `h3_finite_n50_d4.jl`, and `h3_finite_n100_d4.jl` form the validated `D = 4` numerical ladder.
- `h4_symmetry_probe_d5.jl`, `h4_d5_n6.jl`, and `h4_finite_n100_d5.jl` validate the `D = 5` numerical entry gate.
- `h4_finite_n200_d5.jl` is a high-memory numerical entrypoint for depths `50`, `75`, and `100`. It deliberately cannot pass the challenge target until a strict post-certificate is integrated at that scale.

## Retained negative probes

- `h2_hybrid_rdm.jl` shows that RDM8 is inactive at shallow `n_super = 4`; the retained exact-compatible runs show that it becomes useful at active depth.
- `h2_state_optimality.jl`, `h2_psd_optimality.jl`, and `h2_su2_invariance.jl` preserve the three independently valid constraints that showed no residual-resolved tightening at their tested points.
- `h9_optimality_screen_n20.jl` rechecks state-optimality constraints at the active exact-dyadic `N = 20` point and rejects them before strict certification.
- `predecessor_validation.jl` centralizes the H0/H1/H3 provenance checks consumed by later entrypoints.

## Tests

- `test_metadata.jl` checks numerical schema migration without stored results; explicit H0/H1/H3 paths enable full provenance and tamper checks.
- `test_h4.jl` checks higher-`D` ladder invariants without stored results; explicit generated-result paths enable full predecessor replay.
- `test_rg_rdm8.jl` checks exact dyadic compatibility without requiring stored results; optional result paths enable full H3 and certificate replay.
- `test_h11_full_sohs.jl` checks deterministic full-model reconstruction without solving; `--result=<path>` enables exact replay, status gates, and transcript tamper rejection.

Run the complete solver-free suite:

```bash
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_metadata.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h4.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_rg_rdm8.jl
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_h11_full_sohs.jl
```

## Reproduce the RG/RDM8 certificate

The complete `N = 20` chain starts from generated VUMPS tensors and passes every predecessor explicitly:

```bash
run_stamp=$(date +%Y%m%d-%H%M%S)
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/t1_vumps.jl --Ds=3 --seed=20260729
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/h0_finite_n20.jl --run-name="${run_stamp}-h0-n20"
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/symmetry_probe.jl --n-super=6 --mosek-threads=16 --run-name="${run_stamp}-h1-d4-n6"
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/h3_finite_n20_d4.jl --run-name="${run_stamp}-h3-d4-n20" --h0-result="results/${run_stamp}-h0-n20/finite_n20.json" --h1-result="results/${run_stamp}-h1-d4-n6/symmetry_probe.json"
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/rg_rdm8_local.jl --run-name="${run_stamp}-rg-rdm8-local" --mosek-threads=16
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/rg_rdm8_finite_n20.jl --run-name="${run_stamp}-rg-rdm8-n20" --mosek-threads=16 --h3-result="results/${run_stamp}-h3-d4-n20/finite_n20_d4.json"
julia +1.11.5 --startup-file=no --project=. scripts/heisenberg1d/test_rg_rdm8.jl --h3="results/${run_stamp}-h3-d4-n20/finite_n20_d4.json" --local="results/${run_stamp}-rg-rdm8-local/dyadic_d4_n6.json" --finite="results/${run_stamp}-rg-rdm8-n20/finite_n20_rdm8.json"
```

Every result-producing entrypoint refuses an existing result directory. `results/` and `artifacts/coarse_grainers/` are generated locally and ignored by Git. Residual-adjusted quantities are numerical diagnostics, and only runs whose replay validator returns a certified status are described as lower bounds.
