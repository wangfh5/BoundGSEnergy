# Structured-RG size campaign

## Scientific boundary

The accepted fused certificate is currently only for periodic N=20. For each matched `(N, n_super)` point, the general size entrypoint `scripts/heisenberg1d/structured_rg_size_screen.jl` solves three routes with the same D=4 dyadic coarse-grainer: RG+RDM8 without NPA, structured NPA without RG, and their fused relaxation. It performs same-process exact replay for the structured and fused routes, but deliberately records `formal_certificate=false`.

The ratified scientific gate is N=30, `n_super=8`. Run it before treating larger sizes as evidence of scalable fusion. If the user explicitly asks to launch N=50/100/200 in parallel before that gate, do so only as an exploratory campaign and preserve the `formal_certificate=false` label.

## Parameter set

Use the accepted model family unchanged: order 4, `extra=1`, RDM9, PSO3, LSO, 13-bit dyadic D=4 map, and the RDM8 link. Bond dimension is fixed at D=4. The small sweep axis is `n_super`, whose physical support is `2*n_super`. Unnormalized runs use equality scales `(1,16,1)`; normalized recovery runs use the exact power-of-two scales `(1,32,8)` to balance the normalized RG columns without introducing depth-dependent row growth.

Scale-stable recovery runs must set `BGE_OMEGA_NORMALIZATION=dyadic_trace`. This applies the exact positive power-of-two substitution `omega_tilde[k]=omega[k]/alpha[k]`, with the analytic trace bound divided by the same `alpha[k]`; it changes neither the feasible set nor the strict correction. Select `BGE_COEFFICIENT_SCALE` at N=20 from the small set `512, 2048, 16384`, then keep the selected value fixed through N=30 and larger gates. Do not weaken RDM9, PSO3, LSO, solver-status requirements, or exact replay to make a large point terminate.

Recommended points:

| Purpose | N | n_super | Notes |
|---|---:|---:|---|
| Numerical scale control | 20 | 5 | Select the coefficient scale using raw-objective agreement and strict replay loss |
| Scientific gate | 30 | 8 | Must have positive strict improvement before a certification-scale claim |
| N=50 fixed-depth control | 50 | 5 | Separates size growth from RG coverage |
| N=50 proportional coverage | 50 | 13 | Main N=50 point |
| N=100 proportional coverage | 100 | 25 | Exploratory |
| N=200 proportional coverage | 200 | 50 | Exploratory; not a formal certificate |

## Initial resource requests

Use separate single-node jobs because the memory needs differ:

| Point | Partition | CPUs / Mosek threads | Memory | Walltime |
|---|---|---:|---:|---:|
| N=20, k=5 control | `64c512g` | 8 | 48G | 02:00:00 |
| N=30, k=8 | `64c512g` | 12 | 48G | 01:00:00 |
| N=50, k=5 or 13 | `64c512g` | 12 | 96G | 02:00:00 |
| N=100, k=25 | `64c512g` | 32 | 256G | 12:00:00 |
| N=200, k=50 | `64c512g` | 64 | 480G | 7-00:00:00 |

These CPU counts keep every request at or below the partition's 8,000 MB-per-CPU ceiling and favor campaign-level parallelism over ineffective inner-solver oversubscription. MOSEK uses `SLURM_CPUS_PER_TASK`, while Julia model construction and exact replay remain mostly serial. Use shared nodes by default; test `--exclusive` or explicit CPU binding only as a measured performance A/B after the numerical Gate passes. After completion, use `sacct MaxRSS` to right-size later runs. Move N=200 to `huge` only after N=100 passes and a documented OOM or measured projection exceeds one 512 GB node.

## Submission workflow

Create one fresh campaign timestamp locally:

```bash
campaign="$(date +%Y%m%d-%H%M%S)"
```

For each point, construct an exact command of this form:

```bash
run-remote -H siyuan sbatch --test-only \
  --partition=64c512g --qos=normal --nodes=1 --ntasks=1 \
  --cpus-per-task=<cpus> --mem=<memory> --time=<walltime> \
  --job-name=bge-n<N>-k<K> --output=slurm-%x-%j.out \
  --export=ALL,BGE_SYSTEM_SIZE=<N>,BGE_N_SUPER=<K>,BGE_RUN_NAME=<campaign>-n<N>-k<K>,BGE_COEFFICIENT_SCALE=<scale>,BGE_OMEGA_NORMALIZATION=dyadic_trace,BGE_SOLVER_LOG=1,BGE_FEASIBILITY_TOLERANCE=1e-7 \
  .agents/skills/boundgsenergy-hpc-operations/scripts/structured-rg-size-screen.sbatch
```

Run every `--test-only` request first. Show the exact resource table and scheduler response. After explicit authorization, remove `--test-only` and submit. Capture each parsable job id.

Do not combine these points into a Slurm array: a uniform 480 GB request would waste resources for smaller systems, while a smaller request would make N=200 unsafe.

## Acceptance checks

For each completed job:

1. `sacct` state is `COMPLETED`, exit code is zero, and `MaxRSS` is recorded.
2. `results/<run>/size_screen.json` exists and parses.
3. All three solves report `OPTIMAL`.
4. `raw_monotonicity_passed`, `fusion_dominates_rg_only`, `all_raw_objectives_below_finite_reference`, `local_rdm_exactness_passed`, and `physical_relaxation_compatibility_proven` are true.
5. `omega_normalization_exact` is true, and the stored omega scales equal the deterministic minimal dyadic scales.
6. Across scale-equivalent runs at the same `(N, n_super)`, define the robust strict gain as the hybrid strict endpoint minus the highest baseline strict endpoint among those runs. Require this robust gain to be positive; never accept a positive paired gain produced by degrading its own baseline replay.
7. Report raw and robust strict improvement, PSD repair shifts, and raw-versus-repaired moment residuals separately.
8. A negative strict improvement is a valid scientific outcome, not an operational failure.
9. Keep `formal_certificate=false` in every report.

After pulling results, summarize them:

```bash
python3 .agents/skills/boundgsenergy-hpc-operations/scripts/summarize_size_screen.py \
  results/<run1>/size_screen.json \
  results/<run2>/size_screen.json
```

Do not start a formal N=200 certificate claim until the size-specific model construction, independently replayable transcript, and finite-size upper certificate have been generalized beyond N=20.
