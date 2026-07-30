# Structured-RG size campaign

## Scientific boundary

The accepted fused certificate is currently only for periodic N=20. For each matched `(N, n_super)` point, the general size entrypoint `scripts/heisenberg1d/structured_rg_size_screen.jl` solves three routes with the same D=4 dyadic coarse-grainer: RG+RDM8 without NPA, structured NPA without RG, and their fused relaxation. It performs same-process exact replay for the structured and fused routes, but deliberately records `formal_certificate=false`.

The ratified scientific gate is N=30, `n_super=8`. Run it before treating larger sizes as evidence of scalable fusion. If the user explicitly asks to launch N=50/100/200 in parallel before that gate, do so only as an exploratory campaign and preserve the `formal_certificate=false` label.

## Parameter set

Use the accepted model family unchanged: order 4, `extra=1`, RDM9, PSO3, LSO, 13-bit dyadic D=4 map, RDM8 link, and equality scales `(1,16,1)`. Bond dimension is fixed at D=4. The small sweep axis is `n_super`, whose physical support is `2*n_super`.

Recommended points:

| Purpose | N | n_super | Notes |
|---|---:|---:|---|
| Scientific gate | 30 | 8 | Must have positive strict improvement before a certification-scale claim |
| N=50 fixed-depth control | 50 | 5 | Separates size growth from RG coverage |
| N=50 proportional coverage | 50 | 13 | Main N=50 point |
| N=100 proportional coverage | 100 | 25 | Exploratory |
| N=200 proportional coverage | 200 | 50 | Exploratory; not a formal certificate |

## Initial resource requests

Use separate single-node jobs because the memory needs differ:

| Point | Partition | CPUs / Mosek threads | Memory | Walltime |
|---|---|---:|---:|---:|
| N=30, k=8 | `64c512g` | 16 | 64G | 06:00:00 |
| N=50, k=5 or 13 | `64c512g` | 20 | 128G | 12:00:00 |
| N=100, k=25 | `64c512g` | 40 | 256G | 2-00:00:00 |
| N=200, k=50 | `64c512g` | 64 | 480G | 7-00:00:00 |

These CPU counts keep every request below the partition's 8,000 MB-per-CPU ceiling. They are conservative first requests, not measured peaks. After completion, use `sacct MaxRSS` to right-size later runs. Move N=200 to `huge` only after a documented OOM or a measured projection above one 512 GB node.

## Submission workflow

Create one fresh campaign timestamp locally:

```bash
campaign="$(date +%Y%m%d-%H%M%S)"
```

For each point, construct an exact command of this form:

```bash
run-remote -H siyuan sbatch --test-only \
  --partition=64c512g --qos=normal --nodes=1 \
  --cpus-per-task=<cpus> --mem=<memory> --time=<walltime> \
  --job-name=bge-n<N>-k<K> --output=slurm-%x-%j.out \
  --export=ALL,BGE_SYSTEM_SIZE=<N>,BGE_N_SUPER=<K>,BGE_RUN_NAME=<campaign>-n<N>-k<K> \
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
5. Report raw and strict improvement separately.
6. A negative `strict_improvement` is a valid scientific outcome, not an operational failure.
7. Keep `formal_certificate=false` in every report.

After pulling results, summarize them:

```bash
python3 .agents/skills/boundgsenergy-hpc-operations/scripts/summarize_size_screen.py \
  results/<run1>/size_screen.json \
  results/<run2>/size_screen.json
```

Do not start a formal N=200 certificate claim until the size-specific model construction, independently replayable transcript, and finite-size upper certificate have been generalized beyond N=20.
