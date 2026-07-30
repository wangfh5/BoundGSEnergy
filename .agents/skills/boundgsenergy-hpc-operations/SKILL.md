---
name: boundgsenergy-hpc-operations
description: Operate BoundGSEnergy Structured-NPA/RG jobs on the Siyuan Slurm cluster through sync-remote and run-remote. Use for remote bootstrap, Julia/Mosek checks, N=30/50/100/200 size screens, job submission, monitoring, result retrieval, failure diagnosis, or safe resubmission.
---

# BoundGSEnergy HPC Operations

Use `run-remote` as the only remote command channel and `sync-remote` as the only repository transfer channel. Do not replace them with ad-hoc `ssh`, `scp`, `rsync`, or GitHub-mediated deployment.

## Route by request

| Request | Read |
|---|---|
| First connection, remote mirror, Julia 1.11.5, Mosek license, partition choice | [references/siyuan-environment.md](references/siyuan-environment.md) |
| Structured-NPA + RG size campaign, resource matrix, submit commands, result gates | [references/structured-rg-campaign.md](references/structured-rg-campaign.md) |
| Queue status, logs, resource use, result retrieval and summaries | [references/job-monitoring.md](references/job-monitoring.md) |
| OOM, timeout, node failure, solver failure, cancellation or resubmission | [references/failure-recovery.md](references/failure-recovery.md) |

Load only the references needed for the current request.

## Always apply

1. Verify `command -v run-remote && command -v sync-remote`.
2. Run `run-remote -H siyuan precheck` and capture `git status --short --branch`.
3. State the target returned by `run-remote -H siyuan whoami` before any mutation.
4. Treat status and file inspection as read-only. Before `sync-remote`, license transfer, `sbatch`, `scancel`, deletion, or overwrite: show the exact action and obtain explicit authorization.
5. Preview repository transfer with `sync-remote -n`. Preview scheduler feasibility with `run-remote sbatch --test-only ...` using the exact resource request.
6. Refuse to ship uncommitted tracked changes silently. For the first mirror use reviewed `handoff`; after the mirror exists prefer `git-push`, which is commit-only and fast-forward-only.
7. Never put `mosek.lic` in the repository, a result directory, a Slurm log, or command output.
8. Use fresh run names. Never overwrite or reuse a generated result directory.
9. A scheduler success is not scientific success. Fetch and inspect `size_screen.json`; require `OPTIMAL`, physical compatibility, exact coefficient audit, and monotonicity before accepting an exploratory result.
10. `structured_rg_size_screen.jl` sets `formal_certificate=false`. Never describe its endpoints as certified finite-size bounds, even when same-process strict replay is positive.

## Operational contract

- Siyuan commands run on the login alias `siyuan`; do not run Slurm commands through the data host.
- The canonical remote mirror is `~/code/BoundGSEnergy`, matching the local path below `$HOME`.
- Use independent jobs rather than one array for heterogeneous N=50/100/200 memory requests.
- Keep each job on one node. Mosek is shared-memory; multi-node allocation does not distribute one SDP.
- Match `--mosek-threads` to `SLURM_CPUS_PER_TASK` through the bundled batch entrypoint.
- Set BLAS/OpenMP helper threads to one so they do not oversubscribe the Mosek allocation.

## Bundled helpers

- `scripts/structured-rg-size-screen.sbatch`: compute-node entrypoint driven by `BGE_SYSTEM_SIZE`, `BGE_N_SUPER`, and `BGE_RUN_NAME`.
- `scripts/mosek_smoke.jl`: tiny license and solver smoke test.
- `scripts/summarize_size_screen.py`: local summary and structural-gate checker for fetched `size_screen.json` files.
