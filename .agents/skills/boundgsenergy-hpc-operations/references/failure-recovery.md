# Failure recovery

Diagnose first. Never cancel, delete, edit remote files, or resubmit without explicit authorization.

## Triage

```bash
run-remote -H siyuan sacct -j <jobid> --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList,Reason -P
run-remote -H siyuan scontrol show job <jobid>
run-remote -H siyuan tail -100 slurm-<job-name>-<jobid>.out
```

Classify the outcome:

| Class | Evidence | Next action |
|---|---|---|
| Environment | Julia is not 1.11.5, package import fails, license error | Fix bootstrap; rerun the smoke job |
| OOM | `OUT_OF_MEMORY`, memory cgroup kill, MaxRSS near request | Increase memory on `64c512g`; use `huge` only if one node is insufficient |
| Timeout | `TIMEOUT` with continuing progress | Increase walltime within partition limit or reduce the parameter point |
| Node/system | `NODE_FAIL`, filesystem or network error | Retry once; exclude a node only after repeated correlated evidence |
| Solver status | `SLOW_PROGRESS`, infeasible, unknown, or non-OPTIMAL | Treat as scientific/numerical failure; do not call it a bound |
| Strict gate | Operational success but nonpositive strict improvement | Preserve result; stop certification-scale interpretation and revisit residual control |
| Model error | Support, coefficient audit, monotonicity, or finite-ceiling assertion fails | Do not retry unchanged; diagnose code/model provenance |

## Safe cancellation

Show the exact job first:

```bash
run-remote -H siyuan squeue -j <jobid> -o "%.18i %.32j %.2t %.10M %.6D %R"
```

After explicit authorization:

```bash
run-remote -H siyuan scancel <jobid>
```

Verify the final state with `sacct`.

## Resubmission

- Always choose a fresh `BGE_RUN_NAME`; never reuse a partially created result directory.
- Preserve the failed log and any partial result directory.
- Change only the resource or numerical parameter justified by the diagnosis.
- Repeat `sbatch --test-only`, show the delta from the prior request, and obtain authorization.
- After submission, verify `WorkDir`, command, resources and exported N/`n_super` with `scontrol show job`.

For a transient node failure, retry the same scientific point. For OOM or timeout, record both attempts in the final campaign table. For non-OPTIMAL solver status or failed mathematical gates, do not hide the failure behind an operational retry.
