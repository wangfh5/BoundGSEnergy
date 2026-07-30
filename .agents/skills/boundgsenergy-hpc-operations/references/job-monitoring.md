# Job monitoring

## Queue snapshot

Use `run-remote` and lead with anomalies:

```bash
run-remote -H siyuan squeue --me -o "%.18i %.9P %.32j %.2t %.10M %.10l %.6D %R"
run-remote -H siyuan squeue --me -t R -o "%.18i %.32j %.10M %.6D %R"
```

For a specific job:

```bash
run-remote -H siyuan scontrol show job <jobid>
run-remote -H siyuan tail -50 slurm-<job-name>-<jobid>.out
```

Confirm `WorkDir` is the canonical BoundGSEnergy mirror and the requested CPU, memory, partition, QOS and environment variables match the intended point.

## Progress interpretation

The model construction and exact replay phases can be quiet. A quiet log alone does not prove a hang. Check:

- log modification time and most recent phase message;
- scheduler elapsed time and node;
- live memory/CPU with `sstat`;
- comparison with a smaller point in the same campaign.

```bash
run-remote -H siyuan sstat -j <jobid>.batch --format=JobID,AveCPU,MaxRSS,AveRSS,MaxVMSize
```

Do not estimate completion time until at least one comparable phase has completed.

## Completion and retrieval

Use accounting only after a job leaves the queue:

```bash
run-remote -H siyuan sacct -j <jobid> --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList -P
```

Inspect the remote result before pulling:

```bash
run-remote -H siyuan test -f results/<run>/size_screen.json
run-remote -H siyuan du -sh results/<run>
```

Pull generated results without deleting either side:

```bash
sync-remote -n -m copy-pull --include-only results
sync-remote -m copy-pull --include-only results
```

Because `results/` is gitignored and the remote campaign mirror starts without historical results, pulling the top-level results directory is the reliable `sync-remote` path. Do not use destructive `pull` mode.

Fetch selected Slurm logs separately if needed:

```bash
sync-remote -m copy-pull --include-only "slurm-*.out"
```

Then run `scripts/summarize_size_screen.py` locally and report one row per point with N, `n_super`, job state, peak RSS, elapsed time, raw improvement, strict improvement and all structural gates.
