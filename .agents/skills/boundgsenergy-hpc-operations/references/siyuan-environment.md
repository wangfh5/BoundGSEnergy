# Siyuan environment

## Known baseline

The following facts were probed on 2026-07-30 and must be rechecked before submission:

- SSH alias: `siyuan`
- Account observed through SSH: `phyxxy-wfh`
- Remote home: `/dssg/home/acct-phyxxy/phyxxy-wfh`
- Canonical mirror: `~/code/BoundGSEnergy`
- Scheduler: Slurm
- Primary partition: `64c512g`, x86_64, 64 CPU and 512 GB per node, seven-day limit
- Smoke partition: `debug64c512g`, one-hour limit
- High-memory fallback: `huge`, about 2.9 TB per node, two-day limit
- QOS used by the account: `normal`; use `debug` only for the debug partition
- Installed Julia default: 1.11.7 through juliaup
- Required Julia: exactly 1.11.5

Queue state is dynamic. Re-probe it:

```bash
run-remote -H siyuan sinfo -h -o "%P|%a|%l|%D|%t|%c|%m|%G"
run-remote -H siyuan squeue --me -o "%.18i %.9P %.32j %.2t %.10M %.10l %.6D %R"
```

## Bootstrap the remote mirror

The first transfer uses `handoff` because it carries the complete Git worktree and can create an absent mirror. Review both tracked and untracked files first:

```bash
git status --short --branch
sync-remote -n -H siyuan -m handoff --suffix boundgsenergy
```

Only after explicit authorization:

```bash
sync-remote -H siyuan -m handoff --suffix boundgsenergy
run-remote -H siyuan git status --short --branch
run-remote -H siyuan git rev-parse HEAD
```

Later updates use committed history:

```bash
sync-remote -n -H siyuan -m git-push
sync-remote -H siyuan -m git-push
```

`git-push` refuses a missing mirror, diverged history, dirty tracked files, or non-fast-forward updates. Do not bypass those checks with destructive synchronization.

## Pin Julia and instantiate

The entrypoints reject every Julia version other than 1.11.5:

```bash
run-remote -H siyuan juliaup status
run-remote -H siyuan juliaup add 1.11.5
run-remote -H siyuan julia +1.11.5 --version
run-remote -H siyuan julia +1.11.5 --startup-file=no --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

The last three commands mutate the remote user environment and require authorization. Run package installation on the login node, not inside every compute job.

## Install the Mosek license securely

The local license is expected at `$HOME/mosek/mosek.lic`. Check only its presence and permissions; never print it:

```bash
test -r "$HOME/mosek/mosek.lic"
```

After the remote mirror exists and after explicit authorization, transfer it over `run-remote` stdin:

```bash
run-remote -H siyuan mkdir -p /dssg/home/acct-phyxxy/phyxxy-wfh/mosek
run-remote -H siyuan tee /dssg/home/acct-phyxxy/phyxxy-wfh/mosek/mosek.lic >/dev/null < "$HOME/mosek/mosek.lic"
run-remote -H siyuan chmod 600 /dssg/home/acct-phyxxy/phyxxy-wfh/mosek/mosek.lic
```

Do not use `sync-remote` for the license because repository synchronization must never carry credentials.

## Compute-node smoke test

Test the exact runtime, project environment and license on a compute node before a campaign:

```bash
run-remote -H siyuan sbatch --test-only --partition=debug64c512g --qos=debug --nodes=1 --cpus-per-task=4 --mem=16G --time=00:10:00 --job-name=bge-mosek-smoke --output=slurm-%x-%j.out --wrap="julia +1.11.5 --startup-file=no --project=. .agents/skills/boundgsenergy-hpc-operations/scripts/mosek_smoke.jl"
```

After the feasibility preview and explicit authorization, repeat without `--test-only`. Capture the job id, wait for completion, then require `OPTIMAL` in its log.

## Connection gotcha

The local SSH configuration may attempt a remote forward on port 45698 and print a warning when that port is already occupied. The warning is harmless when `run-remote ... echo ok` exits zero. Remote GitHub access may inherit a broken localhost proxy; the campaign does not depend on it because code moves directly through `sync-remote`.
