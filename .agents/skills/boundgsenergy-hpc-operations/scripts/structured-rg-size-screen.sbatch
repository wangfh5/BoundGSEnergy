#!/usr/bin/env bash
#SBATCH --job-name=bge-structured-rg
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

: "${BGE_SYSTEM_SIZE:?Set BGE_SYSTEM_SIZE}"
: "${BGE_N_SUPER:?Set BGE_N_SUPER}"
: "${BGE_RUN_NAME:?Set BGE_RUN_NAME}"
: "${SLURM_SUBMIT_DIR:?Run this entrypoint through sbatch}"
: "${SLURM_CPUS_PER_TASK:?Request --cpus-per-task}"

[[ "$BGE_SYSTEM_SIZE" =~ ^[0-9]+$ ]] || {
    echo "BGE_SYSTEM_SIZE must be an integer" >&2
    exit 2
}
[[ "$BGE_N_SUPER" =~ ^[0-9]+$ ]] || {
    echo "BGE_N_SUPER must be an integer" >&2
    exit 2
}
[[ "$BGE_RUN_NAME" != */* ]] || {
    echo "BGE_RUN_NAME must be one directory name" >&2
    exit 2
}

cd "$SLURM_SUBMIT_DIR"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "commit=$(git rev-parse HEAD)"
echo "host=$(hostname)"
echo "julia=$(julia +1.11.5 --version)"
echo "N=$BGE_SYSTEM_SIZE n_super=$BGE_N_SUPER mosek_threads=$SLURM_CPUS_PER_TASK"
echo "started=$(date -u +%Y-%m-%dT%H:%M:%SZ)"

julia +1.11.5 --startup-file=no --project=. \
    scripts/heisenberg1d/structured_rg_size_screen.jl \
    "--system-size=$BGE_SYSTEM_SIZE" \
    "--n-super=$BGE_N_SUPER" \
    "--run-name=$BGE_RUN_NAME" \
    "--mosek-threads=$SLURM_CPUS_PER_TASK" \
    --feasibility-tolerance=1e-7

echo "finished=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
