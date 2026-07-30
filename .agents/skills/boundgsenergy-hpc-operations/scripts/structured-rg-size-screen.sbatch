#!/usr/bin/env bash
#SBATCH --job-name=bge-structured-rg
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

: "${BGE_SYSTEM_SIZE:?Set BGE_SYSTEM_SIZE}"
: "${BGE_N_SUPER:?Set BGE_N_SUPER}"
: "${BGE_RUN_NAME:?Set BGE_RUN_NAME}"
: "${SLURM_SUBMIT_DIR:?Run this entrypoint through sbatch}"
: "${SLURM_CPUS_PER_TASK:?Request --cpus-per-task}"

BGE_COEFFICIENT_SCALE="${BGE_COEFFICIENT_SCALE:-16384}"
BGE_SOLVER_LOG="${BGE_SOLVER_LOG:-0}"
BGE_OMEGA_NORMALIZATION="${BGE_OMEGA_NORMALIZATION:-none}"
BGE_FEASIBILITY_TOLERANCE="${BGE_FEASIBILITY_TOLERANCE:-1e-7}"

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
[[ "$BGE_COEFFICIENT_SCALE" =~ ^[0-9]+$ ]] || {
    echo "BGE_COEFFICIENT_SCALE must be an integer" >&2
    exit 2
}
[[ "$BGE_SOLVER_LOG" =~ ^(0|1)$ ]] || {
    echo "BGE_SOLVER_LOG must be zero or one" >&2
    exit 2
}
[[ "$BGE_OMEGA_NORMALIZATION" =~ ^(none|dyadic_trace)$ ]] || {
    echo "BGE_OMEGA_NORMALIZATION must be none or dyadic_trace" >&2
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
echo "N=$BGE_SYSTEM_SIZE n_super=$BGE_N_SUPER mosek_threads=$SLURM_CPUS_PER_TASK coefficient_scale=$BGE_COEFFICIENT_SCALE solver_log=$BGE_SOLVER_LOG omega_normalization=$BGE_OMEGA_NORMALIZATION feasibility_tolerance=$BGE_FEASIBILITY_TOLERANCE"
echo "started=$(date -u +%Y-%m-%dT%H:%M:%SZ)"

solver_log_arg=()
if [[ "$BGE_SOLVER_LOG" == 1 ]]; then
    solver_log_arg=(--solver-log)
fi

julia +1.11.5 --startup-file=no --project=. \
    scripts/heisenberg1d/structured_rg_size_screen.jl \
    "--system-size=$BGE_SYSTEM_SIZE" \
    "--n-super=$BGE_N_SUPER" \
    "--run-name=$BGE_RUN_NAME" \
    "--mosek-threads=$SLURM_CPUS_PER_TASK" \
    "--coefficient-scale=$BGE_COEFFICIENT_SCALE" \
    "--omega-normalization=$BGE_OMEGA_NORMALIZATION" \
    "${solver_log_arg[@]}" \
    "--feasibility-tolerance=$BGE_FEASIBILITY_TOLERANCE"

echo "finished=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
