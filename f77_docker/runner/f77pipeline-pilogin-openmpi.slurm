#!/usr/bin/env bash

#SBATCH --job-name=f77pipeline-pmi2
#SBATCH --partition=cpu
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=f77pipeline-pmi2-%j.out
#SBATCH --error=f77pipeline-pmi2-%j.err

set -euo pipefail

export SLURM_CPUS_PER_TASK=1
export SLURM_TRES_PER_TASK=cpu=1

RUNNER_DIR="${F77_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
export F77_RUNNER_DIR="${RUNNER_DIR}"
export F77_RUNNER_ENV_FILE="${F77_RUNNER_ENV_FILE:-${RUNNER_DIR}/f77pipeline.env}"

exec "${RUNNER_DIR}/f77pipeline.slurm" "$@"
