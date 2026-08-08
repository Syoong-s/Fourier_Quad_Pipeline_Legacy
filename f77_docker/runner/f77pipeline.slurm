#!/usr/bin/env bash

#SBATCH --job-name=f77pipeline
#SBATCH --partition=cpu
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=f77pipeline-%j.out
#SBATCH --error=f77pipeline-%j.err

set -euo pipefail

RUNNER_DIR="${F77_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
ENV_FILE="${F77_RUNNER_ENV_FILE:-${RUNNER_DIR}/f77pipeline.env}"
RUNNER="${RUNNER_DIR}/run-apptainer.sh"

# ==========================================
# Function: Stop a batch job with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Load the host modules requested by the HPC configuration
# Method: Validate the Bash module array and load it before MPI or Apptainer use
# ==========================================
load_host_modules() {
    local module_declaration

    module_declaration="$(declare -p HPC_MODULES 2>/dev/null)" ||
        die "HPC_MODULES must be defined as a Bash array in ${ENV_FILE}"
    [[ "${module_declaration}" == "declare -a "* ]] ||
        die "HPC_MODULES must be a Bash indexed array in ${ENV_FILE}"

    if ((${#HPC_MODULES[@]} == 0)); then
        return
    fi

    command -v module >/dev/null 2>&1 ||
        die "the module command is unavailable but HPC_MODULES is not empty"
    module load "${HPC_MODULES[@]}" ||
        die "failed to load one or more modules from HPC_MODULES"
}

# ==========================================
# Function: Load user paths and cluster module initialization
# Method: Source trusted configuration, optional site setup, and configured modules
# ==========================================
load_job_environment() {
    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    [[ -x "${RUNNER}" ]] || die "Apptainer runner is not executable: ${RUNNER}"

    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    if [[ -n "${SITE_ENV_SCRIPT:-}" ]]; then
        [[ -r "${SITE_ENV_SCRIPT}" ]] ||
            die "SITE_ENV_SCRIPT is not readable: ${SITE_ENV_SCRIPT}"
        # shellcheck disable=SC1090
        source "${SITE_ENV_SCRIPT}"
    fi

    load_host_modules

    [[ -n "${F77_SOURCE_CONTAINER:-}" &&
       "${F77_SOURCE_CONTAINER}" == /* ]] ||
        die "F77_SOURCE_CONTAINER must be an absolute container path"
    [[ -n "${F77_EXPO_LIST_CONTAINER:-}" &&
       "${F77_EXPO_LIST_CONTAINER}" == /* ]] ||
        die "F77_EXPO_LIST_CONTAINER must be an absolute container path"
    [[ -n "${F77_EXECUTABLE:-}" && "${F77_EXECUTABLE}" == /* ]] ||
        die "F77_EXECUTABLE must be an absolute container path"
}

# ==========================================
# Function: Compile the currently mounted f77pipeline source exactly once
# Method: Run make on the batch node before starting any multi-rank job step
# ==========================================
compile_pipeline() {
    [[ "${F77_BUILD_JOBS:-}" =~ ^[1-9][0-9]*$ ]] ||
        die "F77_BUILD_JOBS must be a positive integer"
    [[ "${F77_MAKE_CLEAN:-}" =~ ^[01]$ ]] ||
        die "F77_MAKE_CLEAN must be 0 or 1"

    if [[ "${F77_MAKE_CLEAN}" == "1" ]]; then
        F77_RUNNER_ENV_FILE="${ENV_FILE}" \
            "${RUNNER}" --pwd "${F77_SOURCE_CONTAINER}" -- \
            make -C "${F77_SOURCE_CONTAINER}" clean
    fi

    F77_RUNNER_ENV_FILE="${ENV_FILE}" \
        "${RUNNER}" --pwd "${F77_SOURCE_CONTAINER}" -- \
        make -C "${F77_SOURCE_CONTAINER}" -j"${F77_BUILD_JOBS}"
}

# ==========================================
# Function: Launch one container process for every Slurm MPI rank
# Method: Select host MPICH hybrid mode or an explicitly validated srun PMI mode
# ==========================================
launch_pipeline() {
    local pipeline_arguments=("$@")
    local container_command

    [[ -n "${SLURM_NTASKS:-}" ]] ||
        die "SLURM_NTASKS is unavailable; submit this file with sbatch"
    [[ -n "${F77_EXECUTABLE:-}" && "${F77_EXECUTABLE}" == /* ]] ||
        die "F77_EXECUTABLE must be an absolute container path"

    if ((${#pipeline_arguments[@]} == 0)); then
        pipeline_arguments=("${F77_EXPO_LIST_CONTAINER}")
    fi

    container_command=(
        "${RUNNER}"
        --env-file "${ENV_FILE}"
        --pwd "${PROCESS_DATA_CONTAINER}"
        --
        "${F77_EXECUTABLE}"
        "${pipeline_arguments[@]}"
    )

    case "${MPI_LAUNCH_MODE:-}" in
        mpiexec)
            [[ -n "${MPI_LAUNCHER:-}" ]] ||
                die "MPI_LAUNCHER is empty for mpiexec mode"
            command -v "${MPI_LAUNCHER}" >/dev/null 2>&1 ||
                die "host MPICH launcher is unavailable: ${MPI_LAUNCHER}"
            "${MPI_LAUNCHER}" -n "${SLURM_NTASKS}" \
                "${container_command[@]}"
            ;;
        srun)
            command -v srun >/dev/null 2>&1 || die "srun is unavailable"
            [[ -n "${SLURM_MPI_TYPE:-}" ]] ||
                die "SLURM_MPI_TYPE is empty for srun mode"
            srun --kill-on-bad-exit=1 --ntasks="${SLURM_NTASKS}" \
                --mpi="${SLURM_MPI_TYPE}" "${container_command[@]}"
            ;;
        *)
            die "MPI_LAUNCH_MODE must be mpiexec or srun"
            ;;
    esac
}

# ==========================================
# Function: Validate the image, build the pipeline, and run the requested stage
# Method: Sequence all single-process preparation before the MPI launch barrier
# ==========================================
main() {
    load_job_environment

    printf 'Slurm job %s: %s nodes, %s ranks\n' \
        "${SLURM_JOB_ID:-unknown}" \
        "${SLURM_JOB_NUM_NODES:-unknown}" \
        "${SLURM_NTASKS:-unknown}"

    command -v mpichversion >/dev/null 2>&1 && mpichversion || true
    F77_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" --check
    compile_pipeline
    F77_RUNNER_ENV_FILE="${ENV_FILE}" \
        "${RUNNER}" --pwd "${F77_SOURCE_CONTAINER}" -- \
        test -x "${F77_EXECUTABLE}"
    launch_pipeline "$@"
}

main "$@"
