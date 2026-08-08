#!/usr/bin/env bash

#SBATCH --job-name=f77-mpi-smoke
#SBATCH --partition=cpu
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00
#SBATCH --output=f77-mpi-smoke-%j.out
#SBATCH --error=f77-mpi-smoke-%j.err

set -euo pipefail

RUNNER_DIR="${F77_RUNNER_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
ENV_FILE="${F77_RUNNER_ENV_FILE:-${RUNNER_DIR}/f77pipeline.env}"
RUNNER="${RUNNER_DIR}/run-apptainer.sh"
TEST_SOURCE_DIR="${RUNNER_DIR}/tests"
SMOKE_DIR=

# ==========================================
# Function: Stop a smoke test with a concise diagnostic
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
# Function: Remove only the unique shared smoke directory created by this job
# Method: Validate the mktemp result and recursively remove that exact path
# ==========================================
cleanup_smoke_directory() {
    if [[ -n "${SMOKE_DIR}" && -d "${SMOKE_DIR}" ]]; then
        rm -rf -- "${SMOKE_DIR}"
    fi
}

# ==========================================
# Function: Load paths and cluster modules for the smoke test
# Method: Source trusted configuration, optional site setup, and configured modules
# ==========================================
load_job_environment() {
    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    [[ -x "${RUNNER}" ]] || die "Apptainer runner is not executable: ${RUNNER}"
    [[ -r "${TEST_SOURCE_DIR}/mpi_identity.f" ]] ||
        die "MPI identity test source is missing"
    [[ -r "${TEST_SOURCE_DIR}/mpi_lapack_smoke.f" ]] ||
        die "MPI/LAPACK test source is missing"

    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    if [[ -n "${SITE_ENV_SCRIPT:-}" ]]; then
        [[ -r "${SITE_ENV_SCRIPT}" ]] ||
            die "SITE_ENV_SCRIPT is not readable: ${SITE_ENV_SCRIPT}"
        # shellcheck disable=SC1090
        source "${SITE_ENV_SCRIPT}"
    fi

    load_host_modules
}

# ==========================================
# Function: Compile MPI identity and LAPACK tests into shared temporary storage
# Method: Bind read-only test sources and one unique writable job directory
# ==========================================
compile_smoke_programs() {
    [[ -n "${HPC_SHARED_SCRATCH_HOST:-}" ]] ||
        die "HPC_SHARED_SCRATCH_HOST is empty"
    [[ -d "${HPC_SHARED_SCRATCH_HOST}" && -w "${HPC_SHARED_SCRATCH_HOST}" ]] ||
        die "shared smoke parent is not writable: ${HPC_SHARED_SCRATCH_HOST}"

    SMOKE_DIR="$(mktemp -d \
        "${HPC_SHARED_SCRATCH_HOST%/}/.f77pipeline-smoke.${SLURM_JOB_ID}.XXXXXX")"
    trap cleanup_smoke_directory EXIT

    F77_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" \
        --bind "${TEST_SOURCE_DIR}:/tests:ro" \
        --bind "${SMOKE_DIR}:/smoke:rw" \
        --pwd /smoke -- \
        mpif77 -O2 /tests/mpi_identity.f -o /smoke/mpi_identity

    F77_RUNNER_ENV_FILE="${ENV_FILE}" "${RUNNER}" \
        --bind "${TEST_SOURCE_DIR}:/tests:ro" \
        --bind "${SMOKE_DIR}:/smoke:rw" \
        --pwd /smoke -- \
        mpif77 -O2 /tests/mpi_lapack_smoke.f \
        -L/opt/f77stack/lib -llapack -lblas \
        -o /smoke/mpi_lapack_smoke
}

# ==========================================
# Function: Launch one shared smoke executable across the allocation
# Method: Reuse the configured host mpiexec or validated srun PMI mode
# ==========================================
launch_smoke_program() {
    local executable_name="$1"
    local container_command=(
        "${RUNNER}"
        --env-file "${ENV_FILE}"
        --bind "${SMOKE_DIR}:/smoke:rw"
        --pwd /smoke
        --
        "/smoke/${executable_name}"
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
# Function: Build and run multi-node MPI validation programs
# Method: Print rank placement first, then verify MPI and reference LAPACK linkage
# ==========================================
main() {
    [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_NTASKS:-}" ]] ||
        die "submit this file with sbatch"
    load_job_environment
    compile_smoke_programs
    launch_smoke_program mpi_identity
    launch_smoke_program mpi_lapack_smoke
    printf 'Multi-node MPI smoke test completed successfully.\n'
}

main "$@"
