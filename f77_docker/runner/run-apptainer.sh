#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${F77_RUNNER_ENV_FILE:-${SCRIPT_DIR}/f77pipeline.env}"
CONTAINER_PWD=
CHECK_ONLY=0
EXTRA_BINDS=()

# ==========================================
# Function: Display the Apptainer runner command-line contract
# Method: Print supported configuration, bind, working-directory, and check options
# ==========================================
usage() {
    printf '%s\n' \
        "Usage: $0 [--env-file FILE] [--pwd CONTAINER_DIR]" \
        "          [--bind SRC:DEST[:ro|rw]] [--check] [--] [COMMAND [ARG...]]"
}

# ==========================================
# Function: Stop execution with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Require one nonempty shell variable
# Method: Resolve the variable indirectly and reject missing configuration
# ==========================================
require_variable() {
    local variable_name="$1"
    [[ -n "${!variable_name:-}" ]] ||
        die "required variable ${variable_name} is empty in ${ENV_FILE}"
}

# ==========================================
# Function: Validate a host bind source directory
# Method: Require an absolute path to an existing directory visible on this node
# ==========================================
require_host_directory() {
    local variable_name="$1"
    local directory_path="${!variable_name}"

    [[ "${directory_path}" == /* ]] ||
        die "${variable_name} must be an absolute host path"
    [[ -d "${directory_path}" ]] ||
        die "${variable_name} does not exist on $(hostname): ${directory_path}"
}

# ==========================================
# Function: Validate a container destination path
# Method: Require an absolute path and reject bind-delimiter characters
# ==========================================
require_container_path() {
    local variable_name="$1"
    local container_path="${!variable_name}"

    [[ "${container_path}" == /* ]] ||
        die "${variable_name} must be an absolute container path"
    [[ "${container_path}" != *:* && "${container_path}" != *,* ]] ||
        die "${variable_name} contains an unsupported bind delimiter"
}

# ==========================================
# Function: Locate the HPC container runtime
# Method: Honor APPTAINER_BIN, then try Apptainer and Singularity in order
# ==========================================
select_runtime() {
    if [[ -n "${APPTAINER_BIN:-}" ]]; then
        command -v "${APPTAINER_BIN}" >/dev/null 2>&1 ||
            die "container runtime is not available: ${APPTAINER_BIN}"
        printf '%s\n' "${APPTAINER_BIN}"
        return
    fi

    if command -v apptainer >/dev/null 2>&1; then
        printf '%s\n' apptainer
    elif command -v singularity >/dev/null 2>&1; then
        printf '%s\n' singularity
    else
        die "neither apptainer nor singularity is available"
    fi
}

# ==========================================
# Function: Load and validate all shared-filesystem and container paths
# Method: Source the trusted env file and check every required runtime input
# ==========================================
load_configuration() {
    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"

    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    local required_variables=(
        F77_SIF
        F77_SOURCE_HOST
        F77_SOURCE_CONTAINER
        ASTROMETRY_CAT_HOST
        ASTROMETRY_CAT_CONTAINER
        SOURCE_CAT_HOST
        SOURCE_CAT_CONTAINER
        FLAT_PATH_HOST
        FLAT_PATH_CONTAINER
        PROCESS_DATA_HOST
        PROCESS_DATA_CONTAINER
    )
    local variable_name

    for variable_name in "${required_variables[@]}"; do
        require_variable "${variable_name}"
    done

    [[ "${F77_SIF}" == /* ]] || die "F77_SIF must be an absolute path"
    [[ -r "${F77_SIF}" ]] || die "SIF image is not readable: ${F77_SIF}"
    [[ -n "${CONTAINER_PWD}" ]] ||
        CONTAINER_PWD="${F77_SOURCE_CONTAINER}"

    for variable_name in \
        F77_SOURCE_HOST \
        ASTROMETRY_CAT_HOST \
        SOURCE_CAT_HOST \
        FLAT_PATH_HOST \
        PROCESS_DATA_HOST; do
        require_host_directory "${variable_name}"
    done

    for variable_name in \
        F77_SOURCE_CONTAINER \
        ASTROMETRY_CAT_CONTAINER \
        SOURCE_CAT_CONTAINER \
        FLAT_PATH_CONTAINER \
        PROCESS_DATA_CONTAINER; do
        require_container_path "${variable_name}"
    done
}

# ==========================================
# Function: Remove host OpenMPI transport selections from the container command
# Method: Unset only OpenMPI MCA overrides while preserving Slurm and PMI variables
# ==========================================
scrub_host_openmpi_environment() {
    [[ "${HPC_SCRUB_OPENMPI_ENV:-}" =~ ^[01]$ ]] ||
        die "HPC_SCRUB_OPENMPI_ENV must be 0 or 1"

    if [[ "${HPC_SCRUB_OPENMPI_ENV}" == "1" ]]; then
        unset OMPI_MCA_mtl
        unset OMPI_MCA_osc
    fi
}

# ==========================================
# Function: Parse runner options without consuming the application command
# Method: Shift known flags and retain all trailing arguments for exec
# ==========================================
parse_arguments() {
    while (($#)); do
        case "$1" in
            --env-file)
                (($# >= 2)) || die "--env-file requires a path"
                ENV_FILE="$2"
                shift 2
                ;;
            --pwd)
                (($# >= 2)) || die "--pwd requires a container path"
                CONTAINER_PWD="$2"
                shift 2
                ;;
            --bind)
                (($# >= 2)) || die "--bind requires SRC:DEST[:ro|rw]"
                EXTRA_BINDS+=("$2")
                shift 2
                ;;
            --check)
                CHECK_ONLY=1
                shift
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                die "unknown option: $1"
                ;;
            *)
                break
                ;;
        esac
    done

    RUN_COMMAND=("$@")
}

# ==========================================
# Function: Execute a command in the immutable SIF with Compose-equivalent binds
# Method: Build isolated bind arguments and replace the wrapper with the runtime
# ==========================================
run_container() {
    local runtime_command="$1"
    local runtime_arguments=(
        exec
        --no-home
        --pwd "${CONTAINER_PWD}"
        --bind "${F77_SOURCE_HOST}:${F77_SOURCE_CONTAINER}:rw"
        --bind "${ASTROMETRY_CAT_HOST}:${ASTROMETRY_CAT_CONTAINER}:ro"
        --bind "${SOURCE_CAT_HOST}:${SOURCE_CAT_CONTAINER}:ro"
        --bind "${FLAT_PATH_HOST}:${FLAT_PATH_CONTAINER}:ro"
        --bind "${PROCESS_DATA_HOST}:${PROCESS_DATA_CONTAINER}:rw"
        --env "F77_SOURCE_CONTAINER=${F77_SOURCE_CONTAINER}"
        --env "PROCESS_DATA_CONTAINER=${PROCESS_DATA_CONTAINER}"
    )
    local bind_specification

    for bind_specification in "${EXTRA_BINDS[@]}"; do
        runtime_arguments+=(--bind "${bind_specification}")
    done

    if [[ -n "${HPC_EXTRA_BINDS:-}" ]]; then
        runtime_arguments+=(--bind "${HPC_EXTRA_BINDS}")
    fi
    if [[ -n "${FI_PROVIDER:-}" ]]; then
        runtime_arguments+=(--env "FI_PROVIDER=${FI_PROVIDER}")
    fi
    if [[ -n "${FI_PROVIDER_PATH:-}" ]]; then
        runtime_arguments+=(--env "FI_PROVIDER_PATH=${FI_PROVIDER_PATH}")
    fi

    if ((CHECK_ONLY)); then
        RUN_COMMAND=(
            bash -lc
            'set -euo pipefail
             test "$(gcc -dumpversion)" = "4.8.5"
             test "$(gfortran -dumpversion)" = "4.8.5"
             mpichversion
             mpichversion | grep -F "MPICH Version:" | grep -F "4.1.2"
             mpiexec -info |
                 grep -F "Launchers available:" | grep -F "slurm"
             mpiexec -info |
                 grep -F "Resource management kernels available:" |
                 grep -F "slurm"
             grep -F "Version: 4.3.1" /opt/f77stack/lib/pkgconfig/cfitsio.pc
             test -d "${F77_SOURCE_CONTAINER}"
             test -d "${PROCESS_DATA_CONTAINER}"
             printf "Apptainer image and bind checks passed on %s.\n" "$(hostname)"'
        )
    elif ((${#RUN_COMMAND[@]} == 0)); then
        RUN_COMMAND=(bash)
    fi

    exec "${runtime_command}" "${runtime_arguments[@]}" \
        "${F77_SIF}" "${RUN_COMMAND[@]}"
}

# ==========================================
# Function: Coordinate argument parsing, configuration checks, and container exec
# Method: Resolve all inputs before replacing the wrapper process
# ==========================================
main() {
    local runtime_command

    parse_arguments "$@"
    load_configuration
    require_container_path CONTAINER_PWD
    scrub_host_openmpi_environment
    runtime_command="$(select_runtime)"
    run_container "${runtime_command}"
}

main "$@"
