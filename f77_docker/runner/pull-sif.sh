#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${F77_RUNNER_ENV_FILE:-${SCRIPT_DIR}/f77pipeline.env}"

# ==========================================
# Function: Stop SIF acquisition with a concise diagnostic
# Method: Write one error message to stderr and return a nonzero status
# ==========================================
die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

# ==========================================
# Function: Locate the user-space HPC container runtime
# Method: Honor APPTAINER_BIN, then try Apptainer and Singularity
# ==========================================
select_runtime() {
    if [[ -n "${APPTAINER_BIN:-}" ]]; then
        command -v "${APPTAINER_BIN}" >/dev/null 2>&1 ||
            die "container runtime is not available: ${APPTAINER_BIN}"
        printf '%s\n' "${APPTAINER_BIN}"
    elif command -v apptainer >/dev/null 2>&1; then
        printf '%s\n' apptainer
    elif command -v singularity >/dev/null 2>&1; then
        printf '%s\n' singularity
    else
        die "neither apptainer nor singularity is available"
    fi
}

# ==========================================
# Function: Pull one immutable SIF without overwriting an existing image
# Method: Read the OCI reference from configuration and invoke the runtime pull
# ==========================================
main() {
    local runtime_command
    local output_directory

    [[ -r "${ENV_FILE}" ]] || die "configuration file is not readable: ${ENV_FILE}"
    # shellcheck disable=SC1090
    source "${ENV_FILE}"

    [[ -n "${OCI_IMAGE_URI:-}" ]] || die "OCI_IMAGE_URI is empty in ${ENV_FILE}"
    [[ -n "${F77_SIF:-}" && "${F77_SIF}" == /* ]] ||
        die "F77_SIF must be an absolute path"
    [[ ! -e "${F77_SIF}" ]] ||
        die "refusing to overwrite existing SIF: ${F77_SIF}"

    output_directory="$(dirname "${F77_SIF}")"
    [[ -d "${output_directory}" && -w "${output_directory}" ]] ||
        die "SIF output directory is not writable: ${output_directory}"

    runtime_command="$(select_runtime)"
    "${runtime_command}" pull "${F77_SIF}" "docker://${OCI_IMAGE_URI}"
    [[ -s "${F77_SIF}" ]] ||
        die "runtime returned without creating a nonempty SIF: ${F77_SIF}"
    printf 'Created SIF: %s\n' "${F77_SIF}"
}

main "$@"
