#!/usr/bin/env bash

set -euo pipefail

REPOSITORY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ==========================================
# Function: Reject local source archives and oversized repository files
# Method: Inspect tracked project candidates without modifying the repository
# ==========================================
check_repository_payload() {
    if find "${REPOSITORY_DIR}" -type f \
        \( -name '*.tar' -o -name '*.tar.gz' -o -name '*.tar.bz2' \
           -o -name '*.tar.xz' \) -print -quit | grep -q .; then
        printf 'Source archive found in public repository.\n' >&2
        return 1
    fi

    if find "${REPOSITORY_DIR}" -type f -size +95M -print -quit | grep -q .; then
        printf 'A repository file exceeds the 95 MiB safety limit.\n' >&2
        return 1
    fi
}

# ==========================================
# Function: Reject private cluster paths and local source COPY directives
# Method: Search the public project for known private path markers
# ==========================================
check_public_content() {
    local private_home_marker="/home/""alatrion"
    local private_cluster_marker="/lustre/home/""acct-"

    if grep -R -n -F \
        -e "${private_home_marker}" \
        -e "${private_cluster_marker}" \
        --exclude=.env \
        --exclude=f77pipeline.env \
        --exclude=site-env.sh \
        --exclude-dir=.git "${REPOSITORY_DIR}"; then
        printf 'Private path marker found in public repository.\n' >&2
        return 1
    fi

    if grep -n -E '^[[:space:]]*COPY[[:space:]]+sources/' \
        "${REPOSITORY_DIR}/Dockerfile"; then
        printf 'Dockerfile copies a local source directory.\n' >&2
        return 1
    fi
}

# ==========================================
# Function: Keep Docker Compose and HPC mount defaults synchronized
# Method: Compare each shared variable directly between the two example env files
# ==========================================
check_mount_defaults() {
    local docker_environment="${REPOSITORY_DIR}/.env.example"
    local hpc_environment
    local variable_name
    local docker_value
    local hpc_value
    local hpc_environments=(
        "${REPOSITORY_DIR}/runner/f77pipeline.env.example"
        "${REPOSITORY_DIR}/runner/f77pipeline.pilogin-openmpi.env.example"
    )
    local shared_variables=(
        F77_SOURCE_HOST
        ASTROMETRY_CAT_HOST
        ASTROMETRY_CAT_CONTAINER
        SOURCE_CAT_HOST
        SOURCE_CAT_CONTAINER
        FLAT_PATH_HOST
        FLAT_PATH_CONTAINER
        PROCESS_DATA_HOST
        PROCESS_DATA_CONTAINER
    )

    for hpc_environment in "${hpc_environments[@]}"; do
        for variable_name in "${shared_variables[@]}"; do
            docker_value="$(
                sed -n "s/^${variable_name}=//p" "${docker_environment}"
            )"
            hpc_value="$(
                sed -n "s/^${variable_name}=//p" "${hpc_environment}"
            )"
            if [[ -z "${docker_value}" ||
                  "${docker_value}" != "${hpc_value}" ]]; then
                printf 'Mount default differs for %s in %s.\n' \
                    "${variable_name}" "${hpc_environment}" >&2
                return 1
            fi
        done
    done
}

# ==========================================
# Function: Enforce the single-file HPC runtime configuration contract
# Method: Source both examples and reject duplicated path or module literals
# ==========================================
check_hpc_configuration_contract() {
    local hpc_environment
    local hpc_environments=(
        "${REPOSITORY_DIR}/runner/f77pipeline.env.example"
        "${REPOSITORY_DIR}/runner/f77pipeline.pilogin-openmpi.env.example"
    )
    local runtime_scripts=(
        "${REPOSITORY_DIR}/runner/f77pipeline.slurm"
        "${REPOSITORY_DIR}/runner/mpi-smoke-test.slurm"
        "${REPOSITORY_DIR}/runner/run-apptainer.sh"
    )

    for hpc_environment in "${hpc_environments[@]}"; do
        if ! bash -eu -o pipefail -c '
            source "$1"
            module_declaration="$(declare -p HPC_MODULES)"
            [[ "${module_declaration}" == "declare -a "* ]]
            [[ "${F77_SOURCE_CONTAINER}" == /* ]]
            [[ "${F77_EXPO_LIST_CONTAINER}" == "${PROCESS_DATA_CONTAINER%/}/expo_list.list" ]]
            [[ "${F77_EXECUTABLE}" == "${F77_SOURCE_CONTAINER%/}/Fourier_Quad_Pipe" ]]
            [[ "${F77_BUILD_JOBS}" =~ ^[1-9][0-9]*$ ]]
            [[ "${F77_MAKE_CLEAN}" =~ ^[01]$ ]]
            [[ "${HPC_SCRUB_OPENMPI_ENV}" =~ ^[01]$ ]]
            [[ "${MPI_LAUNCH_MODE}" == "mpiexec" || "${MPI_LAUNCH_MODE}" == "srun" ]]
        ' _ "${hpc_environment}"; then
            printf 'Invalid HPC configuration contract in %s.\n' \
                "${hpc_environment}" >&2
            return 1
        fi
    done

    if grep -n -F '/workspace/f77' "${runtime_scripts[@]}"; then
        printf 'Source container path is duplicated outside the env examples.\n' \
            >&2
        return 1
    fi

    if grep -n -E '^[[:space:]]*module[[:space:]]+load[[:space:]]+gcc/' \
        "${REPOSITORY_DIR}/runner/f77pipeline-pilogin-openmpi.slurm" \
        "${REPOSITORY_DIR}/runner/mpi-smoke-test-pilogin-openmpi.slurm"; then
        printf 'Host module versions are duplicated outside f77pipeline.env.\n' \
            >&2
        return 1
    fi
}

# ==========================================
# Function: Enforce the standalone runner directory layout
# Method: Require bundled smoke sources and reject legacy directory contracts
# ==========================================
check_runner_layout_contract() {
    local legacy_hpc_directory="${REPOSITORY_DIR}/""hpc"
    local legacy_hpc_marker="F77_""HPC_"
    local legacy_project_marker="F77_""PROJECT_DIR"
    local legacy_source_marker="Code_""auto"
    local layout_files=(
        "${REPOSITORY_DIR}/README.md"
        "${REPOSITORY_DIR}/README-CN.md"
        "${REPOSITORY_DIR}/runner"
    )

    [[ ! -e "${legacy_hpc_directory}" ]] || {
        printf 'Legacy hpc directory still exists.\n' >&2
        return 1
    }

    [[ -r "${REPOSITORY_DIR}/runner/tests/mpi_identity.f" &&
       -r "${REPOSITORY_DIR}/runner/tests/mpi_lapack_smoke.f" ]] || {
        printf 'Standalone runner smoke-test sources are incomplete.\n' >&2
        return 1
    }

    if grep -R -n -E \
        -e "${legacy_hpc_marker}(DIR|ENV_FILE)" \
        -e "${legacy_project_marker}" \
        -e "${legacy_source_marker}" \
        -- "${layout_files[@]}"; then
        printf 'Legacy HPC layout marker found in runtime assets.\n' >&2
        return 1
    fi
}

# ==========================================
# Function: Validate portable text and executable metadata for HPC assets
# Method: Reject CRLF, overlong fixed-form Fortran, or non-executable runners
# ==========================================
check_hpc_assets() {
    local executable_file
    local executable_files=(
        "${REPOSITORY_DIR}/runner/f77pipeline.slurm"
        "${REPOSITORY_DIR}/runner/f77pipeline-pilogin-openmpi.slurm"
        "${REPOSITORY_DIR}/runner/inspect-cluster-mpi.sh"
        "${REPOSITORY_DIR}/runner/mpi-smoke-test.slurm"
        "${REPOSITORY_DIR}/runner/mpi-smoke-test-pilogin-openmpi.slurm"
        "${REPOSITORY_DIR}/runner/pull-sif.sh"
        "${REPOSITORY_DIR}/runner/run-apptainer.sh"
    )

    if grep -R -I -l $'\r$' \
        --exclude-dir=.git "${REPOSITORY_DIR}" | grep -q .; then
        printf 'CRLF line ending found in public repository text.\n' >&2
        return 1
    fi

    if awk 'length($0) > 72 { exit 1 }' \
        "${REPOSITORY_DIR}/runner/tests/mpi_identity.f" \
        "${REPOSITORY_DIR}/runner/tests/mpi_lapack_smoke.f"; then
        :
    else
        printf 'Fixed-form Fortran line exceeds 72 columns.\n' >&2
        return 1
    fi

    for executable_file in "${executable_files[@]}"; do
        [[ -x "${executable_file}" ]] || {
            printf 'HPC runner is not executable: %s\n' \
                "${executable_file}" >&2
            return 1
        }
    done
}

# ==========================================
# Function: Run all repository publication checks
# Method: Execute payload, content, mount, and HPC checks in deterministic order
# ==========================================
main() {
    check_repository_payload
    check_public_content
    check_mount_defaults
    check_hpc_configuration_contract
    check_runner_layout_contract
    check_hpc_assets
    printf 'Public repository checks passed.\n'
}

main "$@"
