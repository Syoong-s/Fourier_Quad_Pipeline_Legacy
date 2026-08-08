#!/usr/bin/env bash

set -u

# ==========================================
# Function: Print one clearly delimited inspection section
# Method: Format a heading without writing any report file
# ==========================================
print_heading() {
    printf '\n== %s ==\n' "$1"
}

# ==========================================
# Function: Run a read-only command when it is available
# Method: Resolve the executable first and report absence without failing the audit
# ==========================================
run_if_available() {
    local command_name="$1"
    shift

    if command -v "${command_name}" >/dev/null 2>&1; then
        "$@"
    else
        printf '%s is not available in the current environment.\n' "${command_name}"
    fi
}

# ==========================================
# Function: Inspect architecture, scheduler, MPI, fabric, and container runtime
# Method: Execute version and configuration queries only; never create files
# ==========================================
main() {
    print_heading "Host"
    printf 'hostname: %s\n' "$(hostname)"
    printf 'architecture: %s\n' "$(uname -m)"
    printf 'kernel: %s\n' "$(uname -r)"
    getconf GNU_LIBC_VERSION 2>/dev/null || true

    print_heading "Commands"
    command -v sbatch srun mpiexec mpicc mpifort \
        apptainer singularity fi_info ibv_devinfo 2>/dev/null || true

    print_heading "Slurm"
    run_if_available srun srun --mpi=list
    if command -v scontrol >/dev/null 2>&1; then
        scontrol show config 2>/dev/null |
            grep -E '^(SlurmctldVersion|MpiDefault|MpiParams)' || true
    else
        printf 'scontrol is not available in the current environment.\n'
    fi

    print_heading "MPICH"
    run_if_available mpichversion mpichversion
    if command -v mpichversion >/dev/null 2>&1; then
        run_if_available mpicc mpicc -show
        run_if_available mpifort mpifort -show
        run_if_available mpiexec mpiexec -info
    fi

    print_heading "OpenMPI"
    run_if_available ompi_info ompi_info --version
    if command -v ompi_info >/dev/null 2>&1; then
        run_if_available mpicc mpicc --showme
        run_if_available mpifort mpifort --showme
        ompi_info --parsable --all 2>/dev/null |
            grep -E '^(ompi:version:full|config:cli|compiler:.*:version)' ||
            true
    fi

    print_heading "Fabric and RDMA"
    run_if_available fi_info fi_info --version
    run_if_available fi_info fi_info -l
    run_if_available ibv_devinfo ibv_devinfo -l
    printf 'locked-memory limit: '
    ulimit -l

    print_heading "Apptainer or Singularity"
    if command -v apptainer >/dev/null 2>&1; then
        apptainer version
    elif command -v singularity >/dev/null 2>&1; then
        singularity version
    else
        printf 'No Apptainer or Singularity command is currently available.\n'
    fi

    print_heading "Audit result"
    printf '%s\n' \
        "Read-only inspection completed." \
        "Preserve this terminal output when selecting mpiexec or srun mode."
}

main "$@"
