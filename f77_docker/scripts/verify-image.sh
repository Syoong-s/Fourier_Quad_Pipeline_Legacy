#!/usr/bin/env bash

set -euo pipefail

IMAGE_NAME="${1:-f77pipeline-dev:gnu4.8.5}"
REPOSITORY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ==========================================
# Function: Fail early when Docker is not usable by the current user
# Method: Query daemon metadata without changing daemon state
# ==========================================
require_docker() {
    command -v docker >/dev/null 2>&1
    docker info >/dev/null
}

# ==========================================
# Function: Verify compiler and library versions inside the image
# Method: Run version commands and assert exact required markers
# ==========================================
verify_versions() {
    docker run --rm "${IMAGE_NAME}" bash -lc '
        test "$(gcc -dumpversion)" = "4.8.5"
        test "$(gfortran -dumpversion)" = "4.8.5"
        test "$(getconf GNU_LIBC_VERSION)" = "glibc 2.28"
        mpichversion | grep -F "MPICH Version:" | grep -F "4.1.2"
        mpichversion | grep -F "MPICH ABI:" | grep -F "15:1:3"
        mpichversion | grep -F "MPICH Device:" | grep -F "ch4:ofi"
        mpiexec -info |
            grep -F "Launchers available:" | grep -F "slurm"
        mpiexec -info |
            grep -F "Resource management kernels available:" |
            grep -F "slurm"
        grep -F "Version: 4.3.1" /opt/f77stack/lib/pkgconfig/cfitsio.pc
    '
}

# ==========================================
# Function: Verify the unified dependency directory and runtime loader
# Method: Check required files, SONAMEs, dependencies, and a two-rank launch
# ==========================================
verify_stack() {
    docker run --rm "${IMAGE_NAME}" bash -lc '
        test -e /opt/f77stack/lib/libmpi.so.12
        test -e /opt/f77stack/lib/libmpifort.so.12
        test -e /opt/f77stack/lib/libcfitsio.so.10
        test -s /opt/f77stack/lib/liblapack.a
        test -s /opt/f77stack/lib/libblas.a
        test -d /workspace/f77
        test -d /data/catalogs/AstroDir
        test -d /data/catalogs/ExtSrcDir
        test -d /data/calib/FlatDir
        test -d /data/DataProcess
        readelf -d /opt/f77stack/lib/libcfitsio.so.10 |
            grep -F "Library soname: [libcfitsio.so.10]"
        ldd /opt/f77stack/lib/libmpifort.so.12 |
            grep -F "libgfortran.so.3"
        ldd /opt/f77stack/lib/libmpi.so.12 | grep -F "librdmacm.so.1"
        ldd /opt/f77stack/lib/libmpi.so.12 | grep -F "libibverbs.so.1"
        ! ldd /opt/f77stack/lib/libmpifort.so.12 | grep -F "not found"
        ! ldd /opt/f77stack/lib/libmpi.so.12 | grep -F "not found"
        mpiexec -n 2 /bin/true
    '
}

# ==========================================
# Function: Verify compiler wrappers and scientific-library linkage
# Method: Build and run small MPI/LAPACK and CFITSIO programs in a disposable container
# ==========================================
verify_compiler_linkage() {
    docker run --rm \
        --volume "${REPOSITORY_DIR}/tests:/tests:ro" \
        --volume "${REPOSITORY_DIR}/runner/tests:/runner-tests:ro" \
        "${IMAGE_NAME}" bash -lc '
        mkdir -p /tmp/f77pipeline-smoke
        mpif77 -O2 /runner-tests/mpi_lapack_smoke.f \
            -L/opt/f77stack/lib -llapack -lblas \
            -o /tmp/f77pipeline-smoke/mpi_lapack_smoke
        mpif77 -O2 /runner-tests/mpi_identity.f \
            -o /tmp/f77pipeline-smoke/mpi_identity
        gcc -O2 /tests/cfitsio_smoke.c \
            -I/opt/f77stack/include -L/opt/f77stack/lib -lcfitsio \
            -o /tmp/f77pipeline-smoke/cfitsio_smoke
        mpiexec -n 2 /tmp/f77pipeline-smoke/mpi_identity
        mpiexec -n 2 /tmp/f77pipeline-smoke/mpi_lapack_smoke
        /tmp/f77pipeline-smoke/cfitsio_smoke
    '
}

# ==========================================
# Function: Verify that pipeline and third-party build sources were excluded
# Method: Search final runtime locations for forbidden source and archive names
# ==========================================
verify_source_absence() {
    docker run --rm "${IMAGE_NAME}" bash -lc '
        test ! -e /workspace/f77/para.inc
        test ! -e /workspace/f77/cust_para.inc
        test ! -e /workspace/f77/Fourier_Quad_Pipe
        test -z "$(find /opt /workspace -xdev -type f \
            \( -name "para.inc" -o -name "cust_para.inc" \
               -o -name "Fourier_Quad_Pipe" -o -name "*.tar.gz" \
               -o -name "*.tar.bz2" \) -print -quit)"
    '
}

# ==========================================
# Function: Run the complete non-destructive image verification suite
# Method: Execute independent version, library, MPI, and source checks
# ==========================================
main() {
    require_docker
    verify_versions
    verify_stack
    verify_compiler_linkage
    verify_source_absence
    printf 'Image verification passed: %s\n' "${IMAGE_NAME}"
}

main "$@"
