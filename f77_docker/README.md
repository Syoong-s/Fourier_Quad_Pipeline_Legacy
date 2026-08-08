# f77pipeline build environment

This project builds a development image matching the pilogin toolchain used by
f77pipeline:

- GNU GCC, G++, and GFortran 4.8.5
- MPICH 4.1.2, built from source with the `ch4:ofi` device
- CFITSIO 4.3.1, built from source
- LAPACK 3.8.0 and its reference BLAS, built from source
- glibc 2.28 from the pinned Rocky Linux 8.10 base image

The image deliberately contains no f77pipeline source code or observation data.
The source tree, catalogues, calibration files, and processing directory are
bind-mounted when the container starts.

For non-root, multi-node Slurm clusters, keep this Dockerfile as the image
build source and replace Docker Compose at runtime with Apptainer/Singularity.
See the [English HPC reference](runner/README.md) or the
[beginner-oriented Chinese tutorial](runner/README-CN.md).

## Reproducibility and supported platform

The Dockerfile downloads all source archives from the official upstream
locations listed in [SOURCES.md](SOURCES.md). Every archive is verified against
a committed SHA-256 digest before extraction. No `sources/` directory is used
by the build.

The validated target is `linux/amd64`. On an x86-64 Linux host with Docker,
build the image with:

```text
docker build --platform linux/amd64 --build-arg BUILD_JOBS=4 -t f77pipeline-dev:gnu4.8.5 .
```

The initial build compiles GCC and the complete scientific stack from source,
so it can take a substantial amount of time and needs several gigabytes of free
Docker storage. Later builds can reuse Docker's layer cache.

## Verify the image

```text
bash scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

The verification checks exact compiler and library versions, MPI's `ch4:ofi`
device, the shared-library loader, a two-rank MPI launch, real MPI
Fortran/LAPACK and CFITSIO compile-and-run tests, and the absence of pipeline or
dependency source files in the final image.

## Mount f77pipeline and its data

Copy the example environment file and edit the host paths:

```text
cp .env.example .env
```

The three container catalogue paths must be identical to the values compiled
into `para.inc`:

- `ASTROMETRY_CAT_CONTAINER` must match `ASTROMETRY_CAT`.
- `SOURCE_CAT_CONTAINER` must match `SOURCE_CAT`.
- `FLAT_PATH_CONTAINER` must match `FLAT_PATH`.

`PROCESS_DATA_HOST` is the processing directory containing `EXPO_LIST` and the
FITS images. Its writable container path is `PROCESS_DATA_CONTAINER`.

Build and enter the development container:

```text
docker compose build
docker compose run --rm FourierQuad-F77
```

Inside the container, f77pipeline is available at `/workspace/f77`. Recompile it
there whenever a processing stage changes. All dependency libraries are in one
directory:

```text
/opt/f77stack/lib
```

Headers are in `/opt/f77stack/include`. The environment also exports
`LIBRARY_PATH`, `CPATH`, `LAPACK_LIB_DIR`, and `CFITSIO_LIB_DIR`, all pointing
to the unified installation as appropriate. Makefile link flags can therefore
use `-L/opt/f77stack/lib`.

## Publication check

Before committing the directory to GitHub, run:

```text
bash scripts/check-public-repo.sh
```

This rejects accidentally added source archives, oversized files, private
cluster paths, and Dockerfile instructions that copy a local source directory.

## Publish the OCI image for HPC

The repository workflow at
[`../.github/workflows/f77pipeline-container.yml`](../.github/workflows/f77pipeline-container.yml)
builds the pinned `linux/amd64` image and publishes it to GHCR when manually
dispatched or when a `v*` tag is pushed. A normal HPC user can then create an immutable SIF with
`runner/pull-sif.sh`; no Docker daemon or root privilege is required at runtime.

Before selecting the MPI launcher, run the read-only
`runner/inspect-cluster-mpi.sh` audit. MPI version equality alone does not prove
compatibility with the site's Slurm PMI plugin or high-speed interconnect.
pilogin also has validated wrappers that read the configured GCC 12.3.0 and
OpenMPI 4.1.6 modules from `runner/f77pipeline.env`, then use `srun --mpi=pmi2`
to launch the MPICH 4.1.2 application inside the SIF. The host OpenMPI
`mpirun` is not mixed with the MPICH application.
See [runner/PILOGIN-AUDIT.md](runner/PILOGIN-AUDIT.md).

## License

Repository-authored files are licensed under the MIT License. Downloaded
dependencies and the GCC compatibility patch remain subject to their upstream
licenses; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
