# Fourier_Quad F77 Pipeline

MPI-parallel Fortran implementation of the Fourier_Quad weak-lensing shear
measurement pipeline. It processes CCD exposures from preprocessing and
astrometric calibration through PSF modeling, shear estimation, and final
catalog generation.

> 中文使用说明见 [F77_GUIDE_CN.md](F77_GUIDE_CN.md). The complete English
> reference is [F77_GUIDE.md](F77_GUIDE.md).

## Choose a version

| Version | Recommended for |
|---|---|
| `f77/` | Full pipeline, including optional multi-scale PCA PSF reconstruction |
| `f77_Lite/` | Smaller frozen production variant when PCA reconstruction is not needed |

New users who do not require `PSF_Ms=1` can start with `f77_Lite`. Both
versions build the executable `Fourier_Quad_Pipe` and use the same command-line
interface.

Clone the repository with the URL shown by GitHub's **Code** button, or download
the selected source and Docker archives from the repository's **Releases**
page. The Docker quick start requires both the source tree and `f77_docker/`.

## Quick start with Docker

The published image is the easiest way to obtain the legacy GNU Fortran, MPICH,
CFITSIO, and LAPACK stack. The image contains the toolchain only; your selected
source directory and data are mounted when the container starts.

Requirements: Docker with Compose support and an x86-64 Linux host.

### 1. Prepare an isolated processing directory

Do not run the pipeline directly against an irreplaceable raw-data tree. Output
paths are derived from the FITS paths, and existing intermediate files may be
replaced. Create a writable processing tree and place or link the input FITS and
DQ-mask files below it.

You may use two scripts to prepare the input data:
- [`Decompose Images`](Tools/gen_cat_mpi_ver2.6.py)
- [`Decompose Masks`](Tools/decom_mask_fz_v2.6.py)

After preparation, the input data must be organized into a directory tree with 
the top-level exposure list contains one exposure-list path and chip count per
line:

```text
"/data/DataProcess/run/g2019/stamps/exposure_001.list" 5
```

Each per-exposure list contains one Science FITS path per line:

```text
/data/DataProcess/run/g2019/science/exposure_001_1.fits
/data/DataProcess/run/g2019/science/exposure_001_2.fits
```

When using a container, every path inside these lists must be a container path,
not its host equivalent. See the [HPC tutorial](f77_docker/runner/README-CN.md)
for a complete safe data-tree example.

### 2. Configure the pipeline

In the selected source directory, edit `para.inc` before compiling:

- set `ASTROMETRY_CAT`, `SOURCE_CAT`, and `FLAT_PATH` to the paths visible
  inside the container;
- select stages with `PROCESS_stage`;
- review the image geometry, mask/flat switches, catalog columns, and PSF
  parameters.

The default `PROCESS_stage = 2*3*5*7*11*13*17*19*23` runs all nine stages.
Remove a stage's prime factor to skip that stage; the mapping is documented in
[F77_GUIDE.md](F77_GUIDE.md#pipeline-stages).

### 3. Configure the container mounts

From the repository root:

```bash
cd f77_docker
cp .env.example .env
```

Edit `.env` and set:

- `IMAGE_NAME=ghcr.io/syoong-s/fourier_quad_pipeline_legacy:latest`;
- `F77_SOURCE_HOST` to the absolute host path of `f77/` or `f77_Lite/`;
- the astrometry catalog, source catalog, flat-field, and processing-data host
  paths;
- `HOST_UID` and `HOST_GID` to your numeric user and group IDs when necessary.

The `*_CONTAINER` catalog paths in `.env` must exactly match the strings in
`para.inc`.

### 4. Compile and run

Pull the published image and enter the container:

```bash
docker compose pull
docker compose run --rm FourierQuad-F77
```

Inside the container:

```bash
make clean
make LAPACK_LIB_DIR=/opt/f77stack/lib \
  CFITSIO_LIB_DIR=/opt/f77stack/lib \
  CFITSIO_LIB=/opt/f77stack/lib/libcfitsio.so
mpiexec -n 4 ./Fourier_Quad_Pipe \
  /data/DataProcess/expo_list.list
```

Adjust the rank count for the number of exposures and available resources. The
pipeline dynamically assigns exposures to MPI ranks; chips within one exposure
are processed sequentially.

### 5. Find the results

For each dataset, the main products are written below the directory derived
from its FITS paths:

- `result/<EXPOSURE>_all.cat`: combined per-exposure shear catalog;
- `expo_info.dat`: exposure-level diagnostics;
- `astrometry/`, `stamps/`, `rescale/`, and related directories: intermediate
  products used by later stages.

Check the program output for missing catalog tiles, invalid chips, and PSF
quality failures before using the shear catalogs scientifically.

## Build from source

Use this route when a compatible scientific stack is already installed.

Required software:

- GNU Fortran and MPI with `mpif77`;
- CFITSIO;
- LAPACK and BLAS;
- GNU Make.

After configuring `para.inc`, build and run from the repository root:

```bash
export SCIENCE_PREFIX=/path/to/scientific-stack
make -C f77 \
  LAPACK_LIB_DIR="$SCIENCE_PREFIX/lib" \
  CFITSIO_LIB_DIR="$SCIENCE_PREFIX/lib"
mpirun -np 4 ./f77/Fourier_Quad_Pipe \
  /path/to/expo_list.list
```

Replace `f77` with `f77_Lite` for the Lite variant. If CFITSIO does not provide
`libcfitsio.so` at that location, pass its full path as `CFITSIO_LIB`.

## Run on a Slurm cluster

An Apptainer/Singularity runner is included for systems where Docker is not
available on compute nodes:

```bash
cd f77_docker/runner
cp f77pipeline.env.example f77pipeline.env
```

Edit `f77pipeline.env` with the GHCR image, SIF destination, source/data mounts,
MPI launch mode, and site modules. Then validate in order:

*Note:* GHCR image is available as `ghcr.io/syoong-s/fourier_quad_pipeline_legacy:latest`.

```bash
bash pull-sif.sh
bash inspect-cluster-mpi.sh
bash run-apptainer.sh --check
sbatch mpi-smoke-test.slurm
sbatch f77pipeline.slurm
```

Pin `OCI_IMAGE_URI` by digest for production runs. Do not combine a host
OpenMPI launcher with the MPICH application unless compatibility has been
validated for the site. See [the runner guide](f77_docker/runner/README.md) for
configuration and multi-node launch details.

## Manual for AI

Here also provides a manual wraped as a codex/claude code plugin, which gives agent abilities to help you build environment, switch parameters and run pipeline. For more details, please refer to the [Codex/Claude Code Plugin](https://github.com/Syoong-s/FQLegacyAIManual).

## Documentation

- [F77 source, parameters, Docker, and runner guide](F77_GUIDE.md)
- [中文完整指南](F77_GUIDE_CN.md)
- [Container environment](f77_docker/README.md)
- [HPC runner](f77_docker/runner/README.md)

## License

Distributed under the [MIT License](LICENSE).
