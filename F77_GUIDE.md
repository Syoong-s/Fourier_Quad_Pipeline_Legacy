# Fourier_Quad F77 Pipeline Guide

Source structure, configuration, building, Docker environment, and HPC runner
for the Fortran 77 (`f77` / `f77_Lite`) pipeline. For the project overview and
the C++ pipeline see [`README.md`](README.md) and [`CPP_GUIDE.md`](CPP_GUIDE.md).

> **中文文档**：请参阅 [F77_GUIDE_CN.md](F77_GUIDE_CN.md)

---

## Source Structure

### Source directories

#### `f77/` — Full Fortran 77 pipeline

| File | Description |
|---|---|
| `main.f` | Main program entry point. Initializes MPI, reads the exposure list, and dispatches pipeline stages. |
| `para.inc` | Master parameter file. Controls image dimensions, stage selection (`PROCESS_stage`), catalog paths, PSF order, stamp size, thresholds, and catalog column indices. |
| `path_layout.inc` | WFST-compatible product-directory contract used by the full `f77/` pipeline. |
| `cust_para.inc` | Custom parameters for CCD geometry and PCA PSF decomposition. |
| `sig_para.inc` | F6 mode-bar noise-plane estimator parameters. |
| `pre_process.f` | **Stage 1**: flat-field, mask handling, background/noise estimation. |
| `proc_astrometry.f` | **Stage 2**: astrometric calibration using Gaia reference catalog. |
| `proc_source.f` | **Stage 3**: source detection and stamp extraction. |
| `proc_FFT_st1.f` | **Stage 4**: first-stage Fourier transform of galaxy stamps. |
| `proc_PSF.f` | **Stage 5**: PSF modeling from stellar stamps (local polynomial fit). |
| `proc_psfreconsV2.f` | PSF PCA reconstruction (multi-scale, `PSF_Ms=1` only). |
| `00_psf_module.f` | Fortran 90 module for global PSF PCA storage. |
| `proc_FFT_st2.f` | **Stage 6**: second-stage Fourier transform. |
| `proc_shear.f` | **Stage 7**: Fourier\_Quad shear estimation. |
| `proc_info.f` | **Stage 8**: per-exposure statistics collection. |
| `proc_combine_shear_catalog.f` | **Stage 9**: shear catalog combination and calibration. |
| `astrometry_calib.f` | Astrometric calibration utility routines. |
| `ex_star.f` | Star extraction and classification. |
| `rw_fits_image.f` | FITS image read/write routines (via CFITSIO). |
| `universal.f` | General-purpose utility functions. |
| `univ_imag_proc.f` | Image processing utilities. |
| `press.f` | Numerical Recipes routines (sorting, statistics, interpolation). |
| `mpi_routines.f` | MPI distribution helper (`mpi_distribute`). |
| `FFTPACK.f` | FFTPACK 5.1 FFT library. |
| `Makefile` | Build file. Uses `mpif77`, links against LAPACK, BLAS, and CFITSIO. |

#### `f77_Lite/` — Simplified Fortran 77 pipeline

Closely follows `f77/`, but `00_psf_module.f` is absent (PCA PSF
reconstruction is removed). All eight compile-time switches are frozen to
production values and dead-code branches are physically removed. The
WFST-compatible I/O migration documented below applies to `f77/`; `f77_Lite/`
continues to use its existing legacy layout.


## Pipeline Stages

The pipeline consists of 9 stages. Stage execution is controlled by the
`PROCESS_stage` parameter (defined in `para.inc` / `LensingConfig.hpp`), which is
the product of prime factors. A stage runs when `PROCESS_stage` is divisible by
its prime. The default value `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23` enables
all stages.

| Stage | Prime | Function | Description |
|---|---|---|---|
| 1 | 2 | `pre_process` / `PreProcess` | Read FITS images, apply flat-field and mask corrections, estimate background noise (F6 mode-bar estimator). |
| 2 | 3 | `proc_astrometry` / `Astrometry` | Astrometric calibration using Gaia reference catalog; WCS fitting. |
| 3 | 5 | `proc_source` / `SourceExtractor` | Source detection, deblending, and stamp extraction. |
| 4 | 7 | `proc_FFT_st1` / `FourierTransformSt1` | First-stage Fourier transform of galaxy stamps. |
| 5 | 11 | `proc_PSF` / `PSFModel` | PSF modeling from stellar stamps via local polynomial fitting. Optional PCA reconstruction (`PSF_Ms=1`). |
| 6 | 13 | `proc_FFT_st2` / `FourierTransformSt2` | Second-stage Fourier transform. |
| 7 | 17 | `proc_shear` / `ShearMeasurement` | Fourier\_Quad shear estimation from fourth-order Fourier moments. |
| 8 | 19 | `proc_info` / `ExposureInfo` | Collect per-exposure statistics (PSF FWHM, star count, etc.). |
| 9 | 23 | `proc_combine_shear_catalog` / `CatalogCombiner` | Combine shear catalogs across exposures and apply calibration corrections. |

To disable a stage, divide `PROCESS_stage` by its prime factor. For example,
setting `PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19` (omit 23) skips
catalog combination.

---


## Configuration

### Fortran 77 (`para.inc`, `cust_para.inc`, `sig_para.inc`)

All compile-time parameters live in three include files:

- **`para.inc`** — Image dimensions (`npx`, `npy`), stage control
  (`PROCESS_stage`), catalog paths, PSF parameters, stamp size, detection
  thresholds, and catalog column indices.
- **`cust_para.inc`** — CCD geometry (`chipnx`, `chipny`, `Camera_ccd_num`) and
  PCA PSF decomposition parameters.
- **`sig_para.inc`** — F6 mode-bar noise-plane estimator parameters.

> Paths such as `ASTROMETRY_CAT`, `SOURCE_CAT`, and `FLAT_PATH` must match the
> container mount paths defined in the Docker `.env` files.

`SOURCE_CAT_TILE_PREFIX` controls the external source-catalog tile basename.
Its default is `extern_`, matching the WFST initializer and `process_extcat`.


## Dataset Structure

```text
<dataset>/
├── science/<exposure>/<chip>.fits
├── dqmask/<exposure>/<exposure>_<CCDNUM>.fits
├── expolists/<exposure>.list
├── astrometry/{Head,dat_Chk,dat_Astro}/...
├── stamps/
│   ├── {Norm,cat_Orig,dat_StarCanInfo,fits_StarCan}/<exposure>/...
│   ├── {fits_StarCanN,fits_StarCanP,dat_SrcInfo,fits_Src}/<exposure>/...
│   ├── {fits_Noise,fits_SrcP,dat_PsfFit,fits_PsfLocal}/<exposure>/...
│   ├── {dat_Shear,dat_StarXY,fits_PsfResi}/<exposure>/...
│   └── {dat_StarInfo,fits_StarP,fits_PsfSrc,dat_ExpoInfo}/...
└── result/<exposure>_all.cat
```


## Building from Source

### Fortran 77

Prerequisites: `gfortran`, `mpif77` (MPICH), CFITSIO, LAPACK/BLAS.

```bash
cd f77          # or f77_Lite
# Edit para.inc to set your catalog paths and parameters
make
# Executable: ./Fourier_Quad_Pipe
```

The Makefile accepts override variables:

```bash
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
```

Local verification for this I/O migration used GNU Fortran 15.2.0 through
Open MPI 5.0.10, CFITSIO 4.6.4, and LAPACK/BLAS 3.11.0. Modern GNU Fortran
requires the compatibility override below for pre-existing FFTPACK implicit
argument mismatches:

```bash
export SCIENCE_PREFIX=/path/to/scientific-stack
make FC=mpifort \
     FFLAGS='-mcmodel=medium -w -fallow-argument-mismatch' \
     LAPACK_LIB_DIR="$SCIENCE_PREFIX/lib" \
     CFITSIO_LIB_DIR="$SCIENCE_PREFIX/lib"
```

On the production cluster, load the site MPI/GNU Fortran, CFITSIO, and
LAPACK/BLAS modules (or use the documented container stack) and keep the same
relative build and run commands.


### Running the pipeline

```bash
mpirun -np <N> ./Fourier_Quad_Pipe <EXPO_LIST>             
```

---


## Docker Environment


The Docker environments provide a reproducible build toolchain without requiring
manual installation of compilers and scientific libraries.


### Docker directories

#### `f77_docker/`

Builds a development image matching the pilogin cluster toolchain:

| Component | Version |
|---|---|
| Base image | Rocky Linux 8.10 (glibc 2.28) |
| GCC / G++ / GFortran | 4.8.5 (compiled from source) |
| MPICH | 4.1.2 (ch4:ofi device, Slurm launcher) |
| CFITSIO | 4.3.1 |
| LAPACK + reference BLAS | 3.8.0 |

Key files: `Dockerfile`, `compose.yaml`, `.env.example`, `scripts/verify-image.sh`,
`scripts/check-public-repo.sh`, `runner/` (HPC deployment scripts), `config/`,
`patches/`, `checksums.sha256`, `SOURCES.md`, `THIRD_PARTY_NOTICES.md`.


### Quick start

```bash
cd f77_docker
cp .env.example .env
# Edit .env: set F77_SOURCE_HOST, catalog paths, and PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-F77
```

Inside the container:

```bash
cd /workspace/f77
make
mpirun -np 4 ./Fourier_Quad_Pipe /data/DataProcess/expo_list.list
```


### Verifying an image

Each Docker directory includes a verification script:

```bash
bash f77_docker/scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

For detailed Docker environment documentation, see:
- [`f77_docker/README.md`](f77_docker/README.md) / [`f77_docker/README-CN.md`](f77_docker/README-CN.md)
---


## HPC Deployment

Docker environments include `runner/` directories with Slurm/Apptainer
deployment scripts. The typical workflow is:

1. Pull the GHCR image and convert to a SIF:
   ```bash
   bash f77_docker/runner/pull-sif.sh    
   ```

2. Configure the environment:
   ```bash
   cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
   # Edit paths for your cluster
   ```

3. Audit MPI compatibility (read-only):
   ```bash
   bash f77_docker/runner/inspect-cluster-mpi.sh
   ```

4. Submit the pipeline:
   ```bash
   sbatch f77_docker/runner/f77pipeline.slurm
   ```

For detailed HPC runner documentation, see:
- [`f77_docker/runner/README.md`](f77_docker/runner/README.md) / [`f77_docker/runner/README-CN.md`](f77_docker/runner/README-CN.md)
---
