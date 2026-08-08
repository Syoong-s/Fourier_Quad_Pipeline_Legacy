# f77pipeline on a non-root Slurm HPC cluster

This directory replaces Docker Compose at HPC runtime. The Dockerfile remains
the single toolchain build source, while Apptainer or Singularity converts the
published OCI image into an immutable SIF and Slurm launches one container
process per MPI rank.

For a command-by-command introduction to Docker, Apptainer, Slurm, safe data
layout, pilogin-specific settings, and troubleshooting, see the
[complete Chinese tutorial](README-CN.md).

The production target remains:

- Linux x86-64
- GNU GCC and GFortran 4.8.5
- MPICH 4.1.2
- CFITSIO 4.3.1
- LAPACK 3.8.0 and reference BLAS
- Bash 4 or newer for indexed-array configuration
- Apptainer or Singularity supplied by the cluster
- Slurm with a process manager compatible with the cluster MPICH build

No f77pipeline source or observation data is stored in the image.

## Runtime model

Docker Compose is retained for local development only. On the cluster the
equivalent flow is:

1. Build and publish the OCI image outside the cluster.
2. Pull an immutable SIF as an unprivileged cluster user.
3. Bind the source, catalogues, calibration files, and processing data.
4. Compile the mounted source once on the allocated batch node.
5. Use host MPICH `mpiexec` or a validated Slurm `srun` PMI mode to launch one
   `apptainer exec` process per MPI rank.

Apptainer runs with the submitting user's real UID and GID. Docker
`HOST_UID`, `HOST_GID`, `USER`, service names, and container names do not apply
to an HPC job.

## Files

- `f77pipeline.env.example`: paths, host modules, default EXPO_LIST, and MPI selection.
- `f77pipeline.pilogin-openmpi.env.example`: validated pilogin module/PMI2 settings.
- `site-env.example.sh`: optional advanced initialization before module loading.
- `inspect-cluster-mpi.sh`: read-only cluster environment audit.
- `PILOGIN-AUDIT.md`: recorded pilogin/container compatibility comparison.
- `pull-sif.sh`: non-root OCI-to-SIF acquisition without overwrite.
- `run-apptainer.sh`: Compose-equivalent binds and container execution.
- `mpi-smoke-test.slurm`: two-node rank placement and MPI/LAPACK test.
- `mpi-smoke-test-pilogin-openmpi.slurm`: pilogin resource/PMI2 wrapper.
- `f77pipeline.slurm`: compile-once and multi-rank pipeline job template.
- `f77pipeline-pilogin-openmpi.slurm`: pilogin resource/PMI2 pipeline wrapper.

## 1. Inspect the cluster without changing it

Load the site's intended compiler/MPI and Apptainer modules, if required, then
run:

```bash
bash inspect-cluster-mpi.sh
```

Keep the complete terminal output. In particular, compare:

- `mpichversion`, including configure options and the CH4 device;
- `ompi_info` and the OpenMPI wrapper configuration when OpenMPI is loaded;
- `mpicc -show` and `mpifort -show`;
- `srun --mpi=list` and `MpiDefault`;
- libfabric providers and InfiniBand devices;
- Apptainer or Singularity version.

The pilogin read-only audit performed for this project found GNU Fortran 4.8.5,
glibc 2.28, MPICH 4.1.2 ABI 15:1:3 with `ch4:ofi`, and Hydra support for the
Slurm launcher and resource manager. The image has the same MPICH version,
ABI, device, process manager, Hydra launchers, and resource managers. It has
subsequently passed a two-node, four-rank smoke test and a complete
f77pipeline data run on pilogin. Other clusters still require their own
compute-node validation because version and build-option agreement alone does
not prove high-speed-network compatibility.

A second pilogin validation loaded GCC 12.3.0 and OpenMPI 4.1.6, then used
Slurm `srun --mpi=pmi2` to start the MPICH 4.1.2 application in the SIF. A
two-node four-rank smoke test and a complete pipeline run succeeded. This does
not make OpenMPI and MPICH interchangeable: host OpenMPI `mpirun` must not
launch the MPICH-linked application.

The audit also found KOS5 `libibverbs` and `librdmacm` 37.2 on pilogin, while
the pinned Rocky 8.10 runtime currently installs the same SONAMEs from
rdma-core 48.0. Do not blindly replace the whole container `/lib64`. If a
compute-node test requires host provider libraries or configuration, add only
the administrator-recommended read-only mappings to `HPC_EXTRA_BINDS`.
`FI_PROVIDER` and `FI_PROVIDER_PATH` are available for controlled diagnostics.

## 2. Configure shared paths

On the cluster, upload this complete `runner/` directory and enter it before
running any command:

```bash
cd /shared/project/f77pipeline/runner
cp f77pipeline.env.example f77pipeline.env
```

When upgrading an existing `f77pipeline.env`, add
`F77_SOURCE_CONTAINER`, `F77_EXPO_LIST_CONTAINER`, and the Bash array
`HPC_MODULES`; the jobs intentionally reject an incomplete older file instead
of silently restoring duplicated script defaults.

Edit every host path. All of the following must be visible at the same absolute
path from every allocated node:

- `F77_SIF`
- `F77_SOURCE_HOST`
- `ASTROMETRY_CAT_HOST`
- `SOURCE_CAT_HOST`
- `FLAT_PATH_HOST`
- `PROCESS_DATA_HOST`
- the `runner/` directory itself

`f77pipeline.env` is the single routine HPC runtime configuration. It also
defines `F77_SOURCE_CONTAINER`, the catalogue/calibration/data container
destinations, `F77_EXPO_LIST_CONTAINER`, `F77_EXECUTABLE`, `HPC_MODULES`, and
the MPI launch settings. The default container paths intentionally match
`.env.example` and `compose.yaml`:

```text
/workspace/f77
/data/catalogs/AstroDir
/data/catalogs/ExtSrcDir
/data/calib/FlatDir
/data/DataProcess
```

`ASTROMETRY_CAT_CONTAINER`, `SOURCE_CAT_CONTAINER`, and
`FLAT_PATH_CONTAINER` must exactly match the corresponding strings compiled
into `para.inc`. The processing host directory must contain `EXPO_LIST` and
the FITS images. `F77_EXPO_LIST_CONTAINER` is used only when the batch script
receives no explicit application argument.

For the validated pilogin OpenMPI-module mode, copy
`f77pipeline.pilogin-openmpi.env.example` instead. Its `srun`/PMI2 selection
decouples process launch from the host OpenMPI runtime. Both generic job
scripts load the `HPC_MODULES` array after reading the env file and before
using MPI or Apptainer. Leave it as `HPC_MODULES=()` when no modules are
needed.

Slurm allocation directives remain in the `.slurm` files by design. Edit those
files, or override them on the `sbatch` command line, when changing partition,
node/task counts, CPU, memory, time, account, QoS, or log paths.

## 3. Publish and acquire the SIF

The repository GitHub Actions workflow in
[`../../.github/workflows/f77pipeline-container.yml`](../../.github/workflows/f77pipeline-container.yml)
publishes the Dockerfile result to GHCR when manually dispatched or when a
`v*` tag is pushed. Update `OCI_IMAGE_URI` to the published package.

For a production run, prefer an OCI digest over a mutable tag:

```text
OCI_IMAGE_URI=ghcr.io/OWNER/REPOSITORY@sha256:IMAGE_DIGEST
```

Create the SIF as a normal user:

```bash
bash pull-sif.sh
```

Standalone `pull-sif.sh` and `run-apptainer.sh --check` commands do not load
`HPC_MODULES`, because `run-apptainer.sh` is also executed once per MPI rank.
If the container runtime itself is module-provided, load that runtime module
before these standalone commands or set `APPTAINER_BIN` to its executable.
Scheduled pipeline and smoke jobs always load `HPC_MODULES` themselves.

The script refuses to overwrite an existing SIF. If compute nodes cannot
access the internet, create the SIF on a permitted Linux host and transfer the
single SIF file to the configured shared path.

## 4. Verify one container instance

Run a version and bind check without a scheduler allocation:

```bash
bash run-apptainer.sh --check
```

This check is read-only except for normal runtime temporary files. It verifies
GNU 4.8.5, MPICH 4.1.2, CFITSIO 4.3.1, and the configured bind destinations.

## 5. Validate MPI before real data

Submit the two-node smoke test from inside `runner/`:

```bash
sbatch mpi-smoke-test.slurm
```

The smoke job creates one unique temporary directory below
`HPC_SHARED_SCRATCH_HOST`, compiles two small programs there, launches them
across the allocation, and removes that exact temporary directory on exit.
Its log must show multiple ranks distributed across the requested hostnames.
It loads the same `HPC_MODULES`, paths, runtime, and MPI settings as the
pipeline job.

The initial recommendation is:

```text
MPI_LAUNCH_MODE=mpiexec
MPI_LAUNCHER=mpiexec
```

This is Apptainer's hybrid model: cluster MPICH launches one container command
per rank, and the application uses compatible MPICH libraries inside the SIF.
The host and container MPI builds should be as close as possible.

Use `MPI_LAUNCH_MODE=srun` only after the audit confirms that the container
MPICH process management interface matches `SLURM_MPI_TYPE`. An incompatible
PMI setting can produce hangs, communication failures, or many independent
rank-zero processes.

On pilogin this test has now passed with:

```text
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
HPC_SCRUB_OPENMPI_ENV=1
```

Submit `mpi-smoke-test-pilogin-openmpi.slurm` when using the current GCC
12.3.0/OpenMPI 4.1.6 module stack. It preserves Slurm/PMI variables while
removing host OpenMPI transport overrides before entering the MPICH image. The
module versions come from `HPC_MODULES` in `f77pipeline.env`, not from the
Slurm wrapper.

## 6. Submit f77pipeline

The template defaults to the `cpu` partition, two nodes, eight total ranks,
four ranks per node, one CPU per rank, and 2 GiB per rank. Override all
dependent Slurm resource values together or edit a site-specific copy:

On pilogin with the current modules, submit
`f77pipeline-pilogin-openmpi.slurm`. It includes `--exclusive`, sets
`SLURM_CPUS_PER_TASK=1` and `SLURM_TRES_PER_TASK=cpu=1` consistently, and
delegates execution to the generic pipeline template. The generic template
loads `HPC_MODULES` and uses the MPI mode configured in `f77pipeline.env`. The
[Chinese tutorial](README-CN.md) records the tested procedure.

```bash
sbatch --nodes=4 --ntasks=64 --ntasks-per-node=16 \
    f77pipeline.slurm
```

With no batch-script arguments, the executable receives
`F77_EXPO_LIST_CONTAINER`. Explicit arguments override that configured default
and are preserved:

```bash
sbatch f77pipeline.slurm
sbatch f77pipeline.slurm EXPO_LIST
```

The job performs these operations in order:

1. verifies the image and binds;
2. optionally runs `make clean`;
3. compiles `F77_SOURCE_CONTAINER` once;
4. changes the container working directory to `PROCESS_DATA_CONTAINER`;
5. launches `F77_EXECUTABLE` with all allocated ranks.

Set `F77_MAKE_CLEAN=0` to retain existing objects. The source directory is
writable, so compilation products remain on the shared host filesystem.
Avoid concurrent jobs that clean or rebuild the same source tree. Use a
separate source/build directory for each concurrently compiled pipeline stage.

## High-speed interconnect notes

Apptainer uses the host network namespace by default; do not add an isolated
container network. For InfiniBand or another RDMA fabric, the container must
also see compatible devices, libfabric providers, and user-space driver
libraries. Exact requirements depend on the cluster audit and may require
additional read-only binds through `HPC_EXTRA_BINDS` or an image rebuild.
Keep it empty for the first smoke test, then change one provider-related
setting at a time based on compute-node evidence.

Passing a local two-rank Docker test is not sufficient. On a new cluster,
treat the image as multi-node ready only after:

1. single-rank execution;
2. multiple ranks on one node;
3. one rank on each of two nodes;
4. multiple ranks across two nodes;
5. fabric/provider and representative performance validation.

## Compile and run references

The image compiler wrappers are `mpif77` and `mpifort`. Scientific libraries
and headers remain unified at:

```text
/opt/f77stack/lib
/opt/f77stack/include
```

The wrapper scripts were syntax-tested locally with Bash 5.3.9. Production HPC
jobs require Bash 4 or newer, Slurm, and the modules listed in
`HPC_MODULES`; compilation and scientific libraries still come from the SIF.

The local development command remains `docker compose run --rm
FourierQuad-F77`. The cluster production command is `sbatch
f77pipeline.slurm`; no Windows or WSL wrapper belongs in a cluster job.
