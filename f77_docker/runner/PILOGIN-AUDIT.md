# pilogin compatibility audit

Audit date: 2026-07-25

The inspection was read-only. It queried versions and runtime configuration on
the login service and did not create files, submit jobs, or change cluster
state. Private installation prefixes are intentionally omitted from this
public document.

## Host and scheduler

| Component | pilogin result |
| --- | --- |
| Architecture | x86-64 |
| Kernel | 4.18.0 KOS5 |
| glibc API version | 2.28 |
| Slurm | 25.11.2 |
| `MpiDefault` | unset |
| `MpiParams` | unset |
| Apptainer | 1.5.2 |
| Singularity compatibility command | Apptainer 1.5.2 |
| Locked-memory limit on login nodes | 64 KiB |

With the current module stack loaded, `srun --mpi=list true` reports `pmi2`,
`none`, and `cray_shasta`. It does not report a Slurm PMIx plugin.

## Compiler and MPICH

| Property | pilogin | Container |
| --- | --- | --- |
| GNU Fortran | 4.8.5 | 4.8.5 |
| MPICH | 4.1.2 | 4.1.2 |
| MPICH ABI | 15:1:3 | 15:1:3 |
| Device | `ch4:ofi` | `ch4:ofi` |
| Configure options | prefix only | prefix only |
| C compiler | `gcc -std=gnu99 -O2` | `gcc -std=gnu99 -O2` |
| C++ compiler | `g++ -O2` | `g++ -O2` |
| F77/FC compiler | `gfortran -O2` | `gfortran -O2` |
| Process manager | PMI | PMI |
| Hydra Slurm launcher | available | available |
| Hydra Slurm resource manager | available | available |
| Topology library | embedded hwloc | embedded hwloc |
| Network layer | embedded libfabric | embedded libfabric |

The container MPICH configure command was reduced to `--prefix` only after
this audit so its automatically selected compiler flags match pilogin. Prefix
and source-tree paths necessarily differ and do not affect the ABI.

## Current module environment

The site-standard module command tested for this project is:

```text
module load gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0
```

It provides GCC/GFortran 12.3.0 and OpenMPI 4.1.6. The OpenMPI build uses
internal PMIx and Slurm support, but was configured `--without-pmi`. The
Slurm installation exposes PMI2 but no PMIx plugin.

Host compiler selection does not replace the image compiler. Apptainer still
uses GFortran 4.8.5 and MPICH 4.1.2 from the SIF for pipeline compilation and
runtime linkage.

## RDMA user-space packages

| Package | pilogin | Container runtime |
| --- | --- | --- |
| `libibverbs` | KOS5 37.2 | Rocky 8 rdma-core 48.0 |
| `librdmacm` | KOS5 37.2 | Rocky 8 rdma-core 48.0 |
| `libnl3` | KOS5 3.5.0 | Rocky 8 3.7.0 |

Both MPI libraries resolve the same required SONAMEs, including
`libibverbs.so.1`, `librdmacm.so.1`, `libnl-3.so.200`, and
`libnl-route-3.so.200`. This is necessary but not sufficient to guarantee the
compute-node fabric provider.

The first multi-node smoke test should use the image defaults. If it fails in
provider discovery or RDMA initialization, use administrator-confirmed
read-only paths in `HPC_EXTRA_BINDS`; do not replace the entire container
`/lib64`. `FI_PROVIDER` and `FI_PROVIDER_PATH` are exposed only for controlled
diagnostics.

## Deployment decisions

The original validated launch model is host-MPICH hybrid mode:

1. Slurm allocates the nodes and ranks.
2. The cluster MPICH 4.1.2 `mpiexec` uses its Hydra Slurm integration.
3. Hydra starts `run-apptainer.sh` once per rank.
4. Each containerized application rank initializes the matching MPICH 4.1.2
   library from the SIF.

The current pilogin module environment uses the second validated model:

1. The batch script first loads GCC 12.3.0 and OpenMPI 4.1.6.
2. It does not use the host OpenMPI `mpirun` for the MPICH application.
3. Slurm `srun --mpi=pmi2` starts one Apptainer command per rank.
4. The SIF's MPICH 4.1.2 library initializes through Slurm PMI2.
5. The runner removes only host `OMPI_MCA_mtl` and `OMPI_MCA_osc` overrides,
   preserving `PMI_*` and `SLURM_*`.

Host OpenMPI and container MPICH are not ABI-compatible MPI implementations.
The successful PMI2 mode must not be interpreted as permission to mix host
OpenMPI `mpirun` with a MPICH-linked executable.

## Compute-node validation result

The post-audit pilogin validation completed successfully:

- a two-node allocation launched four ranks, with two ranks on each node;
- MPI identity and MPI/LAPACK smoke programs completed;
- a representative f77pipeline run processed one exposure and five valid
  chips through all stages from preprocessing to catalogue combination;
- the job exited successfully with empty stderr and generated the expected
  result files.

The GCC 12.3.0/OpenMPI 4.1.6 module path was separately validated:

- native OpenMPI launched four ranks across two nodes;
- Slurm PMI2 launched the container MPICH identity test across two nodes;
- the complete containerized pipeline finished in 3 minutes 52 seconds of
  MPI runtime with exit code zero;
- all stages completed and the merged catalogue was regenerated;
- the MPI step used approximately 5.4 GiB peak memory.

The successful run also established two pilogin-specific requirements for
low-rank testing: request the multi-node allocation with `--exclusive`, and
set `SLURM_CPUS_PER_TASK=1` and `SLURM_TRES_PER_TASK=cpu=1` consistently before
the MPI launcher starts. See [README-CN.md](README-CN.md) for the exact
workflow and the dedicated pilogin wrappers.

This result applies to pilogin. A different cluster must repeat its own
compute-node smoke test because devices, provider selection, CPU affinity, and
cross-node networking are site-specific.
