# F77 dataset initializer

`init_program.py` replaces the former `gencat.py`, `decom_mask.py`, and
`server.slurm`. It is one Python entry point with embedded `#SBATCH` directives,
so the same file supports local serial runs, MPI runs, and direct Slurm batch
submission.

The program reads the original science and DQ `.fits.fz` repositories in place.
It never copies, moves, or deletes those archives. Outputs use the current
`cpp_Standard`/`cpp_Lite` and `f77`/`f77_Lite` contract:

```text
<output-root>/
├── <target>/
│   ├── science/<exposure>/<exposure>_<sequence>.fits
│   ├── dqmask/<exposure>/<exposure>_<CCDNUM>.fits
│   ├── expolists/<exposure>.list
│   ├── result/
│   ├── stamps/<product>/[<exposure>/]
│   └── astrometry/<product>/[<exposure>/]
├── expo_<target>.list
├── fits_<target>.list
└── init_<target>_manifest.json
```

Science outputs retain `CCDNUM` in their FITS headers; their basename suffix is
the C++ two-dimensional-HDU occurrence number. DQ outputs use physical
`CCDNUM`, and DQ exposure stems map `ood` to `ooi`.

## Runtime requirements

- Python 3.10 or newer
- Astropy
- mpi4py linked to the MPI used by the launcher
- NumPy (installed as an Astropy dependency)
- Slurm `srun` only for direct `sbatch` execution

The verified WSL environment uses Python 3.12.13, Astropy 8.0.0, mpi4py 4.1.2,
NumPy 2.5.0, and Open MPI 5.0.10. On the cluster, load the site's Python/MPI
modules or activate an equivalent virtual environment before submitting.

## Portable commands

Inspect all options without opening data:

```bash
python init_program.py --help
```

Run serially (MPI size one):

```bash
python init_program.py \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work \
  --target gband --prefix c4d_ --contains v1 \
  --existing fail
```

Run under MPI:

```bash
mpirun -np 16 python init_program.py \
  --science-root /data/archive/science \
  --dq-root /data/archive/dq \
  --output-root /data/work \
  --target gband --prefix c4d_ --contains v1 \
  --existing resume
```

Submit the same file to Slurm. The batch process automatically starts one
`srun` step using the allocation requested by the directives at the top of the
Python file:

```bash
sbatch init_program.py \
  --science-root /lustre/archive/science \
  --dq-root /lustre/archive/dq \
  --output-root /lustre/work \
  --target gband --prefix c4d_ --contains v1 \
  --existing resume
```

Edit the embedded `#SBATCH` lines when another partition, task count, node
layout, log path, or time limit is required. `--no-srun` intentionally keeps a
batch allocation serial for diagnosis.

## Existing-output modes

- `fail`: reject an archive if any planned chip output already exists.
- `resume`: validate and reuse good outputs, then regenerate missing or invalid
  chips.
- `overwrite`: regenerate every planned chip and atomically replace it.

Each chip is first written below `<target>/.fq_init_tmp/<run>/` and then moved
into place on the same filesystem. Lists and the manifest are also atomically
replaced. Any failed archive produces a nonzero program exit and is recorded in
the manifest; successfully initialized exposures remain usable for diagnosis or
a later `--existing resume` run.
