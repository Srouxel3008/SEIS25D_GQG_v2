# Building and running on an HPC cluster

This workflow is intended for the Marmousi VSP example and other cases that
need more memory than a workstation. Cluster module names, account,
partition, CPU count, memory policy and time limits vary between institutions
and must be adapted locally.

The Marmousi example requires more than 170 GB of available RAM. Verify the
memory actually available to the scheduled job and retain additional headroom
for the operating system, MPI and numerical libraries. A sufficiently equipped
high-memory workstation can also run it, but a multi-node cluster provides more
CPU cores for the expensive Fréchet-derivative calculations. The supplied
six-node allocation targets a runtime below the 48-hour production limit on
the reference cluster.

Before preparing a cluster case, read [Input files](INPUT_FORMATS.md) for the
main configuration, model, grid/topography, acquisition-geometry and external
spectral-data formats.

## Files

- `src/Makefile.hpc` supplies visible HPC build defaults and includes the
  common `src/Makefile` so the source list and dependencies remain identical.
- `scripts/run_slurm.sh` is a sanitized, single-case Slurm example.

## 1. Load the cluster toolchain

List the modules available on your system, then load the Intel compiler, Intel
MPI and oneMKL modules. For example only:

```bash
module purge
module load intel
module load intel-mpi
module load mkl
# module load mumps  # uncomment only when building with USE_MUMPS=1
```

The exact names must be replaced with those used by the target cluster.

## 2. Build

### Release build

From the repository root, build the optimized version used for normal
production runs:

```bash
make -C src -f Makefile.hpc CONFIG=release
```

`release` is the default, so omitting `CONFIG=release` produces the same build:

```bash
make -C src -f Makefile.hpc
```

The default build uses PARDISO, matching the published Marmousi input. The
executable is written to `bin/Seis2D_GQG_FWI`.

### Debug build

For troubleshooting, compile with runtime checks and debugging information:

```bash
make -C src -f Makefile.hpc CONFIG=debug
```

Debug builds are slower and should not be used for normal Marmousi production
runs. Before changing between release and debug modes, remove the objects from
the previous build and then rebuild:

```bash
make -C src -f Makefile.hpc clean
make -C src -f Makefile.hpc CONFIG=debug
```

The Slurm launcher does not compile the program or choose a build mode. It runs
the executable currently present in `bin/Seis2D_GQG_FWI`, so compile the
intended configuration before submitting the job.

### Optional MUMPS build

To try MUMPS, load the cluster's MUMPS module and build with:

```bash
make -C src -f Makefile.hpc USE_MUMPS=1
```

MUMPS can also be combined with a debug build:

```bash
make -C src -f Makefile.hpc CONFIG=debug USE_MUMPS=1
```

The solver must then be selected in `MainInput.inp`. The MUMPS module must
expose its headers and libraries. If it does not define `MUMPS_DIR`, pass the
installation root to `make`, and override `MUMPS_LIBS` when the cluster uses
different dependency or library names. The editable settings and link pathway
are shown at the top of `src/Makefile.hpc`.

## 3. Configure Slurm

Open `scripts/run_slurm.sh` and edit the first two settings sections:

- uncomment and set `--partition` and `--account` when required;
- choose a time limit accepted by the cluster;
- retain six nodes for the published Marmousi runtime target, or adapt after
  benchmarking on the destination cluster;
- ensure that only one MPI rank per node is used;
- set `--cpus-per-task` to the physical/logical CPUs assigned to each rank;
- request a node with enough memory for Marmousi;
- replace the generic module names with the local names.

The supplied launcher selects `examples/Marmousi_VSP/MainInput_C33.inp` by
default. To run another published inversion, change the visible `MainInput`
filename near the top of the launcher. This single assignment can later be
replaced by a `SLURM_ARRAY_TASK_ID` mapping when several configurations must be
submitted together.

The example requests all available memory on each allocated node with
`--mem=0`; some clusters require an explicit per-node value instead.

## 4. Submit

From the repository root:

```bash
sbatch scripts/run_slurm.sh
```

The script locates the repository automatically and selects the C33
configuration under `examples/Marmousi_VSP`. Results are written to:

```text
runs/Marmousi_VSP/slurm_JOB_ID/
```

The result directory contains `environment.txt`, `configuration_summary.txt`,
and all program outputs. The configuration summary preserves the selected
MainInput filename, checksum and complete contents. Slurm stdout and stderr are
copied into the directory every five minutes and again when the job exits.

Monitor with the commands provided by the cluster, commonly:

```bash
squeue -u "$USER"
tail -f slurm-JOB_ID.out
```

Cancel a job with:

```bash
scancel JOB_ID
```

## Other datasets and filenames

The Slurm launcher keeps the five input roles together in one visible block.
Edit `INPUT_DIR` and the required filenames directly. Optional roles may be
assigned empty strings while their positional slots remain intact.

See [Input files](INPUT_FORMATS.md) for the required contents, units, indexing
and supported combinations of these files.

## Before relying on a new cluster

Perform a short Small TTI build/run first. Confirm module compatibility,
thread placement, memory reporting, Slurm integration and numerical output
before allocating resources to Marmousi.
