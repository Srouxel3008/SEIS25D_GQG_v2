# Building and running on a workstation

This is the complete workflow for the bundled Small TTI example on native
Linux or Ubuntu under Windows Subsystem for Linux (WSL). The workstation build
uses one MPI rank and Intel oneMKL PARDISO. MUMPS is not required.
Small TTI requires a workstation with at least 16 GB of installed RAM and more
than 10 GB currently available to Linux/WSL. It has been run successfully up to
22 Hz on a 16 GB Lenovo Legion 5. The Marmousi example requires more than 170 GB
of available RAM. It can run on a sufficiently equipped high-memory
workstation, but a multi-node cluster is preferred because additional CPU cores
substantially accelerate the Fréchet-derivative calculations.
HPC build and run are covered separately in
[Building and running on an HPC cluster](HPC.md).

Before preparing a cluster case, read [Input files](INPUT_FORMATS.md) for the
main configuration, model, grid/topography, acquisition-geometry and external
spectral-data formats.

## 1. Install the required tools

### Windows: install WSL 2 and Ubuntu

Native Linux users can skip this subsection. On Windows 11, follow the
[official Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).
From an Administrator PowerShell window, the normal installation command is:

```powershell
wsl --install
```

Restart Windows when requested, open Ubuntu, and complete the initial Linux
username and password setup. An
[illustrated Windows 11 walkthrough](https://contabo.com/blog/how-to-install-wsl2-on-windows-11/)
is available as a supplementary guide.

Inside Ubuntu, install the standard build utilities:

```bash
sudo apt update
sudo apt install -y build-essential make
```

### Linux/WSL: install Intel oneAPI

Follow Intel's
[oneAPI installation guide for Linux](https://www.intel.com/content/www/us/en/docs/oneapi-toolkit/installation-guide-linux/latest/install-oneapi-toolkit-with-installer.html).
Allow approximately 14 GB for a complete installation; the exact requirement
depends on the oneAPI release and selected components. If the installer offers
component selection, install at least:

- [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html), providing `ifx`;
- [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html);
- [Intel oneAPI Threading Building Blocks (oneTBB)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb-download.html);
- [Intel oneAPI Collective Communications Library (oneCCL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/oneccl-download.html);
- [Intel MPI Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library-download.html), providing `mpiifx` and `mpirun`.

The documented workflow uses this Intel oneAPI toolchain because it is the
configuration currently tested. Other compiler and library combinations are
not yet supported.

Native Windows compilation is not supported. Windows users should build and
run from Ubuntu in WSL.

## 2. Open Linux and initialize oneAPI

On Windows, open the standalone Ubuntu application or a WSL terminal. A
repository downloaded at:

```text
C:\Users\path\SEIS25D_GQG_v2
```

is visible inside WSL at:

```text
/mnt/c/Users/path/SEIS25D_GQG_v2
```

These paths refer to the same files. Initialize Intel oneAPI after opening each
new Linux terminal:

```bash
source /opt/intel/oneapi/setvars.sh
```

## 3. Build

Move to the repository root, then build the normal optimized executable:

```bash
make -C src CONFIG=release USE_MUMPS=0
```

The executable is written to `bin/Seis2D_GQG_FWI`. Objects and generated
Fortran module files are written under `build/`.

For troubleshooting source changes, use the slower debug configuration:

```bash
make -C src clean
make -C src CONFIG=debug USE_MUMPS=0
```

The debug configuration enables compiler warnings, runtime checks,
floating-point traps and tracebacks. Clean before switching between debug and
release so that objects from different configurations are not mixed.

The available build settings are:

| Option | Values | Default | Meaning |
| --- | --- | --- | --- |
| `CONFIG` | `release`, `debug` | `release` | Optimized or diagnostic build |
| `GF_STORAGE` | `sp`, `dp` | `sp` | Green's-function storage precision |
| `USE_MUMPS` | `0`, `1` | `0` | Compile optional MUMPS support |
| `FC` | compiler path/name | `mpiifx` | MPI Fortran compiler wrapper |

MUMPS has not been validated for the workstation workflow. Its optional HPC
configuration is documented in [the HPC guide](HPC.md).

## 4. Run the Small TTI example

Compilation and execution are separate. The launcher runs the executable
currently in `bin/Seis2D_GQG_FWI` and does not rebuild it.

From the repository root:

```bash
./scripts/run_linux.sh
```

If execution permission was lost when downloading or copying the repository:

```bash
chmod +x scripts/run_linux.sh
```

The launcher:

- locates the repository without a user-specific path;
- uses the files in `examples/Small_TTI`;
- fixes MPI to one rank;
- detects the logical CPUs visible to Linux;
- sets the OpenMP, oneMKL and program thread limits consistently;
- detects total and currently available Linux/WSL memory;
- warns, but does not refuse to run, when available memory is low;
- creates a timestamped result directory;
- records normal output, errors and the execution environment separately.

The numerical solver is selected in `MainInput.inp`, not in the launcher.

## 5. Inputs and overrides

The executable receives five positional file roles:

1. main job configuration;
2. working elastic model;
3. working attenuation/Q model;
4. true/reference elastic model;
5. true/reference attenuation/Q model.

The working elastic model is mandatory. For pure forward modelling, the model
to simulate uses this same slot. Optional roles remain positional, so the
launcher passes explicit empty arguments when a middle role is unused.

The bundled examples use:

```text
MainInput.inp
InitInput.inp
InitInputQ.inp
TrueInput.inp
TrueInputQ.inp
```

To use different filenames, edit the visible input block near the top of
`scripts/run_linux.sh`:

```text
INPUT_DIR="${PROJ_ROOT}/examples/Small_TTI"
MainInput="${INPUT_DIR}/MainInput.inp"
InitInput="${INPUT_DIR}/InitInput.inp"
```

For a synthetic elastic inversion without Q files, keep all five positional
roles but assign an empty string to the unused Q roles:

```text
InitInputQ=""
TrueInput="${INPUT_DIR}/TrueInput.inp"
TrueInputQ=""
```

See [Input files](INPUT_FORMATS.md) for the complete formats and scientific
roles.

## 6. CPU and memory settings

By default the launcher uses all logical CPUs visible to Linux or WSL. A lower
thread count can be requested for testing:

```bash
OMP_NUM_THREADS=8 ./scripts/run_linux.sh
```

Do not request more threads than the CPUs available to WSL; oversubscription
can substantially reduce performance. MPI remains fixed to one workstation
rank.

The launcher exports currently available memory to the program as `MEM_MAX`
and records it in `environment.txt`. 

Monitor memory from another terminal with:

```bash
watch -n 2 'free -h'
```

## 7. Results and logs

Each run creates a directory such as:

```text
runs/Small_TTI/output_YYYYMMDD_HHMMSS/
```

| File | Contents |
| --- | --- |
| `run.log` | Normal program and timing output |
| `error.log` | MPI/Fortran errors and available kernel OOM evidence |
| `environment.txt` | Compiler, MPI, CPU, memory and thread settings |
| `configuration_summary.txt` | Selected MainInput filename, SHA-256 checksum and complete contents |
| `out_diag.txt` | Per-iteration inversion diagnostics |
| `out_resid.txt` | Data and residual information |
| `out_linesearch.txt` | Line-search history |
| `out_lbfgs_hist.txt` | L-BFGS monitoring history |

Model, gradient, Green's-function and calculated-data products are written in
the same directory. An empty or nearly empty `error.log` is normal. After an
abnormal exit, inspect both `error.log` and the end of `run.log`.

## 8. Protect long runs with tmux

Run long inversions from the standalone Ubuntu terminal inside `tmux` so that
a VS Code or terminal disconnection does not stop the process:

```bash
tmux new -s seis2d
./scripts/run_linux.sh
```

Detach by pressing `Ctrl+B`, releasing both keys, and then pressing `D`.
Reconnect later with:

```bash
tmux attach -t seis2d
```

## 9. Common build problem

If `mpiifx` or `mpirun` is not found, initialize oneAPI in the current terminal:

```bash
source /opt/intel/oneapi/setvars.sh
```

Building under `/mnt/c` is supported and convenient. Keeping the repository in
WSL's native Linux filesystem may improve compilation and I/O performance, but
is not required for Small TTI.
