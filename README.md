# SEIS25D_GQG_v2

SEIS25D_GQG _v2 is a Fortran90 research code for 2.5-D frequency-domain seismic
forward modeling, and full-waveform inversion (FWI)
in elastic to viscoelastic TTI media.

The repository currently provides a small three-component tilted transverse
isotropy (TTI) example that can be compiled and run on a Linux workstation or
on Windows through Ubuntu in Windows Subsystem for Linux (WSL).

> **Project status:** the code and documentation are being prepared for their
> first public release.

## Capabilities

- 2D/2.5D frequency-domain forward modeling forelastic, viscoelastic and viscoelastic TTI media;
- 2D/2.5D frequency-domain forward FWI for elastic, viscoelastic and viscoelastic TTI media; ;

## Supported workstation platform

- native Linux, or Ubuntu under WSL on Windows;

Native Windows compilation is not currently supported. MUMPS is optional and
is disabled in the standard workstation build.

The bundled Small TTI example requires a minimum of 10 GB available RAM but has been run
successfully up to 22 Hz on a Lenovo Legion 5 laptop. The Marmousi_VSP example requires over 170 GB of available RAM. It can run on a
sufficiently equipped high-memory workstation, although a multi-node cluster
is the preferred option.

See [Building and running on a workstation](docs/WORKSTATION.md) or [Building and running on HPC](docs/HPC.md) for the full
installation, compilation, execution and troubleshooting workflow.

## Get the complete code

The repository must be downloaded as one complete folder because the source,
launchers, examples and documentation use the directory structure shown
below. Do not download individual Fortran or input files separately.

To download without Git:

1. Open the project page on GitHub.
2. Select **Code**, then **Download ZIP**.
3. Extract the ZIP file to a location with sufficient free disk space.
4. Open the extracted top-level folder—the repository root—which contains
   `README.md` together with `src`, `scripts`, `examples` and `docs`.

Windows users should keep the complete extracted folder together and access it
from Ubuntu/WSL. For example:

```text
C:\Users\path\SEIS25D_GQG_v2
```

is available inside WSL as:

```text
/mnt/c/Users/path/SEIS25D_GQG_v2
```

Users familiar with Git may instead clone the complete repository and run the
same commands from the cloned repository root.

## Quick start

Open a Linux or Ubuntu/WSL terminal and initialize Intel oneAPI:

```bash
source /opt/intel/oneapi/setvars.sh
```

Move to the downloaded repository, then compile a release build without
MUMPS:

```bash
make -C src CONFIG=release USE_MUMPS=0
```

Run the bundled Small TTI example:

```bash
./scripts/run_linux.sh
```

If the downloaded script is not executable, enable it once:

```bash
chmod +x scripts/run_linux.sh
```

Compilation and execution are deliberately separate. Recompile after changing
the Fortran source; changing only an input file does not require recompilation.

## What the launcher does

The workstation launcher:

- locates the repository automatically, so no fixed project path is required;
- runs one MPI rank;
- detects the available logical CPUs and selects an OpenMP thread count;
- records total and available WSL/Linux memory and warns when memory is low;
- creates a timestamped directory under `runs/Small_TTI`;
- writes `run.log`, `error.log`, and `environment.txt`;

Thread counts and input/output locations can be overridden when needed. See
the [workstation guide](docs/WORKSTATION.md) for settings and monitoring.

## Input files

The executable currently accepts five positional input roles:

1. main job configuration;
2. working elastic model (always required);
3. working attenuation model, when attenuation is enabled;
4. true/reference elastic model for synthetic inversion;
5. true/reference attenuation model for synthetic viscoelastic inversion.

For forward modeling, the model to simulate is supplied in the working elastic
model position. Pure elastic runs do not require the Q files. Because the
current interface is positional, unused middle positions must remain explicit
empty arguments when a later file is required; the supplied launcher handles
the bundled combination.

See [Input files](docs/INPUT_FORMATS.md) for file layouts, units, indexing,
geometry, acquisition patterns, and external spectral data.

## Repository layout

```text
src/                 Fortran source and Makefile
scripts/             Workstation launcher and cluster examples
examples/Small_TTI/  Small runnable TTI input set
docs/                Build, running, and input documentation
bin/                 Compiled executable (generated)
build/               Objects and module files (generated)
runs/                Run outputs and logs (generated, not committed)
tests/                Local validation outputs (not committed)
```

Large generated results are intentionally excluded from Git. Preserve the
input set, `environment.txt`, and relevant logs with any result that must be
reproduced.

## Documentation

- [Standalone Word user manual](SEIS25D_GQG_v2_Manual.docx)
- [Build and run on a workstation](docs/WORKSTATION.md)
- [Build and run on Slurm](docs/HPC.md)
- [Prepare input files](docs/INPUT_FORMATS.md)
- [Understand output files](docs/OUTPUTS.md)
- [Source organization](docs/SOURCE_CLASSIFICATION.md)
- [Authoritative source-copy verification](docs/SOURCE_COPY_VERIFICATION.md)

## Previous presentation

Preliminary results produced with this code were presented at the AGU25 Annual Meeting.
This repository is being prepared to accompany the first full publication of
the software and methodology.

Rouxel, S., and Bouchaala, F. (2025). *Advancing 2.5D Frequency-Domain Full
Waveform Inversion in Viscoelastic Anisotropic Media for VSP Applications*.
AGU25 Annual Meeting, presentation S33B-0258, abstract 1940271.
[AGU abstract](https://agu.confex.com/agu/agu25/meetingapp.cgi/Paper/1940271)

## Contact

S. Rouxel PhD in Earth Sciences, Khalifa University of Science and Technology  
Institutional email: [100061882@ku.ac.ae](mailto:100061882@ku.ac.ae)  
Alternative email: [sedrouxelian@hotmail.com](mailto:sedrouxelian@hotmail.com)

## License and citation

This project is distributed under the [GNU General Public License version 3](LICENSE).
The preferred scientific citation will be added before the first public
release.
