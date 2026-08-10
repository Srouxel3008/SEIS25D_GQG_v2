# Output files

This guide describes the files produced by SEIS25D_GQG_v2 and by the supplied
workstation and Slurm launchers. The exact set depends on the selected run mode,
active parameters, frequencies, number of iterations, preprocessing options
and diagnostic-output settings. A forward-only run therefore produces fewer
files than a full inversion.

See [Building and running on a workstation](WORKSTATION.md) or [Building and
running on an HPC cluster](HPC.md) for result-directory creation and launch
instructions.

## Plotting scripts

The publication package provides two equivalent plotting scripts under
`scripts/` for displaying the principal model, gradient and spectral-data
outputs described below:

- `scripts/seis25d_plot.py` for Python;
- `scripts/seis25d_plot.m` for MATLAB.

The Python version uses freely available libraries and is the recommended
option for users without MATLAB. Both scripts keep the case-dependent plotting
settings near the top of the file so that users can readily adjust paths,
parameter names, well position and plotting-depth limits. These scripts are
intended for common QC; specialist diagnostic files and publication quality
may require additional plotting code.

## Result directory

The launchers change into a new result directory before starting the program,
so the program's relative output filenames are collected in one place.

Workstation example:

```text
runs/Small_TTI/output_YYYYMMDD_HHMMSS/
```

Slurm example:

```text
runs/Marmousi_VSP/slurm_JOB_ID/
```

The `runs/` tree is ignored by Git because a scientific run can produce many
gigabytes. Copy or archive important results separately.

## Launcher records

These files are created by the launchers rather than by the Fortran program.

| File | Purpose |
| --- | --- |
| `run.log` | Standard program output, progress messages and timings |
| `error.log` | Standard error and available out-of-memory evidence |
| `environment.txt` | Executable, inputs, host, compiler/MPI information, CPU, memory and thread settings |
| `configuration_summary.txt` | Selected MainInput filename, SHA-256 checksum and complete file contents |
| `slurm-JOB_ID.out` | Slurm scheduler standard output, copied into the HPC result directory |
| `slurm-JOB_ID.err` | Slurm scheduler standard error, copied into the HPC result directory |

An empty or nearly empty `error.log` is normal. For a failed run, inspect the
end of both `run.log` and `error.log`. On HPC, also inspect the original Slurm
logs in the submission directory if copying was interrupted.

## Geometry and model-grid exports

| Pattern | Meaning |
| --- | --- |
| `m_Grid.dat` | Numbered computational-grid point coordinates |
| `Seis2D_GQG.txt` | Legacy coordinate export, including source/receiver interpolation information |
| `m_Top.dat` | Topography and related plotting coordinates |
| `m_SR.dat` | Source/receiver geometry prepared for plotting and QC |
| `data_mapping.txt` | Mapping between internal spectral-data indices and acquisition entries, when written |
| `out.Gmask.txt` | Gradient mask before the optional taper |
| `out_GmaskTap.txt` | Tapered gradient mask, when tapering is enabled |

`m_Grid.dat` and `Seis2D_GQG.txt` contain a point index followed by the two
model coordinates. Geometry coordinates use the same coordinate convention and
units as the corresponding input files.

## Gridded field format

Model, gradient, updated-model, Green-function and many Fréchet files use a
plain-text rectangular grid designed for direct plotting:

```text
NZ  z(1)  z(2) ... z(NZ)
x(1)  value(1,1) ... value(1,NZ)
x(2)  value(2,1) ... value(2,NZ)
...
```

The first number on the first row is the number of vertical samples. The
remaining values on that row are the vertical coordinates. Each later row
starts with one horizontal coordinate followed by the field values at all
vertical coordinates. Some interpolation routines represent air cells with a
large sentinel value (currently approximately `1.0E20`); plotting scripts
should mask these values rather than display them as physical properties.

The values retain the physical or derived units of the named parameter. Model
exports are not automatically converted to plotting units such as GPa or
g/cm³.

## Initial and reference models

| Pattern | Meaning |
| --- | --- |
| `mI_{PARAM}.dat` | Initial or working model for `PARAM` |
| `mT_{PARAM}.dat` | True/reference model for a synthetic case |

`PARAM` is the internal parameter label, for example `C33`, `rho`, `the` or a
Q parameter. Only parameters available in the selected physical model are
written. A real-data run without a true/reference model will not produce the
corresponding `mT_*` files.

## Spectral data and calculated responses

### Trace files

`GT0_IFQ_ITER.txt` contains the observed or synthetic-reference complex
spectrum. `GO_IFQ_ITER.txt` contains the response calculated from the current
model. `IFQ` is the one-based frequency index and `ITER` is the iteration
number; reference data are commonly written with iteration zero.

Each row has 14 columns:

| Column | Meaning |
| ---: | --- |
| 1 | Frequency index `IFQ` |
| 2 | Iteration |
| 3 | Data index |
| 4 | Source index |
| 5 | Source component index |
| 6 | Receiver index |
| 7 | Receiver component index |
| 8–9 | Source X and Z coordinates |
| 10–11 | Receiver X and Z coordinates |
| 12–13 | Real and imaginary spectral values |
| 14 | Complex amplitude/magnitude |

Calculated `GO_*` files are QC-dependent and may not be retained for every
iteration in a reduced-output run.

### Forward-model export

For forward-only runs, `OBS_FREQ_FREQUENCY.txt` contains:

```text
frequency  data_index  real_part  imaginary_part
```

Despite the historical `OBS` prefix, this is the spectral dataset exported by
the forward-only path and can be used as external spectral input for a later
run after confirming the expected input naming and ordering.

## Residual and inversion histories

### `out_resid.txt`

When acquisition geometry is available, each row contains:

```text
iteration  data_index  source_index  source_component
receiver_index  receiver_component  source_x  source_z
receiver_x  receiver_z  residual_real  residual_imag  acquisition_weight
```

Rows are appended during the run. Their order follows the internal data
mapping, not necessarily a simple receiver-only sequence.

### `out_diag.txt`

One row is written for each active parameter and iteration:

```text
frequency  iteration  previous_cost  scaled_gradient_norm  residual_RMS
iteration_seconds  frequency_band  frequency_index  parameter_order
parameter_index  active_parameter_index  global_iteration
```

Use the frequency, iteration and parameter columns together when comparing
rows from a multi-frequency or multi-parameter inversion.

### `out_linesearch.txt`

The header written by the program defines:

```text
FREQ  ITER  trial  TAG  ALPHA  trial_cost  reference_cost
```

There may be several rows for one inversion iteration when a step is reduced
or retried.

### `out_lbfgs_hist.txt`

This is the detailed L-BFGS diagnostic history. Its header names the gradient,
direction, step and curvature quantities, accepted/skipped history columns,
and diagnostic reason codes. Treat the header in the generated file as the
authoritative column definition because the optimizer diagnostics may evolve.

## Gradient files

Gradient files are generated for each inverted parameter, frequency and
iteration:

| Pattern | Meaning |
| --- | --- |
| `GRADr_{PARAM}_FREQUENCY_ITNN.dat` | Masked physical/data gradient in parameter-scaled model space |
| `GRADrsc_{PARAM}_FREQUENCY_ITNN.dat` | Scaled and parameter-balanced gradient |
| `GRADsm_{PARAM}_FREQUENCY_ITNN.dat` | Smoothed, preconditioned gradient, when smoothing is enabled |
| `GRADrpm_{PARAM}_FREQUENCY_ITNN.dat` | Preconditioned gradient, when preconditioning is enabled |

`GRADrpm` is intentionally absent when preconditioning is disabled. The Small
TTI configuration used during workstation validation therefore produced
`GRADr` and `GRADrsc` but no `GRADrpm`; the Marmousi VSP configuration enables
and exports it.

## Updated model files

The current historical naming pattern is:

```text
PCBB_{{PARAM}_FREQUENCY_ITNN}.dat
PCBA_{{PARAM}_FREQUENCY_ITNN}.dat
```

`PCBB` is written from the current accepted model state. `PCBA` is an
additional updated/trial-model diagnostic written only when the corresponding
debug output is enabled. The doubled braces are part of the current filename,
not Markdown notation. Users should match these files by pattern rather than
constructing names manually.

## Fréchet-derivative files

When detailed Fréchet export is enabled, filenames include the data/source
index, parameter and frequency. Current prefixes include:

```text
F_RE...
F_IM...
F_RE_SRC1_SUM...
F_IM_SRC1_SUM...
F_ABS_SRC1_SUM...
F_RMS_SRC1_SUM...
```

These are diagnostic fields for the real part, imaginary part, and selected
source summaries. They can be numerous and very large. Their absence does not
mean that Fréchet derivatives were not calculated internally; file export is a
separate diagnostic choice.

## Green-function files

Optional wavefield exports use real/imaginary pairs with names containing the
field prefix, source/component indices and frequency, for example patterns of
the form:

```text
PREFIX_ReG{SC}_FREQUENCY.dat
PREFIX_ImG{SC}_FREQUENCY.dat
```

Green functions are normally held in memory for modelling and derivative
calculation. Disk export is intended for QC and can consume substantial space.

## Which files should be preserved?

For a reproducible scientific result, retain at minimum:

- every input file used by the run;
- `environment.txt`, `run.log` and `error.log`;
- `out_diag.txt`, `out_resid.txt`, `out_linesearch.txt` and
  `out_lbfgs_hist.txt` when applicable;
- the final accepted `PCBB_*` model for every inverted parameter;
- the initial and true/reference model exports used in figures;
- the exact plotting script and any manually selected plotting limits.

Intermediate wavefields, Fréchet fields, gradients and trial models may be
discarded after validation if storage is limited, provided they can be
regenerated from the archived inputs and software version.
