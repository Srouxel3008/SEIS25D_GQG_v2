# Input files

This document describes the input files accepted by SEIS25D_GQG_v2 and how they
fit together. It is based on the current Fortran reader and the bundled
`examples/Small_TTI` example.

> **Current transition:** the grid is still embedded in `MainInput.inp`.
> Topography and acquisition geometry may be embedded or read from external files, but
> their selection remains inside `MainInput.inp` rather than using named
> command-line options.

## Files passed to the executable

The executable currently uses five positional file slots:

| Position | Role | When required |
|---:|---|---|
| 1 | `MainInput` job configuration | Always |
| 2 | Working elastic model (`InitInput`) | Always |
| 3 | Working attenuation model (`InitInputQ`) | Viscoelastic runs only |
| 4 | True/reference elastic model (`TrueInput`) | Synthetic inversion |
| 5 | True/reference attenuation model (`TrueInputQ`) | Synthetic viscoelastic inversion |

For forward modeling, position 2 contains the model to simulate even though
its historical name is `InitInput`. For pure elastic work, positions 3 and 5
are unused.

The interface is positional. If an unused slot occurs before a used one, pass
an empty quoted argument (`""`) so that later files do not shift into the
wrong roles. An `InitInputQ` without an elastic `InitInput` is not a valid
combination.

The supplied Linux launcher contains the recommended argument combinations.

## `MainInput.inp`

`MainInput.inp` controls the job, frequencies, physics, solver, optimization,
grid, geometry, and acquisition pattern. Comment-like heading lines are part
of the current ordered format: retain the structure and change the values
beneath the headings.

### Job and mode

```text
--JOB AND MODE CONTROL INV -- I25D--
1 1
```

The first integer is `INV`:

| Value | Operation |
|---:|---|
| 0 | Synthetic forward modeling |
| 1 | Synthetic full-waveform inversion |
| 2 | Inversion of externally supplied observations |


`I25D` selects the 2.5-D formulation (`1`) or 2D mode (`0`).


### Frequency bands and observed data

```text
--FREQUENCY BANDS--
2
3 7.0 12.0 22.0
2 37.0 50.0
--FREQUENCY DATA FILES--
-
```

The first value is the number of frequency bands. Each following band line
starts with the number of frequencies in that band, followed by those
frequencies in hertz.

A single `-` indicates that no external spectra are supplied and requests
internally generated synthetic data. 

## External spectral data
For external observations, supply one spectral-data filename per frequency.
The filename order must match the frequency order.

External-data runs require one text file per frequency. A typical file is:

```text
Freq ND Real Imag
10.0 1 10.94 13.05
10.0 2 11.20 12.74
```

Each row contains frequency, acquisition-row index, and the real and imaginary
parts of the complex spectrum. The row order and count must match the selected
acquisition pattern. A short text header is accepted. 

### Physics and inverted parameters

```text
--ANISOTROPIC MODEL (IANISO ITHOM IVISCO CMIN CMAX)--
7 0 1 1.0 2.5
'rho', 'c11', 'c13', 'c33', 'c44', 'c66', 'the', 'Q11', ...
0 0 0 1 0 0 0 0 0 0 0 0
```

The current reader consumes `IANISO`, `ITHOM`, `IVISCO`, `CMIN`, and `CMAX`
from the value line. 
`IANISO` selects TTI (`7`), VTI (`6`) or isotropic (`3`)
elastic parameterization. Small TTI uses `rho, C11, C13, C33, C44, C66,
theta`; Marmousi is VTI and uses `rho, C11, C13, C33, C44, C66`.

`ITHOM` enables (`1`) or disables (`0`) Thomsen parameterization using `rho`,
`Vp`, `Vs`, `epsilon`, `delta`, and `gamma`. Viscoelastic inversion is not yet
available in that parameterization. 
`IVISCO` enables (`1`) or disables (`0`)
attenuation.
`CMIN` and `CMAX` are the expected minimum Vs and maximum Vp used
by the preliminary grid checks.

The next line names the active elastic and Q parameters in their file order.
The following integers assign inversion groups: `0` keeps a parameter fixed;
a positive integer assigns it to an inversion group/stage. The number of
entries must match the selected parameterization. Use one assignment line for
the entire inversion, or one tailored line per active frequency band. For two
frequency bands:
```text
'rho', 'c11', 'c13', 'c33', 'c44', 'c66', 'the', 'Q11', ...
0 0 0 1 0 0 0 0 0 0 0 0
0 0 0 1 1 0 0 0 0 2 0 0
```

### Solver and inversion controls

```text
--SOLVER CHOICE--
2
--REGULARIZATION (REG lambda)--
0 0.0
--LBFGS SETTINGS (MML LBFGS_TYPE)--
5 0
--OPTIMIZATION CONTROLS (MAXITER cost_conv_tol misfit_conv_tol)--
10 1e-3 1e-4
--PRECONDITIONING (active, beta value)--
0 0.001
--GRADIENT TAPER AND SMOOTHING (active, Z1, Z2, sigmax, sigmaz) --
0 0.0 0.0 0 0
```

Solver choices are:

| Value | Solver |
|---:|---|
| 1 | LU legacy solver |
| 2 | Intel oneMKL PARDISO |
| 3 | MUMPS |

A build made without MUMPS will reject solver `3` with an explicit message.
The remaining lines control regularization, L-BFGS history/type, maximum
iterations and convergence tolerances, preconditioning strength, and gradient
tapering/smoothing. Begin with the supplied example values. Regularization is
not fully implemented and should remain set to `0 0`.



### Grid sampling

```text
--GLOBAL GRID (NORD DX DZ)--
3 5.0 5.0
```

- `NORD` is the GQG interpolation order, typically `3` or `5`;
- `DX` is the horizontal GQG subdomain size in metres;
- `DZ` is the vertical GQG subdomain size in metres.

### Frequency-dependent grids (optional)

For a multiband run, each frequency band may use its own grid order and
sampling. Place the optional block immediately after the global grid line:

```text
--GLOBAL GRID (NORD DX DZ)--
NORD(NFBAND) DX(NFBAND) DZ(NFBAND)
--Band meshing
NORD(1) DX(1) DZ(1)
NORD(2) DX(2) DZ(2)
...
NORD(NFBAND) DX(NFBAND) DZ(NFBAND)
```

The heading must contain the case-sensitive text `Band meshing`. Exactly one
value row is required for each frequency band, in the same order as the bands
listed under `--FREQUENCY BANDS--`. Each row contains:

- `NORD`: GQG interpolation order for that band, from 2 to 10;
- `DX`: horizontal GQG subdomain size in metres, greater than zero;
- `DZ`: vertical GQG subdomain size in metres, greater than zero.

For example, a two-band inversion can use:

```text
--GLOBAL GRID (NORD DX DZ)--
5 50.0 50.0
--Band meshing
5 60.0 60.0
5 50.0 50.0
```

The first row applies to frequency band 1 and the second to band 2. This makes
it possible to use a coarser, less expensive grid for lower frequencies and a
finer grid for higher frequencies.

The global `NORD DX DZ` line remains mandatory even when this optional block
is present and should describe the finest desired grid. Every selected grid is
checked against `CMIN`, `CMAX`, and its band's frequency range before modeling
begins.

## Grid block

The grid currently follows the global-grid settings inside `MainInput.inp`:

```text
MX
MZ(1) MZ(2) ... MZ(MX)
--GRID COORDINATES--
x(1) z(1,1) ... z(1,MZ(1))
x(2) z(2,1) ... z(2,MZ(2))
...
```

`MX` is the number of horizontal grid columns. `MZ(i)` is the number of
vertical points in column `i`. Each coordinate row begins with its horizontal
coordinate and contains that column's vertical coordinates. Coordinates are
in metres; the example uses negative depth below `z = 0`. Rows may have
different lengths when the topography is irregular.

Every elastic and Q model must use the same `MX`, `MZ`, and row ordering.
`grid_file.txt` in the example is currently a reference copy, not a separate
command-line input.

## Topography

`MainInput.inp` first supplies the topography size and surface flag:

```text
NTO IS0
```

This is followed either by `NTO` coordinate rows embedded directly in the
main input, or by the path to an external topography file. The external file
repeats the header:

```text
NTO IS0
x(1) z(1)
...
x(NTO) z(NTO)
```

`NTO` is the number of surface points, followed by `NTO` x/z coordinate pairs
in metres. `IS0` selects a free surface (`0`) or top absorbing boundary (`1`).
Small TTI embeds its flat surface; Marmousi supplies `topo.txt` externally.

## Elastic model files

`InitInput` and `TrueInput` contain one block per elastic parameter, in the
exact order selected by `IANISO` and `ITHOM`. Each block has a label and `MX`
data rows:

```text
--PARAMETER NAME--
x(1) value(1,1) ... value(1,MZ(1))
...
x(MX) value(MX,1) ... value(MX,MZ(MX))
```

For the bundled TTI case, the blocks are `rho`, `C11`, `C13`, `C33`, `C44`,
`C66`, and `theta`. Density is in kg/m³, stiffness values are supplied in GPa,
and angle is in degrees. `InitInput` is the working/starting model;
`TrueInput` is the reference model used to generate synthetic observations.

For isotropic and VTI parameterizations, elastic properties may instead be
provided as `rho`, `Vp`, and `Vs`, with velocities in km/s. Thomsen VTI also
uses `epsilon`, `delta`, and `gamma` and their associated Q values; viscoelastic
FWI is not currently available in that parameterization.

## Attenuation model files

`InitInputQ` and `TrueInputQ` use the same grid-shaped block format. The Small
TTI case contains `Q11`, `Q13`, `Q33`, `Q44`, and `Q66`, which are
dimensionless. These files are used only when `IVISCO = 1`.

## Source and receiver geometry

The main input defines the total number of source/receiver points (`NSR`), the
number of sources (`NSS`), and the retained legacy setting `ISR90`. It accepts
either a compact regular-geometry generator or individually listed points.

### Regular source/receiver generator

The generator is selected automatically when the first non-blank line after
`NSR NSS ISR90` is a short numeric component-setting line. Its block is:

```text
NSR NSS ISR90
NCOMPS NCOMPR
source_x0 source_dx source_z0 source_dz NSS
top_x0 top_dx top_z0 top_dz N_top
left_x0 left_dx left_z0 left_dz N_left
right_x0 right_dx right_z0 right_dz N_right
bottom_x0 bottom_dx bottom_z0 bottom_dz N_bottom
```

`NCOMPS` and `NCOMPR` are the numbers of source and receiver components
(1 to 3). Only these first two values are used on that line; a trailing legacy
value in older inputs is ignored. For point `j`, the generator uses
`x = x0 + (j - 1) dx` and `z = z0 + (j - 1) dz`; all generated `y` coordinates
are zero. Set a receiver-side count to `0` to omit that side.

Points are ordered as sources, then top, left, right and bottom receivers. The
program calculates `NRR = N_top + N_left + N_right + N_bottom` and recomputes
`NSR = NSS + NRR`. Keep the written `NSR` consistent with that total, and keep
all generated coordinates inside the model bounds. The bundled Small TTI input
is the working example.

### Individually listed geometry

The first `NSS` rows are sources and the remaining `NSR - NSS` rows are
receivers. Rows may be embedded in `MainInput.inp` or supplied through an
external `geom.txt` path. Each row is:

```text
x y z ICSR dir1_x dir1_y dir1_z [dir2_x dir2_y dir2_z ...]
```

Coordinates are in metres. `ICSR` is the number of components (1 to 3),
followed by that many direction vectors. Component conventions are X, Y, Z;
a one-component case is Z, while two components are X and Z.

The current reader expects exactly `NSR` data rows and does **not** consume a
count header from `geom.txt`.

## Acquisition pattern

After the grid, the data-acquisition section contains an optional acquisition
weight filename  (`-` for none), followed by the pattern selection.

A value of `1` generates every valid source/receiver/component pairing.
Alternatively, the list can be embedded or supplied in an external
`acq_pattern.txt` file:

```text
ND
row_id source_id source_component receiver_id receiver_component
...
```

The row ID may be omitted, leaving four integers per row. Component IDs are
`1 = X`, `2 = Y`, and `3 = Z`. Source IDs run from `1` to `NSS`. Receiver IDs
use the global geometry numbering, from `NSS + 1` through `NSR`. `ND` is the
number of acquisition rows.

Acquisition weighting is supported, but the legacy weight-file layout still
needs a dedicated reader audit before it is presented as a stable public
format.

The input format is intentionally documented as it works today.
