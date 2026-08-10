# Marmousi VSP test results

This directory contains the final results selected from the 11 independent HPC
inversions displayed in the Marmousi VTI composite figure. There is no
tilt-angle inversion in this VTI test.

Unlike the Small TTI test, the Marmousi inversions use parameter-specific beta
values. The exact numbered `MainInput` files used for these runs are supplied
under `examples/Marmousi_VSP`.

## Contents

- `results/mI_{PARAM}.dat`: initial model for each parameter;
- `results/mT_{PARAM}.dat`: true model for each parameter;
- `results/PCBB_{{PARAM}_FREQUENCY_ITNN}.dat`: final accepted model selected
  from each inversion;
- `results/m_Top.dat` and `results/m_SR.dat`: plotting geometry prepared by
  the program;
- `diagnostics/out_diag_PARAM.txt`: inversion diagnostic history;
- `logs/hpc_PARAM.log`: sanitized HPC standard output. Personal account
  identifiers and home-directory paths were replaced by neutral placeholders;
- `figures/`: the recovered-model, true-model, profile and normalized-cost
  composite figures.

## Parameter configurations and final models

| Input ID | Parameter | Beta | Final frequency | Final iteration |
| ---: | --- | ---: | ---: | ---: |
| 1 | rho | 0.005 | 8.10 Hz | 5 |
| 2 | C11 | 0.005 | 9.00 Hz | 4 |
| 3 | C13 | 0.005 | 9.00 Hz | 7 |
| 4 | C33 | 0.001 | 9.00 Hz | 1 |
| 5 | C44 | 0.001 | 9.00 Hz | 7 |
| 6 | C66 | 0.001 | 9.00 Hz | 1 |
| 7 | Q11 | 0.007 | 9.00 Hz | 6 |
| 8 | Q13 | 0.005 | 9.00 Hz | 7 |
| 9 | Q33 | 0.001 | 7.40 Hz | 1 |
| 10 | Q44 | 0.001 | 6.80 Hz | 7 |
| 11 | Q66 | 0.005 | 8.85 Hz | 3 |

For example, input ID 4 corresponds to
`examples/Marmousi_VSP/MainInput_C33.inp`. The other configurations use the
same `MainInput_PARAMETER.inp` naming pattern.

## Deliberate exclusions

These are selected final scientific test outputs, not complete raw run
directories. Intermediate model states, residual histories, solver temporary
files and duplicate plotting products were omitted. Several raw Marmousi
residual files are larger than GitHub's normal 100 MiB per-file limit and are
not required to reproduce the final model composite.

The supplied Python and MATLAB plotting scripts can inspect the individual
model and diagnostic files, but they do not recreate the exact 11-panel
publication layout. The historical composite is included as a visual reference.
