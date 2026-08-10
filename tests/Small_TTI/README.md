# Small TTI test results

This directory contains the final results selected from the 12 independent HPC
inversions used for the Small TTI example. Each inversion recovers one model
parameter while the other parameters remain fixed.

## Contents

- `results/mI_{PARAM}.dat`: initial model for each parameter;
- `results/mT_{PARAM}.dat`: true model for each parameter;
- `results/PCBB_{{PARAM}_FREQUENCY_ITNN}.dat`: final accepted model selected
  from each inversion;
- `results/m_Top.dat`: topography prepared by the program;
- `results/m_SR.dat`: source and receiver geometry prepared by the program;
- `figures/Final_Inverted_Models_AllParams.png`: composite recovered-model
  figure for all 12 parameters;
- `figures/True_Models_AllParams.png`: corresponding true models;
- `figures/Final_Profiles_AllParams_x100.png`: final profile comparison.
- `logs/hpc_PARAM.log`: sanitized HPC standard output for each of the 12
  parameter inversions. Personal account identifiers and home-directory paths
  have been replaced by neutral placeholders.

The final accepted files are:

| Parameter | Final frequency | Final iteration |
| --- | ---: | ---: |
| C11 | 50 Hz | 2 |
| C13 | 50 Hz | 3 |
| C33 | 50 Hz | 1 |
| C44 | 50 Hz | 2 |
| C66 | 37 Hz | 5 |
| Q11 | 50 Hz | 4 |
| Q13 | 50 Hz | 3 |
| Q33 | 50 Hz | 4 |
| Q44 | 50 Hz | 1 |
| Q66 | 37 Hz | 5 |
| rho | 50 Hz | 10 |
| theta | 50 Hz | 4 |

The composite images were created with the historical plotting workflow used
for the study. The supplied `scripts/seis25d_plot.py` and
`scripts/seis25d_plot.m` scripts can be used to inspect and display the
individual numerical model files, but they do not recreate the exact
12-panel publication layout.

These are selected final scientific outputs, not complete run directories.
Intermediate models, residual histories, temporary files and duplicate
plotting products were intentionally omitted from the repository. The original
scheduler filenames were normalized to parameter-based names; job timings and
scientific program output were retained.
