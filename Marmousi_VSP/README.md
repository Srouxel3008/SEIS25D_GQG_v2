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
- `logs/hpc_PARAM.log`: sanitized HPC standard output. 
- `figures/`: the recovered-model, true-model, profile and normalized-cost
  composite figures.

