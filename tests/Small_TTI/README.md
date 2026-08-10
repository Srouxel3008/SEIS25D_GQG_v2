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
- `logs/hpc_PARAM.log`: sanitized HPC standard output for each of the 12
  parameter inversions. Personal account identifiers and home-directory paths
  have been replaced by neutral placeholders.

