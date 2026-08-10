! file: vea_constants.f90
MODULE constant_mod
  USE iso_fortran_env, ONLY: dp => real64
  IMPLICIT NONE
  PRIVATE
  ! Core circle constants
  REAL(dp), PARAMETER, PUBLIC :: pi      = ACOS(-1.0_dp)
  ! Angle conversion
  REAL(dp), PARAMETER, PUBLIC :: deg2rad = pi/180.0_dp
  REAL(dp), PARAMETER, PUBLIC :: rad2deg = 180.0_dp/pi
  REAL(dp), PARAMETER, PUBLIC :: tiny_dp = 1.0e-50_dp
  REAL(dp), PARAMETER, PUBLIC :: BYTES_GB     = 1024.0_dp**3
REAL(dp), PARAMETER,PUBLIC :: BYTES_CPLX16 = 16.0_dp
END MODULE constant_mod
