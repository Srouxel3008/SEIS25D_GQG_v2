module input_mod
   USE mpi
   USE err_mpi_mod
   USE constant_mod
   use iso_fortran_env, only: dp => real64
   implicit NONE

contains
!***************************************************
!include subroutine to load model parameters depending on selected mode
! if INV=0 synthetic forward modeling only, if NSR=NRX=1 also output Frechet derivatives
! only need to load one starting model, +source and receiver patterns
! if INV=1 synthetic modeling FWI
! require ture model and starting model +source and receiver patterns
! if INV=2 real case inversion
! requires only a starting model + source and receiver locations + real data
!****************************************************
   SUBROUTINE Read_Job_Setup( &
      INV, I25D, NFBAND, NFQ, MAXITER, REG, FREQN, FREQ1, FREQ2, &
      IANISO, ITHOM, IVISCO, CMIN, CMAX, NPAR, NPAR_INV, INVP_ORD, INVP_ORD_BAND, &
      NORD, DX, DZ, GZSMTH, GZ1, GZ2, PRECOND, SOLVER_KIND, MODEL, &
      FREQ_FILE, LAMBDA, BETA, EXTERNAL_DATA, &
      LBFGS_TYPE, MML, USE_GR_SMOOTH, SIGMA_X_GRAD, SIGMA_Z_GRAD, &
      NFQ_PER_BAND, NORD_BAND, DX_BAND, DZ_BAND, &
      cost_conv_tol, model_conv_tol, my_rank, DEBUG_OUTPUT)

      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32
      IMPLICIT NONE

      ! ==== Args ====
      INTEGER, INTENT(OUT) :: INV, I25D
      INTEGER, INTENT(OUT) :: NFBAND, NFQ, MAXITER
      INTEGER, INTENT(OUT) :: REG, SOLVER_KIND, PRECOND
      REAL(KIND=8), INTENT(OUT) :: FREQ1, FREQ2
      REAL(KIND=8), ALLOCATABLE, INTENT(OUT) :: FREQN(:)

      ! anisotropy / physics
      INTEGER, INTENT(OUT) :: IANISO, ITHOM, IVISCO
      REAL(KIND=8), INTENT(OUT) :: CMIN, CMAX
      INTEGER, INTENT(OUT) :: NPAR, NPAR_INV
      INTEGER, ALLOCATABLE, INTENT(OUT) :: INVP_ORD(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: INVP_ORD_BAND(:, :)
      CHARACTER(LEN=3), INTENT(OUT) :: MODEL(22)

      ! grid / gradient shaping
      INTEGER, INTENT(OUT) :: NORD
      REAL(KIND=8), INTENT(OUT) :: DX, DZ
      INTEGER, INTENT(OUT) :: GZSMTH          ! taper on/off (0/1)
      REAL(KIND=8), INTENT(OUT) :: GZ1, GZ2        ! taper depths [m]

      ! frequency data files
      CHARACTER(LEN=256), ALLOCATABLE, INTENT(OUT) :: FREQ_FILE(:)

      ! IO / diag
      INTEGER, INTENT(IN)  :: my_rank
      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT

      ! ---- new outputs ----
      REAL(KIND=8), INTENT(OUT) :: LAMBDA, BETA
      INTEGER, OPTIONAL, INTENT(OUT) :: LBFGS_TYPE
      INTEGER, INTENT(OUT) :: MML
      LOGICAL, INTENT(OUT) :: USE_GR_SMOOTH, EXTERNAL_DATA
      REAL(KIND=8), INTENT(OUT) :: SIGMA_X_GRAD, SIGMA_Z_GRAD
      INTEGER, ALLOCATABLE, INTENT(OUT) :: NFQ_PER_BAND(:)
      INTEGER, OPTIONAL, ALLOCATABLE, INTENT(OUT) :: NORD_BAND(:)
      REAL(KIND=8), OPTIONAL, ALLOCATABLE, INTENT(OUT) :: DX_BAND(:), DZ_BAND(:)
      REAL(KIND=8), INTENT(OUT) :: cost_conv_tol, model_conv_tol

      ! ==== locals ====
      INTEGER, PARAMETER :: unit_main = 1
      CHARACTER(LEN=160) :: LAB
      INTEGER :: ios, ios2
      LOGICAL :: dbg

      ! frequency bands
      INTEGER                        :: ib, j, total_nf, NFQ_B, k
      CHARACTER(LEN=512)             :: bandline, invp_band_line, line
      REAL(KIND=8), ALLOCATABLE      :: tmpfreq(:)
      INTEGER, ALLOCATABLE           :: invp_band_tmp(:, :)

      ! LBFGS controls (local, then propagate if PRESENT)
      INTEGER :: lbfgs_type_loc

      ! band sampling helpers
      INTEGER, ALLOCATABLE :: band_start(:), band_end(:)
      INTEGER :: sumq, i1b, i2b, NORD_use
      REAL(KIND=8) :: DX_use, DZ_use, FREQ_min_b, FREQ_max_b, WLMAX_b, WLMIN_b
      INTEGER :: nloc
      REAL(KIND=8) :: dxloc, dzloc

      ! misc
      REAL(KIND=8) :: sigma_x_grad_in, sigma_z_grad_in

50    FORMAT(160A)

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      !--------------------------------------------------------------------
      ! (1) Job definition INV  I25D
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) INV, I25D
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading INV, I25D')
      IF (INV < 0 .OR. INV > 3) CALL fail('Read_Job_Setup: INV must be 0,1,2,3  ')
      IF (I25D < 0 .OR. I25D > 1) CALL fail('Read_Job_Setup: I25D must be 0 or 1')
      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, *) 'Job type INV =', INV, ' I25D =', I25D
      END IF

      !--------------------------------------------------------------------
      ! (2) Frequency bands
      ! label
      !   NFBAND
      ! then NFBAND lines:
      !   NFQ_B  f1  f2  ...  fNFQ_B
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) NFBAND
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading NFBAND')
      IF (NFBAND < 1) CALL fail('Read_Job_Setup: NFBAND must be >= 1')
      IF (dbg) write (*, *) 'NFBAND = ', NFBAND

      ALLOCATE (NFQ_PER_BAND(NFBAND))
      ALLOCATE (FREQN(0))
      total_nf = 0

      DO ib = 1, NFBAND
         READ (unit_main, '(A)', IOSTAT=ios) bandline
         IF (ios /= 0) CALL fail('Read_Job_Setup: error reading frequency band line #'//itoa(ib))

         READ (bandline, *, IOSTAT=ios2) NFQ_B
         IF (ios2 /= 0 .OR. NFQ_B < 1) CALL fail('Read_Job_Setup: invalid NFQ_B on band #'//itoa(ib))
         NFQ_PER_BAND(ib) = NFQ_B

         ALLOCATE (tmpfreq(NFQ_B))
         READ (bandline, *, IOSTAT=ios2) NFQ_B, (tmpfreq(j), j=1, NFQ_B)
         IF (ios2 /= 0) CALL fail('Read_Job_Setup: error parsing frequencies on band #'//itoa(ib))

         CALL append_freqs(FREQN, tmpfreq, NFQ_B)
         DEALLOCATE (tmpfreq)
         total_nf = total_nf + NFQ_B
      END DO

      NFQ = total_nf
      IF (NFQ < 1) CALL fail('Read_Job_Setup: total NFQ < 1')
      FREQ1 = MINVAL(FREQN)
      FREQ2 = MAXVAL(FREQN)

      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, *) 'Bands =', NFBAND, '  Total NFQ =', NFQ, '  fmin =', FREQ1, '  fmax =', FREQ2
      END IF

      !--------------------------------------------------------------------
      ! (3) Per-frequency data files
      ! label
      ! then NFQ lines: one path per frequency (or 'N/A' for synthetic)
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      ! IF (INV == 2) THEN
      !    ALLOCATE (FREQ_FILE(NFQ))
      !    DO k = 1, NFQ
      !       READ (unit_main, '(A)', IOSTAT=ios) FREQ_FILE(k)
      !       IF (ios /= 0) CALL fail('Read_Job_Setup: error reading FREQ_FILE #'//itoa(k))
      !       FREQ_FILE(k) = ADJUSTL(TRIM(FREQ_FILE(k)))
      !       if (my_rank == 0 .and. dbg) WRITE (*, *) 'FREQ_FILE entry #', k, ' = ', FREQ_FILE(k)
      !    END DO
      ! ELSE
      !    ! For forward modeling (0) or integrated synthetic inversion (1): skip lines & keep consistent file structure
      !    ALLOCATE (FREQ_FILE(1))
      !    READ (unit_main, '(A)', IOSTAT=ios) FREQ_FILE(1)     ! consume line but ignore it
      ! END IFREAD(unit_main, 50) LAB

      READ (unit_main, '(A)', IOSTAT=ios) line
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading first FREQ_FILE line')

      line = ADJUSTL(TRIM(line))

      IF (TRIM(line) == '-' .OR. LEN_TRIM(line) == 0) THEN
         EXTERNAL_DATA = .FALSE.
         ALLOCATE (FREQ_FILE(1))
         FREQ_FILE(1) = ''
      ELSE
         EXTERNAL_DATA = .TRUE.
         ALLOCATE (FREQ_FILE(NFQ))
         FREQ_FILE(1) = line

         DO k = 2, NFQ
            READ (unit_main, '(A)', IOSTAT=ios) FREQ_FILE(k)
            IF (ios /= 0) CALL fail('Read_Job_Setup: error reading FREQ_FILE #'//itoa(k))
            FREQ_FILE(k) = ADJUSTL(TRIM(FREQ_FILE(k)))
         END DO
      END IF

      !--------------------------------------------------------------------
      ! (4) Anisotropic model header
      ! label
      !   IANISO  ITHOM   IVISCO  CMIN  CMAX
      ! then MODEL names, then INVP_ORD
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) IANISO, ITHOM, IVISCO, CMIN, CMAX
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading IANISO, ITHOM, IVISCO, CMIN, CMAX')

      IF (IVISCO < 0 .OR. IVISCO > 1) CALL fail('Read_Job_Setup: IVISCO must be 0/1')
      IF (IANISO < 0 .OR. IANISO > 22) CALL fail('Read_Job_Setup: IANISO out of range')
      IF (ITHOM < 0 .OR. ITHOM > 1) CALL fail('Read_Job_Setup: ITHOM must be 0/1')
      IF (CMAX < CMIN .OR. CMIN < 0.0_8) CALL fail('Read_Job_Setup: invalid CMIN/CMAX')

      ! convert km/s -> m/s if desired (keeps previous behaviour)
      CMIN = 1.0e3_8*CMIN
      CMAX = 1.0e3_8*CMAX

      ! NPAR logic as before
      IF (IVISCO == 0) THEN
         NPAR = IANISO
      ELSE
         IF (IANISO /= 7) THEN
            NPAR = 2*(IANISO - 1) + 1
         ELSE
            NPAR = 2*(IANISO - 2) + 2
         END IF
      END IF
      IF (NPAR > SIZE(MODEL)) CALL fail('Read_Job_Setup: NPAR exceeds MODEL array size')

      READ (unit_main, *, IOSTAT=ios) (MODEL(j), j=1, NPAR)
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading MODEL names')

      ALLOCATE (INVP_ORD(NPAR))
      READ (unit_main, *, IOSTAT=ios) (INVP_ORD(j), j=1, NPAR)
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading INVP_ORD')
      NPAR_INV = COUNT(INVP_ORD(1:NPAR) /= 0)

      !--------------------------------------------------------------------
      ! Optional per-band inversion-order overrides.
      !
      ! Supported layouts after the MODEL names:
      !   1) Legacy global-only:
      !        INVP_ORD(1:NPAR)
      !
      !   2) Preferred simplified band mode:
      !        band_1_order(1:NPAR)   ! also used as global fallback
      !        ...
      !        band_NFBAND_order(1:NPAR)
      !
      !   3) Older explicit-global + band overrides:
      !        INVP_ORD(1:NPAR)
      !        INVP_ORD_BAND(1,1:NPAR)
      !        ...
      !        INVP_ORD_BAND(NFBAND,1:NPAR)
      !
      ! The first row is always read into INVP_ORD first. If additional
      ! integer rows follow, detect whether they represent the simplified
      ! band table or the older explicit-global layout.
      !--------------------------------------------------------------------
      READ (unit_main, '(A)', IOSTAT=ios2) invp_band_line
      IF (ios2 == 0) THEN
         IF (NFBAND == 1) THEN
            ALLOCATE (INVP_ORD_BAND(1, NPAR))
            READ (invp_band_line, *, IOSTAT=ios) (INVP_ORD_BAND(1, j), j=1, NPAR)
            IF (ios /= 0) THEN
               DEALLOCATE (INVP_ORD_BAND)
               BACKSPACE (unit_main)
            ELSE
               NPAR_INV = 0
               DO j = 1, NPAR
                  IF (INVP_ORD_BAND(1, j) /= 0 .OR. INVP_ORD(j) /= 0) NPAR_INV = NPAR_INV + 1
               END DO
            END IF
         ELSE
            ALLOCATE (invp_band_tmp(NFBAND, NPAR))
            invp_band_tmp = 0

            READ (invp_band_line, *, IOSTAT=ios) (invp_band_tmp(1, j), j=1, NPAR)
            IF (ios == 0) THEN
               DO ib = 2, NFBAND - 1
                  READ (unit_main, *, IOSTAT=ios) (invp_band_tmp(ib, j), j=1, NPAR)
                  IF (ios /= 0) CALL fail('Read_Job_Setup: error reading INVP_ORD_BAND row #'//itoa(ib))
               END DO

               READ (unit_main, '(A)', IOSTAT=ios2) invp_band_line
               IF (ios2 == 0) THEN
                  READ (invp_band_line, *, IOSTAT=ios) (invp_band_tmp(NFBAND, j), j=1, NPAR)
                  IF (ios == 0) THEN
                     ! Older explicit-global layout:
                     ! global row already in INVP_ORD, rows read here are band 1..NFBAND.
                     ALLOCATE (INVP_ORD_BAND(NFBAND, NPAR))
                     INVP_ORD_BAND(:, :) = invp_band_tmp(:, :)
                  ELSE
                     ! Preferred simplified layout:
                     ! first row already stored in INVP_ORD and also becomes band 1.
                     BACKSPACE (unit_main)
                     ALLOCATE (INVP_ORD_BAND(NFBAND, NPAR))
                     INVP_ORD_BAND(1, :) = INVP_ORD(:)
                     INVP_ORD_BAND(2:NFBAND, :) = invp_band_tmp(1:NFBAND - 1, :)
                  END IF
               ELSE
                  ! EOF after simplified band rows.
                  ALLOCATE (INVP_ORD_BAND(NFBAND, NPAR))
                  INVP_ORD_BAND(1, :) = INVP_ORD(:)
                  INVP_ORD_BAND(2:NFBAND, :) = invp_band_tmp(1:NFBAND - 1, :)
               END IF

               NPAR_INV = 0
               DO j = 1, NPAR
                  IF (ANY(INVP_ORD_BAND(:, j) /= 0) .OR. INVP_ORD(j) /= 0) NPAR_INV = NPAR_INV + 1
               END DO
            ELSE
               DEALLOCATE (invp_band_tmp)
               BACKSPACE (unit_main)
            END IF

            IF (ALLOCATED(invp_band_tmp)) DEALLOCATE (invp_band_tmp)
         END IF
      END IF

      !--------------------------------------------------------------------
      ! (5) Solver choice
      ! label
      !   SOLVER_KIND
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) SOLVER_KIND
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading SOLVER_KIND')
      IF (SOLVER_KIND < 1 .OR. SOLVER_KIND > 3) CALL fail('Read_Job_Setup: SOLVER_KIND must be 1..3')

      !--------------------------------------------------------------------
      ! (6) Regularization
      ! label
      !   REG   lambda
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) REG, LAMBDA
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading REG, lambda')
      IF (REG < 0 .OR. REG > 1) CALL fail('Read_Job_Setup: REG must be 0/1')

      !--------------------------------------------------------------------
      ! (7) LBFGS controls
      ! label
      !   MML   LBFGS_TYPE
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) MML, lbfgs_type_loc
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading MML, LBFGS_TYPE')
      IF (MML <= 0) CALL fail('Read_Job_Setup: MML must be > 0')

      IF (PRESENT(LBFGS_TYPE)) LBFGS_TYPE = lbfgs_type_loc

      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, *) 'LBFGS: MML=', MML, ' type=', lbfgs_type_loc
      END IF

      !--------------------------------------------------------------------
      ! (8) Convergence controls
      ! label
      !   MAXITER   cost_conv_tol   model_conv_tol
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) MAXITER, cost_conv_tol, model_conv_tol
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading MAXITER, cost_conv_tol, model_conv_tol')
      IF (MAXITER <= 0) CALL fail('Read_Job_Setup: MAXITER must be > 0')

      !--------------------------------------------------------------------
      ! (9) Preconditioning (single switch + beta)
      ! label
      !   PRECOND   beta
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) PRECOND, BETA
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading PRECOND, beta')
      IF (PRECOND < 0 .OR. PRECOND > 1) CALL fail('Read_Job_Setup: PRECOND must be 0/1')

      !--------------------------------------------------------------------
      ! (10) Gradient taper + smoothing
      ! label
      !   GZSMTH   z1[m]   z2[m]   SIGMAX   SIGMAZ
      ! z1,z2 are physical depths; SIGMAX/SIGMAZ are Gaussian sigma values in grid points
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) GZSMTH, GZ1, GZ2, sigma_x_grad_in, sigma_z_grad_in
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading GZSMTH, z1, z2, SIGMAX, SIGMAZ')
      IF (GZSMTH < 0 .OR. GZSMTH > 1) CALL fail('Read_Job_Setup: GZSMTH must be 0/1')
      IF (sigma_x_grad_in < 0.0_8) CALL fail('Read_Job_Setup: SIGMAX must be >= 0')
      IF (sigma_z_grad_in < 0.0_8) CALL fail('Read_Job_Setup: SIGMAZ must be >= 0')
      SIGMA_X_GRAD = sigma_x_grad_in
      SIGMA_Z_GRAD = sigma_z_grad_in
      USE_GR_SMOOTH = (GZSMTH == 1 .AND. (SIGMA_X_GRAD > 0.0_8 .OR. SIGMA_Z_GRAD > 0.0_8))

      !--------------------------------------------------------------------
      ! (11) Global grid: NORD DX DZ
      ! label
      !   NORD   DX   DZ
      !--------------------------------------------------------------------
      READ (unit_main, 50) LAB
      READ (unit_main, *, IOSTAT=ios) NORD, DX, DZ
      IF (ios /= 0) CALL fail('Read_Job_Setup: error reading NORD, DX, DZ')
      IF (NORD < 2 .OR. NORD > 10) CALL fail('Read_Job_Setup: NORD out of range')
      IF (DX <= 0.0_8 .OR. DZ <= 0.0_8) CALL fail('Read_Job_Setup: DX/DZ must be > 0')

      !--------------------------------------------------------------------
      ! (11b) Optional per-band meshing block:
      ! label line containing 'Band meshing'
      ! then NFBAND lines:
      !   NORD_i   DX_i   DZ_i
      !--------------------------------------------------------------------
      READ (unit_main, 50, IOSTAT=ios2) LAB

      IF (ios2 == 0) THEN
         IF (INDEX(ADJUSTL(LAB), 'Band meshing') > 0) THEN
            ! We really have a "Band meshing" block: allocate and read it
            IF (PRESENT(NORD_BAND)) THEN
               IF (.NOT. ALLOCATED(NORD_BAND)) ALLOCATE (NORD_BAND(NFBAND))
            END IF
            IF (PRESENT(DX_BAND)) THEN
               IF (.NOT. ALLOCATED(DX_BAND)) ALLOCATE (DX_BAND(NFBAND))
            END IF
            IF (PRESENT(DZ_BAND)) THEN
               IF (.NOT. ALLOCATED(DZ_BAND)) ALLOCATE (DZ_BAND(NFBAND))
            END IF

            DO ib = 1, NFBAND
               READ (unit_main, *, IOSTAT=ios) nloc, dxloc, dzloc
               IF (ios /= 0) CALL fail('Read_Job_Setup: error reading band meshing line #'//itoa(ib))
               IF (PRESENT(NORD_BAND)) NORD_BAND(ib) = nloc
               IF (PRESENT(DX_BAND)) DX_BAND(ib) = dxloc
               IF (PRESENT(DZ_BAND)) DZ_BAND(ib) = dzloc
            END DO

         ELSE
            ! We read a line that actually belongs to the NEXT section (e.g. topography).
            ! Push the record pointer back so the next routine sees this label.
            BACKSPACE (unit_main)
         END IF

      ELSE
         ! ios2 /= 0: probably EOF or error; no band meshing and nothing to backspace.
         ! You may want to handle ios2 < 0 or > 0 explicitly if needed.
      END IF

      !--------------------------------------------------------------------
      ! (12) Per-band quadrature sampling early check
      !--------------------------------------------------------------------
      IF (.NOT. ALLOCATED(NFQ_PER_BAND)) THEN
         CALL fail('Read_Job_Setup: NFQ_PER_BAND not allocated; banded setup required.')
      END IF
      IF (NFBAND <= 0) CALL fail('Read_Job_Setup: NFBAND must be >= 1')

      ALLOCATE (band_start(NFBAND), band_end(NFBAND))
      sumq = 0
      DO ib = 1, NFBAND
         band_start(ib) = sumq + 1
         band_end(ib) = sumq + NFQ_PER_BAND(ib)
         sumq = sumq + NFQ_PER_BAND(ib)
      END DO
      IF (sumq /= NFQ) CALL fail('Read_Job_Setup: Sum(NFQ_PER_BAND) != NFQ')

      DO ib = 1, NFBAND
         i1b = band_start(ib)
         i2b = band_end(ib)
         IF (i2b < i1b) CALL fail('Read_Job_Setup: empty band encountered')

         FREQ_min_b = MINVAL(FREQN(i1b:i2b))
         FREQ_max_b = MAXVAL(FREQN(i1b:i2b))

         IF (PRESENT(NORD_BAND) .AND. ALLOCATED(NORD_BAND)) THEN
            NORD_use = NORD_BAND(ib)
         ELSE
            NORD_use = NORD
         END IF

         IF (PRESENT(DX_BAND) .AND. ALLOCATED(DX_BAND)) THEN
            DX_use = DX_BAND(ib)
         ELSE
            DX_use = DX
         END IF

         IF (PRESENT(DZ_BAND) .AND. ALLOCATED(DZ_BAND)) THEN
            DZ_use = DZ_BAND(ib)
         ELSE
            DZ_use = DZ
         END IF

         IF (NORD_use < 2 .OR. NORD_use > 10) CALL fail('Read_Job_Setup: band '//TRIM(ADJUSTL(itoa(ib)))//' invalid NORD')
         IF (DX_use <= 0.0_8 .OR. DZ_use <= 0.0_8) CALL fail('Read_Job_Setup: band '//TRIM(ADJUSTL(itoa(ib)))//' invalid DX/DZ')

         CALL CheckQuadratureSampling(CMIN, CMAX, FREQ_min_b, FREQ_max_b, &
                                      DX_use, DZ_use, NORD_use, WLMAX_b, WLMIN_b, my_rank)

         IF (dbg .AND. my_rank == 0) THEN
            WRITE (*, *) 'Band ', ib, '  fmin=', FREQ_min_b, '  fmax=', FREQ_max_b, &
               '  DX=', DX_use, '  DZ=', DZ_use, '  NORD=', NORD_use
         END IF
      END DO

      DEALLOCATE (band_start, band_end)

      IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Read_Job_Setup: done.'

      RETURN

   CONTAINS

      PURE FUNCTION itoa(i) RESULT(s)
         INTEGER, INTENT(IN) :: i
         CHARACTER(LEN=32)   :: s
         WRITE (s, '(I0)') i
      END FUNCTION itoa

      SUBROUTINE append_freqs(vec, newvals, nnew)
         REAL(KIND=8), ALLOCATABLE, INTENT(INOUT) :: vec(:)
         REAL(KIND=8), INTENT(IN)    :: newvals(:)
         INTEGER, INTENT(IN)    :: nnew
         REAL(KIND=8), ALLOCATABLE :: tmp(:)
         INTEGER :: oldn

         oldn = SIZE(vec)
         ALLOCATE (tmp(oldn + nnew))
         IF (oldn > 0) tmp(1:oldn) = vec
         tmp(oldn + 1:oldn + nnew) = newvals(1:nnew)
         CALL MOVE_ALLOC(tmp, vec)
      END SUBROUTINE append_freqs

   END SUBROUTINE Read_Job_Setup

!*****************************************************************************
!Topography inputs

!-----------------------------------------------------------------------
!
!     ReadTopographyData reads topography data from an input file and
!     computes the minimum and maximum coordinates for the model domain.
!
!     Inputs:
!       NTO............... Number of topography points
!       IS0............... Absorbing zone type (0 = free surface, 1= absorbing top and bottom
!
!     Outputs:
!       XTO(NTO).......... x-coordinates of the topography points
!       ZTO(NTO).......... z-coordinates of the topography points
!       XMIN, XMAX........ Minimum and maximum x-coordinates
!       YMIN, YMAX........ Minimum and maximum y-coordinates (defaulted to 0)
!       ZMIN, ZMAX........ Minimum and maximum z-coordinates
!
!-----------------------------------------------------------------------
   SUBROUTINE ReadTopographyData(NTO, IS0, XTO, ZTO, &
                                 XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, &
                                 my_rank, DEBUG_OUTPUT)
!-----------------------------------------------------------------------
! Reads topography points and computes domain bounds.
! All prints except the final bounds line are conditional on DEBUG_OUTPUT.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! -- arguments
      INTEGER, INTENT(OUT)    :: NTO, IS0
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XTO(:), ZTO(:)
      REAL(dp), INTENT(INOUT)  :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      INTEGER, INTENT(IN)     :: my_rank
      LOGICAL, OPTIONAL, INTENT(IN)     :: DEBUG_OUTPUT
      ! -- locals
      CHARACTER(LEN=80) :: LAB
      CHARACTER(LEN=256) :: topo_file
      CHARACTER(LEN=256) :: line, trimmed, next_line, first_point_line
      INTEGER           :: I, ios, ios2, unit_topo, nto_tmp, is0_tmp
      LOGICAL           :: dbg
      LOGICAL           :: use_external
      LOGICAL           :: have_prefetched_point, file_exists

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT
      have_prefetched_point = .FALSE.

50    FORMAT(80A)

      ! Header
      READ (1, 50, IOSTAT=ios) LAB
      if (my_rank == 0) WRITE (*, *) LAB
      IF (ios /= 0) THEN
         WRITE (*, *) 'Error reading topography header; IOSTAT=', ios
         STOP
      END IF
      IF (dbg .AND. my_rank == 0) WRITE (*, *) TRIM(LAB)

      call next_data_line(1, line, ios)
      IF (ios /= 0) THEN
         WRITE (*, *) 'Error reading topography spec line; IOSTAT=', ios
         STOP
      END IF

      trimmed = ADJUSTL(line)
      CALL parse_topo_spec(trimmed, NTO, IS0, ios2)
      use_external = .FALSE.

      IF (ios2 == 0) THEN
         ! Support both:
         !   1) inline points after "NTO IS0"
         !   2) external file path after "NTO IS0"
         call next_data_line(1, next_line, ios)
         IF (ios /= 0) THEN
            WRITE (*, *) 'Error reading first topography line after NTO/IS0; IOSTAT=', ios
            STOP
         END IF

         trimmed = ADJUSTL(next_line)
         topo_file = trimmed
         CALL strip_quotes_inplace(topo_file)
         INQUIRE (FILE=TRIM(topo_file), EXIST=file_exists)

         IF (file_exists) THEN
            use_external = .TRUE.
            IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography external file:', TRIM(topo_file)

            OPEN (NEWUNIT=unit_topo, FILE=TRIM(topo_file), STATUS='OLD', ACTION='READ', IOSTAT=ios)
            IF (ios /= 0) THEN
               WRITE (*, *) 'Error opening topography file: ', TRIM(topo_file), ' IOSTAT=', ios
               STOP
            END IF

            call next_data_line(unit_topo, line, ios)
            IF (ios /= 0) THEN
               WRITE (*, *) 'Error reading topography file header; IOSTAT=', ios
               STOP
            END IF

            CALL parse_topo_spec(line, nto_tmp, is0_tmp, ios2)
            IF (ios2 == 0) THEN
               NTO = nto_tmp
               IS0 = is0_tmp
               IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography external spec line:', TRIM(line)
            ELSE
               first_point_line = line
               have_prefetched_point = .TRUE.
               IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography uses main-file NTO/IS0 with external points file.'
            END IF
         ELSE
            unit_topo = 1
            first_point_line = next_line
            have_prefetched_point = .TRUE.
         END IF
      ELSE
         ! Path-only external file: the external file must provide its own "NTO IS0" header.
         topo_file = trimmed
         CALL strip_quotes_inplace(topo_file)
         IF (LEN_TRIM(topo_file) == 0) THEN
            WRITE (*, *) 'Error: empty topography filename'
            STOP
         END IF
         IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography external file:', TRIM(topo_file)

         OPEN (NEWUNIT=unit_topo, FILE=TRIM(topo_file), STATUS='OLD', ACTION='READ', IOSTAT=ios)
         IF (ios /= 0) THEN
            WRITE (*, *) 'Error opening topography file: ', TRIM(topo_file), ' IOSTAT=', ios
            STOP
         END IF

         call next_data_line(unit_topo, line, ios)
         IF (ios /= 0) THEN
            WRITE (*, *) 'Error reading topography file header; IOSTAT=', ios
            STOP
         END IF

         CALL parse_topo_spec(line, NTO, IS0, ios2)
         IF (ios2 /= 0) THEN
            WRITE (*, *) 'Error parsing topography header/spec: ', TRIM(line)
            IF (use_external) CLOSE (unit_topo)
            STOP
         END IF

         use_external = .TRUE.
         IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography external spec line:', TRIM(line)
      END IF
      IF ((NTO <= 0) .OR. (IS0 < 0)) THEN
         IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography invalid state from line:', TRIM(line)
         WRITE (*, *) 'Error with the topography data! (NTO<=0 or IS0<0)'
         IF (use_external) CLOSE (unit_topo)
         STOP
      END IF

      ! Absorbing flag normalization (0 => free surface => 1; otherwise => 2)
      IF (IS0 == 0) THEN
         IS0 = 1
      ELSE
         IS0 = 2
      END IF
      IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Topography points:', NTO, '  IS0:', IS0

      ! Allocate & init bounds
      ALLOCATE (XTO(NTO), ZTO(NTO))
      XMIN = +1.0e10_dp; XMAX = -1.0e10_dp
      YMIN = +1.0e10_dp; YMAX = -1.0e10_dp
      ZMIN = +1.0e10_dp; ZMAX = -1.0e10_dp

      ! Read points & update bounds
      DO I = 1, NTO
         IF (have_prefetched_point) THEN
            line = first_point_line
            have_prefetched_point = .FALSE.
         ELSE
            call next_data_line(unit_topo, line, ios)
            IF (ios /= 0) THEN
               WRITE (*, *) 'Error reading topography point line #', I, '; IOSTAT=', ios
               IF (use_external) CLOSE (unit_topo)
               STOP
            END IF
         END IF
         READ (line, *, IOSTAT=ios2) XTO(I), ZTO(I)
         IF (ios2 /= 0) THEN
            WRITE (*, *) 'Error parsing topography point #', I, ': ', TRIM(line), ' IOSTAT=', ios2
            IF (use_external) CLOSE (unit_topo)
            STOP
         END IF

         IF (dbg .AND. my_rank == 0) THEN
            WRITE (*, '(A,I6,2(1X,A,F12.3))') 'Topography point:', I, 'X=', XTO(I), 'Z=', ZTO(I)
         END IF

         IF (XTO(I) < XMIN) XMIN = XTO(I)
         IF (XTO(I) > XMAX) XMAX = XTO(I)
         IF (ZTO(I) < ZMIN) ZMIN = ZTO(I)
         IF (ZTO(I) > ZMAX) ZMAX = ZTO(I)
      END DO

      ! Final bounds line (always printed on rank 0)
      IF (my_rank == 0) THEN
         WRITE (*, *)
         WRITE (*, '(A,F10.2,1X,A,F10.2,1X,A,F10.2,1X,A,F10.2)') &
            'XMIN =', XMIN, 'XMAX =', XMAX, 'ZMIN =', ZMIN, 'ZMAX =', ZMAX
      END IF

      IF (use_external) CLOSE (unit_topo)

   CONTAINS

      SUBROUTINE next_data_line(unit_in, line_out, ios_out)
         INTEGER, INTENT(IN) :: unit_in
         CHARACTER(LEN=*), INTENT(OUT) :: line_out
         INTEGER, INTENT(OUT) :: ios_out
         CHARACTER(LEN=LEN(line_out)) :: raw

         DO
            READ (unit_in, '(A)', IOSTAT=ios_out) raw
            IF (ios_out /= 0) THEN
               line_out = ''
               RETURN
            END IF
            raw = strip_comment(raw)
            IF (LEN_TRIM(ADJUSTL(raw)) == 0) CYCLE
            line_out = raw
            RETURN
         END DO
      END SUBROUTINE next_data_line

      SUBROUTINE parse_topo_spec(spec_line, nto_out, is0_out, ios_out)
         CHARACTER(LEN=*), INTENT(IN) :: spec_line
         INTEGER, INTENT(OUT) :: nto_out, is0_out, ios_out

         nto_out = 0
         is0_out = 0
         READ (spec_line, *, IOSTAT=ios_out) nto_out, is0_out
         IF (ios_out /= 0) THEN
            nto_out = 0
            is0_out = 0
            READ (spec_line, *, IOSTAT=ios_out) nto_out
         END IF
      END SUBROUTINE parse_topo_spec

      FUNCTION strip_comment(s) RESULT(out)
         CHARACTER(LEN=*), INTENT(IN) :: s
         CHARACTER(LEN=LEN(s)) :: out
         INTEGER :: bang

         out = s
         bang = INDEX(out, '!')
         IF (bang == 1) THEN
            out = ''
         ELSEIF (bang > 1) THEN
            out = out(1:bang - 1)
         END IF
      END FUNCTION strip_comment

      SUBROUTINE strip_quotes_inplace(s)
         CHARACTER(LEN=*), INTENT(INOUT) :: s
         INTEGER :: n

         s = ADJUSTL(TRIM(s))
         n = LEN_TRIM(s)
         IF (n >= 2) THEN
            IF ((s(1:1) == '''' .AND. s(n:n) == '''') .OR. &
                (s(1:1) == '"' .AND. s(n:n) == '"')) THEN
               s = s(2:n - 1)
            END IF
         END IF
         s = ADJUSTL(TRIM(s))
      END SUBROUTINE strip_quotes_inplace
   END SUBROUTINE ReadTopographyData

!**********************************************************
! unified SR input: auto-detect synth-pattern vs real listing/file
! (no extra switch on NSR/NSS/ISR90 line)
   SUBROUTINE input_SR(NSR, NSS, NRR, XSR, YSR, ZSR, ICSR, VSR, &
                       ISR90, NCOMPS, NCOMPR, NCOMP, &
                       XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, &
                       my_rank, DEBUG_OUTPUT)

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: NSR, NSS, NRR, ISR90, NCOMPS, NCOMPR, NCOMP
      REAL(dp), INTENT(IN)  :: XMIN, XMAX, ZMIN, ZMAX
      REAL(dp), INTENT(INOUT) :: YMIN, YMAX
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XSR(:), YSR(:), ZSR(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: ICSR(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: VSR(:, :, :)
      INTEGER, INTENT(IN)  :: my_rank
      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT

      CHARACTER(LEN=256) :: header_line, nsr_line, nsr_parse, probe_line, probe_trim
      INTEGER :: ios, nsr_probe, nss_probe, isr90_probe
      INTEGER :: nprobe, ib, bang
      LOGICAL :: use_synth, dbg

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT
      use_synth = .FALSE.

      READ (1, '(A)', IOSTAT=ios) header_line
      IF (ios /= 0) STOP 'input_SR: failed reading SR header line'

      READ (1, '(A)', IOSTAT=ios) nsr_line
      IF (ios /= 0) STOP 'input_SR: failed reading NSR/NSS/ISR90 line'

      nsr_parse = nsr_line
      bang = INDEX(nsr_parse, '!')
      IF (bang == 1) THEN
         nsr_parse = ''
      ELSEIF (bang > 1) THEN
         nsr_parse = nsr_parse(1:bang - 1)
      END IF

      READ (nsr_parse, *, IOSTAT=ios) nsr_probe, nss_probe, isr90_probe
      IF (ios /= 0) STOP 'input_SR: cannot parse NSR/NSS/ISR90 line'

      nprobe = 0
      DO
         READ (1, '(A)', IOSTAT=ios) probe_line
         IF (ios /= 0) STOP 'input_SR: failed probing SR format line'
         nprobe = nprobe + 1
         probe_trim = ADJUSTL(probe_line)
         IF (LEN_TRIM(probe_trim) == 0) CYCLE
         IF (probe_trim(1:1) == '!') CYCLE
         EXIT
      END DO
      use_synth = is_synth_pattern_line(probe_trim)

      DO ib = 1, nprobe
         BACKSPACE (1)
      END DO

      BACKSPACE (1)
      BACKSPACE (1)

      IF (my_rank == 0 .AND. dbg) THEN
         WRITE (*, '(A,L1)') 'input_SR: use_synth=', use_synth
      END IF

      IF (use_synth) THEN
         CALL input_SR_synth(NSR, NSS, NRR, XSR, YSR, ZSR, ICSR, VSR, &
                             ISR90, NCOMPS, NCOMPR, NCOMP, &
                             XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, &
                             my_rank, DEBUG_OUTPUT)
      ELSE
         CALL input_SR_real(NSR, NSS, NRR, XSR, YSR, ZSR, ICSR, VSR, &
                            XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, ISR90, &
                            NCOMPS, NCOMPR, NCOMP, my_rank, DEBUG_OUTPUT)
      END IF

   END SUBROUTINE input_SR

   LOGICAL FUNCTION is_synth_pattern_line(line)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: line
      CHARACTER(LEN=LEN(line)) :: t
      CHARACTER(1) :: c1
      INTEGER :: ntok

      t = ADJUSTL(line)
      IF (LEN_TRIM(t) == 0) THEN
         is_synth_pattern_line = .TRUE.
         RETURN
      END IF

      c1 = t(1:1)

      IF ((c1 < '0' .OR. c1 > '9') .AND. c1 /= '-' .AND. c1 /= '+' .AND. c1 /= '.') THEN
         is_synth_pattern_line = .FALSE.
         RETURN
      END IF

      ntok = count_tokens_before_comment(t)
      is_synth_pattern_line = (ntok < 7)
   END FUNCTION is_synth_pattern_line

   INTEGER FUNCTION count_tokens_before_comment(line)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: line
      CHARACTER(LEN=LEN(line)) :: t
      INTEGER :: i, l, bang
      LOGICAL :: in_token

      t = line
      bang = INDEX(t, '!')
      IF (bang == 1) THEN
         t = ''
      ELSEIF (bang > 1) THEN
         t = t(1:bang - 1)
      END IF

      l = LEN_TRIM(t)
      count_tokens_before_comment = 0
      in_token = .FALSE.

      DO i = 1, l
         IF (t(i:i) /= ' ' .AND. t(i:i) /= ',' .AND. t(i:i) /= CHAR(9)) THEN
            IF (.NOT. in_token) count_tokens_before_comment = count_tokens_before_comment + 1
            in_token = .TRUE.
         ELSE
            in_token = .FALSE.
         END IF
      END DO
   END FUNCTION count_tokens_before_comment

!**********************************************************
!input for synthetic FWI

   SUBROUTINE input_SR_synth(NSR, NSS, NRR, XSR, YSR, ZSR, ICSR, VSR, &
                             ISR90, NCOMPS, NCOMPR, NCOMP, &
                             XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, &
                             my_rank, DEBUG_OUTPUT)
!-----------------------------------------------------------------------
!
!     input_SR_synth reads and computes the source and receiver
!     positions for synthetic forward modeling.
!
!     Inputs:
!       LAB............... Label for the input data
!       ISX_0, ISDX....... Initial x-coordinate and spacing for sources
!       ISZ_0, ISDZ....... Initial z-coordinate and spacing for sources
!       NSS............... Number of sources
!       IRX_B_0, IRDXB.... Initial x-coordinate and spacing for bottom receivers
!       IRZ_B_0, IRDZB.... Initial z-coordinate and spacing for bottom receivers
!       N_BR.............. Number of bottom receivers
!       IRX_T_0, IRDXT.... Initial x-coordinate and spacing for top receivers
!       IRZ_T_0, IRDZT.... Initial z-coordinate and spacing for top receivers
!       N_TR.............. Number of top receivers
!       IRX_L_0, IRDXL.... Initial x-coordinate and spacing for left receivers
!       IRZ_L_0, IRDZL.... Initial z-coordinate and spacing for left receivers
!       N_LR.............. Number of left receivers
!       IRX_R_0, IRDXR.... Initial x-coordinate and spacing for right receivers
!       IRZ_R_0, IRDZR.... Initial z-coordinate and spacing for right receivers
!       N_RR.............. Number of right receivers
!       XMIN, XMAX........ Lateral bounds
!       YMIN, YMAX........ Lateral bounds (y) not used
!       ZMIN, ZMAX........ Depth bounds (z)
!       my_rank........... MPI rank (only rank 0 prints output)
!
!     Outputs:
!       NSR............... Total number of sources and receivers
!       XSR(NSR).......... x-coordinates of sources and receivers
!       YSR(NSR).......... y-coordinates of sources and receivers
!       ZSR(NSR).......... z-coordinates of sources and receivers
!       ICSR(NSR)......... Source/receiver type indices
!       VSR(NSR,3,3)...... Source/receiver vectors
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      !-------------------- arguments --------------------
      INTEGER, INTENT(OUT) :: NSR, NSS, NRR, ISR90, NCOMPR, NCOMPS, NCOMP
      REAL(dp), INTENT(IN)  :: XMIN, XMAX, ZMIN, ZMAX
      REAL(dp), INTENT(INOUT) :: YMIN, YMAX
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XSR(:), YSR(:), ZSR(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: ICSR(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: VSR(:, :, :)
      INTEGER, INTENT(IN)  :: my_rank
      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT

      !-------------------- locals -----------------------
      CHARACTER(LEN=80) :: LAB
      INTEGER :: I, J, K, ios
      LOGICAL :: dbg
      ! source pattern
      REAL(dp) :: ISX_0, ISDX, ISZ_0, ISDZ
      ! receivers per side
      REAL(dp) :: IRX_B_0, IRDXB, IRZ_B_0, IRDZB
      REAL(dp) :: IRX_T_0, IRDXT, IRZ_T_0, IRDZT
      REAL(dp) :: IRX_L_0, IRDXL, IRZ_L_0, IRDZL
      REAL(dp) :: IRX_R_0, IRDXR, IRZ_R_0, IRDZR
      INTEGER :: N_BR, N_TR, N_LR, N_RR_side
      CHARACTER(1) :: mark

50    FORMAT(80A)

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      ! Read a label line from input file (unit 1)
      READ (1, 50, IOSTAT=ios) LAB; IF (ios /= 0) STOP 'input_SR_synth: read header label failed'
      IF (my_rank == 0 .AND. dbg) WRITE (*, *) LAB
      READ (1, *, IOSTAT=ios) NSR, NSS, ISR90; IF (ios /= 0) STOP 'input_SR_synth: read NSR/NSS/ISR90 failed'
      if (my_rank == 0 .AND. dbg) write (*, *) 'NSR, NSS', NSR, NSS, ISR90
      READ (1, *, IOSTAT=ios) NCOMPS, NCOMPR; IF (ios /= 0) STOP 'input_SR_synth: read NCOMPS/NCOMPR failed'
      if (my_rank == 0 .AND. dbg) write (*, *) 'NCOMP', NCOMPS, NCOMPR
      NCOMP = NCOMPS*NCOMPR
      IF (NSS .GE. NSR) THEN
         if (my_rank == 0) WRITE (*, *) '  WRONG SOURCE-RECEIVER NUMBER !'
         STOP
      END IF

      ! Read source parameters
      READ (1, *, IOSTAT=ios) ISX_0, ISDX, ISZ_0, ISDZ, NSS; IF (ios /= 0) STOP 'input_SR_synth: read source pattern failed'
      if (my_rank == 0 .AND. dbg) write (*, *) 'source setup', ISX_0, ISDX, ISZ_0, ISDZ, NSS
      ! Read receiver positions for different locations
      READ (1, *, IOSTAT=ios) IRX_T_0, IRDXT, IRZ_T_0, IRDZT, N_TR; IF (ios /= 0) STOP 'input_SR_synth: read 1d recvs failed'
      if (my_rank == 0 .AND. dbg) write (*, *) '2d', IRX_T_0, IRDXT, IRZ_T_0, IRDZT, N_TR
      READ (1, *, IOSTAT=ios) IRX_L_0, IRDXL, IRZ_L_0, IRDZL, N_LR; IF (ios /= 0) STOP 'input_SR_synth: read 2d recvs failed'
      if (my_rank == 0 .AND. dbg) write (*, *) '3rd', IRX_L_0, IRDXL, IRZ_L_0, IRDZL, N_LR
      READ (1, *, IOSTAT=ios) IRX_R_0, IRDXR, IRZ_R_0, IRDZR, N_RR_side; IF (ios /= 0) STOP 'input_SR_synth: 3t right recvs failed'
      if (my_rank == 0 .AND. dbg) write (*, *) '4th', IRX_R_0, IRDXR, IRZ_R_0, IRDZR, N_RR_side
      READ (1, *, IOSTAT=ios) IRX_B_0, IRDXB, IRZ_B_0, IRDZB, N_BR; IF (ios /= 0) STOP 'input_SR_synth: read 4s recvs failed'
      if (my_rank == 0 .AND. dbg) write (*, *) '1st', IRX_B_0, IRDXB, IRZ_B_0, IRDZB, N_BR
      ! Compute total number of receivers and sources+receivers
      NRR = N_BR + N_TR + N_LR + N_RR_side
      NSR = NRR + NSS
      if (NSR /= NRR + NSS) THEN
         if (my_rank == 0) write (*, *) 'NSR inconsistent'
         STOP
      END IF
      IF (my_rank == 0 .AND. dbg) WRITE (*, *) 'Src nbr: ', NSS, '  RCV nbr: ', NRR, '  NSR total: ', NSR, 'NCOMPS=', NCOMPS, 'NCOMPR=', NCOMPR

      ! Allocate arrays to store source and receiver positions
      ALLOCATE (XSR(NSR), YSR(NSR), ZSR(NSR), ICSR(NSR), VSR(NSR, 3, 3))
      XSR = 0.0_dp
      YSR = 0.0_dp
      ZSR = 0.0_dp
      ICSR = 0
      VSR = 0.0_dp
      ! y-coordinates are zero for all sources and receivers in this setup
      YSR(:) = 0.0_dp

      ! Print the total number of sources and receivers

      DO I = 1, NSR
      IF (I <= NSS) THEN
         XSR(I) = ISDX*(I - 1) + ISX_0
         ZSR(I) = ISDZ*(I - 1) + ISZ_0
         ICSR(I) = NCOMPS
      ELSE IF ((N_TR > 0) .AND. (I >= NSS + 1) .AND. (I <= (NSS + N_TR))) THEN
         XSR(I) = IRDXT*(I - NSS - 1) + IRX_T_0
         ZSR(I) = IRDZT*(I - NSS - 1) + IRZ_T_0
         ICSR(I) = NCOMPR
      ELSE IF ((N_LR > 0) .AND. (I >= (NSS + N_TR + 1)) .AND. (I <= (NSS + N_TR + N_LR))) THEN
         XSR(I) = IRDXL*(I - NSS - N_TR - 1) + IRX_L_0
         ZSR(I) = IRDZL*(I - NSS - N_TR - 1) + IRZ_L_0
         ICSR(I) = NCOMPR
      ELSE IF ((N_RR_side > 0) .AND. (I >= (NSS + N_TR + N_LR + 1)) .AND. (I <= (NSS + N_TR + N_LR + N_RR_side))) THEN
         XSR(I) = IRDXR*(I - NSS - N_TR - N_LR - 1) + IRX_R_0
         ZSR(I) = IRDZR*(I - NSS - N_TR - N_LR - 1) + IRZ_R_0
         ICSR(I) = NCOMPR
      ELSE IF ((N_BR > 0) .AND. (I >= (NSS + N_TR + N_LR + N_RR_side + 1))) THEN
         XSR(I) = IRDXB*(I - NSS - N_TR - N_LR - N_RR_side - 1) + IRX_B_0
         ZSR(I) = IRDZB*(I - NSS - N_TR - N_LR - N_RR_side - 1) + IRZ_B_0
         ICSR(I) = NCOMPR
      END IF

      VSR(I, :, :) = 0.0_dp

      SELECT CASE (ICSR(I))
      CASE (1)
         ! z-only
         VSR(i, 1, 1:3) = (/0.0_dp, 0.0_dp, 0.0_dp/)   ! x row empty
         VSR(i, 2, 1:3) = (/0.0_dp, 0.0_dp, 0.0_dp/)   ! y row empty
         VSR(i, 3, 1:3) = (/0.0_dp, 0.0_dp, 1.0_dp/)   ! z row
      CASE (2)
         ! x and z (no y)
         VSR(i, 1, 1:3) = (/1.0_dp, 0.0_dp, 0.0_dp/)   ! x
         VSR(i, 2, 1:3) = (/0.0_dp, 0.0_dp, 0.0_dp/)   ! y empty
         VSR(i, 3, 1:3) = (/0.0_dp, 0.0_dp, 1.0_dp/)   ! z
      CASE (3)
         ! full 3C
         VSR(i, 1, 1:3) = (/1.0_dp, 0.0_dp, 0.0_dp/)   ! x
         VSR(i, 2, 1:3) = (/0.0_dp, 1.0_dp, 0.0_dp/)   ! y
         VSR(i, 3, 1:3) = (/0.0_dp, 0.0_dp, 1.0_dp/)   ! z
      CASE DEFAULT
         STOP 'ICSR must be 1..3'
      END SELECT

      IF (YSR(I) < YMIN) YMIN = YSR(I)
      IF (YSR(I) > YMAX) YMAX = YSR(I)

      IF (XSR(I) < XMIN) THEN
         IF (my_rank == 0) WRITE (*, *) 'FAIL: XSR(', I, ') < XMIN'
      END IF
      IF (XSR(I) > XMAX) THEN
         IF (my_rank == 0) WRITE (*, *) 'FAIL: XSR(', I, ') > XMAX'
      END IF
      IF (ICSR(I) < 1) THEN
         IF (my_rank == 0) WRITE (*, *) 'FAIL: ICSR(', I, ') < 1'
      END IF
      IF (NSS > NSR) THEN
         IF (my_rank == 0) WRITE (*, *) 'FAIL: NSS > NSR'
      END IF

      ! If any failed, stop

      IF ((XSR(I) < XMIN) .OR. (XSR(I) > XMAX) .OR. (ICSR(I) < 1)) THEN
         if (my_rank == 0) WRITE (*, *) 'WRONG SOURCE/GEOPHONE LOCATIONS!'
         STOP
      END IF
      END DO

      IF (my_rank == 0 .AND. dbg) THEN
         DO I = 1, NSR
            WRITE (*, '(A, F12.6, 2X, A, F12.6)') 'XSR: ', XSR(I), 'ZSR: ', ZSR(I)
         END DO
         WRITE (*, '(A,F10.2,1X,A,F10.2,1X,A,F10.2,1X,A,F10.2)') &
            'XMIN =', XMIN, '  XMAX =', XMAX, '  ZMIN =', ZMIN, 'ZMAX =', ZMAX
      END IF

   END SUBROUTINE input_SR_synth

!-----------------------------------------------------------------------
!
!     input_SR_real reads real (predefined) source and receiver positions
!     from a file and initializes source/receiver vectors and bounds.
!
!     Inputs:
!       NSR................ Total number of sources and receivers
!       my_rank............ MPI rank (only rank 0 prints output)
!
!     Outputs:
!       NSS................ Number of sources
!       NRR................ Number of receivers
!       XSR(NSR)........... x-coordinates of sources and receivers
!       YSR(NSR)........... y-coordinates of sources and receivers
!       ZSR(NSR)........... z-coordinates of sources and receivers
!       ICSR(NSR).......... Number of independent vectors per SR point
!       VSR(NSR,3,3)....... Source/receiver vector definitions
!       XMIN, XMAX......... Lateral bounds
!       YMIN, YMAX......... Depth bounds (y)
!       ZMIN, ZMAX......... Depth bounds (z)
!
!-----------------------------------------------------------------------

   SUBROUTINE input_SR_real(NSR, NSS, NRR, XSR, YSR, ZSR, ICSR, VSR, &
                            XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, ISR90, &
                            NCOMPS, NCOMPR, NCOMP, my_rank, DEBUG_OUTPUT)

      IMPLICIT NONE

      !-------------------- arguments --------------------
      INTEGER, INTENT(OUT)   :: NSR, NSS, NRR, ISR90
      INTEGER, INTENT(OUT)   :: NCOMPS, NCOMPR, NCOMP
      REAL(dp), INTENT(IN)    :: XMIN, XMAX, ZMIN, ZMAX
      REAL(dp), INTENT(INOUT) :: YMIN, YMAX
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: XSR(:), YSR(:), ZSR(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: ICSR(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: VSR(:, :, :)     ! (NSR, 3, 3)
      INTEGER, INTENT(IN)    :: my_rank
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      !-------------------- locals -----------------------
      CHARACTER(LEN=80)   :: LAB
      CHARACTER(LEN=256)  :: line, SR_FILE
      INTEGER             :: i, j, ios, ios2, r
      REAL(dp)            :: dir(3, 3)  ! temp buffer (max 3 rows)
      LOGICAL             :: dbg
      CHARACTER(1)        :: mark
      INTEGER             :: unit_sr
      CHARACTER(LEN=256)  ::  trimmed
      CHARACTER(LEN=1)    :: c1
      LOGICAL             :: inline_mode, from_external

50    FORMAT(80A)

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT
      ! ---- Read header label on unit 1 ----
      READ (1, 50, IOSTAT=ios) LAB
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'input_SR_real: header label read failed, IOSTAT=', ios
         STOP
      END IF
      IF (my_rank == 0 .AND. dbg) WRITE (*, '(A)') TRIM(LAB)

      ! ---- NSR / NSS / ISR90 ----
      READ (1, *, IOSTAT=ios) NSR, NSS, ISR90
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'input_SR_real: NSR/NSS/ISR90 read failed, IOSTAT=', ios
         STOP
      END IF
      IF (my_rank == 0) WRITE (*, '(A,I0,2X,A,I0,2X,A,I0)') ' Total SR:', NSR, 'Sources:', NSS, 'ISR90:', ISR90

      IF (NSS >= NSR) THEN
         IF (my_rank == 0) WRITE (*, *) 'input_SR_real: WRONG SOURCE-RECEIVER NUMBER (NSS >= NSR)'
         STOP
      END IF

      ! ---- Allocate storage ----
      ALLOCATE (XSR(NSR), YSR(NSR), ZSR(NSR), ICSR(NSR), VSR(NSR, 3, 3))
      VSR = 0.0_dp

      !==========================================================
      ! Decide: inline SR block vs external SR file (NO internal READ)
      !==========================================================
      READ (1, '(A)', IOSTAT=ios) line
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'input_SR_real: failed to read SR data/file line. IOSTAT=', ios
         STOP
      END IF

      trimmed = ADJUSTL(line)

      ! Skip blank lines if needed
      IF (LEN_TRIM(trimmed) == 0) THEN
         ! could loop here until non-blank; for now assume not blank
      END IF

      c1 = trimmed(1:1)

      ! Heuristic: if first non-blank char is numeric-ish, assume inline SR row
      inline_mode = (c1 >= '0' .AND. c1 <= '9') .OR. c1 == '-' .OR. c1 == '+' .OR. c1 == '.'

      IF (my_rank == 0) THEN
         WRITE (*, '(A,A)') 'input_SR_real: SR line = "', TRIM(trimmed), '"'
         WRITE (*, '(A,L1)') 'input_SR_real: inline_mode? ', inline_mode
      END IF

      IF (inline_mode) THEN
         !-------------------------------------------------------
         ! INLINE MODE: line we just read is the first SR row
         !-------------------------------------------------------
         from_external = .FALSE.
         unit_sr = 1

         ! Put the line back so we can read it normally with list-directed READ
         BACKSPACE (unit_sr)

         DO i = 1, NSR
            READ (unit_sr, *, IOSTAT=ios) XSR(i), YSR(i), ZSR(i), ICSR(i), &
               (dir(1, j), dir(2, j), dir(3, j), j=1, ICSR(i))
              IF (my_rank == 0.AND. dbg) WRITE (*, *) 'DEBUG: read SR point ', i, ' XSR=', XSR(i), ' YSR=', YSR(i), ' ZSR=', ZSR(i), ' ICSR=', ICSR(i)
            IF (ios /= 0) THEN
               IF (my_rank == 0) WRITE (*, *) 'input_SR_real: SR read failed at i=', i, ' IOSTAT=', ios
               STOP
            END IF

            IF (ICSR(i) < 1 .OR. ICSR(i) > 3) THEN
               IF (my_rank == 0) WRITE (*, *) 'input_SR_real: ICSR(', i, ') out of range [1..3]: ', ICSR(i)
               STOP
            END IF

            VSR(i, :, :) = 0.0_dp
            SELECT CASE (ICSR(i))
            CASE (1)
               IF (.NOT. ANY(NINT(dir(1, 1:1)) == 0 .AND. NINT(dir(2, 1:1)) == 0 .AND. NINT(dir(3, 1:1)) == 1)) THEN
                  IF (my_rank == 0) WRITE (*, *) 'input_SR_real: ICSR=1 requires Z only at i=', i
                  STOP
               END IF
               VSR(i, 3, :) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            CASE (2)
               IF (.NOT. (ANY(NINT(dir(1, 1:2)) == 1 .AND. NINT(dir(2, 1:2)) == 0 .AND. NINT(dir(3, 1:2)) == 0) .AND. &
                          ANY(NINT(dir(1, 1:2)) == 0 .AND. NINT(dir(2, 1:2)) == 0 .AND. NINT(dir(3, 1:2)) == 1))) THEN
                  IF (my_rank == 0) WRITE (*, *) 'input_SR_real: ICSR=2 must be {X,Z} at i=', i
                  STOP
               END IF
               VSR(i, 1, :) = (/1.0_dp, 0.0_dp, 0.0_dp/)
               VSR(i, 3, :) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            CASE (3)
               VSR(i, 1, :) = (/1.0_dp, 0.0_dp, 0.0_dp/)
               VSR(i, 2, :) = (/0.0_dp, 1.0_dp, 0.0_dp/)
               VSR(i, 3, :) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            END SELECT

            ! Y-span and bounds checks as before
            IF (YSR(i) < YMIN) YMIN = YSR(i)
            IF (YSR(i) > YMAX) YMAX = YSR(i)

            IF (XSR(i) < XMIN) THEN
               IF (my_rank == 0) WRITE (*, *) 'FAIL: XSR(', i, ') < XMIN'
               STOP
            END IF
            IF (XSR(i) > XMAX) THEN
               IF (my_rank == 0) WRITE (*, *) 'FAIL: XSR(', i, ') > XMAX'
               STOP
            END IF
            IF (ZSR(i) > ZMAX) THEN
               IF (my_rank == 0) WRITE (*, *) 'FAIL: ZSR(', i, ') > ZMAX'
               STOP
            END IF
         END DO

      ELSE
         !-------------------------------------------------------
         ! EXTERNAL FILE MODE: "trimmed" is a path to SR file
         !-------------------------------------------------------
         from_external = .TRUE.
         SR_FILE = trimmed

         ! Strip quotes if present
         IF (LEN_TRIM(SR_FILE) >= 2) THEN
            IF (SR_FILE(1:1) == '''' .AND. SR_FILE(LEN_TRIM(SR_FILE):LEN_TRIM(SR_FILE)) == '''') THEN
               SR_FILE = SR_FILE(2:LEN_TRIM(SR_FILE) - 1)
            ELSEIF (SR_FILE(1:1) == '"' .AND. SR_FILE(LEN_TRIM(SR_FILE):LEN_TRIM(SR_FILE)) == '"') THEN
               SR_FILE = SR_FILE(2:LEN_TRIM(SR_FILE) - 1)
            END IF
         END IF

         IF (my_rank == 0) WRITE (*, '(A,A)') 'input_SR_real: using external SR file: ', TRIM(SR_FILE)

         OPEN (NEWUNIT=unit_sr, FILE=TRIM(SR_FILE), STATUS='OLD', ACTION='READ', IOSTAT=ios)
         IF (ios /= 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'input_SR_real: cannot open SR file. IOSTAT=', ios
            STOP
         END IF

         DO i = 1, NSR
            READ (unit_sr, *, IOSTAT=ios) XSR(i), YSR(i), ZSR(i), ICSR(i), &
               (dir(1, j), dir(2, j), dir(3, j), j=1, ICSR(i))
            IF (ios /= 0) THEN
               IF (my_rank == 0) WRITE (*, *) 'input_SR_real: SR read failed at i=', i, ' from SR file; IOSTAT=', ios
               STOP
            END IF

            IF (ICSR(i) < 1 .OR. ICSR(i) > 3) THEN
               IF (my_rank == 0) WRITE (*, *) 'input_SR_real: ICSR(', i, ') out of range [1..3]: ', ICSR(i)
               STOP
            END IF

            VSR(i, :, :) = 0.0_dp
            SELECT CASE (ICSR(i))
            CASE (1)
               IF (.NOT. ANY(NINT(dir(1, 1:1)) == 0 .AND. NINT(dir(2, 1:1)) == 0 .AND. NINT(dir(3, 1:1)) == 1)) THEN
                  IF (my_rank == 0) WRITE (*, *) 'input_SR_real: ICSR=1 requires Z only at i=', i
                  STOP
               END IF
               VSR(i, 3, :) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            CASE (2)
               IF (.NOT. (ANY(NINT(dir(1, 1:2)) == 1 .AND. NINT(dir(2, 1:2)) == 0 .AND. NINT(dir(3, 1:2)) == 0) .AND. &
                          ANY(NINT(dir(1, 1:2)) == 0 .AND. NINT(dir(2, 1:2)) == 0 .AND. NINT(dir(3, 1:2)) == 1))) THEN
                  IF (my_rank == 0) WRITE (*, *) 'input_SR_real: ICSR=2 must be {X,Z} at i=', i
                  STOP
               END IF
               VSR(i, 1, :) = (/1.0_dp, 0.0_dp, 0.0_dp/)
               VSR(i, 3, :) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            CASE (3)
               VSR(i, 1, :) = (/1.0_dp, 0.0_dp, 0.0_dp/)
               VSR(i, 2, :) = (/0.0_dp, 1.0_dp, 0.0_dp/)
               VSR(i, 3, :) = (/0.0_dp, 0.0_dp, 1.0_dp/)
            END SELECT

            ! same YMIN/YMAX and bounds checks here
            IF (YSR(i) < YMIN) YMIN = YSR(i)
            IF (YSR(i) > YMAX) YMAX = YSR(i)

            IF (XSR(i) < XMIN) THEN
               IF (my_rank == 0) WRITE (*, *) 'FAIL: XSR(', i, ') < XMIN'
               STOP
            END IF
            IF (XSR(i) > XMAX) THEN
               IF (my_rank == 0) WRITE (*, *) 'FAIL: XSR(', i, ') > XMAX'
               STOP
            END IF
            IF (ZSR(i) > ZMAX) THEN
               IF (my_rank == 0) WRITE (*, *) 'FAIL: ZSR(', i, ') > ZMAX'
               STOP
            END IF
         END DO

         CLOSE (unit_sr)

      END IF

      ! Receiver count & component summary
      NRR = NSR - NSS
      NCOMPS = MAXVAL(ICSR(1:NSS))
      NCOMPR = MAXVAL(ICSR(NSS + 1:NSR))
      NCOMP = NCOMPS*NCOMPR

      ! Optional full VSR dump
      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, '(A)') '--- Full VSR dump (all rows; * = active row, - = padded) ---'
         DO i = 1, NSR
            WRITE (*, '(A,I0,A,I0)') 'SR ', i, '  ICSR=', ICSR(i)
            DO r = 1, 3
               mark = MERGE('*', '-', r <= ICSR(i))
               WRITE (*, '(A,1X,3(I1,1X))') mark, NINT(VSR(i, r, 1)), NINT(VSR(i, r, 2)), NINT(VSR(i, r, 3))
            END DO
         END DO
      END IF

   END SUBROUTINE input_SR_real

!----------------------------------------------------------------------
!  SUBROUTINE Read_Grid_Geometry
!----------------------------------------------------------------------
!  Purpose:
!    Read geometry grid coordinates (XM(:), ZM(:)) from the main input,
!    given MX and MZ(:). Updates XMIN/XMAX/ZMIN/ZMAX.
!
!  Entries:
!    unit_main [in]   : file unit already positioned at the line AFTER the "coordinates" label
!    MX        [in]   : number of x-columns
!    MZ(MX)    [in]   : vertical counts per column
!
!    XMIN,XMAX,ZMIN,ZMAX [inout] : domain extents
!    XM(MX)    [out]  : x coordinates (1D, column headers)
!    ZM(:)     [out]  : stacked z coordinates (sum(MZ))
!
!  Return:
!    ZM is filled as [ZM(1:MZ(1)), ZM(MZ(1)+1: MZ(1)+MZ(2)), ...]
!----------------------------------------------------------------------

   SUBROUTINE Read_Grid_Geometry(MX, MZ, MZT, XM, ZM, XMIN, XMAX, ZMIN, ZMAX, my_rank, comm, DEBUG_OUTPUT)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: my_rank, comm
      REAL(dp), INTENT(INOUT) :: XMIN, XMAX, ZMIN, ZMAX
      INTEGER, INTENT(OUT), ALLOCATABLE :: MZ(:)
      REAL(dp), INTENT(OUT), ALLOCATABLE :: ZM(:), XM(:)
      INTEGER, INTENT(OUT) :: MX, MZT
      LOGICAL, OPTIONAL, INTENT(IN)     :: DEBUG_OUTPUT
      ! -- locals
      REAL(dp) :: X0, ZK
      INTEGER :: I, K, I0, ios, mpierr, unit_main
      CHARACTER(LEN=80) :: LAB
      LOGICAL           :: dbg

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

50    FORMAT(80A)
      unit_main = 1  ! main input file unit

      READ (unit_main, 50, IOSTAT=ios) LAB
      IF (ios /= 0) CALL abort_rank('read coords label @main', 701)
      IF (dbg) WRITE (*, *) TRIM(LAB)
      READ (unit_main, *, IOSTAT=ios) MX; IF (ios /= 0) CALL abort_rank('fail read MX @main', 701)
      if (dbg .and. my_rank == 0) WRITE (*, *) ' Grid MX = ', MX

      ALLOCATE (MZ(MX))
      READ (unit_main, *, IOSTAT=ios) (MZ(I), I=1, MX); IF (ios /= 0) CALL fail('read fail MZ @main')

      IF (dbg) WRITE (*, *) ' Grid MZ = ', (MZ(I), I=1, MX)
!Coordinates from MAIN
      READ (unit_main, 50, IOSTAT=ios) LAB
      if (dbg .and. my_rank == 0) WRITE (*, *) TRIM(LAB)

      MZT = SUM(MZ)
      ALLOCATE (XM(MX), ZM(MZT))

      I0 = 0
      DO I = 1, MX
         READ (unit_main, *, IOSTAT=ios) XM(I), (ZM(I0 + K), K=1, MZ(I))
         if (dbg .and. my_rank == 0) WRITE (*, *) ' Debug: read column ', I, ' XM=', XM(I)
         IF (ios /= 0) CALL fail('fail read coords @main')

         if (dbg .and. my_rank == 0) WRITE (*, *) ' XM, MZ ', XM(I), MZ(I)
         IF (XM(I) .LT. XMIN) XMIN = XM(I)
         IF (XM(I) .GT. XMAX) XMAX = XM(I)
         DO K = 1, MZ(I)
            ZK = ZM(I0 + K)
            IF (ZK .LT. ZMIN) ZMIN = ZK
            IF (ZK .GT. ZMAX) ZMAX = ZK
         END DO
         I0 = I0 + MZ(I)
      END DO

      RETURN
   END SUBROUTINE Read_Grid_Geometry
!**********************************************************
!----------------------------------------------------------------------
!  SUBROUTINE Read_Model_Parameters
!----------------------------------------------------------------------
!  Purpose:
!    Read CR0/CRT (and CI0/CIT if IVISCO=1) according to anisotropy mode.
!    Applies unit scaling rules. Copies TRUE->INIT when appropriate.
!
!  Entries:
!    File units must be positioned at the first parameter label line:
!      unit_init  : initial (MAIN) parameter section
!      unit_true  : true parameter section (optional)
!      unit_q_init: initial Q section (optional)
!      unit_q_true: true Q section (optional)
!
!  Entries / Globals:
!    IANISO, ITHOM, IVISCO, NSR, INVP_ORD(:), MODEL(:)
!    MX, MZ(:) and stacked size MZ = SUM(MZ)
!
!  Return:
!    CR0(IANISO, MZ), CRT(IANISO, MZ), CMZ(IANISO, MZ), CIT(IANISO, MZ)
!----------------------------------------------------------------------

   SUBROUTINE Read_Model_Parameters(MX, MZ, have_true, have_true_q, IANISO, ITHOM, IVISCO, NSR, INVP_ORD, MODEL, &
                                    CR0, CRT, CI0, CIT, my_rank, comm, DEBUG_OUTPUT)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)  :: have_true, have_true_q
      INTEGER, INTENT(IN) ::  IANISO, ITHOM, IVISCO, NSR, my_rank, comm
      INTEGER, INTENT(IN) :: MZ(:), INVP_ORD(:), MX
      CHARACTER(LEN=3), INTENT(IN) :: MODEL(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: CRT(:, :), CIT(:, :), CR0(:, :), CI0(:, :)

      LOGICAL, INTENT(IN) :: DEBUG_OUTPUT

      ! locals
      CHARACTER(LEN=80) :: LAB
      REAL(dp) :: X0
      INTEGER :: I, K, IP, I0, ios, mpierr
      LOGICAL :: dbg
      INTEGER :: MZT, II
      INTEGER, PARAMETER :: unit_init = 2, unit_q_init = 3, unit_true = 4, unit_q_true = 5

50    FORMAT(80A)
      dbg = DEBUG_OUTPUT
      MZT = SUM(MZ)

      ! Param branches classification (same as before, but parenthesized)
      ! II=1: isotropic (rho, alp, bet); II=2: Thomsen; II=3: VTI/TTI in Cij; II=4: general >7
      if (my_rank == 0) WRITE (*, *) '      (4) anisotropic model'

      II = 0
      IF (IANISO == 3) II = 1
      IF (((IANISO == 6) .OR. (IANISO == 7)) .AND. (ITHOM == 1)) II = 2
      IF (((IANISO == 6) .OR. (IANISO == 7)) .AND. (ITHOM == 0)) II = 3
      IF (IANISO > 7) II = 4

      ALLOCATE (CRT(IANISO, MZT), CIT(IANISO, MZT))    ! Complex Real/Imag (true model)
      ALLOCATE (CR0(IANISO, MZT), CI0(IANISO, MZT))

      ! Zero-out targets
      CR0 = 0.0_dp; CRT = 0.0_dp; CI0 = 0.0_dp; CIT = 0.0_dp

      DO IP = 1, IANISO
         SELECT CASE (II)
         CASE (1)  ! isotropic (rho, alp, bet) + optional Q
            READ (unit_init, 50)
            if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading MAIN label'
            IF (have_true) THEN
               READ (unit_true, 50)    ! TRUE label (keep files aligned)
               if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading TRUE label'
            END IF

            I0 = 0
            DO I = 1, MX
               ! Initial from MAIN -> CR0
               READ (unit_init, *, IOSTAT=ios) X0, (CR0(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('fail read CR0 @InitC')
               if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CR0(IP, I0 + K), K=1, MZ(I))
               IF (have_true) THEN !TRUE -> CRT
                  READ (unit_true, *, IOSTAT=ios) X0, (CRT(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('fail read CRT @TrueC')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CRT(IP, I0 + K), K=1, MZ(I))
               END IF
               I0 = I0 + MZ(I)
            END DO
            IF (have_true .AND. (INVP_ORD(IP) .EQ. 0 .OR. NSR == 2)) THEN
               CRT(IP, 1:I0) = CR0(IP, 1:I0)
               ! CR0(IP, 1:I0) = CRT(IP, 1:I0)
            END IF

            ! Scale velocities (not density) if in km/s
            IF (IP > 1 .AND. CR0(IP, 1) < 10.0_dp) THEN
               CR0(IP, :) = 1.0e3_dp*CR0(IP, :)
               CRT(IP, :) = 1.0e3_dp*CRT(IP, :)
            END IF

            ! Imag part (Q) if visco
            IF (IVISCO /= 0 .AND. IP > 1) THEN
               READ (unit_q_init, 50) LAB
               if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading INIT label', LAB
               IF (have_true_q) THEN
                  READ (unit_q_true, 50) LAB
                  if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading TRUE label', LAB
               END IF

               I0 = 0
               DO I = 1, MX
                  READ (unit_q_init, *, IOSTAT=ios) X0, (CI0(IP, I0 + K), K=1, MZ(I)); 
                  IF (ios /= 0) CALL fail('read CI0 @InitnQ')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CI0(IP, I0 + K), K=1, MZ(I))

                  IF (have_true_q) THEN
                     READ (unit_q_true, *, IOSTAT=ios) X0, (CIT(IP, I0 + K), K=1, MZ(I)); 
                     IF (ios /= 0) CALL fail('read CIT @TrueQ')
                     if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CIT(IP, I0 + K), K=1, MZ(I))
                  END IF
                  I0 = I0 + MZ(I)
               END DO
               IF (have_true_q .AND. (INVP_ORD(IP + IANISO - 1) .EQ. 0 .OR. NSR == 2)) THEN
                  ! CIT(IP, 1:I0) = CI0(IP, 1:I0)
                  CI0(IP, 1:I0) = CIT(IP, 1:I0)
               END IF
            END IF

         CASE (2)  ! VTI/TTI: Thomsen or Cij branch (MAIN -> CR0, TRUE -> CRT)
            IF (ITHOM == 1) THEN
               !-------------------------
               ! Real part: MAIN / TRUE
               !-------------------------
               READ (unit_init, 50) LAB
               if (dbg .AND. my_rank == 0) WRITE (*, *) LAB
               IF (have_true) THEN
                  READ (unit_true, 50) LAB
                  if (dbg .AND. my_rank == 0) WRITE (*, *) LAB   ! TRUE label (keep files aligned)
               END IF

               I0 = 0
               DO I = 1, MX
                  ! Initial from init -> CR0
                  READ (unit_init, *, IOSTAT=ios) X0, (CR0(IP, I0 + K), K=1, MZ(I))
                  IF (ios /= 0) CALL fail('fail to read CR0  @init')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CR0(IP, I0 + K), K=1, MZ(I))

                  ! TRUE -> CRT
                  IF (have_true) THEN
                     READ (unit_true, *, IOSTAT=ios) X0, (CRT(IP, I0 + K), K=1, MZ(I))
                     IF (ios /= 0) CALL fail('fail to read CRT  @true')
                     if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CRT(IP, I0 + K), K=1, MZ(I))
                  END IF
                  I0 = I0 + MZ(I)
               END DO

               ! Copy rule (same logic as CASE 3)
               IF (have_true .AND. ((INVP_ORD(IP) .EQ. 0))) THEN
                  ! CRT(IP, 1:I0) = CR0(IP, 1:I0)
                  CR0(IP, 1:I0) = CRT(IP, 1:I0)   ! keeps non inverted parameter to true values
               else if (have_true .AND. (NSR == 2)) THEN
                  CRT(IP, 1:I0) = CR0(IP, 1:I0)

               END IF

               !---------------------------------------
               ! Thomsen scaling (keep original logic)
               !---------------------------------------
               IF (IP == 2 .OR. IP == 3 .AND. CR0(IP, 1) < 10D0) THEN
                  CR0(IP, :) = 1.0D3*CR0(IP, :)
                  CRT(IP, :) = 1.0D3*CRT(IP, :)
               END IF

               !-------------------------
               ! Imag part (Q) if visco
               !-------------------------
               IF (IVISCO == 1 .AND. ((IP .GT. 1) .AND. (IP .LT. 7))) THEN
                  READ (unit_q_init, 50) LAB
                  if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading MAIN label', LAB
                  IF (have_true_q) THEN
                     READ (unit_q_true, 50) LAB
                     if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading TRUE label', LAB
                  END IF

                  I0 = 0
                  DO I = 1, MX
                     ! initial Q -> CI0 from MAIN Q
                     READ (unit_q_init, *, IOSTAT=ios) X0, (CI0(IP, I0 + K), K=1, MZ(I))
                     IF (ios /= 0) CALL fail('read CI0 @InitQ')
                     if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CI0(IP, I0 + K), K=1, MZ(I))

                     ! true Q -> CIT from TRUE Q
                     IF (have_true_q) THEN
                        READ (unit_q_true, *, IOSTAT=ios) X0, (CIT(IP, I0 + K), K=1, MZ(I))
                        IF (ios /= 0) CALL fail('read CIT @TrueQ')
                        if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CIT(IP, I0 + K), K=1, MZ(I))
                     END IF
                     I0 = I0 + MZ(I)
                  END DO

                  IF (have_true_q .AND. (INVP_ORD(IP + IANISO - 1) .EQ. 0 .OR. NSR == 2)) THEN
                     ! CIT(IP, 1:I0) = CI0(IP, 1:I0)
                     CI0(IP, 1:I0) = CIT(IP, 1:I0)
                  END IF
               END IF
            END IF
            ! Cij: MAIN->CR0 (Pa), INIT->CRT (Pa)

         CASE (3) !VTI/TTI
            READ (unit_init, 50) LAB
            if (dbg .AND. my_rank == 0) WRITE (*, *) LAB
            IF (have_true) THEN
               READ (unit_true, 50) LAB
               if (dbg .AND. my_rank == 0) WRITE (*, *) LAB    ! TRUE label (keep files aligned)
            END IF
            I0 = 0
            DO I = 1, MX
               ! Initial from MAIN -> CR0
               READ (unit_init, *, IOSTAT=ios) X0, (CR0(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CR0 @init')
               if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CR0(IP, I0 + K), K=1, MZ(I))
               IF (have_true) THEN
                  READ (unit_true, *, IOSTAT=ios) X0, (CRT(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CRT @true')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CRT(IP, I0 + K), K=1, MZ(I))
               END IF
               I0 = I0 + MZ(I)
            END DO
            IF (have_true .AND. ((INVP_ORD(IP) .EQ. 0) .OR. (NSR == 2))) THEN
               ! CRT(IP, 1:I0) = CR0(IP, 1:I0)!freeze non inverted parameter to initial values
               CR0(IP, 1:I0) = CRT(IP, 1:I0)   ! keeps non inverted parameter to true values
            END IF

            ! Scale moduli (not density) if in km/s

            IF (IP > 1 .AND. IP < 7) THEN
               CR0(IP, :) = (10.0_dp**9)*CR0(IP, :)
               CRT(IP, :) = (10.0_dp**9)*CRT(IP, :)
            END IF

            ! Imag part (Q) if visco
            IF (IVISCO == 1 .AND. ((IP .GT. 1) .AND. (IP .LT. 7))) THEN
               READ (unit_q_init, 50) LAB
               if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading MAIN label', LAB
               IF (have_true) THEN
                  READ (unit_q_true, 50) LAB
                  if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading TRUE label', LAB

               END IF
               I0 = 0
               DO I = 1, MX
                  ! initial Q -> CI0 from MAIN Q
                  READ (unit_q_init, *, IOSTAT=ios) X0, (CI0(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CI0 @mainQ')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CI0(IP, I0 + K), K=1, MZ(I))
                  IF (have_true_q) THEN
                     READ (unit_q_true, *, IOSTAT=ios) X0, (CIT(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CIT @initQ')
                     if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CIT(IP, I0 + K), K=1, MZ(I))
                  END IF
                  I0 = I0 + MZ(I)
               END DO
               IF (have_true_q .AND. (INVP_ORD(IP + IANISO - 1) .EQ. 0 .OR. NSR == 2)) THEN
                  ! CIT(IP, 1:I0) = CI0(IP, 1:I0)
                  CI0(IP, 1:I0) = CIT(IP, 1:I0)
               END IF
            END IF

!===============================================================
! II = 3: General anisotropy
!===============================================================
         CASE (4)
            READ (unit_init, 50) LAB
            if (dbg .AND. my_rank == 0) WRITE (*, *) LAB
            IF (have_true) THEN
               READ (unit_true, 50) LAB
               if (dbg .AND. my_rank == 0) WRITE (*, *) LAB    ! TRUE label (keep files aligned)
            END IF
            I0 = 0
            DO I = 1, MX
               ! Initial from MAIN -> CR0
               READ (unit_init, *, IOSTAT=ios) X0, (CR0(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CR0 @main')
               if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CR0(IP, I0 + K), K=1, MZ(I))
               IF (have_true) THEN
                  READ (unit_true, *, IOSTAT=ios) X0, (CRT(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CRT @init')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CRT(IP, I0 + K), K=1, MZ(I))
               END IF
               I0 = I0 + MZ(I)
            END DO
            IF (have_true .AND. (INVP_ORD(IP) .EQ. 0 .OR. NSR == 2)) THEN
               CRT(IP, 1:I0) = CR0(IP, 1:I0)
               ! CR0(IP, 1:I0) = CRT(IP, 1:I0)
            END IF

            ! Scale moduli (not density) if in km/s

            IF (IP > 1 .AND. IP < 7) THEN
               CR0(IP, :) = 1.0D9*CR0(IP, :)
               CRT(IP, :) = 1.0D9*CRT(IP, :)
            END IF

            ! Imag part (Q) if visco
            IF (IVISCO == 1 .AND. ((IP .GT. 1) .AND. (IP .LT. 7))) THEN
               READ (unit_q_init, 50) LAB
               if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading MAIN label', LAB
               IF (have_true) THEN
                  READ (unit_q_true, 50) LAB
                  if (dbg .AND. my_rank == 0) WRITE (*, *) 'Reading TRUE label', LAB

               END IF
               I0 = 0
               DO I = 1, MX
                  ! initial Q -> CI0 from MAIN Q
                  READ (unit_q_init, *, IOSTAT=ios) X0, (CI0(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CI0 @mainQ')
                  if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CI0(IP, I0 + K), K=1, MZ(I))
                  IF (have_true_q) THEN
                     READ (unit_q_true, *, IOSTAT=ios) X0, (CIT(IP, I0 + K), K=1, MZ(I)); IF (ios /= 0) CALL fail('read CIT @initQ')
                     if (dbg .AND. my_rank == 0) write (*, *) I, X0, (CIT(IP, I0 + K), K=1, MZ(I))
                  END IF
                  I0 = I0 + MZ(I)
               END DO
               IF (have_true_q .AND. (INVP_ORD(IP + IANISO - 1) .EQ. 0 .OR. NSR == 2)) THEN
                  CIT(IP, 1:I0) = CI0(IP, 1:I0)
                  ! CI0(IP, 1:I0) = CIT(IP, 1:I0)
               END IF
            END IF
         CASE DEFAULT
            CALL abort_rank('Unsupported parameterization II', 200)

         END SELECT
      END DO

      RETURN
      !   CRT(IP, :) = CR0(IP, :)
!
! Effect:
!   The reference/observed model is forced to match the initial model
!   for frozen parameters.
!
! Consequence:
!   Frozen parameters generate no data residual contribution since:
!
!      d_obs = F(m_active_true , m_frozen_initial)
!      d_syn = F(m_active_initial , m_frozen_initial)
!
!   Therefore the residual is only driven by the active parameters.
!
! Usage:
!   Useful for controlled sensitivity/isolation tests where only the
!   selected parameter(s) should contribute to the misfit.
      !   CR0(IP, :) = CRT(IP, :)
!
! Effect:
!   The initial/inversion model is forced to match the true model
!   for frozen parameters.
!
! Consequence:
!   Frozen parameters are assumed perfectly known during inversion:
!
!      d_obs = F(m_true)
!      d_syn = F(m_active_initial , m_frozen_true)
!
!   Only the active parameters remain mismatched, which strongly reduces:
!      - parameter cross-talk
!      - leakage from frozen parameters
!      - line-search instability
!      - artificial compensation effects
!
! Usage:
!   Commonly used in single-parameter or staged multiparameter
!   inversion studies to evaluate recoverability under controlled
!   conditions.

   END SUBROUTINE Read_Model_Parameters

!**********************************************************
   SUBROUTINE Read_acquisition(NFQ, NSS, NRR, ICSR, VSR, ACQ_WEIGHT_FILE, &
                               NDATA, NS, NSV, NR, NRV, NCOMPS, NCOMPR, COMP, &
                               my_rank, DEBUG_OUTPUT)

      IMPLICIT NONE

      !----------------- Arguments -----------------
      INTEGER, INTENT(IN)    :: NFQ, NSS, NRR
      INTEGER, INTENT(IN)    :: my_rank, NCOMPS, NCOMPR
      INTEGER, INTENT(INOUT) :: ICSR(:)
      INTEGER, INTENT(INOUT) :: NDATA(:)
      INTEGER, INTENT(INOUT) :: NS(:), NSV(:), NR(:), NRV(:), COMP(:)
      REAL(dp), INTENT(IN)    :: VSR(:, :, :)
      CHARACTER(len=256), INTENT(OUT)   :: ACQ_WEIGHT_FILE
      LOGICAL, OPTIONAL, INTENT(IN)    :: DEBUG_OUTPUT

      !----------------- Locals -----------------
      CHARACTER(80)      :: LAB
      CHARACTER(256)     :: line
      CHARACTER(256)     :: ACQ_FILE
      INTEGER            :: ios, ios2
      INTEGER            :: IFQ, J, ND, nd_file
      INTEGER            :: shotid, ics, recid, icr
      INTEGER            :: sc, rc, IS, IR, scomp, rcomp, IDD
      LOGICAL            :: dbg
      REAL(dp)           :: Svec(3), Gvec(3)
      INTEGER            :: unit_acq          ! unit used for acquisition mapping
      LOGICAL            :: from_external_acq

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

50    FORMAT(80A)

      !----------------------------------------------------------
      ! (0) Header + acquisition weight file path (on unit 1)
      !----------------------------------------------------------
      READ (1, 50, IOSTAT=ios) LAB
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: header read failed (unit=1). IOSTAT=', ios
         STOP
      END IF
      IF (dbg .AND. my_rank == 0) WRITE (*, *) TRIM(LAB)

      ! Next line: acquisition weight file path (may be '-' or empty)
      READ (1, '(A)', IOSTAT=ios) line
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: failed to read ACQ_WEIGHT_FILE line. IOSTAT=', ios
         STOP
      END IF

      ACQ_WEIGHT_FILE = ADJUSTL(line)
      IF (TRIM(ACQ_WEIGHT_FILE) == '-' .OR. LEN_TRIM(ACQ_WEIGHT_FILE) == 0) THEN
         ACQ_WEIGHT_FILE = ''   ! force empty string as "no weight file"
      END IF

      !----------------------------------------------------------
      ! (1) Acquisition format dispatch based on the input line
      !----------------------------------------------------------
      READ (1, '(A)', IOSTAT=ios) line
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: failed to read acquisition ID line. IOSTAT=', ios
         STOP
      END IF

      READ (line, *, IOSTAT=ios2) IDD
      IF (ios2 == 0 .AND. IDD == 1) THEN
         ! Generated synthetic-style mapping
         IF (dbg .AND. my_rank == 0) WRITE (*, '(A)') 'Read_acquisition: using generated mapping pattern from ID line = 1'

         ND = NSS*NRR*NCOMPS*NCOMPR

         IF (SIZE(NS) < ND .OR. SIZE(NSV) < ND .OR. &
             SIZE(NR) < ND .OR. SIZE(NRV) < ND .OR. &
             SIZE(COMP) < ND) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: arrays too small for ND=', ND
            STOP
         END IF

         IF (SIZE(NS) > ND .OR. SIZE(NSV) > ND .OR. &
             SIZE(NR) > ND .OR. SIZE(NRV) > ND .OR. &
             SIZE(COMP) > ND) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: arrays too big for ND=', ND
            STOP
         END IF

         J = 0

         DO IS = 1, NSS
            DO sc = 1, NCOMPS
               ! map logical source component index
               scomp = MERGE(3, MERGE(2*sc - 1, sc, NCOMPS == 2), NCOMPS == 1)

               DO rc = 1, NCOMPR
                  ! map logical receiver component index
                  rcomp = MERGE(3, MERGE(2*rc - 1, rc, NCOMPR == 2), NCOMPR == 1)

                  DO IR = 1, NRR
                     J = J + 1

                     NS(J) = IS
                     NR(J) = NSS + IR
                     NSV(J) = scomp
                     NRV(J) = rcomp

                     ! cross-component gating rule
                     IF ((NSV(J) == 2 .AND. NRV(J) == 1) .OR. &
                         (NSV(J) == 1 .AND. NRV(J) == 2) .OR. &
                         (NSV(J) == 2 .AND. NRV(J) == 3) .OR. &
                         (NSV(J) == 3 .AND. NRV(J) == 2)) THEN
                        COMP(J) = 0
                     ELSE
                        COMP(J) = 1
                     END IF

                  END DO  ! IR
               END DO     ! rc
            END DO        ! sc
         END DO           ! IS

         IF (J /= ND) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: mapping size mismatch, J=', J, ' ND=', ND
            STOP
         END IF

         ! NS must be 1..NSS
         IF (MAXVAL(NS(1:ND)) > NSS .OR. MINVAL(NS(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NS out of range [1,NSS]'
            STOP 'Read_acquisition: NS index out of range'
         END IF

         ! NR must be NSS+1 .. NSS+NRR
         IF (MAXVAL(NR(1:ND)) > NSS + NRR .OR. MINVAL(NR(1:ND)) < NSS + 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NR out of range [NSS+1, NSS+NRR]'
            STOP 'Read_acquisition: NR index out of range'
         END IF

         ! NSV/NRV must be in {1,2,3}
         IF (MAXVAL(NSV(1:ND)) > 3 .OR. MINVAL(NSV(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NSV outside [1,3]'
            STOP 'Read_acquisition: NSV index out of range'
         END IF
         IF (MAXVAL(NRV(1:ND)) > 3 .OR. MINVAL(NRV(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NRV outside [1,3]'
            STOP 'Read_acquisition: NRV index out of range'
         END IF

         ! COMP must be 0 or 1
         IF (MAXVAL(COMP(1:ND)) > 1 .OR. MINVAL(COMP(1:ND)) < 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: COMP not in {0,1}'
            STOP 'Read_acquisition: invalid COMP flag'
         END IF

         ! VSR range consistency (VSR comes from elsewhere)
         IF (NSS + NRR > SIZE(VSR, 1)) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NSS+NRR=', NSS + NRR, &
               ' > SIZE(VSR,1)=', SIZE(VSR, 1)
            STOP 'Read_acquisition: VSR first dimension too small'
         END IF

         ! All NS, NR must be valid indices into VSR(:,*,*)
         IF (MAXVAL(NS(1:ND)) > SIZE(VSR, 1) .OR. MAXVAL(NR(1:ND)) > SIZE(VSR, 1)) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NS or NR exceeds SIZE(VSR,1)'
            STOP 'Read_acquisition: NS/NR vs VSR size mismatch'
         END IF
         IF (SIZE(VSR, 2) < 3 .OR. SIZE(VSR, 3) < 3) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: VSR second/third dims < 3'
            STOP 'Read_acquisition: VSR must be (N_s+N_r,3,3)'
         END IF
      ELSE IF (ios2 == 0) THEN
         ! Inline acquisition rows follow on the main input file
         ND = IDD
         from_external_acq = .FALSE.
         unit_acq = 1

         IF (dbg .AND. my_rank == 0) WRITE (*, '(A,I0)') 'Read_acquisition: using inline mapping with ND = ', ND

         IF (SIZE(NS) < ND .OR. SIZE(NSV) < ND .OR. &
             SIZE(NR) < ND .OR. SIZE(NRV) < ND .OR. &
             SIZE(COMP) < ND) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: arrays too small for ND=', ND
            STOP
         END IF

         IF (dbg .AND. my_rank == 0) WRITE (*, *) ' ND check', ND, NSS*NRR

         DO J = 1, ND
            READ (unit_acq, '(A)', IOSTAT=ios) line
            IF (ios /= 0) THEN
               IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: line read failed at row', J, ' IOSTAT=', ios
               STOP
            END IF

            ! try 5-int form: "ND Shot ICS Rec ICR"
            READ (line, *, IOSTAT=ios) nd_file, shotid, ics, recid, icr
            IF (ios == 0) THEN
               NS(J) = shotid
               NSV(J) = ics
               NR(J) = recid
               NRV(J) = icr
            ELSE
               ! fallback 4-int form: "Shot ICS Rec ICR"
               READ (line, *, IOSTAT=ios) NS(J), NSV(J), NR(J), NRV(J)
               IF (ios /= 0) THEN
                  IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: bad row ', J, ':', TRIM(line)
                  STOP
               END IF
            END IF

            IF (dbg .AND. my_rank == 0) WRITE (*, *) NS(J), NSV(J), NR(J), NRV(J)

            ! 2.5D keep/zero rule
            IF ((NSV(J) == 2 .AND. NRV(J) == 1) .OR. &
                (NSV(J) == 1 .AND. NRV(J) == 2) .OR. &
                (NSV(J) == 2 .AND. NRV(J) == 3) .OR. &
                (NSV(J) == 3 .AND. NRV(J) == 2)) THEN
               COMP(J) = 0
            ELSE
               COMP(J) = 1
            END IF
         END DO

         IF (MAXVAL(NS(1:ND)) > NSS .OR. MINVAL(NS(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NS out of range [1,NSS]'
            STOP 'Read_acquisition: NS index out of range'
         END IF

         IF (MAXVAL(NR(1:ND)) > NSS + NRR .OR. MINVAL(NR(1:ND)) < NSS + 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NR out of range [NSS+1,NSS+NRR]'
            STOP 'Read_acquisition: NR index out of range'
         END IF

         IF (MAXVAL(NSV(1:ND)) > 3 .OR. MINVAL(NSV(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NSV outside [1,3]'
            STOP 'Read_acquisition: NSV index out of range'
         END IF

         IF (MAXVAL(NRV(1:ND)) > 3 .OR. MINVAL(NRV(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NRV outside [1,3]'
            STOP 'Read_acquisition: NRV index out of range'
         END IF

         IF (MAXVAL(COMP(1:ND)) > 1 .OR. MINVAL(COMP(1:ND)) < 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: COMP not in {0,1}'
            STOP 'Read_acquisition: invalid COMP flag'
         END IF

         IF (NSS + NRR > SIZE(VSR, 1)) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NSS+NRR=', NSS + NRR, &
               ' > SIZE(VSR,1)=', SIZE(VSR, 1)
            STOP 'Read_acquisition: VSR first dimension too small'
         END IF

         IF (MAXVAL(NS(1:ND)) > SIZE(VSR, 1) .OR. MAXVAL(NR(1:ND)) > SIZE(VSR, 1)) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NS or NR exceeds SIZE(VSR,1)'
            STOP 'Read_acquisition: NS/NR vs VSR size mismatch'
         END IF

         IF (SIZE(VSR, 2) < 3 .OR. SIZE(VSR, 3) < 3) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: VSR second/third dims < 3'
            STOP 'Read_acquisition: VSR must be (N_s+N_r,3,3)'
         END IF
      ELSE
         ! External acquisition file
         from_external_acq = .TRUE.
         ACQ_FILE = ADJUSTL(line)
         IF (my_rank == 0) WRITE (*, '(A,A)') 'Read_acquisition: using external ACQ file: ', TRIM(ACQ_FILE)

         OPEN (NEWUNIT=unit_acq, FILE=TRIM(ACQ_FILE), STATUS='OLD', &
               ACTION='READ', IOSTAT=ios)
         IF (ios /= 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: cannot open ACQ file. IOSTAT=', ios
            STOP
         END IF

         ! First record in ACQ_FILE must be ND
         READ (unit_acq, *, IOSTAT=ios) ND
         IF (ios /= 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: ND read failed from ACQ file. IOSTAT=', ios
            STOP
         END IF

         IF (SIZE(NS) < ND .OR. SIZE(NSV) < ND .OR. &
             SIZE(NR) < ND .OR. SIZE(NRV) < ND .OR. &
             SIZE(COMP) < ND) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: arrays too small for ND=', ND
            STOP
         END IF

         IF (dbg .AND. my_rank == 0) WRITE (*, *) ' ND check', ND, NSS*NRR

         DO J = 1, ND
            READ (unit_acq, '(A)', IOSTAT=ios) line
            IF (ios /= 0) THEN
               IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: line read failed at row', J, ' IOSTAT=', ios
               STOP
            END IF

            ! try 5-int form: "ND Shot ICS Rec ICR"
            READ (line, *, IOSTAT=ios) nd_file, shotid, ics, recid, icr
            IF (ios == 0) THEN
               NS(J) = shotid
               NSV(J) = ics
               NR(J) = recid
               NRV(J) = icr
            ELSE
               ! fallback 4-int form: "Shot ICS Rec ICR"
               READ (line, *, IOSTAT=ios) NS(J), NSV(J), NR(J), NRV(J)
               IF (ios /= 0) THEN
                  IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: bad row ', J, ':', TRIM(line)
                  STOP
               END IF
            END IF

            IF (dbg .AND. my_rank == 0) WRITE (*, *) NS(J), NSV(J), NR(J), NRV(J)

            ! 2.5D keep/zero rule
            IF ((NSV(J) == 2 .AND. NRV(J) == 1) .OR. &
                (NSV(J) == 1 .AND. NRV(J) == 2) .OR. &
                (NSV(J) == 2 .AND. NRV(J) == 3) .OR. &
                (NSV(J) == 3 .AND. NRV(J) == 2)) THEN
               COMP(J) = 0
            ELSE
               COMP(J) = 1
            END IF
         END DO

         CLOSE (unit_acq)

         IF (MAXVAL(NS(1:ND)) > NSS .OR. MINVAL(NS(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NS out of range [1,NSS]'
            STOP 'Read_acquisition: NS index out of range'
         END IF

         IF (MAXVAL(NR(1:ND)) > NSS + NRR .OR. MINVAL(NR(1:ND)) < NSS + 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NR out of range [NSS+1,NSS+NRR]'
            STOP 'Read_acquisition: NR index out of range'
         END IF

         IF (MAXVAL(NSV(1:ND)) > 3 .OR. MINVAL(NSV(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NSV outside [1,3]'
            STOP 'Read_acquisition: NSV index out of range'
         END IF

         IF (MAXVAL(NRV(1:ND)) > 3 .OR. MINVAL(NRV(1:ND)) < 1) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NRV outside [1,3]'
            STOP 'Read_acquisition: NRV index out of range'
         END IF

         IF (MAXVAL(COMP(1:ND)) > 1 .OR. MINVAL(COMP(1:ND)) < 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: COMP not in {0,1}'
            STOP 'Read_acquisition: invalid COMP flag'
         END IF

         IF (NSS + NRR > SIZE(VSR, 1)) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NSS+NRR=', NSS + NRR, &
               ' > SIZE(VSR,1)=', SIZE(VSR, 1)
            STOP 'Read_acquisition: VSR first dimension too small'
         END IF

         IF (MAXVAL(NS(1:ND)) > SIZE(VSR, 1) .OR. MAXVAL(NR(1:ND)) > SIZE(VSR, 1)) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: NS or NR exceeds SIZE(VSR,1)'
            STOP 'Read_acquisition: NS/NR vs VSR size mismatch'
         END IF

         IF (SIZE(VSR, 2) < 3 .OR. SIZE(VSR, 3) < 3) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_acquisition: VSR second/third dims < 3'
            STOP 'Read_acquisition: VSR must be (N_s+N_r,3,3)'
         END IF
      END IF

      ! same ND for all frequencies
      NDATA(1:NFQ) = ND

      ! -------- optional debug export --------
      IF (dbg .AND. my_rank == 0) THEN
         OPEN (UNIT=61, FILE='data_mapping.txt', STATUS='REPLACE', ACTION='WRITE')
         WRITE (61, '(A,I0)') 'ND = ', ND
         WRITE (61, '(A)') ' J   NS  NSV  NR  NRV   COMP   Sx Sy Sz   Gx Gy Gz'

         DO J = 1, ND
            Svec(1:3) = VSR(NS(J), NSV(J), 1:3)
            Gvec(1:3) = VSR(NR(J), NRV(J), 1:3)

            WRITE (61, '(I4,1X, I3,1X,I1,2X, I3,1X,I1,2X, I3,2X, 3(I1,1X),2X, 3(I1,1X))') &
               J, NS(J), NSV(J), NR(J), NRV(J), COMP(J), &
               NINT(Svec(1)), NINT(Svec(2)), NINT(Svec(3)), &
               NINT(Gvec(1)), NINT(Gvec(2)), NINT(Gvec(3))
         END DO

         CLOSE (61)
      END IF

   END SUBROUTINE Read_acquisition

   SUBROUTINE Read_seismic_spec(FREQ_FILE, IFQ, FREQ, ND, GT0, my_rank, DEBUG_OUTPUT)
      USE constant_mod, ONLY: PI   ! <<-- use shared constant
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ND, my_rank, IFQ
      REAL(dp), INTENT(IN)  :: FREQ
      COMPLEX(dp), INTENT(OUT) :: GT0(:)
      CHARACTER(len=*), INTENT(IN)  :: FREQ_FILE(:)
      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT

      ! ==== Locals ====
      LOGICAL            :: dbg, first_is_header
      CHARACTER(len=256) :: line
      REAL(dp)           :: fcol, re, im
      INTEGER            :: idx, i, ios, unit
      COMPLEX(dp)        :: jomega, vval

      dbg = .FALSE.
      IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

      !------------------ Open file safely ------------------
      OPEN (NEWUNIT=unit, FILE=TRIM(FREQ_FILE(IFQ)), STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'Read_seismic_spec: cannot open file ', TRIM(FREQ_FILE(IFQ)), ' IOSTAT=', ios
         STOP
      ELSEIF (dbg .AND. my_rank == 0) THEN
         WRITE (*, *) 'Opened seismic file ', TRIM(FREQ_FILE(IFQ))
      END IF

      !------------------ Frequency guard ------------------
      IF (ABS(FREQ) <= 1.0D-12) THEN
         IF (my_rank == 0) WRITE (*, *) 'Read_seismic_spec: FREQ is zero or too small.'
         STOP
      END IF
      jomega = CMPLX(0.0_dp, 2.0_dp*PI*FREQ, KIND=dp)

      !------------------ Detect header line ------------------
      READ (unit, '(A)', IOSTAT=ios) line
      IF (ios /= 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'Read_seismic_spec: cannot read first line; IOSTAT=', ios
         STOP
      END IF

      READ (line, *, IOSTAT=ios) fcol, idx, re, im
      IF (ios == 0) THEN
         first_is_header = .FALSE.
         BACKSPACE (unit)
      ELSE
         first_is_header = .TRUE.
         IF (dbg .AND. my_rank == 0) WRITE (*, *) 'Header detected in file ', TRIM(FREQ_FILE(IFQ))
      END IF

      !------------------ Read ND rows ------------------
      DO i = 1, ND
         READ (unit, *, IOSTAT=ios) fcol, idx, re, im
         IF (ios /= 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'Read_seismic_spec: premature EOF at row ', i
            STOP
         END IF
         vval = CMPLX(re, im, KIND=dp)
         GT0(i) = vval/jomega  ! convert velocity -> displacement

         IF (dbg .AND. my_rank == 0) WRITE (*, '(A,I5,2(A,F10.5))') 'Row=', i, ' re=', re, ' im=', im
      END DO

      CLOSE (unit)
      IF (my_rank == 0) WRITE (*, *) 'Data imported for frequency', FREQ

   CONTAINS

      SUBROUTINE strip_inplace(s)
         CHARACTER(len=*), INTENT(INOUT) :: s
         INTEGER :: i, j, n, c
         n = LEN_TRIM(s); j = 1
         DO i = 1, n
            c = IACHAR(s(i:i))
            IF (c == 9 .OR. c == 13 .OR. c == 160 .OR. c < 32 .OR. c == 0) CYCLE  ! TAB, CR, NBSP, control, NUL
            s(j:j) = s(i:i); j = j + 1
         END DO
         IF (j <= LEN(s)) s(j:) = ' '
      END SUBROUTINE strip_inplace

      SUBROUTINE dump_bytes(label, s)
         CHARACTER(len=*), INTENT(IN) :: label
         CHARACTER(len=*), INTENT(IN) :: s
         INTEGER :: i, n
         n = LEN_TRIM(s)
         WRITE (*, '(A,I0)') TRIM(label)//' len_trim=', n
         WRITE (*, '(A)', ADVANCE='NO') '  bytes:'
         DO i = 1, n
            WRITE (*, '(1X,I3)', ADVANCE='NO') IACHAR(s(i:i))
         END DO
         WRITE (*, *)
      END SUBROUTINE dump_bytes

   END SUBROUTINE Read_seismic_spec

   SUBROUTINE CheckQuadratureSampling(CMIN, CMAX, FREQ1, FREQ2, DX, DZ, NORD, WLMAX, WLMIN, my_rank)

   IMPLICIT NONE

   REAL(dp), INTENT(IN)  :: CMIN, CMAX, FREQ1, FREQ2, DX, DZ
   INTEGER,  INTENT(IN)  :: NORD, my_rank
   REAL(dp), INTENT(OUT) :: WLMAX, WLMIN

   REAL(dp) :: epw_x, epw_z, epw_min
   REAL(dp) :: qpw_x, qpw_z, qpw_min
   REAL(dp) :: dx_req_epw, dz_req_epw
   REAL(dp) :: dx_req_qpw, dz_req_qpw
   REAL(dp) :: dx_req, dz_req

   REAL(dp), PARAMETER :: REQ_EPW = 2.0_dp   ! minimum elements/subdomains per wavelength
   REAL(dp), PARAMETER :: REQ_QPW = 5.0_dp   ! minimum GQG points per wavelength

   ! --- Band-edge wavelengths ---
   WLMAX = CMAX / FREQ1
   WLMIN = CMIN / FREQ2

   ! --- Elements/subdomains per minimum wavelength ---
   epw_x = WLMIN / DX
   epw_z = WLMIN / DZ
   epw_min = MIN(epw_x, epw_z)

   ! --- GQG points per minimum wavelength ---
   ! Criterion uses the number of quadrature points per element, not nodal intervals.
   qpw_x = REAL(NORD, dp) * epw_x
   qpw_z = REAL(NORD, dp) * epw_z
   qpw_min = MIN(qpw_x, qpw_z)

   IF (my_rank == 0) THEN
      WRITE (*,'(A,F10.4)')       'lambda_min:       ', WLMIN
      WRITE (*,'(A,F10.4)')       'lambda_max:       ', WLMAX
      WRITE (*,'(A,F8.2,2X,A,F8.2)') 'Element spacing:  DX=', DX, 'DZ=', DZ
      WRITE (*,'(A,I2)')          'NORD (GQG points/elem)=', NORD

      WRITE (*,'(A,F8.3,2X,A,F8.3)') 'Elements/wavelength:  X=', epw_x, 'Z=', epw_z
      WRITE (*,'(A,F8.3,2X,A,F8.3)') 'GQG points/wavelength: X=', qpw_x, 'Z=', qpw_z

      WRITE (*,'(A,F4.1)')        'Required elements/wavelength = ', REQ_EPW
      WRITE (*,'(A,F4.1)')        'Required GQG points/wavelength = ', REQ_QPW
   END IF

   IF (epw_min < REQ_EPW .OR. qpw_min < REQ_QPW) THEN

      ! --- Maximum admissible element size from each criterion ---
      dx_req_epw = WLMIN / REQ_EPW
      dz_req_epw = WLMIN / REQ_EPW

      dx_req_qpw = REAL(NORD, dp) * WLMIN / REQ_QPW
      dz_req_qpw = REAL(NORD, dp) * WLMIN / REQ_QPW

      ! --- Need to satisfy both criteria, so take the smaller allowed size ---
      dx_req = MIN(dx_req_epw, dx_req_qpw)
      dz_req = MIN(dz_req_epw, dz_req_qpw)

      IF (my_rank == 0) THEN
         WRITE (*,*) 'NOT ENOUGH SAMPLING.'
         WRITE (*,'(A,F8.3)') '  Minimum elements/wavelength found: ', epw_min
         WRITE (*,'(A,F8.3)') '  Minimum GQG points/wavelength found: ', qpw_min

         WRITE (*,'(A,F10.4)') '  Maximum DX/DZ allowed by element criterion ≈ ', dx_req_epw
         WRITE (*,'(A,F10.4)') '  Maximum DX/DZ allowed by GQG-point criterion ≈ ', dx_req_qpw
         WRITE (*,'(A,F10.4)') '  Required maximum DX/DZ satisfying both ≈ ', dx_req

         WRITE (*,'(A,I2,A)') '  For NORD=', NORD, ', suggested element DX/DZ:'
         WRITE (*,'(A,F8.1,2X,A,F8.1)') '   - snapped to 5 m :  DX=', &
            floor_to_step(dx_req, 5.0_dp), 'DZ=', floor_to_step(dz_req, 5.0_dp)
         WRITE (*,'(A,F8.1,2X,A,F8.1)') '   - snapped to 10 m:  DX=', &
            floor_to_step(dx_req, 10.0_dp), 'DZ=', floor_to_step(dz_req, 10.0_dp)
      END IF

      STOP

   ELSE

      IF (my_rank == 0) WRITE (*,*) 'Sampling OK.'

   END IF

CONTAINS

   REAL(dp) FUNCTION floor_to_step(val, step)
      REAL(dp), INTENT(IN) :: val, step
      floor_to_step = MAX(step, FLOOR(val / step) * step)
   END FUNCTION floor_to_step

END SUBROUTINE CheckQuadratureSampling

   SUBROUTINE ReadWaveletFile(wavelet, NFQ)
      !----------------------------------------------------------------------
      !  ReadWaveletFile reads a complex-valued wavelet from a text file.
      !
      !  Input:
      !    NFQ           - Number of frequency points
      !
      !  Output:
      !    wavelet(NFQ)  - COMPLEX(dp) array containing wavelet spectrum
      !
      !----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NFQ
      COMPLEX(dp), INTENT(OUT) :: wavelet(NFQ)

      INTEGER :: I, istat
      real(dp) :: fq, reval, imagval

      OPEN (6, FILE='../wavelet.txt', STATUS='OLD', ACTION='READ', IOSTAT=istat)

      IF (istat /= 0) THEN
         WRITE (*, *) 'Error opening wavelet.txt, IOSTAT = ', istat
         STOP
      END IF

      DO I = 1, NFQ
         READ (6, *) fq, reval, imagval
         wavelet(I) = DCMPLX(reval, imagval)
      END DO

      CLOSE (6)
   END SUBROUTINE ReadWaveletFile
   SUBROUTINE safe_deallocate_vec(arr, label)
      real(dp), ALLOCATABLE, INTENT(INOUT) :: arr(:)
      CHARACTER(*), INTENT(IN) :: label
      INTEGER :: ierr
      IF (ALLOCATED(arr)) THEN
         DEALLOCATE (arr, STAT=ierr)
         IF (ierr /= 0) WRITE (*, '(A,": DEALLOCATE failed (STAT=",I0,")")') TRIM(label), ierr
      END IF
   END SUBROUTINE safe_deallocate_vec

   SUBROUTINE safe_deallocate_mat(arr, label)
      real(dp), ALLOCATABLE, INTENT(INOUT) :: arr(:, :)
      CHARACTER(*), INTENT(IN) :: label
      INTEGER :: ierr
      IF (ALLOCATED(arr)) THEN
         DEALLOCATE (arr, STAT=ierr)
         IF (ierr /= 0) WRITE (*, '(A,": DEALLOCATE failed (STAT=",I0,")")') TRIM(label), ierr
      END IF
   END SUBROUTINE safe_deallocate_mat

   SUBROUTINE safe_deallocate_cmat(arr, label)
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT) :: arr(:, :)
      CHARACTER(*), INTENT(IN) :: label
      INTEGER :: ierr
      IF (ALLOCATED(arr)) THEN
         DEALLOCATE (arr, STAT=ierr)
         IF (ierr /= 0) WRITE (*, '(A,": DEALLOCATE failed (STAT=",I0,")")') TRIM(label), ierr
      END IF
   END SUBROUTINE safe_deallocate_cmat
END MODULE input_mod
