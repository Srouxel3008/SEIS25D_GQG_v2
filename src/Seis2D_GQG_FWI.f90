program Seis2D_GQG_FWI

   !-----------------------------------------------------------------------
   USE OMP_LIB
   USE MPI

   USE shared_mod
   USE err_mpi_mod
   USE input_mod
   USE gridtype_mod
   USE output_mod
   USE hardware_mod     ! max memory and thread
   USE ky_sampling_mod         ! Y-dimension wavenumber functions
   USE grid_mod         ! GQG grid related functions
   USE boundaries_mod   ! absorbing boundaries
   USE stiffness_assembly_mod ! Stiffness matrix assembly routines
   USE solvers_mod   ! choice of solver
   USE F_modeling       ! forward modeling returning Green's functions
   USE Frechet_mod  ! Frechét  derivative computation
   USE gradient_mod     ! gradient computation
   USE lbfgs_mod        ! LBFGS-related subroutines
   USE lbfgs_lsearch
   USE constant_mod
   USE scalers_mod
   USE misfit_mod

   use iso_fortran_env, only: dp => real64, sp => real32

   IMPLICIT NONE

   ! ======= constants / simple chars =======

   CHARACTER(len=80)  :: LAB, msg
   CHARACTER(len=15)  :: FNAME
   CHARACTER(len=30)  :: FNAME20
   CHARACTER(len=1)   :: T
   CHARACTER(len=2)   :: IT
   CHARACTER(len=150) :: MainInput, InitInput, InitInputQ, TrueInput, TrueInputQ

! ======= timers / mpi =======
   INTEGER :: iTimes1, iTimes2, rate, clock_max, istat
   INTEGER :: ierr, my_rank, n_process, comm
   INTEGER :: iter_clock_start, iter_clock_end, iter_tick_delta
   REAL(dp) :: t_global_begin, t_global_end
   LOGICAL ::  converged, exit_iteration_loop, recompute_model_scaling
   LOGICAL :: DEBUG_OUTPUT, DEBUG_INPUT, DEBUG_DBG, DEBUG_QC
   LOGICAL :: active_q_stage, active_nonq_stage, active_q_only_stage
! ======= job / inversion controls =======
   INTEGER :: INV, I25D, NFBAND, NFQ, MAXITER, ITER, IFQ, total_iter_count
   INTEGER :: REG, GZSMTH, NORD, SOLVER_KIND, PRECOND
   REAL(dp) :: GZ1, GZ2
   REAL(dp) :: DX, FREQ, FREQ1, FREQ2, WLMAX, WLMIN, DZ_USE, DZ, DX_use
   REAL(dp) :: iter_elapsed, GRAD_scaled_norm_conv, grad0_norm
   REAL(dp) :: factr, cost_conv_tol, rel_cost_drop, rel_cost_denom
   REAL(dp) :: model_conv_tol, m_coarse0_norm, norm_diff
   LOGICAL  :: have_true, have_true_q, have_init_q, stage_has_new_param, EXTERNAL_DATA

! NEW: per-frequency observation files from input
   CHARACTER(len=256), ALLOCATABLE :: FREQ_FILE(:)
   INTEGER, ALLOCATABLE :: NFQ_PER_BAND(:)
   INTEGER, ALLOCATABLE :: band_of_ifq(:)    ! size NFQ, maps IFQ -> band index (1..NFBAND)
   INTEGER, ALLOCATABLE :: band_start(:), band_end(:)
   REAL(dp), ALLOCATABLE :: band_fmax(:), band_fmin(:)
   INTEGER :: ib, last_band, f1b, f2b, sumq, NORD_USE
   INTEGER :: j
! ======= acquisition =======
   INTEGER :: NSS, NRR, NSR, NSSt, IZ, ISR90
   INTEGER :: IDD, ND, IS00, NCOMPS, NCOMPR, NCOMP
   INTEGER, ALLOCATABLE :: MSR(:), MSR1(:, :), COMP(:)
   INTEGER, ALLOCATABLE :: NS(:), NR(:), NSV(:), NRV(:), ICSR(:), NDATA(:)
   REAL(dp), ALLOCATABLE :: XSR(:), YSR(:), ZSR(:), VSR(:, :, :)
   REAL(dp), ALLOCATABLE :: FSR(:, :)
   CHARACTER(len=256):: ACQ_WEIGHT_FILE

! ======= grid =======
   INTEGER, ALLOCATABLE :: MZ(:)
   REAL(dp), ALLOCATABLE :: XTO(:), ZTO(:), FREQN(:)
   REAL(dp) :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
   REAL(dp), ALLOCATABLE :: X(:), Z(:), XM(:), ZM(:), XP(:), ZP(:)
   TYPE(InversionGridType) :: IG
   LOGICAL, ALLOCATABLE :: vPML(:)
   INTEGER :: MX, IS0, NTO, IE0, MZT
   INTEGER :: NX, NZ, NNX, NNZ, NN, NPT, NBLOCK, NM, NMM
   TYPE(PhysicalState) :: canon
   REAL(dp) :: FREQ_MESH, DX_REF, DZ_REF, CREF
   REAL(dp) :: src_burial_max, zsurf_src
   LOGICAL  :: REUSE_MESH_FROM_HIGHEST, MESH_BUILT
   ! REAL(dp), PARAMETER :: WIDTH = 1.0_dpREAL(dp) :: WIDTH
   REAL(dp) :: WIDTH
   CHARACTER(len=32) :: WIDTH_ENV, WIDTH_ENV_TRIM
   INTEGER :: WIDTH_ENV_STAT
! OPTIONAL per-band meshing (filled if present in input)
   INTEGER, ALLOCATABLE :: NORD_BAND(:)
   REAL(dp), ALLOCATABLE :: DX_BAND(:), DZ_BAND(:)

! ======= material / parameters =======
   INTEGER :: IVISCO, IANISO, ITHOM, NPAR, NPAR_INV, NGRP_INV, NGRP_INV_BAND
   INTEGER, ALLOCATABLE :: INVP_ORD(:), INVP(:)
   INTEGER, ALLOCATABLE :: INVP_ORD_BAND(:, :)
   CHARACTER(len=3) :: MODEL(22)
   REAL(dp), ALLOCATABLE :: CRT(:, :), CIT(:, :), CR(:, :), CI(:, :)
   REAL(dp), ALLOCATABLE :: CR0(:, :), CI0(:, :), CRR0(:, :), CII0(:, :)
   REAL(dp), ALLOCATABLE :: CRREG(:, :), CIREG(:, :)
   REAL(dp), ALLOCATABLE :: CRRU(:, :), CIIU(:, :)

! ======= forward modeling =======
   INTEGER :: NK
   REAL(dp), ALLOCATABLE :: FKY(:), WTK(:), AS(:), WT(:), SIG(:), FKYI(:), WTKY(:)
   REAL(dp), ALLOCATABLE :: AMPT0(:), AMP0(:)
   COMPLEX(sp), ALLOCATABLE :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
   COMPLEX(sp), ALLOCATABLE :: GFXT(:, :, :, :), GFYT(:, :, :, :), GFZT(:, :, :, :)
   COMPLEX(dp), ALLOCATABLE ::G0(:), GT0(:), WAVELET(:)
   COMPLEX(dp), ALLOCATABLE :: S(:), G(:), SourceScaler(:)

! ======= frechet / gradient =======
   COMPLEX(dp), ALLOCATABLE :: FRECHET(:, :)
   REAL(dp), ALLOCATABLE :: GRADr(:), GRAD_scaled(:),GRAD_smooth(:), GRADsc_precond(:), GRADr_NORM(:), GRAD_scaled_NORM(:), GRADsc_precond_NORM(:)

   REAL(dp), ALLOCATABLE :: GR1(:), GR2(:), PRE(:), PREE(:), HESS(:), HESSDI_scaled(:), HESSDI_scaled_norm(:)
   REAL(dp), ALLOCATABLE :: GMASK(:), GTAPPER(:), GMASK_FINAL(:)
   INTEGER :: NZ_ACT
   LOGICAL :: USE_GR_SMOOTH
   REAL(dp) :: sigma_grad_x, sigma_grad_z

! ======= LBFGS =======
   INTEGER :: IFLAG, USE_LBFGS_TYPE, MML
   CHARACTER(len=5) :: PARAM(22)
   REAL(dp) :: g2, denom, eps_h, beta, s_norm
   REAL(dp), ALLOCATABLE :: p_k(:), p_kf(:), grad_prev(:), BF_grad_res(:, :), BF_s_hist(:, :)
   REAL(dp), ALLOCATABLE :: m_coarse_prev(:), m_fine(:), m_coarse(:)
   REAL(dp), ALLOCATABLE :: BF_a(:), BF_p(:), BF_q(:), BF_r(:)

   LOGICAL :: USE_PRECOND, USE_RMS, USE_GROUP

   !====SCALERS============
   REAL(dp), ALLOCATABLE :: PAR_SCALE(:), PAR_SCALE_REG(:), SCALER(:)
   REAL(dp), ALLOCATABLE :: F_SCALE_L2(:), F_SCALE(:)
   REAL(dp), ALLOCATABLE :: BALANCE_SCALE(:)
   REAL(dp) :: BALANCE_PCTL, BALANCE_WMIN, BALANCE_WMAX
   INTEGER :: balance_ref_param
   ! NEW: regularization lambda from input + LBFGS controls from input
   REAL(dp) :: LAMBDA
   INTEGER  :: LBFGS_TYPE
   LOGICAL :: balance_initialized

! ======= cost =======
   REAL(dp) :: FCOST0, FCOST0_amp, FCOSTMIN, lmbd, FCOST0_REGm, W_Q, FCOST_prev
   REAL(dp), ALLOCATABLE :: RESID_REAL(:), RESID_amp(:), m_coarse_reg(:), WD_acq(:)
   COMPLEX(dp), ALLOCATABLE :: RESID(:), RESID_m(:)
   REAL(dp) :: REG_LAMBDA
   real(dp), ALLOCATABLE :: diff_model(:)
   real(dp) :: CMIN, CMAX, RESID_L2, RESID_RMS, WD_amp, WQS
   REAL(dp) :: DZ0, Z0, ZMIN0
   INTEGER, PARAMETER:: MIN_ITER_CHECK = 3

   !-------------Loops-----------------------------------
   INTEGER :: ISO, NSIG, MAXDA, III
   INTEGER :: I, IP, IM, II, IA, IK, cs, ce, ps, pe, i1, i2, col, PA
   real(dp) :: ALPHA00_prev, ALPHA00, ALPHA, ALPHA_prev, pnorm, dot_g_p      ! loop/index
   !----------------------------- MPI Init -----------------------------
   real(dp) ::FCOST0_ref, dot_ref
   real(dp), allocatable :: p_k_ref(:), GRAD_scaled_ref(:)
! ======= MPI init =======
   CALL MPI_Init(ierr)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   CALL MPI_Comm_size(MPI_COMM_WORLD, n_process, ierr)
   CALL err_init(MPI_COMM_WORLD, my_rank)
   comm = MPI_COMM_WORLD
   ! CALL initialize_hardware()
   balance_initialized = .FALSE.
   ALPHA00_prev = 1.0_dp
   IF (my_rank == 0) WRITE (*, *) 'n_process', n_process
   if (my_rank == 0) WRITE (*, '(A)') ' Build mode   : '//BUILD_MODE

   IF (my_rank == 0) THEN
      WRITE (*, *) '                                                      '
      WRITE (*, *) '      ************************************************'
      WRITE (*, *) '      *                                              *'
      WRITE (*, *) '      * 2D/2.5D frequency-domain seismic wave model- *'
      WRITE (*, *) '      * ing  in  arbitrary  anisotropic viscoelastic *'
      WRITE (*, *) '      * media with GQG approach and multiparameter FWI *'
      WRITE (*, *) '      ************************************************'
   END IF
   WIDTH = 1.2_dp
   BALANCE_PCTL = 0.95_dp
   BALANCE_WMIN = 0.005_dp
   BALANCE_WMAX = 1.0e1_dp
! ======= Input / Output =======
   CALL SYSTEM_CLOCK(count_rate=rate, count_max=clock_max)
   t_global_begin = MPI_Wtime()

   CALL getarg(1, MainInput)
   CALL getarg(2, InitInput)
   CALL getarg(3, InitInputQ)
   CALL getarg(4, TrueInput)
   CALL getarg(5, TrueInputQ)
   IF (LEN_TRIM(MainInput) == 0) THEN
      PRINT *, 'Error: No input file provided as argument.'
      STOP
   END IF

   IF (my_rank == 0) THEN
      WRITE (*, *) 'MainInput file', TRIM(MainInput)
      WRITE (*, *) 'InitInput file', TRIM(InitInput)

   END IF

   OPEN (1, FILE=TRIM(MainInput), STATUS='OLD', ACTION='READ')
   total_iter_count = 0
   DEBUG_INPUT = .FALSE. !Turn on write off for checking input read
   DEBUG_DBG = .FALSE.
   DEBUG_QC = .TRUE.

! ======= (1) Job definition t =======
   IF (my_rank == 0) WRITE (*, *) '---- (1) working-job definition ----------------'

   CALL Read_Job_Setup( &
      INV, I25D, NFBAND, NFQ, MAXITER, REG, FREQN, FREQ1, FREQ2, &
      IANISO, ITHOM, IVISCO, CMIN, CMAX, NPAR, NPAR_INV, INVP_ORD, INVP_ORD_BAND, &
      NORD, DX, DZ, GZSMTH, GZ1, GZ2, PRECOND, SOLVER_KIND, MODEL, &
      FREQ_FILE, LAMBDA, BETA, EXTERNAL_DATA, &
      LBFGS_TYPE, MML, USE_GR_SMOOTH, sigma_grad_x, sigma_grad_z, &
      NFQ_PER_BAND, NORD_BAND, DX_BAND, DZ_BAND, &
      cost_conv_tol, model_conv_tol, my_rank, DEBUG_OUTPUT=DEBUG_INPUT)
   IF (NPAR > 0) THEN
      NGRP_INV = MAXVAL(INVP_ORD(1:NPAR))
   ELSE
      NGRP_INV = 0
   END IF

   REG_LAMBDA = 0.0_dp
   ! ======= (2) Model loading =======

   OPEN (2, FILE=TRIM(InitInput), STATUS='OLD', ACTION='READ')

   have_init_q = (LEN_TRIM(InitInputQ) > 0)
   if (my_rank == 0) WRITE (*, *) 'InitInputQ file', TRIM(InitInputQ)
   IF (have_init_q) OPEN (3, FILE=TRIM(InitInputQ), STATUS='OLD', ACTION='READ')

   have_true = (LEN_TRIM(TrueInput) > 0 .AND. (INV == 1 .or. INV == 0 .or. INV == 3))
   have_true_q = (LEN_TRIM(TrueInputQ) > 0 .AND. (INV == 1 .or. INV == 0. .or. INV == 3))
   IF (have_true) OPEN (4, FILE=TRIM(TrueInput), STATUS='OLD', ACTION='READ')
   IF (have_true_q) OPEN (5, FILE=TRIM(TrueInputQ), STATUS='OLD', ACTION='READ')

   !---- (3)topography curves ---------------------

   CALL ReadTopographyData(NTO, IS0, XTO, ZTO, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, my_rank, DEBUG_OUTPUT=DEBUG_INPUT)

!---- (2) sources & geophones'---------------------
   CALL input_SR(NSR, NSS, NRR, XSR, YSR, ZSR, ICSR, VSR, &
                 ISR90, NCOMPS, NCOMPR, NCOMP, &
                 XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, &
                 my_rank, DEBUG_OUTPUT=DEBUG_INPUT)

   IF (my_rank == 0) WRITE (*, *) '  '
   IF (my_rank == 0 .and. NSR == 2 .or. INV == 3) WRITE (*, *) 'Frechet derivatives/Gradient computations only'
   IF (my_rank == 0 .and. INV == 0) WRITE (*, *) 'Forward modeling only'

   !------grid geometry-----------------------------

   CALL Read_Grid_Geometry(MX, MZ, MZT, XM, ZM, XMIN, XMAX, ZMIN, ZMAX, my_rank, mpi_comm_world, DEBUG_OUTPUT=DEBUG_INPUT)

   !---- anisotropic model ---------------------

   CALL Read_Model_Parameters(MX, MZ, have_true, have_true_q, IANISO, ITHOM, IVISCO, NSR, INVP_ORD, MODEL, &
                              CR0, CRT, CI0, CIT, my_rank, mpi_comm_world, DEBUG_OUTPUT=DEBUG_INPUT)

   if (my_rank == 0) WRITE (*, *) 'MODEL PARAMETERS LOADED'

   !---- data acquisition ---------------------------

   IF (my_rank == 0) WRITE (*, *) '      (5) data acquisition mode'

   ! Allocate arrays for storing data acquisition parameters

   ALLOCATE (NDATA(NFQ))
   MAXDA = MAX(NSS*NRR*NCOMPS*NCOMPR, 1)
   ALLOCATE (NS(MAXDA), NSV(MAXDA), NR(MAXDA), NRV(MAXDA), COMP(MAXDA))

   CALL Read_acquisition(NFQ, NSS, NRR, ICSR, VSR, ACQ_WEIGHT_FILE, &
                         NDATA, NS, NSV, NR, NRV, NCOMPS, NCOMPR, COMP, &
                         my_rank, DEBUG_OUTPUT=DEBUG_INPUT)

   ALLOCATE (wavelet(NFQ))
   ! Default to a unit source spectrum unless an explicit wavelet file is read.
   wavelet = CMPLX(1.0_dp, 0.0_dp, dp)
   ! wavelet = CMPLX(100000000.0_dp, 0.0_dp, dp)
   if (my_rank == 0) WRITE (*, *) 'beta for hess', beta
   ! CALL ReadWaveletFile(wavelet, NFQ)
   CLOSE (5)
   CLOSE (4)
   CLOSE (3)
   CLOSE (2)

   CLOSE (1)
   IF (my_rank == 0) THEN !create debug and output files
      OPEN (UNIT=70, FILE='Seis2D_GQG.txt', STATUS='UNKNOWN', POSITION='APPEND')
      OPEN (UNIT=64, FILE='out_resid.txt', STATUS='UNKNOWN', POSITION='APPEND')
      OPEN (UNIT=66, FILE='out_diag.txt', STATUS='UNKNOWN', POSITION='APPEND')
      ! OPEN (UNIT=67, FILE='out_descent.txt', STATUS='UNKNOWN', POSITION='APPEND')
      ! OPEN (UNIT=69, FILE='out_descentC.txt', STATUS='UNKNOWN', POSITION='APPEND')
   END IF

   CALL SET_PARAMETERS(IANISO, IVISCO, ITHOM, PARAM, III)

   IF (my_rank == 0) WRITE (*, *) 'INITIAL INPUT LOADING COMPLETE'

   CALL Freeze_Phys_Grid(XMIN, XMAX, ZMIN, ZMAX, &
                         NTO, XTO, ZTO, XSR, YSR, ZSR, NSR, &
                         XM, ZM, MX, MZ, MZT, canon, msg, my_rank, DEBUG_OUTPUT=DEBUG_DBG)
   CREF = 0.5_dp*(CMIN + CMAX)    ! fix: average of [CMIN,CMAX]

   ALLOCATE (PAR_SCALE(NPAR), F_SCALE_L2(NPAR), F_SCALE(NPAR), SCALER(NPAR), BALANCE_SCALE(NPAR), PAR_SCALE_REG(NPAR))
   ALLOCATE (GRADr_NORM(NPAR), GRAD_scaled_NORM(NPAR), GRADsc_precond_NORM(NPAR))
   ALLOCATE (INVP(NPAR + 1))
   INVP = 0
   SCALER = 1
   CALL ComputeModelScales(CR0, CI0, NPAR, IANISO, ITHOM, PAR_SCALE, my_rank, &
                           USE_RMS=.FALSE., DEBUG_OUTPUT=DEBUG_QC)

   USE_PRECOND = (PRECOND /= 0)
   PAR_SCALE = 1.0_dp

! cumulative offsets per band
   ALLOCATE (band_of_ifq(NFQ), band_start(NFBAND), band_end(NFBAND), band_fmax(NFBAND), band_fmin(NFBAND))

   sumq = 0
   DO ib = 1, NFBAND
      band_start(ib) = sumq + 1
      band_end(ib) = sumq + NFQ_PER_BAND(ib)
      sumq = sumq + NFQ_PER_BAND(ib)
   END DO

   DO ib = 1, NFBAND
      DO j = band_start(ib), band_end(ib)
         band_of_ifq(j) = ib
      END DO
      band_fmax(ib) = MAXVAL(FREQN(band_start(ib):band_end(ib)))
      band_fmin(ib) = MINVAL(FREQN(band_start(ib):band_end(ib)))
   END DO

   REUSE_MESH_FROM_HIGHEST = .NOT. (ALLOCATED(NORD_BAND) .OR. ALLOCATED(DX_BAND) .OR. ALLOCATED(DZ_BAND))
   MESH_BUILT = .FALSE.

! BAND loop (mesh per band), then FREQ loop
!========================
   band_loop: DO ib = 1, NFBAND

      IF (ALLOCATED(INVP_ORD_BAND)) THEN
         NGRP_INV_BAND = MAXVAL(INVP_ORD_BAND(ib, 1:NPAR))
      ELSE
         NGRP_INV_BAND = NGRP_INV
      END IF

      IF (ALLOCATED(NORD_BAND)) THEN
         NORD_use = NORD_BAND(ib)
      ELSE
         NORD_use = NORD
      END IF
      IF (ALLOCATED(DX_BAND)) THEN
         DX_use = DX_BAND(ib)
      ELSE
         DX_use = DX
      END IF
      IF (ALLOCATED(DZ_BAND)) THEN
         DZ_use = DZ_BAND(ib)
      ELSE
         DZ_use = DZ
      END IF

      f1b = band_start(ib)
      f2b = band_end(ib)
      src_burial_max = 0.0_dp
      DO j = 1, NSS
         zsurf_src = ZH(canon%NTOC, canon%XTOC, canon%ZTOC, canon%XSRC(j))
         src_burial_max = MAX(src_burial_max, MAX(0.0_dp, zsurf_src - canon%ZSRC(j)))
      END DO
      ! Keep at least one excluded cell around the source region so surface sources
      ! do not turn off the near-source mask entirely.
      IZ = MAX(1, CEILING(src_burial_max/DZ_use))
      IF (REUSE_MESH_FROM_HIGHEST) THEN
         IF (.NOT. MESH_BUILT) THEN
            ! No band-meshing block: build one global mesh once and reuse it.
            FREQ_MESH = MAXVAL(FREQN)
            WLMAX = CMAX/MINVAL(FREQN)
            IF (my_rank == 0) THEN
               WRITE (*, *) '--- Building single global mesh for all bands', &
                  ' (NORD=', NORD_use, ', DX=', DX_use, ', DZ=', DZ_use, ')  FREQ_MESH=', FREQ_MESH
               WRITE (*, *) 'Max source burial for GQG mask (m) = ', src_burial_max
               WRITE (*, *) 'IZ for GQG mesh = ', IZ
            END IF
            CALL BuildMeshForBand_FromCanon( &
               FREQ_MESH, ib, WLMAX, WIDTH, IS0, NORD_use, DX_use, DZ_use, &
               IANISO, ITHOM, IVISCO, PARAM, NPAR, &
               CMIN, CMAX, FREQ_MESH, NSIG, ISO, &
               NSR, NSS, INV, my_rank, NTO, XTO, ZTO, &
               canon, CR0, CI0, &
               IE0, NX, NZ, NZ_ACT, NNX, NNZ, NPT, NBLOCK, DZ0, Z0, &
               X, Z, XP, ZP, AS, WT, &
               XSR, YSR, ZSR, MSR, MSR1, FSR, ICSR, VSR, &
               CRR0, CII0, CRREG, CIREG, CRRU, CIIU, SIG, XMIN, XMAX, ZMIN, ZMAX, IG, &
               msg, CRT, CIT, CR, CI, DEBUG_OUTPUT=DEBUG_DBG)
            MESH_BUILT = .TRUE.
         ELSE
            IF (my_rank == 0) THEN
               WRITE (*, *) '--- Reusing global mesh for band', ib, ' (no band meshing override present)'
            END IF
         END IF
      ELSE
         ! Band-meshing block present: rebuild once per band only if band meshing active
         FREQ_MESH = band_fmax(ib)
         WLMAX = CMAX/band_fmin(ib)
         IF (my_rank == 0) THEN
            WRITE (*, *) '--- Building mesh for band', ib, &
               ' (NORD=', NORD_use, ', DX=', DX_use, ', DZ=', DZ_use, ')  FREQ_MESH=', FREQ_MESH
            WRITE (*, *) 'Max source burial for GQG mask (m) = ', src_burial_max
            WRITE (*, *) 'IZ for GQG mesh = ', IZ
         END IF

         CALL BuildMeshForBand_FromCanon( &
            FREQ_MESH, ib, WLMAX, WIDTH, IS0, NORD_use, DX_use, DZ_use, &
            IANISO, ITHOM, IVISCO, PARAM, NPAR, &
            CMIN, CMAX, FREQ_MESH, NSIG, ISO, &
            NSR, NSS, INV, my_rank, NTO, XTO, ZTO, &
            canon, CR0, CI0, &
            IE0, NX, NZ, NZ_ACT, NNX, NNZ, NPT, NBLOCK, DZ0, Z0, &
            X, Z, XP, ZP, AS, WT, &
            XSR, YSR, ZSR, MSR, MSR1, FSR, ICSR, VSR, &
            CRR0, CII0, CRREG, CIREG, CRRU, CIIU, SIG, XMIN, XMAX, ZMIN, ZMAX, IG, &
            msg, CRT, CIT, CR, CI, DEBUG_OUTPUT=DEBUG_DBG)
      END IF

      if (my_rank == 0) write (*, '(A,1X,I0,1X,I0)') 'IE0, IZ', IE0, IZ

!--------------------------------------------------------------------------------------
      freq_loop: DO IFQ = f1b, f2b

         FREQ = FREQN(IFQ)
         ND = NDATA(IFQ)
         IF (my_rank == 0) WRITE (*, *) '--------- Frequency', FREQ, ' (band ', ib, ')'

         IF (ALLOCATED(WD_acq)) DEALLOCATE (WD_acq)
         ALLOCATE (WD_acq(ND))
         CALL ComputeAcqWeight(NFQ, ND, NSS, NS, NR, NSV, NRV, XSR, ZSR, &
                               WD_acq, ierr, my_rank, ACQ_WEIGHT_FILE, DEBUG_OUTPUT=DEBUG_DBG)
         ! WRITE (*, '(A,1X,I0)') 'Data weighting computed for frequency', SUM(WD_acq)

         CALL Count_NKY(FREQ, XSR(NSR), XSR(1), ZSR(1), ZSR(NSR), INV, I25D, NSIG, SIG, NK, FKY, WTK, NPT, my_rank)

         CALL freq_loop_vars(ND, NSS, .TRUE., GT0, G0, SourceScaler, RESID, RESID_m, ampT0, amp0, RESID_amp)

         ITER = 0
         IF (EXTERNAL_DATA) THEN
            CALL Read_seismic_spec(FREQ_FILE, IFQ, FREQ, ND, GT0, my_rank, DEBUG_OUTPUT=DEBUG_INPUT)

            ampT0 = ABS(GT0)
            if (my_rank == 0) WRITE (*, '(A,1X,G0)') 'Data imported for frequency', FREQ

         ELSE
            ! Synthetic path: compute GF for the true model and then G0

            ALLOCATE (GFXT(NSS, 3, NK, NPT), GFYT(NSS, 3, NK, NPT), GFZT(NSS, 3, NK, NPT))

            CALL SYSTEM_CLOCK(iTimes1)

            CALL GF(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD_USE, &
                    AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                    IE0, IS0, DZ0, GFXT, GFYT, GFZT, NBLOCK, IG, &
                    my_rank, n_process, MPI_COMM_WORLD, ITER, SOLVER_KIND, DEBUG_OUTPUT=.TRUE.)

            CALL SYSTEM_CLOCK(iTimes2)
            IF (my_rank == 0) WRITE (*, '(A,1X,G0,1X,A)') 'TIME for 1st forward modeling = ', &
               REAL(iTimes2 - iTimes1)/REAL(rate)/60.D0, ' min'

            IF ((my_rank == 0) .AND. (NSS == 1)) THEN ! .AND. (NRR <= 3)
               CALL Process_GF(ICSR, IS0, NPT, NK, FKY, WTK, YSR, XP, ZP, NNX, NNZ, NTO, XTO, ZTO, &
                               GFXT, GFYT, GFZT, VSR, FREQ, 'T')

            END IF

            CALL MPI_Barrier(mpi_comm_world, ierr)
            CALL SYSTEM_CLOCK(iTimes1)

            CALL Compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, MSR, MSR1, FSR, GFXT, GFYT, GFZT, &
                            YSR, FREQ, FKY, WTK, GT0, WAVELET, SourceScaler, IFQ, my_rank, ampT0, DEBUG_OUTPUT=DEBUG_DBG)
            CALL SYSTEM_CLOCK(iTimes2)
            IF (my_rank == 0) WRITE (*, '(A,1X,G0,1X,A)') 'TIME for GT0 = ', &
               REAL(iTimes2 - iTimes1)/REAL(rate)/60.D0, ' min'

            DEALLOCATE (GFXT, GFYT, GFZT)
         END IF

         ! Output for QC
         CALL WriteObservedDataFile(IFQ, FREQ, ND, GT0, my_rank, INV, DEBUG_OUTPUT=.TRUE.)
         CALL WriteTraceFile('GT0', IFQ, ITER, ND, NS, NR, NSV, NRV, XSR, ZSR, GT0, AMPT0, my_rank, DEBUG_OUTPUT=.TRUE.)

         ! Skip inversion loop for INV == 0 (export only)
         IF (INV == 0) THEN
            CALL freq_loop_vars(ND, NSS, .FALSE., GT0, G0, SourceScaler, RESID, RESID_m, ampT0, amp0, RESID_amp)

            IF (ALLOCATED(FKY)) DEALLOCATE (FKY)
            IF (ALLOCATED(WTK)) DEALLOCATE (WTK)
            IF (ALLOCATED(WD_acq)) DEALLOCATE (WD_acq)

            CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
            CYCLE freq_loop
         END IF

         CALL ComputeDataScalerW(ND, GT0, WD_acq, WD_amp, INV)

         ! if (my_rank == 0) WRITE (*, '(A,1X,G0)') 'Data scaler W for L2 norm:', WD_amp

!-----------------------------------------------------------C PARAMETERS LOOP
         param_loop: DO IA = 1, NGRP_INV_BAND

            GRADr_NORM = 1.0D0
            GRAD_scaled_NORM = 1.0D0
            GRADsc_precond_NORM = 1.0D0

            IM = 0
            stage_has_new_param = .FALSE.
            DO IP = 1, NPAR
               IF (ALLOCATED(INVP_ORD_BAND)) THEN
                  IF (INVP_ORD_BAND(ib, IP) == IA) THEN
                     INVP(IP) = 1
                     IM = IM + 1
                     stage_has_new_param = .TRUE.
                  ELSE
                     INVP(IP) = 0
                  END IF
               ELSEIF (INVP_ORD(IP) == IA) THEN
                  INVP(IP) = 1
                  IM = IM + 1
                  stage_has_new_param = .TRUE.
               ELSE
                  INVP(IP) = 0
               END IF
            END DO
            INVP(NPAR + 1) = 0

            NM = IM*NBLOCK !arrayslenght for block computation versions (frechet based)
            NMM = IM*NPT !arrays lenght for point based computation (forward modeling)

            if (my_rank == 0) then

               WRITE (*, '(A,I10)') 'Total number of active inversion parameters for current sequential stage (IM): ', IM
               WRITE (*, '(A,I10,A,I10)') 'Total model blocks points (NM = IM * NBLOCK): ', NM, ' on grid (NMM = IM * NPT):', NMM
            end if

            IF (.NOT. stage_has_new_param) THEN
               IF (my_rank == 0) THEN
    WRITE (*, '(A,I10,A)') 'Skipping sequential stage IA=', IA, ' because no parameters are assigned to this inversion-order label.'
               END IF
               CYCLE param_loop
            END IF

            ! Reset the balance scale for each parameter stage.
            ! Use a default Q scaling when any Q parameters are active.
            active_q_stage = .FALSE.
            active_nonq_stage = .FALSE.
            DO II = 1, NPAR
               IF (INVP(II) /= 1) CYCLE
               IF (is_q_parameter_name(PARAM(II))) THEN
                  active_q_stage = .TRUE.
               ELSE
                  active_nonq_stage = .TRUE.
               END IF
            END DO
            active_q_only_stage = active_q_stage .AND. .NOT. active_nonq_stage

            BALANCE_SCALE(1:NPAR) = 1.0_dp
            IF (active_q_stage .AND. IANISO > 0) THEN
               BALANCE_SCALE(NPAR - IANISO + 1:NPAR) = 0.005_dp
            END IF
            balance_initialized = .FALSE.
            IF (my_rank == 0) THEN
               PRINT *, 'Gradient balance (initial):', BALANCE_SCALE
               IF (active_q_stage) WRITE (*, '(A)') 'Active stage includes Q parameters; default Q balance applied.'
            END IF

            CALL param_loop_vars(.TRUE., NPAR, ND, NM, NMM, mml, FRECHET, &
                                 m_coarse, m_coarse_prev, m_fine, m_coarse_reg, grad_prev, diff_model, &
                                 BF_grad_res, BF_s_hist, &
                                 GRADr, GRAD_scaled, GRAD_smooth, GRADsc_precond, GMASK, GMASK_FINAL, &
                                 p_k, p_kf, HESS, HESSDI_scaled, HESSDI_scaled_norm, &
                                 GZSMTH, GTAPPER)
            ALLOCATE (p_k_ref(NM), GRAD_scaled_ref(NM))

            if (my_rank == 0) write (*, *) 'IE0, IZ', IE0, IZ
            CALL gradient_mask(GMASK, NX, NZ, IE0, IS0, IZ, DZ_USE, IG, NTO, XTO, ZTO, &
                               my_rank, DEBUG_OUTPUT=DEBUG_DBG)

! Default: no additional taper
            GMASK_FINAL(:) = GMASK(:)

            IF (GZSMTH /= 0) THEN

               IF (my_rank == 0) THEN
                  WRITE (*, *) 'IE0, IZ, GZ1, GZ2 = ', IE0, IZ, GZ1, GZ2
               END IF

               CALL gradient_mask_taper(NX, NZ, IE0, IS0, my_rank, NM, NZ_ACT, IZ, &
                                        GZ1, GZ2, DZ_USE, IG, NTO, XTO, ZTO, GTAPPER, &
                                        DEBUG_OUTPUT=DEBUG_DBG)

               GMASK_FINAL(:) = GMASK(:)*GTAPPER(:)

            END IF

            FCOST_prev = 0.0
            grad0_norm = 0.0
            recompute_model_scaling = .TRUE.

            balance_ref_param = 0

            iter_loop: DO ITER = 1, MAXITER!------------------------------------------------------START ITERATION

               total_iter_count = total_iter_count + 1
               CALL SYSTEM_CLOCK(iter_clock_start)

               IF (INV >= 1) NSSt = NSR

               if (my_rank == 0) write (*, *) ' '
               IF (my_rank == 0) WRITE (*, *) 'Frequency:', IFQ, 'Iteration:', ITER

               ALLOCATE (GFX(NSSt, 3, NK, NPT), GFY(NSSt, 3, NK, NPT), GFZ(NSSt, 3, NK, NPT))

               CALL SYSTEM_CLOCK(iTimes1)

               CALL GF(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CRR0, CII0, NORD_USE, &
                       AS, WT, FREQ, NK, FKY, NSR, NSSt, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                       IE0, IS0, DZ0, GFX, GFY, GFZ, NBLOCK, IG, &
                       my_rank, n_process, mpi_comm_world, ITER, SOLVER_KIND, DEBUG_OUTPUT=.TRUE.)

               CALL SYSTEM_CLOCK(iTimes2)

            IF (my_rank == 0 .AND. ITER == 1) WRITE (*, '(A,1X,G0,1X,A)') 'TIME for 2d forward modeling_i calc = ', REAL(iTimes2 - iTimes1)/REAL(rate)/60.D0, 'min'

               if ((my_rank == 0) .AND. (NSS == 1) .AND. (NRR <= 3)) then
                  CALL Process_GF(ICSR, IS0, NPT, NK, FKY, WTK, YSR, XP, ZP, NNX, NNZ, NTO, XTO, ZTO, &
                                  GFX, GFY, GFZ, VSR, FREQ, 'I')
               end if

               !     (5) pick up the synthetic observed data G0(NDATA,NFQ) at geophones         C

               IF ((INV == 2)) THEN
                  IF (IA == 1 .AND. ITER == 1) THEN
                     call compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, &
                                     MSR, MSR1, FSR, GFX, GFY, GFZ, &
                                     YSR, FREQ, FKY, WTK, G0, WAVELET, SourceScaler, IFQ, my_rank, amp0, DEBUG_OUTPUT=DEBUG_DBG)

                     ! CALL Calibrate_Source(IFQ, ND, GT0, G0, &
                     !                       SourceScaler, Wd_Acq, &
                     !                       NS, NSS, &
                     !                       my_rank, DEBUG_OUTPUT=.TRUE.)

                     CALL Calibrate_PerSource(IFQ, ND, GT0, G0, &
                                              SourceScaler, Wd_Acq, &
                                              NS, NSS, &
                                              my_rank, DEBUG_OUTPUT=.TRUE.)
                     amp0 = ABS(G0)

                  ELSE
                     call compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, &
                                     MSR, MSR1, FSR, GFX, GFY, GFZ, &
                                     YSR, FREQ, FKY, WTK, G0, WAVELET, SourceScaler, IFQ, my_rank, amp0, DEBUG_OUTPUT=DEBUG_DBG)
                  END IF

               ELSE !(synthetic case)
                  SourceScaler(:) = (1.0_dp, 0.0_dp)
                  call compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, &
                                  MSR, MSR1, FSR, GFX, GFY, GFZ, &
                                  YSR, FREQ, FKY, WTK, G0, WAVELET, SourceScaler, IFQ, my_rank, amp0, DEBUG_OUTPUT=DEBUG_DBG)
               END IF

               CALL WriteTraceFile('GO', IFQ, ITER, ND, NS, NR, NSV, NRV, XSR, ZSR, G0, AMP0, my_rank, DEBUG_OUTPUT=DEBUG_QC)

! !----------------- Frechét Derivative Computation -------------------

               CALL SYSTEM_CLOCK(iTimes1)

               CALL QFRECHET(ITER, PARAM, IFQ, FREQ, NSR, ND, NX, NZ, NNX, NNZ, X, IE0, DZ0, NPT, NBLOCK, &
                             IANISO, ITHOM, IVISCO, NPAR, INVP, IS0, CRR0, CII0, &
                             NS, NSV, NR, NRV, NTO, XTO, ZTO, YSR, VSR, &
                             AS, WT, NK, FKY, WTK, GFX, GFY, GFZ, &
                             NM, FRECHET, WAVELET, SourceScaler, NORD_USE, NCOMP, NCOMPS, IG, &
                             my_rank, n_process, MPI_COMM_WORLD, DEBUG_OUTPUT=DEBUG_DBG)

               CALL SYSTEM_CLOCK(iTimes2)
               IF (my_rank == 0) THEN
                  WRITE (*, *)
                  WRITE (*, '(A,1X,G0,1X,A)') 'TIME for Frechet calc = ', REAL(iTimes2 - iTimes1)/REAL(rate)/60.D0, 'min'
                  WRITE (*, *)
               END IF

               ! Build Frechet-based column scalers once (first iteration only),
               ! then keep them fixed for sThubsequent iterations at this (IFQ, IA) loop.
               IF (ITER == 1 .AND. IFQ == 1) THEN
                  CALL ComputeFrechetScalers(FRECHET, ND, NBLOCK, NM, NPAR, INVP, &
                                             WD_amp, F_SCALE, F_SCALE_L2, my_rank, DEBUG_OUTPUT=.TRUE.)
                  ! F_SCALE=1.0_dp
               end IF

               ! Scale the model-side vectors once after the first Frechet

               IF (recompute_model_scaling) THEN
                  CALL ComputeModelScaling(CRR0, CII0, &
                                           NPT, NBLOCK, IG, &
                                           NNX, NNZ, NX, NZ, NORD, &
                                           NTO, XTO, ZTO, &
                                           NPAR, IANISO, INVP, &
                                           PAR_SCALE, GMASK, &
                                           m_fine, m_coarse_prev, &
                                           m_coarse, m_coarse_reg, &
                                           NM, IFQ, m_coarse0_norm, &
                                           my_rank, DEBUG_OUTPUT=.TRUE.)
                  recompute_model_scaling = .FALSE.
               END IF

!             !----------------------------------------------------------------------C
               IF (my_rank == 0) WRITE (*, *)
               IF (my_rank == 0) WRITE (*, *) ' compute the residual and the cost function'

               CALL ComputeFCOST0(ND, G0, GT0, WD_acq, WD_amp, ITER, my_rank, &
                                  NM, m_coarse, m_coarse_reg, REG_LAMBDA, &
                                  RESID, RESID_L2, RESID_RMS, FCOST0, &
                                  NS, NR, NSV, NRV, XSR, ZSR, DEBUG_OUTPUT)

               DEALLOCATE (GFX, GFY, GFZ)
               CALL MPI_Barrier(mpi_comm_world, ierr)
               IF (NSR == 2) EXIT ITER_LOOP

! !----------------------------------------------------------------------C
!              Computate Hessian diagonal inverse approximation and gradient            C

               IF (ITER == 1) then
                  IF (USE_PRECOND) THEN
                     CALL ComputeHessDI(FRECHET, ND, IM, NBLOCK, NM, NPAR, INVP, &
                                        GMASK, F_SCALE, PAR_SCALE, &
                                        HESS, HESSDI_scaled, REG_LAMBDA, &
                                        WD_acq, WD_amp, BETA, &
                                        NX, NZ, IG, FREQ, CREF, my_rank, &
                                        DEBUG_OUTPUT=.TRUE.)
                  ELSE
                     ! HESS(1:NM) = 1.0_dp
                     HESSDI_scaled(1:NM) = 1.0_dp
                  END IF
               end if

      CALL compute_gradient_scaled(FRECHET, RESID, NPAR, F_SCALE, PAR_SCALE, BALANCE_SCALE, GMASK, ND, NM, &
                                            m_coarse, m_coarse_reg, REG_LAMBDA, &
                                            GRADr, GRAD_scaled, GRADr_NORM, GRAD_scaled_NORM, &
                                            NBLOCK, INVP, USE_PRECOND=USE_PRECOND, my_rank=my_rank, &
                                            HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                                            GRADsc_Precond_NORM=GRADsc_precond_NORM, DEBUG_OUTPUT=.TRUE., &
                                            PARAM=Param, FREQ=freq, ITER=iter, &
                                            NX=nx, NZ=nz, IG=IG, NTO=nto, XTO=xto, ZTO=zto, IE0=IE0, IS0=IS0, &
                                            XMINC=canon%XMINC, XMAXC=canon%XMAXC, ZMINC=canon%ZMINC, ZMAXC=canon%ZMAXC)

       IF (.NOT. balance_initialized .AND. ITER == 1 .AND. COUNT(INVP == 1) > 1 .AND. .NOT. active_q_only_stage) THEN
                  balance_ref_param = 0
                  DO II = 1, NPAR
                     IF (INVP(II) /= 1) CYCLE
                     IF (INDEX(TRIM(ADJUSTL(Param(II))), 'C33') > 0) THEN
                        balance_ref_param = II
                        EXIT
                     END IF
                  END DO
                  IF (balance_ref_param == 0) THEN
                     DO II = 1, NPAR
                        IF (INVP(II) == 1) THEN
                           balance_ref_param = II
                           EXIT
                        END IF
                     END DO
                  END IF

                  CALL ComputeParamBalanceScale(GRADr, GMASK, NPAR, NBLOCK, INVP, &
                                                balance_ref_param, BALANCE_PCTL, BALANCE_WMIN, BALANCE_WMAX, &
                                                BALANCE_SCALE, my_rank, PARAM=Param, DEBUG_OUTPUT=.TRUE.)

                  IF (active_q_stage .AND. IANISO > 0) THEN
                     BALANCE_SCALE(NPAR - IANISO + 1:NPAR) = 0.005_dp
                  END IF

                  balance_initialized = .TRUE.
                  IF (my_rank == 0) THEN
                     WRITE (*, '(A)') 'Automatic parameter balance initialized; recomputing gradient with fixed BALANCE_SCALE.'
                  END IF
                  CALL compute_gradient_scaled(FRECHET, RESID, NPAR, F_SCALE, PAR_SCALE, BALANCE_SCALE, GMASK, ND, NM, &
                                               m_coarse, m_coarse_reg, REG_LAMBDA, &
                                               GRADr, GRAD_scaled, GRADr_NORM, GRAD_scaled_NORM, &
                                               NBLOCK, INVP, USE_PRECOND=USE_PRECOND, my_rank=my_rank, &
                                               HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                                               GRADsc_Precond_NORM=GRADsc_precond_NORM, DEBUG_OUTPUT=.TRUE., &
                                               PARAM=Param, FREQ=freq, ITER=iter, &
                                               NX=nx, NZ=nz, IG=IG, NTO=nto, XTO=xto, ZTO=zto, IE0=IE0, IS0=IS0, &
                                               XMINC=canon%XMINC, XMAXC=canon%XMAXC, ZMINC=canon%ZMINC, ZMAXC=canon%ZMAXC)
               END IF

               IF (USE_GR_SMOOTH) THEN
                  CALL smooth_gradient_blocks( &
                     GRAD_in=GRAD_scaled, GMASK=GMASK_FINAL, &
                     NBLOCK=NBLOCK, INVP=INVP, NPAR=NPAR, &
                     NX=NX, NZ=NZ, IG=IG, &
                     FREQ=FREQ, CREF=CREF, my_rank=my_rank, &
                     USE_GR_SMOOTH=.TRUE., USE_PRECOND=USE_PRECOND, &
                     HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                     GRADsc_Precond_NORM=GRADsc_precond_NORM, &
                     PARAM=Param, ITER=iter, NTO=nto, XTO=xto, ZTO=zto, IE0=IE0, IS0=IS0, &
                     XMINC=canon%XMINC, XMAXC=canon%XMAXC, ZMINC=canon%ZMINC, ZMAXC=canon%ZMAXC, &
                     DEBUG_OUTPUT=.FALSE., &
                     SIGMA_X=sigma_grad_x, SIGMA_Z=sigma_grad_z, &
                     GRAD_out=GRAD_smooth) !
                  GRAD_scaled(1:NM) = GRAD_smooth(1:NM)
                  CALL recompute_masked_norms(GRAD_scaled, GRAD_scaled_NORM, INVP, NBLOCK)
               END IF

               if (my_rank == 0) WRITE (*, '(A)') 'Gradient finished'

! !---------------- LBFGS METHOD --------------------------------------

               ! For INV == 0 (forward only with gradient output), skip inversion and move to next frequency
               IF (INV == 0) THEN
                  exit_iteration_loop = .TRUE.
               ELSE
                  exit_iteration_loop = .FALSE.
                  CALL LBFGS_update(exit_iteration_loop)
               END IF
               ! CALL TNGN_update(exit_iteration_loop)
               IF (exit_iteration_loop) EXIT iter_loop

!!---------------- Cleanup Per Iteration ------------------------------
               CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
            END DO iter_loop   ! End ITERATION loop

!!---------------- Cleanup Per Parameter ------------------------------
            ! write (*, *) ' calling param_loop_vars'

            CALL param_loop_vars(.FALSE., NPAR, ND, NM, NMM, mml, FRECHET, &
                                 m_coarse, m_coarse_prev, m_fine, m_coarse_reg, grad_prev, diff_model, &
                                 BF_grad_res, BF_s_hist, &
                                 GRADr, GRAD_scaled, GRAD_smooth, GRADsc_precond, GMASK, GMASK_FINAL, &
                                 p_k, p_kf, HESS, HESSDI_scaled, HESSDI_scaled_norm, &
                                 GZSMTH, GTAPPER)
            DEALLOCATE (p_k_ref, GRAD_scaled_ref)
            IF (NSR == 2) EXIT PARAM_LOOP
            CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
         END DO param_loop
!---------------- Cleanup Per Frequency ------------------------------

         CALL freq_loop_vars(ND, NSS, .FALSE., GT0, G0, SourceScaler, RESID, RESID_m, ampT0, amp0, RESID_amp)
         IF (ALLOCATED(FKY)) DEALLOCATE (FKY)
         IF (ALLOCATED(WTK)) DEALLOCATE (WTK)
         IF (NSR == 2) CYCLE freq_loop
         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

      END DO freq_loop  ! End FREQUENCY loop

!---------------- Cleanup Per Frequency band variable related to purely ND------------------------------
      if (my_rank == 0) write (*, *) ' projecting back to physical grid for next frequency'
      CALL Project_CRR0_2_CR0( &
         NNX, NNZ, NPT, XP, ZP, &
         NTO, XTO, ZTO, IS0, &
         IANISO, CRR0, &
         canon, &
         CR0)
      CALL Project_CRR0_2_CR0( &
         NNX, NNZ, NPT, XP, ZP, &
         NTO, XTO, ZTO, IS0, &
         IANISO, CII0, &
         canon, &
         CI0)

      IF (my_rank == 0) THEN
         write (*, *) ' end of  frequency ', IFQ
         WRITE (*, *) 'CR0 shape', SIZE(CR0, 1), SIZE(CR0, 2), SIZE(CR0)
         WRITE (*, *) 'C0 shape', SIZE(CI0, 1), SIZE(CI0, 2), SIZE(CI0)
      END IF

   END DO band_loop

! !---------------- Final Cleanup ------------------------------
! !
   DEALLOCATE (XP, ZP, CRR0, CII0, SIG, wavelet, PAR_SCALE, COMP, CRRU, CIIU, F_SCALE, F_SCALE_L2)
   DEALLOCATE (XSR, YSR, ZSR, ICSR, VSR, ZTO, XTO, MSR, MSR1, FSR, X, Z, NSV, NRV, NS, NR, NDATA, INVP, INVP_ORD)
   IF (ALLOCATED(CR)) DEALLOCATE (CR)
   IF (ALLOCATED(CI)) DEALLOCATE (CI)

   DEALLOCATE (AS, WT, CRREG, CIREG, PAR_SCALE_REG)

   IF (my_rank == 0) THEN
      CLOSE (70)
      CLOSE (62)
      CLOSE (63)
      CLOSE (64)
      CLOSE (65)
      CLOSE (66)
      ! CLOSE (67)
      WRITE (*, *) 'RESULT FILES CLOSED'
   END IF
   t_global_end = MPI_Wtime()
   IF (my_rank == 0) THEN
      PRINT *, 'Total iterations executed = ', total_iter_count
      PRINT *, 'Time of operation was ', (t_global_end - t_global_begin)/60.D0, 'min'
   END IF
   CALL MPI_FINALIZE(ierr)

contains

   SUBROUTINE freq_loop_vars(ND_in, NSS_in, DO_ALLOC, GT0_in, G0_in, SourceScaler_in, RESID_in, RESID_m_in, ampT0_in, amp0_in, RESID_amp_in)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ND_in, NSS_in
      LOGICAL, INTENT(IN) :: DO_ALLOC
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT) :: GT0_in(:), RESID_in(:), RESID_m_in(:), G0_in(:), SourceScaler_in(:)
      real(dp), ALLOCATABLE, INTENT(INOUT) :: ampT0_in(:), amp0_in(:), RESID_amp_in(:)

      IF (DO_ALLOC) THEN
         ! Complex arrays
         ALLOCATE (GT0_in(ND_in)); GT0_in = (0.0D0, 0.0D0)
         ALLOCATE (RESID_in(ND_in)); RESID_in = (0.0D0, 0.0D0)
         ALLOCATE (RESID_m_in(ND_in)); RESID_m_in = (0.0D0, 0.0D0)
         ALLOCATE (G0_in(ND_in)); G0_in = (0.0D0, 0.0D0)
         ALLOCATE (SourceScaler_in(NSS_in)); SourceScaler_in = (1.0D0, 0.0D0)

         ! Real arrays
         ALLOCATE (ampT0_in(ND_in)); ampT0_in = 0.0D0
         ALLOCATE (amp0_in(ND_in)); amp0_in = 0.0D0
         ALLOCATE (RESID_amp_in(ND_in)); RESID_amp_in = 0.0D0

      ELSE
         IF (ALLOCATED(GT0_in)) DEALLOCATE (GT0_in)
         IF (ALLOCATED(RESID_in)) DEALLOCATE (RESID_in)
         IF (ALLOCATED(RESID_m_in)) DEALLOCATE (RESID_m_in)
         IF (ALLOCATED(G0_in)) DEALLOCATE (G0_in)
         IF (ALLOCATED(ampT0_in)) DEALLOCATE (ampT0_in)
         IF (ALLOCATED(amp0_in)) DEALLOCATE (amp0_in)
         IF (ALLOCATED(RESID_amp_in)) DEALLOCATE (RESID_amp_in)
         IF (ALLOCATED(SourceScaler_in)) DEALLOCATE (SourceScaler_in)
      END IF
   END SUBROUTINE freq_loop_vars

   SUBROUTINE param_loop_vars(DO_ALLOC, NPAR, ND_in, NM, NMM, mml, FRECHET, &
                              m_coarse, m_coarse_prev, m_fine, m_coarse_reg, grad_prev, diff_model, &
                              BF_grad_res, BF_s_hist, &
                              GRADr, GRAD_scaled, GRAD_smooth, GRADsc_precond, GMASK, GMASK_FINAL, &
                              p_k, p_kf, HESS, HESSDI_scaled, HESSDI_scaled_norm, &
                              GZSMTH, GTAPPER)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ND_in, NM, NMM, mml, NPAR
      INTEGER, INTENT(IN) :: GZSMTH
      LOGICAL, INTENT(IN) :: DO_ALLOC

      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: m_coarse(:), m_coarse_prev(:), m_fine(:), m_coarse_reg(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: grad_prev(:), diff_model(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: BF_grad_res(:, :), BF_s_hist(:, :)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: GRADr(:), GRAD_scaled(:), GRADsc_precond(:), GMASK(:), GRAD_smooth(:), GMASK_FINAL(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: p_k(:), p_kf(:), HESS(:), HESSDI_scaled(:), HESSDI_scaled_norm(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: GTAPPER(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT) :: FRECHET(:, :)

      IF (DO_ALLOC) THEN

         ALLOCATE (GRADr(NM)); GRADr = 0.0_dp
         ALLOCATE (GRAD_scaled(NM)); GRAD_scaled = 0.0_dp
         ALLOCATE (GRADsc_precond(NM)); GRADsc_precond = 0.0_dp
         ALLOCATE (GRAD_smooth(NM)); GRAD_smooth = 0.0_dp
         ALLOCATE (GMASK(NM)); GMASK = 1.0_dp
         ALLOCATE (GMASK_FINAL(NM)); GMASK_FINAL = 1.0_dp

! ---- L-BFGS (coarse / NM) ----
         ALLOCATE (m_coarse_prev(NM)); m_coarse_prev = 0.0_dp
         ALLOCATE (m_coarse(NM)); m_coarse = 0.0_dp
         ALLOCATE (m_coarse_reg(NM)); m_coarse_reg = 0.0_dp

         ALLOCATE (grad_prev(NM)); grad_prev = 0.0_dp
         ALLOCATE (BF_grad_res(NM, mml)); BF_grad_res = 0.0_dp
         ALLOCATE (BF_s_hist(NM, mml)); BF_s_hist = 0.0_dp
         ALLOCATE (p_k(NM)); p_k = 0.0_dp
         ALLOCATE (HESS(NM)); HESS = 0.0_dp
         ALLOCATE (HESSDI_scaled(NM)); HESSDI_scaled = 0.0_dp
         ALLOCATE (HESSDI_scaled_norm(NPAR)); HESSDI_scaled_norm = 0.0_dp

! ---- Grid-space (NMM = NPT * #active params) ----
         ALLOCATE (m_fine(NMM)); m_fine = 0.0_dp
         ALLOCATE (p_kf(NMM)); p_kf = 0.0_dp
         ALLOCATE (diff_model(NMM)); diff_model = 0.0_dp

         ! 2D complex
         ALLOCATE (FRECHET(ND_in, NM)); FRECHET = (0.0_dp, 0.0_dp)

         ! Vertical taper buffer (same size as NM)
         IF (GZSMTH /= 0) THEN
            IF (ALLOCATED(GTAPPER)) DEALLOCATE (GTAPPER)
            ALLOCATE (GTAPPER(NM))
            GTAPPER = 1.0_dp
         END IF

      ELSE

         IF (ALLOCATED(m_coarse_prev)) DEALLOCATE (m_coarse_prev)
         IF (ALLOCATED(m_coarse)) DEALLOCATE (m_coarse)
         IF (ALLOCATED(m_fine)) DEALLOCATE (m_fine)
         IF (ALLOCATED(m_coarse_reg)) DEALLOCATE (m_coarse_reg)
         IF (ALLOCATED(grad_prev)) DEALLOCATE (grad_prev)
         IF (ALLOCATED(diff_model)) DEALLOCATE (diff_model)
         IF (ALLOCATED(BF_grad_res)) DEALLOCATE (BF_grad_res)
         IF (ALLOCATED(BF_s_hist)) DEALLOCATE (BF_s_hist)
         IF (ALLOCATED(FRECHET)) DEALLOCATE (FRECHET)
         IF (ALLOCATED(GRADr)) DEALLOCATE (GRADr)
         IF (ALLOCATED(GRAD_scaled)) DEALLOCATE (GRAD_scaled)
         IF (ALLOCATED(GRAD_smooth)) DEALLOCATE (GRAD_smooth)
         IF (ALLOCATED(GRADsc_precond)) DEALLOCATE (GRADsc_precond)
         IF (ALLOCATED(GMASK)) DEALLOCATE (GMASK)
         IF (ALLOCATED(GMASK_FINAL)) DEALLOCATE (GMASK_FINAL)
         IF (ALLOCATED(p_k)) DEALLOCATE (p_k)
         IF (ALLOCATED(p_kf)) DEALLOCATE (p_kf)
         IF (ALLOCATED(HESS)) DEALLOCATE (HESS)
         IF (ALLOCATED(HESSDI_scaled)) DEALLOCATE (HESSDI_scaled)
         IF (ALLOCATED(HESSDI_scaled_norm)) DEALLOCATE (HESSDI_scaled_norm)
         IF (ALLOCATED(GTAPPER)) DEALLOCATE (GTAPPER)

      END IF

   END SUBROUTINE param_loop_vars

   SUBROUTINE LBFGS_update(exit_iter_loop)
      IMPLICIT NONE
      INTEGER :: ia_loc, im_loc, cs_loc, ce_loc
      REAL(dp), ALLOCATABLE :: pk_up(:)
      REAL(dp) :: pk_min, pk_max, pk_l2
      LOGICAL, INTENT(OUT) :: exit_iter_loop

      exit_iter_loop = .FALSE.

      IF (my_rank == 0) WRITE (*, '(A)') 'Start LBFGS '
      CALL SYSTEM_CLOCK(iTimes1)

      CALL lbfgs_basic(iter, NM, mml, INVP, NPAR, NBLOCK, &
                       GRADsc_Precond, GRADsc_Precond_NORM, HESSDI_scaled, &
                       p_k, grad_prev, BF_grad_res, BF_s_hist, &
                       m_coarse, m_coarse_prev, FCOST0, FCOST_prev, USE_LBFGS_TYPE, USE_PRECOND, my_rank)

      grad_prev(1:NM) = GRADsc_Precond(1:NM)
      FCOST_prev = FCOST0
      GRAD_scaled_ref(1:NM) = GRAD_scaled(1:NM)
      p_k_ref(1:NM) = p_k(1:NM)
      FCOST0_ref = FCOST0
      dot_ref = DOT_PRODUCT(GRAD_scaled_ref(1:NM), p_k_ref(1:NM))

      ! Diagnostics in physical parameter scale: unscale p_k blockwise by SCALER(IA)
      ALLOCATE (pk_up(NBLOCK))
      im_loc = 0
      IF (my_rank == 0) WRITE (*, '(A)') 'LBFGS p_k (unscaled) stats per active parameter:'
      DO ia_loc = 1, NPAR
         IF (INVP(ia_loc) /= 1) CYCLE
         im_loc = im_loc + 1
         cs_loc = (im_loc - 1)*NBLOCK + 1
         ce_loc = im_loc*NBLOCK

         pk_up(1:NBLOCK) = p_k(cs_loc:ce_loc)*PAR_SCALE(ia_loc)*F_SCALE(ia_loc) ! unscale to physical parameter units for diagnostics
         ! pk_up(1:NBLOCK) = p_k(cs_loc:ce_loc)
         pk_min = MINVAL(pk_up(1:NBLOCK))
         pk_max = MAXVAL(pk_up(1:NBLOCK))
         pk_l2 = SQRT(SUM(pk_up(1:NBLOCK)*pk_up(1:NBLOCK)))

         IF (my_rank == 0) THEN
            WRITE (*, '(A,I3,2X,A,A,2X,A,1PE12.5,2X,A,1PE12.5,2X,A,1PE12.5)') &
               '  IA=', ia_loc, 'PAR=', TRIM(PARAM(ia_loc)), &
               'min=', pk_min, 'max=', pk_max, 'L2=', pk_l2
         END IF
      END DO
      DEALLOCATE (pk_up)

      CALL lbfgs_project_direction(NPAR, INVP, NBLOCK, NPT, NNX, NNZ, NX, NZ, NTO, XTO, ZTO, &
                                   IG, NORD, p_k, p_kf, ITER, my_rank, &
                                   unit_block=69, unit_point=67, DEBUG_OUTPUT=DEBUG_DBG, &
                                   pnorm_out=pnorm)

      CALL SYSTEM_CLOCK(iTimes2)
      CALL MPI_Barrier(mpi_comm_world, ierr)

      IF (my_rank == 0) WRITE (*, '(A)') 'Start linesearch'

      ! IF ((ITER == 1 .AND. INV < 2) .OR. (ITER <= 3 .AND. INV == 2)) THEN
      IF (ITER == 1) THEN
         ALPHA00 = 1.0_dp
         ALPHA00_prev = 1.0_dp
         CALL LINESEARCH_BACKTRACK(p_k, p_kf, GRAD_scaled, m_fine, FCOST0, SCALER, F_SCALE, PAR_SCALE, INVP, &
                                   CRR0, CII0, ICSR, VSR, FSR, m_coarse, m_coarse_reg, REG_LAMBDA, &
                                   NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD_USE, X, DZ0, ALPHA00, &
                                   NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, &
                                   MSR, MSR1, COMP, FREQ, GT0, WAVELET, SourceScaler, WD_acq, WD_amp, IFQ, NK, FKY, WTK, AS, WT, &
                                   I25D, IS0, IE0, ITER, ALPHA, NBLOCK, IG, &
                                   my_rank, n_process, mpi_comm_world, IFLAG, FCOSTMIN, SOLVER_KIND, DEBUG_OUTPUT=.TRUE.)

      ELSE
         ALPHA00 = MIN(1.0_dp, ALPHA00_prev)

      CALL LINESEARCH_MT(ALPHA00, p_k, p_kf, GRAD_scaled, m_fine, SCALER, F_SCALE, PAR_SCALE, BALANCE_SCALE, INVP, CRR0, CII0, ICSR, VSR, FSR, &
                            m_coarse, m_coarse_reg, REG_LAMBDA, NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD_USE, X, DZ0, &
                            NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, MSR, MSR1, COMP, &
                            FREQ, GT0, WAVELET, SourceScaler, WD_amp, WD_acq, IFQ, NK, FKY, WTK, AS, WT, I25D, IS0, IE0, ITER, &
                   FRECHET, GMASK, GMASK_FINAL, GRAD_scaled_NORM, GRADr_NORM, CREF, HESSDI_scaled, GRADr, GRAD_smooth, sigma_grad_x, sigma_grad_z, &
                            GRADsc_precond, GRADsc_precond_NORM, USE_PRECOND, NBLOCK, IG, USE_GR_SMOOTH, &
                            ALPHA, FCOST0, FCOSTMIN, IFLAG, SOLVER_KIND, my_rank, n_process, mpi_comm_world, &
                            PARAM, ITHOM, IVISCO, NCOMP, NCOMPS, FCOST0_ref, GRAD_scaled_ref, p_k_ref, dot_ref, &
                            canon%XMINC, canon%XMAXC, canon%ZMINC, canon%ZMAXC, DEBUG_OUTPUT=.TRUE.)
      END IF

      IF (ALPHA > 0.0_dp .AND. ALPHA == ALPHA) ALPHA00_prev = ALPHA*2.0_dp ! increase the initial trial step for the next iteration if current linesearch succeeded with a positive step

      ! CALL LINESEARCH_STD(p_k, p_kf, GRAD_scaled, m_fine, SCALER, F_SCALE, PAR_SCALE, BALANCE_SCALE, INVP, CRR0, CII0, ICSR, VSR, FSR, &
      !                     NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD_USE, X, DZ0, &
      !                     NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSSt, NS, NR, NSV, NRV, ND, MSR, MSR1, COMP, &
      !                     FREQ, GT0, WAVELET, SourceScaler, WD_amp, WD_acq, IFQ, NK, FKY, WTK, AS, WT, I25D, IS0, IE0, ITER, &
      !                     FRECHET, GMASK, GRAD_scaled_NORM, GRADr_NORM, CREF, HESSDI_scaled, GRADr, GRAD_smooth, &
      !                     GRADsc_precond, GRADsc_precond_NORM, USE_PRECOND, NBLOCK, IG, USE_GR_SMOOTH, &
      !                     m_coarse, m_coarse_reg, REG_LAMBDA, ALPHA, FCOST0, FCOSTMIN, IFLAG, SOLVER_KIND, my_rank, n_process, mpi_comm_world, &
      !                     PARAM, ITHOM, IVISCO, NCOMP, NCOMPS, FCOST0_ref, GRAD_scaled_ref, p_k_ref, dot_ref, DEBUG_OUTPUT=.TRUE., &
      !                     sigma_grad_x, sigma_grad_z, canon%XMINC, canon%XMAXC, canon%ZMINC, canon%ZMAXC)

      SELECT CASE (IFLAG)
      CASE (2)
         IF (my_rank == 0) WRITE (*, '(A)') 'Convergence reached, updating parameter grid'

         m_coarse_prev(1:NM) = m_coarse(1:NM)
         ALPHA_prev = ALPHA

         IF (my_rank == 0) WRITE (*, '(A,1X,G0)') 'ALPHA from line search = ', ALPHA
         IF (my_rank == 0) WRITE (*, *) 'FCOST0 after linesearch', FCOSTMIN

         m_coarse(1:NM) = m_coarse_prev(1:NM) + ALPHA*p_k(1:NM) ! update in coarse model space (scaled)

         s_norm = DNRM2(NM, m_coarse(1:NM) - m_coarse_prev(1:NM), 1)

         IF (my_rank == 0) WRITE (*, '(A,1X,G0)') 'Norm of s (s_norm) = ', s_norm
         ! IF (s_norm < 1.0e-08_dp) THEN
         !    WRITE (*, *) 'ERROR: accepted step but s_norm ~ 0. Model not updated.'
         !    exit_iter_loop = .TRUE.
         !    RETURN
         ! END IF

         GRAD_scaled_norm_conv = DNRM2(NM, GRAD_scaled(1:NM), 1)
         CALL GRID2D_UPDATE_MODEL(p_kf, m_fine, CRR0, CII0, INVP, IANISO, ALPHA, NPT, PAR_SCALE, NPAR, F_SCALE)

         IF (ITER <= 2) THEN
            CALL GRID2D_UPDATE_LBFG(p_kf, CRRU, CIIU, INVP, IANISO, ALPHA, &
                                    NPT, PAR_SCALE, NPAR)
         END IF

         CALL output_update(my_rank, NPAR, IANISO, INVP, CRR0, CRRU, CII0, CIIU, PARAM, &
                            FREQ, ITER, NPT, NBLOCK, NNX, NNZ, XP, ZP, NTO, XTO, ZTO, &
                            IE0, IS0, canon%XMINC, canon%XMAXC, canon%ZMINC, canon%ZMAXC, &
                            DEBUG_OUTPUT=.false.)
         WRITE (*, *) ' '

         CALL SYSTEM_CLOCK(iter_clock_end)
         iter_tick_delta = iter_clock_end - iter_clock_start
         IF (iter_tick_delta < 0) THEN
            IF (clock_max > 0) THEN
               iter_tick_delta = iter_tick_delta + clock_max
            ELSE
               iter_tick_delta = -iter_tick_delta
            END IF
         END IF
         iter_elapsed = DBLE(iter_tick_delta)/DBLE(rate)/60.D0

         IF (my_rank == 0) THEN
            DO I = 1, NPAR
               IF (INVP(I) == 1) THEN
                  WRITE (66, '(F7.2,"   ",I4,"   ",ES17.10,"   ",ES17.10,"   ",ES17.10,"   ",F12.3,"   ",I4,"   ",I4,"   ",I4,"   ",I4,"   ",I4,"   ",I8)') &
                     FREQ, ITER, FCOST_prev, GRAD_scaled_NORM(I), RESID_RMS, iter_elapsed, ib, IFQ, IA, I, IM, total_iter_count
               END IF
            END DO
         END IF

         CALL check_convergence(ITER, MAXITER, NM, &
                                m_coarse_prev, m_coarse, &
                                FCOST_prev, FCOSTMIN, &
                                m_coarse0_norm, GRAD_scaled_norm_conv, grad0_norm, &
                                cost_conv_tol, model_conv_tol, &
                                MIN_ITER_CHECK, &
                                my_rank, converged)

         IF (converged) THEN
            CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
            exit_iter_loop = .TRUE.
            RETURN
         END IF

      CASE (11)
         IF (my_rank == 0) WRITE (*, *) 'Line search failed (IFLAG=', IFLAG, '), inversion diverged.'

         CALL SYSTEM_CLOCK(iter_clock_end)
         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
         exit_iter_loop = .TRUE.
         RETURN

      CASE DEFAULT
         IF (my_rank == 0) WRITE (*, *) 'Unknown IFLAG code from LINESEARCH_MT:', IFLAG
         STOP
      END SELECT

      CALL SYSTEM_CLOCK(iTimes2)
      IF (my_rank == 0) WRITE (*, '(A,1X,G0,1X,A)') 'Time for LBFGS = ', REAL(iTimes2 - iTimes1)/REAL(rate)/60.D0, 'min'
   END SUBROUTINE LBFGS_update

   ! SUBROUTINE TNGN_update(exit_iter_loop)
   !    IMPLICIT NONE
   !    LOGICAL, INTENT(OUT) :: exit_iter_loop
   !    ! Local temporaries
   !    REAL(dp), ALLOCATABLE :: dm(:), p_k_local(:), p_kf_local(:)
   !    REAL(dp), ALLOCATABLE :: CRR0_old(:, :), CII0_old(:, :), m_fine_old(:)
   !    INTEGER :: i, ntries, it_used
   !    REAL(dp) :: alpha, alpha_gn, relres, tol_rel
   !    INTEGER :: maxit
   !    REAL(dp) :: FCOST_try, FCOST_saved
   !    INTEGER :: NMM
   !    LOGICAL :: accepted_step

   !    exit_iter_loop = .FALSE.

   !    IF (my_rank == 0) WRITE (*, '(A)') 'Start TNGN'
   !    CALL SYSTEM_CLOCK(iTimes1)

   !    ! Configure GN solver parameters (tunable)
   !    tol_rel = 1.0e-3_dp
   !    maxit = 50
   !    ntries = 6

   !    ! sizes
   !    NMM = SIZE(m_fine)

   !    ALLOCATE (dm(NM)); dm = 0.0_dp
   !    ALLOCATE (p_k_local(NM)); p_k_local = 0.0_dp
   !    ALLOCATE (p_kf_local(NMM)); p_kf_local = 0.0_dp

   !    ! save current model to allow rollback on rejection
   !    ALLOCATE (CRR0_old(SIZE(CRR0, 1), SIZE(CRR0, 2)))
   !    ALLOCATE (CII0_old(SIZE(CII0, 1), SIZE(CII0, 2)))
   !    ALLOCATE (m_fine_old(SIZE(m_fine)))
   !    CRR0_old = CRR0
   !    CII0_old = CII0
   !    m_fine_old = m_fine

   !    FCOST_prev = FCOST0
   !    FCOST_saved = FCOST0

   !    ! Call GN wrapper to compute proposed step (dm -> p_k and p_kf)
   !    CALL GaussNewton_step(ND, NM, FRECHET, GRAD_scaled, HESSDI_scaled, USE_PRECOND, lambda, &
   !                          maxit, tol_rel, dm, p_k_local, p_kf_local, NPAR, INVP, NBLOCK, NPT, &
   !                          SCALER, IANISO, CRR0, CII0, m_fine, NNX, NNZ, NX, NZ, NTO, XTO, ZTO, IG, NORD, &
   !                          alpha_gn, it_used, relres, my_rank)

   !    p_k(1:NM) = p_k_local(1:NM)
   !    p_kf(1:NMM) = p_kf_local(1:NMM)
   !    alpha = alpha_gn
   !    accepted_step = .FALSE.

   !    DO i = 1, ntries
   !       ! apply step to model (grid arrays), do not change global m_coarse here
   !       CALL GRID2D_UPDATE_MODEL(p_kf_local, m_fine, CRR0, CII0, INVP, IANISO, alpha, NPT, SCALER, NPAR)

   !       ! recompute forward model and cost for updated model
   !       ALLOCATE (GFX(NSSt, 3, NK, NPT), GFY(NSSt, 3, NK, NPT), GFZ(NSSt, 3, NK, NPT))
   !       CALL GF(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CRR0, CII0, NORD_USE, &
   !               AS, WT, FREQ, NK, FKY, NSR, NSSt, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
   !               IE0, IS0, DZ0, GFX, GFY, GFZ, NBLOCK, IG, &
   !               my_rank, n_process, mpi_comm_world, ITER, SOLVER_KIND, DEBUG_OUTPUT=.FALSE.)

   !       CALL Compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, MSR, MSR1, FSR, GFX, GFY, GFZ, &
   !                       YSR, FREQ, FKY, WTK, G0, WAVELET, SourceScaler, IFQ, my_rank, amp0, DEBUG_OUTPUT=.FALSE.)

   !       CALL ComputeFCOST0(ND, G0, GT0, WD_acq, WD_amp, ITER, my_rank, RESID, RESID_L2, RESID_RMS, FCOST_try, &
   !                          NS, NR, NSV, NRV, XSR, ZSR, DEBUG_OUTPUT=.FALSE.)

   !       DEALLOCATE (GFX, GFY, GFZ)

   !       IF (my_rank == 0) WRITE (*, '(A,1X,G0)') 'TNGN try alpha=', alpha, ' FCOST_try=', FCOST_try

   !       IF (FCOST_try <= FCOST_saved) THEN
   !          ! accept step
   !          accepted_step = .TRUE.
   !          ALPHA_prev = alpha
   !          m_coarse_prev(1:NM) = m_coarse(1:NM)
   !          m_coarse(1:NM) = m_coarse_prev(1:NM) + alpha*p_k_local(1:NM)
   !          s_norm = DNRM2(NM, m_coarse(1:NM) - m_coarse_prev(1:NM), 1)
   !          GRAD_scaled_norm_conv = DNRM2(NM, GRAD_scaled(1:NM), 1)

   !          FCOST0 = FCOST_try
   !          CALL GRID2D_UPDATE_LBFG(p_kf_local, CRRU, CIIU, INVP, IANISO, alpha, NPT, SCALER, NPAR)

   !          IF (my_rank == 0) THEN
   !             WRITE (*, '(A,1X,I0,1X,A,1X,G0)') 'TNGN accepted after try', i, 'alpha =', alpha
   !             WRITE (*, '(A,1X,I0,1X,A,1X,G0)') 'TNGN inner iterations =', it_used, 'relres =', relres
   !             WRITE (*, '(A,1X,G0)') 'Norm of s (s_norm) = ', s_norm
   !          END IF

   !          CALL output_update(my_rank, NPAR, IANISO, INVP, CRR0, CRRU, CII0, CIIU, PARAM, &
   !                             FREQ, ITER, NPT, NBLOCK, NNX, NNZ, XP, ZP, NTO, XTO, ZTO, &
   !                             IE0, IS0, canon%XMINC, canon%XMAXC, canon%ZMINC, canon%ZMAXC, &
   !                             DEBUG_OUTPUT=.FALSE.)

   !          CALL SYSTEM_CLOCK(iter_clock_end)
   !          iter_tick_delta = iter_clock_end - iter_clock_start
   !          IF (iter_tick_delta < 0) THEN
   !             IF (clock_max > 0) THEN
   !                iter_tick_delta = iter_tick_delta + clock_max
   !             ELSE
   !                iter_tick_delta = -iter_tick_delta
   !             END IF
   !          END IF
   !          iter_elapsed = DBLE(iter_tick_delta)/DBLE(rate)/60.D0

   !          IF (my_rank == 0) THEN
   !             DO I = 1, NPAR
   !                IF (INVP(I) == 1) THEN
   !                   WRITE (66, '(F7.2,"   ",I4,"   ",ES17.10,"   ",ES17.10,"   ",ES17.10,"   ", F12.3)') &
   !                      FREQ, ITER, FCOST_prev, GRAD_scaled_NORM(I), RESID_RMS, iter_elapsed
   !                END IF
   !             END DO
   !          END IF

   !          CALL check_convergence(ITER, MAXITER, NM, &
   !                                 m_coarse_prev, m_coarse, &
   !                                 FCOST_prev, FCOST0, &
   !                                 m_coarse0_norm, GRAD_scaled_norm_conv, grad0_norm, &
   !                                 cost_conv_tol, model_conv_tol, &
   !                                 MIN_ITER_CHECK, &
   !                                 my_rank, converged)

   !          IF (converged) THEN
   !             CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
   !             exit_iter_loop = .TRUE.
   !          END IF
   !          EXIT
   !       ELSE
   !          ! reject: rollback and reduce alpha
   !          CRR0 = CRR0_old
   !          CII0 = CII0_old
   !          m_fine = m_fine_old
   !          alpha = alpha*0.5_dp
   !          IF (my_rank == 0) WRITE (*, '(A)') 'TNGN: reject, reducing alpha'
   !       END IF
   !    END DO

   !    IF (.NOT. accepted_step) THEN
   !       IF (my_rank == 0) WRITE (*, '(A)') 'TNGN: all trial steps rejected, stopping iteration loop'
   !       exit_iter_loop = .TRUE.
   !    END IF

   !    CALL SYSTEM_CLOCK(iTimes2)
   !    IF (my_rank == 0) WRITE (*, '(A,1X,G0,1X,A)') 'Time for TNGN = ', REAL(iTimes2 - iTimes1)/REAL(rate)/60.D0, 'min'

   !    ! cleanup
   !    IF (ALLOCATED(dm)) DEALLOCATE (dm)
   !    IF (ALLOCATED(p_k_local)) DEALLOCATE (p_k_local)
   !    IF (ALLOCATED(p_kf_local)) DEALLOCATE (p_kf_local)
   !    IF (ALLOCATED(CRR0_old)) DEALLOCATE (CRR0_old)
   !    IF (ALLOCATED(CII0_old)) DEALLOCATE (CII0_old)
   !    IF (ALLOCATED(m_fine_old)) DEALLOCATE (m_fine_old)

   ! END SUBROUTINE TNGN_update

end program
