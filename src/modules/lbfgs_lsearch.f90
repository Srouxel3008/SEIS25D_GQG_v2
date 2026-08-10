module lbfgs_lsearch
   USE omp_lib
   use MPI
   use output_mod
   USE shared_mod
   use hardware_mod
   use F_modeling
   USE gradient_mod
   USE gridtype_mod
   USE misfit_mod
   USE constant_mod
   use iso_fortran_env, only: dp => real64, sp => real32
   implicit NONE
   !-----------------------------------------------------------------------
   !This module provides the line search routines for the L-BFGS optimization
   !algorithm, including backtracking line search and cubic interpolation.
   !Included subroutines
   ! BACKTRACK_LINESEARCH_MT : Maure Thuente line search for multiparameter
   !                             seismic inversion (Nocedal, 2006)
   ! BACKTRACK_LINESEARCH : Armijo only line search with quadratic interpolation
   ! ALPHA_INTERP : Quadratic interpolation to refine step length for backtracking
   !subfcostcal : Computes the cost function and residuals for linesearch
   ! UPDATEPARAMETERS : Updates model parameters from the search direction
   !polint: Quadratic interpolation subroutine
   !cubic_interp_fcost : Cubic interpolation for step length refinement for MT linesearch

CONTAINS

   !---------------------------------------------------------------

! !     BACKTRACK_LINESEARCH performs a backtracking line search to
! !     determine the step length (ALPHA) !     This routine iteratively tests a decreasing sequence of ALPHA values
! !     by updating the model parameters along the search direction `p_k`,
! !     and evaluating the cost function. If a sufficient decrease (Armijo)
! !     is not met, the step is reduced. Optionally, quadratic interpolation
! !     (alpha2 refinement) is used to refine the step.
! !
! !     Inputs:
! !       ALPHA................ Initial estimate for step length (modified in-place)
! !       FCOST0............... Cost function at current model
! !       p_k(NMM)............. Search direction vector (unscaled)
! !       CRR0(IANISO,NPT)....... Real-valued model parameters (elastic moduli)
! !       CII0(IANISO,NPT)....... Imaginary-valued model parameters (attenuation)
! !       m_fine(NMM)..... Previous model vector (unscaled)
! !       SCALER(NPAR)...... Per-parameter scaling factors (optional use)
! !       INVP(NPAR)........... Inversion mask (1=update, 0=fixed)
! !       X(NX)................ Horizontal grid coordinates
! !       XTO(NTO), ZTO(NTO)... Receiver coordinates
! !       XSR(NSR), ZSR(NSR)... Source coordinates
! !       YSR(NSR)............. Source Y-coordinates (for 2.5D)
! !       MSR(NSS), MSR1(NSS,*) Receiver map for each source
! !       FSR(NSS,*)........... Interpolation weights
! !       VSR(NSS,3,3)......... Source amplitude vector components
! !       GT0(ND).............. Observed data vector (complex)
! !       WAVELET(NK).......... Source wavelet (frequency domain)
! !       FREQ................. Current frequency
! !       FKY(NK), WTK(NK)..... Wavenumber and weights for y-direction integration
! !       AS(NORD), WT(NORD)... Quadrature abscissas and weights (x-z)
! !       DZ0.................  Nominal grid spacing in z-direction
! !       NX, NZ............... Number of primary grid points in X and Z
! !       NNX, NNZ............. Number of fine grid points in X and Z
! !       NPT.................. Total number of model grid points (NNX*NNZ)
! !       IANISO............... Number of elastic parameters
! !       NPAR................. Total number of model parameters (elastic + attenuation)
! !       NK................... Number of quadrature points in ky-direction
! !       NMM.................. Length of model vector (sum active params x NPT)
! !       ND................... Total number of data points
! !       NS(ND), NR(ND)....... Source and receiver index per data point
! !       NSV(ND), NRV(ND)..... Source/receiver vector index per data point
! !       IE0, IS0............. Absorbing layer size and surface index
! !       I25D................. Flag for 2.5D
! !       IFQ.................. Flag for Q attenuation
! !       COMP(ND)............. Component mask (1 to retain, 0 to ignore)
! !       ITER................. Outer iteration number (for logging)
! !
! !     Outputs:
! !       ALPHA................ Optimal step length (satisfying Armijo)
! !       FCOSTMIN............. Final cost value for selected ALPHA
! !
! !     Dependencies:
! !       - subfcostcal         (computes cost from forward modeling)
! !       - updateparameter     (builds CR1/CI1 from vector representation)
! !       - ALPHA_INTERP        (optional quadratic interpolation)
! !
! !-----------------------------------------------------------------------

   SUBROUTINE LINESEARCH_BACKTRACK(p_k, p_kf, GRAD_Scaled, m_fine, FCOST0, SCALER, F_SCALE, PAR_SCALE, INVP, &
                                   CRR0, CII0, ICSR, VSR, FSR, m_coarse, m_coarse_reg, REG_LAMBDA, &
                                   NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, X, DZ0, ALPHA00, &
                                   NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, &
                                   MSR, MSR1, COMP, FREQ, GT0, WAVELET, SourceScaler, WD_acq, WD_amp, IFQ, NK, FKY, WTK, AS, WT, &
                                   I25D, IS0, IE0, ITER, ALPHA, NBLOCK, IG, &
                                   my_rank, n_process, comm, IFLAG, FCOSTMIN, SOLVER_KIND, &
                                   DEBUG_OUTPUT)

      USE iso_fortran_env, ONLY: dp => real64, sp => real32
      IMPLICIT NONE

      !------------------------------
      ! Inputs
      !------------------------------
      INTEGER, INTENT(IN)  :: NPAR, ND, NSR, NSS, IE0, IFQ, NMM, NBLOCK, ITER, NM
      INTEGER, INTENT(IN)  :: I25D, IS0, NTO, NX, NNX, NZ, NNZ, NPT, IANISO, NORD, NK
      INTEGER, INTENT(IN)  :: comm, my_rank, n_process, SOLVER_KIND
      INTEGER, INTENT(IN)  :: NS(:), NR(:), NSV(:), NRV(:), INVP(:), COMP(:)
      INTEGER, INTENT(IN)  :: MSR(:), MSR1(:, :), ICSR(:)

      REAL(dp), INTENT(IN) :: CRR0(:, :), CII0(:, :), VSR(:, :, :), FSR(:, :)
      REAL(dp), INTENT(IN) :: m_coarse(:), m_coarse_reg(:), REG_LAMBDA
      REAL(dp), INTENT(IN) :: FREQ, DZ0, X(:), XTO(:), ZTO(:), XSR(:), ZSR(:), YSR(:)
      REAL(dp), INTENT(IN) :: AS(:), WT(:), FKY(:), WTK(:), WD_acq(:), WD_amp
      REAL(dp), INTENT(IN) :: p_k(:), p_kf(:), SCALER(:), m_fine(:), F_SCALE(:), PAR_SCALE(:)
      REAL(dp), INTENT(IN) :: GRAD_Scaled(:)
      REAL(dp), INTENT(IN) :: FCOST0, ALPHA00

      COMPLEX(dp), INTENT(IN) :: GT0(:), WAVELET(:), SourceScaler(:)
      TYPE(InversionGridType), INTENT(IN) :: IG
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT

      !------------------------------
      ! Outputs
      !------------------------------
      REAL(dp), INTENT(OUT) :: ALPHA, FCOSTMIN
      INTEGER, INTENT(OUT)  :: IFLAG

      !------------------------------
      ! Parameters
      !------------------------------
      REAL(dp), PARAMETER :: C1 = 1.0e-4_dp
      REAL(dp), PARAMETER :: TAU = 0.5_dp
      REAL(dp), PARAMETER :: MIN_STEP = 1.0e-5_dp
      REAL(dp), PARAMETER :: MAX_STEP = 10.0_dp
      INTEGER, PARAMETER :: MAX_IT = 15
      REAL(dp), PARAMETER :: TAU_COST = 1.0e-5_dp

      ! Relative-update cap parameters
      REAL(dp), PARAMETER :: eps_r = 1.0e-30_dp
      REAL(dp), PARAMETER :: r_cap = 15.0e-2_dp   ! % relative max update cap

      CHARACTER(len=80), PARAMETER :: FMT_ALPHA = &
                                      '(F7.2, "   ", I5, " ", I5, "  ", A, " ", ES13.6, "   ", 2(ES17.10, "   "))'

      !------------------------------
      ! Locals
      !------------------------------
      REAL(dp) :: grad_dot_p0, FCOST1, alpha_rhs
      REAL(dp) :: alpha_cap, alpha_trial, alpha_probe, rmax1
      REAL(dp), ALLOCATABLE :: V_STEP(:), PAR_STEP(:)
      REAL(dp), ALLOCATABLE :: CR1(:, :), CI1(:, :)
      COMPLEX(sp), ALLOCATABLE :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE :: RESID(:)

      INTEGER :: it, ierr, ios
      LOGICAL :: dbg, ok_eval, is_open
      CHARACTER(len=200) :: iomsg

      LOGICAL, SAVE :: ls_ready = .FALSE.
      INTEGER, SAVE :: ls_unit = 65

      !------------------------------
      ! Init
      !------------------------------
      dbg = .FALSE.
      IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

      IFLAG = 11
      ALPHA = 0.0_dp
      FCOSTMIN = FCOST0

      ALLOCATE (CR1(IANISO, NPT), CI1(IANISO, NPT), V_STEP(NMM), PAR_STEP(NMM), RESID(ND))
      ALLOCATE (GFX(NSS, 3, NK, NPT), GFY(NSS, 3, NK, NPT), GFZ(NSS, 3, NK, NPT))

      !------------------------------
      ! Log file open once
      !------------------------------
      IF (my_rank == 0 .AND. .NOT. ls_ready) THEN
         INQUIRE (UNIT=ls_unit, OPENED=is_open)
         IF (.NOT. is_open) THEN
            OPEN (UNIT=ls_unit, FILE='out_linesearch.txt', STATUS='UNKNOWN', POSITION='APPEND', &
                  ACTION='WRITE', IOSTAT=ios, IOMSG=iomsg)
            IF (ios == 0) THEN
               WRITE (ls_unit, '(A)') '# FREQ   ITER   it   TAG         ALPHA           FCOST1            FCOST0'
               FLUSH (ls_unit)
            ELSE
               WRITE (*, *) 'LINESEARCH_BACKTRACK: cannot open out_linesearch.txt: ', TRIM(iomsg)
            END IF
         END IF
         ls_ready = .TRUE.
      END IF

      !------------------------------
      ! Descent check at alpha = 0
      !------------------------------
      grad_dot_p0 = DOT_PRODUCT(GRAD_Scaled(1:NM), p_k(1:NM))

      IF (my_rank == 0 .AND. dbg) THEN
         WRITE (*, '(A,1PE12.4,2X,A,1PE12.4)') 'BT entry: phi(0)=', FCOST0, 'g·p(0)=', grad_dot_p0
      END IF

      IF (grad_dot_p0 >= 0.0_dp) THEN
         IF (my_rank == 0) WRITE (*, *) 'Backtrack: non-descent at base; abort.'
         GOTO 900
      END IF

      !------------------------------
      ! Initial trial alpha from caller
      !------------------------------
      alpha_trial = ALPHA00
      IF ((alpha_trial /= alpha_trial) .OR. (alpha_trial <= 0.0_dp)) alpha_trial = 1.0_dp
      alpha_trial = MAX(MIN_STEP, MIN(MAX_STEP, alpha_trial))
      alpha_cap = MAX_STEP

      !------------------------------
      ! Compute relative-update cap ONCE per call
      !------------------------------
      alpha_probe = 1.0_dp

      CALL UPDATEPARAMETER(alpha_probe, NMM, NPT, NPAR, IANISO, my_rank, &
                           p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                           CRR0, CII0, CR1, CI1, F_SCALE)

!       ! CALL compute_rmax_rel(CR1, CI1, CRR0, CII0, INVP, NPAR, IANISO, NPT, eps_r, rmax1)
      CALL compute_rmax_rel(CR1, CI1, CRR0, CII0, INVP, NPAR, IANISO, NPT, eps_r, rmax1, my_rank, dbg)
      alpha_cap = 1.0_dp
      IF (rmax1 > r_cap) alpha_cap = r_cap/rmax1
      IF ((alpha_cap /= alpha_cap) .OR. (alpha_cap <= 0.0_dp)) alpha_cap = MIN_STEP
      alpha_cap = MAX(MIN_STEP, MIN(MAX_STEP, alpha_cap))

      IF (my_rank == 0 .AND. dbg) THEN
         WRITE (*, '(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') &
            'BT cap probe: rmax1=', rmax1, 'alpha_cap=', alpha_cap, 'alpha00=', alpha_trial
      END IF

      ! Apply cap once to the initial guess
      alpha_trial = MIN(alpha_trial, alpha_cap)
      alpha_trial = MAX(MIN_STEP, MIN(MAX_STEP, alpha_trial))
! alpha_trial = 1.0_dp ! disable cap for now, as it can be too aggressive and cause failure to find a step in early iterations when the gradient is very steep. We keep the cap logic here for future use and testing.
      !------------------------------
      ! Backtracking loop
      !------------------------------
      DO it = 1, MAX_IT

         ALPHA = alpha_trial

         CALL UPDATEPARAMETER(ALPHA, NMM, NPT, NPAR, IANISO, my_rank, &
                              p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                              CRR0, CII0, CR1, CI1, F_SCALE)

         CALL MPI_BARRIER(comm, ierr)

         CALL subfcostcal(GFX, GFY, GFZ, I25D, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR1, CI1, NORD, &
                          NM, AS, WT, FREQ, NK, FKY, WTK, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                          IE0, IS0, DZ0, GT0, ND, NS, NR, NSV, NRV, YSR, NPT, WAVELET, SourceScaler, IFQ, WD_amp, WD_acq, &
                          NBLOCK, IG, FCOST1, RESID, m_coarse, m_coarse_reg, REG_LAMBDA, &
                          my_rank, n_process, comm, ITER, COMP, SOLVER_KIND, DEBUG_OUTPUT=dbg)

         ok_eval = .TRUE.
         IF (.NOT. (FCOST1 == FCOST1)) ok_eval = .FALSE.
         IF (.NOT. (ABS(FCOST1) < HUGE(1.0_dp))) ok_eval = .FALSE.

         alpha_rhs = FCOST0 + C1*ALPHA*grad_dot_p0

         IF (my_rank == 0 .AND. ls_ready) THEN
            WRITE (ls_unit, FMT_ALPHA) FREQ, ITER, it, 'ALPHA BT TRY', ALPHA, FCOST1, FCOST0
            FLUSH (ls_unit)
            IF (dbg) WRITE (*, '(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') 'BT: alpha=', ALPHA, 'phi=', FCOST1, 'rhs=', alpha_rhs
         END IF
!--------------------------------------------------------
! EARLY STOP: cost stagnation (avoid useless tiny alpha)
!--------------------------------------------------------
         IF (ok_eval .AND. FCOST1 > FCOST0) THEN
            IF ((FCOST1 - FCOST0)/MAX(FCOST0, tiny_dp) < TAU_COST) THEN

               IF (dbg .AND. my_rank == 0) THEN
                  WRITE (*, '(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') &
                     'BT STOP (stagnation): alpha=', ALPHA, &
                     'FCOST1=', FCOST1, 'FCOST0=', FCOST0
               END IF

               EXIT   ! exit backtracking loop early

            END IF
         END IF

         IF (ok_eval .AND. (FCOST1 <= alpha_rhs)) THEN
            FCOSTMIN = FCOST1
            IFLAG = 2
            GOTO 900
         END IF

         ! Reduce step, but never exceed cap
         alpha_trial = MAX(MIN_STEP, TAU*ALPHA)
         alpha_trial = MIN(alpha_trial, alpha_cap)

         ! No point continuing if already pinned at MIN_STEP
         IF (alpha_trial <= MIN_STEP) EXIT

      END DO

      !------------------------------
      ! Final tiny attempt at MIN_STEP
      !------------------------------
      ALPHA = MIN_STEP

      CALL UPDATEPARAMETER(ALPHA, NMM, NPT, NPAR, IANISO, my_rank, &
                           p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                           CRR0, CII0, CR1, CI1, F_SCALE)

      CALL subfcostcal(GFX, GFY, GFZ, I25D, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR1, CI1, NORD, &
                       NM, AS, WT, FREQ, NK, FKY, WTK, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                       IE0, IS0, DZ0, GT0, ND, NS, NR, NSV, NRV, YSR, NPT, WAVELET, SourceScaler, IFQ, WD_amp, WD_acq, &
                       NBLOCK, IG, FCOST1, RESID, m_coarse, m_coarse_reg, REG_LAMBDA, &
                       my_rank, n_process, comm, ITER, COMP, SOLVER_KIND, DEBUG_OUTPUT=dbg)

      alpha_rhs = FCOST0 + C1*ALPHA*grad_dot_p0

      IF ((FCOST1 == FCOST1) .AND. (ABS(FCOST1) < HUGE(1.0_dp)) .AND. (FCOST1 <= alpha_rhs)) THEN
         FCOSTMIN = FCOST1
         IFLAG = 2
      ELSE
         FCOSTMIN = FCOST0
         IFLAG = 11
         IF (my_rank == 0) WRITE (*, *) 'Backtrack: failed after MAX_IT'
      END IF

900   CONTINUE

      IF (ALLOCATED(V_STEP)) DEALLOCATE (V_STEP)
      IF (ALLOCATED(PAR_STEP)) DEALLOCATE (PAR_STEP)
      IF (ALLOCATED(CR1)) DEALLOCATE (CR1)
      IF (ALLOCATED(CI1)) DEALLOCATE (CI1)
      IF (ALLOCATED(RESID)) DEALLOCATE (RESID)
      IF (ALLOCATED(GFX)) DEALLOCATE (GFX)
      IF (ALLOCATED(GFY)) DEALLOCATE (GFY)
      IF (ALLOCATED(GFZ)) DEALLOCATE (GFZ)

      RETURN

   END SUBROUTINE LINESEARCH_BACKTRACK

SUBROUTINE LINESEARCH_MT(ALPHA00_prev, p_k, p_kf, GRAD_scaled, m_fine, SCALER, F_SCALE, PAR_SCALE, BALANCE_SCALE, INVP, CRR0, CII0, ICSR, VSR, FSR,  m_coarse, m_coarse_reg,  REG_LAMBDA,&
                            NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, X, DZ0, &
                            NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, MSR, MSR1, COMP, &
                            FREQ, GT0, WAVELET, SourceScaler, WD_amp, WD_acq, IFQ, NK, FKY, WTK, AS, WT, I25D, IS0, IE0, ITER, &
                   FRECHET, GMASK, GMASK_FINAL, GRAD_scaled_NORM, GRADr_NORM, CREF, HESSDI_scaled, GRADr, GRAD_smooth, sigma_grad_x, sigma_grad_z, &
                            GRADsc_precond, GRADsc_Precond_NORM, USE_PRECOND, &
                            NBLOCK, IG, USE_GR_SMOOTH, &
                            ALPHA, FCOST0, FCOSTMIN, IFLAG, SOLVER_KIND, my_rank, n_process, comm, &
                            PARAM, ITHOM, IVISCO, NCOMP, NCOMPS, FCOST0_ref, GRAD_scaled_ref, p_k_ref, dot_ref, &
                            XMINC, XMAXC, ZMINC, ZMAXC, DEBUG_OUTPUT)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! ---------- INTENTs ----------
      INTEGER, INTENT(IN) :: NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, ND
      INTEGER, INTENT(IN) :: NTO, NK, IFQ, IE0, IS0, ITER, I25D, NBLOCK, SOLVER_KIND
      INTEGER, INTENT(IN) :: NSR, NSS, comm, my_rank, n_process
      INTEGER, INTENT(IN) :: ITHOM, IVISCO, NCOMP, NCOMPS
      INTEGER, INTENT(IN) :: NS(:), NR(:), NSV(:), NRV(:), INVP(:), COMP(:), MSR(:), MSR1(:, :), ICSR(:)

      REAL(dp), INTENT(IN) :: FREQ, DZ0, X(:), XTO(:), ZTO(:), XSR(:), ZSR(:), YSR(:)
      REAL(dp), INTENT(IN) :: AS(:), WT(:), FKY(:), WTK(:), CREF, F_SCALE(:), PAR_SCALE(:), BALANCE_SCALE(:)
      REAL(dp), INTENT(IN) :: m_coarse(:), m_coarse_reg(:), REG_LAMBDA
      REAL(dp), INTENT(IN) :: XMINC, XMAXC, ZMINC, ZMAXC
      REAL(dp), INTENT(IN) :: p_k(:), p_kf(:), SCALER(:), m_fine(:), FCOST0, ALPHA00_prev
      REAL(dp), INTENT(IN) :: CRR0(:, :), CII0(:, :), VSR(:, :, :), FSR(:, :)
      REAL(dp), INTENT(IN) :: GMASK(:), GMASK_FINAL(:), HESSDI_scaled(:)
      REAL(dp), INTENT(INOUT) :: GRADr(:), GRAD_scaled(:), GRAD_smooth(:)
      REAL(dp), INTENT(INOUT) :: GRADsc_precond(:)
      REAL(dp), INTENT(INOUT) :: GRADr_NORM(:), GRAD_scaled_NORM(:), GRADsc_Precond_NORM(:)
      REAL(dp), INTENT(IN)    :: WD_amp, WD_acq(:)
      COMPLEX(dp), INTENT(IN) :: GT0(:), WAVELET(:), SourceScaler(:)
      CHARACTER(LEN=5), INTENT(IN) :: PARAM(22)
      REAL(dp), INTENT(IN) :: FCOST0_ref, dot_ref
      REAL(dp), INTENT(IN) :: p_k_ref(NM), GRAD_scaled_ref(NM)

      TYPE(InversionGridType), INTENT(IN) :: IG
      LOGICAL, INTENT(IN) :: USE_PRECOND, USE_GR_SMOOTH
      REAL(dp), INTENT(IN), OPTIONAL :: sigma_grad_x, sigma_grad_z
      COMPLEX(dp), INTENT(INOUT) :: FRECHET(ND, NM)

      REAL(dp), INTENT(OUT) :: ALPHA, FCOSTMIN
      INTEGER, INTENT(OUT)  :: IFLAG
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT

      ! ---------- params (More-Thuente via dcsrch) ----------
      REAL(dp), PARAMETER :: C1 = 1.0e-4_dp
      REAL(dp), PARAMETER :: C2 = 0.9_dp
      REAL(dp), PARAMETER :: MIN_STEP = 1.0e-5_dp
      REAL(dp), PARAMETER :: MAX_STEP = 10.0_dp
      INTEGER, PARAMETER :: MAX_ITER = 10
      REAL(dp), PARAMETER :: XTOL = 1.0e-12_dp

      ! Relative-update cap parameters
      REAL(dp), PARAMETER :: TAU_COST = 1.0e-4_dp
      REAL(dp), PARAMETER :: eps_r = 1.0e-30_dp
      REAL(dp), PARAMETER :: r_cap = 15.0e-2_dp   ! % relative max update cap

      CHARACTER(len=80), PARAMETER :: FMT_ALPHA = &
                                      '(F7.2, "   ", I5, " ", I5, "  ", A, " ", ES13.6, "   ", 2(ES17.10, "   "))'

      ! ---------- locals ----------
      REAL(dp), ALLOCATABLE :: V_STEP(:), PAR_STEP(:), CR1(:, :), CI1(:, :)
      COMPLEX(dp), ALLOCATABLE :: RESID(:)
      COMPLEX(sp), ALLOCATABLE :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)

      REAL(dp) :: FCOST1, grad_dot_p0, grad_dot_p, armijo_rhs
      REAL(dp) :: stp, f_ls, g_ls, stpmax_cap
      REAL(dp) :: alpha_cap, alpha_probe, rmax1
      INTEGER  :: ITER_LS, ierr, NSSt
      LOGICAL  :: dbg, FOUND, use_fallback, hit_minstep, armijo_ok, curvature_ok
      CHARACTER(len=60) :: task
      INTEGER :: isave(2)
      REAL(dp) :: dsave(13)
      INTEGER :: clk_start, clk_end, clk_rate
      REAL(dp) :: ls_iter_time

      LOGICAL, SAVE :: ls_ready = .FALSE.
      INTEGER, SAVE  :: ls_unit = 65
      LOGICAL  :: is_open
      INTEGER  :: ios
      CHARACTER(len=200) :: iomsg

      dbg = .FALSE.
      IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      NSSt = NSR
      FOUND = .FALSE.
      use_fallback = .FALSE.
      hit_minstep = .FALSE.

      IFLAG = 11
      ALPHA = 0.0_dp
      FCOSTMIN = FCOST0

      ! ---------- allocs ----------
      ALLOCATE (CR1(IANISO, NPT), CI1(IANISO, NPT), V_STEP(NMM), PAR_STEP(NMM), RESID(ND))
      ALLOCATE (GFX(NSSt, 3, NK, NPT), GFY(NSSt, 3, NK, NPT), GFZ(NSSt, 3, NK, NPT))

      ! ---------- log file ----------
      IF (my_rank == 0 .AND. .NOT. ls_ready) THEN
         INQUIRE (UNIT=ls_unit, OPENED=is_open)
         IF (.NOT. is_open) THEN
            OPEN (UNIT=ls_unit, FILE='out_linesearch.txt', STATUS='UNKNOWN', POSITION='APPEND', &
                  ACTION='WRITE', IOSTAT=ios, IOMSG=iomsg)
            IF (ios == 0) THEN
               WRITE (ls_unit, '(A)') '# FREQ   ITER   it   TAG         ALPHA           FCOST1            FCOST0'
               FLUSH (ls_unit)
            ELSE IF (dbg .AND. my_rank == 0) THEN
               WRITE (*, *) 'LINESEARCH_MT: cannot open out_linesearch.txt: ', TRIM(iomsg)
            END IF
         END IF
         ls_ready = .TRUE.
      END IF

      ! Keep start checks unchanged.
      IF (FCOST0 <= 0.0_dp) THEN
         IFLAG = 11
         ALPHA = 0.0_dp
         GOTO 900
      END IF

      grad_dot_p0 = DOT_PRODUCT(GRAD_scaled_ref(1:NM), p_k(1:NM))
      IF (grad_dot_p0 >= 0.0_dp) THEN
         IFLAG = 11
         ALPHA = 0.0_dp
         GOTO 900
      END IF

      stp = ALPHA00_prev

      stp = MAX(MIN_STEP, MIN(MAX_STEP, stp))

      IF (my_rank == 0) THEN
         WRITE (*, '(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') &
            'MT entry: phi(0)=', FCOST0, 'dg0=', grad_dot_p0, 'alpha0=', stp
      END IF

!       ! ------------------------------
!       ! Compute relative-update cap ONCE per call (matching BT logic)
!       ! ------------------------------
!       alpha_probe = 1.0_dp

!       CALL UPDATEPARAMETER(alpha_probe, NMM, NPT, NPAR, IANISO, my_rank, &
!                            p_kf, m_fine, PAR_STEP, V_STEP, INVP, SCALER, &
!                            CRR0, CII0, CR1, CI1,F_SCALE)

!       ! CALL compute_rmax_rel(CR1, CI1, CRR0, CII0, INVP, NPAR, IANISO, NPT, eps_r, rmax1)
! CALL compute_rmax_rel(CR1, CI1, CRR0, CII0, INVP, NPAR, IANISO, NPT, eps_r, rmax1, my_rank, dbg)
!       alpha_cap = 1.0_dp
!       IF (rmax1 > r_cap) alpha_cap = r_cap / rmax1
!       IF ((alpha_cap /= alpha_cap) .OR. (alpha_cap <= 0.0_dp)) alpha_cap = MIN_STEP
!       alpha_cap = MAX(MIN_STEP, MIN(MAX_STEP, alpha_cap))

!       IF (my_rank == 0 .AND. dbg) THEN
!          WRITE(*,'(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') &
!             'MT cap probe: rmax1=', rmax1, 'alpha_cap=', alpha_cap, 'alpha00=', stp
!       END IF

!       stp = MIN(stp, alpha_cap)
!       stp = MAX(MIN_STEP, MIN(MAX_STEP, stp))

!       stpmax_cap = MIN(MAX_STEP, alpha_cap)
      stpmax_cap = stp

      f_ls = FCOST0
      g_ls = grad_dot_p0
      task = 'START'
      ITER_LS = 0

      DO
         CALL dcsrch(f_ls, g_ls, stp, C1, C2, XTOL, MIN_STEP, stpmax_cap, task, isave, dsave)

         IF (task(1:2) == 'FG') THEN
            ITER_LS = ITER_LS + 1
            IF (ITER_LS > MAX_ITER) THEN
               task = 'WARNING: MAX ITER REACHED'
               use_fallback = .TRUE.
               EXIT
            END IF

            CALL SYSTEM_CLOCK(clk_start, clk_rate)

            CALL UPDATEPARAMETER(stp, NMM, NPT, NPAR, IANISO, my_rank, &
                                 p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                                 CRR0, CII0, CR1, CI1, F_SCALE)

            CALL MPI_Barrier(comm, ierr)

            CALL subfcostcal(GFX, GFY, GFZ, I25D, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR1, CI1, NORD, &
                             NM, AS, WT, FREQ, NK, FKY, WTK, NSR, NSSt, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                             IE0, IS0, DZ0, GT0, ND, NS, NR, NSV, NRV, YSR, NPT, WAVELET, SourceScaler, IFQ, &
                             WD_amp, WD_acq, NBLOCK, IG, FCOST1, RESID, m_coarse, m_coarse_reg, REG_LAMBDA, &
                             my_rank, n_process, comm, ITER, COMP, SOLVER_KIND, DEBUG_OUTPUT=dbg)

            IF (.NOT. (FCOST1 == FCOST1) .OR. ABS(FCOST1) >= HUGE(1.0_dp)) FCOST1 = HUGE(1.0_dp)

            IF (my_rank == 0 .AND. ls_ready) THEN
               WRITE (ls_unit, FMT_ALPHA) FREQ, ITER, ITER_LS, 'ALPHA TRY  ', stp, FCOST1, FCOST0
               FLUSH (ls_unit)
            END IF

!--------------------------------------------------------
! EARLY STOP: cost stagnation (avoid useless tiny alpha)
!--------------------------------------------------------
            IF (FCOST1 > FCOST0) THEN
               IF ((FCOST1 - FCOST0)/MAX(FCOST0, tiny_dp) < TAU_COST) THEN

                  IF (dbg .AND. my_rank == 0) THEN
                     WRITE (*, '(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') &
                        'LS STOP (stagnation): alpha=', stp, &
                        'FCOST1=', FCOST1, 'FCOST0=', FCOST0
                  END IF

                  use_fallback = .TRUE.
                  EXIT   ! exit dcsrch loop early

               END IF
            END IF
            CALL QFRECHET(ITER, PARAM, IFQ, FREQ, NSR, ND, NX, NZ, NNX, NNZ, X, IE0, DZ0, NPT, NBLOCK, &
                          IANISO, ITHOM, IVISCO, NPAR, INVP, IS0, CR1, CI1, &
                          NS, NSV, NR, NRV, NTO, XTO, ZTO, YSR, VSR, &
                          AS, WT, NK, FKY, WTK, GFX, GFY, GFZ, &
                          NM, FRECHET, WAVELET, SourceScaler, NORD, NCOMP, NCOMPS, &
                          IG, my_rank, n_process, comm, DEBUG_OUTPUT=.FALSE.)

            CALL compute_gradient_scaled(FRECHET, RESID, NPAR, F_SCALE, PAR_SCALE, BALANCE_SCALE, GMASK, ND, NM, &
                                         m_coarse, m_coarse_reg, REG_LAMBDA, &
                                         GRADr, GRAD_scaled, GRADr_NORM, GRAD_scaled_NORM, &
                                         NBLOCK, INVP, USE_PRECOND=USE_PRECOND, my_rank=my_rank, &
                                         HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                                         GRADsc_Precond_NORM=GRADsc_Precond_NORM, DEBUG_OUTPUT=.FALSE., &
                                         PARAM=PARAM, FREQ=FREQ, ITER=ITER, &
                                         NX=NX, NZ=NZ, IG=IG, NTO=NTO, XTO=XTO, ZTO=ZTO, IE0=IE0, IS0=IS0, &
                                         XMINC=XMINC, XMAXC=XMAXC, ZMINC=ZMINC, ZMAXC=ZMAXC)

            IF (USE_GR_SMOOTH) THEN
               CALL smooth_gradient_blocks( &
                  GRAD_in=GRAD_scaled, GMASK=GMASK_FINAL, &
                  NBLOCK=NBLOCK, INVP=INVP, NPAR=NPAR, &
                  NX=NX, NZ=NZ, IG=IG, &
                  FREQ=FREQ, CREF=CREF, my_rank=my_rank, &
                  USE_GR_SMOOTH=.TRUE., USE_PRECOND=USE_PRECOND, &
                  HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                  GRADsc_Precond_NORM=GRADsc_Precond_NORM, &
                  PARAM=Param, ITER=iter, NTO=nto, XTO=xto, ZTO=zto, IE0=IE0, IS0=IS0, &
                  XMINC=XMINC, XMAXC=XMAXC, ZMINC=ZMINC, ZMAXC=ZMAXC, &
                  DEBUG_OUTPUT=.FALSE., &
                  SIGMA_X=sigma_grad_x, SIGMA_Z=sigma_grad_z, &
                  GRAD_out=GRAD_smooth) !
               GRAD_scaled(1:NM) = GRAD_smooth(1:NM)
               CALL recompute_masked_norms(GRAD_scaled, GRAD_scaled_NORM, INVP, NBLOCK)
            END IF

            grad_dot_p = DOT_PRODUCT(GRAD_scaled(1:NM), p_k(1:NM))
            armijo_rhs = FCOST0 + C1*stp*grad_dot_p0
            armijo_ok = (FCOST1 <= armijo_rhs)
            curvature_ok = (ABS(grad_dot_p) <= C2*ABS(grad_dot_p0))

            f_ls = FCOST1
            g_ls = grad_dot_p
            IF (.NOT. (g_ls == g_ls) .OR. ABS(g_ls) >= HUGE(1.0_dp)) g_ls = 0.0_dp

            IF (my_rank == 0) THEN
               WRITE (*, '(A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4,2X,A,1PE12.4)') &
                  'MT try: alpha=', stp, 'phi=', FCOST1, 'dg=', grad_dot_p, 'dg0=', grad_dot_p0
               WRITE (*, '(A,1PE12.4,2X,A,L1,2X,A,L1,2X,2A)') &
                  'MT try: armijo_rhs=', armijo_rhs, 'armijo=', armijo_ok, 'curvature=', curvature_ok, 'task=', TRIM(task)
               CALL SYSTEM_CLOCK(clk_end)
               ls_iter_time = DBLE(clk_end - clk_start)/DBLE(clk_rate)
               WRITE (*, '(A,I3,2X,A,F8.3,A)') 'MT linesearch iter=', ITER_LS, ' time=', ls_iter_time, ' s'
            END IF

            CYCLE
         ELSE IF (task(1:4) == 'CONV') THEN
            FOUND = .TRUE.
            IFLAG = 2
            ALPHA = stp
            FCOSTMIN = f_ls
            EXIT
         ELSE IF (task(1:4) == 'WARN') THEN
            use_fallback = .TRUE.
            IF ((task(1:21) == 'WARNING: STP = STPMIN') .OR. &
                (stp <= MIN_STEP*(1.0_dp + 1.0e-12_dp))) hit_minstep = .TRUE.
            EXIT
         ELSE IF (task(1:5) == 'ERROR') THEN
            IFLAG = 11
            ALPHA = 0.0_dp
            FCOSTMIN = FCOST0
            use_fallback = .TRUE.
            EXIT
         ELSE
            use_fallback = .TRUE.
            EXIT
         END IF
      END DO

      IF (.NOT. FOUND) THEN
         IF (use_fallback .and. ITER < 5) THEN

            IF (dbg .AND. my_rank == 0) THEN
               WRITE (*, *) 'LINESEARCH_MT(dcsrch): fallback to backtrack, task=', TRIM(task), &
                  ' stp=', stp, ' hit_minstep=', hit_minstep
            END IF
            CALL LINESEARCH_BACKTRACK(p_k, p_kf, GRAD_scaled_ref, m_fine, FCOST0, SCALER, F_SCALE, PAR_SCALE, INVP, &
                                      CRR0, CII0, ICSR, VSR, FSR, m_coarse, m_coarse_reg, REG_LAMBDA, &
                                      NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, X, DZ0, 0.0_dp, &
                                      NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, &
                                     MSR, MSR1, COMP, FREQ, GT0, WAVELET, SourceScaler, WD_acq, WD_amp, IFQ, NK, FKY, WTK, AS, WT, &
                                      I25D, IS0, IE0, ITER, ALPHA, NBLOCK, IG, &
                                      my_rank, n_process, comm, IFLAG, FCOSTMIN, SOLVER_KIND, DEBUG_OUTPUT=dbg)
         ELSE
            IFLAG = 11
            ALPHA = 0.0_dp
            FCOSTMIN = FCOST0
         END IF
      END IF

900   CONTINUE
      IF (ALLOCATED(V_STEP)) DEALLOCATE (V_STEP)
      IF (ALLOCATED(PAR_STEP)) DEALLOCATE (PAR_STEP)
      IF (ALLOCATED(CR1)) DEALLOCATE (CR1)
      IF (ALLOCATED(CI1)) DEALLOCATE (CI1)
      IF (ALLOCATED(RESID)) DEALLOCATE (RESID)
      IF (ALLOCATED(GFX)) DEALLOCATE (GFX)
      IF (ALLOCATED(GFY)) DEALLOCATE (GFY)
      IF (ALLOCATED(GFZ)) DEALLOCATE (GFZ)

   END SUBROUTINE LINESEARCH_MT

!*******************************************************************************
!>
!  This subroutine computes a safeguarded step for a search
!  procedure and updates an interval that contains a step that
!  satisfies a sufficient decrease and a curvature condition.
!
!  The parameter `stx` contains the step with the least function
!  value. If `brackt` is set to .true. then a minimizer has
!  been bracketed in an interval with endpoints `stx` and `sty`.
!  The parameter `stp` contains the current step.
!  The subroutine assumes that if `brackt` is set to .true. then
!
!  `min(stx,sty) < stp < max(stx,sty)`
!
!  and that the derivative at `stx` is negative in the direction
!  of the step.
!
!### Credits
!###github
!  * MINPACK-1 Project. June 1983
!    Argonne National Laboratory.
!    Jorge J. More' and David J. Thuente.
!  * MINPACK-2 Project. October 1993.
!    Argonne National Laboratory and University of Minnesota.
!    Brett M. Averick and Jorge J. More'.

   subroutine dcstep(Stx, Fx, Dx, Sty, Fy, Dy, Stp, Fp, Dpval, Brackt, Stpmin, &
                     Stpmax)
      use iso_fortran_env, only: dp => real64
      implicit none

      logical, intent(inout) :: Brackt  !! On entry `brackt` specifies if a minimizer has been bracketed.
                                       !! Initially `brackt` must be set to .false.
                                       !! On exit `brackt` specifies if a minimizer has been bracketed.
                                       !! When a minimizer is bracketed `brackt` is set to .true.
      real(dp), intent(inout) :: Stx !! On entry `stx` is the best step obtained so far and is an
                                    !! endpoint of the interval that contains the minimizer.
                                    !! On exit `stx is the updated best step.
      real(dp), intent(inout) :: Fx !! On entry `fx` is the function at `stx`.
                                   !! On exit `fx` is the function at `stx`.
      real(dp), intent(inout) :: Dx !! On entry `dx` is the derivative of the function at
                                   !! `stx`. The derivative must be negative in the direction of
                                   !! the step, that is, `dx` and `stp - stx` must have opposite
                                   !! signs.
                                   !! On exit `dx` is the derivative of the function at `stx`.
      real(dp), intent(inout) :: Sty !! On entry `sty` is the second endpoint of the interval that contains the minimizer.
                                    !! On exit `sty` is the updated endpoint of the interval that contains the minimizer.
      real(dp), intent(inout) :: Fy !! On entry `fy` is the function at `sty`.
                                   !! On exit `fy` is the function at `sty`.
      real(dp), intent(inout) :: Dy !! On entry `dy` is the derivative of the function at `sty`.
                                   !! On exit `dy` is the derivative of the function at the exit `sty`.
      real(dp), intent(inout) :: Stp !! On entry `stp` is the current step. If `brackt` is set to .true.
                                    !! then on input `stp` must be between `stx` and `sty`.
                                    !! On exit `stp` is a new trial step.
      real(dp), intent(in) :: Fp !! the function at `stp`.
      real(dp), intent(in) :: Dpval !! the derivative of the function at `stp`.
      real(dp), intent(in) :: Stpmin !! a lower bound for the step.
      real(dp), intent(in) :: Stpmax !! an upper bound for the step.

      real(dp), parameter :: p66 = 0.66_dp

      real(dp) :: gamma, p, q, r, s, sgnd, stpc, stpf, &
                  stpq, theta

      sgnd = Dpval*(Dx/abs(Dx))

      ! First case: A higher function value. The minimum is bracketed.
      ! If the cubic step is closer to stx than the quadratic step, the
      ! cubic step is taken, otherwise the average of the cubic and
      ! quadratic steps is taken.

      if (Fp > Fx) then
         theta = 3.0_dp*(Fx - Fp)/(Stp - Stx) + Dx + Dpval
         s = max(abs(theta), abs(Dx), abs(Dpval))
         gamma = s*sqrt((theta/s)**2 - (Dx/s)*(Dpval/s))
         if (Stp < Stx) gamma = -gamma
         p = (gamma - Dx) + theta
         q = ((gamma - Dx) + gamma) + Dpval
         r = p/q
         stpc = Stx + r*(Stp - Stx)
         stpq = Stx + ((Dx/((Fx - Fp)/(Stp - Stx) + Dx))/2.0_dp)*(Stp - Stx)
         if (abs(stpc - Stx) < abs(stpq - Stx)) then
            stpf = stpc
         else
            stpf = stpc + (stpq - stpc)/2.0_dp
         end if
         Brackt = .true.

         ! Second case: A lower function value and derivatives of opposite
         ! sign. The minimum is bracketed. If the cubic step is farther from
         ! stp than the secant step, the cubic step is taken, otherwise the
         ! secant step is taken.

      else if (sgnd < 0.0_dp) then
         theta = 3.0_dp*(Fx - Fp)/(Stp - Stx) + Dx + Dpval
         s = max(abs(theta), abs(Dx), abs(Dpval))
         gamma = s*sqrt((theta/s)**2 - (Dx/s)*(Dpval/s))
         if (Stp > Stx) gamma = -gamma
         p = (gamma - Dpval) + theta
         q = ((gamma - Dpval) + gamma) + Dx
         r = p/q
         stpc = Stp + r*(Stx - Stp)
         stpq = Stp + (Dpval/(Dpval - Dx))*(Stx - Stp)
         if (abs(stpc - Stp) > abs(stpq - Stp)) then
            stpf = stpc
         else
            stpf = stpq
         end if
         Brackt = .true.

         ! Third case: A lower function value, derivatives of the same sign,
         ! and the magnitude of the derivative decreases.

      else if (abs(Dpval) < abs(Dx)) then

         ! The cubic step is computed only if the cubic tends to infinity
         ! in the direction of the step or if the minimum of the cubic
         ! is beyond stp. Otherwise the cubic step is defined to be the
         ! secant step.

         theta = 3.0_dp*(Fx - Fp)/(Stp - Stx) + Dx + Dpval
         s = max(abs(theta), abs(Dx), abs(Dpval))

         ! The case gamma = 0 only arises if the cubic does not tend
         ! to infinity in the direction of the step.

         gamma = s*sqrt(max(0.0_dp, (theta/s)**2 - (Dx/s)*(Dpval/s)))
         if (Stp > Stx) gamma = -gamma
         p = (gamma - Dpval) + theta
         q = (gamma + (Dx - Dpval)) + gamma
         r = p/q
         if (r < 0.0_dp .and. gamma /= 0.0_dp) then
            stpc = Stp + r*(Stx - Stp)
         else if (Stp > Stx) then
            stpc = Stpmax
         else
            stpc = Stpmin
         end if
         stpq = Stp + (Dpval/(Dpval - Dx))*(Stx - Stp)

         if (Brackt) then

            ! A minimizer has been bracketed. If the cubic step is
            ! closer to stp than the secant step, the cubic step is
            ! taken, otherwise the secant step is taken.

            if (abs(stpc - Stp) < abs(stpq - Stp)) then
               stpf = stpc
            else
               stpf = stpq
            end if
            if (Stp > Stx) then
               stpf = min(Stp + p66*(Sty - Stp), stpf)
            else
               stpf = max(Stp + p66*(Sty - Stp), stpf)
            end if
         else

            ! A minimizer has not been bracketed. If the cubic step is
            ! farther from stp than the secant step, the cubic step is
            ! taken, otherwise the secant step is taken.

            if (abs(stpc - Stp) > abs(stpq - Stp)) then
               stpf = stpc
            else
               stpf = stpq
            end if
            stpf = min(Stpmax, stpf)
            stpf = max(Stpmin, stpf)
         end if

         ! Fourth case: A lower function value, derivatives of the same sign,
         ! and the magnitude of the derivative does not decrease. If the
         ! minimum is not bracketed, the step is either stpmin or stpmax,
         ! otherwise the cubic step is taken.

      else
         if (Brackt) then
            theta = 3.0_dp*(Fp - Fy)/(Sty - Stp) + Dy + Dpval
            s = max(abs(theta), abs(Dy), abs(Dpval))
            gamma = s*sqrt((theta/s)**2 - (Dy/s)*(Dpval/s))
            if (Stp > Sty) gamma = -gamma
            p = (gamma - Dpval) + theta
            q = ((gamma - Dpval) + gamma) + Dy
            r = p/q
            stpc = Stp + r*(Sty - Stp)
            stpf = stpc
         else if (Stp > Stx) then
            stpf = Stpmax
         else
            stpf = Stpmin
         end if
      end if

      ! Update the interval which contains a minimizer.

      if (Fp > Fx) then
         Sty = Stp
         Fy = Fp
         Dy = Dpval
      else
         if (sgnd < 0.0_dp) then
            Sty = Stx
            Fy = Fx
            Dy = Dx
         end if
         Stx = Stp
         Fx = Fp
         Dx = Dpval
      end if

      ! Compute the new step.

      Stp = stpf

   end subroutine dcstep
!>
!  This subroutine finds a step that satisfies a sufficient
!  decrease condition and a curvature condition.
!
!  Each call of the subroutine updates an interval with
!  endpoints stx and sty. The interval is initially chosen
!  so that it contains a minimizer of the modified function
!
!  `psi(stp) = f(stp) - f(0) - ftol*stp*f'(0)`.
!
!  If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
!  interval is chosen so that it contains a minimizer of f.
!
!  The algorithm is designed to find a step that satisfies
!  the sufficient decrease condition
!
!  `f(stp) <= f(0) + ftol*stp*f'(0)`,
!
!  and the curvature condition
!
!  `abs(f'(stp)) <= gtol*abs(f'(0))`.
!
!  If ftol is less than gtol and if, for example, the function
!  is bounded below, then there is always a step which satisfies
!  both conditions.
!
!  If no step can be found that satisfies both conditions, then
!  the algorithm stops with a warning. In this case stp only
!  satisfies the sufficient decrease condition.
!
!  A typical invocation of dcsrch has the following outline:
!
!```fortran
!     task = 'START'
!     main : block
!       call dcsrch( ... )
!       if (task == 'FG') then
!          ! Evaluate the function and the gradient at stp
!          cycle main
!       end if
!```
!
!  Note: The user must not alter work arrays between calls.
!
!### Credits
!
!  * MINPACK-1 Project. June 1983.
!    Argonne National Laboratory.
!    Jorge J. More' and David J. Thuente.
!  * MINPACK-2 Project. October 1993.
!    Argonne National Laboratory and University of Minnesota.
!    Brett M. Averick, Richard G. Carter, and Jorge J. More'.

   subroutine dcsrch(f, g, Stp, Ftol, Gtol, Xtol, Stpmin, Stpmax, Task, Isave, Dsave)
      use iso_fortran_env, only: dp => real64
      implicit none

      character(len=*), intent(inout) :: Task !! `task` is a character variable of length at least 60:
                                             !!
                                             !!  * On initial entry `task` must be set to 'START'.
                                             !!  * On exit `task` indicates the required action:
                                             !!     * If `task(1:2) = 'FG'` then evaluate the function and
                                             !!       derivative at stp and call dcsrch again.
                                             !!     * If `task(1:4) = 'CONV'` then the search is successful.
                                             !!     * If `task(1:4) = 'WARN'` then the subroutine is not able
                                             !!       to satisfy the convergence conditions. The exit value of
                                             !!       `stp` contains the best point found during the search.
                                             !!     * If `task(1:5) = 'ERROR'` then there is an error in the
                                             !!       input arguments.
                                             !!  * On exit with convergence, a warning or an error, the
                                             !!    variable task contains additional information.
      real(dp), intent(inout) :: f !! * On initial entry `f` is the value of the function at 0.
                                  !!   On subsequent entries `f` is the value of the
                                  !!   function at `stp`.
                                  !! * On exit `f` is the value of the function at `stp`.
      real(dp), intent(inout) :: g !! * On initial entry `g` is the derivative of the function at 0.
                                  !!   On subsequent entries `g` is the  of the
                                  !!   function at `stp`.
                                  !! * On exit `g` is the derivative of the function at `stp`.
      real(dp), intent(inout) :: Stp !! * On entry `stp` is the current estimate of a satisfactory
                                    !!   step. On initial entry, a positive initial estimate
                                    !!   must be provided.
                                    !! * On exit `stp` is the current estimate of a satisfactory step
                                    !!   if `task = 'FG'`. If `task = 'CONV'` then `stp` satisfies
                                    !!   the sufficient decrease and curvature condition.
      real(dp), intent(in) :: Ftol !! `ftol` specifies a nonnegative tolerance for the
                                  !! sufficient decrease condition.
      real(dp), intent(in) :: Gtol !! `gtol` specifies a nonnegative tolerance for the curvature condition.
      real(dp), intent(in) :: Xtol !! `xtol` specifies a nonnegative relative tolerance
                                  !! for an acceptable step. The subroutine exits with a
                                  !! warning if the relative difference between `sty` and `stx`
                                  !! is less than `xtol`.
      real(dp), intent(in) :: Stpmin !! a nonnegative lower bound for the step.
      real(dp), intent(in) :: Stpmax !! a nonnegative upper bound for the step.
      integer :: Isave(2) !! integer work array
      real(dp) :: Dsave(13) !! real work array

      real(dp), parameter :: p5 = 0.5_dp
      real(dp), parameter :: p66 = 0.66_dp
      real(dp), parameter :: xtrapl = 1.1_dp
      real(dp), parameter :: xtrapu = 4.0_dp

      logical :: brackt
      integer :: stage
      real(dp) :: finit, ftest, fm, fx, fxm, fy, fym, &
                  ginit, gtest, gm, gx, gxm, gy, gym, stx, &
                  sty, stmin, stmax, width, width1

      ! Initialization block.

      if (Task(1:5) == 'START') then

         ! Check the input arguments for errors.

         if (Stp < Stpmin) Task = 'ERROR: STP < STPMIN'
         if (Stp > Stpmax) Task = 'ERROR: STP > STPMAX'
         if (g >= 0.0_dp) Task = 'ERROR: INITIAL G >= 0.0_dp'
         if (Ftol < 0.0_dp) Task = 'ERROR: FTOL < 0.0_dp'
         if (Gtol < 0.0_dp) Task = 'ERROR: GTOL < 0.0_dp'
         if (Xtol < 0.0_dp) Task = 'ERROR: XTOL < 0.0_dp'
         if (Stpmin < 0.0_dp) Task = 'ERROR: STPMIN < 0.0_dp'
         if (Stpmax < Stpmin) Task = 'ERROR: STPMAX < STPMIN'

         ! Exit if there are errors on input.

         if (Task(1:5) == 'ERROR') return

         ! Initialize local variables.

         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = Ftol*ginit
         width = Stpmax - Stpmin
         width1 = width/p5

         ! The variables stx, fx, gx contain the values of the step,
         ! function, and derivative at the best step.
         ! The variables sty, fy, gy contain the value of the step,
         ! function, and derivative at sty.
         ! The variables stp, f, g contain the values of the step,
         ! function, and derivative at stp.

         stx = 0.0_dp
         fx = finit
         gx = ginit
         sty = 0.0_dp
         fy = finit
         gy = ginit
         stmin = 0.0_dp
         stmax = Stp + xtrapu*Stp
         Task = 'FG'

         call save_locals()
         return

      else

         ! Restore local variables.

         if (Isave(1) == 1) then
            brackt = .true.
         else
            brackt = .false.
         end if
         stage = Isave(2)
         ginit = Dsave(1)
         gtest = Dsave(2)
         gx = Dsave(3)
         gy = Dsave(4)
         finit = Dsave(5)
         fx = Dsave(6)
         fy = Dsave(7)
         stx = Dsave(8)
         sty = Dsave(9)
         stmin = Dsave(10)
         stmax = Dsave(11)
         width = Dsave(12)
         width1 = Dsave(13)

      end if

      ! If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
      ! algorithm enters the second stage.

      ftest = finit + Stp*gtest
      if (stage == 1 .and. f <= ftest .and. g >= 0.0_dp) stage = 2

      ! Test for warnings.

      if (brackt .and. (Stp <= stmin .or. Stp >= stmax)) &
         Task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if (brackt .and. stmax - stmin <= Xtol*stmax) &
         Task = 'WARNING: XTOL TEST SATISFIED'
      if (Stp == Stpmax .and. f <= ftest .and. g <= gtest) &
         Task = 'WARNING: STP = STPMAX'
      if (Stp == Stpmin .and. (f > ftest .or. g >= gtest)) &
         Task = 'WARNING: STP = STPMIN'

      ! Test for convergence.

      if (f <= ftest .and. abs(g) <= Gtol*(-ginit)) Task = 'CONVERGENCE'

      ! Test for termination.

      if (Task(1:4) == 'WARN' .or. Task(1:4) == 'CONV') then
         call save_locals()
         return
      end if

      ! A modified function is used to predict the step during the
      ! first stage if a lower function value has been obtained but
      ! the decrease is not sufficient.

      if (stage == 1 .and. f <= fx .and. f > ftest) then

         ! Define the modified function and derivative values.

         fm = f - Stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

         ! Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx, fxm, gxm, sty, fym, gym, Stp, fm, gm, brackt, stmin, &
                     stmax)

         ! Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else

         ! Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx, fx, gx, sty, fy, gy, Stp, f, g, brackt, stmin, stmax)

      end if

      ! Decide if a bisection step is needed.

      if (brackt) then
         if (abs(sty - stx) >= p66*width1) Stp = stx + p5*(sty - stx)
         width1 = width
         width = abs(sty - stx)
      end if

      ! Set the minimum and maximum steps allowed for stp.

      if (brackt) then
         stmin = min(stx, sty)
         stmax = max(stx, sty)
      else
         stmin = Stp + xtrapl*(Stp - stx)
         stmax = Stp + xtrapu*(Stp - stx)
      end if

      ! Force the step to be within the bounds stpmax and stpmin.

      Stp = max(Stp, Stpmin)
      Stp = min(Stp, Stpmax)

      ! If further progress is not possible, let stp be the best
      ! point obtained during the search.

      if (brackt .and. (Stp <= stmin .or. Stp >= stmax) .or. &
         & (brackt .and. stmax - stmin <= Xtol*stmax)) Stp = stx

      ! Obtain another function and derivative.

      Task = 'FG'

      call save_locals()

   contains

      subroutine save_locals()

         !! Save local variables.

         if (brackt) then
            Isave(1) = 1
         else
            Isave(1) = 0
         end if
         Isave(2) = stage
         Dsave(1) = ginit
         Dsave(2) = gtest
         Dsave(3) = gx
         Dsave(4) = gy
         Dsave(5) = finit
         Dsave(6) = fx
         Dsave(7) = fy
         Dsave(8) = stx
         Dsave(9) = sty
         Dsave(10) = stmin
         Dsave(11) = stmax
         Dsave(12) = width
         Dsave(13) = width1

      end subroutine save_locals

   end subroutine dcsrch
   SUBROUTINE LINESEARCH_STD(p_k, p_kf, GRAD_scaled, m_fine, SCALER, F_SCALE, PAR_SCALE, BALANCE_SCALE, INVP, CRR0, CII0, ICSR, VSR, FSR, &
                             NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, X, DZ0, &
                             NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, MSR, MSR1, COMP, &
                             FREQ, GT0, WAVELET, SourceScaler, WD_amp, WD_acq, IFQ, NK, FKY, WTK, AS, WT, I25D, IS0, IE0, ITER, &
                             FRECHET, GMASK, GRAD_scaled_NORM, GRADr_NORM, CREF, HESSDI_scaled, GRADr, GRAD_smooth, &
                             GRADsc_precond, GRADsc_Precond_NORM, USE_PRECOND, &
                             NBLOCK, IG, USE_GR_SMOOTH, m_coarse, m_coarse_reg, REG_LAMBDA, &
                             ALPHA, FCOST0, FCOSTA, IFLAG, SOLVER_KIND, my_rank, n_process, comm, &
                             PARAM, ITHOM, IVISCO, NCOMP, NCOMPS, FCOST0_ref, GRAD_scaled_ref, p_k_ref, dot_ref, DEBUG_OUTPUT, &
                             sigma_grad_x, sigma_grad_z, XMINC, XMAXC, ZMINC, ZMAXC)

      IMPLICIT NONE

      ! ---------- INTENTs (match MT/BT) ----------
      INTEGER, INTENT(IN) :: NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, ND
      INTEGER, INTENT(IN) :: NTO, NK, IFQ, IE0, IS0, ITER, I25D, NBLOCK, SOLVER_KIND
      INTEGER, INTENT(IN) :: NSR, NSS, comm, my_rank, n_process
      INTEGER, INTENT(IN) :: ITHOM, IVISCO, NCOMP, NCOMPS
      INTEGER, INTENT(IN) :: NS(:), NR(:), NSV(:), NRV(:), INVP(:), COMP(:), MSR(:), MSR1(:, :), ICSR(:)

      REAL(dp), INTENT(IN) :: FREQ, DZ0, X(:), XTO(:), ZTO(:), XSR(:), ZSR(:), YSR(:)
      REAL(dp), INTENT(IN) :: AS(:), WT(:), FKY(:), WTK(:), CREF, PAR_SCALE(:), BALANCE_SCALE(:)
      REAL(dp), INTENT(IN) :: m_coarse(:), m_coarse_reg(:), REG_LAMBDA
      REAL(dp), INTENT(IN) :: XMINC, XMAXC, ZMINC, ZMAXC
      REAL(dp), INTENT(IN) :: p_k(:), p_kf(:), SCALER(:), F_SCALE(:), m_fine(:), FCOST0
      REAL(dp), INTENT(IN) :: CRR0(:, :), CII0(:, :), VSR(:, :, :), FSR(:, :)
      REAL(dp), INTENT(IN) :: GMASK(:), HESSDI_scaled(:)
      REAL(dp), INTENT(INOUT) :: GRADr(:), GRAD_scaled(:), GRAD_smooth(:)
      REAL(dp), INTENT(INOUT) :: GRADsc_precond(:)
      REAL(dp), INTENT(INOUT) :: GRADr_NORM(:), GRAD_scaled_NORM(:), GRADsc_Precond_NORM(:)
      REAL(dp), INTENT(IN)    :: WD_amp, WD_acq(:)
      COMPLEX(dp), INTENT(IN) :: GT0(:), WAVELET(:), SourceScaler(:)
      CHARACTER(LEN=5), INTENT(IN) :: PARAM(22)
      REAL(dp), INTENT(IN) :: FCOST0_ref, dot_ref
      REAL(dp), INTENT(IN) :: p_k_ref(NM), GRAD_scaled_ref(NM)

      TYPE(InversionGridType), INTENT(IN) :: IG
      LOGICAL, INTENT(IN) :: USE_PRECOND, USE_GR_SMOOTH
      REAL(dp), INTENT(IN), OPTIONAL :: sigma_grad_x, sigma_grad_z
      COMPLEX(dp), INTENT(INOUT) :: FRECHET(ND, NM)

      REAL(dp), INTENT(OUT) :: ALPHA, FCOSTA
      INTEGER, INTENT(OUT)  :: IFLAG
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT

      ! ---------- params ----------
      REAL(dp), PARAMETER :: M1 = 1.0e-4_dp
      REAL(dp), PARAMETER :: M2 = 0.9_dp
      REAL(dp), PARAMETER :: MIN_STEP = 1.0e-5_dp
      REAL(dp), PARAMETER :: MAX_STEP = 10.0_dp
      REAL(dp), PARAMETER :: MULT_FACTOR = 2.0_dp
      INTEGER, PARAMETER :: NLS_MAX = 16
      CHARACTER(len=80), PARAMETER :: FMT_ALPHA = &
                                      '(F7.2, "   ", I5, " ", I5, "  ", A, " ", ES13.6, "   ", 2(ES17.10, "   "))'

      ! ---------- locals ----------
      REAL(dp), ALLOCATABLE :: V_STEP(:), PAR_STEP(:), CR1(:, :), CI1(:, :)
      COMPLEX(dp), ALLOCATABLE :: RESID(:)
      COMPLEX(sp), ALLOCATABLE :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      REAL(dp), ALLOCATABLE :: p_km_scaled(:)

      REAL(dp) :: q0, q, alpha_L, alpha_R, new_alpha
      REAL(dp) :: fcost, phi_rhs
      REAL(dp) :: norm_pk, norm_pkm, scale_pk_to_pkm
      REAL(dp), PARAMETER :: eps_r = 1.0e-30_dp
      REAL(dp), PARAMETER :: r_cap = 2.0e-1_dp
      REAL(dp) :: rmax1, alpha_cap, alpha_probe
      INTEGER  :: cpt_ls, ierr, NSSt
      LOGICAL  :: dbg
      LOGICAL  :: cap_ready = .FALSE.
      LOGICAL, SAVE :: ls_ready = .FALSE.
      INTEGER, SAVE :: ls_unit = 65
      LOGICAL  :: is_open
      INTEGER  :: ios
      CHARACTER(len=200) :: iomsg

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      NSSt = NSR

      ! ---------- allocs (same as MT/BT) ----------
      ALLOCATE (CR1(IANISO, NPT), CI1(IANISO, NPT), V_STEP(NMM), PAR_STEP(NMM), RESID(ND))
      ALLOCATE (GFX(NSSt, 3, NK, NPT), GFY(NSSt, 3, NK, NPT), GFZ(NSSt, 3, NK, NPT))
      ALLOCATE (p_km_scaled(NMM))

      ! ---------- log file (only once, rank 0) ----------
      IF (my_rank == 0 .AND. .NOT. ls_ready) THEN
         INQUIRE (UNIT=ls_unit, OPENED=is_open)
         IF (.NOT. is_open) THEN
            OPEN (UNIT=ls_unit, FILE='out_linesearch.txt', STATUS='UNKNOWN', POSITION='APPEND', &
                  ACTION='WRITE', IOSTAT=ios, IOMSG=iomsg)
            IF (ios == 0) THEN
               WRITE (ls_unit, '(A)') '# FREQ   ITER   it   TAG         ALPHA           FCOST1            FCOST0'
               FLUSH (ls_unit)
            ELSE IF (dbg .AND. my_rank == 0) THEN
               WRITE (*, *) 'LINESEARCH_STD: cannot open out_linesearch.txt: ', TRIM(iomsg)
            END IF
         END IF
         ls_ready = .TRUE.
      END IF

      ! ---- init / scaling (match MT/BT) ----
      FCOSTA = FCOST0
      !--------------------------------------------------------
      ! (A) Compute cap once using a probe step (alpha=1.0)
      !     so alpha_cap = r_cap / rmax( alpha=1.0 )
      !--------------------------------------------------------
      IF (.NOT. cap_ready) THEN
         alpha_probe = 1.0_dp

         CALL UPDATEPARAMETER(alpha_probe, NMM, NPT, NPAR, IANISO, my_rank, &
                              p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                              CRR0, CII0, CR1, CI1, F_SCALE)

         ! CALL compute_rmax_rel(CR1, CI1, CRR0, CII0, INVP, NPAR, IANISO, NPT, eps_r, rmax1)
         CALL compute_rmax_rel(CR1, CI1, CRR0, CII0, INVP, NPAR, IANISO, NPT, eps_r, rmax1, my_rank, dbg)
         alpha_cap = 1.0_dp
         IF (rmax1 > r_cap) alpha_cap = r_cap/rmax1

         alpha_cap = MAX(MIN_STEP, MIN(MAX_STEP, alpha_cap))
         cap_ready = .TRUE.

         IF (my_rank == 0 .AND. dbg) THEN
            WRITE (*, '(A,1PE12.4,2X,A,1PE12.4)') 'STD cap probe: rmax1=', rmax1, 'alpha_cap=', alpha_cap
         END IF
      END IF

      !--------------------------------------------------------
      ! (B) Apply the cap to the current trial ALPHA
      !--------------------------------------------------------
      ALPHA = MIN(ALPHA, alpha_cap)
      ALPHA = MAX(MIN_STEP, MIN(MAX_STEP, ALPHA))

      q0 = DOT_PRODUCT(GRAD_scaled(1:NM), p_k(1:NM))
      IF (q0 >= 0.0_dp) THEN
         IF (my_rank == 0) WRITE (*, *) 'STD-LS: non-descent base direction, g·p(0)=', q0
         IFLAG = 11; ALPHA = 0.0_dp
         GOTO 900
      END IF

      ALPHA = 1.0_dp
      ALPHA = MIN(MAX(ALPHA, MIN_STEP), MAX_STEP)
      alpha_L = 0.0_dp
      alpha_R = 0.0_dp
      cpt_ls = 0

      ! ---- main loop  ----
      DO
         CALL UPDATEPARAMETER(ALPHA, NMM, NPT, NPAR, IANISO, my_rank, &
                              p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                              CRR0, CII0, CR1, CI1, F_SCALE)
         CALL MPI_Barrier(comm, ierr)

         CALL subfcostcal(GFX, GFY, GFZ, I25D, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR1, CI1, NORD, &
                          NM, AS, WT, FREQ, NK, FKY, WTK, NSR, NSSt, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                          IE0, IS0, DZ0, GT0, ND, NS, NR, NSV, NRV, YSR, NPT, WAVELET, SourceScaler, IFQ, WD_amp, WD_acq, &
  NBLOCK, IG, fcost, RESID, m_coarse, m_coarse_reg, REG_LAMBDA, my_rank, n_process, comm, ITER, COMP, SOLVER_KIND, DEBUG_OUTPUT=dbg)

         IF (my_rank == 0 .AND. ls_ready) THEN
            WRITE (ls_unit, FMT_ALPHA) FREQ, ITER, cpt_ls + 1, 'ALPHA STD TRY', ALPHA, fcost, FCOST0
            FLUSH (ls_unit)
         END IF

         IF (.NOT. (fcost == fcost)) THEN
            alpha_R = ALPHA
            new_alpha = 0.5_dp*(alpha_L + alpha_R)
            ALPHA = MAX(MIN_STEP, new_alpha)
            cpt_ls = cpt_ls + 1
            IF (cpt_ls >= NLS_MAX) EXIT
            CYCLE
         END IF
         phi_rhs = FCOST0 + M1*ALPHA*q0

         ! Armijo precheck: skip expensive gradient/Frechet if sufficient decrease fails.
         IF (fcost > phi_rhs) THEN
            alpha_R = ALPHA
            new_alpha = 0.5_dp*(alpha_L + alpha_R)
            ALPHA = MAX(MIN_STEP, new_alpha)
         ELSE
            CALL QFRECHET(ITER, PARAM, IFQ, FREQ, NSR, ND, NX, NZ, NNX, NNZ, X, IE0, DZ0, NPT, NBLOCK, &
                          IANISO, ITHOM, IVISCO, NPAR, INVP, IS0, CR1, CI1, &
                          NS, NSV, NR, NRV, NTO, XTO, ZTO, YSR, VSR, &
                          AS, WT, NK, FKY, WTK, GFX, GFY, GFZ, &
                          NM, FRECHET, WAVELET, SourceScaler, NORD, NCOMP, NCOMPS, &
                          IG, my_rank, n_process, comm, DEBUG_OUTPUT=.FALSE.)

            CALL compute_gradient_scaled(FRECHET, RESID, NPAR, F_SCALE, PAR_SCALE, BALANCE_SCALE, GMASK, ND, NM, &
                                         m_coarse, m_coarse_reg, REG_LAMBDA, &
                                         GRADr, GRAD_scaled, GRADr_NORM, GRAD_scaled_NORM, &
                                         NBLOCK, INVP, USE_PRECOND=USE_PRECOND, my_rank=my_rank, &
                                         HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                                         GRADsc_Precond_NORM=GRADsc_Precond_NORM, DEBUG_OUTPUT=.FALSE., &
                                         PARAM=Param, FREQ=freq, ITER=iter, &
                                         NX=nx, NZ=nz, IG=IG, NTO=nto, XTO=xto, ZTO=zto, IE0=IE0, IS0=IS0, &
                                         XMINC=XMINC, XMAXC=XMAXC, ZMINC=ZMINC, ZMAXC=ZMAXC)

            IF (USE_GR_SMOOTH) THEN
               CALL smooth_gradient_blocks( &
                  GRAD_in=GRAD_scaled, GMASK=GMASK, &
                  NBLOCK=NBLOCK, INVP=INVP, NPAR=NPAR, &
                  NX=NX, NZ=NZ, IG=IG, &
                  FREQ=FREQ, CREF=CREF, my_rank=my_rank, &
                  USE_GR_SMOOTH=.TRUE., USE_PRECOND=USE_PRECOND, &
                  HESSDI_scaled=HESSDI_scaled, GRADsc_precond=GRADsc_precond, &
                  GRADsc_Precond_NORM=GRADsc_Precond_NORM, &
                  PARAM=Param, ITER=iter, NTO=nto, XTO=xto, ZTO=zto, IE0=IE0, IS0=IS0, &
                  XMINC=XMINC, XMAXC=XMAXC, ZMINC=ZMINC, ZMAXC=ZMAXC, &
                  DEBUG_OUTPUT=.FALSE., &
                  SIGMA_X=sigma_grad_x, SIGMA_Z=sigma_grad_z, &
                  GRAD_out=GRAD_smooth) !
               GRAD_scaled(1:NM) = GRAD_smooth(1:NM)
               CALL recompute_masked_norms(GRAD_scaled, GRAD_scaled_NORM, INVP, NBLOCK)
            END IF

            q = DOT_PRODUCT(GRAD_scaled(1:NM), p_k(1:NM))

            IF (q >= M2*q0) THEN
               FCOSTA = fcost
               IFLAG = 2
               EXIT
            ELSE
               alpha_L = ALPHA
               IF (alpha_R > 0.0_dp) THEN
                  new_alpha = 0.5_dp*(alpha_L + alpha_R)
               ELSE
                  new_alpha = MIN(MAX_STEP, MULT_FACTOR*ALPHA)
               END IF
               ALPHA = MAX(MIN_STEP, new_alpha)
            END IF
         END IF

         ALPHA = MIN(MAX(ALPHA, MIN_STEP), MAX_STEP)
         cpt_ls = cpt_ls + 1
         IF (cpt_ls >= NLS_MAX) EXIT
         IF (ALPHA <= MIN_STEP*(1.0_dp + 1.0e-12_dp)) EXIT
         IF (ALPHA >= MAX_STEP*(1.0_dp - 1.0e-12_dp) .AND. alpha_R > 0.0_dp) EXIT
      END DO

      IF (IFLAG /= 2) THEN
         IF ((cpt_ls >= NLS_MAX) .AND. (fcost < FCOST0)) THEN
            FCOSTA = fcost
            IFLAG = 2
         ELSE
            IFLAG = 11
            ALPHA = 0.0_dp
            FCOSTA = FCOST0
         END IF
      END IF

900   CONTINUE
      DEALLOCATE (V_STEP, PAR_STEP, CR1, CI1, RESID, GFX, GFY, GFZ, p_km_scaled)
      RETURN
   END SUBROUTINE LINESEARCH_STD

   SUBROUTINE subfcostcal(GFX, GFY, GFZ, I25D, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CRR0, CII0, NORD, &
                          NM, AS, WT, FREQ, NK, FKY, WTK, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                          IE0, IS0, DZ0, GT0, ND, NS, NR, NSV, NRV, YSR, NPT, WAVELET, SourceScaler, IFQ, WD_amp, WD_acq, &
                          NBLOCK, IG, fcost, RESID, m_coarse, m_coarse_reg, REG_LAMBDA, &
                          my_rank, n_process, comm, ITER, COMP, SOLVER_KIND, DEBUG_OUTPUT)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I25D, NTO, NX, NZ, IANISO, NORD, NK, NNX, NNZ, NBLOCK, ITER
      INTEGER, INTENT(IN) :: NSR, NSS, IE0, IS0, ND, NPT, IFQ, SOLVER_KIND, NM
      INTEGER, INTENT(IN) :: my_rank, n_process, comm
      INTEGER, INTENT(IN) :: NS(:), NR(:), NSV(:), NRV(:), COMP(:)
      INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
      REAL(dp), INTENT(IN) :: FREQ, DZ0, WD_amp
      REAL(dp), INTENT(IN) :: WD_acq(:)
      REAL(dp), INTENT(IN) :: m_coarse_reg(:), m_coarse(:), REG_LAMBDA
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:)
      REAL(dp), INTENT(IN) :: AS(:), WT(:), FKY(:), WTK(:)
      REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), YSR(:), FSR(:, :), VSR(:, :, :)
      REAL(dp), INTENT(IN) :: CRR0(:, :), CII0(:, :)

      COMPLEX(dp), INTENT(IN) :: WAVELET(:), GT0(:), SourceScaler(:)
      TYPE(InversionGridType), INTENT(IN) :: IG
      COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)

      COMPLEX(dp), INTENT(INOUT) :: RESID(:)
      REAL(dp), INTENT(OUT) :: fcost
      REAL(dp) ::  FCOST0_amp, RESID_L2, RESID_RMS, FCOST0
      COMPLEX(dp), ALLOCATABLE :: G0(:)
      REAL(dp), ALLOCATABLE :: amp0(:), AMPT0(:)
      LOGICAL, OPTIONAL :: DEBUG_OUTPUT
      fcost = 0.0_dp

      CALL GF(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CRR0, CII0, NORD, &
              AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
              IE0, IS0, DZ0, GFX, GFY, GFZ, NBLOCK, IG, &
              my_rank, n_process, comm, ITER, SOLVER_KIND)

      ALLOCATE (G0(ND), amp0(ND), AMPT0(ND))

      G0 = (0.D0, 0.D0)
      amp0 = 0.D0
      AMPT0 = ABS(GT0)
      CALL Compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, &
                      MSR, MSR1, FSR, GFX, GFY, GFZ, &
                      YSR, FREQ, FKY, WTK, G0, WAVELET, SourceScaler, IFQ, my_rank, amp0, DEBUG_OUTPUT=.TRUE.)

      !   call ComputeFCOST0(ND, G0, GT0, WD_acq, WD_amp, ITER, my_rank, &
      !                       RESID, RESID_L2, RESID_RMS, FCOST0, &
      !                       NS, NR, NSV, NRV, XSR, ZSR, DEBUG_OUTPUT)
      call ComputeFCOST0(ND, G0, GT0, WD_acq, WD_amp, ITER, my_rank, &
                         NM, m_coarse, m_coarse_reg, REG_LAMBDA, &
                         RESID, RESID_L2, RESID_RMS, FCOST0, &
                         NS, NR, NSV, NRV, XSR, ZSR, DEBUG_OUTPUT)
      fcost = FCOST0
      RETURN
   END SUBROUTINE subfcostcal

   SUBROUTINE UPDATEPARAMETER(ALPHA, NMM, NPT, NPAR, IANISO, my_rank, &
                              p_kf, m_fine, PAR_STEP, V_STEP, INVP, PAR_SCALE, &
                              CRR0, CII0, CR1, CI1, F_SCALE)

!-------------------------------------------------------------------------------
! SUBROUTINE: UPDATEPARAMETER
!
! PURPOSE:
!    Computes updated model parameters from a previous model vector and search
!    direction `p_k`, then constructs the updated parameter grids `CR1`, `CI1`
!    for use in alpha trial forward modeling .
!
! ENTRIES:
!    1. ALPHA        - Step length multiplier for search direction.
!    2. NMM          - Total number of model variables (sum_active NPT).
!    3. NPT          - Number of grid points per parameter (e.g., NNX x NNZ).
!    4. NPAR         - Total number of parameters (elastic + attenuation).
!    5. IANISO       - Number of elastic parameters (first IANISO entries).
!    6. p_k(NMM)     - Search direction vector (unscaled).
!    7. m_fine  - Previous model vector (unscaled), length NMM.
!    8. SCALER(:) - Per-parameter scaling factor (usually 1/max(abs(Param))).
!    9. CRR0(IANISO,:) - Current elastic parameter grid (real).
!   10. CII0(IANISO,:) - Current attenuation parameter grid (imag).
!   11. INVP(:)      - Inversion mask per parameter (1 = active, 0 = fixed).
!
! RETURN:
!    1. PAR_STEP(NMM) - Updated unscaled model vector (m_fine + alpha dot p_k).
!    2. V_STEP(NMM)   - Same as PAR_STEP safegard
!    3. CR1(:,:)      - Updated real-valued elastic parameters for modeling.
!    4. CI1(:,:)      - Updated complex-valued attenuation parameters.
!
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NMM, NPT, NPAR, IANISO, my_rank
      REAL(dp), INTENT(IN) :: ALPHA
      REAL(dp), INTENT(IN) :: p_kf(:), m_fine(:), PAR_SCALE(:), F_SCALE(:)
      REAL(dp), INTENT(IN) :: CRR0(:, :), CII0(:, :)
      INTEGER, INTENT(IN) :: INVP(:)
      REAL(dp), INTENT(INOUT) :: PAR_STEP(:), V_STEP(:)
      REAL(dp), INTENT(INOUT) :: CR1(:, :), CI1(:, :)

      INTEGER :: IM, K, ps, pe   ! <-- indices defined here
      CR1 = 0.0_dp
      CI1 = 0.0_dp
      V_STEP = 0.0_dp
      PAR_STEP = 0.0_dp

      !-- Pre-fill trial grids (ensures inactive elastic & Q are correct)
      CR1(:, :) = CRR0(:, :)
      CI1(:, :) = CII0(:, :)

      !-- Trial stacked vector in active-family space scaled
      IF (ALL(PAR_SCALE == 1.0_dp)) THEN
         if (my_rank == 0) WRITE (*, *) 'UPDATEPARAMETER: f-scale only'
         IM = 0
         DO K = 1, NPAR
            IF (INVP(K) /= 1) CYCLE
            IM = IM + 1
            ps = (IM - 1)*NPT + 1
            pe = IM*NPT

            V_STEP(ps:pe) = m_fine(ps:pe) + ALPHA*p_kf(ps:pe)*F_SCALE(K)
            PAR_STEP(ps:pe) = V_STEP(ps:pe)

            IF (K <= IANISO) THEN
               ! Elastic family K -> overwrite real grid row K
               CR1(K, 1:NPT) = PAR_STEP(ps:pe)
            ELSE
               ! Q family: map to corresponding elastic row (no Q for rho)
               CI1(K - (IANISO - 1), 1:NPT) = PAR_STEP(ps:pe)
            END IF
         END DO
      ELSE
         if (my_rank == 0) WRITE (*, *) 'UPDATEPARAMETER: par_scale only or combined scheme, applying par-scale *F-scale'

         IM = 0
         DO K = 1, NPAR
            IF (INVP(K) /= 1) CYCLE
            IM = IM + 1
            ps = (IM - 1)*NPT + 1
            pe = IM*NPT
            PAR_STEP(1:NMM) = m_fine(1:NMM) + ALPHA*p_kf(1:NMM)*F_SCALE(K)
            V_STEP(1:NMM) = PAR_STEP(1:NMM)
            IF (K <= IANISO) THEN
               ! Elastic family K -> overwrite real grid row K
               V_STEP(1:NMM) = PAR_STEP(1:NMM)*PAR_SCALE(K)  ! <-- apply par-scale inverse here
               CR1(K, 1:NPT) = PAR_STEP(ps:pe)*PAR_SCALE(K)
            ELSE
               ! Q family: map to corresponding elastic row (no Q for rho)
               CI1(K - (IANISO - 1), 1:NPT) = PAR_STEP(ps:pe)*PAR_SCALE(K)
            END IF
         END DO
      END IF

      IF (my_rank == 0) THEN
         WRITE (*, '(A,2(1X,ES17.10))') 'UPDATEPARAMETER PAR_STEP min/max, V_STEP min/max =', MINVAL(PAR_STEP), MAXVAL(PAR_STEP),MINVAL(V_STEP), MAXVAL(V_STEP)
      END IF
   END SUBROUTINE UPDATEPARAMETER

   SUBROUTINE ALPHA_INTERP(ALPHA0, ALPHA1, ALPHA2, FCOST0, FCOST1, FCOST2, ALPHAMIN, FCOSTMIN, MY_RANK)
      IMPLICIT real(dp) (A - H, O - Z)

      INTEGER, INTENT(IN) :: MY_RANK
      real(dp), INTENT(IN) :: ALPHA0, ALPHA1, ALPHA2
      real(dp), INTENT(IN) :: FCOST0, FCOST1, FCOST2
      real(dp), INTENT(OUT) :: ALPHAMIN, FCOSTMIN

      real(dp) :: XA(3), YA(3), XX(100), YY(100), OY, DALPHA
      INTEGER :: J, J2, NA, NS, N, M
      real(dp) :: DIF, DIFT, HO, HP, W, DEN, DY
      real(dp), ALLOCATABLE :: C(:), D(:)

      XA(1) = ALPHA0
      XA(2) = ALPHA1
      XA(3) = ALPHA2

      YA(1) = FCOST0
      YA(2) = FCOST1
      YA(3) = FCOST2

      DALPHA = (ALPHA2 - ALPHA0)/50.0D0
      IF (MY_RANK == 0) WRITE (*, *) 'DALPHA:', DALPHA

      NA = 3
      FCOSTMIN = FCOST0
      ALPHAMIN = ALPHA0

      ALLOCATE (C(NA), D(NA))

      DO J = 1, 100
         XX(J) = ALPHA0 + (J - 1)*DALPHA

         ! ==== POLINT ====
         N = NA
         NS = 1
         DIF = DABS(XX(J) - XA(1))

         DO M = 1, N
            DIFT = DABS(XX(J) - XA(M))
            IF (DIFT < DIF) THEN
               NS = M
               DIF = DIFT
            END IF
            C(M) = YA(M)
            D(M) = YA(M)
         END DO

         OY = YA(NS)
         NS = NS - 1

         DO M = 1, N - 1
            DO J2 = 1, N - M
               HO = XA(J2) - XX(J)
               HP = XA(J2 + M) - XX(J)
               W = C(J2 + 1) - D(J2)
               DEN = HO - HP
               IF (DEN == 0.0D0) THEN
                  WRITE (*, *) 'POLINT FAILURE'
                  EXIT
               END IF
               DEN = W/DEN
               D(J2) = HP*DEN
               C(J2) = HO*DEN
            END DO
            IF (2*NS < (N - M)) THEN
               DY = C(NS + 1)
            ELSE
               DY = D(NS)
               NS = NS - 1
            END IF
            OY = OY + DY
         END DO

         YY(J) = OY

         IF (FCOSTMIN > YY(J)) THEN
            FCOSTMIN = YY(J)
            ALPHAMIN = XX(J)
         END IF
      END DO

      DEALLOCATE (C, D)

   END SUBROUTINE ALPHA_INTERP

   !----------------------------------------------------------------------
   !----------------------------------------------------------------------
   SUBROUTINE polint(xa, ya, n, x, y)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: xa(:), ya(:), x
      REAL(dp), INTENT(OUT) :: y

      INTEGER :: j, m, ns
      REAL(dp) :: dif, dift, ho, hp, w, den, dy
      REAL(dp), ALLOCATABLE :: c(:), d(:)

      allocate (d(n), c(n))
      ns = 1
      dif = ABS(x - xa(1))
      do j = 1, n
         dift = ABS(x - xa(j))
         if (dift < dif) then
            ns = j
            dif = dift
         end if
         c(j) = ya(j)
         d(j) = ya(j)
      end do

      y = ya(ns)
      ns = ns - 1
      do m = 1, (n - 1)
         do j = 1, n - m
            ho = xa(j) - x
            hp = xa(j + m) - x
            w = c(j + 1) - d(j)
            den = ho - hp
            if (den == 0.0_dp) then
               write (*, *) 'polint failure'
               exit
            end if
            den = w/den
            d(j) = hp*den
            c(j) = ho*den
         end do
         if ((2*ns) < (n - m)) then
            dy = c(ns + 1)
         else
            dy = d(ns)
            ns = ns - 1
         end if

         y = y + dy
      end do
      deallocate (d, c)
      return
   end

   SUBROUTINE cubic_interp_fcost(alpha_lo, fcost_lo, dfcost_lo, alpha_hi, fcost_hi, dfcost_hi, alpha_j)
      IMPLICIT NONE

      REAL(dp), INTENT(IN) :: alpha_lo, fcost_lo, dfcost_lo
      REAL(dp), INTENT(IN) :: alpha_hi, fcost_hi, dfcost_hi
      REAL(dp), INTENT(OUT) :: alpha_j

      REAL(dp) :: d1, d2, delta_alpha, radicand, numerator, denominator

      delta_alpha = alpha_hi - alpha_lo
      d1 = dfcost_lo + dfcost_hi - 3.0_dp*((fcost_lo - fcost_hi)/delta_alpha)
      radicand = d1*d1 - dfcost_lo*dfcost_hi

      IF (radicand < 0.0_dp) THEN
         alpha_j = 0.5_dp*(alpha_lo + alpha_hi)
         RETURN
      END IF

      d2 = SQRT(radicand)
      IF (alpha_hi > alpha_lo) THEN
         numerator = d2 - dfcost_lo
      ELSE
         numerator = -d2 - dfcost_lo
      END IF

      denominator = dfcost_hi - dfcost_lo + 2.0D0*d2

      IF (denominator == 0.0_dp) THEN
         alpha_j = 0.5_dp*(alpha_lo + alpha_hi)
      ELSE
         alpha_j = alpha_lo + (numerator/denominator)*delta_alpha
      END IF

      IF ((alpha_j < MIN(alpha_lo, alpha_hi)) .OR. (alpha_j > MAX(alpha_lo, alpha_hi))) THEN
         alpha_j = 0.5_dp*(alpha_lo + alpha_hi)
      END IF
   END SUBROUTINE cubic_interp_fcost

   !===========================================================
! Helper: max relative change over active parameters
! - Assumes elastic rows: CR0(1:IANISO,:) are valid
! - Assumes Q rows:      CI0(2:IANISO,:) are valid (row 1 reserved/unused)
! - Uses a per-row robust denominator to avoid blow-ups when CI0 ~ 0
!===========================================================
   ! SUBROUTINE compute_rmax_rel(CR1, CI1, CR0, CI0, INVP, NPAR, IANISO, NPT, eps, rmax)
   !    USE iso_fortran_env, ONLY: dp => real64
   !    IMPLICIT NONE
   !    INTEGER, INTENT(IN) :: NPAR, IANISO, NPT
   !    INTEGER, INTENT(IN) :: INVP(:)
   !    REAL(dp), INTENT(IN) :: CR1(:, :), CI1(:, :)
   !    REAL(dp), INTENT(IN) :: CR0(:, :), CI0(:, :)
   !    REAL(dp), INTENT(IN) :: eps
   !    REAL(dp), INTENT(OUT):: rmax

   !    INTEGER  :: k, itp, row
   !    REAL(dp) :: rloc, den, row_scale

   !    rmax = 0.0_dp

   !    DO k = 1, NPAR
   !       IF (INVP(k) /= 1) CYCLE

   !       IF (k <= IANISO) THEN
   !          ! -------- Elastic: CR rows are 1..IANISO --------
   !          row = k
   !          IF (row < 1 .OR. row > SIZE(CR0, 1)) THEN
   !             WRITE (*, *) 'compute_rmax_rel: CR row OOB: k=', k, ' row=', row, ' size=', SIZE(CR0, 1)
   !             STOP
   !          END IF

   !          ! robust row scale (avoid den ~ 0)
   !          row_scale = MAXVAL(ABS(CR0(row, 1:NPT)))
   !          DO itp = 1, NPT
   !             den = MAX(ABS(CR0(row, itp)), eps*MAX(1.0_dp, row_scale))
   !             rloc = ABS(CR1(row, itp) - CR0(row, itp))/den
   !             rmax = MAX(rmax, rloc)
   !          END DO

   !       ELSE
   !          ! -------- Q: CI rows are 2..IANISO (row 1 reserved) --------
   !          row = k - (IANISO - 1)     ! gives 2..IANISO when k=IANISO+1..2*IANISO-1
   !          IF (row < 2 .OR. row > SIZE(CI0, 1)) THEN
   !             WRITE (*, *) 'compute_rmax_rel: CI row OOB: k=', k, ' row=', row, ' size=', SIZE(CI0, 1)
   !             STOP
   !          END IF

   !          row_scale = MAXVAL(ABS(CI0(row, 1:NPT)))
   !          DO itp = 1, NPT
   !             den = MAX(ABS(CI0(row, itp)), eps*MAX(1.0_dp, row_scale))
   !             rloc = ABS(CI1(row, itp) - CI0(row, itp))/den
   !             rmax = MAX(rmax, rloc)
   !          END DO
   !       END IF
   !    END DO
   ! END SUBROUTINE compute_rmax_rel
   SUBROUTINE compute_rmax_rel(CR1, CI1, CR0, CI0, INVP, NPAR, IANISO, NPT, eps, rmax, my_rank, dbg)
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPAR, IANISO, NPT, my_rank
      INTEGER, INTENT(IN) :: INVP(:)
      REAL(dp), INTENT(IN) :: CR1(:, :), CI1(:, :)
      REAL(dp), INTENT(IN) :: CR0(:, :), CI0(:, :)
      REAL(dp), INTENT(IN) :: eps
      REAL(dp), INTENT(OUT):: rmax
      LOGICAL, INTENT(IN) :: dbg

      INTEGER  :: k, itp, row
      REAL(dp) :: rloc, den, row_scale
      REAL(dp) :: rmax_k

      rmax = 0.0_dp

      DO k = 1, NPAR
         IF (INVP(k) /= 1) CYCLE

         rmax_k = 0.0_dp

         IF (k <= IANISO) THEN
            row = k
            row_scale = MAXVAL(ABS(CR0(row, 1:NPT)))

            DO itp = 1, NPT
               den = MAX(ABS(CR0(row, itp)), eps*MAX(1.0_dp, row_scale))
               rloc = ABS(CR1(row, itp) - CR0(row, itp))/den
               rmax_k = MAX(rmax_k, rloc)
            END DO

            IF (my_rank == 0 .AND. dbg) THEN
               WRITE (*, '(A,I3,A,I3,A,1PE12.4,A,1PE12.4)') &
                  'REL STEP elastic: k=', k, ' row=', row, &
                  ' rmax_k=', rmax_k, ' row_scale=', row_scale
            END IF

         ELSE
            row = k - (IANISO - 1)
            row_scale = MAXVAL(ABS(CI0(row, 1:NPT)))

            DO itp = 1, NPT
               den = MAX(ABS(CI0(row, itp)), eps*MAX(1.0_dp, row_scale))
               rloc = ABS(CI1(row, itp) - CI0(row, itp))/den
               rmax_k = MAX(rmax_k, rloc)
            END DO

            IF (my_rank == 0 .AND. dbg) THEN
               WRITE (*, '(A,I3,A,I3,A,1PE12.4,A,1PE12.4)') &
                  'REL STEP Q      : k=', k, ' row=', row, &
                  ' rmax_k=', rmax_k, ' row_scale=', row_scale
            END IF
         END IF

         rmax = MAX(rmax, rmax_k)
      END DO

   END SUBROUTINE compute_rmax_rel

!       SUBROUTINE LINESEARCH_STD(p_k, p_kf, GRAD_scaled, m_fine, SCALER, INVP, CRR0, CII0, ICSR, VSR, FSR, &
!                              NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, X, DZ0, &
!                              NTO, XTO, ZTO, XSR, ZSR, YSR, NSR, NSS, NS, NR, NSV, NRV, ND, MSR, MSR1, COMP, &
!                              FREQ, GT0, WAVELET, SourceScaler, IFQ, NK, FKY, WTK, AS, WT, I25D, IS0, IE0, ITER, &
!                               FRECHET, WD_amp, WD_acq, &
!                              GMASK, GRAD_scaled_NORM, GRADr_NORM, &
!                              NBLOCK, IG, &
!                              ALPHA, FCOST0, FCOSTMIN, IFLAG, SOLVER_KIND, my_rank, n_process, comm, &
!                              DEBUG_OUTPUT)

!       IMPLICIT real(dp)(A - H, O - Z)

!       ! ---- INTENTs (mirrors MT/BT) ----
!       INTEGER, INTENT(IN) :: NM, NMM, NPT, NPAR, IANISO, NX, NZ, NNX, NNZ, NORD, ND
!       INTEGER, INTENT(IN) :: NTO, NK, IFQ, IE0, IS0, ITER, I25D, NBLOCK, SOLVER_KIND
!       INTEGER, INTENT(IN) :: NSR, NSS, comm, my_rank, n_process
!       INTEGER, INTENT(IN) :: NS(:), NR(:), NSV(:), NRV(:), INVP(:), COMP(:), MSR(:), MSR1(:, :), N0_BLOCK(:), ICSR(:)
!       real(dp), INTENT(IN) :: FREQ, DZ0, X(:), XTO(:), ZTO(:), XSR(:), ZSR(:), YSR(:)
!       real(dp), INTENT(IN) :: AS(:), WT(:), FKY(:), WTK(:)
!       real(dp), INTENT(IN) :: p_k(:), p_kf(:), SCALER(:), m_fine(:), FCOST0
!       real(dp), INTENT(IN) :: CRR0(:, :), CII0(:, :), VSR(:, :, :), FSR(:, :)
!       real(dp), INTENT(IN) :: DX_BLOCK(:), C1_BLOCK(:)
!       real(dp), INTENT(IN) :: Z1_OUT(:, :), Z2_OUT(:, :), T1_OUT(:, :), T2_OUT(:, :)
!       real(dp), INTENT(IN) :: GRADr_NORM(:), GRAD_scaled_NORM(:), GMASK(:)
!       real(dp), INTENT(IN) ::  WD_amp, WD_acq(:)
!       COMPLEX(dp), INTENT(IN)       :: GT0(:), WAVELET(:), FRECHET(ND, NM), SourceScaler(:)
!       real(dp), INTENT(IN) :: GRAD_scaled(:)

!       real(dp), INTENT(OUT):: ALPHA, FCOSTMIN
!       INTEGER, INTENT(OUT):: IFLAG
!       LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT

!       ! ---- Parameters (Wolfe + bounds + limits) ----
!       real(dp), PARAMETER :: M1 = 1.0D-4        ! Armijo
!       real(dp), PARAMETER :: M2 = 0.9D0         ! Curvature (sign-aware test, as in your snippet)
!       real(dp), PARAMETER :: MIN_STEP = 1.0D-08
!       real(dp), PARAMETER :: MAX_STEP = 10.0D0
!       real(dp), PARAMETER :: MULT_FACTOR = 2.0D0
!       INTEGER, PARAMETER :: NLS_MAX = 20

!       ! ---- Locals ----
!       real(dp), ALLOCATABLE :: V_STEP(:), PAR_STEP(:), CR1(:, :), CI1(:, :)
!       COMPLEX(dp), ALLOCATABLE   :: RESID(:), GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
!       real(dp), ALLOCATABLE :: GRADr_loc(:), GRAD_scaled_loc(:)
!       real(dp), ALLOCATABLE :: GRADr_NORM_loc(:), GRAD_scaled_NORM_loc(:)

!       real(dp) :: q0, q, alpha_L, alpha_R, new_alpha
!       real(dp) :: fcost, phi_rhs
!       INTEGER :: cpt_ls, ierr
!       LOGICAL :: dbg
!       LOGICAL, SAVE :: ls_ready = .FALSE.
!       INTEGER, SAVE :: ls_unit = 65
!       LOGICAL :: is_open
!       INTEGER :: ios
!       CHARACTER(len=200) :: iomsg
!       CHARACTER(len=80), PARAMETER :: FMT_ALPHA = '(F7.2, "   ", I5, " ", I5, "  ", A, " ", ES13.6, "   ", 2(ES17.10, "   "))'
!       dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

!       ! ---- allocs ----
!       ALLOCATE (CR1(IANISO, NPT), CI1(IANISO, NPT), V_STEP(NMM), PAR_STEP(NMM), RESID(ND))
!       ALLOCATE (GFX(NSS, 3, NK, NPT), GFY(NSS, 3, NK, NPT), GFZ(NSS, 3, NK, NPT))
!       ALLOCATE (GRADr_loc(NM), GRAD_scaled_loc(NM))
!       ALLOCATE (GRADr_NORM_loc(NPAR), GRAD_scaled_NORM_loc(NPAR))

!       IF (my_rank == 0 .AND. .NOT. ls_ready) THEN
!          INQUIRE (UNIT=ls_unit, OPENED=is_open)
!          IF (.NOT. is_open) THEN
!             OPEN (UNIT=ls_unit, FILE='out_linesearch.txt', STATUS='UNKNOWN', POSITION='APPEND', &
!                   ACTION='WRITE', IOSTAT=ios, IOMSG=iomsg)
!             IF (ios /= 0) THEN
!                WRITE (*, *) 'LINESEARCH_BACKTRACK: cannot open out_linesearch.txt: ', TRIM(iomsg)
!             ELSE
!                ! Optional: header (only once)
!                WRITE (ls_unit, '(A)') '# FREQ   ITER   it   TAG         ALPHA           FCOST1            FCOST0'
!                FLUSH (ls_unit)
!             END IF
!          END IF
!          ls_ready = .TRUE.
!       END IF
!       ! ---- init (coarse) ----
!       FCOSTMIN = FCOST0
!       q0 = DOT_PRODUCT(GRAD_scaled, p_k)     ! g(0)·p  (must be < 0 for descent)
!       IF (q0 >= 0.0D0) THEN
!          IF (my_rank == 0) WRITE (*, *) 'STD-LS: non-descent base direction, g·p(0)=', q0
!          IFLAG = 11; ALPHA = 0.0D0
!          GOTO 900
!       END IF

!       ! start from a safe α
!       ALPHA = 1.0D0
!       ALPHA = MIN(MAX(ALPHA, MIN_STEP), MAX_STEP)

!       alpha_L = 0.0D0
!       alpha_R = 0.0D0
!       cpt_ls = 0

!       IF (my_rank == 0 .AND. dbg) THEN
!          WRITE (*, *) 'STD-LS init: phi(0)=', FCOST0, '  g·p(0)=', q0, '  alpha=', ALPHA
!       END IF

!       ! ---- main loop ----
!       DO
!          ! Build trial params at current alpha (coarse step -> physical CR1/CI1)
!          CALL UPDATEPARAMETER(ALPHA, NMM, NPT, NPAR, IANISO, my_rank, &
!                               p_k, m_fine, PAR_STEP, V_STEP, INVP, SCALER, &
!                               CRR0, CII0, CR1, CI1)
!          CALL MPI_Barrier(comm, ierr)

!          ! Evaluate cost and residual
!          CALL subfcostcal(GFX, GFY, GFZ, I25D, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR1, CI1, NORD, &
!                           AS, WT, FREQ, NK, FKY, WTK, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
!                           IE0, IS0, DZ0, GT0, ND, NS, NR, NSV, NRV, YSR, NPT, WAVELET, SourceScaler, IFQ, WD_amp, WD_acq, &
!                           NBLOCK, DX_BLOCK, C1_BLOCK, N0_BLOCK, Z1_OUT, Z2_OUT, IG,T1_OUT, T2_OUT, fcost, RESID, &
!                           my_rank, n_process, comm, ITER, COMP, SOLVER_KIND, DEBUG_OUTPUT=dbg)

!          IF (.NOT. (fcost == fcost)) THEN
!             ! NaN cost: shrink
!             alpha_R = ALPHA
!             new_alpha = 0.5D0*(alpha_L + alpha_R)
!             ALPHA = MAX(MIN_STEP, new_alpha)
!             cpt_ls = cpt_ls + 1
!             IF (cpt_ls >= NLS_MAX) EXIT
!             CYCLE
!          END IF

!          ! Compute coarse gradient at alpha for curvature test (but Armijo uses q0)
!          CALL compute_gradient_scaled(FRECHET, RESID, NPAR, SCALER, GMASK, ND, NM, &
!                                       GRADr_loc, GRAD_scaled_loc, GRADr_NORM_loc, GRAD_scaled_NORM_loc, &
!                                       NBLOCK, INVP, NX, NZ,  NTO, XTO, ZTO, ITER, FREQ, my_rank, WD_acq)

!          q = DOT_PRODUCT(GRAD_scaled_loc, p_k)     ! g(α)·p (coarse)

!          IF (my_rank == 0 .AND. dbg) THEN
!             WRITE (*, *) 'STD-LS try: α=', ALPHA, '  φ(α)=', fcost, '  g(α)·p=', q
!          END IF
!          IF (my_rank == 0 .AND. ls_ready) THEN
!             WRITE (*, *) 'STD-LS try alpha=', ALPHA, '  phi(alpha)=', FCOST, '  Armijo rhs=', q
!             WRITE (65, FMT_ALPHA) FREQ, ITER, it, 'ALPHA BT TRY', ALPHA, FCOST, FCOST0
!             FLUSH (65)
!          END IF
!          ! ---- Wolfe tests (as in your snippet) ----
!          phi_rhs = FCOST0 + M1*ALPHA*q0

!          IF ((fcost <= phi_rhs) .AND. (q >= M2*q0)) THEN
!             ! Success
!             FCOSTMIN = fcost
!             IFLAG = 2
!             IF (my_rank == 0 .AND. dbg) WRITE (*, *) 'STD-LS: Wolfe satisfied, accept α=', ALPHA
!             EXIT
!          ELSEIF (fcost > phi_rhs) THEN
!             ! Failure 1: Armijo violated -> shrink (bracket right)
!             alpha_R = ALPHA
!             new_alpha = 0.5D0*(alpha_L + alpha_R)
!             ALPHA = MAX(MIN_STEP, new_alpha)
!          ELSE
!             ! Failure 2: Armijo OK but curvature not -> expand if not bracketed
!             alpha_L = ALPHA
!             IF (alpha_R > 0.0D0) THEN
!                new_alpha = 0.5D0*(alpha_L + alpha_R)
!             ELSE
!                new_alpha = MIN(MAX_STEP, MULT_FACTOR*ALPHA)
!             END IF
!             ALPHA = MAX(MIN_STEP, new_alpha)
!          END IF

!          ! Clamp and termination checks
!          ALPHA = MIN(MAX(ALPHA, MIN_STEP), MAX_STEP)
!          cpt_ls = cpt_ls + 1
!          IF (cpt_ls >= NLS_MAX) EXIT
!          IF (ALPHA <= MIN_STEP*(1.0D0 + 1.0D-12)) EXIT
!          IF (ALPHA >= MAX_STEP*(1.0D0 - 1.0D-12) .AND. alpha_R > 0.0D0) EXIT

!       END DO

!       ! ---- Fallback / accept-if-decrease ----
!       IF (IFLAG /= 2) THEN
!          IF ((cpt_ls >= NLS_MAX) .AND. (fcost < FCOST0)) THEN
!             ! accept decrease even if Wolfe not met (matches your snippet policy)
!             FCOSTMIN = fcost
!             IFLAG = 2
!             IF (my_rank == 0 .AND. dbg) WRITE (*, *) 'STD-LS: max iters, but φ decreased; accept α=', ALPHA
!          ELSE
!             ! fail
!             IFLAG = 11
!             ALPHA = 0.0D0
!             FCOSTMIN = FCOST0
!             IF (my_rank == 0) WRITE (*, *) 'STD-LS: failed to satisfy conditions'
!          END IF
!       END IF

! 900   CONTINUE
!    DEALLOCATE (V_STEP, PAR_STEP, CR1, CI1, RESID, GFX, GFY, GFZ, GRADr_loc, GRAD_scaled_loc, GRADr_NORM_loc, GRAD_scaled_NORM_loc)
!       RETURN
!    END SUBROUTINE LINESEARCH_STD

end module

