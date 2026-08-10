module misfit_mod
   use shared_mod
   use grid_mod
   use output_mod
   USE gridtype_mod
   use constant_mod
   use iso_fortran_env, ONLY: dp => real64
   implicit none

contains

SUBROUTINE ComputeFCOST0(ND, G0, GT0, WD_acq, WD_amp, ITER, my_rank, &
                         NM, MODEL_CURR, MODEL_REF, REG_LAMBDA, &
                         RESID, RESID_L2, RESID_RMS, FCOST0,  &
                         NS, NR, NSV, NRV, XSR, ZSR, DEBUG_OUTPUT)

   IMPLICIT NONE

   ! -------- inputs --------
   INTEGER, INTENT(IN) :: ND, ITER, my_rank
   INTEGER, INTENT(IN) :: NM
   COMPLEX(dp), INTENT(IN) :: G0(:), GT0(:)          ! complex model/target
   REAL(dp),    INTENT(IN) :: WD_acq(:)              ! acquisition weights (>=0)
   REAL(dp),    INTENT(IN) :: WD_amp                 ! scalar amplitude weight

   REAL(dp),    INTENT(IN) :: MODEL_CURR(:)          ! current model vector (size NM)
   REAL(dp),    INTENT(IN) :: MODEL_REF(:)           ! reference model vector (size NM)
   REAL(dp),    INTENT(IN) :: REG_LAMBDA             ! 0th-order Tikhonov weight

   ! -------- outputs --------
   COMPLEX(dp), INTENT(OUT) :: RESID(:)              ! for gradient:
                                                     ! RESID = WD_acq * (WD_amp^2 * (G0 - GT0))
   REAL(dp),    INTENT(OUT) :: RESID_L2              ! weighted L2 of amp-scaled residual
   REAL(dp),    INTENT(OUT) :: RESID_RMS             ! weighted RMS of amp-scaled residual
   REAL(dp),    INTENT(OUT) :: FCOST0                ! total cost = data + regularization
   ! REAL(dp),    INTENT(OUT) :: FCOST_DATA            ! data misfit part only
   ! REAL(dp),    INTENT(OUT) :: FCOST_REG             ! regularization part only

   ! -------- optionals (for debug dump) --------
   LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT
   INTEGER, OPTIONAL, INTENT(IN) :: NS(:), NR(:), NSV(:), NRV(:)
   REAL(dp), OPTIONAL, INTENT(IN) :: XSR(:), ZSR(:)

   ! -------- locals --------
   INTEGER  :: I
   REAL(dp):: FCOST_REG                ! total cost = data + regularization
   REAL(dp) :: FCOST_DATA            ! data misfit part only
   REAL(dp) :: sum_wr2, sum_w, sum_reg, dm
   LOGICAL  :: dbg, is_open
   COMPLEX(dp) :: r0(ND), r_sc

   dbg = .FALSE.
   IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

   sum_wr2 = 0.0_dp
   sum_w   = 0.0_dp
   sum_reg = 0.0_dp

   !-------------------------------------------------------
   ! DATA TERM
   !
   ! 1) r0(i)      = G0(i) - GT0(i)           (raw residual)
   ! 2) r_sc(i)    = WD_amp * r0(i)           (amplitude-scaled)
   ! 3) Phi_data   = 0.5 * SUM WD_acq(i)*|r_sc(i)|^2
   ! 4) 4) RESID(i) = WD_acq(i) * WD_amp**2 * r0(i)
   !
   ! => effective W_d(i) = WD_amp^2 * WD_acq(i)
   !    Phi_data = 0.5 * r^T W_d r
   !    g_data   = J^T W_d r = J^T [WD_acq * r_sc]
   !-------------------------------------------------------
   DO I = 1, ND

      ! raw residual
      r0(I) = G0(I) - GT0(I)

      ! amplitude-scaled residual
      r_sc = WD_amp * r0(I)

      ! cost accumulation using acquisition weights (NOT squared)
      sum_wr2 = sum_wr2 + WD_acq(I) * (REAL(r_sc, dp)**2 + AIMAG(r_sc)**2)
      sum_w   = sum_w   + WD_acq(I)

      ! scaled residual for gradient
      RESID(I) = CMPLX(WD_amp * WD_amp * WD_acq(I), 0.0_dp, dp) * r0(I)

   END DO

   FCOST_DATA = 0.5_dp * sum_wr2

   !-------------------------------------------------------
   ! 0th-ORDER TIKHONOV REGULARIZATION
   !
   ! Phi_reg = 0.5 * REG_LAMBDA * SUM_j ( MODEL_CURR(j) - MODEL_REF(j) )^2
   !-------------------------------------------------------
   IF (REG_LAMBDA > 0.0_dp) THEN
      DO I = 1, NM
         dm = MODEL_CURR(I) - MODEL_REF(I) !sqquer tg\he difference
         sum_reg = sum_reg + dm * dm
      END DO
      FCOST_REG = 0.5_dp * REG_LAMBDA * sum_reg / NM  ! note: dividing by NM here to make the regularization term more comparable in magnitude to the
   ELSE
      FCOST_REG = 0.0_dp
   END IF

   ! total cost
   FCOST0 = FCOST_DATA + FCOST_REG

   ! Weighted L2 and RMS of the amp-scaled residual (data term only)
   RESID_L2 = SQRT(MAX(sum_wr2, 0.0_dp))
   IF (sum_w > 0.0_dp) THEN
      RESID_RMS = SQRT(MAX(sum_wr2 / sum_w, 0.0_dp))
   ELSE
      RESID_RMS = 0.0_dp
   END IF

   IF (my_rank == 0) THEN
      WRITE (*, '(A,1PE12.5,3(A,1PE12.5))') &
         ' WD_amp=', WD_amp, '  RESID_RMS=', RESID_RMS, &
         '  FCOST_DATA=', FCOST_DATA, '  FCOST0=', FCOST0

      IF (REG_LAMBDA > 0.0_dp) THEN
         WRITE (*, '(A,1PE12.5,A,1PE12.5)') &
            ' REG_LAMBDA=', REG_LAMBDA, '  FCOST_REG=', FCOST_REG
      END IF
   END IF

   ! ----- optional early debug trace dump using r0 & WD_acq -----
   IF (dbg .AND. my_rank == 0 .AND. ITER < 3) THEN
      INQUIRE (UNIT=64, OPENED=is_open)
      IF (is_open) THEN
         DO I = 1, ND
            IF (PRESENT(NS) .AND. PRESENT(NRV) .AND. PRESENT(NR) .AND. PRESENT(NSV) .AND. &
                PRESENT(XSR) .AND. PRESENT(ZSR)) THEN
               WRITE (64, '(I5," ", I5,"  ",I5,"   ",I5,"   ",I5,"   ",I5,"  ", &
      &                   F10.3,"   ",F10.3,"   ",F10.3,"   ",F10.3,"   ",ES22.14," ",ES22.14,"  ",1PE11.3)') &
                      ITER, I, NS(I), NSV(I), NR(I), NRV(I), &
                      XSR(NS(I)), ZSR(NS(I)), XSR(NR(I)), ZSR(NR(I)), &
                      REAL(r0(I), dp), AIMAG(r0(I)), WD_acq(I)
            ELSE
               WRITE (64, '(I5," ",I5,"  ",ES22.14," ",ES22.14,"  ",1PE11.3,ES22.14," ",ES22.14)') &
                  ITER, I, REAL(r0(I), dp), AIMAG(r0(I)), WD_acq(I),REAL(RESID(I), dp), AIMAG(RESID(I))
            END IF
         END DO
      END IF
   END IF

END SUBROUTINE ComputeFCOST0
 


   SUBROUTINE build_data_weights(NFQ, ND, NSS, NS, NR, NSV, NRV, XSR, ZSR, &
                                 n_data_weight, W_D, ierr, my_rank)

      IMPLICIT real(dp) (A - H, O - Z)

      !====================================================================
      ! build_data_weights
      !--------------------------------------------------------------------
      ! Purpose:
      !   Build per-trace data weights W_D(J) using vertical-distance bins
      !   read from a simple text file. Falls back to uniform weights with
      !   clear error messages and codes when anything goes wrong.
      !
      ! Error codes (ierr):
      !    0  OK
      !   10  Could not OPEN weight file -> using uniform weights
      !   11  Failed to READ bin sizes (dx_bin,dz_bin) -> uniform
      !   12  Invalid nbin (<=0) from n_data_weight -> uniform
      !   13  Failed to READ weight_table -> uniform
      !   14  Failed to READ component weights (c1..c3) -> set to 1, continue
      !   15  Non-positive dz_bin -> uniform
      ! 1001  ALLOCATE(W_D) failed
      ! 1002  ALLOCATE(weight_table) failed -> uniform
      ! 2001  Index out of bounds for NS(J) or NR(J) -> set W_D(J)=1 for that J
      !
      ! Notes:
      !   - Does NOT change the interface.
      !   - Prints detailed diagnostics when my_rank==0.
      !====================================================================

      ! ---- inputs ----
      INTEGER, INTENT(IN)  :: NFQ, ND, NSS, n_data_weight, my_rank
      INTEGER, INTENT(IN)  :: NS(:), NR(:), NSV(:), NRV(:)
      real(dp), INTENT(IN)  :: XSR(:), ZSR(:)

      ! ---- outputs ----
      real(dp), ALLOCATABLE, INTENT(OUT) :: W_D(:)
      INTEGER, INTENT(OUT) :: ierr

      ! ---- locals ----
      INTEGER :: J, ios, idist, nbin, rc, stat_alloc
      real(dp) :: xs, zs, xr, zr, dz, dist
      real(dp) :: dx_bin, dz_bin, bin_w, w_geom, w_comp
      real(dp), ALLOCATABLE :: weight_table(:)
      real(dp) :: c_weight(3)
      LOGICAL :: verbose, dump_file
      INTEGER :: unit_file
      CHARACTER(len=256) :: dbg_name

      CHARACTER(len=*), PARAMETER :: data_weight_file = &
                                     '/home/kunet.ae/100061882/Models/TTI/data_weight_file.txt'

      ierr = 0
      nbin = n_data_weight
      verbose = (my_rank == 0)
      dump_file = .FALSE.      ! set to .TRUE. if you also want a file dump

      ! --- allocate W_D ---
      ALLOCATE (W_D(ND), STAT=stat_alloc)
      IF (stat_alloc /= 0) THEN
         ierr = 1001
         IF (verbose) WRITE (*, '(A,I0)') '[build_data_weights] ERROR: ALLOCATE(W_D) failed, STAT=', stat_alloc
         RETURN
      END IF

      ! --- open weight file ---
      OPEN (45, FILE=TRIM(data_weight_file), STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
         ierr = 10
     IF (verbose) WRITE(*,'(A,1X,A,1X,A,I0)') '[build_data_weights] WARNING:cannot open weight file:', TRIM(data_weight_file), 'IOSTAT=', ios
         W_D(1:ND) = 1D0
         RETURN
      END IF

      ! --- read line 1: dx_bin, dz_bin ---
      READ (45, *, IOSTAT=ios) dx_bin, dz_bin
      IF (ios /= 0) THEN
         ierr = 11
         IF (verbose) WRITE (*, '(A,I0)') '[build_data_weights] WARNING: failed to read (dx_bin,dz_bin), IOSTAT=', ios
         CLOSE (45)
         W_D(1:ND) = 1D0
         RETURN
      END IF

      ! --- validate nbin ---
      IF (nbin <= 0) THEN
         ierr = 12
         IF (verbose) WRITE (*, '(A,I0)') '[build_data_weights] WARNING: invalid nbin from n_data_weight =', nbin
         CLOSE (45)
         W_D(1:ND) = 1D0
         RETURN
      END IF

      ! --- read line 2: weight_table(1..nbin) ---
      ALLOCATE (weight_table(nbin), STAT=stat_alloc)
      IF (stat_alloc /= 0) THEN
         ierr = 1002
         IF (verbose) WRITE (*, '(A,I0)') '[build_data_weights] ERROR: ALLOCATE(weight_table) failed, STAT=', stat_alloc
         CLOSE (45)
         W_D(1:ND) = 1D0
         RETURN
      END IF

      READ (45, *, IOSTAT=ios) (weight_table(J), J=1, nbin)
      IF (ios /= 0) THEN
         ierr = 13
         IF (verbose) WRITE (*, '(A,I0)') '[build_data_weights] WARNING: failed to read weight table, IOSTAT=', ios
         DEALLOCATE (weight_table)
         CLOSE (45)
         W_D(1:ND) = 1D0
         RETURN
      END IF

      ! --- read line 3: component weights (c1..c3); fallback to 1 if missing ---
      c_weight = 1D0
      READ (45, *, IOSTAT=ios) c_weight(1), c_weight(2), c_weight(3)
      IF (ios /= 0) THEN
         ierr = 14   ! non-fatal; continue
         c_weight = 1D0
         IF (verbose) WRITE (*, '(A,I0)') '[build_data_weights] NOTICE: component weights not found; using 1, IOSTAT=', ios
      END IF

      CLOSE (45)

      ! --- choose binning axis (fixed: vertical distance) ---
      bin_w = dz_bin
      IF (bin_w <= 0D0) THEN
         ierr = 15
 IF (verbose) WRITE(*,'(A,1X,ES12.4)') '[build_data_weights] WARNING: non-positive dz_bin -> using uniform weights; dz_bin=', dz_bin
         DEALLOCATE (weight_table)
         W_D(1:ND) = 1D0
         RETURN
      END IF

      ! --- optional per-rank dump file ---
      IF (dump_file) THEN
         WRITE (dbg_name, '(A,I0,A)') 'W_D_debug_rank', my_rank, '.txt'
         OPEN (UNIT=99, FILE=TRIM(dbg_name), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
         IF (ios == 0) THEN
            unit_file = 99
            IF (verbose) WRITE (*, '(A,1X,A)') '[build_data_weights] dumping W_D to', TRIM(dbg_name)
            WRITE (unit_file, '(A)') '# J   idist   w_geom    rc   w_comp     W_D(J)'
         ELSE
            unit_file = -1
            IF (verbose) WRITE (*, '(A)') '[build_data_weights] WARNING: cannot open debug dump file; continuing without file dump.'
         END IF
      ELSE
         unit_file = -1
      END IF

      ! --- build per-trace weights ---
      DO J = 1, ND

         ! bounds check for indices used to access XSR/ZSR
         IF (J < 1 .OR. J > SIZE(NS)) THEN
            ierr = MAX(ierr, 2001)
        IF (verbose) WRITE (*, '(A,I0,A,I0)') '[build_data_weights] ERROR: J out of bounds for NS(:): J=', J, ' SIZE(NS)=', SIZE(NS)
            W_D(J) = 1D0
            CYCLE
         END IF
         IF (J < 1 .OR. J > SIZE(NR)) THEN
            ierr = MAX(ierr, 2001)
        IF (verbose) WRITE (*, '(A,I0,A,I0)') '[build_data_weights] ERROR: J out of bounds for NR(:): J=', J, ' SIZE(NR)=', SIZE(NR)
            W_D(J) = 1D0
            CYCLE
         END IF

         IF (NS(J) < 1 .OR. NS(J) > SIZE(XSR) .OR. NS(J) > SIZE(ZSR)) THEN
            ierr = MAX(ierr, 2001)
    IF (verbose) WRITE (*, '(A,I0,A,I0)') '[build_data_weights] ERROR: NS(J) out of bounds: NS(J)=', NS(J), ' SIZE(XSR)=', SIZE(XSR)
            W_D(J) = 1D0
            CYCLE
         END IF
         IF (NR(J) < 1 .OR. NR(J) > SIZE(XSR) .OR. NR(J) > SIZE(ZSR)) THEN
            ierr = MAX(ierr, 2001)
    IF (verbose) WRITE (*, '(A,I0,A,I0)') '[build_data_weights] ERROR: NR(J) out of bounds: NR(J)=', NR(J), ' SIZE(XSR)=', SIZE(XSR)
            W_D(J) = 1D0
            CYCLE
         END IF

         xs = XSR(NS(J)); zs = ZSR(NS(J))
         xr = XSR(NR(J)); zr = ZSR(NR(J))

         dz = DABS(zs - zr)
         dist = dz

         idist = INT(dist/bin_w) + 1
         IF (idist < 1) idist = 1
         IF (idist > nbin) idist = nbin

         w_geom = weight_table(idist)

         rc = 1
         IF (J >= 1 .AND. J <= SIZE(NRV)) THEN
            rc = MAX(1, MIN(3, NRV(J)))
         ELSE
            ierr = MAX(ierr, 2001)
   IF (verbose) WRITE (*, '(A,I0,A,I0)') '[build_data_weights] WARNING: J out of bounds for NRV(:): J=', J, ' SIZE(NRV)=', SIZE(NRV)
         END IF
         w_comp = c_weight(rc)

         W_D(J) = w_geom*w_comp

         ! --- per-trace prints ---
         IF (verbose) THEN
            WRITE (*, '(A,I0,2X,A,ES12.4,2X,A,ES12.4,2X,A,I0,2X,A,ES12.4,2X,A,I0,2X,A,ES12.4)') &
               'J=', J, 'dist=', dist, 'bin_w=', bin_w, 'idist=', idist, 'w_geom=', w_geom, &
               'rc=', rc, 'w_comp=', w_comp
            WRITE (*, '(A,I0,A,ES12.4)') 'W_D(', J, ')=', W_D(J)
         END IF
         IF (unit_file > 0) THEN
            WRITE (unit_file, '(I6,2X,I6,2X,ES12.4,2X,I2,2X,ES12.4,2X,ES12.4)') J, idist, w_geom, rc, w_comp, W_D(J)
         END IF
      END DO

      IF (ALLOCATED(weight_table)) DEALLOCATE (weight_table)
      IF (unit_file > 0) CLOSE (unit_file)

   END SUBROUTINE build_data_weights

   SUBROUTINE check_convergence(iter, maxiter, NM, &
                                m_coarse_prev, m_coarse, &
                                FCOST_prev, FCOST_curr, &
                                norm_init, grad_norm, grad0_norm, &
                                cost_conv_tol, model_conv_tol, &
                                MIN_ITER_BEFORE_CHECK, &
                                my_rank, converged)

      IMPLICIT NONE

      ! Kinds / constants
      INTEGER, PARAMETER :: dp = KIND(1.0D0)
      REAL(dp), PARAMETER :: tiny_dp = 1.0e-30_dp
      REAL(dp), PARAMETER :: grad_conv_tol = 1.0e-3_dp   ! <<< internal gradient tol

      ! Inputs
      INTEGER, INTENT(IN) :: iter, maxiter, NM, my_rank
      REAL(dp), INTENT(IN) :: m_coarse_prev(NM), m_coarse(NM)
      REAL(dp), INTENT(IN) :: FCOST_prev, FCOST_curr
      REAL(dp), INTENT(IN) :: grad_norm
      REAL(dp), INTENT(IN) :: cost_conv_tol, model_conv_tol
      INTEGER, INTENT(IN)  :: MIN_ITER_BEFORE_CHECK

      ! Input / output
      REAL(dp), INTENT(INOUT) :: norm_init
      REAL(dp), INTENT(INOUT) :: grad0_norm

      ! Output
      LOGICAL, INTENT(OUT) :: converged

      ! Locals
      REAL(dp) :: rel_cost_drop
      REAL(dp) :: norm_diff
      REAL(dp) :: grad_rel
      REAL(dp), ALLOCATABLE :: diff_model(:)
      LOGICAL :: cost_converged, model_converged, grad_converged, cost_increased

      !--------------------------------------------------------------------
      converged = .FALSE.
      cost_converged = .FALSE.
      model_converged = .FALSE.
      grad_converged = .FALSE.
      cost_increased = .FALSE.
      grad_rel = 1.0_dp
      norm_diff = 0.0_dp

      ! Initialize reference gradient norm (for gradient criterion)
      IF (grad0_norm <= tiny_dp .AND. grad_norm > tiny_dp) THEN
         grad0_norm = grad_norm
      END IF
      IF (grad0_norm <= tiny_dp) grad0_norm = 1.0_dp

      ! If too early in iterations, or trivial maxiter, do not check yet
      IF (iter < MIN_ITER_BEFORE_CHECK) RETURN

      !------------------------- Cost criterion ----------------------------
      rel_cost_drop = (FCOST_prev - FCOST_curr)/ &
                      MAX(ABS(FCOST_prev), tiny_dp)
                            rel_cost_drop = (FCOST_prev - FCOST_curr)/ ABS(FCOST_prev)
      cost_increased = (rel_cost_drop < 0.0_dp)

      IF (.NOT. cost_increased .AND. rel_cost_drop <= cost_conv_tol) THEN
         cost_converged = .TRUE.
      END IF

      !----------------------- Model update criterion ----------------------
      ! ALLOCATE (diff_model(NM))
      ! diff_model(1:NM) = m_coarse(1:NM) - m_coarse_prev(1:NM)

      ! norm_diff = DNRM2(NM, diff_model, 1)/MAX(norm_init, tiny_dp)
      ! IF (norm_diff <= model_conv_tol) THEN
      !    model_converged = .TRUE.
      ! END IF

      ! DEALLOCATE (diff_model)

      !---------------------- Gradient criterion ---------------------------
      IF (grad_norm > tiny_dp) THEN
         grad_rel = grad_norm/MAX(grad0_norm, tiny_dp)
         IF (grad_rel <= grad_conv_tol) THEN
            grad_converged = .TRUE.
         END IF
      END IF

      !---------------------- Combine and report ---------------------------
      IF (my_rank == 0) THEN
         WRITE (*, '(A,I4)') 'Convergence check at iteration ', iter
         WRITE (*, '(A,1PE11.4)') '  Relative cost decrease  = ', rel_cost_drop
         WRITE (*, '(A,1PE11.4)') '  Relative model update   = ', norm_diff
         WRITE (*, '(A,1PE11.4)') '  Relative grad. norm     = ', grad_rel
         IF (cost_increased) WRITE (*, '(A)') '  Cost increased: convergence is disabled for this iterate.'
      END IF

      ! Policy: declare convergence if ANY criterion is satisfied
      IF ((.NOT. cost_increased) .AND. (cost_converged .OR. model_converged .OR. grad_converged)) THEN
         converged = .TRUE.
         IF (my_rank == 0) THEN
            WRITE (*, *) '*** Convergence conditions satisfied:'
            IF (cost_converged) WRITE (*, *) '    - Cost change below tolerance.'
            IF (model_converged) WRITE (*, *) '    - Model update below tolerance.'
            IF (grad_converged) WRITE (*, *) '    - Gradient norm below tolerance.'
            WRITE (*, *) ' '
         END IF
      END IF

   END SUBROUTINE check_convergence

end module misfit_mod
