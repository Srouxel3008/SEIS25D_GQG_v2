module gradient_mod
   use Frechet_mod
   use misfit_mod
   use shared_mod
   use grid_mod
   use output_mod
   use scalers_mod
   USE gridtype_mod
   use iso_fortran_env, ONLY: dp => real64

   implicit none

contains

   SUBROUTINE compute_gradient_scaled( &
      FRECHET, RESID, NPAR, F_SCALE, PAR_SCALE, BALANCE_SCALE, GMASK, ND, NM, &
      MODEL_CURR, MODEL_REF, REG_LAMBDA, &
      GRADr, GRAD_scaled, &
      GRADr_NORM, GRAD_scaled_NORM, &
      NBLOCK, INVP, &
      USE_PRECOND, my_rank, &
      HESSDI_scaled, GRADsc_precond, &
      GRADsc_Precond_NORM, DEBUG_OUTPUT, &
      PARAM, FREQ, ITER, NX, NZ, IG, NTO, XTO, ZTO, IE0, IS0, &
      XMINC, XMAXC, ZMINC, ZMAXC)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! ---- required ----
      INTEGER, INTENT(IN) :: ND, NM, NPAR, NBLOCK
      COMPLEX(dp), INTENT(IN) :: FRECHET(ND, NM)     ! (ND, NM)
      COMPLEX(dp), INTENT(IN) :: RESID(:)            ! length ND (already weighted)
      REAL(dp), INTENT(IN) :: F_SCALE(:), PAR_SCALE(:), BALANCE_SCALE(:)  ! length NPAR (for optional precond, not used here)
      REAL(dp), INTENT(IN) :: GMASK(:)            ! length NM (hardzero+taper)
      INTEGER, INTENT(IN) :: INVP(NPAR)

      ! ---- regularization inputs ----
      REAL(dp), INTENT(IN) :: MODEL_CURR(:)       ! length NM
      REAL(dp), INTENT(IN) :: MODEL_REF(:)        ! length NM
      REAL(dp), INTENT(IN) :: REG_LAMBDA

      ! ---- outputs ----
      REAL(dp), INTENT(OUT) :: GRADr(:)              ! masked
      REAL(dp), INTENT(OUT) :: GRAD_scaled(:)        ! masked
      REAL(dp), INTENT(OUT) :: GRADr_NORM(:)         ! masked norms
      REAL(dp), INTENT(OUT) :: GRAD_scaled_NORM(:)   ! masked norms

      ! ---- optionals ----
      LOGICAL, INTENT(IN), OPTIONAL :: USE_PRECOND
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT
      REAL(dp), INTENT(IN), OPTIONAL :: HESSDI_scaled(:)         ! length NM
      REAL(dp), INTENT(OUT), OPTIONAL :: GRADsc_precond(:)       ! masked
      REAL(dp), INTENT(OUT), OPTIONAL :: GRADsc_Precond_NORM(:)  ! masked norms

      ! ---- required/optional dummy arguments for export ----
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: PARAM(:)
      REAL(dp), INTENT(IN)           :: FREQ
      INTEGER, INTENT(IN)           :: ITER
      INTEGER, INTENT(IN)           :: NX, NZ
      INTEGER, INTENT(IN)           :: NTO, IE0, IS0
      REAL(dp), INTENT(IN), OPTIONAL :: XMINC, XMAXC, ZMINC, ZMAXC

      TYPE(InversionGridType), INTENT(IN) :: IG
      REAL(dp), INTENT(IN)           :: XTO(:), ZTO(:)
      INTEGER, INTENT(IN)           :: my_rank

      ! ---- locals ----
      LOGICAL :: dbg, use_precond_opt, can_precond
      INTEGER :: II, IM, cs, ce
      INTEGER :: NNX_OUT, NNZ_OUT
      REAL(dp) :: sc_ia
      REAL(dp) :: DX_OUT
      REAL(dp) :: raw_min, raw_max, raw_l2

      REAL(dp), PARAMETER :: tiny_dp = TINY(1.0_dp)

      REAL(dp) :: DNRM2      ! BLAS
      CHARACTER(LEN=256) :: FNAME   ! local scratch for filenames
      REAL(dp), ALLOCATABLE :: GRADr_raw(:)

      ! ---- options ----
      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      use_precond_opt = .FALSE.
      IF (PRESENT(USE_PRECOND)) use_precond_opt = USE_PRECOND
      can_precond = use_precond_opt .AND. PRESENT(HESSDI_scaled) .AND. PRESENT(GRADsc_precond)

      ! ---- minimal sanity checks (only if dbg) ----
      IF (dbg) THEN
         IF (SIZE(RESID) /= ND) THEN
            WRITE (*, *) 'compute_gradient_scaled: SIZE(RESID)=', SIZE(RESID), ' ND=', ND
            STOP 'compute_gradient_scaled: RESID size mismatch'
         END IF

         IF (SIZE(FRECHET, 1) /= ND .OR. SIZE(FRECHET, 2) /= NM) THEN
            WRITE (*, *) 'compute_gradient_scaled: FRECHET shape=', &
               SIZE(FRECHET, 1), SIZE(FRECHET, 2), ' expected=', ND, NM
            STOP 'compute_gradient_scaled: FRECHET shape mismatch'
         END IF

         IF (SIZE(MODEL_CURR) /= NM) THEN
            WRITE (*, *) 'compute_gradient_scaled: SIZE(MODEL_CURR)=', SIZE(MODEL_CURR), ' NM=', NM
            STOP 'compute_gradient_scaled: MODEL_CURR size mismatch'
         END IF

         IF (SIZE(MODEL_REF) /= NM) THEN
            WRITE (*, *) 'compute_gradient_scaled: SIZE(MODEL_REF)=', SIZE(MODEL_REF), ' NM=', NM
            STOP 'compute_gradient_scaled: MODEL_REF size mismatch'
         END IF

         IF (NBLOCK*COUNT(INVP == 1) /= NM) THEN
            WRITE (*, *) 'compute_gradient_scaled: NBLOCK*COUNT(INVP==1)=', &
               NBLOCK*COUNT(INVP == 1), ' but NM=', NM
            STOP 'compute_gradient_scaled: NBLOCK/INVP/NM inconsistency'
         END IF

         IF (PRESENT(HESSDI_scaled)) THEN
            IF (SIZE(HESSDI_scaled) /= NM) THEN
               WRITE (*, *) 'compute_gradient_scaled: SIZE(HESSDI_scaled)=', SIZE(HESSDI_scaled), ' NM=', NM
               STOP 'compute_gradient_scaled: HESSDI_scaled size mismatch'
            END IF
         END IF

         IF (PRESENT(GRADsc_precond)) THEN
            IF (SIZE(GRADsc_precond) /= NM) THEN
               WRITE (*, *) 'compute_gradient_scaled: SIZE(GRADsc_precond)=', SIZE(GRADsc_precond), ' NM=', NM
               STOP 'compute_gradient_scaled: GRADsc_precond size mismatch'
            END IF
         END IF

         IF (PRESENT(GRADsc_Precond_NORM)) THEN
            IF (SIZE(GRADsc_Precond_NORM) /= NPAR) THEN
               WRITE (*, *) 'compute_gradient_scaled: SIZE(GRADsc_Precond_NORM)=', &
                  SIZE(GRADsc_Precond_NORM), ' NPAR=', NPAR
               STOP 'compute_gradient_scaled: GRADsc_Precond_NORM size mismatch'
            END IF
         END IF
      END IF

      ! ---- init outputs ----
      GRADr(:) = 0.0_dp
      GRAD_scaled(:) = 0.0_dp
      GRADr_NORM(:) = 0.0_dp
      GRAD_scaled_NORM(:) = 0.0_dp
      IF (PRESENT(GRADsc_precond)) GRADsc_precond(:) = 0.0_dp
      IF (PRESENT(GRADsc_Precond_NORM)) GRADsc_Precond_NORM(:) = 0.0_dp
      ALLOCATE (GRADr_raw(NM))
      GRADr_raw(:) = 0.0_dp

      ! ---- raw gradient (UNMASKED here): g_data = Re{ J^H RESID } ----
      GRADr_raw = REAL(MATMUL(CONJG(TRANSPOSE(FRECHET)), RESID), dp)
      GRADr(:) = GRADr_raw(:)

      IM = 0
      DO II = 1, NPAR
         IF (INVP(II) /= 1) CYCLE
         IM = IM + 1
         cs = (IM - 1)*NBLOCK + 1
         ce = IM*NBLOCK

         ! 1) convert physical data gradient to PAR-scaled model space
         GRADr(cs:ce) = GRADr_raw(cs:ce)*PAR_SCALE(II)

         ! 2) add Tikhonov gradient already in PAR-scaled model space
         IF (REG_LAMBDA > 0.0_dp) THEN
            GRADr(cs:ce) = GRADr(cs:ce) &
                           + REG_LAMBDA/REAL(NM, dp)*(MODEL_CURR(cs:ce) - MODEL_REF(cs:ce))
         END IF

         ! 3) apply Frechet sensitivity scaling and Q/parameter balance
         GRAD_scaled(cs:ce) = GRADr(cs:ce)*F_SCALE(II)*BALANCE_SCALE(II)

         ! 4) inverse Hessian preconditioning, if requested
         ! IF (can_precond) THEN
            GRADsc_precond(cs:ce) = HESSDI_scaled(cs:ce)*GRAD_scaled(cs:ce)
         ! END IF
      END DO

      ! ---- apply GMASK once (hardzero+taper) ----
      GRADr(:) = GRADr(:)*GMASK(:)
      GRAD_scaled(:) = GRAD_scaled(:)*GMASK(:)
      GRADsc_precond(:) = GRADsc_precond(:)*GMASK(:)

      ! ---- masked norms per parameter-block ----
      IM = 0
      DO II = 1, NPAR
         IF (INVP(II) /= 1) CYCLE
         IM = IM + 1
         cs = (IM - 1)*NBLOCK + 1
         ce = IM*NBLOCK

         GRADr_NORM(II) = DNRM2(NBLOCK, GRADr(cs:ce), 1)
         GRAD_scaled_NORM(II) = DNRM2(NBLOCK, GRAD_scaled(cs:ce), 1)
         if (my_rank == 0) then
            WRITE (*, '(A, I3, A, 1ES12.4)') 'Grad block ', II, ' GRAD_SCALED MORM = ', GRAD_scaled_NORM(II)
         END IF

         ! IF (can_precond .AND. PRESENT(GRADsc_Precond_NORM)) THEN
            GRADsc_Precond_NORM(II) = DNRM2(NBLOCK, GRADsc_precond(cs:ce), 1)
         ! END IF
      END DO

      ! 6) optional exports (always masked exactly once)
      !    When canonical physical bounds are provided, export to a fixed grid
      !    (same extent and point locations across bands). Otherwise, fall back
      !    to the legacy grid-resolved export.
      IF (dbg .AND. my_rank == 0 .AND. PRESENT(PARAM) .AND. ITER < 3) THEN
         IM = 0
         IF (PRESENT(XMINC) .AND. PRESENT(XMAXC) .AND. PRESENT(ZMINC) .AND. PRESENT(ZMAXC)) THEN
            NNX_OUT = NTO
            IF (NNX_OUT < 2) RETURN
            DX_OUT = (XMAXC - XMINC)/REAL(NNX_OUT - 1, dp)
            NNZ_OUT = 1 + NINT((ZMAXC - ZMINC)/DX_OUT)
            IF (NNZ_OUT < 2) NNZ_OUT = 2
         ELSE
            NNX_OUT = 0
            NNZ_OUT = 0
         END IF

         DO II = 1, NPAR
            IF (INVP(II) /= 1) CYCLE
            IM = IM + 1
            cs = (IM - 1)*NBLOCK + 1
            ce = IM*NBLOCK

            raw_min = MINVAL(GRADr_raw(cs:ce))
            raw_max = MAXVAL(GRADr_raw(cs:ce))
            raw_l2 = DNRM2(NBLOCK, GRADr_raw(cs:ce), 1)
            WRITE (*, '(A,I3,2X,A,A,2X,A,1ES12.4,2X,A,1ES12.4,2X,A,1ES12.4)') &
               'Raw grad block', II, 'PAR=', TRIM(PARAM(II)), &
               'min=', raw_min, 'max=', raw_max, 'L2=', raw_l2

            ! CALL CFNAME_GRADIENT('GRADraw_', PARAM(II), FREQ, ITER, '.dat', FNAME)
            ! IF (NNX_OUT > 0) THEN
            !    CALL GRID2D_OUT_FIXED(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADr_raw(cs:ce), &
            !                          NTO, XTO, ZTO, XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
            ! ELSE
            !    CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADr_raw(cs:ce), NTO, XTO, ZTO, IE0, IS0)
            ! END IF

            CALL CFNAME_GRADIENT('GRADr_', PARAM(II), FREQ, ITER, '.dat', FNAME)
            IF (NNX_OUT > 0) THEN
               CALL GRID2D_OUT_FIXED(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADr(cs:ce), &
                                     NTO, XTO, ZTO, XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
            ELSE
               CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADr(cs:ce), NTO, XTO, ZTO, IE0, IS0)
            END IF

            CALL CFNAME_GRADIENT('GRADrsc_', PARAM(II), FREQ, ITER, '.dat', FNAME)
            IF (NNX_OUT > 0) THEN
               CALL GRID2D_OUT_FIXED(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRAD_scaled(cs:ce), &
                                     NTO, XTO, ZTO, XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
            ELSE
               CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRAD_scaled(cs:ce), NTO, XTO, ZTO, IE0, IS0)
            END IF

            IF (can_precond) THEN
               CALL CFNAME_GRADIENT('GRADrpm_', PARAM(II), FREQ, ITER, '.dat', FNAME)
               IF (NNX_OUT > 0) THEN
                  CALL GRID2D_OUT_FIXED(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADsc_precond(cs:ce), &
                                        NTO, XTO, ZTO, XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
               ELSE
                  CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADsc_precond(cs:ce), NTO, XTO, ZTO, IE0, IS0)
               END IF
            END IF
         END DO
      END IF

      IF (ALLOCATED(GRADr_raw)) DEALLOCATE (GRADr_raw)

   END SUBROUTINE compute_gradient_scaled

   SUBROUTINE ComputeHessDI(FRECHET, ND, IM, NBLOCK, NM, NPAR, INVP, &
                            GMASK, F_SCALE, PAR_SCALE, &
                            HESS, HESSDI_scaled, REG_LAMBDA, &
                            WD_acq, WD_amp, BETA, &
                            NX, NZ, IG, FREQ, CREF, my_rank, &
                            DEBUG_OUTPUT, USE_LM, LAMBDA, GLOBAL_MAXP)

      IMPLICIT NONE

      ! ---- required ----
      INTEGER, INTENT(IN)        :: ND, IM, NBLOCK, NM, NPAR
      INTEGER, INTENT(IN)        :: NX, NZ, my_rank
      INTEGER, INTENT(IN)        :: INVP(:)
      COMPLEX(dp), INTENT(IN)    :: FRECHET(ND, NM)
      REAL(dp), INTENT(IN)    :: REG_LAMBDA
      REAL(dp), INTENT(IN)    :: GMASK(:)        ! (NM)
      REAL(dp), INTENT(IN)    :: F_SCALE(:)    ! per-IA (NPAR)
      REAL(dp), INTENT(IN)    :: PAR_SCALE(:)    ! per-IA (NPAR)
      REAL(dp), INTENT(IN)    :: FREQ, CREF
      TYPE(InversionGridType), INTENT(IN) :: IG
      REAL(dp), INTENT(OUT)   :: HESS(:)         ! diag(Jz^H Wd Jz) (scaled/weighted)
      REAL(dp), INTENT(OUT)   :: HESSDI_scaled(:) ! (diag + damp)^{-1} ⊙ GMASK
      REAL(dp), INTENT(IN)    :: WD_acq(:)       ! per-trace weights (ND)
      REAL(dp), INTENT(IN)    :: WD_amp          ! global amplitude scaler
      REAL(dp), INTENT(IN)    :: BETA            ! damping strength (β·maxP)

      ! ---- optional ----
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT
      LOGICAL, OPTIONAL, INTENT(IN) :: USE_LM
      REAL(dp), OPTIONAL, INTENT(IN):: LAMBDA
      LOGICAL, OPTIONAL, INTENT(IN) :: GLOBAL_MAXP
      REAL(dp) :: Hphys, Hphys_F2
      ! ---- locals ----
      LOGICAL :: dbg, do_lm, do_global
      INTEGER :: ii, jj, ip, i, j, np, nvals, idx
      REAL(dp) :: wamp2, s_ia, tsum, maxP_all, maxP_param, denom
      REAL(dp) :: Hmin_F1, Hmax_F1, Hmin_F2, Hmax_F2
      REAL(dp), ALLOCATABLE :: H_work(:)
      ! ---- size checks ----
      IF (SIZE(FRECHET, 2) /= NM) STOP 'ComputeHessDI: FRECHET second dim must be NM'
      IF (SIZE(FRECHET, 1) < ND) STOP 'ComputeHessDI: FRECHET first dim < ND'
      IF (SIZE(GMASK) /= NM) STOP 'ComputeHessDI: GMASK size must be NM'
      IF (SIZE(HESS) /= NM) STOP 'ComputeHessDI: HESS size must be NM'
      IF (SIZE(HESSDI_scaled) /= NM) STOP 'ComputeHessDI: HESSDI_scaled size must be NM'
      IF (SIZE(WD_acq) < ND) STOP 'ComputeHessDI: WD_acq size < ND'

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      do_lm = .FALSE.; IF (PRESENT(USE_LM)) do_lm = USE_LM
      do_global = .TRUE.; IF (PRESENT(GLOBAL_MAXP)) do_global = GLOBAL_MAXP

      ! effective amplitude scaling factor squared

      ! initialise outputs
      HESS(:) = 0.0_dp
      HESSDI_scaled(:) = 0.0_dp

      ! -----------------------------------------------------------
! Pass 1: HESS = diag(Jz^H Wd Jz)
!
! One-F0 version returned:
!   HESS = PAR_SCALE^2 * F_SCALE * Hphys
!
! F0^2 version only printed for comparison:
!   HESS_F2 = PAR_SCALE^2 * F_SCALE^2 * Hphys
! -----------------------------------------------------------
      ip = 0
      DO ii = 1, NPAR

         IF (INVP(ii) /= 1) CYCLE

         ip = ip + 1
         np = (ip - 1)*NBLOCK

         s_ia = MAX(PAR_SCALE(ii), tiny_dp)

         Hmin_F1 = HUGE(1.0_dp)
         Hmax_F1 = -HUGE(1.0_dp)
         Hmin_F2 = HUGE(1.0_dp)
         Hmax_F2 = -HUGE(1.0_dp)

         DO j = 1, NBLOCK

            tsum = 0.0_dp

            !$OMP SIMD REDUCTION(+:tsum)
            DO i = 1, ND
               tsum = tsum + WD_acq(i)*Wd_amp**2* &
                      (DREAL(FRECHET(i, np + j))**2 + &
                       DIMAG(FRECHET(i, np + j))**2) !
            END DO

            ! Physical weighted Hessian
            Hphys = tsum

            ! Returned version: single-F0 conditioner
            ! HESS(np + j) = s_ia*s_ia*F_SCALE(ii)*F_SCALE(ii)*Hphys + &
            !                REG_LAMBDA/REAL(IM*NBLOCK, dp)
            HESS(np + j) = s_ia*s_ia*F_SCALE(ii)*F_SCALE(ii)*Hphys + &
                           REG_LAMBDA/REAL(IM*NBLOCK, dp)

            ! Diagnostic only: mathematically consistent F0^2
            Hphys_F2 = s_ia*s_ia*F_SCALE(ii)*F_SCALE(ii)*Hphys + &
                       REG_LAMBDA/REAL(IM*NBLOCK, dp)

            Hmin_F1 = MIN(Hmin_F1, HESS(np + j))
            Hmax_F1 = MAX(Hmax_F1, HESS(np + j))

            Hmin_F2 = MIN(Hmin_F2, Hphys_F2)
            Hmax_F2 = MAX(Hmax_F2, Hphys_F2)

         END DO

         IF (dbg.and.my_rank == 0) THEN
            WRITE (*, '(A,I3,2X,A,1PE12.5,2X,A,1PE12.5)') &
               'IA=', ii, 'H(single F0) min=', Hmin_F1, 'max=', Hmax_F1

            WRITE (*, '(A,I3,2X,A,1PE12.5,2X,A,1PE12.5)') &
               'IA=', ii, 'H(F0^2)     min=', Hmin_F2, 'max=', Hmax_F2
         END IF
      END DO

      ! -----------------------------------------------------------
      ! Pass 2: damping and inverse -> HESSDI_scaled
      ! -----------------------------------------------------------
      IF (.NOT. do_lm) THEN
         IF (do_global) THEN
            maxP_all = MAXVAL(HESS*GMASK)
            IF (maxP_all < tiny_dp) maxP_all = 0.0_dp
         ELSE
            maxP_all = 0.0_dp
         END IF
      ELSE
         maxP_all = 0.0_dp
      END IF

      ip = 0
      DO ii = 1, NPAR
         IF (INVP(ii) /= 1) CYCLE
         ip = ip + 1
         np = (ip - 1)*NBLOCK

         IF (.NOT. do_lm .AND. .NOT. do_global) THEN
            ALLOCATE (H_work(NBLOCK))
            H_work(1:NBLOCK) = HESS(np + 1:np + NBLOCK)!*GMASK(np + 1:np + NBLOCK)
            nvals = COUNT(H_work > tiny_dp)
            IF (nvals > 0) THEN
               H_work = PACK(H_work, H_work > tiny_dp)
               CALL quicksort_real_dp(H_work, 1, nvals)
               idx = MAX(1, CEILING(0.95_dp*REAL(nvals, dp)))
               maxP_param = H_work(idx)
            ELSE
               maxP_param = 0.0_dp
            END IF
            IF (maxP_param < tiny_dp) maxP_param = 0.0_dp
            WRITE (*, '(A,1PE12.5,1PE12.5)') 'ComputeHessDI: using param-wise 95th pct H=', maxP_param, BETA*maxP_param
            DEALLOCATE (H_work)
         ELSE
            maxP_param = maxP_all
         END IF

         DO j = 1, NBLOCK
            IF (do_lm) THEN
               denom = HESS(np + j) + MAX(0.0_dp, MERGE(LAMBDA, 0.0_dp, PRESENT(LAMBDA)))
            ELSE
               denom = HESS(np + j) + MAX(0.0_dp, BETA)*maxP_param
            END IF

            HESSDI_scaled(np + j) = (1.0_dp/denom)*GMASK(np + j)
         END DO

         IF (dbg.and.my_rank == 0) THEN
            WRITE (*, '(A,I4,2X,A,1PE12.5,2X,A,1PE12.5)') 'HESSD: ia=', ii, &
               'min=', MINVAL(HESSDI_scaled(np + 1:np + NBLOCK)), &
               'max=', MAXVAL(HESSDI_scaled(np + 1:np + NBLOCK))
         END IF
      END DO
   END SUBROUTINE ComputeHessDI

   SUBROUTINE smooth_gradient_blocks( &
      GRAD_in, GMASK, &
      NBLOCK, INVP, NPAR, &
      NX, NZ, IG, &
      FREQ, CREF, my_rank, &
      USE_GR_SMOOTH, USE_PRECOND, HESSDI_scaled, GRADsc_precond, GRADsc_Precond_NORM, &
      PARAM, ITER, NTO, XTO, ZTO, IE0, IS0, &
      XMINC, XMAXC, ZMINC, ZMAXC, DEBUG_OUTPUT, &
      SIGMA_X, SIGMA_Z, &
      GRAD_out)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! ---- required ----
      INTEGER, INTENT(IN) :: NBLOCK, NPAR
      INTEGER, INTENT(IN) :: NX, NZ, my_rank
      REAL(dp), INTENT(IN) :: GRAD_in(:)       ! length NM
      REAL(dp), INTENT(IN) :: GMASK(:)         ! length NM
      INTEGER, INTENT(IN) :: INVP(NPAR)
      REAL(dp), INTENT(IN) :: FREQ, CREF
      TYPE(InversionGridType), INTENT(IN) :: IG
      LOGICAL, INTENT(IN) :: USE_GR_SMOOTH
      LOGICAL, INTENT(IN), OPTIONAL :: USE_PRECOND
      REAL(dp), INTENT(IN), OPTIONAL :: HESSDI_scaled(:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: GRADsc_precond(:)
      REAL(dp), INTENT(OUT), OPTIONAL :: GRADsc_Precond_NORM(:)
      REAL(dp), INTENT(IN), OPTIONAL :: SIGMA_X, SIGMA_Z
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: PARAM(:)
      INTEGER, INTENT(IN), OPTIONAL :: ITER
      INTEGER, INTENT(IN), OPTIONAL :: NTO, IE0, IS0
      REAL(dp), INTENT(IN), OPTIONAL :: XTO(:), ZTO(:)
      REAL(dp), INTENT(IN), OPTIONAL :: XMINC, XMAXC, ZMINC, ZMAXC
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT

      ! ---- output ----
      REAL(dp), INTENT(OUT) :: GRAD_out(:)     ! same size as GRAD_in, masked

      ! ---- locals ----
      INTEGER :: II, IM, cs, ce
      INTEGER :: NNX_OUT, NNZ_OUT
      REAL(dp), ALLOCATABLE :: GZ_work(:)
      REAL(dp) :: DX_OUT
      CHARACTER(LEN=256) :: FNAME
      LOGICAL :: dbg, use_precond_opt, can_precond
      REAL(dp) :: DNRM2

      ! ---- init ----
      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      use_precond_opt = .FALSE.
      IF (PRESENT(USE_PRECOND)) use_precond_opt = USE_PRECOND
      can_precond = use_precond_opt .AND. PRESENT(HESSDI_scaled) .AND. PRESENT(GRADsc_precond)
      GRAD_out(:) = GRAD_in(:)
      IF (PRESENT(GRADsc_Precond_NORM)) GRADsc_Precond_NORM(:) = 0.0_dp

      IF (dbg) THEN
         IF (PRESENT(HESSDI_scaled)) THEN
            IF (SIZE(HESSDI_scaled) /= SIZE(GRAD_in)) THEN
               WRITE (*, *) 'smooth_gradient_blocks: SIZE(HESSDI_scaled)=', SIZE(HESSDI_scaled), ' NM=', SIZE(GRAD_in)
               STOP 'smooth_gradient_blocks: HESSDI_scaled size mismatch'
            END IF
         END IF

         IF (PRESENT(GRADsc_precond)) THEN
            IF (SIZE(GRADsc_precond) /= SIZE(GRAD_in)) THEN
               WRITE (*, *) 'smooth_gradient_blocks: SIZE(GRADsc_precond)=', SIZE(GRADsc_precond), ' NM=', SIZE(GRAD_in)
               STOP 'smooth_gradient_blocks: GRADsc_precond size mismatch'
            END IF
         END IF

         IF (PRESENT(GRADsc_Precond_NORM)) THEN
            IF (SIZE(GRADsc_Precond_NORM) /= NPAR) THEN
               WRITE (*, *) 'smooth_gradient_blocks: SIZE(GRADsc_Precond_NORM)=', SIZE(GRADsc_Precond_NORM), ' NPAR=', NPAR
               STOP 'smooth_gradient_blocks: GRADsc_Precond_NORM size mismatch'
            END IF
         END IF
      END IF

      IF (USE_GR_SMOOTH) THEN
         ALLOCATE (GZ_work(NBLOCK))

         IM = 0
         DO II = 1, NPAR
            IF (INVP(II) /= 1) CYCLE
            IM = IM + 1
            cs = (IM - 1)*NBLOCK + 1
            ce = IM*NBLOCK

            GZ_work(1:NBLOCK) = GRAD_in(cs:ce)

            CALL gradient_smooth_gauss_block( &
               G=GZ_work, &
               GMASK=GMASK(cs:ce), &  ! boundary/taper guidance
               NX=NX, &
               NZ=NZ, &
               DX_BLOCK=IG%DX_BLOCK, &
               ZBC=IG%ZBC, &
               FREQ=FREQ, &
               CREF=CREF, &
               my_rank=my_rank, &
               SIGMA_X=SIGMA_X, SIGMA_Z=SIGMA_Z)

            GRAD_out(cs:ce) = GZ_work(1:NBLOCK)
         END DO

         IF (can_precond) THEN
            IM = 0
            DO II = 1, NPAR
               IF (INVP(II) /= 1) CYCLE
               IM = IM + 1
               cs = (IM - 1)*NBLOCK + 1
               ce = IM*NBLOCK

               CALL gradient_smooth_gauss_block( &
                  G=GRADsc_precond(cs:ce), &
                  GMASK=GMASK(cs:ce), &  ! boundary/taper guidance
                  NX=NX, &
                  NZ=NZ, &
                  DX_BLOCK=IG%DX_BLOCK, &
                  ZBC=IG%ZBC, &
                  FREQ=FREQ, &
                  CREF=CREF, &
                  my_rank=my_rank, &
                  SIGMA_X=SIGMA_X, SIGMA_Z=SIGMA_Z)
            END DO
         END IF

         DEALLOCATE (GZ_work)
      END IF

      ! Apply GMASK once (to enforce hard zero + taper)
      GRAD_out(:) = GRAD_out(:)*GMASK(:)

      IF (can_precond) THEN
         GRADsc_precond(:) = GRADsc_precond(:)*GMASK(:)
         IF (PRESENT(GRADsc_Precond_NORM)) THEN
            IM = 0
            DO II = 1, NPAR
               IF (INVP(II) /= 1) CYCLE
               IM = IM + 1
               cs = (IM - 1)*NBLOCK + 1
               ce = IM*NBLOCK

               GRADsc_Precond_NORM(II) = DNRM2(NBLOCK, GRADsc_precond(cs:ce), 1)
               ! GRADsc_precond(cs:ce) = GRADsc_precond(cs:ce)/MAX(GRADsc_Precond_NORM(II), tiny_dp)
            END DO
         END IF
      END IF

      IF (dbg .AND. my_rank == 0 .AND. PRESENT(PARAM) .AND. PRESENT(ITER) .AND. &
          PRESENT(NTO) .AND. PRESENT(XTO) .AND. PRESENT(ZTO) .AND. &
          PRESENT(IE0) .AND. PRESENT(IS0)) THEN
         IF (ITER < 4) THEN
            IF (PRESENT(XMINC) .AND. PRESENT(XMAXC) .AND. PRESENT(ZMINC) .AND. PRESENT(ZMAXC)) THEN
               NNX_OUT = NTO
               IF (NNX_OUT < 2) RETURN
               DX_OUT = (XMAXC - XMINC)/REAL(NNX_OUT - 1, dp)
               NNZ_OUT = 1 + NINT((ZMAXC - ZMINC)/DX_OUT)
               IF (NNZ_OUT < 2) NNZ_OUT = 2
            ELSE
               NNX_OUT = 0
               NNZ_OUT = 0
            END IF

            IM = 0
            DO II = 1, NPAR
               IF (INVP(II) /= 1) CYCLE
               IM = IM + 1
               cs = (IM - 1)*NBLOCK + 1
               ce = IM*NBLOCK

               CALL CFNAME_GRADIENT('GRADsm_', PARAM(II), FREQ, ITER, '.dat', FNAME)
               IF (NNX_OUT > 0) THEN
                  CALL GRID2D_OUT_FIXED(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRAD_out(cs:ce), &
                                        NTO, XTO, ZTO, XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
               ELSE
                  CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRAD_out(cs:ce), NTO, XTO, ZTO, IE0, IS0)
               END IF

               IF (can_precond) THEN
                  CALL CFNAME_GRADIENT('GRADrpm_', PARAM(II), FREQ, ITER, '.dat', FNAME)
                  IF (NNX_OUT > 0) THEN
                     CALL GRID2D_OUT_FIXED(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADsc_precond(cs:ce), &
                                           NTO, XTO, ZTO, XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
                  ELSE
                     CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, GRADsc_precond(cs:ce), NTO, XTO, ZTO, IE0, IS0)
                  END IF
               END IF
            END DO
         END IF
      END IF
   END SUBROUTINE smooth_gradient_blocks

   SUBROUTINE gradient_mask(GMASK, NX, NZ, IE0, IS0, IZ, DZ, IG, NTO, XTO, ZTO, my_rank, DEBUG_OUTPUT)

!----------------------------------------------------------------------C
!                                                                      C
!     gradient_mask generates a binary box-shaped mask for the         C
!     gradient vector. The mask is applied to suppress values in the   C
!     side and bottom padding zones.               C
!                                                                      C
!     Entries:                                                         C
!       GMASK(*)...............Flattened Gradient scaling array ;     C
!       NX,NZ...................Coarse Grid size in X and Z directions;C
!       IE0.....................Number of grid points for side PML;    C
!       IS0.....................Number of absorbing layer (1 or  2);   C
!       IZ...................... Z-padding around the source depth;    C
!       my_rank.................MPI rank (used for optional logging);  C
!                                                                      C
!     Return:                                                          C
!       GMASK(*)...............Updated with 1 in active zone, 0 else;
!                 ISO=1                    ISO=2
!         +--------------------+     +--------------------+
!                                    | 0 |0 0 0 0 0| 0 |←Top PML
!         | 0 |0 0 0 0 0| 0 |        | 0 |0 0 0 0 0| 0 |← Source zone
!         | 0 |1 1 1 1 1| 0 |        | 0 |1 1 1 1 1| 0 |
!         | 0 |1 1 1 1 1| 0 |        | 0 |1 1 1 1 1| 0 |
!         | 0 |1 1 1 1 1| 0 |        | 0 |1 1 1 1 1| 0 |
!         | 0 |1 1 1 1 1| 0 |        | 0 |1 1 1 1 1| 0 |
!         | 0 |0 0 0 0 0| 0 |        | 0 |0 0 0 0 0| 0 |← Bottom PML
!         +--------------------+     +--------------------+
!                                                                      C
!----------------------------------------------------------------------C
      IMPLICIT NONE
      TYPE(InversionGridType), INTENT(IN) :: IG
      INTEGER, INTENT(IN) :: NX, NZ, IE0, IS0, IZ, NTO, my_rank
      REAL(dp), INTENT(IN) :: DZ
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:)
      REAL(dp), INTENT(INOUT) :: GMASK(:)    ! flattened mask (length (NX-1)*(NZ-1))
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT

      INTEGER :: K, L, LK, unit_out
      INTEGER :: block_count, NBLOCK, nstack, istack, ps, pe
      LOGICAL :: dbg
      REAL(dp) :: zsurf

      block_count = 0
      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      NBLOCK = (NX - 1)*(NZ - 1)
      IF (NBLOCK <= 0) RETURN
      IF (MOD(SIZE(GMASK), NBLOCK) /= 0) THEN
         WRITE (*, *) 'gradient_mask: SIZE(GMASK)=', SIZE(GMASK), ' incompatible with NBLOCK=', NBLOCK
         STOP 'gradient_mask: invalid stacked mask size'
      END IF
      nstack = SIZE(GMASK)/NBLOCK

      DO L = 1, NX - 1
         zsurf = ZH(NTO, XTO, ZTO, IG%XBC(L))
         DO K = 1, NZ - 1
            block_count = block_count + 1
            LK = (L - 1)*(NZ - 1) + K
            IF ((K > IE0) .AND. (K < NZ - 1 - (IS0 - 1)*IE0 - IZ) .AND. &
                (L > IE0) .AND. (L <= NX - 1 - IE0) .AND. &
                (.NOT. IG%vPML(LK)) .AND. (IG%ZBC(LK) <= zsurf)) THEN
               GMASK(LK) = 1.0_dp
            ELSE
               GMASK(LK) = 0.0_dp
            END IF
         END DO
      END DO

      DO istack = 2, nstack
         ps = (istack - 1)*NBLOCK + 1
         pe = istack*NBLOCK
         GMASK(ps:pe) = GMASK(1:NBLOCK)
      END DO

      IF (dbg .AND. my_rank == 0) THEN
         OPEN (NEWUNIT=unit_out, FILE='out.Gmask.txt', STATUS='REPLACE', ACTION='WRITE')
         DO K = 1, NZ - 1
            DO L = 1, NX - 1
               LK = (L - 1)*(NZ - 1) + K
               WRITE (unit_out, '(F6.2)', ADVANCE='NO') GMASK(LK)
               IF (L < NX - 1) WRITE (unit_out, '(A)', ADVANCE='NO') ' '
            END DO
            WRITE (unit_out, *)
         END DO
         CLOSE (unit_out)
      END IF

   END SUBROUTINE gradient_mask

   !----------------------------------------------------------------------
   !
   !     gradient_mask_taper applies a Gaussian-like taper vertically in
   !     the active inversion zone. Zones outside this are zeroed.
   !
   !     The taper is applied from the top of the active zone downward.
   !     The masking follows the same logic as gradient_mask.
   !
   !     Entries:
   !       GMASKSM(*)...............1D Gradient scaling array (output)
   !       NX,NZ...................Grid dimensions in X and Z
   !       IE0.....................PML width
   !       IS0.....................Absorption scheme (1 or 2)
   !       IZ......................Depth taper at bottom
   !       MZ0.....................Number of points for top taper zone
   !       DZ......................Grid spacing in Z (used by taper)
   !       n_z1, n_z2..................Taper shape parameters passed to Z_taper
   !       my_rank.................MPI rank (for debug output)
   !       DEBUG_OUTPUT............Logical flag: if true, writes mask to file
   !
   !     Returns:
   !       GMASKSM..................Smooth mask (tapered at top, 1 inside, 0 outside)
   !     Visual layout (IE0 = 1, MZ0 = 4):
   !                 IS0 = 1                    IS0 = 2
   !         +--------------------+      +--------------------+
   !         | 0 |p0 p0 p0 p0 p0| 0 |    | 0 |0  0  0  0  0 | 0 |← Top PML
   !         | 0 |p1 p1 p1 p1 p1| 0 |    | 0 |p0 p0 p0 p0 p0| 0 |
   !         | 0 |p2 p2 p2 p2 p2| 0 |    | 0 |p1 p1 p1 p1 p1| 0 |
   !         | 0 |1  1  1  1  1 | 0 |    | 0 |p2 p2 p2 p2 p2| 0 |
   !         | 0 |0  0  0  0  0 | 0 |    | 0 |1  1  1  1  1 | 0 |
   !                                     | 0 |0  0  0  0  0 | 0 | ← Bottom PML
   !         +--------------------+      +--------------------+
   !     Note: 'p' represents taper values between 0 and 1 from PRE(:)
   !
   !----------------------------------------------------------------------
   SUBROUTINE gradient_mask_taper(NX, NZ, IE0, IS0, my_rank, NM, MZ0, IZ, &
                                  GZ1, GZ2, DZ, IG, NTO, XTO, ZTO, GMASK, DEBUG_OUTPUT)
      IMPLICIT NONE

      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG_OUTPUT
      TYPE(InversionGridType), INTENT(IN) :: IG
      INTEGER, INTENT(IN) :: NX, NZ, IE0, IS0, my_rank, NM, MZ0, IZ, NTO
      REAL(dp), INTENT(IN) :: GZ1, GZ2      ! <-- depths in meters
      REAL(dp), INTENT(IN) :: DZ
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:)
      REAL(dp), INTENT(INOUT) :: GMASK(:)

      LOGICAL :: dbg
      INTEGER :: K, L, LK, unit_out
      INTEGER :: n_z1, n_z2         ! <-- internal block indices
      INTEGER :: k_from_top
      INTEGER :: NBLOCK, nstack, istack, ps, pe
      REAL(dp), ALLOCATABLE :: PRE(:)
      REAL(dp) :: zsurf

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

      ! Convert depths in meters to block indices (1-based)
      n_z1 = NINT(GZ1/DZ)
      n_z2 = NINT(GZ2/DZ)

      ! (Optional) simple safety guards
      IF (n_z1 < 1) n_z1 = 1
      IF (n_z2 <= n_z1) n_z2 = n_z1 + 1
      IF (n_z2 > MZ0) n_z2 = MZ0

      ! Derived dimensions
      NBLOCK = (NX - 1)*(NZ - 1)
      IF (NBLOCK <= 0) RETURN
      IF (MOD(NM, NBLOCK) /= 0) THEN
         WRITE (*, *) 'gradient_mask_taper: NM=', NM, ' incompatible with NBLOCK=', NBLOCK
         STOP 'gradient_mask_taper: invalid stacked mask size'
      END IF
      nstack = NM/NBLOCK
      ALLOCATE (PRE(MZ0))
      GMASK(1:NM) = 0.0_dp

      ! Build a taper indexed by distance below the topographic surface.
      CALL Z_taper(MZ0, DZ, n_z1, n_z2, PRE)

      DO L = 1, NX - 1
         zsurf = ZH(NTO, XTO, ZTO, IG%XBC(L))
         DO K = 1, NZ - 1
            LK = (L - 1)*(NZ - 1) + K
            IF ((K > IE0) .AND. (K < NZ - 1 - (IS0 - 1)*IE0 - IZ) .AND. &
                (L > IE0) .AND. (L <= NX - 1 - IE0) .AND. &
                (.NOT. IG%vPML(LK)) .AND. (IG%ZBC(LK) <= zsurf)) THEN
               k_from_top = NZ - 1 - (IS0 - 1)*IE0 - IZ - K
               IF (k_from_top >= 1 .AND. k_from_top <= MZ0) THEN
                  GMASK(LK) = PRE(k_from_top)
               ELSE
                  GMASK(LK) = 1.0_dp
               END IF
            ELSE
               GMASK(LK) = 0.0_dp
            END IF
         END DO
      END DO

      DO istack = 2, nstack
         ps = (istack - 1)*NBLOCK + 1
         pe = istack*NBLOCK
         GMASK(ps:pe) = GMASK(1:NBLOCK)
      END DO

      IF (dbg .AND. my_rank == 0) THEN
         OPEN (NEWUNIT=unit_out, FILE='out_GmaskTap.txt', STATUS='REPLACE', ACTION='WRITE')
         DO K = 1, NZ - 1
            DO L = 1, NX - 1
               LK = (L - 1)*(NZ - 1) + K
               WRITE (unit_out, '(F6.2)', ADVANCE='NO') GMASK(LK)
               IF (L < NX - 1) WRITE (unit_out, '(A)', ADVANCE='NO') ' '
            END DO
            WRITE (unit_out, *)
         END DO
         CLOSE (unit_out)
      END IF

      DEALLOCATE (PRE)

   END SUBROUTINE gradient_mask_taper

   SUBROUTINE Z_taper(nz, dz, n_z1, n_z2, P)
!----------------------------------------------------------------------
      !
!     Z_taper defines a 1D vertical tapering function P(z)
!     over a depth vector defined by nz and dz.
!     suppress artifacts near the sources, and noisy near surface
!     Create smoother gradient transitions.
!
!     It divides the vertical domain into 3 regions:
!     1. For z ≤ z1 = n_z1 * dz:
!        - P(z) = 0.0 (fully suppressed, typically source zone and water layer)
!     2. For z1 < z ≤ z2 = n_z2 * dz:
!        - Gaussian-like transition: downweight the impact of near surface
!          P(z) = exp(-0.5 * ((2 * n_z1 * (z - z1 - delta1) / delta1)^2))
!     3. For z > z2:
!        - P(z) increases linearly: P(z) = 0.5 + z/z2
!
!     Entries:
!       nz..................Number of vertical grid points to be tapered
!       dz..................Vertical grid spacing
!       n_z1.........Controls steepness of Gaussian taper
!       n_z2.........Defines width and exponent of taper region
!
!     Return:
!       P(nz)...............Taper vector (output)
!
!----------------------------------------------------------------------
      IMPLICIT real(dp) (A - H, O - Z)

      INTEGER, INTENT(IN)           :: nz, n_z1, n_z2
      real(dp), INTENT(IN)          :: dz
      real(dp), INTENT(OUT)         :: P(nz)

      INTEGER                       :: iz
      real(dp)                      :: a, delta1, zs1, zs2, zg, temp

      a = 3.0_dp
      delta1 = (n_z2 - n_z1)*dz
      zs1 = n_z1*dz
      zs2 = n_z2*dz

      DO iz = 1, nz
         zg = iz*dz
         IF ((zg .gt. 0.0D0) .AND. (zg .le. zs1)) THEN
            P(iz) = 0.0D0
         ELSE IF ((zg .gt. zs1) .AND. (zg .le. zs2)) THEN
            temp = (zg - zs1)/delta1
            P(iz) = EXP(-0.5D0*(a*(1.0D0 - temp))**2)
         ELSE IF (zg .GT. zs2) THEN
            P(iz) = 1.0D0
         END IF
      END DO

      RETURN
   END SUBROUTINE Z_taper

!----------------------------------------------------------------------C
!  gradient_smooth_gauss_block                                        C
!                                                                      C
!  Frequency-dependent Gaussian smoothing of the gradient on the      C
!  (NX-1) x (NZ-1) block grid.                                        C
!                                                                      C
!  - Uses DX_BLOCK(:) and IG%ZBC(:) to estimate effective horizontal      C
!    and vertical block spacings.                                      C
!  - Computes SPANX/SPANZ from wavelength: λ_ref = CREF / FREQ.        C
!  - Applies separable 1D Gaussian smoothing (Z then X) with           C
!    local renormalisation and GMASK-based exclusion of PML / water.   C
!                                                                      C
!  Entries:                                                            C
!    G(*)...........Gradient vector on block grid (length (NX-1)*(NZ-1))C
!    GMASK(*).......Mask/taper on same grid (0 in PML/water, 1 inside) C
!    NX,NZ..........Coarse grid size in X and Z (blocks + 1)           C
!    DX_BLOCK(*)....Per-block horizontal size (from ComputeBlockCenters)C
!    IG%ZBC(*).........Per-block depth of block centres (length NBLOCK)   C
!    FREQ...........Current frequency (Hz)                             C
!    CREF...........Reference velocity for wavelength (e.g. max Vp)    C
!    my_rank........MPI rank (for optional debug)                      C
!                                                                      C
!  Return:                                                             C
!    G(*)...........Smoothed gradient (in-place)                       C
!----------------------------------------------------------------------C
   SUBROUTINE gradient_smooth_gauss_block(G, GMASK, NX, NZ, &
                                          DX_BLOCK, ZBC, FREQ, CREF, &
                                          my_rank, &
                                          SIGMA_X, SIGMA_Z, SMOOTH_SCALE, SPAN_MULT, NPASS)
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! ---------------------
      ! Arguments
      ! ---------------------
      INTEGER, INTENT(IN) :: NX, NZ, my_rank
      REAL(dp), INTENT(IN) :: DX_BLOCK(:), ZBC(:)
      REAL(dp), INTENT(IN) :: FREQ, CREF
      REAL(dp), INTENT(IN) :: GMASK(:)
      REAL(dp), INTENT(INOUT) :: G(:)

      ! Optional controls (backwards-compatible)
      REAL(dp), INTENT(IN), OPTIONAL :: SIGMA_X, SIGMA_Z
      REAL(dp), INTENT(IN), OPTIONAL :: SMOOTH_SCALE ! fraction of wavelength -> sigma_m = SMOOTH_SCALE * lambda
      REAL(dp), INTENT(IN), OPTIONAL :: SPAN_MULT    ! multiplier for truncation (default 3.0)
      INTEGER, INTENT(IN), OPTIONAL    :: NPASS      ! number of full separable passes (default 1)

      ! ---------------------
      ! Locals
      ! ---------------------
      INTEGER :: NX1, NZ1, NBLOCK
      REAL(dp) :: DXeff_block, DZeff_block
      REAL(dp) :: LAMBDA_REF, sigmaX, sigmaZ
      REAL(dp) :: Rx, Rz
      REAL(dp) :: smooth_scale_local
      REAL(dp) :: span_mult_local
      INTEGER :: SPANX, SPANZ
      INTEGER :: SPASS, pass
      INTEGER :: L, K, LK, offset, idx, row0
      REAL(dp) :: wsum, gsum, wloc
      REAL(dp) :: dsx, dsz   ! for possible diagnostics (kept)
      INTEGER :: offset0

      ! --- cached work buffers / metadata (preserve across calls) ---
      LOGICAL, SAVE :: kernels_valid = .FALSE.
      REAL(dp), ALLOCATABLE, SAVE :: KX_ws(:), KZ_ws(:), Gtmp_ws(:)
      INTEGER, SAVE :: cap_nblock = 0, cap_nx1 = 0, cap_nz1 = 0, spanx_ws = 0, spanz_ws = 0
      REAL(dp), SAVE :: sigx_ws = -1.0D0, sigz_ws = -1.0D0, spanmult_ws = -1.0D0
      INTEGER, SAVE :: npass_ws = 0

      ! ---------------------
      ! Constants and tiny guards
      ! ---------------------
      REAL(dp), PARAMETER :: tiny_dp = 1.0e-50_dp
      REAL(dp), PARAMETER :: default_smooth_scale = 0.15_dp  ! safer default than prior 0.5
      REAL(dp), PARAMETER :: default_span_mult = 2.0_dp
      INTEGER, PARAMETER :: default_npass = 1
      REAL(dp) :: sigma_z_pts_default

      ! ---------------------
      ! Begin
      ! ---------------------
      NX1 = NX - 1
      NZ1 = NZ - 1
      NBLOCK = NX1*NZ1
      if (PRESENT(SIGMA_Z) .AND. SIGMA_Z > 0.0_dp) then
         sigma_z_pts_default = 1.25_dp
      else
         sigma_z_pts_default = 1.25_dp
      end if
      ! Effective block spacings (representative)
      IF (SIZE(DX_BLOCK) >= 1) THEN
         DXeff_block = DX_BLOCK(1)
      ELSE
         DXeff_block = 1.0_dp
      END IF

      IF (SIZE(ZBC) >= 2) THEN
         DZeff_block = ABS(ZBC(2) - ZBC(1))
      ELSE IF (SIZE(ZBC) == 1) THEN
         DZeff_block = ZBC(1)
      ELSE
         DZeff_block = 1.0_dp
      END IF

      ! Determine span multiplier and npass
      span_mult_local = default_span_mult
      IF (PRESENT(SPAN_MULT)) THEN
         IF (SPAN_MULT > 0.0_dp) span_mult_local = SPAN_MULT
      END IF

      SPASS = default_npass
      IF (PRESENT(NPASS)) THEN
         IF (NPASS >= 1) SPASS = NPASS
      END IF

      ! Reference wavelength (guarded)
      IF (FREQ > 0.0_dp .AND. CREF > 0.0_dp) THEN
         LAMBDA_REF = CREF/FREQ
      ELSE
         LAMBDA_REF = 0.0_dp
      END IF

      ! Compute sigma in grid points (sigmaX, sigmaZ)
      IF (PRESENT(SIGMA_X) .AND. SIGMA_X > 0.0_dp) THEN
         sigmaX = MAX(SIGMA_X, tiny_dp)
      ELSE IF (PRESENT(SIGMA_Z) .AND. SIGMA_Z > 0.0_dp) THEN
         sigmaX = MAX(1.0_dp, tiny_dp)
      ELSE
         ! Use smooth_scale (fraction of wavelength) -> convert to grid points
         smooth_scale_local = default_smooth_scale
         IF (PRESENT(SMOOTH_SCALE)) THEN
            IF (SMOOTH_SCALE > 0.0_dp) smooth_scale_local = SMOOTH_SCALE
         END IF

         IF (LAMBDA_REF > 0.0_dp .AND. DXeff_block > 0.0_dp) THEN
            Rx = smooth_scale_local*LAMBDA_REF
            sigmaX = Rx/DXeff_block
         ELSE
            sigmaX = 0.0_dp
         END IF
      END IF

      IF (PRESENT(SIGMA_Z) .AND. SIGMA_Z > 0.0_dp) THEN
         sigmaZ = MAX(SIGMA_Z, tiny_dp)
      ELSE IF (PRESENT(SIGMA_X) .AND. SIGMA_X > 0.0_dp) THEN
         sigmaZ = sigma_z_pts_default
      ELSE
         IF (LAMBDA_REF > 0.0_dp .AND. DZeff_block > 0.0_dp) THEN
            Rz = smooth_scale_local*LAMBDA_REF
            sigmaZ = Rz/DZeff_block
         ELSE
            sigmaZ = 0.0_dp
         END IF
      END IF

      ! Guard small values
      IF (sigmaX <= 0.0_dp) sigmaX = tiny_dp
      IF (sigmaZ <= 0.0_dp) sigmaZ = tiny_dp

      ! Convert to integer kernel half-widths (truncation radius)
      SPANX = MAX(0, NINT(span_mult_local*sigmaX))
      SPANZ = MAX(0, NINT(span_mult_local*sigmaZ))

      ! Clamp to domain size
      IF (SPANX > NX1 - 1) SPANX = NX1 - 1
      IF (SPANZ > NZ1 - 1) SPANZ = NZ1 - 1

      ! ----------------------------------------------------------------------
      ! (Re)allocate / rebuild kernels if needed
      ! ----------------------------------------------------------------------
      IF (.NOT. kernels_valid .OR. &
          NX1 /= cap_nx1 .OR. NZ1 /= cap_nz1 .OR. NBLOCK /= cap_nblock .OR. &
          SPANX /= spanx_ws .OR. SPANZ /= spanz_ws .OR. &
          ABS(sigmaX - sigx_ws) > 1.0D-6 .OR. ABS(sigmaZ - sigz_ws) > 1.0D-6 .OR. &
          ABS(span_mult_local - spanmult_ws) > 1.0D-9 .OR. SPASS /= npass_ws) THEN

         ! Deallocate old if allocated
         IF (ALLOCATED(KX_ws)) DEALLOCATE (KX_ws)
         IF (ALLOCATED(KZ_ws)) DEALLOCATE (KZ_ws)
         IF (ALLOCATED(Gtmp_ws)) DEALLOCATE (Gtmp_ws)

         ! Allocate new workspace arrays
         ALLOCATE (KX_ws(-SPANX:SPANX))
         ALLOCATE (KZ_ws(-SPANZ:SPANZ))
         ALLOCATE (Gtmp_ws(NBLOCK))

         ! Build kernels (K*_ws normalized)
         CALL build_gaussian_kernel_sigma(SPANX, MAX(sigmaX, tiny_dp), KX_ws)
         CALL build_gaussian_kernel_sigma(SPANZ, MAX(sigmaZ, tiny_dp), KZ_ws)

         ! Update cache metadata
         cap_nx1 = NX1
         cap_nz1 = NZ1
         cap_nblock = NBLOCK
         spanx_ws = SPANX
         spanz_ws = SPANZ
         sigx_ws = sigmaX
         sigz_ws = sigmaZ
         spanmult_ws = span_mult_local
         npass_ws = SPASS
         kernels_valid = .TRUE.
      END IF

      ! ----------------------------------------------------------------------
      ! Separable Gaussian smoothing with mask weighting:
      !   For NPASS > 1, repeat the full separable pass NPASS times (variance adds).
      !   Each pass performs:
      !     1) horizontal (X) pass: G -> Gtmp_ws
      !     2) vertical (Z) pass:   Gtmp_ws -> G
      ! ----------------------------------------------------------------------
      DO pass = 1, SPASS

         ! ---- horizontal (X) pass ----
         DO K = 1, NZ1
            DO L = 1, NX1
               wsum = 0.0_dp
               gsum = 0.0_dp
               DO idx = -SPANX, SPANX
                  row0 = reflect_index(L + idx, NX1)
                  offset = (row0 - 1)*NZ1 + K
                  wloc = KX_ws(idx)
                  IF (GMASK(offset) > 0.0_dp) THEN
                     wsum = wsum + wloc*GMASK(offset)
                     gsum = gsum + wloc*GMASK(offset)*G(offset)
                  END IF
               END DO
               LK = (L - 1)*NZ1 + K
               IF (wsum > 0.0_dp) THEN
                  Gtmp_ws(LK) = gsum/wsum
               ELSE
                  Gtmp_ws(LK) = 0.0_dp
               END IF
            END DO
         END DO

         ! ---- vertical (Z) pass ----
         DO L = 1, NX1
            DO K = 1, NZ1
               wsum = 0.0_dp
               gsum = 0.0_dp
               DO idx = -SPANZ, SPANZ
                  row0 = reflect_index(K + idx, NZ1)
                  offset = (L - 1)*NZ1 + row0
                  wloc = KZ_ws(idx)
                  IF (GMASK(offset) > 0.0_dp) THEN
                     wsum = wsum + wloc*GMASK(offset)
                     gsum = gsum + wloc*GMASK(offset)*Gtmp_ws(offset)
                  END IF
               END DO
               LK = (L - 1)*NZ1 + K
               IF (wsum > 0.0_dp) THEN
                  G(LK) = gsum/wsum
               ELSE
                  G(LK) = 0.0_dp
               End IF
            END DO
         END DO

      END DO ! pass loop

   CONTAINS

      SUBROUTINE build_gaussian_kernel_sigma(SPAN, SIGMA, KERNEL)
         INTEGER, INTENT(IN) :: SPAN
         REAL(dp), INTENT(IN) :: SIGMA
         REAL(dp), INTENT(OUT) :: KERNEL(-SPAN:SPAN)
         INTEGER :: m
         REAL(dp) :: wsum
         IF (SPAN <= 0) THEN
            KERNEL(0) = 1.0_dp
            RETURN
         END IF
         wsum = 0.0_dp
         DO m = -SPAN, SPAN
            KERNEL(m) = EXP(-0.5_dp*(REAL(m, dp)/SIGMA)**2)
            wsum = wsum + KERNEL(m)
         END DO
         DO m = -SPAN, SPAN
            KERNEL(m) = KERNEL(m)/wsum
         END DO
      END SUBROUTINE build_gaussian_kernel_sigma

      PURE INTEGER FUNCTION reflect_index(k, n) RESULT(r)
         INTEGER, INTENT(IN) :: k, n
         INTEGER :: period, m
         IF (n <= 1) THEN
            r = 1
            RETURN
         END IF
         period = 2*n - 2
         m = MOD(k - 1, period)
         IF (m < 0) m = m + period
         IF (m < n) THEN
            r = m + 1
         ELSE
            r = 2*n - m - 1
         END IF
      END FUNCTION reflect_index

      PURE INTEGER FUNCTION reflect_inline(k, n) RESULT(r)
         INTEGER, INTENT(IN) :: k, n
         r = reflect_index(k, n)
      END FUNCTION reflect_inline

   END SUBROUTINE gradient_smooth_gauss_block

end module gradient_mod
