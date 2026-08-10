module lbfgs_mod
   use gradient_mod
   use Frechet_mod
   use shared_mod
   use grid_mod
   use gridtype_mod
   use iso_fortran_env, only: dp => real64, real64
   use constant_mod, only: tiny_dp
   implicit None

   INTEGER, SAVE :: lbfgs_hist_unit = -1
   LOGICAL, SAVE :: lbfgs_hist_ready = .FALSE.
   LOGICAL, SAVE :: lbfgs_hist_pending = .FALSE.
   INTEGER, SAVE :: lbfgs_hist_iter = 0
   INTEGER, SAVE :: lbfgs_hist_m = 0
   INTEGER, SAVE :: lbfgs_hist_curv_ok = 0
   REAL(real64), SAVE :: lbfgs_hist_gnorm = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_pnorm = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_gdotp = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_snorm = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_ynorm = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_sty = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_desc_cos = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_cos_sy = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_gamma0 = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_sty_min = 0.0_real64
   INTEGER, SAVE :: lbfgs_hist_used_pairs = 0
   INTEGER, SAVE :: lbfgs_hist_skipped_pairs = 0
   INTEGER, SAVE :: lbfgs_hist_bad_pair_col = 0
   INTEGER, SAVE :: lbfgs_hist_bad_reason = 0
   REAL(real64), SAVE :: lbfgs_hist_bad_bp = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_bad_snorm = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_bad_ynorm = 0.0_real64
   REAL(real64), SAVE :: lbfgs_hist_bad_cos = 0.0_real64

contains
   ! Module for L-BFGS optimization
   ! This module contains the lbfgs_step subroutine and other related functions:
   !1. lbfgs_step,  performs the L-BFGS update step.
   !2. ComputeHess, inverse  Hessian approximation diagonal.
   !3. d

   SUBROUTINE lbfgs_basic(iter, NM, mml, INVP, NPAR, NBLOCK, &
                          GRAD_scaled, GRAD_scaled_norm, HESSDI_scaled, &
                          p_k, grad_prev, BF_grad_res, BF_s_hist, &
                          m_coarse, m_coarse_prev, FCOST0, FCOST_prev, USE_LBFGS_TYPE, USE_PRECOND, my_rank)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! --------- dummy args ---------
      INTEGER, INTENT(IN)              :: iter, NM, mml, my_rank, NPAR, NBLOCK
      INTEGER, INTENT(IN)              :: INVP(:)
      REAL(dp), INTENT(IN)             :: GRAD_scaled(:)          ! gradient in scaled (LBFGS) space (RAW, not precond)
      REAL(dp), INTENT(IN)             :: GRAD_scaled_norm(:)     ! block norms (kept for possible scaling)
      REAL(dp), INTENT(IN)             :: HESSDI_scaled(:)         ! diag approx of inverse Hessian in scaled space
      REAL(dp), INTENT(INOUT)          :: p_k(:)                  ! search direction (scaled space)
      REAL(dp), INTENT(IN)             :: grad_prev(:)             ! gradient at previous iterate
      REAL(dp), INTENT(INOUT)          :: BF_grad_res(:, :)       ! L-BFGS y history (Δg)
      REAL(dp), INTENT(INOUT)          :: BF_s_hist(:, :)          ! L-BFGS s history (Δm)
      REAL(dp), INTENT(INOUT)          :: m_coarse(:), m_coarse_prev(:)
      REAL(dp), INTENT(IN)             :: FCOST0
      REAL(dp), INTENT(IN)             :: FCOST_prev
      INTEGER, INTENT(IN)              :: USE_LBFGS_TYPE
      LOGICAL, INTENT(IN)              :: USE_PRECOND             ! kept for interface compatibility

      ! --------- locals ----------
      INTEGER  :: iiter, m_LBFGS, im, jcol, IA, i1, i2, col
      INTEGER  :: nActive, IMAP
      REAL(dp) :: BP, BA, BB, BF_b
      REAL(dp) :: dot_g_p, gnorm, pnorm
      REAL(dp) :: sTy, sTs, yTy, sTy_min
      LOGICAL  :: restart_now

      REAL(dp) :: BF_a(mml), BF_p(mml)
      REAL(dp) :: BF_q(NM), BF_r(NM)
      LOGICAL  :: use_pair(mml)
      REAL(dp) :: BP_min, pscaler

      REAL(dp), PARAMETER :: ZERO = 0.0_dp
      REAL(dp), PARAMETER :: C_CURV = 1.0e-3_dp !0.1
      REAL(dp), PARAMETER :: TINYD = TINY(1.0_dp)

      INTERFACE
         FUNCTION DNRM2(N, X, INCX) RESULT(res)
            USE iso_fortran_env, ONLY: dp => real64
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N, INCX
            REAL(dp), INTENT(IN) :: X(*)
            REAL(dp) :: res
         END FUNCTION DNRM2
      END INTERFACE

      ! -------- init local arrays ----------
      BF_a = 0.0_dp
      BF_p = 0.0_dp
      BF_q = 0.0_dp
      BF_r = 0.0_dp
      use_pair = .FALSE.

      ! absolute floor for s^T y usage in two-loop (scaled space)
      BP_min = 1.0e-12_dp

      ! --------- sanity checks ----------
      IF (SIZE(GRAD_scaled) /= NM) THEN
         WRITE (*, *) 'LBFGS: SIZE(GRAD_scaled)=', SIZE(GRAD_scaled), ' NM=', NM
         STOP 'LBFGS: GRAD_scaled size mismatch'
      END IF

      IF (SIZE(HESSDI_scaled) /= NM) THEN
         WRITE (*, *) 'LBFGS: SIZE(HESSDI_scaled)=', SIZE(HESSDI_scaled), ' NM=', NM
         STOP 'LBFGS: HESSDI_scaled size mismatch'
      END IF

      IF (SIZE(p_k) < NM) THEN
         WRITE (*, *) 'LBFGS: SIZE(p_k)=', SIZE(p_k), ' NM=', NM
         STOP 'LBFGS: p_k size mismatch'
      END IF

      IF (SIZE(m_coarse) < NM .OR. SIZE(m_coarse_prev) < NM) THEN
         WRITE (*, *) 'LBFGS: VV sizes: final=', SIZE(m_coarse), ' prev=', SIZE(m_coarse_prev), ' NM=', NM
         STOP 'LBFGS: BF_VV_* size mismatch'
      END IF

      IF (SIZE(BF_grad_res, 1) < NM .OR. SIZE(BF_s_hist, 1) < NM) THEN
         WRITE (*, *) 'LBFGS: BF_*_res dim1: grad=', SIZE(BF_grad_res, 1), ' V=', SIZE(BF_s_hist, 1), ' NM=', NM
         STOP 'LBFGS: BF_*_res dim1 mismatch'
      END IF

      IF (SIZE(BF_grad_res, 2) < mml .OR. SIZE(BF_s_hist, 2) < mml) THEN
         WRITE (*, *) 'LBFGS: BF_*_res dim2: grad=', SIZE(BF_grad_res, 2), ' V=', SIZE(BF_s_hist, 2), ' mml=', mml
         STOP 'LBFGS: BF_*_res dim2 mismatch'
      END IF

      ! quick exit if no active parameters
      nActive = COUNT(INVP == 1)
      IF (nActive == 0 .OR. NM <= 0) THEN
         p_k(1:NM) = 0.0_dp
         RETURN
      END IF

      ! -------- diagnostics ----------
      gnorm = DNRM2(NM, GRAD_scaled, 1)
      p_k(1:NM) = 0.0_dp

      CALL lbfgs_hist_init(my_rank)
      CALL lbfgs_hist_reset_pairs()

      iiter = iter - 1

      ! ==========================
      ! Iter 0: preconditioned steepest direction
      ! ==========================
      IF (iiter == 0) THEN

         ! Preconditioned steepest descent: p = -H0 * g  with H0 = HESSDI_scaled
         p_k(1:NM) = 0.0_dp

         call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                            GRAD_scaled, HESSDI_scaled, &
                            USE_PRECOND, my_rank, p_k)

         dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)

         CALL lbfgs_hist_stage(iter, 0, gnorm, pnorm, dot_g_p, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1)

         RETURN
      END IF

      ! ==========================
      ! iiter > 0: L-BFGS two-loop
      ! ==========================
      m_LBFGS = MIN(iiter, mml)

      ! ----------- Warm-up: 1 <= iiter <= mml -----------
      IF (iiter > 0 .AND. iiter <= mml) THEN

         col = mml + 1 - iiter

         ! y_k = g_k - g_{k-1} in scaled space
         BF_grad_res(1:NM, col) = GRAD_scaled(1:NM) - grad_prev(1:NM)
         ! s_k = m_k - m_{k-1} in scaled model space
         BF_s_hist(1:NM, col) = m_coarse(1:NM) - m_coarse_prev(1:NM)

         IF (USE_LBFGS_TYPE == 1) THEN
            CALL lbfgs_cost_damp_pair(NM, FCOST0, FCOST_prev, grad_prev(1:NM), GRAD_scaled(1:NM), &
                                      BF_s_hist(1:NM, col), BF_grad_res(1:NM, col))
         END IF

         sTy = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_grad_res(1:NM, col))
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_s_hist(1:NM, col))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
         sTy_min = C_CURV*SQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
         restart_now = ((sTy <= sTy_min) .OR. .NOT. (sTy == sTy) .OR. (sTy == 0.0_dp))

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (warm): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

1000     CONTINUE
         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0_dp
            BF_s_hist(1:NM, 1:mml) = 0.0_dp

            ! restart with preconditioned steepest descent

            call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                               GRAD_scaled, HESSDI_scaled, &
                               USE_PRECOND, my_rank, p_k)
            dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
            pnorm = DNRM2(NM, p_k, 1)

            IF (dot_g_p >= ZERO) THEN
               WRITE (*, '(A,1X,G0)') 'LBFGS ERROR: warm restart produced non-descent direction, g·p =', dot_g_p
               STOP
            END IF

            IF (my_rank == 0) THEN
               WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm
            END IF

            CALL lbfgs_hist_stage(iter, 0, gnorm, pnorm, dot_g_p, SQRT(MAX(sTs, 0.0_dp)), &
                                  SQRT(MAX(yTy, 0.0_dp)), sTy, sTy_min, 0)

            RETURN
         END IF

         BF_a(1:m_LBFGS) = 0.0_dp
         BF_p(1:m_LBFGS) = 0.0_dp
         BF_q(1:NM) = GRAD_scaled(1:NM)
         use_pair(1:m_LBFGS) = .FALSE.

         ! ---- first loop (backward application of stored updates) ----
         DO im = 1, m_LBFGS
            jcol = mml - iiter + im

            BP = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_s_hist(1:NM, jcol))   ! s^T y
            IF (.NOT. (BP == BP)) THEN
               CALL lbfgs_hist_note_pair(jcol, BP, &
                                         DNRM2(NM, BF_s_hist(1:NM, jcol), 1), &
                                         DNRM2(NM, BF_grad_res(1:NM, jcol), 1), &
                                         .FALSE., 1)
               CYCLE
            END IF
            IF (ABS(BP) <= BP_min) THEN
               CALL lbfgs_hist_note_pair(jcol, BP, &
                                         DNRM2(NM, BF_s_hist(1:NM, jcol), 1), &
                                         DNRM2(NM, BF_grad_res(1:NM, jcol), 1), &
                                         .FALSE., 2)
               CYCLE
            END IF

            BA = DOT_PRODUCT(BF_s_hist(1:NM, jcol), BF_q(1:NM))                ! s^T q

            BF_p(im) = 1.0_dp/BP
            BF_a(im) = BF_p(im)*BA
            BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, jcol)

            use_pair(im) = .TRUE.
            CALL lbfgs_hist_note_pair(jcol, BP, &
                                      DNRM2(NM, BF_s_hist(1:NM, jcol), 1), &
                                      DNRM2(NM, BF_grad_res(1:NM, jcol), 1), &
                                      .TRUE., 0)
         END DO

         CALL lbfgs_hist_report_pairs(my_rank)

         ! ---- H0 application: diagonal inverse Hessian in scaled space ----
         ! DO i1 = 1, NM
         !    IF (HESSDI_scaled(i1) > 0.0_dp) THEN
         !       BF_r(i1) = HESSDI_scaled(i1)*BF_q(i1)
         !    ELSE
         !       BF_r(i1) = BF_q(i1)               ! fallback: identity if H0 entry invalid
         !    END IF
         ! END DO
  BF_r(1:NM) = BF_q(1:NM)
         ! ---- second loop ----
         DO im = m_LBFGS, 1, -1
            IF (.NOT. use_pair(im)) CYCLE

            jcol = mml - iiter + im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_r(1:NM))
            BF_b = BF_p(im)*BB
            BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, jcol)*(BF_a(im) - BF_b)
         END DO

         p_k(1:NM) = -BF_r(1:NM)
         dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)

         IF (dot_g_p >= ZERO) THEN
            WRITE (*, *) 'LBFGS ERROR: warm two-loop produced non-descent direction, g·p =', dot_g_p
            STOP
         END IF

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step: g·p =', dot_g_p, '||p||=', pnorm
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step: min(p)=', MINVAL(p_k), &
               'max(p)=', MAXVAL(p_k)
         END IF

         CALL lbfgs_hist_stage(iter, m_LBFGS, gnorm, pnorm, dot_g_p, SQRT(MAX(sTs, 0.0_dp)), &
                               SQRT(MAX(yTy, 0.0_dp)), sTy, sTy_min, 1)

         RETURN
      END IF

      ! ----------- Full memory: iiter > mml -----------
      IF (iiter > mml) THEN

         ! shift history down
         DO im = mml, 2, -1
            BF_grad_res(1:NM, im) = BF_grad_res(1:NM, im - 1)
            BF_s_hist(1:NM, im) = BF_s_hist(1:NM, im - 1)
         END DO
         BF_grad_res(1:NM, 1) = GRAD_scaled(1:NM) - grad_prev(1:NM)
         BF_s_hist(1:NM, 1) = m_coarse(1:NM) - m_coarse_prev(1:NM)

         ! IF (USE_LBFGS_TYPE == 1) THEN
         !    CALL lbfgs_cost_damp_pair(NM, FCOST0, FCOST_prev, grad_prev(1:NM), GRAD_scaled(1:NM), &
         !                              BF_s_hist(1:NM, 1), BF_grad_res(1:NM, 1))
         ! END IF

         sTy = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_grad_res(1:NM, 1))
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_s_hist(1:NM, 1))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
         sTy_min = C_CURV*SQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
         restart_now = ((sTy <= sTy_min) .OR. .NOT. (sTy == sTy) .OR. (sTy == 0.0_dp))

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (full): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

2000     CONTINUE
         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0_dp
            BF_s_hist(1:NM, 1:mml) = 0.0_dp

            ! restart with preconditioned steepest descent

            call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                               GRAD_scaled, HESSDI_scaled, &
                               USE_PRECOND, my_rank, p_k)
            dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
            pnorm = DNRM2(NM, p_k, 1)

            IF (dot_g_p >= ZERO) THEN
               WRITE (*, *) 'LBFGS ERROR: full restart produced non-descent direction, g·p =', dot_g_p
               STOP
            END IF

            IF (my_rank == 0) THEN
               WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm
            END IF

            CALL lbfgs_hist_stage(iter, 0, gnorm, pnorm, dot_g_p, SQRT(MAX(sTs, 0.0_dp)), &
                                  SQRT(MAX(yTy, 0.0_dp)), sTy, sTy_min, 0)

            RETURN
         END IF

         m_LBFGS = mml
         BF_a(1:m_LBFGS) = 0.0_dp
         BF_p(1:m_LBFGS) = 0.0_dp
         BF_q(1:NM) = GRAD_scaled(1:NM)
         use_pair(1:m_LBFGS) = .FALSE.

         ! ---- first loop ----
         DO im = 1, m_LBFGS
            col = im

            BP = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))     ! s^T y
            IF (.NOT. (BP == BP)) THEN
               CALL lbfgs_hist_note_pair(col, BP, &
                                         DNRM2(NM, BF_s_hist(1:NM, col), 1), &
                                         DNRM2(NM, BF_grad_res(1:NM, col), 1), &
                                         .FALSE., 1)
               CYCLE
            END IF
            IF (ABS(BP) <= BP_min) THEN
               CALL lbfgs_hist_note_pair(col, BP, &
                                         DNRM2(NM, BF_s_hist(1:NM, col), 1), &
                                         DNRM2(NM, BF_grad_res(1:NM, col), 1), &
                                         .FALSE., 2)
               CYCLE
            END IF

            BA = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_q(1:NM))                ! s^T q

            BF_p(im) = 1.0_dp/BP
            BF_a(im) = BF_p(im)*BA
            BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, col)

            use_pair(im) = .TRUE.
            CALL lbfgs_hist_note_pair(col, BP, &
                                      DNRM2(NM, BF_s_hist(1:NM, col), 1), &
                                      DNRM2(NM, BF_grad_res(1:NM, col), 1), &
                                      .TRUE., 0)
         END DO

         CALL lbfgs_hist_report_pairs(my_rank)

         ! ---- H0 application: diag inverse Hessian ----
         ! DO i1 = 1, NM
         !    IF (HESSDI_scaled(i1) > 0.0_dp) THEN
         !       BF_r(i1) = HESSDI_scaled(i1)*BF_q(i1)
         !    ELSE
         !       BF_r(i1) = BF_q(i1)
         !    END IF
         ! END DO
 BF_r(1:NM) = BF_q(1:NM)
         ! ---- second loop ----
         DO im = m_LBFGS, 1, -1
            IF (.NOT. use_pair(im)) CYCLE

            col = im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_r(1:NM))
            BF_b = BF_p(im)*BB
            BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, col)*(BF_a(im) - BF_b)
         END DO

         p_k(1:NM) = -BF_r(1:NM)
         dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)

         IF (dot_g_p >= ZERO) THEN
            WRITE (*, *) 'LBFGS ERROR: full two-loop produced non-descent direction, g·p =', dot_g_p
            STOP
         END IF

         IF (my_rank == 0) THEN
            WRITE (*, *)
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step(full): g·p =', dot_g_p, '||p||=', pnorm
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step(full): min(p)=', MINVAL(p_k), &
               'max(p)=', MAXVAL(p_k)
         END IF

         CALL lbfgs_hist_stage(iter, m_LBFGS, gnorm, pnorm, dot_g_p, SQRT(MAX(sTs, 0.0_dp)), &
                               SQRT(MAX(yTy, 0.0_dp)), sTy, sTy_min, 1)

         RETURN
      END IF

   END SUBROUTINE lbfgs_basic

   SUBROUTINE lbfgs_cost_damp_pair(NM, FCOST0, FCOST_prev, grad_prev, grad_curr, svec, yvec)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NM
      REAL(dp), INTENT(IN) :: FCOST0, FCOST_prev
      REAL(dp), INTENT(IN) :: grad_prev(:), grad_curr(:), svec(:)
      REAL(dp), INTENT(INOUT) :: yvec(:)

      REAL(dp) :: s2, theta, ys_hat
      REAL(dp), PARAMETER :: EPS_CURV = 1.0e-5_dp
      REAL(dp), PARAMETER :: EPS_THETA = 1.0e-30_dp
      REAL(dp), PARAMETER :: TINYD = 1.0e-37_dp

      s2 = DOT_PRODUCT(svec(1:NM), svec(1:NM))
      theta = 6.0_dp*(FCOST0 - FCOST_prev) + &
              3.0_dp*DOT_PRODUCT(grad_prev(1:NM) + grad_curr(1:NM), svec(1:NM))

      IF (s2 > TINYD) THEN
         IF (ABS(theta) > EPS_THETA) yvec(1:NM) = yvec(1:NM) + (theta/s2)*svec(1:NM)
      END IF

      ys_hat = DOT_PRODUCT(yvec(1:NM), svec(1:NM))
      IF (ys_hat <= EPS_CURV*MAX(s2, 1.0_dp)) THEN
         yvec(1:NM) = yvec(1:NM) + ((EPS_CURV*MAX(s2, 1.0_dp) - ys_hat)/MAX(s2, 1.0_dp))*svec(1:NM)
      END IF
   END SUBROUTINE lbfgs_cost_damp_pair

   SUBROUTINE steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                            GRAD_scaled, HESSDI_scaled, &
                            USE_PRECOND, &
                            my_rank, p_k)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      !------------- dummy args -------------
      INTEGER, INTENT(IN)  :: iter, NM, NPAR, NBLOCK
      INTEGER, INTENT(IN)  :: INVP(:)
      REAL(dp), INTENT(IN) :: GRAD_scaled(:)
      REAL(dp), INTENT(IN) :: HESSDI_scaled(:)
      LOGICAL, INTENT(IN)  :: USE_PRECOND
      INTEGER, INTENT(IN)  :: my_rank
      REAL(dp), INTENT(INOUT) :: p_k(:)

      !------------- locals -------------
      INTEGER  :: IA, IMAP, i1, i2, ik
      REAL(dp) :: pscaler
      REAL(dp) :: dot_g_p, pnorm
      REAL(dp), PARAMETER :: ZERO = 0.0_dp
      REAL(dp), PARAMETER :: TINYD = TINY(1.0_dp)

      INTERFACE
         FUNCTION DNRM2(N, X, INCX) RESULT(res)
            USE iso_fortran_env, ONLY: dp => real64
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N, INCX
            REAL(dp), INTENT(IN) :: X(*)
            REAL(dp) :: res
         END FUNCTION DNRM2
      END INTERFACE

      !===========================================================
      ! Build steepest-descent direction
      ! If no precond: HESSDI_scaled is conceptually identity (1)
      !===========================================================

      p_k(1:NM) = 0.0_dp

      IMAP = 0
      DO IA = 1, NPAR
         IF (INVP(IA) == 1) THEN
            IMAP = IMAP + 1
            i1 = (IMAP - 1)*NBLOCK + 1
            i2 = IMAP*NBLOCK
             p_k(i1:i2) = -GRAD_scaled(i1:i2)! *HESSDI_scaled(i1:i2) 
             if (use_precond == .FALSE.) then
               pscaler = 1.0_dp/DNRM2(NBLOCK, p_k(i1:i2), 1)
               ! pscaler = 1.0_dp
               p_k(i1:i2) = p_k(i1:i2)*pscaler
             end if

         END IF
      END DO

      !===========================
      ! Descent check and logging
      !===========================
      dot_g_p = DOT_PRODUCT(GRAD_scaled(1:NM), p_k(1:NM))
      pnorm = DNRM2(NM, p_k, 1)

      IF (dot_g_p >= ZERO) THEN
         WRITE (*, *) 'LBFGS ERROR:  non-descent direction, g·p =', dot_g_p
         STOP
      END IF

      IF (my_rank == 0) THEN
         WRITE (*, '(A,1X,I6)') 'LBFGS steepest step , iter=', iter
         WRITE (*, '(A,ES12.4,2X,A,ES12.4)') '  g·p =', dot_g_p, '||p|| =', pnorm
         WRITE (*, '(A,ES12.4,2X,A,ES12.4)') '  min(p)=', MINVAL(p_k(1:NM)), 'max(p)=', MAXVAL(p_k(1:NM))
      END IF

   END SUBROUTINE steep_descent

   SUBROUTINE lbfgs_step_loop(iter, NM, mml, INVP, NPAR, NBLOCK, &
                              GRAD_Scaled, GRAD_scaled_norm, HESSDI_scaled, &
                              p_k, grad_prev, BF_grad_res, &
                              BF_s_hist, m_coarse, m_coarse_prev, FCOST0, FCOST_prev, &
                              USE_LBFGS_TYPE, USE_PRECOND, my_rank)

      USE, INTRINSIC :: ieee_arithmetic            ! [E] finite checks
      IMPLICIT NONE

      ! --------- dummy args ---------
      INTEGER, INTENT(IN)    :: iter, NM, mml, my_rank, INVP(:), NPAR, NBLOCK
      REAL(dp), INTENT(IN)   :: GRAD_Scaled(:)                 ! scaled gradient g_k
      REAL(dp), INTENT(IN)   :: GRAD_Scaled_Norm(:)            ! per-parameter norms
      REAL(dp), INTENT(IN)   :: HESSDI_scaled(:)                ! diagonal preconditioner in scaled space
      REAL(dp), INTENT(INOUT):: p_k(:)
      REAL(dp), INTENT(IN)   :: m_coarse(:), m_coarse_prev(:)
      REAL(dp), INTENT(IN)   :: grad_prev(:)                    ! previous scaled gradient g_{k-1}
      REAL(dp), INTENT(INOUT):: BF_grad_res(:, :)              ! Y history (y = g_k - g_{k-1})
      REAL(dp), INTENT(INOUT):: BF_s_hist(:, :)                 ! S history (s = x_k - x_{k-1})
      REAL(dp), INTENT(IN)   :: FCOST0
      REAL(dp), INTENT(IN)   :: FCOST_prev
      INTEGER, INTENT(IN)    :: USE_LBFGS_TYPE
      LOGICAL, INTENT(IN)    :: USE_PRECOND

      ! --------- locals ----------
      REAL(dp) :: BP, BA, BPK, BB, BF_b
      INTEGER  :: iiter, ik, m_LBFGS, im, jcol, IA, i1, i2, col, IMAP
      REAL(dp) :: dot_g_p, gnorm, pnorm, gHg, minH, maxH, posMaxH
      REAL(dp) :: BF_a(mml), BF_p(mml)
      REAL(dp) :: BF_q(NM), BF_r(NM), pscaler, precond_grad_norm, precond_grad(NM)
      REAL(dp) :: sTy, sTs, yTy, sTy_min
      REAL(dp) :: denom, tau, yy, gamma
      REAL(dp) :: pnorm_blk, gblk_norm
      LOGICAL  :: restart_now
      REAL(dp), PARAMETER :: ZERO = 0.0_dp
      REAL(dp), PARAMETER :: C_CURV = 1.0e-10_dp
      REAL(dp), PARAMETER :: TINYD = 1.0e-37_dp

      ! --------- cost-damping / trust-region ----------
      REAL(dp) :: s2, theta, ys_raw, ys_hat, yy_damp
      REAL(dp) :: ytmp(NM), svec(NM), yhat(NM)
      REAL(dp), PARAMETER :: EPS_CURV = 1.0e-5_dp
      REAL(dp), PARAMETER :: DELTA_TR = 1.0_dp
      REAL(dp) :: pnorm_tr, scale_tr

      REAL(dp), PARAMETER :: EPS_THETA = 1.0e-30_dp

      restart_now = .FALSE.
      BF_a(1:mml) = 0.0_dp
      BF_p(1:mml) = 0.0_dp
      BF_q(1:NM) = 0.0_dp
      BF_r(1:NM) = 0.0_dp
      precond_grad(1:NM) = 0.0_dp
      ytmp(1:NM) = 0.0_dp
      svec(1:NM) = 0.0_dp
      yhat(1:NM) = 0.0_dp
      ! ========= diagnostics =========
      gnorm = DNRM2(NM, GRAD_Scaled, 1)
      ! IF (USE_PRECOND) THEN

      !    gHg = 0.0_dp
      !    DO ik = 1, NM
      !       IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) CYCLE
      !       gHg = gHg + GRAD_Scaled(ik)*(HESSDI_scaled(ik)*GRAD_Scaled(ik))
      !    END DO
      !    minH = MINVAL(HESSDI_scaled)
      !    maxH = MAXVAL(HESSDI_scaled)
      !    posMaxH = MAXVAL(HESSDI_scaled, MASK=(HESSDI_scaled > 0.0_dp))
      ! ELSE
      !    gHg = DOT_PRODUCT(GRAD_Scaled, GRAD_Scaled)
      !    minH = 1.0_dp
      !    maxH = 1.0_dp
      !    posMaxH = 1.0_dp
      ! END IF

      ! IF (my_rank == 0) THEN
      !    WRITE (*, *) '  '
      !    WRITE (*, '(A,I6,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
      !       'LBFGS: iter=', iter, '||g||=', gnorm, 'min(H)=', minH, 'max(H)=', maxH, 'max+(H)=', posMaxH
      !    WRITE (*, '(A,ES12.4)') 'LBFGS: g^T (H ⊙ g) =', gHg
      ! END IF

      iiter = iter - 1

      ! ==========================
      ! Iter 0: steepest direction
      ! ==========================
      IF (iiter == 0) THEN

         call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                            GRAD_scaled, HESSDI_scaled, &
                            USE_PRECOND, my_rank, p_k)

         dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)
         IF (my_rank == 0) WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'iter0: min(p)=', MINVAL(p_k), 'max(p)=', MAXVAL(p_k)
         IF (.NOT. ieee_is_finite(dot_g_p) .OR. dot_g_p >= ZERO) THEN
            IF (my_rank == 0) WRITE (*, *) 'WARNING: Non-descent/NaN at iter=0 → stop.'
            STOP
         END IF

         IF (USE_LBFGS_TYPE >= 2) THEN
            pnorm_tr = DNRM2(NM, p_k, 1)
            IF (pnorm_tr > DELTA_TR) THEN
               scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
               p_k(1:NM) = p_k(1:NM)*scale_tr
               IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(iter0): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
            END IF
         END IF

         RETURN
      END IF

      ! ==========================
      ! iiter > 0: L-BFGS two-loop
      ! ==========================
      m_LBFGS = MIN(iiter, mml)

      ! -------- Warm-up window: 1 <= iiter <= mml --------
      IF (iiter > 0 .AND. iiter <= mml) THEN
         col = mml + 1 - iiter

         BF_grad_res(1:NM, col) = GRAD_Scaled(1:NM) - grad_prev(1:NM)     ! y_k
         BF_s_hist(1:NM, col) = m_coarse(1:NM) - m_coarse_prev(1:NM)   ! s_k

         ! ---- cost-damped y (type>=1) ----
         IF (USE_LBFGS_TYPE >= 1) THEN
            ytmp(1:NM) = BF_grad_res(1:NM, col)
            svec(1:NM) = BF_s_hist(1:NM, col)
            s2 = DOT_PRODUCT(svec, svec)
            theta = 6.0_dp*(FCOST0 - FCOST_prev) + &
                    3.0_dp*DOT_PRODUCT(grad_prev(1:NM) + GRAD_Scaled(1:NM), svec(1:NM))

            IF (s2 > TINYD) THEN
               yhat(1:NM) = ytmp(1:NM)
               IF (DABS(theta) > EPS_THETA) yhat(1:NM) = yhat(1:NM) + (theta/s2)*svec(1:NM)
            ELSE
               yhat(1:NM) = ytmp(1:NM)
            END IF

            ys_raw = DOT_PRODUCT(ytmp, svec)
            ys_hat = DOT_PRODUCT(yhat, svec)
            yy_damp = DOT_PRODUCT(yhat, yhat)

            IF (ys_hat <= EPS_CURV*MAX(s2, 1.0_dp)) THEN
               yhat(1:NM) = yhat(1:NM) + ((EPS_CURV*MAX(s2, 1.0_dp) - ys_hat)/MAX(s2, 1.0_dp))*svec(1:NM)
               ys_hat = DOT_PRODUCT(yhat, svec)
            END IF

            BF_grad_res(1:NM, col) = yhat(1:NM)

            IF (my_rank == 0) THEN
               WRITE (*, '(A)') 'LBFGS[damp]: cost-based y update'
               WRITE (*, '(2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
                  'ΔJ=', (FCOST0 - FCOST_prev), 'theta=', theta, 'y·s (raw)=', ys_raw, 'y·s (damped)=', ys_hat
            END IF
         END IF

         ! curvature gate
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_s_hist(1:NM, col))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
         sTy = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_grad_res(1:NM, col))
         IF (.NOT. ieee_is_finite(sTs) .OR. .NOT. ieee_is_finite(yTy) .OR. &
             .NOT. ieee_is_finite(sTy)) THEN
            restart_now = .TRUE.
         ELSE
            sTy_min = C_CURV*DSQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
            restart_now = (sTy <= sTy_min)
         END IF
         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (warm): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0_dp
            BF_s_hist(1:NM, 1:mml) = 0.0_dp

            ! fallback steepest

            call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                               GRAD_scaled, HESSDI_scaled, &
                               USE_PRECOND, my_rank, p_k)

            dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
            pnorm = DNRM2(NM, p_k, 1)
            IF (my_rank == 0) WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm

            IF (USE_LBFGS_TYPE >= 2) THEN
               pnorm_tr = DNRM2(NM, p_k, 1)
               IF (pnorm_tr > DELTA_TR) THEN
                  scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
                  p_k(1:NM) = p_k(1:NM)*scale_tr
                  IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(restart): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
               END IF
            END IF

            RETURN
         END IF

         ! two-loop prep
         BF_a(1:m_LBFGS) = 0.0_dp
         BF_p(1:m_LBFGS) = 0.0_dp
         BF_q(1:NM) = GRAD_Scaled(1:NM)

         ! backward loop
         DO im = 1, m_LBFGS
            jcol = mml - iiter + im
            IF (jcol < 1 .OR. jcol > SIZE(BF_grad_res, 2) .OR. jcol > SIZE(BF_s_hist, 2)) THEN
               WRITE (*, *) 'LBFGS_STEP_LOOP: bad jcol on rank', my_rank, &
                  ' jcol=', jcol, ' size2(BF_grad_res)=', SIZE(BF_grad_res, 2), &
                  ' size2(BF_s_hist)=', SIZE(BF_s_hist, 2)
               STOP 'LBFGS: jcol out of bounds'
            END IF
            BP = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_s_hist(1:NM, jcol))   ! y·s
            BA = DOT_PRODUCT(BF_s_hist(1:NM, jcol), BF_q(1:NM))                ! s·q

            IF (.NOT. ieee_is_finite(BP) .OR. DABS(BP) < TINYD) THEN
               BF_p(im) = 0.0_dp
               BF_a(im) = 0.0_dp
            ELSE
               BF_p(im) = 1.0_dp/BP
               IF (.NOT. ieee_is_finite(BA)) BA = 0.0_dp
               BF_a(im) = BF_p(im)*BA
            END IF

            IF (BF_a(im) /= 0.0_dp) THEN
               BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, jcol)
            END IF
         END DO

         BPK = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))     ! y·s

         ! IF (USE_PRECOND) THEN
         !    ! Diagonal preconditioning: H₀ = τ·diag(H)
         !    denom = 0.0_dp
         !    DO ik = 1, NM
         !       IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) CYCLE
         !       denom = denom + BF_grad_res(ik, col)*(HESSDI_scaled(ik)*BF_grad_res(ik, col))
         !    END DO
         !    denom = MAX(denom, TINYD)
         !    IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
         !    tau = BPK/denom
         !    IF (.NOT. ieee_is_finite(tau) .OR. tau <= 0.0_dp) tau = 1.0_dp

         !    DO ik = 1, NM
         !       IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
         !          BF_r(ik) = BF_q(ik)
         !       ELSE
         !          BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
         !       END IF
         !    END DO
         ! ELSE
         !    ! Scalar initial Hessian: H₀ = γ·I
            yy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
            yy = MAX(yy, TINYD)
            IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
            gamma = BPK/yy
            IF (.NOT. ieee_is_finite(gamma) .OR. gamma <= 0.0_dp) gamma = 1.0_dp
            BF_r(1:NM) = gamma*BF_q(1:NM)
         ! END IF

         ! forward loop
         DO im = m_LBFGS, 1, -1
            jcol = mml - iiter + im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_r(1:NM))   ! y·r
            IF (.NOT. ieee_is_finite(BF_p(im))) CYCLE
            IF (.NOT. ieee_is_finite(BB)) CYCLE
            BF_b = BF_p(im)*BB
            IF (.NOT. ieee_is_finite(BF_b)) CYCLE
            BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, jcol)*(BF_a(im) - BF_b)
         END DO

         p_k(1:NM) = -BF_r(1:NM)

         IF (USE_LBFGS_TYPE >= 2) THEN
            pnorm_tr = DNRM2(NM, p_k, 1)
            IF (pnorm_tr > DELTA_TR) THEN
               scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
               p_k(1:NM) = p_k(1:NM)*scale_tr
               IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(warm): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
            END IF
         END IF

         dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)
         IF (.NOT. ieee_is_finite(dot_g_p) .OR. .NOT. ieee_is_finite(pnorm)) THEN
            IF (my_rank == 0) WRITE (*, *) 'LBFGS produced non-finite step → fallback steepest.'

            call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                               GRAD_scaled, HESSDI_scaled, &
                               USE_PRECOND, my_rank, p_k)
         END IF

         RETURN
      END IF

      ! -------- Full memory: iiter > mml --------
      IF (iiter > mml) THEN
         DO im = mml, 2, -1
            BF_grad_res(1:NM, im) = BF_grad_res(1:NM, im - 1)
            BF_s_hist(1:NM, im) = BF_s_hist(1:NM, im - 1)
         END DO
         BF_grad_res(1:NM, 1) = GRAD_Scaled(1:NM) - grad_prev(1:NM)
         BF_s_hist(1:NM, 1) = m_coarse(1:NM) - m_coarse_prev(1:NM)

         IF (USE_LBFGS_TYPE >= 1) THEN
            ytmp(1:NM) = BF_grad_res(1:NM, 1)
            svec(1:NM) = BF_s_hist(1:NM, 1)
            s2 = DOT_PRODUCT(svec, svec)
            theta = 6.0_dp*(FCOST0 - FCOST_prev) + &
                    3.0_dp*DOT_PRODUCT(GRAD_Scaled(1:NM) + grad_prev(1:NM), svec(1:NM))
            IF (s2 > TINYD) THEN
               yhat(1:NM) = ytmp(1:NM) + (theta/s2)*svec(1:NM)
            ELSE
               yhat(1:NM) = ytmp(1:NM)
            END IF
            ys_raw = DOT_PRODUCT(ytmp, svec)
            ys_hat = DOT_PRODUCT(yhat, svec)
            yy_damp = DOT_PRODUCT(yhat, yhat)
            IF (ys_hat <= EPS_CURV*MAX(s2, 1.0_dp)) THEN
               yhat(1:NM) = yhat(1:NM) + ((EPS_CURV*MAX(s2, 1.0_dp) - ys_hat)/MAX(s2, 1.0_dp))*svec(1:NM)
               ys_hat = DOT_PRODUCT(yhat, svec)
            END IF
            BF_grad_res(1:NM, 1) = yhat(1:NM)
            IF (my_rank == 0) THEN
               WRITE (*, '(A)') 'LBFGS[damp]: cost-based y update'
               WRITE (*, '(2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
                  'ΔJ=', (FCOST0 - FCOST_prev), 'theta=', theta, 'y·s (raw)=', ys_raw, 'y·s (damped)=', ys_hat
            END IF
         END IF

         ! curvature gate on newest
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_s_hist(1:NM, 1))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
         sTy = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_grad_res(1:NM, 1))
         sTy_min = C_CURV*DSQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
         restart_now = (.NOT. ieee_is_finite(sTy)) .OR. (sTy <= sTy_min)

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (full): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0_dp
            BF_s_hist(1:NM, 1:mml) = 0.0_dp

            call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                               GRAD_scaled, HESSDI_scaled, &
                               USE_PRECOND, my_rank, p_k)

            dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
            pnorm = DNRM2(NM, p_k, 1)
            IF (my_rank == 0) WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm

            IF (USE_LBFGS_TYPE >= 2) THEN
               pnorm_tr = DNRM2(NM, p_k, 1)
               IF (pnorm_tr > DELTA_TR) THEN
                  scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
                  p_k(1:NM) = p_k(1:NM)*scale_tr
                  IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(restart): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
               END IF
            END IF

            RETURN
         END IF

         m_LBFGS = mml
         BF_a(1:m_LBFGS) = 0.0_dp
         BF_p(1:m_LBFGS) = 0.0_dp
         BF_q(1:NM) = GRAD_Scaled(1:NM)

         ! backward loop (full)
         DO im = 1, m_LBFGS
            col = im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))   ! y·s
            BP = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_q(1:NM))               ! s·q

            IF (.NOT. ieee_is_finite(BB) .OR. DABS(BB) < TINYD) THEN
               BF_p(im) = 0.0_dp
               BF_a(im) = 0.0_dp
            ELSE
               BF_p(im) = 1.0_dp/BB
               IF (.NOT. ieee_is_finite(BP)) BP = 0.0_dp
               BF_a(im) = BF_p(im)*BP
            END IF

            IF (BF_a(im) /= 0.0_dp) THEN
               BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, col)
            END IF
         END DO

         ! H0 from newest (col=1)
         BPK = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_s_hist(1:NM, 1))     ! y·s
         ! IF (USE_PRECOND) THEN
         !    denom = 0.0_dp
         !    DO ik = 1, NM
         !       IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) CYCLE
         !       denom = denom + BF_grad_res(ik, 1)*(HESSDI_scaled(ik)*BF_grad_res(ik, 1))
         !    END DO
         !    denom = MAX(denom, TINYD)
         !    IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
         !    tau = BPK/denom
         !    IF (.NOT. ieee_is_finite(tau) .OR. tau <= 0.0_dp) tau = 1.0_dp
         !    DO ik = 1, NM
         !       IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
         !          BF_r(ik) = BF_q(ik)
         !       ELSE
         !          BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
         !       END IF
         !    END DO
         ! ELSE
            yy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
            yy = MAX(yy, TINYD)
            IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
            gamma = BPK/yy
            IF (.NOT. ieee_is_finite(gamma) .OR. gamma <= 0.0_dp) gamma = 1.0_dp
            BF_r(1:NM) = gamma*BF_q(1:NM)
         ! END IF

         ! forward loop (full)
         DO im = m_LBFGS, 1, -1
            col = im
            BP = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_r(1:NM))       ! y·r
            IF (.NOT. ieee_is_finite(BF_p(im))) CYCLE
            IF (.NOT. ieee_is_finite(BP)) CYCLE
            BF_b = BF_p(im)*BP
            IF (.NOT. ieee_is_finite(BF_b)) CYCLE
            BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, col)*(BF_a(im) - BF_b)
         END DO

         p_k(1:NM) = -BF_r(1:NM)

         IF (USE_LBFGS_TYPE >= 2) THEN
            pnorm_tr = DNRM2(NM, p_k, 1)
            IF (pnorm_tr > DELTA_TR) THEN
               scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
               p_k(1:NM) = p_k(1:NM)*scale_tr
               IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(full): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
            END IF
         END IF

         dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)
         IF (.NOT. ieee_is_finite(dot_g_p) .OR. .NOT. ieee_is_finite(pnorm)) THEN
            IF (my_rank == 0) WRITE (*, *) 'LBFGS produced non-finite step → fallback steepest.'

            call steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                               GRAD_scaled, HESSDI_scaled, &
                               USE_PRECOND, my_rank, p_k)
         END IF
         RETURN
      END IF

   END SUBROUTINE lbfgs_step_loop

   SUBROUTINE GRID2D_UPDATE_MODEL(p_kf, m_fine, CR, CI, INVP, IANISO, ALPHA, &
                                  NPT, PAR_SCALE, NPAR, F_SCALE)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NPT, NPAR, IANISO
      INTEGER, INTENT(IN)    :: INVP(:)
      REAL(dp), INTENT(IN)    :: ALPHA
      REAL(dp), INTENT(IN)    :: PAR_SCALE(:), p_kf(:), F_SCALE(:)
      REAL(dp), INTENT(INOUT) :: CR(:, :), CI(:, :), m_fine(:)

      INTEGER :: I, IM, ps, pe

      IF (ALL(PAR_SCALE == 1.0_dp)) THEN
         IM = 0
         DO I = 1, NPAR
            IF (INVP(I) /= 1) CYCLE
            IM = IM + 1
            ps = (IM - 1)*NPT + 1
            pe = IM*NPT

            ! update grid baseline stack
            m_fine(ps:pe) = m_fine(ps:pe) + ALPHA*p_kf(ps:pe)*F_SCALE(I) !if F_SCALE**2 in Hessian

            ! map back to grid arrays for next forward step
            IF (I <= IANISO) THEN
               CR(I, 1:NPT) = m_fine(ps:pe)
            ELSE
               CI(I - (IANISO - 1), 1:NPT) = m_fine(ps:pe)
            END IF
         END DO
      ELSE
         IM = 0
         DO I = 1, NPAR
            IF (INVP(I) /= 1) CYCLE
            IM = IM + 1
            ps = (IM - 1)*NPT + 1
            pe = IM*NPT

            ! update grid baseline stack
            m_fine(ps:pe) = m_fine(ps:pe) + ALPHA*p_kf(ps:pe)*F_SCALE(I) !if F_SCALE**2 in Hessian

            ! map back to grid arrays for next forward step
            IF (I <= IANISO) THEN
               CR(I, 1:NPT) = m_fine(ps:pe)*PAR_SCALE(I)
            ELSE
               CI(I - (IANISO - 1), 1:NPT) = m_fine(ps:pe)*PAR_SCALE(I)
            END IF
         END DO
      END IF

      RETURN
   END SUBROUTINE

   SUBROUTINE GRID2D_UPDATE_LBFG(p_kf, CR, CI, INVP, IANISO, ALPHA, &
                                 NPT, SCALER, NPAR)

      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NPT, NPAR, IANISO
      INTEGER, INTENT(IN)    :: INVP(:)
      REAL(dp), INTENT(IN)    :: ALPHA
      REAL(dp), INTENT(IN)    :: SCALER(:), p_kf(:)
      REAL(dp), INTENT(INOUT) :: CR(:, :), CI(:, :)

      INTEGER :: I, IM, ps, pe

      IM = 0
      DO I = 1, NPAR
         IF (INVP(I) .EQ. 1) THEN
            IM = IM + 1
            ps = (IM - 1)*NPT + 1
            pe = IM*NPT

            ! map to grids for next iteration
            IF (I .LE. IANISO) THEN
               CR(I, 1:NPT) = ALPHA*p_kf(ps:pe)*SCALER(I)
            ELSEIF (I .GT. IANISO) THEN
               CI(I - (IANISO - 1), 1:NPT) = ALPHA*p_kf(ps:pe)*SCALER(I)
            END IF
         END IF
      END DO

      RETURN
   END SUBROUTINE GRID2D_UPDATE_LBFG

   SUBROUTINE lbfgs_project_direction(NPAR, INVP, NBLOCK, NPT, NNX, NNZ, NX, NZ, NTO, XTO, ZTO, &
                                      IG, NORD, p_k, p_kf, ITER, my_rank, unit_block, unit_point, &
                                      DEBUG_OUTPUT, pnorm_out)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPAR, NBLOCK, NPT, NNX, NNZ, NX, NZ, NTO, ITER, my_rank
      INTEGER, INTENT(IN) :: INVP(:), NORD
      TYPE(InversionGridType), INTENT(IN) :: IG
      REAL(dp), INTENT(INOUT) :: p_k(:)
      REAL(dp), INTENT(INOUT) :: p_kf(:)
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:)
      INTEGER, OPTIONAL, INTENT(IN) :: unit_block, unit_point
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT
      REAL(dp), OPTIONAL, INTENT(OUT) :: pnorm_out

      INTEGER :: II, IM, cs, ce, ps, pe, u_blk, u_pt, nmm
      LOGICAL :: dbg
      REAL(dp) :: pnorm_local, proj_scaler

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      u_blk = -1; IF (PRESENT(unit_block)) u_blk = unit_block
      u_pt = -1; IF (PRESENT(unit_point)) u_pt = unit_point

      p_kf = 0.0_dp

      CALL map_block2npt_all(NPAR, INVP, NPT, NBLOCK, NORD, NNZ, IG%N0_BLOCK, &
                             p_k, p_kf)

      nmm = SIZE(p_kf)
      pnorm_local = DNRM2(nmm, p_kf, 1)
      proj_scaler = DNRM2(SIZE(p_k), p_k, 1)/MAX(pnorm_local, TINY(1.0_dp))
      IF (PRESENT(pnorm_out)) pnorm_out = pnorm_local
      if (my_rank==0) WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') 'min(p_k)=', MINVAL(p_k), 'max(p_k)=', MAXVAL(p_k), '||p_k||=', DNRM2(SIZE(p_k), p_k, 1)
      if (my_rank==0) WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') 'min(p_kf)=', MINVAL(p_kf), 'max(p_kf)=', MAXVAL(p_kf), '||p_kf||=', pnorm_local
      if (my_rank == 0) WRITE (*, '(A,ES12.4)') 'p_k/p_kf  = ', DNRM2(SIZE(p_k), p_k, 1)/MAX(TINY(1.0_dp), pnorm_local)
      ! p_kf = p_kf*proj_scaler

      pnorm_local = DNRM2(nmm, p_kf, 1)
      CALL lbfgs_hist_flush_with_p_kf(my_rank, pnorm_local)

      IF (dbg .AND. my_rank == 0 .AND. ITER <= 5) THEN
         IF (u_blk > 0) THEN
            DO II = 1, SIZE(p_k)
               IF (p_k(II) /= 0.0_dp) WRITE (u_blk, '(I3, 2X, I10, 2X, E16.8)') ITER, II, p_k(II)
            END DO
         END IF
         IF (u_pt > 0) THEN
            DO II = 1, SIZE(p_kf)
               IF (p_kf(II) /= 0.0_dp) WRITE (u_pt, '(I3, 2X, I10, 2X, E16.8)') ITER, II, p_kf(II)
            END DO
         END IF
      END IF
   END SUBROUTINE lbfgs_project_direction

   SUBROUTINE cg_basic(iter, NM, GRAD_scaled, GRAD_scaled_norm, &
                       HESS_scaled, p_k, grad_prev, USE_PRECOND, &
                       my_rank)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: dp => real64
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)    :: iter       ! iteration index (>=1)
      INTEGER, INTENT(IN)    :: NM         ! number of unknowns
      REAL(dp), INTENT(IN)    :: GRAD_scaled(:)   ! current scaled gradient g_k
      REAL(dp), INTENT(OUT)   :: GRAD_scaled_norm ! ||g_eff|| (see below)
      REAL(dp), INTENT(IN)    :: HESS_scaled(:)   ! diag precond (inverse Hessian)
      REAL(dp), INTENT(INOUT) :: p_k(:)           ! search direction
      REAL(dp), INTENT(INOUT) :: grad_prev(:)      ! previous effective gradient
      LOGICAL, INTENT(IN)    :: USE_PRECOND      ! use HESS_scaled as precond
      INTEGER, INTENT(IN)    :: my_rank

      ! Locals
      INTEGER  :: i
      REAL(dp) :: gnorm2, gnorm2_old, beta
      REAL(dp) :: g_dot_p
      REAL(dp), ALLOCATABLE :: g_eff(:)   ! effective gradient (precond or not)

      INTEGER, PARAMETER :: restart_freq = 10
      REAL(dp), PARAMETER :: beta_min = 0.0_dp
      REAL(dp), PARAMETER :: beta_max = 1.0_dp

      ! Safety checks (optional)
      IF (SIZE(GRAD_scaled) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic: GRAD_scaled size mismatch: ', SIZE(GRAD_scaled), NM
      END IF
      IF (SIZE(p_k) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic: p_k size mismatch: ', SIZE(p_k), NM
      END IF
      IF (SIZE(grad_prev) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic: grad_prev size mismatch: ', SIZE(grad_prev), NM
      END IF
      IF (SIZE(HESS_scaled) /= NM .AND. USE_PRECOND .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic: HESS_scaled size mismatch: ', SIZE(HESS_scaled), NM
      END IF

      ! Allocate effective gradient
      ALLOCATE (g_eff(NM))

      ! Build effective gradient:
      !  - If USE_PRECOND: g_eff = M g = HESS_scaled * GRAD_scaled (inverse Hess diag)
      !  - Else:           g_eff = g
      IF (USE_PRECOND) THEN
         DO i = 1, NM
            g_eff(i) = HESS_scaled(i)*GRAD_scaled(i)
         END DO
      ELSE
         DO i = 1, NM
            g_eff(i) = GRAD_scaled(i)
         END DO
      END IF

      ! Compute ||g_eff||^2
      gnorm2 = 0.0_dp
      DO i = 1, NM
         gnorm2 = gnorm2 + g_eff(i)*g_eff(i)
      END DO
      GRAD_scaled_norm = SQRT(gnorm2)   ! norm of effective gradient

      ! --- First iteration: steepest descent in effective metric ---
      IF (iter <= 1) THEN

         DO i = 1, NM
            p_k(i) = -g_eff(i)        ! p_1 = -g_eff
            grad_prev(i) = g_eff(i)        ! store effective gradient as "previous"
         END DO

         DEALLOCATE (g_eff)
         RETURN
      END IF

      ! --- Fletcher–Reeves beta using effective gradients ---
      gnorm2_old = 0.0_dp
      DO i = 1, NM
         gnorm2_old = gnorm2_old + grad_prev(i)*grad_prev(i)
      END DO

      IF (gnorm2_old <= 0.0_dp) THEN
         beta = 0.0_dp
      ELSE
         beta = gnorm2/gnorm2_old
      END IF

      ! Periodic restart & safety bounds on beta
      IF (MOD(iter, restart_freq) == 0) beta = 0.0_dp
      IF (beta < beta_min .OR. beta > beta_max) beta = 0.0_dp

      ! Update search direction: p_k = -g_eff + beta * p_{k-1}
      DO i = 1, NM
         p_k(i) = -g_eff(i) + beta*p_k(i)
      END DO

      ! Descent check in true gradient metric: g^T p should be < 0
      g_dot_p = 0.0_dp
      DO i = 1, NM
         g_dot_p = g_dot_p + GRAD_scaled(i)*p_k(i)
      END DO

      IF (g_dot_p >= 0.0_dp) THEN
         IF (my_rank == 0) WRITE (*, *) 'cg_basic: g·p >= 0, restart to -g_eff.'
         DO i = 1, NM
            p_k(i) = -g_eff(i)
         END DO
         beta = 0.0_dp
      END IF

      ! Store current effective gradient for next iteration
      DO i = 1, NM
         grad_prev(i) = g_eff(i)
      END DO

      DEALLOCATE (g_eff)

   END SUBROUTINE cg_basic

   SUBROUTINE cg_basicpr(iter, NM, GRAD_scaled, GRAD_scaled_norm, &
                         HESS_scaled, p_k, grad_prev, g_prev, USE_PRECOND, &
                         my_rank)
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: dp => real64
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)     :: iter
      INTEGER, INTENT(IN)     :: NM
      REAL(dp), INTENT(IN)    :: GRAD_scaled(:)       ! g_k (scaled gradient)
      REAL(dp), INTENT(OUT)   :: GRAD_scaled_norm     ! ||z_k|| where z_k = M g_k (or g_k)
      REAL(dp), INTENT(IN)    :: HESS_scaled(:)       ! M (diag precond): approx inverse Hessian
      REAL(dp), INTENT(INOUT) :: p_k(:)               ! search direction (stores p_{k-1} on entry)
      REAL(dp), INTENT(INOUT) :: grad_prev(:)         ! stores z_{k-1} on entry; overwritten with z_k
      REAL(dp), INTENT(INOUT) :: g_prev(:)            ! stores g_{k-1} on entry; overwritten with g_k
      LOGICAL, INTENT(IN)     :: USE_PRECOND
      INTEGER, INTENT(IN)     :: my_rank

      ! Locals
      INTEGER  :: i
      REAL(dp) :: beta, denom, numer
      REAL(dp) :: g_dot_p
      REAL(dp) :: z_norm2
      REAL(dp), ALLOCATABLE :: z(:)     ! z_k = M g_k  (or g_k)
      REAL(dp), ALLOCATABLE :: dg(:)    ! g_k - g_{k-1}

      INTEGER, PARAMETER :: restart_freq = 10
      REAL(dp), PARAMETER :: tiny_denom = 1.0e-30_dp

      ! Safety checks (optional)
      IF (SIZE(GRAD_scaled) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic(PR+): GRAD_scaled size mismatch: ', SIZE(GRAD_scaled), NM
      END IF
      IF (SIZE(p_k) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic(PR+): p_k size mismatch: ', SIZE(p_k), NM
      END IF
      IF (SIZE(grad_prev) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic(PR+): grad_prev size mismatch: ', SIZE(grad_prev), NM
      END IF
      IF (SIZE(g_prev) /= NM .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic(PR+): g_prev size mismatch: ', SIZE(g_prev), NM
      END IF
      IF (SIZE(HESS_scaled) /= NM .AND. USE_PRECOND .AND. my_rank == 0) THEN
         WRITE (*, *) 'cg_basic(PR+): HESS_scaled size mismatch: ', SIZE(HESS_scaled), NM
      END IF

      ALLOCATE (z(NM))

      ! Build preconditioned (effective) gradient z_k
      IF (USE_PRECOND) THEN
         DO i = 1, NM
            z(i) = HESS_scaled(i)*GRAD_scaled(i)
         END DO
      ELSE
         DO i = 1, NM
            z(i) = GRAD_scaled(i)
         END DO
      END IF

      ! Norm of effective gradient (for logging)
      z_norm2 = 0.0_dp
      DO i = 1, NM
         z_norm2 = z_norm2 + z(i)*z(i)
      END DO
      GRAD_scaled_norm = SQRT(z_norm2)

      ! --- First iteration: steepest descent (preconditioned) ---
      IF (iter <= 1) THEN
         DO i = 1, NM
            p_k(i) = -z(i)
            grad_prev(i) = z(i)            ! store z_1
            g_prev(i) = GRAD_scaled(i)  ! store g_1
         END DO
         DEALLOCATE (z)
         RETURN
      END IF

      ! --- PR+ beta (preconditioned) ---
      ! beta_PR = ( z_k^T (g_k - g_{k-1}) ) / ( z_{k-1}^T g_{k-1} )
      ALLOCATE (dg(NM))
      DO i = 1, NM
         dg(i) = GRAD_scaled(i) - g_prev(i)
      END DO

      denom = 0.0_dp
      DO i = 1, NM
         denom = denom + grad_prev(i)*g_prev(i)   ! z_{k-1}^T g_{k-1}
      END DO

      IF (ABS(denom) <= tiny_denom) THEN
         beta = 0.0_dp
      ELSE
         numer = 0.0_dp
         DO i = 1, NM
            numer = numer + z(i)*dg(i)            ! z_k^T (g_k - g_{k-1})
         END DO
         beta = numer/denom
         IF (beta < 0.0_dp) beta = 0.0_dp           ! PR+
      END IF

      DEALLOCATE (dg)

      ! Periodic restart
      IF (MOD(iter, restart_freq) == 0) beta = 0.0_dp

      ! Update direction: p_k = -z_k + beta * p_{k-1}
      DO i = 1, NM
         p_k(i) = -z(i) + beta*p_k(i)
      END DO

      ! Descent check in true gradient metric (scaled gradient)
      g_dot_p = 0.0_dp
      DO i = 1, NM
         g_dot_p = g_dot_p + GRAD_scaled(i)*p_k(i)
      END DO

      IF (g_dot_p >= 0.0_dp) THEN
         IF (my_rank == 0) WRITE (*, *) 'cg_basic(PR+): g·p >= 0, restart to -z.'
         DO i = 1, NM
            p_k(i) = -z(i)
         END DO
         beta = 0.0_dp
      END IF

      ! Store for next iteration
      DO i = 1, NM
         grad_prev(i) = z(i)            ! z_k
         g_prev(i) = GRAD_scaled(i)  ! g_k
      END DO

      DEALLOCATE (z)

   END SUBROUTINE cg_basicpr

! subroutine CG_solve_GN_basic(ND, NM, FRECHET, b, x, maxit, tol, it_used, rnorm)
!   use iso_fortran_env, only: dp => real64
!   implicit none

!   integer, intent(in) :: ND, NM, maxit
!   complex(dp), intent(in) :: FRECHET(ND, NM)
!   real(dp), intent(in) :: b(NM)          ! RHS (for GN: b = -g)
!   real(dp), intent(inout) :: x(NM)       ! solution (init to 0)
!   real(dp), intent(in) :: tol
!   integer, intent(out) :: it_used
!   real(dp), intent(out) :: rnorm

!   real(dp), allocatable :: r(:), p(:), Ap(:)
!   complex(dp), allocatable :: y(:)
!   real(dp) :: rr, rr0, rr_new, alpha, beta, denom
!   integer :: k

!   allocate(r(NM), p(NM), Ap(NM))
!   allocate(y(ND))

!   ! x assumed provided; for GN simplest: x=0
!   call GN_apply_Ap_basic(ND, NM, FRECHET, x, Ap, y)
!   r = b - Ap
!   p = r

!   rr  = dot_product(r, r)
!   rr0 = rr
!   if (rr0 <= 1.0E-300_dp) then
!      it_used = 0
!      rnorm = 0.0_dp
!      deallocate(r,p,Ap,y)
!      return
!   end if

!   do k = 1, maxit
!      call GN_apply_Ap_basic(ND, NM, FRECHET, p, Ap, y)

!      denom = dot_product(p, Ap)
!      if (abs(denom) <= 1.0E-300_dp) exit

!      alpha = rr / denom
!      x = x + alpha * p
!      r = r - alpha * Ap

!      rr_new = dot_product(r, r)
!      if (rr_new <= (tol*tol) * rr0) exit

!      beta = rr_new / rr
!      p = r + beta * p
!      rr = rr_new
!   end do

!   it_used = k
!   rnorm = sqrt(rr_new / rr0)

!   deallocate(r,p,Ap,y)
! end subroutine CG_solve_GN_basic

!=======================================================================
!  CG_solve_GN
!-----------------------------------------------------------------------
!  Purpose:
!    Solve the (damped/undamped) Gauss–Newton normal equation system
!      A * dm = b
!    using (P)CG with matrix-free operator products from explicit FRECHET.
!
!    BASIC version: NO MPI, NO OpenMP, NO line search.
!    Operator (undamped): A(p) = Re( J^H * (J * p) )
!    Optional damping   : A(p) = Re( J^H * (J * p) ) + lambda * p
!    Optional diag precond: z = M^{-1} r using HESS_scaled(:) provided
!      (interpreted as an inverse-diagonal preconditioner).
!
!  Inputs:
!    ND, NM        : data/model sizes
!    FRECHET       : J (ND x NM) complex(dp)
!    b             : RHS (NM) real(dp)   (for GN: b = -g)
!    HESS_scaled   : (NM) real(dp) diagonal preconditioner (inv diag) if used
!    USE_PRECOND   : enable preconditioning
!    lambda        : damping (>=0). Set 0 for pure GN.
!    maxit         : maximum CG/PCG iterations (e.g., 10-30)
!    tol_rel       : relative tolerance on preconditioned residual norm
!                    stopping when sqrt(rz/rz0) <= tol_rel
!
!  In/Out:
!    dm            : (NM) real(dp)
!                    On entry: initial guess (usually zeros)
!                    On exit : solution estimate
!
!  Outputs:
!    it_used       : iterations used
!    relres        : achieved relative residual (sqrt(rz/rz0))
!=======================================================================

   subroutine CG_solve_GN(ND, NM, FRECHET, b, HESS_scaled, USE_PRECOND, lambda, &
                          dm, maxit, tol_rel, it_used, relres)

      use iso_fortran_env, only: dp => real64
      implicit none

      integer, intent(in) :: ND, NM, maxit
      complex(dp), intent(in) :: FRECHET(ND, NM)
      real(dp), intent(in) :: b(NM)
      real(dp), intent(in) :: HESS_scaled(NM)
      logical, intent(in) :: USE_PRECOND
      real(dp), intent(in) :: lambda
      real(dp), intent(inout) :: dm(NM)
      real(dp), intent(in) :: tol_rel
      integer, intent(out) :: it_used
      real(dp), intent(out) :: relres

      ! Locals
      real(dp), allocatable :: r(:), z(:), p(:), Ap(:)
      complex(dp), allocatable :: y(:)
      real(dp) :: rz, rz0, rz_new, alpha, beta, denom
      integer :: k
      real(dp), parameter :: tiny = 1.0e-300_dp

      allocate (r(NM), z(NM), p(NM), Ap(NM))
      allocate (y(ND))

      !---------------------------------------------------------------
      ! r = b - A*dm
      !---------------------------------------------------------------
      call GN_vec_product(ND, NM, FRECHET, dm, Ap, y)   ! Ap = Re(J^H J dm)
      if (lambda > 0.0_dp) Ap = Ap + lambda*dm
      r = b - Ap

      !---------------------------------------------------------------
      ! z = M^{-1} r  (or z=r)
      !---------------------------------------------------------------
      if (USE_PRECOND) then
         z = HESS_scaled*r
      else
         z = r
      end if

      p = z
      rz = dot_product(r, z)
      rz0 = rz

      if (abs(rz0) <= tiny) then
         it_used = 0
         relres = 0.0_dp
         deallocate (r, z, p, Ap, y)
         return
      end if

      relres = 1.0_dp

      !---------------------------------------------------------------
      ! PCG iterations
      !---------------------------------------------------------------
      do k = 1, maxit

         ! Ap = A * p
         call GN_vec_product(ND, NM, FRECHET, p, Ap, y)   ! Ap = Re(J^H J p)
         if (lambda > 0.0_dp) Ap = Ap + lambda*p

         denom = dot_product(p, Ap)
         if (abs(denom) <= tiny) exit

         alpha = rz/denom

         dm = dm + alpha*p
         r = r - alpha*Ap

         ! z = M^{-1} r
         if (USE_PRECOND) then
            z = HESS_scaled*r
         else
            z = r
         end if

         rz_new = dot_product(r, z)

         relres = sqrt(max(rz_new, 0.0_dp)/rz0)
         if (relres <= tol_rel) then
            rz = rz_new
            exit
         end if

         beta = rz_new/rz
         p = z + beta*p

         rz = rz_new

      end do

      it_used = k

      deallocate (r, z, p, Ap, y)

   end subroutine CG_solve_GN
   SUBROUTINE lbfgs_hist_init(my_rank)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: my_rank
      INTEGER :: ios
      CHARACTER(len=200) :: iomsg

      IF (my_rank /= 0 .OR. lbfgs_hist_ready) RETURN

      OPEN (NEWUNIT=lbfgs_hist_unit, FILE='out_lbfgs_hist.txt', STATUS='UNKNOWN', &
            POSITION='APPEND', ACTION='WRITE', IOSTAT=ios, IOMSG=iomsg)
      IF (ios == 0) THEN
         WRITE (lbfgs_hist_unit, '(A)') '# iter  m  ||g||        ||p||        g·p          '// &
            '||s||       ||y||       s·y         desc_cos    cos_sy      gamma0      '// &
            'sTy_min     curv_ok  hist_used  hist_skip  bad_col  bad_reason  '// &
            'bad_s·y     bad_cos_sy  bad_||s||   bad_||y||   ||p_kf||    ||p_kf||/||p||'
         FLUSH (lbfgs_hist_unit)
         lbfgs_hist_ready = .TRUE.
      ELSE
         WRITE (*, *) 'LBFGS: cannot open out_lbfgs_hist.txt: ', TRIM(iomsg)
      END IF
   END SUBROUTINE lbfgs_hist_init

   SUBROUTINE lbfgs_hist_reset_pairs()
      IMPLICIT NONE

      lbfgs_hist_used_pairs = 0
      lbfgs_hist_skipped_pairs = 0
      lbfgs_hist_bad_pair_col = 0
      lbfgs_hist_bad_reason = 0
      lbfgs_hist_bad_bp = 0.0_real64
      lbfgs_hist_bad_snorm = 0.0_real64
      lbfgs_hist_bad_ynorm = 0.0_real64
      lbfgs_hist_bad_cos = 0.0_real64
   END SUBROUTINE lbfgs_hist_reset_pairs

   SUBROUTINE lbfgs_hist_note_pair(pair_col, bp, snorm, ynorm, used_pair, reason_code)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: pair_col, reason_code
      REAL(real64), INTENT(IN) :: bp, snorm, ynorm
      LOGICAL, INTENT(IN) :: used_pair
      REAL(real64) :: denom_sy, cos_sy

      denom_sy = MAX(snorm*ynorm, tiny_dp)
      cos_sy = bp/denom_sy

      IF (used_pair) THEN
         lbfgs_hist_used_pairs = lbfgs_hist_used_pairs + 1
         RETURN
      END IF

      lbfgs_hist_skipped_pairs = lbfgs_hist_skipped_pairs + 1
      IF (lbfgs_hist_bad_pair_col == 0 .OR. ABS(bp) < ABS(lbfgs_hist_bad_bp)) THEN
         lbfgs_hist_bad_pair_col = pair_col
         lbfgs_hist_bad_reason = reason_code
         lbfgs_hist_bad_bp = bp
         lbfgs_hist_bad_snorm = snorm
         lbfgs_hist_bad_ynorm = ynorm
         lbfgs_hist_bad_cos = cos_sy
      END IF
   END SUBROUTINE lbfgs_hist_note_pair

   SUBROUTINE lbfgs_hist_report_pairs(my_rank)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: my_rank

      IF (my_rank /= 0) RETURN
      IF (lbfgs_hist_skipped_pairs <= 0) RETURN

      WRITE (*, '(A,I2,2X,A,I2,2X,A,I2,2X,A,1PE12.4,2X,A,1PE12.4)') &
         'LBFGS history: used=', lbfgs_hist_used_pairs, &
         'skipped=', lbfgs_hist_skipped_pairs, &
         'worst_col=', lbfgs_hist_bad_pair_col, &
         'bad_s·y=', lbfgs_hist_bad_bp, &
         'bad_cos_sy=', lbfgs_hist_bad_cos
   END SUBROUTINE lbfgs_hist_report_pairs

   SUBROUTINE lbfgs_hist_stage(iter, m_hist, gnorm, pnorm, gdotp, snorm, ynorm, sty, sty_min, curv_ok)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iter, m_hist, curv_ok
      REAL(real64), INTENT(IN) :: gnorm, pnorm, gdotp, snorm, ynorm, sty, sty_min
      REAL(real64) :: denom_gp, denom_sy, ysq

      lbfgs_hist_iter = iter
      lbfgs_hist_m = m_hist
      lbfgs_hist_gnorm = gnorm
      lbfgs_hist_pnorm = pnorm
      lbfgs_hist_gdotp = gdotp
      lbfgs_hist_snorm = snorm
      lbfgs_hist_ynorm = ynorm
      lbfgs_hist_sty = sty
      lbfgs_hist_sty_min = sty_min
      lbfgs_hist_curv_ok = curv_ok

      denom_gp = MAX(gnorm*pnorm, tiny_dp)
      lbfgs_hist_desc_cos = -gdotp/denom_gp

      denom_sy = MAX(snorm*ynorm, tiny_dp)
      lbfgs_hist_cos_sy = sty/denom_sy

      ysq = ynorm*ynorm
      IF (ysq > tiny_dp) THEN
         lbfgs_hist_gamma0 = sty/ysq
      ELSE
         lbfgs_hist_gamma0 = 0.0_real64
      END IF

      lbfgs_hist_pending = .TRUE.
   END SUBROUTINE lbfgs_hist_stage

   SUBROUTINE lbfgs_hist_flush_with_p_kf(my_rank, pnorm_kf)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: my_rank
      REAL(real64), INTENT(IN) :: pnorm_kf
      REAL(real64) :: ratio

      IF (my_rank /= 0) RETURN
      IF (.NOT. lbfgs_hist_pending) RETURN
      CALL lbfgs_hist_init(my_rank)
      IF (.NOT. lbfgs_hist_ready) RETURN

      IF (lbfgs_hist_pnorm > tiny_dp) THEN
         ratio = pnorm_kf/lbfgs_hist_pnorm
      ELSE
         ratio = 0.0_real64
      END IF

      WRITE (lbfgs_hist_unit, '(I5,1X,I2,1P,1X,10E13.5,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2,1X,4E13.5,1X,2E13.5)') &
         lbfgs_hist_iter, lbfgs_hist_m, &
         lbfgs_hist_gnorm, lbfgs_hist_pnorm, lbfgs_hist_gdotp, &
         lbfgs_hist_snorm, lbfgs_hist_ynorm, lbfgs_hist_sty, &
         lbfgs_hist_desc_cos, lbfgs_hist_cos_sy, lbfgs_hist_gamma0, &
         lbfgs_hist_sty_min, lbfgs_hist_curv_ok, &
         lbfgs_hist_used_pairs, lbfgs_hist_skipped_pairs, &
         lbfgs_hist_bad_pair_col, lbfgs_hist_bad_reason, &
         lbfgs_hist_bad_bp, lbfgs_hist_bad_cos, lbfgs_hist_bad_snorm, lbfgs_hist_bad_ynorm, &
         pnorm_kf, ratio
      FLUSH (lbfgs_hist_unit)
      lbfgs_hist_pending = .FALSE.
   END SUBROUTINE lbfgs_hist_flush_with_p_kf
!=======================================================================
!  GN_vec_product
!-----------------------------------------------------------------------
!  Purpose:
!    Compute the Gauss–Newton normal-operator vector product
!
!      Hp = Re( J^H * (J * p_k) )
!
!    where:
!      J      = FRECHET (ND x NM) complex(dp)
!      p_k    = real(dp) model-space vector (NM)
!      y      = complex(dp) workspace vector (ND) : y = J p_k
!      Hp     = real(dp) output vector (NM)
!
!  Notes:
!    - Basic loop implementation (no BLAS, no parallel).
!    - A BLAS ZGEMV-based implementation is included below (commented out).
!=======================================================================

   subroutine GN_vec_product(ND, NM, FRECHET, p_k, Hp, y)
      use iso_fortran_env, only: dp => real64
      implicit none

      integer, intent(in) :: ND, NM
      complex(dp), intent(in)  :: FRECHET(ND, NM)  ! J
      real(dp), intent(in)  :: p_k(NM)          ! input vector
      real(dp), intent(out) :: Hp(NM)           ! Hp = Re(J^H J p_k)
      complex(dp), intent(out) :: y(ND)            ! workspace: y = J p_k

      integer :: i, j
      complex(dp) :: acc
      ! print *, '||p|| = ', dnrm2(NM, p_k, 1)
      !----------------------------
      ! 1) y = J * p_k
      !----------------------------
      do i = 1, ND
         acc = (0.0_dp, 0.0_dp)
         do j = 1, NM
            acc = acc + FRECHET(i, j)*p_k(j)
         end do
         y(i) = acc
      end do
!  -- debug prints removed for production runs
      !----------------------------
      ! 2) Hp = Re( J^H * y )
      !----------------------------
      do j = 1, NM
         acc = (0.0_dp, 0.0_dp)
         do i = 1, ND
            acc = acc + conjg(FRECHET(i, j))*y(i)
         end do
         Hp(j) = real(acc, dp)
      end do
!  -- debug prints removed for production runs
      return

      !=====================================================================
      !  BLAS version (ZGEMV) - COMMENTED OUT
      !
      !  This avoids explicit double loops and is usually much faster with MKL.
      !  It does two matvecs:
      !    y   = J * p_k
      !    tmp = J^H * y
      !    Hp  = Re(tmp)
      !
      !  Important:
      !    - ZGEMV requires complex input vector; we promote p_k to complex.
      !    - tmp is complex(NM) workspace.
      !    - Ensure you have a BLAS interface (MKL) available in your build.
      !=====================================================================

      !use mkl_blas95, only: gemv   ! (Option A: BLAS95)
      ! ! or declare external ZGEMV for F77 BLAS (Option B)
      ! ! external zgemv

      !complex(dp), allocatable :: pk_c(:), tmp(:)
      !complex(dp) :: alpha, beta
      !integer :: j

      !alpha = (1.0_dp, 0.0_dp)
      !beta  = (0.0_dp, 0.0_dp)

      !allocate(pk_c(NM), tmp(NM))
      !pk_c(:) = cmplx(p_k(:), 0.0_dp, kind=dp)

      ! ! y = J * pk_c
      !call zgemv('N', ND, NM, alpha, FRECHET, ND, pk_c, 1, beta, y, 1)

      ! ! tmp = J^H * y
      !call zgemv('C', ND, NM, alpha, FRECHET, ND, y, 1, beta, tmp, 1)

      !do j = 1, NM
      !   Hp(j) = real(tmp(j), dp)
      !end do

      !deallocate(pk_c, tmp)

   end subroutine GN_vec_product

!-----------------------------------------------------------------------
! Gauss-Newton step wrapper
! - Builds RHS b = -g, solves A dm = b with CG_solve_GN (matrix-free)
! - Maps dm (block space) to point space via existing projector
! - Applies a conservative damping/backtracking on alpha to avoid huge steps
! Note: This wrapper does NOT perform forward-model cost evaluation; it
!       returns a recommended step (`alpha`) and the mapped step `p_kf`.
!-----------------------------------------------------------------------
   subroutine GaussNewton_step(ND, NM, FRECHET, GRAD_scaled, HESSDI_scaled, USE_PRECOND, lambda, maxit, tol_rel, &
                               dm, p_k, p_kf, NPAR, INVP, NBLOCK, NPT, SCALER, IANISO, CRR0, CII0, m_fine, &
                               NNX, NNZ, NX, NZ, NTO, XTO, ZTO, IG, NORD, alpha_out, it_used, relres, my_rank)

      use iso_fortran_env, only: dp => real64
      implicit none

      integer, intent(in) :: ND, NM, maxit, NPAR, NBLOCK, NPT, IANISO, my_rank
      integer, intent(in) :: NNX, NNZ, NX, NZ, NTO, NORD
      real(dp), intent(in) :: XTO(:), ZTO(:)
      TYPE(InversionGridType), INTENT(IN) :: IG
      complex(dp), intent(in) :: FRECHET(ND, NM)
      real(dp), intent(in) :: GRAD_scaled(NM), HESSDI_scaled(NM)
      logical, intent(in) :: USE_PRECOND
      real(dp), intent(in) :: lambda, tol_rel

      ! outputs / inout
      real(dp), intent(inout) :: dm(NM)         ! initial guess on entry, solution on exit
      real(dp), intent(inout) :: p_k(NM)        ! block-space step (will be set to dm)
      real(dp), intent(out) :: p_kf(:)          ! mapped point-space step
      real(dp), intent(in) :: SCALER(:)
      real(dp), intent(inout) :: m_fine(:)
      real(dp), intent(inout) :: CRR0(:, :), CII0(:, :)

      real(dp), intent(out) :: alpha_out
      integer, intent(out) :: it_used
      real(dp), intent(out) :: relres
      integer, intent(in) :: INVP(:)

      real(dp), allocatable :: b_local(:)
      real(dp) :: alpha
      real(dp) :: max_abs, max_allowed
      integer :: i
      integer :: iter_local, unit_block, unit_point
      real(dp) :: pnorm_out_local
      logical :: debug_local

      ! Build RHS: b = -g
      allocate (b_local(NM))
      do i = 1, NM
         b_local(i) = -GRAD_scaled(i)
      end do

      ! Initialize dm if not already
      dm = 0.0_dp

      ! Solve GN normal equations (matrix-free PCG)
      call CG_solve_GN(ND, NM, FRECHET, b_local, HESSDI_scaled, USE_PRECOND, lambda, dm, maxit, tol_rel, it_used, relres)

      ! Use dm as block-space step
      p_k(1:NM) = dm(1:NM)
      p_kf = 0.0_dp

      ! Map to point space via existing projector. Caller must allocate p_kf sized correctly.
      ! build locals for projector optional args
      iter_local = 0
      unit_block = -1
      unit_point = -1
      debug_local = .FALSE.
      pnorm_out_local = 0.0_dp

      call lbfgs_project_direction(NPAR, INVP, NBLOCK, NPT, NNX, NNZ, NX, NZ, NTO, XTO, ZTO, &
                                   IG, NORD, p_k, p_kf, ITER=iter_local, my_rank=my_rank, &
                                   unit_block=unit_block, unit_point=unit_point, &
                                   DEBUG_OUTPUT=debug_local, pnorm_out=pnorm_out_local)

      ! Conservative backtracking on alpha based on max absolute mapped step
      alpha = 1.0_dp
      max_allowed = 1.0_dp   ! default max absolute step (tunable)
      do
         max_abs = 0.0_dp
         do i = 1, size(p_kf)
            max_abs = max(max_abs, abs(p_kf(i)*alpha))
         end do
         if (max_abs <= max_allowed .or. alpha <= 1.0e-4_dp) exit
         alpha = alpha*0.5_dp
      end do

      alpha_out = alpha

      ! Do not actually update CR/CI here; caller may accept step and call GRID2D_UPDATE_MODEL/GRID2D_UPDATE_LBFG
      if (allocated(b_local)) deallocate (b_local)
      return
   end subroutine GaussNewton_step

!------------------------------------------------------------------------------
!  Subroutine: UpdateAdaptiveTR
!------------------------------------------------------------------------------
!  Purpose:
!    Maintain EMAs of ||g|| and accepted ||p|| and map them to an adaptive
!    trust-region radius for clamping the L-BFGS search direction.
!
!  Strategy:
!    - Gradient-driven radius:   DELTA_TR_G = BASE * K_G / max(EMA(||g||), eps)
!    - Step-driven radius:       DELTA_TR_P = K_P * EMA(||p||_accepted)
!    - Combined recommendation:  DELTA_TR   = min(DELTA_TR_G, DELTA_TR_P) with
!                                   optional [CLAMP_MIN, CLAMP_MAX] bounds.
!
!  Notes:
!    * Persistence: EMAs are kept with SAVE so values persist across calls.
!    * Call twice per iteration if you want both updates:
!        (a) before line-search (have_acc_step=.FALSE.) to get a proposed TR
!        (b) after line-search accept (have_acc_step=.TRUE., pnorm_acc=...)
!    * If you want to reset between frequency batches, pass RESET=.TRUE.
!
!  Entries:
!    gnorm          [in]  : Euclidean norm of current gradient ||g|| (scaled space)
!    have_acc_step  [in]  : .TRUE. iff a step was accepted and pnorm_acc is valid
!    pnorm_acc      [in]  : Norm of the accepted step ||p|| (ignored if have_acc_step=.FALSE.)
!    delta_tr_out   [out] : Recommended trust-region radius for clamping ||p||
!
!  Optional tuning (all OPTIONAL; reasonable defaults provided):
!    reset          [in]  : .TRUE. to clear EMAs this call
!    alpha_g_in     [in]  : EMA weight for gradient norm (default 0.20D0)
!    alpha_p_in     [in]  : EMA weight for step norm     (default 0.30D0)
!    k_g_in         [in]  : Gain for gradient-driven TR  (default 1.0D0)
!    k_p_in         [in]  : Gain for step-driven TR      (default 1.25D0)
!    base_in        [in]  : Baseline TR scale (scaled space) (default 1.0D0)
!    clamp_min_in   [in]  : Lower bound on DELTA_TR (default 0.0D0 => no lower clamp)
!    clamp_max_in   [in]  : Upper bound on DELTA_TR (default huge => no upper clamp)
!
!  Return:
!    delta_tr_out         : Final recommended radius after min-combine and clamping
!

!------------------------------------------------------------------------------
! SUBROUTINE UpdateAdaptiveTR(gnorm, have_acc_step, pnorm_acc, delta_tr_out, &
!                             reset, alpha_g_in, alpha_p_in, k_g_in, k_p_in,  &
!                             base_in, clamp_min_in, clamp_max_in)

!   IMPLICIT real(dp) (A-H,O-Z)

!   ! ---- required ----
!   real(dp), INTENT(IN)  :: gnorm
!   LOGICAL,          INTENT(IN)  :: have_acc_step
!   real(dp), INTENT(IN)  :: pnorm_acc
!   real(dp), INTENT(OUT) :: delta_tr_out

!   ! ---- optional ----
!   LOGICAL,          INTENT(IN),  OPTIONAL :: reset
!   real(dp), INTENT(IN),  OPTIONAL :: alpha_g_in, alpha_p_in
!   real(dp), INTENT(IN),  OPTIONAL :: k_g_in, k_p_in, base_in
!   real(dp), INTENT(IN),  OPTIONAL :: clamp_min_in, clamp_max_in

!   ! ---- persistent EMAs ----
!   real(dp), SAVE :: gnorm_ema = 0.0D0
!   real(dp), SAVE :: pnorm_ema = 0.0D0

!   ! ---- locals / params ----
!   real(dp) :: ALPHA_G, ALPHA_P, K_G, K_P, BASE_TR
!   real(dp) :: CLAMP_MIN, CLAMP_MAX
!   real(dp), PARAMETER :: G_EPS = 1.0D-16
!   real(dp), PARAMETER :: P_EPS = 1.0D-16
!   real(dp) :: delta_tr_g, delta_tr_p, delta_tr_comb
!   LOGICAL :: do_reset

!   ! ---- defaults ----
!   ALPHA_G   = 0.20D0
!   ALPHA_P   = 0.30D0
!   K_G       = 1.00D0
!   K_P       = 1.25D0
!   BASE_TR   = 1.00D0
!   CLAMP_MIN = 0.00D0
!   CLAMP_MAX = HUGE(1.0D0)

!   IF (PRESENT(alpha_g_in))   ALPHA_G   = alpha_g_in
!   IF (PRESENT(alpha_p_in))   ALPHA_P   = alpha_p_in
!   IF (PRESENT(k_g_in))       K_G       = k_g_in
!   IF (PRESENT(k_p_in))       K_P       = k_p_in
!   IF (PRESENT(base_in))      BASE_TR   = base_in
!   IF (PRESENT(clamp_min_in)) CLAMP_MIN = clamp_min_in
!   IF (PRESENT(clamp_max_in)) CLAMP_MAX = clamp_max_in
!   do_reset = .FALSE.; IF (PRESENT(reset)) do_reset = reset

!   ! ---- reset EMAs if requested ----
!   IF (do_reset) THEN
!      gnorm_ema = 0.0D0
!      pnorm_ema = 0.0D0
!   END IF

!   ! ---- update EMA(||g||) ----
!   IF (gnorm_ema == 0.0D0) THEN
!      gnorm_ema = gnorm
!   ELSE
!      gnorm_ema = (1.0D0 - ALPHA_G) * gnorm_ema + ALPHA_G * gnorm
!   END IF

!   ! ---- optionally update EMA(||p||_accepted) ----
!   IF (have_acc_step) THEN
!      IF (pnorm_ema == 0.0D0) THEN
!         pnorm_ema = pnorm_acc
!      ELSE
!         pnorm_ema = (1.0D0 - ALPHA_P) * pnorm_ema + ALPHA_P * pnorm_acc
!      END IF
!   END IF

!   ! ---- map EMAs to radii ----
!   delta_tr_g = BASE_TR * K_G / MAX(gnorm_ema, G_EPS)
!   delta_tr_p = MAX(K_P     * pnorm_ema,      P_EPS)

!   ! ---- combine & clamp ----
!   delta_tr_comb = MIN(delta_tr_g, delta_tr_p)
!   delta_tr_out  = MIN(MAX(delta_tr_comb, CLAMP_MIN), CLAMP_MAX)

! END SUBROUTINE UpdateAdaptiveTR

end module lbfgs_mod
