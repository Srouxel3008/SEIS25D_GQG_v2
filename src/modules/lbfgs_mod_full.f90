module lbfgs_mod
   use gradient_mod
   use Frechet_mod
   use shared_mod
   use grid_mod
   use iso_fortran_env, only: dp => real64
   use constant_mod, only: tiny_dp
   implicit None

contains
   ! Module for L-BFGS optimization
   ! This module contains the lbfgs_step subroutine and other related functions:
   !1. lbfgs_step,  performs the L-BFGS update step.
   !2. ComputeHess, inverse  Hessian approximation diagonal.
   !3. d
!

   SUBROUTINE lbfgs_step_loop(iter, NM, mml, INVP, NPAR, NBLOCK, &
                              GRAD_Scaled, GRAD_scaled_norm, HESSDI_scaled, &
                              p_k, grad_prev, BF_grad_res, &
                              BF_s_hist, m_coarse, m_coarse_prev, FCOST0, FCOST_prev, &
                              USE_SCALAR_H0, USE_LBFGS_TYPE, my_rank)

      USE, INTRINSIC :: ieee_arithmetic            ! [E] finite checks
      IMPLICIT NONE

      ! --------- dummy args ---------
      INTEGER, INTENT(IN)    :: iter, NM, mml, my_rank, INVP(:), NPAR, NBLOCK
      REAL(dp), INTENT(IN)    :: GRAD_Scaled(:)                  ! scaled gradient g_k
      REAL(dp), INTENT(IN)    :: GRAD_Scaled_Norm(:)             ! per-parameter norms
      REAL(dp), INTENT(IN)    :: HESSDI_scaled(:)          ! diagonal preconditioner in scaled space
      REAL(dp), INTENT(INOUT) :: p_k(:)
      REAL(dp), INTENT(IN) :: m_coarse(:), m_coarse_prev(:)
      REAL(dp), INTENT(IN)    :: grad_prev(:)              ! previous scaled gradient g_{k-1}
      REAL(dp), INTENT(INOUT) :: BF_grad_res(:, :)        ! Y history (y = g_k - g_{k-1})
      REAL(dp), INTENT(INOUT) :: BF_s_hist(:, :)           ! S history (s = x_k - x_{k-1})
      REAL(dp), INTENT(IN)    :: FCOST0
      REAL(dp), INTENT(IN)    :: FCOST_prev
      LOGICAL, INTENT(IN)    :: USE_SCALAR_H0
      INTEGER, INTENT(IN)    :: USE_LBFGS_TYPE

      ! --------- locals ----------
      REAL(dp) :: BP, BA, BPK, BB, BF_b
      INTEGER          :: iiter, ik, m_LBFGS, im, jcol, IA, i1, i2, col, IMAP
      REAL(dp) :: dot_g_p, gnorm, pnorm, gHg
      REAL(dp) :: BF_a(mml), BF_p(mml)
      REAL(dp) :: BF_q(NM), BF_r(NM), pscaler, precond_grad_norm, precond_grad(NM)
      REAL(dp) :: sTy, sTs, yTy, sTy_min
      REAL(dp) :: denom, tau, yy, gamma
      LOGICAL          :: restart_now
      REAL(dp), PARAMETER :: ZERO = 0.0_dp
      REAL(dp), PARAMETER :: C_CURV = 1.0e-8_dp
      REAL(dp), PARAMETER :: TINYD = 1.0e-36_dp

      ! --------- cost-damping / trust-region ----------
      REAL(dp) :: s2, theta, ys_raw, ys_hat, yy_damp
      REAL(dp) :: ytmp(NM), svec(NM), yhat(NM)
      REAL(dp), PARAMETER :: EPS_CURV = 1.0e-14_dp
      REAL(dp), PARAMETER :: DELTA_TR = 1.0_dp
      REAL(dp) :: pnorm_tr, scale_tr

      ! [B] denominator/perturbation floors
      REAL(dp), PARAMETER :: EPS_DEN = 1.0e-30_dp
      REAL(dp), PARAMETER :: EPS_THETA = 1.0e-30_dp

      ! INTERFACE
      !    REAL(dp) FUNCTION DNRM2(N, X, INCX)
      !       INTEGER :: N, INCX
      !       REAL(dp) :: X(*)
      !    END FUNCTION DNRM2
      ! END INTERFACE

      ! ! [E] Early finite checks on inputs (fallback if bad)
      ! IF (.NOT.all(ieee_is_finite(GRAD_Scaled))) THEN
      !    IF (my_rank == 0) WRITE(*,*) 'LBFGS: non-finite GRAD_Scaled → fallback steepest.'
      !    p_k(1:NM) = 0.0D0
      !    DO ik=1,NM
      !       p_k(ik) = -MERGE(HESSDI_scaled(ik)*GRAD_Scaled(ik), -GRAD_Scaled(ik), &
      !                        .NOT.ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0)
      !    END DO
      !    RETURN
      ! END IF

      gnorm = DNRM2(NM, GRAD_Scaled, 1)
      IF (.NOT. USE_SCALAR_H0) THEN

         gHg = 0.0D0
         DO ik = 1, NM
            IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) CYCLE
            gHg = gHg + GRAD_Scaled(ik)*(HESSDI_scaled(ik)*GRAD_Scaled(ik))
         END DO
      ELSE
         gHg = 0.0D0
      END IF

      IF (my_rank == 0) THEN
         WRITE (*, *) ' '
         WRITE (*, '(A,I6,2X,A,ES12.4)') 'LBFGS: iter=', iter, '||g||=', gnorm
         IF (.NOT. USE_SCALAR_H0) WRITE (*, '(A,ES12.4)') 'LBFGS: g^T (H_precond ⊙ g) =', gHg
      END IF

      iiter = iter - 1

      ! ==========================
      ! Iter 0: steepest direction
      ! ==========================
      IF (iiter == 0) THEN
         IF (USE_SCALAR_H0) THEN
            IMAP = 0
            DO IA = 1, NPAR
               IF (INVP(IA) == 1) THEN
                  IMAP = IMAP + 1
                  i1 = (IMAP - 1)*NBLOCK + 1
                  i2 = IMAP*NBLOCK
                  ! [E] guard GRAD_Scaled_Norm(IA)
                  pscaler = 1.0D0/MAX(DABS(GRAD_Scaled_Norm(IA)), EPS_DEN)
                  p_k(i1:i2) = -GRAD_Scaled(i1:i2)*pscaler
                  ! p_k(i1:i2) = -GRAD_Scaled(i1:i2)
                  IF (my_rank == 0) WRITE (*, '(A,I4,2X,A,ES12.4)') '  steepest (no-precond) IA=', IA, 'scale=', pscaler
               END IF
            END DO
         ELSE
            IMAP = 0
            DO IA = 1, NPAR
               IF (INVP(IA) == 1) THEN
                  IMAP = IMAP + 1
                  i1 = (IMAP - 1)*NBLOCK + 1
                  i2 = IMAP*NBLOCK
                  ! [E] sanitize preconditioner on the fly
                  DO ik = i1, i2
                     IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                        precond_grad(ik) = GRAD_Scaled(ik)               ! fallback to identity
                     ELSE
                        precond_grad(ik) = HESSDI_scaled(ik)*GRAD_Scaled(ik)
                     END IF
                  END DO
                  precond_grad_norm = MAX(DNRM2(NBLOCK, precond_grad(i1), 1), EPS_DEN)   ! [B]
                  pscaler = 1/precond_grad_norm     ! [B,E]MAX(GRAD_Scaled_Norm(IA), EPS_DEN)
                  p_k(i1:i2) = -pscaler*precond_grad(i1:i2)
                  IF (my_rank == 0) WRITE (*, '(A,I4,2X,A,ES12.4,2X,A,ES12.4)') &
                     '  steepest (precond)  IA=', IA, '||H_precond*g||_blk=', precond_grad_norm, 'full_scale=', pscaler
               END IF
            END DO
         END IF

         dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)
      IF (my_rank == 0) WRITE (*,'(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') 'iter0: min(p)=', MINVAL(p_k), 'max(p)=', MAXVAL(p_k), '||p||=', pnorm
         IF (.NOT. ieee_is_finite(dot_g_p) .OR. dot_g_p >= ZERO) THEN     ! [E]
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

         BF_grad_res(1:NM, col) = GRAD_Scaled(1:NM) - grad_prev(1:NM)   ! y_k
         BF_s_hist(1:NM, col) = m_coarse(1:NM) - m_coarse_prev(1:NM) ! s_k

         ! ---- cost-damped y (type>=1) ----
         IF (USE_LBFGS_TYPE >= 1) THEN
            ytmp(1:NM) = BF_grad_res(1:NM, col)
            svec(1:NM) = BF_s_hist(1:NM, col)
            s2 = DOT_PRODUCT(svec, svec)
            theta = 6.0D0*(FCOST_prev - FCOST0) + 3.0D0*DOT_PRODUCT(grad_prev(1:NM) + GRAD_Scaled(1:NM), svec(1:NM))
            ! [B] robust damping term
            IF (s2 > EPS_DEN) THEN
               yhat(1:NM) = ytmp(1:NM)
               IF (DABS(theta) > EPS_THETA) yhat(1:NM) = yhat(1:NM) + (theta/s2)*svec(1:NM)
            ELSE
               yhat(1:NM) = ytmp(1:NM)
            END IF

            ys_raw = DOT_PRODUCT(ytmp, svec)
            ys_hat = DOT_PRODUCT(yhat, svec)
            yy_damp = DOT_PRODUCT(yhat, yhat)

            IF (ys_hat <= EPS_CURV*MAX(s2, 1.0D0)) THEN
               yhat(1:NM) = yhat(1:NM) + ((EPS_CURV*MAX(s2, 1.0D0) - ys_hat)/MAX(s2, 1.0D0))*svec(1:NM)
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
         sTy_min = C_CURV*DSQRT(MAX(sTs, 0D0)*MAX(yTy, 0D0))
         restart_now = (.NOT. ieee_is_finite(sTy)) .OR. (sTy <= sTy_min)    ! [D]

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (warm): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0D0     ! [A] clear history on restart
            BF_s_hist(1:NM, 1:mml) = 0.0D0

            ! fallback steepest
            IF (USE_SCALAR_H0) THEN
               IMAP = 0
               DO IA = 1, NPAR
                  IF (INVP(IA) == 1) THEN
                     IMAP = IMAP + 1
                     i1 = (IMAP - 1)*NBLOCK + 1
                     i2 = IMAP*NBLOCK
                     pscaler = 1.0D0/MAX(DABS(GRAD_Scaled_Norm(IA)), EPS_DEN)     ! [B,E]

                     p_k(i1:i2) = -GRAD_Scaled(i1:i2)*pscaler
                     ! p_k(i1:i2) = -GRAD_Scaled(i1:i2)
                     IF (my_rank == 0) WRITE (*, '(A,I4,2X,A,ES12.4)') '  restart→steepest (no-precond) IA=', IA, 'scale=', pscaler
                  END IF
               END DO
            ELSE
               IMAP = 0
               DO IA = 1, NPAR
                  IF (INVP(IA) == 1) THEN
                     IMAP = IMAP + 1
                     i1 = (IMAP - 1)*NBLOCK + 1
                     i2 = IMAP*NBLOCK
                     DO ik = i1, i2
                        IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                           precond_grad(ik) = GRAD_Scaled(ik)
                        ELSE
                           precond_grad(ik) = HESSDI_scaled(ik)*GRAD_Scaled(ik)
                        END IF
                     END DO
                     precond_grad_norm = MAX(DNRM2(NBLOCK, precond_grad(i1), 1), EPS_DEN)
                     pscaler = 1/precond_grad_norm !MAX(GRAD_Scaled_Norm(IA), EPS_DEN)
                     p_k(i1:i2) = -pscaler*precond_grad(i1:i2)
                     IF (my_rank == 0) WRITE (*, '(A,I4,2X,A,ES12.4,2X,A,ES12.4)') &
                        '  restart→steepest (precond)  IA=', IA, '||H_precond*g||_blk=', precond_grad_norm, 'scale=', pscaler
                  END IF
               END DO
            END IF

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
         BF_a(1:m_LBFGS) = 0.0D0
         BF_p(1:m_LBFGS) = 0.0D0
         BF_q(1:NM) = GRAD_Scaled(1:NM)

         ! backward loop  (unchanged indices; guarded inverses)  [B,D]
         DO im = 1, m_LBFGS
            jcol = mml - iiter + im
            BP = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_s_hist(1:NM, jcol))   ! y·s
            BA = DOT_PRODUCT(BF_s_hist(1:NM, jcol), BF_q(1:NM))                ! s·q

            IF (.NOT. ieee_is_finite(BP) .OR. DABS(BP) < EPS_DEN) THEN
               BF_p(im) = 0.0D0
               BF_a(im) = 0.0D0
            ELSE
               BF_p(im) = 1.0D0/BP
               IF (.NOT. ieee_is_finite(BA)) BA = 0.0D0
               BF_a(im) = BF_p(im)*BA
            END IF

            IF (BF_a(im) .NE. 0.0D0) THEN
               BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, jcol)
            END IF
         END DO

         ! H0  [C,E]
         BPK = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))       ! y·s
         IF (USE_SCALAR_H0) THEN
            yy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))  ! y·y
            yy = MAX(yy, EPS_DEN)
            IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0D0
            gamma = BPK/yy
            IF (.NOT. ieee_is_finite(gamma) .OR. gamma <= 0.0D0) gamma = 1.0D0
            BF_r(1:NM) = gamma*BF_q(1:NM)
         ELSE
            denom = 0.0D0
            DO ik = 1, NM
               IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) CYCLE
               denom = denom + BF_grad_res(ik, col)*(HESSDI_scaled(ik)*BF_grad_res(ik, col))
            END DO
            denom = MAX(denom, EPS_DEN)
            IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0D0
            tau = BPK/denom
            IF (.NOT. ieee_is_finite(tau) .OR. tau <= 0.0D0) tau = 1.0D0
            DO ik = 1, NM
               IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                  BF_r(ik) = BF_q(ik)      ! identity fallback
               ELSE
                  BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
               END IF
            END DO
         END IF

         ! forward loop (guarded)  [B,D]
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
         IF (.NOT. ieee_is_finite(dot_g_p) .OR. .NOT. ieee_is_finite(pnorm)) THEN   ! [E]
            IF (my_rank == 0) WRITE (*, *) 'LBFGS produced non-finite step → fallback steepest.'
            ! fallback to preconditioned steepest
            IF (USE_SCALAR_H0) THEN
               p_k(1:NM) = -GRAD_Scaled(1:NM)
            ELSE
               DO ik = 1, NM
                  IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                     p_k(ik) = -GRAD_Scaled(ik)
                  ELSE
                     p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
                  END IF
               END DO
            END IF
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
            theta = 6.0D0*(FCOST0 - FCOST_prev) + 3.0D0*DOT_PRODUCT(GRAD_Scaled(1:NM) + grad_prev(1:NM), svec(1:NM))
            IF (s2 > EPS_DEN) THEN                                       ! [B]
               yhat(1:NM) = ytmp(1:NM) + (theta/s2)*svec(1:NM)
            ELSE
               yhat(1:NM) = ytmp(1:NM)
            END IF
            ys_raw = DOT_PRODUCT(ytmp, svec)
            ys_hat = DOT_PRODUCT(yhat, svec)
            yy_damp = DOT_PRODUCT(yhat, yhat)
            IF (ys_hat <= EPS_CURV*MAX(s2, 1.0D0)) THEN
               yhat(1:NM) = yhat(1:NM) + ((EPS_CURV*MAX(s2, 1.0D0) - ys_hat)/MAX(s2, 1.0D0))*svec(1:NM)
               ys_hat = DOT_PRODUCT(yhat, svec)
            END IF
            BF_grad_res(1:NM, 1) = yhat(1:NM)
            IF (my_rank == 0) THEN
               WRITE (*, '(A)') 'LBFGS[damp]: cost-based y update'
               WRITE (*, '(2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
                  'ΔJ=', (FCOST0 - FCOST_prev), 'theta=', theta, 'y·s (raw)=', ys_raw, 'y·s (damped)=', ys_hat
            END IF
         END IF

         ! curvature gate on newest  [D]
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_s_hist(1:NM, 1))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
         sTy = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_grad_res(1:NM, 1))
         sTy_min = C_CURV*DSQRT(MAX(sTs, 0D0)*MAX(yTy, 0D0))
         restart_now = (.NOT. ieee_is_finite(sTy)) .OR. (sTy <= sTy_min)

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (full): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0D0     ! [A]
            BF_s_hist(1:NM, 1:mml) = 0.0D0

            IF (USE_SCALAR_H0) THEN
               IMAP = 0
               DO IA = 1, NPAR
                  IF (INVP(IA) == 1) THEN
                     IMAP = IMAP + 1
                     i1 = (IMAP - 1)*NBLOCK + 1
                     i2 = IMAP*NBLOCK
                     pscaler = 1.0D0/MAX(DABS(GRAD_Scaled_Norm(IA)), EPS_DEN)      ! [B,E]
                     p_k(i1:i2) = -GRAD_Scaled(i1:i2)*pscaler
                     IF (my_rank == 0) WRITE (*, '(A,I4,2X,A,ES12.4)') '  restart→steepest (no-precond) IA=', IA, 'scale=', pscaler
                  END IF
               END DO
            ELSE
               IMAP = 0
               DO IA = 1, NPAR
                  IF (INVP(IA) == 1) THEN
                     IMAP = IMAP + 1
                     i1 = (IMAP - 1)*NBLOCK + 1
                     i2 = IMAP*NBLOCK
                     DO ik = i1, i2
                        IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                           precond_grad(ik) = GRAD_Scaled(ik)
                        ELSE
                           precond_grad(ik) = HESSDI_scaled(ik)*GRAD_Scaled(ik)
                        END IF
                     END DO
                     precond_grad_norm = MAX(DNRM2(NBLOCK, precond_grad(i1), 1), EPS_DEN)
                     pscaler = MAX(GRAD_Scaled_Norm(IA), EPS_DEN)/precond_grad_norm
                     p_k(i1:i2) = -pscaler*precond_grad(i1:i2)
                     IF (my_rank == 0) WRITE (*, '(A,I4,2X,A,ES12.4,2X,A,ES12.4)') &
                        '  restart→steepest (precond)  IA=', IA, '||H_precond*g||_blk=', precond_grad_norm, 'scale=', pscaler
                  END IF
               END DO
            END IF

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
         BF_a(1:m_LBFGS) = 0.0D0
         BF_p(1:m_LBFGS) = 0.0D0
         BF_q(1:NM) = GRAD_Scaled(1:NM)

         ! backward loop (full)  [B,D]
         DO im = 1, m_LBFGS
            col = im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))   ! y·s
            BP = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_q(1:NM))               ! s·q

            IF (.NOT. ieee_is_finite(BB) .OR. DABS(BB) < EPS_DEN) THEN
               BF_p(im) = 0.0D0
               BF_a(im) = 0.0D0
            ELSE
               BF_p(im) = 1.0D0/BB
               IF (.NOT. ieee_is_finite(BP)) BP = 0.0D0
               BF_a(im) = BF_p(im)*BP
            END IF

            IF (BF_a(im) .NE. 0.0D0) THEN
               BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, col)
            END IF
         END DO

         ! H0 from newest (col=1) [C,E]
         BPK = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_s_hist(1:NM, 1))     ! y·s
         IF (USE_SCALAR_H0) THEN
            yy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
            yy = MAX(yy, EPS_DEN)
            IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0D0
            gamma = BPK/yy
            IF (.NOT. ieee_is_finite(gamma) .OR. gamma <= 0.0D0) gamma = 1.0D0
            BF_r(1:NM) = gamma*BF_q(1:NM)
         ELSE
            denom = 0.0D0
            DO ik = 1, NM
               IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) CYCLE
               denom = denom + BF_grad_res(ik, 1)*(HESSDI_scaled(ik)*BF_grad_res(ik, 1))
            END DO
            denom = MAX(denom, EPS_DEN)
            IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0D0
            tau = BPK/denom
            IF (.NOT. ieee_is_finite(tau) .OR. tau <= 0.0D0) tau = 1.0D0
            DO ik = 1, NM
               IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                  BF_r(ik) = BF_q(ik)
               ELSE
                  BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
               END IF
            END DO
         END IF

         ! forward loop (full)  [B,D]
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
         IF (.NOT. ieee_is_finite(dot_g_p) .OR. .NOT. ieee_is_finite(pnorm)) THEN   ! [E]
            IF (my_rank == 0) WRITE (*, *) 'LBFGS produced non-finite step → fallback steepest.'
            IF (USE_SCALAR_H0) THEN
               p_k(1:NM) = -GRAD_Scaled(1:NM)
            ELSE
               DO ik = 1, NM
                  IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0D0) THEN
                     p_k(ik) = -GRAD_Scaled(ik)
                  ELSE
                     p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
                  END IF
               END DO
            END IF
         END IF
         RETURN
      END IF

   END SUBROUTINE lbfgs_step_loop

   SUBROUTINE lbfgs_basic(iter, NM, mml, INVP, NPAR, NBLOCK, &
                          GRAD_scaled, GRAD_scaled_norm, HESSDI_scaled, &
                          p_k, grad_prev, BF_grad_res, BF_s_hist, &
                          m_coarse, m_coarse_prev, FCOST0, USE_PRECOND, my_rank)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! --------- dummy args ---------
      INTEGER, INTENT(IN)              :: iter, NM, mml, my_rank, NPAR, NBLOCK
      INTEGER, INTENT(IN)              :: INVP(:)
      REAL(dp), INTENT(IN)             :: GRAD_scaled(:)
      REAL(dp), INTENT(IN)             :: GRAD_scaled_norm(:)
      REAL(dp), INTENT(IN)             :: HESSDI_scaled(:)
      REAL(dp), INTENT(INOUT)          :: p_k(:)
      REAL(dp), INTENT(IN)             :: grad_prev(:)
      REAL(dp), INTENT(INOUT)          :: BF_grad_res(:, :)
      REAL(dp), INTENT(INOUT)          :: BF_s_hist(:, :)
      REAL(dp), INTENT(INOUT)          :: m_coarse(:), m_coarse_prev(:)
      REAL(dp), INTENT(IN)             :: FCOST0
      LOGICAL, INTENT(IN)              :: USE_PRECOND

      ! --------- locals ----------
      REAL(dp) :: BP = 0.0_dp
      REAL(dp) :: BA = 0.0_dp
      REAL(dp) :: BPK = 0.0_dp
      REAL(dp) :: BB = 0.0_dp
      REAL(dp) :: BF_b = 0.0_dp
      INTEGER  :: iiter = 0
      INTEGER  :: ik = 0
      INTEGER  :: m_LBFGS = 0
      INTEGER  :: im = 0
      INTEGER  :: jcol = 0
      INTEGER  :: IA = 0
      INTEGER  :: i1 = 0
      INTEGER  :: i2 = 0
      INTEGER  :: col = 0

      REAL(dp) :: dot_g_p = 0.0_dp
      REAL(dp) :: gnorm = 0.0_dp
      REAL(dp) :: pnorm = 0.0_dp
      REAL(dp) :: minH = 0.0_dp
      REAL(dp) :: maxH = 0.0_dp
      REAL(dp) :: posMaxH = 0.0_dp
      REAL(dp) :: gHg = 0.0_dp

      REAL(dp) :: BF_a(mml)
      REAL(dp) :: BF_p(mml)
      REAL(dp) :: BF_q(NM)
      REAL(dp) :: BF_r(NM)
      REAL(dp) :: pscaler = 0.0_dp

      REAL(dp) :: sTy = 0.0_dp
      REAL(dp) :: sTs = 0.0_dp
      REAL(dp) :: yTy = 0.0_dp
      REAL(dp) :: sTy_min = 0.0_dp
      REAL(dp) :: denom = 0.0_dp
      REAL(dp) :: tau = 0.0_dp
      REAL(dp) :: yy = 0.0_dp
      REAL(dp) :: gamma = 1.0_dp

      LOGICAL  :: restart_now = .FALSE.

      REAL(dp), PARAMETER :: ZERO = 0.0_dp
      REAL(dp), PARAMETER :: C_CURV = 1.0e-8_dp
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

      ! --------- LBFGS history logging ----------
      INTEGER, SAVE :: lbfgs_unit = -1
      LOGICAL, SAVE :: lbfgs_ready = .FALSE.
      INTEGER        :: ios
      CHARACTER(len=200) :: iomsg

      ! -------- init local arrays ----------
      BF_a(:) = 0.0_dp
      BF_p(:) = 0.0_dp
      BF_q(:) = 0.0_dp
      BF_r(:) = 0.0_dp

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

      ! quick exit if no active parameters (protect against uninitialised blocks)
      im = COUNT(INVP == 1)
      IF (im == 0 .OR. NM <= 0) THEN
         p_k(1:NM) = 0.0_dp
         RETURN
      END IF

      ! -------- diagnostics ----------
      gnorm = DNRM2(NM, GRAD_scaled, 1)
      p_k(1:NM) = 0.0_dp
      IF (USE_PRECOND) THEN
         minH = MINVAL(HESSDI_scaled)
         maxH = MAXVAL(HESSDI_scaled)
         posMaxH = MAXVAL(HESSDI_scaled, MASK=(HESSDI_scaled > 0.0_dp))
         gHg = DOT_PRODUCT(GRAD_scaled, HESSDI_scaled*GRAD_scaled)
      ELSE
         minH = 1.0_dp
         maxH = 1.0_dp
         posMaxH = 1.0_dp
         gHg = 0.0_dp
      END IF

      ! ----- open LBFGS history log (rank 0 only) -----
      IF (my_rank == 0 .AND. .NOT. lbfgs_ready) THEN
         OPEN (NEWUNIT=lbfgs_unit, FILE='out_lbfgs_hist..txt', STATUS='UNKNOWN', &
               POSITION='APPEND', ACTION='WRITE', IOSTAT=ios, IOMSG=iomsg)
         IF (ios == 0) THEN
            WRITE (lbfgs_unit, '(A)') '# iter  m  ||g||        ||p||        g·p          '// &
               '||s||       ||y||       s·y'
            FLUSH (lbfgs_unit)
            lbfgs_ready = .TRUE.
         ELSE
            WRITE (*, *) 'LBFGS: cannot open out_lbfgs_hist..txt: ', TRIM(iomsg)
         END IF
      END IF

      IF (my_rank == 0) THEN
         WRITE (*, *)
         ! WRITE (*, '(A,I6,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
         !    'LBFGS: iter=', iter, '||g||=', gnorm, &
         !    'min(H)=', minH, 'max(H)=', maxH, 'max+(H)=', posMaxH
         IF (USE_PRECOND) WRITE (*, '(A,ES12.4)') 'LBFGS: g^T (H ⊙ g) =', gHg
      END IF

      iiter = iter - 1

      ! ==========================
      ! Iter 0: steepest direction
      ! ==========================
      IF (iiter == 0) THEN

         ! IF (.NOT. USE_PRECOND) THEN
         !    im = 0
         !    DO IA = 1, NPAR
         !       IF (INVP(IA) == 1) THEN
         !          im = im + 1
         !          i1 = (im - 1)*NBLOCK + 1
         !          i2 = im*NBLOCK
         !           pscaler = 1.0_dp/MAX(TINYD, DABS(GRAD_scaled_norm(IA)))
         !          p_k(i1:i2) = -GRAD_scaled(i1:i2)*pscaler
         !       END IF
         !    END DO
         ! ELSE
         !    DO ik = 1, NM
         !       IF (HESSDI_scaled(ik) <= 0.0_dp) THEN
         !          p_k(ik) = -GRAD_scaled(ik)
         !       ELSE
         !          p_k(ik) = -HESSDI_scaled(ik)*GRAD_scaled(ik)
         !       END IF
         !    END DO
         !    ! Normalize full search direction in preconditioned case; avoid stale block indices.
         !    ! NOTE: For multi-parameter problems this global norm may not match block-wise scaling;
         !    ! revisit if parameter-wise norms are required.
         !    pscaler = 1.0_dp /MAX(DNRM2(NM, p_k, 1), TINYD)
         !    p_k(1:NM) = pscaler*p_k(1:NM)
         ! END IF
         IF (.NOT. USE_PRECOND) THEN
            im = 0
            DO IA = 1, NPAR
               IF (INVP(IA) == 1) THEN
                  im = im + 1
                  i1 = (im - 1)*NBLOCK + 1
                  i2 = im*NBLOCK

                  pscaler = 1.0_dp/MAX(TINYD, DABS(GRAD_scaled_norm(IA)))
                  p_k(i1:i2) = -GRAD_scaled(i1:i2)*pscaler
               END IF
            END DO

         ELSE
            DO ik = 1, NM
               IF (HESSDI_scaled(ik) <= 0.0_dp) THEN
                  p_k(ik) = -GRAD_scaled(ik)
               ELSE
                  p_k(ik) = -HESSDI_scaled(ik)*GRAD_scaled(ik)
               END IF
            END DO

            ! Normalize search direction PER active parameter-block (size NBLOCK each).
            ! This avoids a single global norm being dominated by one block.
            im = 0
            DO IA = 1, NPAR
               IF (INVP(IA) == 1) THEN
                  im = im + 1
                  i1 = (im - 1)*NBLOCK + 1
                  i2 = im*NBLOCK

                  pscaler = 1.0_dp/MAX(DNRM2(NBLOCK, p_k(i1:i2), 1), TINYD)
                  p_k(i1:i2) = pscaler*p_k(i1:i2)
               END IF
            END DO
         END IF

         dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)

         IF (dot_g_p >= ZERO) THEN
            WRITE (*, *) 'LBFGS ERROR: iter=0 non-descent direction, g·p =', dot_g_p
            STOP
         END IF

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'iter0: min(p)=', MINVAL(p_k), 'max(p)=', MAXVAL(p_k), '||p||=', pnorm
         END IF

         IF (my_rank == 0 .AND. lbfgs_ready) THEN
            WRITE (lbfgs_unit, '(I5,1X,I2,1P,1X,7E13.5)') iter, 0, &
               gnorm, pnorm, dot_g_p, 0.0_dp, 0.0_dp, 0.0_dp
            FLUSH (lbfgs_unit)
         END IF

         RETURN
      END IF

      ! ==========================
      ! iiter > 0: L-BFGS two-loop
      ! ==========================
      m_LBFGS = MIN(iiter, mml)

      ! ----------- Warm-up: 1 <= iiter <= mml -----------
      IF (iiter > 0 .AND. iiter <= mml) THEN

         col = mml + 1 - iiter

         BF_grad_res(1:NM, col) = GRAD_scaled(1:NM) - grad_prev(1:NM)
         BF_s_hist(1:NM, col) = m_coarse(1:NM) - m_coarse_prev(1:NM)

         sTy = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_grad_res(1:NM, col))
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_s_hist(1:NM, col))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
         sTy_min = C_CURV*SQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
         restart_now = ((sTy <= sTy_min) .OR. .NOT. (sTy == sTy))

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (warm): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

         IF (sTy == 0.0_dp) THEN
            restart_now = .TRUE.
         END IF

1000     CONTINUE
         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0_dp
            BF_s_hist(1:NM, 1:mml) = 0.0_dp

            ! steepest according to USE_PRECOND
            IF (.NOT. USE_PRECOND) THEN
               im = 0
               DO IA = 1, NPAR
                  IF (INVP(IA) == 1) THEN
                     im = im + 1
                     i1 = (im - 1)*NBLOCK + 1
                     i2 = im*NBLOCK
                     pscaler = 1.0_dp/MAX(TINYD, DABS(GRAD_scaled_norm(IA)))
                     p_k(i1:i2) = -GRAD_scaled(i1:i2)*pscaler
                  END IF
               END DO
            ELSE
               DO ik = 1, NM
                  IF (HESSDI_scaled(ik) <= 0.0_dp) THEN
                     p_k(ik) = -GRAD_scaled(ik)
                  ELSE
                     p_k(ik) = -HESSDI_scaled(ik)*GRAD_scaled(ik)
                  END IF
               END DO
            END IF

            dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
            pnorm = DNRM2(NM, p_k, 1)

            IF (dot_g_p >= ZERO) THEN
               WRITE (*, *) 'LBFGS ERROR: warm restart produced non-descent direction, g·p =', dot_g_p
               STOP
            END IF

            IF (my_rank == 0) THEN
               WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm
            END IF

            IF (my_rank == 0 .AND. lbfgs_ready) THEN
               WRITE (lbfgs_unit, '(I5,1X,I2,1P,1X,7E13.5)') iter, 0, &
                  gnorm, pnorm, dot_g_p, 0.0_dp, 0.0_dp, 0.0_dp
               FLUSH (lbfgs_unit)
            END IF

            RETURN
         END IF

         BF_a(1:m_LBFGS) = 0.0_dp
         BF_p(1:m_LBFGS) = 0.0_dp
         BF_q(1:NM) = GRAD_scaled(1:NM)

         DO im = 1, m_LBFGS
            jcol = mml - iiter + im
            BP = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_s_hist(1:NM, jcol))
            BA = DOT_PRODUCT(BF_s_hist(1:NM, jcol), BF_q(1:NM))
            IF (BP == 0.0_dp) THEN
               restart_now = .TRUE.
               EXIT
            END IF
            BF_p(im) = 1.0_dp/BP
            BF_a(im) = BF_p(im)*BA
            BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, jcol)
         END DO
         IF (restart_now) THEN
            GOTO 1000
         END IF

         BPK = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))

         IF (.NOT. USE_PRECOND) THEN
            yy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
            IF (yy > 0.0_dp) THEN
               gamma = BPK/yy
            ELSE
               gamma = 1.0_dp
            END IF
            BF_r(1:NM) = gamma*BF_q(1:NM)
         ELSE
            denom = 0.0_dp
            DO ik = 1, NM
               denom = denom + BF_grad_res(ik, col)*(HESSDI_scaled(ik)*BF_grad_res(ik, col))
            END DO
            IF (denom > 0.0_dp) THEN
               tau = BPK/denom
            ELSE
               tau = 1.0_dp
            END IF
            DO ik = 1, NM
               BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
            END DO
         END IF

         DO im = m_LBFGS, 1, -1
            jcol = mml - iiter + im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_r(1:NM))
            BF_b = BF_p(im)*BB
            BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, jcol)*(BF_a(im) - BF_b)
         END DO

         p_k(1:NM) = -BF_r(1:NM)
         dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
         pnorm = DNRM2(NM, p_k, 1)
         p_k = p_k/pnorm  ! normalize step direction
         IF (dot_g_p >= ZERO) THEN
            WRITE (*, *) 'LBFGS ERROR: warm two-loop produced non-descent direction, g·p =', dot_g_p
            STOP
         END IF

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step: g·p =', dot_g_p, '||p||=', pnorm
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step: min(p)=', MINVAL(p_k), 'max(p)=', MAXVAL(p_k)
         END IF

         IF (my_rank == 0 .AND. lbfgs_ready) THEN
            WRITE (lbfgs_unit, '(I5,1X,I2,1P,1X,7E13.5)') iter, m_LBFGS, &
               gnorm, pnorm, dot_g_p, SQRT(MAX(sTs, 0.0_dp)), SQRT(MAX(yTy, 0.0_dp)), sTy
            FLUSH (lbfgs_unit)
         END IF

         RETURN
      END IF

      ! ----------- Full memory: iiter > mml -----------
      IF (iiter > mml) THEN

         DO im = mml, 2, -1
            BF_grad_res(1:NM, im) = BF_grad_res(1:NM, im - 1)
            BF_s_hist(1:NM, im) = BF_s_hist(1:NM, im - 1)
         END DO
         BF_grad_res(1:NM, 1) = GRAD_scaled(1:NM) - grad_prev(1:NM)
         BF_s_hist(1:NM, 1) = m_coarse(1:NM) - m_coarse_prev(1:NM)

         sTy = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_grad_res(1:NM, 1))
         sTs = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_s_hist(1:NM, 1))
         yTy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
         sTy_min = C_CURV*SQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
         restart_now = ((sTy <= sTy_min) .OR. .NOT. (sTy == sTy))

         IF (my_rank == 0) THEN
            WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
               'curvature (full): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
            IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
         END IF

2000     CONTINUE
         IF (restart_now) THEN
            BF_grad_res(1:NM, 1:mml) = 0.0_dp
            BF_s_hist(1:NM, 1:mml) = 0.0_dp

            IF (.NOT. USE_PRECOND) THEN
               im = 0
               DO IA = 1, NPAR
                  IF (INVP(IA) == 1) THEN
                     im = im + 1
                     i1 = (im - 1)*NBLOCK + 1
                     i2 = im*NBLOCK
                     pscaler = 1.0_dp/MAX(TINYD, DABS(GRAD_scaled_norm(IA)))
                     p_k(i1:i2) = -GRAD_scaled(i1:i2)*pscaler
                  END IF
               END DO
            ELSE
               DO ik = 1, NM
                  IF (HESSDI_scaled(ik) <= 0.0_dp) THEN
                     p_k(ik) = -GRAD_scaled(ik)
                  ELSE
                     p_k(ik) = -HESSDI_scaled(ik)*GRAD_scaled(ik)
                  END IF
               END DO
            END IF

            dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
            pnorm = DNRM2(NM, p_k, 1)
            p_k = p_k/pnorm  ! normalize step direction
            IF (dot_g_p >= ZERO) THEN
               WRITE (*, *) 'LBFGS ERROR: full restart produced non-descent direction, g·p =', dot_g_p
               STOP
            END IF

            IF (my_rank == 0) THEN
               WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm
            END IF

            IF (my_rank == 0 .AND. lbfgs_ready) THEN
               WRITE (lbfgs_unit, '(I5,1X,I2,1P,1X,7E13.5)') iter, 0, &
                  gnorm, pnorm, dot_g_p, 0.0_dp, 0.0_dp, 0.0_dp
               FLUSH (lbfgs_unit)
            END IF

            RETURN
         END IF

         m_LBFGS = mml
         BF_a(1:m_LBFGS) = 0.0_dp
         BF_p(1:m_LBFGS) = 0.0_dp
         BF_q(1:NM) = GRAD_scaled(1:NM)

         DO im = 1, m_LBFGS
            col = im
            BB = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))
            BP = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_q(1:NM))
            IF (BB == 0.0_dp) THEN
               restart_now = .TRUE.
               EXIT
            END IF
            BF_p(im) = 1.0_dp/BB
            BF_a(im) = BF_p(im)*BP
            BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, col)
         END DO
         IF (restart_now) THEN
            GOTO 2000
         END IF

         BPK = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_s_hist(1:NM, 1))

         IF (.NOT. USE_PRECOND) THEN
            yy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
            IF (yy > 0.0_dp) THEN
               gamma = BPK/yy
            ELSE
               gamma = 1.0_dp
            END IF
            BF_r(1:NM) = gamma*BF_q(1:NM)
         ELSE
            denom = 0.0_dp
            DO ik = 1, NM
               denom = denom + BF_grad_res(ik, 1)*(HESSDI_scaled(ik)*BF_grad_res(ik, 1))
            END DO
            IF (denom > 0.0_dp) THEN
               tau = BPK/denom
            ELSE
               tau = 1.0_dp
            END IF
            DO ik = 1, NM
               BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
            END DO
         END IF

         DO im = m_LBFGS, 1, -1
            col = im
            BP = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_r(1:NM))
            BF_b = BF_p(im)*BP
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
            WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'LBFGS step(full): min(p)=', MINVAL(p_k), 'max(p)=', MAXVAL(p_k)
         END IF

         IF (my_rank == 0 .AND. lbfgs_ready) THEN
            WRITE (lbfgs_unit, '(I5,1X,I2,1P,1X,7E13.5)') iter, m_LBFGS, &
               gnorm, pnorm, dot_g_p, SQRT(MAX(sTs, 0.0_dp)), SQRT(MAX(yTy, 0.0_dp)), sTy
            FLUSH (lbfgs_unit)
         END IF

         RETURN
      END IF

   END SUBROUTINE lbfgs_basic

   SUBROUTINE steep_descent(iter, NM, NPAR, NBLOCK, INVP, &
                            GRAD_scaled, GRAD_scaled_norm, HESSDI_scaled, &
                            USE_PRECOND, NORMALIZE_PRECOND, &
                            my_rank, label, &
                            p_k)

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      !------------- dummy args -------------
      INTEGER, INTENT(IN)  :: iter, NM, NPAR, NBLOCK
      INTEGER, INTENT(IN)  :: INVP(:)
      REAL(dp), INTENT(IN)  :: GRAD_scaled(:)
      REAL(dp), INTENT(IN)  :: GRAD_scaled_norm(:)
      REAL(dp), INTENT(IN)  :: HESSDI_scaled(:)
      LOGICAL, INTENT(IN)  :: USE_PRECOND
      LOGICAL, INTENT(IN)  :: NORMALIZE_PRECOND   ! .TRUE. for iter=0, .FALSE. for restarts
      INTEGER, INTENT(IN)  :: my_rank
      CHARACTER(LEN=*), INTENT(IN)  :: label               ! e.g. 'iter0', 'warm restart', ...
      REAL(dp), INTENT(INOUT) :: p_k(:)            ! search direction (output)

      !------------- locals -------------
      INTEGER  :: IA, im, i1, i2, ik
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
      ! Build steepest-descent direction according to USE_PRECOND
      !===========================================================

      p_k(1:NM) = 0.0_dp
      IF (.NOT. USE_PRECOND) THEN
         ! Block-wise scaled steepest descent
         im = 0
         DO IA = 1, NPAR
            IF (INVP(IA) == 1) THEN
               im = im + 1
               i1 = (im - 1)*NBLOCK + 1
               i2 = im*NBLOCK

               pscaler = 1.0_dp/MAX(TINYD, ABS(GRAD_scaled_norm(IA)))
               p_k(i1:i2) = -GRAD_scaled(i1:i2)*pscaler
            END IF
         END DO

      ELSE
         ! Diagonal-preconditioned steepest descent
         DO ik = 1, NM
            IF (HESSDI_scaled(ik) <= 0.0_dp) THEN
               p_k(ik) = -GRAD_scaled(ik)
            ELSE
               p_k(ik) = -HESSDI_scaled(ik)*GRAD_scaled(ik)
            END IF
         END DO

         IF (NORMALIZE_PRECOND) THEN
            ! Match iter=0 behavior in lbfgs_basic: normalize ||p|| = 1
            pscaler = 1.0_dp/MAX(DNRM2(NM, p_k, 1), TINYD)
            p_k(1:NM) = pscaler*p_k(1:NM)
         END IF
      END IF

      !===========================
      ! Descent check and logging
      !===========================
      dot_g_p = DOT_PRODUCT(GRAD_scaled, p_k)
      pnorm = DNRM2(NM, p_k, 1)

      IF (dot_g_p >= ZERO) THEN
         WRITE (*, *) 'LBFGS ERROR: ', TRIM(label), &
            ' produced non-descent direction, g·p =', dot_g_p
         STOP
      END IF

      IF (my_rank == 0) THEN
         WRITE (*, '(A,1X,A,1X,A,I6)') 'LBFGS steepest step [', &
            TRIM(label), '], iter=', iter
         WRITE (*, '(A,ES12.4,2X,A,ES12.4)') '  g·p =', dot_g_p, '||p|| =', pnorm
         WRITE (*, '(A,ES12.4,2X,A,ES12.4)') '  min(p)=', MINVAL(p_k), &
            'max(p)=', MAXVAL(p_k)
      END IF

   END SUBROUTINE steep_descent

   ! SUBROUTINE lbfgs_step_loop(iter, NM, mml, INVP, NPAR, NBLOCK, &
   !                            GRAD_Scaled, GRAD_scaled_norm, HESSDI_scaled, &
   !                            p_k, grad_prev, BF_grad_res, &
   !                            BF_s_hist, m_coarse, m_coarse_prev, FCOST0, FCOST_prev, &
   !                            USE_LBFGS_TYPE, USE_PRECOND, my_rank)

   !    USE, INTRINSIC :: ieee_arithmetic            ! [E] finite checks
   !    IMPLICIT NONE

   !    ! --------- dummy args ---------
   !    INTEGER, INTENT(IN)    :: iter, NM, mml, my_rank, INVP(:), NPAR, NBLOCK
   !    REAL(dp), INTENT(IN)   :: GRAD_Scaled(:)                 ! scaled gradient g_k
   !    REAL(dp), INTENT(IN)   :: GRAD_Scaled_Norm(:)            ! per-parameter norms
   !    REAL(dp), INTENT(IN)   :: HESSDI_scaled(:)                ! diagonal preconditioner in scaled space
   !    REAL(dp), INTENT(INOUT):: p_k(:)
   !    REAL(dp), INTENT(IN)   :: m_coarse(:), m_coarse_prev(:)
   !    REAL(dp), INTENT(IN)   :: grad_prev(:)                    ! previous scaled gradient g_{k-1}
   !    REAL(dp), INTENT(INOUT):: BF_grad_res(:, :)              ! Y history (y = g_k - g_{k-1})
   !    REAL(dp), INTENT(INOUT):: BF_s_hist(:, :)                 ! S history (s = x_k - x_{k-1})
   !    REAL(dp), INTENT(IN)   :: FCOST0
   !    REAL(dp), INTENT(IN)   :: FCOST_prev
   !    INTEGER, INTENT(IN)    :: USE_LBFGS_TYPE
   !    LOGICAL, INTENT(IN)    :: USE_PRECOND

   !    ! --------- locals ----------
   !    REAL(dp) :: BP, BA, BPK, BB, BF_b
   !    INTEGER  :: iiter, ik, m_LBFGS, im, jcol, IA, i1, i2, col, IMAP
   !    REAL(dp) :: dot_g_p, gnorm, pnorm, gHg, minH, maxH, posMaxH
   !    REAL(dp) :: BF_a(mml), BF_p(mml)
   !    REAL(dp) :: BF_q(NM), BF_r(NM), pscaler, precond_grad_norm, precond_grad(NM)
   !    REAL(dp) :: sTy, sTs, yTy, sTy_min
   !    REAL(dp) :: denom, tau, yy, gamma
   !    REAL(dp) :: pnorm_blk, gblk_norm
   !    LOGICAL  :: restart_now
   !    REAL(dp), PARAMETER :: ZERO = 0.0_dp
   !    REAL(dp), PARAMETER :: C_CURV = 1.0e-10_dp
   !    REAL(dp), PARAMETER :: TINYD = 1.0e-37_dp

   !    ! --------- cost-damping / trust-region ----------
   !    REAL(dp) :: s2, theta, ys_raw, ys_hat, yy_damp
   !    REAL(dp) :: ytmp(NM), svec(NM), yhat(NM)
   !    REAL(dp), PARAMETER :: EPS_CURV = 1.0e-14_dp
   !    REAL(dp), PARAMETER :: DELTA_TR = 1.0_dp
   !    REAL(dp) :: pnorm_tr, scale_tr

   !    REAL(dp), PARAMETER :: EPS_THETA = 1.0e-30_dp

   !    restart_now = .FALSE.
   !    BF_a(1:mml) = 0.0_dp
   !    BF_p(1:mml) = 0.0_dp
   !    BF_q(1:NM) = 0.0_dp
   !    BF_r(1:NM) = 0.0_dp
   !    precond_grad(1:NM) = 0.0_dp
   !    ytmp(1:NM) = 0.0_dp
   !    svec(1:NM) = 0.0_dp
   !    yhat(1:NM) = 0.0_dp
   !    ! ========= diagnostics =========
   !    gnorm = DNRM2(NM, GRAD_Scaled, 1)
   !    IF (USE_PRECOND) THEN

   !       gHg = 0.0_dp
   !       DO ik = 1, NM
   !          IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) CYCLE
   !          gHg = gHg + GRAD_Scaled(ik)*(HESSDI_scaled(ik)*GRAD_Scaled(ik))
   !       END DO
   !       minH = MINVAL(HESSDI_scaled)
   !       maxH = MAXVAL(HESSDI_scaled)
   !       posMaxH = MAXVAL(HESSDI_scaled, MASK=(HESSDI_scaled > 0.0_dp))
   !    ELSE
   !       gHg = DOT_PRODUCT(GRAD_Scaled, GRAD_Scaled)
   !       minH = 1.0_dp
   !       maxH = 1.0_dp
   !       posMaxH = 1.0_dp
   !    END IF

   !    IF (my_rank == 0) THEN
   !       WRITE (*, *) '  '
   !       WRITE (*, '(A,I6,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
   !          'LBFGS: iter=', iter, '||g||=', gnorm, 'min(H)=', minH, 'max(H)=', maxH, 'max+(H)=', posMaxH
   !       WRITE (*, '(A,ES12.4)') 'LBFGS: g^T (H ⊙ g) =', gHg
   !    END IF

   !    iiter = iter - 1

   !    ! ==========================
   !    ! Iter 0: steepest direction
   !    ! ==========================
   !    IF (iiter == 0) THEN
   !       IMAP = 0
   !       DO IA = 1, NPAR
   !          IF (INVP(IA) == 1) THEN
   !             IMAP = IMAP + 1
   !             i1 = (IMAP - 1)*NBLOCK + 1
   !             i2 = IMAP*NBLOCK

   !             IF (USE_PRECOND) THEN
   !                ! [PC] Pure preconditioned steepest: p = -HESSDI_scaled ⊙ g (no block normalisation)
   !                DO ik = i1, i2
   !                   IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                      p_k(ik) = -GRAD_Scaled(ik)
   !                   ELSE
   !                      p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
   !                   END IF
   !                END DO
   !                pnorm_blk = SQRT(SUM(p_k(i1:i2)*p_k(i1:i2)))
   !                IF (my_rank == 0) WRITE (*, '(A,I4,2X,2(A,ES12.4,2X))') &
   !                   '  steepest (precond)  IA=', IA, '||g_blk||=', GRAD_Scaled_Norm(IA), &
   !                   '||p_blk||=', pnorm_blk
   !             ELSE
   !                ! No-precond: scalar steepest using block norm
   !                gblk_norm = MAX(DABS(GRAD_Scaled_Norm(IA)), TINYD)
   !                pscaler = 1.0_dp/gblk_norm
   !                p_k(i1:i2) = -GRAD_Scaled(i1:i2)*pscaler
   !                pnorm_blk = SQRT(SUM(p_k(i1:i2)*p_k(i1:i2)))
   !                IF (my_rank == 0) WRITE (*, '(A,I4,2X,3(A,ES12.4,2X))') &
   !                   '  steepest (no-precond) IA=', IA, '||g_blk||=', gblk_norm, &
   !                   'scale=', pscaler, '||p_blk||=', pnorm_blk
   !             END IF
   !          END IF
   !       END DO

   !       dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
   !       pnorm = DNRM2(NM, p_k, 1)
   !       IF (my_rank == 0) WRITE (*, '(A,ES12.4,2X,A,ES12.4)') 'iter0: min(p)=', MINVAL(p_k), 'max(p)=', MAXVAL(p_k)
   !       IF (.NOT. ieee_is_finite(dot_g_p) .OR. dot_g_p >= ZERO) THEN
   !          IF (my_rank == 0) WRITE (*, *) 'WARNING: Non-descent/NaN at iter=0 → stop.'
   !          STOP
   !       END IF

   !       IF (USE_LBFGS_TYPE >= 2) THEN
   !          pnorm_tr = DNRM2(NM, p_k, 1)
   !          IF (pnorm_tr > DELTA_TR) THEN
   !             scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
   !             p_k(1:NM) = p_k(1:NM)*scale_tr
   !             IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(iter0): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
   !          END IF
   !       END IF

   !       RETURN
   !    END IF

   !    ! ==========================
   !    ! iiter > 0: L-BFGS two-loop
   !    ! ==========================
   !    m_LBFGS = MIN(iiter, mml)

   !    ! -------- Warm-up window: 1 <= iiter <= mml --------
   !    IF (iiter > 0 .AND. iiter <= mml) THEN
   !       col = mml + 1 - iiter

   !       BF_grad_res(1:NM, col) = GRAD_Scaled(1:NM) - grad_prev(1:NM)     ! y_k
   !       BF_s_hist(1:NM, col) = m_coarse(1:NM) - m_coarse_prev(1:NM)   ! s_k

   !       ! ---- cost-damped y (type>=1) ----
   !       IF (USE_LBFGS_TYPE >= 1) THEN
   !          ytmp(1:NM) = BF_grad_res(1:NM, col)
   !          svec(1:NM) = BF_s_hist(1:NM, col)
   !          s2 = DOT_PRODUCT(svec, svec)
   !          theta = 6.0_dp*(FCOST0 - FCOST_prev) + &
   !                  3.0_dp*DOT_PRODUCT(grad_prev(1:NM) + GRAD_Scaled(1:NM), svec(1:NM))

   !          IF (s2 > TINYD) THEN
   !             yhat(1:NM) = ytmp(1:NM)
   !             IF (DABS(theta) > EPS_THETA) yhat(1:NM) = yhat(1:NM) + (theta/s2)*svec(1:NM)
   !          ELSE
   !             yhat(1:NM) = ytmp(1:NM)
   !          END IF

   !          ys_raw = DOT_PRODUCT(ytmp, svec)
   !          ys_hat = DOT_PRODUCT(yhat, svec)
   !          yy_damp = DOT_PRODUCT(yhat, yhat)

   !          IF (ys_hat <= EPS_CURV*MAX(s2, 1.0_dp)) THEN
   !             yhat(1:NM) = yhat(1:NM) + ((EPS_CURV*MAX(s2, 1.0_dp) - ys_hat)/MAX(s2, 1.0_dp))*svec(1:NM)
   !             ys_hat = DOT_PRODUCT(yhat, svec)
   !          END IF

   !          BF_grad_res(1:NM, col) = yhat(1:NM)

   !          IF (my_rank == 0) THEN
   !             WRITE (*, '(A)') 'LBFGS[damp]: cost-based y update'
   !             WRITE (*, '(2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
   !                'ΔJ=', (FCOST0 - FCOST_prev), 'theta=', theta, 'y·s (raw)=', ys_raw, 'y·s (damped)=', ys_hat
   !          END IF
   !       END IF

   !       ! curvature gate
   !       sTs = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_s_hist(1:NM, col))
   !       yTy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
   !       sTy = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_grad_res(1:NM, col))
   !       IF (.NOT. ieee_is_finite(sTs) .OR. .NOT. ieee_is_finite(yTy) .OR. &
   !           .NOT. ieee_is_finite(sTy)) THEN
   !          restart_now = .TRUE.
   !       ELSE
   !          sTy_min = C_CURV*DSQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
   !          restart_now = (sTy <= sTy_min)
   !       END IF
   !       IF (my_rank == 0) THEN
   !          WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
   !             'curvature (warm): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
   !          IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
   !       END IF

   !       IF (restart_now) THEN
   !          BF_grad_res(1:NM, 1:mml) = 0.0_dp
   !          BF_s_hist(1:NM, 1:mml) = 0.0_dp

   !          ! fallback steepest
   !          IMAP = 0
   !          DO IA = 1, NPAR
   !             IF (INVP(IA) == 1) THEN
   !                IMAP = IMAP + 1
   !                i1 = (IMAP - 1)*NBLOCK + 1
   !                i2 = IMAP*NBLOCK

   !                IF (USE_PRECOND) THEN
   !                   ! [PC] restart → pure preconditioned steepest
   !                   DO ik = i1, i2
   !                      IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                         p_k(ik) = -GRAD_Scaled(ik)
   !                      ELSE
   !                         p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
   !                      END IF
   !                   END DO
   !                   pnorm_blk = SQRT(SUM(p_k(i1:i2)*p_k(i1:i2)))
   !                   IF (my_rank == 0) WRITE (*, '(A,I4,2X,2(A,ES12.4,2X))') &
   !                      '  restart→steepest (precond)  IA=', IA, '||g_blk||=', GRAD_Scaled_Norm(IA), &
   !                      '||p_blk||=', pnorm_blk
   !                ELSE
   !                   gblk_norm = MAX(DABS(GRAD_Scaled_Norm(IA)), TINYD)
   !                   pscaler = 1.0_dp/gblk_norm
   !                   p_k(i1:i2) = -GRAD_Scaled(i1:i2)*pscaler
   !                   pnorm_blk = SQRT(SUM(p_k(i1:i2)*p_k(i1:i2)))
   !                   IF (my_rank == 0) WRITE (*, '(A,I4,2X,3(A,ES12.4,2X))') &
   !                      '  restart→steepest (no-precond) IA=', IA, '||g_blk||=', gblk_norm, &
   !                      'scale=', pscaler, '||p_blk||=', pnorm_blk
   !                END IF
   !             END IF
   !          END DO

   !          dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
   !          pnorm = DNRM2(NM, p_k, 1)
   !          IF (my_rank == 0) WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm

   !          IF (USE_LBFGS_TYPE >= 2) THEN
   !             pnorm_tr = DNRM2(NM, p_k, 1)
   !             IF (pnorm_tr > DELTA_TR) THEN
   !                scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
   !                p_k(1:NM) = p_k(1:NM)*scale_tr
   !                IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(restart): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
   !             END IF
   !          END IF

   !          RETURN
   !       END IF

   !       ! two-loop prep
   !       BF_a(1:m_LBFGS) = 0.0_dp
   !       BF_p(1:m_LBFGS) = 0.0_dp
   !       BF_q(1:NM) = GRAD_Scaled(1:NM)

   !       ! backward loop
   !       DO im = 1, m_LBFGS
   !          jcol = mml - iiter + im
   !          IF (jcol < 1 .OR. jcol > SIZE(BF_grad_res, 2) .OR. jcol > SIZE(BF_s_hist, 2)) THEN
   !             WRITE (*, *) 'LBFGS_STEP_LOOP: bad jcol on rank', my_rank, &
   !                ' jcol=', jcol, ' size2(BF_grad_res)=', SIZE(BF_grad_res, 2), &
   !                ' size2(BF_s_hist)=', SIZE(BF_s_hist, 2)
   !             STOP 'LBFGS: jcol out of bounds'
   !          END IF
   !          BP = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_s_hist(1:NM, jcol))   ! y·s
   !          BA = DOT_PRODUCT(BF_s_hist(1:NM, jcol), BF_q(1:NM))                ! s·q

   !          IF (.NOT. ieee_is_finite(BP) .OR. DABS(BP) < TINYD) THEN
   !             BF_p(im) = 0.0_dp
   !             BF_a(im) = 0.0_dp
   !          ELSE
   !             BF_p(im) = 1.0_dp/BP
   !             IF (.NOT. ieee_is_finite(BA)) BA = 0.0_dp
   !             BF_a(im) = BF_p(im)*BA
   !          END IF

   !          IF (BF_a(im) /= 0.0_dp) THEN
   !             BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, jcol)
   !          END IF
   !       END DO

   !       ! H0
   !       ! H0 - RESPECT USE_PRECOND!
   !       BPK = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))     ! y·s

   !       IF (USE_PRECOND) THEN
   !          ! Diagonal preconditioning: H₀ = τ·diag(H)
   !          denom = 0.0_dp
   !          DO ik = 1, NM
   !             IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) CYCLE
   !             denom = denom + BF_grad_res(ik, col)*(HESSDI_scaled(ik)*BF_grad_res(ik, col))
   !          END DO
   !          denom = MAX(denom, TINYD)
   !          IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
   !          tau = BPK/denom
   !          IF (.NOT. ieee_is_finite(tau) .OR. tau <= 0.0_dp) tau = 1.0_dp

   !          DO ik = 1, NM
   !             IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                BF_r(ik) = BF_q(ik)
   !             ELSE
   !                BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
   !             END IF
   !          END DO
   !       ELSE
   !          ! Scalar initial Hessian: H₀ = γ·I
   !          yy = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_grad_res(1:NM, col))
   !          yy = MAX(yy, TINYD)
   !          IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
   !          gamma = BPK/yy
   !          IF (.NOT. ieee_is_finite(gamma) .OR. gamma <= 0.0_dp) gamma = 1.0_dp
   !          BF_r(1:NM) = gamma*BF_q(1:NM)
   !       END IF

   !       ! forward loop
   !       DO im = m_LBFGS, 1, -1
   !          jcol = mml - iiter + im
   !          BB = DOT_PRODUCT(BF_grad_res(1:NM, jcol), BF_r(1:NM))   ! y·r
   !          IF (.NOT. ieee_is_finite(BF_p(im))) CYCLE
   !          IF (.NOT. ieee_is_finite(BB)) CYCLE
   !          BF_b = BF_p(im)*BB
   !          IF (.NOT. ieee_is_finite(BF_b)) CYCLE
   !          BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, jcol)*(BF_a(im) - BF_b)
   !       END DO

   !       p_k(1:NM) = -BF_r(1:NM)

   !       IF (USE_LBFGS_TYPE >= 2) THEN
   !          pnorm_tr = DNRM2(NM, p_k, 1)
   !          IF (pnorm_tr > DELTA_TR) THEN
   !             scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
   !             p_k(1:NM) = p_k(1:NM)*scale_tr
   !             IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(warm): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
   !          END IF
   !       END IF

   !       dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
   !       pnorm = DNRM2(NM, p_k, 1)
   !       IF (.NOT. ieee_is_finite(dot_g_p) .OR. .NOT. ieee_is_finite(pnorm)) THEN
   !          IF (my_rank == 0) WRITE (*, *) 'LBFGS produced non-finite step → fallback steepest.'
   !          IF (USE_PRECOND) THEN
   !             DO ik = 1, NM
   !                IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                   p_k(ik) = -GRAD_Scaled(ik)
   !                ELSE
   !                   p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
   !                END IF
   !             END DO
   !          ELSE
   !             p_k(1:NM) = -GRAD_Scaled(1:NM)
   !          END IF
   !       END IF
   !       RETURN
   !    END IF

   !    ! -------- Full memory: iiter > mml --------
   !    IF (iiter > mml) THEN
   !       DO im = mml, 2, -1
   !          BF_grad_res(1:NM, im) = BF_grad_res(1:NM, im - 1)
   !          BF_s_hist(1:NM, im) = BF_s_hist(1:NM, im - 1)
   !       END DO
   !       BF_grad_res(1:NM, 1) = GRAD_Scaled(1:NM) - grad_prev(1:NM)
   !       BF_s_hist(1:NM, 1) = m_coarse(1:NM) - m_coarse_prev(1:NM)

   !       IF (USE_LBFGS_TYPE >= 1) THEN
   !          ytmp(1:NM) = BF_grad_res(1:NM, 1)
   !          svec(1:NM) = BF_s_hist(1:NM, 1)
   !          s2 = DOT_PRODUCT(svec, svec)
   !          theta = 6.0_dp*(FCOST0 - FCOST_prev) + &
   !                  3.0_dp*DOT_PRODUCT(GRAD_Scaled(1:NM) + grad_prev(1:NM), svec(1:NM))
   !          IF (s2 > TINYD) THEN
   !             yhat(1:NM) = ytmp(1:NM) + (theta/s2)*svec(1:NM)
   !          ELSE
   !             yhat(1:NM) = ytmp(1:NM)
   !          END IF
   !          ys_raw = DOT_PRODUCT(ytmp, svec)
   !          ys_hat = DOT_PRODUCT(yhat, svec)
   !          yy_damp = DOT_PRODUCT(yhat, yhat)
   !          IF (ys_hat <= EPS_CURV*MAX(s2, 1.0_dp)) THEN
   !             yhat(1:NM) = yhat(1:NM) + ((EPS_CURV*MAX(s2, 1.0_dp) - ys_hat)/MAX(s2, 1.0_dp))*svec(1:NM)
   !             ys_hat = DOT_PRODUCT(yhat, svec)
   !          END IF
   !          BF_grad_res(1:NM, 1) = yhat(1:NM)
   !          IF (my_rank == 0) THEN
   !             WRITE (*, '(A)') 'LBFGS[damp]: cost-based y update'
   !             WRITE (*, '(2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
   !                'ΔJ=', (FCOST0 - FCOST_prev), 'theta=', theta, 'y·s (raw)=', ys_raw, 'y·s (damped)=', ys_hat
   !          END IF
   !       END IF

   !       ! curvature gate on newest
   !       sTs = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_s_hist(1:NM, 1))
   !       yTy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
   !       sTy = DOT_PRODUCT(BF_s_hist(1:NM, 1), BF_grad_res(1:NM, 1))
   !       sTy_min = C_CURV*DSQRT(MAX(sTs, 0.0_dp)*MAX(yTy, 0.0_dp))
   !       restart_now = (.NOT. ieee_is_finite(sTy)) .OR. (sTy <= sTy_min)

   !       IF (my_rank == 0) THEN
   !          WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
   !             'curvature (full): s·y=', sTy, '||s||^2=', sTs, '||y||^2=', yTy, 'sTy_min=', sTy_min
   !          IF (restart_now) WRITE (*, *) '  -> restart: insufficient/invalid curvature'
   !       END IF

   !       IF (restart_now) THEN
   !          BF_grad_res(1:NM, 1:mml) = 0.0_dp
   !          BF_s_hist(1:NM, 1:mml) = 0.0_dp

   !          IMAP = 0
   !          DO IA = 1, NPAR
   !             IF (INVP(IA) == 1) THEN
   !                IMAP = IMAP + 1
   !                i1 = (IMAP - 1)*NBLOCK + 1
   !                i2 = IMAP*NBLOCK

   !                IF (USE_PRECOND) THEN
   !                   ! [PC] restart→steepest (precond) full-memory
   !                   DO ik = i1, i2
   !                      IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                         p_k(ik) = -GRAD_Scaled(ik)
   !                      ELSE
   !                         p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
   !                      END IF
   !                   END DO
   !                   pnorm_blk = SQRT(SUM(p_k(i1:i2)*p_k(i1:i2)))
   !                   IF (my_rank == 0) WRITE (*, '(A,I4,2X,2(A,ES12.4,2X))') &
   !                      '  restart→steepest (precond)  IA=', IA, '||g_blk||=', GRAD_Scaled_Norm(IA), &
   !                      '||p_blk||=', pnorm_blk
   !                ELSE
   !                   gblk_norm = MAX(DABS(GRAD_Scaled_Norm(IA)), TINYD)
   !                   pscaler = 1.0_dp/gblk_norm
   !                   p_k(i1:i2) = -GRAD_Scaled(i1:i2)*pscaler
   !                   pnorm_blk = SQRT(SUM(p_k(i1:i2)*p_k(i1:i2)))
   !                   IF (my_rank == 0) WRITE (*, '(A,I4,2X,3(A,ES12.4,2X))') &
   !                      '  restart→steepest (no-precond)  IA=', IA, '||g_blk||=', gblk_norm, &
   !                      'scale=', pscaler, '||p_blk||=', pnorm_blk
   !                END IF
   !             END IF
   !          END DO

   !          dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
   !          pnorm = DNRM2(NM, p_k, 1)
   !          IF (my_rank == 0) WRITE (*, '("LBFGS restart: g·p=",ES12.4,"  ||p||=",ES12.4)') dot_g_p, pnorm

   !          IF (USE_LBFGS_TYPE >= 2) THEN
   !             pnorm_tr = DNRM2(NM, p_k, 1)
   !             IF (pnorm_tr > DELTA_TR) THEN
   !                scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
   !                p_k(1:NM) = p_k(1:NM)*scale_tr
   !                IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(restart): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
   !             END IF
   !          END IF

   !          RETURN
   !       END IF

   !       m_LBFGS = mml
   !       BF_a(1:m_LBFGS) = 0.0_dp
   !       BF_p(1:m_LBFGS) = 0.0_dp
   !       BF_q(1:NM) = GRAD_Scaled(1:NM)

   !       ! backward loop (full)
   !       DO im = 1, m_LBFGS
   !          col = im
   !          BB = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_s_hist(1:NM, col))   ! y·s
   !          BP = DOT_PRODUCT(BF_s_hist(1:NM, col), BF_q(1:NM))               ! s·q

   !          IF (.NOT. ieee_is_finite(BB) .OR. DABS(BB) < TINYD) THEN
   !             BF_p(im) = 0.0_dp
   !             BF_a(im) = 0.0_dp
   !          ELSE
   !             BF_p(im) = 1.0_dp/BB
   !             IF (.NOT. ieee_is_finite(BP)) BP = 0.0_dp
   !             BF_a(im) = BF_p(im)*BP
   !          END IF

   !          IF (BF_a(im) /= 0.0_dp) THEN
   !             BF_q(1:NM) = BF_q(1:NM) - BF_a(im)*BF_grad_res(1:NM, col)
   !          END IF
   !       END DO

   !       ! H0 from newest (col=1)
   !       BPK = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_s_hist(1:NM, 1))     ! y·s
   !       IF (USE_PRECOND) THEN
   !          denom = 0.0_dp
   !          DO ik = 1, NM
   !             IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) CYCLE
   !             denom = denom + BF_grad_res(ik, 1)*(HESSDI_scaled(ik)*BF_grad_res(ik, 1))
   !          END DO
   !          denom = MAX(denom, TINYD)
   !          IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
   !          tau = BPK/denom
   !          IF (.NOT. ieee_is_finite(tau) .OR. tau <= 0.0_dp) tau = 1.0_dp
   !          DO ik = 1, NM
   !             IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                BF_r(ik) = BF_q(ik)
   !             ELSE
   !                BF_r(ik) = tau*HESSDI_scaled(ik)*BF_q(ik)
   !             END IF
   !          END DO
   !       ELSE
   !          yy = DOT_PRODUCT(BF_grad_res(1:NM, 1), BF_grad_res(1:NM, 1))
   !          yy = MAX(yy, TINYD)
   !          IF (.NOT. ieee_is_finite(BPK)) BPK = 0.0_dp
   !          gamma = BPK/yy
   !          IF (.NOT. ieee_is_finite(gamma) .OR. gamma <= 0.0_dp) gamma = 1.0_dp
   !          BF_r(1:NM) = gamma*BF_q(1:NM)
   !       END IF

   !       ! forward loop (full)
   !       DO im = m_LBFGS, 1, -1
   !          col = im
   !          BP = DOT_PRODUCT(BF_grad_res(1:NM, col), BF_r(1:NM))       ! y·r
   !          IF (.NOT. ieee_is_finite(BF_p(im))) CYCLE
   !          IF (.NOT. ieee_is_finite(BP)) CYCLE
   !          BF_b = BF_p(im)*BP
   !          IF (.NOT. ieee_is_finite(BF_b)) CYCLE
   !          BF_r(1:NM) = BF_r(1:NM) + BF_s_hist(1:NM, col)*(BF_a(im) - BF_b)
   !       END DO

   !       p_k(1:NM) = -BF_r(1:NM)

   !       IF (USE_LBFGS_TYPE >= 2) THEN
   !          pnorm_tr = DNRM2(NM, p_k, 1)
   !          IF (pnorm_tr > DELTA_TR) THEN
   !             scale_tr = DELTA_TR/MAX(pnorm_tr, TINYD)
   !             p_k(1:NM) = p_k(1:NM)*scale_tr
   !             IF (my_rank == 0) WRITE (*, '(A,ES12.4,A,ES12.4)') 'TR clamp(full): ||p|| ', pnorm_tr, ' -> ', DELTA_TR
   !          END IF
   !       END IF

   !       dot_g_p = DOT_PRODUCT(GRAD_Scaled, p_k)
   !       pnorm = DNRM2(NM, p_k, 1)
   !       IF (.NOT. ieee_is_finite(dot_g_p) .OR. .NOT. ieee_is_finite(pnorm)) THEN
   !          IF (my_rank == 0) WRITE (*, *) 'LBFGS produced non-finite step → fallback steepest.'
   !          IF (USE_PRECOND) THEN
   !             DO ik = 1, NM
   !                IF (.NOT. ieee_is_finite(HESSDI_scaled(ik)) .OR. HESSDI_scaled(ik) <= 0.0_dp) THEN
   !                   p_k(ik) = -GRAD_Scaled(ik)
   !                ELSE
   !                   p_k(ik) = -HESSDI_scaled(ik)*GRAD_Scaled(ik)
   !                END IF
   !             END DO
   !          ELSE
   !             p_k(1:NM) = -GRAD_Scaled(1:NM)
   !          END IF
   !       END IF
   !       RETURN
   !    END IF

   ! END SUBROUTINE lbfgs_step_loop

   SUBROUTINE GRID2D_UPDATE_MODEL(p_kf, m_fine, CR, CI, INVP, IANISO, ALPHA, &
                                  NPT, SCALER, NPAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NPT, NPAR, IANISO
      INTEGER, INTENT(IN)    :: INVP(:)
      REAL(dp), INTENT(IN)    :: ALPHA
      REAL(dp), INTENT(IN)    :: SCALER(:), p_kf(:)
      REAL(dp), INTENT(INOUT) :: CR(:, :), CI(:, :), m_fine(:)

      INTEGER :: I, IM, ps, pe

      IM = 0
      DO I = 1, NPAR
         IF (INVP(I) /= 1) CYCLE
         IM = IM + 1
         ps = (IM - 1)*NPT + 1
         pe = IM*NPT

         ! update grid baseline stack
         m_fine(ps:pe) = m_fine(ps:pe) + ALPHA*p_kf(ps:pe)

         ! map back to grid arrays for next forward step
         IF (I <= IANISO) THEN
            CR(I, 1:NPT) = m_fine(ps:pe)*SCALER(I)
         ELSE
            CI(I - (IANISO - 1), 1:NPT) = m_fine(ps:pe)*SCALER(I)
         END IF
      END DO
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

   PURE real(dp) FUNCTION alpha0_from_hess_precond(g, h, a_min, a_max) RESULT(a0)
      IMPLICIT real(dp) (A - H, O - Z)
      real(dp), INTENT(IN) :: g(:), h(:)
      real(dp), INTENT(IN), OPTIONAL :: a_min, a_max
      real(dp) :: amin, amax, g2, denom
      INTEGER :: i, n, npos

      amin = MERGE(a_min, 1.0D-3, PRESENT(a_min))   ! <- pick a sensible floor; 1e-3 works well
      amax = MERGE(a_max, 1.0D+3, PRESENT(a_max))   ! <- avoid absurdly large α0

      n = SIZE(g)
      IF (SIZE(h) /= n) THEN
         a0 = 1.0D0
         RETURN
      END IF

      g2 = 0.0D0
      denom = 0.0D0
      npos = 0
      DO i = 1, n
         g2 = g2 + g(i)*g(i)
         IF (h(i) > 0.0D0) THEN
            denom = denom + (g(i)*g(i))/h(i)
            npos = npos + 1
         END IF
      END DO

      ! Not enough positive entries? Fall back.
      IF (npos < MAX(10, n/100)) THEN
         a0 = 1.0D0
         RETURN
      END IF

      IF (denom > 0.0D0 .AND. g2 > 0.0D0) THEN
         a0 = g2/denom
         a0 = MAX(amin, MIN(amax, a0))
      ELSE
         a0 = 1.0D0
      END IF
   END FUNCTION alpha0_from_hess_precond

   SUBROUTINE lbfgs_project_direction(NPAR, INVP, NBLOCK, NPT, NNX, NNZ, NX, NZ, NTO, XTO, ZTO, &
                                      XBC, ZBC, p_k, p_kf, ITER, my_rank, unit_block, unit_point, &
                                      DEBUG_OUTPUT, pnorm_out)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPAR, NBLOCK, NPT, NNX, NNZ, NX, NZ, NTO, ITER, my_rank
      INTEGER, INTENT(IN) :: INVP(:)
      REAL(dp), INTENT(IN) :: XBC(:), ZBC(:), XTO(:), ZTO(:)
      REAL(dp), INTENT(INOUT) :: p_k(:)
      REAL(dp), INTENT(INOUT) :: p_kf(:)
      INTEGER, OPTIONAL, INTENT(IN) :: unit_block, unit_point
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT
      REAL(dp), OPTIONAL, INTENT(OUT) :: pnorm_out

      INTEGER :: II, IM, cs, ce, ps, pe, u_blk, u_pt, nmm
      LOGICAL :: dbg
      REAL(dp) :: pnorm_local

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      u_blk = -1; IF (PRESENT(unit_block)) u_blk = unit_block
      u_pt = -1; IF (PRESENT(unit_point)) u_pt = unit_point

      p_kf = 0.0_dp
      IM = 0
      DO II = 1, NPAR
         IF (INVP(II) == 1) THEN
            IM = IM + 1
            cs = (IM - 1)*NBLOCK + 1
            ce = IM*NBLOCK
            ps = (IM - 1)*NPT + 1
            pe = IM*NPT

            CALL GRID_TRANSFORM_BLOCK2NPT(NNX, NNZ, NX - 1, NZ - 1, XBC, ZBC, &
                                          p_k(cs:ce), NTO, XTO, ZTO, p_kf(ps:pe))
         END IF
      END DO

      nmm = SIZE(p_kf)
      pnorm_local = DNRM2(nmm, p_kf, 1)
      IF (PRESENT(pnorm_out)) pnorm_out = pnorm_local
      if (my_rank==0) WRITE (*, '(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') 'min(p_kf)=', MINVAL(p_kf), 'max(p_kf)=', MAXVAL(p_kf), '||p_kf||=', pnorm_local
      ! p_kf = p_kf/MAX(pnorm_local, 1.0D0)

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
!  Author: you
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
