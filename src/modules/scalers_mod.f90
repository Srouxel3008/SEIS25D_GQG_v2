module scalers_mod
   use constant_mod
   use iso_fortran_env, only: dp => real64
   use grid_mod
   use gridtype_mod

   IMPLICIT NONE

CONTAINS

   SUBROUTINE ComputeModelScales(CR, CI, NPAR, IANISO, ITHOM, param_scale_raw, my_rank, &
                                 USE_RMS, DEBUG_OUTPUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPAR, IANISO, ITHOM, my_rank
      REAL(dp), INTENT(IN)  :: CR(:, :), CI(:, :)
      REAL(dp), INTENT(OUT) :: param_scale_raw(:)          ! length NPAR
      LOGICAL, OPTIONAL, INTENT(IN) :: USE_RMS, DEBUG_OUTPUT

      LOGICAL :: use_rms_local, dbg
      INTEGER :: IA, irow, npts, ncr, nci
      REAL(dp), PARAMETER :: FLOOR_REAL = 1.0_dp
      REAL(dp), PARAMETER :: FLOOR_Q = 1.0_dp
      REAL(dp) :: one_rad
      LOGICAL :: has_theta

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

! Default: RMS scaling
      use_rms_local = .TRUE.
      IF (PRESENT(USE_RMS)) use_rms_local = USE_RMS
      npts = SIZE(CR, 2)
      ncr = SIZE(CR, 1)
      nci = SIZE(CI, 1)

      one_rad = 1.0_dp
      has_theta = (ITHOM == 0 .AND. IANISO == 7)

      IF (SIZE(param_scale_raw) < NPAR) THEN
         IF (my_rank == 0) WRITE (*, *) 'Compute RAWParameterScales: param_scale_raw too small.'
         STOP
      END IF

      DO IA = 1, NPAR
         IF (IA <= IANISO) THEN
            irow = IA
            IF (irow >= 1 .AND. irow <= ncr) THEN
               param_scale_raw(IA) = stat_val(ABS(CR(irow, 1:npts)), use_rms_local)
            ELSE
               param_scale_raw(IA) = FLOOR_REAL
            END IF
            ! real parameters: keep a safe floor
            param_scale_raw(IA) = MAX(param_scale_raw(IA), FLOOR_REAL)
         ELSE
            irow = IA - (IANISO - 1)
            IF (irow >= 1 .AND. irow <= nci) THEN
               param_scale_raw(IA) = stat_val(ABS(CI(irow, 1:npts)), use_rms_local)
            ELSE
               param_scale_raw(IA) = FLOOR_Q
            END IF
            ! imag/Q-like: also keep floor
            param_scale_raw(IA) = MAX(param_scale_raw(IA), FLOOR_Q)
         END IF
      END DO

      ! Force theta scale if applicable (stiffness IANISO=7) or Thomsen theta
      IF (ITHOM == 0 .AND. has_theta) THEN
         IF (NPAR >= 7) param_scale_raw(7) = param_scale_raw(7)
      ELSEIF (ITHOM == 1) THEN
         IF (NPAR >= 7) param_scale_raw(7) = param_scale_raw(7)
      END IF

      IF (dbg .AND. my_rank == 0) THEN
         DO IA = 1, NPAR
            WRITE (*, '(A,I3,A,1PE12.5)') ' param_scale_raw(IA=', IA, ') = ', param_scale_raw(IA)
         END DO
      END IF

   CONTAINS
      FUNCTION stat_val(vec, use_rms_local) RESULT(s)
         REAL(dp), INTENT(IN) :: vec(:)
         LOGICAL, INTENT(IN)  :: use_rms_local
         REAL(dp) :: s
         IF (SIZE(vec) <= 0) THEN
            s = 0.0_dp
         ELSEIF (use_rms_local) THEN
            s = SQRT(SUM(vec*vec)/REAL(SIZE(vec), dp))
         ELSE
            s = MAXVAL(vec)
         END IF
      END FUNCTION stat_val

   END SUBROUTINE ComputeModelScales

   SUBROUTINE ComputeFrechetScalers(FRECHET, ND, NBLOCK, NM, NPAR, INVP, &
                                    WD_amp, F_SCALE, F_SCALE_L2, my_rank, DEBUG_OUTPUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ND, NBLOCK, NM, NPAR, my_rank
      INTEGER, INTENT(IN) :: INVP(:)
      COMPLEX(dp), INTENT(IN) :: FRECHET(ND, NM)
      REAL(dp), INTENT(IN) :: WD_amp

      REAL(dp), INTENT(OUT) :: F_SCALE(:)      ! 1 / max |FRECHET|
      REAL(dp) :: F_SCALE_L2(:)   ! 1 / RMS(|FRECHET|)

      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      LOGICAL :: dbg
      INTEGER :: IA, IP, NP, iBlk, iDat
      INTEGER :: nsamp
      REAL(dp) :: vmax_f, l2sum_f, rms_f, val_abs
      REAL(dp), PARAMETER :: EPS = 1.0e-300_dp

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

      IF (SIZE(F_SCALE) < NPAR .OR. SIZE(F_SCALE_L2) < NPAR) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeFrechetScalers: output arrays too small.'
         STOP
      END IF
      IF (SIZE(INVP) < NPAR) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeFrechetScalers: INVP size < NPAR.'
         STOP
      END IF

      F_SCALE(:) = 1.0_dp
      F_SCALE_L2(:) = 1.0_dp
      nsamp = MAX(1, ND*NBLOCK)

      IP = 0
      DO IA = 1, NPAR
         IF (INVP(IA) /= 1) CYCLE
         IP = IP + 1
         NP = (IP - 1)*NBLOCK

         vmax_f = 0.0_dp
         l2sum_f = 0.0_dp

         DO iBlk = 1, NBLOCK
            DO iDat = 1, ND
               val_abs = ABS(FRECHET(iDat, NP + iBlk)*Wd_amp)
               vmax_f = MAX(vmax_f, val_abs)
               l2sum_f = l2sum_f + val_abs*val_abs
            END DO
         END DO

         rms_f = SQRT(MAX(l2sum_f/REAL(nsamp, dp), 0.0_dp))

         ! For active parameters, use the scalers to normalize the Frechet derivatives
         F_SCALE(IA) = 1.0_dp/MAX(vmax_f, EPS)
         F_SCALE_L2(IA) = 1.0_dp/MAX(rms_f, EPS)

         IF (dbg .AND. my_rank == 0) THEN
            WRITE (*, '(A,I3,2X,A,1PE12.5,2X,A,1PE12.5,2(2X,A,1PE12.5))') &
               'IA=', IA, 'max|F|=', vmax_f, 'rms|F|=', rms_f, &
               'F_SCALE(Inf)=', F_SCALE(IA), 'F_SCALE_L2=', F_SCALE_L2(IA)
         END IF
      END DO

   END SUBROUTINE ComputeFrechetScalers

   SUBROUTINE ComputeGroupScales(param_scale_raw, NPAR, IANISO, ITHOM, INVP, PAR_SCALE, my_rank, &
                                 USE_GROUP, DEBUG_OUTPUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPAR, IANISO, ITHOM, my_rank
      INTEGER, INTENT(IN) :: INVP(:)          ! length NPAR
      REAL(dp), INTENT(IN)  :: param_scale_raw(:)       ! length NPAR, from ComputePureParameterScales
      REAL(dp), INTENT(OUT) :: PAR_SCALE(:)   ! length NPAR
      LOGICAL, OPTIONAL, INTENT(IN) :: USE_GROUP, DEBUG_OUTPUT

      LOGICAL :: use_grp, dbg
      INTEGER :: IA
      REAL(dp), PARAMETER :: FLOOR_REAL = 1.0_dp
      REAL(dp), PARAMETER :: FLOOR_Q = 1.0_dp
      REAL(dp) :: one_rad
      LOGICAL :: has_theta

      REAL(dp) :: s_rho, s_cij, s_q, s_vel, s_aniso

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT
      use_grp = .TRUE.; IF (PRESENT(USE_GROUP)) use_grp = USE_GROUP

      one_rad = 1.0_dp
      has_theta = (ITHOM == 0 .AND. IANISO == 7)

      IF (SIZE(PAR_SCALE) < NPAR) THEN
         IF (my_rank == 0) WRITE (*, *) 'BuildGroupedParameterScales: PAR_SCALE too small.'
         STOP
      END IF

      ! Default: inactive -> 1 (harmless)
      PAR_SCALE(:) = FLOOR_REAL

      ! -----------------------------------------
      ! No grouping: active IA uses param_scale_raw(IA)
      ! -----------------------------------------
      IF (.NOT. use_grp) THEN
         DO IA = 1, NPAR
            IF (INVP(IA) /= 1) CYCLE
            PAR_SCALE(IA) = MAX(param_scale_raw(IA), FLOOR_REAL)
            IF (IA > IANISO) PAR_SCALE(IA) = MAX(PAR_SCALE(IA), FLOOR_Q)
         END DO
         IF (dbg .AND. my_rank == 0) THEN
            DO IA = 1, NPAR
               WRITE (*, '(A,I3,A,1PE12.5)') 'INDIV PAR_SCALE(IA=', IA, ') = ', PAR_SCALE(IA)
            END DO
         END IF
         RETURN
      END IF

      ! -----------------------------------------
      ! Grouping: use ONLY ACTIVE IA to build buckets
      ! -----------------------------------------
      IF (ITHOM == 0) THEN
         ! stiffness mode
         s_rho = FLOOR_REAL
         s_cij = FLOOR_REAL
         s_q = FLOOR_Q

         ! rho bucket (only if rho is active)
         IF (NPAR >= 1) THEN
            IF (INVP(1) == 1) s_rho = MAX(s_rho, param_scale_raw(1))
         END IF

         ! stiffness bucket: IA=2..min(IANISO,6), active only
         DO IA = 2, MIN(IANISO, 6)
            IF (INVP(IA) /= 1) CYCLE
            s_cij = MAX(s_cij, param_scale_raw(IA))
         END DO

         ! Q bucket: IA>IANISO, active only
         DO IA = IANISO + 1, NPAR
            IF (INVP(IA) /= 1) CYCLE
            s_q = MAX(s_q, param_scale_raw(IA))
         END DO
         s_q = MAX(s_q, FLOOR_Q)

         DO IA = 1, NPAR
            IF (INVP(IA) /= 1) CYCLE
            IF (IA == 1) THEN
               PAR_SCALE(IA) = MAX(s_rho, FLOOR_REAL)
            ELSEIF (IA >= 2 .AND. IA <= MIN(IANISO, 6)) THEN
               PAR_SCALE(IA) = MAX(s_cij, FLOOR_REAL)
            ELSEIF (has_theta .AND. IA == 7) THEN
               PAR_SCALE(IA) = one_rad
            ELSEIF (IA > IANISO) THEN
               PAR_SCALE(IA) = MAX(s_q, FLOOR_Q)
            ELSE
               ! Any other active real parameter outside 2..6:
               ! safest default: treat like stiffness magnitude (NOT like Q)
               PAR_SCALE(IA) = MAX(s_cij, FLOOR_REAL)
            END IF
         END DO

      ELSE
         ! Thomsen mode
         s_rho = FLOOR_REAL
         s_vel = FLOOR_REAL
         s_aniso = 1.0_dp
         s_q = FLOOR_Q

         IF (NPAR >= 1) THEN
            IF (INVP(1) == 1) s_rho = MAX(s_rho, param_scale_raw(1))
         END IF

         ! velocities 2:3, active only
         DO IA = 2, MIN(3, NPAR)
            IF (INVP(IA) /= 1) CYCLE
            s_vel = MAX(s_vel, param_scale_raw(IA))
         END DO

         ! anisotropy 4:6 kept O(1) by design (your new bucket logic)
         s_aniso = 1.0_dp

         ! Q bucket
         DO IA = IANISO + 1, NPAR
            IF (INVP(IA) /= 1) CYCLE
            s_q = MAX(s_q, param_scale_raw(IA))
         END DO
         s_q = MAX(s_q, FLOOR_Q)

         DO IA = 1, NPAR
            IF (INVP(IA) /= 1) CYCLE
            SELECT CASE (IA)
            CASE (1)
               PAR_SCALE(IA) = MAX(s_rho, FLOOR_REAL)
            CASE (2:3)
               PAR_SCALE(IA) = MAX(s_vel, FLOOR_REAL)
            CASE (4:6)
               PAR_SCALE(IA) = s_aniso
            CASE (7)
               PAR_SCALE(IA) = one_rad
            CASE DEFAULT
               IF (IA > IANISO) THEN
                  PAR_SCALE(IA) = MAX(s_q, FLOOR_Q)
               ELSE
                  PAR_SCALE(IA) = MAX(s_vel, FLOOR_REAL)
               END IF
            END SELECT
         END DO
      END IF

      IF (dbg .AND. my_rank == 0) THEN
         DO IA = 1, NPAR
            WRITE (*, '(A,I3,A,1PE12.5,A,I2)') 'GROUP PAR_SCALE(IA=', IA, ') = ', PAR_SCALE(IA), &
               '  inv=', INVP(IA)
         END DO
      END IF

   END SUBROUTINE ComputeGroupScales

   SUBROUTINE ComputeModelScaling(CR, CI, &
                                  NPT, NBLOCK, IG, &
                                  NNX, NNZ, NX, NZ, NORD, &
                                  NTO, XTO, ZTO, &
                                  NPAR, IANISO, INVP, &
                                  PAR_SCALE, GMASK, &
                                  m_fine, m_coarse_prev, &
                                  m_coarse, m_coarse_reg, &
                                  NM, IFQ, norm_init, &
                                  my_rank, DEBUG_OUTPUT)

      USE gridtype_mod, ONLY: InversionGridType
      IMPLICIT NONE

      ! ---- geometry / sizes ----
      INTEGER, INTENT(IN) :: NPT, NBLOCK, NNX, NNZ, NX, NZ, NTO, NORD
      INTEGER, INTENT(IN) :: NPAR, IANISO, my_rank, NM, IFQ
      INTEGER, INTENT(IN) :: INVP(:)
      ! ---- coordinates (coarse/block & topo/GQG) ----
      TYPE(InversionGridType), INTENT(IN)  :: IG
      REAL(dp), INTENT(IN)  :: XTO(:), ZTO(:)   ! size NTO
      ! ---- model (value-space, per-parameter rows) ----
      REAL(dp), INTENT(IN)  :: CR(:, :), CI(:, :)
      REAL(dp), INTENT(IN)  :: PAR_SCALE(:)     ! per-IA scale used to DIVIDE model

      ! ---- outputs ----
      REAL(dp), INTENT(INOUT) :: m_fine(:)   ! stacked per-IA vectors on NPT (NPT × nActive)
      REAL(dp), INTENT(INOUT) :: m_coarse_prev(:)    ! stacked per-IA vectors on NBLOCK (NBLOCK × nActive)
      REAL(dp), INTENT(INOUT) :: m_coarse(:)
      REAL(dp), INTENT(INOUT) :: m_coarse_reg(:), GMASK(:)
      REAL(dp), INTENT(OUT)   :: norm_init

      ! ---- optional ----
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      ! ---- locals ----
      INTEGER :: IA, IM, irow, i1, i2, actN
      INTEGER :: II, cs, ce, ps, pe
      INTEGER :: nprint
      REAL(dp) :: sc_par
      LOGICAL :: dbg
      REAL(dp), PARAMETER :: EPS = 1.0e-50_dp  ! guard for divisions

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      ! -------------------------
      ! Basic size checks
      ! -------------------------
      actN = COUNT(INVP == 1)
      IF (SIZE(m_fine) < actN*NPT) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeModelScaling: m_fine too small. Have ', &
            SIZE(m_fine), ' need ', actN*NPT
         STOP
      END IF
      IF (SIZE(m_coarse_prev) < actN*NBLOCK) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeModelScaling: m_coarse_prev too small. Have ', &
            SIZE(m_coarse_prev), ' need ', actN*NBLOCK
         STOP
      END IF
      IF (SIZE(PAR_SCALE) < NPAR) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeModelScaling: PAR_SCALE too small. Need ', NPAR
         STOP
      END IF

      ! -------------------------
      ! Build per-parameter model on NPT and divide by SCALER
      ! -------------------------
      IM = 0
      DO IA = 1, NPAR
         IF (INVP(IA) /= 1) CYCLE
         IM = IM + 1
         i1 = (IM - 1)*NPT + 1
         i2 = IM*NPT

         IF (IA <= IANISO) THEN
            irow = IA
            m_fine(i1:i2) = CR(irow, 1:NPT)
         ELSE
            irow = IA - (IANISO - 1)
            m_fine(i1:i2) = CI(irow, 1:NPT)
         END IF

         sc_par = MAX(PAR_SCALE(IA), EPS)            ! guard
         m_fine(i1:i2) = m_fine(i1:i2)/sc_par   ! m_s = m / PAR_SCALE

      END DO

! --- new fine→coarse mapping on inversion grid ---
      CALL map_npt2block_all(NPAR, INVP, NPT, NBLOCK, NORD, NNZ, IG%N0_BLOCK, &
                             m_fine, m_coarse_prev)

      IF (dbg .and. my_rank == 0) THEN
         WRITE (*, '(A,1X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,I8)') &
            'm_fine:', &
            'min=', MINVAL(m_fine), &
            'max=', MAXVAL(m_fine), &
            'mean=', SUM(m_fine)/REAL(SIZE(m_fine), dp), &
            '||.||2=', DNRM2(SIZE(m_fine), m_fine, 1), &
            'size=', SIZE(m_fine)

         WRITE (*, '(A,1X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,I8)') &
            'm_coarse_prev:', &
            'min=', MINVAL(m_coarse_prev), &
            'max=', MAXVAL(m_coarse_prev), &
            'mean=', SUM(m_coarse_prev)/REAL(SIZE(m_coarse_prev), dp), &
            '||.||2=', DNRM2(SIZE(m_coarse_prev), m_coarse_prev, 1), &
            'size=', SIZE(m_coarse_prev)

      END IF
      m_coarse_prev = m_coarse_prev*GMASK
      IF (IFQ == 1) m_coarse_reg = m_coarse_prev
      m_coarse(1:NM) = m_coarse_prev(1:NM)
      norm_init = DNRM2(NM, m_coarse_prev(1:NM), 1)

   END SUBROUTINE ComputeModelScaling

   SUBROUTINE ComputeModelScaling2D(CR, CI, &
                                    NPT, NBLOCK, IG, &
                                    NNX, NNZ, NX, NZ, NORD, &
                                    NTO, XTO, ZTO, &
                                    NPAR, IANISO, INVP, &
                                    PAR_SCALE, GMASK, &
                                    m_fine, m_coarse_prev, &
                                    m_coarse, m_coarse_reg, &
                                    NM, IFQ, norm_init, &
                                    my_rank, DEBUG_OUTPUT)

      USE gridtype_mod, ONLY: InversionGridType
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPT, NBLOCK, NNX, NNZ, NX, NZ, NTO, NORD
      INTEGER, INTENT(IN) :: NPAR, IANISO, my_rank, NM, IFQ
      INTEGER, INTENT(IN) :: INVP(:)
      TYPE(InversionGridType), INTENT(IN)  :: IG
      REAL(dp), INTENT(IN)  :: XTO(:), ZTO(:)
      REAL(dp), INTENT(IN)  :: CR(:, :), CI(:, :)
      REAL(dp), INTENT(IN)  :: PAR_SCALE(:)

      REAL(dp), INTENT(INOUT) :: m_fine(:)
      REAL(dp), INTENT(INOUT) :: m_coarse_prev(:)
      REAL(dp), INTENT(INOUT) :: m_coarse(:)
      REAL(dp), INTENT(INOUT) :: m_coarse_reg(:), GMASK(:)
      REAL(dp), INTENT(OUT)   :: norm_init

      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      INTEGER :: IA, IM, irow, i1, i2, actN
      INTEGER :: cs, ce
      REAL(dp) :: sc_par
      LOGICAL :: dbg
      REAL(dp), PARAMETER :: EPS = 1.0e-50_dp

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      actN = COUNT(INVP == 1)
      IF (SIZE(m_fine) < actN*NPT) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeModelScaling2D: m_fine too small. Have ', &
            SIZE(m_fine), ' need ', actN*NPT
         STOP
      END IF
      IF (SIZE(m_coarse_prev) < actN*NBLOCK) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeModelScaling2D: m_coarse_prev too small. Have ', &
            SIZE(m_coarse_prev), ' need ', actN*NBLOCK
         STOP
      END IF
      IF (SIZE(PAR_SCALE) < NPAR) THEN
         IF (my_rank == 0) WRITE (*, *) 'ComputeModelScaling2D: PAR_SCALE too small. Need ', NPAR
         STOP
      END IF

      ! Build the raw fine-grid parameter stack first.
      IM = 0
      DO IA = 1, NPAR
         IF (INVP(IA) /= 1) CYCLE
         IM = IM + 1
         i1 = (IM - 1)*NPT + 1
         i2 = IM*NPT

         IF (IA <= IANISO) THEN
            irow = IA
            m_fine(i1:i2) = CR(irow, 1:NPT)
         ELSE
            irow = IA - (IANISO - 1)
            m_fine(i1:i2) = CI(irow, 1:NPT)
         END IF
      END DO

      ! Build the coarse model from the raw fine-grid stack.
      CALL map_npt2block_all(NPAR, INVP, NPT, NBLOCK, NORD, NNZ, IG%N0_BLOCK, &
                             m_fine, m_coarse_prev)

      ! Apply parameter scaling only after the raw coarse grid has been created.
      IM = 0
      DO IA = 1, NPAR
         IF (INVP(IA) /= 1) CYCLE
         IM = IM + 1
         cs = (IM - 1)*NBLOCK + 1
         ce = IM*NBLOCK

         sc_par = MAX(PAR_SCALE(IA), EPS)
         m_coarse_prev(cs:ce) = m_coarse_prev(cs:ce)/sc_par
      END DO

      IF (dbg .and. my_rank == 0) THEN
         WRITE (*, '(A,1X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,I8)') &
            'm_fine:', &
            'min=', MINVAL(m_fine), &
            'max=', MAXVAL(m_fine), &
            'mean=', SUM(m_fine)/REAL(SIZE(m_fine), dp), &
            '||.||2=', DNRM2(SIZE(m_fine), m_fine, 1), &
            'size=', SIZE(m_fine)

         WRITE (*, '(A,1X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,I8)') &
            'm_coarse_prev:', &
            'min=', MINVAL(m_coarse_prev), &
            'max=', MAXVAL(m_coarse_prev), &
            'mean=', SUM(m_coarse_prev)/REAL(SIZE(m_coarse_prev), dp), &
            '||.||2=', DNRM2(SIZE(m_coarse_prev), m_coarse_prev, 1), &
            'size=', SIZE(m_coarse_prev)
      END IF

      m_coarse_prev = m_coarse_prev*GMASK
      IF (IFQ == 1) m_coarse_reg = m_coarse_prev
      m_coarse(1:NM) = m_coarse_prev(1:NM)
      norm_init = DNRM2(NM, m_coarse_prev(1:NM), 1)

   END SUBROUTINE ComputeModelScaling2D

   SUBROUTINE ComputeDataScalerW(ND, GT0, WD_acq, WD_amp, INV)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ND, INV
      COMPLEX(dp), INTENT(IN) :: GT0(:)
      REAL(dp), INTENT(IN) :: WD_acq(:)
      REAL(dp), INTENT(OUT) :: WD_amp
      REAL(dp) :: sum_wd2, sum_w, data_rms
      INTEGER :: i

      sum_wd2 = 0.0_dp
      sum_w = 0.0_dp

      ! if (INV == 2) then
      !    WD_amp = 1.0D0
      !    RETURN
      ! else

      DO i = 1, ND
         sum_wd2 = sum_wd2 + WD_acq(i)*(REAL(GT0(i), dp)**2 + AIMAG(GT0(i))**2)
         sum_w = sum_w + WD_acq(i)
      END DO

      IF (sum_w > 0.0_dp) THEN
         data_rms = SQRT(MAX(sum_wd2/sum_w, tiny_dp))
      ELSE
         data_rms = 1.0_dp
      END IF

      WD_amp = 1.0_dp/MAX(data_rms, tiny_dp)
      ! end if
   END SUBROUTINE

   SUBROUTINE ComputeAcqWeight(NFQ, ND, NSS, NS, NR, NSV, NRV, XSR, ZSR, &
                               WD_acq, ierr, my_rank, acq_weight_file, DEBUG_OUTPUT)

      IMPLICIT NONE

      ! ---- inputs ----
      INTEGER, INTENT(IN)  :: NFQ, ND, NSS, my_rank
      INTEGER, INTENT(IN)  :: NS(:), NR(:), NSV(:), NRV(:)
      REAL(dp), INTENT(IN)  :: XSR(:), ZSR(:)
      CHARACTER(len=*), INTENT(IN)  :: acq_weight_file
      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT

      ! ---- outputs ----
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: WD_acq(:)
      INTEGER, INTENT(OUT)     :: ierr

      ! ---- locals ----
      INTEGER  :: j, ios, idist, nbin, rc, stat_alloc
      REAL(dp) :: xs, zs, xr, zr, dz, dist
      REAL(dp) :: dx_bin, dz_bin, bin_w, w_geom, w_comp
      REAL(dp), ALLOCATABLE :: weight_table(:)
      REAL(dp) :: c_weight(3)
      LOGICAL  :: verbose, dbg, dump_file
      INTEGER  :: unit_file
      CHARACTER(len=256) :: dbg_name
      REAL(dp) :: sum_w

      ! ------------------------------------------------------------------
      ! init
      ! ------------------------------------------------------------------
      ierr = 0
      verbose = (my_rank == 0)
      dbg = .FALSE.

      IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT .AND. verbose
      dump_file = dbg

      ! basic consistency checks before any allocation
      IF (ND <= 0 .OR. NSS <= 0) THEN
         ierr = 900
         IF (verbose) WRITE (*, '(A,2I8)') '[ComputeAcqWeight] invalid ND/NSS', ND, NSS
         RETURN
      END IF
      IF (ND > SIZE(NS) .OR. ND > SIZE(NR) .OR. ND > SIZE(NSV) .OR. ND > SIZE(NRV)) THEN
         ierr = 901
         IF (verbose) WRITE (*, '(A,4I8)') '[ComputeAcqWeight] ND exceeds array sizes: ', &
            SIZE(NS), SIZE(NR), SIZE(NSV), SIZE(NRV)
         RETURN
      END IF
      IF (NSS > SIZE(XSR) .OR. NSS > SIZE(ZSR)) THEN
         ierr = 902
         IF (verbose) WRITE (*, '(A,2I8)') '[ComputeAcqWeight] NSS exceeds XSR/ZSR sizes: ', SIZE(XSR), SIZE(ZSR)
         RETURN
      END IF

      ! allocate output weights
      ALLOCATE (WD_acq(ND), STAT=stat_alloc)

      IF (stat_alloc /= 0) THEN
         ierr = 1001
         IF (verbose) WRITE (*, '(A,I0)') &
            '[ComputeAcqWeight] ERROR: ALLOCATE(WD_acq) failed, STAT=', stat_alloc
         RETURN
      END IF

      ! --- if no file path given, fallback to 1.0 and exit ---
      IF (LEN_TRIM(acq_weight_file) == 0) THEN
         IF (verbose) WRITE (*, '(A)') '[ComputeAcqWeight] empty acq_weight_file, using uniform weights = 1.0'
         WD_acq(:) = 1.0_dp/REAL(ND, dp)
         RETURN
      END IF

      ! --- open acquisition weight file ---
      OPEN (NEWUNIT=unit_file, FILE=TRIM(acq_weight_file), STATUS='OLD', &
            ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
         ierr = 10
         IF (verbose) WRITE (*, '(A,1X,A,1X,A,I0)') &
            '[ComputeAcqWeight] WARNING: cannot open acq weight file:', &
            TRIM(acq_weight_file), 'IOSTAT=', ios
         WD_acq(:) = 1.0_dp/REAL(ND, dp)
         RETURN
      END IF

      ! --- line 1: dx_bin, dz_bin ---
      READ (unit_file, *, IOSTAT=ios) dx_bin, dz_bin
      IF (ios /= 0) THEN
         ierr = 11
         IF (verbose) WRITE (*, '(A)') '[ComputeAcqWeight] ERROR reading dx_bin, dz_bin'
         WD_acq(:) = 1.0_dp
         CLOSE (unit_file)
         RETURN
      END IF
      ! --- line 2: nbin ---
      READ (unit_file, *, IOSTAT=ios) nbin
      IF (ios /= 0 .OR. nbin <= 0) THEN
         ierr = 12
         IF (verbose) WRITE (*, '(A,I0)') &
            '[build_acq_weight] WARNING: invalid nbin read, IOSTAT=', ios
         CLOSE (unit_file)
         WD_acq(1:ND) = 1.0_dp
         RETURN
      END IF

      ! --- line 3: weight_table(1..nbin) ---
      ALLOCATE (weight_table(nbin), STAT=stat_alloc)
      IF (stat_alloc /= 0) THEN
         ierr = 1002
         IF (verbose) WRITE (*, '(A,I0)') &
            '[build_acq_weight] ERROR: ALLOCATE(weight_table) failed, STAT=', stat_alloc
         CLOSE (unit_file)
         WD_acq(1:ND) = 1.0_dp
         RETURN
      END IF

      READ (unit_file, *, IOSTAT=ios) (weight_table(J), J=1, nbin)
      IF (ios /= 0) THEN
         ierr = 13
         IF (verbose) WRITE (*, '(A,I0)') &
            '[build_acq_weight] WARNING: failed to read weight table, IOSTAT=', ios
         DEALLOCATE (weight_table)
         CLOSE (unit_file)
         WD_acq(1:ND) = 1.0_dp
         RETURN
      END IF

      ! --- line 4: component weights ---
      c_weight = 1.0_dp
      READ (unit_file, *, IOSTAT=ios) c_weight(1), c_weight(2), c_weight(3)
      IF (ios /= 0) THEN
         ierr = 14
         c_weight = 1.0_dp
         IF (verbose) WRITE (*, '(A,I0)') &
            '[build_acq_weight] NOTICE: component weights missing; using defaults.'
      END IF

      CLOSE (unit_file)

      ! --- binning axis (fixed: vertical distance) ---
      bin_w = dx_bin
      IF (bin_w <= 0.0_dp) THEN
         ierr = 15
         IF (verbose) WRITE (*, '(A,ES12.4)') &
            '[build_acq_weight] WARNING: non-positive dz_bin, using uniform weights; dz_bin=', dz_bin
         DEALLOCATE (weight_table)
         WD_acq(1:ND) = 1.0_dp
         RETURN
      END IF

      ! --- optional debug dump file ---
      IF (dump_file) THEN
         WRITE (dbg_name, '(A,I0,A)') 'WD_acq_debug_rank', my_rank, '.txt'
         OPEN (UNIT=99, FILE=TRIM(dbg_name), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
         IF (ios == 0) THEN
            unit_file = 99
            WRITE (unit_file, '(A)') '# J   idist   w_geom    rc   w_comp     WD_acq(J)'
         ELSE
            unit_file = -1
         END IF
      ELSE
         unit_file = -1
      END IF

      ! ------------------------------------------------------------------
      ! BUILD THE ACQUISITION WEIGHT
      ! ------------------------------------------------------------------
      DO J = 1, ND

         ! --- geometry bounds checking ---
         IF (NS(J) < 1 .OR. NS(J) > SIZE(XSR) .OR. NR(J) < 1 .OR. NR(J) > SIZE(XSR)) THEN
            ierr = MAX(ierr, 2001)
            WD_acq(J) = 1.0_dp
            CYCLE
         END IF
         IF (NSV(J) < 1 .OR. NSV(J) > SIZE(NSV)) THEN
            ierr = MAX(ierr, 2002)
            WD_acq(J) = 1.0_dp
            CYCLE
         END IF
         IF (NRV(J) < 1 .OR. NRV(J) > SIZE(NRV)) THEN
            ierr = MAX(ierr, 2003)
            WD_acq(J) = 1.0_dp
            CYCLE
         END IF

         xs = XSR(NS(J)); zs = ZSR(NS(J))
         xr = XSR(NR(J)); zr = ZSR(NR(J))

         dist = SQRT((xs - xr)**2 + (zs - zr)**2)
         idist = INT(dist/bin_w) + 1
         idist = MAX(1, MIN(nbin, idist))

         w_geom = weight_table(idist)

         ! component weight
         rc = MAX(1, MIN(3, NRV(J)))
         w_comp = c_weight(rc)

         ! final acquisition weight
         WD_acq(J) = w_geom*w_comp

         IF (dbg) THEN
            WRITE (*, '("J=",I5," dist=",ES12.4," id=",I3," w_geom=",ES10.3," comp=",I2," w_comp=",ES10.3," WD_acq=",ES10.3)') &
               J, dist, idist, w_geom, rc, w_comp, WD_acq(J)
         END IF

         IF (unit_file > 0) THEN
            WRITE (unit_file, '(I6,2X,I6,2X,ES12.4,2X,I2,2X,ES12.4,2X,ES12.4)') &
               J, idist, w_geom, rc, w_comp, WD_acq(J)
         END IF
      END DO

      ! cleanup
      IF (ALLOCATED(weight_table)) DEALLOCATE (weight_table)
      IF (unit_file > 0) CLOSE (unit_file)

! ------------------------------------------------------------------
! NORMALIZATION
! Normalize acquisition weights so that SUM(WD_acq)=1.
! This makes the data misfit an average weighted misfit rather than
! a quantity that scales with ND.
! ------------------------------------------------------------------
      sum_w = SUM(WD_acq)

      IF (sum_w > tiny_dp) THEN
         WD_acq = WD_acq/sum_w
      ELSE
         WD_acq = 1.0_dp/REAL(ND, dp)
      END IF

   END SUBROUTINE ComputeAcqWeight

   SUBROUTINE ApplyQFrechetScaleCorrection(F_SCALE_INOUT, F_SCALE_L2_INOUT, PARAM_IN, INVP_IN, NPAR_IN, IM_IN, &
                                           Q_CORR, my_rank_in, DEBUG_OUTPUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPAR_IN, IM_IN, my_rank_in
      INTEGER, INTENT(IN) :: INVP_IN(:)
      REAL(dp), INTENT(INOUT) :: F_SCALE_INOUT(:), F_SCALE_L2_INOUT(:)
      REAL(dp), INTENT(IN) :: Q_CORR
      CHARACTER(LEN=*), INTENT(IN) :: PARAM_IN(:)
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      LOGICAL :: dbg
      INTEGER :: ia_loc

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      IF (Q_CORR == 1.0_dp) RETURN
      IF (IM_IN <= 1) RETURN
      IF (SIZE(F_SCALE_INOUT) < NPAR_IN .OR. SIZE(F_SCALE_L2_INOUT) < NPAR_IN) RETURN
      IF (SIZE(INVP_IN) < NPAR_IN .OR. SIZE(PARAM_IN) < NPAR_IN) RETURN

      DO ia_loc = 1, NPAR_IN
         IF (INVP_IN(ia_loc) /= 1) CYCLE
         IF (.NOT. is_q_parameter_name(PARAM_IN(ia_loc))) CYCLE
         F_SCALE_INOUT(ia_loc) = F_SCALE_INOUT(ia_loc)*Q_CORR
         F_SCALE_L2_INOUT(ia_loc) = F_SCALE_L2_INOUT(ia_loc)*Q_CORR
         IF (dbg .AND. my_rank_in == 0) THEN
            WRITE (*, '(A,I3,2X,A,2X,A,1PE12.5)') 'Applied Q F_SCALE correction to IA=', ia_loc, &
               TRIM(PARAM_IN(ia_loc)), 'factor=', Q_CORR
         END IF
      END DO
   END SUBROUTINE ApplyQFrechetScaleCorrection
   LOGICAL FUNCTION is_q_parameter_name(param_name)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: param_name
      CHARACTER(LEN=LEN(param_name)) :: param_trimmed

      param_trimmed = ADJUSTL(TRIM(param_name))
      is_q_parameter_name = INDEX(param_trimmed, '{Q') == 1 .OR. INDEX(param_trimmed, '{q') == 1
   END FUNCTION is_q_parameter_name

   REAL(dp) FUNCTION percentile_masked_abs_real_dp(values, mask, percentile)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: values(:)
      REAL(dp), INTENT(IN) :: mask(:)
      REAL(dp), INTENT(IN) :: percentile

      INTEGER :: nvals, idx
      REAL(dp), PARAMETER :: tiny_dp_val = TINY(1.0_dp)
      REAL(dp), ALLOCATABLE :: valid(:)

      IF (SIZE(values) /= SIZE(mask)) THEN
         STOP 'percentile_masked_abs_real_dp: values and mask length mismatch'
      END IF

      nvals = COUNT(mask > tiny_dp_val)
      IF (nvals <= 0) THEN
         percentile_masked_abs_real_dp = 0.0_dp
         RETURN
      END IF

      valid = PACK(ABS(values), mask > tiny_dp_val)
      CALL quicksort_real_dp(valid, 1, nvals)
      idx = MAX(1, MIN(nvals, CEILING(percentile*REAL(nvals, dp))))
      percentile_masked_abs_real_dp = valid(idx)
      DEALLOCATE (valid)
   END FUNCTION percentile_masked_abs_real_dp

   SUBROUTINE ComputeParamBalanceScale( &
      GRAD_FOR_BALANCE, GMASK, NPAR, NBLOCK, INVP, &
      REF_PARAM, PCTL, W_MIN, W_MAX, &
      BALANCE_SCALE, my_rank, PARAM, DEBUG_OUTPUT)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPAR, NBLOCK, REF_PARAM, my_rank
      INTEGER, INTENT(IN) :: INVP(:)
      REAL(dp), INTENT(IN) :: GRAD_FOR_BALANCE(:), GMASK(:)
      REAL(dp), INTENT(IN) :: PCTL, W_MIN, W_MAX
      REAL(dp), INTENT(INOUT) :: BALANCE_SCALE(:)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: PARAM(:)
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      LOGICAL :: dbg
      INTEGER :: II, active_count, active_idx, cs, ce, ref_idx
      REAL(dp), PARAMETER :: tiny_dp_val = TINY(1.0_dp)
      REAL(dp) :: amp_ref, amp_ii, w
      REAL(dp), ALLOCATABLE :: amp(:)
      CHARACTER(LEN=32) :: pname

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT

      IF (SIZE(GRAD_FOR_BALANCE) /= SIZE(GMASK)) THEN
         STOP 'ComputeParamBalanceScale: GRAD_FOR_BALANCE and GMASK length mismatch'
      END IF
      IF (SIZE(BALANCE_SCALE) < NPAR) THEN
         STOP 'ComputeParamBalanceScale: BALANCE_SCALE too small'
      END IF
      IF (SIZE(INVP) < NPAR) THEN
         STOP 'ComputeParamBalanceScale: INVP too small'
      END IF

      ALLOCATE (amp(NPAR))
      amp = 0.0_dp
      active_count = COUNT(INVP == 1)

      IF (active_count <= 1) THEN
         BALANCE_SCALE(1:NPAR) = 1.0_dp
         DEALLOCATE (amp)
         RETURN
      END IF

      active_idx = 0
      DO II = 1, NPAR
         IF (INVP(II) /= 1) CYCLE
         active_idx = active_idx + 1
         cs = (active_idx - 1)*NBLOCK + 1
         ce = active_idx*NBLOCK
         IF (cs < 1 .OR. ce > SIZE(GRAD_FOR_BALANCE) .OR. ce > SIZE(GMASK)) THEN
            WRITE (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
               'ComputeParamBalanceScale: invalid balance slice ', II, ' cs=', cs, ' ce=', ce, &
               ' size(grads)=', SIZE(GRAD_FOR_BALANCE), ' size(gmask)=', SIZE(GMASK)
            STOP 'ComputeParamBalanceScale: invalid slice bounds'
         END IF
         ! amp(II) = percentile_masked_abs_real_dp(GRAD_FOR_BALANCE(cs:ce), GMASK(cs:ce), PCTL)
         amp(II) = maxval(ABS(GRAD_FOR_BALANCE(cs:ce)*GMASK(cs:ce)))
      END DO

      ref_idx = 0
      IF (REF_PARAM >= 1 .AND. REF_PARAM <= NPAR) THEN
         IF (INVP(REF_PARAM) == 1 .AND. amp(REF_PARAM) > tiny_dp_val) THEN
            ref_idx = REF_PARAM
         END IF
      END IF
      IF (ref_idx == 0) THEN
         DO II = 1, NPAR
            IF (INVP(II) /= 1) CYCLE
            IF (amp(II) > tiny_dp_val) THEN
               ref_idx = II
               EXIT
            END IF
         END DO
      END IF

      IF (ref_idx == 0) THEN
         BALANCE_SCALE(1:NPAR) = 1.0_dp
         DEALLOCATE (amp)
         RETURN
      END IF

      amp_ref = amp(ref_idx)
      DO II = 1, NPAR
         IF (INVP(II) /= 1) THEN
            BALANCE_SCALE(II) = 1.0_dp
            CYCLE
         END IF
         amp_ii = amp(II)
         w = amp_ref/MAX(amp_ii, tiny_dp_val)
         w = MAX(W_MIN, MIN(W_MAX, w))
         BALANCE_SCALE(II) = round_to_one_significant_digit(w)
      END DO
      BALANCE_SCALE(ref_idx) = 1.0_dp

      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, '(A)') 'ComputeParamBalanceScale: balance factors:'
         DO II = 1, NPAR
            IF (INVP(II) /= 1) CYCLE
            IF (PRESENT(PARAM)) THEN
               pname = TRIM(ADJUSTL(PARAM(II)))
            ELSE
               WRITE (pname, '(I0)') II
            END IF
            WRITE (*, '(A,I3,A,A,A,1PE12.5,A,1PE12.5)') &
               ' IA=', II, ' PARAM=', TRIM(pname), ' amp=', amp(II), ' scale=', BALANCE_SCALE(II)
         END DO
      END IF

      DEALLOCATE (amp)
   CONTAINS
      PURE FUNCTION round_to_one_significant_digit(x) RESULT(y)
         REAL(dp), INTENT(IN) :: x
         REAL(dp) :: y, exponent, scale

         IF (x == 0.0_dp) THEN
            y = 0.0_dp
            RETURN
         END IF
         exponent = FLOOR(LOG10(ABS(x)))
         scale = 10.0_dp**exponent
         y = NINT(x/scale)*scale
      END FUNCTION round_to_one_significant_digit
   END SUBROUTINE ComputeParamBalanceScale

   RECURSIVE SUBROUTINE quicksort_real_dp(A, L, R)
      REAL(dp), INTENT(INOUT) :: A(:)
      INTEGER, INTENT(IN) :: L, R
      INTEGER :: i, j
      REAL(dp) :: pivot, tmp

      i = L
      j = R
      pivot = A((L + R)/2)

      DO
         DO WHILE (A(i) < pivot)
            i = i + 1
         END DO

         DO WHILE (pivot < A(j))
            j = j - 1
         END DO

         IF (i <= j) THEN
            tmp = A(i)
            A(i) = A(j)
            A(j) = tmp
            i = i + 1
            j = j - 1
         END IF

         IF (i > j) EXIT
      END DO

      IF (L < j) CALL quicksort_real_dp(A, L, j)
      IF (i < R) CALL quicksort_real_dp(A, i, R)
   END SUBROUTINE quicksort_real_dp
end module scalers_mod
