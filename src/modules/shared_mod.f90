module shared_mod
   use iso_fortran_env, only: dp => real64
   use constant_mod, only: BYTES_GB, BYTES_CPLX16
   implicit none

   public :: DNRM2   ! optional if you set PRIVATE by default

   interface
      real(dp) function DNRM2(n, x, incx)
         import dp
         integer, intent(in) :: n, incx
         real(dp), intent(in) :: x(*)
      end function DNRM2
   end interface
   !--------------------------------------
   !subroutines used in multiple modules
   !FFT_GL(main, modeling, Frechet lbfgs)
   !CLDI (grid_mod, F_modeling,)
   !ZH
   !--------------------------------------

contains

! SUBROUTINE FFT_GL(IS,IR,N,X,W,Y,R1,R2,F1,F2)
   !-----------------------------------------------------------------------
!     FFT_GL computes a cosine or sine-weighted quadrature transform
!     over real-valued input arrays. Used in 2.5D modeling  where
!    source-receiver symmetry, sign changes are handled via parity-based weighting.
!
!     Behavior:
!     - Applies either cosine or sine kernel based on (-1)^(IS+IR).
!     - Operates over Gauss-Legendre or similar quadrature points X,
!       using associated weights W.
!
!
!     Inputs:
!       IS.................... Source index
!       IR.................... Receiver index
!       N..................... Number of quadrature points
!       X(N).................. Quadrature abscissae (e.g., Gauss-Legendre)
!       W(N).................. Quadrature weights
!       Y..................... Wavenumber or phase scaling coefficient
!       R1(N)................. Real-valued array (real part)
!       R2(N)................. Real-valued array (imag part)
!
!     Outputs:
!       F1.................... Accumulated real result
!       F2.................... Accumulated imaginary result
!
!-----------------------------------------------------------------------

   SUBROUTINE FFT_GL(IS, IR, N, X, W, Y, R1, R2, F1, F2)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IS, IR, N
  REAL(dp), CONTIGUOUS, INTENT(IN) :: X(:), W(:), R1(:), R2(:)
  REAL(dp), INTENT(IN) :: Y
  REAL(dp), INTENT(OUT) :: F1, F2
  INTEGER :: I
  REAL(dp) :: yy, t, c, s
  LOGICAL :: even_parity

  even_parity = MOD(IS + IR, 2) == 0
  yy = Y
  F1 = 0.0_dp
  F2 = 0.0_dp

  ! Preserve the legacy 2D/NK=1 semantics used in the old code.
  IF (N == 1) THEN
     IF (even_parity) THEN
        F1 = R1(1)
        F2 = R2(1)
     END IF
     RETURN
  END IF

  IF (even_parity) THEN
     !$omp simd reduction(+:F1,F2)
     DO I = 1, N
        t = yy*X(I)
        c = COS(t)
        F1 = F1 + W(I)*R1(I)*c
        F2 = F2 + W(I)*R2(I)*c
     END DO
  ELSE
     !$omp simd reduction(+:F1,F2)
     DO I = 1, N
        t = yy*X(I)
        s = SIN(t)
        F2 = F2 + W(I)*R1(I)*s
        F1 = F1 - W(I)*R2(I)*s
     END DO
  END IF
END SUBROUTINE FFT_GL

SUBROUTINE FFT_GL_PRECOMP(IS, IR, N, X, W, Y, WC, WS, even_parity)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IS, IR, N
  REAL(dp), CONTIGUOUS, INTENT(IN) :: X(:), W(:)
  REAL(dp), INTENT(IN) :: Y
  REAL(dp), CONTIGUOUS, INTENT(OUT) :: WC(:), WS(:)
  LOGICAL, INTENT(OUT) :: even_parity
  INTEGER :: I
  REAL(dp) :: t

  even_parity = MOD(IS + IR, 2) == 0

  IF (even_parity) THEN
     !$omp simd
     DO I = 1, N
        t = Y * X(I)
        WC(I) = W(I) * COS(t)
     END DO
  ELSE
     !$omp simd
     DO I = 1, N
        t = Y * X(I)
        WS(I) = W(I) * SIN(t)
     END DO
  END IF
END SUBROUTINE


SUBROUTINE FFT_GL_W(IS, IR, N, WC, WS, even_parity, R1, R2, F1, F2)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IS, IR, N
  LOGICAL, INTENT(IN) :: even_parity
  REAL(dp), CONTIGUOUS, INTENT(IN) :: WC(:), WS(:), R1(:), R2(:)
  REAL(dp), INTENT(OUT) :: F1, F2
  INTEGER :: I

  F1 = 0.0_dp
  F2 = 0.0_dp

  ! Match FFT_GL for the 2D/NK=1 case regardless of precomputed weights.
  IF (N == 1) THEN
     IF (even_parity) THEN
        F1 = R1(1)
        F2 = R2(1)
     END IF
     RETURN
  END IF

  IF (even_parity) THEN
     !$omp simd reduction(+:F1,F2)
     DO I = 1, N
        F1 = F1 + WC(I) * R1(I)
        F2 = F2 + WC(I) * R2(I)
     END DO
  ELSE
     !$omp simd reduction(+:F1,F2)
     DO I = 1, N
        F2 = F2 + WS(I) * R1(I)
        F1 = F1 - WS(I) * R2(I)
     END DO
  END IF
END SUBROUTINE

   SUBROUTINE split_work(total_items, my_rank, n_process, start_idx, end_idx, task_count)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: total_items, my_rank, n_process
      INTEGER, INTENT(OUT) :: start_idx, end_idx, task_count
      INTEGER :: per_proc, remainder

      start_idx = 0
      end_idx = 0
      task_count = 0

      ! nothing to do
      IF (total_items <= 0) RETURN

      ! single process or serial
      IF (n_process <= 1) THEN
         start_idx = 1
         end_idx = total_items
         task_count = total_items
         RETURN
      END IF

      ! ranks beyond pool get nothing
      IF (my_rank < 0 .OR. my_rank >= n_process) THEN
         RETURN
      END IF

      ! more ranks than work: give first total_items ranks exactly 1 item
      IF (n_process > total_items) THEN
         IF (my_rank < total_items) THEN
            start_idx = my_rank + 1
            end_idx = start_idx
            task_count = 1
         END IF
         RETURN
      END IF

      ! balanced partition with remainder
      per_proc = total_items/n_process
      remainder = MOD(total_items, n_process)

      IF (my_rank < remainder) THEN
         start_idx = my_rank*(per_proc + 1) + 1
         end_idx = start_idx + per_proc
      ELSE
         start_idx = my_rank*per_proc + remainder + 1
         end_idx = start_idx + per_proc - 1
      END IF

      ! clamp & count
      IF (start_idx > total_items) THEN
         start_idx = 0
         end_idx = 0
         task_count = 0
      ELSE
         end_idx = MIN(end_idx, total_items)
         task_count = end_idx - start_idx + 1
      END IF
   END SUBROUTINE split_work

   !----------------------------------------------------------C
   !                                                          C
   !     derivatives of Lagrange function                     C
   !                                                          C
   !     Entries:                                             C
   !                                                          C
   !            X(N)..................N points for Lagrange;  C
   !            Xi.......................interpolated point;  C
   !                                                          C
   !     Returns:                                             C
   !                                                          C
   !            DLI(N).........Lagrange functions at X point. C
   !                                                          C
   !----------------------------------------------------------C
   SUBROUTINE CDLI(XI, N, X, DLI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(dp), INTENT(IN) :: XI
      REAL(dp), INTENT(IN) :: X(:)
      REAL(dp), INTENT(OUT) :: DLI(:)
      INTEGER, ALLOCATABLE :: ID(:)
      INTEGER :: I, J, K, II, I1, J1
      REAL(dp) :: F1, F2

      ALLOCATE (ID(N - 1))

      DO I = 1, N

         F1 = 1.0_dp
         II = 0
         DO J = 1, N
            IF (J == I) CYCLE
            F1 = F1*(X(I) - X(J))
            II = II + 1
            ID(II) = J
         END DO

         DLI(I) = 0.0_dp

         DO K = 1, II
            I1 = ID(K)
            F2 = 1.0_dp
            DO J = 1, II
               J1 = ID(J)
               IF (J1 == I1) CYCLE
               F2 = F2*(XI - X(J1))
            END DO
            DLI(I) = DLI(I) + F2/F1
         END DO

      END DO

      DEALLOCATE (ID)
      RETURN
   END SUBROUTINE CDLI

   PURE FUNCTION itoa(i) RESULT(s)
      INTEGER, INTENT(IN) :: i
      CHARACTER(LEN=32)   :: s
      WRITE (s, '(I0)') i
   END FUNCTION itoa

   REAL(dp) FUNCTION ZH(NTO, XTO, ZTO, XI)
!-----------------------------------------------------------------------
!
!     ZH performs linear interpolation to compute the z-coordinate
!     corresponding to a given x-coordinate (XI) based on the provided
!     topography data.
!
!     Inputs:
!       NTO............... Number of topography points
!       XTO(NTO).......... x-coordinates of the topography points
!       ZTO(NTO).......... z-coordinates of the topography points
!       XI................ x-coordinate for which the z-coordinate is needed
!
!     Outputs:
!       ZH................ Interpolated z-coordinate corresponding to XI
!
!     Behavior:
!       - If XI is less than the smallest x-coordinate, returns the
!         corresponding smallest z-coordinate.
!       - If XI is greater than the largest x-coordinate, returns the
!         corresponding largest z-coordinate.
!       - Otherwise, performs linear interpolation between the two
!         nearest points.
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTO
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), XI
      REAL(dp) :: X1, X2, ZH1, ZH2
      INTEGER :: I

      IF (XI <= XTO(1)) THEN
         ZH = ZTO(1)
         RETURN
      END IF

      IF (XI >= XTO(NTO)) THEN
         ZH = ZTO(NTO)
         RETURN
      END IF

      DO I = 1, NTO - 1
         X1 = XTO(I)
         X2 = XTO(I + 1)
         IF ((XI >= X1) .AND. (XI <= X2)) THEN
            ZH1 = ZTO(I)
            ZH2 = ZTO(I + 1)
            ZH = (ZH2 - ZH1)*(XI - X1)/(X2 - X1) + ZH1
            RETURN
         END IF
      END DO

      ! fallback (should not hit if inputs are consistent)
      ZH = ZTO(NTO)
      RETURN
   END FUNCTION ZH

   SUBROUTINE ConvertThomsenModuli(CR, CI, NPT, reverse)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPT
      LOGICAL, INTENT(IN), OPTIONAL :: reverse
      REAL(dp), INTENT(INOUT) :: CR(:, :), CI(:, :)

      REAL(dp) :: RHO, ALP, BET, EPS, DEL, GAM
      REAL(dp) :: C11, C13, C33, C44, C66
      REAL(dp) :: QP, QS, QEPS, QDEL, QGAM
      REAL(dp) :: Q11, Q13, Q33, Q44, Q66
      REAL(dp) :: A, B
      INTEGER :: IP
      LOGICAL :: backward

      backward = .FALSE.
      IF (PRESENT(reverse)) backward = reverse

      DO IP = 1, NPT
         IF (.NOT. backward) THEN
            ! ---------------------------
            ! Thomsen → Moduli direction
            ! ---------------------------
            RHO = CR(1, IP)
            ALP = CR(2, IP)
            BET = CR(3, IP)
            EPS = CR(4, IP)
            DEL = CR(5, IP)
            GAM = CR(6, IP)

            C11 = RHO*(1.0_dp + 2.0_dp*EPS)*ALP**2
            C13 = RHO*SQRT(ALP**4*DEL + (ALP**2 - BET**2)*((1.0_dp + EPS)*ALP**2 - BET**2)) - RHO*BET**2
            C33 = RHO*ALP**2
            C44 = RHO*BET**2
            C66 = RHO*(1.0_dp + 2.0_dp*GAM)*BET**2

            CR(2, IP) = C11
            CR(3, IP) = C13
            CR(4, IP) = C33
            CR(5, IP) = C44
            CR(6, IP) = C66

            ! Q conversion
            QP = CI(2, IP)
            QS = CI(3, IP)
            QEPS = CI(4, IP)
            QDEL = CI(5, IP)
            QGAM = CI(6, IP)

            Q33 = QP
            Q44 = QS
            Q11 = Q33/(1.0_dp + QEPS)
            Q66 = Q44/(1.0_dp + QGAM)

            A = C33*(C33 - C44)/(2.0_dp*C13*(C13 + C44))
            B = C44*(C13 + C33)**2/(2.0_dp*C13*(C13 + C44)*(C33 - C44))
            Q13 = 1.0_dp/((1.0_dp + A/QDEL + B)/Q33 - B/Q44)

            CI(2, IP) = Q11
            CI(3, IP) = Q13
            CI(4, IP) = Q33
            CI(5, IP) = Q44
            CI(6, IP) = Q66

         ELSE
            ! ---------------------------
            ! Moduli → Thomsen direction
            ! ---------------------------
            RHO = CR(1, IP)
            C11 = CR(2, IP)
            C13 = CR(3, IP)
            C33 = CR(4, IP)
            C44 = CR(5, IP)
            C66 = CR(6, IP)

            ALP = SQRT(C33/RHO)
            BET = SQRT(C44/RHO)
            EPS = 0.5_dp*(C11 - C33)/C33
            DEL = (2.0_dp*(C13 + C44)**2 - (C33 - C44)*(C11 + C33 - 2.0_dp*C44))/(2.0_dp*C33**2)
            GAM = 0.5_dp*(C66 - C44)/C44

            CR(2, IP) = ALP
            CR(3, IP) = BET
            CR(4, IP) = EPS
            CR(5, IP) = DEL
            CR(6, IP) = GAM

            ! Inverse Q conversion
            Q11 = CI(2, IP)
            Q13 = CI(3, IP)
            Q33 = CI(4, IP)
            Q44 = CI(5, IP)
            Q66 = CI(6, IP)

            QEPS = Q33/Q11 - 1.0_dp
            QGAM = Q44/Q66 - 1.0_dp

            A = C33*(C33 - C44)/(2.0_dp*C13*(C13 + C44))
            B = C44*(C13 + C33)**2/(2.0_dp*C13*(C13 + C44)*(C33 - C44))
            QDEL = A/(1.0_dp/Q13 - B*(1.0_dp/Q33 - 1.0_dp/Q44)) - B

            CI(2, IP) = Q33
            CI(3, IP) = Q44
            CI(4, IP) = QEPS
            CI(5, IP) = QDEL
            CI(6, IP) = QGAM

         END IF
      END DO

   END SUBROUTINE ConvertThomsenModuli

   !-----Take average
   REAL(dp) FUNCTION aver(tem, t)
      IMPLICIT NONE
      integer, intent(in) ::t
      real(dp), intent(in) :: tem(:)
      aver = SUM(tem)/real(t, dp)
      return
   end function aver

   !----------------------------------------------------------------------C
   !                                                                      C
   !     DFFT: Discrete Fast Fourier Transform.                           C
   !                                                                      C                                   C
   !     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N), (forward)           C
   !                                                                      C
   !     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N),  (inverse)   C
   !                                                                      C
   !     Entries:                                                         C
   !                                                                      C
   !     (1) x...real array containing real parts of transform sequence;  C
   !     (2) y...real array containing imag parts of transform sequence;  C
   !     (3) n.........integer length of transform:=64,128,256,512,1024;  C
   !     (4) m.........integer such that n = 2**m: = 6,  7,  8,  9,  10;  C
   !     (5)itype.....................=1 forward, =-1 inverse transform;  C
   !                                                                      C
   !     Returns:                                                         C
   !                                                                      C
   !     (1) x(*)..............itype=1,real part/itype=-1,real sequence;  C
   !     (2) y(*).................itype=1,imag.part/itype=-1,zero array;  C
   !                                                                      C
   !                                                                      C
   !----------------------------------------------------------------------C
   SUBROUTINE FFT(X, Y, N, M, ITYPE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, M, ITYPE
      REAL(dp), INTENT(INOUT) :: X(:), Y(:)
      REAL(dp), PARAMETER :: TWOPI = 6.2831853071795864769_dp
      REAL(dp) :: A, E, A3, CC1, CC3, SS1, SS3
      REAL(dp) :: R1, R2, S1, S2, S3, XT
      INTEGER :: I0, I1, I2, I3, IS, ID, I, J, K, N1, N2, N4
      IF (N .EQ. 1) RETURN
      IF (ITYPE .EQ. -1) THEN
      DO 1, I = 1, N
         Y(I) = -Y(I)
1        CONTINUE
         END IF

         N2 = 2*N
         DO 10, K = 1, M - 1
            N2 = N2/2
            N4 = N2/4
            E = TWOPI/N2
            A = 0.0_dp
            DO 20, J = 1, N4
               A3 = 3*A
               CC1 = COS(A)
               SS1 = SIN(A)
               CC3 = COS(A3)
               SS3 = SIN(A3)
               A = J*E
               IS = J
               ID = 2*N2
40             DO 30, I0 = IS, N - 1, ID
                  I1 = I0 + N4
                  I2 = I1 + N4
                  I3 = I2 + N4
                  R1 = X(I0) - X(I2)
                  X(I0) = X(I0) + X(I2)
                  R2 = X(I1) - X(I3)
                  X(I1) = X(I1) + X(I3)
                  S1 = Y(I0) - Y(I2)
                  Y(I0) = Y(I0) + Y(I2)
                  S2 = Y(I1) - Y(I3)
                  Y(I1) = Y(I1) + Y(I3)
                  S3 = R1 - S2
                  R1 = R1 + S2
                  S2 = R2 - S1
                  R2 = R2 + S1
                  X(I2) = R1*CC1 - S2*SS1
                  Y(I2) = -S2*CC1 - R1*SS1
                  X(I3) = S3*CC3 + R2*SS3
                  Y(I3) = R2*CC3 - S3*SS3
30                CONTINUE
                  IS = 2*ID - N2 + J
                  ID = 4*ID
                  IF (IS .LT. N) GOTO 40
20                CONTINUE
10                CONTINUE

                  IS = 1
                  ID = 4
50                DO 60, I0 = IS, N, ID
                     I1 = I0 + 1
                     R1 = X(I0)
                     X(I0) = R1 + X(I1)
                     X(I1) = R1 - X(I1)
                     R1 = Y(I0)
                     Y(I0) = R1 + Y(I1)
                     Y(I1) = R1 - Y(I1)
60                   CONTINUE
                     IS = 2*ID - 1
                     ID = 4*ID
                     IF (IS .LT. N) GOTO 50

100                  J = 1
                     N1 = N - 1
                     DO 104, I = 1, N1
                        IF (I .GE. J) GOTO 101
                        XT = X(J)
                        X(J) = X(I)
                        X(I) = XT
                        XT = Y(J)
                        Y(J) = Y(I)
                        Y(I) = XT
101                     K = N/2
102                     IF (K .GE. J) GOTO 103
                        J = J - K
                        K = K/2
                        GOTO 102
103                     J = J + K
104                     CONTINUE

                        IF (ITYPE .EQ. -1) THEN
                        DO 2, I = 1, N
                           X(I) = X(I)/N
                           Y(I) = -Y(I)/N
2                          CONTINUE
                           END IF

                           RETURN
                        END
                        subroutine log_rss(tag, unit, my_rank)
                           use iso_fortran_env, only: int64, real64
                           implicit none
                           character(len=*), intent(in) :: tag
                           integer, intent(in) :: unit, my_rank

                           integer :: iu, ios
                           character(len=256) :: line
                           integer(int64) :: rss_kb, hwm_kb
                           integer(int64) :: size_pages, rss_pages

                           rss_kb = -1_int64
                           hwm_kb = -1_int64

                           open (newunit=iu, file="/proc/self/status", status="old", action="read", iostat=ios)
                           if (ios /= 0) then
                              if (my_rank == 0) write (unit, '(A)') trim(tag)//" | /proc/self/status unavailable"
                              return
                           end if

                           do
                              read (iu, '(A)', iostat=ios) line
                              if (ios /= 0) exit
                              if (index(line, "VmRSS:") > 0) then
                                 read (line(index(line,"VmRSS:")+6:), *, iostat=ios) rss_kb
                              else if (index(line, "VmHWM:") > 0) then
                                 read (line(index(line,"VmHWM:")+6:), *, iostat=ios) hwm_kb
                              end if
                           end do
                           close (iu)

                           ! Fallback: /proc/self/statm (RSS pages) if VmRSS not found
                           if (rss_kb < 0_int64) then
                              size_pages = -1_int64
                              rss_pages  = -1_int64
                              open (newunit=iu, file="/proc/self/statm", status="old", action="read", iostat=ios)
                              if (ios == 0) then
                                 read (iu, *, iostat=ios) size_pages, rss_pages
                                 close (iu)
                                 if (ios == 0 .and. rss_pages > 0_int64) then
                                    ! statm reports pages; assume 4KB pages if no sysconf available
                                    rss_kb = rss_pages * 4_int64
                                 end if
                              end if
                           end if

                           if (hwm_kb < 0_int64 .and. rss_kb >= 0_int64) then
                              hwm_kb = rss_kb
                           end if

                           if (my_rank == 0) then
                              if (rss_kb < 0_int64) rss_kb = 0_int64
                              if (hwm_kb < 0_int64) hwm_kb = 0_int64
                              write (unit, '(A,2(A,F8.3),A)') trim(tag), &
                                 " | VmRSS=", real(rss_kb, real64)/1024.0_real64/1024.0_real64, &
                                 " GB  VmHWM=", real(hwm_kb, real64)/1024.0_real64/1024.0_real64, " GB"
                           end if
                        end subroutine log_rss

                        SUBROUTINE FFT0(FKC, N, R1, R2, F1, F2)
!-----------------------------------------------------------------------
!     FFT0 computes a single-frequency weighted sum over real-valued
!     sequences,
!
!     Behavior:
!     - For N = 1, the function returns the input value as-is.
!     - For N > 1, a normalization based on π and N is applied.
!     - The function assumes precomputed real (R1) and imaginary (R2)
!       components of a signal or spatial field.
!
!     Note: This is not a true FFT, but a simplified summation routine.
!
!     Inputs:
!       FKC.................. Scaling constant for normalization
!       N.................... Number of samples
!       R1(N)................ Real-valued array (real part)
!       R2(N)................ Real-valued array (imag part)
!
!     Outputs:
!       F1................... Accumulated real result (projection)
!       F2................... Accumulated imaginary result (projection)
!
!-----------------------------------------------------------------------
                           IMPLICIT NONE
                           REAL(dp), INTENT(IN) :: FKC
                           INTEGER, INTENT(IN) :: N
                           REAL(dp), INTENT(IN) :: R1(:), R2(:)
                           REAL(dp), INTENT(OUT) :: F1, F2
                           REAL(dp) :: FACT
                           INTEGER :: I
                           IF (N .EQ. 1) THEN
                              F1 = R1(N)
                              F2 = R2(N)
                           ELSE
                              FACT = FKC/(3.1415926_dp*REAL(N - 1, dp))
                              F1 = 0.0_dp
                              F2 = 0.0_dp
                              DO I = 1, N
                                 F1 = F1 + FACT*R1(I)
                                 F2 = F2 + FACT*R2(I)
                              END DO
                           END IF
                           RETURN
                        END

                        PURE REAL(dp) FUNCTION safe_div(a, b, eps)
                           implicit none
                           REAL(dp), INTENT(IN) :: a, b, eps
                           safe_div = a/MAX(b, eps)
                        END FUNCTION safe_div

                        SUBROUTINE recompute_masked_norms(GRAD_vec, GRAD_norm_vec, INVP, NBLOCK)

                           USE iso_fortran_env, ONLY: dp => real64
                           IMPLICIT NONE
                           REAL(dp), INTENT(IN) :: GRAD_vec(:)
                           REAL(dp), INTENT(OUT) :: GRAD_norm_vec(:)
                           INTEGER, INTENT(IN) :: INVP(:)
                           INTEGER, INTENT(IN) :: NBLOCK
                           INTEGER :: II, IM, cs, ce

                           IM = 0
                           DO II = 1, SIZE(INVP)
                              IF (INVP(II) /= 1) THEN
                                 GRAD_norm_vec(II) = 0.0_dp
                                 CYCLE
                              END IF
                              IM = IM + 1
                              cs = (IM - 1)*NBLOCK + 1
                              ce = IM*NBLOCK
                              GRAD_norm_vec(II) = DNRM2(NBLOCK, GRAD_vec(cs:ce), 1)
                           END DO
                        END SUBROUTINE recompute_masked_norms

                        ! SUBROUTINE FFT_GLCall(D,IS,IR NK, FKY, WTK, YY, R1, R2, F1, F2, R,)
                        ! FFT_GL(IS,IR,N,X,W,Y,R1,R2,F1,F2)
                        !     IMPLICIT real(dp) (A-H, O-Z)
                        !     INTEGER, INTENT(IN) :: NK
                        !     COMPLEX(dp), INTENT(IN) :: D(:)
                        !     real(dp), INTENT(IN) :: FKY(:), WTK(:), YY
                        !     real(dp), INTENT(INOUT) :: R1(:), R2(:)
                        !     real(dp) :: F1, F2
                        !     real(dp), INTENT(OUT) :: R(:,:)

                        !     R1(1:NK) = DBLE(REAL(D(1:NK)))
                        !     R2(1:NK) = DBLE(AIMAG(D(1:NK)))
                        !     CALL FFT_GL(1, 1, NK, FKY, WTK, YY, R1, R2, F1, F2)
                        !     R(IS,IR)= DCMPLX(F1,F2)
                        !   END SUBROUTINE FFT_GLCall

!     !------------------ build_active_map.f90 ------------------
! SUBROUTINE BuildActiveMap(INVP, NPAR, NBLOCK, NPT,             &
!                           act_ids, ip_of_K,                    &
!                           pt_s_of_ip, pt_e_of_ip,              &
!                           blk_s_of_ip, blk_e_of_ip,            &
!                           IM, NMM, NM)
!   IMPLICIT NONE
!   INTEGER, INTENT(IN)  :: NPAR, NBLOCK, NPT
!   INTEGER, INTENT(IN)  :: INVP(:)            ! 1..NPAR
!   INTEGER, INTENT(OUT) :: act_ids(:)         ! 1..IM  (global param IDs K)
!   INTEGER, INTENT(OUT) :: ip_of_K(:)         ! 1..NPAR (0 if inactive, else ip_idx)
!   INTEGER, INTENT(OUT) :: pt_s_of_ip(:), pt_e_of_ip(:)     ! 1..IM
!   INTEGER, INTENT(OUT) :: blk_s_of_ip(:), blk_e_of_ip(:)   ! 1..IM
!   INTEGER, INTENT(OUT) :: IM, NMM, NM

!   INTEGER :: K, ip

!   ! initialize
!   ip_of_K = 0
!   ip      = 0

!   ! enumerate active families in increasing K
!   DO K = 1, NPAR
!     IF (INVP(K) == 1) THEN
!       ip = ip + 1
!       act_ids(ip) = K
!       ip_of_K(K)  = ip
!     END IF
!   END DO

!   IM  = ip
!   NMM = IM * NPT
!   NM  = IM * NBLOCK

!   DO ip = 1, IM
!     pt_s_of_ip(ip)  = (ip-1)*NPT   + 1
!     pt_e_of_ip(ip)  =  ip   *NPT
!     blk_s_of_ip(ip) = (ip-1)*NBLOCK + 1
!     blk_e_of_ip(ip) =  ip   *NBLOCK
!   END DO
!   !use
!   ! after you set INVP(:) for this stage:
! ! allocate temps once per iter (sizes: act_ids/ip_of_K = NPAR, others = <= NPAR)
! ALLOCATE(act_ids(NPAR), ip_of_K(NPAR))
! ALLOCATE(pt_s_of_ip(NPAR), pt_e_of_ip(NPAR))
! ALLOCATE(blk_s_of_ip(NPAR), blk_e_of_ip(NPAR))

! CALL BuildActiveMap(INVP, NPAR, NBLOCK, NPT,                        &
!                     act_ids, ip_of_K,                               &
!                     pt_s_of_ip, pt_e_of_ip,                         &
!                     blk_s_of_ip, blk_e_of_ip,                       &
!                     IM, NMM, NM)

! ! Now reuse these in every subroutine this iteration:
! ! - gradient assembly
! ! - UPDATEPARAMETER / GRID2D_*_UPDATE
! ! - exports (CFNAME_* + GRID2D_OUT2)
!                     ! Example: write updated CR1/CI1 for all active families
! DO K = 1, NPAR
!   IF (INVP(K) == 1) THEN
!     ip   = ip_of_K(K)
!     ps   = pt_s_of_ip(ip)
!     pe   = pt_e_of_ip(ip)
!     ! unscaled:
!     CR1(K,1:NPT) = PAR_STEP(ps:pe)
!     CI1(K,1:NPT) = CI(K,1:NPT)
!   ELSE IF (K <= IANISO) THEN
!     CR1(K,1:NPT) = CR(K,1:NPT)
!     CI1(K,1:NPT) = CI(K,1:NPT)
!   ELSE
!     ! Q-part indexing same rule you already use
!     CI1(K-(IANISO-1),1:NPT) = CI(K-(IANISO-1),1:NPT)
!   END IF
! END DO

! ip = ip_of_K(K)
! bs = blk_s_of_ip(ip)
! be = blk_e_of_ip(ip)
! ! ! e.g., export GRADr(bs:be) for family K

! END SUBROUTINE BuildActiveMap

                        end module shared_mod
