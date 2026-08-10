module stiffness_assembly_mod
    USE omp_lib          ! OpenMP runtime library
   USE mkl_service      ! Intel MKL runtime services
   USE MPI              ! MPI (message passing)
   USE hardware_mod     ! Hardware configuration (threads, memory)
   USE shared_mod       ! Shared variables/constants across program
   USE boundaries_mod   ! Absorbing boundary condition logic
   USE output_mod       ! Output file and logging utilities
   USE gridtype_mod
   USE err_mpi_mod     ! MPI error handling routines
   USE constant_mod   ! Physical and numerical constants
   use iso_fortran_env, only: dp => real64
    implicit none
    contains
    
   SUBROUTINE CA_core(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                      CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, &   ! you can pass these even if some modes don’t need them
                      ADD_ENTRY)   ! callback for storage

      IMPLICIT NONE
      ! ==== arguments ====
      INTEGER, INTENT(IN) :: N, KL, KU, LDAB
      INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
      REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0
      TYPE(InversionGridType), INTENT(IN) :: IG
      ! ==== allocatable temporaries ====
      REAL(dp), ALLOCATABLE :: Z1(:), Z2(:), T1(:), T2(:)
      REAL(dp), ALLOCATABLE :: NOX(:), NOZ(:), DLX(:), DLZ(:)
      REAL(dp), ALLOCATABLE :: DNX(:), DNZ(:), XP(:)
      INTEGER, ALLOCATABLE          :: INDX(:), INDZ(:)
      COMPLEX(dp), ALLOCATABLE       :: P(:), Q(:), CJ(:)

      ! ==== scalars ====
      INTEGER :: I, J, K, L, N0, MM, I0, II, JJ, NBLOCK, IM
      INTEGER :: IX, IZ, ID, K1, L1, IP, IQ, JW
      REAL(dp) :: Z0, OMIGA, X1, X2, DX, C1, XI, DM, S, AK, DZ, C2, C3, BL
      COMPLEX(dp) :: OMIGA2, RHO, WIJ, AIJ

      ! Dummy procedure interface for the callback
      INTERFACE
         SUBROUTINE ADD_ENTRY(row, col, val)
            USE iso_fortran_env, ONLY: dp => real64
            IMPLICIT NONE
            INTEGER, INTENT(IN)    :: row, col
            COMPLEX(dp), INTENT(IN) :: val
         END SUBROUTINE
      END INTERFACE

! write(*,*)"---- initialisation --------------"
      MM = 2*(NORD - 1) + 1
      I0 = KL + KU + 1
      Z0 = IE0*DZ0
      OMIGA = 2.0_dp*3.1415926_dp*FREQ
      OMIGA2 = CMPLX(OMIGA*OMIGA, 0.0_dp, dp)

      ALLOCATE (Z1(NORD), Z2(NORD), T1(NORD), T2(NORD))
      ALLOCATE (NOX(NORD), NOZ(NORD), INDX(MM), INDZ(NORD))
      ALLOCATE (DLX(NORD), DLZ(NORD), DNX(MM), DNZ(NORD))
      ALLOCATE (XP(NORD), P(21), Q(81), CJ(27))

      NBLOCK = SIZE(IG%DX_BLOCK)

      DO I = 1, NX - 1
         DO J = 1, NZ - 1

            ! ----------------------------------------
            ! Block index in IG (same ordering as before)
            ! ----------------------------------------
            IM = (I - 1)*(NZ - 1) + J

            ! --- geometry for this block from IG ---
            DX = IG%DX_BLOCK(IM)
            C1 = IG%C1_BLOCK(IM)
            N0 = IG%N0_BLOCK(IM)

            DO K = 1, NORD
               Z1(K) = IG%Z1_OUT(K, IM)
               Z2(K) = IG%Z2_OUT(K, IM)
               T1(K) = IG%T1_OUT(K, IM)
               T2(K) = IG%T2_OUT(K, IM)
            END DO

            ! ---- GQ over (K,L) ----
            DO K = 1, NORD
               AK = AS(K)
               DZ = Z2(K) - Z1(K)
               C3 = 2.D0/DZ
               CALL CDLI(AK, NORD, AS, DLX)

               DO L = 1, NORD
                  ID = N0 + (K - 1)*NNZ + (L - 1)

                  CALL C21(ID, IANISO, CR, CI, RHO, P)
                  CALL Q81_NewGSRM(FREQ, NX - 1, NZ - 1, I, J, K, L, NORD, IE0, IS0, RHO, P, Q)

                  WIJ = CMPLX(0.25_dp*DX*DZ*WT(K)*WT(L), 0.0_dp, dp)
                  BL = AS(L)
                  C2 = -((T2(K) - T1(K))*BL + (T1(K) + T2(K)))/DZ
                  CALL CDLI(BL, NORD, AS, DLZ)

                  ! build INDX/DNX and INDZ/DNZ
                  DO K1 = 1, NORD
                     NOX(K1) = (N0 + (L - 1)) + (K1 - 1)*NNZ
                     NOZ(K1) = (N0 + (K - 1)*NNZ) + (K1 - 1)
                  END DO

                  IX = 0
                  DO K1 = 1, K - 1
                     IX = IX + 1
                     INDX(IX) = NOX(K1)
                     DNX(IX) = C1*DLX(K1)
                  END DO
                  DO L1 = 1, L - 1
                     IX = IX + 1; INDX(IX) = NOZ(L1); DNX(IX) = C2*DLZ(L1)
                  END DO
                  IX = IX + 1
                  INDX(IX) = NOX(K)                ! NOX(K) = NOZ(L)
                  DNX(IX) = C1*DLX(K) + C2*DLZ(L)
                  DO L1 = L + 1, NORD
                     IX = IX + 1; INDX(IX) = NOZ(L1); DNX(IX) = C2*DLZ(L1)
                  END DO
                  DO K1 = K + 1, NORD
                     IX = IX + 1; INDX(IX) = NOX(K1); DNX(IX) = C1*DLX(K1)
                  END DO

                  DO L1 = 1, NORD
                     INDZ(L1) = NOZ(L1)
                     DNZ(L1) = C3*DLZ(L1)
                  END DO
                  IZ = NORD

                  ! sanity (your checks)
                  IF (NNX*NNZ > 0) THEN
                     DO K1 = 1, IX
                        IF (INDX(K1) < 1 .OR. INDX(K1) > N) THEN
                           WRITE (*, '(A,3I12)') 'CA: bad INDX', K1, INDX(K1), N
                           STOP
                        END IF
                     END DO
                     DO L1 = 1, IZ
                        IF (INDZ(L1) < 1 .OR. INDZ(L1) > N) THEN
                           WRITE (*, '(A,3I12)') 'CA: bad INDZ', L1, INDZ(L1), N
                           STOP
                        END IF
                     END DO
                  END IF

                  ! ===== weights loop (unchanged math) =====
                  DO JW = 1, 3
                     CALL QCJ(JW, Q, CJ)

                     ! ---- A-block: c1j11..,c1j13..,c3j11..,c3j13.. ----
                     !c1j11-term: (DxDx)Gx
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 2
                           AIJ = WIJ*CJ(1)*DCMPLX(DNX(IP)*DNX(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 2
                           AIJ = WIJ*CJ(2)*DCMPLX(DNX(IP)*DNZ(IQ), 0.D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 2
                           AIJ = WIJ*CJ(3)*DCMPLX(DNZ(IP)*DNX(IQ), 0.D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 2
                           AIJ = WIJ*CJ(4)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     JJ = 3*ID - 2
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = DCMPLX(0D0, -FK)*CJ(5)*WIJ*DCMPLX(DNX(IP), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = DCMPLX(0D0, -FK)*CJ(6)*WIJ*DCMPLX(DNZ(IP), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ) - 2
                        AIJ = -DCMPLX(0D0, -FK)*WIJ*CJ(7)*DCMPLX(DNX(IQ), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ) - 2
                        AIJ = -DCMPLX(0D0, -FK)*WIJ*CJ(8)*DCMPLX(DNZ(IQ), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     JJ = 3*ID - 2
                     IF (II == JJ) THEN
                        AIJ = DCMPLX(FK*FK, 0D0)*WIJ*CJ(9)
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     IF (JW == 1) THEN
                        AIJ = -RHO*OMIGA2*WIJ
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     ! ---- B-block (Gy) ----
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 1
                           AIJ = WIJ*CJ(10)*DCMPLX(DNX(IP)*DNX(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 1
                           AIJ = WIJ*CJ(11)*DCMPLX(DNX(IP)*DNZ(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 1
                           AIJ = WIJ*CJ(12)*DCMPLX(DNZ(IP)*DNX(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 1
                           AIJ = WIJ*CJ(13)*DCMPLX(DNZ(IP)*DNZ(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     JJ = 3*ID - 1
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = DCMPLX(0D0, -FK)*CJ(14)*WIJ*DCMPLX(DNX(IP), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = DCMPLX(0D0, -FK)*CJ(15)*WIJ*DCMPLX(DNZ(IP), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ) - 1
                        AIJ = -DCMPLX(0D0, -FK)*CJ(16)*WIJ*DCMPLX(DNX(IQ), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ) - 1
                        AIJ = -DCMPLX(0D0, -FK)*CJ(17)*WIJ*DCMPLX(DNZ(IQ), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     JJ = 3*ID - 1
                     IF (JJ == II) THEN
                        AIJ = WIJ*DCMPLX(FK*FK, 0D0)*CJ(18)
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     IF (JW == 2) THEN
                        AIJ = -WIJ*RHO*OMIGA2
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     ! ---- C-block (Gz) ----
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ)
                           AIJ = WIJ*CJ(19)*DCMPLX(DNX(IP)*DNX(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ)
                           AIJ = WIJ*CJ(20)*DCMPLX(DNX(IP)*DNZ(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ)
                           AIJ = WIJ*CJ(21)*DCMPLX(DNZ(IP)*DNX(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ)
                           AIJ = WIJ*CJ(22)*DCMPLX(DNZ(IP)*DNZ(IQ), 0D0)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     JJ = 3*ID
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = WIJ*DCMPLX(0D0, -FK)*CJ(23)*DCMPLX(DNX(IP), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = WIJ*DCMPLX(0D0, -FK)*CJ(24)*DCMPLX(DNZ(IP), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ)
                        AIJ = -WIJ*DCMPLX(0D0, -FK)*CJ(25)*DCMPLX(DNX(IQ), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ)
                        AIJ = -WIJ*DCMPLX(0D0, -FK)*CJ(26)*DCMPLX(DNZ(IQ), 0D0)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     JJ = 3*ID
                     IF (JJ == II) THEN
                        AIJ = WIJ*CJ(27)*DCMPLX(FK*FK, 0D0)
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     IF (JW == 3) THEN
                        AIJ = -WIJ*RHO*OMIGA2
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                  END DO  ! JW
               END DO    ! L
            END DO      ! K
         END DO        ! J
      END DO          ! I
      DEALLOCATE (Z1, Z2, T1, T2, DLX, DLZ, DNX, DNZ)
      DEALLOCATE (XP, P, NOX, NOZ, INDX, INDZ, Q, CJ)
   END SUBROUTINE CA_core

!============================= BANDED (LU) =============================
   SUBROUTINE CA_core_ori(  FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                 CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU,   &   ! you can pass these even if some modes don’t need them
            ADD_ENTRY )   ! callback for storage

  IMPLICIT NONE
  ! ==== arguments ====
  INTEGER, INTENT(IN) :: N, KL, KU, LDAB
  INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ
  DOUBLE PRECISION, INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
  DOUBLE PRECISION, INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0

  ! ==== allocatable temporaries ====
  DOUBLE PRECISION, ALLOCATABLE :: Z1(:), Z2(:), T1(:), T2(:)
  DOUBLE PRECISION, ALLOCATABLE :: NOX(:), NOZ(:), DLX(:), DLZ(:)
  DOUBLE PRECISION, ALLOCATABLE :: DNX(:), DNZ(:), XP(:)
  INTEGER, ALLOCATABLE          :: INDX(:), INDZ(:)
  COMPLEX*16, ALLOCATABLE       :: P(:), Q(:), CJ(:)

  ! ==== scalars ====
  INTEGER :: I, J, K, L, N0, MM, I0, II, JJ
  INTEGER :: IX, IZ, ID, K1, L1, IP, IQ, JW
  DOUBLE PRECISION :: Z0, OMIGA, X1, X2, DX, C1, XI, DM, S, AK, DZ, C2, C3, BL
  COMPLEX*16 :: OMIGA2, RHO, WIJ, AIJ

  ! Dummy procedure interface for the callback
  INTERFACE
    SUBROUTINE ADD_ENTRY(row, col, val)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: row, col
      COMPLEX*16, INTENT(IN) :: val
    END SUBROUTINE
  END INTERFACE

! write(*,*)"---- initialisation --------------"
      MM = 2*(NORD - 1) + 1
      I0 = KL + KU + 1
      Z0 = IE0*DZ0
      OMIGA = 2.D0*3.1415926D0*FREQ
      OMIGA2 = DCMPLX(OMIGA*OMIGA, 0.D0)

      ALLOCATE (Z1(NORD), Z2(NORD), T1(NORD), T2(NORD))
      ALLOCATE (NOX(NORD), NOZ(NORD), INDX(MM), INDZ(NORD))
      ALLOCATE (DLX(NORD), DLZ(NORD), DNX(MM), DNZ(NORD))
      ALLOCATE (XP(NORD), P(21), Q(81), CJ(27))

    DO I = 1, NX - 1
    X1 = X(I); X2 = X(I+1); DX = X2 - X1
    C1 = 2.D0 / DX

    DO J = 1, NZ - 1
      N0 = (I - 1)*(NORD - 1)*NNZ + (J - 1)*(NORD - 1) + 1

      ! --- z slabs: Z1/Z2 and their derivatives ---
      IF (J == 1) THEN
        DO K = 1, NORD
          Z1(K) = 0.D0
          Z2(K) = Z1(K) + DZ0
          T1(K) = 0.D0
          T2(K) = 0.D0
        END DO
      ELSEIF (J <= IE0) THEN
        DO K = 1, NORD
          Z1(K) = Z2(K)
          Z2(K) = Z1(K) + DZ0
          T1(K) = 0.D0
          T2(K) = 0.D0
        END DO
      ELSE
        DO K = 1, NORD
          XI = 0.5D0*(X2 - X1)*AS(K) + 0.5D0*(X1 + X2)
          XP(K) = XI
          DM = (ZH(NTO, XTO, ZTO, XI) - Z0)/DBLE(FLOAT((NZ - 1) - IE0))
          Z1(K) = Z2(K)
          T1(K) = T2(K)
          Z2(K) = Z1(K) + DM
        END DO
        DO K = 1, NORD
          XI = XP(K)
          CALL CDLI(XI, NORD, XP, DLX)
          S = 0.D0
          DO L = 1, NORD
            S = S + DLX(L)*Z2(L)
          END DO
          T2(K) = S
        END DO
      END IF

      ! ---- GQ over (K,L) ----
      DO K = 1, NORD
        AK = AS(K)
        DZ = Z2(K) - Z1(K)
        C3 = 2.D0 / DZ
        CALL CDLI(AK, NORD, AS, DLX)

        DO L = 1, NORD
          ID = N0 + (K - 1)*NNZ + (L - 1)

          CALL C21(ID, IANISO, CR, CI, RHO, P)
          CALL Q81_NewGSRM(FREQ, NX - 1, NZ - 1, I, J, K, L, NORD, IE0, IS0, RHO, P, Q)

          WIJ = DCMPLX(0.25D0*DX*DZ*WT(K)*WT(L), 0.D0)
          BL  = AS(L)
          C2 = -((T2(K) - T1(K))*BL + (T1(K) + T2(K))) / DZ
          CALL CDLI(BL, NORD, AS, DLZ)         

          ! build INDX/DNX and INDZ/DNZ
          DO K1 = 1, NORD
            NOX(K1) = (N0 + (L - 1)) + (K1 - 1)*NNZ
            NOZ(K1) = (N0 + (K - 1)*NNZ) + (K1 - 1)
          END DO

          IX = 0
          DO K1 = 1, K - 1
            IX = IX + 1
             INDX(IX) = NOX(K1)
              DNX(IX) = C1*DLX(K1)
          END DO
          DO L1 = 1, L - 1
            IX = IX + 1; INDX(IX) = NOZ(L1); DNX(IX) = C2*DLZ(L1)
          END DO
          IX = IX + 1
          INDX(IX) = NOX(K)                ! NOX(K) = NOZ(L)
          DNX(IX)  = C1*DLX(K) + C2*DLZ(L)
          DO L1 = L + 1, NORD
            IX = IX + 1; INDX(IX) = NOZ(L1); DNX(IX) = C2*DLZ(L1)
          END DO
          DO K1 = K + 1, NORD
            IX = IX + 1; INDX(IX) = NOX(K1); DNX(IX) = C1*DLX(K1)
          END DO

          DO L1 = 1, NORD
            INDZ(L1) = NOZ(L1)
            DNZ(L1)  = C3*DLZ(L1)
          END DO
          IZ = NORD

          ! sanity (your checks)
          IF (NNX*NNZ > 0) THEN
            DO K1 = 1, IX
              IF (INDX(K1) < 1 .OR. INDX(K1) > N) THEN
                WRITE(*,'(A,3I12)') 'CA: bad INDX', K1, INDX(K1), N
                STOP
              END IF
            END DO
            DO L1 = 1, IZ
              IF (INDZ(L1) < 1 .OR. INDZ(L1) > N) THEN
                WRITE(*,'(A,3I12)') 'CA: bad INDZ', L1, INDZ(L1), N
                STOP
              END IF
            END DO
          END IF

          ! ===== weights loop (unchanged math) =====
          DO JW = 1, 3
            CALL QCJ(JW, Q, CJ)

            ! ---- A-block: c1j11..,c1j13..,c3j11..,c3j13.. ----
            !c1j11-term: (DxDx)Gx
            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              DO IQ = 1, IX
                JJ = 3*INDX(IQ) - 2
                AIJ = WIJ*CJ(1)*DCMPLX(DNX(IP)*DNX(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              DO IQ = 1, IZ
                JJ = 3*INDZ(IQ) - 2
                AIJ = WIJ*CJ(2)*DCMPLX(DNX(IP)*DNZ(IQ),0.D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              DO IQ = 1, IX
                JJ = 3*INDX(IQ) - 2
                AIJ = WIJ*CJ(3)*DCMPLX(DNZ(IP)*DNX(IQ),0.D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              DO IQ = 1, IZ
                JJ = 3*INDZ(IQ) - 2
                AIJ = WIJ*CJ(4)*DCMPLX(DNZ(IP)*DNZ(IQ),0.D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            JJ = 3*ID - 2
            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              AIJ = DCMPLX(0D0,-FK)*CJ(5)*WIJ*DCMPLX(DNX(IP),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              AIJ = DCMPLX(0D0,-FK)*CJ(6)*WIJ*DCMPLX(DNZ(IP),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            II = 3*ID - (3 - JW)
            DO IQ = 1, IX
              JJ = 3*INDX(IQ) - 2
              AIJ = -DCMPLX(0D0,-FK)*WIJ*CJ(7)*DCMPLX(DNX(IQ),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            DO IQ = 1, IZ
              JJ = 3*INDZ(IQ) - 2
              AIJ = -DCMPLX(0D0,-FK)*WIJ*CJ(8)*DCMPLX(DNZ(IQ),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            JJ = 3*ID - 2
            IF (II == JJ) THEN
              AIJ = DCMPLX(FK*FK,0D0)*WIJ*CJ(9)
              CALL ADD_ENTRY(JJ, JJ, AIJ)
            END IF

            IF (JW == 1) THEN
              AIJ = -RHO*OMIGA2*WIJ
              CALL ADD_ENTRY(JJ, JJ, AIJ)
            END IF

            ! ---- B-block (Gy) ----
            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              DO IQ = 1, IX
                JJ = 3*INDX(IQ) - 1
                AIJ = WIJ*CJ(10)*DCMPLX(DNX(IP)*DNX(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              DO IQ = 1, IZ
                JJ = 3*INDZ(IQ) - 1
                AIJ = WIJ*CJ(11)*DCMPLX(DNX(IP)*DNZ(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              DO IQ = 1, IX
                JJ = 3*INDX(IQ) - 1
                AIJ = WIJ*CJ(12)*DCMPLX(DNZ(IP)*DNX(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              DO IQ = 1, IZ
                JJ = 3*INDZ(IQ) - 1
                AIJ = WIJ*CJ(13)*DCMPLX(DNZ(IP)*DNZ(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            JJ = 3*ID - 1
            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              AIJ = DCMPLX(0D0,-FK)*CJ(14)*WIJ*DCMPLX(DNX(IP),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              AIJ = DCMPLX(0D0,-FK)*CJ(15)*WIJ*DCMPLX(DNZ(IP),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            II = 3*ID - (3 - JW)
            DO IQ = 1, IX
              JJ = 3*INDX(IQ) - 1
              AIJ = -DCMPLX(0D0,-FK)*CJ(16)*WIJ*DCMPLX(DNX(IQ),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            DO IQ = 1, IZ
              JJ = 3*INDZ(IQ) - 1
              AIJ = -DCMPLX(0D0,-FK)*CJ(17)*WIJ*DCMPLX(DNZ(IQ),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            JJ = 3*ID - 1
            IF (JJ == II) THEN
              AIJ = WIJ*DCMPLX(FK*FK,0D0)*CJ(18)
              CALL ADD_ENTRY(JJ, JJ, AIJ)
            END IF

            IF (JW == 2) THEN
              AIJ = -WIJ*RHO*OMIGA2
              CALL ADD_ENTRY(JJ, JJ, AIJ)
            END IF

            ! ---- C-block (Gz) ----
            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              DO IQ = 1, IX
                JJ = 3*INDX(IQ)
                AIJ = WIJ*CJ(19)*DCMPLX(DNX(IP)*DNX(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              DO IQ = 1, IZ
                JJ = 3*INDZ(IQ)
                AIJ = WIJ*CJ(20)*DCMPLX(DNX(IP)*DNZ(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              DO IQ = 1, IX
                JJ = 3*INDX(IQ)
                AIJ = WIJ*CJ(21)*DCMPLX(DNZ(IP)*DNX(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              DO IQ = 1, IZ
                JJ = 3*INDZ(IQ)
                AIJ = WIJ*CJ(22)*DCMPLX(DNZ(IP)*DNZ(IQ),0D0)
                CALL ADD_ENTRY(II, JJ, AIJ)
              END DO
            END DO

            JJ = 3*ID
            DO IP = 1, IX
              II = 3*INDX(IP) - (3 - JW)
              AIJ = WIJ*DCMPLX(0D0,-FK)*CJ(23)*DCMPLX(DNX(IP),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            DO IP = 1, IZ
              II = 3*INDZ(IP) - (3 - JW)
              AIJ = WIJ*DCMPLX(0D0,-FK)*CJ(24)*DCMPLX(DNZ(IP),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            II = 3*ID - (3 - JW)
            DO IQ = 1, IX
              JJ = 3*INDX(IQ)
              AIJ = -WIJ*DCMPLX(0D0,-FK)*CJ(25)*DCMPLX(DNX(IQ),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            DO IQ = 1, IZ
              JJ = 3*INDZ(IQ)
              AIJ = -WIJ*DCMPLX(0D0,-FK)*CJ(26)*DCMPLX(DNZ(IQ),0D0)
              CALL ADD_ENTRY(II, JJ, AIJ)
            END DO

            JJ = 3*ID
            IF (JJ == II) THEN
              AIJ = WIJ*CJ(27)*DCMPLX(FK*FK,0D0)
              CALL ADD_ENTRY(JJ, JJ, AIJ)
            END IF

            IF (JW == 3) THEN
              AIJ = -WIJ*RHO*OMIGA2
              CALL ADD_ENTRY(JJ, JJ, AIJ)
            END IF

          END DO  ! JW
        END DO    ! L
      END DO      ! K
    END DO        ! J
  END DO          ! I
      DEALLOCATE (Z1, Z2, T1, T2, DLX, DLZ, DNX, DNZ)
      DEALLOCATE (XP, P, NOX, NOZ, INDX, INDZ, Q, CJ)
END SUBROUTINE CA_core_ori

SUBROUTINE CAGB_ori( FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                 CI, NORD, AS, WT, IE0, IS0, DZ0, N, KL, KU, LDAB, AGB )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, KL, KU, LDAB
    INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ

      DOUBLE PRECISION, INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
      DOUBLE PRECISION, INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0
      INTEGER :: I, J, K, L, N0, MM, I0, II, KK
      COMPLEX*16 :: OMIGA2, RHO, WIJ, AIJ
      COMPLEX*16, INTENT(INOUT) :: AGB(LDAB, N)

  AGB(:,:) = (0.D0,0.D0)

  CALL CA_core_ori( FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                 CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU,   &   
            ADD_ENTRY )

    CONTAINS
    SUBROUTINE ADD_ENTRY(row, col, val)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: row, col
      COMPLEX*16, INTENT(IN) :: val
      INTEGER :: ib
      IF (val==(0D0,0D0)) RETURN
      ib = KL + KU + 1 + row - col
      IF (ib>=1 .AND. ib<=LDAB) AGB(ib, col) = AGB(ib, col) + val
    END SUBROUTINE ADD_ENTRY

END SUBROUTINE CAGB_ori

   SUBROUTINE CAGB(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                   CI, NORD, AS, WT, IE0, IS0, DZ0, N, KL, KU, LDAB, IG, AGB)
      USE shared_mod, ONLY: dp
      USE err_mpi_mod, ONLY: fail
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, KL, KU, LDAB
      INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ
      TYPE(InversionGridType), INTENT(IN) :: IG
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
      REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0
      COMPLEX(dp), INTENT(INOUT) :: AGB(LDAB, N)

      REAL(dp), PARAMETER :: ZERO_TOL = 0.0_dp !ZERO_TOL = 1.0e-30_dp

      AGB(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)

      CALL CA_core(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                   CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)

   CONTAINS
      SUBROUTINE ADD_ENTRY(row, col, val)
         USE shared_mod, ONLY: dp
         IMPLICIT NONE
         INTEGER, INTENT(IN)     :: row, col
         COMPLEX(dp), INTENT(IN) :: val
         INTEGER :: ib
         IF (ABS(val) < ZERO_TOL) RETURN
         ib = KL + KU + 1 + row - col
         IF (ib < 1 .OR. ib > LDAB) RETURN              ! out of band → ignore
         IF (col < 1 .OR. col > N .OR. row < 1 .OR. row > N) THEN
            CALL fail('CAGB: ADD_ENTRY indices out of range')
         END IF
         AGB(ib, col) = AGB(ib, col) + val
      END SUBROUTINE ADD_ENTRY
   END SUBROUTINE CAGB

   !============================= COO (MUMPS/future) ======================
   SUBROUTINE CA_COO(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                     CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IRN, JCN, ACOO, NNZ_OUT, IG)
      USE shared_mod, ONLY: dp
      USE err_mpi_mod, ONLY: fail
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, KL, KU, LDAB
      INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ
      TYPE(InversionGridType), INTENT(IN) :: IG
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
      REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0

      INTEGER, ALLOCATABLE, INTENT(OUT) :: IRN(:), JCN(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: ACOO(:)
      INTEGER, INTENT(OUT) :: NNZ_OUT

      INTEGER, ALLOCATABLE :: irn_buf(:), jcn_buf(:)
      COMPLEX(dp), ALLOCATABLE :: acoo_buf(:)
      INTEGER :: cap, nnz0, ierr
      REAL(dp), PARAMETER :: ZERO_TOL = 1.0e-30_dp

      cap = MAX(200000, 10*N)   ! heuristic; adjust if needed
      nnz0 = 0
      ALLOCATE (irn_buf(cap), jcn_buf(cap), acoo_buf(cap), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO: initial buffer alloc failed')

      CALL CA_core(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                   CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)

      NNZ_OUT = nnz0
      ALLOCATE (IRN(MAX(1, nnz0)), JCN(MAX(1, nnz0)), ACOO(MAX(1, nnz0)), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO: final alloc failed')

      IF (nnz0 > 0) THEN
         IRN = irn_buf(1:nnz0)
         JCN = jcn_buf(1:nnz0)
         ACOO = acoo_buf(1:nnz0)
      END IF

      DEALLOCATE (irn_buf, jcn_buf, acoo_buf)
      RETURN

   CONTAINS
      SUBROUTINE GROW()
         INTEGER :: newcap, ierr
         INTEGER, ALLOCATABLE :: i2(:), j2(:)
         COMPLEX(dp), ALLOCATABLE :: a2(:)
         newcap = INT(1.5_dp*cap) + 100000
         ALLOCATE (i2(newcap), j2(newcap), a2(newcap), STAT=ierr)
         IF (ierr /= 0) CALL fail('CA_COO: GROW alloc failed')
         IF (nnz0 > 0) THEN
            i2(1:nnz0) = irn_buf(1:nnz0)
            j2(1:nnz0) = jcn_buf(1:nnz0)
            a2(1:nnz0) = acoo_buf(1:nnz0)
         END IF
         CALL MOVE_ALLOC(i2, irn_buf)
         CALL MOVE_ALLOC(j2, jcn_buf)
         CALL MOVE_ALLOC(a2, acoo_buf)
         cap = newcap
      END SUBROUTINE GROW

      SUBROUTINE ADD_ENTRY(row, col, val)
         USE shared_mod, ONLY: dp
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: row, col
         COMPLEX(dp), INTENT(IN) :: val
         IF (ABS(val) < ZERO_TOL) RETURN
         IF (row < 1 .OR. row > N .OR. col < 1 .OR. col > N) THEN
            CALL fail('CA_COO: ADD_ENTRY indices out of range')
         END IF
         nnz0 = nnz0 + 1
         IF (nnz0 > cap) CALL GROW()
         irn_buf(nnz0) = row
         jcn_buf(nnz0) = col
         acoo_buf(nnz0) = val
      END SUBROUTINE ADD_ENTRY
   END SUBROUTINE CA_COO

   !============================= CSR (PARDISO) ===========================
   SUBROUTINE CA_CSR(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                     AS, WT, IE0, IS0, DZ0, N, IA, JA, AV, NBLOCK, LDAB, KL, KU, IG)
      USE shared_mod, ONLY: dp
      USE err_mpi_mod, ONLY: fail
      IMPLICIT NONE
      ! ---- inputs ----
      INTEGER, INTENT(IN) :: N, NTO, NX, NZ, NNX, NNZ, NORD, IANISO, IE0, IS0, NBLOCK, LDAB, KL, KU
      REAL(dp), INTENT(IN) :: FREQ, FK, DZ0, XTO(:), ZTO(:), X(:), AS(:), WT(:)
      REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :)
      TYPE(InversionGridType), INTENT(IN) :: IG
      ! ---- outputs (1-based CSR) ----
      INTEGER, ALLOCATABLE, INTENT(OUT) :: IA(:), JA(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: AV(:)

      ! ---- locals ----
      INTEGER, ALLOCATABLE :: BUCKET(:), IA0(:), POS(:), JV_BUF(:)
      COMPLEX(dp), ALLOCATABLE :: AV_BUF(:)
      INTEGER :: i, nnz_raw, nnz_final, ierr
      LOGICAL :: is_count_pass
      REAL(dp), PARAMETER :: ZERO_TOL = 1.0e-30_dp

      ! -------- PASS 1: count per-row entries --------
      ALLOCATE (BUCKET(N), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc BUCKET')
      BUCKET = 0
      is_count_pass = .TRUE.

      CALL CA_core(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                   AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)

      ! prefix sum (raw IA0)
      ALLOCATE (IA0(N + 1), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc IA0')
      IA0(1) = 1
      DO i = 1, N
         IA0(i + 1) = IA0(i) + BUCKET(i)
      END DO
      nnz_raw = IA0(N + 1) - 1
      IF (nnz_raw < 0) nnz_raw = 0

      ! -------- PASS 2: fill JV/AV with in-row dedup --------
      ALLOCATE (JV_BUF(MAX(1, nnz_raw)), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc JV_BUF')
      ALLOCATE (AV_BUF(MAX(1, nnz_raw)), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc AV_BUF')
      JV_BUF = 0
      AV_BUF = CMPLX(0.0_dp, 0.0_dp, dp)

      ALLOCATE (POS(N), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc POS')
      POS = IA0(1:N)

      is_count_pass = .FALSE.

      CALL CA_core(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                   AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)

      ! compute final row sizes
      DO i = 1, N
         BUCKET(i) = POS(i) - IA0(i)
      END DO

      ! build final IA
      ALLOCATE (IA(N + 1), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc IA')
      IA(1) = 1
      DO i = 1, N
         IA(i + 1) = IA(i) + BUCKET(i)
      END DO

      nnz_final = IA(N + 1) - 1
      ALLOCATE (JA(MAX(1, nnz_final)), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc JA')
      ALLOCATE (AV(MAX(1, nnz_final)), STAT=ierr); IF (ierr /= 0) CALL fail('CA_CSR: alloc AV')

      ! pack each row and sort
      DO i = 1, N
         IF (BUCKET(i) > 0) THEN
            JA(IA(i):IA(i + 1) - 1) = JV_BUF(IA0(i):IA0(i) + BUCKET(i) - 1)
            AV(IA(i):IA(i + 1) - 1) = AV_BUF(IA0(i):IA0(i) + BUCKET(i) - 1)
            CALL SORT_ROW(JA(IA(i):IA(i + 1) - 1), AV(IA(i):IA(i + 1) - 1))
         END IF
      END DO

      DEALLOCATE (BUCKET, IA0, POS, JV_BUF, AV_BUF)
      RETURN

   CONTAINS
      SUBROUTINE ADD_ENTRY(row, col, val)
         USE shared_mod, ONLY: dp
         IMPLICIT NONE
         INTEGER, INTENT(IN)     :: row, col
         COMPLEX(dp), INTENT(IN) :: val
         INTEGER :: p, k_left, k_right
         LOGICAL :: found

         IF (ABS(val) < ZERO_TOL) RETURN
         IF (row < 1 .OR. row > N .OR. col < 1 .OR. col > N) THEN
            CALL fail('CA_CSR: ADD_ENTRY indices out of range')
         END IF

         ! optional band mask
         IF ((row - col) > KL .OR. (col - row) > KU) RETURN

         IF (is_count_pass) THEN
            BUCKET(row) = BUCKET(row) + 1
         ELSE
            found = .FALSE.
            k_left = IA0(row)
            k_right = POS(row) - 1
            DO p = k_left, k_right
               IF (JV_BUF(p) == col) THEN
                  AV_BUF(p) = AV_BUF(p) + val
                  found = .TRUE.; EXIT
               END IF
            END DO
            IF (.NOT. found) THEN
               IF (POS(row) > SIZE(JV_BUF)) CALL fail('CA_CSR: buffer overflow')
               JV_BUF(POS(row)) = col
               AV_BUF(POS(row)) = val
               POS(row) = POS(row) + 1
            END IF
         END IF
      END SUBROUTINE ADD_ENTRY

      SUBROUTINE SORT_ROW(jv, zv)
         IMPLICIT NONE
         INTEGER, INTENT(INOUT) :: jv(:)
         COMPLEX(dp), INTENT(INOUT) :: zv(:)
         INTEGER :: lb, ub, idx, j, keyj
         COMPLEX(dp) :: keyz
         lb = LBOUND(jv, 1); ub = UBOUND(jv, 1)
         IF (ub - lb <= 0) RETURN
         DO idx = lb + 1, ub
            keyj = jv(idx); keyz = zv(idx); j = idx - 1
            DO WHILE (j >= lb)
               IF (jv(j) <= keyj) EXIT
               jv(j + 1) = jv(j)
               zv(j + 1) = zv(j)
               j = j - 1
            END DO
            jv(j + 1) = keyj
            zv(j + 1) = keyz
         END DO
      END SUBROUTINE SORT_ROW
   END SUBROUTINE CA_CSR

   SUBROUTINE CA_core_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                   CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, &
                   ADD_ENTRY)   ! callback for storage

   USE omp_lib
   IMPLICIT NONE
   ! ==== arguments ====
   INTEGER, INTENT(IN) :: N, KL, KU, LDAB
   INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ
   REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
   REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0
   TYPE(InversionGridType), INTENT(IN) :: IG

   ! ==== allocatable temporaries (will be allocated per-thread) ====
   REAL(dp), ALLOCATABLE :: Z1(:), Z2(:), T1(:), T2(:)
   REAL(dp), ALLOCATABLE :: NOX(:), NOZ(:), DLX(:), DLZ(:)
   REAL(dp), ALLOCATABLE :: DNX(:), DNZ(:), XP(:)
   INTEGER,  ALLOCATABLE :: INDX(:), INDZ(:)
   COMPLEX(dp), ALLOCATABLE :: P(:), Q(:), CJ(:)

   ! ==== scalars ====
   INTEGER :: I, J, K, L, N0, MM, I0, NBLOCK, IM
   INTEGER :: IX, IZ, ID, K1, L1, IP, IQ, JW, II, JJ
   REAL(dp) :: Z0, OMIGA, DX, C1, AK, DZ, C2, C3, BL
   COMPLEX(dp) :: OMIGA2, RHO, WIJ, AIJ

   ! Dummy procedure interface for the callback
   INTERFACE
      SUBROUTINE ADD_ENTRY(row, col, val)
         USE iso_fortran_env, ONLY: dp => real64
         IMPLICIT NONE
         INTEGER, INTENT(IN)     :: row, col
         COMPLEX(dp), INTENT(IN) :: val
      END SUBROUTINE
   END INTERFACE

   ! ---- initialisation (shared constants) ----
   MM    = 2*(NORD - 1) + 1
   I0    = KL + KU + 1
   Z0    = IE0*DZ0
   OMIGA = 2.0_dp*3.1415926_dp*FREQ
   OMIGA2 = CMPLX(OMIGA*OMIGA, 0.0_dp, dp)

   NBLOCK = SIZE(IG%DX_BLOCK)
   ! write(*,*) 'CA_core_par nthreads = ', omp_get_max_threads()

!$omp parallel default(shared) &
!$omp private( I,J,K,L,K1,L1,IM,N0,DX,C1,AK,DZ,C2,C3,BL,ID,IX,IZ,IP,IQ,JW, &
!$omp          II,JJ,WIJ,RHO,AIJ,                                     &
!$omp          Z1,Z2,T1,T2,NOX,NOZ,DLX,DLZ,DNX,DNZ,XP,INDX,INDZ,P,Q,CJ )

      ! per-thread allocations
      ALLOCATE (Z1(NORD), Z2(NORD), T1(NORD), T2(NORD))
      ALLOCATE (NOX(NORD), NOZ(NORD), INDZ(NORD))
      ALLOCATE (DLX(NORD), DLZ(NORD), DNZ(NORD))
      ALLOCATE (INDX(MM), DNX(MM))
      ALLOCATE (XP(NORD), P(21), Q(81), CJ(27))


!$omp do collapse(2) schedule(static)
      DO I = 1, NX - 1
         DO J = 1, NZ - 1

            ! ----------------------------------------
            ! Block index in IG (same ordering as before)
            ! ----------------------------------------
            IM = (I - 1)*(NZ - 1) + J

            ! --- geometry for this block from IG ---
            DX = IG%DX_BLOCK(IM)
            C1 = IG%C1_BLOCK(IM)
            N0 = IG%N0_BLOCK(IM)

            DO K = 1, NORD
               Z1(K) = IG%Z1_OUT(K, IM)
               Z2(K) = IG%Z2_OUT(K, IM)
               T1(K) = IG%T1_OUT(K, IM)
               T2(K) = IG%T2_OUT(K, IM)
            END DO

            ! ---- GQ over (K,L) ----
            DO K = 1, NORD
               AK = AS(K)
               DZ = Z2(K) - Z1(K)
               C3 = 2.0_dp / DZ
               CALL CDLI(AK, NORD, AS, DLX)

               DO L = 1, NORD
                  ID = N0 + (K - 1)*NNZ + (L - 1)

                  CALL C21(ID, IANISO, CR, CI, RHO, P)
                  CALL Q81_NewGSRM(FREQ, NX - 1, NZ - 1, I, J, K, L, NORD, IE0, IS0, RHO, P, Q)

                  WIJ = CMPLX(0.25_dp*DX*DZ*WT(K)*WT(L), 0.0_dp, dp)
                  BL  = AS(L)
                  C2  = -((T2(K) - T1(K))*BL + (T1(K) + T2(K))) / DZ
                  CALL CDLI(BL, NORD, AS, DLZ)

                  ! build INDX/DNX and INDZ/DNZ
                  DO K1 = 1, NORD
                     NOX(K1) = (N0 + (L - 1)) + (K1 - 1)*NNZ
                     NOZ(K1) = (N0 + (K - 1)*NNZ) + (K1 - 1)
                  END DO

                  IX = 0
                  DO K1 = 1, K - 1
                     IX = IX + 1
                     INDX(IX) = NOX(K1)
                     DNX(IX)  = C1*DLX(K1)
                  END DO
                  DO L1 = 1, L - 1
                     IX = IX + 1
                     INDX(IX) = NOZ(L1)
                     DNX(IX)  = C2*DLZ(L1)
                  END DO
                  IX = IX + 1
                  INDX(IX) = NOX(K)                ! NOX(K) = NOZ(L)
                  DNX(IX)  = C1*DLX(K) + C2*DLZ(L)
                  DO L1 = L + 1, NORD
                     IX = IX + 1
                     INDX(IX) = NOZ(L1)
                     DNX(IX)  = C2*DLZ(L1)
                  END DO
                  DO K1 = K + 1, NORD
                     IX = IX + 1
                     INDX(IX) = NOX(K1)
                     DNX(IX)  = C1*DLX(K1)
                  END DO

                  DO L1 = 1, NORD
                     INDZ(L1) = NOZ(L1)
                     DNZ(L1)  = C3*DLZ(L1)
                  END DO
                  IZ = NORD

                  ! ! sanity (your checks)
                  ! IF (NNX*NNZ > 0) THEN
                  !    DO K1 = 1, IX
                  !       IF (INDX(K1) < 1 .OR. INDX(K1) > N) THEN
                  !          WRITE (*, '(A,3I12)') 'CA: bad INDX', K1, INDX(K1), N
                  !          STOP
                  !       END IF
                  !    END DO
                  !    DO L1 = 1, IZ
                  !       IF (INDZ(L1) < 1 .OR. INDZ(L1) > N) THEN
                  !          WRITE (*, '(A,3I12)') 'CA: bad INDZ', L1, INDZ(L1), N
                  !          STOP
                  !       END IF
                  !    END DO
                  ! END IF

                  ! ===== weights loop (unchanged math) =====
                  DO JW = 1, 3
                     CALL QCJ(JW, Q, CJ)

                     ! ---- A-block: c1j11..,c1j13..,c3j11..,c3j13.. ----
                     ! c1j11-term: (DxDx)Gx
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 2
                           AIJ = WIJ*CJ(1)*DCMPLX(DNX(IP)*DNX(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 2
                           AIJ = WIJ*CJ(2)*DCMPLX(DNX(IP)*DNZ(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 2
                           AIJ = WIJ*CJ(3)*DCMPLX(DNZ(IP)*DNX(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 2
                           AIJ = WIJ*CJ(4)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     JJ = 3*ID - 2
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = DCMPLX(0.0_dp, -FK)*CJ(5)*WIJ*DCMPLX(DNX(IP), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = DCMPLX(0.0_dp, -FK)*CJ(6)*WIJ*DCMPLX(DNZ(IP), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ) - 2
                        AIJ = -DCMPLX(0.0_dp, -FK)*WIJ*CJ(7)*DCMPLX(DNX(IQ), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ) - 2
                        AIJ = -DCMPLX(0.0_dp, -FK)*WIJ*CJ(8)*DCMPLX(DNZ(IQ), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     JJ = 3*ID - 2
                     IF (II == JJ) THEN
                        AIJ = DCMPLX(FK*FK, 0.0_dp)*WIJ*CJ(9)
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     IF (JW == 1) THEN
                        AIJ = -RHO*OMIGA2*WIJ
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     ! ---- B-block (Gy) ----
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 1
                           AIJ = WIJ*CJ(10)*DCMPLX(DNX(IP)*DNX(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 1
                           AIJ = WIJ*CJ(11)*DCMPLX(DNX(IP)*DNZ(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 1
                           AIJ = WIJ*CJ(12)*DCMPLX(DNZ(IP)*DNX(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 1
                           AIJ = WIJ*CJ(13)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     JJ = 3*ID - 1
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = DCMPLX(0.0_dp, -FK)*CJ(14)*WIJ*DCMPLX(DNX(IP), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = DCMPLX(0.0_dp, -FK)*CJ(15)*WIJ*DCMPLX(DNZ(IP), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ) - 1
                        AIJ = -DCMPLX(0.0_dp, -FK)*CJ(16)*WIJ*DCMPLX(DNX(IQ), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ) - 1
                        AIJ = -DCMPLX(0.0_dp, -FK)*CJ(17)*WIJ*DCMPLX(DNZ(IQ), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     JJ = 3*ID - 1
                     IF (JJ == II) THEN
                        AIJ = WIJ*DCMPLX(FK*FK, 0.0_dp)*CJ(18)
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     IF (JW == 2) THEN
                        AIJ = -WIJ*RHO*OMIGA2
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     ! ---- C-block (Gz) ----
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ)
                           AIJ = WIJ*CJ(19)*DCMPLX(DNX(IP)*DNX(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ)
                           AIJ = WIJ*CJ(20)*DCMPLX(DNX(IP)*DNZ(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ)
                           AIJ = WIJ*CJ(21)*DCMPLX(DNZ(IP)*DNX(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ)
                           AIJ = WIJ*CJ(22)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.0_dp)
                           CALL ADD_ENTRY(II, JJ, AIJ)
                        END DO
                     END DO

                     JJ = 3*ID
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = WIJ*DCMPLX(0.0_dp, -FK)*CJ(23)*DCMPLX(DNX(IP), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = WIJ*DCMPLX(0.0_dp, -FK)*CJ(24)*DCMPLX(DNZ(IP), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ)
                        AIJ = -WIJ*DCMPLX(0.0_dp, -FK)*CJ(25)*DCMPLX(DNX(IQ), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ)
                        AIJ = -WIJ*DCMPLX(0.0_dp, -FK)*CJ(26)*DCMPLX(DNZ(IQ), 0.0_dp)
                        CALL ADD_ENTRY(II, JJ, AIJ)
                     END DO

                     JJ = 3*ID
                     IF (JJ == II) THEN
                        AIJ = WIJ*CJ(27)*DCMPLX(FK*FK, 0.0_dp)
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                     IF (JW == 3) THEN
                        AIJ = -WIJ*RHO*OMIGA2
                        CALL ADD_ENTRY(JJ, JJ, AIJ)
                     END IF

                  END DO  ! JW
               END DO    ! L
            END DO      ! K

         END DO        ! J
      END DO          ! I
!$omp end do

      ! per-thread cleanup
      DEALLOCATE (Z1, Z2, T1, T2, DLX, DLZ, DNX, DNZ)
      DEALLOCATE (XP, P, NOX, NOZ, INDX, INDZ, Q, CJ)

!$omp end parallel

END SUBROUTINE CA_core_par

! SUBROUTINE CA_CSR_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
!                       AS, WT, IE0, IS0, DZ0, N, IA, JA, AV, NBLOCK, LDAB, KL, KU, IG)
!    USE shared_mod,  ONLY: dp
!    USE err_mpi_mod, ONLY: fail
!    USE omp_lib
!    IMPLICIT NONE
!    ! ---- inputs ----
!    INTEGER, INTENT(IN) :: N, NTO, NX, NZ, NNX, NNZ, NORD, IANISO, IE0, IS0, NBLOCK, LDAB, KL, KU
!    REAL(dp), INTENT(IN) :: FREQ, FK, DZ0, XTO(:), ZTO(:), X(:), AS(:), WT(:)
!    REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :)
!    TYPE(InversionGridType), INTENT(IN) :: IG
!    ! ---- outputs (1-based CSR) ----
!    INTEGER,          ALLOCATABLE, INTENT(OUT) :: IA(:), JA(:)
!    COMPLEX(dp),      ALLOCATABLE, INTENT(OUT) :: AV(:)

!    ! ---- per-thread COO buffers ----
!    TYPE ThreadBuf
!       INTEGER,     ALLOCATABLE :: row(:), col(:)
!       COMPLEX(dp), ALLOCATABLE :: val(:)
!       INTEGER,     ALLOCATABLE :: row_count(:)  ! per-row nnz per thread
!       INTEGER                   :: nnz_used
!       INTEGER                   :: cap
!    END TYPE ThreadBuf

!    TYPE(ThreadBuf), ALLOCATABLE :: buf(:)

!    ! ---- global CSR helpers ----
!    INTEGER,     ALLOCATABLE :: row_nnz(:)                ! total nnz per row (sum over threads)

!    INTEGER,     ALLOCATABLE :: BUCKET(:)                 ! final nnz per row after dedup
!    INTEGER,     ALLOCATABLE :: IA_new(:), JA_new(:)
!    COMPLEX(dp), ALLOCATABLE :: AV_new(:)
!    INTEGER, ALLOCATABLE :: row_off(:)      ! size N, running offset inside each row
! INTEGER, ALLOCATABLE :: t_row_start(:)  ! size N, start offset for current thread t
! INTEGER, ALLOCATABLE :: row_next(:)     ! size N, cursor for current thread t


!    INTEGER :: ierr, nthreads, t, init_cap
!    INTEGER :: nnz_total, nnz_final
!    INTEGER :: i, k, r, c, offset
!    INTEGER :: lb, ub, w, kk, cur_col
!    COMPLEX(dp) :: cur_val
!    REAL(dp) :: t_core_start, t_core_end
!    REAL(dp), PARAMETER :: ZERO_TOL = 1.0e-30_dp

!    ! ------------------------------------------------------------------
!    ! 1. Determine number of threads and allocate per-thread buffers
!    ! ------------------------------------------------------------------

! #ifdef _OPENMP
!    nthreads = omp_get_max_threads()
! #else
!    nthreads = 1
! #endif
!    IF (nthreads <= 0) nthreads = 1
! !  WRITE(*,*) 'CA_CSR: nthreads =', nthreads
!    ALLOCATE(buf(nthreads), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc buf(:) failed')

!    ! heuristic initial capacity per thread (can be tuned)
!    init_cap = MAX(200000, 10*N / MAX(1, nthreads))

!    DO t = 1, nthreads
!       buf(t)%cap      = init_cap
!       buf(t)%nnz_used = 0
!       ALLOCATE(buf(t)%row(buf(t)%cap), buf(t)%col(buf(t)%cap), buf(t)%val(buf(t)%cap), STAT=ierr)
!       IF (ierr /= 0) CALL fail('CA_CSR_par: per-thread buffer alloc failed')
!       ALLOCATE(buf(t)%row_count(N), STAT=ierr)
!       IF (ierr /= 0) CALL fail('CA_CSR_par: per-thread row_count alloc failed')
!       buf(t)%row_count = 0
!    END DO

!    ! ------------------------------------------------------------------
!    ! 2. Run CA_core_par; ADD_ENTRY fills per-thread COO buffers + histograms
!    ! ------------------------------------------------------------------

!    CALL CA_core_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
!                     AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)


!    ! ------------------------------------------------------------------
!    ! 3. Build global row_nnz and IA (for raw entries)
!    ! ------------------------------------------------------------------
!    ALLOCATE(row_nnz(N), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc row_nnz failed')
!    row_nnz = 0

!    DO t = 1, nthreads
!       DO i = 1, N
!          row_nnz(i) = row_nnz(i) + buf(t)%row_count(i)
!       END DO
!    END DO

!    ALLOCATE(IA(N+1), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc IA failed')

!    IA(1) = 1
!    DO i = 1, N
!       IA(i+1) = IA(i) + row_nnz(i)
!    END DO

!    nnz_total = IA(N+1) - 1
!    IF (nnz_total < 0) nnz_total = 0

!    ALLOCATE(JA(MAX(1, nnz_total)), AV(MAX(1, nnz_total)), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc JA/AV failed')
!    ! ------------------------------------------------------------------
!    ! 4-5. Compute per-thread row offsets WITHOUT row_start/row_next 2D arrays
!    !      and scatter per-thread COO into global CSR (unsorted within row)
!    ! ------------------------------------------------------------------

!    ALLOCATE(row_off(N), t_row_start(N), row_next(N), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc row_off/t_row_start/row_next failed')

!    ! row_off(i) accumulates how many slots in row i have been reserved by previous threads
!    row_off = 0

!    DO t = 1, nthreads

!       ! For current thread t, compute its start offset in each row and init cursor
!       DO i = 1, N
!          t_row_start(i) = IA(i) + row_off(i)
!          row_next(i)    = t_row_start(i)
!          row_off(i)     = row_off(i) + buf(t)%row_count(i)
!       END DO

!       ! Scatter thread t COO entries into the slice assigned to thread t within each row
!       DO k = 1, buf(t)%nnz_used
!          r = buf(t)%row(k)
!          c = buf(t)%col(k)

!          ! safety: row range
!          IF (r < 1 .OR. r > N) THEN
!             CALL fail('CA_CSR_par: row index out of range in scatter')
!          END IF

!          offset = row_next(r)

!          ! offset must remain inside this thread's reserved slice for row r:
!          ! [t_row_start(r), t_row_start(r) + row_count(r) - 1]
!          IF (offset < t_row_start(r) .OR. offset >= t_row_start(r) + buf(t)%row_count(r)) THEN
!             CALL fail('CA_CSR_par: offset out of thread row slice')
!          END IF

!          JA(offset) = c
!          AV(offset) = buf(t)%val(k)
!          row_next(r) = offset + 1
!       END DO

!    END DO  ! t loop

!    DEALLOCATE(row_off, t_row_start, row_next)

!    ! ! ------------------------------------------------------------------
!    ! ! 4. Compute per-thread row_start/row_next offsets in JA/AV
!    ! ! ------------------------------------------------------------------
!    ! ALLOCATE(row_start(nthreads, N), row_next(nthreads, N), STAT=ierr)
!    ! IF (ierr /= 0) CALL fail('CA_CSR_par: alloc row_start/row_next failed')

!    ! DO i = 1, N
!    !    offset = IA(i)
!    !    DO t = 1, nthreads
!    !       row_start(t,i) = offset
!    !       row_next(t,i)  = offset
!    !       offset = offset + buf(t)%row_count(i)
!    !    END DO
!    ! END DO

!    ! ! ------------------------------------------------------------------
!    ! ! 5. Scatter per-thread COO into global CSR (unsorted within row)
!    ! ! ------------------------------------------------------------------
!    ! DO t = 1, nthreads
!    !    DO k = 1, buf(t)%nnz_used
!    !       r = buf(t)%row(k)
!    !       c = buf(t)%col(k)

!    !       ! safety: row range
!    !       IF (r < 1 .OR. r > N) THEN
!    !          CALL fail('CA_CSR_par: row index out of range in scatter')
!    !       END IF

!    !       offset = row_next(t, r)
!    !       IF (offset < IA(r) .OR. offset > IA(r+1)) THEN
!    !          CALL fail('CA_CSR_par: offset out of row bounds')
!    !       END IF

!    !       JA(offset) = c
!    !       AV(offset) = buf(t)%val(k)
!    !       row_next(t, r) = offset + 1
!    !    END DO
!    ! END DO

!    ! ------------------------------------------------------------------
!    ! 6. Per-row sort, deduplicate (sum duplicates) WITHOUT zero-drop
!    ! ------------------------------------------------------------------
!    ALLOCATE(BUCKET(N), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc BUCKET failed')
!    BUCKET = 0
!    t_core_start = omp_get_wtime()
!    !$omp parallel do default(shared) private(i,lb,ub,w,kk,cur_col,cur_val) schedule(static)
!    DO i = 1, N
!       lb = IA(i)
!       ub = IA(i+1) - 1
!       IF (ub < lb) CYCLE

!       ! sort row by column index
!       CALL SORT_ROW_SEGMENT(JA, AV, lb, ub)

!       ! dedup only (keep zeros) in-place in [lb,ub]
!       w = lb
!       cur_col = JA(lb)
!       cur_val = AV(lb)

!       DO kk = lb + 1, ub
!          IF (JA(kk) == cur_col) THEN
!             cur_val = cur_val + AV(kk)
!          ELSE
!             JA(w) = cur_col
!             AV(w) = cur_val
!             w = w + 1

!             cur_col = JA(kk)
!             cur_val = AV(kk)
!          END IF
!       END DO

!       ! flush last accumulated
!       JA(w) = cur_col
!       AV(w) = cur_val
!       w = w + 1

!       BUCKET(i) = w - lb
!    END DO
! !$omp end parallel do
!    ! ------------------------------------------------------------------
!    ! 7. Build final IA_new and compact JA_new/AV_new to exact nnz_final
!    ! ------------------------------------------------------------------
!    ALLOCATE(IA_new(N+1), STAT=ierr)
!    IF (ierr /= 0) CALL fail('CA_CSR_par: alloc IA_new failed')
!    IA_new(1) = 1
!    DO i = 1, N
!       IA_new(i+1) = IA_new(i) + BUCKET(i)
!    END DO

!    nnz_final = IA_new(N+1) - 1
!    IF (nnz_final < 0) nnz_final = 0

!    IF (nnz_final <= 0) THEN
!       ! empty matrix case
!       DEALLOCATE(JA, AV)
!       ALLOCATE(JA(1), AV(1), STAT=ierr)
!       IF (ierr /= 0) CALL fail('CA_CSR_par: alloc empty JA/AV failed')
!       JA(1) = 0
!       AV(1) = CMPLX(0.0_dp, 0.0_dp, dp)

!       DEALLOCATE(IA)
!       CALL MOVE_ALLOC(IA_new, IA)
!    ELSE
!       ALLOCATE(JA_new(nnz_final), AV_new(nnz_final), STAT=ierr)
!       IF (ierr /= 0) CALL fail('CA_CSR_par: alloc JA_new/AV_new failed')

!       DO i = 1, N
!          lb = IA(i)
!          ub = lb + BUCKET(i) - 1
!          IF (BUCKET(i) > 0) THEN
!             JA_new(IA_new(i):IA_new(i+1)-1) = JA(lb:ub)
!             AV_new(IA_new(i):IA_new(i+1)-1) = AV(lb:ub)
!          END IF
!       END DO

!       DEALLOCATE(JA, AV)
!       DEALLOCATE(IA)
!       CALL MOVE_ALLOC(IA_new, IA)
!       CALL MOVE_ALLOC(JA_new, JA)
!       CALL MOVE_ALLOC(AV_new, AV)
!    END IF
!    t_core_end = omp_get_wtime()
!    ! WRITE(*,*) 'CA_CSR_par: CA_core_par elapsed (s)=', t_core_end - t_core_start
!    ! ------------------------------------------------------------------
!    ! 8. Cleanup per-thread buffers and helpers
!    ! ------------------------------------------------------------------
!    DO t = 1, nthreads
!       IF (ALLOCATED(buf(t)%row))       DEALLOCATE(buf(t)%row)
!       IF (ALLOCATED(buf(t)%col))       DEALLOCATE(buf(t)%col)
!       IF (ALLOCATED(buf(t)%val))       DEALLOCATE(buf(t)%val)
!       IF (ALLOCATED(buf(t)%row_count)) DEALLOCATE(buf(t)%row_count)
!    END DO
!    DEALLOCATE(buf)
!    DEALLOCATE(row_nnz, BUCKET) !row_next,

!    RETURN

! CONTAINS

!    SUBROUTINE GROW_LOCAL(tid)
!       INTEGER, INTENT(IN) :: tid
!       INTEGER :: newcap, ierr
!       INTEGER,     ALLOCATABLE :: r2(:), c2(:)
!       COMPLEX(dp), ALLOCATABLE :: v2(:)

!       newcap = INT(1.5_dp*buf(tid)%cap) + 100000

!       ALLOCATE(r2(newcap), c2(newcap), v2(newcap), STAT=ierr)
!       IF (ierr /= 0) CALL fail('CA_CSR_par: GROW_LOCAL alloc failed')

!       IF (buf(tid)%nnz_used > 0) THEN
!          r2(1:buf(tid)%nnz_used) = buf(tid)%row(1:buf(tid)%nnz_used)
!          c2(1:buf(tid)%nnz_used) = buf(tid)%col(1:buf(tid)%nnz_used)
!          v2(1:buf(tid)%nnz_used) = buf(tid)%val(1:buf(tid)%nnz_used)
!       END IF

! #ifdef _OPENMP
! !$omp critical (grow_local_critical)
! #endif
!       CALL MOVE_ALLOC(r2, buf(tid)%row)
!       CALL MOVE_ALLOC(c2, buf(tid)%col)
!       CALL MOVE_ALLOC(v2, buf(tid)%val)
!       buf(tid)%cap = newcap
! #ifdef _OPENMP
! !$omp end critical (grow_local_critical)
! #endif
!    END SUBROUTINE GROW_LOCAL

!    SUBROUTINE ADD_ENTRY(row, col, val)
!       USE shared_mod, ONLY: dp
!       IMPLICIT NONE
!       INTEGER,     INTENT(IN) :: row, col
!       COMPLEX(dp), INTENT(IN) :: val
!       INTEGER :: tid, pos

!       IF (ABS(val) < ZERO_TOL) RETURN
!       IF (row < 1 .OR. row > N .OR. col < 1 .OR. col > N) THEN
!          CALL fail('CA_CSR_par: ADD_ENTRY indices out of range')
!       END IF

!       ! band mask
!       IF ((row - col) > KL .OR. (col - row) > KU) RETURN

! #ifdef _OPENMP
!       tid = omp_get_thread_num() + 1
! #else
!       tid = 1
! #endif

!       pos = buf(tid)%nnz_used + 1
!       IF (pos > buf(tid)%cap) CALL GROW_LOCAL(tid)

!       buf(tid)%row(pos) = row
!       buf(tid)%col(pos) = col
!       buf(tid)%val(pos) = val
!       buf(tid)%nnz_used = pos

!       ! per-thread row histogram
!       buf(tid)%row_count(row) = buf(tid)%row_count(row) + 1
!    END SUBROUTINE ADD_ENTRY

!    ! Sort row segment JA(lb:ub) and permute AV consistently (insertion sort)
!    SUBROUTINE SORT_ROW_SEGMENT(JA, AV, lb, ub)
!       IMPLICIT NONE
!       INTEGER, INTENT(INOUT) :: JA(:)
!       COMPLEX(dp), INTENT(INOUT) :: AV(:)
!       INTEGER, INTENT(IN) :: lb, ub
!       INTEGER :: idx, j, keyj
!       COMPLEX(dp) :: keyz

!       IF (ub - lb <= 0) RETURN

!       DO idx = lb + 1, ub
!          keyj = JA(idx); keyz = AV(idx); j = idx - 1
!          DO WHILE (j >= lb)
!             IF (JA(j) <= keyj) EXIT
!             JA(j+1) = JA(j)
!             AV(j+1) = AV(j)
!             j = j - 1
!          END DO
!          JA(j+1) = keyj
!          AV(j+1) = keyz
!       END DO
!    END SUBROUTINE SORT_ROW_SEGMENT

! END SUBROUTINE CA_CSR_par

SUBROUTINE CA_CSR_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                      AS, WT, IE0, IS0, DZ0, N, IA, JA, AV, NBLOCK, LDAB, KL, KU, IG,my_rank, DEBUG_OUTPUT)
   USE shared_mod,  ONLY: dp
   USE err_mpi_mod, ONLY: fail
   USE omp_lib
   IMPLICIT NONE
   ! ---- inputs ----
   INTEGER, INTENT(IN) :: N, NTO, NX, NZ, NNX, NNZ, NORD, IANISO, IE0, IS0, NBLOCK, LDAB, KL, KU, my_rank
   REAL(dp), INTENT(IN) :: FREQ, FK, DZ0, XTO(:), ZTO(:), X(:), AS(:), WT(:)
   REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :)
   TYPE(InversionGridType), INTENT(IN) :: IG
   LOGICAL, INTENT(IN) :: DEBUG_OUTPUT
   ! ---- outputs (1-based CSR) ----
   INTEGER,          ALLOCATABLE, INTENT(OUT) :: IA(:), JA(:)
   COMPLEX(dp),      ALLOCATABLE, INTENT(OUT) :: AV(:)

   ! ---- per-thread COO buffers ----
   TYPE ThreadBuf
      INTEGER,     ALLOCATABLE :: row(:), col(:)
      COMPLEX(dp), ALLOCATABLE :: val(:)
      INTEGER,     ALLOCATABLE :: row_count(:)  ! per-row nnz per thread
      INTEGER                   :: nnz_used
      INTEGER                   :: cap
   END TYPE ThreadBuf

   TYPE(ThreadBuf), ALLOCATABLE :: buf(:)

   ! ---- global CSR helpers ----
   INTEGER,     ALLOCATABLE :: row_nnz(:)                ! total nnz per row (sum over threads)

   INTEGER,     ALLOCATABLE :: BUCKET(:)                 ! final nnz per row after dedup
   INTEGER,     ALLOCATABLE :: IA_new(:), JA_new(:)
   COMPLEX(dp), ALLOCATABLE :: AV_new(:)
   INTEGER, ALLOCATABLE :: row_off(:)      ! size N, running offset inside each row
   INTEGER, ALLOCATABLE :: t_row_start(:)  ! size N, start offset for current thread t
   INTEGER, ALLOCATABLE :: row_next(:)     ! size N, cursor for current thread t

   INTEGER :: ierr, nthreads, t, init_cap
   INTEGER :: nnz_total, nnz_final
   INTEGER :: i, k, r, c, offset
   INTEGER :: lb, ub, w, kk, cur_col
   COMPLEX(dp) :: cur_val
   REAL(dp) :: t_core_start, t_core_end

   ! --- DEBUG: CSR/memory accounting ---
   INTEGER :: SZ_INT_BYTES, SZ_CPLX_BYTES
   REAL(dp) :: bytes_A_raw, bytes_IA_raw, bytes_JA_raw, GB_raw
   REAL(dp) :: bytes_A_final, bytes_IA_final, bytes_JA_final, GB_final
   REAL(dp) :: bytes_coo_tot, GB_coo
   ! ------------------------------------

   REAL(dp), PARAMETER :: ZERO_TOL = 1.0e-30_dp

   ! ------------------------------------------------------------------
   ! 1. Determine number of threads and allocate per-thread buffers
   ! ------------------------------------------------------------------

#ifdef _OPENMP
   nthreads = omp_get_max_threads()
#else
   nthreads = 1
#endif
   IF (nthreads <= 0) nthreads = 1

   ! element sizes in bytes (portable)
   SZ_INT_BYTES  = STORAGE_SIZE(0) / 8
   SZ_CPLX_BYTES = STORAGE_SIZE(CMPLX(0.0_dp, 0.0_dp, dp)) / 8

   ALLOCATE(buf(nthreads), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc buf(:) failed')

   ! heuristic initial capacity per thread (can be tuned)
   init_cap = MAX(200000, 10*N / MAX(1, nthreads))

   DO t = 1, nthreads
      buf(t)%cap      = init_cap
      buf(t)%nnz_used = 0
      ALLOCATE(buf(t)%row(buf(t)%cap), buf(t)%col(buf(t)%cap), buf(t)%val(buf(t)%cap), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_CSR_par: per-thread buffer alloc failed')
      ALLOCATE(buf(t)%row_count(N), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_CSR_par: per-thread row_count alloc failed')
      buf(t)%row_count = 0
   END DO

   ! Optional: report per-thread COO scratch capacity
   IF (DEBUG_OUTPUT) THEN
      bytes_coo_tot = 0.0_dp
      DO t = 1, nthreads
         bytes_coo_tot = bytes_coo_tot + REAL(buf(t)%cap, dp) * &
                        ( REAL(2*SZ_INT_BYTES + SZ_CPLX_BYTES, dp) )
      END DO
      GB_coo = bytes_coo_tot / (1024._dp**3)
      if (my_rank==0) WRITE(*,'(A,1X,F8.3,1X,A,1X,I4,1X,A,1X,I10,1X,A,F8.3)') &
           'CA_CSR_par: thread COO buffers  FREQ=', FREQ, &
           'nthreads=', nthreads, 'init_cap=', init_cap, 'GB≈', GB_coo
   END IF

   ! ------------------------------------------------------------------
   ! 2. Run CA_core_par; ADD_ENTRY fills per-thread COO buffers + histograms
   ! ------------------------------------------------------------------

   CALL CA_core_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                    AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)

   ! ------------------------------------------------------------------
   ! 3. Build global row_nnz and IA (for raw entries)
   ! ------------------------------------------------------------------
   ALLOCATE(row_nnz(N), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc row_nnz failed')
   row_nnz = 0

   DO t = 1, nthreads
      DO i = 1, N
         row_nnz(i) = row_nnz(i) + buf(t)%row_count(i)
      END DO
   END DO

   ALLOCATE(IA(N+1), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc IA failed')

   IA(1) = 1
   DO i = 1, N
      IA(i+1) = IA(i) + row_nnz(i)
   END DO

   nnz_total = IA(N+1) - 1
   IF (nnz_total < 0) nnz_total = 0

   ALLOCATE(JA(MAX(1, nnz_total)), AV(MAX(1, nnz_total)), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc JA/AV failed')

   ! --- DEBUG: raw CSR size (before dedup) ---
   bytes_A_raw  = REAL(nnz_total, dp) * REAL(SZ_CPLX_BYTES, dp)
   bytes_IA_raw = REAL(N+1,      dp) * REAL(SZ_INT_BYTES,  dp)
   bytes_JA_raw = REAL(nnz_total, dp) * REAL(SZ_INT_BYTES, dp)
   GB_raw = (bytes_A_raw + bytes_IA_raw + bytes_JA_raw) / (1024._dp**3)

   IF (DEBUG_OUTPUT.and.my_rank==0) THEN
      WRITE(*,'(A,1X,F8.3,1X,A,1X,I10,1X,A,1X,I12,1X,A,F8.3)') &
           'CA_CSR_par: RAW CSR  FREQ=', FREQ, &
           'N=', N, 'nnz_total=', nnz_total, 'GB≈', GB_raw
   END IF
   ! ------------------------------------------

   ! ------------------------------------------------------------------
   ! 4-5. Compute per-thread row offsets and scatter per-thread COO
   !      into global CSR (unsorted within row)
   ! ------------------------------------------------------------------

   ALLOCATE(row_off(N), t_row_start(N), row_next(N), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc row_off/t_row_start/row_next failed')

   ! row_off(i) accumulates how many slots in row i have been reserved by previous threads
   row_off = 0

   DO t = 1, nthreads

      ! For current thread t, compute its start offset in each row and init cursor
      DO i = 1, N
         t_row_start(i) = IA(i) + row_off(i)
         row_next(i)    = t_row_start(i)
         row_off(i)     = row_off(i) + buf(t)%row_count(i)
      END DO

      ! Scatter thread t COO entries into the slice assigned to thread t within each row
      DO k = 1, buf(t)%nnz_used
         r = buf(t)%row(k)
         c = buf(t)%col(k)

         ! safety: row range
         IF (r < 1 .OR. r > N) THEN
            CALL fail('CA_CSR_par: row index out of range in scatter')
         END IF

         offset = row_next(r)

         ! offset must remain inside this thread's reserved slice for row r:
         ! [t_row_start(r), t_row_start(r) + row_count(r) - 1]
         IF (offset < t_row_start(r) .OR. offset >= t_row_start(r) + buf(t)%row_count(r)) THEN
            CALL fail('CA_CSR_par: offset out of thread row slice')
         END IF

         JA(offset) = c
         AV(offset) = buf(t)%val(k)
         row_next(r) = offset + 1
      END DO

   END DO  ! t loop

   DEALLOCATE(row_off, t_row_start, row_next)

   ! ------------------------------------------------------------------
   ! 6. Per-row sort, deduplicate (sum duplicates) WITHOUT zero-drop
   ! ------------------------------------------------------------------
   ALLOCATE(BUCKET(N), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc BUCKET failed')
   BUCKET = 0

   t_core_start = omp_get_wtime()
!$omp parallel do default(shared) private(i,lb,ub,w,kk,cur_col,cur_val) schedule(static)
   DO i = 1, N
      lb = IA(i)
      ub = IA(i+1) - 1
      IF (ub < lb) CYCLE

      ! sort row by column index
      CALL SORT_ROW_SEGMENT(JA, AV, lb, ub)

      ! dedup only (keep zeros) in-place in [lb,ub]
      w = lb
      cur_col = JA(lb)
      cur_val = AV(lb)

      DO kk = lb + 1, ub
         IF (JA(kk) == cur_col) THEN
            cur_val = cur_val + AV(kk)
         ELSE
            JA(w) = cur_col
            AV(w) = cur_val
            w = w + 1

            cur_col = JA(kk)
            cur_val = AV(kk)
         END IF
      END DO

      ! flush last accumulated
      JA(w) = cur_col
      AV(w) = cur_val
      w = w + 1

      BUCKET(i) = w - lb
   END DO
!$omp end parallel do
   t_core_end = omp_get_wtime()

   IF (DEBUG_OUTPUT.and.my_rank==0) THEN
      WRITE(*,'(A,1X,F8.3,1X,A,F10.3)') 'CA_CSR_par: sort/dedup time (s)=', &
           t_core_end - t_core_start, ' FREQ=', FREQ
   END IF

   ! ------------------------------------------------------------------
   ! 7. Build final IA_new and compact JA_new/AV_new to exact nnz_final
   ! ------------------------------------------------------------------
   ALLOCATE(IA_new(N+1), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_CSR_par: alloc IA_new failed')
   IA_new(1) = 1
   DO i = 1, N
      IA_new(i+1) = IA_new(i) + BUCKET(i)
   END DO

   nnz_final = IA_new(N+1) - 1
   IF (nnz_final < 0) nnz_final = 0

   ! compute final CSR memory (always printed)
   bytes_A_final  = REAL(MAX(nnz_final, 0), dp) * REAL(SZ_CPLX_BYTES, dp)
   bytes_IA_final = REAL(N+1,               dp) * REAL(SZ_INT_BYTES,  dp)
   bytes_JA_final = REAL(MAX(nnz_final, 0), dp) * REAL(SZ_INT_BYTES,  dp)
   GB_final = (bytes_A_final + bytes_IA_final + bytes_JA_final) / (1024._dp**3)

   IF (nnz_final <= 0) THEN
      ! empty matrix case
      DEALLOCATE(JA, AV)
      ALLOCATE(JA(1), AV(1), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_CSR_par: alloc empty JA/AV failed')
      JA(1) = 0
      AV(1) = CMPLX(0.0_dp, 0.0_dp, dp)

      DEALLOCATE(IA)
      CALL MOVE_ALLOC(IA_new, IA)
   ELSE
      ALLOCATE(JA_new(nnz_final), AV_new(nnz_final), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_CSR_par: alloc JA_new/AV_new failed')

      !$omp parallel do default(shared) private(i,lb,ub) schedule(static)
      DO i = 1, N
         lb = IA(i)
         ub = lb + BUCKET(i) - 1
         IF (BUCKET(i) > 0) THEN
            JA_new(IA_new(i):IA_new(i+1)-1) = JA(lb:ub)
            AV_new(IA_new(i):IA_new(i+1)-1) = AV(lb:ub)
         END IF
      END DO
      !$omp end parallel do

      DEALLOCATE(JA, AV)
      DEALLOCATE(IA)
      CALL MOVE_ALLOC(IA_new, IA)
      CALL MOVE_ALLOC(JA_new, JA)
      CALL MOVE_ALLOC(AV_new, AV)
   END IF
   ! --- FINAL size print (always) ---
!    IF (DEBUG_OUTPUT.and.my_rank==0) THEN 
!       WRITE(*,'(A,1X,F8.3,1X,A,1X,I10,1X,A,1X,I12,1X,A,F8.3)') &
!         'CA_CSR_par: FINAL CSR FREQ=', FREQ, &
!         'N=', N, 'nnz_final=', nnz_final, 'A_CSR (GB)≈', GB_final
! END IF
   ! ---------------------------------

   ! ------------------------------------------------------------------
   ! 8. Cleanup per-thread buffers and helpers
   ! ------------------------------------------------------------------
   DO t = 1, nthreads
      IF (ALLOCATED(buf(t)%row))       DEALLOCATE(buf(t)%row)
      IF (ALLOCATED(buf(t)%col))       DEALLOCATE(buf(t)%col)
      IF (ALLOCATED(buf(t)%val))       DEALLOCATE(buf(t)%val)
      IF (ALLOCATED(buf(t)%row_count)) DEALLOCATE(buf(t)%row_count)
   END DO
   DEALLOCATE(buf)
   DEALLOCATE(row_nnz, BUCKET)

   RETURN

CONTAINS

   SUBROUTINE GROW_LOCAL(tid)
      INTEGER, INTENT(IN) :: tid
      INTEGER :: newcap, ierr
      INTEGER,     ALLOCATABLE :: r2(:), c2(:)
      COMPLEX(dp), ALLOCATABLE :: v2(:)

      newcap = INT(1.5_dp*buf(tid)%cap) + 100000

      ALLOCATE(r2(newcap), c2(newcap), v2(newcap), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_CSR_par: GROW_LOCAL alloc failed')

      IF (buf(tid)%nnz_used > 0) THEN
         r2(1:buf(tid)%nnz_used) = buf(tid)%row(1:buf(tid)%nnz_used)
         c2(1:buf(tid)%nnz_used) = buf(tid)%col(1:buf(tid)%nnz_used)
         v2(1:buf(tid)%nnz_used) = buf(tid)%val(1:buf(tid)%nnz_used)
      END IF

#ifdef _OPENMP
!$omp critical (grow_local_critical)
#endif
      CALL MOVE_ALLOC(r2, buf(tid)%row)
      CALL MOVE_ALLOC(c2, buf(tid)%col)
      CALL MOVE_ALLOC(v2, buf(tid)%val)
      buf(tid)%cap = newcap
#ifdef _OPENMP
!$omp end critical (grow_local_critical)
#endif
   END SUBROUTINE GROW_LOCAL

   SUBROUTINE ADD_ENTRY(row, col, val)
      USE shared_mod, ONLY: dp
      IMPLICIT NONE
      INTEGER,     INTENT(IN) :: row, col
      COMPLEX(dp), INTENT(IN) :: val
      INTEGER :: tid, pos

      IF (ABS(val) < ZERO_TOL) RETURN
      IF (row < 1 .OR. row > N .OR. col < 1 .OR. col > N) THEN
         CALL fail('CA_CSR_par: ADD_ENTRY indices out of range')
      END IF

      ! band mask
      IF ((row - col) > KL .OR. (col - row) > KU) RETURN

#ifdef _OPENMP
      tid = omp_get_thread_num() + 1
#else
      tid = 1
#endif

      pos = buf(tid)%nnz_used + 1
      IF (pos > buf(tid)%cap) CALL GROW_LOCAL(tid)

      buf(tid)%row(pos) = row
      buf(tid)%col(pos) = col
      buf(tid)%val(pos) = val
      buf(tid)%nnz_used = pos

      ! per-thread row histogram
      buf(tid)%row_count(row) = buf(tid)%row_count(row) + 1
   END SUBROUTINE ADD_ENTRY

   ! Sort row segment JA(lb:ub) and permute AV consistently (insertion sort)
   SUBROUTINE SORT_ROW_SEGMENT(JA, AV, lb, ub)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: JA(:)
      COMPLEX(dp), INTENT(INOUT) :: AV(:)
      INTEGER, INTENT(IN) :: lb, ub
      INTEGER :: idx, j, keyj
      COMPLEX(dp) :: keyz

      IF (ub - lb <= 0) RETURN

      DO idx = lb + 1, ub
         keyj = JA(idx); keyz = AV(idx); j = idx - 1
         DO WHILE (j >= lb)
            IF (JA(j) <= keyj) EXIT
            JA(j+1) = JA(j)
            AV(j+1) = AV(j)
            j = j - 1
         END DO
         JA(j+1) = keyj
         AV(j+1) = keyz
      END DO
   END SUBROUTINE SORT_ROW_SEGMENT

END SUBROUTINE CA_CSR_par

SUBROUTINE CA_COO_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                      CI, NORD, AS, WT, IE0, IS0, DZ0, N,IRN, JCN, ACOO,NBLOCK, LDAB, KL, KU,          &
                       NNZ_OUT, IG, my_rank)
   USE shared_mod,  ONLY: dp
   USE err_mpi_mod, ONLY: fail
   USE omp_lib
   IMPLICIT NONE

   ! ---- inputs ----
   INTEGER, INTENT(IN) :: N, KL, KU, LDAB,NBLOCK,my_rank
   INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ
   TYPE(InversionGridType), INTENT(IN) :: IG
   REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
   REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0

   ! ---- outputs (COO) ----
   INTEGER,          ALLOCATABLE, INTENT(OUT) :: IRN(:), JCN(:)
   COMPLEX(dp),      ALLOCATABLE, INTENT(OUT) :: ACOO(:)
   INTEGER,          INTENT(OUT)              :: NNZ_OUT

   ! ---- per-thread buffer type ----
   TYPE ThreadBuf
      INTEGER,     ALLOCATABLE :: irn(:), jcn(:)
      COMPLEX(dp), ALLOCATABLE :: a(:)
      INTEGER                   :: nnz_used
      INTEGER                   :: cap
   END TYPE ThreadBuf

   TYPE(ThreadBuf), ALLOCATABLE :: buf(:)

   INTEGER :: ierr, nthreads, t, init_cap
   INTEGER :: nnz_total, offset
   REAL(dp), PARAMETER :: ZERO_TOL = 1.0e-30_dp

   ! ---- dedup helpers (CSR) ----
   INTEGER,          ALLOCATABLE :: IA_csr(:), JA_csr(:), rowcnt(:), nxt(:)
   COMPLEX(dp),      ALLOCATABLE :: AV_csr(:)

   INTEGER,          ALLOCATABLE :: BUCKET(:), IA_new(:), JA_new(:)
   COMPLEX(dp),      ALLOCATABLE :: AV_new(:)

   INTEGER :: i, k, lb, ub, p, nnz_final
   INTEGER :: w, kk, cur_col
   COMPLEX(dp) :: cur_val
   REAL(dp) :: mean_i, mean_j, mean_re, mean_im, mean_abs
   LOGICAL :: dbg_coo
   INTEGER :: SZ_INT_BYTES, SZ_CPLX_BYTES
   REAL(dp) :: bytes_A_final, bytes_IRN_final, bytes_JCN_final, GB_final

   ! ------------------------------------------------------------------
   ! 1. Allocate per-thread COO buffers
   ! ------------------------------------------------------------------
#ifdef _OPENMP
   nthreads = omp_get_max_threads()
#else
   nthreads = 1
#endif
   IF (nthreads <= 0) nthreads = 1
   dbg_coo = (my_rank == 0)
   SZ_INT_BYTES  = STORAGE_SIZE(0) / 8
   SZ_CPLX_BYTES = STORAGE_SIZE(CMPLX(0.0_dp, 0.0_dp, dp)) / 8

   ALLOCATE(buf(nthreads), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_COO_par: alloc buf(:) failed')

   init_cap = MAX(200000, 10*N / MAX(1, nthreads))

   DO t = 1, nthreads
      buf(t)%cap      = init_cap
      buf(t)%nnz_used = 0
      ALLOCATE(buf(t)%irn(buf(t)%cap), buf(t)%jcn(buf(t)%cap), buf(t)%a(buf(t)%cap), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO_par: per-thread buffer alloc failed')
   END DO

   ! ------------------------------------------------------------------
   ! 2. Run CA_core_par in parallel; ADD_ENTRY fills per-thread buffers
   ! ------------------------------------------------------------------
   CALL CA_core_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                    CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, IG, ADD_ENTRY)

   ! ------------------------------------------------------------------
   ! 3. Gather per-thread buffers into a raw COO (may have duplicates)
   ! ------------------------------------------------------------------
   nnz_total = 0
   DO t = 1, nthreads
      nnz_total = nnz_total + buf(t)%nnz_used
   END DO

   NNZ_OUT = nnz_total

   ALLOCATE(IRN(MAX(1, nnz_total)), JCN(MAX(1, nnz_total)), ACOO(MAX(1, nnz_total)), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_COO_par: final COO alloc failed')

   IF (nnz_total > 0) THEN
      offset = 0
      DO t = 1, nthreads
         IF (buf(t)%nnz_used > 0) THEN
            IRN(offset+1:offset+buf(t)%nnz_used)  = buf(t)%irn(1:buf(t)%nnz_used)
            JCN(offset+1:offset+buf(t)%nnz_used)  = buf(t)%jcn(1:buf(t)%nnz_used)
            ACOO(offset+1:offset+buf(t)%nnz_used) = buf(t)%a(1:buf(t)%nnz_used)
            offset = offset + buf(t)%nnz_used
         END IF
      END DO
   END IF

   IF (dbg_coo) THEN
      IF (nnz_total > 0) THEN
         mean_i = SUM(REAL(IRN(1:nnz_total), dp))/REAL(nnz_total, dp)
         mean_j = SUM(REAL(JCN(1:nnz_total), dp))/REAL(nnz_total, dp)
         mean_re = SUM(REAL(ACOO(1:nnz_total), dp))/REAL(nnz_total, dp)
         mean_im = SUM(AIMAG(ACOO(1:nnz_total)))/REAL(nnz_total, dp)
         mean_abs = SUM(ABS(ACOO(1:nnz_total)))/REAL(nnz_total, dp)
      ELSE
         mean_i = 0.0_dp
         mean_j = 0.0_dp
         mean_re = 0.0_dp
         mean_im = 0.0_dp
         mean_abs = 0.0_dp
      END IF
!      WRITE (*, '(A,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
!         'CA_COO_par RAW:', 'nnz=', nnz_total, 'size(IA)=', SIZE(IRN), 'size(JA)=', SIZE(JCN), 'size(A)=', SIZE(ACOO)
!      WRITE (*, '(A,1X,A,ES13.5,1X,A,ES13.5,1X,A,ES13.5,1X,A,ES13.5,1X,A,ES13.5)') &
!         'CA_COO_par RAW mean:', 'i=', mean_i, 'j=', mean_j, 'Re(A)=', mean_re, 'Im(A)=', mean_im, 'Abs(A)=', mean_abs
!      CALL flush(6)
   END IF

   ! ------------------------------------------------------------------
   ! 4. Cleanup thread buffers (done early; dedup works on gathered COO)
   ! ------------------------------------------------------------------
   DO t = 1, nthreads
      IF (ALLOCATED(buf(t)%irn)) DEALLOCATE (buf(t)%irn)
      IF (ALLOCATED(buf(t)%jcn)) DEALLOCATE (buf(t)%jcn)
      IF (ALLOCATED(buf(t)%a))   DEALLOCATE (buf(t)%a)
   END DO
   DEALLOCATE(buf)

   ! ------------------------------------------------------------------
   ! 5. COO -> CSR (row-bucket) to enable per-row sort+dedup+sum
   !     NOTE: NO ZERO-DROP (matches your CSR_par behavior)
   ! ------------------------------------------------------------------
   IF (nnz_total <= 0) THEN
      ! keep as a valid empty COO
      NNZ_OUT = 0
      RETURN
   END IF

   ALLOCATE(rowcnt(N), IA_csr(N+1), nxt(N), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_COO_par: alloc rowcnt/IA_csr/nxt failed')
   rowcnt = 0

   DO k = 1, nnz_total
      i = IRN(k)
      IF (i < 1 .OR. i > N) CALL fail('CA_COO_par: IRN out of range')
      rowcnt(i) = rowcnt(i) + 1
   END DO

   IA_csr(1) = 1
   DO i = 1, N
      IA_csr(i+1) = IA_csr(i) + rowcnt(i)
   END DO

   ALLOCATE(JA_csr(nnz_total), AV_csr(nnz_total), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_COO_par: alloc JA_csr/AV_csr failed')

   nxt = IA_csr(1:N)
   DO k = 1, nnz_total
      i = IRN(k)
      p = nxt(i)
      JA_csr(p) = JCN(k)
      AV_csr(p) = ACOO(k)
      nxt(i) = p + 1
   END DO

   ! we will rebuild output COO from deduped CSR, so drop raw COO now
   DEALLOCATE(IRN, JCN, ACOO)

   ! ------------------------------------------------------------------
   ! 6. Per-row sort, deduplicate (sum duplicates) WITHOUT zero-drop
   ! ------------------------------------------------------------------
   ALLOCATE(BUCKET(N), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_COO_par: alloc BUCKET failed')
   BUCKET = 0

!$omp parallel do default(shared) private(i,lb,ub,w,kk,cur_col,cur_val) schedule(static)
   DO i = 1, N
      lb = IA_csr(i)
      ub = IA_csr(i+1) - 1
      IF (ub < lb) CYCLE

      ! sort row by column index
      CALL SORT_ROW_SEGMENT(JA_csr, AV_csr, lb, ub)

      ! dedup only (keep zeros) in-place in [lb,ub]
      w = lb
      cur_col = JA_csr(lb)
      cur_val = AV_csr(lb)

      DO kk = lb + 1, ub
         IF (JA_csr(kk) == cur_col) THEN
            cur_val = cur_val + AV_csr(kk)
         ELSE
            JA_csr(w) = cur_col
            AV_csr(w) = cur_val
            w = w + 1
            cur_col = JA_csr(kk)
            cur_val = AV_csr(kk)
         END IF
      END DO

      ! flush last accumulated
      JA_csr(w) = cur_col
      AV_csr(w) = cur_val
      w = w + 1

      BUCKET(i) = w - lb
   END DO
!$omp end parallel do

   ! ------------------------------------------------------------------
   ! 7. Build final IA_new and compact JA_new/AV_new to exact nnz_final
   ! ------------------------------------------------------------------
   ALLOCATE(IA_new(N+1), STAT=ierr)
   IF (ierr /= 0) CALL fail('CA_COO_par: alloc IA_new failed')
   IA_new(1) = 1
   DO i = 1, N
      IA_new(i+1) = IA_new(i) + BUCKET(i)
   END DO

   nnz_final = IA_new(N+1) - 1
   IF (nnz_final < 0) nnz_final = 0
   bytes_A_final   = REAL(MAX(nnz_final, 0), dp) * REAL(SZ_CPLX_BYTES, dp)
   bytes_IRN_final = REAL(MAX(nnz_final, 0), dp) * REAL(SZ_INT_BYTES,  dp)
   bytes_JCN_final = REAL(MAX(nnz_final, 0), dp) * REAL(SZ_INT_BYTES,  dp)
   GB_final = (bytes_A_final + bytes_IRN_final + bytes_JCN_final) / (1024._dp**3)

   IF (nnz_final <= 0) THEN
      ! empty matrix case (still return valid COO arrays of size 1)
      NNZ_OUT = 0
      ALLOCATE(IRN(1), JCN(1), ACOO(1), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO_par: alloc empty COO failed')
      IRN(1)  = 0
      JCN(1)  = 0
      ACOO(1) = CMPLX(0.0_dp, 0.0_dp, dp)
   ELSE
      ALLOCATE(JA_new(nnz_final), AV_new(nnz_final), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO_par: alloc JA_new/AV_new failed')

      DO i = 1, N
         lb = IA_csr(i)
         ub = lb + BUCKET(i) - 1
         IF (BUCKET(i) > 0) THEN
            JA_new(IA_new(i):IA_new(i+1)-1) = JA_csr(lb:ub)
            AV_new(IA_new(i):IA_new(i+1)-1) = AV_csr(lb:ub)
         END IF
      END DO

      ! ----------------------------------------------------------------
      ! 8. CSR -> COO (deduped) for MUMPS: IRN/JCN/ACOO
      ! ----------------------------------------------------------------
      ALLOCATE(IRN(nnz_final), JCN(nnz_final), ACOO(nnz_final), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO_par: alloc final COO failed')

      p = 0
      DO i = 1, N
         lb = IA_new(i)
         ub = IA_new(i+1) - 1
         DO k = lb, ub
            p = p + 1
            IRN(p)  = i
            JCN(p)  = JA_new(k)
            ACOO(p) = AV_new(k)
         END DO
      END DO

      NNZ_OUT = nnz_final
   END IF
   ! if(my_rank==0) WRITE(*,'(A,1X,F8.3,1X,A,1X,I10,1X,A,1X,I12,1X,A,F8.3)') &
   !      'CA_COO_par: FINAL COO FREQ=', FREQ, &
   !      'N=', N, 'nnz_final=', nnz_final, 'A_COO (GB)≈', GB_final

   IF (dbg_coo) THEN
      IF (NNZ_OUT > 0) THEN
         mean_i = SUM(REAL(IRN(1:NNZ_OUT), dp))/REAL(NNZ_OUT, dp)
         mean_j = SUM(REAL(JCN(1:NNZ_OUT), dp))/REAL(NNZ_OUT, dp)
         mean_re = SUM(REAL(ACOO(1:NNZ_OUT), dp))/REAL(NNZ_OUT, dp)
         mean_im = SUM(AIMAG(ACOO(1:NNZ_OUT)))/REAL(NNZ_OUT, dp)
         mean_abs = SUM(ABS(ACOO(1:NNZ_OUT)))/REAL(NNZ_OUT, dp)
      ELSE
         mean_i = 0.0_dp
         mean_j = 0.0_dp
         mean_re = 0.0_dp
         mean_im = 0.0_dp
         mean_abs = 0.0_dp
      END IF
!      WRITE (*, '(A,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
!         'CA_COO_par FINAL:', 'nnz=', NNZ_OUT, 'size(IA)=', SIZE(IRN), 'size(JA)=', SIZE(JCN), 'size(A)=', SIZE(ACOO)
!      WRITE (*, '(A,1X,A,ES13.5,1X,A,ES13.5,1X,A,ES13.5,1X,A,ES13.5,1X,A,ES13.5)') &
!         'CA_COO_par FINAL mean:', 'i=', mean_i, 'j=', mean_j, 'Re(A)=', mean_re, 'Im(A)=', mean_im, 'Abs(A)=', mean_abs
!      CALL flush(6)
   END IF

   ! ------------------------------------------------------------------
   ! 9. Cleanup CSR temporaries
   ! ------------------------------------------------------------------
   IF (ALLOCATED(rowcnt))  DEALLOCATE(rowcnt)
   IF (ALLOCATED(nxt))     DEALLOCATE(nxt)
   IF (ALLOCATED(BUCKET))  DEALLOCATE(BUCKET)

   IF (ALLOCATED(IA_csr))  DEALLOCATE(IA_csr)
   IF (ALLOCATED(JA_csr))  DEALLOCATE(JA_csr)
   IF (ALLOCATED(AV_csr))  DEALLOCATE(AV_csr)

   IF (ALLOCATED(IA_new))  DEALLOCATE(IA_new)
   IF (ALLOCATED(JA_new))  DEALLOCATE(JA_new)
   IF (ALLOCATED(AV_new))  DEALLOCATE(AV_new)

   RETURN

CONTAINS

   SUBROUTINE GROW_LOCAL(tid)
      INTEGER, INTENT(IN) :: tid
      INTEGER :: newcap, ierr
      INTEGER,     ALLOCATABLE :: irn2(:), jcn2(:)
      COMPLEX(dp), ALLOCATABLE :: a2(:)

      newcap = INT(1.5_dp*buf(tid)%cap) + 100000

      ALLOCATE(irn2(newcap), jcn2(newcap), a2(newcap), STAT=ierr)
      IF (ierr /= 0) CALL fail('CA_COO_par: GROW_LOCAL alloc failed')

      IF (buf(tid)%nnz_used > 0) THEN
         irn2(1:buf(tid)%nnz_used) = buf(tid)%irn(1:buf(tid)%nnz_used)
         jcn2(1:buf(tid)%nnz_used) = buf(tid)%jcn(1:buf(tid)%nnz_used)
         a2(1:buf(tid)%nnz_used)   = buf(tid)%a(1:buf(tid)%nnz_used)
      END IF

      CALL MOVE_ALLOC(irn2, buf(tid)%irn)
      CALL MOVE_ALLOC(jcn2, buf(tid)%jcn)
      CALL MOVE_ALLOC(a2,   buf(tid)%a)

      buf(tid)%cap = newcap
   END SUBROUTINE GROW_LOCAL

   SUBROUTINE ADD_ENTRY(row, col, val)
      USE shared_mod, ONLY: dp
      IMPLICIT NONE
      INTEGER,     INTENT(IN) :: row, col
      COMPLEX(dp), INTENT(IN) :: val
      INTEGER :: tid, pos

      IF (ABS(val) < ZERO_TOL) RETURN
      IF (row < 1 .OR. row > N .OR. col < 1 .OR. col > N) THEN
         CALL fail('CA_COO_par: ADD_ENTRY indices out of range')
      END IF

#ifdef _OPENMP
      tid = omp_get_thread_num() + 1      ! 1-based
#else
      tid = 1
#endif

      pos = buf(tid)%nnz_used + 1
      IF (pos > buf(tid)%cap) CALL GROW_LOCAL(tid)

      buf(tid)%irn(pos) = row
      buf(tid)%jcn(pos) = col
      buf(tid)%a(pos)   = val
      buf(tid)%nnz_used = pos
   END SUBROUTINE ADD_ENTRY

   ! ---------------------------------------------------------------
   ! Sort JA(lb:ub) ascending, applying same permutation to AV.
   ! You can replace this with your existing sorter if you have one.
   ! ---------------------------------------------------------------
   SUBROUTINE SORT_ROW_SEGMENT(JA, AV, lb, ub)
      USE shared_mod, ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: JA(:)
      COMPLEX(dp), INTENT(INOUT) :: AV(:)
      INTEGER, INTENT(IN) :: lb, ub
      CALL QSORT_INT_CPLX(JA, AV, lb, ub)
   END SUBROUTINE SORT_ROW_SEGMENT

   RECURSIVE SUBROUTINE QSORT_INT_CPLX(JA, AV, l, r)
      USE shared_mod, ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: JA(:)
      COMPLEX(dp), INTENT(INOUT) :: AV(:)
      INTEGER, INTENT(IN) :: l, r
      INTEGER :: i, j, piv, tmpi
      COMPLEX(dp) :: tmpv

      IF (l >= r) RETURN
      piv = JA((l+r)/2)
      i = l
      j = r
      DO
         DO WHILE (JA(i) < piv); i = i + 1; END DO
         DO WHILE (JA(j) > piv); j = j - 1; END DO
         IF (i <= j) THEN
            tmpi = JA(i); JA(i) = JA(j); JA(j) = tmpi
            tmpv = AV(i); AV(i) = AV(j); AV(j) = tmpv
            i = i + 1
            j = j - 1
         END IF
         IF (i > j) EXIT
      END DO
      IF (l < j) CALL QSORT_INT_CPLX(JA, AV, l, j)
      IF (i < r) CALL QSORT_INT_CPLX(JA, AV, i, r)
   END SUBROUTINE QSORT_INT_CPLX

END SUBROUTINE CA_COO_par

!-----------------------------------------------------------------------
! Dump CSR as triplet text file: first line "N nnz", then lines: row col re im
! Uses 1-based indices as produced by CA_CSR / CA_CSR_par.
!-----------------------------------------------------------------------
SUBROUTINE DUMP_CSR_TRIPLET(filename, N, IA, JA, AV)
  USE iso_fortran_env, ONLY: dp => real64
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: IA(:), JA(:)
  COMPLEX(dp), INTENT(IN) :: AV(:)
  INTEGER :: i, p, row_start, row_end
  INTEGER :: unit, ierr
  REAL(dp) :: re, im
  INTEGER :: nnz

  nnz = IA(N+1) - 1
  OPEN(NEWUNIT=unit, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierr)
  IF (ierr /= 0) THEN
     WRITE(*,*) 'DUMP_CSR_TRIPLET: failed to open ', TRIM(filename), ' IOSTAT=', ierr
     RETURN
  END IF

  ! header
  WRITE(unit, '(I12,1X,I12)') N, nnz

  DO i = 1, N
     row_start = IA(i)
     row_end   = IA(i+1) - 1
     IF (row_end >= row_start) THEN
        DO p = row_start, row_end
           re = REAL(AV(p))
           im = AIMAG(AV(p))
           WRITE(unit, '(I8,1X,I8,1X,1PE25.16,1X,1PE25.16)') i, JA(p), re, im
        END DO
     END IF
  END DO

  CLOSE(unit)
END SUBROUTINE DUMP_CSR_TRIPLET

SUBROUTINE COMPARE_CSR(N, IA1, JA1, A1, IA2, JA2, A2)
  USE iso_fortran_env, ONLY: dp => real64
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: IA1(:), JA1(:), IA2(:), JA2(:)
  COMPLEX(dp), INTENT(IN) :: A1(:), A2(:)

  INTEGER :: nnz1, nnz2
  INTEGER :: i, p1, p2, s1, e1, s2, e2
  INTEGER :: c1, c2, bigint
  COMPLEX(dp) :: v1, v2
  REAL(dp) :: diff, max_diff
  REAL(dp), PARAMETER :: tol_abs = 1.0e-12_dp, tol_rel = 1.0e-10_dp

  INTEGER :: count_only1, count_only2, count_diff
  INTEGER, PARAMETER :: max_report = 20
  INTEGER :: report_only1, report_only2, report_diff
  INTEGER :: row_missing_diag, row_tiny_diag
  REAL(dp) :: mag1, mag2, mag_diag

  LOGICAL :: ok_ia

  nnz1 = IA1(N+1) - 1
  nnz2 = IA2(N+1) - 1

  WRITE(*,'(A,I12)') 'COMPARE_CSR: nnz (serial)   = ', nnz1
  WRITE(*,'(A,I12)') 'COMPARE_CSR: nnz (parallel) = ', nnz2

  ! Basic IA monotonicity checks
  ok_ia = .TRUE.
  DO i = 1, N
     IF (IA1(i+1) < IA1(i)) THEN
        WRITE(*,'(A,I6)') 'COMPARE_CSR: IA1 not monotone at row ', i
        ok_ia = .FALSE.
     END IF
     IF (IA2(i+1) < IA2(i)) THEN
        WRITE(*,'(A,I6)') 'COMPARE_CSR: IA2 not monotone at row ', i
        ok_ia = .FALSE.
     END IF
  END DO
  IF (.NOT. ok_ia) THEN
     WRITE(*,*) 'COMPARE_CSR: IA monotonicity failed, stopping comparison.'
     RETURN
  END IF

  bigint = HUGE(1)

  count_only1  = 0
  count_only2  = 0
  count_diff   = 0
  report_only1 = 0
  report_only2 = 0
  report_diff  = 0
  row_missing_diag = 0
  row_tiny_diag    = 0
  max_diff = 0.0_dp

  WRITE(*,*) 'COMPARE_CSR: scanning rows...'

  DO i = 1, N
     s1 = IA1(i)
     e1 = IA1(i+1) - 1
     s2 = IA2(i)
     e2 = IA2(i+1) - 1

     p1 = s1
     p2 = s2

     ! track diagonal existence/magnitude in parallel matrix
     mag_diag = -1.0_dp

     DO
        IF (p1 > e1 .AND. p2 > e2) EXIT

        IF (p1 <= e1) THEN
           c1 = JA1(p1)
        ELSE
           c1 = bigint
        END IF

        IF (p2 <= e2) THEN
           c2 = JA2(p2)
        ELSE
           c2 = bigint
        END IF

        IF (c1 == c2) THEN
           v1 = A1(p1)
           v2 = A2(p2)

           IF (c1 == i) THEN
              mag_diag = ABS(v2)
           END IF

           ! check for NaNs / Infs
           IF ( (REAL(v1) /= REAL(v1)) .OR. (AIMAG(v1) /= AIMAG(v1)) ) THEN
              WRITE(*,'(A,I6,A,I6,A)') 'WARNING: NaN in serial at row ', i, ', col ', c1, ''
           END IF
           IF ( (REAL(v2) /= REAL(v2)) .OR. (AIMAG(v2) /= AIMAG(v2)) ) THEN
              WRITE(*,'(A,I6,A,I6,A)') 'WARNING: NaN in parallel at row ', i, ', col ', c2, ''
           END IF

           mag1 = ABS(v1)
           mag2 = ABS(v2)
           diff = ABS(v1 - v2)

           IF (diff > max_diff) max_diff = diff

           IF (diff > tol_abs .AND. diff > tol_rel*MAX(mag1, mag2)) THEN
              count_diff = count_diff + 1
              IF (report_diff < max_report) THEN
                 report_diff = report_diff + 1
                 WRITE(*,'(A,I6,A,I6,A,1PE12.4,A,2(1PE12.4,1X,1PE12.4))') &
                      'DIFF: row=', i, ' col=', c1, ' diff=', diff, '  v1=', REAL(v1), AIMAG(v1), &
                      '  v2=', REAL(v2), AIMAG(v2)
              END IF
           END IF

           p1 = p1 + 1
           p2 = p2 + 1

        ELSE IF (c1 < c2) THEN
           ! entry only in serial
           count_only1 = count_only1 + 1
           IF (report_only1 < max_report) THEN
              report_only1 = report_only1 + 1
              WRITE(*,'(A,I6,A,I6)') 'ONLY_SERIAL: row=', i, ' col=', c1
           END IF
           p1 = p1 + 1

        ELSE
           ! entry only in parallel
           count_only2 = count_only2 + 1
           IF (report_only2 < max_report) THEN
              report_only2 = report_only2 + 1
              WRITE(*,'(A,I6,A,I6)') 'ONLY_PAR: row=', i, ' col=', c2
           END IF
           IF (c2 == i) THEN
              mag_diag = ABS(A2(p2))
           END IF
           p2 = p2 + 1
        END IF

     END DO

     IF (mag_diag < 0.0_dp) THEN
        row_missing_diag = row_missing_diag + 1
        IF (row_missing_diag <= max_report) THEN
           WRITE(*,'(A,I6)') 'MISSING_DIAG_PAR: row=', i
        END IF
     ELSE IF (mag_diag < 1.0e-18_dp) THEN
        row_tiny_diag = row_tiny_diag + 1
        IF (row_tiny_diag <= max_report) THEN
           WRITE(*,'(A,I6,1X,A,1PE12.4)') 'TINY_DIAG_PAR: row=', i, ' |diag|=', mag_diag
        END IF
     END IF

  END DO

  WRITE(*,'(A,I12)') 'COMPARE_CSR: entries only in serial   = ', count_only1
  WRITE(*,'(A,I12)') 'COMPARE_CSR: entries only in parallel = ', count_only2
  WRITE(*,'(A,I12)') 'COMPARE_CSR: entries with significant diff = ', count_diff
  WRITE(*,'(A,1PE12.4)') 'COMPARE_CSR: max absolute difference = ', max_diff
  WRITE(*,'(A,I12)') 'COMPARE_CSR: rows missing diagonal (parallel) = ', row_missing_diag
  WRITE(*,'(A,I12)') 'COMPARE_CSR: rows with tiny diagonal (parallel) = ', row_tiny_diag

END SUBROUTINE COMPARE_CSR
SUBROUTINE COMPARE_CSR_PATTERN(N, ia_ref, ja_ref, nnz_ref, ia_new, ja_new, nnz_new, my_rank, IK)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N, nnz_ref, nnz_new, my_rank, IK
   INTEGER, INTENT(IN) :: ia_ref(:), ja_ref(:)
   INTEGER, INTENT(IN) :: ia_new(:), ja_new(:)
   INTEGER :: i, k

   ! quick nnz check
   IF (nnz_ref /= nnz_new) THEN
      IF (my_rank == 0) THEN
         WRITE(*,*) 'CSR pattern change at IK=', IK, ' : nnz_ref=', nnz_ref, ' nnz_new=', nnz_new
      END IF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, i)
   END IF

   ! check row pointers
   DO i = 1, N+1
      IF (ia_ref(i) /= ia_new(i)) THEN
         IF (my_rank == 0) THEN
            WRITE(*,*) 'CSR row pointer mismatch at IK=', IK, ' row=', i, &
                       ' ref=', ia_ref(i), ' new=', ia_new(i)
         END IF
         CALL MPI_Abort(MPI_COMM_WORLD, 1, i)
      END IF
   END DO

   ! check column indices
   DO k = 1, nnz_ref
      IF (ja_ref(k) /= ja_new(k)) THEN
         IF (my_rank == 0) THEN
            WRITE(*,*) 'CSR column index mismatch at IK=', IK, ' pos=', k, &
                       ' ref=', ja_ref(k), ' new=', ja_new(k)
         END IF
         CALL MPI_Abort(MPI_COMM_WORLD, 1, i)
      END IF
   END DO

END SUBROUTINE COMPARE_CSR_PATTERN

!original CA subroutine 

   SUBROUTINE CA(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                 CI, NORD, AS, WT, IE0, IS0, DZ0, N, LDAB, KL, KU, A)

      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER, INTENT(IN) :: KL, KU, LDAB, N
      INTEGER, INTENT(IN) :: NTO, NX, NZ, NORD, IANISO, IE0, IS0, NNX, NNZ

      real(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FK
      real(dp), INTENT(IN) :: CR(:, :), CI(:, :), FREQ, DZ0
      COMPLEX(dp), INTENT(INOUT) :: A(:, :)

      real(dp), ALLOCATABLE :: Z1(:), Z2(:), T1(:), T2(:)
      real(dp), ALLOCATABLE :: NOX(:), NOZ(:), DLX(:), DLZ(:)
      real(dp), ALLOCATABLE :: DNX(:), DNZ(:), XP(:)
      INTEGER, ALLOCATABLE :: INDX(:), INDZ(:)
      COMPLEX(dp), ALLOCATABLE :: P(:), Q(:), CJ(:)
      INTEGER :: I, J, K, L, N0, MM, I0, II, KK, L1, IX, IZ, ID, IP, IQ, JW
      INTEGER :: K1, NPT
      real(dp) :: Z0, OMIGA, X1, X2, DX, C1, XI, DM, S, AK, DZ, C2, C3, BL
      COMPLEX(dp) :: OMIGA2, RHO, WIJ, AIJ, JJ
! write(*,*)"---- initialisation --------------"
      MM = 2*(NORD - 1) + 1
      I0 = KL + KU + 1
      Z0 = IE0*DZ0
      OMIGA = 2_dp*PI*FREQ
      OMIGA2 = DCMPLX(OMIGA*OMIGA, 0.D0)

      ALLOCATE (Z1(NORD), Z2(NORD), T1(NORD), T2(NORD))
      ALLOCATE (NOX(NORD), NOZ(NORD), INDX(MM), INDZ(NORD))
      ALLOCATE (DLX(NORD), DLZ(NORD), DNX(MM), DNZ(NORD))
      ALLOCATE (XP(NORD), P(21), Q(81), CJ(27))

      A(1:LDAB, 1:N) = (0.D0, 0.D0)

      DO I = 1, NX - 1
         X1 = X(I)
         X2 = X(I + 1)
         DX = X2 - X1
         C1 = 2.D0/DX

!loop for [Zj,Zj+1]
         DO J = 1, NZ - 1
            N0 = (I - 1)*(NORD - 1)*NNZ + (J - 1)*(NORD - 1) + 1

            !find Z1(x) & Z2(x)
            IF (J .EQ. 1) THEN
            DO K = 1, NORD
               Z1(K) = 0.D0
               Z2(K) = Z1(K) + DZ0
               T1(K) = 0.D0
               T2(K) = 0.D0
            END DO
            GO TO 2
            END IF

            IF ((J .GT. 1) .AND. (J .LE. IE0)) THEN
            DO K = 1, NORD
               Z1(K) = Z2(K)
               Z2(K) = Z1(K) + DZ0
               T1(K) = 0.D0
               T2(K) = 0.D0
            END DO
            GO TO 2
            END IF

!find Z1'(x) & Z2'(x)
            DO K = 1, NORD
               XI = 0.5D0*(X2 - X1)*AS(K) + 0.5D0*(X1 + X2)
               XP(K) = XI
               DM = (ZH(NTO, XTO, ZTO, XI) - Z0)/DBLE(FLOAT((NZ - 1) - IE0))
               Z1(K) = Z2(K)
               T1(K) = T2(K)
               Z2(K) = Z1(K) + DM
            END DO

            DO K = 1, NORD
               XI = XP(K)
               CALL CDLI(XI, NORD, XP, DLX)
               S = 0.D0
               DO L = 1, NORD
                  S = S + DLX(L)*Z2(L)
               END DO
               T2(K) = S
            END DO

!---- Loop for GQ abscissa ---------
2           DO K = 1, NORD
               AK = AS(K)
               DZ = Z2(K) - Z1(K)
               C3 = 2.D0/DZ
               CALL CDLI(AK, NORD, AS, DLX)

               DO L = 1, NORD
                  ID = N0 + (K - 1)*NNZ + (L - 1)

                  !(CR,CI) => RHO,P(21)
                  CALL C21(ID, IANISO, CR, CI, RHO, P)
                  !P(21) => Q(81) (complex values)
                  CALL Q81_NewGSRM(FREQ, NX - 1, NZ - 1, I, J, K, L, NORD, IE0, IS0, RHO, P, Q)

                  !GQ weights at the point
                  WIJ = DCMPLX(0.25D0*DX*DZ*WT(K)*WT(L), 0.D0)
                  BL = AS(L)
                  C2 = -((T2(K) - T1(K))*BL + (T1(K) + T2(K)))/DZ
                  CALL CDLI(BL, NORD, AS, DLZ)

                  !INDX(*) & DNX(*)
                  DO K1 = 1, NORD
                     NOX(K1) = (N0 + (L - 1)) + (K1 - 1)*NNZ ! index for DLX(*)
                     NOZ(K1) = (N0 + (K - 1)*NNZ) + (K1 - 1) ! index for DLZ(*)
                  END DO

                  IX = 0
                  DO K1 = 1, K - 1
                     IX = IX + 1
                     INDX(IX) = NOX(K1)
                     DNX(IX) = C1*DLX(K1)
                  END DO

                  DO L1 = 1, L - 1
                     IX = IX + 1
                     INDX(IX) = NOZ(L1)
                     DNX(IX) = C2*DLZ(L1)
                  END DO

                  IX = IX + 1
                  INDX(IX) = NOX(K) !NOX(K)=NOZ(L)
                  DNX(IX) = C1*DLX(K) + C2*DLZ(L)

                  DO L1 = L + 1, NORD
                     IX = IX + 1
                     INDX(IX) = NOZ(L1)
                     DNX(IX) = C2*DLZ(L1)
                  END DO

                  DO K1 = K + 1, NORD
                     IX = IX + 1
                     INDX(IX) = NOX(K1)
                     DNX(IX) = C1*DLX(K1)
                  END DO

                  !INDZ(*) & DNZ(*)
                  DO L1 = 1, NORD
                     INDZ(L1) = NOZ(L1)
                     DNZ(L1) = C3*DLZ(L1)
                  END DO
                  IZ = NORD
! After INDX/DNX built:
                  NPT = NNX*NNZ
                  do k1 = 1, IX
                     if (INDX(k1) < 1 .or. INDX(k1) > NPT) then
                        write (*, '(A,3I12)') 'CA: bad INDX value:', k1, INDX(k1), NPT
                        stop
                     end if
                  end do
                  do l1 = 1, IZ
                     if (INDZ(l1) < 1 .or. INDZ(l1) > NPT) then
                        write (*, '(A,3I12)') 'CA: bad INDZ value:', l1, INDZ(l1), NPT
                        stop
                     end if
                  end do

                  !---- Wj-loop for 3 weights vectors ---
                  DO JW = 1, 3
                     CALL QCJ(JW, Q, CJ) ! Q(81) => Cj(27)

                     !---- matrix: Aj(*,*) -----------------
                     !c1j11-term: (DxDx)Gx
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 2
                           AIJ = WIJ*CJ(1)*DCMPLX(DNX(IP)*DNX(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c1j13-term: (DxDz)Gx
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 2
                           AIJ = WIJ*CJ(2)*DCMPLX(DNX(IP)*DNZ(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c3j11-term: (DzDx)Gx
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 2
                           AIJ = WIJ*CJ(3)*DCMPLX(DNZ(IP)*DNX(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c3j13-term: (DzDz)Gx
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 2
                           AIJ = WIJ*CJ(4)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c1j12-term: (DxDelta)Gx
                     JJ = 3*ID - 2
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = DCMPLX(0.D0, -FK)*CJ(5)*WIJ*DCMPLX(DNX(IP), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c3j12-term: (DzDetla)Gx
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = DCMPLX(0.D0, -FK)*CJ(6)*WIJ*DCMPLX(DNZ(IP), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j11-term: (DeltaDx)Gx
                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ) - 2
                        AIJ = -DCMPLX(0.D0, -FK)*WIJ*CJ(7)*DCMPLX(DNX(IQ), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j13-term: (DeltaDz)Gx
                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ) - 2
                        AIJ = -DCMPLX(0.D0, -FK)*WIJ*CJ(8)*DCMPLX(DNZ(IQ), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j12-term: (DeltaDelta)Gx
                     JJ = 3*ID - 2
                     IF (II .EQ. JJ) THEN
                        AIJ = DCMPLX(FK*FK, 0.D0)*WIJ*CJ(9)
                        A(I0, JJ) = A(I0, JJ) + AIJ
                     END IF

                     !rho-term:
                     IF (JW .EQ. 1) THEN
                        AIJ = -RHO*OMIGA2*WIJ
                        A(I0, JJ) = A(I0, JJ) + AIJ
                     END IF

                     !---- matrix: Bj(*,*) --------------------
                     !c1j21-term: (DxDx)Gy
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 1
                           AIJ = WIJ*CJ(10)*DCMPLX(DNX(IP)*DNX(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c1j23-term: (DxDz)Gy
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 1
                           AIJ = WIJ*CJ(11)*DCMPLX(DNX(IP)*DNZ(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c3j21-term: (DzDx)Gy
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ) - 1
                           AIJ = WIJ*CJ(12)*DCMPLX(DNZ(IP)*DNX(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c3j23-term: (DzDz)Gy
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ) - 1
                           AIJ = WIJ*CJ(13)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c1j22-term: (DxDelta)Gy
                     JJ = 3*ID - 1
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = DCMPLX(0.D0, -FK)*CJ(14)*WIJ*DCMPLX(DNX(IP), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c3j22-term: (DzDelta)Gy
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = DCMPLX(0.D0, -FK)*CJ(15)*WIJ*DCMPLX(DNZ(IP), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j21-term: (DeltaDx)Gy
                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ) - 1
                        AIJ = -DCMPLX(0.D0, -FK)*CJ(16)*WIJ*DCMPLX(DNX(IQ), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j23-term: (DeltaDz)Gy
                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ) - 1
                        AIJ = -DCMPLX(0.D0, -FK)*CJ(17)*WIJ*DCMPLX(DNZ(IQ), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j22-term: (DeltaDelta)Gy
                     JJ = 3*ID - 1
                     IF (JJ .EQ. II) THEN
                        AIJ = WIJ*DCMPLX(FK*FK, 0.D0)*CJ(18)
                        A(I0, JJ) = A(I0, JJ) + AIJ
                     END IF

                     !RHO-term:
                     IF (JW .EQ. 2) THEN
                        AIJ = -WIJ*RHO*OMIGA2
                        A(I0, JJ) = A(I0, JJ) + AIJ
                     END IF

                     !---- matrix: Cj(*,*) ----------
                     !c1j31-term: (DxDx)Gz
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ)
                           AIJ = WIJ*CJ(19)*DCMPLX(DNX(IP)*DNX(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c1j33-term: (DxDz)Gz
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ)
                           AIJ = WIJ*CJ(20)*DCMPLX(DNX(IP)*DNZ(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c3j31-term: (DzDx)Gz
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IX
                           JJ = 3*INDX(IQ)
                           AIJ = WIJ*CJ(21)*DCMPLX(DNZ(IP)*DNX(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c3j33-term: (DzDz)Gz
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        DO IQ = 1, IZ
                           JJ = 3*INDZ(IQ)
                           AIJ = WIJ*CJ(22)*DCMPLX(DNZ(IP)*DNZ(IQ), 0.D0)
                           A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                        END DO
                     END DO

                     !c1j32-term: (DxDelta)Gz
                     JJ = 3*ID
                     DO IP = 1, IX
                        II = 3*INDX(IP) - (3 - JW)
                        AIJ = WIJ*DCMPLX(0.D0, -FK)*CJ(23)*DCMPLX(DNX(IP), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c3j32-term: (DzDelta)Gz
                     DO IP = 1, IZ
                        II = 3*INDZ(IP) - (3 - JW)
                        AIJ = WIJ*DCMPLX(0.D0, -FK)*CJ(24)*DCMPLX(DNZ(IP), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j31-term: (DeltaDx)Gz
                     II = 3*ID - (3 - JW)
                     DO IQ = 1, IX
                        JJ = 3*INDX(IQ)
                        AIJ = -WIJ*DCMPLX(0.D0, -FK)*CJ(25)*DCMPLX(DNX(IQ), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !c2j33-term: (DeltaDz)Gz
                     DO IQ = 1, IZ
                        JJ = 3*INDZ(IQ)
                        AIJ = -WIJ*DCMPLX(0.D0, -FK)*CJ(26)*DCMPLX(DNZ(IQ), 0.D0)
                        A(I0 + II - JJ, JJ) = A(I0 + II - JJ, JJ) + AIJ
                     END DO

                     !(c2j32,rho)-term: (DeltaDelta)Gz
                     JJ = 3*ID
                     IF (JJ .EQ. II) THEN
                        AIJ = WIJ*CJ(27)*DCMPLX(FK*FK, 0.D0)
                        A(I0, JJ) = A(I0, JJ) + AIJ
                     END IF

                     !RHO-term:
                     IF (JW .EQ. 3) THEN
                        AIJ = -WIJ*RHO*OMIGA2
                        A(I0, JJ) = A(I0, JJ) + AIJ
                     END IF

                  END DO ! (weights'loop end)
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (Z1, Z2, T1, T2, DLX, DLZ, DNX, DNZ)
      DEALLOCATE (XP, P, NOX, NOZ, INDX, INDZ, Q, CJ)

      RETURN
   END SUBROUTINE CA
    
    end module stiffness_assembly_mod
