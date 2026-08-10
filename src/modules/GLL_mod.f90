module GLL_mod
   ! Contains all Gauss-Lobatto-Legendre nodes and quadrature weights computations subroutine
   use iso_fortran_env, only: dp => real64
   implicit none
   private

   public :: GLL

contains
!-----------------------------------------------------------------------
!  GAMMAF computes the Gamma function for specific values of X.
!
!  Inputs:
!    X ............... Input value for which the Gamma function is computed
!
!  Output:
!    GAMMAF .......... Value of the Gamma function at X
!
!  Notes:
!    - Provides fixed known values for standard arguments used by
!      Gauss–Jacobi routines.
!    - This is a BLAS-standard compatibility function; logic unchanged.
!-----------------------------------------------------------------------
   FUNCTION GAMMAF(X) RESULT(g)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 REAL(dp), INTENT(IN) :: X
                                 REAL(dp) :: g

                                 ! Local constants
                                 REAL(dp), PARAMETER :: ZERO = 0.0_dp
                                 REAL(dp), PARAMETER :: HALF = 0.5_dp
                                 REAL(dp), PARAMETER :: ONE = 1.0_dp
                                 REAL(dp), PARAMETER :: TWO = 2.0_dp
                                 REAL(dp), PARAMETER :: FOUR = 4.0_dp
                                 REAL(dp) :: PI

                                 PI = FOUR*ATAN(ONE)
                                 g = ONE

                                 IF (X == -HALF) g = -TWO*SQRT(PI)
                                 IF (X == HALF) g = SQRT(PI)
                                 IF (X == ONE) g = ONE
                                 IF (X == TWO) g = ONE
                                 IF (X == 1.5_dp) g = SQRT(PI)/2.0_dp
                                 IF (X == 2.5_dp) g = 1.5_dp*SQRT(PI)/2.0_dp
                                 IF (X == 3.5_dp) g = 2.5_dp*1.5_dp*SQRT(PI)/2.0_dp
                                 IF (X == 3.0_dp) g = 2.0_dp
                                 IF (X == 4.0_dp) g = 6.0_dp
                                 IF (X == 5.0_dp) g = 24.0_dp
                                 IF (X == 6.0_dp) g = 120.0_dp
                              END FUNCTION GAMMAF

!-----------------------------------------------------------------------
!
!     GLL calculates the Gauss-Lobatto-Legendre points and weights
!     for integration over the interval [-1, 1].
!
!     Inputs:
!       NP................ Number of points (NP-1 is the polynomial degree)
!
!     Outputs:
!       Z(NP)............. Gauss-Lobatto-Legendre points
!       W(NP)............. Weights for integration

!     this subroutine involves the following routines:
!            ZWGJD(), JACOBF(), ENDW1(), ENDW2(),
!            GAMMAF(), JACG(), PNORMJ().
!
!-----------------------------------------------------------------------
                              SUBROUTINE GLL(NP, Z, W)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN)  :: NP
                                 REAL(dp), INTENT(OUT) :: Z(NP), W(NP)

                                 ! Locals
                                 INTEGER  :: N, NM1, I
                                 REAL(dp) :: ALPHA, BETA, ALPG, BETG
                                 REAL(dp) :: ONE, TWO
                                 REAL(dp) :: P, PD, PM1, PDM1, PM2, PDM2

                                 ONE = 1.0_dp
                                 TWO = 2.0_dp
                                 ALPHA = 0.0_dp
                                 BETA = 0.0_dp

                                 ! Polynomial degree
                                 N = NP - 1
                                 NM1 = N - 1

                                 ! Basic validation
                                 IF (NP < 2) ERROR STOP 'GLL: minimum number of Gauss-Lobatto points is 2.'
                                 IF (ALPHA <= -ONE .OR. BETA <= -ONE) ERROR STOP 'GLL: alpha,beta must be > -1.'

                                 ! Interior nodes/weights for Jacobi(alpha+1,beta+1)
                                 IF (NM1 > 0) THEN
                                    ALPG = ALPHA + ONE
                                    BETG = BETA + ONE
                                    CALL ZWGJD(Z(2), W(2), NM1, ALPG, BETG)
                                 END IF

                                 ! Endpoints
                                 Z(1) = -ONE
                                 Z(NP) = ONE

                                 ! Convert interior weights to Lobatto weights
                                 DO I = 2, NP - 1
                                    W(I) = W(I)/(ONE - Z(I)*Z(I))
                                 END DO

                                 ! Endpoint weights via derivative of P_n at ±1
                                 CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(1))
                                 W(1) = ENDW1(N, ALPHA, BETA)/(TWO*PD)

                                 CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(NP))
                                 W(NP) = ENDW2(N, ALPHA, BETA)/(TWO*PD)
                              END SUBROUTINE GLL

                              !-----------------------------------------------------------------------
                              !
                              !     ZWGJD calculates the points and weights for the Gauss-Jacobi
                              !     quadrature over the interval [-1, 1].
                              !
                              !     Inputs:
                              !       NP................ Number of points (NP-1 is the polynomial degree)
                              !       ALPHA............. Alpha parameter for the Jacobi polynomial
                              !       BETA.............. Beta parameter for the Jacobi polynomial
                              !
                              !     Outputs:
                              !       Z(NP)............. Gauss-Jacobi points
                              !       W(NP)............. Weights for integration
                              !
                              !-----------------------------------------------------------------------
                              SUBROUTINE ZWGJD(Z, W, NP, ALPHA, BETA)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN)  :: NP
                                 REAL(dp), INTENT(IN)  :: ALPHA, BETA
                                 REAL(dp), INTENT(OUT) :: Z(NP), W(NP)

                                 ! Locals
                                 INTEGER  :: N, NP1, NP2, I
                                 REAL(dp) :: ONE, TWO
                                 REAL(dp) :: DN, DNP1, DNP2
                                 REAL(dp) :: APB, FAC1, FAC2, FAC3
                                 REAL(dp) :: FNORM, RCOEF
                                 REAL(dp) :: P, PD, PM1, PDM1, PM2, PDM2

                                 ONE = 1.0_dp
                                 TWO = 2.0_dp

                                 ! Basic validation
                                 IF (NP < 1) ERROR STOP 'ZWGJD: minimum number of Gauss-Jacobi points is 1.'
                                 IF (ALPHA <= -ONE .OR. BETA <= -ONE) ERROR STOP 'ZWGJD: alpha,beta must be > -1.'

                                 N = NP - 1
                                 DN = REAL(N, dp)
                                 APB = ALPHA + BETA

                                 ! Special case: single point
                                 IF (NP == 1) THEN
                                    Z(1) = (BETA - ALPHA)/(APB + TWO)
                                    W(1) = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + ONE)/GAMMAF(APB + TWO)*TWO**(APB + ONE)
                                    RETURN
                                 END IF

                                 ! Compute Jacobi roots (nodes)
                                 CALL JACG(Z, NP, ALPHA, BETA)

                                 NP1 = N + 1
                                 NP2 = N + 2
                                 DNP1 = REAL(NP1, dp)
                                 DNP2 = REAL(NP2, dp)

                                 FAC1 = DNP1 + ALPHA + BETA + ONE
                                 FAC2 = FAC1 + DNP1
                                 FAC3 = FAC2 + ONE
                                 FNORM = PNORMJ(NP1, ALPHA, BETA)

                                 RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)

                                 ! Compute weights
                                 DO I = 1, NP
                                    CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, NP2, ALPHA, BETA, Z(I))
                                    W(I) = -RCOEF/(P*PDM1)
                                 END DO
                              END SUBROUTINE ZWGJD

!-----------------------------------------------------------------------
!  JACOBF computes the values of Jacobi polynomials and their derivatives
!  at a given point.
!
!  Inputs:
!    N ............... Degree of the Jacobi polynomial
!    ALP, BET ........ Alpha and Beta parameters for the Jacobi polynomial
!    X ............... Point at which the polynomial is evaluated
!
!  Outputs:
!    POLY ............ Value of the Jacobi polynomial at X
!    PDER ............ Derivative of the Jacobi polynomial at X
!    POLYM1, PDERM1 .. Values and derivatives of the previous polynomial
!    POLYM2, PDERM2 .. Values and derivatives of the second previous polynomial
!
!  Notes:
!    - Based on the standard three-term recurrence relation for Jacobi
!      polynomials (see Abramowitz & Stegun, §22.7).
!    - Handles N = 0 and N = 1 explicitly for numerical stability.
!-----------------------------------------------------------------------
                              SUBROUTINE JACOBF(POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2, &
                                                N, ALP, BET, X)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE

                                 ! Arguments
                                 INTEGER, INTENT(IN)    :: N
                                 REAL(dp), INTENT(IN)    :: ALP, BET, X
                                 REAL(dp), INTENT(OUT)   :: POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2

                                 ! Locals
                                 REAL(dp) :: APB, POLYL, PDERL, POLYN, PDERN, PSAVE, PDSAVE
                                 REAL(dp) :: DK, A1, A2, A3, A4, B3
                                 INTEGER  :: K
                                 REAL(dp), PARAMETER :: ONE = 1.0_dp, TWO = 2.0_dp

                                 !--------------------------------------------------------------------
                                 ! Initialization and base cases
                                 !--------------------------------------------------------------------
                                 APB = ALP + BET
                                 POLY = ONE
                                 PDER = 0.0_dp

                                 IF (N == 0) THEN
                                    POLYM1 = POLY
                                    PDERM1 = PDER
                                    POLYM2 = 0.0_dp
                                    PDERM2 = 0.0_dp
                                    RETURN
                                 END IF

                                 POLYL = POLY
                                 PDERL = PDER

                                 ! First-degree polynomial
                                 POLY = (ALP - BET + (APB + TWO)*X)/TWO
                                 PDER = (APB + TWO)/TWO

                                 IF (N == 1) THEN
                                    POLYM1 = POLYL
                                    PDERM1 = PDERL
                                    POLYM2 = 0.0_dp
                                    PDERM2 = 0.0_dp
                                    RETURN
                                 END IF

                                 !--------------------------------------------------------------------
                                 ! Recurrence for N >= 2
                                 !--------------------------------------------------------------------
                                 DO K = 2, N
                                    DK = REAL(K, dp)

                                    A1 = TWO*DK*(DK + APB)*(TWO*DK + APB - TWO)
                                    A2 = (TWO*DK + APB - ONE)*(ALP**2 - BET**2)
                                    B3 = (TWO*DK + APB - TWO)
                                    A3 = B3*(B3 + ONE)*(B3 + TWO)
                                    A4 = TWO*(DK + ALP - ONE)*(DK + BET - ONE)*(TWO*DK + APB)

                                    POLYN = ((A2 + A3*X)*POLY - A4*POLYL)/A1
                                    PDERN = ((A2 + A3*X)*PDER - A4*PDERL + A3*POLY)/A1

                                    ! Shift previous polynomial states
                                    PSAVE = POLYL
                                    PDSAVE = PDERL
                                    POLYL = POLY
                                    POLY = POLYN
                                    PDERL = PDER
                                    PDER = PDERN
                                 END DO

                                 !--------------------------------------------------------------------
                                 ! Return previous two polynomials and their derivatives
                                 !--------------------------------------------------------------------
                                 POLYM1 = POLYL
                                 PDERM1 = PDERL
                                 POLYM2 = PSAVE
                                 PDERM2 = PDSAVE
                              END SUBROUTINE JACOBF

!-----------------------------------------------------------------------
!  ENDW1 computes the weight for the first Gauss–Jacobi quadrature point.
!
!  Inputs:
!    N ............... Degree of the Jacobi polynomial
!    ALPHA, BETA ..... Alpha and Beta parameters for the Jacobi polynomial
!
!  Output:
!    ENDW1 ........... Weight for the first quadrature point
!
!  Notes:
!    - Uses recurrence formulas for stability with increasing N.
!    - References: Abramowitz & Stegun, §22.7
!-----------------------------------------------------------------------
                              FUNCTION ENDW1(N, ALPHA, BETA) RESULT(w)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN) :: N
                                 REAL(dp), INTENT(IN) :: ALPHA, BETA
                                 REAL(dp)             :: w

                                 ! Uses module procedure GAMMAF (defined below); no local declaration

                                 ! Locals
                                 INTEGER :: I
                                 REAL(dp) :: ZERO, ONE, TWO, THREE, FOUR
                                 REAL(dp) :: APB, F1, F2, F3, FINT1, FINT2
                                 REAL(dp) :: A1, A2, A3, DI, ABN, ABNN

                                 ZERO = 0.0_dp
                                 ONE = 1.0_dp
                                 TWO = 2.0_dp
                                 THREE = 3.0_dp
                                 FOUR = 4.0_dp

                                 APB = ALPHA + BETA

                                 IF (N == 0) THEN
                                    w = ZERO
                                    RETURN
                                 END IF

                                 F1 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + ONE)/GAMMAF(APB + THREE)
                                 F1 = F1*(APB + TWO)*TWO**(APB + TWO)/TWO
                                 IF (N == 1) THEN
                                    w = F1
                                    RETURN
                                 END IF

                                 FINT1 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + ONE)/GAMMAF(APB + THREE)
                                 FINT1 = FINT1*TWO**(APB + TWO)
                                 FINT2 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + TWO)/GAMMAF(APB + FOUR)
                                 FINT2 = FINT2*TWO**(APB + THREE)
                                 F2 = (-TWO*(BETA + TWO)*FINT1 + (APB + FOUR)*FINT2)*(APB + THREE)/FOUR
                                 IF (N == 2) THEN
                                    w = F2
                                    RETURN
                                 END IF

                                 DO I = 3, N
                                    DI = REAL(I - 1, dp)
                                    ABN = ALPHA + BETA + DI
                                    ABNN = ABN + DI
                                    A1 = -TWO*(DI + ALPHA)*(DI + BETA)/(ABN*ABNN*(ABNN + ONE))
                                    A2 = TWO*(ALPHA - BETA)/(ABNN*(ABNN + TWO))
                                    A3 = TWO*(ABN + ONE)/((ABNN + TWO)*(ABNN + ONE))
                                    F3 = -(A2*F2 + A1*F1)/A3
                                    F1 = F2
                                    F2 = F3
                                 END DO

                                 w = F3
                              END FUNCTION ENDW1

!-----------------------------------------------------------------------
!  ENDW2 computes the weight for the last Gauss–Jacobi quadrature point.
!
!  Inputs:
!    N ............... Degree of the Jacobi polynomial
!    ALPHA, BETA ..... Alpha and Beta parameters for the Jacobi polynomial
!
!  Output:
!    ENDW2 ........... Weight for the last quadrature point
!
!  Notes:
!    - Same recurrence as ENDW1 but with reversed symmetry (Alpha ↔ Beta).
!    - References: Abramowitz & Stegun, §22.7
!-----------------------------------------------------------------------
                              FUNCTION ENDW2(N, ALPHA, BETA) RESULT(w)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN) :: N
                                 REAL(dp), INTENT(IN) :: ALPHA, BETA
                                 REAL(dp)             :: w

                                 ! Uses module procedure GAMMAF (defined below); no local declaration

                                 ! Locals
                                 INTEGER :: I
                                 REAL(dp) :: ZERO, ONE, TWO, THREE, FOUR
                                 REAL(dp) :: APB, F1, F2, F3, FINT1, FINT2
                                 REAL(dp) :: A1, A2, A3, DI, ABN, ABNN

                                 ZERO = 0.0_dp
                                 ONE = 1.0_dp
                                 TWO = 2.0_dp
                                 THREE = 3.0_dp
                                 FOUR = 4.0_dp

                                 APB = ALPHA + BETA

                                 IF (N == 0) THEN
                                    w = ZERO
                                    RETURN
                                 END IF

                                 F1 = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + TWO)/GAMMAF(APB + THREE)
                                 F1 = F1*(APB + TWO)*TWO**(APB + TWO)/TWO
                                 IF (N == 1) THEN
                                    w = F1
                                    RETURN
                                 END IF

                                 FINT1 = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + TWO)/GAMMAF(APB + THREE)
                                 FINT1 = FINT1*TWO**(APB + TWO)
                                 FINT2 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + TWO)/GAMMAF(APB + FOUR)
                                 FINT2 = FINT2*TWO**(APB + THREE)
                                 F2 = (TWO*(ALPHA + TWO)*FINT1 - (APB + FOUR)*FINT2)*(APB + THREE)/FOUR
                                 IF (N == 2) THEN
                                    w = F2
                                    RETURN
                                 END IF

                                 DO I = 3, N
                                    DI = REAL(I - 1, dp)
                                    ABN = ALPHA + BETA + DI
                                    ABNN = ABN + DI
                                    A1 = -TWO*(DI + ALPHA)*(DI + BETA)/(ABN*ABNN*(ABNN + ONE))
                                    A2 = TWO*(ALPHA - BETA)/(ABNN*(ABNN + TWO))
                                    A3 = TWO*(ABN + ONE)/((ABNN + TWO)*(ABNN + ONE))
                                    F3 = -(A2*F2 + A1*F1)/A3
                                    F1 = F2
                                    F2 = F3
                                 END DO

                                 w = F3
                              END FUNCTION ENDW2

!-----------------------------------------------------------------------
!  JACG computes the roots of the Jacobi polynomial, which are
!  used as the Gauss–Jacobi quadrature points.
!
!  Inputs:
!    NP ............... Number of points (NP-1 is the polynomial degree)
!    ALPHA, BETA ...... Alpha and Beta parameters for the Jacobi polynomial
!
!  Outputs:
!    Z(NP) ............ Roots of the Jacobi polynomial
!-----------------------------------------------------------------------
                              SUBROUTINE JACG(Z, NP, ALPHA, BETA)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN)  :: NP
                                 REAL(dp), INTENT(IN)  :: ALPHA, BETA
                                 REAL(dp), INTENT(OUT) :: Z(NP)

                                 ! Local variables
                                 INTEGER  :: KSTOP, N, J, K, JM, I, JMIN
                                 REAL(dp) :: EPS, DTH, X, X1, X2, XLAST, DELX, RECSUM
                                 REAL(dp) :: P, PD, PM1, PDM1, PM2, PDM2
                                 REAL(dp) :: XMIN, SWAP
                                 REAL(dp), PARAMETER :: PI = 4.0_dp*ATAN(1.0_dp)

                                 KSTOP = 10
                                 EPS = 1.0e-12_dp
                                 N = NP - 1
                                 DTH = PI/(2.0_dp*REAL(N, dp) + 2.0_dp)

                                 DO J = 1, NP
                                    IF (J == 1) THEN
                                       X = COS((2.0_dp*(REAL(J, dp) - 1.0_dp) + 1.0_dp)*DTH)
                                    ELSE
                                       X1 = COS((2.0_dp*(REAL(J, dp) - 1.0_dp) + 1.0_dp)*DTH)
                                       X2 = XLAST
                                       X = 0.5_dp*(X1 + X2)
                                    END IF

                                    DO K = 1, KSTOP
                                       CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, X)
                                       RECSUM = 0.0_dp
                                       JM = J - 1
                                       IF (JM > 0) THEN
                                          DO I = 1, JM
                                             RECSUM = RECSUM + 1.0_dp/(X - Z(NP - I + 1))
                                          END DO
                                       END IF
                                       DELX = -P/(PD - RECSUM*P)
                                       X = X + DELX
                                       IF (ABS(DELX) < EPS) EXIT
                                    END DO

                                    Z(NP - J + 1) = X
                                    XLAST = X
                                 END DO

                                 ! Sort Z in ascending order
                                 DO I = 1, NP
                                    XMIN = 2.0_dp
                                    JMIN = I
                                    DO J = I, NP
                                       IF (Z(J) < XMIN) THEN
                                          XMIN = Z(J)
                                          JMIN = J
                                       END IF
                                    END DO
                                    IF (JMIN /= I) THEN
                                       SWAP = Z(I)
                                       Z(I) = Z(JMIN)
                                       Z(JMIN) = SWAP
                                    END IF
                                 END DO
                              END SUBROUTINE JACG

!-----------------------------------------------------------------------
!  PNORMJ computes the normalization constant for the Jacobi
!  polynomial of degree N.
!
!  Inputs:
!    N ............... Degree of the Jacobi polynomial
!    ALPHA, BETA ..... Alpha and Beta parameters for the Jacobi polynomial
!
!  Output:
!    PNORMJ .......... Normalization constant for the Jacobi polynomial
!-----------------------------------------------------------------------
                              FUNCTION PNORMJ(N, ALPHA, BETA) RESULT(pnorm)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN) :: N
                                 REAL(dp), INTENT(IN) :: ALPHA, BETA

                                 REAL(dp)             :: pnorm
                                 ! Uses module procedure GAMMAF (defined below); no local declaration
                                 ! Locals
                                 INTEGER  :: I
                                 REAL(dp) :: ONE, TWO, DN, CONST, PROD, DINDX, FRAC

                                 ONE = 1.0_dp
                                 TWO = 2.0_dp
                                 DN = REAL(N, dp)
                                 CONST = ALPHA + BETA + ONE

                                 IF (N <= 1) THEN
                                    PROD = GAMMAF(DN + ALPHA)*GAMMAF(DN + BETA)
                                    PROD = PROD/(GAMMAF(DN)*GAMMAF(DN + ALPHA + BETA))
                                    pnorm = PROD*TWO**CONST/(TWO*DN + CONST)
                                    RETURN
                                 END IF

                                 PROD = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + ONE)
                                 PROD = PROD/(TWO*(ONE + CONST)*GAMMAF(CONST + ONE))
                                 PROD = PROD*(ONE + ALPHA)*(TWO + ALPHA)
                                 PROD = PROD*(ONE + BETA)*(TWO + BETA)

                                 DO I = 3, N
                                    DINDX = REAL(I, dp)
                                    FRAC = (DINDX + ALPHA)*(DINDX + BETA)/(DINDX*(DINDX + ALPHA + BETA))
                                    PROD = PROD*FRAC
                                 END DO

                                 pnorm = PROD*TWO**CONST/(TWO*DN + CONST)
                              END FUNCTION PNORMJ

end module GLL_mod

