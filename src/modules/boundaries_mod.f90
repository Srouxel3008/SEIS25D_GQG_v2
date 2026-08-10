module boundaries_mod
   use iso_fortran_env, only: dp => real64
   USE constant_mod,  ONLY: pi, deg2rad
   implicit NONE

contains
  !!----------------------------------------------------------------------C
   !This module contains all subroutines related to the absorbing boundaries
   !of the model domain.
   !1.  Q81_NewGSRM: Convert 21 elastic moduli into 81 COMPLEX(dp) elastic moduli
   !2.  QCJ: Extract a directional stiffness submatrix from the full stiffness tensor
   !4.  C21: Convert 21 elastic moduli into 81 COMPLEX(dp) elastic moduli

   SUBROUTINE QCJ(J, Q, CJ)
      !----------------------------------------------------------------------
      !
      !     SUBROUTINE QCJ extracts a directional stiffness submatrix CJ(27)
      !     from the full complex stiffness tensor Q(81) for use in forward
      !     modeling of wave propagation in general anisotropic media.
      !
      !     In 3D elastic modeling, the stiffness tensor c_{ijkl} (81 components
      !     in Voigt notation) is directionally projected depending on the wave
      !     propagation direction J:
      !
      !       J = 1   → Wave propagates along x-direction
      !       J = 2   → Wave propagates along y-direction
      !       J = 3   → Wave propagates along z-direction
      !
      !     This subroutine selects and arranges 27 complex stiffness elements
      !     from Q into CJ, corresponding to the stress-strain relationship
      !     for the given propagation direction.
      !
      !     Entries:
      !       J......................Wave propagation direction (1: x, 2: y, 3: z)
      !       Q(81)..................Full complex stiffness tensor in Voigt format
      !
      !     Return:
      !       CJ(27).................Directional stiffness submatrix used to solve
      !                              the wave equation along the selected direction
      !
      !     Notes:
      !       - The mapping follows internal conventions specific to the
      !         numerical formulation used in this modeling code.
      !       - If J is not 1, 2, or 3, the subroutine prints an error and stops.
      !
      !----------------------------------------------------------------------

      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: J
      COMPLEX(dp), INTENT(IN) :: Q(81)
      COMPLEX(dp), INTENT(OUT) :: CJ(27)

      !---- J=1: (W,0,0) -----------
      IF (J .EQ. 1) THEN
         CJ(1) = Q(1)   !c1111=q( 1)
         CJ(2) = Q(10)  !c1113=q(10)
         CJ(3) = Q(13)  !c3111=q(13)
         CJ(4) = Q(68)  !c3113=q(68)
         CJ(5) = Q(14)  !c1112=q(14)
         CJ(6) = Q(74)  !c3112=q(74)
         CJ(7) = Q(17)  !c2111=q(17)
         CJ(8) = Q(76)  !c2113=q(76))
         CJ(9) = Q(80)  !c2112=q(80)

         CJ(10) = Q(15) !c1121=q(15)
         CJ(11) = Q(6) !c1123=q( 6)
         CJ(12) = Q(77) !c3121=q(77)
         CJ(13) = Q(54) !c3123=q(54)
         CJ(14) = Q(2) !c1122=q( 2)
         CJ(15) = Q(28) !c3122=q(28)
         CJ(16) = Q(81) !c2121=q(81)
         CJ(17) = Q(60) !c2123=q(60)
         CJ(18) = Q(30) !c2122=q(30)

         CJ(19) = Q(11) !c1131=q(11)
         CJ(20) = Q(4) !c1133=q( 4)
         CJ(21) = Q(69) !c3131=q(69)
         CJ(22) = Q(39) !c3133=q(39)
         CJ(23) = Q(7) !c1132=q( 7)
         CJ(24) = Q(55) !c3132=q(55)
         CJ(25) = Q(75) !c2131=q(75)
         CJ(26) = Q(43) !c2133=q(43)
         CJ(27) = Q(61) !c2132=q(61)
         RETURN
      END IF

      !---- J=2: (0,W,0) ------------
      IF (J .EQ. 2) THEN
         CJ(1) = Q(16)  !c1211=q(16)
         CJ(2) = Q(70)  !c1213=q(70)
         CJ(3) = Q(9)  !c3211=q( 9)
         CJ(4) = Q(56)  !c3213=q(56)
         CJ(5) = Q(78)  !c1212=q(78)
         CJ(6) = Q(64)  !c3212=q(64)
         CJ(7) = Q(3)  !c2211=q( 3)
         CJ(8) = Q(26)  !c2213=q(26)
         CJ(9) = Q(31)  !c2212=q(31)

         CJ(10) = Q(79) !c1221=q(79)
         CJ(11) = Q(58) !c1223=q(58)
         CJ(12) = Q(65) !c3221=q(65)
         CJ(13) = Q(48) !c3223=q(48)
         CJ(14) = Q(29) !c1222=q(29)
         CJ(15) = Q(24) !c3222=q(24)
         CJ(16) = Q(32) !c2221=q(32)
         CJ(17) = Q(21) !c2223=q(21)
         CJ(18) = Q(18) !c2222=q(18)

         CJ(19) = Q(71) !c1231=q(71)
         CJ(20) = Q(42) !c1233=q(42)
         CJ(21) = Q(57) !c3231=q(57)
         CJ(22) = Q(35) !c3233=q(35)
         CJ(23) = Q(59) !c1232=q(59)
         CJ(24) = Q(49) !c3232=q(49)
         CJ(25) = Q(27) !c2231=q(27)
         CJ(26) = Q(19) !c2233=q(19)
         CJ(27) = Q(22) !c2232=q(22)
         RETURN
      END IF

      !---- J=3: (0.0,W) -----------
      IF (J .EQ. 3) THEN
         CJ(1) = Q(12)   !c1311=q(12)
         CJ(2) = Q(66)   !c1313=q(66)
         CJ(3) = Q(5)   !c3311=q( 5)
         CJ(4) = Q(40)   !c3313=q(40)
         CJ(5) = Q(72)   !c1312=q(72)
         CJ(6) = Q(44)   !c3312=q(44)
         CJ(7) = Q(8)   !c2311=q( 8)
         CJ(8) = Q(52)   !c2313=q(52)
         CJ(9) = Q(62)   !c2312=q(62)

         CJ(10) = Q(73)  !c1321=q(73)
         CJ(11) = Q(50)  !c1323=q(50)
         CJ(12) = Q(45)  !c3321=q(45)
         CJ(13) = Q(36)  !c3323=q(36)
         CJ(14) = Q(25)  !c1322=q(25)
         CJ(15) = Q(20)  !c3322=q(20)
         CJ(16) = Q(63)  !c2321=q(63)
         CJ(17) = Q(46)  !c2323=q(49)
         CJ(18) = Q(23)  !c2322=q(23)

         CJ(19) = Q(67)  !c1331=q(67)
         CJ(20) = Q(38)  !c1333=q(38)
         CJ(21) = Q(41)  !c3331=q(41)
         CJ(22) = Q(33)  !c3333=q(33)
         CJ(23) = Q(51)  !c1332=q(51)
         CJ(24) = Q(37)  !c3332=q(37)
         CJ(25) = Q(53)  !c2331=q(53)
         CJ(26) = Q(34)  !c2333=q(34)
         CJ(27) = Q(47)  !c2332=q(47)
         RETURN
      END IF
      WRITE (*, *) '     something wrong with subroutine CJ !'
      STOP
   END

   ! !----------------------------------------------------------------------C
   !                                                                      C
   !     Convert 21 elastic moduli into 81 COMPLEX(dp) elastic moduli      C
   !     for absorbing boundaries using GSRM damping                                                        C
   !                                                                      C
   !      (1) NDX,NDZ.........how many subdomains in x- and z-direction;  C
   !      (2) I,J...........................considered sub-domain index;  C
   !      (3) IE0...........number of the sub-domains in absorbing zone;  C
   !      (4) BI & BJ.................absorbing functions b(x) and b(z);  C
   !      (5) IS0=1 or 2.......No absorbing PML at free-surface or dose;  C
   !      (5) RHO...............................................density;  C
   !      (5) P(21)...............21 complex-valued elastic moduli: Cij;
   !                                                                      C
   !      Returns:                                                        C
   !                                                                      C
   !      (1) RHO...........................................PML-density;  C
   !      (1) Q(81)....................81 COMPLEX(dp) PML-elatsic moduli;  C
   !                                                                      C
   !          Q( 1)=c1111,  Q(21)=C2223,  Q(41)=C3331, Q(61)=C2132,       C
   !          Q( 2)=C1122,  Q(22)=C2232,  Q(42)=C1233, Q(62)=C2312,       C
   !          Q( 3)=C2211,  Q(23)=C2322,  Q(43)=C2133, Q(63)=C2321,       C
   !          Q( 4)=C1133,  Q(24)=C3222,  Q(44)=C3312, Q(64)=C3212,       C
   !          Q( 5)=C3311,  Q(25)=C1322,  Q(45)=C3321, Q(65)=C3221,       C
   !          Q( 6)=C1123,  Q(26)=C2213,  Q(46)=C2323, Q(66)=C1313,       C
   !          Q( 7)=C1132,  Q(27)=C2231,  Q(47)=C2332, Q(67)=C1331,       C
   !          Q( 8)=C2311,  Q(28)=C3122,  Q(48)=C3223, Q(68)=C3113,       C
   !          Q( 9)=C3211,  Q(29)=C1222,  Q(49)=C3232, Q(69)=C3131,       C
   !          Q(10)=C1113,  Q(30)=C2122,  Q(50)=C1323, Q(70)=C1213,       C
   !          Q(11)=C1131,  Q(31)=C2212,  Q(51)=C1332, Q(71)=C1231,       C
   !          Q(12)=C1311,  Q(32)=C2221,  Q(52)=C2313, Q(72)=C1312,       C
   !          Q(13)=C3111,  Q(33)=C3333,  Q(53)=C2331, Q(73)=C1321,       C
   !          Q(14)=C1112,  Q(34)=C2333,  Q(54)=C3123, Q(74)=C3112,       C
   !          Q(15)=C1121,  Q(35)=C3233,  Q(55)=C3132, Q(75)=C2131,       C
   !          Q(16)=C1211,  Q(36)=C3323,  Q(56)=C3231, Q(76)=C2113,       C
   !          Q(17)=C2111,  Q(37)=C3332,  Q(57)=C3213, Q(77)=C3121,       C
   !          Q(18)=C2222,  Q(38)=C1333,  Q(58)=C1223, Q(78)=C1212,       C
   !          Q(19)=C2233,  Q(39)=C3133,  Q(59)=C1232, Q(79)=C1221,       C
   !          Q(20)=C3322,  Q(40)=C3313,  Q(60)=C2123, Q(80)=C2112,       C
   !                                                            Q(81)=C2121,       C
   !                                                                      C
   !----------------------------------------------------------------------C

   SUBROUTINE Q81_NewGSRM(FREQ, NDX, NDZ, I, J, K, L, NORD, IE0, IS0, RHO, P, Q)
   USE iso_fortran_env, ONLY: dp => real64
   IMPLICIT NONE
   REAL(dp),    INTENT(IN)    :: FREQ
   INTEGER,     INTENT(IN)    :: NDX, NDZ, I, J, K, L, NORD, IE0, IS0
   COMPLEX(dp), INTENT(INOUT) :: RHO
   COMPLEX(dp), INTENT(IN)    :: P(21)
   COMPLEX(dp), INTENT(OUT)   :: Q(81)

   ! scalars
   REAL(dp) :: A0, E, DI, DJ, OMIGA, DAMP, REDU, REDUi
   INTEGER  :: MPL, I2, J2, II, JJ
   COMPLEX(dp) :: ONE

   A0    = 2.0_dp
   E     = 1.0e-16_dp
   DI    = 0.0_dp  ! d(x)=0
   DJ    = 0.0_dp  ! d(z)=0
   ONE   = CMPLX(1.0_dp, 0.0_dp, kind=dp)
   OMIGA = 2.0_dp*3.1415926_dp*FREQ

   !---- for d(x) & d(z) -----------------
   MPL = IE0*(NORD - 1) + 1

   ! left side
   IF (I <= IE0) THEN
      II = (I - 1)*(NORD - 1) + K
      DI = A0 * ( REAL(II - MPL, dp) / REAL(1 - MPL, dp) )**2
   END IF

   ! right side
   I2 = NDX - IE0 + 1
   IF (I >= I2) THEN
      II = (I - I2)*(NORD - 1) + K
      DI = A0 * ( REAL(II - 1, dp) / REAL(MPL - 1, dp) )**2
   END IF

   ! bottom
   IF (J <= IE0) THEN
      JJ = (J - 1)*(NORD - 1) + L
      DJ = A0 * ( REAL(JJ - MPL, dp) / REAL(1 - MPL, dp) )**2
   END IF

   ! top (surface), optional
   IF (IS0 == 2) THEN
      J2 = NDZ - IE0 + 1
      IF (J >= J2) THEN
         JJ = (J - J2)*(NORD - 1) + L
         DJ = A0 * ( REAL(JJ - 1, dp) / REAL(MPL - 1, dp) )**2
      END IF
   END IF

   ! damping factors
   DAMP  = 1000.0_dp * SQRT(DI*DI + DJ*DJ)
   REDU  = EXP( -SQRT(DI*DI + DJ*DJ) )
   REDUi = -EXP(  SQRT(DI*DI + DJ*DJ) )
   RHO   = RHO * ( ONE + CMPLX(0.0_dp, -DAMP/OMIGA, kind=dp) )

   !---- complex*16 elastic moduli mapping (logic unchanged) -------
   ! c1jk1:
   Q( 1) = CMPLX(REDU*REAL(P( 1),dp), REDUi*AIMAG(P( 1)), kind=dp) ! C1111
   Q(15) = CMPLX(REDU*REAL(P( 6),dp), REDUi*AIMAG(P( 6)), kind=dp) ! C1121
   Q(11) = CMPLX(REDU*REAL(P( 5),dp), REDUi*AIMAG(P( 5)), kind=dp) ! C1131
   Q(16) = CMPLX(REDU*REAL(P( 6),dp), REDUi*AIMAG(P( 6)), kind=dp) ! C1211
   Q(79) = CMPLX(REDU*REAL(P(21),dp), REDUi*AIMAG(P(21)), kind=dp) ! C1221
   Q(71) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C1231
   Q(12) = CMPLX(REDU*REAL(P( 5),dp), REDUi*AIMAG(P( 5)), kind=dp) ! C1311
   Q(73) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C1321
   Q(67) = CMPLX(REDU*REAL(P(19),dp), REDUi*AIMAG(P(19)), kind=dp) ! C1331

   ! c1jk2:
   Q(14) = CMPLX(REDU*REAL(P( 6),dp), REDUi*AIMAG(P( 6)), kind=dp) ! C1112
   Q( 2) = CMPLX(REDU*REAL(P( 2),dp), REDUi*AIMAG(P( 2)), kind=dp) ! C1122
   Q( 7) = CMPLX(REDU*REAL(P( 4),dp), REDUi*AIMAG(P( 4)), kind=dp) ! C1132
   Q(78) = CMPLX(REDU*REAL(P(21),dp), REDUi*AIMAG(P(21)), kind=dp) ! C1212
   Q(29) = CMPLX(REDU*REAL(P(11),dp), REDUi*AIMAG(P(11)), kind=dp) ! C1222
   Q(59) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C1232
   Q(72) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C1312
   Q(25) = CMPLX(REDU*REAL(P(10),dp), REDUi*AIMAG(P(10)), kind=dp) ! C1322
   Q(51) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C1332

   ! c1jk3:
   Q(10) = CMPLX(REDU*REAL(P( 5),dp), REDUi*AIMAG(P( 5)), kind=dp) ! C1113
   Q( 6) = CMPLX(REDU*REAL(P( 4),dp), REDUi*AIMAG(P( 4)), kind=dp) ! C1123
   Q( 4) = CMPLX(REDU*REAL(P( 3),dp), REDUi*AIMAG(P( 3)), kind=dp) ! C1133
   Q(70) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C1213
   Q(58) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C1223
   Q(42) = CMPLX(REDU*REAL(P(15),dp), REDUi*AIMAG(P(15)), kind=dp) ! C1233
   Q(66) = CMPLX(REDU*REAL(P(19),dp), REDUi*AIMAG(P(19)), kind=dp) ! C1313
   Q(50) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C1323
   Q(38) = CMPLX(REDU*REAL(P(14),dp), REDUi*AIMAG(P(14)), kind=dp) ! C1333

   ! c2jk1:
   Q(17) = CMPLX(REDU*REAL(P( 6),dp), REDUi*AIMAG(P( 6)), kind=dp) ! C2111
   Q(81) = CMPLX(REDU*REAL(P(21),dp), REDUi*AIMAG(P(21)), kind=dp) ! C2121
   Q(75) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C2131
   Q( 3) = CMPLX(REDU*REAL(P( 2),dp), REDUi*AIMAG(P( 2)), kind=dp) ! C2211
   Q(32) = CMPLX(REDU*REAL(P(11),dp), REDUi*AIMAG(P(11)), kind=dp) ! C2221
   Q(27) = CMPLX(REDU*REAL(P(10),dp), REDUi*AIMAG(P(10)), kind=dp) ! C2231
   Q( 8) = CMPLX(REDU*REAL(P( 4),dp), REDUi*AIMAG(P( 4)), kind=dp) ! C2311
   Q(63) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C2321
   Q(53) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C2331

   ! c2jk2:
   Q(80) = CMPLX(REDU*REAL(P(21),dp), REDUi*AIMAG(P(21)), kind=dp) ! C2112
   Q(30) = CMPLX(REDU*REAL(P(11),dp), REDUi*AIMAG(P(11)), kind=dp) ! C2122
   Q(61) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C2132
   Q(31) = CMPLX(REDU*REAL(P(11),dp), REDUi*AIMAG(P(11)), kind=dp) ! C2212
   Q(18) = CMPLX(REDU*REAL(P( 7),dp), REDUi*AIMAG(P( 7)), kind=dp) ! C2222
   Q(22) = CMPLX(REDU*REAL(P( 9),dp), REDUi*AIMAG(P( 9)), kind=dp) ! C2232
   Q(62) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C2312
   Q(23) = CMPLX(REDU*REAL(P( 9),dp), REDUi*AIMAG(P( 9)), kind=dp) ! C2322
   Q(47) = CMPLX(REDU*REAL(P(16),dp), REDUi*AIMAG(P(16)), kind=dp) ! C2332

   ! c2jk3:
   Q(76) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C2113
   Q(60) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C2123
   Q(43) = CMPLX(REDU*REAL(P(15),dp), REDUi*AIMAG(P(15)), kind=dp) ! C2133
   Q(26) = CMPLX(REDU*REAL(P(10),dp), REDUi*AIMAG(P(10)), kind=dp) ! C2213
   Q(21) = CMPLX(REDU*REAL(P( 9),dp), REDUi*AIMAG(P( 9)), kind=dp) ! C2223
   Q(19) = CMPLX(REDU*REAL(P( 8),dp), REDUi*AIMAG(P( 8)), kind=dp) ! C2233
   Q(52) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C2313
   Q(46) = CMPLX(REDU*REAL(P(16),dp), REDUi*AIMAG(P(16)), kind=dp) ! C2323
   Q(34) = CMPLX(REDU*REAL(P(13),dp), REDUi*AIMAG(P(13)), kind=dp) ! C2333

   ! c3jk1:
   Q(13) = CMPLX(REDU*REAL(P( 5),dp), REDUi*AIMAG(P( 5)), kind=dp) ! C3111
   Q(77) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C3121
   Q(69) = CMPLX(REDU*REAL(P(19),dp), REDUi*AIMAG(P(19)), kind=dp) ! C3131
   Q( 9) = CMPLX(REDU*REAL(P( 4),dp), REDUi*AIMAG(P( 4)), kind=dp) ! C3211
   Q(65) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C3221
   Q(57) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C3231
   Q( 5) = CMPLX(REDU*REAL(P( 3),dp), REDUi*AIMAG(P( 3)), kind=dp) ! C3311
   Q(45) = CMPLX(REDU*REAL(P(15),dp), REDUi*AIMAG(P(15)), kind=dp) ! C3321
   Q(41) = CMPLX(REDU*REAL(P(14),dp), REDUi*AIMAG(P(14)), kind=dp) ! C3331

   ! c3jk2:
   Q(74) = CMPLX(REDU*REAL(P(20),dp), REDUi*AIMAG(P(20)), kind=dp) ! C3112
   Q(28) = CMPLX(REDU*REAL(P(10),dp), REDUi*AIMAG(P(10)), kind=dp) ! C3122
   Q(55) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C3132
   Q(64) = CMPLX(REDU*REAL(P(18),dp), REDUi*AIMAG(P(18)), kind=dp) ! C3212
   Q(24) = CMPLX(REDU*REAL(P( 9),dp), REDUi*AIMAG(P( 9)), kind=dp) ! C3222
   Q(49) = CMPLX(REDU*REAL(P(16),dp), REDUi*AIMAG(P(16)), kind=dp) ! C3232
   Q(44) = CMPLX(REDU*REAL(P(15),dp), REDUi*AIMAG(P(15)), kind=dp) ! C3312
   Q(20) = CMPLX(REDU*REAL(P( 8),dp), REDUi*AIMAG(P( 8)), kind=dp) ! C3322
   Q(37) = CMPLX(REDU*REAL(P(13),dp), REDUi*AIMAG(P(13)), kind=dp) ! C3332

   ! c3jk3:
   Q(68) = CMPLX(REDU*REAL(P(19),dp), REDUi*AIMAG(P(19)), kind=dp) ! C3113
   Q(54) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C3123
   Q(39) = CMPLX(REDU*REAL(P(14),dp), REDUi*AIMAG(P(14)), kind=dp) ! C3133
   Q(56) = CMPLX(REDU*REAL(P(17),dp), REDUi*AIMAG(P(17)), kind=dp) ! C3213 (note: comment fixed)
   Q(48) = CMPLX(REDU*REAL(P(16),dp), REDUi*AIMAG(P(16)), kind=dp) ! C3223
   Q(35) = CMPLX(REDU*REAL(P(13),dp), REDUi*AIMAG(P(13)), kind=dp) ! C3233
   Q(40) = CMPLX(REDU*REAL(P(14),dp), REDUi*AIMAG(P(14)), kind=dp) ! C3313
   Q(36) = CMPLX(REDU*REAL(P(13),dp), REDUi*AIMAG(P(13)), kind=dp) ! C3323
   Q(33) = CMPLX(REDU*REAL(P(12),dp), REDUi*AIMAG(P(12)), kind=dp) ! C3333

   RETURN
END SUBROUTINE Q81_NewGSRM


!----------------------------------------------------------------------C
!     Converts model parameters defined by IANISO into a              C
!     full 21-component complex stiffness tensor P(1:21) and          C
!     complex density RHO. The conversion depends on the              C
!     anisotropy type indicated by IANISO, and supports isotropic,    C
!     VTI, TTI, and general cases (both 2D and 3D).                   C
!                                                                      C
!     Entries:                                                         C
!       IP.................Grid index (model point in GQG)            C
!       IANISO.............Number of independent moduli + 1           C
!                          (includes density)                         C
!       CR(IANISO,*).......Real part of input parameters              C
!       CI(IANISO,*).......Imaginary part of input parameters         C
!                                                                      C
!     Returns:                                                         C
!       RHO.................COMPLEX(dp) density value at IP            C
!       P(21)...............COMPLEX(dp) stiffness tensor Cij           C
!                          stored in Voigt notation                   C
!                                                                      C
!     Supported anisotropy types:                                     C
!       IANISO = 3   -> Isotropic (VTI acoustic)                      C
!       IANISO = 6   -> 2D VTI                                        C
!       IANISO = 7   -> 2D TTI                                        C
!       IANISO = 8   -> 3D TTI                                        C
!       IANISO = 10  -> General symmetry with 10 moduli               C
!       IANISO = 22  -> Full 21 independent moduli                    C
!                                                                      C
!----------------------------------------------------------------------C
      SUBROUTINE C21(IP, IANISO, CR, CI, RHO, P)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP, IANISO
      REAL(dp), INTENT(IN) ::  CR(:, :), CI(:, :)
      COMPLEX(dp), INTENT(OUT) :: RHO, P(21)

      INTEGER :: I, J, K, L, I1, J1, K1, L1
      INTEGER :: II, JJ, KK, LL
      INTEGER :: IDELT(3, 3)
      REAL(dp) :: Vp, Vs
      REAL(dp) :: CIM(22)
      REAL(dp) :: ZT0, PH0
      COMPLEX(dp) :: lambda_star, mu_star, Vp_star, Vs_star
      COMPLEX(dp) :: C11, C13, C33, C44, C66, COSZT0, SINZT0, COSZT2, SINZT2
      COMPLEX(dp) :: C1(6, 6), AN(3, 3), Q(6, 6),AA
      


      DATA IDELT(1, 1:3)/1, 0, 0/
      DATA IDELT(2, 1:3)/0, 1, 0/
      DATA IDELT(3, 1:3)/0, 0, 1/

      !write(*,*)'entry variable in C21 IP,IANISO',IP,IANISO

      !---- Pick-up RHO & Setting Im(cij) ------
      RHO = CMPLX(CR(1, IP), CI(1, IP), kind=dp)
      !write(*,*)'C21 RHO return',RHO

      DO I = 2, IANISO
         IF (CI(I, IP) == 0.0_dp) THEN
            CIM(I) = 0.0_dp
         ELSE
            CIM(I) = -CR(I, IP)/CI(I, IP)
         END IF
      END DO
      !---- Isotropic viscoelastic medium: IANISO = 3 -----------------------
!  Uses the RHO and CIM(*) already set above:
!     RHO      = DCMPLX(CR(1,IP), CI(1,IP))
!     CIM(I)   = 0              if CI(I,IP) = 0
!              = -CR(I,IP)/CI(I,IP) otherwise
!
!  Two possible input parameterisations, selected by TH:
!    TH = 0  :  CR(2,IP) = λ_R  (Lamé parameter, Pa or GPa)
!               CR(3,IP) = μ_R  (shear modulus, Pa or GPa)
!               CIM(2)   = λ_I  (e.g. from Q_λ)
!               CIM(3)   = μ_I  (e.g. from Q_μ)
!
!    TH ≠ 0 :  CR(2,IP) = Vp   (P-wave velocity, m/s)
!               CR(3,IP) = Vs   (S-wave velocity, m/s)
!               CIM(2)   ≈ -Vp/Qp
!               CIM(3)   ≈ -Vs/Qs
!               so Vp* = Vp + i*CIM(2), Vs* = Vs + i*CIM(3)
!
!  Complex Lamé parameters:
!      μ*     = μ_R + i μ_I
!      λ*     = λ_R + i λ_I
!
!  Isotropic stiffness (Voigt):
!      c11 = λ* + 2 μ*,   c12 = c13 = λ*,
!      c44 = c55 = c66 = μ*, others = 0.
!----------------------------------------------------------------------
      SELECT CASE (IANISO)
      CASE (3)

         !    IF (TH .EQ. 0) THEN
         !       ! ---- Case 1: input is λ, μ in stiffness units (Pa or GPa) ----
         !       lambda_star = DCMPLX(CR(2,IP), CIM(2))
         !       mu_star     = DCMPLX(CR(3,IP), CIM(3))

         !    ELSE
         ! ---- Case 2: input is Vp, Vs (m/s), Qp, Qs via CIM(2:3) ----
         Vp = CR(2, IP)           ! m/s
         Vs = CR(3, IP)           ! m/s

         Vp_star = CMPLX(Vp, CIM(2), kind=dp)   ! ≈ Vp (1 - i/Qp)
         Vs_star = CMPLX(Vs, CIM(3), kind=dp)   ! ≈ Vs (1 - i/Qs)

         ! Use complex density RHO coming from your original assignment
         mu_star     = RHO*Vs_star*Vs_star
         lambda_star = RHO*Vp_star*Vp_star - 2.0_dp*mu_star
         !    END IF

         ! Build isotropic stiffness tensor
         C11 = lambda_star + 2.0_dp*mu_star
         C44 = mu_star

         P(1)  = C11            ! c11
         P(2)  = lambda_star    ! c12
         P(3)  = lambda_star    ! c13
         P(4)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c14
         P(5)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c15
         P(6)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c16

         P(7)  = C11            ! c22
         P(8)  = lambda_star    ! c23
         P(9)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c24
         P(10) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c25
         P(11) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c26

         P(12) = C11            ! c33
         P(13) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c34
         P(14) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c35
         P(15) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c36

         P(16) = C44            ! c44
         P(17) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c45
         P(18) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c46
         P(19) = C44            ! c55
         P(20) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c56
         P(21) = C44            ! c66

      CASE (6)
         P(1)  = CMPLX(CR(2, IP), CIM(2),  kind=dp)                    ! c11
         P(2)  = P(1) - 2.0_dp*CMPLX(CR(6, IP), CIM(6), kind=dp)       ! c12
         P(3)  = CMPLX(CR(3, IP), CIM(3),  kind=dp)                    ! c13
         P(4)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c14
         P(5)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c15
         P(6)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c16
         P(7)  = P(1)                                                 ! c22
         P(8)  = P(3)                                                 ! c23
         P(9)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c24
         P(10) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c25
         P(11) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c26
         P(12) = CMPLX(CR(4, IP), CIM(4),  kind=dp)                    ! c33
         P(13) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c34
         P(14) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c35
         P(15) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c36
         P(16) = CMPLX(CR(5, IP), CIM(5),  kind=dp)                    ! c44
         P(17) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c45
         P(18) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c46
         P(19) = P(16)                                                ! c55
         P(20) = CMPLX(0.0_dp, 0.0_dp, kind=dp)                        ! c56
         P(21) = CMPLX(CR(6, IP), CIM(6),  kind=dp)                    ! c66

      CASE (7)

         !---- 2D TTI-medium ---------
         IF (IANISO == 7) THEN
            ! Complex VTI stiffnesses in the local symmetry-axis frame
            C11 = CMPLX(CR(2, IP), CIM(2), kind=dp)
            C13 = CMPLX(CR(3, IP), CIM(3), kind=dp)
            C33 = CMPLX(CR(4, IP), CIM(4), kind=dp)
            C44 = CMPLX(CR(5, IP), CIM(5), kind=dp)
            C66 = CMPLX(CR(6, IP), CIM(6), kind=dp)
            ZT0    = REAL(CR(7, IP))*3.1415926D0/180.D0  ! Tilt angle: CR(7,IP) is in degrees, convert to radians
            COSZT0 = CMPLX(COS(ZT0), 0.0_dp, kind=dp)
            SINZT0 = CMPLX(SIN(ZT0), 0.0_dp, kind=dp)
            COSZT2 = CMPLX(COS(2.0_dp*ZT0), 0.0_dp, kind=dp)
            SINZT2 = CMPLX(SIN(2.0_dp*ZT0), 0.0_dp, kind=dp)

            ! c11 : even in tilt sign (depends on SINZT2**2)
            P(1) = C11*COSZT0**4 + C33*SINZT0**4 + &
                   0.5_dp*(C13 + 2.0_dp*C44)*SINZT2**2
            P(2) = (C11 - 2.0_dp*C66)*COSZT0**2 + C13*SINZT0**2   ! c12 : even in tilt sign
            ! c13 : even in tilt sign
            P(3) = 0.25_dp*(C11 + C33 - 4.0_dp*C44)*SINZT2**2 + &
                   C13*(COSZT0**4 + SINZT0**4)
            P(4) = CMPLX(0.0_dp, 0.0_dp, kind=dp)   ! c14
            ! c15 : odd in tilt sign (∝ SINZT0)
            P(5) = -(C11 - C13 - 2.0_dp*C44)*COSZT0**3*SINZT0 - &
                   (C13 - C33 + 2.0_dp*C44)*SINZT0**3*COSZT0
            P(6) = CMPLX(0.0_dp, 0.0_dp, kind=dp)   ! c16
            ! c22 : remains C11
            P(7) = C11           ! c22
            ! c23 : even in tilt sign
            P(8) = (C11 - 2.0_dp*C66)*SINZT0**2 + C13*COSZT0**2
            P(9) = CMPLX(0.0_dp, 0.0_dp, kind=dp)   ! c24
            ! c25 : odd in tilt sign (∝ SINZT0*COSZT0)
            P(10) = -(C11 - 2.0_dp*C66 - C13)*COSZT0*SINZT0
            P(11) = CMPLX(0.0_dp, 0.0_dp, kind=dp)   ! c26
            ! c33 : even in tilt sign
            P(12) = C11*SINZT0**4 + C33*COSZT0**4 + &
                    0.5_dp*(C13 + 2.0_dp*C44)*SINZT2**2
            P(13) = CMPLX(0.0_dp, 0.0_dp, kind=dp)   ! c34
            P(14) = (C13 + 2.0_dp*C44 - C11)*SINZT0**3*COSZT0 + &      ! c35 : odd in tilt sign
                    (C33 - C13 - 2.0_dp*C44)*COSZT0**3*SINZT0
            P(15) = CMPLX(0.0_dp, 0.0_dp, kind=dp)   ! c36

            ! c44, c55, c66 : even in tilt sign
            P(16) = C44*COSZT0**2 + C66*SINZT0**2   ! c44
            P(17) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c45
            P(18) = 0.5_dp*(C44 - C66)*SINZT2       ! c46 : odd in tilt sign (∝ SINZT2)
            P(19) = 0.25_dp*(C11 + C33 - 2.0_dp*C13)*SINZT2**2 + &
                    C44*COSZT2**2                   ! c55
            P(20) = CMPLX(0.0_dp, 0.0_dp, kind=dp)  ! c56
            P(21) = C44*SINZT0**2 + C66*COSZT0**2   ! c66
         END IF
      CASE (8)
         C11 = CMPLX(CR(2, IP), CIM(2), kind=dp)
         C13 = CMPLX(CR(3, IP), CIM(3), kind=dp)
         C33 = CMPLX(CR(4, IP), CIM(4), kind=dp)
         C44 = CMPLX(CR(5, IP), CIM(5), kind=dp)
         C66 = CMPLX(CR(6, IP), CIM(6), kind=dp)
         ZT0 = CR(7, IP)*deg2rad
         PH0 = CR(8, IP)*deg2rad

         AN(1, 1) = CMPLX(COS(ZT0)*COS(PH0), 0.0_dp, kind=dp)
         AN(1, 2) = CMPLX(COS(ZT0)*SIN(PH0), 0.0_dp, kind=dp)
         AN(1, 3) = CMPLX(-SIN(ZT0),          0.0_dp, kind=dp)
         AN(2, 1) = CMPLX(-SIN(PH0),          0.0_dp, kind=dp)
         AN(2, 2) = CMPLX(COS(PH0),           0.0_dp, kind=dp)
         AN(2, 3) = CMPLX(0.0_dp,             0.0_dp, kind=dp)
         AN(3, 1) = CMPLX(SIN(ZT0)*COS(PH0),  0.0_dp, kind=dp)
         AN(3, 2) = CMPLX(SIN(ZT0)*SIN(PH0),  0.0_dp, kind=dp)
         AN(3, 3) = CMPLX(COS(ZT0),           0.0_dp, kind=dp)

         DO I = 1, 3
            DO J = 1, 3
               II = I*IDELT(I, J) + (9 - I - J)*(1 - IDELT(I, J))
               DO K = 1, 3
                  DO L = 1, 3
                     JJ = K*IDELT(K, L) + (9 - K - L)*(1 - IDELT(K, L))
                     Q(II, JJ) = CMPLX(0.0_dp, 0.0_dp, kind=dp)
                     DO I1 = 1, 3
                        DO J1 = 1, 3
                           KK = I1*IDELT(I1, J1) + (9 - I1 - J1)*(1 - IDELT(I1, J1))
                           DO K1 = 1, 3
                              DO L1 = 1, 3
                                 LL = K1*IDELT(K1, L1) + (9 - K1 - L1)*(1 - IDELT(K1, L1))
                                 Q(II, JJ) = Q(II, JJ) + C1(KK, LL)*AN(I1, I)*AN(J1, J)*AN(K1, K)*AN(L1, L)
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         
         
         
         
         K= 0
         DO I = 1, 6
            DO J = I, 6
               K = K + 1
               P(K) = Q(I, J)
            END DO
         END DO

      CASE (10)
         P(1)  = CMPLX(CR(2,  IP), CIM(2),  kind=dp)   !c11
         P(2)  = CMPLX(CR(3,  IP), CIM(3),  kind=dp)   !c12
         P(3)  = CMPLX(CR(4,  IP), CIM(4),  kind=dp)   !c13
         P(4)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c14
         P(5)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c15
         P(6)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c16
         P(7)  = CMPLX(CR(5,  IP), CIM(5),  kind=dp)   !c22
         P(8)  = CMPLX(CR(6,  IP), CIM(6),  kind=dp)   !c23
         P(9)  = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c24
         P(10) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c35
         P(11) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c26
         P(12) = CMPLX(CR(7,  IP), CIM(7),  kind=dp)   !c33
         P(13) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c34
         P(14) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c35
         P(15) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c36
         P(16) = CMPLX(CR(8,  IP), CIM(8),  kind=dp)   !c44
         P(17) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c45
         P(18) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c46
         P(19) = CMPLX(CR(9,  IP), CIM(9),  kind=dp)   !c55
         P(20) = CMPLX(0.0_dp, 0.0_dp, kind=dp)        !c56
         P(21) = CMPLX(CR(10, IP), CIM(10), kind=dp)   !c66
      CASE (22)
         DO I = 1, 21
            P(I) = CMPLX(CR(I + 1, IP), CIM(I + 1), kind=dp)
         END DO
      END SELECT

      RETURN
END SUBROUTINE C21

END MODULE boundaries_mod
   ! SUBROUTINE Q81_NewGSRM(FREQ, NDX, NDZ, I, J, K, L, NORD, IE0, IS0, RHO, P, Q)
   !    USE iso_fortran_env, ONLY: dp => real64
   !    IMPLICIT NONE
   !    REAL(dp), INTENT(IN) :: FREQ
   !    INTEGER, INTENT(IN) :: NDX, NDZ, I, J, K, L, NORD, IE0, IS0
   !    COMPLEX(dp), INTENT(INOUT) :: RHO
   !    COMPLEX(dp), INTENT(IN) :: P(21)
   !    COMPLEX(dp), INTENT(OUT) :: Q(81)

   !    REAL(dp) :: A0, E, DI, DJ, OMIGA, DAMP, REDU, REDUi
   !    REAL(dp) :: MPL, I2, J2, II, JJ
   !    COMPLEX(dp) :: ONE

   !    A0 = 2.D0
   !    E = 1.D-16
   !    DI = 0.D0  !d(x)=0
   !    DJ = 0.D0  !d(z)=0
   !    ONE = DCMPLX(1.D0, 0.D0)
   !    OMIGA = 2.D0*3.1415926D0*FREQ

   !    !---- for d(x) & d(z) -----------------
   !    ! Compute damping coordinates DI (x-dir) and DJ (z-dir)
   !    MPL = IE0*(NORD - 1) + 1

   !    !DI        spatial damping coordinate in the x-direction, computed from subdomain I and local Gauss point K
   !    !DJ        spatial damping coordinate in the z-direction, computed from subdomain J and local Gauss point L

   !    IF (I .LE. IE0) THEN ! for left side
   !       II = (I - 1)*(NORD - 1) + K
   !       DI = A0*(DBLE(FLOAT(II - MPL))/DBLE(FLOAT(1 - MPL)))**2
   !    END IF

   !    I2 = NDX - IE0 + 1
   !    IF (I .GE. I2) THEN  ! for right side
   !       II = (I - I2)*(NORD - 1) + K
   !       DI = A0*(DBLE(FLOAT(II - 1))/DBLE(FLOAT(MPL - 1)))**2
   !    END IF

   !    IF (J .LE. IE0) THEN ! for bottom
   !       JJ = (J - 1)*(NORD - 1) + L
   !       DJ = A0*(DBLE(FLOAT(JJ - MPL))/DBLE(FLOAT(1 - MPL)))**2
   !    END IF

   !    IF (IS0 .EQ. 2) THEN        ! Optional damping at the top if IS0 = 2
   !       J2 = NDZ - IE0 + 1
   !       IF (J .GE. J2) THEN  ! for top (surface)
   !          JJ = (J - J2)*(NORD - 1) + L
   !          DJ = A0*(DBLE(FLOAT(JJ - 1))/DBLE(FLOAT(MPL - 1)))**2
   !       END IF
   !    END IF

   !    !Apply exponential damping factors to moduli for the MRM
   !    ! DSQRT(DI**2 + DJ**2)        Euclidean magnitude of damping (distance in damping coordinate space)
   !    ! REDU        Real exponential decay factor used to scale the real parts of elastic moduli
   !    ! REDUi        Negative exponential growth factor used to scale the imaginary parts of moduli
   !    DAMP = 1000.0D0*DSQRT(DI*DI + DJ*DJ) !d(x,z)=sqrt(dx)^2+dz^2)
   !    REDU = DEXP(-DSQRT(DI*DI + DJ*DJ))  !e^(-d(x,z))
   !    REDUi = -DEXP(DSQRT(DI*DI + DJ*DJ))
   !    RHO = RHO*(ONE + DCMPLX(0.D0, -DAMP/OMIGA))

   !    !---- for COMPLEX(dp)*16 elastic moduli -------

   !    ! c1jk1:
   !    Q(1) = DCMPLX(REDU*REAL(P(1)), REDUi*AIMAG(P(1)))!DCMPLX(0.D0,E) ! Q( 1)=C1111, c11=P( 1)
   !    Q(15) = DCMPLX(REDU*REAL(P(6)), REDUi*AIMAG(P(6)))!DCMPLX(0.D0,E) ! Q(15)=C1121, c16=P( 6)
   !    Q(11) = DCMPLX(REDU*REAL(P(5)), REDUi*AIMAG(P(5)))!DCMPLX(0.D0,E) ! Q(11)=C1131, c15=P( 5)
   !    Q(16) = DCMPLX(REDU*REAL(P(6)), REDUi*AIMAG(P(6)))!DCMPLX(0.D0,E) ! Q(16)=C1211, c61=P( 6)
   !    Q(79) = DCMPLX(REDU*REAL(P(21)), REDUi*AIMAG(P(21)))!DCMPLX(0.D0,E) ! Q(79)=C1221, c66=P(21)
   !    Q(71) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(71)=C1231, c65=P(20)
   !    Q(12) = DCMPLX(REDU*REAL(P(5)), REDUi*AIMAG(P(5)))!DCMPLX(0.D0,E) ! Q(12)=C1311, c51=P( 5)
   !    Q(73) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(73)=C1321, c56=P(20)
   !    Q(67) = DCMPLX(REDU*REAL(P(19)), REDUi*AIMAG(P(19)))!DCMPLX(0.D0,E) ! Q(67)=C1331, c55=P(19)

   !    ! c1jk2:
   !    Q(14) = DCMPLX(REDU*REAL(P(6)), REDUi*AIMAG(P(6)))!DCMPLX(0.D0,E) ! Q(14)=C1112, c16=P( 6)
   !    Q(2) = DCMPLX(REDU*REAL(P(2)), REDUi*AIMAG(P(2)))!DCMPLX(0.D0,E) ! Q( 2)=C1122, c12=P( 2)
   !    Q(7) = DCMPLX(REDU*REAL(P(4)), REDUi*AIMAG(P(4)))!DCMPLX(0.D0,E) ! Q( 7)=C1132, c14=P( 4)
   !    Q(78) = DCMPLX(REDU*REAL(P(21)), REDUi*AIMAG(P(21)))!DCMPLX(0.D0,E) ! Q(78)=C1212, c66=P(21)
   !    Q(29) = DCMPLX(REDU*REAL(P(11)), REDUi*AIMAG(P(11)))!DCMPLX(0.D0,E) ! Q(29)=C1222, c62=P(11)
   !    Q(59) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(59)=C1232, c64=P(18)
   !    Q(72) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(72)=C1312, c56=P(20)
   !    Q(25) = DCMPLX(REDU*REAL(P(10)), REDUi*AIMAG(P(10)))!DCMPLX(0.D0,E) ! Q(25)=C1322, c52=P(10)
   !    Q(51) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(51)=C1332, c54=P(17)

   !    ! c1jk3:
   !    Q(10) = DCMPLX(REDU*REAL(P(5)), REDUi*AIMAG(P(5)))!DCMPLX(0.D0,E) !Q(10)=C1113, c15=P( 5)
   !    Q(6) = DCMPLX(REDU*REAL(P(4)), REDUi*AIMAG(P(4)))!DCMPLX(0.D0,E) !Q( 6)=C1123, c14=P( 4)
   !    Q(4) = DCMPLX(REDU*REAL(P(3)), REDUi*AIMAG(P(3)))!DCMPLX(0.D0,E) !Q( 4)=C1133, c13=P( 3)
   !    Q(70) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) !Q(70)=C1213, c56=P(20)
   !    Q(58) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) !Q(58)=C1223, c46=P(18)
   !    Q(42) = DCMPLX(REDU*REAL(P(15)), REDUi*AIMAG(P(15)))!DCMPLX(0.D0,E) !Q(42)=C1233, c36=P(15)
   !    Q(66) = DCMPLX(REDU*REAL(P(19)), REDUi*AIMAG(P(19)))!DCMPLX(0.D0,E) !Q(66)=C1313, c55=P(19)
   !    Q(50) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) !Q(50)=C1323, c45=P(17)
   !    Q(38) = DCMPLX(REDU*REAL(P(14)), REDUi*AIMAG(P(14)))!DCMPLX(0.D0,E) !Q(38)=C1333, c35=P(14)

   !    ! c2jk1:
   !    Q(17) = DCMPLX(REDU*REAL(P(6)), REDUi*AIMAG(P(6)))!DCMPLX(0.D0,E) ! Q(17)=C2111, c61=P( 6)
   !    Q(81) = DCMPLX(REDU*REAL(P(21)), REDUi*AIMAG(P(21)))!DCMPLX(0.D0,E) ! Q(81)=C2121, c66=P(21)
   !    Q(75) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(75)=C2131, c65=P(20)
   !    Q(3) = DCMPLX(REDU*REAL(P(2)), REDUi*AIMAG(P(2)))!DCMPLX(0.D0,E) ! Q( 3)=C2211, c21=P( 2)
   !    Q(32) = DCMPLX(REDU*REAL(P(11)), REDUi*AIMAG(P(11)))!DCMPLX(0.D0,E) ! Q(32)=C2221, c26=P(11)
   !    Q(27) = DCMPLX(REDU*REAL(P(10)), REDUi*AIMAG(P(10)))!DCMPLX(0.D0,E) ! Q(27)=C2231, c25=P(10)
   !    Q(8) = DCMPLX(REDU*REAL(P(4)), REDUi*AIMAG(P(4)))!DCMPLX(0.D0,E) ! Q( 8)=C2311, c41=P( 4)
   !    Q(63) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(63)=C2321, c46=P(18)
   !    Q(53) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(53)=C2331, c45=P(17)

   !    ! c2jk2:
   !    Q(80) = DCMPLX(REDU*REAL(P(21)), REDUi*AIMAG(P(21)))!DCMPLX(0.D0,E) ! Q(80)=C2112, c66=P(21)
   !    Q(30) = DCMPLX(REDU*REAL(P(11)), REDUi*AIMAG(P(11)))!DCMPLX(0.D0,E) ! Q(30)=C2122, c62=P(11)
   !    Q(61) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(61)=C2132, c64=P(18)
   !    Q(31) = DCMPLX(REDU*REAL(P(11)), REDUi*AIMAG(P(11)))!DCMPLX(0.D0,E) ! Q(31)=C2212, c26=P(11)
   !    Q(18) = DCMPLX(REDU*REAL(P(7)), REDUi*AIMAG(P(7)))!DCMPLX(0.D0,E) ! Q(18)=C2222, c22=P( 7)
   !    Q(22) = DCMPLX(REDU*REAL(P(9)), REDUi*AIMAG(P(9)))!DCMPLX(0.D0,E) ! Q(22)=C2232, c24=P( 9)
   !    Q(62) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(62)=C2312, c46=P(18)
   !    Q(23) = DCMPLX(REDU*REAL(P(9)), REDUi*AIMAG(P(9)))!DCMPLX(0.D0,E) ! Q(23)=C2322, c42=P( 9)
   !    Q(47) = DCMPLX(REDU*REAL(P(16)), REDUi*AIMAG(P(16)))!DCMPLX(0.D0,E) ! Q(47)=C2332, c44=P(16)

   !    !c2jk3:
   !    Q(76) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(76)=C2113, c65=P(20)
   !    Q(60) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(60)=C2123, c64=P(18)
   !    Q(43) = DCMPLX(REDU*REAL(P(15)), REDUi*AIMAG(P(15)))!DCMPLX(0.D0,E) ! Q(43)=C2133, c63=P(15)
   !    Q(26) = DCMPLX(REDU*REAL(P(10)), REDUi*AIMAG(P(10)))!DCMPLX(0.D0,E) ! Q(26)=C2213, c25=P(10)
   !    Q(21) = DCMPLX(REDU*REAL(P(9)), REDUi*AIMAG(P(9)))!DCMPLX(0.D0,E) ! Q(21)=C2223, c24=P( 9)
   !    Q(19) = DCMPLX(REDU*REAL(P(8)), REDUi*AIMAG(P(8)))!DCMPLX(0.D0,E) ! Q(19)=C2233, c23=P( 8)
   !    Q(52) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(52)=C2313, c45=P(17)
   !    Q(46) = DCMPLX(REDU*REAL(P(16)), REDUi*AIMAG(P(16)))!DCMPLX(0.D0,E) ! Q(46)=C2323, c44=P(16)
   !    Q(34) = DCMPLX(REDU*REAL(P(13)), REDUi*AIMAG(P(13)))!DCMPLX(0.D0,E) ! Q(34)=C2333, c43=P(13)

   !    ! c3jk1:
   !    Q(13) = DCMPLX(REDU*REAL(P(5)), REDUi*AIMAG(P(5)))!DCMPLX(0.D0,E) ! Q(13)=C3111, c15=P( 5)
   !    Q(77) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(77)=C3121, c56=P(20)
   !    Q(69) = DCMPLX(REDU*REAL(P(19)), REDUi*AIMAG(P(19)))!DCMPLX(0.D0,E) ! Q(69)=C3131, c55=P(19)
   !    Q(9) = DCMPLX(REDU*REAL(P(4)), REDUi*AIMAG(P(4)))!DCMPLX(0.D0,E) ! Q( 9)=C3211, c14=P( 4)
   !    Q(65) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(65)=C3221, c46=P(18)
   !    Q(57) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(57)=C3231, c45=P(17)
   !    Q(5) = DCMPLX(REDU*REAL(P(3)), REDUi*AIMAG(P(3)))!DCMPLX(0.D0,E) ! Q( 5)=C3311, c13=P( 3)
   !    Q(45) = DCMPLX(REDU*REAL(P(15)), REDUi*AIMAG(P(15)))!DCMPLX(0.D0,E) ! Q(45)=C3321, c36=P(15)
   !    Q(41) = DCMPLX(REDU*REAL(P(14)), REDUi*AIMAG(P(14)))!DCMPLX(0.D0,E) ! Q(41)=C3331, c35=P(14)

   !    ! c3jk2:
   !    Q(74) = DCMPLX(REDU*REAL(P(20)), REDUi*AIMAG(P(20)))!DCMPLX(0.D0,E) ! Q(74)=C3112, c56=P(20)
   !    Q(28) = DCMPLX(REDU*REAL(P(10)), REDUi*AIMAG(P(10)))!DCMPLX(0.D0,E) ! Q(28)=C3122, c52=P(10)
   !    Q(55) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(55)=C3132, c54=P(17)
   !    Q(64) = DCMPLX(REDU*REAL(P(18)), REDUi*AIMAG(P(18)))!DCMPLX(0.D0,E) ! Q(64)=C3212, c46=P(18)
   !    Q(24) = DCMPLX(REDU*REAL(P(9)), REDUi*AIMAG(P(9)))!DCMPLX(0.D0,E) ! Q(24)=C3222, c42=P( 9)
   !    Q(49) = DCMPLX(REDU*REAL(P(16)), REDUi*AIMAG(P(16)))!DCMPLX(0.D0,E) ! Q(49)=C3232, c44=P(16)
   !    Q(44) = DCMPLX(REDU*REAL(P(15)), REDUi*AIMAG(P(15)))!DCMPLX(0.D0,E) ! Q(44)=C3312, c36=P(15)
   !    Q(20) = DCMPLX(REDU*REAL(P(8)), REDUi*AIMAG(P(8)))!DCMPLX(0.D0,E) ! Q(20)=C3322, c32=P( 8)
   !    Q(37) = DCMPLX(REDU*REAL(P(13)), REDUi*AIMAG(P(13)))!DCMPLX(0.D0,E) ! Q(37)=C3332, c34=P(13)

   !    ! c3jk3:
   !    Q(68) = DCMPLX(REDU*REAL(P(19)), REDUi*AIMAG(P(19)))!DCMPLX(0.D0,E) ! Q(68)=C3113, c55=P(19)
   !    Q(54) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(54)=C3123, c54=P(17)
   !    Q(39) = DCMPLX(REDU*REAL(P(14)), REDUi*AIMAG(P(14)))!DCMPLX(0.D0,E) ! Q(39)=C3133, c53=P(14)
   !    Q(56) = DCMPLX(REDU*REAL(P(17)), REDUi*AIMAG(P(17)))!DCMPLX(0.D0,E) ! Q(57)=C3213, c45=P(17)
   !    Q(48) = DCMPLX(REDU*REAL(P(16)), REDUi*AIMAG(P(16)))!DCMPLX(0.D0,E) ! Q(48)=C3223, c44=P(16)
   !    Q(35) = DCMPLX(REDU*REAL(P(13)), REDUi*AIMAG(P(13)))!DCMPLX(0.D0,E) ! Q(35)=C3233, c43=P(13)
   !    Q(40) = DCMPLX(REDU*REAL(P(14)), REDUi*AIMAG(P(14)))!DCMPLX(0.D0,E) ! Q(40)=C3313, c35=P(14)
   !    Q(36) = DCMPLX(REDU*REAL(P(13)), REDUi*AIMAG(P(13)))!DCMPLX(0.D0,E) ! Q(36)=C3323, c34=P(13)
   !    Q(33) = DCMPLX(REDU*REAL(P(12)), REDUi*AIMAG(P(12)))!DCMPLX(0.D0,E) ! Q(33)=C3333, c33=P(12)

   !    RETURN
   ! END

         !---------------------------------------------------------------------
         !     SUBROUTINE Q81_GSRM                                                                    C
         !     Convert 21 elastic moduli into 81 COMPLEX(dp) elastic moduli
         !     for SRM:
         !
         !      (1) NDX,NDZ.........how many subdomains in x- and z-direction;
         !      (2) I,J...........................considered sub-domain index;
         !      (3) IE0...........number of the sub-domains in absorbing zone;
         !      (4) BI & BJ.................absorbing functions b(x) and b(z);
         !      (5) IS0=1 or 2.......No absorbing PML at free-surface or dose;
         !      (5) RHO...............................................density;
         !      (5) P(21)...............21 complex-valued elastic moduli: Cij;
         !
         !      Returns:
         !
         !      (1) RHO...........................................PML-density;
         !      (1) Q(81)....................81 COMPLEX(dp) PML-elatsic moduli;
         !
         !          Q( 1)=c1111,  Q(21)=C2223,  Q(41)=C3331, Q(61)=C2132,
         !          Q( 2)=C1122,  Q(22)=C2232,  Q(42)=C1233, Q(62)=C2312,
         !          Q( 3)=C2211,  Q(23)=C2322,  Q(43)=C2133, Q(63)=C2321,
         !          Q( 4)=C1133,  Q(24)=C3222,  Q(44)=C3312, Q(64)=C3212,       C
         !          Q( 5)=C3311,  Q(25)=C1322,  Q(45)=C3321, Q(65)=C3221,       C
         !          Q( 6)=C1123,  Q(26)=C2213,  Q(46)=C2323, Q(66)=C1313,       C
         !          Q( 7)=C1132,  Q(27)=C2231,  Q(47)=C2332, Q(67)=C1331,       C
         !          Q( 8)=C2311,  Q(28)=C3122,  Q(48)=C3223, Q(68)=C3113,       C
         !          Q( 9)=C3211,  Q(29)=C1222,  Q(49)=C3232, Q(69)=C3131,       C
         !          Q(10)=C1113,  Q(30)=C2122,  Q(50)=C1323, Q(70)=C1213,       C
         !          Q(11)=C1131,  Q(31)=C2212,  Q(51)=C1332, Q(71)=C1231,       C
         !          Q(12)=C1311,  Q(32)=C2221,  Q(52)=C2313, Q(72)=C1312,       C
         !          Q(13)=C3111,  Q(33)=C3333,  Q(53)=C2331, Q(73)=C1321,       C
         !          Q(14)=C1112,  Q(34)=C2333,  Q(54)=C3123, Q(74)=C3112,       C
         !          Q(15)=C1121,  Q(35)=C3233,  Q(55)=C3132, Q(75)=C2131,       C
         !          Q(16)=C1211,  Q(36)=C3323,  Q(56)=C3231, Q(76)=C2113,       C
         !          Q(17)=C2111,  Q(37)=C3332,  Q(57)=C3213, Q(77)=C3121,       C
         !          Q(18)=C2222,  Q(38)=C1333,  Q(58)=C1223, Q(78)=C1212,       C
         !          Q(19)=C2233,  Q(39)=C3133,  Q(59)=C1232, Q(79)=C1221,       C
         !          Q(20)=C3322,  Q(40)=C3313,  Q(60)=C2123, Q(80)=C2112,       C
         !                                                            Q(81)=C2121,
         !                                                                      C
         !----------------------------------------------------------------------C
         !     SUBROUTINE Q81_GSRM(FREQ,NDX,NDZ,I,J,K,L,NORD,IE0,IS0,RHO,P,Q)
         !       IMPLICIT real(dp) (A-H,O-Z)
         !       COMPLEX(dp) P(21)
         !       COMPLEX(dp) RHO,Q(81)

         !       A0=2.D0
         !       E=1.D-16
         !       DI=0.D0  !d(x)=0
         !       DJ=0.D0  !d(z)=0
         !       ONE=DCMPLX(1.D0,0.D0)
         !       OMIGA =2.D0*3.1415926D0*FREQ

         ! !---- for d(x) & d(z) -----------------
         !       MPL=IE0*(NORD-1)+1

         !       IF(I.LE.IE0)THEN ! for left side
         !       II=(I-1)*(NORD-1)+K
         !       DI=A0*(DBLE(FLOAT(II-MPL))/DBLE(FLOAT(1-MPL)))**2
         !       ENDIF

         !       I2=NDX-IE0+1
         !       IF(I.GE.I2)THEN  ! for right side
         !       II=(I-I2)*(NORD-1)+K
         !       DI=A0*(DBLE(FLOAT(II-1))/DBLE(FLOAT(MPL-1)))**2
         !       ENDIF

         !       IF(J.LE.IE0)THEN ! for bottom
         !       JJ=(J-1)*(NORD-1)+L
         !       DJ=A0*(DBLE(FLOAT(JJ-MPL))/DBLE(FLOAT(1-MPL)))**2
         !       ENDIF

         !       IF(IS0.EQ.2)THEN
         !       J2=NDZ-IE0+1
         !       IF(J.GE.J2)THEN  ! for top (surface)
         !       JJ=(J-J2)*(NORD-1)+L
         !       DJ=A0*(DBLE(FLOAT(JJ-1))/DBLE(FLOAT(MPL-1)))**2
         !       ENDIF
         !       ENDIF

         !       !for the MRM
         !       DAMP=1000.0D0*DSQRT(DI*DI+DJ*DJ) !d(x,z)=sqrt(dx)^2+dz^2)
         !       REDU=DEXP(-DSQRT(DI*DI+DJ*DJ))    !e^(-d(x,z))
         !       REDUi=-DEXP(DSQRT(DI*DI+DJ*DJ))!e^(+d(x,z))
         !       RHO=RHO*(ONE+DCMPLX(0.D0,-DAMP/OMIGA))

         ! !---- for COMPLEX(dp)*16 elastic moduli -------
         !       ! c1jk1:
         !       Q( 1)=DCMPLX(REDU*REAL(P( 1)),REDUi*AIMAG(P( 1)))+DCMPLX(0.D0,E) ! Q( 1)=C1111, c11=P( 1)
         !       Q(15)=DCMPLX(REDU*REAL(P( 6)),REDUi*AIMAG(P( 6)))+DCMPLX(0.D0,E) ! Q(15)=C1121, c16=P( 6)
         !       Q(11)=DCMPLX(REDU*REAL(P( 5)),REDUi*AIMAG(P( 5)))+DCMPLX(0.D0,E) ! Q(11)=C1131, c15=P( 5)
         !       Q(16)=DCMPLX(REDU*REAL(P( 6)),REDUi*AIMAG(P( 6)))+DCMPLX(0.D0,E) ! Q(16)=C1211, c61=P( 6)
         !       Q(79)=DCMPLX(REDU*REAL(P(21)),REDUi*AIMAG(P(21)))+DCMPLX(0.D0,E) ! Q(79)=C1221, c66=P(21)
         !       Q(71)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(71)=C1231, c65=P(20)
         !       Q(12)=DCMPLX(REDU*REAL(P( 5)),REDUi*AIMAG(P( 5)))+DCMPLX(0.D0,E) ! Q(12)=C1311, c51=P( 5)
         !       Q(73)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(73)=C1321, c56=P(20)
         !       Q(67)=DCMPLX(REDU*REAL(P(19)),REDUi*AIMAG(P(19)))+DCMPLX(0.D0,E) ! Q(67)=C1331, c55=P(19)

         !       ! c1jk2:
         !       Q(14)=DCMPLX(REDU*REAL(P( 6)),REDUi*AIMAG(P( 6)))+DCMPLX(0.D0,E) ! Q(14)=C1112, c16=P( 6)
         !       Q( 2)=DCMPLX(REDU*REAL(P( 2)),REDUi*AIMAG(P( 2)))+DCMPLX(0.D0,E) ! Q( 2)=C1122, c12=P( 2)
         !       Q( 7)=DCMPLX(REDU*REAL(P( 4)),REDUi*AIMAG(P( 4)))+DCMPLX(0.D0,E) ! Q( 7)=C1132, c14=P( 4)
         !       Q(78)=DCMPLX(REDU*REAL(P(21)),REDUi*AIMAG(P(21)))+DCMPLX(0.D0,E) ! Q(78)=C1212, c66=P(21)
         !       Q(29)=DCMPLX(REDU*REAL(P(11)),REDUi*AIMAG(P(11)))+DCMPLX(0.D0,E) ! Q(29)=C1222, c62=P(11)
         !       Q(59)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(59)=C1232, c64=P(18)
         !       Q(72)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(72)=C1312, c56=P(20)
         !       Q(25)=DCMPLX(REDU*REAL(P(10)),REDUi*AIMAG(P(10)))+DCMPLX(0.D0,E) ! Q(25)=C1322, c52=P(10)
         !       Q(51)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(51)=C1332, c54=P(17)

         !       ! c1jk3:
         !       Q(10)=DCMPLX(REDU*REAL(P( 5)),REDUi*AIMAG(P( 5)))+DCMPLX(0.D0,E) !Q(10)=C1113, c15=P( 5)
         !       Q( 6)=DCMPLX(REDU*REAL(P( 4)),REDUi*AIMAG(P( 4)))+DCMPLX(0.D0,E) !Q( 6)=C1123, c14=P( 4)
         !       Q( 4)=DCMPLX(REDU*REAL(P( 3)),REDUi*AIMAG(P( 3)))+DCMPLX(0.D0,E) !Q( 4)=C1133, c13=P( 3)
         !       Q(70)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) !Q(70)=C1213, c56=P(20)
         !       Q(58)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) !Q(58)=C1223, c46=P(18)
         !       Q(42)=DCMPLX(REDU*REAL(P(15)),REDUi*AIMAG(P(15)))+DCMPLX(0.D0,E) !Q(42)=C1233, c36=P(15)
         !       Q(66)=DCMPLX(REDU*REAL(P(19)),REDUi*AIMAG(P(19)))+DCMPLX(0.D0,E) !Q(66)=C1313, c55=P(19)
         !       Q(50)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) !Q(50)=C1323, c45=P(17)
         !       Q(38)=DCMPLX(REDU*REAL(P(14)),REDUi*AIMAG(P(14)))+DCMPLX(0.D0,E) !Q(38)=C1333, c35=P(14)

         !       ! c2jk1:
         !       Q(17)=DCMPLX(REDU*REAL(P( 6)),REDUi*AIMAG(P( 6)))+DCMPLX(0.D0,E) ! Q(17)=C2111, c61=P( 6)
         !       Q(81)=DCMPLX(REDU*REAL(P(21)),REDUi*AIMAG(P(21)))+DCMPLX(0.D0,E) ! Q(81)=C2121, c66=P(21)
         !       Q(75)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(75)=C2131, c65=P(20)
         !       Q( 3)=DCMPLX(REDU*REAL(P( 2)),REDUi*AIMAG(P( 2)))+DCMPLX(0.D0,E) ! Q( 3)=C2211, c21=P( 2)
         !       Q(32)=DCMPLX(REDU*REAL(P(11)),REDUi*AIMAG(P(11)))+DCMPLX(0.D0,E) ! Q(32)=C2221, c26=P(11)
         !       Q(27)=DCMPLX(REDU*REAL(P(10)),REDUi*AIMAG(P(10)))+DCMPLX(0.D0,E) ! Q(27)=C2231, c25=P(10)
         !       Q( 8)=DCMPLX(REDU*REAL(P( 4)),REDUi*AIMAG(P( 4)))+DCMPLX(0.D0,E) ! Q( 8)=C2311, c41=P( 4)
         !       Q(63)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(63)=C2321, c46=P(18)
         !       Q(53)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(53)=C2331, c45=P(17)

         !       ! c2jk2:
         !       Q(80)=DCMPLX(REDU*REAL(P(21)),REDUi*AIMAG(P(21)))+DCMPLX(0.D0,E) ! Q(80)=C2112, c66=P(21)
         !       Q(30)=DCMPLX(REDU*REAL(P(11)),REDUi*AIMAG(P(11)))+DCMPLX(0.D0,E) ! Q(30)=C2122, c62=P(11)
         !       Q(61)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(61)=C2132, c64=P(18)
         !       Q(31)=DCMPLX(REDU*REAL(P(11)),REDUi*AIMAG(P(11)))+DCMPLX(0.D0,E) ! Q(31)=C2212, c26=P(11)
         !       Q(18)=DCMPLX(REDU*REAL(P( 7)),REDUi*AIMAG(P( 7)))+DCMPLX(0.D0,E) ! Q(18)=C2222, c22=P( 7)
         !       Q(22)=DCMPLX(REDU*REAL(P( 9)),REDUi*AIMAG(P( 9)))+DCMPLX(0.D0,E) ! Q(22)=C2232, c24=P( 9)
         !       Q(62)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(62)=C2312, c46=P(18)
         !       Q(23)=DCMPLX(REDU*REAL(P( 9)),REDUi*AIMAG(P( 9)))+DCMPLX(0.D0,E) ! Q(23)=C2322, c42=P( 9)
         !       Q(47)=DCMPLX(REDU*REAL(P(16)),REDUi*AIMAG(P(16)))+DCMPLX(0.D0,E) ! Q(47)=C2332, c44=P(16)

         !       !c2jk3:
         !       Q(76)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(76)=C2113, c65=P(20)
         !       Q(60)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(60)=C2123, c64=P(18)
         !       Q(43)=DCMPLX(REDU*REAL(P(15)),REDUi*AIMAG(P(15)))+DCMPLX(0.D0,E) ! Q(43)=C2133, c63=P(15)
         !       Q(26)=DCMPLX(REDU*REAL(P(10)),REDUi*AIMAG(P(10)))+DCMPLX(0.D0,E) ! Q(26)=C2213, c25=P(10)
         !       Q(21)=DCMPLX(REDU*REAL(P( 9)),REDUi*AIMAG(P( 9)))+DCMPLX(0.D0,E) ! Q(21)=C2223, c24=P( 9)
         !       Q(19)=DCMPLX(REDU*REAL(P( 8)),REDUi*AIMAG(P( 8)))+DCMPLX(0.D0,E) ! Q(19)=C2233, c23=P( 8)
         !       Q(52)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(52)=C2313, c45=P(17)
         !       Q(46)=DCMPLX(REDU*REAL(P(16)),REDUi*AIMAG(P(16)))+DCMPLX(0.D0,E) ! Q(46)=C2323, c44=P(16)
         !       Q(34)=DCMPLX(REDU*REAL(P(13)),REDUi*AIMAG(P(13)))+DCMPLX(0.D0,E) ! Q(34)=C2333, c43=P(13)

         !       ! c3jk1:
         !       Q(13)=DCMPLX(REDU*REAL(P( 5)),REDUi*AIMAG(P( 5)))+DCMPLX(0.D0,E) ! Q(13)=C3111, c15=P( 5)
         !       Q(77)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(77)=C3121, c56=P(20)
         !       Q(69)=DCMPLX(REDU*REAL(P(19)),REDUi*AIMAG(P(19)))+DCMPLX(0.D0,E) ! Q(69)=C3131, c55=P(19)
         !       Q( 9)=DCMPLX(REDU*REAL(P( 4)),REDUi*AIMAG(P( 4)))+DCMPLX(0.D0,E) ! Q( 9)=C3211, c14=P( 4)
         !       Q(65)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(65)=C3221, c46=P(18)
         !       Q(57)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(57)=C3231, c45=P(17)
         !       Q( 5)=DCMPLX(REDU*REAL(P( 3)),REDUi*AIMAG(P( 3)))+DCMPLX(0.D0,E) ! Q( 5)=C3311, c13=P( 3)
         !       Q(45)=DCMPLX(REDU*REAL(P(15)),REDUi*AIMAG(P(15)))+DCMPLX(0.D0,E) ! Q(45)=C3321, c36=P(15)
         !       Q(41)=DCMPLX(REDU*REAL(P(14)),REDUi*AIMAG(P(14)))+DCMPLX(0.D0,E) ! Q(41)=C3331, c35=P(14)

         !       ! c3jk2:
         !       Q(74)=DCMPLX(REDU*REAL(P(20)),REDUi*AIMAG(P(20)))+DCMPLX(0.D0,E) ! Q(74)=C3112, c56=P(20)
         !       Q(28)=DCMPLX(REDU*REAL(P(10)),REDUi*AIMAG(P(10)))+DCMPLX(0.D0,E) ! Q(28)=C3122, c52=P(10)
         !       Q(55)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(55)=C3132, c54=P(17)
         !       Q(64)=DCMPLX(REDU*REAL(P(18)),REDUi*AIMAG(P(18)))+DCMPLX(0.D0,E) ! Q(64)=C3212, c46=P(18)
         !       Q(24)=DCMPLX(REDU*REAL(P( 9)),REDUi*AIMAG(P( 9)))+DCMPLX(0.D0,E) ! Q(24)=C3222, c42=P( 9)
         !       Q(49)=DCMPLX(REDU*REAL(P(16)),REDUi*AIMAG(P(16)))+DCMPLX(0.D0,E) ! Q(49)=C3232, c44=P(16)
         !       Q(44)=DCMPLX(REDU*REAL(P(15)),REDUi*AIMAG(P(15)))+DCMPLX(0.D0,E) ! Q(44)=C3312, c36=P(15)
         !       Q(20)=DCMPLX(REDU*REAL(P( 8)),REDUi*AIMAG(P( 8)))+DCMPLX(0.D0,E) ! Q(20)=C3322, c32=P( 8)
         !       Q(37)=DCMPLX(REDU*REAL(P(13)),REDUi*AIMAG(P(13)))+DCMPLX(0.D0,E) ! Q(37)=C3332, c34=P(13)

         !       ! c3jk3:
         !       Q(68)=DCMPLX(REDU*REAL(P(19)),REDUi*AIMAG(P(19)))+DCMPLX(0.D0,E) ! Q(68)=C3113, c55=P(19)
         !       Q(54)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(54)=C3123, c54=P(17)
         !       Q(39)=DCMPLX(REDU*REAL(P(14)),REDUi*AIMAG(P(14)))+DCMPLX(0.D0,E) ! Q(39)=C3133, c53=P(14)
         !       Q(56)=DCMPLX(REDU*REAL(P(17)),REDUi*AIMAG(P(17)))+DCMPLX(0.D0,E) ! Q(57)=C3213, c45=P(17)
         !       Q(48)=DCMPLX(REDU*REAL(P(16)),REDUi*AIMAG(P(16)))+DCMPLX(0.D0,E) ! Q(48)=C3223, c44=P(16)
         !       Q(35)=DCMPLX(REDU*REAL(P(13)),REDUi*AIMAG(P(13)))+DCMPLX(0.D0,E) ! Q(35)=C3233, c43=P(13)
         !       Q(40)=DCMPLX(REDU*REAL(P(14)),REDUi*AIMAG(P(14)))+DCMPLX(0.D0,E) ! Q(40)=C3313, c35=P(14)
         !       Q(36)=DCMPLX(REDU*REAL(P(13)),REDUi*AIMAG(P(13)))+DCMPLX(0.D0,E) ! Q(36)=C3323, c34=P(13)
         !       Q(33)=DCMPLX(REDU*REAL(P(12)),REDUi*AIMAG(P(12)))+DCMPLX(0.D0,E) ! Q(33)=C3333, c33=P(12)

         !       RETURN
         !       END

         !----------------------------------------------------------------------C
         !                                                                      C
         !     Convert 21 elastic moduli into 81 COMPLEX(dp) elastic moduli      C
         !     for SRM:                                                         C
         !                                                                      C
         !      (1) NDX,NDZ.........how many subdomains in x- and z-direction;  C
         !      (2) I,J...........................considered sub-domain index;  C
         !      (3) IE0...........number of the sub-domains in absorbing zone;  C
         !      (4) BI & BJ.................absorbing functions b(x) and b(z);  C
         !      (5) IS0=1 or 2.......No absorbing PML at free-surface or dose;  C
         !      (5) RHO...............................................density;  C
         !      (5) P(21)...............21 complex-valued elastic moduli: Cij;  C
         !                                                                      C
         !      Returns:                                                        C
         !                                                                      C
         !      (1) RHO...........................................PML-density;  C
         !      (1) Q(81)....................81 COMPLEX(dp) PML-elatsic moduli;  C
         !                                                                      C
         !          Q( 1)=c1111,  Q(21)=C2223,  Q(41)=C3331, Q(61)=C2132,       C
         !          Q( 2)=C1122,  Q(22)=C2232,  Q(42)=C1233, Q(62)=C2312,       C
         !          Q( 3)=C2211,  Q(23)=C2322,  Q(43)=C2133, Q(63)=C2321,       C
         !          Q( 4)=C1133,  Q(24)=C3222,  Q(44)=C3312, Q(64)=C3212,       C
         !          Q( 5)=C3311,  Q(25)=C1322,  Q(45)=C3321, Q(65)=C3221,       C
         !          Q( 6)=C1123,  Q(26)=C2213,  Q(46)=C2323, Q(66)=C1313,       C
         !          Q( 7)=C1132,  Q(27)=C2231,  Q(47)=C2332, Q(67)=C1331,       C
         !          Q( 8)=C2311,  Q(28)=C3122,  Q(48)=C3223, Q(68)=C3113,       C
         !          Q( 9)=C3211,  Q(29)=C1222,  Q(49)=C3232, Q(69)=C3131,       C
         !          Q(10)=C1113,  Q(30)=C2122,  Q(50)=C1323, Q(70)=C1213,       C
         !          Q(11)=C1131,  Q(31)=C2212,  Q(51)=C1332, Q(71)=C1231,       C
         !          Q(12)=C1311,  Q(32)=C2221,  Q(52)=C2313, Q(72)=C1312,       C
         !          Q(13)=C3111,  Q(33)=C3333,  Q(53)=C2331, Q(73)=C1321,       C
         !          Q(14)=C1112,  Q(34)=C2333,  Q(54)=C3123, Q(74)=C3112,       C
         !          Q(15)=C1121,  Q(35)=C3233,  Q(55)=C3132, Q(75)=C2131,       C
         !          Q(16)=C1211,  Q(36)=C3323,  Q(56)=C3231, Q(76)=C2113,       C
         !          Q(17)=C2111,  Q(37)=C3332,  Q(57)=C3213, Q(77)=C3121,       C
         !          Q(18)=C2222,  Q(38)=C1333,  Q(58)=C1223, Q(78)=C1212,       C
         !          Q(19)=C2233,  Q(39)=C3133,  Q(59)=C1232, Q(79)=C1221,       C
         !          Q(20)=C3322,  Q(40)=C3313,  Q(60)=C2123, Q(80)=C2112,       C
         !                                                            Q(81)=C2121,       C
         !                                                                      C
         !----------------------------------------------------------------------C
         ! SUBROUTINE Q81_SRM(FREQ,NDX,NDZ,I,J,K,L,NORD,IE0,IS0,RHO,P,Q)
         !       IMPLICIT real(dp) (A-H,O-Z)
         !       COMPLEX(dp) P(21)
         !       COMPLEX(dp) RHO,Q(81),FACT

         !       A0=2.D0
         !       E=1.D-16
         !       DI=0.D0  !d(x)=0
         !       DJ=0.D0  !d(z)=0
         !       ONE=DCMPLX(1.D0,0.D0)
         !       OMIGA =2.D0*3.1415926D0*FREQ

         ! !---- for d(x) & d(z) -----------------
         !       MPL=IE0*(NORD-1)+1

         !       IF(I.LE.IE0)THEN ! for left side
         !       II=(I-1)*(NORD-1)+K
         !       DI=A0*(DBLE(FLOAT(II-MPL))/DBLE(FLOAT(1-MPL)))**2
         !       ENDIF

         !       I2=NDX-IE0+1
         !       IF(I.GE.I2)THEN  ! for right side
         !       II=(I-I2)*(NORD-1)+K
         !       DI=A0*(DBLE(FLOAT(II-1))/DBLE(FLOAT(MPL-1)))**2
         !       ENDIF

         !       IF(J.LE.IE0)THEN ! for bottom
         !       JJ=(J-1)*(NORD-1)+L
         !       DJ=A0*(DBLE(FLOAT(JJ-MPL))/DBLE(FLOAT(1-MPL)))**2
         !       ENDIF

         !       IF(IS0.EQ.2)THEN
         !       J2=NDZ-IE0+1
         !       IF(J.GE.J2)THEN  ! for top (surface)
         !       JJ=(J-J2)*(NORD-1)+L
         !       DJ=A0*(DBLE(FLOAT(JJ-1))/DBLE(FLOAT(MPL-1)))**2
         !       ENDIF
         !       ENDIF

         !       !for the MRM
         !       DAMP=1000.0D0*DSQRT(DI*DI+DJ*DJ) !d(x,z)=sqrt(dx)^2+dz^2)
         !       REDU=DEXP(-DSQRT(DI*DI+DJ*DJ))  !e^(-d(x,z))
         !       FACT=DCMPLX(REDU,0.D0)
         !       RHO=RHO*(ONE+DCMPLX(0.D0,-DAMP/OMIGA))

         ! !---- for COMPLEX(dp)*16 elastic moduli -------
         !       ! c1jk1:
         !       Q( 1)=FACT*P( 1)+DCMPLX(0.D0,E) ! Q( 1)=C1111, c11=P( 1)
         !       Q(15)=FACT*P( 6)+DCMPLX(0.D0,E) ! Q(15)=C1121, c16=P( 6)
         !       Q(11)=FACT*P( 5)+DCMPLX(0.D0,E) ! Q(11)=C1131, c15=P( 5)
         !       Q(16)=FACT*P( 6)+DCMPLX(0.D0,E) ! Q(16)=C1211, c61=P( 6)
         !       Q(79)=FACT*P(21)+DCMPLX(0.D0,E) ! Q(79)=C1221, c66=P(21)
         !       Q(71)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(71)=C1231, c65=P(20)
         !       Q(12)=FACT*P( 5)+DCMPLX(0.D0,E) ! Q(12)=C1311, c51=P( 5)
         !       Q(73)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(73)=C1321, c56=P(20)
         !       Q(67)=FACT*P(19)+DCMPLX(0.D0,E) ! Q(67)=C1331, c55=P(19)

         !       ! c1jk2:
         !       Q(14)=FACT*P( 6)+DCMPLX(0.D0,E) ! Q(14)=C1112, c16=P( 6)
         !       Q( 2)=FACT*P( 2)+DCMPLX(0.D0,E) ! Q( 2)=C1122, c12=P( 2)
         !       Q( 7)=FACT*P( 4)+DCMPLX(0.D0,E) ! Q( 7)=C1132, c14=P( 4)
         !       Q(78)=FACT*P(21)+DCMPLX(0.D0,E) ! Q(78)=C1212, c66=P(21)
         !       Q(29)=FACT*P(11)+DCMPLX(0.D0,E) ! Q(29)=C1222, c62=P(11)
         !       Q(59)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(59)=C1232, c64=P(18)
         !       Q(72)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(72)=C1312, c56=P(20)
         !       Q(25)=FACT*P(10)+DCMPLX(0.D0,E) ! Q(25)=C1322, c52=P(10)
         !       Q(51)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(51)=C1332, c54=P(17)

         !       ! c1jk3:
         !       Q(10)=FACT*P( 5)+DCMPLX(0.D0,E) !Q(10)=C1113, c15=P( 5)
         !       Q( 6)=FACT*P( 4)+DCMPLX(0.D0,E) !Q( 6)=C1123, c14=P( 4)
         !       Q( 4)=FACT*P( 3)+DCMPLX(0.D0,E) !Q( 4)=C1133, c13=P( 3)
         !       Q(70)=FACT*P(20)+DCMPLX(0.D0,E) !Q(70)=C1213, c56=P(20)
         !       Q(58)=FACT*P(18)+DCMPLX(0.D0,E) !Q(58)=C1223, c46=P(18)
         !       Q(42)=FACT*P(15)+DCMPLX(0.D0,E) !Q(42)=C1233, c36=P(15)
         !       Q(66)=FACT*P(19)+DCMPLX(0.D0,E) !Q(66)=C1313, c55=P(19)
         !       Q(50)=FACT*P(17)+DCMPLX(0.D0,E) !Q(50)=C1323, c45=P(17)
         !       Q(38)=FACT*P(14)+DCMPLX(0.D0,E) !Q(38)=C1333, c35=P(14)

         !       ! c2jk1:
         !       Q(17)=FACT*P( 6)+DCMPLX(0.D0,E) ! Q(17)=C2111, c61=P( 6)
         !       Q(81)=FACT*P(21)+DCMPLX(0.D0,E) ! Q(81)=C2121, c66=P(21)
         !       Q(75)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(75)=C2131, c65=P(20)
         !       Q( 3)=FACT*P( 2)+DCMPLX(0.D0,E) ! Q( 3)=C2211, c21=P( 2)
         !       Q(32)=FACT*P(11)+DCMPLX(0.D0,E) ! Q(32)=C2221, c26=P(11)
         !       Q(27)=FACT*P(10)+DCMPLX(0.D0,E) ! Q(27)=C2231, c25=P(10)
         !       Q( 8)=FACT*P( 4)+DCMPLX(0.D0,E) ! Q( 8)=C2311, c41=P( 4)
         !       Q(63)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(63)=C2321, c46=P(18)
         !       Q(53)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(53)=C2331, c45=P(17)

         !       ! c2jk2:
         !       Q(80)=FACT*P(21)+DCMPLX(0.D0,E) ! Q(80)=C2112, c66=P(21)
         !       Q(30)=FACT*P(11)+DCMPLX(0.D0,E) ! Q(30)=C2122, c62=P(11)
         !       Q(61)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(61)=C2132, c64=P(18)
         !       Q(31)=FACT*P(11)+DCMPLX(0.D0,E) ! Q(31)=C2212, c26=P(11)
         !       Q(18)=FACT*P( 7)+DCMPLX(0.D0,E) ! Q(18)=C2222, c22=P( 7)
         !       Q(22)=FACT*P( 9)+DCMPLX(0.D0,E) ! Q(22)=C2232, c24=P( 9)
         !       Q(62)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(62)=C2312, c46=P(18)
         !       Q(23)=FACT*P( 9)+DCMPLX(0.D0,E) ! Q(23)=C2322, c42=P( 9)
         !       Q(47)=FACT*P(16)+DCMPLX(0.D0,E) ! Q(47)=C2332, c44=P(16)

         !       !c2jk3:
         !       Q(76)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(76)=C2113, c65=P(20)
         !       Q(60)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(60)=C2123, c64=P(18)
         !       Q(43)=FACT*P(15)+DCMPLX(0.D0,E) ! Q(43)=C2133, c63=P(15)
         !       Q(26)=FACT*P(10)+DCMPLX(0.D0,E) ! Q(26)=C2213, c25=P(10)
         !       Q(21)=FACT*P( 9)+DCMPLX(0.D0,E) ! Q(21)=C2223, c24=P( 9)
         !       Q(19)=FACT*P( 8)+DCMPLX(0.D0,E) ! Q(19)=C2233, c23=P( 8)
         !       Q(52)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(52)=C2313, c45=P(17)
         !       Q(46)=FACT*P(16)+DCMPLX(0.D0,E) ! Q(46)=C2323, c44=P(16)
         !       Q(34)=FACT*P(13)+DCMPLX(0.D0,E) ! Q(34)=C2333, c43=P(13)

         !       ! c3jk1:
         !       Q(13)=FACT*P( 5)+DCMPLX(0.D0,E) ! Q(13)=C3111, c15=P( 5)
         !       Q(77)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(77)=C3121, c56=P(20)
         !       Q(69)=FACT*P(19)+DCMPLX(0.D0,E) ! Q(69)=C3131, c55=P(19)
         !       Q( 9)=FACT*P( 4)+DCMPLX(0.D0,E) ! Q( 9)=C3211, c14=P( 4)
         !       Q(65)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(65)=C3221, c46=P(18)
         !       Q(57)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(57)=C3231, c45=P(17)
         !       Q( 5)=FACT*P( 3)+DCMPLX(0.D0,E) ! Q( 5)=C3311, c13=P( 3)
         !       Q(45)=FACT*P(15)+DCMPLX(0.D0,E) ! Q(45)=C3321, c36=P(15)
         !       Q(41)=FACT*P(14)+DCMPLX(0.D0,E) ! Q(41)=C3331, c35=P(14)

         !       ! c3jk2:
         !       Q(74)=FACT*P(20)+DCMPLX(0.D0,E) ! Q(74)=C3112, c56=P(20)
         !       Q(28)=FACT*P(10)+DCMPLX(0.D0,E) ! Q(28)=C3122, c52=P(10)
         !       Q(55)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(55)=C3132, c54=P(17)
         !       Q(64)=FACT*P(18)+DCMPLX(0.D0,E) ! Q(64)=C3212, c46=P(18)
         !       Q(24)=FACT*P( 9)+DCMPLX(0.D0,E) ! Q(24)=C3222, c42=P( 9)
         !       Q(49)=FACT*P(16)+DCMPLX(0.D0,E) ! Q(49)=C3232, c44=P(16)
         !       Q(44)=FACT*P(15)+DCMPLX(0.D0,E) ! Q(44)=C3312, c36=P(15)
         !       Q(20)=FACT*P( 8)+DCMPLX(0.D0,E) ! Q(20)=C3322, c32=P( 8)
         !       Q(37)=FACT*P(13)+DCMPLX(0.D0,E) ! Q(37)=C3332, c34=P(13)

         !       ! c3jk3:
         !       Q(68)=FACT*P(19)+DCMPLX(0.D0,E) ! Q(68)=C3113, c55=P(19)
         !       Q(54)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(54)=C3123, c54=P(17)
         !       Q(39)=FACT*P(14)+DCMPLX(0.D0,E) ! Q(39)=C3133, c53=P(14)
         !       Q(56)=FACT*P(17)+DCMPLX(0.D0,E) ! Q(57)=C3213, c45=P(17)
         !       Q(48)=FACT*P(16)+DCMPLX(0.D0,E) ! Q(48)=C3223, c44=P(16)
         !       Q(35)=FACT*P(13)+DCMPLX(0.D0,E) ! Q(35)=C3233, c43=P(13)
         !       Q(40)=FACT*P(14)+DCMPLX(0.D0,E) ! Q(40)=C3313, c35=P(14)
         !       Q(36)=FACT*P(13)+DCMPLX(0.D0,E) ! Q(36)=C3323, c34=P(13)
         !       Q(33)=FACT*P(12)+DCMPLX(0.D0,E) ! Q(33)=C3333, c33=P(12)

         !       RETURN
         !       END
         !----------------------------------------------------------------------C
         !                                                                      C
         !     Convert 21 elastic moduli into 81 COMPLEX(dp) elastic moduli      C
         !     for PML:                                                         C
         !                                                                      C
         !      (1) NDX,NDZ.........how many subdomains in x- and z-direction;  C
         !      (2) I,J...........................considered sub-domain index;  C
         !      (3) IE0...........number of the sub-domains in absorbing zone;  C
         !      (4) BI & BJ.................absorbing functions b(x) and b(z);  C
         !      (5) IS0=1 or 2.......No absorbing PML at free-surface or dose;  C
         !      (5) R0................................................density;  C
         !      (5) P(21)....................21 complex-valued elastic moduli;  C
         !                                                                      C
         !      Returns:                                                        C
         !                                                                      C
         !      (1) RHO...........................................PML-density;  C
         !      (1) Q(81)....................81 COMPLEX(dp) PML-elatsic moduli;  C
         !                                                                      C
         !          Q( 1)=c1111,  Q(21)=C2223,  Q(41)=C3331, Q(61)=C2132,       C
         !          Q( 2)=C1122,  Q(22)=C2232,  Q(42)=C1233, Q(62)=C2312,       C
         !          Q( 3)=C2211,  Q(23)=C2322,  Q(43)=C2133, Q(63)=C2321,       C
         !          Q( 4)=C1133,  Q(24)=C3222,  Q(44)=C3312, Q(64)=C3212,       C
         !          Q( 5)=C3311,  Q(25)=C1322,  Q(45)=C3321, Q(65)=C3221,       C
         !          Q( 6)=C1123,  Q(26)=C2213,  Q(46)=C2323, Q(66)=C1313,       C
         !          Q( 7)=C1132,  Q(27)=C2231,  Q(47)=C2332, Q(67)=C1331,       C
         !          Q( 8)=C2311,  Q(28)=C3122,  Q(48)=C3223, Q(68)=C3113,       C
         !          Q( 9)=C3211,  Q(29)=C1222,  Q(49)=C3232, Q(69)=C3131,       C
         !          Q(10)=C1113,  Q(30)=C2122,  Q(50)=C1323, Q(70)=C1213,       C
         !          Q(11)=C1131,  Q(31)=C2212,  Q(51)=C1332, Q(71)=C1231,       C
         !          Q(12)=C1311,  Q(32)=C2221,  Q(52)=C2313, Q(72)=C1312,       C
         !          Q(13)=C3111,  Q(33)=C3333,  Q(53)=C2331, Q(73)=C1321,       C
         !          Q(14)=C1112,  Q(34)=C2333,  Q(54)=C3123, Q(74)=C3112,       C
         !          Q(15)=C1121,  Q(35)=C3233,  Q(55)=C3132, Q(75)=C2131,       C
         !          Q(16)=C1211,  Q(36)=C3323,  Q(56)=C3231, Q(76)=C2113,       C
         !          Q(17)=C2111,  Q(37)=C3332,  Q(57)=C3213, Q(77)=C3121,       C
         !          Q(18)=C2222,  Q(38)=C1333,  Q(58)=C1223, Q(78)=C1212,       C
         !          Q(19)=C2233,  Q(39)=C3133,  Q(59)=C1232, Q(79)=C1221,       C
         !          Q(20)=C3322,  Q(40)=C3313,  Q(60)=C2123, Q(80)=C2112,       C
         !                                                            Q(81)=C2121,       C
         !                                                                      C
         !----------------------------------------------------------------------C
         !       SUBROUTINE Q81_PML(NDX,NDZ,I,J,K,L,NORD,IE0,IS0,RHO,P,Q)
         !       IMPLICIT real(dp) (A-H,O-Z)
         !       COMPLEX(dp) P(21)
         !       COMPLEX(dp) RHO,Q(81),FACT1,FACT3

         !       A0=2.D0
         !       E=1.D-16
         !       BI=0.D0
         !       BJ=0.D0

         ! !---- for b(x) & b(z) -----------------
         !       MPL=IE0*(NORD-1)+1

         !       IF(I.LE.IE0)THEN ! for left side
         !       II=(I-1)*(NORD-1)+K
         !       BI=A0*DBLE(FLOAT(II-MPL))/DBLE(FLOAT(1-MPL))
         !       ENDIF

         !       I2=NDX-IE0+1
         !       IF(I.GE.I2)THEN  ! for right side
         !       II=(I-I2)*(NORD-1)+K
         !       BI=A0*DBLE(FLOAT(II-1))/DBLE(FLOAT(MPL-1))
         !       ENDIF

         !       IF(J.LE.IE0)THEN ! for bottom
         !       JJ=(J-1)*(NORD-1)+L
         !       BJ=A0*DBLE(FLOAT(JJ-MPL))/DBLE(FLOAT(1-MPL))
         !       ENDIF

         !       IF(IS0.EQ.2)THEN
         !       J2=NDZ-IE0+1
         !       IF(J.GE.J2)THEN  ! for top (surface)
         !       JJ=(J-J2)*(NORD-1)+L
         !       BJ=A0*DBLE(FLOAT(JJ-1))/DBLE(FLOAT(MPL-1))
         !       ENDIF
         !       ENDIF

         !       FACT1=DCMPLX(1.D0,-BI)
         !       FACT3=DCMPLX(1.D0,-BJ)
         !       RHO=(FACT1*FACT3)*RHO

         ! !---- for COMPLEX(dp)*16 elastic moduli -------
         !       ! c1jk1:
         !       Q( 1)=(P( 1)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q( 1)=C1111, c11=P( 1)
         !       Q(15)=(P( 6)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(15)=C1121, c16=P( 6)
         !       Q(11)=(P( 5)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(11)=C1131, c15=P( 5)
         !       Q(16)=(P( 6)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(16)=C1211, c61=P( 6)
         !       Q(79)=(P(21)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(79)=C1221, c66=P(21)
         !       Q(71)=(P(20)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(71)=C1231, c65=P(20)
         !       Q(12)=(P( 5)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(12)=C1311, c51=P( 5)
         !       Q(73)=(P(20)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(73)=C1321, c56=P(20)
         !       Q(67)=(P(19)+DCMPLX(0.0,E))*FACT3/FACT1 ! Q(67)=C1331, c55=P(19)

         !       ! c1jk2:
         !       Q(14)=(P( 6)+DCMPLX(0.0,E))*FACT3       ! Q(14)=C1112, c16=P( 6)
         !       Q( 2)=(P( 2)+DCMPLX(0.0,E))*FACT3       ! Q( 2)=C1122, c12=P( 2)
         !       Q( 7)=(P( 4)+DCMPLX(0.0,E))*FACT3       ! Q( 7)=C1132, c14=P( 4)
         !       Q(78)=(P(21)+DCMPLX(0.0,E))*FACT3       ! Q(78)=C1212, c66=P(21)
         !       Q(29)=(P(11)+DCMPLX(0.0,E))*FACT3       ! Q(29)=C1222, c62=P(11)
         !       Q(59)=(P(18)+DCMPLX(0.0,E))*FACT3       ! Q(59)=C1232, c64=P(18)
         !       Q(72)=(P(20)+DCMPLX(0.0,E))*FACT3       ! Q(72)=C1312, c56=P(20)
         !       Q(25)=(P(10)+DCMPLX(0.0,E))*FACT3       ! Q(25)=C1322, c52=P(10)
         !       Q(51)=(P(17)+DCMPLX(0.0,E))*FACT3       ! Q(51)=C1332, c54=P(17)

         !       ! c1jk3:
         !       Q(10)=P( 5)+DCMPLX(0.0,E)               !Q(10)=C1113, c15=P( 5)
         !       Q( 6)=P( 4)+DCMPLX(0.0,E)               !Q( 6)=C1123, c14=P( 4)
         !       Q( 4)=P( 3)+DCMPLX(0.0,E)               !Q( 4)=C1133, c13=P( 3)
         !       Q(70)=P(20)+DCMPLX(0.0,E)               !Q(70)=C1213, c56=P(20)
         !       Q(58)=P(18)+DCMPLX(0.0,E)               !Q(58)=C1223, c46=P(18)
         !       Q(42)=P(15)+DCMPLX(0.0,E)               !Q(42)=C1233, c36=P(15)
         !       Q(66)=P(19)+DCMPLX(0.0,E)               !Q(66)=C1313, c55=P(19)
         !       Q(50)=P(17)+DCMPLX(0.0,E)               !Q(50)=C1323, c45=P(17)
         !       Q(38)=P(14)+DCMPLX(0.0,E)               !Q(38)=C1333, c35=P(14)

         !       ! c2jk1:
         !       Q(17)=(P( 6)+DCMPLX(0.0,E))*FACT3       ! Q(17)=C2111, c61=P( 6)
         !       Q(81)=(P(21)+DCMPLX(0.0,E))*FACT3       ! Q(81)=C2121, c66=P(21)
         !       Q(75)=(P(20)+DCMPLX(0.0,E))*FACT3       ! Q(75)=C2131, c65=P(20)
         !       Q( 3)=(P( 2)+DCMPLX(0.0,E))*FACT3       ! Q( 3)=C2211, c21=P( 2)
         !       Q(32)=(P(11)+DCMPLX(0.0,E))*FACT3       ! Q(32)=C2221, c26=P(11)
         !       Q(27)=(P(10)+DCMPLX(0.0,E))*FACT3       ! Q(27)=C2231, c25=P(10)
         !       Q( 8)=(P( 4)+DCMPLX(0.0,E))*FACT3       ! Q( 8)=C2311, c41=P( 4)
         !       Q(63)=(P(18)+DCMPLX(0.0,E))*FACT3       ! Q(63)=C2321, c46=P(18)
         !       Q(53)=(P(17)+DCMPLX(0.0,E))*FACT3       ! Q(53)=C2331, c45=P(17)

         !       ! c2jk2:
         !       Q(80)=(P(21)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(80)=C2112, c66=P(21)
         !       Q(30)=(P(11)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(30)=C2122, c62=P(11)
         !       Q(61)=(P(18)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(61)=C2132, c64=P(18)
         !       Q(31)=(P(11)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(31)=C2212, c26=P(11)
         !       Q(18)=(P( 7)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(18)=C2222, c22=P( 7)
         !       Q(22)=(P( 9)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(22)=C2232, c24=P( 9)
         !       Q(62)=(P(18)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(62)=C2312, c46=P(18)
         !       Q(23)=(P( 9)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(23)=C2322, c42=P( 9)
         !       Q(47)=(P(16)+DCMPLX(0.0,E))*FACT1*FACT3 ! Q(47)=C2332, c44=P(16)

         !       !c2jk3:
         !       Q(76)=(P(20)+DCMPLX(0.0,E))*FACT1       ! Q(76)=C2113, c65=P(20)
         !       Q(60)=(P(18)+DCMPLX(0.0,E))*FACT1       ! Q(60)=C2123, c64=P(18)
         !       Q(43)=(P(15)+DCMPLX(0.0,E))*FACT1       ! Q(43)=C2133, c63=P(15)
         !       Q(26)=(P(10)+DCMPLX(0.0,E))*FACT1       ! Q(26)=C2213, c25=P(10)
         !       Q(21)=(P( 9)+DCMPLX(0.0,E))*FACT1       ! Q(21)=C2223, c24=P( 9)
         !       Q(19)=(P( 8)+DCMPLX(0.0,E))*FACT1       ! Q(19)=C2233, c23=P( 8)
         !       Q(52)=(P(17)+DCMPLX(0.0,E))*FACT1       ! Q(52)=C2313, c45=P(17)
         !       Q(46)=(P(16)+DCMPLX(0.0,E))*FACT1       ! Q(46)=C2323, c44=P(16)
         !       Q(34)=(P(13)+DCMPLX(0.0,E))*FACT1       ! Q(34)=C2333, c43=P(13)

         !       ! c3jk1:
         !       Q(13)=P( 5)+DCMPLX(0.0,E)               ! Q(13)=C3111, c15=P( 5)
         !       Q(77)=P(20)+DCMPLX(0.0,E)               ! Q(77)=C3121, c56=P(20)
         !       Q(69)=P(19)+DCMPLX(0.0,E)               ! Q(69)=C3131, c55=P(19)
         !       Q( 9)=P( 4)+DCMPLX(0.0,E)               ! Q( 9)=C3211, c14=P( 4)
         !       Q(65)=P(18)+DCMPLX(0.0,E)               ! Q(65)=C3221, c46=P(18)
         !       Q(57)=P(17)+DCMPLX(0.0,E)               ! Q(57)=C3231, c45=P(17)
         !       Q( 5)=P( 3)+DCMPLX(0.0,E)               ! Q( 5)=C3311, c13=P( 3)
         !       Q(45)=P(15)+DCMPLX(0.0,E)               ! Q(45)=C3321, c36=P(15)
         !       Q(41)=P(14)+DCMPLX(0.0,E)               ! Q(41)=C3331, c35=P(14)

         !       ! c3jk2:
         !       Q(74)=(P(20)+DCMPLX(0.0,E))*FACT1       ! Q(74)=C3112, c56=P(20)
         !       Q(28)=(P(10)+DCMPLX(0.0,E))*FACT1       ! Q(28)=C3122, c52=P(10)
         !       Q(55)=(P(17)+DCMPLX(0.0,E))*FACT1       ! Q(55)=C3132, c54=P(17)
         !       Q(64)=(P(18)+DCMPLX(0.0,E))*FACT1       ! Q(64)=C3212, c46=P(18)
         !       Q(24)=(P( 9)+DCMPLX(0.0,E))*FACT1       ! Q(24)=C3222, c42=P( 9)
         !       Q(49)=(P(16)+DCMPLX(0.0,E))*FACT1       ! Q(49)=C3232, c44=P(16)
         !       Q(44)=(P(15)+DCMPLX(0.0,E))*FACT1       ! Q(44)=C3312, c36=P(15)
         !       Q(20)=(P( 8)+DCMPLX(0.0,E))*FACT1       ! Q(20)=C3322, c32=P( 8)
         !       Q(37)=(P(13)+DCMPLX(0.0,E))*FACT1       ! Q(37)=C3332, c34=P(13)

         !       ! c3jk3:
         !       Q(68)=(P(19)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(68)=C3113, c55=P(19)
         !       Q(54)=(P(17)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(54)=C3123, c54=P(17)
         !       Q(39)=(P(14)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(39)=C3133, c53=P(14)
         !       Q(56)=(P(17)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(57)=C3213, c45=P(17)
         !       Q(48)=(P(16)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(48)=C3223, c44=P(16)
         !       Q(35)=(P(13)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(35)=C3233, c43=P(13)
         !       Q(40)=(P(14)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(40)=C3313, c35=P(14)
         !       Q(36)=(P(13)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(36)=C3323, c34=P(13)
         !       Q(33)=(P(12)+DCMPLX(0.0,E))*FACT1/FACT3 ! Q(33)=C3333, c33=P(12)

         !       RETURN
         !       END


