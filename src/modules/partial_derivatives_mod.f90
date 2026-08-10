module partial_derivatives_mod
use iso_fortran_env, only: dp => real64
!-------------------------------
!contains the subroutines and function to compute
!    the partial FRECHET derivatives
!-----------------------------------

contains

!----------------------------------------------------------------------C
   !                                                                      C
   !     functions for the derivatives: d_Gij/d_rho, d_Gij/d_lmd          C
   !                                    d_Gij/d_mul, d_Gij/d_alp          C
   !                                    d_Gij/d_bet.                      C
   !                                    see equations (29) to (30)        C
   !----------------------------------------------------------------------C

   FUNCTION DRHO0(GS, GR)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: GS(3), GR(3)
      COMPLEX(dp)             :: DRHO0
      DRHO0 = GS(1)*GR(1) + GS(2)*GR(2) + GS(3)*GR(3)
   END FUNCTION DRHO0

   !----------------------------------------------------------------------C
   !     DRHO computes the Frechét derivative dGij/drho                  C
   !     based on stiffness coefficients and directional Green's terms. C
   !     Used in the calculation of the gradient with respect to RHO.   C
   !----------------------------------------------------------------------C
   FUNCTION DCRHO(P, DGRHO, DG11, DG13, DG33, DG44, DG66)
      IMPLICIT NONE
      ! Inputs
      COMPLEX(dp), INTENT(IN) :: P(:)                 ! [rho, alp, bet, eps, del, gam] (complex OK)
      COMPLEX(dp), INTENT(IN) :: DGRHO, DG11, DG13, DG33, DG44, DG66
      ! Return
      COMPLEX(dp)             :: DCRHO

      ! Locals
      COMPLEX(dp) :: RHO, ALP, BET, EPS, DEL, GAM
      COMPLEX(dp) :: DRHO11, DRHO13, DRHO33, DRHO44, DRHO66
      COMPLEX(dp) :: ONE

      ! Unpack parameters
      RHO = P(1)
      ALP = P(2)
      BET = P(3)
      EPS = P(4)
      DEL = P(5)
      GAM = P(6)
      ONE = (1.0D0, 0.0D0)

      ! Component sensitivities
      DRHO11 = (2.0D0*EPS + ONE)*ALP*ALP
      DRHO13 = SQRT(DEL*ALP**4 + (ALP*ALP - BET*BET)*(ALP*ALP*(EPS + ONE) - BET*BET)) - BET*BET
      DRHO33 = ALP*ALP
      DRHO44 = BET*BET
      DRHO66 = (2.0D0*GAM + ONE)*BET*BET

      ! Assemble derivative
      DCRHO = DG11*DRHO11 + DG13*DRHO13 + DG33*DRHO33 + DG44*DRHO44 + DG66*DRHO66 + DGRHO
   END FUNCTION DCRHO

   !----------------------------------------------------------
   ! DLMD: Frechét derivative w.r.t. Lamé parameter λ
   !----------------------------------------------------------
   FUNCTION DLMD(FK, GS, GR, DXGS, DZGS, DXGR, DZGR)
      IMPLICIT NONE
      REAL(dp), INTENT(IN)     :: FK
      COMPLEX(dp), INTENT(IN)  :: GS(3), DXGS(3), DZGS(3)
      COMPLEX(dp), INTENT(IN)  :: GR(3), DXGR(3), DZGR(3)
      COMPLEX(dp)              :: DGS, DGR
      COMPLEX(dp)              :: DLMD

      DGS = DXGS(1) + DZGS(3)
      DGR = DXGR(1) + DZGR(3)

      DLMD = -(DGS*DGR + DCMPLX(FK*FK, 0.D0)*GS(2)*GR(2) &
               + DCMPLX(0.D0, FK)*(GS(2)*DGR - GR(2)*DGS))  ! MY CODE

      RETURN
   END

   !----------------------------------------------------------
   ! DLMDQ: Frechét derivative w.r.t. Qλ (attenuation of λ)
   !----------------------------------------------------------
   FUNCTION DLMDQ(FK, GS, GR, DXGS, DZGS, DXGR, DZGR, IVISCO, QLMD)
      IMPLICIT NONE
      REAL(dp), INTENT(IN)     :: FK, QLMD
      INTEGER, INTENT(IN)      :: IVISCO
      COMPLEX(dp), INTENT(IN)  :: GS(3), DXGS(3), DZGS(3)
      COMPLEX(dp), INTENT(IN)  :: GR(3), DXGR(3), DZGR(3)
      COMPLEX(dp)              :: DGS, DGR
      COMPLEX(dp)              :: DLMDQ

      DGS = DXGS(1) + DZGS(3)
      DGR = DXGR(1) + DZGR(3)

      DLMDQ = -(DGS*DGR + DCMPLX(FK*FK, 0.D0)*GS(2)*GR(2) &
                + DCMPLX(0.D0, FK)*(GS(2)*DGR - GR(2)*DGS))  ! MY CODE

      IF (IVISCO .EQ. 1) DLMDQ = DLMDQ*DCMPLX(1.D0, -1/QLMD)

      RETURN
   END

   !----------------------------------------------------------
   ! DMUL: Frechét derivative w.r.t. Lamé parameter μ
   !----------------------------------------------------------
   FUNCTION DMUL(FK, GS, GR, DXGS, DZGS, DXGR, DZGR)
      IMPLICIT NONE
      REAL(dp), INTENT(IN)     :: FK
      COMPLEX(dp), INTENT(IN)  :: GS(3), DXGS(3), DZGS(3)
      COMPLEX(dp), INTENT(IN)  :: GR(3), DXGR(3), DZGR(3)
      COMPLEX(dp)              :: T1, T2, T3, T4
      COMPLEX(dp)              :: DMUL

      T1 = DXGS(1)*DXGR(1) + DZGS(1)*DXGR(3) &
           + DXGS(3)*DZGR(1) + DZGS(3)*DZGR(3)

      T2 = DXGS(1)*DXGR(1) + DZGS(1)*DZGR(1) &
           + DXGS(2)*DXGR(2) + DZGS(2)*DZGR(2) &
           + DXGS(3)*DXGR(3) + DZGS(3)*DZGR(3)

      T3 = DCMPLX(FK*FK, 0.D0)*(GS(1)*GR(1) + GS(2)*GR(2) + GS(3)*GR(3))

      T4 = DCMPLX(0.D0, FK)*(GS(1)*DXGR(2) + GS(3)*DZGR(2) &
                             - GR(1)*DXGS(2) - GR(3)*DZGS(2))

      DMUL = -(T1 + T2 + T3 + T4)  ! MY CODE

      RETURN
   END

   !----------------------------------------------------------
   ! DMULQ: Frechét derivative w.r.t. Qμ (attenuation of μ)
   !----------------------------------------------------------
   FUNCTION DMULQ(FK, GS, GR, DXGS, DZGS, DXGR, DZGR, IVISCO, QMUL)
      IMPLICIT NONE
      REAL(dp), INTENT(IN)     :: FK, QMUL
      INTEGER, INTENT(IN)      :: IVISCO
      COMPLEX(dp), INTENT(IN)  :: GS(3), DXGS(3), DZGS(3)
      COMPLEX(dp), INTENT(IN)  :: GR(3), DXGR(3), DZGR(3)
      COMPLEX(dp)              :: T1, T2, T3, T4
      COMPLEX(dp)              :: DMULQ

      T1 = DXGS(1)*DXGR(1) + DZGS(1)*DXGR(3) &
           + DXGS(3)*DZGR(1) + DZGS(3)*DZGR(3)

      T2 = DXGS(1)*DXGR(1) + DZGS(1)*DZGR(1) &
           + DXGS(2)*DXGR(2) + DZGS(2)*DZGR(2) &
           + DXGS(3)*DXGR(3) + DZGS(3)*DZGR(3)

      T3 = DCMPLX(FK*FK, 0.D0)*(GS(1)*GR(1) + GS(2)*GR(2) + GS(3)*GR(3))

      T4 = DCMPLX(0.D0, FK)*(GS(1)*DXGR(2) + GS(3)*DZGR(2) &
                             - GR(1)*DXGS(2) - GR(3)*DZGS(2))

      DMULQ = -(T1 + T2 + T3 + T4)
      IF (IVISCO .EQ. 1) DMULQ = DMULQ*DCMPLX(1.D0, -1/QMUL)

      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equations (39 & B3) for the Frechet derivative: dGij/dalp        C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DALP(P, DG11, DG13, DG33)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: P(:), DG11, DG13, DG33
      COMPLEX(dp)             :: DALP
      COMPLEX(dp)             :: RHO, ALP, BET, EPS, DEL
      COMPLEX(dp), PARAMETER  :: ONE = (1.0D0, 0.0D0), TWO = (2.0D0, 0.0D0)
      COMPLEX(dp)             :: DALP11, DALP13, DALP33

      RHO = P(1)
      ALP = P(2)
      BET = P(3)
      EPS = P(4)
      DEL = P(5)

      DALP11 = 2.0D0*RHO*ALP*(2.0D0*EPS + ONE)
      DALP13 = RHO*ALP*(2.0D0*ALP*ALP*(DEL + EPS + ONE) - BET*BET*(EPS + TWO)) &
               / SQRT(DEL*ALP**4 + (ALP*ALP - BET*BET)*(ALP*ALP*(EPS + ONE) &
                                                       - BET*BET))
      DALP33 = 2.0D0*RHO*ALP
      DALP = DG11*DALP11 + DG13*DALP13 + DG33*DALP33
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equations (39 & B3) for the Frechet derivative: dGij/dbet        C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DBET(P, DG13, DG44, DG66)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: P(:), DG13, DG44, DG66
      COMPLEX(dp)             :: DBET
      COMPLEX(dp)             :: RHO, ALP, BET, EPS, DEL, GAM
      COMPLEX(dp)             :: DBET13, DBET44, DBET66
      COMPLEX(dp), PARAMETER  :: ONE = (1.0D0, 0.0D0), TWO = (2.0D0, 0.0D0)

      ! Extract Thomsen parameters from P(:)
      RHO = P(1)
      ALP = P(2)
      BET = P(3)
      EPS = P(4)
      DEL = P(5)
      GAM = P(6)

      ! Derivatives of β with respect to cij
      DBET13 = ((TWO*BET*BET - ALP*ALP*(EPS + TWO)) / SQRT(DEL*ALP**4 + (ALP*ALP - BET*BET)*(ALP*ALP*(EPS + ONE) - BET*BET)) - TWO) * RHO * BET
      DBET44 = TWO*RHO*BET
      DBET66 = TWO*RHO*BET*(TWO*GAM + ONE)

      ! Weighted sum of kernel sensitivities
      DBET = DG13*DBET13 + DG44*DBET44 + DG66*DBET66
      RETURN
   END

   SUBROUTINE DVPS(P, DRHO_OUT, DLMD_OUT, DMUL_OUT, DALP_OUT, DBET_OUT)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN)    :: P(:)
      COMPLEX(dp), INTENT(INOUT) :: DRHO_OUT, DLMD_OUT, DMUL_OUT, DALP_OUT, DBET_OUT
      COMPLEX(dp)                :: RHO, ALP, BET

      RHO = P(1)
      ALP = P(2)
      BET = P(3)

      ! Gradient w.r.t. density
      DRHO_OUT = DRHO_OUT + (ALP**2 - 2.0D0*BET**2)*DLMD_OUT + BET**2*DMUL_OUT

      ! Gradient w.r.t. α (P-wave speed)
      DALP_OUT = 2.0D0*RHO*ALP*DLMD_OUT

      ! Gradient w.r.t. β (S-wave speed)
      DBET_OUT = 2.0D0*RHO*BET*(DMUL_OUT - 2.0D0*DLMD_OUT)
      RETURN
   END

   SUBROUTINE DVPSQ(P, DRHO_OUT, DLMD_OUT, DMUL_OUT, DALP_OUT, DBET_OUT, DQLMD, DQMUL, DQVP, DQVS, IVISCO)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN)    :: P(:)
      COMPLEX(dp), INTENT(INOUT) :: DRHO_OUT, DLMD_OUT, DMUL_OUT, DALP_OUT, DBET_OUT, DQLMD, DQMUL, DQVP, DQVS
      INTEGER, INTENT(IN)        :: IVISCO
      COMPLEX(dp)                :: RHO, ALP, BET
      REAL(dp)                   :: RHOR, ALPR, BETR, QVP, QVS
      REAL(dp)                   :: DLMDDRHO, DMULDRHO, DLMDDALP, DQLMDDALP
      REAL(dp)                   :: DLMDDBET, DMULDBET, DQLMDDBET
      REAL(dp)                   :: DLMDDQVP, DQLMDDQVP, DLMDDQVS, DQLMDDQVS, DMULDQVS, DQMULDQVS

      IF (IVISCO .EQ. 0) THEN
         RHO = REAL(P(1))
         ALP = REAL(P(2))
         BET = REAL(P(3))
         DRHO_OUT = DRHO_OUT + (ALP*ALP - 2.*BET*BET)*DLMD_OUT + BET*BET*DMUL_OUT!-------DLAM??? DLMD_OUT
         DALP_OUT = 2.*RHO*ALP*DLMD_OUT
         DBET_OUT = 2.*RHO*BET*(DMUL_OUT - 2.*DLMD_OUT)

      ELSEIF (IVISCO .EQ. 1) THEN
         RHOR = REAL(P(1))
         ALPR = REAL(P(2)); QVP = AIMAG(P(2))
         BETR = REAL(P(3)); QVS = AIMAG(P(3))
         DLMDDRHO = ALPR**2*(1 - (1/QVP)**2) - 2*BETR**2*(1 - (1/QVS)**2)
         DMULDRHO = BETR**2*(1 - (1/QVS)**2)
         DRHO_OUT = DRHO_OUT + DLMD_OUT*DLMDDRHO + DMUL_OUT*DMULDRHO
         !-------------------------------------------
         DLMDDALP = 2*ALPR*RHOR*(1 - (1/QVP)**2)
         DQLMDDALP = 2*ALPR*BETR**2*(1/QVP - 1/QVS**2/QVP - 1/QVS + 1/QVP**2/QVS) &
                     /(ALPR**2/QVP - 2*BETR**2/QVS)**2
         DALP_OUT = DLMD_OUT*DLMDDALP + DQLMD*DQLMDDALP
         !---------------------------------------------
         DLMDDBET = -4*RHOR*BETR*(1 - 1/(QVS**2))
         DMULDBET = 2*RHOR*BETR*(1 - 1/(QVS**2))
         DQLMDDBET = 2*ALPR**2*BETR*(1/QVP - 1/QVS**2/QVP - 1/QVS + 1/QVP**2/QVS) &
                     /(ALPR**2/QVP - 2*BETR**2/QVS)**2
         DBET_OUT = DLMD_OUT*DLMDDBET + DMUL_OUT*DMULDBET + DQLMD*DQLMDDBET
         !---------------------------------------------
         DLMDDQVP = 2*RHOR*ALPR**2/(QVP**3)            !dLam/dQvp
         DQLMDDQVP = 0.5*(ALPR/QVP)**2*(ALPR**2*(1 + (1/QVP)**2) - 2*BETR**2*(1 - (1/QVS)**2 + 2/QVP/QVS)) &
                     /(ALPR**2/QVP - 2*BETR**2/QVS)**2  !dQlam/dQvp
         DQVP = DLMD_OUT*DLMDDQVP + DQLMD*DQLMDDQVP
         !---------------------------------------------
         DLMDDQVS = -4*RHOR*BETR**2/(QVS**3)
         DQLMDDQVS = -BETR**2/(QVS**2)*(ALPR**2*(1 - (1/QVP)**2 + 2/QVP/QVS) - 2*BETR**2*(1 + (1/QVS)**2)) &
                     /(ALPR**2/QVP - 2*BETR**2/QVS)**2
         DMULDQVS = 2*RHOR*BETR**2/(QVS**3)
         DQMULDQVS = (1/QVS**2 + 1)/2
         DQVS = DLMD_OUT*DLMDDQVS + DQLMD*DQLMDDQVS + DMUL_OUT*DMULDQVS + DQMUL*DQMULDQVS
      END IF
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equation (34) for the Frechet derivative: dGij/dc11              C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DC11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)

      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), P(7), DC11
      COMPLEX(dp) FACT1, FACT2, FACT3, FACT4, FACT5, FACT6, FACT7, FACT8, FACT9

      FACT1 = DCMPLX(DCOS(ZTA)**4, 0.D0)
      FACT2 = DCMPLX(0.25*DSIN(2.*ZTA)**2, 0.D0)
      FACT3 = -DCMPLX(DSIN(ZTA)*DCOS(ZTA)**3, 0.D0)
      FACT4 = DCMPLX(DSIN(ZTA)**4, 0.D0)
      FACT5 = -DCMPLX(DCOS(ZTA)*DSIN(ZTA)**3, 0.D0)
      FACT6 = DCMPLX(FK*FK, 0.D0)
      FACT7 = DCMPLX(0.D0, FK)*DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT8 = DCMPLX(0.D0, FK)*DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT9 = -DCMPLX(0.D0, FK)*DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)

      DC11 = -(FACT1*DXGS(1)*DXGR(1) &
               + FACT2*(DZGS(3)*DXGR(1) + DXGS(1)*DZGR(3)) &
               + FACT3*(DZGS(1)*DXGR(1) + DXGS(3)*DXGR(1) + DXGS(1)*DZGR(1) + DXGS(1)*DXGR(3)) &
               + FACT4*DZGS(3)*DZGR(3) &
               + FACT5*(DZGS(1)*DZGR(3) + DXGS(3)*DZGR(3) + DZGS(3)*DXGR(3) + DZGS(3)*DZGR(1)) &
               + FACT2*(DXGS(3)*DXGR(3) + DZGS(1)*DXGR(3) + DXGS(3)*DZGR(1) + DZGS(1)*DZGR(1)) &
               + FACT6*GS(2)*GR(2) &
               + FACT7*(GS(2)*DXGR(1) - GR(2)*DXGS(1)) &
               + FACT8*(GS(2)*DZGR(3) - GR(2)*DZGS(3)) &
               + FACT9*(GS(2)*DXGR(3) - GR(2)*DXGS(3) + GS(2)*DZGR(1) - GR(2)*DZGS(1)))
      IF (IVISCO .EQ. 1) DC11 = DC11*DCMPLX(1.D0, -1/AIMAG(P(2)))
      RETURN
   END
   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equation (35) for the Frechet derivative: dGij/dc13              C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DC13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), P(7), DC13
      COMPLEX(dp) FACT1, FACT2, FACT3, FACT4, FACT5, FACT6, FACT7, FACT8, FACT9

         !   FACT1=DCMPLX(0.5*DSIN(ZTA)**2,0.D0)! ZHOU CODE
      FACT1 = DCMPLX(0.5*DSIN(2.*ZTA)**2, 0.D0)! MY CODE
      FACT2 = DCMPLX(DCOS(ZTA)**4 + DSIN(ZTA)**4, 0.D0)
      FACT3 = DCMPLX(0.25*DSIN(4.*ZTA), 0.D0)
      FACT4 = DCMPLX(0.5*DSIN(2.*ZTA)**2, 0.D0)
      FACT5 = -DCMPLX(0.25*DSIN(4.*ZTA), 0.D0)
         !   FACT6=-DCMPLX(0.25*DSIN(2.*ZTA)**2,0.D0)! ZHOU CODE
      FACT6 = -DCMPLX(0.5*DSIN(2.*ZTA)**2, 0.D0)! MY CODE
      FACT7 = DCMPLX(0.D0, FK)*DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT8 = DCMPLX(0.D0, FK)*DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT9 = DCMPLX(0.D0, FK)*DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)

      DC13 = -(FACT1*DXGS(1)*DXGR(1) &
               + FACT2*(DZGS(3)*DXGR(1) + DXGS(1)*DZGR(3)) &
               + FACT3*(DXGS(3)*DXGR(1) + DXGS(1)*DXGR(3) &
                        + DZGS(1)*DXGR(1) + DXGS(1)*DZGR(1)) &
               + FACT4*DZGS(3)*DZGR(3) &
               + FACT5*(DZGS(3)*DXGR(3) + DXGS(3)*DZGR(3) &
                        + DZGS(3)*DZGR(1) + DZGS(1)*DZGR(3)) &
               + FACT6*(DXGS(3)*DXGR(3) + DZGS(1)*DXGR(3) &
                        + DXGS(3)*DZGR(1) + DZGS(1)*DZGR(1)) &
               + FACT7*(GS(2)*DXGR(1) - GR(2)*DXGS(1)) &
               + FACT8*(GS(2)*DZGR(3) - GR(2)*DZGS(3)) &
               + FACT9*(GS(2)*DXGR(3) - GR(2)*DXGS(3) &
                        + GS(2)*DZGR(1) - GR(2)*DZGS(1)))
      IF (IVISCO .EQ. 1) DC13 = DC13*DCMPLX(1.D0, -1/AIMAG(P(3)))
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equation (36) for the Frechet derivative: dGij/dc33              C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DC33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) DXGS(3), DXGR(3), DZGS(3), DZGR(3), P(7)
      COMPLEX(dp) FACT1, FACT2, FACT3, FACT4, FACT5, FACT6, DC33

      FACT1 = DCMPLX(DSIN(ZTA)**4, 0.D0)
      FACT2 = DCMPLX(0.25*DSIN(2.*ZTA)**2, 0.D0)
      FACT3 = DCMPLX(DCOS(ZTA)*DSIN(ZTA)**3, 0.D0)
      FACT4 = DCMPLX(DCOS(ZTA)**4, 0.D0)
      FACT5 = DCMPLX(DSIN(ZTA)*DCOS(ZTA)**3, 0.D0)
      FACT6 = DCMPLX(0.25*DSIN(2.*ZTA)**2, 0.D0)

      DC33 = -(FACT1*DXGS(1)*DXGR(1) &
               + FACT2*(DZGS(3)*DXGR(1) + DXGS(1)*DZGR(3)) &
               + FACT3*(DXGS(3)*DXGR(1) + DXGS(1)*DXGR(3) &
                        + DZGS(1)*DXGR(1) + DXGS(1)*DZGR(1)) &
               + FACT4*DZGS(3)*DZGR(3) &
               + FACT5*(DZGS(3)*DXGR(3) + DXGS(3)*DZGR(3) &
                        + DZGS(3)*DZGR(1) + DZGS(1)*DZGR(3)) &
               + FACT6*(DXGS(3)*DXGR(3) + DZGS(1)*DXGR(3) &
                        + DXGS(3)*DZGR(1) + DZGS(1)*DZGR(1)))
      IF (IVISCO .EQ. 1) DC33 = DC33*DCMPLX(1.D0, -1/AIMAG(P(4)))
      RETURN
   END

   FUNCTION DQ11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      !COMPLEX(dp) GS(3),GR(3),DXGS(3),DXGR(3),DZGS(3),DZGR(3),DQ11,DC11,P(7)
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), DQ11, P(7)
      DQ11 = 1/DCMPLX(1.D0, -1/AIMAG(P(2)))*DCMPLX(0.D0, -REAL(P(2))/AIMAG(P(2))**2)* &
             DC11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      RETURN
   END

   ! !----------------------------------------------------------------------C
   !                                                                      C
   !     Equation (37) for the Frechet derivative: dGij/dc44              C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DC44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), P(7)
      COMPLEX(dp) FACT01, FACT02, FACT03, FACT04, FACT05, FACT06
      COMPLEX(dp) FACT07, FACT08, FACT09, FACT10, FACT11, FACT12
      COMPLEX(dp) FACT13, FACT14, FACT15, DC44

      FACT01 = DCMPLX(DSIN(2.*ZTA)**2, 0.D0)
      FACT02 = -DCMPLX(DSIN(2.*ZTA)**2, 0.D0)
      FACT03 = DCMPLX(0.5*DSIN(4.*ZTA), 0.D0)
      FACT04 = DCMPLX(DSIN(2.*ZTA)**2, 0.D0)
      FACT05 = -DCMPLX(0.5*DSIN(4.*ZTA), 0.D0)
      FACT06 = DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT07 = DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)
      FACT08 = DCMPLX(DCOS(2.*ZTA)**2, 0.D0)
      FACT09 = DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT10 = DCMPLX(FK*FK, 0.D0)*DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT11 = DCMPLX(FK*FK, 0.D0)*DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)
      FACT12 = DCMPLX(FK*FK, 0.D0)*DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT13 = DCMPLX(0.D0, FK)*DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT14 = DCMPLX(0.D0, FK)*DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)
      FACT15 = DCMPLX(0.D0, FK)*DCMPLX(DSIN(ZTA)**2, 0.D0)

      DC44 = -(FACT01*DXGS(1)*DXGR(1) &
               + FACT02*(DZGS(3)*DXGR(1) + DXGS(1)*DZGR(3)) &
               + FACT03*(DXGS(3)*DXGR(1) + DXGS(1)*DXGR(3) &
                         + DZGS(1)*DXGR(1) + DXGS(1)*DZGR(1)) &
               + FACT04*DZGS(3)*DZGR(3) &
               + FACT05*(DZGS(3)*DXGR(3) + DXGS(3)*DZGR(3) &
                         + DZGS(3)*DZGR(1) + DZGS(1)*DZGR(3)) &
               + FACT06*DZGS(2)*DZGR(2) &
               + FACT07*(DZGS(2)*DXGR(2) + DXGS(2)*DZGR(2)) &
               + FACT08*(DXGS(3)*DXGR(3) + DZGS(1)*DXGR(3) &
                         + DXGS(3)*DZGR(1) + DZGS(1)*DZGR(1)) &
               + FACT09*DXGS(2)*DXGR(2) &
               + FACT10*GS(3)*GR(3) &
               + FACT11*(GS(3)*GR(1) + GS(1)*GR(3)) &
               + FACT12*GS(1)*GR(1) &
               + FACT13*(GS(3)*DZGR(2) - GR(3)*DZGS(2)) &
               + FACT14*(GS(3)*DXGR(2) - GR(3)*DXGS(2) &
                         + GS(1)*DZGR(2) - GR(1)*DZGS(2)) &
               + FACT15*(GS(1)*DXGR(2) - GR(1)*DXGS(2)))
      IF (IVISCO .EQ. 1) DC44 = DC44*DCMPLX(1.D0, -1/AIMAG(P(5)))
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equation (38) for the Frechet derivative: dGij/dc66              C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DC66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), P(7)
      COMPLEX(dp) FACT01, FACT02, FACT03, FACT04, FACT05, FACT06
      COMPLEX(dp) FACT07, FACT08, FACT09, FACT10, FACT11, FACT12
      COMPLEX(dp) DC66

      FACT01 = DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT02 = -DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)
      FACT03 = DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT04 = DCMPLX(FK*FK, 0.D0)*DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT05 = -DCMPLX(FK*FK, 0.D0)*DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)
      FACT06 = DCMPLX(FK*FK, 0.D0)*DCMPLX(DCOS(ZTA)**2, 0.D0)
      FACT07 = -DCMPLX(0.D0, FK)*DCMPLX(2.*DCOS(ZTA)**2, 0.D0)
      FACT08 = -DCMPLX(0.D0, FK)*DCMPLX(2.*DSIN(ZTA)**2, 0.D0)
      FACT09 = DCMPLX(0.D0, FK)*DCMPLX(DSIN(2.*ZTA), 0.D0)
      FACT10 = DCMPLX(0.D0, FK)*DCMPLX(DSIN(ZTA)**2, 0.D0)
      FACT11 = -DCMPLX(0.D0, FK)*DCMPLX(0.5*DSIN(2.*ZTA), 0.D0)
      FACT12 = DCMPLX(0.D0, FK)*DCMPLX(DCOS(ZTA)**2, 0.D0)

      DC66 = -(FACT01*DZGS(2)*DZGR(2) &
               + FACT02*(DZGS(2)*DXGR(2) + DXGS(2)*DZGR(2)) &
               + FACT03*DXGS(2)*DXGR(2) &
               + FACT04*GS(3)*GR(3) &
               + FACT05*(GS(3)*GR(1) + GS(1)*GR(3)) &
               + FACT06*GS(1)*GR(1) &
               + FACT07*(GS(2)*DXGR(1) - GR(2)*DXGS(1)) &
               + FACT08*(GS(2)*DZGR(3) - GR(2)*DZGS(3)) &
               + FACT09*(GS(2)*DXGR(3) - GR(2)*DXGS(3) &
                         + GS(2)*DZGR(1) - GR(2)*DZGS(1)) &
               + FACT10*(GS(3)*DZGR(2) - GR(3)*DZGS(2)) &
               + FACT11*(GS(3)*DXGR(2) - GR(3)*DXGS(2) &
                         + GS(1)*DZGR(2) - GR(1)*DZGS(2)) &
               + FACT12*(GS(1)*DXGR(2) - GR(1)*DXGS(2)))
      IF (IVISCO .EQ. 1) DC66 = DC66*DCMPLX(1.D0, -1/AIMAG(P(6)))
      RETURN
   END
   ! !----------------------------------------------------------------------C
   !                                                                      C
   !     Equation (A1-A3) & (B-1) for the Frechet derivative: dGij/dzta0  C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DCZTA(P, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) P(7), SINZTA, COSZTA, SINZT2, COSZT2, SINZT4, COSZT4
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3)
      COMPLEX(dp) FACT1, FACT2, FACT3, FACT4, FACT5, FACT6, FACT7
      COMPLEX(dp) FACT8, FACT9, P1, P2, P3, DCZTA
      COMPLEX(dp) C11, C13, C33, C44, C66
      IF (IVISCO .EQ. 0) THEN
         C11 = DCMPLX(REAL(P(2)), 0.D0)
         C13 = DCMPLX(REAL(P(3)), 0.D0)
         C33 = DCMPLX(REAL(P(4)), 0.D0)
         C44 = DCMPLX(REAL(P(5)), 0.D0)
         C66 = DCMPLX(REAL(P(6)), 0.D0)
      ELSEIF (IVISCO .EQ. 1) THEN
         C11 = DCMPLX(REAL(P(2)), -REAL(P(2))/AIMAG(P(2)))
         C13 = DCMPLX(REAL(P(3)), -REAL(P(3))/AIMAG(P(3)))
         C33 = DCMPLX(REAL(P(4)), -REAL(P(4))/AIMAG(P(4)))
         C44 = DCMPLX(REAL(P(5)), -REAL(P(5))/AIMAG(P(5)))
         C66 = DCMPLX(REAL(P(6)), -REAL(P(6))/AIMAG(P(6)))
      END IF
      ZTA = REAL(P(7))*3.1415926D0/180.D0
      SINZTA = DCMPLX(DSIN(ZTA), 0.D0)
      COSZTA = DCMPLX(DCOS(ZTA), 0.D0)
      SINZT2 = DCMPLX(DSIN(2*ZTA), 0.D0)
      COSZT2 = DCMPLX(DCOS(2*ZTA), 0.D0)
      SINZT4 = DCMPLX(DSIN(4*ZTA), 0.D0)
      COSZT4 = DCMPLX(DCOS(4*ZTA), 0.D0)

      !---- 1st tensor-product part: eq.(A1) ---------
      !d_C11/d_zta
      FACT1 = 2.*SINZT2*(C33*SINZTA**2 - C11*COSZTA**2) &
              + (C13 + 2.*C44)*SINZT4
      !d_C13/d_zta
      FACT2 = (0.5*(C11 + C33 - 4.*C44) - C13)*SINZT4
      !d_C15/d_zta
      FACT3 = (2.*C44 + C13 - C11)*(COSZTA**4 - (3./4.) &
                                    *SINZT2**2) + (C33 - 2.*C44 - C13)*((3./4.) &
                                                                        *SINZT2**2 - SINZTA**4)
      !d_C33/d_zta
      FACT4 = 2.*SINZT2*(C11*SINZTA**2 - C33*COSZTA**2) &
              + (C13 + 2.*C44)*SINZT4
      !d_C35/d_zta
      FACT5 = (C13 + 2.*C44 - C11)*((3./4.)*SINZT2**2 &
                                    - SINZTA**4) + (C33 - C13 - 2.*C44)*(COSZTA**4 &
                                                                         - (3./4.)*SINZT2**2)
                                                                         
      !d_C44/d_zta
      FACT6 = (C66 - C44)*SINZT2
      !d_C46/d_zta
      FACT7 = (C44 - C66)*COSZT2
      !d_C55/d_zta
           FACT8=(0.5*(C11+C33-2.*C13)-2.*C44)*SINZT4!ZHOU CODE
      ! FACT8 = 0.5*(C11 + C33 - 2.*C13 - 2.*C44)*SINZT4! MY CODE
      !d_C66/d_zta
      FACT9 = (C44 - C66)*SINZT2

      P1 = FACT1*DXGS(1)*DXGR(1) &
           + FACT2*(DZGS(3)*DXGR(1) + DXGS(1)*DZGR(3)) &
           + FACT3*(DXGS(3)*DXGR(1) + DXGS(1)*DXGR(3) &
                    + DZGS(1)*DXGR(1) + DXGS(1)*DZGR(1)) &
           + FACT4*DZGS(3)*DZGR(3) &
           + FACT5*(DZGS(3)*DXGR(3) + DXGS(3)*DZGR(3) &
                    + DZGS(3)*DZGR(1) + DZGS(1)*DZGR(3)) &
           + FACT6*DZGS(2)*DZGR(2) &
           + FACT7*(DZGS(2)*DXGR(2) + DXGS(2)*DZGR(2)) &
           + FACT8*(DXGS(3)*DXGR(3) + DZGS(1)*DXGR(3) &
                    + DXGS(3)*DZGR(1) + DZGS(1)*DZGR(1)) &
           + FACT9*DXGS(2)*DXGR(2)

      !---- 2nd tensor-product part: eq.(A2) ------------
      !d_C44/d_zta
      FACT1 = DCMPLX(FK*FK, 0.D0)*(C66 - C44)*SINZT2
      !d_C46/d_zta
      FACT2 = DCMPLX(FK*FK, 0.D0)*(C44 - C66)*COSZT2
      !d_C66/d_zta
      FACT3 = DCMPLX(FK*FK, 0.D0)*(C44 - C66)*SINZT2

      P2 = FACT1*GS(3)*GR(3) + FACT2*(GS(3)*GR(1) + GS(1)*GR(3)) &
           + FACT3*GS(1)*GR(1)

      !---- 3rd tensor-product part: eq.(A3) ------------
      !d_C12/d_zta
      FACT1 = DCMPLX(0.D0, FK)*(C13 - C11 + 2.*C66)*SINZT2
      !d_C23/d_zta
      FACT2 = DCMPLX(0.D0, FK)*(C11 - 2*C66 - C13)*SINZT2
      !d_C25/d_zta
           FACT3=DCMPLX(0.D0,FK)*(2.*C66+C13-C11)*COSZT2!  ZHOU CODE
      ! FACT3 = DCMPLX(0.D0, FK)*(2.*C66 + C13 - C11)*2.*COSZT2! MY CODE
      !d_C44/d_zta
      FACT4 = DCMPLX(0.D0, FK)*(C66 - C44)*SINZT2
      !d_C46/d_zta
      FACT5 = DCMPLX(0.D0, FK)*(C44 - C66)*COSZT2
      !d_C66/d_zta
      FACT6 = DCMPLX(0.D0, FK)*(C44 - C66)*SINZT2

      P3 = FACT1*(GS(2)*DXGR(1) - GR(2)*DXGS(1)) &
           + FACT2*(GS(2)*DZGR(3) - GR(2)*DZGS(3)) &
           + FACT3*(GS(2)*DXGR(3) - GR(2)*DXGS(3) &
                    + GS(2)*DZGR(1) - GR(2)*DZGS(1)) &
           + FACT4*(GS(3)*DZGR(2) - GR(3)*DZGS(2)) &
           + FACT5*(GS(3)*DXGR(2) - GR(3)*DXGS(2) &
                    + GS(1)*DZGR(2) - GR(1)*DZGS(2)) &
           + FACT6*(GS(1)*DXGR(2) - GR(1)*DXGS(2))

      DCZTA = -(P1 + P2 + P3)
      RETURN
   END

   !----------------------------------------------------------------------C
   !     DQ13 computes the Frechét derivative dGij/dQ13                  C
   !    
   !----------------------------------------------------------------------C
   FUNCTION DQ13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), DQ13, P(7)

      DQ13 = 1/DCMPLX(1.D0, -1/AIMAG(P(3)))*DCMPLX(0.D0, -REAL(P(3))/AIMAG(P(3))**2)* &
             DC13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      RETURN
   END

   !----------------------------------------------------------------------C
   !     DQ33 computes the Frechét derivative dGij/dQ33                  C
   !     based on the chain rule and viscoelastic perturbation theory.  C
   !----------------------------------------------------------------------C
   FUNCTION DQ33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) DXGS(3), DXGR(3), DZGS(3), DZGR(3), DQ33, P(7)

      DQ33 = 1/DCMPLX(1.D0, -1/AIMAG(P(4)))*DCMPLX(0.D0, -REAL(P(4))/AIMAG(P(4))**2)* &
             DC33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      RETURN
   END

   !----------------------------------------------------------------------C
   !     DQ44 computes the Frechét derivative dGij/dQ44                  C
   !     from directional Green's function components and Q perturb.    C
   !----------------------------------------------------------------------C
   FUNCTION DQ44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT real(dp) (A - H, O - Z)
      INTEGER :: IVISCO
      COMPLEX(dp) GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), DQ44, P(7)

      DQ44 = 1/DCMPLX(1.D0, -1/AIMAG(P(5)))*DCMPLX(0.D0, -REAL(P(5))/AIMAG(P(5))**2)* &
             DC44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      RETURN
   END

   !----------------------------------------------------------------------C
   !     DQ66 computes the Frechét derivative dGij/dQ66                  C
   !     applying the complex chain rule and viscoelastic DC66 formula. C
   !----------------------------------------------------------------------C
   FUNCTION DQ66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      IMPLICIT NONE
      REAL(dp), INTENT(IN)     :: ZTA, FK
      INTEGER, INTENT(IN)      :: IVISCO
      COMPLEX(dp), INTENT(IN)  :: GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3), P(7)
      COMPLEX(dp)              :: DQ66

      DQ66 = 1/DCMPLX(1.D0, -1/AIMAG(P(6)))*DCMPLX(0.D0, -REAL(P(6))/AIMAG(P(6))**2)* &
             DC66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equations (39 & B3) for the Frechet derivative: dGij/deps        C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DEPS(P, DG11, DG13)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: P(:), DG11, DG13
      COMPLEX(dp)             :: DEPS
      COMPLEX(dp)             :: RHO, ALP, BET, EPS, DEL
      COMPLEX(dp), PARAMETER  :: ONE = (1.0D0, 0.0D0)
      COMPLEX(dp)             :: DEPS11, DEPS13
      RHO = P(1)
      ALP = P(2)
      BET = P(3)
      EPS = P(4)
      DEL = P(5)
      DEPS11 = 2.0D0*RHO*ALP*ALP
      DEPS13 = RHO*ALP*ALP*(ALP*ALP - BET*BET) &
               /(2.0D0*SQRT(DEL*ALP**4 + (ALP*ALP - BET*BET) &
                         *(ALP*ALP*(EPS + ONE) - BET*BET)))
      DEPS = DG11*DEPS11 + DG13*DEPS13
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equations (39 & B3) for the Frechet derivative: dGij/ddel        C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DDEL(P, DG13)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: P(:), DG13
      COMPLEX(dp)             :: DDEL
      COMPLEX(dp)             :: RHO, ALP, BET, EPS, DEL
      COMPLEX(dp), PARAMETER  :: ONE = (1.0D0, 0.0D0)
      COMPLEX(dp)             :: DDEL13
      RHO = P(1)
      ALP = P(2)
      BET = P(3)
      EPS = P(4)
      DEL = P(5)
      DDEL13 = RHO*ALP**4/(2.0D0*SQRT(DEL*ALP**4 &
                                   + (ALP*ALP - BET*BET)*(ALP*ALP*(EPS + ONE) - BET*BET)))
      DDEL = DG13*DDEL13
      RETURN
   END

   !----------------------------------------------------------------------C
   !                                                                      C
   !     Equations (39 & B3) for the Frechet derivative: dGij/dgam        C
   !                                                                      C
   !----------------------------------------------------------------------C
   FUNCTION DGAM(P, DG66)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(IN) :: P(:), DG66
      COMPLEX(dp)             :: DGAM
      COMPLEX(dp)             :: RHO, BET
      RHO = P(1)
      BET = P(3)
      DGAM = DG66*2.0D0*RHO*BET*BET
      RETURN
   END

   ! !----------------------------------------------------------------------C
   ! !                                                                      C
   ! !     Equations (39 & B3) for the Frechet derivative: dGij/dQalp        C
   ! !                                                                      C
   ! !----------------------------------------------------------------------C
   ! FUNCTION DQALP(DG11, DG33, DG13, QALP, QBET, QEPS, QDEL, C13, C33, C44)
   !   IMPLICIT real(dp) (A-H,O-Z)
   !   COMPLEX(dp) DG11, DG33, DG13, DQALP
   !   real(dp) QALP, QBET, QEPS, QDEL, C13, C33, C44

   !   ! Compute intermediate parameters
   !   real(dp) A, B
   !   A = C33*(C33-C44)/(2.*C13*(C13+C44))
   !   B = C44*(C13+C33)**2/(2.*C13*(C13+C44)*(C33-C44))

   !   COMPLEX(dp) Q13, DQ11_DQALP, DQ13_DQALP
   !   Q13 = 1.0 / ((1.0 + A/QDEL + B)/QALP - B/QBET)

   !   DQ11_DQALP = 1.0 / (1.0 + QEPS)
   !   DQ13_DQALP = (Q13**2) * (1.0 + A/QDEL + B) / (QALP**2)

   !   DQALP = DG33 + DG11*DQ11_DQALP + DG13*DQ13_DQALP
   !   RETURN
   !   END
   !   !----------------------------------------------------------------------C
   !   !                                                                      C
   !   !     Equations (39 & B3) for the Frechet derivative: dGij/dQbet        C
   !   !                                                                      C
   !   !----------------------------------------------------------------------C
   ! FUNCTION DQBET(DG44, DG66, DG13, QALP, QBET, QGAM, QDEL, C13, C33, C44)
   !   IMPLICIT real(dp) (A-H,O-Z)
   !   COMPLEX(dp) DG44, DG66, DG13, DQBET
   !   real(dp) QALP, QBET, QGAM, QDEL, C13, C33, C44

   !   ! Compute intermediate parameters
   !   real(dp) A, B
   !   A = C33*(C33-C44)/(2.*C13*(C13+C44))
   !   B = C44*(C13+C33)**2/(2.*C13*(C13+C44)*(C33-C44))

   !   COMPLEX(dp) Q13, DQ44_DQBET, DQ66_DQBET, DQ13_DQBET
   !   Q13 = 1.0 / ((1.0 + A/QDEL + B)/QALP - B/QBET)

   !   DQ44_DQBET = 1.0
   !   DQ66_DQBET = -1.0 / (1.0 + QGAM)**2
   !   DQ13_DQBET = (Q13**2) * B / (QBET**2)

   !   DQBET = DG44*DQ44_DQBET + DG66*DQ66_DQBET + DG13*DQ13_DQBET
   !   RETURN
   !   END

   ! !----------------------------------------------------------------------C
   ! !                                                                      C
   ! !     Equations (39 & B3) for the Frechet derivative: dGij/dQeps        C
   ! !                                                                      C
   ! !----------------------------------------------------------------------C
   ! FUNCTION DQEPS(DG11, QEPS)
   !   IMPLICIT real(dp) (A-H,O-Z)
   !   COMPLEX(dp) DG11, DQEPS
   !   real(dp) QEPS

   !   COMPLEX(dp) DQ11_DQEPS
   !   DQ11_DQEPS = -QALP/(1.0 + QEPS)**2

   !   DQEPS = DG11 * DQ11_DQEPS
   !   RETURN
   !   END

   !   !----------------------------------------------------------------------C
   !   !                                                                      C
   !   !     Equations (39 & B3) for the Frechet derivative: dGij/dQdel        C
   !   !                                                                      C
   !   !----------------------------------------------------------------------C
   ! FUNCTION DQDEL(DG13, QALP, QBET, QDEL, C13, C33, C44)
   !   IMPLICIT real(dp) (A-H,O-Z)
   !   COMPLEX(dp) DG13, DQDEL
   !   real(dp) QALP, QBET, QDEL, C13, C33, C44

   !   ! Compute intermediate parameters
   !   real(dp) A, B
   !   A = C33*(C33-C44)/(2.*C13*(C13+C44))
   !   B = C44*(C13+C33)**2/(2.*C13*(C13+C44)*(C33-C44))

   !   COMPLEX(dp) Q13, DQ13_DQDEL
   !   Q13 = 1.0 / ((1.0 + A/QDEL + B)/QALP - B/QBET)

   !   DQ13_DQDEL = (Q13**2) * A / (QDEL**2 * QALP)

   !   DQDEL = DG13 * DQ13_DQDEL
   !   RETURN
   !   END

   !   !----------------------------------------------------------------------C
   !   !                                                                      C
   !   !     Equations (39 & B3) for the Frechet derivative: dGij/dgQam        C
   !   !                                                                      C
   !   !----------------------------------------------------------------------C
   ! FUNCTION DQGAM(DG66, QBET, QGAM)
   !   IMPLICIT real(dp) (A-H,O-Z)
   !   COMPLEX(dp) DG66, DQGAM
   !   real(dp) QBET, QGAM

   !   COMPLEX(dp) DQ66_DQGAM
   !   DQ66_DQGAM = -QBET/(1.0 + QGAM)**2

   !   DQGAM = DG66 * DQ66_DQGAM
   !   RETURN
   !   END

end module partial_derivatives_mod

!original
!   ! For viscoelastic case, extract physical quantities from complex parameters
!   IF (IVISCO .EQ. 1) THEN
!     RHOr = REAL(P(1))            ! Density
!     ALP  = REAL(P(2))            ! P-wave velocity
!     QVP  = AIMAG(P(2))            ! Qp
!     BET  = REAL(P(3))            ! S-wave velocity
!     QVS  = AIMAG(P(3))            ! Qs

!     ! Compute Lamé parameters and corresponding Q values
!     CLMD = (ALP**2 * (1 - 1/QVP**2) - 2 * BET**2 * (1 - 1/QVS**2)) * RHOr
!     CMUL = BET**2 * (1 - 1/QVS**2) * RHOr

!     QLMD = (ALP**2 * (1 - 1/QVP**2) - 2 * BET**2 * (1 - 1/QVS**2)) / &
!            (2 * ALP**2 / QVP - 4 * BET**2 / QVS)
!     QMUL = (1 - (1/QVS)**2) / (2/QVS)
!   ENDIF

!   ! Compute partial derivatives
!   DGLMD = DLMDQ(FK, GS, GR, DXGS, DZGS, DXGR, DZGR, IVISCO, QLMD)
!   DGMUL = DMULQ(FK, GS, GR, DXGS, DZGS, DXGR, DZGR, IVISCO, QMUL)

!   ! Derivatives w.r.t. QLMD and QMUL (complex-valued)
!   IF (IVISCO .EQ. 1) THEN
!     DGQLMD = DGLMD / DCMPLX(1.D0, -1/QLMD) * DCMPLX(0.D0, CLMD / QLMD**2)
!     DGQMUL = DGMUL / DCMPLX(1.D0, -1/QMUL) * DCMPLX(0.D0, CMUL / QMUL**2)
!   ENDIF

!   ! Compute derivatives for ALP and BET from Lamé parameters
!   CALL DVPSQ(P, DGRHO, DGLMD, DGMUL, DGALP, DGBET, DGQLMD, DGQMUL, DQVP, DQVS, IVISCO)

!   ! Store frequency-domain derivatives
!   D1(IK) = DGRHO
!   D2(IK) = DGALP
!   D3(IK) = DGBET
!   IF (IVISCO .EQ. 1) THEN
!     D8(IK) = DQVP
!     D9(IK) = DQVS
!   ENDIF
!   GO TO 12  ! Skip anisotropic treatment
! ENDIF
