module ky_sampling_mod
  use iso_fortran_env, ONLY : dp => real64
  use constant_mod, only : pi
  implicit none
    
contains
!-----------------------------------------------------------------------
!this module contains all the subroutine related to wavenumber singularities
!and wavenumber sampling:
! 1. CSIG: find the singularities of the wavenumbers
! 2.ind the different values in array A(*) and put them into B(*)
!    (used in CSIG)
! 3. ODER: put the elements of array B(*) in normal order
!   (used in CSIG)
! 4. CKY1: irregular sampling of wavenumbers
! 5. CKY2: regular sampling of wavenumbers
! 6. GL: Gauss-Legendre points and weights
! 7. Count_NKY: compute NK
!---------------------------------------------


    

SUBROUTINE CSIG(IANISO,NPT,CR,DX,DZ,CMIN,CMAX,FREQ2,NSIG,SIG,my_rank)
  IMPLICIT NONE
  !-----------------------------------------------------------------------
  !     CSIG computes the set of unique singular values (velocities) 
  !     relevant to wavenumber integration in 2.5D modeling or inversion,
  !     based on the input elastic moduli (CR). These singularities are
  !     used to define wavenumber sampling domains provide check for 
  !     subdomain adequacy .
  !
  !     Inputs:
  !       IANISO.............. Number of model parameters (including density)
  !       NPT................. Number of points in the model (GQG grid)
  !       CR(IANISO,NPT)...... Real-valued model parameters (moduli)
  !       CMIN................ Minimum expected velocity (for lower bound)
  !       CMAX................ Maximum expected velocity (for upper bound)
  !       FREQ2............... Maximum frequency (used to compute WLMIN)
  !       DX.................. Grid spacing in horizontal direction (X)
  !       DZ.................. Grid spacing in vertical direction (Z)
  !
  !     Outputs:
  !       NSIG................ Number of distinct singularities
  !       SIG(NSIG)........... Array of singularity values 
  !
  !     Uses:
  !       - CB, ODER.
  !
  !-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: IANISO, NPT, my_rank
  REAL(dp), INTENT(IN) :: CR(:,:)
  REAL(dp), INTENT(IN) :: CMIN, CMAX, FREQ2
  INTEGER, INTENT(OUT) :: NSIG
  REAL(dp), INTENT(OUT) :: SIG(:)
  REAL(dp), INTENT(IN) :: DX, DZ

  INTEGER, ALLOCATABLE :: NO(:)
  INTEGER :: I, I0, II, IP, J, K, MP, NP
  REAL(dp), ALLOCATABLE :: A(:), B(:)
  REAL(dp) :: DE, DD, E, WLMIN, X1, X2, new_DE

  ALLOCATE (A(NPT),B(NPT),NO(6))
  E = 5.0_dp   !(5%)
  II = IANISO-1    ! exclude density

  SELECT CASE (II)
  CASE (2)
    NO(1)=2 ! Vp
    NO(2)=3 ! Vs
    MP=2
  CASE (5,6,7)
    NO(1)=2 ! C11
    NO(2)=4 ! C33
    NO(3)=5 ! C44
    NO(4)=6 ! C66
    MP=4
  CASE DEFAULT
    NO(1)=2   !C11
    NO(2)=8   !C22
    NO(3)=13  !C33
    NO(4)=17  !C44
    NO(5)=20  !C55
    NO(6)=22  !C66
    MP=6
  END SELECT
! IF (my_rank == 0) THEN
!   WRITE(*,*) 'CSIG DEBUG: IANISO=', IANISO, ' II=', II, ' MP=', MP
!   WRITE(*,*) 'CSIG DEBUG: NO = ', (NO(i), i=1,MP)
! END IF
  IP=0
  DO I=1,MP
      !pick-up the diagonal elements
    IF (MP.EQ.2) THEN
      A(1:NPT)=CR(NO(I),1:NPT)
    ELSE
      A(1:NPT)=sqrt(CR(NO(I),1:NPT)/CR(1,1:NPT))
    ENDIF

    CALL CB(NPT,A,NP,B)
! IF (my_rank == 0) THEN
!   WRITE(*,*) 'CSIG DEBUG: I=', I, ' NO(I)=', NO(I), ' NP=', NP
!   WRITE(*,*) 'CSIG DEBUG: B(1:NP) from this diag:'
!   WRITE(*,'(1P,5E16.8)') (B(j), j=1,NP)
! END IF

    DO J=1,NP
      I0=0
      DO K=1,IP
        IF (B(J).EQ.SIG(K)) I0=1
      ENDDO
      IF (I0.EQ.0) THEN
        IP=IP+1
        SIG(IP)=B(J)
      ENDIF
    ENDDO
  ENDDO
! IF (my_rank == 0) THEN
!   WRITE(*,*) 'CSIG DEBUG: after merging, IP=', IP
!   WRITE(*,*) 'CSIG DEBUG: SIG(1:IP) before sorting:'
!   WRITE(*,'(1P,5E16.8)') (SIG(k), k=1,IP)
! END IF

  NSIG=IP
  CALL ODER(NSIG,SIG)
! IF (my_rank == 0) THEN
!   WRITE(*,*) 'CSIG DEBUG: after ODER, NSIG=', NSIG
!   WRITE(*,*) 'CSIG DEBUG: SIG sorted:'
!   WRITE(*,'(1P,5E16.8)') (SIG(k), k=1,NSIG)
! END IF

  IP=0
  X1=0.0_dp
  DO I=1,NSIG
    X2=SIG(I)
    DD=100.0_dp*(X2-X1)/X2
    IF(DD.GE.E) THEN
      IP=IP+1
      B(IP)=SIG(I)
      X1=X2
    END IF
  END DO

  NSIG=IP
  SIG(1:NSIG)=B(1:NSIG)
!   IF (my_rank == 0) THEN
!   WRITE(*,*) 'CSIG DEBUG: after 5% pruning, NSIG=', NSIG
!   WRITE(*,*) 'CSIG DEBUG: SIG(1:NSIG) after pruning:'
!   WRITE(*,'(1P,5E16.8)') (SIG(k), k=1,NSIG)
! END IF
  DEALLOCATE (A,B,NO)

  SIG(1)   = MIN(CMIN,SIG(1))         !min.speed
  SIG(NSIG)= MAX(CMAX,SIG(NSIG))      !max.speed

!   IF (my_rank == 0) THEN
!   WRITE(*,*) 'CSIG DEBUG: final SIG(1:NSIG):'
!   WRITE(*,'(1P,5E16.8)') (SIG(k), k=1,NSIG)
! END IF

  IF (my_rank==0) WRITE(*,'(A,1X,G0,1X,A,1X,G0)') 'SIG(1)', SIG(1), 'FREQ2', FREQ2
  WLMIN=SIG(1)/FREQ2
  DE=MIN(DX,DZ)

IF (my_rank == 0)  then
  WRITE(*,84) WLMIN
84 FORMAT('       Min.wavel.:', F12.2, '  Min.subd.=', F12.2)
END IF
IF (DE .GT. 0.5_dp * WLMIN) THEN
  IF (my_rank == 0) then
  WRITE(*,'(A)') '      Subdomain (block) size does not match the wavelength:'
   new_DE = 5.0_dp * FLOOR((0.5_dp * WLMIN) / 5.0_dp)
   WRITE(*,'(A,1X,G0)') '      CSIG: Subdomain size too large for the minimum wavelength. Adjusting DX and DZ :', new_DE

  END IF
  STOP 
END IF
  RETURN
END
  

SUBROUTINE Count_NKY(FREQ, XMAX, XMIN, ZMAX, ZMIN, INV, I25D, NSIG, SIG, NK, FKY, WTK, NPT, my_rank)
!-----------------------------------------------------------------------
!
!     Count_NKY determines the number of ky components (NK) and sets up
!     the corresponding quadrature points and weights depending on
!     whether the modeling is 2D or 2.5D, and on the chosen sampling strategy.
!
!     Inputs:
!       FREQ.............. Current frequency
!       XMAX, XMIN........ Lateral model extent
!       ZMAX, ZMIN........ Vertical model extent
!       INV............... Inversion flag (0 = modeling, >0 = inversion)
!       I25D.............. 2D/2.5D switch (0 = 2D, 1 = 2.5D)
!       NSIG.............. Number of singularities
!       SIG(NSIG)......... List of singularities
!       NPT............... Max length for working arrays FKYI, WTKY
!
!     Outputs:
!       NK................ Number of ky samples
!       FKY(NK)........... Output ky wavenumber array (allocated inside)
!       WTK(NK)........... Output weight array (allocated inside)
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: INV, I25D, NSIG, NPT, my_rank
  REAL(dp), INTENT(IN) :: FREQ, XMAX, XMIN, ZMAX, ZMIN
  REAL(dp), INTENT(INOUT) :: SIG(:)
  INTEGER, INTENT(OUT) :: NK
  REAL(dp), ALLOCATABLE, INTENT(OUT) :: FKY(:), WTK(:)

  REAL(dp) :: XX, ZZ, RMAX
  REAL(dp), ALLOCATABLE :: FKYI(:), WTKY(:)
  INTEGER :: II, KYSAM
  
  ALLOCATE (FKYI(NPT),WTKY(NPT))

  XX = (XMAX - XMIN)
  ZZ = (ZMAX - ZMIN)
  !   XX = 0.5_dp * (XMAX - XMIN)
  ! ZZ = 0.5_dp * (ZMAX - ZMIN)
  RMAX = sqrt(XX * XX + ZZ * ZZ)
  if (my_rank==0) write(*,'(A,1X,G0,1X,G0,1X,G0,1X,G0)') 'XSRmin, XSRmax, ZSRmin, ZSRMax', XMIN, XMAX, ZMIN, ZMAX


  ! IF (INV == 0) THEN
  !   KYSAM = 1  ! Irregular sampling
  ! ELSE
    KYSAM = 2  ! Regular sampling
  ! END IF

  IF (I25D .EQ. 0) THEN
    NK = 1
    FKYI(NK) = 0.0_dp
    WTKY(NK) = 1.0_dp
  ELSE
    IF (KYSAM .EQ. 1) THEN
      CALL CKY1(FREQ, RMAX, NSIG, SIG, NK, FKYI, WTKY)
    ELSE
      CALL CKY2(FREQ, RMAX, NSIG, SIG, NK, FKYI, WTKY, my_rank)
    END IF
  END IF

  ALLOCATE(FKY(NK), WTK(NK))
  FKY(1:NK) = FKYI(1:NK)
  WTK(1:NK) = WTKY(1:NK)

  IF (my_rank == 0) THEN
    WRITE(*,'()')
    WRITE(*,'(A, I5)') '  Number of wavenumber samples (NK) for FREQ =', NK
    WRITE(70,87)
    WRITE(70,88) (SIG(II), II = 1, NSIG)
    WRITE(70,*) '                 '
    WRITE(70,*) '  NK =', NK
    WRITE(70,88) (FKY(II), II = 1, NK)
    write(70,*) '                 '
    WRITE(*,'(A,F8.2)') '  RMAX =', RMAX
  END IF

87 FORMAT('--------------- singularities & wavenumbers ------------')
88 FORMAT(5(F18.10, 1X))

DEALLOCATE (FKYI,WTKY)
RETURN

END SUBROUTINE Count_NKY
          
!-------------------------------------------------------C
!                                                       C
!     this subroutine is to find the different values   C
!     in array A(*) and put them into B(*).             C
!                                                       C
!-------------------------------------------------------C      
    SUBROUTINE CB(NPT,A,NP,B)
    USE iso_fortran_env, ONLY : dp => real64
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NPT
    REAL(dp), INTENT(IN) :: A(:)
    REAL(dp), INTENT(OUT) :: B(:)
    INTEGER, INTENT(OUT) :: NP
    INTEGER :: I, J, IEQ

    NP = 1
    B(NP) = A(1)
    DO I = 2, NPT
      IEQ = 0
      DO J = 1, NP
        IF (A(I) .EQ. B(J)) IEQ = 1
      ENDDO
      IF (IEQ .EQ. 1) CYCLE
      NP = NP + 1
      B(NP) = A(I)
    ENDDO
    RETURN
    END
    
!------------------------------------------------------C
!                                                      C
!     put the elements of array B(*) in normal order   C
!                                                      C
!------------------------------------------------------C      
    SUBROUTINE ODER(NP,B)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP
      REAL(dp), INTENT(INOUT) :: B(:)

      REAL(dp), ALLOCATABLE :: C(:)
      REAL(dp) :: B0
      INTEGER :: I, I0, IP

      ALLOCATE (C(NP))

      DO IP = 1, NP
        B0 = huge(1.0_dp)
        I0 = 1
        DO I = 1, NP
          IF (B(I) .LT. B0) THEN
            B0 = B(I)
            I0 = I
          ENDIF
        ENDDO
        C(IP) = B0
        B(I0) = huge(1.0_dp)
      ENDDO

      B(1:NP) = C(1:NP)
      DEALLOCATE (C)
      RETURN
    END


    SUBROUTINE CKY1(FREQ,RMAX,NSIG,SIG,NK,FKYI,WTKY)
!-----------------------------------------------------------------------
!     CKY1 calculates the ky sampling points and their quadrature weights
!     using an irregular Gauss–Legendre strategy based on input singularities.
!     It constructs wavenumber segments between singularities and integrates
!     over each segment using Gauss–Legendre quadrature. The singularities
!     are internally converted to cutoffs in wavenumber space, and the 
!     integration is scaled using the geometry (RMAX) and frequency (FREQ).
!   WARNING: cannot be used for inversion
!
!     Inputs:
!       FREQ............. Frequency in Hz
!       RMAX............. Maximum source-receiver offset [m]
!       NSIG............. Number of singularities
!       SIG(NSIG)........ Array of input singularities
!
!     Outputs:
!       NK............... Number of ky samples computed
!       FKYI(NK)......... Irregularly sampled ky values
!       WTKY(NK)......... Associated quadrature weights
!
!     Internal:
!       Z, W............. Gauss–Legendre points and weights
!       A(NSIG).......... Converted cutoff wavenumbers
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NSIG
      REAL(dp), INTENT(IN) :: FREQ, RMAX
      REAL(dp), INTENT(INOUT) :: SIG(NSIG)
      INTEGER, INTENT(OUT) :: NK
      REAL(dp), INTENT(OUT) :: FKYI(:), WTKY(:)

      REAL(dp), ALLOCATABLE :: A(:), W(:), Z(:)
      REAL(dp) :: F, PIr, WLMAX, X1, X2
      INTEGER :: I, II, J, JJ

      ALLOCATE (Z(1024),W(1024),A(NSIG))

      PIr = acos(-1.0_dp)
      DO I = 1, NSIG
        A(NSIG+1-I)=2.0_dp*PIr*FREQ/SIG(I)
      ENDDO

      X2=0.0_dp
      NK=0
      WLMAX=2.0_dp*PI/A(1)

      F=2.0_dp*RMAX/WLMAX
      DO I=1,NSIG
        II=INT(2.0_dp*F/(2.0_dp**REAL(NSIG-I,dp)))
        JJ=MAX(5,II)
        CALL  GL(JJ,Z,W)
        X1=X2
        X2=A(I)
        DO J=1,JJ
          FKYI(NK+J)=0.5_dp*(X2-X1)*Z(J)+0.5_dp*(X2+X1)
          WTKY(NK+J)=0.5_dp*(X2-X1)*W(J)/PI
        ENDDO
        NK=NK+JJ
      ENDDO

      SIG(1:NSIG)=A(1:NSIG)
      DEALLOCATE (Z,W,A)
      RETURN
    END

    SUBROUTINE CKY2(FREQ,RMAX,NSIG,SIG,NK,FKYI,WTKY,my_rank)
        !-----------------------------------------------------------------------
!     CKY2 calculates the ky sampling points and their quadrature weights
!     using a regular (uniform) sampling strategy based on input singularities.
!     The number of ky points (NK) is computed based on geometry and 
!     maximum wavelength, capped at 512. Each ky point is spaced 
!     uniformly from 0 to cutoff (FKYC).
!        COMPATIBLE WITH INVERSION
!
!     Inputs:
!       FREQ............. Frequency in Hz
!       RMAX............. Maximum source-receiver offset [m] (proxy  model size half-diagonal)
!       NSIG............. Number of singularities
!       SIG(NSIG)........ Array of input singularities (used to compute CMIN/CMAX)
!
!     Returns:
!       NK............... Number of ky samples 
!       FKYI(NK)......... Uniform ky sampling points
!       WTKY(NK)......... Associated quadrature weights
!
!-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NSIG, my_rank
    REAL(dp), INTENT(IN) :: FREQ, RMAX
    REAL(dp), INTENT(INOUT) :: SIG(:)
    INTEGER, INTENT(OUT) :: NK
    REAL(dp), INTENT(OUT) :: FKYI(:), WTKY(:)

    REAL(dp) :: CMIN, CMAX, FKYC, FKYCL, WLMAX
    INTEGER :: II, N0, I
    
    ! PI = 3.1415926_dp
    CMIN = SIG(1)
    CMAX = SIG(NSIG)
    FKYC  = 2.0_dp * PI * FREQ / CMIN
      if(my_rank==0)write(*,'(A)') '1 FKYC'
    FKYC = 1.2_dp * FKYC  ! add 20% buffer to the cutoff wavenumber to ensure coverage of evanescent wavefield
    if(my_rank==0)write(*,'(A)') '1.2 FKYC'
    WLMAX = CMAX / FREQ
    if(my_rank==0)write(*,'(A,1X,I0,1X,G0,1X,G0,1X,G0)') 'my_rank, FKYC,WLMAX, RMAX', my_rank, FKYC, WLMAX, RMAX
    
    !optimal scheme
    II = INT(2.0_dp * RMAX / WLMAX)
    N0 = MIN(2*NSIG,10)*MAX(5,II) !(Nse from paper)

    !normal scheme: NK=64,128,256,512
    NK = MIN(N0,512)
    NK=INT(1.2_dp*NK) 
    if(my_rank==0)write(*,'(A)') '1.1 NK'
      if (my_rank==0) write(*,'(A,1X,I0)') 'NK in cky2=', NK
    DO I=1,NK
      FKYI(I) = REAL(I,dp) * FKYC / REAL(NK,dp)
      WTKY(I) = FKYC / (PI * REAL(NK-1,dp))
    ENDDO
    
    RETURN
    END

    !--------------------------------------------------------------------C
  !                                                                    C
  !     This subroutine is to find the abscissas X(N) and the          C
  !     weighting W(N) in (-1,1) interval for the Gauss-Legendra       C
  !     points                                                         C
  !                                                                    C
  !     Entries:  N.....number of points,                              C
  !                                                                    C
  !     Retruns:  X(*), W(*)...Gauss_Legendra points and weights       C
  !                                                                    C
  !--------------------------------------------------------------------C
    SUBROUTINE GL(N,X,W)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL(dp), INTENT(OUT) :: X(:), W(:)

      REAL(dp) :: EPS, XL, XM, X1, X2, Z, Z1, P1, P2, P3, PP
      INTEGER :: I, J, M

      X1=-1.0_dp
      X2= 1.0_dp
      EPS=1.0e-10_dp
      M=(N+1)/2
      XM=0.5_dp*(X2+X1)
      XL=0.5_dp*(X2-X1)

      DO I=1,M          ! 2nd half points
        Z=cos(PI*(REAL(I,dp)-0.25_dp)/(REAL(N,dp)+0.5_dp))

        DO
          P1= 1.0_dp
          P2= 0.0_dp
          DO J=1,N
            P3=P2
            P2=P1
            P1=((2.0_dp*REAL(J,dp)-1.0_dp)*Z*P2-(REAL(J,dp)-1.0_dp)*P3)/REAL(J,dp)
          ENDDO

          PP=REAL(N,dp)*(Z*P1-P2)/(Z*Z-1.0_dp)
          Z1=Z
          Z=Z1-P1/PP
          IF(ABS(Z-Z1).LE.EPS) EXIT
        END DO

        X(N+1-I)=XM+XL*Z
        W(N+1-I)=2.0_dp*XL/((1.0_dp-Z*Z)*PP*PP)
        X(I)=-X(N+1-I)
        W(I)=W(N+1-I)
      END DO

      RETURN
    END      
           
  

end module ky_sampling_mod
