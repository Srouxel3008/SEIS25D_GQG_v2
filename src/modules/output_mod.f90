module output_mod
  use, intrinsic :: iso_fortran_env, only : dp => real64, output_unit, error_unit
  use shared_mod
  implicit none
  contains
!----------------------------------------------------------------------
! Output utilities for writing gridded and trace data to files.
!
! Key responsibilities:
! 1. File naming helpers for Green's functions, Frechet derivatives, and
!    model updates (CFNAME_* routines).
! 2. Field export helpers:
!    - Output_GF writes real/imag Green's functions on a gridded output.
!    - GRID2D_OUT and GRID2D_OUT2 interpolate scattered (X,Z) data onto a
!      regular grid and write it in a simple columnar format suitable for
!      plotting. GRID2D_OUT optionally trims PML regions.
!    - GRID2D_OUT_FIXED writes to an explicit physical extent and resolution
!      (common grid across bands), while interpolating from shifted coordinates.
!----------------------------------------------------------------------



    SUBROUTINE WriteTraceFile(TAG, IFQ, ITER, ND, NS, NR, NSV, NRV, XSR, ZSR, DATA_C, AMP, my_rank, DEBUG_OUTPUT)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: TAG
        INTEGER, INTENT(IN) :: IFQ, ITER, ND, my_rank
        INTEGER, INTENT(IN) :: NS(:), NR(:), NSV(:), NRV(:)
        REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), AMP(:)
        COMPLEX(dp), INTENT(IN) :: DATA_C(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

        CHARACTER(LEN=30) :: FNAME20
        INTEGER :: ID, unit_out, ios
        LOGICAL :: enabled

        enabled = .TRUE.
        IF (PRESENT(DEBUG_OUTPUT)) enabled = DEBUG_OUTPUT
        IF (my_rank /= 0 .OR. .NOT. enabled) RETURN

        WRITE (FNAME20, '(A,"_",I0,"_",I0,".txt")') TRIM(TAG), IFQ, ITER
        OPEN (NEWUNIT=unit_out, FILE=TRIM(FNAME20), STATUS='UNKNOWN', POSITION='APPEND', ACTION='WRITE', IOSTAT=ios)
        IF (ios /= 0) RETURN

        DO ID = 1, ND
           WRITE (unit_out, '(7(I5,1X), 4(F10.3,1X), 2(ES22.14,1X), ES22.14)') &
              IFQ, ITER, ID, NS(ID), NSV(ID), NR(ID), NRV(ID), &
              XSR(NS(ID)), ZSR(NS(ID)), XSR(NR(ID)), ZSR(NR(ID)), &
              REAL(DATA_C(ID)), AIMAG(DATA_C(ID)), AMP(ID)
        END DO
        CLOSE (unit_out)
    END SUBROUTINE WriteTraceFile

SUBROUTINE WriteObservedDataFile(IFQ, FREQ, ND, DATA_C, my_rank, INV, DEBUG_OUTPUT)
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: IFQ, ND, my_rank, INV
   REAL(dp), INTENT(IN) :: FREQ
   COMPLEX(dp), INTENT(IN) :: DATA_C(:)
   LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

   CHARACTER(LEN=256) :: FNAME
   CHARACTER(LEN=16)  :: freq_str
   INTEGER :: I, unit_out, ios
   LOGICAL :: enabled

   enabled = .TRUE.
   IF (PRESENT(DEBUG_OUTPUT)) enabled = DEBUG_OUTPUT
   IF (my_rank /= 0 .OR. .NOT. enabled) RETURN
   IF (INV /= 0) RETURN

   WRITE(freq_str,'(F6.2)') FREQ
   freq_str = ADJUSTL(freq_str)

   IF (FREQ >= 0.0_dp .AND. FREQ < 10.0_dp .AND. freq_str(1:1) /= '-') THEN
      freq_str = '0' // TRIM(freq_str)
   END IF

   FNAME = 'OBS_FREQ_' // TRIM(freq_str) // '.txt'

   OPEN(NEWUNIT=unit_out, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
   IF (ios /= 0) RETURN

   DO I = 1, ND
      WRITE(unit_out,'(F10.4,1X,I8,1X,2(ES22.14,1X))') &
         FREQ, I, REAL(DATA_C(I)), AIMAG(DATA_C(I))
   END DO

   CLOSE(unit_out)
END SUBROUTINE WriteObservedDataFile

SUBROUTINE Output_GF(FREQ,PREFIX, IC, COMPO, NNX, NNZ, XP, ZP, FF1, FF2, NTO, XTO, ZTO)
  IMPLICIT NONE
  ! arguments
  CHARACTER(LEN=*), INTENT(IN) :: PREFIX
  INTEGER, INTENT(IN) :: IC, COMPO, NNX, NNZ, NTO
  REAL(dp), INTENT(IN) :: XP(:), ZP(:), FF1(:), FF2(:)
  REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), FREQ

  ! locals
  CHARACTER(LEN=32) :: freq_str
  CHARACTER(LEN=256) :: tmp_name
  CHARACTER(LEN=20) :: FNAME
  INTEGER :: i

  ! add leading zero for [0,10) for nicer filenames.
  WRITE(freq_str,'(F6.2)') FREQ
  freq_str = ADJUSTL(freq_str)
  ! pad leading zero for positive frequencies < 10.0 and not negative
  IF (FREQ >= 0.0_dp .AND. FREQ < 10.0_dp .AND. freq_str(1:1) /= '-') THEN
    tmp_name = '0' // TRIM(freq_str)
  ELSE
    tmp_name = TRIM(freq_str)
  END IF
  ! remove blanks and keep a short string suitable for filename
  freq_str = ADJUSTL(ADJUSTL(tmp_name))

  ! real/imaginary parts
  CALL CFNAME_GREEN(TRIM(PREFIX) // '_ReG{', IC, COMPO, '}_'// TRIM(freq_str)//'.dat', FNAME)
  CALL GRID2D_OUT0(FNAME, NNX, NNZ, XP, ZP, FF1, NTO, XTO, ZTO)

  CALL CFNAME_GREEN(TRIM(PREFIX) // '_ImG{', IC, COMPO, '}_'// TRIM(freq_str)//'.dat', FNAME)
  CALL GRID2D_OUT0(FNAME, NNX, NNZ, XP, ZP, FF2, NTO, XTO, ZTO)

END SUBROUTINE Output_GF


SUBROUTINE CFNAME_GRADIENT(PREFIX, PARAM, FREQ, ITER, EXT, FNAME)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: PREFIX, PARAM, EXT
  REAL(dp),         INTENT(IN) :: FREQ
  INTEGER,          INTENT(IN) :: ITER
  CHARACTER(LEN=*), INTENT(OUT):: FNAME

  CHARACTER(LEN=16)  :: freq_str
  CHARACTER(LEN=16)  :: iter_str
  CHARACTER(LEN=256) :: tmp_name

  WRITE(freq_str,'(F6.2)') FREQ
  freq_str = ADJUSTL(freq_str)
  IF (FREQ >= 0.0_dp .AND. FREQ < 10.0_dp .AND. freq_str(1:1) /= '-') THEN
    freq_str = '0' // TRIM(freq_str)
  END IF

  WRITE(iter_str,'(I0)') ITER
  IF (ITER >= 0 .AND. ITER < 10) THEN
    iter_str = '0' // TRIM(iter_str)
  END IF

  tmp_name = TRIM(PREFIX) // TRIM(PARAM) // '_' // TRIM(freq_str) // '_IT' // &
             TRIM(iter_str) // TRIM(EXT)

  IF (LEN(FNAME) >= LEN_TRIM(tmp_name)) THEN
    FNAME = TRIM(tmp_name)
  ELSE
    FNAME = tmp_name(1:LEN(FNAME))
  END IF
END SUBROUTINE CFNAME_GRADIENT


    !     !***************************************************************
    !CFNAME1 creates filenames for the Green Function (wavefield export) export
    SUBROUTINE CFNAME_GREEN(PREFIX, NO1, NO2, EXT, FNAME)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN)  :: PREFIX, EXT
        INTEGER,      INTENT(IN)  :: NO1, NO2
        CHARACTER(*), INTENT(OUT) :: FNAME

        CHARACTER(1)   :: idx1_char, idx2_char
        INTEGER        :: idx1_value, idx2_value
        CHARACTER(256) :: tmp_name

        idx1_value = MAX(1, MIN(NO1, 9))
        idx2_value = MAX(1, MIN(NO2, 9))
        WRITE(idx1_char, '(I1)') idx1_value
        WRITE(idx2_char, '(I1)') idx2_value

        tmp_name = TRIM(PREFIX) // idx1_char // idx2_char // TRIM(EXT)
        IF (LEN(FNAME) >= LEN_TRIM(tmp_name)) THEN
            FNAME = TRIM(tmp_name)
        ELSE
            FNAME = tmp_name(1:LEN(FNAME))
        END IF
    END SUBROUTINE CFNAME_GREEN
    !     !***************************************************************
!CFNAME_FRECHET creates filenames for the Frechet derivatives export
        SUBROUTINE CFNAME_FRECHET(PREFIX, NO, PARAM, FREQ, EXT, FNAME)
          IMPLICIT NONE
          CHARACTER(*), INTENT(IN)  :: PREFIX, PARAM, EXT
          INTEGER,      INTENT(IN)  :: NO
          REAL(dp),     INTENT(IN) :: FREQ
          CHARACTER(*), INTENT(OUT) :: FNAME

          CHARACTER(LEN=2)  :: idx_str
          CHARACTER(LEN=16) :: freq_str
          INTEGER           :: freq_pos, idx_value
          CHARACTER(LEN=256):: tmp_name

          idx_value = MAX(0, MIN(NO, 99))
          WRITE(idx_str, '(I2.2)') idx_value

          WRITE(freq_str, '(F6.2)') FREQ
          freq_str = ADJUSTL(freq_str)
          IF (FREQ >= 0.0_dp .AND. FREQ < 10.0_dp .AND. freq_str(1:1) /= '-') THEN
            freq_str = '0' // TRIM(freq_str)
          END IF
          freq_pos = INDEX(freq_str, '.')
          IF (freq_pos > 0) freq_str(freq_pos:freq_pos) = '_'

          tmp_name = TRIM(PREFIX) // idx_str // TRIM(PARAM) // '_' // TRIM(freq_str) // &
                     TRIM(EXT)
          IF (LEN(FNAME) >= LEN_TRIM(tmp_name)) THEN
            FNAME = TRIM(tmp_name)
          ELSE
            FNAME = tmp_name(1:LEN(FNAME))
          END IF
        END SUBROUTINE CFNAME_FRECHET
   
!********************************************
SUBROUTINE CFNAME_MODELS(PREFIX, PARAM, FREQ, ITER, EXT, FNAME20)
  IMPLICIT NONE
  ! inputs
  CHARACTER(*), INTENT(IN) :: PREFIX, PARAM, EXT
  REAL(dp),     INTENT(IN) :: FREQ
  INTEGER,      INTENT(IN) :: ITER
  ! output
  CHARACTER(*), INTENT(OUT) :: FNAME20

  ! locals
  CHARACTER(LEN=16) :: freq_str
  CHARACTER(LEN=16) :: iter_str
  CHARACTER(LEN=256):: tmp_name

  WRITE(freq_str,'(F8.2)') FREQ
  freq_str = ADJUSTL(freq_str)
  IF (FREQ >= 0.0_dp .AND. FREQ < 10.0_dp .AND. freq_str(1:1) /= '-') THEN
    freq_str = '0' // TRIM(freq_str)
  END IF

  WRITE(iter_str,'(I0)') ITER   ! minimal-width integer; Intel compilers support I0
  IF (ITER >= 0 .AND. ITER < 10) THEN
    iter_str = '0' // TRIM(iter_str)
  END IF

  tmp_name = TRIM(PREFIX)//TRIM(PARAM)//'_'//TRIM(freq_str)//'_IT'// &
             TRIM(iter_str)//TRIM(EXT)

  IF (LEN(FNAME20) >= LEN_TRIM(tmp_name)) THEN
    FNAME20 = TRIM(tmp_name)
  ELSE
    FNAME20 = tmp_name(1:LEN(FNAME20))
  END IF

  RETURN
END SUBROUTINE CFNAME_MODELS

    SUBROUTINE GRID2D_OUT0(FNAME,NX,NZ,X,Z,F,NTO,XTO,ZTO)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: FNAME
        INTEGER, INTENT(IN) :: NX, NZ, NTO
        REAL(dp), INTENT(IN) ::X(:),Z(:),F(:),XTO(:),ZTO(:)
        REAL(dp),ALLOCATABLE :: A(:),ZM(:)
        CHARACTER (len=256) :: errmsg, iomsg
        INTEGER :: istat, unit_out
        INTEGER :: MAXPT, NNX, NNZ, IP, I, J, II, JJ, I01, I02, J01, J02, J03, J04
        REAL(dp) :: FAIR, EPS, EEE, XA, XB, ZA, ZB, DM, XL, ZL, RL, DX, DZ, B0
        REAL(dp) :: XI, ZI, Z0, RX, XX, ZJ, ZZ, RR, W1, W2, W3, W4, RZ
    
        ! write(*,*)'inside grid2DOut'
        MAXPT=300      ! <=500 points
        FAIR =1.0e20_dp   ! Air treatment
        EPS=1.0e-6_dp     ! Mmin. distance
        EEE=1.0e-50_dp    ! Equivalent to Zero 
        
    !---- determine z-range --------      
        XA=+1.0e10_dp
        XB=-1.0e10_dp
        ZA=+1.0e10_dp
        ZB=-1.0e10_dp
        DM=0.0_dp
        
        IP=0
        DO I=1,NX
        XI=X(I)
        
        IF(I.LT.NX)THEN
        DX=X(I+1)-X(I)
        IF(DX.GT.DM)DM=DX
        ENDIF
        
        IF(XI.LT.XA)XA=XI
        IF(XI.GT.XB)XB=XI
        DO J=1,NZ
        IP=IP+1
        ZI=Z(IP)
        
        IF(J.LT.NZ)THEN
        DZ=Z(IP+1)-Z(IP)
        IF(DZ.GT.DM)DM=DZ
        ENDIF
        
        IF(ZI.LT.ZA)ZA=ZI
        IF(ZI.GT.ZB)ZB=ZI
        ENDDO
        ENDDO
        
        XL=XB-XA
        ZL=ZB-ZA
        ! WRITE(*,*) 'DBG local extents: XA=', XA, ' XB=', XB, ' ZA=', ZA, ' ZB=', ZB, ' DM=', DM
        RL=MAX(XL,ZL)
        NNX=INT((XL/RL)*REAL(MAXPT,dp))+1
        NNZ=INT((ZL/RL)*REAL(MAXPT,dp))+1      
        DX=(XB-XA)/REAL(NNX-1,dp)
        DZ=(ZB-ZA)/REAL(NNZ-1,dp)
        DM=0.5_dp*(DX+DZ)
        ! WRITE(*,*) 'DBG local grid: XA=', XA, ' XB=', XB, ' ZA=', ZA, ' ZB=', ZB
        ! WRITE(*,*) 'DBG local grid: RL=', RL, ' NNX=', NNX, ' NNZ=', NNZ, ' DX=', DX, ' DZ=', DZ, ' DM=', DM
        
        ALLOCATE (A(NNZ),ZM(NNZ))
        OPEN(newunit=unit_out, FILE=FNAME, STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=istat, IOMSG=iomsg)
        IF (istat /= 0) THEN
            SELECT CASE (istat)
                CASE (1001)
                    errmsg = "File not found."
                CASE (1002)
                    errmsg = "Permission denied."
                CASE DEFAULT
                    errmsg = TRIM(iomsg)
            END SELECT
            WRITE(error_unit,*) "Error opening file: ", TRIM(FNAME)
            WRITE(error_unit,*) TRIM(errmsg)
            STOP "File open failure"
        END IF
        
    !---- grid the domain ----------
        DO J=1,NNZ
        ZM(J)=ZA+REAL(J-1,dp)*DZ
        ENDDO
        
        B0=REAL(NNZ,dp)
        WRITE(unit_out,12)B0,(ZM(J),J=1,NNZ)
        
        DO 10 I=1,NNX
        XI=XA+REAL(I-1,dp)*DX
        Z0=ZH(NTO,XTO,ZTO,XI)
        
        RX=1.0e10_dp
        DO 2 II=1,NX
        XX=ABS(XI-X(II))
        IF(XX.GT.RX)GO TO 2
        RX=XX
        I01=II
     2  CONTINUE
        
        RX=1.0e10_dp
        DO 3 II=1,NX
        IF(II.EQ.I01)GO TO 3
        XX=ABS(XI-X(II))
        IF(XX.GT.RX)GO TO 3
        RX=XX
        I02=II
     3  CONTINUE
        
        DO 8 J=1,NNZ
        ZJ=ZM(J)
        A(J)=FAIR   ! for air treatment
     
        IF(ZJ.GT.Z0)GO TO 8
        RZ=1.0e10_dp
        DO 4 JJ=1,NZ
        IP=(I01-1)*NZ+JJ
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 4
        RZ=ZZ
        J01=IP
     4  CONTINUE
     
        RZ=1.0e10_dp
        DO 5 JJ=1,NZ
        IP=(I01-1)*NZ+JJ
        IF(IP.EQ.J01)GO TO 5
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 5
        RZ=ZZ
        J02=IP
     5  CONTINUE
        
        RZ=1.0e10_dp
        DO 6 JJ=1,NZ
        IP=(I02-1)*NZ+JJ
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 6
        RZ=ZZ
        J03=IP
     6  CONTINUE
     
        RZ=1.0e10_dp
        DO 7 JJ=1,NZ
        IP=(I02-1)*NZ+JJ 
        IF(IP.EQ.J03)GO TO 7
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 7
        RZ=ZZ
        J04=IP
     7  CONTINUE
        
        XX=XI-X(I01)
        ZZ=ZJ-Z(J01)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(J)=F(J01)
        GO TO 8
        ELSE
        W1=1.0_dp/RR
        ENDIF
        
        ZZ=ZJ-Z(J02)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(J)=F(J02)
        GO TO 8
        ELSE
        W2=1.0_dp/RR
        ENDIF
        
        XX=XI-X(I02)
        ZZ=ZJ-Z(J03)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(J)=F(J03)
        GO TO 8
        ELSE
        W3=1.0_dp/RR
        ENDIF
        
        ZZ=ZJ-Z(J04)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(J)=F(J04)
        GO TO 8
        ELSE
        W4=1.0_dp/RR
        ENDIF
        
        A(J)=(W1*F(J01)+W2*F(J02)+W3*F(J03)+W4*F(J04))/(W1+W2+W3+W4)
        IF(ABS(A(J)).LE.EEE)A(J)=0.0_dp
     8  CONTINUE
     
        WRITE(unit_out,12)XI,(A(J),J=1,NNZ)     

    10  CONTINUE
        CLOSE(unit_out)
        
    12  FORMAT(502(E13.6,4X))
        DEALLOCATE (A,ZM)
        RETURN
        END     SUBROUTINE GRID2D_OUT0
    
      !----------------------------------------------------------------------C
      !                                                                      C
      !     grid 2D scattered data for viewing the results                   C
      !                                                                      C
      !     Entries:                                                         C
      !                                                                      C
      !      (1) X(NX),Z(NX*NZ),F(NX*NZ)......2D scattered points & values;  C
      !      (2) XTO(NTO),ZTO(NTO)......2D topography curve on the surface;
      !      (3) F......................data to be plotted;
      !                                                                      C
      !----------------------------------------------------------------------C
    SUBROUTINE GRID2D_OUT(FNAME,NX,NZ,X,Z,F,NTO,XTO,ZTO,IE0,IS0)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: FNAME
        INTEGER, INTENT(IN) :: NX, NZ, NTO
        INTEGER, OPTIONAL, INTENT(IN) :: IE0, IS0
        REAL(dp), INTENT(IN) ::X(:),Z(:),F(:),XTO(:),ZTO(:)
        REAL(dp),ALLOCATABLE :: A(:),ZM(:)
        CHARACTER (len=256) :: errmsg, iomsg
        INTEGER :: istat, unit_out
        INTEGER :: MAXPT, NNX, NNZ, IP, I, J, II, JJ, I01, I02, J01, J02, J03, J04
        INTEGER :: JP1, JP2, IP1, IP2, TRIM_NNX, TRIM_NNZ, JOUT
        REAL(dp) :: FAIR, EPS, EEE, XA, XB, ZA, ZB, DM, XL, ZL, RL, DX, DZ, B0
        REAL(dp) :: XI, ZI, Z0, RX, XX, ZJ, ZZ, RR, W1, W2, W3, W4, RZ
        LOGICAL :: do_trim
    
        ! write(*,*)'inside grid2DOut'
        MAXPT=300         ! Target max grid size for output (fixed cap).
        FAIR =1.0e20_dp   ! Air treatment value used above topography.
        EPS=1.0e-6_dp     ! Minimum distance to avoid divide-by-zero.
        EEE=1.0e-50_dp    ! Equivalent-to-zero threshold.
        
    !---- determine data bounds and spacing (scattered input) --------
        XA=+1.0e10_dp
        XB=-1.0e10_dp
        ZA=+1.0e10_dp
        ZB=-1.0e10_dp
        DM=0.0_dp
        
        IP=0
        DO I=1,NX
        XI=X(I)
        
        IF(I.LT.NX)THEN
        DX=X(I+1)-X(I)
        IF(DX.GT.DM)DM=DX
        ENDIF
        
        IF(XI.LT.XA)XA=XI
        IF(XI.GT.XB)XB=XI
        DO J=1,NZ
        IP=IP+1
        ZI=Z(IP)
        
        IF(J.LT.NZ)THEN
        DZ=Z(IP+1)-Z(IP)
        IF(DZ.GT.DM)DM=DZ
        ENDIF
        
        IF(ZI.LT.ZA)ZA=ZI
        IF(ZI.GT.ZB)ZB=ZI
        ENDDO
        ENDDO
        
        XL=XB-XA
        ZL=ZB-ZA
        ! WRITE(*,*) 'DBG local extents: XA=', XA, ' XB=', XB, ' ZA=', ZA, ' ZB=', ZB, ' DM=', DM
        RL=MAX(XL,ZL)
        ! Scale output grid size to the larger axis, capped by MAXPT.
        NNX=INT((XL/RL)*REAL(MAXPT,dp))+1
        NNZ=INT((ZL/RL)*REAL(MAXPT,dp))+1      
        DX=(XB-XA)/REAL(NNX-1,dp)
        DZ=(ZB-ZA)/REAL(NNZ-1,dp)
        DM=0.5_dp*(DX+DZ)
        ! WRITE(*,*) 'DBG local grid: XA=', XA, ' XB=', XB, ' ZA=', ZA, ' ZB=', ZB
        ! WRITE(*,*) 'DBG local grid: RL=', RL, ' NNX=', NNX, ' NNZ=', NNZ, ' DX=', DX, ' DZ=', DZ, ' DM=', DM

        ! Optional trimming of PML/absorbing zones:
        ! IS0 = 2 -> trim top and bottom; IS0 = 1 -> trim bottom only.
        ! When not provided, output the full grid.
        do_trim = PRESENT(IE0) .AND. PRESENT(IS0)
        IF (do_trim) THEN
          IF (IS0 == 2) THEN
            JP1 = NINT(IE0 * REAL(NNZ, dp) / REAL(NZ, dp)) + 1
            JP2 = NNZ - (JP1 - 1)
          ELSEIF (IS0 == 1) THEN
            JP1 = 1
            JP2 = NNZ - NINT(IE0 * REAL(NNZ, dp) / REAL(NZ, dp))
          ELSE
            JP1 = 1
            JP2 = NNZ
          ENDIF
          TRIM_NNZ = JP2 - JP1 + 1

          IP1 = NINT(IE0 * REAL(NNX, dp) / REAL(NZ, dp)) + 1
          IP2 = NNX - (IP1 - 1)
          TRIM_NNX = IP2 - IP1 + 1
        ELSE
          JP1 = 1
          JP2 = NNZ
          TRIM_NNZ = NNZ
          IP1 = 1
          IP2 = NNX
          TRIM_NNX = NNX
        ENDIF

        ! Allocate trimmed vertical output buffers.
        ALLOCATE (A(TRIM_NNZ),ZM(TRIM_NNZ))
        OPEN(newunit=unit_out, FILE=FNAME, STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=istat, IOMSG=iomsg)
        IF (istat /= 0) THEN
            SELECT CASE (istat)
                CASE (1001)
                    errmsg = "File not found."
                CASE (1002)
                    errmsg = "Permission denied."
                CASE DEFAULT
                    errmsg = TRIM(iomsg)
            END SELECT
            WRITE(error_unit,*) "Error opening file: ", TRIM(FNAME)
            WRITE(error_unit,*) TRIM(errmsg)
            STOP "File open failure"
        END IF
        
    !---- build output vertical coordinate array ----------
        DO J=JP1,JP2
        ZM(J-JP1+1)=ZA+REAL(J-1,dp)*DZ
        ENDDO
        
        B0=REAL(TRIM_NNZ,dp)
        WRITE(unit_out,12)B0,(ZM(J),J=1,TRIM_NNZ)
        
        !---- interpolate per output X position ----------
        DO 10 I=IP1,IP2
        XI=XA+REAL(I-1,dp)*DX
        ! Topography height at XI (values above get FAIR).
        Z0=ZH(NTO,XTO,ZTO,XI)
        
        ! Find two nearest X columns in the input grid.
        RX=1.0e10_dp
        DO 2 II=1,NX
        XX=ABS(XI-X(II))
        IF(XX.GT.RX)GO TO 2
        RX=XX
        I01=II
     2  CONTINUE
        
        RX=1.0e10_dp
        DO 3 II=1,NX
        IF(II.EQ.I01)GO TO 3
        XX=ABS(XI-X(II))
        IF(XX.GT.RX)GO TO 3
        RX=XX
        I02=II
     3  CONTINUE
        
        DO 8 J=JP1,JP2
        ZJ=ZA+REAL(J-1,dp)*DZ
        JOUT=J-JP1+1
        A(JOUT)=FAIR   ! for air treatment
     
        ! Skip points above the topography surface.
        IF(ZJ.GT.Z0)GO TO 8
        ! Find the two nearest Z points in each of the two nearest X columns.
        RZ=1.0e10_dp
        DO 4 JJ=1,NZ
        IP=(I01-1)*NZ+JJ
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 4
        RZ=ZZ
        J01=IP
     4  CONTINUE
     
        RZ=1.0e10_dp
        DO 5 JJ=1,NZ
        IP=(I01-1)*NZ+JJ
        IF(IP.EQ.J01)GO TO 5
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 5
        RZ=ZZ
        J02=IP
     5  CONTINUE
        
        RZ=1.0e10_dp
        DO 6 JJ=1,NZ
        IP=(I02-1)*NZ+JJ
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 6
        RZ=ZZ
        J03=IP
     6  CONTINUE
     
        RZ=1.0e10_dp
        DO 7 JJ=1,NZ
        IP=(I02-1)*NZ+JJ 
        IF(IP.EQ.J03)GO TO 7
        ZZ=ABS(ZJ-Z(IP))
        IF(ZZ.GT.RZ)GO TO 7
        RZ=ZZ
        J04=IP
     7  CONTINUE
        
        ! Inverse-distance weighting using the four nearest neighbors.
        XX=XI-X(I01)
        ZZ=ZJ-Z(J01)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(JOUT)=F(J01)
        GO TO 8
        ELSE
        W1=1.0_dp/RR
        ENDIF
        
        ZZ=ZJ-Z(J02)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(JOUT)=F(J02)
        GO TO 8
        ELSE
        W2=1.0_dp/RR
        ENDIF
        
        XX=XI-X(I02)
        ZZ=ZJ-Z(J03)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(JOUT)=F(J03)
        GO TO 8
        ELSE
        W3=1.0_dp/RR
        ENDIF
        
        ZZ=ZJ-Z(J04)
        RR=SQRT(XX*XX+ZZ*ZZ)
        IF(RR.LE.EPS)THEN
        A(JOUT)=F(J04)
        GO TO 8
        ELSE
        W4=1.0_dp/RR
        ENDIF
        
        A(JOUT)=(W1*F(J01)+W2*F(J02)+W3*F(J03)+W4*F(J04))/(W1+W2+W3+W4)
        IF(ABS(A(JOUT)).LE.EEE)A(JOUT)=0.0_dp
     8  CONTINUE
     
        ! Write one output row: X value followed by vertical profile.
        WRITE(unit_out,12)XI,(A(J),J=1,TRIM_NNZ)     

    10  CONTINUE
        CLOSE(unit_out)
        
    12  FORMAT(502(E13.6,4X))
        DEALLOCATE (A,ZM)
        RETURN
        END     SUBROUTINE GRID2D_OUT

    
     !----------------------------------------------------------------------C
      !                                                                      C
      !     grid 2D scattered data for viewing the results                   C
      !                                                                      C
      !     Entries:                                                         C
      !                                                                      C
      !      (1) X(NX),Z(NX*NZ),F(NX*NZ)......2D scattered points & values;  C
      !      (2) XTO(NTO),ZTO(NTO)......2D topography curve on the surface;  C
      !                                                                      C
      !----------------------------------------------------------------------C
        SUBROUTINE GRID2D_OUT2(FNAME, NX, NZ, X, Z, F, NTO, XTO, ZTO)
          use iso_fortran_env, only : dp => real64
          IMPLICIT NONE
          CHARACTER(LEN=*), INTENT(IN) :: FNAME
          INTEGER,      INTENT(IN) :: NX, NZ, NTO
          REAL(dp),     INTENT(IN) :: X(:), Z(:), F(:), XTO(:), ZTO(:)

          REAL(dp),    ALLOCATABLE :: A(:), ZM(:)
          INTEGER :: MAXPT, NNX, NNZ, unit_out2, istat
          CHARACTER(len=256) :: iomsg
          REAL(dp) :: FAIR, EPS, EEE, XA, XB, ZA, ZB, DM
          REAL(dp) :: XL, ZL, RL, DX, DZ, B0, XI, ZI, Z0, RX,RZ, XX, ZJ, ZZ, RR
          REAL(dp) :: W1, W2, W3, W4
          INTEGER :: IP, I, J, II, JJ, I01, I02, J01, J02, J03, J04
            
          MAXPT = MAX(NX - 1, NZ - 1)
          FAIR  = aver(F, NX * NZ)   ! Air treatment
          EPS   = 1.0e-6_dp          ! Min. distance
          EEE   = 1.0e-50_dp         ! Equivalent to zero 
            
          ! Determine z-range
          XA = 1.0e10_dp
          XB = -1.0e10_dp
          ZA = 1.0e10_dp
          ZB = -1.0e10_dp
          DM = 0.0_dp
            
          IP = 0
          DO I = 1, NX
              XI = X(I)
              IF (I < NX) THEN
                  DX = X(I+1) - X(I)
                  IF (DX > DM) DM = DX
              END IF
            
              IF (XI < XA) XA = XI
              IF (XI > XB) XB = XI
              DO J = 1, NZ
                  IP = IP + 1
                  ZI = Z(IP)
            
                  IF (J < NZ) THEN
                      DZ = Z(IP+1) - Z(IP)
                      IF (DZ > DM) DM = DZ
                  END IF
            
                  IF (ZI < ZA) ZA = ZI
                  IF (ZI > ZB) ZB = ZI
              END DO
          END DO
            

          ! Proceed with spacing and resolution
          XL = XB - XA
          ZL = ZB - ZA
          RL = MAX(XL, ZL)
              
          NNX = INT((XL / RL) * REAL(MAXPT, dp)) + 1
          NNZ = INT((ZL / RL) * REAL(MAXPT, dp)) + 1      
          DX = (XB - XA) / REAL(NNX - 1, dp)
          DZ = (ZB - ZA) / REAL(NNZ - 1, dp)
          DM = 0.5_dp * (DX + DZ)

          ALLOCATE (A(NNZ), ZM(NNZ))

          OPEN(newunit=unit_out2,FILE=FNAME,STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=istat, IOMSG=iomsg)
          IF (istat /= 0) THEN
            WRITE(error_unit,*) "Error opening file: ", TRIM(FNAME)
            WRITE(error_unit,*) TRIM(iomsg)
            RETURN
          END IF
            
        !---- grid the domain ----------
          DO J = 1, NNZ
            ZM(J) = ZA + REAL(J-1, dp) * DZ
          ENDDO
            
          B0 = REAL(NNZ, dp)
          WRITE(unit_out2,12) B0, (ZM(J), J=1, NNZ)
            
          DO 10 I = 1, NNX
            XI = XA + REAL(I-1, dp) * DX
            Z0 = ZH(NTO, XTO, ZTO, XI)
            
            RX = 1.0e10_dp
            DO 2 II=1,NX
              XX = ABS(XI - X(II))
              IF (XX > RX) GO TO 2
              RX  = XX
              I01 = II
         2  CONTINUE
            
            RX = 1.0e10_dp
            DO 3 II=1,NX
              IF (II == I01) GO TO 3
              XX = ABS(XI - X(II))
              IF (XX > RX) GO TO 3
              RX  = XX
              I02 = II
         3  CONTINUE
            
            DO 8 J=1,NNZ
              ZJ = ZM(J)
              A(J) = FAIR   ! for air treatment
         
              IF (ZJ > Z0) GO TO 8
              RZ = 1.0e10_dp
              DO 4 JJ=1,NZ
                IP = (I01-1)*NZ + JJ
                ZZ = ABS(ZJ - Z(IP))
                IF (ZZ > RZ) GO TO 4
                RZ  = ZZ
                J01 = IP
         4    CONTINUE
         
              RZ = 1.0e10_dp
              DO 5 JJ=1,NZ
                IP = (I01-1)*NZ + JJ
                IF (IP == J01) GO TO 5
                ZZ = ABS(ZJ - Z(IP))
                IF (ZZ > RZ) GO TO 5
                RZ  = ZZ
                J02 = IP
         5    CONTINUE
            
              RZ = 1.0e10_dp
              DO 6 JJ=1,NZ
                IP = (I02-1)*NZ + JJ
                ZZ = ABS(ZJ - Z(IP))
                IF (ZZ > RZ) GO TO 6
                RZ  = ZZ
                J03 = IP
         6    CONTINUE
         
              RZ = 1.0e10_dp
              DO 7 JJ=1,NZ
                IP = (I02-1)*NZ + JJ 
                IF (IP == J03) GO TO 7
                ZZ = ABS(ZJ - Z(IP))
                IF (ZZ > RZ) GO TO 7
                RZ  = ZZ
                J04 = IP
         7    CONTINUE
            
              XX = XI - X(I01)
              ZZ = ZJ - Z(J01)
              RR = SQRT(XX*XX + ZZ*ZZ)
              IF (RR <= EPS) THEN
                A(J) = F(J01)
                GO TO 8
              ELSE
                W1 = 1.0_dp / RR
              ENDIF
            
              ZZ = ZJ - Z(J02)
              RR = SQRT(XX*XX + ZZ*ZZ)
              IF (RR <= EPS) THEN
                A(J) = F(J02)
                GO TO 8
              ELSE
                W2 = 1.0_dp / RR
              ENDIF
            
              XX = XI - X(I02)
              ZZ = ZJ - Z(J03)
              RR = SQRT(XX*XX + ZZ*ZZ)
              IF (RR <= EPS) THEN
                A(J) = F(J03)
                GO TO 8
              ELSE
                W3 = 1.0_dp / RR
              ENDIF
            
              ZZ = ZJ - Z(J04)
              RR = SQRT(XX*XX + ZZ*ZZ)
              IF (RR <= EPS) THEN
                A(J) = F(J04)
                GO TO 8
              ELSE
                W4 = 1.0_dp / RR
              ENDIF
            
              A(J)=(W1*F(J01)+W2*F(J02)+W3*F(J03)+W4*F(J04))/(W1+W2+W3+W4)
              IF (ABS(A(J)) <= EEE) A(J)=0.0_dp
         8  CONTINUE
         
            WRITE(unit_out2,12) XI,(A(J),J=1,NNZ)     
    10  CONTINUE
          CLOSE(unit_out2)
            
        12  FORMAT(502(E13.6,4X))
          DEALLOCATE (A,ZM)
          RETURN
        END SUBROUTINE GRID2D_OUT2     

!-----------------------------------------------------------------------
!
!     GRID2D_OUT_FIXED interpolates scattered (X,Z) data onto a user-defined
!     physical grid (XMINC..XMAXC, ZMINC..ZMAXC with NNX_OUT/NNZ_OUT points).
!     The input mesh is assumed to be in shifted coordinates. To avoid mixing
!     physical points with PML-shifted depths, the routine first infers the
!     interior (non-PML) shifted window, then maps physical targets into that
!     window before interpolation.
!
!     Inputs:
!       FNAME............... Output filename
!       NX, NZ.............. Input mesh sizes for X/Z columns
!       X(NX), Z(NX*NZ)...... Shifted input coordinates (flattened Z)
!       F(NX*NZ)............ Field values at input grid points
!       NTO, XTO, ZTO........Shifted topography for ZH() masking
!       XMINC, XMAXC.........Physical output bounds in X
!       ZMINC, ZMAXC.........Physical output bounds in Z
!       NNX_OUT, NNZ_OUT.....Output grid point counts in X and Z
!       IE0, IS0.............PML settings (IS0=1 bottom only, IS0=2 top+bottom)
!
!     Output:
!       File FNAME written with physical coordinates; values above topography
!       are set to FAIR (average value).
!
!-----------------------------------------------------------------------
  SUBROUTINE GRID2D_OUT_FIXED(FNAME, NX, NZ, X, Z, F, NTO, XTO, ZTO, &
                              XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    INTEGER, INTENT(IN) :: NX, NZ, NTO, NNX_OUT, NNZ_OUT, IE0, IS0
    REAL(dp), INTENT(IN) :: X(:), Z(:), F(:), XTO(:), ZTO(:)
    REAL(dp), INTENT(IN) :: XMINC, XMAXC, ZMINC, ZMAXC

    REAL(dp), ALLOCATABLE :: A(:), ZM(:)
    INTEGER :: unit_out, istat
    CHARACTER(len=256) :: iomsg
    REAL(dp) :: FAIR, EPS, EEE, DX, DZ
    REAL(dp) :: XA, XB, ZA, ZB, XL_phys, ZL_phys, XL_shift, ZL_shift
    REAL(dp) :: excess_x, excess_z, XA_phys_shift, ZA_phys_shift
    REAL(dp) :: XI_phys, ZJ_phys, XI, ZJ, Z0
    REAL(dp) :: RX, RZ, XX, ZZ, RR, W1, W2, W3, W4
    INTEGER :: I, J, II, JJ, IP, I01, I02, J01, J02, J03, J04

    IF (NNX_OUT < 2 .OR. NNZ_OUT < 2) RETURN
    IF (XMAXC <= XMINC .OR. ZMAXC <= ZMINC) RETURN

    ! FAIR = aver(F, NX * NZ)   ! Air treatment value
    FAIR =1.0e20_dp
    EPS  = 1.0e-6_dp
    EEE  = 1.0e-50_dp

    DX = (XMAXC - XMINC) / REAL(NNX_OUT - 1, dp)
    DZ = (ZMAXC - ZMINC) / REAL(NNZ_OUT - 1, dp)

    ! Infer shifted input extents and interior (non-PML) shifted window.
    XA = MINVAL(X(1:NX)); XB = MAXVAL(X(1:NX))
    ZA = MINVAL(Z(1:NX*NZ)); ZB = MAXVAL(Z(1:NX*NZ))
    XL_phys = XMAXC - XMINC
    ZL_phys = ZMAXC - ZMINC
    XL_shift = XB - XA
    ZL_shift = ZB - ZA
    excess_x = MAX(0.0_dp, XL_shift - XL_phys)
    excess_z = MAX(0.0_dp, ZL_shift - ZL_phys)

    ! Side PML is symmetric in x.
    XA_phys_shift = XA + 0.5_dp*excess_x

    ! Vertical crop depends on absorption mode.
    IF (IS0 == 1) THEN
      ! Bottom PML only: physical window is the top part of shifted Z range.
      ZA_phys_shift = ZB - ZL_phys
    ELSEIF (IS0 == 2) THEN
      ! Top+bottom PML: physical window centered in shifted Z range.
      ZA_phys_shift = ZA + 0.5_dp*excess_z
    ELSE
      ZA_phys_shift = ZA + 0.5_dp*excess_z
    END IF
    ! Clamp for numerical safety.
    XA_phys_shift = MIN(MAX(XA_phys_shift, XA), XB - XL_phys)
    ZA_phys_shift = MIN(MAX(ZA_phys_shift, ZA), ZB - ZL_phys)

    ALLOCATE (A(NNZ_OUT), ZM(NNZ_OUT))
    OPEN(newunit=unit_out, FILE=FNAME, STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=istat, IOMSG=iomsg)
    IF (istat /= 0) THEN
      WRITE(error_unit,*) "Error opening file: ", TRIM(FNAME)
      WRITE(error_unit,*) TRIM(iomsg)
      DEALLOCATE(A, ZM)
      RETURN
    END IF

    DO J = 1, NNZ_OUT
      ZM(J) = ZMINC + REAL(J - 1, dp) * DZ
    END DO

    WRITE(unit_out,12) REAL(NNZ_OUT, dp), (ZM(J), J=1, NNZ_OUT)

    DO I = 1, NNX_OUT
      XI_phys = XMINC + REAL(I - 1, dp) * DX
      XI = XA_phys_shift + (XI_phys - XMINC)
      Z0 = ZH(NTO, XTO, ZTO, XI)

      ! Find two nearest X columns in the input grid (shifted coords).
      RX = 1.0e10_dp
      DO II = 1, NX
        XX = ABS(XI - X(II))
        IF (XX < RX) THEN
          RX  = XX
          I01 = II
        END IF
      END DO

      RX = 1.0e10_dp
      DO II = 1, NX
        IF (II == I01) CYCLE
        XX = ABS(XI - X(II))
        IF (XX < RX) THEN
          RX  = XX
          I02 = II
        END IF
      END DO

      DO J = 1, NNZ_OUT
        ZJ_phys = ZM(J)
        ZJ = ZA_phys_shift + (ZJ_phys - ZMINC)
        A(J) = FAIR
        IF (ZJ > Z0) CYCLE

        ! Find the two nearest Z points in each of the two nearest X columns.
        RZ = 1.0e10_dp
        DO JJ = 1, NZ
          IP = (I01 - 1) * NZ + JJ
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ  = ZZ
            J01 = IP
          END IF
        END DO

        RZ = 1.0e10_dp
        DO JJ = 1, NZ
          IP = (I01 - 1) * NZ + JJ
          IF (IP == J01) CYCLE
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ  = ZZ
            J02 = IP
          END IF
        END DO

        RZ = 1.0e10_dp
        DO JJ = 1, NZ
          IP = (I02 - 1) * NZ + JJ
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ  = ZZ
            J03 = IP
          END IF
        END DO

        RZ = 1.0e10_dp
        DO JJ = 1, NZ
          IP = (I02 - 1) * NZ + JJ
          IF (IP == J03) CYCLE
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ  = ZZ
            J04 = IP
          END IF
        END DO

        XX = XI - X(I01)
        ZZ = ZJ - Z(J01)
        RR = SQRT(XX*XX + ZZ*ZZ)
        IF (RR <= EPS) THEN
          A(J) = F(J01)
        ELSE
          W1 = 1.0_dp / RR
          ZZ = ZJ - Z(J02)
          W2 = 1.0_dp / SQRT(XX*XX + ZZ*ZZ)
          XX = XI - X(I02)
          ZZ = ZJ - Z(J03)
          W3 = 1.0_dp / SQRT(XX*XX + ZZ*ZZ)
          ZZ = ZJ - Z(J04)
          W4 = 1.0_dp / SQRT(XX*XX + ZZ*ZZ)
          A(J) = (W1 * F(J01) + W2 * F(J02) + W3 * F(J03) + W4 * F(J04)) / (W1 + W2 + W3 + W4)
          IF (ABS(A(J)) <= EEE) A(J) = 0.0_dp
        END IF
      END DO

      WRITE(unit_out,12) XI_phys, (A(J), J=1, NNZ_OUT)
    END DO

    CLOSE(unit_out)
12  FORMAT(502(E13.6,4X))
    DEALLOCATE(A, ZM)
  END SUBROUTINE GRID2D_OUT_FIXED

!-----------------------------------------------------------------------
!
!     GRID2D_OUT3 interpolates a scattered 2D field defined on irregular GQG
!     (X,Z) coordinates onto a regular vertical grid, trimming out the 
!     PML (absorbing) region based on the input IE0 and topography (XTO,ZTO).
!
!     Writes the output field to a file `FNAME` in a format suitable for
!     visualization tools. Interpolation is done using inverse-distance
!     weighting over the four nearest neighbors in the (X,Z) grid.
!       - The vertical extent is trimmed by IE0 to remove absorbing zones.
!       - The number of vertical grid points (NNZ) is auto-scaled from NZ.
!       - The interpolation uses 4-point inverse-distance weighting.
!       - The routine uses the ZH() function to determine surface height.
!
!     Inputs:
!       FNAME............... Output filename (CHARACTER*(*))
!       NX.................. Number of horizontal X locations
!       NZ.................. Number of vertical Z values per X location
!       X(NX)............... X-coordinates of the input grid
!       Z(NX*NZ)............ Z-coordinates of the input grid (flattened)
!       F(NX*NZ)............ Field values at input grid points
!       NTO................. Number of topography control points
!       XTO(NTO), ZTO(NTO).. Topography curve coordinates
!       IE0................. Number of vertical cells in the absorbing zone
!
!     Outputs:
!       File FNAME is created with gridded interpolated values. The first
!       row contains the vertical coordinate array. Each subsequent row
!       contains one vertical profile along X.
!

!
!     Internally Called Subroutines/Functions:
!       ZH.................. Topography height function at given X
!       AVER................ Computes average of input array F
!
!-----------------------------------------------------------------------
  SUBROUTINE GRID2D_OUT3(FNAME, NX, NZ, X, Z, F, NTO, XTO, ZTO, IE0, IS0)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    INTEGER, INTENT(IN) :: NX, NZ, NTO, IE0, IS0
    REAL(dp), INTENT(IN) :: X(:), Z(:), F(:), XTO(:), ZTO(:)

    REAL(dp), ALLOCATABLE :: A(:), ZM(:)
    INTEGER :: I, J, K, II, JJ, IP, I01, I02, J01, J02, J03, J04
    INTEGER :: TRIM_NNX, TRIM_NNZ, NNX, NNZ, JP1, JP2, IP1, IP2, MAXPT, unit_out3, istat
    INTEGER :: my_rank
    CHARACTER(len=256) :: iomsg
    REAL(dp) :: EPS, EEE, FAKEVAL, XL, ZL, RL, DX, DZ
    REAL(dp) :: XA, XB, ZA, ZB, XI, ZI, ZJ, XX, ZZ, Z0, RX, RZ, RR
    REAL(dp) :: W1, W2, W3, W4

    MAXPT = MAX(NX - 1, NZ - 1)
    FAKEVAL = AVER(F, NX * NZ)
    EPS = 1.0e-6_dp
    EEE = 1.0e-50_dp
    my_rank = 0

    ! Determine bounding box
    XA = 1.0e10_dp
    XB = -1.0e10_dp
    ZA = 1.0e10_dp
    ZB = -1.0e10_dp
    IP = 0
    DO I = 1, NX
      XI = X(I)
      XA = MIN(XA, XI)
      XB = MAX(XB, XI)
      DO J = 1, NZ
        IP = IP + 1
        ZI = Z(IP)
        ZA = MIN(ZA, ZI)
        ZB = MAX(ZB, ZI)
      ENDDO
    ENDDO

    ! Regular grid sizing
    XL = XB - XA
    ZL = ZB - ZA
    RL = MAX(XL, ZL)

    NNX = INT((XL / RL) * REAL(MAXPT, dp)) + 1
    NNZ = INT((ZL / RL) * REAL(MAXPT, dp)) + 1
    DX = XL / REAL(NNX - 1, dp)
    DZ = ZL / REAL(NNZ - 1, dp)

    ! Determine vertical crop based on IS0
    IF (IS0 == 2) THEN
      JP1 = NINT(IE0 * REAL(NNZ, dp) / REAL(NZ, dp)) + 1     ! Remove top PML
      JP2 = NNZ - (JP1 - 1)                          ! Remove bottom PML
    ELSEIF (IS0 == 1) THEN
      JP1 = 1                                        ! Keep top
      JP2 = NNZ - NINT(IE0 * REAL(NNZ, dp) / REAL(NZ, dp))   ! Remove bottom only
    ELSE
      JP1 = 1
      JP2 = NNZ
    ENDIF
    TRIM_NNZ = JP2 - JP1 + 1

    ! Trim horizontal (side) PMLs always
    IP1 = NINT(IE0 * REAL(NNX, dp) / REAL(NZ, dp)) + 1
    IP2 = NNX - (IP1 - 1)
    TRIM_NNX = IP2 - IP1 + 1

    ALLOCATE(A(TRIM_NNZ), ZM(TRIM_NNZ))
    unit_out3 = -1
    IF (my_rank == 0) THEN
      OPEN(newunit=unit_out3, FILE=FNAME, STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=istat, IOMSG=iomsg)
      IF (istat /= 0) THEN
        WRITE(error_unit,*) "Error opening file: ", TRIM(FNAME)
        WRITE(error_unit,*) TRIM(iomsg)
        DEALLOCATE(A, ZM)
        RETURN
      END IF
    END IF

    DO J = JP1, JP2
      ZM(J - JP1 + 1) = ZA + REAL(J - 1, dp) * DZ
    ENDDO
    IF (my_rank == 0)WRITE(unit_out3, 12) REAL(TRIM_NNZ, dp), (ZM(J), J=1, TRIM_NNZ)

    DO I = IP1, IP2
      XI = XA + REAL(I - 1, dp) * DX
      Z0 = ZH(NTO, XTO, ZTO, XI)

      ! Find nearest X points
      RX = 1.0e10_dp
      DO II = 1, NX
        XX = ABS(XI - X(II))
        IF (XX < RX) THEN
          RX = XX
          I01 = II
        ENDIF
      ENDDO

      RX = 1.0e10_dp
      DO II = 1, NX
        IF (II == I01) CYCLE
        XX = ABS(XI - X(II))
        IF (XX < RX) THEN
          RX = XX
          I02 = II
        ENDIF
      ENDDO

      DO JJ = JP1, JP2
        J = JJ
        ZJ = ZA + REAL(J - 1, dp) * DZ
        A(J - JP1 + 1) = FAKEVAL
        IF (ZJ > Z0) CYCLE

        ! Find interpolation points
        RZ = 1.0e10_dp
        DO K = 1, NZ
          IP = (I01 - 1) * NZ + K
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ = ZZ
            J01 = IP
          ENDIF
        ENDDO

        RZ = 1.0e10_dp
        DO K = 1, NZ
          IP = (I01 - 1) * NZ + K
          IF (IP == J01) CYCLE
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ = ZZ
            J02 = IP
          ENDIF
        ENDDO

        RZ = 1.0e10_dp
        DO K = 1, NZ
          IP = (I02 - 1) * NZ + K
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ = ZZ
            J03 = IP
          ENDIF
        ENDDO

        RZ = 1.0e10_dp
        DO K = 1, NZ
          IP = (I02 - 1) * NZ + K
          IF (IP == J03) CYCLE
          ZZ = ABS(ZJ - Z(IP))
          IF (ZZ < RZ) THEN
            RZ = ZZ
            J04 = IP
          ENDIF
        ENDDO

        ! Interpolation
        XX = XI - X(I01)
        ZZ = ZJ - Z(J01)
        RR = SQRT(XX**2 + ZZ**2)
        IF (RR <= EPS) THEN
          A(J - JP1 + 1) = F(J01)
        ELSE
          W1 = 1.0_dp / RR
          ZZ = ZJ - Z(J02)
          W2 = 1.0_dp / SQRT(XX**2 + ZZ**2)
          XX = XI - X(I02)
          ZZ = ZJ - Z(J03)
          W3 = 1.0_dp / SQRT(XX**2 + ZZ**2)
          ZZ = ZJ - Z(J04)
          W4 = 1.0_dp / SQRT(XX**2 + ZZ**2)
          A(J - JP1 + 1) = (W1 * F(J01) + W2 * F(J02) + W3 * F(J03) + W4 * F(J04)) / (W1 + W2 + W3 + W4)
          IF (ABS(A(J - JP1 + 1)) <= EEE) A(J - JP1 + 1) = 0.0_dp
        ENDIF
      ENDDO

     IF (my_rank == 0) WRITE(unit_out3, 12) XI, (A(J), J = 1, TRIM_NNZ)
    ENDDO

  IF (my_rank == 0 .AND. unit_out3 > 0)  CLOSE(unit_out3)
    DEALLOCATE(A, ZM)

12 FORMAT(502(E13.6, 4X))
RETURN
END SUBROUTINE GRID2D_OUT3

SUBROUTINE output_update(my_rank, NPAR, IANISO, INVP, CRR0, CRRU, CII0, CIIU, PARAM, &
                            FREQ, ITER, NPT, NBLOCK, NNX, NNZ, XP, ZP, NTO, XTO, ZTO, &
                            IE0, IS0, XMINC, XMAXC, ZMINC, ZMAXC, DEBUG_OUTPUT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: my_rank, NPAR, IANISO, ITER, NPT, NBLOCK, NNX, NNZ, NTO
      INTEGER, INTENT(IN) :: IE0, IS0
      REAL(dp), INTENT(IN) :: XMINC, XMAXC, ZMINC, ZMAXC
      INTEGER, INTENT(IN) :: INVP(:)
      REAL(dp), INTENT(IN) :: CRR0(:, :), CRRU(:, :), CII0(:, :), CIIU(:, :)
      REAL(dp), INTENT(IN) :: XP(:), ZP(:), XTO(:), ZTO(:)
      REAL(dp), INTENT(IN) :: FREQ
      CHARACTER(len=*), INTENT(IN) :: PARAM(:)
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      REAL(dp), ALLOCATABLE :: GR2(:), GR1(:)
      INTEGER :: I, IM
      INTEGER :: NNX_OUT, NNZ_OUT
      REAL(dp) :: DX_OUT
      CHARACTER(len=30) :: FNAME20
      LOGICAL :: dbg

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      IF (my_rank /= 0) RETURN

      NNX_OUT = NTO
      IF (NNX_OUT < 2) RETURN
      DX_OUT = (XMAXC - XMINC) / REAL(NNX_OUT - 1, dp)
      NNZ_OUT = 1 + NINT((ZMAXC - ZMINC)/DX_OUT + 1.0e-10_dp)

      IF (NNZ_OUT < 2) NNZ_OUT = 2

      ALLOCATE (GR2(NPT), GR1(NBLOCK))
      IM = 0

      DO I = 1, NPAR
         IF ((I <= IANISO) .AND. (INVP(I) == 1)) THEN
            IM = IM + 1
            GR2(1:NPT) = CRR0(I, 1:NPT)
            CALL CFNAME_MODELS('PCBB_{', PARAM(I), FREQ, ITER, '}.dat', FNAME20)
            IF (dbg) WRITE (*, *) 'Exported ', FNAME20
            CALL GRID2D_OUT_FIXED(FNAME20, NNX, NNZ, XP, ZP, GR2, NTO, XTO, ZTO, &
                                  XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
            IF (dbg) THEN
               GR2(1:NPT) = CRRU(I, 1:NPT)
               CALL CFNAME_MODELS('PCBA_{', PARAM(I), FREQ, ITER, '}.dat', FNAME20)
               CALL GRID2D_OUT2(FNAME20, NNX, NNZ, XP, ZP, GR2, NTO, XTO, ZTO)
            END IF

         ELSEIF ((I > IANISO) .AND. (INVP(I) == 1)) THEN
            IM = IM + 1
            GR2(1:NPT) = CII0(I - (IANISO - 1), 1:NPT)
            CALL CFNAME_MODELS('PCBB_{', PARAM(I), FREQ, ITER, '}.dat', FNAME20)
            IF (dbg) WRITE (*, *) 'Exported ', FNAME20
            CALL GRID2D_OUT_FIXED(FNAME20, NNX, NNZ, XP, ZP, GR2, NTO, XTO, ZTO, &
                                  XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
            IF (dbg) THEN
               GR2(1:NPT) = CIIU(I - (IANISO - 1), 1:NPT)
               CALL CFNAME_MODELS('PCBA_{', PARAM(I), FREQ, ITER, '}.dat', FNAME20)
               CALL GRID2D_OUT2(FNAME20, NNX, NNZ, XP, ZP, GR2, NTO, XTO, ZTO)
            END IF

         END IF
      END DO

      IF (ALLOCATED(GR2)) DEALLOCATE (GR2)
      IF (ALLOCATED(GR1)) DEALLOCATE (GR1)
   END SUBROUTINE output_update
   
SUBROUTINE STRIP_LEADING_SPACES(STR)
  IMPLICIT NONE
  CHARACTER*(*) STR
  INTEGER I

  DO I = 1, LEN(STR)
      IF (STR(I:I) .NE. ' ') THEN
          STR = STR(I:)
          RETURN
      END IF
  END DO
END SUBROUTINE STRIP_LEADING_SPACES



 SUBROUTINE SET_PARAMETERS(IANISO, IVISCO, ITHOM, PARAM, III)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IANISO, IVISCO, ITHOM
  CHARACTER(LEN=*), INTENT(OUT) :: PARAM(:)
  INTEGER, INTENT(OUT) :: III
            
            III = 0
            PARAM(1) = '{rho}'
          
            IF (IANISO .EQ. 3) THEN
              PARAM(2) = '{Vp}'
              PARAM(3) = '{Vs}'
              III = 1
              IF (IVISCO .EQ. 1) THEN
                PARAM(4) = '{Qvp}'
                PARAM(5) = '{Qvs}'
              ENDIF
            ENDIF
          
            IF (IANISO .EQ. 6) THEN
              IF (ITHOM .EQ. 1) THEN
                PARAM(2) = '{alp}'
                PARAM(3) = '{bet}'
                PARAM(4) = '{eps}'
                PARAM(5) = '{del}'
                PARAM(6) = '{gam}'
                III = 1
                IF (IVISCO .EQ. 1) THEN
                  PARAM(7) = '{Qal}'
                  PARAM(8) = '{Qbe}'
                  PARAM(9) = '{Qep}'
                  PARAM(10) = '{Qde}'
                  PARAM(11) = '{Qga}'
                ENDIF
              ELSE
                PARAM(2) = '{C11}'
                PARAM(3) = '{C13}'
                PARAM(4) = '{C33}'
                PARAM(5) = '{C44}'
                PARAM(6) = '{C66}'
                III = 1
                IF (IVISCO .EQ. 1) THEN
                  PARAM(7) = '{Q11}'
                  PARAM(8) = '{Q13}'
                  PARAM(9) = '{Q33}'
                  PARAM(10) = '{Q44}'
                  PARAM(11) = '{Q66}'
                ENDIF
              ENDIF
            ENDIF
          
            IF (IANISO .EQ. 7) THEN
              PARAM(7) = '{the}'
              IF (ITHOM .EQ. 1) THEN
                PARAM(2) = '{alp}'
                PARAM(3) = '{bet}'
                PARAM(4) = '{eps}'
                PARAM(5) = '{del}'
                PARAM(6) = '{gam}'
                III = 1
                PARAM(8) = '{Q11}'
                PARAM(9) = '{Q13}'
                PARAM(10) = '{Q33}'
                PARAM(11) = '{Q44}'
                PARAM(12) = '{Q66}'
              ELSE
                PARAM(2) = '{C11}'
                PARAM(3) = '{C13}'
                PARAM(4) = '{C33}'
                PARAM(5) = '{C44}'
                PARAM(6) = '{C66}'
                III = 1
                IF (IVISCO .EQ. 1) THEN
                  PARAM(8) = '{Q11}'
                  PARAM(9) = '{Q13}'
                  PARAM(10) = '{Q33}'
                  PARAM(11) = '{Q44}'
                  PARAM(12) = '{Q66}'
                ENDIF
              ENDIF
            ENDIF
          
            IF (IANISO .GT. 7) THEN
              PARAM(2) = '{c11}'
              PARAM(3) = '{c12}'
              PARAM(4) = '{c13}'
              PARAM(5) = '{c14}'
              PARAM(6) = '{c15}'
              PARAM(7) = '{c16}'
              PARAM(8) = '{c22}'
              PARAM(9) = '{c23}'
              PARAM(10) = '{c24}'
              PARAM(11) = '{c25}'
              PARAM(12) = '{c26}'
              PARAM(13) = '{c33}'
              PARAM(14) = '{c34}'
              PARAM(15) = '{c35}'
              PARAM(16) = '{c36}'
              PARAM(17) = '{c44}'
              PARAM(18) = '{c45}'
              PARAM(19) = '{c46}'
              PARAM(20) = '{c55}'
              PARAM(21) = '{c56}'
              PARAM(22) = '{c66}'
              III = 1
            ENDIF
          
            IF (III .EQ. 0) THEN
              WRITE(*,*) '              Wrong value of IANISO !!' 
              STOP
            ENDIF
          END SUBROUTINE SET_PARAMETERS
         
    end module output_mod
