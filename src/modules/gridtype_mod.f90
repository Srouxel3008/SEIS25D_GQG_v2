MODULE gridtype_mod
   USE shared_mod, ONLY: ZH, CDLI
   use iso_fortran_env, only: dp => real64

   implicit NONE

   TYPE :: PhysicalState
      ! ---- status ----
      LOGICAL :: is_frozen = .FALSE.

      ! ---- immutable bounds ----
      real(dp) :: XMINC = 0D0, XMAXC = 0D0
      real(dp) :: ZMINC = 0D0, ZMAXC = 0D0

      ! ---- immutable sizes / counters ----
      INTEGER :: MXC = 0          ! number of horizontal columns (SIZE(XM))
      INTEGER :: MZTC = 0          ! total number of ZM entries (SIZE(ZM))
      INTEGER :: NTOC = 0          ! number of topo points
      INTEGER :: NSRC = 0          ! number of SR points

      ! ---- immutable geometry/topography/source-receiver copies ----
      real(dp), ALLOCATABLE :: XMC(:)      ! length MXC
      INTEGER, ALLOCATABLE :: MZC(:)      ! per-column counts, length MXC
      real(dp), ALLOCATABLE :: ZMC(:)      ! concatenated Z depths, length MZC

      real(dp), ALLOCATABLE :: XTOC(:), ZTOC(:)  ! topo
      real(dp), ALLOCATABLE :: XSRC(:), YSRC(:), ZSRC(:)  ! src/rec

   END TYPE PhysicalState

   !-----------------------------------------------------------------
   ! Inversion grid / block geometry built from a specific GQG mesh
   !-----------------------------------------------------------------
   TYPE :: InversionGridType
      real(dp), ALLOCATABLE :: XBC(:)        ! (NX-1)
      real(dp), ALLOCATABLE :: ZBC(:)        ! (NBLOCK)
      real(dp), ALLOCATABLE :: DX_BLOCK(:)   ! (NBLOCK)
      real(dp), ALLOCATABLE :: C1_BLOCK(:)   ! (NBLOCK)
      INTEGER, ALLOCATABLE :: N0_BLOCK(:)   ! (NBLOCK)
      LOGICAL, ALLOCATABLE :: vPML(:)       ! (NBLOCK)
      real(dp), ALLOCATABLE :: Z1_OUT(:, :)  ! (NORD, NBLOCK)
      real(dp), ALLOCATABLE :: Z2_OUT(:, :)  ! (NORD, NBLOCK)
      real(dp), ALLOCATABLE :: T1_OUT(:, :)  ! (NORD, NBLOCK)
      real(dp), ALLOCATABLE :: T2_OUT(:, :)  ! (NORD, NBLOCK)
   END TYPE InversionGridType

CONTAINS

   SUBROUTINE Freeze_Phys_Grid(XMIN, XMAX, ZMIN, ZMAX, &
                               NTO, XTO, ZTO, XSR, YSR, ZSR, NSR, &
                               XM, ZM, MX, MZ, MZT, canon, msg, my_rank, DEBUG_OUTPUT)
      IMPLICIT NONE

      real(dp), INTENT(IN) :: XMIN, XMAX, ZMIN, ZMAX
      real(dp), INTENT(IN) :: XM(:), ZM(:), XTO(:), ZTO(:), XSR(:), YSR(:), ZSR(:)
      INTEGER, INTENT(IN) :: NTO, NSR, MZ(:), my_rank, MX, MZT
      TYPE(PhysicalState), INTENT(INOUT) :: canon
      CHARACTER(*), INTENT(OUT) :: msg
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      LOGICAL :: dbg

      dbg = .FALSE.; IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      msg = 'Freeze: input validation pending'

      IF ((SIZE(ZTO) /= NTO) .OR. (SIZE(XTO) /= NTO)) THEN
         msg = 'Freeze: XTO/ZTO mismatch'; RETURN
      END IF
      IF ((SIZE(ZSR) /= NSR) .OR. (SIZE(XSR) /= NSR) .OR. (SIZE(YSR) /= NSR)) THEN
         msg = 'Freeze: XSR/YSR/ZSR mismatch'; RETURN
      END IF
      IF (SIZE(XM) /= MX) THEN
         msg = 'Freeze: SIZE(XM) /= MX'; RETURN
      END IF
      IF (SIZE(MZ) /= MX) THEN
         msg = 'Freeze: SIZE(MZ) /= MX'; RETURN
      END IF
      IF (SIZE(ZM) /= MZT) THEN
         msg = 'Freeze: SIZE(ZM) /= MZT'; RETURN
      END IF
      IF (MZT /= SUM(MZ)) THEN
         msg = 'Freeze: MZT /= SUM(MZ)'; RETURN
      END IF

      ! ---- Store scalars ----
      canon%XMINC = XMIN; canon%XMAXC = XMAX
      canon%ZMINC = ZMIN; canon%ZMAXC = ZMAX
      canon%NTOC = NTO
      canon%MXC = MX
      canon%MZTC = MZT
      canon%NSRC = NSR

      ! ---- Deep-copy immutable originals (dealloc first if needed) ----
      IF (ALLOCATED(canon%XTOC)) DEALLOCATE (canon%XTOC)
      IF (ALLOCATED(canon%ZTOC)) DEALLOCATE (canon%ZTOC)
      IF (ALLOCATED(canon%XSRC)) DEALLOCATE (canon%XSRC)
      IF (ALLOCATED(canon%YSRC)) DEALLOCATE (canon%YSRC)
      IF (ALLOCATED(canon%ZSRC)) DEALLOCATE (canon%ZSRC)
      IF (ALLOCATED(canon%XMC)) DEALLOCATE (canon%XMC)
      IF (ALLOCATED(canon%MZC)) DEALLOCATE (canon%MZC)
      IF (ALLOCATED(canon%ZMC)) DEALLOCATE (canon%ZMC)

      ALLOCATE (canon%XTOC(NTO)); canon%XTOC = XTO
      ALLOCATE (canon%ZTOC(NTO)); canon%ZTOC = ZTO
      ALLOCATE (canon%XSRC(NSR)); canon%XSRC = XSR
      ALLOCATE (canon%YSRC(NSR)); canon%YSRC = YSR
      ALLOCATE (canon%ZSRC(NSR)); canon%ZSRC = ZSR

      ALLOCATE (canon%XMC(MX)); canon%XMC = XM
      ALLOCATE (canon%MZC(MX)); canon%MZC = MZ
      ALLOCATE (canon%ZMC(canon%MZTC)); canon%ZMC = ZM(1:canon%MZTC)

      canon%is_frozen = .TRUE.

      msg = 'Freeze: canonical geometry & bounds recorded'
      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, '(A,2F12.3,2F12.3,4I8)') 'Freeze: X/Z bounds=', &
            XMIN, XMAX, ZMIN, ZMAX, MX, SIZE(ZM), NTO, NSR
      END IF

   END SUBROUTINE Freeze_Phys_Grid

   SUBROUTINE Allocate_InversionGrid(IG, nxm1, NBLOCK, NORD)
      TYPE(InversionGridType), INTENT(INOUT) :: IG
      INTEGER, INTENT(IN) :: nxm1, NBLOCK, NORD

      IF (ALLOCATED(IG%XBC)) DEALLOCATE (IG%XBC)
      IF (ALLOCATED(IG%ZBC)) DEALLOCATE (IG%ZBC)
      IF (ALLOCATED(IG%DX_BLOCK)) DEALLOCATE (IG%DX_BLOCK)
      IF (ALLOCATED(IG%C1_BLOCK)) DEALLOCATE (IG%C1_BLOCK)
      IF (ALLOCATED(IG%N0_BLOCK)) DEALLOCATE (IG%N0_BLOCK)
      IF (ALLOCATED(IG%vPML)) DEALLOCATE (IG%vPML)
      IF (ALLOCATED(IG%Z1_OUT)) DEALLOCATE (IG%Z1_OUT)
      IF (ALLOCATED(IG%Z2_OUT)) DEALLOCATE (IG%Z2_OUT)
      IF (ALLOCATED(IG%T1_OUT)) DEALLOCATE (IG%T1_OUT)
      IF (ALLOCATED(IG%T2_OUT)) DEALLOCATE (IG%T2_OUT)

      ALLOCATE (IG%XBC(nxm1))
      ALLOCATE (IG%ZBC(NBLOCK), IG%DX_BLOCK(NBLOCK), IG%C1_BLOCK(NBLOCK))
      ALLOCATE (IG%N0_BLOCK(NBLOCK), IG%vPML(NBLOCK))
      ALLOCATE (IG%Z1_OUT(NORD, NBLOCK), IG%Z2_OUT(NORD, NBLOCK))
      ALLOCATE (IG%T1_OUT(NORD, NBLOCK), IG%T2_OUT(NORD, NBLOCK))
      IG%XBC = 0.0_dp
      IG%ZBC = 0.0_dp
      IG%DX_BLOCK = 0.0_dp
      IG%C1_BLOCK = 0.0_dp
      IG%N0_BLOCK = 0
      IG%vPML = .FALSE.
      IG%Z1_OUT = 0.0_dp
      IG%Z2_OUT = 0.0_dp
      IG%T1_OUT = 0.0_dp
      IG%T2_OUT = 0.0_dp
   END SUBROUTINE Allocate_InversionGrid

   SUBROUTINE Precompute_InversionGrid(IG, NX, NZ, NNZ, NBLOCK, NORD, IE0, DZ0, AS, &
                                       NTO, XTO, ZTO, IS0, X, my_rank)

      implicit none

      TYPE(InversionGridType), INTENT(INOUT) :: IG
      INTEGER, INTENT(IN) :: NX, NZ, NNZ, NBLOCK, NORD, IE0, NTO, IS0, my_rank
      real(dp), INTENT(IN) :: DZ0, AS(:), X(:), XTO(:), ZTO(:)

      INTEGER :: nxm1, nzm1
      INTEGER :: I, J, K, L, IM, MID, N0
      real(dp) :: X1, X2, DX, C1, XI, Z0, DM, SS
      real(dp) :: Z1(NORD), Z2(NORD), T1(NORD), T2(NORD)
      real(dp) :: XP_loc(NORD), DLX(NORD)
      INTEGER :: dL, dR, dB, dT
      LOGICAL :: top_on, in_left, in_right, in_bottom, in_top

      nxm1 = NX - 1
      nzm1 = NZ - 1

      Z0 = IE0*DZ0
      IF (SIZE(X) < NX) THEN
         IF (my_rank == 0) WRITE (*, *) 'ERROR Precompute_InversionGrid: SIZE(X) < NX'
         STOP
      END IF
      IF (SIZE(AS) < NORD) THEN
         IF (my_rank == 0) WRITE (*, *) 'ERROR Precompute_InversionGrid: SIZE(AS) < NORD'
         STOP
      END IF
      IF (SIZE(XTO) < NTO .OR. SIZE(ZTO) < NTO) THEN
         IF (my_rank == 0) WRITE (*, *) 'ERROR Precompute_InversionGrid: topography arrays too small for NTO'
         STOP
      END IF
      IF (NX < 2 .OR. NZ < 2 .OR. NORD < 1) THEN
         IF (my_rank == 0) WRITE (*, *) 'ERROR Precompute_InversionGrid: invalid grid/order sizes', NX, NZ, NORD
         STOP
      END IF
      IF (NNZ /= (NZ - 1)*(NORD - 1) + 1) THEN
         IF (my_rank == 0) THEN
            WRITE (*, *) 'ERROR Precompute_InversionGrid: NNZ mismatch. NNZ=', NNZ, &
               ' expected=', (NZ - 1)*(NORD - 1) + 1
         END IF
         STOP
      END IF

      ! Optional consistency check on NBLOCK
      IF (NBLOCK /= nxm1*nzm1) THEN
         IF (my_rank == 0) THEN
            WRITE (*, *) 'ERROR Precompute_InversionGrid: NBLOCK mismatch. NBLOCK=', NBLOCK, &
               ' expected=', nxm1*nzm1
         END IF
         STOP
      END IF

      IF (IE0 < 0 .OR. IE0 > MIN(nxm1, nzm1)) THEN
         IF (my_rank == 0) THEN
            WRITE (*, *) 'ERROR in Precompute_InversionGrid: IE0 out of range. IE0=', IE0, &
               ' NX-1=', nxm1, ' NZ-1=', nzm1
         END IF
         STOP
      END IF

      CALL Allocate_InversionGrid(IG, nxm1, NBLOCK, NORD)

      !-----------------------------------------------------------------
      ! OpenMP-ready structure (parallel over I, serial over J).
      ! Uncomment the !$OMP lines if you want threading here.
      !-----------------------------------------------------------------
      ! !$OMP PARALLEL DO DEFAULT(shared) &
      ! !  PRIVATE(I,J,IM,X1,X2,DX,C1,XI,DM,SS,Z1,Z2,T1,T2,XP_loc,DLX, &
      ! !          dL,dR,dB,dT,top_on,in_left,in_right,in_bottom,in_top, &
      ! !          MID,N0)
      DO I = 1, nxm1

         X1 = X(I)
         X2 = X(I + 1)
         DX = X2 - X1
         C1 = 2.D0/DX

         IG%XBC(I) = 0.5D0*(X1 + X2)

         ! Reset vertical arrays at the start of each column
         DO K = 1, NORD
            Z1(K) = 0.D0
            Z2(K) = 0.D0
            T1(K) = 0.D0
            T2(K) = 0.D0
         END DO

         DO J = 1, nzm1

            ! Analytic block index (safe for parallel over I)
            IM = (I - 1)*nzm1 + J

            ! ---- PML mask: bottom+side always, top if IS0=2 ----
            dL = I - 1
            dR = nxm1 - I
            dB = J - 1
            dT = nzm1 - J

            in_left = (dL < IE0)
            in_right = (dR < IE0)
            in_bottom = (dB < IE0)
            in_top = (dT < IE0)

            top_on = (IS0 == 2)   ! 1: bottom only, 2: top+bottom

            IG%vPML(IM) = (IE0 > 0) .AND. &
                          (in_left .OR. in_right .OR. in_bottom .OR. (top_on .AND. in_top))

            ! ---- vertical geometry (same cases as in old GQG_2D) ----
            IF (J .EQ. 1) THEN

               DO K = 1, NORD
                  Z1(K) = 0.D0
                  Z2(K) = Z1(K) + DZ0
                  T1(K) = 0.D0
                  T2(K) = 0.D0
               END DO

            ELSE IF ((J .GT. 1) .AND. (J .LE. IE0)) THEN

               DO K = 1, NORD
                  Z1(K) = Z2(K)
                  Z2(K) = Z1(K) + DZ0
                  T1(K) = 0.D0
                  T2(K) = 0.D0
               END DO

            ELSE

               DO K = 1, NORD
                  XI = 0.5D0*(X2 - X1)*AS(K) + 0.5D0*(X1 + X2)
                  XP_loc(K) = XI
                  DM = (ZH(NTO, XTO, ZTO, XI) - Z0)/DBLE((NZ - 1) - IE0)
                  Z1(K) = Z2(K)
                  T1(K) = T2(K)
                  Z2(K) = Z1(K) + DM
               END DO

               DO K = 1, NORD
                  XI = XP_loc(K)
                  CALL CDLI(XI, NORD, XP_loc, DLX)
                  SS = 0.D0
                  DO L = 1, NORD
                     SS = SS + DLX(L)*Z2(L)
                  END DO
                  T2(K) = SS
               END DO

            END IF

            MID = INT(DBLE(NORD)/2.D0) + 1
            IG%ZBC(IM) = 0.5D0*(Z1(MID) + Z2(MID))

            N0 = (I - 1)*(NORD - 1)*NNZ + (J - 1)*(NORD - 1) + 1
            IG%N0_BLOCK(IM) = N0
            IG%DX_BLOCK(IM) = DX
            IG%C1_BLOCK(IM) = C1

            IG%Z1_OUT(:, IM) = Z1(:)
            IG%Z2_OUT(:, IM) = Z2(:)
            IG%T1_OUT(:, IM) = T1(:)
            IG%T2_OUT(:, IM) = T2(:)

         END DO
      END DO
      ! !$OMP END PARALLEL DO

   END SUBROUTINE Precompute_InversionGrid

   ! SUBROUTINE Destroy_Canonical_State(canon)
   !   IMPLICIT real(dp) (A-H,O-Z)
   !   TYPE(PhysicalState), INTENT(INOUT) :: canon
   !   IF (ALLOCATED(canon%XMC))  DEALLOCATE(canon%XMC)
   !   IF (ALLOCATED(canon%ZMC))  DEALLOCATE(canon%ZMC)
   !   IF (ALLOCATED(canon%XTOC)) DEALLOCATE(canon%XTOC)
   !   IF (ALLOCATED(canon%ZTOC)) DEALLOCATE(canon%ZTOC)
   !   IF (ALLOCATED(canon%XSRC)) DEALLOCATE(canon%XSRC)
   !   IF (ALLOCATED(canon%ZSRC)) DEALLOCATE(canon%ZSRC)
   !   canon%is_frozen = .FALSE.
   !   canon%XMINC=0D0; canon%XMAXC=0D0; canon%ZMINC=0D0; canon%ZMAXC=0D0
   !   canon%NX_phys=0; canon%NZ_phys=0; canon%NTO_phys=0; canon%NSR_phys=0; canon%IANISO=0
   ! END SUBROUTINE Destroy_Canonical_State

!======================================================================
!  SUBROUTINE FreqBand_Wrapper
!----------------------------------------------------------------------
!  Purpose:
!    Dispatcher for multi-scale: select band IB and call
!    BuildMeshForBand_Plain with SCALAR settings taken from arrays
!    you read from the input file(s).
!
!  Note:
!    - For monoband, pass arrays of length 1 (IB=1).
!    - No automatic DX/DZ anywhere—arrays are from your inputs.
!======================================================================
! SUBROUTINE FreqBand_Wrapper( IB, &
!      FREQB, WLMAXB, WIDTHB, IS0B, NORDB, DXB, DZB, &
!      IANISO, XMINB, XMAXB, ZMINB, ZMAXB, &
!      XM, ZM, XTO, ZTO, CR0, CI0, &
!      XSR, ZSR, NSR, my_rank, &
!      ! sizes (out)
!      IE0, NX, NZ, NNX, NNZ, NN, NPT, NBLOCK, DZ0, Z0, &
!      ! work arrays (inout)
!      X, Z, XP, ZP, AS, WT, XBC, ZBC, MSR, MSR1, FSR, CRR0, CII0, SIG, &
!      DEBUG_OUTPUT )

!   IMPLICIT real(dp) (A - H, O - Z)

!   INTEGER,          INTENT(IN)  :: IB, IANISO, NSR, my_rank
!   real(dp), INTENT(IN)  :: FREQB(:), WLMAXB(:), WIDTHB(:), DXB(:), DZB(:)
!   real(dp), INTENT(IN)  :: XMINB(:), XMAXB(:), ZMINB(:), ZMAXB(:)
!   INTEGER,          INTENT(IN)  :: IS0B(:), NORDB(:)

!   real(dp), INTENT(IN)  :: XM(:), ZM(:), XTO(:), ZTO(:)
!   real(dp), INTENT(IN)  :: CR0(IANISO, :), CI0(IANISO, :)
!   real(dp), INTENT(IN)  :: XSR(:), ZSR(:), YSR(:)

!   INTEGER,          INTENT(OUT) :: IE0, NX, NZ, NNX, NNZ, NN, NPT, NBLOCK
!   real(dp), INTENT(OUT) :: DZ0, Z0

!   real(dp), ALLOCATABLE, INTENT(INOUT) :: X(:), Z(:), XP(:), ZP(:), AS(:), WT(:), XBC(:), ZBC(:)
!   INTEGER,          ALLOCATABLE, INTENT(INOUT) :: MSR(:), MSR1(:,:)
!   real(dp), ALLOCATABLE, INTENT(INOUT) :: FSR(:,:)
!   real(dp), ALLOCATABLE, INTENT(INOUT) :: CRR0(:,:), CII0(:,:)
!   real(dp), ALLOCATABLE, INTENT(INOUT) :: SIG(:)

!   LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

!   CALL BuildMeshForBand_Plain( &
!        FREQB(IB), WLMAXB(IB), WIDTHB(IB), IS0B(IB), NORDB(IB), DXB(IB), DZB(IB), &
!        IANISO, XMINB(IB), XMAXB(IB), ZMINB(IB), ZMAXB(IB), &
!        XM, ZM, XTO, ZTO, CR0, CI0, &
!        XSR,YSR, ZSR, NSR, my_rank, &
!        IE0, NX, NZ, NNX, NNZ, NN, NPT, NBLOCK, DZ0, Z0, &
!        X, Z, XP, ZP, AS, WT, XBC, ZBC, &
!        MSR, MSR1, FSR, CRR0, CII0, SIG, &
!        DEBUG_OUTPUT=DEBUG_OUTPUT )

! END SUBROUTINE FreqBand_Wrapper
! INTEGER :: IB
! IB = 1
! ! arrays of length 1, read from your input file(s)
! FREQB(1)  = FREQN(IFQ)
! WLMAXB(1) = CMAX / FREQB(1)
! WIDTHB(1) = WIDTH
! DXB(1)    = DX_in        ! from file
! DZB(1)    = DZ_in        ! from file
! IS0B(1)   = IS0
! NORDB(1)  = NORD
! XMINB(1)  = XMIN; XMAXB(1) = XMAX
! ZMINB(1)  = ZMIN; ZMAXB(1) = ZMAX

! CALL FreqBand_Wrapper( IB, &
!      FREQB, WLMAXB, WIDTHB, IS0B, NORDB, DXB, DZB, &
!      IANISO, XMINB, XMAXB, ZMINB, ZMAXB, &
!      XM, ZM, XTO, ZTO, CR0, CI0, XSR, ZSR, NSR, my_rank, &
!      IE0, NX, NZ, NNX, NNZ, NN, NPT, NBLOCK, DZ0, Z0, &
!      X, Z, XP, ZP, AS, WT, XBC, ZBC, MSR, MSR1, FSR, CRR0, CII0, SIG, &
!      DEBUG_OUTPUT=.FALSE.)
! DO IB = 1, NBANDS
!    ! FREQB(IB), DXB(IB), DZB(IB), WIDTHB(IB), IS0B(IB), NORDB(IB),
!    ! XMINB/XMAXB/ZMINB/ZMAXB(IB) all come from your input files
!    WLMAXB(IB) = CMAX / FREQB(IB)

!    CALL FreqBand_Wrapper( IB, &
!         FREQB, WLMAXB, WIDTHB, IS0B, NORDB, DXB, DZB, &
!         IANISO, XMINB, XMAXB, ZMINB, ZMAXB, &
!         XM, ZM, XTO, ZTO, CR0, CI0, XSR, ZSR, NSR, my_rank, &
!         IE0, NX, NZ, NNX, NNZ, NN, NPT, NBLOCK, DZ0, Z0, &
!         X, Z, XP, ZP, AS, WT, XBC, ZBC, MSR, MSR1, FSR, CRR0, CII0, SIG, &
!         DEBUG_OUTPUT=.FALSE.)

!    ! Now run your per-band forward/inversion steps on this mesh:
!    ! CALL GF(...), QFRECHET/IND_FRECHET(...), compute_gradient_scaled(...), lbfgs_step_loop(...), etc.
! END DO

END MODULE gridtype_mod
