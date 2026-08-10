module grid_mod
   !----------------------------------------------------------------------
   !  this module contains subroutines for the GQG grid generation
   !  and fill:absorb_zone, GQG_2D, QC +dependencies
   !
   !----------------------------------------------------------------------
   USE omp_lib          ! OpenMP runtime library
   ! USE mkl_service      ! Intel MKL runtime services
   USE hardware_mod     ! Hardware configuration (threads, memory)
   use shared_mod
   ! use GLL_mod
   use GLL_mod, only: GLL
   use output_mod
   use ky_sampling_mod
   USE gridtype_mod
   use iso_fortran_env, ONLY: dp => real64
   IMPLICIT NONE

contains
!======================================================================
!  SUBROUTINE BuildMeshForBand_Plain
!----------------------------------------------------------------------
!  Purpose:
!    Build the mesh for a single band (frequency) using SCALAR inputs
!    you already read from file: DX, DZ, XMIN..ZMAX, WIDTH, etc.
!    No automatic DX/DZ.
!======================================================================
   SUBROUTINE BuildMeshForBand_FromCanon( &
      FREQ, ib, WLMAX, WIDTH, IS0, NORD, DX, DZ, &
      IANISO, ITHOM, IVISCO, PARAM, NPAR, &
      CMIN, CMAX, FREQ2, NSIG, ISO, &
      NSR, NSS, INV, my_rank, NTO, XTO, ZTO, &
      canon, CR0, CI0, &
      IE0, NX, NZ, NZ_ACT, NNX, NNZ, NPT, NBLOCK, DZ0, Z0, &
      X, Z, XP, ZP, AS, WT, &
      XSR, YSR, ZSR, MSR, MSR1, FSR, ICSR, VSR, &
      CRR0, CII0, CRREG, CIREG, CRRU, CIIU, SIG, XMIN, XMAX, ZMIN, ZMAX, IG, &
      msg, CRT, CIT, CR, CI, DEBUG_OUTPUT)

      USE gridtype_mod, ONLY: PhysicalState, InversionGridType
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      ! -------- inputs (scalars) --------
      REAL(dp), INTENT(IN)  :: FREQ, WLMAX, DX, DZ
      REAL(dp), INTENT(IN)  :: CMIN, CMAX, FREQ2, WIDTH
      INTEGER, INTENT(IN)  :: IS0, NORD, IANISO, ITHOM, NPAR, IVISCO
      INTEGER, INTENT(IN)  :: ISO, ib
      INTEGER, INTENT(IN)  :: NSR, NSS, INV, my_rank
      INTEGER, INTENT(IN)  :: NTO
      CHARACTER(LEN=5), INTENT(IN)  :: PARAM(22)

      ! -------- inputs (arrays) --------
      TYPE(PhysicalState), INTENT(IN) :: canon
      REAL(dp), INTENT(IN) :: CR0(:, :), CI0(:, :)
      REAL(dp), OPTIONAL, INTENT(IN) :: CRT(:, :), CIT(:, :)
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      ! -------- outputs (sizes/scalars) --------
      INTEGER, INTENT(OUT) :: IE0, NX, NZ, NNX, NNZ, NPT, NBLOCK, NZ_ACT, NSIG
      REAL(dp), INTENT(OUT) :: DZ0, Z0
      CHARACTER(*), INTENT(OUT) :: msg
      REAL(dp), INTENT(INOUT) :: XMIN, XMAX, ZMIN, ZMAX

      ! -------- inout (allocatable, will (re)allocate here) --------
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: X(:), Z(:), XP(:), ZP(:), AS(:), WT(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: XSR(:), ZSR(:), YSR(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: FSR(:, :)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: MSR(:), MSR1(:, :)
      REAL(dp), ALLOCATABLE, INTENT(INOUT) :: CRR0(:, :), CII0(:, :), CRREG(:, :), CIREG(:, :), CRRU(:, :), CIIU(:, :), SIG(:)

      ! -------- SR mapping helpers --------
      INTEGER, INTENT(INOUT) :: ICSR(:)
      REAL(dp), INTENT(INOUT) :: VSR(:, :, :), XTO(:), ZTO(:)

      ! -------- optional true-model outputs (only if INV=1,0) --------
      REAL(dp), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: CR(:, :), CI(:, :)

      ! -------- locals --------
      LOGICAL :: dbg, have_true_model

      INTEGER :: I
      integer :: NN ! Number of points per subdomain
      REAL(dp), ALLOCATABLE :: XM(:), ZM(:)
      INTEGER, ALLOCATABLE :: MZ(:)

      TYPE(InversionGridType), INTENT(INOUT) :: IG

      ! ---- debug flag ----
      dbg = .FALSE.
      IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      have_true_model = PRESENT(CRT) .AND. PRESENT(CIT) .AND. PRESENT(CR) .AND. PRESENT(CI)

      IF (ALLOCATED(XM)) DEALLOCATE (XM)
      IF (ALLOCATED(ZM)) DEALLOCATE (ZM)
      IF (ALLOCATED(MZ)) DEALLOCATE (MZ)

      ALLOCATE (XM(canon%MXC)); XM = canon%XMC
      ALLOCATE (MZ(canon%MXC)); MZ = canon%MZC
      ALLOCATE (ZM(canon%MZTC)); ZM = canon%ZMC

      XTO = canon%XTOC
      ZTO = canon%ZTOC

      ! ---- source/receiver arrays (as outputs) ----
      IF (ALLOCATED(XSR)) DEALLOCATE (XSR)
      IF (ALLOCATED(ZSR)) DEALLOCATE (ZSR)
      ALLOCATE (XSR(NSR), ZSR(NSR))
      XSR = canon%XSRC(1:NSR)
      ZSR = canon%ZSRC(1:NSR)

      XMIN = canon%XMINC
      XMAX = canon%XMAXC
      ZMIN = canon%ZMINC
      ZMAX = canon%ZMAXC
      NN = 2*NORD*NORD
      !
      IF (my_rank == 0 .AND. dbg) THEN
     WRITE (*, '(A,ES14.6,1X,A,ES14.6,1X,A,ES14.6,1X,A,ES14.6)') 'Original XMIN=', XMIN, 'XMAX=', XMAX, 'ZMIN=', ZMIN, 'ZMAX=', ZMAX
      END IF

      ! ---------------- 1) Absorption & mesh sizes -------------------
      CALL absorb_zone(WIDTH, WLMAX, XMIN, XMAX, ZMIN, ZMAX, DX, DZ, IS0, NORD, &
                       IE0, NX, NZ, NNX, NNZ, NPT, NBLOCK, Z0, DZ0, &
                       NSR, NZ_ACT, my_rank, DEBUG_OUTPUT=dbg)

      ! ---------------- 2) Allocate size-independent arrays ----------
      IF (ALLOCATED(AS)) DEALLOCATE (AS); ALLOCATE (AS(NORD))
      IF (ALLOCATED(WT)) DEALLOCATE (WT); ALLOCATE (WT(NORD))
      IF (ALLOCATED(MSR)) DEALLOCATE (MSR); ALLOCATE (MSR(NSR))
      IF (ALLOCATED(MSR1)) DEALLOCATE (MSR1); ALLOCATE (MSR1(NSR, NN))
      IF (ALLOCATED(FSR)) DEALLOCATE (FSR); ALLOCATE (FSR(NSR, NN))

      ! ---------------- 3) Center / shift coordinates ----------------
      IF (my_rank == 0) WRITE (*, '(A,F6.2)') 'calling ShiftCoordinates for FREQ ', FREQ
      CALL ShiftCoordinates(IS0, NX, canon%MXC, MZ, XM, ZM, IE0, DZ, &
                            XMIN, XMAX, ZMIN, ZMAX, &
                            XTO, ZTO, NTO, XSR, ZSR, NSR, NSS, &
                            my_rank, DEBUG_OUTPUT=dbg)

      ! ---------------- 4) Recompute final shifted mesh sizes --------
      NX = INT((XMAX - XMIN)/DX) + 1
      NZ = INT((ZMAX - ZMIN)/DZ) + 1
      NNX = (NX - 1)*(NORD - 1) + 1
      NNZ = (NZ - 1)*(NORD - 1) + 1
      NPT = NNX*NNZ
      NBLOCK = (NX - 1)*(NZ - 1)
      NZ_ACT = NZ - 1 - IS0*IE0

      IF (my_rank == 0 .AND. dbg) THEN
         WRITE (*, '(A, I8, 1X, I8, A, I8)') 'Shifted mesh NX,NZ = ', NX, NZ, ' NBLOCK = ', NBLOCK
         WRITE (*, '(A, I8, 1X, I8, A, I8)') 'Shifted mesh NNX,NNZ = ', NNX, NNZ, ' NPT = ', NPT
      END IF

      ! ---------------- 5) Allocate size-dependent arrays -----------
      IF (ALLOCATED(X)) DEALLOCATE (X); ALLOCATE (X(NX))
      IF (ALLOCATED(XP)) DEALLOCATE (XP); ALLOCATE (XP(NNX))
      IF (ALLOCATED(ZP)) DEALLOCATE (ZP); ALLOCATE (ZP(NPT))
      IF (ALLOCATED(Z)) DEALLOCATE (Z); ALLOCATE (Z(NBLOCK))

      IF (ALLOCATED(CRR0)) DEALLOCATE (CRR0); ALLOCATE (CRR0(IANISO, NPT))
      IF (ALLOCATED(CII0)) DEALLOCATE (CII0); ALLOCATE (CII0(IANISO, NPT))
      IF (ALLOCATED(CRREG)) DEALLOCATE (CRREG); ALLOCATE (CRREG(IANISO, NPT))
      IF (ALLOCATED(CIREG)) DEALLOCATE (CIREG); ALLOCATE (CIREG(IANISO, NPT))
      IF (ALLOCATED(CRRU)) DEALLOCATE (CRRU); ALLOCATE (CRRU(IANISO, NPT))
      IF (ALLOCATED(CIIU)) DEALLOCATE (CIIU); ALLOCATE (CIIU(IANISO, NPT))
      IF (ALLOCATED(SIG)) DEALLOCATE (SIG); ALLOCATE (SIG(6*NPT))
      X = 0.0_dp
      XP = 0.0_dp
      ZP = 0.0_dp
      Z = 0.0_dp
      CRR0 = 0.0_dp
      CII0 = 0.0_dp
      CRREG = 0.0_dp
      CIREG = 0.0_dp
      CRRU = 0.0_dp
      CIIU = 0.0_dp
      SIG = 0.0_dp

      ! ---------------- 6) GQG + SR mapping + block precompute ------
      IF (my_rank == 0) WRITE (*, '(A,F6.2)') 'calling GQG_2D for FREQ ', FREQ

      CALL GQG_2D(NTO, XTO, ZTO, XMIN, XMAX, ZMIN, ZMAX, DX, DZ, NORD, &
                  NX, NZ, NNX, NNZ, NPT, NBLOCK, &
                  AS, WT, X, XP, ZP, IE0, NSR, NSS, XSR, ZSR, ICSR, VSR, &
                  MSR, MSR1, FSR, my_rank, IS0, IG, DEBUG_OUTPUT=dbg)

      IF (my_rank == 0) THEN
         WRITE (*, *) '--- SHIFT/MAP CHECK (53) ---'
         WRITE (*, '(A,2F14.4)') 'XM first,last = ', XM(1), XM(canon%MXC)
         WRITE (*, '(A)') 'ZM first 10:'
         WRITE (*, '(10(1X,ES14.6))') (ZM(I), I=1, MIN(10, SIZE(ZM)))
         WRITE (*, '(A)') 'XSR,ZSR first 10:'
         DO I = 1, MIN(10, NSR)
            WRITE (*, '(I4,2(1X,F14.4))') I, XSR(I), ZSR(I)
         END DO
         WRITE (*, '(A,I6)') 'MSR(1) = ', MSR(1)
         WRITE (*, '(A)') 'MSR1(1,:) first 10:'
         WRITE (*, '(10(1X,I8))') (MSR1(1, I), I=1, MIN(10, MSR(1)))
         WRITE (*, '(A)') 'FSR(1,:) first 10:'
         WRITE (*, '(10(1X,ES14.6))') (FSR(1, I), I=1, MIN(10, MSR(1)))
      END IF

      ! ---------------- 7) Interpolate base model --------------------
      CALL QC(IANISO, ITHOM, IVISCO, PARAM, NPAR, canon%MXC, XM, MZ, ZM, NNX, NNZ, NPT, XP, ZP, IE0, IS0, NORD, &
              NTO, XTO, ZTO, CR0, CI0, CRR0, CII0, 'STRT', my_rank, &
              canon%XMINC, canon%XMAXC, canon%ZMINC, canon%ZMAXC, ib)

      ! --- 5b) Optional TRUE model interpolation (only if INV=1)
      IF (INV <= 1) THEN
         IF (.NOT. have_true_model) THEN
            IF (my_rank == 0) THEN
               WRITE (*, *) 'BuildMeshForBand_FromCanon: INV=1 requires CRT/CIT/CR/CI to be present.'
            END IF
            STOP 'BuildMeshForBand_FromCanon: missing true-model arrays for INV=1'
         END IF
         IF (have_true_model) THEN
            IF (ALLOCATED(CR)) DEALLOCATE (CR)
            IF (ALLOCATED(CI)) DEALLOCATE (CI)
            ALLOCATE (CR(IANISO, NPT), CI(IANISO, NPT))
            CR = 0.0_dp
            CI = 0.0_dp
            CALL QC(IANISO, ITHOM, IVISCO, PARAM, NPAR, canon%MXC, XM, MZ, ZM, NNX, NNZ, NPT, XP, ZP, IE0, IS0, NORD, &
                    NTO, XTO, ZTO, CRT, CIT, CR, CI, 'TRUE', my_rank, &
                    canon%XMINC, canon%XMAXC, canon%ZMINC, canon%ZMAXC, ib)
         END IF
      END IF

      if (INV == 1 .or. INV == 0) then
         CALL CSIG(IANISO, NPT, CR, DX, DZ, CMIN, CMAX, FREQ2, NSIG, SIG, my_rank)
      else
         CALL CSIG(IANISO, NPT, CRR0, DX, DZ, CMIN, CMAX, FREQ2, NSIG, SIG, my_rank)
      end if

      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, '(A,F6.2,A,3(I8,1X),A,I8)') 'Mesh@', FREQ, ' Hz => NNX,NNZ,NPT= ', &
            NNX, NNZ, NPT, '  NSR=', NSR
      END IF

      msg = 'Mesh built from canonical state'

   END SUBROUTINE BuildMeshForBand_FromCanon

   !-----------------------------------------------------------------------
   !     absorb_zone computes the characteristics of an absorbing boundary
   !     condition based on frequency and wave speed.
   !     It modifies the grid values within a specified zone to reduce
   !     reflections at the boundaries.
   !
   !     Inputs:
   !       WIDTH............. Width of the absorbing zone
   !       CMAX.............. Maximum wave speed
   !       WLMAX............. Maximum wavelength
   !       XMIN, XMAX........ Lateral model extent
   !       ZMIN, ZMAX........ Vertical model extent
   !       DX, DZ............ Grid spacing in x and z directions
   !       IS0............... Absorbing zone type (1 = bottom, 2 = both)
   !       NORD.............. Point per subdomain
   !       NSR............... Number of sources
   !
   !     Outputs:
   !       IE0............... Number of cells in the absorbing zone
   !       NX, NZ............ Grid dimensions in x and z directions
   !       NNX, NNZ.......... Extended grid dimensions
   !       NPT............... Total number of grid points
   !       NBLOCK............ Number of subdomains
   !       Z0................ Depth of the absorbing zone
   !       DZ0............... Adjusted grid spacing in z direction
   !
   !-----------------------------------------------------------------------
   SUBROUTINE absorb_zone(WIDTH, WLMAX, XMIN, XMAX, ZMIN, ZMAX, DX, DZ, IS0, NORD, &
                          IE0, NX, NZ, NNX, NNZ, NPT, NBLOCK, Z0, DZ0, NSR, NZ_ACT, my_rank, DEBUG_OUTPUT)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IS0, NORD, NSR, my_rank
      INTEGER, INTENT(OUT) :: IE0, NX, NZ, NNX, NNZ, NPT, NBLOCK, NZ_ACT
      REAL(dp), INTENT(IN) :: WLMAX, WIDTH, DX, DZ
      REAL(dp), INTENT(INOUT) :: XMIN, XMAX, ZMIN, ZMAX
      REAL(dp), INTENT(OUT) :: Z0, DZ0
      REAL(dp) :: X0
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      ! Local variables
      REAL(dp) :: AIDEA, SLENT, AZONE
      INTEGER :: IS

      ! Calculate ideal absorbing zone
      AIDEA = WIDTH*WLMAX  ! Ideal absorbing zone: 1.5 times the wavelength

      ! Determine maximum domain length (input/physical grid coord)
      SLENT = MAX(XMAX - XMIN, ZMAX - ZMIN)

      ! Calculate actual absorbing zone
      AZONE = MIN(AIDEA, 0.5_dp*SLENT)

      ! Determine number of cells in the absorbing zone
      IE0 = INT(MAX(20.0_dp/REAL(NORD, dp), AZONE/MIN(DX, DZ)))

      ! Absorbing zone adjustments
      IF (my_rank == 0) THEN
         WRITE (*, '(A,1X,F10.3)') 'Original XMIN =', XMIN, ' Original XMAX =', XMAX
         WRITE (*, '(A,1X,F10.3)') 'Original ZMIN =', ZMIN, ' Original ZMAX =', ZMAX
         WRITE (*, *) ' '
      END IF

      DZ0 = DZ
      Z0 = IE0*DZ0
      X0 = IE0*DX
      XMIN = XMIN - X0        ! Extend range for PML
      XMAX = XMAX + X0
      ZMIN = ZMIN - IS0*Z0  ! down if IS0=1, up and down if IS0=2, ZMAX is adjusted in QGQ2D subroutine

      ! Grid dimensions
      NX = INT((XMAX - XMIN)/DX) + 1
      NZ = INT((ZMAX - ZMIN)/DZ) + 1
      NNX = (NX - 1)*(NORD - 1) + 1
      NNZ = (NZ - 1)*(NORD - 1) + 1
      NPT = NNX*NNZ
      NBLOCK = (NX - 1)*(NZ - 1)

      NZ_ACT = NZ - 1 - IS0*IE0 !(FOR  GRADIENT Taper)

      ! IF (PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT .AND. my_rank == 0) then
      IF (my_rank == 0) then
         WRITE (*, *)
         ! WRITE (*, *) 'ABSORB_ZONE OUTPUT'
         ! WRITE (*, '(A, F12.2)') 'Computed AIDEA (ideal absorb. zone for min vle and max freq) = ', AIDEA, ' SLENT (domain length) = ', SLENT
         WRITE (*, '(A, F12.2,A, I4)') 'actual absorb. zone= ', AZONE, 'm  cells= ', IE0
         WRITE (*, '(A,1X,F10.3,2X,A,1X,F10.3)') 'Adjusted XMIN =', XMIN, 'XMAX =', XMAX
         WRITE (*, '(A,1X,F10.3)') 'Adjusted ZMIN =', ZMIN
         WRITE (*, '(A, I8, 1X,I8, 1X,A,1X, I8)') 'Computed NX, NZ  (COARSE grid) = ', NX, NZ, 'NBLOCK (elements) =', NBLOCK
         WRITE (*, '(A, I8, 1X,I8, 1X, A,1X, I8)') 'Computed NNX, NNZ (GQG) = ', NNX, NNZ, 'NPT (grid pts) =', NPT
         WRITE (*, *)
      END IF

      RETURN
   END SUBROUTINE absorb_zone

   SUBROUTINE GQG_2D(NTO, XTO, ZTO, XMIN, XMAX, ZMIN, ZMAX, DX, DZ, NORD, &
                     NX, NZ, NNX, NNZ, NPT, NBLOCK, &
                     AS, WT, X, XP, ZP, IE0, NSR, NSS, XSR, ZSR, ICSR, VSR, &
                     MSR, MSR1, FSR, my_rank, IS0, IG, DEBUG_OUTPUT)
      !----------------------------------------------------------------------
      !
      !     GQG_2D creates a 2D Gaussian quadrature grid for modeling.
      !     It generates grid points, weights, and maps source/receiver
      !     locations to the grid.
      !
      !     Inputs:
      !       NTO............... Number of topography points
      !       XTO(NTO), ZTO(NTO) Topography curve coordinates
      !       XMIN, XMAX........ Lateral model extent
      !       ZMIN, ZMAX........ Vertical model extent
      !       DX, DZ............ Grid spacing in x and z directions
      !       NORD.............. Order of Gaussian quadrature
      !       IE0............... Extensional cell number for PML
      !       NSR............... Number of sources
      !       XSR(NSR), ZSR(NSR) Source coordinates
      !       ICSR(NSR)......... Component indices
      !       VSR(NSR,3,*)...... Source vectors
      !       my_rank........... MPI rank (only rank 0 prints output)
      !       IS0............... Absorbing zone type (1=bottom, 2=top+bottom)
      !       NBLOCK............ Number of blocks ((NX-1)*(NZ-1), from absorb zone)
      !
      !     Outputs:
      !       AS(NORD), WT(NORD) Gaussian points and weights
      !       X(NX)............. Horizontal cell centers
      !       XP(NNX), ZP(NPT).. GQG grid coordinates
      !       MSR(NSR).......... Source mapping counts
      !       MSR1(NSR,*)....... Source mapping indices
      !       FSR(NSR,*)........ Source mapping weights
      !       IG................ InversionGridType (block centers, PML mask, Z1/Z2/T1/T2)
      !
      !-----------------------------------------------------------------------
      USE gridtype_mod, ONLY: InversionGridType, Precompute_InversionGrid

      USE iso_fortran_env, ONLY: dp => real64, output_unit
      IMPLICIT NONE

      TYPE(InversionGridType), INTENT(INOUT) :: IG

      CHARACTER(LEN=80) :: LAB
      INTEGER, INTENT(IN) :: NTO, NORD, IE0, NSR, my_rank, NSS, IS0
      INTEGER, INTENT(INOUT) :: NX, NZ, NNX, NNZ, NPT, NBLOCK

      REAL(dp), INTENT(IN) :: DX, DZ
      REAL(dp), INTENT(IN) :: VSR(:, :, :)
      REAL(dp), INTENT(INOUT) :: XMIN, XMAX, ZMIN, ZMAX

      REAL(dp), INTENT(INOUT) :: AS(:), WT(:)
      INTEGER, INTENT(IN)    :: ICSR(:)
      LOGICAL, OPTIONAL, INTENT(IN) :: DEBUG_OUTPUT

      REAL(dp), INTENT(INOUT) :: XTO(:), ZTO(:), X(:), XP(:), ZP(:)
      REAL(dp), INTENT(INOUT) :: XSR(:), ZSR(:)
      INTEGER, INTENT(INOUT) :: MSR(NSR), MSR1(:, :)
      REAL(dp), INTENT(INOUT) :: FSR(:, :)

      REAL(dp), ALLOCATABLE :: RD(:, :)
      REAL(dp) :: DM, XX, ZZ, RR, S, E, XI, ZI
      INTEGER :: I, J, K, L, II, IP
      INTEGER :: IE, ISR, IS
      INTEGER :: IM, N0, I1, J1, ID
      INTEGER :: unit_grid
      INTEGER :: NPT_local, max_nb
      LOGICAL :: dbg

      dbg = .FALSE.
      IF (PRESENT(DEBUG_OUTPUT)) dbg = DEBUG_OUTPUT
      ! === Generate Gaussian quadrature nodes and weights ===
      CALL GLL(NORD, AS, WT)

      ! === Compute GQG grid size ===
      NX = INT((XMAX - XMIN)/DX) + 1
      NZ = INT((ZMAX - ZMIN)/DZ) + 1
      NNX = (NX - 1)*(NORD - 1) + 1
      NNZ = (NZ - 1)*(NORD - 1) + 1
      NPT = NNX*NNZ
      NBLOCK = (NX - 1)*(NZ - 1)

      IF (SIZE(X) < NX .OR. SIZE(XP) < NNX .OR. SIZE(ZP) < NPT) THEN
         IF (my_rank == 0) THEN
            WRITE (*, *) 'GQG_2D: mesh arrays too small for shifted sizes.'
            WRITE (*, *) '  NX, NZ, NNX, NNZ, NPT = ', NX, NZ, NNX, NNZ, NPT
            WRITE (*, *) '  SIZE(X), SIZE(XP), SIZE(ZP) = ', SIZE(X), SIZE(XP), SIZE(ZP)
         END IF
         STOP 'GQG_2D: mesh arrays not allocated to shifted sizes'
      END IF

      IE = 0
      IF (my_rank == 0 .AND. dbg) THEN
         WRITE (*, *) 'GQG_2D: NX=', NX, ' NZ=', NZ, ' NNX=', NNX, ' NNZ=', NNZ, ' NBLOCK=', NBLOCK
      END IF
      ! === Generate horizontal cell centers ===
      DO I = 1, NX
         X(I) = XMIN + REAL(I - 1, dp)*DX
      END DO

      ! === Build inversion grid geometry + PML index from this mesh ===
      CALL Precompute_InversionGrid(IG, NX, NZ, NNZ, NBLOCK, NORD, IE0, DZ, AS, &
                                    NTO, XTO, ZTO, IS0, X, my_rank)

      ! === Build XP (horizontal GQG abscissae) ===
      DO I = 1, NX - 1
         DO K = 1, NORD
            XI = 0.5_dp*(X(I + 1) - X(I))*AS(K) + 0.5_dp*(X(I + 1) + X(I))
            II = (I - 1)*(NORD - 1) + K
            XP(II) = XI
         END DO
      END DO

      ! === Build ZP using InversionGridType Z1/Z2 ===
      DO I = 1, NX - 1
         DO J = 1, NZ - 1
            IM = (I - 1)*(NZ - 1) + J
            N0 = IG%N0_BLOCK(IM)

            DO I1 = 1, NORD
               DO J1 = 1, NORD
                  ID = N0 + (I1 - 1)*NNZ + (J1 - 1)
                  ZP(ID) = 0.5_dp*(IG%Z2_OUT(I1, IM) - IG%Z1_OUT(I1, IM))*AS(J1) + &
                           0.5_dp*(IG%Z1_OUT(I1, IM) + IG%Z2_OUT(I1, IM))
               END DO
            END DO

         END DO
      END DO

      ! === Save GQG grid to file ===
      IF (my_rank == 0) THEN

         unit_grid = -1
         OPEN (newunit=unit_grid, FILE='m_Grid.dat', STATUS='UNKNOWN')
         LAB = '--------------- GQG_2D: coordinates ----------------'
         WRITE (unit_grid, *) LAB
         IP = 0
         DO I = 1, NNX
            DO J = 1, NNZ
               IP = IP + 1
               WRITE (unit_grid, 29) IP, XP(I), ZP(IP)
               WRITE (70, 29) IP, XP(I), ZP(IP)
29             FORMAT(I8, 2X, 2(F12.4, 4X))
            END DO
         END DO
         CLOSE (unit_grid)
      END IF

      ! === Source-receiver mapping to GQG grid ===
      ALLOCATE (RD(NSR, 2*NORD*NORD))
      DM = 0.5_dp*MIN(DX, DZ)

      ! total GQG grid size
      NPT_local = NNX*NNZ

      ! sanity: MSR1/FSR second dimension must be at least RD second dimension
      max_nb = SIZE(RD, 2)
      IF (SIZE(MSR1, 2) < max_nb .OR. SIZE(FSR, 2) < max_nb) THEN
         IF (my_rank == 0) THEN
            WRITE (*, *) 'GQG_2D: MSR1/FSR second dim too small.'
            WRITE (*, *) '  SIZE(RD,2)=', max_nb, &
               ' SIZE(MSR1,2)=', SIZE(MSR1, 2), &
               ' SIZE(FSR,2)=', SIZE(FSR, 2)
         END IF
         STOP 'GQG_2D: MSR1/FSR arrays too small for RD'
      END IF

      ! initialise mappings defensively
      MSR(1:NSR) = 0
      MSR1(1:NSR, :) = 0
      FSR(1:NSR, :) = 0.0_dp

      DO I = 1, NSR
         XI = XSR(I)
         ZI = ZSR(I)
         IE = 0
         IP = 0

         DO K = 1, NNX
            XX = XI - XP(K)
            DO L = 1, NNZ
               IP = IP + 1
               ZZ = ZI - ZP(IP)
               RR = SQRT(XX*XX + ZZ*ZZ)
               IF (RR .GT. DM) CYCLE

               ! guard against too many neighbours for this source
               IF (IE >= max_nb) THEN
                  IF (my_rank == 0) THEN
                     WRITE (*, *) 'GQG_2D: too many GQG neighbours for source I=', I
                     WRITE (*, *) '  IE would exceed 2*NORD^2 =', max_nb
                     WRITE (*, *) '  DX=', DX, ' DZ=', DZ, ' DM=', DM, ' NORD=', NORD
                  END IF
                  STOP 'GQG_2D: IE exceeds RD/MSR1/FSR bounds'
               END IF

               IE = IE + 1
               RD(I, IE) = RR
               MSR1(I, IE) = IP
            END DO
         END DO

         MSR(I) = IE

         ! extra guard: MSR(I) must be >= 1
         IF (MSR(I) <= 0) THEN
            IF (my_rank == 0) THEN
               WRITE (*, *) 'GQG_2D: source I=', I, ' found zero neighbours; XI,ZI=', XI, ZI
               WRITE (*, *) '  NNX=', NNX, ' NNZ=', NNZ, ' DM=', DM
            END IF
            STOP 'GQG_2D: source with no GQG neighbours'
         END IF

         ! check that all MSR1 indices are in [1, NPT_local]
         IF (MAXVAL(MSR1(I, 1:MSR(I))) > NPT_local .OR. &
             MINVAL(MSR1(I, 1:MSR(I))) < 1) THEN
            IF (my_rank == 0) THEN
               WRITE (*, *) 'GQG_2D: MSR1 out of range for source I=', I
               WRITE (*, *) '  NPT_local=', NPT_local, &
                  ' min(MSR1) =', MINVAL(MSR1(I, 1:MSR(I))), &
                  ' max(MSR1) =', MAXVAL(MSR1(I, 1:MSR(I)))
            END IF
            STOP 'GQG_2D: MSR1 index outside [1, NNX*NNZ]'
         END IF
      END DO

      ! === Normalize interpolation weights for sources/receivers ===
      E = 1.0e-4_dp
      DO I = 1, NSR
         IS = 0
         DO J = 1, MSR(I)
            IF (RD(I, J) .LE. E) IS = J
         END DO

         IF (IS .NE. 0) THEN
            MSR(I) = 1
            MSR1(I, 1) = MSR1(I, IS)
            FSR(I, 1) = 1.0_dp
         ELSE
            S = 0.0_dp
            DO J = 1, MSR(I)
               S = S + 1.0_dp/RD(I, J)
               FSR(I, J) = 1.0_dp/RD(I, J)
            END DO
            IF (S <= 0.0_dp) THEN
               IF (my_rank == 0) THEN
                  WRITE (*, *) 'GQG_2D: zero weight sum for source I=', I
               END IF
               STOP 'GQG_2D: weight sum S <= 0'
            END IF
            DO J = 1, MSR(I)
               FSR(I, J) = FSR(I, J)/S
            END DO
         END IF

         ! final safety: FSR(I,:) should be zero outside 1:MSR(I)
         IF (dbg .AND. my_rank == 0) THEN
            WRITE (*, '("GQG_2D: I=",I6," MSR=",I6," min(MSR1)=",I8," max(MSR1)=",I8)') &
               I, MSR(I), MINVAL(MSR1(I, 1:MSR(I))), MAXVAL(MSR1(I, 1:MSR(I)))
         END IF
      END DO

      DEALLOCATE (RD)

      ! === Output mapping details for each source ===
      IF (my_rank == 0) THEN
         DO I = 1, NSR
            IF (I .LE. NSS) ISR = 1
            WRITE (70, 21) I, XSR(I), ZSR(I)
            WRITE (70, *) ' GQGs No.:'
            WRITE (70, 22) (MSR1(I, J), J=1, MSR(I))
            WRITE (70, *) ' NoteFunc: '
            WRITE (70, 23) (FSR(I, J), J=1, MSR(I))
            WRITE (70, *) ' S_Vector:'
            WRITE (70, 24) ((VSR(I, J, K), K=1, 3), J=1, 3)
         END DO
21       FORMAT('------------ source: ', I4, 2(F12.4, 1X), ' ----------------')
22       FORMAT(3(11X, I12, 1X))
23       FORMAT(11X, E14.7, 1X, E14.7, 1X, E14.7)
24       FORMAT(11X, E14.7, 1X, E14.7, 1X, E14.7, /, &
                11X, E14.7, 1X, E14.7, 1X, E14.7, /, &
                11X, E14.7, 1X, E14.7, 1X, E14.7)
      END IF

      IF (my_rank == 0) WRITE (*, *) 'GQG grid created'
      RETURN
   END SUBROUTINE GQG_2D

   SUBROUTINE SHIFTCoordinates(IS0, NX, MX, MZ, XM, ZM, IE0, DZ, &
                               XMIN, XMAX, ZMIN, ZMAX, &
                               XTO, ZTO, NTO, XSR, ZSR, NSR, NSS, &
                               my_rank, DEBUG_OUTPUT)
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE

      !---- arguments ----
      INTEGER, INTENT(IN)    :: NX, IS0, IE0, NSR, NTO, my_rank, NSS
      INTEGER, INTENT(IN)    :: MX, MZ(:)
      REAL(dp), INTENT(IN)    :: DZ
      REAL(dp), INTENT(INOUT) :: XMIN, XMAX, ZMIN, ZMAX
      REAL(dp), INTENT(INOUT) :: XM(:), ZM(:)
      REAL(dp), INTENT(INOUT) :: XTO(:), ZTO(:)
      REAL(dp), INTENT(INOUT) :: XSR(:), ZSR(:)
      LOGICAL, INTENT(IN)     :: DEBUG_OUTPUT

      !---- locals (declare everything we use) ----
      REAL(dp) :: XC, ZC, Z0
      INTEGER :: I, J, IP, ISR, UNIT_TOP
      CHARACTER(LEN=80) :: LAB

      Z0 = IE0*DZ
      XC = 0.5_dp*(XMIN + XMAX)
      ZC = ZMIN

      ! shift bounds
      XMIN = XMIN - XC
      XMAX = XMAX - XC
      ZMIN = ZMIN - ZC
      ZMAX = ZMAX - ZC

      IF (my_rank == 0) THEN
         WRITE (*, 2) XMIN, XMAX
         WRITE (*, 3) ZMIN, ZMAX
2        FORMAT('       Shifted X_range:', 2(F18.4, 1X))
3        FORMAT('       Shifted Z_range:', 2(F18.4, 1X))
         WRITE (*, *)
      END IF

      ! shift structured grid

      ! shift topography
      DO I = 1, NTO
         XTO(I) = XTO(I) - XC
         IF (IS0 == 1) THEN
            ZTO(I) = ZTO(I) - ZC
         ELSE
            ZTO(I) = ZTO(I) - ZC
         END IF
      END DO

      IF (my_rank == 0) THEN
         LAB = '--------------- topography data shifted --------------------'
         UNIT_TOP = -1
         OPEN (newunit=UNIT_TOP, FILE='m_Top.dat', STATUS='UNKNOWN', ACTION='WRITE')
         DO I = 1, NTO
            WRITE (UNIT_TOP, 4) XTO(I), ZTO(I), XTO(I) + XC, ZTO(I) + ZC
4           FORMAT(4(F12.4, 4X))
         END DO
         CLOSE (UNIT_TOP)
      END IF

      IP = 0
      DO I = 1, MX
         XM(I) = XM(I) - XC
         DO J = 1, MZ(I)
            IP = IP + 1
            IF (IS0 == 1) THEN
               ZM(IP) = ZM(IP) - ZC
            ELSE
               ZM(IP) = ZM(IP) - ZC - Z0
            END IF
         END DO
      END DO

      ! shift sources/receivers

      DO I = 1, NSR
         XSR(I) = XSR(I) - XC
         IF (IS0 == 1) THEN
            ZSR(I) = ZSR(I) - ZC
         ELSE
            ZSR(I) = ZSR(I) - ZC - Z0
         END IF
      END DO
      IF (my_rank == 0) THEN
         OPEN (newunit=UNIT_TOP, FILE='m_SR.dat', STATUS='UNKNOWN', ACTION='WRITE')
         DO I = 1, NSR
            IF (I <= NSS) THEN
               ISR = 1
            ELSE
               ISR = I - NSS + 1
            END IF

            WRITE (UNIT_TOP, 25) ISR, XSR(I), ZSR(I), XSR(I) + XC, ZSR(I) + ZC   ! use (ISR) instead of (I) if you intend remapped indexing
         END DO
25       FORMAT(I12, 3X, 4(F12.3, 3X))
         CLOSE (UNIT_TOP)
      END IF
      IF (DEBUG_OUTPUT .AND. my_rank == 0) THEN
         WRITE (*, '(A,2F12.4)') 'CenterCoordinates: XC,ZC=', XC, ZC
         WRITE (*, '(A,F12.4)') 'CenterCoordinates if top absorption zone: Z0   =', Z0
         WRITE (*, '(A,4F12.4)') 'Bounds after shift: ', XMIN, XMAX, ZMIN, ZMAX
      END IF
   END SUBROUTINE SHIFTCoordinates

   !-----------------------------------------------------------------------
   !
   !     QC interpolates the input model parameters C0(*,*) onto the Gaussian
   !     quadrature grid C(*,*) and transforms Thomsen parameters
   !     into stiffness moduli/Q if needed.
   !
   !     Inputs:
   !       IANISO............ Number of independent elastic moduli + density
   !       ITHOM............. Flag for using Thomsen's parameters (1 = yes, 0 = no)
   !       MODEL(IANISO)..... Names of moduli (e.g., 'rho', 'c11')
   !       MX, XM(MX)........ Number and coordinates of x-points in the model
   !       MZ, ZM(MX*MZ)..... Number and coordinates of z-points in the model
   !       NNX, NNZ.......... Total number of points in x and z directions
   !       XP(NNX), ZP(NNX*NNZ) Coordinates of GQG points
   !       IE0, IS0.......... Extensional blocks and absorbing zone type
   !       NORD.............. Order of Gaussian quadrature
   !       NTO............... Number of topography points
   !       XTO(NTO), ZTO(NTO) Topography curve coordinates
   !       CRE(IANISO,*), CIM(IANISO,*) Input  model parameters
   !          (velocity: Km/s, moduli: (km/s)^2, density: kg/m^3)
   !       TYPE_FLAG......... Type of output (e.g., 'TRUE', 'STRT')
   !       my_rank........... MPI rank (only rank 0 prints output)
   !
   !     Outputs:
   !       CR(IANISO,*), CI(IANISO,*) Interpolated complex-valued moduli
   !       XM(MX), ZM(MZ) modified coordinates for XM and ZM,
   !                     depending on IS0=1 or 2 and topography
   !
   !     Internally Called Subroutines/Functions:
   !       GRID2D_OUT2....... Outputs the grid data to a plotting file
   !
   !-----------------------------------------------------------------------

   SUBROUTINE QC(IANISO, ITHOM, IVISCO, PARAM, NPAR, MX, XM, MZ, ZM, NNX, NNZ, NPT, XP, &
                 ZP, IE0, IS0, NORD, NTO, XTO, ZTO, CRE, CIM, CR, CI, TYPE_FLAG, my_rank, &
                 XMINC, XMAXC, ZMINC, ZMAXC, ib)

      IMPLICIT NONE

      CHARACTER(LEN=2)             :: Prefix
      CHARACTER(LEN=*)             :: PARAM(22)

      CHARACTER(LEN=15)            :: FNAME
      CHARACTER(LEN=*), INTENT(IN) :: TYPE_FLAG

      INTEGER                     :: IM, I, J, K, L, I0, IP
      INTEGER, INTENT(IN) :: IANISO, ITHOM, MX, NORD, IE0, IS0, NTO, NPAR, NNX, NNZ, NPT, my_rank, IVISCO, ib
      INTEGER, INTENT(INOUT) :: MZ(:)
      REAL(dp), INTENT(IN) :: XP(:), ZP(:)
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:)
      REAL(dp), INTENT(IN) :: CRE(:, :), CIM(:, :)
      REAL(dp), INTENT(INOUT) :: XM(:), ZM(:)
      REAL(dp), INTENT(INOUT) :: CR(:, :), CI(:, :)
      REAL(dp), INTENT(IN), OPTIONAL :: XMINC, XMAXC, ZMINC, ZMAXC

      REAL(dp), ALLOCATABLE :: F(:), ZTO1(:)
      REAL(dp) :: ZMAX, XI, XX, X1, X2, ZI, ZZ, RM, RHO, ALP, BET, EPS, DEL, GAM
      REAL(dp) :: C11, C13, C33, C44, C66
      REAL(dp) :: QP, QS, EPSQ, DELQ, GAMQ, Q33, Q44, Q11, Q66, Q13
      REAL(dp) :: A, B, RR
      INTEGER :: I1, I2, J1, J2, ISR, NOT
      INTEGER :: NNX_OUT, NNZ_OUT
      REAL(dp) :: DX_OUT

      !---- four PML points ---
      I1 = (NORD - 1)*IE0 + 1
      I2 = NNX - I1 + 1
      J1 = (NORD - 1)*IE0 + 1
      J2 = NNZ - J1 + 1
      IF (IS0 .EQ. 1) J2 = NNZ
      ZMAX = -1.0e10_dp

      !---- for GQG points ----
      DO I = 1, NNX
         XI = XP(I)

         DO J = 1, NNZ
            IP = (I - 1)*NNZ + J
            ZI = ZP(IP)
            I0 = 1
            IF (ZI .GT. ZMAX) ZMAX = ZI

            !---- xi < xm(1) ---------
            IF (XI .LT. XM(1)) THEN
               XX = XI - XM(1)
               RM = 1.0e10_dp
               DO 2 K = 1, MZ(1)
                  ZZ = ZI - ZM(K)
                  RR = SQRT(XX*XX + ZZ*ZZ)
                  IF (RR .GT. RM) GO TO 2
                  I0 = K
                  RM = RR
2                 CONTINUE
                  GO TO 10
                  END IF

                  !---- xi > xm(mx) --------
                  IF (XI .GT. XM(MX)) THEN
                     IM = 0
                     DO K = 1, MX - 1
                        IM = IM + MZ(K)
                     END DO

                     XX = XI - XM(MX)
                     RM = 1.0e10_dp
                     DO 4 K = 1, MZ(MX)
                        ZZ = ZI - ZM(IM + K)
                        RR = SQRT(XX*XX + ZZ*ZZ)
                        IF (RR .GT. RM) GO TO 4
                        I0 = IM + K
                        RM = RR
4                       CONTINUE
                        GO TO 10
                        END IF

                        !---- x1 <= xi <= x2 -----
                        DO L = 1, MX - 1
                           X1 = XM(L)
                           X2 = XM(L + 1)

                           IF ((XI .GE. X1) .AND. (XI .LE. X2)) THEN
                              IM = 0
                              DO K = 1, L - 1
                                 IM = IM + MZ(K)
                              END DO

                              XX = XI - XM(L)
                              RM = 1.0e10_dp
                              DO 6 K = 1, MZ(L)
                                 ZZ = ZI - ZM(IM + K)
                                 RR = SQRT(XX*XX + ZZ*ZZ)
                                 IF (RR .GT. RM) GO TO 6
                                 RM = RR
                                 I0 = IM + K
6                                CONTINUE

                                 IM = IM + MZ(L)
                                 XX = XI - XM(L + 1)
                                 DO 8 K = 1, MZ(L + 1)
                                    ZZ = ZI - ZM(IM + K)
                                    RR = SQRT(XX*XX + ZZ*ZZ)
                                    IF (RR .GT. RM) GO TO 8
                                    RM = RR
                                    I0 = IM + K
8                                   CONTINUE
                                    GO TO 10
                                    END IF
                                 END DO

10                               DO K = 1, IANISO
                                    CR(K, IP) = CRE(K, I0); 
                                    CI(K, IP) = CIM(K, I0)
                                 END DO

                                 END DO !(J-loop end)
                              END DO !(I-loop end)

                              !---- ouput CR(IANISO,*) & CI(IANISO,*) for first iteration -----
                              ALLOCATE (ZTO1(NTO), F(NPT))

                              ZTO1(1:NTO) = ZTO(1:NTO)
                              IF (IS0 .EQ. 2) ZTO1(1:NTO) = ZMAX

                              IF (my_rank == 0 .AND. ib == 1) THEN
                                 IF (PRESENT(XMINC) .AND. PRESENT(XMAXC) .AND. PRESENT(ZMINC) .AND. PRESENT(ZMAXC)) THEN
                                    NNX_OUT = NTO
                                    IF (NNX_OUT < 2) RETURN
                                    DX_OUT = (XMAXC - XMINC)/REAL(NNX_OUT - 1, dp)
                                    NNZ_OUT = 1 + NINT((ZMAXC - ZMINC)/DX_OUT + 1.0e-10_dp)

                                    IF (NNZ_OUT < 2) NNZ_OUT = 2
                                 ELSE
                                    NNX_OUT = 0
                                    NNZ_OUT = 0
                                 END IF

                                 DO IM = 1, NPAR
                                    IF (TYPE_FLAG == "TRUE") THEN
                                       Prefix = 'mT'
                                    ELSEIF (TYPE_FLAG == 'STRT') THEN
                                       Prefix = 'mI'
                                    ELSE
                                       WRITE (*, *) 'ERROR: Unknown TYPE_FLAG in QC subroutine!'
                                       CYCLE
                                    END IF

                                    IF (IM <= IANISO) THEN
                                       F(1:NPT) = CR(IM, 1:NPT)
                                    ELSE
                                       F(1:NPT) = CI(IM - IANISO + 1, 1:NPT)
                                    END IF

                                    WRITE (FNAME, '(A,A,A,A)') TRIM(Prefix)//'_'//TRIM(PARAM(IM))//'.dat'
                                    IF (NNX_OUT > 0) THEN
                                       CALL GRID2D_OUT_FIXED(FNAME, NNX, NNZ, XP, ZP, F, NTO, XTO, ZTO1, &
                                                             XMINC, XMAXC, ZMINC, ZMAXC, NNX_OUT, NNZ_OUT, IE0, IS0)
                                       !  write(*,*) 'Outputting fixed grid with NNX_OUT=', NNX_OUT, ' NNZ_OUT=', NNZ_OUT
                                    ELSE
                                       CALL GRID2D_OUT(FNAME, NNX, NNZ, XP, ZP, F, NTO, XTO, ZTO1, IE0, IS0)
                                    END IF
                                 END DO
                              END IF

                              DEALLOCATE (F, ZTO1)

                              IF ((IANISO .EQ. 3) .OR. (IANISO .GT. 8)) GO TO 14
                              IF (ITHOM .EQ. 1) THEN
                              DO IP = 1, NPT
                                 RHO = CR(1, IP)
                                 ALP = CR(2, IP)
                                 BET = CR(3, IP)
                                 EPS = CR(4, IP)
                                 DEL = CR(5, IP)
                                 GAM = CR(6, IP)
                                 C11 = RHO*(1.0_dp + 2.0_dp*EPS)*ALP*ALP
                                 C13 = RHO*SQRT(ALP**4*DEL + (ALP*ALP - BET*BET)* &
                                                ((1.0_dp + EPS)*ALP*ALP - BET*BET)) - RHO*BET*BET
                                 C33 = RHO*ALP*ALP
                                 C44 = RHO*BET*BET
                                 C66 = RHO*(1.0_dp + 2.0_dp*GAM)*BET*BET
                                 CR(2, IP) = C11
                                 CR(3, IP) = C13
                                 CR(4, IP) = C33
                                 CR(5, IP) = C44
                                 CR(6, IP) = C66

                                 if (IVISCO == 1) then
                                    QP = CI(2, IP)
                                    QS = CI(3, IP)
                                    EPSQ = CI(4, IP)
                                    DELQ = CI(5, IP)
                                    GAMQ = CI(6, IP)
                                    Q33 = QP
                                    Q44 = QS
                                    Q11 = Q33/(1.0_dp + EPSQ)
                                    Q66 = Q44/(1.0_dp + GAMQ)

                                    A = C33*(C33 - C44)/(2.0_dp*C13*(C13 + C44))
                                    B = C44*(C13 + C33)**2/(2.0_dp*C13*(C13 + C44)*(C33 - C44))
                                    Q13 = 1.0_dp/((1.0_dp + A/DELQ + B)/Q33 - B/Q44)

                                    CI(2, IP) = Q11
                                    CI(3, IP) = Q13
                                    CI(4, IP) = Q33
                                    CI(5, IP) = Q44
                                    CI(6, IP) = Q66
                                 end if

                              END DO
                              END IF

14                            RETURN
                              END subroutine QC

                              !-----------------------------------------------------------------
                              ! Build_Physical_Offsets
                              !
                              ! Purpose:
                              !   Compute 1-based starting offsets for each physical column
                              !   given per-column heights `MZ0` and number of columns `MX0`.
                              !
                              ! Inputs:
                              !   - MZ0(:) : integer array of column heights (number of nodes per column)
                              !   - MX0    : number of columns
                              !
                              ! Outputs:
                              !   - OFF_phys(:) : starting index (1-based) of each column in a
                              !                   flattened storage layout
                              !   - I0_tot      : total number of nodes (last index)
                              !-----------------------------------------------------------------
                              SUBROUTINE Build_Physical_Offsets(MZ0, MX0, OFF_phys, I0_tot)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN)  :: MX0
                                 INTEGER, INTENT(IN)  :: MZ0(:)
                                 INTEGER, INTENT(OUT) :: OFF_phys(:)
                                 INTEGER, INTENT(OUT) :: I0_tot
                                 INTEGER :: k

                                 OFF_phys(1) = 1
                                 DO k = 2, MX0
                                    OFF_phys(k) = OFF_phys(k - 1) + MZ0(k - 1)
                                 END DO
                                 I0_tot = OFF_phys(MX0) + MZ0(MX0) - 1
                              END SUBROUTINE Build_Physical_Offsets

                              !-----------------------------------------------------------------
                              ! Build_Block_NodeDepths
                              !
                              ! Purpose:
                              !   Compute vertical depth positions for each block node in the
                              !   discretized block grid. Near-surface nodes (j <= IE0+1)
                              !   follow uniform spacing `DZ`, while deeper nodes stretch
                              !   to match physical topography.
                              !
                              ! Inputs:
                              !   - X(:)   : horizontal block column centers
                              !   - NX,NZ  : number of block columns and vertical nodes
                              !   - IE0    : number of uniform surface layers
                              !   - DZ     : uniform vertical spacing for the shallow portion
                              !   - NTO,XTO,ZTO : topography table inputs (used by ZH())
                              !
                              ! Output:
                              !   - Znode(:) : per-block-node depth values (size NX*NZ)
                              !-----------------------------------------------------------------
                              SUBROUTINE Build_Block_NodeDepths(X, NX, NZ, IE0, DZ, NTO, XTO, ZTO, Znode)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 REAL(dp), INTENT(IN)  :: X(:), DZ, XTO(:), ZTO(:)
                                 INTEGER, INTENT(IN)  :: NX, NZ, IE0, NTO
                                 REAL(dp), INTENT(OUT) :: Znode(:)

                                 REAL(dp) :: XI, Ztop, Zphys_len
                                 INTEGER :: i, j, ip
                                 REAL(dp) :: denom

                                 DO i = 1, NX
                                    XI = X(i)
                                    Ztop = ZH(NTO, XTO, ZTO, XI)
                                    Zphys_len = MAX(Ztop - IE0*DZ, 0.0_dp)

                                    DO j = 1, NZ
                                       ip = (i - 1)*NZ + j
                                       IF (j <= IE0 + 1) THEN
                                          Znode(ip) = REAL(j - 1, dp)*DZ
                                       ELSE
                                          denom = MAX(REAL(NZ - 1 - IE0, dp), 1.0_dp)
                                          Znode(ip) = IE0*DZ + REAL(j - 1 - IE0, dp)*(Zphys_len/denom)
                                       END IF
                                    END DO
                                 END DO
                              END SUBROUTINE Build_Block_NodeDepths

                              !-----------------------------------------------------------------
                              ! NearestIndexInColumn
                              !
                              ! Purpose:
                              !   Binary-search a sorted column array `ZM0(s:e)` and return the
                              !   index `t` whose value is closest to `z` (ties resolved
                              !   towards the shallower index).
                              !
                              ! Inputs:
                              !   - ZM0(:) : sorted depth values for a physical column
                              !   - s,e    : search interval (inclusive)
                              !   - z      : query depth
                              !
                              ! Output:
                              !   - t : index in [s,e] nearest to `z`
                              !-----------------------------------------------------------------
                              INTEGER FUNCTION NearestIndexInColumn(ZM0, s, e, z) RESULT(t)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 REAL(dp), INTENT(IN) :: ZM0(:), z
                                 INTEGER, INTENT(IN) :: s, e
                                 INTEGER :: lo, hi, mid

                                 lo = s
                                 hi = e
                                 DO WHILE (lo < hi)
                                    mid = (lo + hi)/2
                                    IF (ZM0(mid) <= z) THEN
                                       lo = mid + 1
                                    ELSE
                                       hi = mid
                                    END IF
                                 END DO
                                 IF (lo <= s) THEN
                                    t = s
                                 ELSE IF (lo > e) THEN
                                    t = e
                                 ELSE
                                    IF (ABS(ZM0(lo) - z) < ABS(ZM0(lo - 1) - z)) THEN
                                       t = lo
                                    ELSE
                                       t = lo - 1
                                    END IF
                                 END IF
                              END FUNCTION NearestIndexInColumn

                              !-----------------------------------------------------------------
                              ! Build_Block2Canon_Map
                              !
                              ! Purpose:
                              !   Map each block column center `X(i)` to its corresponding
                              !   canonical column index k (1..MX0) using Voronoi midpoints
                              !   of canonical centers `XM0`.
                              !
                              ! Inputs:
                              !   - X(:)    : block column centers (size NX)
                              !   - XM0(:)  : canonical column centers (size MX0)
                              !   - NX, MX0 : sizes
                              !
                              ! Output:
                              !   - map_i2k(:) : integer mapping of block index i -> canonical k
                              !-----------------------------------------------------------------
                              SUBROUTINE Build_Block2Canon_Map(X, NX, XM0, MX0, map_i2k)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 !---- args ----
                                 REAL(dp), INTENT(IN)  :: X(:), XM0(:)     ! centers
                                 INTEGER, INTENT(IN)  :: NX, MX0
                                 INTEGER, INTENT(OUT) :: map_i2k(:)       ! size NX

                                 !---- locals ----
                                 REAL(dp), ALLOCATABLE :: edge(:)          ! size MX0+1
                                 INTEGER :: i, k, lo, hi, mid

                                 ! guard: trivial case (single canonical column)
                                 IF (MX0 <= 1) THEN
                                    DO i = 1, NX
                                       map_i2k(i) = 1
                                    END DO
                                    RETURN
                                 END IF

                                 ! build Voronoi “borders” from XM0 midpoints
                                 ALLOCATE (edge(MX0 + 1))
                                 edge(1) = XM0(1) - 0.5_dp*(XM0(2) - XM0(1))
                                 DO k = 2, MX0
                                    edge(k) = 0.5_dp*(XM0(k - 1) + XM0(k))
                                 END DO
                                 edge(MX0 + 1) = XM0(MX0) + 0.5_dp*(XM0(MX0) - XM0(MX0 - 1))

                                 ! map each block column center X(i) to canonical cell k with: edge(k) ≤ X(i) < edge(k+1)
                                 DO i = 1, NX
                                    lo = 1
                                    hi = MX0
                                    DO WHILE (lo < hi)
                                       mid = (lo + hi)/2
                                       IF (X(i) >= edge(mid + 1)) THEN
                                          lo = mid + 1
                                       ELSE
                                          hi = mid
                                       END IF
                                    END DO
                                    map_i2k(i) = lo     ! in 1..MX0
                                 END DO

                                 DEALLOCATE (edge)
                              END SUBROUTINE Build_Block2Canon_Map

                              SUBROUTINE GRID_TRANSFORM_BLOCK2NPT(NNX, NNZ, NX, NZ, X, Z, F, NTO, XTO, ZTO, FFF)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN) :: NNX, NNZ, NX, NZ, NTO
                                 REAL(dp), INTENT(IN) :: X(:), Z(:), F(:), XTO(:), ZTO(:)
                                 REAL(dp), INTENT(OUT) :: FFF(:)
                                 REAL(dp), ALLOCATABLE :: A(:), ZM(:), FF(:)
                                 REAL(dp) :: FAIR, EPS, EEE
                                 REAL(dp) :: XA, XB, ZA, ZB, DM, XL, ZL, RL
                                 REAL(dp) :: DX, DZ, XI, ZI, Z0, ZJ, RX, XX, ZZ, RR
                                 REAL(dp) :: RZ, W1, W2, W3, W4
                                 INTEGER :: IP, I, J, II, IPP, IPPP, I01, I02, J01, J02, J03, J04
                                 INTEGER :: NPTI, LI, LFI, LJ, LFJ, JJ, K

                                 FAIR = aver(F, NX*NZ)
                                 EPS = 1.0e-6_dp
                                 EEE = 1.0e-50_dp
                                 FFF(:) = FAIR

                                 XA = 1.0e10_dp
                                 XB = -1.0e10_dp
                                 ZA = 1.0e10_dp
                                 ZB = -1.0e10_dp
                                 DM = 0.0_dp

                                 IP = 0
                                 DO I = 1, NX
                                    XI = X(I)
                                    IF (I .LT. NX) THEN
                                       DX = X(I + 1) - X(I)
                                       IF (DX .GT. DM) DM = DX
                                    END IF
                                    IF (XI .LT. XA) XA = XI
                                    IF (XI .GT. XB) XB = XI
                                    DO J = 1, NZ
                                       IP = IP + 1
                                       ZI = Z(IP)
                                       IF (J .LT. NZ) THEN
                                          DZ = Z(IP + 1) - Z(IP)
                                          IF (DZ .GT. DM) DM = DZ
                                       END IF
                                       IF (ZI .LT. ZA) ZA = ZI
                                       IF (ZI .GT. ZB) ZB = ZI
                                    END DO
                                 END DO

                                 XL = XB - XA
                                 ZL = ZB - ZA
                                 RL = MAX(XL, ZL)

                                 DX = (XB - XA)/REAL(NNX - 2 - 1, dp)
                                 DZ = (ZB - ZA)/REAL(NNZ - 2 - 1, dp)
                                 DM = 0.5_dp*(DX + DZ)
                                 NPTI = (NNX - 2)*(NNZ - 2)
                                 ALLOCATE (FF(NPTI))
                                 ALLOCATE (A(NNZ - 2), ZM(NNZ - 2))

                                 DO J = 1, NNZ - 2
                                    ZM(J) = ZA + REAL(J - 1, dp)*DZ
                                 END DO

                                 IPP = 0
                                 IPPP = 0
                                 DO I = 1, NNX - 2
                                    XI = XA + REAL(I - 1, dp)*DX
                                    Z0 = ZH(NTO, XTO, ZTO, XI)

                                    RX = 1.0e10_dp
                                    DO II = 1, NX
                                       XX = ABS(XI - X(II))
                                       IF (XX .GT. RX) CYCLE
                                       RX = XX
                                       I01 = II
                                    END DO

                                    RX = 1.0e10_dp
                                    DO II = 1, NX
                                       IF (II .EQ. I01) CYCLE
                                       XX = ABS(XI - X(II))
                                       IF (XX .GT. RX) CYCLE
                                       RX = XX
                                       I02 = II
                                    END DO

                                    DO J = 1, NNZ - 2
                                       IPP = (I - 1)*(NNZ - 2) + J
                                       IPPP = (I + 1 - 1)*NNZ + J + 1
                                       ZJ = ZM(J)
                                       A(J) = FAIR
                                       FF(IPP) = A(J)

                                       IF (ZJ .GT. Z0) CYCLE
                                       RZ = 1.0e10_dp
                                       DO JJ = 1, NZ
                                          IP = (I01 - 1)*NZ + JJ
                                          ZZ = ABS(ZJ - Z(IP))
                                          IF (ZZ .GT. RZ) CYCLE
                                          RZ = ZZ
                                          J01 = IP
                                       END DO

                                       RZ = 1.0e10_dp
                                       DO JJ = 1, NZ
                                          IP = (I01 - 1)*NZ + JJ
                                          IF (IP .EQ. J01) CYCLE
                                          ZZ = ABS(ZJ - Z(IP))
                                          IF (ZZ .GT. RZ) CYCLE
                                          RZ = ZZ
                                          J02 = IP
                                       END DO

                                       RZ = 1.0e10_dp
                                       DO JJ = 1, NZ
                                          IP = (I02 - 1)*NZ + JJ
                                          ZZ = ABS(ZJ - Z(IP))
                                          IF (ZZ .GT. RZ) CYCLE
                                          RZ = ZZ
                                          J03 = IP
                                       END DO

                                       RZ = 1.0e10_dp
                                       DO JJ = 1, NZ
                                          IP = (I02 - 1)*NZ + JJ
                                          IF (IP .EQ. J03) CYCLE
                                          ZZ = ABS(ZJ - Z(IP))
                                          IF (ZZ .GT. RZ) CYCLE
                                          RZ = ZZ
                                          J04 = IP
                                       END DO

                                       XX = XI - X(I01)
                                       ZZ = ZJ - Z(J01)
                                       RR = SQRT(XX*XX + ZZ*ZZ)
                                       IF (RR .LE. EPS) THEN
                                          A(J) = F(J01)
                                          FF(IPP) = A(J)
                                          FFF(IPPP) = FF(IPP)
                                          CYCLE
                                       ELSE
                                          W1 = 1.0_dp/RR
                                       END IF

                                       ZZ = ZJ - Z(J02)
                                       RR = SQRT(XX*XX + ZZ*ZZ)
                                       IF (RR .LE. EPS) THEN
                                          A(J) = F(J02)
                                          FF(IPP) = A(J)
                                          FFF(IPPP) = FF(IPP)
                                          CYCLE
                                       ELSE
                                          W2 = 1.0_dp/RR
                                       END IF

                                       XX = XI - X(I02)
                                       ZZ = ZJ - Z(J03)
                                       RR = SQRT(XX*XX + ZZ*ZZ)
                                       IF (RR .LE. EPS) THEN
                                          A(J) = F(J03)
                                          FF(IPP) = A(J)
                                          FFF(IPPP) = FF(IPP)
                                          CYCLE
                                       ELSE
                                          W3 = 1.0_dp/RR
                                       END IF

                                       ZZ = ZJ - Z(J04)
                                       RR = SQRT(XX*XX + ZZ*ZZ)
                                       IF (RR .LE. EPS) THEN
                                          A(J) = F(J04)
                                          FF(IPP) = A(J)
                                          FFF(IPPP) = FF(IPP)
                                          CYCLE
                                       ELSE
                                          W4 = 1.0_dp/RR
                                       END IF

                                       A(J) = (W1*F(J01) + W2*F(J02) + W3*F(J03) + W4*F(J04))/(W1 + W2 + W3 + W4)
                                       IF (ABS(A(J)) .LE. EEE) A(J) = 0.0_dp
                                       FF(IPP) = A(J)
                                       FFF(IPPP) = FF(IPP)
                                    END DO
                                 END DO

                                 J = 0
                                 DO I = 2, NNX - 1
                                    J = J + 1
                                    LI = (I - 1)*NNZ + 1
                                    LFI = (J - 1)*(NNZ - 2) + 1
                                    FFF(LI) = FF(LFI)

                                    LJ = (I - 1)*NNZ + NNZ
                                    LFJ = (J - 1)*(NNZ - 2) + (NNZ - 2)
                                    FFF(LJ) = FF(LFJ)
                                 END DO
                                 DO K = 1, NNZ
                                    FFF(K) = FFF(NNZ + K)
                                 END DO
                                 DO K = 1, NNZ
                                    FFF((NNX - 1)*NNZ + K) = FFF((NNX - 2)*NNZ + K)
                                 END DO

                                 DEALLOCATE (A, ZM, FF)
                                 RETURN
                              END SUBROUTINE GRID_TRANSFORM_BLOCK2NPT

                              SUBROUTINE GRID_TRANSFORM_NPT2BLOCK(NNX, NNZ, NX, NZ, X, Z, FFF, NTO, XTO, ZTO, F)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN) :: NNX, NNZ, NX, NZ, NTO
                                 REAL(dp), INTENT(IN) :: X(:), Z(:), FFF(:), XTO(:), ZTO(:)
                                 REAL(dp), INTENT(OUT) :: F(:)
                                 CHARACTER(LEN=16) :: dbg_env
                                 INTEGER :: dbg_status
                                 LOGICAL :: dbggt
                                 INTEGER :: dbg_count
                                 ! locals
                                 INTEGER :: I, J, ip, ii, jj, idx00, idx10, idx01, idx11
                                 REAL(dp) :: XA, XB, ZA, ZB, XL, ZL, RL, DX, DZ, DM
                                 REAL(dp) :: XI, ZK, Z0, xg, zg, tx, tz, v00, v10, v01, v11, FAIR
                                 INTEGER :: ii_clamp, jj_clamp, nnint, ii0, jj0

                                 REAL(dp) :: tiny
                                 tiny = 1.0e-12_dp

                                 ! ----- derive domain extents from block coordinates (same as forward) -----

                                 XA = +1.0e300_dp; XB = -1.0e300_dp
                                 ZA = +1.0e300_dp; ZB = -1.0e300_dp
                                 DM = 0.0_dp

                                 ip = 0
                                 DO I = 1, NX
                                    XI = X(I)
                                    IF (I < NX) THEN
                                       DX = X(I + 1) - X(I)
                                       IF (DX > DM) DM = DX
                                    END IF
                                    IF (XI < XA) XA = XI
                                    IF (XI > XB) XB = XI
                                    DO J = 1, NZ
                                       ip = ip + 1
                                       ZK = Z(ip)
                                       IF (J < NZ) THEN
                                          DZ = Z(ip + 1) - Z(ip)
                                          IF (DZ > DM) DM = DZ
                                       END IF
                                       IF (ZK < ZA) ZA = ZK
                                       IF (ZK > ZB) ZB = ZK
                                    END DO
                                 END DO

                                 XL = XB - XA
                                 ZL = ZB - ZA
                                 RL = MAX(XL, ZL)

                                 ! Reconstruct the same inner-grid spacing used by GRID_TRANSFORM_BLOCK2NPT:
                                 ! In forward code:
                                 !   DX = (XB - XA) / DBLE(NNX-2-1)   and inner i runs 1..(NNX-2) mapped to fine ii=2..NNX-1
                                 ! => number of inner cells in X is (NNX-2-1)+1 = NNX-2
                                 DX = (XB - XA)/REAL(NNX - 2 - 1, dp)
                                 DZ = (ZB - ZA)/REAL(NNZ - 2 - 1, dp)
                                 DM = 0.5_dp*(DX + DZ)

                                 ! FAIR = mean of interior fine values (ii=2..NNX-1, jj=2..NNZ-1)
                                 FAIR = 0.0D0
                                 nnint = 0
                                 DO ii = 2, NNX - 1
                                    DO jj = 2, NNZ - 1
                                       FAIR = FAIR + FFF((ii - 1)*NNZ + jj)
                                       nnint = nnint + 1
                                    END DO
                                 END DO
                                 IF (nnint > 0) FAIR = FAIR/DBLE(nnint)

                                 ! ----- inverse mapping: for each block node, bilinear from fine grid -----
                                 ip = 0
                                 DO I = 1, NX
                                    XI = X(I)
                                    Z0 = ZH(NTO, XTO, ZTO, XI)

                                    ! find inner fine-grid i-index around XI
                                    ! inner fine-grid node coordinate for ii is: x(ii) = XA + (ii-2)*DX, for ii=2..NNX-1
                                    ! choose ii0 so that x(ii0) <= XI < x(ii0+1)
                                    ii0 = INT((XI - XA)/MAX(DX, tiny)) + 2
                                    IF (ii0 < 2) ii0 = 2
                                    IF (ii0 > NNX - 2) ii0 = NNX - 2

                                    DO J = 1, NZ
                                       ip = ip + 1
                                       ZK = Z(ip)

                                       ! air handling: above topography → FAIR
                                       IF (ZK > Z0) THEN
                                          F(ip) = FAIR
                                          CYCLE
                                       END IF

                                       jj0 = INT((ZK - ZA)/MAX(DZ, tiny)) + 2
                                       IF (jj0 < 2) jj0 = 2
                                       IF (jj0 > NNZ - 2) jj0 = NNZ - 2

                                       ! local coords within the fine cell
                                       xg = XA + DBLE(ii0 - 2)*DX
                                       zg = ZA + DBLE(jj0 - 2)*DZ
                                       tx = (XI - xg)/MAX(DX, tiny)
                                       tz = (ZK - zg)/MAX(DZ, tiny)
                                       IF (tx < 0D0) tx = 0D0
                                       IF (tx > 1D0) tx = 1D0
                                       IF (tz < 0D0) tz = 0D0
                                       IF (tz > 1D0) tz = 1D0

                                       ! neighbors in fine grid (2..NNX-1, 2..NNZ-1)
                                       ii = ii0; jj = jj0
                                       idx00 = (ii - 1)*NNZ + (jj)   ! (ii,   jj)
                                       idx10 = (ii + 1 - 1)*NNZ + (jj)   ! (ii+1, jj)
                                       idx01 = (ii - 1)*NNZ + (jj + 1)   ! (ii,   jj+1)
                                       idx11 = (ii + 1 - 1)*NNZ + (jj + 1)   ! (ii+1, jj+1)

                                       v00 = FFF(idx00)
                                       v10 = FFF(idx10)
                                       v01 = FFF(idx01)
                                       v11 = FFF(idx11)

                                       ! bilinear interpolation
                                       F(ip) = (1D0 - tx)*(1D0 - tz)*v00 + tx*(1D0 - tz)*v10 + (1D0 - tx)*tz*v01 + tx*tz*v11
                                    END DO
                                 END DO
                              END SUBROUTINE GRID_TRANSFORM_NPT2BLOCK

                              SUBROUTINE GRID_TRANSFORM_NPT2BLOCK_ADJ(NNX, NNZ, NX, NZ, X, Z, FFF, NTO, XTO, ZTO, F)

                                 IMPLICIT NONE
                                 INTEGER, INTENT(IN) :: NNX, NNZ, NX, NZ, NTO
                                 REAL(dp), INTENT(IN)  :: X(:), Z(:), FFF(:), XTO(:), ZTO(:)
                                 REAL(dp), INTENT(OUT) :: F(:)

                                 ! locals
                                 INTEGER :: ii, I, J, IP, I01, I02
                                 INTEGER :: K, K1, K2, K3, K4, pcoarse, JF
                                 REAL(dp) :: XA, XB, ZA, ZB, DX, DZ, DM, XI, ZJ, Z0
                                 REAL(dp) :: XL, ZL, RL, FAIR, tiny, xx, zz, rr, w1, w2, w3, w4, val
                                 REAL(dp), ALLOCATABLE :: SUMCOARSE(:), WCOARSE(:)

                                 tiny = 1.0e-12_dp

                                 ! ---- domain extents & FAIR (same as forward) ----
                                 XA = +1.0e300_dp; XB = -1.0e300_dp
                                 ZA = +1.0e300_dp; ZB = -1.0e300_dp
                                 DM = 0.0_dp

                                 IP = 0
                                 DO I = 1, NX
                                    IF (I < NX) THEN
                                       DX = X(I + 1) - X(I)
                                       IF (DX > DM) DM = DX
                                    END IF
                                    IF (X(I) < XA) XA = X(I)
                                    IF (X(I) > XB) XB = X(I)
                                    DO J = 1, NZ
                                       IP = IP + 1
                                       IF (J < NZ) THEN
                                          DZ = Z(IP + 1) - Z(IP)
                                          IF (DZ > DM) DM = DZ
                                       END IF
                                       IF (Z(IP) < ZA) ZA = Z(IP)
                                       IF (Z(IP) > ZB) ZB = Z(IP)
                                    END DO
                                 END DO

                                 XL = XB - XA
                                 ZL = ZB - ZA
                                 RL = MAX(XL, ZL)

                                 ! Same interior spacings used by forward routine:
                                 DX = (XB - XA)/REAL(NNX - 2 - 1, dp)
                                 DZ = (ZB - ZA)/REAL(NNZ - 2 - 1, dp)
                                 DM = 0.5_dp*(DX + DZ)

                                 ! FAIR = mean of interior fine values (optional fallback)
                                 FAIR = 0.0_dp
                                 IP = 0
                                 DO ii = 2, NNX - 1
                                    DO JF = 2, NNZ - 1
                                       IP = (ii - 1)*NNZ + JF
                                       FAIR = FAIR + FFF(IP)
                                    END DO
                                 END DO
                                 IF ((NNX - 2)*(NNZ - 2) > 0) FAIR = FAIR/REAL((NNX - 2)*(NNZ - 2), dp)

                                 ! Accumulators on coarse grid
                                 ALLOCATE (SUMCOARSE(NX*NZ), WCOARSE(NX*NZ))
                                 SUMCOARSE = 0.0_dp
                                 WCOARSE = 0.0_dp

                                 ! ----- adjoint: loop over the SAME interior fine nodes -----
                                 DO ii = 2, NNX - 1
                                    ! forward used I = 1..NNX-2 and mapped ii = I+1, so invert:
                                    I = ii - 1
                                    XI = XA + REAL(I - 1, dp)*DX
                                    Z0 = ZH(NTO, XTO, ZTO, XI)

                                    ! nearest two X coarse nodes: I01, I02
                                    xx = 1.0e300_dp
                                    I01 = 1
                                    DO J = 1, NX
                                       rr = ABS(XI - X(J))
                                       IF (rr < xx) THEN
                                          xx = rr; I01 = J
                                       END IF
                                    END DO
                                    xx = 1.0e300_dp
                                    I02 = I01
                                    DO J = 1, NX
                                       IF (J == I01) CYCLE
                                       rr = ABS(XI - X(J))
                                       IF (rr < xx) THEN
                                          xx = rr; I02 = J
                                       END IF
                                    END DO

                                    DO JF = 2, NNZ - 1
                                       J = JF - 1
                                       ZJ = ZA + REAL(J - 1, dp)*DZ
                                       IF (ZJ > Z0) CYCLE     ! air: contributes nothing back

                                       ! Along column I01: find K1 and K2 (nearest two Z nodes)
                                       zz = 1.0e300_dp
                                       K1 = 1
                                       DO K = 1, NZ
                                          IP = (I01 - 1)*NZ + K
                                          rr = ABS(ZJ - Z(IP))
                                          IF (rr < zz) THEN
                                             zz = rr; K1 = K
                                          END IF
                                       END DO
                                       zz = 1.0e300_dp
                                       K2 = K1
                                       DO K = 1, NZ
                                          IF (K == K1) CYCLE
                                          IP = (I01 - 1)*NZ + K
                                          rr = ABS(ZJ - Z(IP))
                                          IF (rr < zz) THEN
                                             zz = rr; K2 = K
                                          END IF
                                       END DO

                                       ! Along column I02: K3 and K4
                                       zz = 1.0e300_dp
                                       K3 = 1
                                       DO K = 1, NZ
                                          IP = (I02 - 1)*NZ + K
                                          rr = ABS(ZJ - Z(IP))
                                          IF (rr < zz) THEN
                                             zz = rr; K3 = K
                                          END IF
                                       END DO
                                       zz = 1.0e300_dp
                                       K4 = K3
                                       DO K = 1, NZ
                                          IF (K == K3) CYCLE
                                          IP = (I02 - 1)*NZ + K
                                          rr = ABS(ZJ - Z(IP))
                                          IF (rr < zz) THEN
                                             zz = rr; K4 = K
                                          END IF
                                       END DO

                                       ! Same weights as forward: inverse distance to the 4 coarse nodes
                                       xx = XI - X(I01); zz = ZJ - Z((I01 - 1)*NZ + K1); rr = SQRT(xx*xx + zz*zz)
                                       w1 = 1.0_dp/MAX(rr, tiny)

                                       xx = XI - X(I01); zz = ZJ - Z((I01 - 1)*NZ + K2); rr = SQRT(xx*xx + zz*zz)
                                       w2 = 1.0_dp/MAX(rr, tiny)

                                       xx = XI - X(I02); zz = ZJ - Z((I02 - 1)*NZ + K3); rr = SQRT(xx*xx + zz*zz)
                                       w3 = 1.0_dp/MAX(rr, tiny)

                                       xx = XI - X(I02); zz = ZJ - Z((I02 - 1)*NZ + K4); rr = SQRT(xx*xx + zz*zz)
                                       w4 = 1.0_dp/MAX(rr, tiny)

                                       ! normalize weights
                                       val = FFF((ii - 1)*NNZ + JF)
                                       rr = w1 + w2 + w3 + w4
                                       IF (rr <= tiny) CYCLE
                                       w1 = w1/rr; w2 = w2/rr; w3 = w3/rr; w4 = w4/rr

                                       ! scatter-add
                                       pcoarse = (I01 - 1)*NZ + K1
                                       SUMCOARSE(pcoarse) = SUMCOARSE(pcoarse) + w1*val
                                       WCOARSE(pcoarse) = WCOARSE(pcoarse) + w1

                                       pcoarse = (I01 - 1)*NZ + K2
                                       SUMCOARSE(pcoarse) = SUMCOARSE(pcoarse) + w2*val
                                       WCOARSE(pcoarse) = WCOARSE(pcoarse) + w2

                                       pcoarse = (I02 - 1)*NZ + K3
                                       SUMCOARSE(pcoarse) = SUMCOARSE(pcoarse) + w3*val
                                       WCOARSE(pcoarse) = WCOARSE(pcoarse) + w3

                                       pcoarse = (I02 - 1)*NZ + K4
                                       SUMCOARSE(pcoarse) = SUMCOARSE(pcoarse) + w4*val
                                       WCOARSE(pcoarse) = WCOARSE(pcoarse) + w4
                                    END DO
                                 END DO

                                 ! finalize: normalize, fallback to FAIR if no contributions
                                 DO pcoarse = 1, NX*NZ
                                    IF (WCOARSE(pcoarse) > tiny) THEN
                                       F(pcoarse) = SUMCOARSE(pcoarse)/WCOARSE(pcoarse)
                                    ELSE
                                       F(pcoarse) = FAIR
                                    END IF
                                 END DO

                                 DEALLOCATE (SUMCOARSE, WCOARSE)
                              END SUBROUTINE GRID_TRANSFORM_NPT2BLOCK_ADJ

                              SUBROUTINE npt2block_avg_single(NPT, NBLOCK, NORD, NNZ, N0_BLOCK, v_fine, v_block)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE

                                 INTEGER, INTENT(IN)  :: NPT, NBLOCK, NORD, NNZ
                                 INTEGER, INTENT(IN)  :: N0_BLOCK(:)     ! size NBLOCK
                                 REAL(dp), INTENT(IN)  :: v_fine(:)       ! size NPT
                                 REAL(dp), INTENT(OUT) :: v_block(:)      ! size NBLOCK

                                 INTEGER :: NBK, K, L, IP, N0
                                 REAL(dp) :: sum_val, cnt

                                 IF (SIZE(v_fine) < NPT) STOP 'npt2block_avg_single: v_fine too small'
                                 IF (SIZE(v_block) < NBLOCK) STOP 'npt2block_avg_single: v_block too small'

                                 DO NBK = 1, NBLOCK
                                    N0 = N0_BLOCK(NBK)
                                    sum_val = 0.0_dp
                                    cnt = 0.0_dp

                                    DO K = 1, NORD
                                       DO L = 1, NORD
                                          IP = N0 + (K - 1)*NNZ + (L - 1)
                                          IF (IP < 1 .OR. IP > NPT) CYCLE   ! guard
                                          sum_val = sum_val + v_fine(IP)
                                          cnt = cnt + 1.0_dp
                                       END DO
                                    END DO

                                    IF (cnt > 0.0_dp) THEN
                                       v_block(NBK) = sum_val/cnt
                                    ELSE
                                       v_block(NBK) = 0.0_dp
                                    END IF
                                 END DO

                              END SUBROUTINE npt2block_avg_single
                              SUBROUTINE block2npt_const_single(NPT, NBLOCK, NORD, NNZ, N0_BLOCK, v_block, v_fine)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE

                                 INTEGER, INTENT(IN)  :: NPT, NBLOCK, NORD, NNZ
                                 INTEGER, INTENT(IN)  :: N0_BLOCK(:)     ! size NBLOCK
                                 REAL(dp), INTENT(IN)  :: v_block(:)      ! size NBLOCK
                                 REAL(dp), INTENT(OUT) :: v_fine(:)       ! size NPT

                                 INTEGER :: NBK, K, L, IP, N0
                                 REAL(dp) :: val

                                 IF (SIZE(v_block) < NBLOCK) STOP 'block2npt_const_single: v_block too small'
                                 IF (SIZE(v_fine) < NPT) STOP 'block2npt_const_single: v_fine too small'

                                 v_fine(:) = 0.0_dp

                                 DO NBK = 1, NBLOCK
                                    N0 = N0_BLOCK(NBK)
                                    val = v_block(NBK)

                                    DO K = 1, NORD
                                       DO L = 1, NORD
                                          IP = N0 + (K - 1)*NNZ + (L - 1)
                                          IF (IP < 1 .OR. IP > NPT) CYCLE   ! guard
                                          v_fine(IP) = val
                                       END DO
                                    END DO
                                 END DO

                              END SUBROUTINE block2npt_const_single

                              SUBROUTINE map_npt2block_all(NPAR, INVP, NPT, NBLOCK, NORD, NNZ, N0_BLOCK, &
                                                           m_fine, m_coarse_prev)
                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE

                                 INTEGER, INTENT(IN)  :: NPAR, NPT, NBLOCK, NORD, NNZ
                                 INTEGER, INTENT(IN)  :: INVP(:)
                                 INTEGER, INTENT(IN)  :: N0_BLOCK(:)       ! NBLOCK
                                 REAL(dp), INTENT(IN)  :: m_fine(:)    ! actN*NPT
                                 REAL(dp), INTENT(OUT) :: m_coarse_prev(:)     ! actN*NBLOCK

                                 INTEGER :: IM, IA
                                 INTEGER :: ps, pe, cs, ce

                                 IM = 0
                                 DO IA = 1, NPAR
                                    IF (INVP(IA) /= 1) CYCLE

                                    IM = IM + 1
                                    ps = (IM - 1)*NPT + 1
                                    pe = IM*NPT
                                    cs = (IM - 1)*NBLOCK + 1
                                    ce = IM*NBLOCK

                                    CALL npt2block_avg_single(NPT, NBLOCK, NORD, NNZ, N0_BLOCK, &
                                                              m_fine(ps:pe), m_coarse_prev(cs:ce))
                                 END DO

                              END SUBROUTINE map_npt2block_all

                              SUBROUTINE map_block2npt_all(NPAR, INVP, NPT, NBLOCK, NORD, NNZ, N0_BLOCK, &
                                                           p_k, p_kf, &
                                                           NNX, NX, NZ, X, Z, NTO, XTO, ZTO)

                                 USE iso_fortran_env, ONLY: dp => real64
                                 IMPLICIT NONE

                                 INTEGER, INTENT(IN)  :: NPAR, NPT, NBLOCK, NORD, NNZ
                                 INTEGER, INTENT(IN)  :: INVP(:)
                                 INTEGER, INTENT(IN)  :: N0_BLOCK(:)
                                 REAL(dp), INTENT(IN) :: p_k(:)
                                 REAL(dp), INTENT(OUT):: p_kf(:)

                                 ! Optional inputs for GRID_TRANSFORM_BLOCK2NPT
                                 INTEGER, INTENT(IN), OPTIONAL :: NNX, NX, NZ, NTO
                                 REAL(dp), INTENT(IN), OPTIONAL :: X(:), Z(:), XTO(:), ZTO(:)

                                 INTEGER :: IM, IA
                                 INTEGER :: ps, pe, cs, ce
                                 LOGICAL :: use_grid_transform

                                 use_grid_transform = PRESENT(NNX) .AND. PRESENT(NX) .AND. PRESENT(NZ) .AND. &
                                                      PRESENT(X) .AND. PRESENT(Z) .AND. &
                                                      PRESENT(NTO) .AND. PRESENT(XTO) .AND. PRESENT(ZTO)

                                 p_kf(:) = 0.0_dp
                                 IM = 0

                                 DO IA = 1, NPAR
                                    IF (INVP(IA) /= 1) CYCLE

                                    IM = IM + 1

                                    cs = (IM - 1)*NBLOCK + 1
                                    ce = IM*NBLOCK

                                    ps = (IM - 1)*NPT + 1
                                    pe = IM*NPT

                                    IF (use_grid_transform) THEN
                                       CALL GRID_TRANSFORM_BLOCK2NPT(NNX, NNZ, NX, NZ, X, Z, &
                                                                     p_k(cs:ce), NTO, XTO, ZTO, &
                                                                     p_kf(ps:pe))
                                    ELSE
                                       CALL block2npt_const_single(NPT, NBLOCK, NORD, NNZ, N0_BLOCK, &
                                                                   p_k(cs:ce), p_kf(ps:pe))
                                    END IF

                                 END DO

                              END SUBROUTINE map_block2npt_all
                              ! SUBROUTINE map_block2npt_all(NPAR, INVP, NPT, NBLOCK, NORD, NNZ, N0_BLOCK, &
                              !                              p_k, p_kf)
                              !    USE iso_fortran_env, ONLY: dp => real64
                              !    IMPLICIT NONE

                              !    INTEGER, INTENT(IN)  :: NPAR, NPT, NBLOCK, NORD, NNZ
                              !    INTEGER, INTENT(IN)  :: INVP(:)
                              !    INTEGER, INTENT(IN)  :: N0_BLOCK(:)    ! NBLOCK
                              !    REAL(dp), INTENT(IN)  :: p_k(:)         ! actN*NBLOCK
                              !    REAL(dp), INTENT(OUT) :: p_kf(:)        ! actN*NPT

                              !    INTEGER :: IM, IA
                              !    INTEGER :: ps, pe, cs, ce

                              !    p_kf(:) = 0.0_dp
                              !    IM = 0

                              !    DO IA = 1, NPAR
                              !       IF (INVP(IA) /= 1) CYCLE

                              !       IM = IM + 1
                              !       cs = (IM - 1)*NBLOCK + 1
                              !       ce = IM*NBLOCK
                              !       ps = (IM - 1)*NPT + 1
                              !       pe = IM*NPT

                              !       CALL block2npt_const_single(NPT, NBLOCK, NORD, NNZ, N0_BLOCK, &
                              !                                   p_k(cs:ce), p_kf(ps:pe))
                              !    ! CALL GRID_TRANSFORM_BLOCK2NPT(NNX, NNZ, NX, NZ, X, Z, m_coarse_prev(cs:ce), NTO, XTO, ZTO, m_fine(ps:pe))
                              !    END DO

                              ! END SUBROUTINE map_block2npt_all

! !-----------------------------------------------------------------------
! !  GAMMAF computes the Gamma function for specific values of X.
! !
! !  Inputs:
! !    X ............... Input value for which the Gamma function is computed
! !
! !  Output:
! !    GAMMAF .......... Value of the Gamma function at X
! !
! !  Notes:
! !    - Provides fixed known values for standard arguments used by
! !      Gauss–Jacobi routines.
! !    - This is a BLAS-standard compatibility function; logic unchanged.
! !-----------------------------------------------------------------------
!                               FUNCTION GAMMAF(X) RESULT(g)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  REAL(dp), INTENT(IN) :: X
!                                  REAL(dp) :: g

!                                  ! Local constants
!                                  REAL(dp), PARAMETER :: ZERO = 0.0_dp
!                                  REAL(dp), PARAMETER :: HALF = 0.5_dp
!                                  REAL(dp), PARAMETER :: ONE = 1.0_dp
!                                  REAL(dp), PARAMETER :: TWO = 2.0_dp
!                                  REAL(dp), PARAMETER :: FOUR = 4.0_dp
!                                  REAL(dp) :: PI

!                                  PI = FOUR*ATAN(ONE)
!                                  g = ONE

!                                  IF (X == -HALF) g = -TWO*SQRT(PI)
!                                  IF (X == HALF) g = SQRT(PI)
!                                  IF (X == ONE) g = ONE
!                                  IF (X == TWO) g = ONE
!                                  IF (X == 1.5_dp) g = SQRT(PI)/2.0_dp
!                                  IF (X == 2.5_dp) g = 1.5_dp*SQRT(PI)/2.0_dp
!                                  IF (X == 3.5_dp) g = 2.5_dp*1.5_dp*SQRT(PI)/2.0_dp
!                                  IF (X == 3.0_dp) g = 2.0_dp
!                                  IF (X == 4.0_dp) g = 6.0_dp
!                                  IF (X == 5.0_dp) g = 24.0_dp
!                                  IF (X == 6.0_dp) g = 120.0_dp
!                               END FUNCTION GAMMAF

! !-----------------------------------------------------------------------
! !
! !     GLL calculates the Gauss-Lobatto-Legendre points and weights
! !     for integration over the interval [-1, 1].
! !
! !     Inputs:
! !       NP................ Number of points (NP-1 is the polynomial degree)
! !
! !     Outputs:
! !       Z(NP)............. Gauss-Lobatto-Legendre points
! !       W(NP)............. Weights for integration

! !     this subroutine involves the following routines:
! !            ZWGJD(), JACOBF(), ENDW1(), ENDW2(),
! !            GAMMAF(), JACG(), PNORMJ().
! !
! !-----------------------------------------------------------------------
!                               SUBROUTINE GLL(NP, Z, W)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  INTEGER, INTENT(IN)  :: NP
!                                  REAL(dp), INTENT(OUT) :: Z(NP), W(NP)

!                                  ! Locals
!                                  INTEGER  :: N, NM1, I
!                                  REAL(dp) :: ALPHA, BETA, ALPG, BETG
!                                  REAL(dp) :: ONE, TWO
!                                  REAL(dp) :: P, PD, PM1, PDM1, PM2, PDM2

!                                  ONE = 1.0_dp
!                                  TWO = 2.0_dp
!                                  ALPHA = 0.0_dp
!                                  BETA = 0.0_dp

!                                  ! Polynomial degree
!                                  N = NP - 1
!                                  NM1 = N - 1

!                                  ! Basic validation
!                                  IF (NP < 2) ERROR STOP 'GLL: minimum number of Gauss-Lobatto points is 2.'
!                                  IF (ALPHA <= -ONE .OR. BETA <= -ONE) ERROR STOP 'GLL: alpha,beta must be > -1.'

!                                  ! Interior nodes/weights for Jacobi(alpha+1,beta+1)
!                                  IF (NM1 > 0) THEN
!                                     ALPG = ALPHA + ONE
!                                     BETG = BETA + ONE
!                                     CALL ZWGJD(Z(2), W(2), NM1, ALPG, BETG)
!                                  END IF

!                                  ! Endpoints
!                                  Z(1) = -ONE
!                                  Z(NP) = ONE

!                                  ! Convert interior weights to Lobatto weights
!                                  DO I = 2, NP - 1
!                                     W(I) = W(I)/(ONE - Z(I)*Z(I))
!                                  END DO

!                                  ! Endpoint weights via derivative of P_n at ±1
!                                  CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(1))
!                                  W(1) = ENDW1(N, ALPHA, BETA)/(TWO*PD)

!                                  CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(NP))
!                                  W(NP) = ENDW2(N, ALPHA, BETA)/(TWO*PD)
!                               END SUBROUTINE GLL

!                               !-----------------------------------------------------------------------
!                               !
!                               !     ZWGJD calculates the points and weights for the Gauss-Jacobi
!                               !     quadrature over the interval [-1, 1].
!                               !
!                               !     Inputs:
!                               !       NP................ Number of points (NP-1 is the polynomial degree)
!                               !       ALPHA............. Alpha parameter for the Jacobi polynomial
!                               !       BETA.............. Beta parameter for the Jacobi polynomial
!                               !
!                               !     Outputs:
!                               !       Z(NP)............. Gauss-Jacobi points
!                               !       W(NP)............. Weights for integration
!                               !
!                               !-----------------------------------------------------------------------
!                               SUBROUTINE ZWGJD(Z, W, NP, ALPHA, BETA)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  INTEGER, INTENT(IN)  :: NP
!                                  REAL(dp), INTENT(IN)  :: ALPHA, BETA
!                                  REAL(dp), INTENT(OUT) :: Z(NP), W(NP)

!                                  ! Locals
!                                  INTEGER  :: N, NP1, NP2, I
!                                  REAL(dp) :: ONE, TWO
!                                  REAL(dp) :: DN, DNP1, DNP2
!                                  REAL(dp) :: APB, FAC1, FAC2, FAC3
!                                  REAL(dp) :: FNORM, RCOEF
!                                  REAL(dp) :: P, PD, PM1, PDM1, PM2, PDM2

!                                  ONE = 1.0_dp
!                                  TWO = 2.0_dp

!                                  ! Basic validation
!                                  IF (NP < 1) ERROR STOP 'ZWGJD: minimum number of Gauss-Jacobi points is 1.'
!                                  IF (ALPHA <= -ONE .OR. BETA <= -ONE) ERROR STOP 'ZWGJD: alpha,beta must be > -1.'

!                                  N = NP - 1
!                                  DN = REAL(N, dp)
!                                  APB = ALPHA + BETA

!                                  ! Special case: single point
!                                  IF (NP == 1) THEN
!                                     Z(1) = (BETA - ALPHA)/(APB + TWO)
!                                     W(1) = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + ONE)/GAMMAF(APB + TWO)*TWO**(APB + ONE)
!                                     RETURN
!                                  END IF

!                                  ! Compute Jacobi roots (nodes)
!                                  CALL JACG(Z, NP, ALPHA, BETA)

!                                  NP1 = N + 1
!                                  NP2 = N + 2
!                                  DNP1 = REAL(NP1, dp)
!                                  DNP2 = REAL(NP2, dp)

!                                  FAC1 = DNP1 + ALPHA + BETA + ONE
!                                  FAC2 = FAC1 + DNP1
!                                  FAC3 = FAC2 + ONE
!                                  FNORM = PNORMJ(NP1, ALPHA, BETA)

!                                  RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)

!                                  ! Compute weights
!                                  DO I = 1, NP
!                                     CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, NP2, ALPHA, BETA, Z(I))
!                                     W(I) = -RCOEF/(P*PDM1)
!                                  END DO
!                               END SUBROUTINE ZWGJD

! !-----------------------------------------------------------------------
! !  JACOBF computes the values of Jacobi polynomials and their derivatives
! !  at a given point.
! !
! !  Inputs:
! !    N ............... Degree of the Jacobi polynomial
! !    ALP, BET ........ Alpha and Beta parameters for the Jacobi polynomial
! !    X ............... Point at which the polynomial is evaluated
! !
! !  Outputs:
! !    POLY ............ Value of the Jacobi polynomial at X
! !    PDER ............ Derivative of the Jacobi polynomial at X
! !    POLYM1, PDERM1 .. Values and derivatives of the previous polynomial
! !    POLYM2, PDERM2 .. Values and derivatives of the second previous polynomial
! !
! !  Notes:
! !    - Based on the standard three-term recurrence relation for Jacobi
! !      polynomials (see Abramowitz & Stegun, §22.7).
! !    - Handles N = 0 and N = 1 explicitly for numerical stability.
! !-----------------------------------------------------------------------
!                               SUBROUTINE JACOBF(POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2, &
!                                                 N, ALP, BET, X)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE

!                                  ! Arguments
!                                  INTEGER, INTENT(IN)    :: N
!                                  REAL(dp), INTENT(IN)    :: ALP, BET, X
!                                  REAL(dp), INTENT(OUT)   :: POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2

!                                  ! Locals
!                                  REAL(dp) :: APB, POLYL, PDERL, POLYN, PDERN, PSAVE, PDSAVE
!                                  REAL(dp) :: DK, A1, A2, A3, A4, B3
!                                  INTEGER  :: K
!                                  REAL(dp), PARAMETER :: ONE = 1.0_dp, TWO = 2.0_dp

!                                  !--------------------------------------------------------------------
!                                  ! Initialization and base cases
!                                  !--------------------------------------------------------------------
!                                  APB = ALP + BET
!                                  POLY = ONE
!                                  PDER = 0.0_dp

!                                  IF (N == 0) THEN
!                                     POLYM1 = POLY
!                                     PDERM1 = PDER
!                                     POLYM2 = 0.0_dp
!                                     PDERM2 = 0.0_dp
!                                     RETURN
!                                  END IF

!                                  POLYL = POLY
!                                  PDERL = PDER

!                                  ! First-degree polynomial
!                                  POLY = (ALP - BET + (APB + TWO)*X)/TWO
!                                  PDER = (APB + TWO)/TWO

!                                  IF (N == 1) THEN
!                                     POLYM1 = POLYL
!                                     PDERM1 = PDERL
!                                     POLYM2 = 0.0_dp
!                                     PDERM2 = 0.0_dp
!                                     RETURN
!                                  END IF

!                                  !--------------------------------------------------------------------
!                                  ! Recurrence for N >= 2
!                                  !--------------------------------------------------------------------
!                                  DO K = 2, N
!                                     DK = REAL(K, dp)

!                                     A1 = TWO*DK*(DK + APB)*(TWO*DK + APB - TWO)
!                                     A2 = (TWO*DK + APB - ONE)*(ALP**2 - BET**2)
!                                     B3 = (TWO*DK + APB - TWO)
!                                     A3 = B3*(B3 + ONE)*(B3 + TWO)
!                                     A4 = TWO*(DK + ALP - ONE)*(DK + BET - ONE)*(TWO*DK + APB)

!                                     POLYN = ((A2 + A3*X)*POLY - A4*POLYL)/A1
!                                     PDERN = ((A2 + A3*X)*PDER - A4*PDERL + A3*POLY)/A1

!                                     ! Shift previous polynomial states
!                                     PSAVE = POLYL
!                                     PDSAVE = PDERL
!                                     POLYL = POLY
!                                     POLY = POLYN
!                                     PDERL = PDER
!                                     PDER = PDERN
!                                  END DO

!                                  !--------------------------------------------------------------------
!                                  ! Return previous two polynomials and their derivatives
!                                  !--------------------------------------------------------------------
!                                  POLYM1 = POLYL
!                                  PDERM1 = PDERL
!                                  POLYM2 = PSAVE
!                                  PDERM2 = PDSAVE
!                               END SUBROUTINE JACOBF

! !-----------------------------------------------------------------------
! !  ENDW1 computes the weight for the first Gauss–Jacobi quadrature point.
! !
! !  Inputs:
! !    N ............... Degree of the Jacobi polynomial
! !    ALPHA, BETA ..... Alpha and Beta parameters for the Jacobi polynomial
! !
! !  Output:
! !    ENDW1 ........... Weight for the first quadrature point
! !
! !  Notes:
! !    - Uses recurrence formulas for stability with increasing N.
! !    - References: Abramowitz & Stegun, §22.7
! !-----------------------------------------------------------------------
!                               FUNCTION ENDW1(N, ALPHA, BETA) RESULT(w)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  INTEGER, INTENT(IN) :: N
!                                  REAL(dp), INTENT(IN) :: ALPHA, BETA
!                                  REAL(dp)             :: w

!                                  ! Uses module procedure GAMMAF (defined below); no local declaration

!                                  ! Locals
!                                  INTEGER :: I
!                                  REAL(dp) :: ZERO, ONE, TWO, THREE, FOUR
!                                  REAL(dp) :: APB, F1, F2, F3, FINT1, FINT2
!                                  REAL(dp) :: A1, A2, A3, DI, ABN, ABNN

!                                  ZERO = 0.0_dp
!                                  ONE = 1.0_dp
!                                  TWO = 2.0_dp
!                                  THREE = 3.0_dp
!                                  FOUR = 4.0_dp

!                                  APB = ALPHA + BETA

!                                  IF (N == 0) THEN
!                                     w = ZERO
!                                     RETURN
!                                  END IF

!                                  F1 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + ONE)/GAMMAF(APB + THREE)
!                                  F1 = F1*(APB + TWO)*TWO**(APB + TWO)/TWO
!                                  IF (N == 1) THEN
!                                     w = F1
!                                     RETURN
!                                  END IF

!                                  FINT1 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + ONE)/GAMMAF(APB + THREE)
!                                  FINT1 = FINT1*TWO**(APB + TWO)
!                                  FINT2 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + TWO)/GAMMAF(APB + FOUR)
!                                  FINT2 = FINT2*TWO**(APB + THREE)
!                                  F2 = (-TWO*(BETA + TWO)*FINT1 + (APB + FOUR)*FINT2)*(APB + THREE)/FOUR
!                                  IF (N == 2) THEN
!                                     w = F2
!                                     RETURN
!                                  END IF

!                                  DO I = 3, N
!                                     DI = REAL(I - 1, dp)
!                                     ABN = ALPHA + BETA + DI
!                                     ABNN = ABN + DI
!                                     A1 = -TWO*(DI + ALPHA)*(DI + BETA)/(ABN*ABNN*(ABNN + ONE))
!                                     A2 = TWO*(ALPHA - BETA)/(ABNN*(ABNN + TWO))
!                                     A3 = TWO*(ABN + ONE)/((ABNN + TWO)*(ABNN + ONE))
!                                     F3 = -(A2*F2 + A1*F1)/A3
!                                     F1 = F2
!                                     F2 = F3
!                                  END DO

!                                  w = F3
!                               END FUNCTION ENDW1

! !-----------------------------------------------------------------------
! !  ENDW2 computes the weight for the last Gauss–Jacobi quadrature point.
! !
! !  Inputs:
! !    N ............... Degree of the Jacobi polynomial
! !    ALPHA, BETA ..... Alpha and Beta parameters for the Jacobi polynomial
! !
! !  Output:
! !    ENDW2 ........... Weight for the last quadrature point
! !
! !  Notes:
! !    - Same recurrence as ENDW1 but with reversed symmetry (Alpha ↔ Beta).
! !    - References: Abramowitz & Stegun, §22.7
! !-----------------------------------------------------------------------
!                               FUNCTION ENDW2(N, ALPHA, BETA) RESULT(w)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  INTEGER, INTENT(IN) :: N
!                                  REAL(dp), INTENT(IN) :: ALPHA, BETA
!                                  REAL(dp)             :: w

!                                  ! Uses module procedure GAMMAF (defined below); no local declaration

!                                  ! Locals
!                                  INTEGER :: I
!                                  REAL(dp) :: ZERO, ONE, TWO, THREE, FOUR
!                                  REAL(dp) :: APB, F1, F2, F3, FINT1, FINT2
!                                  REAL(dp) :: A1, A2, A3, DI, ABN, ABNN

!                                  ZERO = 0.0_dp
!                                  ONE = 1.0_dp
!                                  TWO = 2.0_dp
!                                  THREE = 3.0_dp
!                                  FOUR = 4.0_dp

!                                  APB = ALPHA + BETA

!                                  IF (N == 0) THEN
!                                     w = ZERO
!                                     RETURN
!                                  END IF

!                                  F1 = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + TWO)/GAMMAF(APB + THREE)
!                                  F1 = F1*(APB + TWO)*TWO**(APB + TWO)/TWO
!                                  IF (N == 1) THEN
!                                     w = F1
!                                     RETURN
!                                  END IF

!                                  FINT1 = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + TWO)/GAMMAF(APB + THREE)
!                                  FINT1 = FINT1*TWO**(APB + TWO)
!                                  FINT2 = GAMMAF(ALPHA + TWO)*GAMMAF(BETA + TWO)/GAMMAF(APB + FOUR)
!                                  FINT2 = FINT2*TWO**(APB + THREE)
!                                  F2 = (TWO*(ALPHA + TWO)*FINT1 - (APB + FOUR)*FINT2)*(APB + THREE)/FOUR
!                                  IF (N == 2) THEN
!                                     w = F2
!                                     RETURN
!                                  END IF

!                                  DO I = 3, N
!                                     DI = REAL(I - 1, dp)
!                                     ABN = ALPHA + BETA + DI
!                                     ABNN = ABN + DI
!                                     A1 = -TWO*(DI + ALPHA)*(DI + BETA)/(ABN*ABNN*(ABNN + ONE))
!                                     A2 = TWO*(ALPHA - BETA)/(ABNN*(ABNN + TWO))
!                                     A3 = TWO*(ABN + ONE)/((ABNN + TWO)*(ABNN + ONE))
!                                     F3 = -(A2*F2 + A1*F1)/A3
!                                     F1 = F2
!                                     F2 = F3
!                                  END DO

!                                  w = F3
!                               END FUNCTION ENDW2

! !-----------------------------------------------------------------------
! !  JACG computes the roots of the Jacobi polynomial, which are
! !  used as the Gauss–Jacobi quadrature points.
! !
! !  Inputs:
! !    NP ............... Number of points (NP-1 is the polynomial degree)
! !    ALPHA, BETA ...... Alpha and Beta parameters for the Jacobi polynomial
! !
! !  Outputs:
! !    Z(NP) ............ Roots of the Jacobi polynomial
! !-----------------------------------------------------------------------
!                               SUBROUTINE JACG(Z, NP, ALPHA, BETA)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  INTEGER, INTENT(IN)  :: NP
!                                  REAL(dp), INTENT(IN)  :: ALPHA, BETA
!                                  REAL(dp), INTENT(OUT) :: Z(NP)

!                                  ! Local variables
!                                  INTEGER  :: KSTOP, N, J, K, JM, I, JMIN
!                                  REAL(dp) :: EPS, DTH, X, X1, X2, XLAST, DELX, RECSUM
!                                  REAL(dp) :: P, PD, PM1, PDM1, PM2, PDM2
!                                  REAL(dp) :: XMIN, SWAP
!                                  REAL(dp), PARAMETER :: PI = 4.0_dp*ATAN(1.0_dp)

!                                  KSTOP = 10
!                                  EPS = 1.0e-12_dp
!                                  N = NP - 1
!                                  DTH = PI/(2.0_dp*REAL(N, dp) + 2.0_dp)

!                                  DO J = 1, NP
!                                     IF (J == 1) THEN
!                                        X = COS((2.0_dp*(REAL(J, dp) - 1.0_dp) + 1.0_dp)*DTH)
!                                     ELSE
!                                        X1 = COS((2.0_dp*(REAL(J, dp) - 1.0_dp) + 1.0_dp)*DTH)
!                                        X2 = XLAST
!                                        X = 0.5_dp*(X1 + X2)
!                                     END IF

!                                     DO K = 1, KSTOP
!                                        CALL JACOBF(P, PD, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, X)
!                                        RECSUM = 0.0_dp
!                                        JM = J - 1
!                                        IF (JM > 0) THEN
!                                           DO I = 1, JM
!                                              RECSUM = RECSUM + 1.0_dp/(X - Z(NP - I + 1))
!                                           END DO
!                                        END IF
!                                        DELX = -P/(PD - RECSUM*P)
!                                        X = X + DELX
!                                        IF (ABS(DELX) < EPS) EXIT
!                                     END DO

!                                     Z(NP - J + 1) = X
!                                     XLAST = X
!                                  END DO

!                                  ! Sort Z in ascending order
!                                  DO I = 1, NP
!                                     XMIN = 2.0_dp
!                                     JMIN = I
!                                     DO J = I, NP
!                                        IF (Z(J) < XMIN) THEN
!                                           XMIN = Z(J)
!                                           JMIN = J
!                                        END IF
!                                     END DO
!                                     IF (JMIN /= I) THEN
!                                        SWAP = Z(I)
!                                        Z(I) = Z(JMIN)
!                                        Z(JMIN) = SWAP
!                                     END IF
!                                  END DO
!                               END SUBROUTINE JACG

! !-----------------------------------------------------------------------
! !  PNORMJ computes the normalization constant for the Jacobi
! !  polynomial of degree N.
! !
! !  Inputs:
! !    N ............... Degree of the Jacobi polynomial
! !    ALPHA, BETA ..... Alpha and Beta parameters for the Jacobi polynomial
! !
! !  Output:
! !    PNORMJ .......... Normalization constant for the Jacobi polynomial
! !-----------------------------------------------------------------------
!                               FUNCTION PNORMJ(N, ALPHA, BETA) RESULT(pnorm)
!                                  USE iso_fortran_env, ONLY: dp => real64
!                                  IMPLICIT NONE
!                                  INTEGER, INTENT(IN) :: N
!                                  REAL(dp), INTENT(IN) :: ALPHA, BETA

!                                  REAL(dp)             :: pnorm
!                                  ! Uses module procedure GAMMAF (defined below); no local declaration
!                                  ! Locals
!                                  INTEGER  :: I
!                                  REAL(dp) :: ONE, TWO, DN, CONST, PROD, DINDX, FRAC

!                                  ONE = 1.0_dp
!                                  TWO = 2.0_dp
!                                  DN = REAL(N, dp)
!                                  CONST = ALPHA + BETA + ONE

!                                  IF (N <= 1) THEN
!                                     PROD = GAMMAF(DN + ALPHA)*GAMMAF(DN + BETA)
!                                     PROD = PROD/(GAMMAF(DN)*GAMMAF(DN + ALPHA + BETA))
!                                     pnorm = PROD*TWO**CONST/(TWO*DN + CONST)
!                                     RETURN
!                                  END IF

!                                  PROD = GAMMAF(ALPHA + ONE)*GAMMAF(BETA + ONE)
!                                  PROD = PROD/(TWO*(ONE + CONST)*GAMMAF(CONST + ONE))
!                                  PROD = PROD*(ONE + ALPHA)*(TWO + ALPHA)
!                                  PROD = PROD*(ONE + BETA)*(TWO + BETA)

!                                  DO I = 3, N
!                                     DINDX = REAL(I, dp)
!                                     FRAC = (DINDX + ALPHA)*(DINDX + BETA)/(DINDX*(DINDX + ALPHA + BETA))
!                                     PROD = PROD*FRAC
!                                  END DO

!                                  pnorm = PROD*TWO**CONST/(TWO*DN + CONST)
!                               END FUNCTION PNORMJ

!======================================================================
!  SUBROUTINE Project_CRR0_2_CR0
!----------------------------------------------------------------------
!  Purpose:
!    Project a full parameter field defined on the current NPT grid
!    directly back to the canonical physical-grid format.
!
!    New logic:
!      1) Treat CRR0_NPT(IA,:) as a field on the old NPT grid
!         of size (NNX,NNZ)
!      2) For each canonical physical point (x,z), interpolate directly
!         from the old NPT grid using the same 4-point IDW logic as
!         GRID2D_OUT_FIXED
!      3) Store the result in packed canonical format CR0(IA,:)
!
!  Entries:
!    1) Old NPT grid geometry:
!       NNX, NNZ, NPT, XP, ZP, NTO, XTO, ZTO, IS0
!    2) Parameter field on old NPT grid:
!       IANISO, CRR0_NPT
!    3) Canonical geometry:
!       canon
!
!  Return:
!    1) CR0(IANISO, I0_tot):
!       canonical physical-grid parameter field
!
!  Notes:
!    - CRR0_NPT is the FULL current model on the old NPT grid.
!    - No block intermediary is used.
!    - This routine follows the interpolation style of GRID2D_OUT_FIXED.
!======================================================================
                              SUBROUTINE Project_CRR0_2_CR0( &
                                 NNX, NNZ, NPT, XP, ZP, &
                                 NTO, XTO, ZTO, IS0, &
                                 IANISO, CRR0_NPT, &
                                 canon, &
                                 CR0)

                                 USE iso_fortran_env, ONLY: dp => real64
                                 USE gridtype_mod, ONLY: PhysicalState
                                 IMPLICIT NONE

                                 INTEGER, INTENT(IN) :: NNX, NNZ, NPT
                                 INTEGER, INTENT(IN) :: NTO, IS0, IANISO
                                 REAL(dp), INTENT(IN) :: XP(:), ZP(:)
                                 REAL(dp), INTENT(IN) :: XTO(:), ZTO(:)
                                 REAL(dp), INTENT(IN) :: CRR0_NPT(:, :)
                                 TYPE(PhysicalState), INTENT(IN) :: canon
                                 REAL(dp), INTENT(INOUT) :: CR0(:, :)

                                 INTEGER :: MX0, I0_tot
                                 INTEGER :: IA, i, j, idx, off
                                 INTEGER :: I01, I02, J01, J02, J03, J04
                                 INTEGER :: II, JJ, IP
                                 REAL(dp) :: xq_phys, zq_phys
                                 REAL(dp) :: xq, zq, zsurf
                                 REAL(dp) :: XA, XB, ZA, ZB
                                 REAL(dp) :: XL_phys, ZL_phys, XL_shift, ZL_shift
                                 REAL(dp) :: excess_x, excess_z
                                 REAL(dp) :: XA_phys_shift, ZA_phys_shift
                                 REAL(dp) :: FAIR, EPS, EEE
                                 REAL(dp) :: RX, RZ, XX, ZZ, RR, W1, W2, W3, W4

                                 INTEGER, ALLOCATABLE :: OFF_phys(:)

                                 EPS = 1.0e-6_dp
                                 EEE = 1.0e-50_dp

                                 MX0 = canon%MXC

                                 ALLOCATE (OFF_phys(MX0))
                                 CALL Build_Physical_Offsets(canon%MZC, MX0, OFF_phys, I0_tot)

                                 IF (SIZE(CR0, 1) /= IANISO) THEN
                                    WRITE (*, *) 'Project_CRR0_2_CR0: inconsistent CR0 first dimension.'
                                    WRITE (*, *) '  SIZE(CR0,1)=', SIZE(CR0, 1), ' IANISO=', IANISO
                                    STOP
                                 END IF

                                 IF (SIZE(CR0, 2) /= I0_tot) THEN
                                    WRITE (*, *) 'Project_CRR0_2_CR0: inconsistent CR0 second dimension.'
                                    WRITE (*, *) '  SIZE(CR0,2)=', SIZE(CR0, 2), ' expected=', I0_tot
                                    STOP
                                 END IF

                                 IF (SIZE(XP) < NNX) THEN
                                    WRITE (*, *) 'Project_CRR0_2_CR0: XP too small for NNX.'
                                    STOP
                                 END IF

                                 IF (SIZE(ZP) < NPT) THEN
                                    WRITE (*, *) 'Project_CRR0_2_CR0: ZP too small for NPT.'
                                    STOP
                                 END IF

                                 IF (SIZE(CRR0_NPT, 1) < IANISO .OR. SIZE(CRR0_NPT, 2) < NPT) THEN
                                    WRITE (*, *) 'Project_CRR0_2_CR0: CRR0_NPT has inconsistent dimensions.'
                                    STOP
                                 END IF

                                 !---------------------------------------------------------------
                                 ! Infer shifted input extents and physical window placement
                                 ! exactly in the same spirit as GRID2D_OUT_FIXED
                                 !---------------------------------------------------------------
                                 XA = MINVAL(XP(1:NNX))
                                 XB = MAXVAL(XP(1:NNX))
                                 ZA = MINVAL(ZP(1:NPT))
                                 ZB = MAXVAL(ZP(1:NPT))

                                 XL_phys = canon%XMAXC - canon%XMINC
                                 ZL_phys = canon%ZMAXC - canon%ZMINC
                                 XL_shift = XB - XA
                                 ZL_shift = ZB - ZA

                                 excess_x = MAX(0.0_dp, XL_shift - XL_phys)
                                 excess_z = MAX(0.0_dp, ZL_shift - ZL_phys)

                                 XA_phys_shift = XA + 0.5_dp*excess_x

                                 IF (IS0 == 1) THEN
                                    ZA_phys_shift = ZB - ZL_phys
                                 ELSEIF (IS0 == 2) THEN
                                    ZA_phys_shift = ZA + 0.5_dp*excess_z
                                 ELSE
                                    ZA_phys_shift = ZA + 0.5_dp*excess_z
                                 END IF

                                 XA_phys_shift = MIN(MAX(XA_phys_shift, XA), XB - XL_phys)
                                 ZA_phys_shift = MIN(MAX(ZA_phys_shift, ZA), ZB - ZL_phys)

                                 !---------------------------------------------------------------
                                 ! Loop over parameters
                                 !---------------------------------------------------------------
                                 DO IA = 1, IANISO

                                    FAIR = aver(CRR0_NPT(IA, 1:NPT), NPT)

                                    off = 0
                                    DO i = 1, MX0

                                       xq_phys = canon%XMC(i)
                                       xq = XA_phys_shift + (xq_phys - canon%XMINC)
                                       zsurf = ZH(NTO, XTO, ZTO, xq)

                                       ! ---- find the two nearest X columns in the old NPT grid ----
                                       RX = 1.0e10_dp
                                       I01 = 1
                                       DO II = 1, NNX
                                          XX = ABS(xq - XP(II))
                                          IF (XX < RX) THEN
                                             RX = XX
                                             I01 = II
                                          END IF
                                       END DO

                                       RX = 1.0e10_dp
                                       I02 = I01
                                       DO II = 1, NNX
                                          IF (II == I01) CYCLE
                                          XX = ABS(xq - XP(II))
                                          IF (XX < RX) THEN
                                             RX = XX
                                             I02 = II
                                          END IF
                                       END DO

                                       DO j = 1, canon%MZC(i)

                                          idx = off + j
                                          zq_phys = canon%ZMC(idx)
                                          zq = ZA_phys_shift + (zq_phys - canon%ZMINC)

                                          CR0(IA, idx) = FAIR
                                          IF (zq > zsurf) CYCLE

                                          ! ---- nearest 2 Z nodes in column I01 ----
                                          RZ = 1.0e10_dp
                                          J01 = (I01 - 1)*NNZ + 1
                                          DO JJ = 1, NNZ
                                             IP = (I01 - 1)*NNZ + JJ
                                             ZZ = ABS(zq - ZP(IP))
                                             IF (ZZ < RZ) THEN
                                                RZ = ZZ
                                                J01 = IP
                                             END IF
                                          END DO

                                          RZ = 1.0e10_dp
                                          J02 = J01
                                          DO JJ = 1, NNZ
                                             IP = (I01 - 1)*NNZ + JJ
                                             IF (IP == J01) CYCLE
                                             ZZ = ABS(zq - ZP(IP))
                                             IF (ZZ < RZ) THEN
                                                RZ = ZZ
                                                J02 = IP
                                             END IF
                                          END DO

                                          ! ---- nearest 2 Z nodes in column I02 ----
                                          RZ = 1.0e10_dp
                                          J03 = (I02 - 1)*NNZ + 1
                                          DO JJ = 1, NNZ
                                             IP = (I02 - 1)*NNZ + JJ
                                             ZZ = ABS(zq - ZP(IP))
                                             IF (ZZ < RZ) THEN
                                                RZ = ZZ
                                                J03 = IP
                                             END IF
                                          END DO

                                          RZ = 1.0e10_dp
                                          J04 = J03
                                          DO JJ = 1, NNZ
                                             IP = (I02 - 1)*NNZ + JJ
                                             IF (IP == J03) CYCLE
                                             ZZ = ABS(zq - ZP(IP))
                                             IF (ZZ < RZ) THEN
                                                RZ = ZZ
                                                J04 = IP
                                             END IF
                                          END DO

                                          ! ---- 4-point IDW interpolation ----
                                          XX = xq - XP(I01)
                                          ZZ = zq - ZP(J01)
                                          RR = SQRT(XX*XX + ZZ*ZZ)

                                          IF (RR <= EPS) THEN
                                             CR0(IA, idx) = CRR0_NPT(IA, J01)
                                          ELSE
                                             W1 = 1.0_dp/RR

                                             ZZ = zq - ZP(J02)
                                             W2 = 1.0_dp/SQRT(XX*XX + ZZ*ZZ)

                                             XX = xq - XP(I02)
                                             ZZ = zq - ZP(J03)
                                             W3 = 1.0_dp/SQRT(XX*XX + ZZ*ZZ)

                                             ZZ = zq - ZP(J04)
                                             W4 = 1.0_dp/SQRT(XX*XX + ZZ*ZZ)

                                             CR0(IA, idx) = ( &
                                                            W1*CRR0_NPT(IA, J01) + &
                                                            W2*CRR0_NPT(IA, J02) + &
                                                            W3*CRR0_NPT(IA, J03) + &
                                                            W4*CRR0_NPT(IA, J04))/(W1 + W2 + W3 + W4)

                                             IF (ABS(CR0(IA, idx)) <= EEE) CR0(IA, idx) = 0.0_dp
                                          END IF

                                       END DO

                                       off = off + canon%MZC(i)
                                    END DO
                                 END DO

                                 DEALLOCATE (OFF_phys)

                              END SUBROUTINE Project_CRR0_2_CR0
                              END module grid_mod
