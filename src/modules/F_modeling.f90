module F_modeling

   USE omp_lib          ! OpenMP runtime library
   USE mkl_service      ! Intel MKL runtime services
   USE MPI              ! MPI (message passing)
   USE hardware_mod     ! Hardware configuration (threads, memory)
   USE shared_mod       ! Shared variables/constants across program
   USE boundaries_mod   ! Absorbing boundary condition logic
   USE solvers_mod
   USE output_mod       ! Output file and logging utilities
   USE gridtype_mod
   USE err_mpi_mod     ! MPI error handling routines
   USE constant_mod   ! Physical and numerical constants
   USE stiffness_assembly_mod ! Stiffness matrix assembly routines
   USE scalers_mod   ! Scaling utilities
   use iso_fortran_env, only: dp => real64, sp => real32
   IMPLICIT NONE

contains
   !----------------------------------------------------------------------C
   !                                                                      C
   !     calculate the Green's function vector in the Frequency domain    C
   !                                                                      C
   !     Entries:                                                         C
   !      (0) I25D =0 OR =1 .....................for 2D & 2.5D modelling; C
   !      (1) XTO(NTO),ZTO(NTO).......................surface topography; C
   !      (2) X(NX),NZ.............................subdomain coordinates; C
   !      (4) CR(IANISO,*),CI(IANISO,*)......Real & imaginary Rho & cij; C
   !      (5) AS(NORD),WT(NORD)..........GQG points & weights in (Dx,Dz); C
   !      (6) FREQ,FKY(NK)....................Frequency & Ky wavenumbers; C
   !      (7) NSS.....sources for Green's tensor (for inversion NSS=NSR); C
   !      (8) XRS(NSR),ZSR(NSR)........coordinates of source & geophones; C
   !      (9) MSR(NSR),MSR1(NSR,*),FSR(NSR,*).......S & G's interpolants. C
   !     (10) IE0,DZ0..........extension-block number & its size for PML; C
   !     (11) IS0=1 or 2.........No absorbing PML at free-surface or has; C
   !                                                                      C
   !     Returns:                                                         C
   !                                                                      C
   !      (1) GFX(NSS,3,NK,*),GFY(NSS,3,NK,*),GFZ(NSS,3,NK,*)..GF tensor. C
   !                                                                      C
   !----------------------------------------------------------------------C
!----------------------------------------------------------------------C

   SUBROUTINE GF(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                 AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
                 IE0, IS0, DZ0, GFX, GFY, GFZ, NBLOCK, IG, &
                 my_rank, n_process, comm, ITER, SOLVER_KIND, DEBUG_OUTPUT)
      IMPLICIT NONE
      ! ----------------- inputs -----------------
      INTEGER, INTENT(IN) :: I25D, IFQ, NTO, NX, NZ, NNX, NNZ, NPT, IANISO, NORD
      INTEGER, INTENT(IN) :: NK, NSR, NSS, IE0, IS0, ITER, NBLOCK
      INTEGER, INTENT(IN) :: my_rank, n_process, comm
      INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
      REAL(dp), INTENT(IN) :: FREQ, DZ0, XTO(:), ZTO(:), X(:), AS(:), WT(:), FKY(:)
      REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
      COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      INTEGER, INTENT(IN) :: SOLVER_KIND   ! 1=LU(auto), 2=PARDISO, 3=MUMPS
      LOGICAL, OPTIONAL :: DEBUG_OUTPUT
      TYPE(InversionGridType), INTENT(IN) :: IG
      ! ----------------- locals -----------------
      INTEGER :: KL, KU, LDAB, N
      INTEGER :: start_index, end_index, modeling_iterations_per_process, remainder, my_tasks
      INTEGER :: NCORE, NCORE_MAX, NCOREI, rate, t0_in, t1_in, t0_all, t1_all
      REAL(dp) :: t_internal_local, t_internal_max
      INTEGER :: ierred_local, ierr, cmp
      REAL(dp) :: t_global_local, t_global_max, sum_check

      REAL(dp) :: AMY, GFMY, SMY, TOL,RHS_MY
      INTEGER :: IK, IS, IC, II
      INTEGER(KIND=8) :: d2_64, chunk_64
      LOGICAL :: dbg


INTEGER  :: n_dof, nrhs_tot, nbatches,MAX_RHS_BATCH

      ! --- simple solver codes
      INTEGER, PARAMETER :: SOLVER_LU = 1, SOLVER_PARDISO = 2, SOLVER_MUMPS = 3
      dbg = .FALSE.
      if (PRESENT(DEBUG_OUTPUT)) dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      d2_64 = 3_8
      ! Target ~128 MiB per chunk (COMPLEX(dp) = 16 bytes)
      ! chunk_64 = (128_8*1024_8*1024_8)/16_8
      chunk_64 = (128_8*1024_8*1024_8)/8_8
      ! --- timers: get tick rate and a global start count

      CALL SYSTEM_CLOCK(COUNT_RATE=rate)
      CALL SYSTEM_CLOCK(t0_all)

      ! zero outputs
      GFX = (0.0_dp, 0.0_dp); GFY = (0.0_dp, 0.0_dp); GFZ = (0.0_dp, 0.0_dp)

      ! sizes for band (LU path needs these)
      N = 3*NPT
      KU = 3*(NORD - 1)*(NNZ + 1) + 3
      KL = KU
      LDAB = 2*KL + KU + 1

      ! ---- split IK among ranks
      CALL split_work(NK, my_rank, n_process, start_index, end_index, my_tasks)

      ! ================== run solver (timed) ==================
      CALL MPI_Barrier(comm, ierr)
      CALL SYSTEM_CLOCK(t0_in)


!---------------- Pre-solver memory preview (approx, GB) ----------------

n_dof    = N                 ! should be 3*NPT
nrhs_tot = 3*NSS   !        ! total RHS to solve for (all sources, all components)
MAX_RHS_BATCH=32  ! heuristic max RHS per batch (solver-agnostic, for memory preview)
nbatches = (nrhs_tot + MAX_RHS_BATCH - 1) / MAX_RHS_BATCH   ! ceil(nrhs_tot / batch)

! Stiffness matrix A (banded complex, approx 2*LDAB rows)
AMY = 2.0_dp * BYTES_CPLX16 * REAL(LDAB,dp) * REAL(n_dof,dp) / BYTES_GB

! Green’s functions: GFX/GFY/GFZ each (NSS,3,NK,NPT)
GFMY = 3.0_dp * BYTES_CPLX16/2 * REAL(NSS,dp) * 3.0_dp * REAL(NK,dp) * REAL(NPT,dp) / BYTES_GB

! Batched RHS + SOL buffers (solver-agnostic preview)
! (Two arrays: rhs + sol) each sized (N, MAX_RHS_BATCH)
RHS_MY = 2.0_dp * BYTES_CPLX16 * REAL(n_dof,dp) * REAL(MAX_RHS_BATCH,dp) / BYTES_GB

! Total preview (keep TOL name)
TOL = AMY + GFMY + RHS_MY

IF (my_rank == 0 .AND. dbg.and.ITER<=1) THEN

   WRITE(*,'(A)') '--- Pre-solver preview (memory + RHS workload) ---'
   WRITE(*,'(A,I10,A,I10,A,I10)') '   GQG dims: NNX=',NNX,'  NNZ=',NNZ,'  NPT=',NPT
   WRITE(*,'(A,I10)')            '   DOF:     N=3*NPT =', n_dof
   WRITE(*,'(A,I10,A,I10,A,I10)')'   Sources: NSR=',NSR,'  NSS=',NSS,'  NK=',NK
   WRITE(*,'(A,I10)')            '   RHS total (per IK): nrhs_tot=3*NSS =', nrhs_tot
   WRITE(*,'(A,I10,A,I10)')      '   Batch:    MAX_RHS_BATCH=', MAX_RHS_BATCH, '  nbatches=', nbatches
   ! WRITE(*,'(A,F10.3,A)') '   AMY   (A banded):        ', AMY,    ' GB'
   WRITE(*,'(A,F10.3,A)') '   GFMY  (GFX+GFY+GFZ):     ', GFMY,   ' GB'
   WRITE(*,'(A,F10.3,A)') '   RHS_MY(rhs+sol buffers): ', RHS_MY, ' GB'
   ! WRITE(*,'(A,F10.3,A)') '   TOL   (A+GF+buffers):    ', TOL,    ' GB'
   WRITE(*,*)
     call flush(6)
   !  call log_rss("1st inside GF", 6, my_rank)
 
END IF


      SELECT CASE (SOLVER_KIND)

      CASE (SOLVER_LU)
         ! --- (optional) memory guide used by LU auto-choice

         ! NCORE: max ky-chunks for LU (from memory model)
         NCORE_MAX = FLOOR((MEM_MAX - TOL)/AMY) - 1
         IF (NCORE_MAX < 1) NCORE_MAX = 1
         NCORE = MIN(NCORE_MAX, NTHREAD_MAX, my_tasks)

         ! NCOREI: max source-chunks (from memory model)
         NCORE_MAX = FLOOR((MEM_MAX - TOL)/SMY) - 1
         IF (NCORE_MAX < 1) NCORE_MAX = 1
         NCOREI = MIN(NCORE_MAX, NSS, NTHREAD_MAX)

         IF ((I25D == 0) .OR. (NCORE <= 0)) NCORE = 1
         IF ((I25D == 0) .OR. (NCOREI <= 0)) NCOREI = 1

         ! --- simplified choice: hybrid vs per-source, no special NSS<=1 case ---
         IF (NCORE > 2*n_process .OR. NCORE > NSS) THEN
         CALL run_LU_hybrid(N, KL, KU, LDAB, start_index, end_index, &
                            I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                            AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                            IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCORE)
         ELSE
            CALL run_LU_mkl_NSS(N, KL, KU, LDAB, start_index, end_index, &
                                I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                                AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                                IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCOREI)
         END IF

      CASE (SOLVER_PARDISO)
         CALL run_ik_loop_pardiso(N, KL, KU, LDAB, start_index, end_index, &
                                  I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                                  AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                                  IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, NBLOCK, &
                                  comm, ITER)

      CASE (SOLVER_MUMPS)
#if HAVE_MUMPS
         CALL run_ik_loop_mumps(N,LDAB, KL, KU, start_index, end_index, &
                                I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                                AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                                IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, NBLOCK,comm, ITER)
#else
         IF (my_rank == 0) WRITE (*, *) &
            'ERROR: solver 3 (MUMPS) was requested, but this executable was built with USE_MUMPS=0.'
         CALL MPI_Abort(comm, 3, ierr)
#endif

      CASE DEFAULT
         IF (my_rank == 0) WRITE (*, *) 'GF_production: unknown SOLVER_KIND =', SOLVER_KIND
         CALL MPI_Abort(comm, 1, ierr)
      END SELECT

      CALL SYSTEM_CLOCK(t1_in)

      ! --- safe elapsed-time computation + reduction
      t_internal_local = 0.0_dp
      IF (rate > 0) t_internal_local = REAL(t1_in - t0_in, dp)/REAL(rate, dp)

      CALL MPI_Barrier(comm, ierr)
   !    IF (dbg) THEN
   !          sum_check = 0.0D0
   !          DO IK = 1, NK
   !             DO IS = 1, NSS
   !                DO IC = 1, 3
   !                   DO II = 1, NPT
   !                      sum_check = sum_check + ABS(GFZ(IS, IC, IK, II))
   !                   END DO
   !                END DO
   !             END DO
   !          END DO
   !          PRINT *, 'Rank', my_rank, 'local GFZ sum before =', sum_check
   !   WRITE (*, '("Chunked_AllReduce_4D args: rank=",I0," comm=",I0," NSS=",I0," d2_64=",I0," NK=",I0," NPT=",I0," chunk_64=",I0)') &
   !          my_rank, comm, NSS, d2_64, NK, NPT, chunk_64
   !    END IF

      CALL Chunked_AllReduce_4D(GFX, NSS, d2_64, NK, NPT, chunk_64, my_rank, comm, ierred_local)
      CALL Chunked_AllReduce_4D(GFY, NSS, d2_64, NK, NPT, chunk_64, my_rank, comm, ierred_local)
      CALL Chunked_AllReduce_4D(GFZ, NSS, d2_64, NK, NPT, chunk_64, my_rank, comm, ierred_local)
      ! if (dbg) THEN
      !    sum_check = 0.0D0
      !    DO IK = 1, NK
      !       DO IS = 1, NSS
      !          DO IC = 1, 3
      !             DO II = 1, NPT
      !                sum_check = sum_check + ABS(GFZ(IS, IC, IK, II))
      !             END DO
      !          END DO
      !       END DO
      !    END DO
      !    PRINT *, 'Rank', my_rank, 'local GFZ sum after =', sum_check
      !    if (my_rank == 0) THEN
      !       WRITE (*, '(A,1ES12.4)') 'max|GFX|=', MAXVAL(ABS(GFX))
      !       WRITE (*, '(A,1ES12.4)') 'max|GFY|=', MAXVAL(ABS(GFY))
      !       WRITE (*, '(A,1ES12.4)') 'max|GFZ|=', MAXVAL(ABS(GFZ))
      !    END IF
      ! end if

      ! --- end-to-end timing
      CALL SYSTEM_CLOCK(t1_all)
      t_global_local = 0.0_dp
      IF (rate > 0) t_global_local = REAL(t1_all - t0_all, dp)/REAL(rate, dp)
      CALL MPI_Reduce(t_global_local, t_global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)

      IF (my_rank == 0 .and. dbg) THEN
         !  WRITE(*,'("t_internal = ",F10.3," s")') t_internal_max
         WRITE (*, '("t_global   = ",F10.3," s")') t_global_max
      END IF

      RETURN
   END SUBROUTINE GF
! !----------------------------------------------------------------------C

   SUBROUTINE Chunked_AllReduce_4D(GFF, NSS, d2, NK, NPT, chunk_size, my_rank, comm, ierred)
      USE iso_fortran_env, ONLY: int64 => int64
      IMPLICIT NONE
      ! COMPLEX(dp), INTENT(INOUT) :: GFF(:, :, :, :)
      COMPLEX(sp), INTENT(INOUT) :: GFF(:, :, :, :)
      INTEGER, INTENT(IN) :: NSS, NK, NPT, my_rank, comm
      INTEGER(int64), INTENT(IN) :: d2, chunk_size
      INTEGER, INTENT(OUT):: ierred

      INTEGER(int64) :: total_size, per_trace, remaining, offset, chunk_len
      INTEGER        :: count32, chunk_start, chunk_end

      INTEGER(int64) :: safe_total, safe_chunk

      ierred = 0

      ! Elements per "trace" (= per fixed 4th index)
      per_trace = INT(NSS, int64)*d2*INT(NK, int64)

#ifdef DEBUG_CHUNK
      safe_total = SIZE(GFF, kind=int64)
#endif

      total_size = per_trace*INT(NPT, int64)

#ifdef DEBUG_CHUNK
      IF (total_size /= safe_total) THEN
         WRITE (*, *) 'Chunked_AllReduce_4D: total_size mismatch: calc=', total_size, ' SIZE(GFF)=', safe_total
         STOP 'Chunked_AllReduce_4D: total_size mismatch'
      END IF
#endif

      ! Branch 1: large array, chunk by NPT so that each chunk <= chunk_size elements
      IF (total_size > chunk_size .OR. total_size > INT(HUGE(count32), int64)) THEN

         remaining = INT(NPT, int64)
         offset = 0_int64

         DO WHILE (remaining > 0_int64)

            ! max traces that fit in chunk_size
            IF (per_trace <= 0_int64) THEN
               STOP 'Chunked_AllReduce_4D: per_trace <= 0'
            END IF

            chunk_len = chunk_size/per_trace
            IF (chunk_len <= 0_int64) chunk_len = 1_int64

            IF (chunk_len > remaining) chunk_len = remaining

            ! 32-bit count safety
            IF (per_trace*chunk_len > INT(HUGE(count32), int64)) THEN
               chunk_len = INT(HUGE(count32), int64)/per_trace
               IF (chunk_len <= 0_int64) THEN
                  STOP 'Chunked_AllReduce_4D: cannot fit even one element into 32-bit count'
               END IF
            END IF

            count32 = INT(per_trace*chunk_len)

            chunk_start = INT(offset) + 1
            chunk_end = chunk_start + INT(chunk_len) - 1

! #ifdef DEBUG_CHUNK
            safe_chunk = SIZE(GFF(:, :, :, chunk_start:chunk_end), kind=int64)
            IF (safe_chunk /= INT(count32, int64)) THEN
               WRITE (*, *) 'Chunked_AllReduce_4D: chunk mismatch'
               WRITE (*, *) '  NSS=', NSS, ' d2=', d2, ' NK=', NK, ' NPT=', NPT
               WRITE (*, *) '  per_trace=', per_trace, ' chunk_len=', chunk_len
               WRITE (*, *) '  chunk_start=', chunk_start, ' chunk_end=', chunk_end
               WRITE (*, *) '  safe_chunk=', safe_chunk, ' count32=', count32
               STOP 'Chunked_AllReduce_4D: chunk SIZE() mismatch'
            END IF
! #endif

            ! CALL MPI_ALLREDUCE(MPI_IN_PLACE, GFF(:, :, :, chunk_start:chunk_end), &
            !                    count32, MPI_COMPLEX16, MPI_SUM, comm, ierred)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, GFF(:, :, :, chunk_start:chunk_end), &
                               count32, MPI_COMPLEX8, MPI_SUM, comm, ierred)

            offset = offset + chunk_len
            remaining = INT(NPT, int64) - offset
         END DO

      ELSE
         ! Branch 2: small enough to reduce in one shot
         IF (total_size > INT(HUGE(count32), int64)) THEN
            STOP 'Chunked_AllReduce_4D: total_size fits chunk_size but not count32'
         END IF

         count32 = INT(total_size)

#ifdef DEBUG_CHUNK
         safe_total = SIZE(GFF, kind=int64)
         IF (safe_total /= INT(count32, int64)) THEN
            WRITE (*, *) 'Chunked_AllReduce_4D: full SIZE mismatch'
            WRITE (*, *) '  safe_total=', safe_total, ' count32=', count32
            STOP 'Chunked_AllReduce_4D: full SIZE mismatch'
         END IF
#endif

         ! CALL MPI_ALLREDUCE(MPI_IN_PLACE, GFF, count32, MPI_COMPLEX16, MPI_SUM, comm, ierred)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, GFF, count32, MPI_COMPLEX8, MPI_SUM, comm, ierred)
      END IF

   END SUBROUTINE Chunked_AllReduce_4D


!    SUBROUTINE GF_test(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                       AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
!                       IE0, IS0, DZ0, GFX, GFY, GFZ, NBLOCK, IG, &
!                       my_rank, n_process, comm, ITER, ND, NS, NR, NSV, NRV, YSR, WTK, WAVELET, SourceScaler)
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: I25D, NTO, NX, NZ, NPT, IANISO, NORD, NK, NNX, NNZ
!       INTEGER, INTENT(IN) :: NSR, NSS, IE0, IS0, ITER, NBLOCK, IFQ
!       INTEGER, INTENT(IN) :: my_rank, n_process, comm
!       INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
!       REAL(dp), INTENT(IN) :: FREQ, DZ0, XTO(:), ZTO(:), X(:), AS(:), WT(:), FKY(:)
!       REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
!       TYPE(InversionGridType), INTENT(IN) :: IG

!       COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
!       INTEGER, INTENT(IN), OPTIONAL :: ND
!       INTEGER, INTENT(IN), OPTIONAL :: NS(:), NR(:), NSV(:), NRV(:)
!       REAL(dp), INTENT(IN), OPTIONAL :: YSR(:), WTK(:)
!       COMPLEX(dp), INTENT(IN), OPTIONAL :: WAVELET(:), SourceScaler(:)

!       ! --- locals (matching GF subroutine variable names and types) ---
!       INTEGER :: KL, KU, LDAB, N
!       INTEGER :: start_index, end_index, my_tasks
!       INTEGER :: NCORE, NCORE_MAX, NCOREI, rate, t0_in, t1_in, t0_all, t1_all
!       REAL(dp) :: t_internal_local, t_internal_max
!       INTEGER :: ierred_local, ierr
!       REAL(dp) :: t_global_local, t_global_max, sum_check
!       REAL(dp) :: AMY, GFMY, SMY, TOL
!       LOGICAL :: use_g0, baseline_ready
!       INTEGER :: ND_loc
!       REAL(dp), PARAMETER :: EPS = 1.0e-14_dp
!       REAL(dp) :: gf_norm_ref, g0_norm_ref, amp_norm_ref
!       COMPLEX(sp), ALLOCATABLE :: GFX_ref(:, :, :, :), GFY_ref(:, :, :, :), GFZ_ref(:, :, :, :)
!       COMPLEX(dp), ALLOCATABLE :: G0_ref(:), G0_curr(:)
!       REAL(dp), ALLOCATABLE :: amp_ref(:), amp_curr(:)
!       INTEGER :: IK, IS, IC, II
!       INTEGER(KIND=8) :: d2_64, chunk_64

!       ! Solver parameters matching GF subroutine
!       INTEGER, PARAMETER :: SOLVER_LU = 1, SOLVER_PARDISO = 2, SOLVER_MUMPS = 3
!       INTEGER, PARAMETER :: SOLVER_LU_NSS_TEST = 1
!       INTEGER, PARAMETER :: SOLVER_LU_HYBRID_TEST = 2
!       INTEGER, PARAMETER :: SOLVER_PARDISO_TEST = 3

!       CALL SYSTEM_CLOCK(COUNT_RATE=rate)
!       CALL SYSTEM_CLOCK(t0_all)

!       ! Setup chunking constants matching GF subroutine
!       d2_64 = 3_8
!       ! Target ~128 MiB per chunk (COMPLEX(dp) = 16 bytes)
!       ! chunk_64 = (128_8*1024_8*1024_8)/16_8
!       chunk_64 = (128_8*1024_8*1024_8)/8_8

!       N = 3*NPT !size of the linear system (number of rows/columns in stiffness matrix): nbr of grid points *3 displacement components
!       KU = 3*(NORD - 1)*(NNZ + 1) + 3
!       KL = KU
!       LDAB = 2*KL + KU + 1

!       AMY = DBLE(FLOAT(2*8*LDAB)/FLOAT(1024))*DBLE(FLOAT(N)/FLOAT(1024*1024))
!       GFMY = DBLE(FLOAT(2*8*NSS*3*NK)/FLOAT(1024))*DBLE(FLOAT(NPT)/FLOAT(1024*1024))*3D0
!       SMY = 16D0*DBLE(N)/(1024D0*1024D0*1024D0)
!       TOL = AMY + GFMY

!       IF (my_rank == 0 .AND. ITER <= 1) THEN
!          WRITE (*, *) 'NSR =', NSR, '   NSS =', NSS
!          ! WRITE (*, '(A, I10, A, I10, A, I10)') '       GQG_dms:', NNX, '    X', NNZ, ' = :', NPT
!     WRITE (*, '(A, I10, A, I10, A, F7.1, A)') '       Max Req_mry A_[m,n] eq banded format:', LDAB, '    X', N, ' =>:', AMY, ' (GB)'
!  WRITE (*, '(A, F10.3, A, F10.3, A, F7.1, A)') '       Total_Req_mry TOL A + 3 GF arrays:', AMY, '    +', GFMY, ' =  ', TOL, ' (GB)'
!          WRITE (*, *)
!       END IF

!       CALL split_work(NK, my_rank, n_process, start_index, end_index, my_tasks)
!       IF (start_index > end_index) RETURN

!       NCORE_MAX = FLOOR((MEM_MAX - (TOL))/TOL) - 1
!       NCORE = MIN(NCORE_MAX, NTHREAD_MAX, my_tasks)
!       NCORE_MAX = FLOOR((MEM_MAX - (TOL))/SMY) - 1
!       NCOREI = MIN(NCORE_MAX, NSS, NTHREAD_MAX)
!       IF ((I25D == 0) .OR. (NCORE <= 0)) NCORE = 1
!       IF ((I25D == 0) .OR. (NCOREI <= 0)) NCOREI = 1

!       use_g0 = PRESENT(ND) .AND. PRESENT(NS) .AND. PRESENT(NR) .AND. PRESENT(NSV) .AND. PRESENT(NRV) &
!                .AND. PRESENT(YSR) .AND. PRESENT(WTK) .AND. PRESENT(WAVELET)
!       baseline_ready = .FALSE.
!       gf_norm_ref = 1.0D0
!       g0_norm_ref = 1.0D0
!       amp_norm_ref = 1.0D0

!       IF (use_g0) THEN
!          ND_loc = ND
!          ALLOCATE (G0_ref(ND_loc), G0_curr(ND_loc), amp_ref(ND_loc), amp_curr(ND_loc))
!       END IF

!       ! Test all solver variants: LU_NSS (baseline), LU_hybrid, and PARDISO
!       CALL RunSolver(SOLVER_LU_NSS_TEST, 'LU_NSS', .TRUE.)
!       CALL RunSolver(SOLVER_LU_HYBRID_TEST, 'LU_hybrid', .FALSE.)
!       CALL RunSolver(SOLVER_PARDISO_TEST, 'PARDISO', .FALSE.)

!       IF (ALLOCATED(GFX_ref)) THEN
!          GFX = GFX_ref
!          GFY = GFY_ref
!          GFZ = GFZ_ref
!          DEALLOCATE (GFX_ref, GFY_ref, GFZ_ref)
!       END IF

!       IF (use_g0) THEN
!          DEALLOCATE (G0_ref, G0_curr, amp_ref, amp_curr)
!       END IF

!       CALL MPI_Barrier(comm, ierr)
!       !CALL mkl_set_dynamic(1)
!       RETURN

!    CONTAINS

!       SUBROUTINE RunSolver(kind, label, is_baseline)
!          INTEGER, INTENT(IN)        :: kind
!          CHARACTER(*), INTENT(IN)   :: label
!          LOGICAL, INTENT(IN)        :: is_baseline

!          REAL(dp) :: t_internal_local_solver, t_internal_max_solver
!          INTEGER :: t_start_solver, t_end_solver
!          REAL(dp) :: diff_gfx_local, diff_gfy_local, diff_gfz_local
!          REAL(dp) :: diff_max_local
!          REAL(dp) :: gfx_local, gfy_local, gfz_local, max_local
!          REAL(dp) :: g0_local, amp_local
!          REAL(dp) :: g0_diff_local, amp_diff_local
!          INTEGER :: ierred_local_test
!          INTEGER :: IK_test, IS_test, IC_test, II_test

!          GFX = (0.0D0, 0.0D0)
!          GFY = (0.0D0, 0.0D0)
!          GFZ = (0.0D0, 0.0D0)

!          CALL MPI_Barrier(comm, ierr)
!          CALL SYSTEM_CLOCK(t_start_solver)

!          IF (my_tasks > 0) THEN
!             SELECT CASE (kind)
!             CASE (SOLVER_LU_NSS_TEST)
!                ! Baseline: per-source LU decomposition (LU_mkl_NSS)
!                CALL run_LU_mkl_NSS(N, KL, KU, LDAB, start_index, end_index, &
!                                    I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                    AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
!                                    IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCOREI)

!             CASE (SOLVER_LU_HYBRID_TEST)
!                ! Hybrid: parallel factorizations across sources (run_LU_hybrid)
!                CALL run_LU_hybrid(N, KL, KU, LDAB, start_index, end_index, &
!                                   I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                   AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
!                                   IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCORE)

!             CASE (SOLVER_PARDISO_TEST)
!                ! PARDISO sparse direct solver
!                CALL run_ik_loop_pardiso(N, KL, KU, LDAB, start_index, end_index, &
!                                         I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                         AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
!                                         IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, NBLOCK, &
!                                         comm, ITER)
!             END SELECT
!          END IF

!          CALL MPI_Barrier(comm, ierr)

!          ! Perform MPI reduction using same chunking as GF subroutine
!          CALL Chunked_AllReduce_4D(GFX, NSS, d2_64, NK, NPT, chunk_64, my_rank, comm, ierred_local_test)
!          CALL Chunked_AllReduce_4D(GFY, NSS, d2_64, NK, NPT, chunk_64, my_rank, comm, ierred_local_test)
!          CALL Chunked_AllReduce_4D(GFZ, NSS, d2_64, NK, NPT, chunk_64, my_rank, comm, ierred_local_test)

!          CALL SYSTEM_CLOCK(t_end_solver)
!          t_internal_local_solver = 0.0_dp
!          IF (rate > 0) t_internal_local_solver = REAL(t_end_solver - t_start_solver, dp)/REAL(rate, dp)

!        CALL MPI_ALLREDUCE(t_internal_local_solver, t_internal_max_solver, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierred_local_test)

!          IF (is_baseline) THEN
!             ! Store baseline results
!             IF (.NOT. ALLOCATED(GFX_ref)) THEN
!                ALLOCATE (GFX_ref(NSS, 3, NK, NPT), GFY_ref(NSS, 3, NK, NPT), GFZ_ref(NSS, 3, NK, NPT))
!             END IF
!             GFX_ref = GFX
!             GFY_ref = GFY
!             GFZ_ref = GFZ

!             gf_norm_ref = MAX(MAXVAL(ABS(GFX_ref)), MAX(MAXVAL(ABS(GFY_ref)), MAXVAL(ABS(GFZ_ref))))
!             IF (gf_norm_ref <= EPS) gf_norm_ref = 1.0D0

!             IF (use_g0) THEN
!                G0_ref = (0.0D0, 0.0D0)
!                amp_ref = 0.0D0

!                CALL Compute_G0(ND_loc, NS, NR, NSV, NRV, VSR, NK, MSR, MSR1, FSR, GFX_ref, GFY_ref, GFZ_ref, &
!                                YSR, FREQ, FKY, WTK, G0_ref, WAVELET, SourceScaler, IFQ, my_rank, amp_ref)
!                g0_norm_ref = MAXVAL(ABS(G0_ref)); IF (g0_norm_ref <= EPS) g0_norm_ref = 1.0D0
!                amp_norm_ref = MAXVAL(ABS(amp_ref)); IF (amp_norm_ref <= EPS) amp_norm_ref = 1.0D0
!             END IF

!             baseline_ready = .TRUE.
!             IF (my_rank == 0) THEN
!                IF (use_g0) THEN
!                   WRITE (*, '("GF_test ",A," runtime=",F10.3," s | max|GF|=",ES12.5," | max|G0|=",ES12.5)') &
!                      TRIM(label), t_internal_max_solver, gf_norm_ref, g0_norm_ref
!                ELSE
!                   WRITE (*, '("GF_test ",A," runtime=",F10.3," s | max|GF|=",ES12.5)') &
!                      TRIM(label), t_internal_max_solver, gf_norm_ref
!                END IF
!             END IF

!          ELSE
!             ! Compare with baseline
!             IF (.NOT. baseline_ready) RETURN

!             gfx_local = MAXVAL(ABS(GFX))
!             gfy_local = MAXVAL(ABS(GFY))
!             gfz_local = MAXVAL(ABS(GFZ))
!             max_local = MAX(gfx_local, MAX(gfy_local, gfz_local))

!             diff_gfx_local = MAXVAL(ABS(GFX - GFX_ref))
!             diff_gfy_local = MAXVAL(ABS(GFY - GFY_ref))
!             diff_gfz_local = MAXVAL(ABS(GFZ - GFZ_ref))
!             diff_max_local = MAX(diff_gfx_local, MAX(diff_gfy_local, diff_gfz_local))

!             g0_local = 0.0D0
!             amp_local = 0.0D0
!             g0_diff_local = 0.0D0
!             amp_diff_local = 0.0D0

!             IF (use_g0) THEN
!                G0_curr = (0.0D0, 0.0D0)
!                amp_curr = 0.0D0
!                CALL Compute_G0(ND_loc, NS, NR, NSV, NRV, VSR, NK, MSR, MSR1, FSR, GFX, GFY, GFZ, &
!                                YSR, FREQ, FKY, WTK, G0_curr, WAVELET, SourceScaler, IFQ, my_rank, amp_curr)
!                g0_diff_local = MAXVAL(ABS(G0_curr - G0_ref))
!                g0_local = MAXVAL(ABS(G0_curr))
!                amp_local = MAXVAL(ABS(amp_curr))
!                amp_diff_local = MAXVAL(ABS(amp_curr - amp_ref))
!             END IF

!             IF (my_rank == 0) THEN
!                IF (use_g0) THEN
!                   WRITE (*, '("GF_test ",A," vs LU_NSS baseline")') TRIM(label)
!                   WRITE (*, '("   runtime        = ",F10.3," s")') t_internal_max_solver
!                   WRITE (*, '("   max|dGF|       = ",ES12.5,"   | rel_diff      = ",ES12.5)') &
!                      max_local, diff_max_local/gf_norm_ref
!                   WRITE (*, '("   max|dG0|       = ",ES12.5,"   | rel_diff_dG0 = ",ES12.5)') &
!                      g0_local, g0_diff_local/g0_norm_ref
!                   WRITE (*, '("   max|damp|      = ",ES12.5,"   | rel_diff_damp = ",ES12.5)') &
!                      amp_local, amp_diff_local/amp_norm_ref
!                ELSE
!                   WRITE (*, '("GF_test ",A," vs LU_NSS baseline")') TRIM(label)
!                   WRITE (*, '("   runtime        = ",F10.3," s")') t_internal_max_solver
!                   WRITE (*, '("   max|dGF|       = ",ES12.5,"   | rel_diff      = ",ES12.5)') &
!                      max_local, diff_max_local/gf_norm_ref
!                END IF
!             END IF
!          END IF

!       END SUBROUTINE RunSolver

!    END SUBROUTINE GF_test

!-----------------------------------------------------------------------
! SUBROUTINE: Compute_G0
!
! PURPOSE:
!   Computes the frequency-domain modeled data vector G0 and the velocity
!   vector V0 for a given source–receiver configuration using the Green’s
!   function components and source–receiver vectors.
!
!   This subroutine performs source–receiver projection, y-direction FFT
!   integration, and wavelet scaling. The output G0 is the modeled data in
!   displacement, and V0 is the corresponding velocity.
!
! ENTRIES:
!   ND       - Total number of data (source-receiver pairs)
!   NS(ND)   - Source indices
!   NR(ND)   - Receiver indices
!   NSV(ND)  - Source vector indices (1–3)
!   NRV(ND)  - Receiver vector indices (1–3)
!   VSR      - Source/receiver polarization vectors
!   NK       - Number of wavenumber samples
!   MSR      - Number of quadrature points per receiver
!   MSR1     - Mapping from receiver to global quadrature index
!   FSR      - Quadrature interpolation weights
!   GFX      - Green’s function component Gij in X
!   GFY      - Green’s function component Gij in Y
!   GFZ      - Green’s function component Gij in Z
!   YSR      - Y-offsets for source/receiver pairs
!   FREQ     - Frequency [Hz]
!   FKY      - Wavenumber sampling array
!   WTK      - Wavenumber weights
!   WAVELET  - Complex-valued source spectrum
!   IFQ      - Frequency index
!   my_rank  - MPI rank (for debug printing)
!
! RETURN:
!   G0(ND)   - Modeled data vector (displacement)
!   amp0(ND) - Amplitude of G0
!
!-----------------------------------------------------------------------

   SUBROUTINE Process_GF(ICSR, IS0, NPT, NK, FKY, WTK, YSR, XP, ZP, NNX, NNZ, NTO, XTO, ZTO, &
                         GFX, GFY, GFZ, VSR, FREQ, PREFIX)

      IMPLICIT NONE

      ! -------- Arguments --------
      INTEGER, INTENT(IN) :: ICSR(:)
      INTEGER, INTENT(IN) :: IS0, NPT, NK, NNX, NNZ, NTO
      REAL(dp), INTENT(IN) :: FKY(:), WTK(:), YSR(:), FREQ
      REAL(dp), INTENT(IN) :: XP(:), ZP(:), XTO(:), ZTO(:)
      COMPLEX(sp), INTENT(IN) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      REAL(dp), INTENT(IN) :: VSR(:, :, :)          
      CHARACTER(1), INTENT(IN) :: PREFIX

      ! -------- Locals --------
      REAL(dp), ALLOCATABLE :: R11(:), R12(:), R21(:), R22(:), R31(:), R32(:)
      REAL(dp), ALLOCATABLE :: FF1(:), FF2(:), ZTO1(:)
      CHARACTER(LEN=80) :: errmsg
      INTEGER :: istat, IC, IP, IK, NSS0, COMPO, IS
      REAL(dp) :: F11, F12, F21, F22, F31, F32, ZMAX
      LOGICAL :: row_active(3)
      REAL(dp) :: eps, YY

      ! -------- Init --------
      NSS0 = 1                  ! processing the first source (as in your original)
      ZMAX = -1.0e10_dp
      eps = 1.0e-12_dp

      ! Quick sanity on ICSR range (for diagnostics only)
      IF ((ICSR(NSS0) < 1) .OR. (ICSR(NSS0) > 3)) THEN
         WRITE (*, *) "ERROR: ICSR(NSS0) has invalid value:", ICSR(NSS0)
         STOP
      END IF

      ! Build active-row mask for the *source* NSS0 using VSR physical rows
      row_active = .FALSE.
      DO IC = 1, 3
         IF (ABS(VSR(NSS0, IC, 1)) > eps .OR. ABS(VSR(NSS0, IC, 2)) > eps .OR. ABS(VSR(NSS0, IC, 3)) > eps) THEN
            row_active(IC) = .TRUE.
         END IF
      END DO

      ! Prepare ZTO1
      ALLOCATE (ZTO1(NTO), STAT=istat, ERRMSG=errmsg); IF (istat /= 0) THEN
         PRINT *, "Error allocating ZTO1:", errmsg; STOP
      END IF
      ZTO1(1:NTO) = ZTO(1:NTO)

      IF (IS0 .EQ. 2) THEN
         DO IP = 1, NPT
            IF (ZP(IP) > ZMAX) ZMAX = ZP(IP)
         END DO
         ZTO1(1:NTO) = ZMAX
      END IF

      ! Allocate transforms
      IF (ALLOCATED(R11)) DEALLOCATE (R11, R12, R21, R22, R31, R32)
      ALLOCATE (R11(NK), R12(NK), R21(NK), R22(NK), R31(NK), R32(NK), STAT=istat, ERRMSG=errmsg)
      IF (istat /= 0) THEN
         PRINT *, "Error allocating R11–R32:", errmsg; STOP
      END IF

      ALLOCATE (FF1(NPT), FF2(NPT), STAT=istat, ERRMSG=errmsg)
      IF (istat /= 0) THEN
         PRINT *, "Error allocating FF1/FF2:", errmsg; STOP
      END IF

      ! -------- Main: loop over physical source rows 1..3, but only if active --------
      DO IS = 1, 3
         IF (.NOT. row_active(IS)) CYCLE
         YY = ABS(YSR(NSS0))

         DO COMPO = 1, 3
            DO IP = 1, NPT
               SELECT CASE (COMPO)
               CASE (1)
                  DO IK = 1, NK
                     R11(IK) = REAL(GFX(NSS0, IS, IK, IP), dp)
                     R12(IK) = AIMAG(GFX(NSS0, IS, IK, IP))
                  END DO

                  CALL FFT_GL(IS, COMPO, NK, FKY, WTK, YY, R11, R12, F11, F12)
                  FF1(IP) = F11
                  FF2(IP) = F12

               CASE (2)
                  DO IK = 1, NK
                     R21(IK) = REAL(GFY(NSS0, IS, IK, IP), dp)
                     R22(IK) = AIMAG(GFY(NSS0, IS, IK, IP))
                  END DO

                  CALL FFT_GL(IS, COMPO, NK, FKY, WTK, YY, R21, R22, F21, F22)
                  FF1(IP) = F21
                  FF2(IP) = F22

               CASE (3)
                  DO IK = 1, NK
                     R31(IK) = REAL(GFZ(NSS0, IS, IK, IP), dp)
                     R32(IK) = AIMAG(GFZ(NSS0, IS, IK, IP))
                  END DO

                  CALL FFT_GL(IS, COMPO, NK, FKY, WTK, YY, R31, R32, F31, F32)
                  FF1(IP) = F31
                  FF2(IP) = F32
               END SELECT
            END DO

            CALL Output_GF(FREQ, PREFIX, IS, COMPO, NNX, NNZ, XP, ZP, FF1, FF2, NTO, XTO, ZTO1)
         END DO
      END DO

      ! -------- Clean up --------
      DEALLOCATE (FF1, FF2, ZTO1, R11, R12, R21, R22, R31, R32)

   END SUBROUTINE Process_GF

   SUBROUTINE Compute_G0(ND, NS, NR, NSV, NRV, VSR, NK, &
                         MSR, MSR1, FSR, GFX, GFY, GFZ, &
                         YSR, FREQ, FKY, WTK, G0, WAVELET, SourceScaler, IFQ, my_rank, amp0, DEBUG_OUTPUT)

      USE OMP_LIB
      IMPLICIT NONE

      !-------------------- Arguments --------------------
      INTEGER, INTENT(IN)             :: ND, NK, IFQ, my_rank
      INTEGER, INTENT(IN)             :: NS(:), NR(:), NSV(:), NRV(:)
      INTEGER, INTENT(IN)             :: MSR(:), MSR1(:, :)
      REAL(dp), INTENT(IN)            :: YSR(:), FSR(:, :), VSR(:, :, :)
      REAL(dp), INTENT(IN)            :: FKY(:), WTK(:), FREQ
      COMPLEX(sp), INTENT(IN)         :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      COMPLEX(dp), INTENT(IN)         :: WAVELET(:)
      COMPLEX(dp), INTENT(IN)         :: SourceScaler(:)

      COMPLEX(dp), INTENT(INOUT)      :: G0(:)
      REAL(dp), INTENT(INOUT)         :: amp0(:)
      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT

      !-------------------- Locals -----------------------
      INTEGER                         :: I, IK, IS, IR, ICS, ICR, K
      INTEGER                         :: NCORE, rate, t_start, t_end
      REAL(dp)                        :: t_elapsed, YY
      COMPLEX(dp)                     :: GX, GY, GZ
      REAL(dp)                        :: F11, F12, F21, F22, F31, F32

      ! Per-thread allocatable work arrays (declared here; allocated inside OMP region)
      COMPLEX(dp), ALLOCATABLE        :: S(:), G(:)
      REAL(dp), ALLOCATABLE           :: R11(:), R12(:), R21(:), R22(:), R31(:), R32(:)
      LOGICAL :: dbg
      INTEGER :: tid, udbg, ierr, ierred
      CHARACTER(LEN=64) :: fname
      ! --- extra locals for debugging ---
      INTEGER            :: IP, KPEEK(3), KMID, KMX
      LOGICAL            :: WATCHIK, WATCHK
      COMPLEX(dp)        :: CSUMX, CSUMY, CSUMZ
      COMPLEX(dp)        :: X11, X12, X13, Y11, Y12, Y13, Z11, Z12, Z13

      REAL(dp) :: WC(NK), WS(NK)
      LOGICAL :: even_parity
      !-------------------- Setup ------------------------
      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT
      ! IF (SIZE(SourceScaler) < MAXVAL(NS)) THEN
      !    IF (my_rank == 0) WRITE (*, *) 'Compute_G0: SourceScaler size < max(NS).'
      !    STOP
      ! END IF

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      !---------------------------------------------
      ! Debug-time sanity checks on sizes / indices
      !---------------------------------------------
      ! IF (dbg .AND. my_rank == 0) THEN
      !    WRITE(*,'("Compute_G0 entry: ND=",I10," NK=",I10)') ND, NK
      !    WRITE(*,'("  SIZE(NS)=",I10," SIZE(NR)=",I10," SIZE(NSV)=",I10," SIZE(NRV)=",I10)') &
      !         SIZE(NS), SIZE(NR), SIZE(NSV), SIZE(NRV)
      !    WRITE(*,'("  SHAPE(VSR)=",3(I10,","))') SHAPE(VSR)
      !    WRITE(*,'("  SHAPE(GFX)=",4(I10,","))') SHAPE(GFX)
      !    WRITE(*,'("  SHAPE(G0) =",I10)') SIZE(G0)
      !    WRITE(*,'("  SHAPE(MSR)=",I10," SHAPE(MSR1)=",2(I10,","))') SIZE(MSR), SHAPE(MSR1)
      !    WRITE(*,'("  SHAPE(FSR)=",2(I10,","))') SHAPE(FSR)
      ! END IF

      IF (dbg) THEN
         ! ND must not exceed mapping arrays
         IF (ND > SIZE(NS) .OR. ND > SIZE(NR) .OR. ND > SIZE(NSV) .OR. ND > SIZE(NRV)) THEN
            WRITE (*, *) 'Compute_G0: ND=', ND, '  SIZE(NS)=', SIZE(NS), ' SIZE(NR)=', SIZE(NR), &
               ' SIZE(NSV)=', SIZE(NSV), ' SIZE(NRV)=', SIZE(NRV)
            STOP 'Compute_G0: ND larger than NS/NR/NSV/NRV'
         END IF

         ! G0, amp0 must be at least ND
         IF (ND > SIZE(G0)) THEN
            WRITE (*, *) 'Compute_G0: ND=', ND, ' SIZE(G0)=', SIZE(G0)
            STOP 'Compute_G0: G0 size mismatch'
         END IF
         IF (ND > SIZE(amp0)) THEN
            WRITE (*, *) 'Compute_G0: ND=', ND, ' SIZE(amp0)=', SIZE(amp0)
            STOP 'Compute_G0: amp0 size mismatch'
         END IF

         ! NK must fit the Green’s function arrays
         IF (NK > SIZE(GFX, 3) .OR. NK > SIZE(GFY, 3) .OR. NK > SIZE(GFZ, 3)) THEN
            WRITE (*, *) 'Compute_G0: NK=', NK, ' SHAPE(GFX)=', SHAPE(GFX)
            STOP 'Compute_G0: NK larger than GFX/GFY/GFZ third dim'
         END IF

         ! SourceScaler / Wavelet lengths
         IF (MAXVAL(NS(1:ND)) > SIZE(SourceScaler)) THEN
            WRITE (*, *) 'Compute_G0: max(NS) =', MAXVAL(NS(1:ND)), ' SIZE(SourceScaler)=', SIZE(SourceScaler)
            STOP 'Compute_G0: SourceScaler too short'
         END IF
         IF (IFQ < 1 .OR. IFQ > SIZE(WAVELET)) THEN
            WRITE (*, *) 'Compute_G0: IFQ=', IFQ, ' SIZE(WAVELET)=', SIZE(WAVELET)
            STOP 'Compute_G0: IFQ out of range for WAVELET'
         END IF

         ! Check NS/NR indices against VSR & MSR
         IF (MAXVAL(NS(1:ND)) > SIZE(VSR, 1) .OR. MINVAL(NS(1:ND)) < 1) THEN
            WRITE (*, *) 'Compute_G0: NS out of [1, SIZE(VSR,1)]'
            STOP 'Compute_G0: NS index out of range'
         END IF
         IF (MAXVAL(NR(1:ND)) > SIZE(VSR, 1) .OR. MINVAL(NR(1:ND)) < 1) THEN
            WRITE (*, *) 'Compute_G0: NR out of [1, SIZE(VSR,1)]'
            STOP 'Compute_G0: NR index out of range'
         END IF

         ! MSR(I) must fit MSR1 and FSR
         IF (MAXVAL(NR(1:ND)) > SIZE(MSR)) THEN
            WRITE (*, *) 'Compute_G0: max(NR)=', MAXVAL(NR(1:ND)), ' SIZE(MSR)=', SIZE(MSR)
            STOP 'Compute_G0: NR exceeds MSR size'
         END IF
         IF (MAXVAL(MSR(1:SIZE(MSR))) > SIZE(MSR1, 2) .OR. MAXVAL(MSR(1:SIZE(MSR))) > SIZE(FSR, 2)) THEN
            WRITE (*, *) 'Compute_G0: max(MSR)=', MAXVAL(MSR), ' SHAPE(MSR1)=', SHAPE(MSR1), ' SHAPE(FSR)=', SHAPE(FSR)
            STOP 'Compute_G0: MSR too large for MSR1/FSR'
         END IF

         ! IP = MSR1(IR,K) must be within GF field size
         IF (ANY(MSR(1:SIZE(MSR)) < 0)) THEN
            WRITE (*, *) 'Compute_G0: negative MSR entry detected'
            STOP 'Compute_G0: MSR negative'
         END IF
      END IF

     
      CALL mkl_set_num_threads(1)

      NCORE = MAX(1, MIN(omp_get_max_threads(), ND))
      CALL omp_set_num_threads(NCORE)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(I,IS,IR,ICS,ICR,IK,K,GX,GY,GZ,YY,F11,F12,F21,F22,F31,F32, &
!$OMP         t_start,t_end,t_elapsed,ierr,S,G,R11,R12,R21,R22,R31,R32,WC,WS,even_parity, udbg, tid, &
!$OMP         IP,KPEEK,KMID,KMX,WATCHIK,WATCHK,CSUMX,CSUMY,CSUMZ, &
!$OMP         X11,X12,X13,Y11,Y12,Y13,Z11,Z12,Z13)
      tid = 0

      ALLOCATE (S(3), G(3), STAT=ierr); IF (ierr /= 0) STOP 'alloc S/G'
      ALLOCATE (R11(NK), R12(NK), R21(NK), R22(NK), R31(NK), R32(NK), STAT=ierr)
      IF (ierr /= 0) STOP 'alloc R**'

      ! IF (dbg) THEN
      !    tid = omp_get_thread_num()
      !    WRITE (fname, '("g0_debug_t",I0,".txt")') tid
      !    OPEN (NEWUNIT=udbg, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')
      ! ELSE
      !    udbg = -1
      ! END IF

!$OMP DO SCHEDULE(static)
      DO I = 1, ND
         CALL SYSTEM_CLOCK(t_start, rate)

         IS = NS(I)
         IR = NR(I)
         ICS = NSV(I)
         ICR = NRV(I)
         YY = ABS(YSR(IS) - YSR(IR))
         CALL FFT_GL_PRECOMP(ICS, ICR, NK, FKY, WTK, YY, WC, WS, even_parity)

          
         ! guards (catch bad mapping early)
         IF (ICS < 1 .OR. ICS > 3 .OR. ICR < 1 .OR. ICR > 3) THEN
            IF (dbg) WRITE (*, *) 'Compute_G0: bad ICS/ICR at I=', I, ' ICS=', ICS, ' ICR=', ICR
            CYCLE
         END IF

         S(1:3) = CMPLX(VSR(IS, ICS, 1:3), 0.0_dp, dp)
         G(1:3) = CMPLX(VSR(IR, ICR, 1:3), 0.0_dp, dp)

         ! IF (dbg) THEN
         !    WRITE (udbg, '(A,I6,2X,A,4I4,2X,A,3I2,2X,A,3I2)') &
         !       'I=', I, ' IS/ICS/IR/ICR=', IS, ICS, IR, ICR, &
         !       ' S=', NINT(REAL(S(1), dp)), NINT(REAL(S(2), dp)), NINT(REAL(S(3), dp)), &
         !       ' G=', NINT(REAL(G(1), dp)), NINT(REAL(G(2), dp)), NINT(REAL(G(3), dp))
         ! END IF

         R11 = 0.0_dp; R12 = 0.0_dp
         R21 = 0.0_dp; R22 = 0.0_dp
         R31 = 0.0_dp; R32 = 0.0_dp

         DO IK = 1, NK
            GX = CMPLX(0.0_dp, 0.0_dp, dp)
            GY = CMPLX(0.0_dp, 0.0_dp, dp)
            GZ = CMPLX(0.0_dp, 0.0_dp, dp)

            ! ---- one-time per-datum summary of MSR/MSR1/FSR and peek a few entries
            ! IF (dbg .AND. IK == 1) THEN
            !    WRITE (udbg, '("I=",I6,"  IR=",I6,"  MSR(IR)=",I6)') I, IR, MSR(IR)
            !    IF (MSR(IR) >= 1) THEN
            !       KMX = MSR(IR)
            !       KMID = MAX(1, KMX/2)
            !       KPEEK(1) = 1
            !       KPEEK(2) = KMID
            !       KPEEK(3) = KMX
            !       WRITE (udbg, '(A,3(I6,1X))') '  K-peek:', KPEEK
            !       WRITE (udbg, '(A)') '  (K, FSR(IR,K), MSR1(IR,K))'
            !       WRITE (udbg, '(I6,1X,ES12.4,1X,I8)') KPEEK(1), FSR(IR, KPEEK(1)), MSR1(IR, KPEEK(1))
            !       WRITE (udbg, '(I6,1X,ES12.4,1X,I8)') KPEEK(2), FSR(IR, KPEEK(2)), MSR1(IR, KPEEK(2))
            !       WRITE (udbg, '(I6,1X,ES12.4,1X,I8)') KPEEK(3), FSR(IR, KPEEK(3)), MSR1(IR, KPEEK(3))
            !    END IF
            ! END IF

            DO K = 1, MSR(IR)
               IP = MSR1(IR, K)
            !    IF (dbg .and. my_rank == 0) THEN
            !       IF (IP < 1 .OR. IP > SIZE(GFX, 4)) THEN
            !          WRITE (*, *) 'Compute_G0: bad IP=', IP, ' for IR=', IR, ' K=', K, ' NPT=', SIZE(GFX, 4)
            !          STOP 'Compute_G0: IP index out of range for GFX/GFY/GFZ'
            !       END IF
            !    END IF
               ! ---- raw tensor entries at this (IS,IK,IP)
               X11 = GFX(IS, 1, IK, IP); X12 = GFX(IS, 2, IK, IP); X13 = GFX(IS, 3, IK, IP)
               Y11 = GFY(IS, 1, IK, IP); Y12 = GFY(IS, 2, IK, IP); Y13 = GFY(IS, 3, IK, IP)
               Z11 = GFZ(IS, 1, IK, IP); Z12 = GFZ(IS, 2, IK, IP); Z13 = GFZ(IS, 3, IK, IP)

               CSUMX = X11*G(1) + X12*G(2) + X13*G(3)
               CSUMY = Y11*G(1) + Y12*G(2) + Y13*G(3)
               CSUMZ = Z11*G(1) + Z12*G(2) + Z13*G(3)

               GX = GX + CMPLX(FSR(IR, K), 0.0_dp, dp)*CSUMX
               GY = GY + CMPLX(FSR(IR, K), 0.0_dp, dp)*CSUMY
               GZ = GZ + CMPLX(FSR(IR, K), 0.0_dp, dp)*CSUMZ

               ! ---- detailed peek: a few K for two IKs
   !             KMX = MSR(IR)
   !             KMID = MAX(1, KMX/2)
   !             WATCHIK = (IK == 1 .OR. IK == NK/2)
   !             WATCHK = (K == 1 .OR. K == KMID .OR. K == KMX)
   !             IF (dbg .AND. WATCHIK .AND. WATCHK) THEN
   !                WRITE (udbg, '("I=",I6," IK=",I6," K=",I6," IP=",I8,"  FSR=",ES12.4)') I, IK, K, IP, FSR(IR, K)
   !                WRITE (udbg, '("  GFX row:  (",ES12.4,",",ES12.4,") (",ES12.4,",",ES12.4,") (",ES12.4,",",ES12.4,")")') &
   !                   REAL(X11, dp), AIMAG(X11), REAL(X12, dp), AIMAG(X12), REAL(X13, dp), AIMAG(X13)
   !                WRITE (udbg, '("  GFY row:  (",ES12.4,",",ES12.4,") (",ES12.4,",",ES12.4,") (",ES12.4,",",ES12.4,")")') &
   !                   REAL(Y11, dp), AIMAG(Y11), REAL(Y12, dp), AIMAG(Y12), REAL(Y13, dp), AIMAG(Y13)
   !                WRITE (udbg, '("  GFZ row:  (",ES12.4,",",ES12.4,") (",ES12.4,",",ES12.4,") (",ES12.4,",",ES12.4,")")') &
   !                   REAL(Z11, dp), AIMAG(Z11), REAL(Z12, dp), AIMAG(Z12), REAL(Z13, dp), AIMAG(Z13)
   !  WRITE(udbg,'("  S(:)=",3("(",ES12.4,",",ES12.4,")",1X))') REAL(S(1),dp),AIMAG(S(1)), REAL(S(2),dp),AIMAG(S(2)), REAL(S(3),dp),AIMAG(S(3))
   !                WRITE (udbg, '("  CSUMX=",2ES12.4,"  CSUMY=",2ES12.4,"  CSUMZ=",2ES12.4)') &
   !                   REAL(CSUMX, dp), AIMAG(CSUMX), REAL(CSUMY, dp), AIMAG(CSUMY), REAL(CSUMZ, dp), AIMAG(CSUMZ)
   !             END IF
            END DO

            R11(IK) = REAL(GX, dp); R12(IK) = AIMAG(GX)
            R21(IK) = REAL(GY, dp); R22(IK) = AIMAG(GY)
            R31(IK) = REAL(GZ, dp); R32(IK) = AIMAG(GZ)

            ! IF (dbg .AND. (IK == NK/5 .OR. IK == NK/2)) THEN
            !    WRITE (udbg, '(A,I6,2X,A,I5,2X,A,2ES12.4,2X,A,2ES12.4,2X,A,2ES12.4)') &
            !       'I=', I, 'IK=', IK, 'GX=', REAL(GX, dp), AIMAG(GX), 'GY=', REAL(GY, dp), AIMAG(GY), 'GZ=', REAL(GZ, dp), AIMAG(GZ)
            ! END IF
         END DO

        
         ! IF (dbg) WRITE (udbg, '(A,I6,2X,A,F10.3)') 'I=', I, 'YY=', YY
         ! CALL FFT_GL_W(ICS, ICR, NK, WC, WS, even_parity, R11, R12, F11, F12)
         ! CALL FFT_GL_W(ICS, ICR, NK, WC, WS, even_parity, R21, R22, F21, F22)
         ! CALL FFT_GL_W(ICS, ICR, NK, WC, WS, even_parity, R31, R32, F31, F32)
         CALL FFT_GL(ICS, ICR, NK, FKY, WTK, YY, R11, R12, F11, F12)
         CALL FFT_GL(ICS, ICR, NK, FKY, WTK, YY, R21, R22, F21, F22)
         CALL FFT_GL(ICS, ICR, NK, FKY, WTK, YY, R31, R32, F31, F32)

         G0(I) = (S(1)*CMPLX(F11, F12, dp) + S(2)*CMPLX(F21, F22, dp) + S(3)*CMPLX(F31, F32, dp))*WAVELET(IFQ)*SourceScaler(IS)
         amp0(I) = ABS(G0(I))

         ! IF (dbg) WRITE (udbg, '(A,I6,2X,A,ES17.10)') 'I=', I, '|G0|=', amp0(I)

         CALL SYSTEM_CLOCK(t_end)
         t_elapsed = REAL(t_end - t_start, dp)/REAL(rate, dp)
      END DO
!$OMP END DO
      ! IF (dbg) CLOSE (udbg)
      DEALLOCATE (R11, R12, R21, R22, R31, R32)
      DEALLOCATE (S, G)
!$OMP END PARALLEL

      CALL omp_set_num_threads(1)


   END SUBROUTINE Compute_G0

   SUBROUTINE Calibrate_Source(IFQ, ND, GT0, G0, &
                               SourceScaler, Wd_Acq, &
                               NS, NSS, &
                               my_rank, DEBUG_OUTPUT)
      IMPLICIT NONE

      !---- inputs
      INTEGER, INTENT(IN)    :: IFQ, ND, NSS, my_rank
      COMPLEX(dp), INTENT(IN)  :: GT0(:)             ! size >= ND
      INTEGER, INTENT(IN)      :: NS(:)              ! (unused for global scaler)
      REAL(dp), INTENT(IN)     :: Wd_Acq(:)           

      !---- in/out
      COMPLEX(dp), INTENT(INOUT) :: G0(:)              ! size >= ND; updated in-place (global scaling)

      !---- outputs
      COMPLEX(dp), INTENT(OUT)   :: SourceScaler(:)   ! size >= NSS; all set to the same global scaler
    
      !---- optional
      LOGICAL, OPTIONAL, INTENT(IN)    :: DEBUG_OUTPUT

      !---- locals
      COMPLEX(dp) :: num, den, alpha
      REAL(dp), PARAMETER :: EPS_DEN = 1.0e-30_dp
      REAL(dp) :: phase_deg, amag
      LOGICAL :: dbg
      INTEGER :: i

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      ! ---- accumulate over ALL traces (GLOBAL calibration) ----
      num = CMPLX(0.0_dp, 0.0_dp, dp)
      den = CMPLX(0.0_dp, 0.0_dp, dp)
DO i = 1, ND
   IF (WD_acq(i) <= 0.0_dp) CYCLE
   num = num + CONJG(G0(i)) * GT0(i) * WD_acq(i)
   den = den + CONJG(G0(i)) * G0(i)  * WD_acq(i)
END DO

      IF (ABS(den) > EPS_DEN) THEN
         alpha = num/den
      ELSE
         alpha = CMPLX(1.0_dp, 0.0_dp, dp)
         IF (dbg .AND. my_rank == 0) WRITE (*, '("IFQ=",I0," global: tiny denominator; scaler=1")') IFQ
      END IF

      ! ---- apply to all G0 ----
      DO i = 1, ND
         G0(i) = G0(i)*alpha
      END DO

      ! ---- set outputs ----
      SourceScaler(1:NSS) = alpha        ! same global scaler copied to all entries
      amag = ABS(alpha)
      ! W_Q_S = 1.0_dp/MAX(EPS_DEN, amag)  ! scalar residual weight

      IF (dbg .AND. my_rank == 0) THEN
         phase_deg = ATAN2(AIMAG(alpha), REAL(alpha, dp))*180.0_dp/PI
         WRITE (*, '("IFQ=",I0,"  GLOBAL  |scaler|=",ES12.4,"  phase(deg)=",F8.3,"  ND=",I0)') &
            IFQ, amag, phase_deg, ND
      END IF

   END SUBROUTINE Calibrate_Source

   SUBROUTINE Calibrate_PerSource(IFQ, ND, GT0, G0, &
                                  SourceScaler, Wd_Acq, &
                                  NS, NSS, &
                                  my_rank, DEBUG_OUTPUT)
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: IFQ, ND, NSS, my_rank
      COMPLEX(dp), INTENT(IN)    :: GT0(:)            ! size ND
      COMPLEX(dp), INTENT(INOUT) :: G0(:)             ! size ND (updated in-place)
      COMPLEX(dp), INTENT(OUT)   :: SourceScaler(:)  ! size NSS (complex per-source)
      REAL(dp), INTENT(IN)      :: Wd_Acq(:)        ! size ND (acquisition weights)
      INTEGER, INTENT(IN)    :: NS(:)         ! size ND, values in 1..NSS
      LOGICAL, OPTIONAL, INTENT(IN)    :: DEBUG_OUTPUT

      COMPLEX(dp) :: num, den
      REAL(dp), PARAMETER :: EPS_DEN = 1.0e-30_dp
      LOGICAL :: dbg
      INTEGER :: s, i, hit_count

      ! ---- basic size checks ----
      ! IF (SIZE(SourceScaler) < NSS) THEN
      !    IF (my_rank == 0) WRITE(*,*) 'Calibrate_PerSource: SourceScaler too small.'
      !    STOP
      ! END IF
      ! IF (SIZE(W_Q_S) < NSS) THEN
      !    IF (my_rank == 0) WRITE(*,*) 'Calibrate_PerSource: W_Q_S too small.'
      !    STOP
      ! END IF
      ! IF (SIZE(NS) < ND) THEN
      !    IF (my_rank == 0) WRITE(*,*) 'Calibrate_PerSource: NS too small.'
      !    STOP
      ! END IF

      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      ! Initialize outputs
      SourceScaler(1:NSS) = CMPLX(1.0_dp, 0.0_dp, dp)
      ! W_Q_S = 1.0_dp

      DO s = 1, NSS
         num = CMPLX(0.0_dp, 0.0_dp, dp)
         den = CMPLX(0.0_dp, 0.0_dp, dp)
         hit_count = 0

         DO i = 1, ND
            IF (NS(i) == s) THEN
                IF (Wd_Acq(i) <= 0.0_dp) CYCLE
              num = num + Wd_Acq(i) * CONJG(G0(i)) * GT0(i)
den = den + Wd_Acq(i) * CONJG(G0(i)) * G0(i)
               hit_count = hit_count + 1
            
            END IF
         END DO

         IF (hit_count == 0) THEN
            SourceScaler(s) = CMPLX(1.0_dp, 0.0_dp, dp)
            IF (dbg .AND. my_rank == 0) WRITE (*, '("Calib src ",I0,": no traces; SourceScaler=1")') s
         ELSEIF (ABS(den) > EPS_DEN) THEN
            SourceScaler(s) = num/den
         ELSE
            SourceScaler(s) = CMPLX(1.0_dp, 0.0_dp, dp)
            IF (dbg .AND. my_rank == 0) WRITE (*, '("Calib src ",I0,": tiny den; SourceScaler=1")') s
         END IF

         ! Apply found scaler to G0 for this source (current frequency)
         DO i = 1, ND
            IF (NS(i) == s) G0(i) = G0(i)*SourceScaler(s)
         END DO

         !  W_Q_S(s) = 1.0D0 / MAX(EPS_DEN, ABS(SourceScaler(s)))

         IF (dbg .AND. my_rank == 0) THEN
            WRITE (*, '("IFQ=",I0,"  src=",I0,"  |SourceScaler|=",ES12.4,"  phase(deg)=",F8.3, &
    &                 "  traces=",I0)') IFQ, s, ABS(SourceScaler(s)), &
                     ATAN2(DIMAG(SourceScaler(s)), DBLE(SourceScaler(s)))*180.0D0/3.141592653589793D0, hit_count
         END IF
      END DO
   END SUBROUTINE Calibrate_PerSource
!    SUBROUTINE GF_ori(I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                      AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, &
!  IE0, IS0, DZ0, GFX, GFY, GFZ, NBLOCK, IG, my_rank, n_process, comm, ITER)

!       IMPLICIT real(dp) (A - H, O - Z)
!       INTEGER, INTENT(IN) :: I25D, NTO, NX, NZ, NPT, IANISO, NORD, NK, NNX, NNZ
!       INTEGER, INTENT(IN) :: NSR, NSS, IE0, IS0, ITER, NBLOCK,IFQ
!       INTEGER, INTENT(IN) :: my_rank, n_process, comm
!       INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
!       real(dp), INTENT(IN) :: FREQ, DZ0, XTO(:), ZTO(:), X(:), AS(:), WT(:), FKY(:)
!       real(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
!       COMPLEX(dp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
!       TYPE(InversionGridType), INTENT(IN) :: IG
!       ! locals
!       INTEGER :: KL, KU, LDAB, N, ierred
!       INTEGER :: start_index, end_index, modeling_iterations_per_process, remainder, my_tasks
!       real(dp) :: AMY, GFMY, SMY, TOL
!       INTEGER :: NCORE, NCORE_MAX, NCOREI, rate, ierr, SOLVER_KIND
!       INTEGER, PARAMETER :: SOLVER_LU = 1, SOLVER_PARDISO = 2, SOLVER_MUMPS = 3
!       INTEGER :: ierred_local

!       CALL SYSTEM_CLOCK(COUNT_RATE=rate)

!       GFX = (0.0D0, 0.0D0); GFY = (0.0D0, 0.0D0); GFZ = (0.0D0, 0.0D0)

!       IF (my_rank == 0) THEN
!          WRITE (*, *) ''
!          WRITE (*, *) '-------------- Forward modeling -----------'
!          WRITE (*, *) 'IFREQ =', IFQ, '   ITERATION =', ITER
!       END IF

!       N = 3*NPT
!       KU = 3*(NORD - 1)*(NNZ + 1) + 3
!       KL = KU
!       LDAB = 2*KL + KU + 1

!       AMY = DBLE(FLOAT(2*8*LDAB)/FLOAT(1024))*DBLE(FLOAT(N)/FLOAT(1024*1024))
!       GFMY = DBLE(FLOAT(2*8*NSS*3*NK)/FLOAT(1024))*DBLE(FLOAT(NPT)/FLOAT(1024*1024))*3
!       SMY = 16D0*DBLE(N)/(1024D0*1024D0*1024D0)
!       TOL = AMY + GFMY
!       IF (my_rank == 0 .AND. ITER <= 1) THEN
!          WRITE (*, *) 'NSR =', NSR, '   NSS =', NSS
!          WRITE (*, '(A, I10, A, I10, A, I10)') '       GQG_dms:', NNX, '    X', NNZ, ' = :', NPT
!          WRITE (*, '(A, I10, A, I10, A, F7.1, A)') '       Max Req_mry A_[m,n]:', LDAB, '    X', N, ' =>:', AMY, ' (GB)'
!          WRITE (*, '(A, F10.3, A, F10.3, A, F7.1, A)') '       Total_Req_mry TOL:', AMY, '    +', GFMY, ' =  ', TOL, ' (GB)'
!          WRITE (*, *)
!       END IF
!       CALL split_work(NK, my_rank, n_process, start_index, end_index, my_tasks)

!       SELECT CASE (SOLVER_KIND)

!       CASE (SOLVER_LU)
!          ! Keep your automatic LU choice here (examples shown)
!          NCORE_MAX = FLOOR((MEM_MAX - TOL)/AMY) - 1
!          IF (NCORE_MAX < 1) NCORE_MAX = 1
!          NCORE = MIN(NCORE_MAX, NTHREAD_MAX, my_tasks)

!          NCORE_MAX = FLOOR((MEM_MAX - TOL)/SMY) - 1
!          IF (NCORE_MAX < 1) NCORE_MAX = 1
!          NCOREI = MIN(NCORE_MAX, NSS, NTHREAD_MAX)

!          IF ((I25D == 0) .OR. (NCORE <= 0)) NCORE = 1
!          IF ((I25D == 0) .OR. (NCOREI <= 0)) NCOREI = 1

!          IF (NCORE > 2*n_process .AND. NSS < 11) THEN
!             CALL run_LU_hybrid(N, KL, KU, LDAB, start_index, end_index, &
!                                I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
!                                IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCORE)

!          ELSE
!             CALL run_LU_mkl_NSS(N, KL, KU, LDAB, start_index, end_index, &
!                                 I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                 AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
!                                 IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCOREI)
!          END IF

!       CASE (SOLVER_PARDISO)
!          CALL run_ik_loop_pardiso(N, KL, KU, LDAB, start_index, end_index, &
!                                   I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                   AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
!                                   IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, NBLOCK, &
!                                   comm, ITER)

!       CASE (SOLVER_MUMPS)
!          ! CALL run_ik_loop_mumps(...)

!       CASE DEFAULT
!          IF (my_rank == 0) WRITE (*, *) 'GF_production: unknown SOLVER_KIND =', SOLVER_KIND
!          CALL MPI_Abort(comm, 1, ierr)
!       END SELECT

!       CALL MPI_Barrier(comm, ierr)

!       CALL Chunked_AllReduce_4D(GFX, NSS, 3_8, NK, NPT,  my_rank, comm, ierred_local)
!       IF (ierred_local /= MPI_SUCCESS) THEN
!    IF (my_rank == 0) WRITE(*,*) 'Error in Chunked_AllReduce_4D for GFX, ier =', ierred_local
!    CALL MPI_Abort(comm, ierred_local, ierr)
! END IF
!       CALL Chunked_AllReduce_4D(GFY, NSS, 3_8, NK, NPT,  my_rank, comm, ierred_local)
!       CALL Chunked_AllReduce_4D(GFZ, NSS, 3_8, NK, NPT,  my_rank, comm, ierred_local)
!       sum_check = 0.0D0

!       CALL MPI_Barrier(comm, ierr)
!       !CALL mkl_set_dynamic(1)
!       RETURN

!    END SUBROUTINE GF_ori

END module F_modeling
