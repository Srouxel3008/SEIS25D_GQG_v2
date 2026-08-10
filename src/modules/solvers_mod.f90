module solvers_mod

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
   USE stiffness_assembly_mod
   use iso_fortran_env, only: dp => real64, sp => real32
   IMPLICIT none

contains


   SUBROUTINE run_LU_hybrid(N, KL, KU, LDAB, i1, i2, &
                            I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                            AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                            IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCORE)

      USE OMP_LIB
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, KL, KU, LDAB, i1, i2, I25D, IFQ, NTO, NX, NZ, NNX, NNZ, NPT, IANISO, NORD
      INTEGER, INTENT(IN) :: NK, NSR, NSS, IE0, IS0, my_rank, comm, ITER, NCORE
      INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FREQ, FKY(:)
      REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
      REAL(dp), INTENT(IN) :: DZ0
      TYPE(InversionGridType), INTENT(IN) :: IG
      COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE :: A(:, :), S(:, :)
      REAL(dp) :: SV(3, 3)
      INTEGER, ALLOCATABLE :: IPIV(:)

      REAL(dp) :: FK, t_elapsed, diff_norm, diff_min, diff_max
      REAL(dp) :: t_ca_elapsed, t_fact_elapsed, t_solve_elapsed
      INTEGER :: IS, IC, II, IK, INFO, ierr_alloc, myid, idx, rate, t_start, t_end
      INTEGER :: t_ca_start, t_ca_end, t_fact_start, t_fact_end, t_solve_start, t_solve_end
      INTEGER :: J, ICS
      LOGICAL :: dbg
      INTEGER :: min_idx_used, max_idx_used
      INTEGER :: min_local, max_local

      CALL mkl_set_num_threads(1)         ! MKL single-threaded inside OMP
      CALL omp_set_num_threads(NCORE)

      IF (my_rank == 0 .AND. ITER <= 1) WRITE (*, *) 'run_LU_hybrid: NCORE =', NCORE
      IF (omp_in_parallel()) THEN
         WRITE (*, *) 'ERROR: run_LU_hybrid entered inside a parallel region'
         STOP
      END IF
      IF (i1 > i2) RETURN

      dbg = (ITER <= 1)
      min_idx_used = HUGE(min_idx_used)
      max_idx_used = -HUGE(max_idx_used)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IK,FK,IPIV,A,SV,S,INFO,ierr_alloc,myid,IS,IC,II,J,&
!$OMP                                  t_start,t_end,t_elapsed,idx,diff_norm,diff_min,diff_max,min_local,max_local,&
!$OMP                                  t_ca_start,t_ca_end,t_fact_start,t_fact_end,t_solve_start,t_solve_end,&
!$OMP                                  t_ca_elapsed,t_fact_elapsed,t_solve_elapsed)
      myid = OMP_GET_THREAD_NUM()
      min_local = HUGE(min_local)
      max_local = -HUGE(max_local)

      ALLOCATE (IPIV(N), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc IPIV'
      ALLOCATE (S(N, 3), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc S'
      ALLOCATE (A(LDAB, N), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc A'

!$OMP DO SCHEDULE(static)
      DO IK = i1, i2
         CALL SYSTEM_CLOCK(t_start, rate)
         FK = FKY(IK)

         ! ===== CA (Assembly) timing =====
         CALL SYSTEM_CLOCK(t_ca_start)
         CALL CAGB(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                   CI, NORD, AS, WT, IE0, IS0, DZ0, N, KL, KU, LDAB, IG, A)

         CALL SYSTEM_CLOCK(t_ca_end)
         t_ca_elapsed = REAL(t_ca_end - t_ca_start, dp)/REAL(rate, dp)

         ! ===== Factorization timing =====
         CALL SYSTEM_CLOCK(t_fact_start)
         CALL CLUD(N, KL, KU, LDAB, IPIV, INFO, A)
         CALL SYSTEM_CLOCK(t_fact_end)
         t_fact_elapsed = REAL(t_fact_end - t_fact_start, dp)/REAL(rate, dp)

         ! Basis tensor (currently identity for 3 components)
         SV(1:3, 1:3) = 0.0_dp
         SV(1, 1) = 1.0_dp
         SV(2, 2) = 1.0_dp
         SV(3, 3) = 1.0_dp

         ! ===== Solve timing (for all sources) =====
         CALL SYSTEM_CLOCK(t_solve_start)
         DO IS = 1, NSS
            ! Zero all 3 RHS at once
            S(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)

            ! Assemble all three RHS columns for this source
            DO II = 1, MSR(IS)
               idx = MSR1(IS, II)   ! global node index
               IF (dbg) THEN
                  IF (idx < 1 .OR. idx > NPT) THEN
                     WRITE (*, *) 'run_LU_hybrid: idx out of range in RHS assembly: IS=', IS, ' II=', II, ' idx=', idx, ' NPT=', NPT, ' IK=', IK, ' rank=', my_rank
                     CALL MPI_Abort(comm, 1, INFO)
                  END IF
               END IF

               ! ICS = 1
               S(3*idx - 2, 1) = S(3*idx - 2, 1) + CMPLX(FSR(IS, II)*SV(1, 1), 0.0_dp, dp)
               S(3*idx - 1, 1) = S(3*idx - 1, 1) + CMPLX(FSR(IS, II)*SV(1, 2), 0.0_dp, dp)
               S(3*idx, 1) = S(3*idx, 1) + CMPLX(FSR(IS, II)*SV(1, 3), 0.0_dp, dp)

               ! ICS = 2
               S(3*idx - 2, 2) = S(3*idx - 2, 2) + CMPLX(FSR(IS, II)*SV(2, 1), 0.0_dp, dp)
               S(3*idx - 1, 2) = S(3*idx - 1, 2) + CMPLX(FSR(IS, II)*SV(2, 2), 0.0_dp, dp)
               S(3*idx, 2) = S(3*idx, 2) + CMPLX(FSR(IS, II)*SV(2, 3), 0.0_dp, dp)

               ! ICS = 3
               S(3*idx - 2, 3) = S(3*idx - 2, 3) + CMPLX(FSR(IS, II)*SV(3, 1), 0.0_dp, dp)
               S(3*idx - 1, 3) = S(3*idx - 1, 3) + CMPLX(FSR(IS, II)*SV(3, 2), 0.0_dp, dp)
               S(3*idx, 3) = S(3*idx, 3) + CMPLX(FSR(IS, II)*SV(3, 3), 0.0_dp, dp)
            END DO

            ! One solve with 3 RHS
             CALL SOLVER('N', N, KL, KU, 3, LDAB, IPIV, S, N, INFO, A)
            ! CALL SOLVER_LEG('No transpose',N,KL,KU,3,LDAB,IPIV,S,N,INFO,A)

            ! Scatter solution back to GFX/GFY/GFZ for all components
            DO II = 1, NPT
               ! ICS = 1
               GFX(IS, 1, IK, II) = S(3*II - 2, 1)
               GFY(IS, 1, IK, II) = S(3*II - 1, 1)
               GFZ(IS, 1, IK, II) = S(3*II, 1)
               ! ICS = 2
               GFX(IS, 2, IK, II) = S(3*II - 2, 2)
               GFY(IS, 2, IK, II) = S(3*II - 1, 2)
               GFZ(IS, 2, IK, II) = S(3*II, 2)
               ! ICS = 3
               GFX(IS, 3, IK, II) = S(3*II - 2, 3)
               GFY(IS, 3, IK, II) = S(3*II - 1, 3)
               GFZ(IS, 3, IK, II) = S(3*II, 3)
               IF (dbg) THEN
                  min_local = MIN(min_local, II)
                  max_local = MAX(max_local, II)
               END IF
            END DO
         END DO   ! IS
         CALL SYSTEM_CLOCK(t_solve_end)
         t_solve_elapsed = REAL(t_solve_end - t_solve_start, dp)/REAL(rate, dp)

         CALL SYSTEM_CLOCK(t_end)

         t_elapsed = REAL(t_end - t_start, dp)/REAL(rate, dp)
         IF (ITER <= 1 .AND. my_rank == 0) THEN
            WRITE (*, '(A,F7.1,A,E11.4,A,I4,A,I4,A,F12.6,A,F8.4,A,F8.4,A,F8.4,A)') &
               'f=', FREQ, ' Hz | ky=', FK, ' | IK=', IK, '/', NK, ' | t=', t_elapsed, &
               ' s | CA=', t_ca_elapsed, ' | fact=', t_fact_elapsed, ' | solve=', t_solve_elapsed
         END IF
      END DO
!$OMP END DO

!$OMP CRITICAL
      min_idx_used = MIN(min_idx_used, min_local)
      max_idx_used = MAX(max_idx_used, max_local)
!$OMP END CRITICAL

      DEALLOCATE (A, S, IPIV)
!$OMP END PARALLEL

      CALL omp_set_num_threads(1)

      IF (dbg .AND. my_rank == 0) THEN
         IF (SIZE(GFX, 1) /= NSS .OR. SIZE(GFX, 2) /= 3 .OR. SIZE(GFX, 3) /= NK .OR. SIZE(GFX, 4) /= NPT) THEN
            WRITE (*, *) 'run_LU_hybrid: GFX shape mismatch: ', SHAPE(GFX), ' expected=', (/NSS, 3, NK, NPT/)
            CALL MPI_Abort(comm, 1, INFO)
         END IF
         WRITE (*, *) 'run_LU_hybrid: idx range used = ', min_idx_used, ' to ', max_idx_used, ' of NPT=', NPT
      END IF

   END SUBROUTINE run_LU_hybrid

   SUBROUTINE run_LU_mkl_NSS(N, KL, KU, LDAB, i1, i2, &
                             I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                             AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                             IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCOREI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, KL, KU, LDAB, i1, i2, I25D, IFQ, NTO, NX, NZ, NNX, NNZ, NPT, IANISO, NORD
      INTEGER, INTENT(IN) :: NK, NSR, NSS, IE0, IS0, my_rank, comm, ITER, NCOREI
      INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FREQ, FKY(:)
      REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
      REAL(dp), INTENT(IN) :: DZ0
      TYPE(InversionGridType), INTENT(IN) :: IG

      COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE :: A(:, :), S(:, :)
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(dp)        :: FK
      REAL(dp)        :: SV(3, 3)
      INTEGER :: IS, IC, II, IK, INFO, ierr_alloc
      INTEGER :: locK, idx
      INTEGER :: rate, t_start, t_end
      INTEGER :: t_ca_start, t_ca_end, t_fact_start, t_fact_end, t_solve_start, t_solve_end
      REAL(dp) :: t_elapsed, t_ca_elapsed, t_fact_elapsed, t_solve_elapsed
      INTEGER :: ncore_eff
      LOGICAL :: use_omp
      CHARACTER(len=128) :: fname_ser, fname_par
      ! effective OMP threads over sources
      ncore_eff = MAX(1, MIN(NCOREI, NSS))
      use_omp = (ncore_eff > 1 .AND. NSS > 1)

      SV = 0.0_dp
      SV(1, 1) = 1.0_dp
      SV(2, 2) = 1.0_dp
      SV(3, 3) = 1.0_dp

      IF (my_rank == 0 .AND. ITER <= 1) THEN
         IF (use_omp) THEN
            WRITE (*, *) 'run_LU_mkl_NSS: OMP over sources, NCORE =', ncore_eff
         ELSE
            WRITE (*, *) 'run_LU_mkl_NSS: serial over sources (NSS=', NSS, ')'
         END IF
      END IF

      ALLOCATE (IPIV(N), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc IPIV'
      ALLOCATE (A(LDAB, N), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc A'

      DO IK = i1, i2
         CALL SYSTEM_CLOCK(t_start, rate)
         FK = FKY(IK)

         ! ===== CA (Assembly) timing =====
         !CALL mkl_set_dynamic(0)
         CALL mkl_set_num_threads(1)
         CALL SYSTEM_CLOCK(t_ca_start)
         CALL CAGB(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
                   CI, NORD, AS, WT, IE0, IS0, DZ0, N, KL, KU, LDAB, IG, A)
         CALL SYSTEM_CLOCK(t_ca_end)
         t_ca_elapsed = REAL(t_ca_end - t_ca_start, dp)/REAL(rate, dp)

         ! ===== Factorization timing =====
         !CALL mkl_set_dynamic(0)
         CALL mkl_set_num_threads(32)
         CALL SYSTEM_CLOCK(t_fact_start)
         CALL CLUD(N, KL, KU, LDAB, IPIV, INFO, A)
         CALL SYSTEM_CLOCK(t_fact_end)
         t_fact_elapsed = REAL(t_fact_end - t_fact_start, dp)/REAL(rate, dp)

         ! ===== Solve timing: MKL single-threaded; optional OMP over sources =====
         CALL mkl_set_num_threads(1)
         CALL omp_set_num_threads(NCOREI)

         CALL SYSTEM_CLOCK(t_solve_start)
!$OMP PARALLEL IF(use_omp) DEFAULT(SHARED) PRIVATE(IS,IC,II,locK,idx,S) FIRSTPRIVATE(SV) NUM_THREADS(ncore_eff)
         ALLOCATE (S(N, 3), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc S'

!$OMP DO SCHEDULE(static)
         DO IS = 1, NSS
            ! basis tensor (identity)

            ! zero all three RHS
            S(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)

            ! build three RHS columns
            DO IC = 1, 3
               DO locK = 1, MSR(IS)
                  idx = MSR1(IS, locK)
                  S(3*idx - 2, IC) = S(3*idx - 2, IC) + CMPLX(FSR(IS, locK)*SV(IC, 1), 0.0_dp, dp)
                  S(3*idx - 1, IC) = S(3*idx - 1, IC) + CMPLX(FSR(IS, locK)*SV(IC, 2), 0.0_dp, dp)
                  S(3*idx, IC) = S(3*idx, IC) + CMPLX(FSR(IS, locK)*SV(IC, 3), 0.0_dp, dp)
               END DO
            END DO

            ! solve 3 RHS at once
            CALL SOLVER('N', N, KL, KU, 3, LDAB, IPIV, S, N, INFO, A)

            ! scatter
            DO II = 1, NPT
               DO IC = 1, 3
                  GFX(IS, IC, IK, II) = S(3*II - 2, IC)
                  GFY(IS, IC, IK, II) = S(3*II - 1, IC)
                  GFZ(IS, IC, IK, II) = S(3*II, IC)
               END DO
            END DO
         END DO
!$OMP END DO

         DEALLOCATE (S)
!$OMP END PARALLEL
         CALL SYSTEM_CLOCK(t_solve_end)
         t_solve_elapsed = REAL(t_solve_end - t_solve_start, dp)/REAL(rate, dp)

         CALL SYSTEM_CLOCK(t_end)
         t_elapsed = REAL(t_end - t_start, dp)/REAL(rate, dp)
         IF (ITER <= 1 .AND. my_rank == 0) THEN
            WRITE (*, '(A,F7.1,A,E11.4,A,I4,A,I4,A,F12.6,A,F8.4,A,F8.4,A,F8.4,A)') &
               'f=', FREQ, ' Hz | ky=', FK, ' | IK=', IK, '/', NK, ' | t=', t_elapsed, &
               ' s | CA=', t_ca_elapsed, ' | fact=', t_fact_elapsed, ' | solve=', t_solve_elapsed
         END IF
      END DO

      DEALLOCATE (A, IPIV)
      CALL omp_set_num_threads(1)

   END SUBROUTINE run_LU_mkl_NSS

   SUBROUTINE run_ik_loop_pardiso(N, KL, KU, LDAB, i1, i2, &
                                  I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                                  AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                                  IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, NBLOCK, &
                                  comm, ITER)
      IMPLICIT NONE
      ! ----------------- inputs -----------------
      INTEGER, INTENT(IN) :: N, KL, KU, LDAB, i1, i2, I25D, IFQ, NTO, NX, NZ, NNX, NNZ, NPT, IANISO, NORD
      INTEGER, INTENT(IN) :: NK, NSR, NSS, IE0, IS0, my_rank, comm, ITER, NBLOCK
      INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
      REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FREQ, FKY(:)
      REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
      REAL(dp), INTENT(IN) :: DZ0
      TYPE(InversionGridType), INTENT(IN) :: IG

      COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)

      ! ----------------- locals -----------------
      REAL(dp) :: SV(3, 3)
      INTEGER :: IK, IS, IC, II, idx
      REAL(dp) :: FK
      INTEGER :: rate, t_start, t_end, t_cas, t_cae
      INTEGER :: t_fact_start, t_fact_end, t_solve_start, t_solve_end
      REAL(dp) :: t_elapsed, t_ca, t_fact_pardiso, t_solve_pardiso
      REAL(dp) :: mem_sym_gb, mem_perm_gb, mem_num_gb, mem_peak_gb
      ! fine-grained timers
      INTEGER :: t_rhs_start, t_rhs_end, t_solve_phase_start, t_solve_phase_end, t_scatter_start, t_scatter_end
      INTEGER :: t_batches_start, t_batches_end
      REAL(dp) :: t_rhs_total, t_solve_phase_total, t_scatter_total, t_batches_total
      ! CSR storage
      INTEGER, ALLOCATABLE :: ia(:), iap(:), ja(:), jap(:)
      COMPLEX(dp), ALLOCATABLE :: A(:), AP(:)
      ! PARDISO controls
      INTEGER(KIND=8) :: pt(64)
      INTEGER         :: iparm(64)
      INTEGER         :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
      INTEGER, ALLOCATABLE :: perm(:)
      COMPLEX(dp), ALLOCATABLE :: rhs(:, :), sol(:, :)
      COMPLEX(dp) :: rhs_dummy(1, 1), sol_dummy(1, 1)
      ! Column map (RHS index -> (IS,IC))
      INTEGER, ALLOCATABLE :: col2is(:), col2ic(:)
      INTEGER :: nrhs_tot, batch_beg, batch_end, batch_size, jcol
      INTEGER, PARAMETER :: MAX_RHS_BATCH = 32
      LOGICAL :: first_ik
      LOGICAL :: did_factorization
      INTEGER :: max_msr1_cols, ierr_abort
      
      INTEGER :: min_idx_used, max_idx_used
      INTEGER :: nnz_par, nnz_ser
      CHARACTER(LEN=256) :: fname_ser, fname_par
      LOGICAL :: dbg
      dbg = (ITER <= 1)

      IF (dbg .AND. my_rank == 0) WRITE (*, *) 'running run_ik_loop_pardiso'
      !CALL mkl_set_dynamic(0)
      CALL mkl_set_num_threads(1)

      ! ---- PARDISO setup  ----
      mtype = 13          ! complex unsymmetric
      maxfct = 1
      mnum = 1
      msglvl = 0
      nrhs = 1
      error = 1
      pt = 0_8
      iparm = 0

      iparm(1) = 1
      iparm(2) = 2
      iparm(3) = 0
      iparm(4) = 0
      iparm(5) = 0
      iparm(9) = 13
      iparm(11) = 0
      iparm(15) = 0      ! peak memory symbolic (output)
      iparm(16) = 0      ! permanent memory symbolic (output)
      iparm(17) = 0      ! memory factor+solve (output)
      iparm(18) = -1     ! report nnz in factors
      iparm(27) = 0
      iparm(34) = 0

      rhs_dummy = CMPLX(0.0_dp, 0.0_dp, dp)
      sol_dummy = CMPLX(0.0_dp, 0.0_dp, dp)

      ! permutation array reused for all IK
      ALLOCATE (perm(N))
      perm = 0

      first_ik = .TRUE.
      did_factorization = .FALSE.
      mem_sym_gb = 0.0_dp
      mem_perm_gb = 0.0_dp
      mem_num_gb = 0.0_dp
      mem_peak_gb = 0.0_dp

      CALL SYSTEM_CLOCK(count_rate=rate)

      max_msr1_cols = SIZE(MSR1, 2)

      SV = 0.0_dp
      SV(1, 1) = 1.0_dp
      SV(2, 2) = 1.0_dp
      SV(3, 3) = 1.0_dp

      ! Fixed-size batch buffers: These stay alive for all IK and all batches.
      nrhs_tot = 3*NSS
      IF (.NOT. ALLOCATED(col2is)) ALLOCATE (col2is(nrhs_tot))
      IF (.NOT. ALLOCATED(col2ic)) ALLOCATE (col2ic(nrhs_tot))

      IF (.NOT. ALLOCATED(rhs)) ALLOCATE (rhs(N, MAX_RHS_BATCH))
      IF (.NOT. ALLOCATED(sol)) ALLOCATE (sol(N, MAX_RHS_BATCH))

      ! Pre-calculate the column map
      jcol = 0
      DO IS = 1, NSS
         DO IC = 1, 3
            jcol = jcol + 1
            col2is(jcol) = IS
            col2ic(jcol) = IC
         END DO
      END DO

      ! Build column map: (IS, IC=1..3) for all sources

      DO IK = i1, i2
         CALL SYSTEM_CLOCK(t_start)
         FK = FKY(IK)
         t_rhs_total = 0.0_dp
         t_solve_phase_total = 0.0_dp
         t_scatter_total = 0.0_dp

! --- run parallel assembler ---

         CALL omp_set_num_threads(MIN(32, NTHREAD_MAX))
         CALL mkl_set_num_threads(1)
         CALL SYSTEM_CLOCK(t_cas)
         if (my_rank == 0 .and. ik == 1 .and.dbg) call log_rss("before CA", 6, my_rank)
         CALL CA_CSR_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                         AS, WT, IE0, IS0, DZ0, N, IA, JA, A, NBLOCK, LDAB, KL, KU, IG, my_rank, DEBUG_OUTPUT=.FALSE.)

         CALL SYSTEM_CLOCK(t_cae)
         t_ca = REAL(t_cae - t_cas, dp)/REAL(rate, dp)

         if (my_rank == 0 .and. ik == 1.and.dbg) call log_rss("after CA", 6, my_rank)

         CALL omp_set_num_threads(1)
         CALL mkl_set_num_threads(MIN(32, NTHREAD_MAX))
         ! ---- analysis + factorization on first IK, then factorization only ----
         IF (first_ik) THEN
            ! IA=IAP
            ! JA=JAP
            ! A=AP
            phase = 12       ! analysis + factorization
         ELSE
            ! CALL COMPARE_CSR(N, IAP,JAP,AP, IA,JA,A)
            phase = 22       ! factorization only (reuse analysis & perm)
         END IF

         ! ===== PARDISO Factorization timing (phase 12 or 22) =====
         CALL SYSTEM_CLOCK(t_fact_start)

         CALL pardiso(pt, maxfct, mnum, mtype, phase, N, A, ia, ja, perm, 1, &
                      iparm, msglvl, rhs_dummy, sol_dummy, error)
         if (my_rank == 0 .and. ik == 1 .and.dbg) call log_rss("after factor", 6, my_rank)

         CALL SYSTEM_CLOCK(t_fact_end)
         t_fact_pardiso = REAL(t_fact_end - t_fact_start, dp)/REAL(rate, dp)

         IF (error /= 0) THEN
            IF (my_rank == 0) WRITE (*, *) 'PARDISO error (phase ', phase, ') =', error
            STOP
         END IF

         did_factorization = .TRUE.
         mem_sym_gb = REAL(iparm(15), dp)/(1024.0_dp*1024.0_dp)
         mem_perm_gb = REAL(iparm(16), dp)/(1024.0_dp*1024.0_dp)
         mem_num_gb = REAL(iparm(17), dp)/(1024.0_dp*1024.0_dp)
         mem_peak_gb = MAX(mem_sym_gb, mem_perm_gb + mem_num_gb)

         first_ik = .FALSE.
         ! IF (ITER <= 1 .AND. my_rank == 0 .AND. did_factorization) THEN
         !    WRITE (*, *) 'PARDISO mem: f=', FREQ, &
         !       ' mem_sym=', mem_sym_gb, ' GB', &
         !       ' mem_perm=', mem_perm_gb, ' GB', &
         !       ' mem_num=', mem_num_gb, ' GB', &
         !       ' mem_peak=', mem_peak_gb, ' GB'
         ! END IF

         ! ---- batched solves ----
         batch_beg = 1

         CALL SYSTEM_CLOCK(t_batches_start)

         DO WHILE (batch_beg <= nrhs_tot)
            batch_end = MIN(nrhs_tot, batch_beg + MAX_RHS_BATCH - 1)
            batch_size = batch_end - batch_beg + 1

            rhs = CMPLX(0.0_dp, 0.0_dp, dp)

            ! ===== PARDISO Solve timing (phase 33) =====

            CALL SYSTEM_CLOCK(t_solve_start)

            ! assemble RHS for each (IS,IC) in this batch using SV (identity)

            CALL SYSTEM_CLOCK(t_rhs_start)
!               !$OMP PARALLEL DO PRIVATE(jcol, IS, IC, II, idx) SCHEDULE(static)

            DO jcol = 1, batch_size
               IS = col2is(batch_beg + jcol - 1)
               IC = col2ic(batch_beg + jcol - 1)

               DO II = 1, MSR(IS)
                  idx = MSR1(IS, II)            ! 1..NPT
                  rhs(3*idx - 2, jcol) = rhs(3*idx - 2, jcol) + CMPLX(FSR(IS, II)*SV(IC, 1), 0.0_dp, dp)
                  rhs(3*idx - 1, jcol) = rhs(3*idx - 1, jcol) + CMPLX(FSR(IS, II)*SV(IC, 2), 0.0_dp, dp)
                  rhs(3*idx, jcol) = rhs(3*idx, jcol) + CMPLX(FSR(IS, II)*SV(IC, 3), 0.0_dp, dp)
               END DO
            END DO
!            !$OMP END PARALLEL DO

            CALL SYSTEM_CLOCK(t_rhs_end)
            t_rhs_total = t_rhs_total + REAL(t_rhs_end - t_rhs_start, dp)/REAL(rate, dp)

            !-- solve---

            phase = 33
            CALL SYSTEM_CLOCK(t_solve_phase_start)

            CALL pardiso(pt, maxfct, mnum, mtype, 33, N, A, ia, ja, perm, batch_size, &
                         iparm, msglvl, rhs(1, 1), sol(1, 1), error)
            CALL SYSTEM_CLOCK(t_solve_phase_end)

            t_solve_phase_total = t_solve_phase_total + REAL(t_solve_phase_end - t_solve_phase_start, dp)/REAL(rate, dp)
            IF (error /= 0) THEN
               IF (my_rank == 0) WRITE (*, *) 'PARDISO error (phase 33) =', error
               STOP
            END IF
            CALL mkl_set_num_threads(1)
            CALL omp_set_num_threads(MIN(32, NTHREAD_MAX))
            ! scatter solutions back GFX to (IS,IC,IK,:)
            CALL SYSTEM_CLOCK(t_scatter_start)

            !$OMP PARALLEL DO PRIVATE(II, jcol, IS, IC) SCHEDULE(static)

            DO jcol = 1, batch_size
               IS = col2is(batch_beg + jcol - 1)
               IC = col2ic(batch_beg + jcol - 1)
               DO II = 1, NPT

                  GFX(IS, IC, IK, II) = sol(3*II - 2, jcol)
                  GFY(IS, IC, IK, II) = sol(3*II - 1, jcol)
                  GFZ(IS, IC, IK, II) = sol(3*II, jcol)
               END DO
            END DO
            !$OMP END PARALLEL DO
            CALL omp_set_num_threads(1)
            CALL mkl_set_num_threads(MIN(32, NTHREAD_MAX))
            CALL SYSTEM_CLOCK(t_scatter_end)

            t_scatter_total = t_scatter_total + REAL(t_scatter_end - t_scatter_start, dp)/REAL(rate, dp)

            batch_beg = batch_end + 1
            ! batch_beg = batch_beg + batch_size
         END DO

         CALL SYSTEM_CLOCK(t_batches_end)     
         t_batches_total = REAL(t_batches_end - t_batches_start, dp)/REAL(rate, dp)

         CALL SYSTEM_CLOCK(t_solve_end)
         t_solve_pardiso = REAL(t_solve_end - t_solve_start, dp)/REAL(rate, dp)

         CALL SYSTEM_CLOCK(t_end)
         t_elapsed = REAL(t_end - t_start, dp)/REAL(rate, dp)

         IF (ITER <= 1 .AND. my_rank == 0) THEN
            WRITE (*, '(A,F10.3,A,E14.6,A,I8,A,I8,A,F12.6,A,  A,F8.4,A,F8.4,A,F8.4,A,F8.4,A,F8.4)') &
               'f=', FREQ, ' Hz | ky=', FK, ' | IK=', IK, '/', NK, ' | t=', t_elapsed, ' s', &
               ' | CA=', t_ca, ' | fact=', t_fact_pardiso, ' | GF=', t_batches_total, &
               ' | solve=', t_solve_phase_total, ' | scatter=', t_scatter_total
            CALL flush (6)
         END IF
         call mkl_set_num_threads(1)
      END DO   ! IK loop

      IF (ITER <= 1 .AND. my_rank == 0 .AND. did_factorization) THEN
         WRITE (*, '(A,F6.2,A,  A,F6.2,A,  A,F6.2,A,  A,F6.2,A)') &
            'PARDISO mem (f=', FREQ, ' Hz): ', &
            'sym=', mem_sym_gb, ' GB  ', &
            'perm=', mem_perm_gb, ' GB  ', &
            'num=', mem_num_gb, ' GB  ', &
            'peak=', mem_peak_gb, ' GB'
      END IF

      ! ---- release PARDISO after all IK are done ----
      IF (did_factorization) THEN
         phase = -1
         nrhs = 0
         CALL pardiso(pt, maxfct, mnum, mtype, phase, 0, A, ia, ja, perm, nrhs, &
                      iparm, msglvl, rhs_dummy, sol_dummy, error)
         IF (ALLOCATED(ia)) DEALLOCATE (ia)
         IF (ALLOCATED(ja)) DEALLOCATE (ja)
         IF (ALLOCATED(A)) DEALLOCATE (A)
      END IF
      IF (ALLOCATED(perm)) DEALLOCATE (perm)
      call mkl_set_num_threads(1)

      IF (dbg .AND. my_rank == 0) THEN
         IF (SIZE(GFX, 1) /= NSS .OR. SIZE(GFX, 2) /= 3 .OR. SIZE(GFX, 3) /= NK .OR. SIZE(GFX, 4) /= NPT) THEN
            WRITE (*, *) 'run_ik_loop_pardiso: GFX shape mismatch: ', SHAPE(GFX), ' expected=', (/NSS, 3, NK, NPT/)
            CALL MPI_Abort(comm, 1, error)
         END IF
      END IF
   END SUBROUTINE run_ik_loop_pardiso



#if HAVE_MUMPS
   ! MUMPS headers, derived types and calls are compiled only when enabled.

SUBROUTINE run_ik_loop_mumps(N,LDAB, KL, KU, i1, i2, &
                            I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
                            AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR, IG, &
                            IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, NBLOCK, comm, ITER)
            
   IMPLICIT NONE


   ! ---- inputs ----
   INTEGER, INTENT(IN) :: N, i1, i2, I25D, IFQ, NTO, NX, NZ, NNX, NNZ, NPT, IANISO, NORD
   INTEGER, INTENT(IN) :: NK, NSR, NSS, IE0, IS0, my_rank, comm, ITER, NBLOCK,LDAB, KL, KU
   INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
   REAL(dp), INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FREQ, FKY(:)
   REAL(dp), INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
   REAL(dp), INTENT(IN) :: DZ0
   TYPE(InversionGridType), INTENT(IN) :: IG
   COMPLEX(sp), CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)

   ! ---- locals ----

   INTEGER :: IK, IS, IC, II, idx
   REAL(dp) :: FK
   INTEGER :: ierr, rate
   INTEGER :: t_start, t_end, t_cas, t_cae, t_fact_s, t_fact_e
   REAL(dp) :: t_elapsed, t_ca, t_fact, t_solve_phase_total, t_scatter_total
   INTEGER, PARAMETER :: MAX_RHS_BATCH = 32
   INTEGER :: nrhs_tot, batch_beg, batch_end, batch_size, jcol
   INTEGER, ALLOCATABLE :: col2is(:), col2ic(:)
   REAL(dp) :: mem_fact_gb, mem_work_gb, mem_tot_gb, mem_peak_gb

   ! ---- COO from assembler ----
   INTEGER, ALLOCATABLE, TARGET :: IA(:), JA(:)
   COMPLEX(dp), ALLOCATABLE, TARGET :: A(:)
   INTEGER :: nnz_out


   ! ---- RHS batching ----
   COMPLEX(dp), ALLOCATABLE, TARGET :: rhs(:, :)
   COMPLEX(dp), POINTER :: rhs_flat(:)
   REAL(dp) :: SV(3,3)

   LOGICAL :: first_ik
   LOGICAL :: dbg
!----MUMPS setup ----
   include 'zmumps_struc.h'
   TYPE(ZMUMPS_STRUC) :: mumps_par 
   ! =========================
   !  MUMPS init (per rank)
   ! =========================
   mumps_par%COMM = MPI_COMM_SELF
   mumps_par%PAR  = 1
   mumps_par%SYM  = 0
   mumps_par%JOB  = -1
   CALL zmumps_call_safe(mumps_par)

   ! Silence
   mumps_par%ICNTL(1) = -1
   mumps_par%ICNTL(2) = -1
   mumps_par%ICNTL(3) = -1
   mumps_par%ICNTL(4) = 2
      dbg = (ITER <= 1)
   CALL SYSTEM_CLOCK(count_rate=rate)
   IF (rate <= 0) rate = 1


   SV = 0.0_dp
   SV(1,1)=1.0_dp; SV(2,2)=1.0_dp; SV(3,3)=1.0_dp

   !prebuild NRHS table
   nrhs_tot = 3*NSS
   ALLOCATE(col2is(nrhs_tot), col2ic(nrhs_tot))
   jcol = 0
   DO IS = 1, NSS
      DO IC = 1, 3
         jcol = jcol + 1
         col2is(jcol) = IS
         col2ic(jcol) = IC
      END DO
   END DO

   ALLOCATE(rhs(N, MAX_RHS_BATCH))
   rhs_flat(1:N*MAX_RHS_BATCH) => rhs


   ! Assembled matrix input
   ! (Most builds use ICNTL(5)=0 for assembled; if your headers differ, adjust.)
   mumps_par%ICNTL(5) = 0

   first_ik = .TRUE.

   DO IK = i1, i2
      CALL SYSTEM_CLOCK(t_start)
      FK = FKY(IK)

      ! ---- assemble COO (parallel) ----
      CALL SYSTEM_CLOCK(t_cas)
         CALL omp_set_num_threads(MIN(32, NTHREAD_MAX))
         CALL mkl_set_num_threads(1)
         
         if (my_rank == 0 .and. ik == 1 .and.dbg) call log_rss("before CA", 6, my_rank)
         CALL CA_COO_par(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, CI, NORD, &
                         AS, WT, IE0, IS0, DZ0, N, IA, JA, A, NBLOCK, LDAB, KL, KU,NNZ_OUT, IG, my_rank)


      CALL SYSTEM_CLOCK(t_cae)
      t_ca = REAL(t_cae - t_cas, dp)/REAL(rate, dp)
         if (my_rank == 0 .and. ik == 1.and.dbg) call log_rss("after CA", 6, my_rank)

CALL omp_set_num_threads(1)
      ! ---- attach matrix to MUMPS ----
      mumps_par%N  = N
      mumps_par%NZ = nnz_out

      mumps_par%IRN => IA
      mumps_par%JCN => JA
      mumps_par%A   => A

      CALL validate_mumps_coo_inputs(N, nnz_out, IA, JA, A, comm, my_rank, IK)

         
         CALL mkl_set_num_threads(MIN(32, NTHREAD_MAX))
      ! =========================
      !  Analysis once (IK=i1)
      ! =========================
      IF (first_ik) THEN
         mumps_par%JOB = 1
         CALL zmumps_call_safe(mumps_par)
         IF (mumps_par%INFOG(1) /= 0) CALL abort_mumps('MUMPS analysis failed', mumps_par, comm)

         first_ik = .FALSE.
      END IF

      ! =========================
      !  Factorization (every IK)
      ! =========================
      CALL SYSTEM_CLOCK(t_fact_s)
      mumps_par%JOB = 2
      CALL zmumps_call_safe(mumps_par)
      CALL SYSTEM_CLOCK(t_fact_e)
      t_fact = REAL(t_fact_e - t_fact_s, dp)/REAL(rate, dp)

      IF (mumps_par%INFOG(1) /= 0) CALL abort_mumps('MUMPS factor failed', mumps_par, comm)


mem_fact_gb = REAL(mumps_par%INFO(16), dp) / 1024.0_dp
mem_work_gb = REAL(mumps_par%INFO(17), dp) / 1024.0_dp
mem_tot_gb  = REAL(mumps_par%INFO(18), dp) / 1024.0_dp
mem_peak_gb = REAL(mumps_par%INFO(22), dp) / 1024.0_dp

IF (my_rank == 0) THEN
   WRITE(*,'(A,F8.3)') 'MUMPS factor memory (GB): ', mem_fact_gb
   WRITE(*,'(A,F8.3)') 'MUMPS work   memory (GB): ', mem_work_gb
   WRITE(*,'(A,F8.3)') 'MUMPS total  memory (GB): ', mem_tot_gb
   WRITE(*,'(A,F8.3)') 'MUMPS peak   memory (GB): ', mem_peak_gb
END IF
      ! =========================
      !  Batched solves (JOB=3)
      ! =========================
      batch_beg = 1
      t_solve_phase_total = 0.0_dp
      t_scatter_total     = 0.0_dp

      DO WHILE (batch_beg <= nrhs_tot)
         batch_end  = MIN(nrhs_tot, batch_beg + MAX_RHS_BATCH - 1)
         batch_size = batch_end - batch_beg + 1

         rhs(:,:) = CMPLX(0.0_dp, 0.0_dp, dp)

         ! build RHS columns (same as yours)
         DO jcol = 1, batch_size
            IS = col2is(batch_beg + jcol - 1)
            IC = col2ic(batch_beg + jcol - 1)
            DO II = 1, MSR(IS)
               idx = MSR1(IS, II)
               rhs(3*idx - 2, jcol) = rhs(3*idx - 2, jcol) + CMPLX(FSR(IS, II)*SV(IC,1), 0.0_dp, dp)
               rhs(3*idx - 1, jcol) = rhs(3*idx - 1, jcol) + CMPLX(FSR(IS, II)*SV(IC,2), 0.0_dp, dp)
               rhs(3*idx,     jcol) = rhs(3*idx,     jcol) + CMPLX(FSR(IS, II)*SV(IC,3), 0.0_dp, dp)
            END DO
         END DO

         ! attach RHS to MUMPS
         mumps_par%NRHS = batch_size
         mumps_par%LRHS = N
         mumps_par%RHS  => rhs_flat

         mumps_par%JOB = 3
         CALL zmumps_call_safe(mumps_par)
         IF (mumps_par%INFOG(1) /= 0) CALL abort_mumps('MUMPS solve failed', mumps_par, comm)

         ! scatter back (same mapping)
         CALL mkl_set_num_threads(1)
         CALL omp_set_num_threads(MIN(32, NTHREAD_MAX))

         !$OMP PARALLEL DO PRIVATE(II, jcol, IS, IC) SCHEDULE(static)
         DO jcol = 1, batch_size
            IS = col2is(batch_beg + jcol - 1)
            IC = col2ic(batch_beg + jcol - 1)
            DO II = 1, NPT
               GFX(IS,IC,IK,II) = rhs(3*II - 2, jcol)
               GFY(IS,IC,IK,II) = rhs(3*II - 1, jcol)
               GFZ(IS,IC,IK,II) = rhs(3*II,     jcol)
            END DO
         END DO
         !$OMP END PARALLEL DO
            CALL omp_set_num_threads(1)
            CALL mkl_set_num_threads(MIN(32, NTHREAD_MAX))
         batch_beg = batch_end + 1
      END DO

      CALL SYSTEM_CLOCK(t_end)
      t_elapsed = REAL(t_end - t_start, dp)/REAL(rate, dp)

      IF (ITER <= 1 .AND. my_rank == 0) THEN
         WRITE(*,'(A,F10.3,A,E14.6,A,I8,A,I8,A,F12.6,A, A,F8.4,A,F8.4)') &
            'f=',FREQ,' Hz | ky=',FK,' | IK=',IK,'/',NK,' | t=',t_elapsed,' s', &
            ' | CA=',t_ca,' | fact=',t_fact
         CALL flush(6)
            call mkl_set_num_threads(1)
      END IF

      ! free per-IK COO buffers
      IF (ALLOCATED(IA))   DEALLOCATE(IA)
      IF (ALLOCATED(JA))   DEALLOCATE(JA)
      IF (ALLOCATED(A))  DEALLOCATE(A)
   END DO

   ! =========================
   !  Finalize MUMPS
   ! =========================
   mumps_par%JOB = -2
   CALL zmumps_call_safe(mumps_par)

   DEALLOCATE(rhs, col2is, col2ic)

END SUBROUTINE run_ik_loop_mumps

   SUBROUTINE abort_mumps(msg, id, comm)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: msg
      INTEGER, INTENT(IN) :: comm
      INCLUDE 'zmumps_struc.h'
      TYPE(ZMUMPS_STRUC), INTENT(IN) :: id
      INTEGER :: ierr_abort, my_rank_local, errcode

      CALL MPI_Comm_rank(comm, my_rank_local, ierr_abort)
      IF (my_rank_local == 0) THEN
         WRITE (*, '(A,2X,A,2X,A,I0,2X,A,I0)') 'MUMPS_ABORT:', TRIM(msg), &
            'INFOG(1)=', id%INFOG(1), 'INFOG(2)=', id%INFOG(2)
      END IF

      errcode = id%INFOG(1)
      IF (errcode == 0) errcode = 1
      CALL MPI_Abort(comm, errcode, ierr_abort)
   END SUBROUTINE abort_mumps

   SUBROUTINE validate_mumps_coo_inputs(N, nz, irn, jcn, aval, comm, my_rank, ik)
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, nz, comm, my_rank, ik
      INTEGER, INTENT(IN) :: irn(:), jcn(:)
      COMPLEX(dp), INTENT(IN) :: aval(:)
      INTEGER :: k, ierr_abort
      REAL(dp) :: ar, ai

      IF (nz < 0) THEN
         IF (my_rank == 0) WRITE (*, *) 'MUMPS input check: negative NZ at IK=', ik, ' NZ=', nz
         CALL MPI_Abort(comm, 97, ierr_abort)
      END IF
      IF (SIZE(irn) < MAX(1, nz) .OR. SIZE(jcn) < MAX(1, nz) .OR. SIZE(aval) < MAX(1, nz)) THEN
         IF (my_rank == 0) WRITE (*, *) 'MUMPS input check: array size < NZ at IK=', ik, &
            ' size(irn)=', SIZE(irn), ' size(jcn)=', SIZE(jcn), ' size(a)=', SIZE(aval), ' NZ=', nz
         CALL MPI_Abort(comm, 98, ierr_abort)
      END IF

      DO k = 1, nz
         IF (irn(k) < 1 .OR. irn(k) > N) THEN
            WRITE (*, '(A,2(I0,A),2I0)') 'MUMPS input check: bad IRN at IK=', ik, ', k=', k, ': ', irn(k), N
            CALL MPI_Abort(comm, 99, ierr_abort)
         END IF
         IF (jcn(k) < 1 .OR. jcn(k) > N) THEN
            WRITE (*, '(A,2(I0,A),2I0)') 'MUMPS input check: bad JCN at IK=', ik, ', k=', k, ': ', jcn(k), N
            CALL MPI_Abort(comm, 100, ierr_abort)
         END IF
         ar = REAL(aval(k), dp)
         ai = AIMAG(aval(k))
         IF ((.NOT. ieee_is_finite(ar)) .OR. (.NOT. ieee_is_finite(ai))) THEN
            WRITE (*, '(A,2(I0,A),2I0,2(A,ES14.6))') 'MUMPS input check: non-finite A at IK=', ik, ', k=', k, &
               ': row=', irn(k), ' col=', jcn(k), ' re=', ar, ' im=', ai
            CALL MPI_Abort(comm, 101, ierr_abort)
         END IF
      END DO
   END SUBROUTINE validate_mumps_coo_inputs

   SUBROUTINE zmumps_call_safe(id)
      USE, INTRINSIC :: ieee_exceptions, ONLY: ieee_get_halting_mode, ieee_set_halting_mode, &
                                               ieee_invalid, ieee_divide_by_zero, ieee_overflow
      IMPLICIT NONE
      INCLUDE 'zmumps_struc.h'
      TYPE(ZMUMPS_STRUC), INTENT(INOUT) :: id
      LOGICAL :: halt_invalid, halt_divzero, halt_overflow

      CALL ieee_get_halting_mode(ieee_invalid, halt_invalid)
      CALL ieee_get_halting_mode(ieee_divide_by_zero, halt_divzero)
      CALL ieee_get_halting_mode(ieee_overflow, halt_overflow)

      CALL ieee_set_halting_mode(ieee_invalid, .FALSE.)
      CALL ieee_set_halting_mode(ieee_divide_by_zero, .FALSE.)
      CALL ieee_set_halting_mode(ieee_overflow, .FALSE.)

      CALL ZMUMPS(id)

      CALL ieee_set_halting_mode(ieee_invalid, halt_invalid)
      CALL ieee_set_halting_mode(ieee_divide_by_zero, halt_divzero)
      CALL ieee_set_halting_mode(ieee_overflow, halt_overflow)
   END SUBROUTINE zmumps_call_safe

#endif
    !-----------------------------------------------------------------------C
    !                                                                       C
   !  CLUD computes the solution to a COMPLEX system of linear equations   C
   !  A * X = B, where A is a band matrix of order N with KL subdiagonals  C
   !  and KU superdiagonals, and X and B are N-by-NRHS matrices.           C
   !                                                                       C
   !  The LU decomposition with partial pivoting and row interchanges is   C
   !  used to factor A as A = L * U, where L is a product of permutation   C
   !  and unit lower triangular matrices with KL subdiagonals, and U is    C
   !  upper triangular with KL+KU superdiagonals.  The factored form of A  C
   !  is then used to solve the system of equations A * X = B.             C
   !                                                                       C
   !  N       (input) INTEGER                                              C
   !          The number of linear equations, i.e., the order of the       C
   !          matrix A.  N >= 0.                                           C
   !                                                                       C
   !  KL      (input) INTEGER                                              C
   !          The number of subdiagonals within the band of A.  KL >= 0.   C
   !                                                                       C
   !  KU      (input) INTEGER                                              C
   !          The number of superdiagonals within the band of A.  KU >= 0. C
   !                                                                       C
   !                                                                       C
   !  AB      (input/output) COMPLEX array, dimension (LDAB,N)             C
   !          On entry, the matrix A in band storage, in rows KL+1 to      C
   !          2*KL+KU+1; rows 1 to KL of the array need not be set.        C
   !          The j-th column of A is stored in the j-th column of the     C
   !          array AB as follows:                                         C
   !          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)   C
   !          On exit, details of the factorization: U is stored as an     C
   !          upper triangular band matrix with KL+KU superdiagonals in    C
   !          rows 1 to KL+KU+1, and the multipliers used during the       C
   !          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.       C
   !          See below for further details.                               C
   !                                                                       C
   !  LDAB    (input) INTEGER                                              C
   !          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   C
   !                                                                       C
   !  IPIV    (output) INTEGER array, dimension (N)                        C
   !          The pivot indices that define the permutation matrix P;      C
   !          row i of the matrix was interchanged with row IPIV(i).       C
   !                                                                       C
   !                                                                       C
   !  INFO    (output) INTEGER                                             C
   !          = 0:  successful exit                                        C
   !          < 0:  if INFO = -i, the i-th argument had an illegal value   C
   !          > 0:  if INFO = i, U(i,i) is exactly zero.The factorization  C
   !                has been completed, but the factor U is exactly        C
   !                singular, and the solution has not been computed.      C
   !                                                                       C
   !  Further Details                                                      C
   !  ===============                                                      C
   !                                                                       C
   !  The band storage scheme is illustrated by the following example,when C
   !  M = N = 6, KL = 2, KU = 1:                                           C
   !                                                                       C
   !  On entry:                       On exit:                             C
   !                                                                       C
   !      *    *    *    +    +    +       *    *    *   u14  u25  u36     C
   !      *    *    +    +    +    +       *    *   u13  u24  u35  u46     C
   !      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56     C
   !     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66     C
   !     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *      C
   !     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *      C
   !                                                                       C
   !  Array elements marked * are not used by the routine; elements marked C
   !  + need not be set on entry, but are required by the routine to store C
   !  elements of U because of fill-in resulting from the row interchanges.C
   !                                                                       C
   !      so the following scheme may be used to assemble the matrix       C
   !      AB(LDAB,N) from A(i,j):                                          C
   !                                                                       C
   !                           __ AB[KU+KL+1-(j-i),j], (j>=i);             C
   !                          |                                            C
   !                   A[i,j]=|                                            C
   !                          |__ AB[KU+KL+1+(i-j),j], (i>j);              C
   !                                                                       C
   !-----------------------------------------------------------------------C
   SUBROUTINE CLUD(N, KL, KU, LDAB, IPIV, INFO, AB)
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE
      !---- arguments ----
      INTEGER, INTENT(IN)    :: N, KL, KU, LDAB
      INTEGER, INTENT(OUT)   :: INFO
      INTEGER, INTENT(INOUT) :: IPIV(:)          ! length >= N
      COMPLEX(dp), INTENT(INOUT) :: AB(LDAB, N)      ! band matrix (LAPACK gb storage)

      !---- externals ----
      EXTERNAL :: XERBLA, ZGBTRF

      !---- body ----
      INFO = 0

      ! input checks (NRHS not applicable here)
      IF (N < 0) THEN; INFO = -1
      ELSE IF (KL < 0) THEN; INFO = -2
      ELSE IF (KU < 0) THEN; INFO = -3
      ELSE IF (LDAB < 2*KL + KU + 1) THEN; INFO = -4
      END IF

      IF (INFO /= 0) THEN
         CALL XERBLA('CLUD  ', -INFO)
         RETURN
      END IF

      ! LU factorization of complex band matrix (gb storage)
      CALL ZGBTRF(N, N, KL, KU, AB, LDAB, IPIV, INFO)

      IF (INFO /= 0) THEN
         WRITE (*, *) 'ZGBTRF failed with INFO=', INFO
         STOP 'Error in ZGBTRF'
      END IF

      RETURN
   END SUBROUTINE CLUD

   SUBROUTINE SOLVER(TRANS, N, KL, KU, NRHS, LDAB, IPIV, B, LDB, INFO, AB)
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN)    :: TRANS
      INTEGER, INTENT(IN)    :: N, KL, KU, NRHS, LDAB, LDB
      INTEGER, INTENT(OUT)   :: INFO
      INTEGER, INTENT(IN)    :: IPIV(:)
      COMPLEX(dp), INTENT(INOUT) :: AB(LDAB, N)      ! LU factors
      COMPLEX(dp), INTENT(INOUT) :: B(LDB, *)        ! RHS/solution

      EXTERNAL :: ZGBTRS

      CALL ZGBTRS(TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
   END SUBROUTINE SOLVER

   !-----------------------------------------------------------------C
   !  SOLVER solves a system of linear equations                     C
   !     A * X = B,  A**T * X = B,  or  A**H * X = B                 C
   !  with a general band matrix A using the LU factorization        C
   !  computed by CGBTRF.                                            C
   !                                                                 C
   !  Arguments                                                      C
   !  =========                                                      C
   !                                                                 C
   !  TRANS (input) CHARACTER*1                                      C
   !        Specifies the form of the system of equations.           C
   !        = 'N':  A * X = B     (No transpose)                     C
   !        = 'T':  A**T * X = B  (Transpose)                        C
   !        = 'C':  A**H * X = B  (Conjugate transpose)              C
   !                                                                 C
   !  N    (input) INTEGER                                           C
   !       The order of the matrix A.  N >= 0.                       C
   !                                                                 C
   !  KL   (input) INTEGER                                           C
   !       The number of subdiagonals within the band of A. KL >= 0  C
   !                                                                 C
   !  KU   (input) INTEGER                                           C
   !       The number of superdiagonals within the band of A.KU >= 0.C
   !                                                                 C
   !  NRHS (input) INTEGER                                           C
   !       The number of right hand sides,eg, the number of columns  C
   !       of the matrix B.  NRHS >= 0.                              C
   !                                                                 C
   !  AB   (input) COMPLEX(dp)*16 array, dimension (LDAB,N)           C
   !       Details of the LU factorization of the band matrix A, as  C
   !       computed by CGBTRF. U stored as an upper triangular band  C
   !       matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,and C
   !       The multipliers used during the factorization are stored  C
   !       in rows KL+KU+2 to 2*KL+KU+1.                             C
   !                                                                 C
   !  LDAB (input) INTEGER                                           C
   !       The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.C
   !                                                                 C
   !  IPIV (input) INTEGER array, dimension (N)                      C
   !       The pivot indices; for 1 <= i <= N, row i of the matrix   C
   !       was interchanged with row IPIV(i).                        C
   !                                                                 C
   !  B    (input/output) COMPLEX(dp)*16 array, dimension (LDB,NRHS)  C
   !       On entry, the right hand side matrix B.                   C
   !       On exit, the solution matrix X.                           C
   !  LDB  (input) INTEGER                                           C
   !       The leading dimension of the array B.  LDB >= max(1,N).   C
   !                                                                 C
   !  INFO (output) INTEGER                                          C
   !       = 0:  successful exit                                     C
   !       < 0: if INFO = -i,the i-th argument had an illegal value  C
   !-----------------------------------------------------------------C

   SUBROUTINE SOLVER_LEG(TRANS, N, KL, KU, NRHS, LDAB, IPIV, B, LDB, INFO, AB)
  USE iso_fortran_env, ONLY: dp => real64
  IMPLICIT NONE
  !---- arguments ----
  CHARACTER(LEN=1), INTENT(IN)    :: TRANS
  INTEGER,          INTENT(IN)    :: N, KL, KU, NRHS, LDAB, LDB
  INTEGER,          INTENT(OUT)   :: INFO
  INTEGER,          INTENT(IN)    :: IPIV(:)
  COMPLEX(dp),      INTENT(IN)    :: AB(LDAB, N)
  COMPLEX(dp),      INTENT(INOUT) :: B(LDB, *)

  !---- locals ----
  COMPLEX(dp), PARAMETER :: ONE = CMPLX(1.0_dp, 0.0_dp, kind=dp)
  LOGICAL :: LNOTI, NOTRAN
  INTEGER :: I, J, KD, L, LM

  ! external LAPACK/BLAS helpers
  LOGICAL :: LSAME
  EXTERNAL :: LSAME, XERBLA
  EXTERNAL :: ZGEMV, ZSWAP, ZGERU, ZLACGV, ZTBSV

  INFO = 0
  NOTRAN = LSAME(TRANS, 'N')

  IF (.NOT. NOTRAN .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT. LSAME(TRANS, 'C')) THEN
     INFO = -1
  ELSE IF (N    < 0) THEN
     INFO = -2
  ELSE IF (KL   < 0) THEN
     INFO = -3
  ELSE IF (KU   < 0) THEN
     INFO = -4
  ELSE IF (NRHS < 0) THEN
     INFO = -5
  ELSE IF (LDAB < (2*KL + KU + 1)) THEN
     INFO = -7
  ELSE IF (LDB  < MAX(1, N)) THEN
     INFO = -10
  END IF
  IF (INFO /= 0) THEN
     CALL XERBLA('SOLVER', -INFO)
     RETURN
  END IF
  IF (N == 0 .OR. NRHS == 0) RETURN

  KD    = KU + KL + 1
  LNOTI = (KL > 0)

  IF (NOTRAN) THEN
     IF (LNOTI) THEN
        DO J = 1, N-1
           LM = MIN(KL, N-J)
           L  = IPIV(J)
           IF (L /= J) CALL ZSWAP(NRHS, B(L,1), LDB, B(J,1), LDB)
           CALL ZGERU(LM, NRHS, -ONE, AB(KD+1,J), 1, B(J,1), LDB, B(J+1,1), LDB)
        END DO
     END IF

     DO I = 1, NRHS
        CALL ZTBSV('Upper','No transpose','Non-unit', N, KL+KU, AB, LDAB, B(1,I), 1)
     END DO

  ELSE IF (LSAME(TRANS,'T')) THEN

     DO I = 1, NRHS
        CALL ZTBSV('Upper','Transpose','Non-unit', N, KL+KU, AB, LDAB, B(1,I), 1)
     END DO

     IF (LNOTI) THEN
        DO J = N-1, 1, -1
           LM = MIN(KL, N-J)
           CALL ZGEMV('Transpose', LM, NRHS, -ONE, B(J+1,1), LDB, AB(KD+1,J), 1, ONE, B(J,1), LDB)
           L = IPIV(J)
           IF (L /= J) CALL ZSWAP(NRHS, B(L,1), LDB, B(J,1), LDB)
        END DO
     END IF

  ELSE
     ! TRANS = 'C' (conjugate transpose)
     DO I = 1, NRHS
        CALL ZTBSV('Upper','Conjugate transpose','Non-unit', N, KL+KU, AB, LDAB, B(1,I), 1)
     END DO

     IF (LNOTI) THEN
        DO J = N-1, 1, -1
           LM = MIN(KL, N-J)
           CALL ZLACGV(NRHS, B(J,1), LDB)
           CALL ZGEMV('Conjugate transpose', LM, NRHS, -ONE, B(J+1,1), LDB, AB(KD+1,J), 1, ONE, B(J,1), LDB)
           CALL ZLACGV(NRHS, B(J,1), LDB)
           L = IPIV(J)
           IF (L /= J) CALL ZSWAP(NRHS, B(L,1), LDB, B(J,1), LDB)
        END DO
     END IF
  END IF

  RETURN
END SUBROUTINE SOLVER_LEG

! SUBROUTINE run_LU_hybrid(N, KL, KU, LDAB, i1, i2, &
!                                    I25D, IFQ, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, NPT, X, IANISO, CR, CI, NORD, &
!                                    AS, WT, FREQ, NK, FKY, NSR, NSS, ICSR, VSR, XSR, ZSR, MSR, MSR1, FSR,IG, &
!                                    IE0, IS0, DZ0, GFX, GFY, GFZ, my_rank, comm, ITER, NCORE)
!       USE OMP_LIB
!       IMPLICIT DOUBLE PRECISION(A - H, O - Z)
!       INTEGER, INTENT(IN) :: N, KL, KU, LDAB, i1, i2, I25D, IFQ, NTO, NX, NZ, NNX, NNZ, NPT, IANISO, NORD
!       INTEGER, INTENT(IN) :: NK, NSR, NSS, IE0, IS0, my_rank, comm, ITER, NCORE
!       INTEGER, INTENT(IN) :: MSR(:), MSR1(:, :), ICSR(:)
!       DOUBLE PRECISION, INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FREQ, FKY(:)
!       DOUBLE PRECISION, INTENT(IN) :: XSR(:), ZSR(:), CR(:, :), CI(:, :), FSR(:, :), VSR(:, :, :)
!       DOUBLE PRECISION, INTENT(IN) :: DZ0
!       COMPLEX*16, CONTIGUOUS, INTENT(INOUT) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
!    TYPE(InversionGridType), INTENT(IN) :: IG
!       COMPLEX*16, ALLOCATABLE :: A(:, :), S(:, :)
!       DOUBLE PRECISION :: SV(3, 3)
!       INTEGER, ALLOCATABLE :: IPIV(:)

!       DOUBLE PRECISION :: FK, t_elapsed
!       INTEGER :: IS, IC, II, IK, INFO, ierr_alloc, myid, idx, rate, t_start, t_end
!       INTEGER :: J

!       CALL mkl_set_dynamic(0)
!       CALL mkl_set_num_threads(1)         ! MKL single-threaded inside OMP
!       CALL omp_set_num_threads(NCORE)

!       IF (my_rank == 0 .AND. ITER <= 1) WRITE (*, *) 'run_LU_hybrid: NCORE =', NCORE
!       IF (omp_in_parallel()) THEN
!          WRITE (*, *) 'ERROR: run_LU_hybrid entered inside a parallel region'
!          STOP
!       END IF

! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(IK,FK,IPIV,A,SV,S,INFO,ierr_alloc,myid,IS,IC,II,J,&
! !$OMP                                  t_start,t_end,t_elapsed,idx,diff_norm,diff_min,diff_max)
!       myid = OMP_GET_THREAD_NUM()

!       ALLOCATE (IPIV(N), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc IPIV'
!       ALLOCATE (S(N, 3), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc S'
!       ALLOCATE (A(LDAB, N), STAT=ierr_alloc); IF (ierr_alloc /= 0) STOP 'alloc A'

! !$OMP DO SCHEDULE(static)
!       DO IK = i1, i2
!          CALL SYSTEM_CLOCK(t_start, rate)
!          FK = FKY(IK)

!          ! Assemble and factor (per thread, per IK)
!          ! CALL CAGB( FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
!          !         CI, NORD, AS, WT, IE0, IS0, DZ0, N, KL, KU, LDAB, A )
!          CALL CAGB(FREQ, FK, NTO, XTO, ZTO, NX, NZ, NNX, NNZ, X, IANISO, CR, &
!                    CI, NORD, AS, WT, IE0, IS0, DZ0, N, KL, KU, LDAB, IG, A)

!          IF (my_rank == 0) THEN
!          END IF
!          CALL CLUD(N, KL, KU, LDAB, IPIV, INFO, A)

!       SV(1:3, 1:3) = 0D0                
!       SV(1,1)=1.D0
!       SV(2,2)=1.D0
!       SV(3,3)=1.D0
!          DO IS = 1, NSS
!             S(:, :) = (0.0D0, 0.0D0)

!             DO II = 1, MSR(IS)
!                idx = MSR1(IS, II)
!                S(3*idx - 2, 1) = S(3*idx - 2, 1) + DCMPLX(FSR(IS, II)*SV(1, 1), 0.0D0)
!                S(3*idx - 1, 1) = S(3*idx - 1, 1) + DCMPLX(FSR(IS, II)*SV(1, 2), 0.0D0)
!                S(3*idx, 1) = S(3*idx, 1) + DCMPLX(FSR(IS, II)*SV(1, 3), 0.0D0)

!                S(3*idx - 2, 2) = S(3*idx - 2, 2) + DCMPLX(FSR(IS, II)*SV(2, 1), 0.0D0)
!                S(3*idx - 1, 2) = S(3*idx - 1, 2) + DCMPLX(FSR(IS, II)*SV(2, 2), 0.0D0)
!                S(3*idx, 2) = S(3*idx, 2) + DCMPLX(FSR(IS, II)*SV(2, 3), 0.0D0)

!                S(3*idx - 2, 3) = S(3*idx - 2, 3) + DCMPLX(FSR(IS, II)*SV(3, 1), 0.0D0)
!                S(3*idx - 1, 3) = S(3*idx - 1, 3) + DCMPLX(FSR(IS, II)*SV(3, 2), 0.0D0)
!                S(3*idx, 3) = S(3*idx, 3) + DCMPLX(FSR(IS, II)*SV(3, 3), 0.0D0)
!             END DO

!             CALL SOLVER('N', N, KL, KU, 3, LDAB, IPIV, S, N, INFO, A)

!             DO II = 1, NPT
!                GFX(IS, 1, IK, II) = S(3*II - 2, 1)
!                GFY(IS, 1, IK, II) = S(3*II - 1, 1)
!                GFZ(IS, 1, IK, II) = S(3*II, 1)

!                GFX(IS, 2, IK, II) = S(3*II - 2, 2)
!                GFY(IS, 2, IK, II) = S(3*II - 1, 2)
!                GFZ(IS, 2, IK, II) = S(3*II, 2)

!                GFX(IS, 3, IK, II) = S(3*II - 2, 3)
!                GFY(IS, 3, IK, II) = S(3*II - 1, 3)
!                GFZ(IS, 3, IK, II) = S(3*II, 3)
!             END DO
!          END DO


!          CALL SYSTEM_CLOCK(t_end)

!          t_elapsed = (t_end - t_start)/DBLE(rate)
!          IF (ITER <= 1 .AND. my_rank == 0) THEN
!    WRITE (*, '(A,F7.1,A,E11.4,A,I4,A,I4,A,F12.6,A)') 'f=', FREQ, ' Hz | ky=', FK, ' | IK=', IK, '/', NK, ' | t=', t_elapsed, ' s'
!          END IF
!       END DO
! !$OMP END DO

!       DEALLOCATE (A, S, IPIV)
! !$OMP END PARALLEL

!       CALL omp_set_num_threads(1)
!    END SUBROUTINE run_LU_hybrid
end module solvers_mod
