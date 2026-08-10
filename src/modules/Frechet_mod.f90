module Frechet_mod

   USE omp_lib            ! OpenMP runtime library
   USE mkl_service        ! Intel MKL runtime services
   use mpi               ! MPI message-passing
   USE hardware_mod       ! Hardware configuration and memory limits
   USE shared_mod         ! Shared constants, parameters, globals
   USE output_mod         ! Output and logging
   USE partial_derivatives_mod  ! Frechét or gradient computations
   USE gridtype_mod
   USE constant_mod, ONLY: pi
   use iso_fortran_env, only: dp => real64, sp => real32
   IMPLICIT none

   !----------------------------------------------------------------------C
   !this moduld contains all subroutines and functions needed to compute
   !the frechet derivatives:
   ! 1. QFRECHET: compute the frechet derivatives of the displacement vectors
   ! 2. IND_FRECHET: compute individual frechet derivatives
   ! 3. CreateDiagonalMatrix: compute the diagonal of the Hessian matrix
   !--------------------------------------------------------------------C

contains
   !----------------------------------------------------------------------C
!                                                                      C
!     This subroutine calculates the Frechet derivatives of the        C
!     displacement vectors with respect to the independent model       C
!     parameters, for isotropic (3),VTI (6), TTI (7), or general anisotropic   C
!     media (21) parameterizations.                                    C
!                                                                      C
!     Entries:                                                         C
!       IFQ, FREQ.................Frequency index and value            C
!       NX, NZ, X(NX).............Grid dimensions and x-coordinates    C
!       IE0, DZ0..................Absorbing boundary settings          C
!       NPT, NBLOCK...............Number of grid points and blocks     C
!       IANISO, ITHOM, IVISCO.....Anisotropy and viscoelastic flags    C
!       NPAR, INVP(*).............Inversion parameter control          C
!       CR, CI....................Real and imaginary model parameters  C
!       NS, NSV, NR, NRV..........Source and receiver components       C
!       NTO, XTO(*), ZTO(*).......Topography and surface coordinates   C
!       YSR(*), VSR(*,3,3)........Source-receiver geometry             C
!       AS(*), WT(*)..............Gaussian quadrature points/weights   C
!       NK, FKY(*), WTK(*)........Wavenumber sampling                  C
!       GFX, GFY, GFZ.............Green's functions                    C
!       ND, NM....................Number of data and model parameters  C
!       WAVELET(*)................Source wavelet                       C
!
!                                                                      C
!     Return:                                                          C
!       FRECHET(ND,NM)............Jacobian matrix for this frequency   C
!----------------------------------------------------------------------C

   SUBROUTINE QFRECHET(ITER, PARAM, IFQ, FREQ, NSR, ND, NX, NZ, NNX, NNZ, X, IE0, DZ0, NPT, NBLOCK, &
                       IANISO, ITHOM, IVISCO, NPAR, INVP, IS0, CR, CI, &
                       NS, NSV, NR, NRV, NTO, XTO, ZTO, YSR, VSR, &
                       AS, WT, NK, FKY, WTK, GFX, GFY, GFZ, &
                       NM, FRECHET, WAVELET, SourceScaler, NORD, NCOMP, NCOMPS, IG, &
                       my_rank, n_process, comm, DEBUG_OUTPUT)

      IMPLICIT real(dp) (A - H, O - Z)

      CHARACTER(LEN=30) :: FNAME

      CHARACTER(LEN=5), INTENT(IN)  :: PARAM(22)
      INTEGER, INTENT(IN) :: ITER, IFQ, NX, NZ, IANISO, ITHOM, IVISCO, NBLOCK, NPT
      INTEGER, INTENT(IN) :: NK, ND, NM, NTO, NPAR, NORD, NNX, NNZ, NCOMP, NCOMPS, IE0, NSR
      INTEGER, INTENT(IN) :: my_rank, n_process, comm
      INTEGER, CONTIGUOUS, INTENT(IN) :: NS(:), NSV(:), NR(:), NRV(:), INVP(:)
      REAL(dp), INTENT(IN) :: FREQ, DZ0
      REAL(dp), CONTIGUOUS, INTENT(IN) :: XTO(:), ZTO(:), X(:), AS(:), WT(:), FKY(:), WTK(:), YSR(:)
      REAL(dp), CONTIGUOUS, INTENT(IN) :: VSR(:, :, :), CR(:, :), CI(:, :)
      COMPLEX(sp), CONTIGUOUS, INTENT(IN) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      COMPLEX(dp), INTENT(INOUT) :: FRECHET(ND, NM)
      COMPLEX(dp), CONTIGUOUS, INTENT(IN) :: WAVELET(:), SourceScaler(:)
      TYPE(InversionGridType), INTENT(IN) :: IG

      LOGICAL, OPTIONAL, INTENT(IN)  :: DEBUG_OUTPUT
      ! locals
      INTEGER :: NP, J, I, III, remainder, NCORE, ierr, my_tasks, NB, IS0
      INTEGER :: IS, ICS, IR, ICR, NBK, myid, ID
      REAL(dp) :: sum_check
      INTEGER :: frechet_iteration_per_process, start_indexf, end_indexf
      INTEGER :: n_ns, n_nsv, n_nr, n_nrv, n_ysr, n_vsr1, n_vsr2
      INTEGER :: ierr_abort
      REAL(dp) :: Z0, OMIG1, Frec_my, YY, IP, IK
      INTEGER :: clk_start, clk_end, clk_rate, ierr_int
      REAL(dp) :: S(3), G(3)
      REAL(dp)  :: elapsed_time
      LOGICAL :: dbg
      REAL(dp), ALLOCATABLE :: R1(:), R2(:)
      REAL(dp), ALLOCATABLE :: RABS(:), RRMS(:)
      INTEGER :: ID1, ID2
      INTEGER :: npar_active
      REAL(dp) :: WC(NK), WS(NK)
      LOGICAL :: parity
      LOGICAL :: do_rho, do_c11, do_c13, do_c33, do_c44, do_c66, do_zta
      LOGICAL :: do_q11, do_q13, do_q33, do_q44, do_q66
!-------------------------------------------
      FRECHET(:, :) = (0.0_dp, 0.0_dp)
      ierr_abort = 0
      n_ns = SIZE(NS)
      n_nsv = SIZE(NSV)
      n_nr = SIZE(NR)
      n_nrv = SIZE(NRV)
      n_ysr = SIZE(YSR)
      n_vsr1 = SIZE(VSR, 1)
      n_vsr2 = SIZE(VSR, 2)

      IF (ND > n_ns .OR. ND > n_nsv .OR. ND > n_nr .OR. ND > n_nrv) THEN
         WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0)') 'QFRECHET: ND exceeds ID-mapping arrays: ND=', ND, &
            ' SIZE(NS)=', n_ns, ' SIZE(NSV)=', n_nsv, ' SIZE(NR)=', n_nr, ' SIZE(NRV)=', n_nrv
         CALL FLUSH (0)
         CALL MPI_Abort(comm, 811, ierr_abort)
      END IF
      IF (my_rank == 0) THEN
         WRITE (*, *) ' '
         WRITE (*, *) '-------------- calculate Frechet derivatives -----------'
         WRITE (*, *) 'IFREQ =', IFQ, '   ITERATION =', ITER
      END IF
      dbg = PRESENT(DEBUG_OUTPUT) .AND. DEBUG_OUTPUT

      do_rho = (NPAR >= 1) .AND. (INVP(1) .EQ. 1)
      do_c11 = (NPAR >= 2) .AND. (INVP(2) .EQ. 1)
      do_c13 = (NPAR >= 3) .AND. (INVP(3) .EQ. 1)
      do_c33 = (NPAR >= 4) .AND. (INVP(4) .EQ. 1)
      do_c44 = (NPAR >= 5) .AND. (INVP(5) .EQ. 1)
      do_c66 = (NPAR >= 6) .AND. (INVP(6) .EQ. 1)

      ! Keep the active-derivative mapping aligned with the accumulation order below.
      ! For VTI (IANISO=6), slot 7 is Q11, not theta/zeta.
      do_zta = (IANISO .EQ. 7) .AND. (NPAR >= 7) .AND. (INVP(7) .EQ. 1)
      do_q11 = (IVISCO .EQ. 1) .AND. ((IANISO .EQ. 6 .AND. NPAR >= 7 .AND. INVP(7) .EQ. 1) .OR. &
                                      (IANISO .EQ. 7 .AND. NPAR >= 8 .AND. INVP(8) .EQ. 1))
      do_q13 = (IVISCO .EQ. 1) .AND. ((IANISO .EQ. 6 .AND. NPAR >= 8 .AND. INVP(8) .EQ. 1) .OR. &
                                      (IANISO .EQ. 7 .AND. NPAR >= 9 .AND. INVP(9) .EQ. 1))
      do_q33 = (IVISCO .EQ. 1) .AND. ((IANISO .EQ. 6 .AND. NPAR >= 9 .AND. INVP(9) .EQ. 1) .OR. &
                                      (IANISO .EQ. 7 .AND. NPAR >= 10 .AND. INVP(10) .EQ. 1))
      do_q44 = (IVISCO .EQ. 1) .AND. ((IANISO .EQ. 6 .AND. NPAR >= 10 .AND. INVP(10) .EQ. 1) .OR. &
                                      (IANISO .EQ. 7 .AND. NPAR >= 11 .AND. INVP(11) .EQ. 1))
      do_q66 = (IVISCO .EQ. 1) .AND. ((IANISO .EQ. 6 .AND. NPAR >= 11 .AND. INVP(11) .EQ. 1) .OR. &
                                      (IANISO .EQ. 7 .AND. NPAR >= 12 .AND. INVP(12) .EQ. 1))

      npar_active = COUNT(INVP(1:NPAR) .EQ. 1)
      IF (npar_active*NBLOCK /= NM) THEN
         WRITE (0, '(A,I0,A,I0,A,I0,A,I0)') 'QFRECHET: active-parameter packing mismatch: active=', npar_active, &
            ' NBLOCK=', NBLOCK, ' active*NBLOCK=', npar_active*NBLOCK, ' NM=', NM
         CALL FLUSH (0)
         CALL MPI_Abort(comm, 812, ierr_abort)
      END IF

      ! IF (my_rank == 0 .AND. ITER <= 1) THEN
      !    WRITE (*, '(A,I0,A,I0,A)') 'QFRECHET active mapping: ', npar_active, ' active parameters packed into ', NM, ' columns'
      !    J0 = 0
      !    DO I = 1, NPAR
      !       IF (INVP(I) /= 1) CYCLE
      !       J0 = J0 + 1
      !       WRITE (*, '(A,I0,A,A,A,I0,A,I0)') '  slot ', J0, ' <- ', TRIM(PARAM(I)), ' (IA=', I, ', cols ', &
      !          (J0 - 1)*NBLOCK + 1, '-', J0*NBLOCK, ')'
      !    END DO
      ! END IF

      Z0 = IE0*DZ0
      OMIG1 = 2_dp*PI*FREQ
      ! WRITE (*, *) 'NPAR, INVP', NPAR, size(INVP)

      Frec_my = REAL(ND, dp)*REAL(NM, dp)*16.0_dp/1073741824.0_dp
      IF (my_rank == 0 .AND. ITER <= 1) WRITE (*, '(A,F10.3,A)') 'Frechet memory size Frec_my: ', Frec_my, ' (Gb)'

      ! ----------- distribute ID over ranks -----------
      ! CALL split_work(ND, my_rank, n_process, start_indexf, end_indexf, my_tasks)

      IF (n_process > ND) THEN
         IF (my_rank < ND) THEN
            start_indexf = my_rank + 1
            end_indexf = my_rank + 1
         ELSE
            RETURN
         END IF
      ELSE
         frechet_iteration_per_process = ND/n_process
         remainder = MOD(ND, n_process)
         IF (my_rank < remainder) THEN
            start_indexf = my_rank*(frechet_iteration_per_process + 1) + 1
            end_indexf = start_indexf + frechet_iteration_per_process
         ELSE
            start_indexf = my_rank*frechet_iteration_per_process + remainder + 1
            end_indexf = start_indexf + frechet_iteration_per_process - 1
         END IF
      END IF

      IF (start_indexf > ND) RETURN
      end_indexf = MIN(end_indexf, ND)
      my_tasks = end_indexf - start_indexf + 1

if (ITER==1) WRITE(*,'(A,I4,A,I6,A,I6,A,I6, A, I6)') 'Rank ', my_rank, ' -> ID range: ', start_indexf, ' to ', end_indexf, ' | Tasks: ', my_tasks
      call flush (6)

      NCORE = MIN(NTHREAD_MAX, my_tasks)
      CALL omp_set_num_threads(NCORE)

      !     ------------- ID + NBLOCK loop --------------

      DO ID = start_indexf, end_indexf
         CALL SYSTEM_CLOCK(clk_start, clk_rate)

         IS = NS(ID); ICS = NSV(ID)
         IR = NR(ID); ICR = NRV(ID)

         YY = DABS(YSR(IS) - YSR(IR))
         CALL FFT_GL_PRECOMP(ICS, ICR, NK, FKY, WTK, YY, WC, WS, parity)
         S(1:3) = VSR(IS, ICS, 1:3)
         G(1:3) = VSR(IR, ICR, 1:3)

         !   WRITE(*,*) 'S =', S
         !   WRITE(*,*) 'G =', G
         ! write (*, *) 'Rank ', my_rank, ' NNZ, NPT, NBLOCK =', NNZ, NPT, NBLOCK
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(NBK,clk_start,clk_end,elapsed_time,myid)
         myid = OMP_GET_THREAD_NUM()
!$OMP DO SCHEDULE(STATIC)
         DO NBK = 1, NBLOCK
            IF (IG%vPML(NBK)) CYCLE

            CALL IND_FRECHET_BLOCK(ID, NBK, NBLOCK, IS, ICS, NNZ, NPT, IR, ICR, &
                                   X, S, G, ND, NCOMPS, IANISO, ITHOM, CR, CI, NORD, AS, WT, &
                                   NK, FKY, WTK, YY, WC, WS, parity, YSR, VSR, IE0, DZ0, GFX, GFY, GFZ, &
                                   NM, INVP, IVISCO, NPAR, FREQ, ITER, &
                                   do_rho, do_c11, do_c13, do_c33, do_c44, do_c66, do_zta, &
                                   do_q11, do_q13, do_q33, do_q44, do_q66, &
                                   IG, FRECHET)
         END DO
!$OMP END DO
!$OMP END PARALLEL

         CALL SYSTEM_CLOCK(clk_end)
         elapsed_time = DBLE(clk_end - clk_start)/DBLE(clk_rate)

         IF (PRESENT(DEBUG_OUTPUT)) THEN
            ! IF (DEBUG_OUTPUT .AND. ITER == 1 .AND. my_rank == 0) THEN
            IF (ITER == 1 .AND. (MOD(ID, 100) == 0 .OR. ID == ND) .AND. my_rank == 0) THEN

               WRITE (*, '(A,F8.1,A,I8,A,I8,A,I8,A,F10.6,A)') 'FREQ: ', FREQ, ' Hz, ITER=', ITER, &
                  '  Frechet Derivative ID=', ID, '/', ND, ' Finished ', elapsed_time, ' s'

               ! call log_rss("Frechet ID="//trim(adjustl(itoa(id))), 6, my_rank)
               call flush (6)
            END IF
         END IF
      END DO
      ! END IF

      ! collapse threads before MPI barrier
      CALL omp_set_num_threads(1)
      CALL MPI_Barrier(comm, ierr)

      if (dbg) then
         sum_check = 0.0D0
         DO IK = 1, ND
            DO NB = 1, NBLOCK
               sum_check = sum_check + ABS(FRECHET(IK, NB))
            END DO
         END DO
         PRINT *, 'Rank', my_rank, 'local Frechet sum before =', sum_check
      end if

      ! ---- chunked MPI Allreduce on FRECHET ----
      CALL Chunked_AllReduce_Frechet(FRECHET, ND, NM, comm, my_rank, dbg=.false., ierr=ierr)

      if (dbg) then
         sum_check = 0.0D0
         DO IK = 1, ND
            DO NB = 1, NBLOCK
               sum_check = sum_check + ABS(FRECHET(IK, NB))
            END DO
         END DO
         PRINT *, 'Rank', my_rank, 'local Frechet sum after =', sum_check
         WRITE (*, '(A,1ES12.4,A,1ES12.4)') '||Frechet||_F=', SUM(ABS(FRECHET)), 'max|Frechet|=', MAXVAL(ABS(FRECHET))
      end if

      DO i = 1, ND
   FRECHET(i, :) = FRECHET(i, :) * WAVELET(IFQ) * SourceScaler(NS(i))
END DO
     if (my_rank == 0) WRITE (*, '(A,1ES12.4,A,1ES12.4)') '||Frechet||_F=', SUM(ABS(FRECHET)), 'max|Frechet|=', MAXVAL(ABS(FRECHET))
      ! --------- optional export of Frechet derivatives ----------
      ! Optional export (unchanged)
      IF (NSR == 2 .AND. my_rank == 0) THEN
         ALLOCATE (R1(NBLOCK), R2(NBLOCK))
         DO ID = 1, NCOMP
            IP = 0
            DO I = 1, NPAR
               IF (INVP(I) == 1) THEN
                  IP = IP + 1
                  NP = (IP - 1)*NBLOCK
                  DO J = 1, NBLOCK
                     R1(J) = DREAL(FRECHET(ID, NP + J))
                     R2(J) = DIMAG(FRECHET(ID, NP + J))
                  END DO
                  CALL CFNAME_FRECHET('F_RE', ID, PARAM(I), FREQ, '.dat', FNAME)
                  CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, R1, NTO, XTO, ZTO, IE0, IS0)
                  CALL CFNAME_FRECHET('F_IM', ID, PARAM(I), FREQ, '.dat', FNAME)
                  CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, R2, NTO, XTO, ZTO, IE0, IS0)
               END IF
            END DO
         END DO
         DEALLOCATE (R1, R2)
      END IF
!Combined receivers
      IF (NSR == 2 .AND. my_rank == 0) THEN
         ALLOCATE (R1(NBLOCK), R2(NBLOCK), RABS(NBLOCK), RRMS(NBLOCK))

         ID1 = 1
         ID2 = NS(1)*NR(1)*NRV(1)

         IP = 0
         DO I = 1, NPAR
            IF (INVP(I) == 1) THEN
               IP = IP + 1
               NP = (IP - 1)*NBLOCK

               R1(:) = 0.0_dp   ! signed real sum
               R2(:) = 0.0_dp   ! signed imag sum
               RABS(:) = 0.0_dp   ! abs sum
               RRMS(:) = 0.0_dp   ! sum of squares, sqrt at end

               DO ID = ID1, ID2
                  DO J = 1, NBLOCK
                     R1(J) = R1(J) + DREAL(FRECHET(ID, NP + J))
                     R2(J) = R2(J) + DIMAG(FRECHET(ID, NP + J))
                     RABS(J) = RABS(J) + ABS(FRECHET(ID, NP + J))
                     RRMS(J) = RRMS(J) + &
                               DREAL(FRECHET(ID, NP + J))**2 + &
                               DIMAG(FRECHET(ID, NP + J))**2
                  END DO
               END DO

               DO J = 1, NBLOCK
                  RRMS(J) = SQRT(RRMS(J))
               END DO

               ! ---- signed complex sum: real part ----
               CALL CFNAME_FRECHET('F_RE_SRC1_SUM', 1, PARAM(I), FREQ, '.dat', FNAME)
               CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, R1, &
                               NTO, XTO, ZTO, IE0, IS0)

               ! ---- signed complex sum: imaginary part ----
               CALL CFNAME_FRECHET('F_IM_SRC1_SUM', 1, PARAM(I), FREQ, '.dat', FNAME)
               CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, R2, &
                               NTO, XTO, ZTO, IE0, IS0)

               ! ---- absolute-value sum ----
               CALL CFNAME_FRECHET('F_ABS_SRC1_SUM', 1, PARAM(I), FREQ, '.dat', FNAME)
               CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, RABS, &
                               NTO, XTO, ZTO, IE0, IS0)

               ! ---- RMS / energy sum ----
               CALL CFNAME_FRECHET('F_RMS_SRC1_SUM', 1, PARAM(I), FREQ, '.dat', FNAME)
               CALL GRID2D_OUT(FNAME, NX - 1, NZ - 1, IG%XBC, IG%ZBC, RRMS, &
                               NTO, XTO, ZTO, IE0, IS0)

            END IF
         END DO

         DEALLOCATE (R1, R2, RABS, RRMS)
      END IF

      RETURN
   END SUBROUTINE QFRECHET

   SUBROUTINE Chunked_AllReduce_Frechet(FRECHET, ND, NM, comm, my_rank, dbg, ierr)
      IMPLICIT NONE
      COMPLEX(dp), INTENT(INOUT) :: FRECHET(ND, NM)
      INTEGER, INTENT(IN) :: ND, NM, comm, my_rank
      LOGICAL, INTENT(IN) :: dbg
      INTEGER, INTENT(OUT) :: ierr
      INTEGER(KIND=8) :: per_trace, total_size, chunk_64, offset, remaining, chunk_len
      INTEGER :: chunk_start, chunk_end, count32, chunk_num

      per_trace = INT8(ND)
      chunk_64 = (128_8*1024_8*1024_8)/16_8   ! Target ~128 MiB per chunk (COMPLEX(dp)=16 bytes)
      total_size = INT8(ND)*INT8(NM)
      ierr = 0

      IF (ND <= 0 .OR. NM <= 0) THEN
         WRITE (*, *) 'Chunked_AllReduce_Frechet: ND or NM non-positive: ', ND, NM
         STOP
      END IF

      IF (dbg .AND. my_rank == 0) THEN
         WRITE (*, '(A,I0,A,I0,A,I0,A)') 'DEBUG: Starting chunked AllReduce: total_size=', total_size, &
            ' chunk_64=', chunk_64, ' HUGE(count32)=', HUGE(count32)
         WRITE (*, '(A,L1)') 'DEBUG: Using chunked AllReduce? ', &
            total_size > chunk_64 .OR. total_size > INT(HUGE(count32), KIND=8)
      END IF

      IF (total_size > chunk_64 .OR. total_size > INT(HUGE(count32), KIND=8)) THEN
         IF (dbg .AND. my_rank == 0) WRITE (*, '(A)') 'DEBUG: Entering chunked loop'
         remaining = INT8(NM)
         offset = 0_8
         chunk_num = 0

         DO WHILE (remaining > 0_8)
            chunk_num = chunk_num + 1
            IF (per_trace <= 0_8) STOP 'Frechet Chunked_AllReduce: per_trace <= 0'

            chunk_len = chunk_64/per_trace
            IF (chunk_len <= 0_8) chunk_len = 1_8
            IF (chunk_len > remaining) chunk_len = remaining

            IF (per_trace*chunk_len > INT(HUGE(count32), KIND=8)) THEN
               chunk_len = INT(HUGE(count32), KIND=8)/per_trace
               IF (chunk_len <= 0_8) STOP 'Frechet Chunked_AllReduce: cannot fit one chunk'
            END IF

            count32 = INT(per_trace*chunk_len)
            chunk_start = INT(offset) + 1
            chunk_end = chunk_start + INT(chunk_len) - 1

            IF (dbg .AND. my_rank == 0) THEN
               WRITE (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I8)') &
                  'DEBUG Chunk ', chunk_num, ': offset=', offset, ' chunk_len=', chunk_len, &
                  ' cols ', chunk_start, '-', chunk_end, ' count32=', count32
            END IF

            CALL MPI_ALLREDUCE(MPI_IN_PLACE, FRECHET(1, chunk_start), count32, &
                               MPI_COMPLEX16, MPI_SUM, comm, ierr)

            IF (ierr /= 0) THEN
               WRITE (*, '(A,I0,A,I0)') 'ERROR: MPI_ALLREDUCE failed on rank ', my_rank, &
                  ' with ierr=', ierr
               STOP 'Frechet Chunked_AllReduce: MPI_ALLREDUCE failed'
            END IF

            offset = offset + chunk_len
            remaining = INT8(NM) - offset

            IF (dbg .AND. my_rank == 0) THEN
               WRITE (*, '(A,I0,A,I0)') 'DEBUG Chunk ', chunk_num, ' complete. Remaining columns: ', remaining
            END IF
         END DO

         IF (dbg .AND. my_rank == 0) WRITE (*, '(A,I0)') 'DEBUG: Chunked AllReduce completed with ', chunk_num, ' chunks'

      ELSE
         IF (dbg .AND. my_rank == 0) WRITE (*, '(A,I8)') 'DEBUG: Single AllReduce call with count32=', INT(total_size)
         count32 = INT(total_size)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, FRECHET, count32, MPI_COMPLEX16, MPI_SUM, comm, ierr)
         IF (ierr /= 0) THEN
            WRITE (*, '(A,I0,A,I0)') 'ERROR: Single MPI_ALLREDUCE failed on rank ', my_rank, &
               ' with ierr=', ierr
            STOP 'Frechet AllReduce: MPI_ALLREDUCE failed'
         END IF
      END IF
   END SUBROUTINE Chunked_AllReduce_Frechet

   SUBROUTINE IND_FRECHET_BLOCK(ID, NBK, NBLOCK, IS, ICS, NNZ, NPT, IR, ICR, &
                                X, S, G, ND, NCOMPS, IANISO, ITHOM, CR, CI, NORD, AS, WT, &
                                NK, FKY, WTK, YY, WC, WS, parity, YSR, VSR, IE0, DZ0, GFX, GFY, GFZ, &
                                NM, INVP, IVISCO, NPAR, FREQ, ITER, &
                                do_rho, do_c11, do_c13, do_c33, do_c44, do_c66, do_zta, &
                                do_q11, do_q13, do_q33, do_q44, do_q66, &
                                IG, FRECHET)
      IMPLICIT NONE

      ! ------------------ arguments ------------------

      INTEGER, INTENT(IN) :: ID, NBK, NBLOCK, NNZ, NPT, NPAR, IANISO, ITHOM, IVISCO, NORD, ITER
      INTEGER, INTENT(IN) :: IS, ICS, IR, ICR, ND, IE0, NK, NM, NCOMPS
      INTEGER, CONTIGUOUS, INTENT(IN) :: INVP(:)
      REAL(dp), INTENT(IN) :: FREQ, DZ0, YY
      REAL(dp), CONTIGUOUS, INTENT(IN) :: X(:), YSR(:)
      ! REAL(dp), CONTIGUOUS, INTENT(IN) :: AS(:), WT(:),, S(:), G(:)
      REAL(dp), CONTIGUOUS, INTENT(IN) :: AS(:), WT(:), WC(:), WS(:), S(:), G(:), FKY(:), WTK(:)
      REAL(dp), INTENT(IN) :: CR(:, :), CI(:, :), VSR(:, :, :)
      TYPE(InversionGridType), INTENT(IN) :: IG
      COMPLEX(sp), CONTIGUOUS, INTENT(IN) :: GFX(:, :, :, :), GFY(:, :, :, :), GFZ(:, :, :, :)
      COMPLEX(dp), INTENT(INOUT) :: FRECHET(ND, NM)
      LOGICAL, INTENT(IN) :: parity
      LOGICAL, INTENT(IN) :: do_rho, do_c11, do_c13, do_c33, do_c44, do_c66, do_zta
      LOGICAL, INTENT(IN) :: do_q11, do_q13, do_q33, do_q44, do_q66

      ! ------------------ locals  ------------------

      INTEGER :: INDX(2*(NORD - 1) + 1), INDZ(NORD)
      INTEGER :: NOX(NORD), NOZ(NORD)
      REAL(dp) :: DLX(NORD), DLZ(NORD)
      REAL(dp) :: DNX(2*(NORD - 1) + 1), DNZ(NORD), XP(NORD)
      COMPLEX(dp) :: cDNX(2*(NORD - 1) + 1), cDNZ(NORD)
      REAL(dp) :: Z1(NORD), Z2(NORD), T1(NORD), T2(NORD)
      COMPLEX(dp) :: P(IANISO + 1)

      COMPLEX(dp) :: DG_RHO(3, 3), DG_ALP(3, 3), DG_BET(3, 3), DG_EPS(3, 3)
      COMPLEX(dp) :: DG_DEL(3, 3), DG_GAM(3, 3), DG_ZTA(3, 3)
      COMPLEX(dp) :: DG_QVP(3, 3), DG_QVS(3, 3)
      COMPLEX(dp) :: DG_Q11(3, 3), DG_Q13(3, 3), DG_Q33(3, 3), DG_Q44(3, 3), DG_Q66(3, 3)

      COMPLEX(dp) :: D1(NK), D2(NK), D3(NK), D4(NK), D5(NK), D6(NK), D7(NK)
      COMPLEX(dp) :: D8(NK), D9(NK), D10(NK), D11(NK), D12(NK)
      REAL(dp) :: R1(NK), R2(NK)

      REAL(dp) :: C1, C2, C3, DX, DZ, Z0, BL, AK, FK
      REAL(dp) :: OMIG1, ZTA
      REAL(dp) :: F1, F2
      COMPLEX(dp) :: FACT, OMIG2, WIJ
      COMPLEX(dp) :: DGRHO, DGZT, DG11, DG13, DG33, DG44, DG66
      COMPLEX(dp) :: DGQ11, DGQ13, DGQ33, DGQ44, DGQ66
      COMPLEX(dp) :: GS(3), GR(3), DXGS(3), DXGR(3), DZGS(3), DZGR(3)
      COMPLEX(dp) :: tmp_cdnx,tmp_cdnz, tmp_gs_x, tmp_gs_y, tmp_gs_z

      INTEGER :: I, J, K, L, K1, L1, II, JJ, LL, J0, IP, IK, IQ, IZ, IX, I0, NO
      INTEGER :: npt_gfx, nsrc_gfx, nk_gfx, ncmp_gfx, npt_cr, ierr_abort
      INTEGER :: npt_gfy, nsrc_gfy, nk_gfy, ncmp_gfy
      INTEGER :: npt_gfz, nsrc_gfz, nk_gfz, ncmp_gfz
      INTEGER :: npt_gf_min, nsrc_gf_min, nk_gf_min, ncmp_gf_min, npt_all_min
      INTEGER :: min_indx, max_indx, min_indz, max_indz, idx, idz
      INTEGER :: IMP(21), I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
      INTEGER :: active_pos(22)
      INTEGER :: N0
      INTEGER :: nactive_invp
      INTEGER :: izta_idx, iqvp_idx, iqvs_idx
      LOGICAL, PARAMETER :: dbg_bounds = .false.
      REAL(dp) :: t_block_start, t_block_end, t_grad, t_fft, t_accum, t_tmp
      INTEGER :: rank_timer, ierr_rank

      ! Set constants
      FACT = CMPLX(1.0_dp, 0.0_dp, kind=dp)

      Z0 = IE0*DZ0                         ! Depth of bottom of physical model
      OMIG1 = 2.D0*PI*FREQ              ! Angular frequency
      OMIG2 = CMPLX(OMIG1*OMIG1, 0.0_dp, kind=dp)   ! Square of angular frequency (used in density term)
      DG_RHO(:, :) = (0.0_dp, 0.0_dp)
      DG_ALP(:, :) = (0.0_dp, 0.0_dp)
      DG_BET(:, :) = (0.0_dp, 0.0_dp)
      DG_EPS(:, :) = (0.0_dp, 0.0_dp)
      DG_DEL(:, :) = (0.0_dp, 0.0_dp)
      DG_GAM(:, :) = (0.0_dp, 0.0_dp)
      DG_ZTA(:, :) = (0.0_dp, 0.0_dp)
      DG_Q11(:, :) = (0.0_dp, 0.0_dp)
      DG_Q13(:, :) = (0.0_dp, 0.0_dp)
      DG_Q33(:, :) = (0.0_dp, 0.0_dp)
      DG_Q44(:, :) = (0.0_dp, 0.0_dp)
      DG_Q66(:, :) = (0.0_dp, 0.0_dp)
      DGQ11 = (0.0_dp, 0.0_dp)
      DGQ13 = (0.0_dp, 0.0_dp)
      DGQ33 = (0.0_dp, 0.0_dp)
      DGQ44 = (0.0_dp, 0.0_dp)
      DGQ66 = (0.0_dp, 0.0_dp)
      IF (IG%vPML(NBK)) RETURN

      C1 = IG%C1_BLOCK(NBK)
      DX = IG%DX_BLOCK(NBK)
      N0 = IG%N0_BLOCK(NBK)

      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank_timer, ierr_rank)
      t_block_start = MPI_Wtime()
      t_grad = 0.0_dp
      t_fft = 0.0_dp
      t_accum = 0.0_dp

      ierr_abort = 0
      npt_gfx = SIZE(GFX, 4)
      nsrc_gfx = SIZE(GFX, 1)
      ncmp_gfx = SIZE(GFX, 2)
      nk_gfx = SIZE(GFX, 3)
      npt_gfy = SIZE(GFY, 4)
      nsrc_gfy = SIZE(GFY, 1)
      ncmp_gfy = SIZE(GFY, 2)
      nk_gfy = SIZE(GFY, 3)
      npt_gfz = SIZE(GFZ, 4)
      nsrc_gfz = SIZE(GFZ, 1)
      ncmp_gfz = SIZE(GFZ, 2)
      nk_gfz = SIZE(GFZ, 3)

      npt_gf_min = MIN(npt_gfx, MIN(npt_gfy, npt_gfz))
      nsrc_gf_min = MIN(nsrc_gfx, MIN(nsrc_gfy, nsrc_gfz))
      ncmp_gf_min = MIN(ncmp_gfx, MIN(ncmp_gfy, ncmp_gfz))
      nk_gf_min = MIN(nk_gfx, MIN(nk_gfy, nk_gfz))

! CR/CI are indexed as CR(IQ, IP): ensure IP fits their 2nd dim too
      npt_cr = MIN(SIZE(CR, 2), SIZE(CI, 2))
      npt_all_min = MIN(npt_gf_min, npt_cr)

      IF (nsrc_gfx /= nsrc_gfy .OR. nsrc_gfx /= nsrc_gfz .OR. &
          ncmp_gfx /= ncmp_gfy .OR. ncmp_gfx /= ncmp_gfz .OR. &
          nk_gfx /= nk_gfy .OR. nk_gfx /= nk_gfz .OR. &
          npt_gfx /= npt_gfy .OR. npt_gfx /= npt_gfz) THEN
         WRITE (0, '(A)') 'IND_FRECHET_BLOCK: GFX/GFY/GFZ shape mismatch'
         WRITE (0, '(A,4(I0,1X))') '  GFX dims=', nsrc_gfx, ncmp_gfx, nk_gfx, npt_gfx
         WRITE (0, '(A,4(I0,1X))') '  GFY dims=', nsrc_gfy, ncmp_gfy, nk_gfy, npt_gfy
         WRITE (0, '(A,4(I0,1X))') '  GFZ dims=', nsrc_gfz, ncmp_gfz, nk_gfz, npt_gfz
         CALL FLUSH (0)
         CALL MPI_Abort(MPI_COMM_WORLD, 900, ierr_abort)
      END IF

      IF (IS < 1 .OR. IS > nsrc_gf_min .OR. IR < 1 .OR. IR > nsrc_gf_min) THEN
         WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0)') &
            'IND_FRECHET_BLOCK OOB SR: ID=', ID, ' NBK=', NBK, ' IS=', IS, ' IR=', IR, &
            ' valid SR max=', nsrc_gf_min
         WRITE (0, '(A,4(I0,1X))') '  GFX dims=', nsrc_gfx, ncmp_gfx, nk_gfx, npt_gfx
         WRITE (0, '(A,4(I0,1X))') '  GFY dims=', nsrc_gfy, ncmp_gfy, nk_gfy, npt_gfy
         WRITE (0, '(A,4(I0,1X))') '  GFZ dims=', nsrc_gfz, ncmp_gfz, nk_gfz, npt_gfz
         CALL FLUSH (0)
         CALL MPI_Abort(MPI_COMM_WORLD, 901, ierr_abort)
      END IF
      IF (ICS < 1 .OR. ICS > ncmp_gf_min .OR. ICR < 1 .OR. ICR > ncmp_gf_min) THEN
         WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0)') &
            'IND_FRECHET_BLOCK OOB CMP: ID=', ID, ' NBK=', NBK, ' ICS=', ICS, ' ICR=', ICR, &
            ' valid CMP max=', ncmp_gf_min
         WRITE (0, '(A,4(I0,1X))') '  GFX dims=', nsrc_gfx, ncmp_gfx, nk_gfx, npt_gfx
         CALL FLUSH (0)
         CALL MPI_Abort(MPI_COMM_WORLD, 902, ierr_abort)
      END IF
      IF (NK < 1 .OR. NK > nk_gf_min) THEN
         WRITE (0, '(A,I0,A,I0,A,I0,A,I0)') &
            'IND_FRECHET_BLOCK OOB NK: ID=', ID, ' NBK=', NBK, ' NK=', NK, ' valid max=', nk_gf_min
         CALL FLUSH (0)
         CALL MPI_Abort(MPI_COMM_WORLD, 903, ierr_abort)
      END IF

      DO K = 1, NORD
         Z1(K) = IG%Z1_OUT(K, NBK)
         Z2(K) = IG%Z2_OUT(K, NBK)
         T1(K) = IG%T1_OUT(K, NBK)
         T2(K) = IG%T2_OUT(K, NBK)
      END DO

      J0 = 0

      DO I0 = 1, NPAR
      IF (INVP(I0) .EQ. 1) THEN
         J0 = J0 + 1
         IMP(J0) = (J0 - 1)*NBLOCK + NBK
         ! IF (IMP(J0) < 1 .OR. IMP(J0) > NM) THEN
         !    WRITE(*,'(A,I6,A,I6,A,I6)') 'IND_FRECHET_BLOCK: IMP out of range', IMP(J0), ' NM=', NM, ' NBK=', NBK
         !    STOP 2
         ! END IF
         ! FRECHET(ID, IMP(J0)) = (0.D0, 0.D0)  ! Set to zero
      END IF
      END DO

      !---- Loop for GQ abscissa ---------
      ! Main quadrature loop over Gauss points inside the block
      DO 16 K = 1, NORD
         AK = AS(K)
         ! write(*,*)'AS',K,AS(K)
         DZ = Z2(K) - Z1(K)
         C3 = 2.D0/DZ                    ! Derivative scaling in Z

         CALL CDLI(AK, NORD, AS, DLX)         ! Compute Lagrange derivatives in X

         DO 17 L = 1, NORD
            IP = N0 + (K - 1)*NNZ + (L - 1)        ! Global GQG index of the (K,L) point
            IF (IP < 1 .OR. IP > npt_all_min) THEN
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                  'IND_FRECHET_BLOCK OOB IP: ID=', ID, ' NBK=', NBK, ' IP=', IP, ' K=', K, ' L=', L, &
                  ' N0=', N0, ' NNZ=', NNZ
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0)') &
                  '  npt_gf_min=', npt_gf_min, ' npt_cr=', npt_cr, ' NPT=', NPT, ' NBLOCK=', NBLOCK
               CALL FLUSH (0)
               CALL MPI_Abort(MPI_COMM_WORLD, 904, ierr_abort)
            END IF

            BL = AS(L)
            C2 = -((T2(K) - T1(K))*BL + &
                   (T1(K) + T2(K)))/DZ
            CALL CDLI(BL, NORD, AS, DLZ)       ! Lagrange derivative in Z
            WIJ = CMPLX(0.25D0*DX*DZ*WT(K)*WT(L), 0.0_dp, kind=dp)

            ! Load elastic/viscoelastic parameters at this Gauss point
            DO IQ = 1, IANISO
               IF (IVISCO == 0) THEN
                  P(IQ) = CMPLX(CR(IQ, IP), 0.0_dp, kind=dp)
               ELSE
                  P(IQ) = CMPLX(CR(IQ, IP), CI(IQ, IP), kind=dp)
               END IF
            END DO

            ! NOX, NOZ give GQG node indices in X and Z
            DO K1 = 1, NORD
               NOX(K1) = (N0 + (L - 1)) + (K1 - 1)*NNZ
               ! write(*,*)'NOX,NBK',NOX
               NOZ(K1) = (N0 + (K - 1)*NNZ) + (K1 - 1)
            END DO

            ! Fill DNX and INDX for x-gradient
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

            INDX(IX) = NOX(K)
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

            ! Fill DNZ and INDZ for z-gradient
            DO L1 = 1, NORD
               INDZ(L1) = NOZ(L1)
               DNZ(L1) = C3*DLZ(L1)
            END DO
            IZ = NORD

            IF (IX < 1 .OR. IX > SIZE(INDX) .OR. IZ < 1 .OR. IZ > SIZE(INDZ)) THEN
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                  'IND_FRECHET_BLOCK OOB IX/IZ: IX=', IX, ' SIZE(INDX)=', SIZE(INDX), &
                  ' IZ=', IZ, ' SIZE(INDZ)=', SIZE(INDZ), ' ID=', ID, ' NBK=', NBK
               CALL FLUSH (0)
               CALL MPI_Abort(MPI_COMM_WORLD, 910, ierr_abort)
            END IF

            min_indx = MINVAL(INDX(1:IX))
            max_indx = MAXVAL(INDX(1:IX))
            min_indz = MINVAL(INDZ(1:IZ))
            max_indz = MAXVAL(INDZ(1:IZ))

            IF (min_indx < 1 .OR. max_indx > npt_gf_min) THEN
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                  'IND_FRECHET_BLOCK OOB INDX: ID=', ID, ' NBK=', NBK, ' K=', K, ' L=', L, &
                  ' min=', min_indx, ' max=', max_indx, ' npt_gf_min=', npt_gf_min
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0)') &
                  '  IP=', IP, ' N0=', N0, ' NNZ=', NNZ, ' NPT=', NPT
               CALL FLUSH (0)
               CALL MPI_Abort(MPI_COMM_WORLD, 911, ierr_abort)
            END IF
            IF (min_indz < 1 .OR. max_indz > npt_gf_min) THEN
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                  'IND_FRECHET_BLOCK OOB INDZ: ID=', ID, ' NBK=', NBK, ' K=', K, ' L=', L, &
                  ' min=', min_indz, ' max=', max_indz, ' npt_gf_min=', npt_gf_min
               WRITE (0, '(A,I0,A,I0,A,I0,A,I0)') &
                  '  IP=', IP, ' N0=', N0, ' NNZ=', NNZ, ' NPT=', NPT
               CALL FLUSH (0)
               CALL MPI_Abort(MPI_COMM_WORLD, 912, ierr_abort)
            END IF

            ! Precompute complex DNX/DNZ to avoid repeated CMPLX conversions
            DO LL = 1, SIZE(DNX)
               cDNX(LL) = CMPLX(DNX(LL), 0.0_dp, kind=dp)
            END DO
            DO LL = 1, SIZE(DNZ)
               cDNZ(LL) = CMPLX(DNZ(LL), 0.0_dp, kind=dp)
            END DO

            !---- picking up Gs & Gg ----------
            SELECT CASE (IANISO)

            CASE (6)

               ! DO II = 1, 3
               !    IF(NCOMPS==1.AND.II<3) CYCLE  ! Skip horizontal components vertical source
               ! II = ICS
               ! JJ = ICR
               ! DO JJ = 1, 3
               II = ICS
               JJ = ICR
               DO IK = 1, NK
                  FK = FKY(IK)

                  ! IP is already validated once per (K,L), and source/component ranges
                  ! are checked once per block. Keep the hot IK loop free of redundant checks.
                  GS(1) = GFX(IS, II, IK, IP)
                  GS(2) = GFY(IS, II, IK, IP)
                  GS(3) = GFZ(IS, II, IK, IP)

                  GR(1) = GFX(IR, JJ, IK, IP)
                  GR(2) = GFY(IR, JJ, IK, IP)
                  GR(3) = GFZ(IR, JJ, IK, IP)

!  for d_Gs/dx, d_Gg/dx
                  ! t_tmp = MPI_Wtime()
                  DXGS(1:3) = (0.D0, 0.D0)
                  DXGR(1:3) = (0.D0, 0.D0)
                  DZGS(1:3) = (0.D0, 0.D0)
                  DZGR(1:3) = (0.D0, 0.D0)
                  ! Set up the Green's functions matrix

                  DO LL = 1, IX !  computing spatial derivatives using the pre-computed DNX and IND
                     idx = INDX(LL)
                     tmp_cdnx = cDNX(LL)
                     IF (idx < 1 .OR. idx > NPT) THEN
                        WRITE (*, *) 'IND_FRECHET_BLOCK: index range out of bounds:', &
                           ' INDX=', idx, 'NBLOCK=', NBLOCK, 'ID', ID
                        CALL MPI_Abort(MPI_COMM_WORLD, 905, ierr_abort)
                     end if

                     DXGS(1) = DXGS(1) + tmp_cdnx*GFX(IS, II, IK, idx)
                     DXGS(2) = DXGS(2) + tmp_cdnx*GFY(IS, II, IK, idx)
                     DXGS(3) = DXGS(3) + tmp_cdnx*GFZ(IS, II, IK, idx)

                     DXGR(1) = DXGR(1) + tmp_cdnx*GFX(IR, JJ, IK, idx)
                     DXGR(2) = DXGR(2) + tmp_cdnx*GFY(IR, JJ, IK, idx)
                     DXGR(3) = DXGR(3) + tmp_cdnx*GFZ(IR, JJ, IK, idx)
                  END DO

                  ! !for d_Gs/dz, d_Gg/dz

                  DO LL = 1, IZ
                     idz = INDZ(LL)
                     tmp_cdnz = cDNZ(LL)
                     IF (idz < 1 .OR. idz > NPT) THEN
                        WRITE (*, *) 'IND_FRECHET_BLOCK: index range out of bounds:', &
                           ' INDZ=', idz, 'NBLOCK=', NBLOCK,'II', II, 'ID', ID,'JJ',JJ
                        CALL MPI_Abort(MPI_COMM_WORLD, 905, ierr_abort)
                     end if
                     DZGS(1) = DZGS(1) + tmp_cdnz*GFX(IS, II, IK, idz)
                     DZGS(2) = DZGS(2) + tmp_cdnz*GFY(IS, II, IK, idz)
                     DZGS(3) = DZGS(3) + tmp_cdnz*GFZ(IS, II, IK, idz)

                     DZGR(1) = DZGR(1) + tmp_cdnz*GFX(IR, JJ, IK, idz)
                     DZGR(2) = DZGR(2) + tmp_cdnz*GFY(IR, JJ, IK, idz)
                     DZGR(3) = DZGR(3) + tmp_cdnz*GFZ(IR, JJ, IK, idz)
                  END DO

                  IF (do_rho) DGRHO = DRHO0(GS, GR)   !d_Gsg/d_rho

                  ! Compute Parameter Derivatives & Store
                  ! Derivatives w.r.t. physical parameters are computed
                  ! using external functions like DRHO0, DLMDQ, DC11, DQ11, etc.
                  ! Results (D1, D2, ..., D12) stored the kernels for ρ, λ, μ, Q, cij, θ

                  !========================================
                  ! VTI case (IANISO = 6): Moduli or Thomsen
                  !========================================
                  ZTA = 0.D0
                  P(7) = DCMPLX(ZTA, 0.D0)

                  IF (do_c11) DG11 = DC11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c13) DG13 = DC13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c33) DG33 = DC33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c44) DG44 = DC44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c66) DG66 = DC66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_zta) DGZT = DCZTA(P, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO)*Pi/180.D0

                  IF (do_q11) DGQ11 = DQ11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q13) DGQ13 = DQ13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q33) DGQ33 = DQ33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q44) DGQ44 = DQ44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q66) DGQ66 = DQ66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)

                  IF (do_rho) D1(IK) = DGRHO
                  IF (do_c11) D2(IK) = DG11
                  IF (do_c13) D3(IK) = DG13
                  IF (do_c33) D4(IK) = DG33
                  IF (do_c44) D5(IK) = DG44
                  IF (do_c66) D6(IK) = DG66
                  IF (do_zta) D7(IK) = DGZT
                  IF (do_q11) D8(IK) = DGQ11
                  IF (do_q13) D9(IK) = DGQ13
                  IF (do_q33) D10(IK) = DGQ33
                  IF (do_q44) D11(IK) = DGQ44
                  IF (do_q66) D12(IK) = DGQ66

                  !=======================
                  ! Thomsen formulation
                  !=======================
                  IF (ITHOM .EQ. 1) THEN
                     ! CALL  ProcessThomsenDerivatives(P, DG11, DG13, DG33, DG44, DG66, &
                     !   DGRHO, DGALP, DGBET, DGEPS, DGDEL, DGGAM, ITHOM, IK, D1, D2, D3, D4, D5, D6)
                  END IF
               END DO

               ! VTI conversion
               IF (do_rho) THEN
                  R1(1:NK) = REAL(D1(1:NK), dp)
                  R2(1:NK) = AIMAG(D1(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_RHO(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_c11) THEN
                  R1(1:NK) = REAL(D2(1:NK), dp)
                  R2(1:NK) = AIMAG(D2(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_ALP(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C11
               END IF
               IF (do_c13) THEN
                  R1(1:NK) = REAL(D3(1:NK), dp)
                  R2(1:NK) = AIMAG(D3(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_BET(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C13
               END IF
               IF (do_c33) THEN
                  R1(1:NK) = REAL(D4(1:NK), dp)
                  R2(1:NK) = AIMAG(D4(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_EPS(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C33
               END IF
               IF (do_c44) THEN
                  R1(1:NK) = REAL(D5(1:NK), dp)
                  R2(1:NK) = AIMAG(D5(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_DEL(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C44
               END IF
               IF (do_c66) THEN
                  R1(1:NK) = REAL(D6(1:NK), dp)
                  R2(1:NK) = AIMAG(D6(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_GAM(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C66
               END IF
               IF (do_zta) THEN
                  R1(1:NK) = REAL(D7(1:NK), dp)
                  R2(1:NK) = AIMAG(D7(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_ZTA(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q11) THEN
                  R1(1:NK) = REAL(D8(1:NK), dp)
                  R2(1:NK) = AIMAG(D8(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q11(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q13) THEN
                  R1(1:NK) = REAL(D9(1:NK), dp)
                  R2(1:NK) = AIMAG(D9(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q13(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q33) THEN
                  R1(1:NK) = REAL(D10(1:NK), dp)
                  R2(1:NK) = AIMAG(D10(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q33(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q44) THEN
                  R1(1:NK) = REAL(D11(1:NK), dp)
                  R2(1:NK) = AIMAG(D11(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q44(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q66) THEN
                  R1(1:NK) = REAL(D12(1:NK), dp)
                  R2(1:NK) = AIMAG(D12(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q66(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               ! END DO
               ! END DO

            CASE (7)
               ! DO II = 1, 3
               !    IF(NCOMPS==1.AND.II<3) CYCLE `
               II = ICS
               JJ = ICR
               ! DO JJ = 1, 3

               DO IK = 1, NK
                  FK = FKY(IK)
                  ! Green's function at this GQ point for shot and receiver
                  GS(1) = GFX(IS, II, IK, IP)
                  GS(2) = GFY(IS, II, IK, IP)
                  GS(3) = GFZ(IS, II, IK, IP)

                  GR(1) = GFX(IR, JJ, IK, IP)
                  GR(2) = GFY(IR, JJ, IK, IP)
                  GR(3) = GFZ(IR, JJ, IK, IP)

!  for d_Gs/dx, d_Gg/dx
                  ! t_tmp = MPI_Wtime()
                  DXGS(1:3) = (0.D0, 0.D0)
                  DXGR(1:3) = (0.D0, 0.D0)
                  DZGS(1:3) = (0.D0, 0.D0)
                  DZGR(1:3) = (0.D0, 0.D0)
                  ! Set up the Green's functions matrix

                  DO LL = 1, IX !  computing gradients
                     idx = INDX(LL)
                     DXGS(1) = DXGS(1) + cDNX(LL)*GFX(IS, II, IK, idx)
                     ! write(*,*) ' DNX, INDX ,IS, II, IK',LL, DNX(LL), INDX(LL),IS, II, IK
                     DXGS(2) = DXGS(2) + cDNX(LL)*GFY(IS, II, IK, idx)
                     DXGS(3) = DXGS(3) + cDNX(LL)*GFZ(IS, II, IK, idx)

                     DXGR(1) = DXGR(1) + cDNX(LL)*GFX(IR, JJ, IK, idx)
                     DXGR(2) = DXGR(2) + cDNX(LL)*GFY(IR, JJ, IK, idx)
                     DXGR(3) = DXGR(3) + cDNX(LL)*GFZ(IR, JJ, IK, idx)
                  END DO

                  ! !for d_Gs/dz, d_Gg/dz

                  DO LL = 1, IZ
                     idz = INDZ(LL)
                     DZGS(1) = DZGS(1) + cDNZ(LL)*GFX(IS, II, IK, idz)
                     DZGS(2) = DZGS(2) + cDNZ(LL)*GFY(IS, II, IK, idz)
                     DZGS(3) = DZGS(3) + cDNZ(LL)*GFZ(IS, II, IK, idz)

                     DZGR(1) = DZGR(1) + cDNZ(LL)*GFX(IR, JJ, IK, idz)
                     DZGR(2) = DZGR(2) + cDNZ(LL)*GFY(IR, JJ, IK, idz)
                     DZGR(3) = DZGR(3) + cDNZ(LL)*GFZ(IR, JJ, IK, idz)
                  END DO

                  IF (do_rho) DGRHO = DRHO0(GS, GR)   !d_Gsg/d_rho

                  ! Compute Parameter Derivatives & Store
                  ! Derivatives w.r.t. physical parameters are computed
                  ! using external functions like DRHO0, DLMDQ, DC11, DQ11, etc.
                  ! Results (D1, D2, ..., D12) stored the kernels for ρ, λ, μ, Q, cij, θ

                  ! TTI case (IANISO = 7)
                  !--------------------
                  ZTA = REAL(P(7))*PI/180.D0

                  IF (do_c11) DG11 = DC11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c13) DG13 = DC13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c33) DG33 = DC33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c44) DG44 = DC44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_c66) DG66 = DC66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_zta) DGZT = DCZTA(P, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO)*PI/180.D0

                  IF (do_q11) DGQ11 = DQ11(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q13) DGQ13 = DQ13(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q33) DGQ33 = DQ33(ZTA, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q44) DGQ44 = DQ44(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)
                  IF (do_q66) DGQ66 = DQ66(ZTA, FK, GS, GR, DXGS, DXGR, DZGS, DZGR, IVISCO, P)

                  IF (do_rho) D1(IK) = DGRHO
                  IF (do_c11) D2(IK) = DG11
                  IF (do_c13) D3(IK) = DG13
                  IF (do_c33) D4(IK) = DG33
                  IF (do_c44) D5(IK) = DG44
                  IF (do_c66) D6(IK) = DG66
                  IF (do_zta) D7(IK) = DGZT
                  IF (do_q11) D8(IK) = DGQ11
                  IF (do_q13) D9(IK) = DGQ13
                  IF (do_q33) D10(IK) = DGQ33
                  IF (do_q44) D11(IK) = DGQ44
                  IF (do_q66) D12(IK) = DGQ66

               END DO

               IF (do_rho) THEN
                  R1(1:NK) = REAL(D1(1:NK), dp)
                  R2(1:NK) = AIMAG(D1(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_RHO(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_c11) THEN
                  R1(1:NK) = REAL(D2(1:NK), dp)
                  R2(1:NK) = AIMAG(D2(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_ALP(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C11
               END IF
               IF (do_c13) THEN
                  R1(1:NK) = REAL(D3(1:NK), dp)
                  R2(1:NK) = AIMAG(D3(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_BET(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C13
               END IF
               IF (do_c33) THEN
                  R1(1:NK) = REAL(D4(1:NK), dp)
                  R2(1:NK) = AIMAG(D4(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_EPS(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C33
               END IF
               IF (do_c44) THEN
                  R1(1:NK) = REAL(D5(1:NK), dp)
                  R2(1:NK) = AIMAG(D5(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_DEL(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C44
               END IF
               IF (do_c66) THEN
                  R1(1:NK) = REAL(D6(1:NK), dp)
                  R2(1:NK) = AIMAG(D6(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_GAM(II, JJ) = CMPLX(F1, F2, kind=dp)  ! C66
               END IF
               IF (do_zta) THEN
                  R1(1:NK) = REAL(D7(1:NK), dp)
                  R2(1:NK) = AIMAG(D7(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_ZTA(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q11) THEN
                  R1(1:NK) = REAL(D8(1:NK), dp)
                  R2(1:NK) = AIMAG(D8(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q11(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q13) THEN
                  R1(1:NK) = REAL(D9(1:NK), dp)
                  R2(1:NK) = AIMAG(D9(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q13(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q33) THEN
                  R1(1:NK) = REAL(D10(1:NK), dp)
                  R2(1:NK) = AIMAG(D10(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q33(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q44) THEN
                  R1(1:NK) = REAL(D11(1:NK), dp)
                  R2(1:NK) = AIMAG(D11(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q44(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               IF (do_q66) THEN
                  R1(1:NK) = REAL(D12(1:NK), dp)
                  R2(1:NK) = AIMAG(D12(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_Q66(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               ! END DO
            CASE (3)
               ! DO II = 1, 3
               !    IF(NCOMPS==1.AND.II<3) CYCLE
               !    DO JJ = 1, 3
               II = ICS
               JJ = ICR
               DO IK = 1, NK
                  FK = FKY(IK)

                  ! IF (ID == 1 .AND. IK == 1) THEN
                  !    min_indx = MINVAL(INDX(1:IX))
                  !    max_indx = MAXVAL(INDX(1:IX))
                  !    min_indz = MINVAL(INDZ(1:IZ))
                  !    max_indz = MAXVAL(INDZ(1:IZ))
                  !    IF (min_indx < 1 .OR. max_indx > npt_gfx .OR. min_indz < 1 .OR. max_indz > npt_gfx) THEN
                  !       WRITE (*, *) 'IND_FRECHET_BLOCK: index range out of bounds:', &
                  !          ' INDX=', min_indx, ' to ', max_indx, &
                  !          ' INDZ=', min_indz, ' to ', max_indz, &
                  !          ' valid=1..', npt_gfx, ' ID=', ID, ' IK=', IK, ' IX=', IX, ' IZ=', IZ
                  !       CALL MPI_Abort(MPI_COMM_WORLD, 905, ierr_abort)
                  !    END IF
                  ! END IF

                  ! Green's function at this GQ point for shot and receiver
                  GS(1) = GFX(IS, II, IK, IP)
                  GS(2) = GFY(IS, II, IK, IP)
                  GS(3) = GFZ(IS, II, IK, IP)

                  GR(1) = GFX(IR, JJ, IK, IP)
                  GR(2) = GFY(IR, JJ, IK, IP)
                  GR(3) = GFZ(IR, JJ, IK, IP)

!  for d_Gs/dx, d_Gg/dx
                  ! t_tmp = MPI_Wtime()
                  DXGS(1:3) = (0.D0, 0.D0)
                  DXGR(1:3) = (0.D0, 0.D0)
                  DZGS(1:3) = (0.D0, 0.D0)
                  DZGR(1:3) = (0.D0, 0.D0)
                  ! Set up the Green's functions matrix

                  DO LL = 1, IX !  computing gradients
                     idx = INDX(LL)
                     DXGS(1) = DXGS(1) + cDNX(LL)*GFX(IS, II, IK, idx)
                     ! write(*,*) ' DNX, INDX ,IS, II, IK',LL, DNX(LL), INDX(LL),IS, II, IK
                     DXGS(2) = DXGS(2) + cDNX(LL)*GFY(IS, II, IK, idx)
                     DXGS(3) = DXGS(3) + cDNX(LL)*GFZ(IS, II, IK, idx)

                     DXGR(1) = DXGR(1) + cDNX(LL)*GFX(IR, JJ, IK, idx)
                     DXGR(2) = DXGR(2) + cDNX(LL)*GFY(IR, JJ, IK, idx)
                     DXGR(3) = DXGR(3) + cDNX(LL)*GFZ(IR, JJ, IK, idx)
                  END DO

                  ! !for d_Gs/dz, d_Gg/dz

                  DO LL = 1, IZ
                     idz = INDZ(LL)
                     DZGS(1) = DZGS(1) + cDNZ(LL)*GFX(IS, II, IK, idz)
                     DZGS(2) = DZGS(2) + cDNZ(LL)*GFY(IS, II, IK, idz)
                     DZGS(3) = DZGS(3) + cDNZ(LL)*GFZ(IS, II, IK, idz)

                     DZGR(1) = DZGR(1) + cDNZ(LL)*GFX(IR, JJ, IK, idz)
                     DZGR(2) = DZGR(2) + cDNZ(LL)*GFY(IR, JJ, IK, idz)
                     DZGR(3) = DZGR(3) + cDNZ(LL)*GFZ(IR, JJ, IK, idz)
                  END DO

                  DGRHO = DRHO0(GS, GR)   !d_Gsg/d_rho
                  D1(IK) = DGRHO
                  D2(IK) = (0.0_dp, 0.0_dp)
                  D3(IK) = (0.0_dp, 0.0_dp)
                  IF (IVISCO == 1) THEN
                     D8(IK) = (0.0_dp, 0.0_dp)
                     D9(IK) = (0.0_dp, 0.0_dp)
                  END IF
               END DO

               ! isotropic conversion
               R1(1:NK) = REAL(D1(1:NK), dp)
               R2(1:NK) = AIMAG(D1(1:NK))
               CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
               DG_RHO(II, JJ) = CMPLX(F1, F2, kind=dp)

               R1(1:NK) = REAL(D2(1:NK), dp)
               R2(1:NK) = AIMAG(D2(1:NK))
               CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
               DG_ALP(II, JJ) = CMPLX(F1, F2, kind=dp)

               R1(1:NK) = REAL(D3(1:NK), dp)
               R2(1:NK) = AIMAG(D3(1:NK))
               CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
               DG_BET(II, JJ) = CMPLX(F1, F2, kind=dp)

               IF (IVISCO == 1) THEN
                  R1(1:NK) = REAL(D8(1:NK), dp)
                  R2(1:NK) = AIMAG(D8(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_QVP(II, JJ) = CMPLX(F1, F2, kind=dp)

                  R1(1:NK) = REAL(D9(1:NK), dp)
                  R2(1:NK) = AIMAG(D9(1:NK))
                  CALL FFT_GL_W(II, JJ, NK, WC, WS, parity, R1, R2, F1, F2)
                  DG_QVS(II, JJ) = CMPLX(F1, F2, kind=dp)
               END IF
               !    END DO
               ! END DO
            CASE DEFAULT
               t_tmp = MPI_Wtime()
               DO II = 1, 3
                  IF (NCOMPS == 1 .AND. II < 3) CYCLE
                  DO JJ = 1, 3
                     ! II = ICS
                     ! JJ = ICR
                     DO IK = 1, NK
                        FK = FKY(IK)

                        ! IF (ID == 1 .AND. IK == 1) THEN
                        !    min_indx = MINVAL(INDX(1:IX))
                        !    max_indx = MAXVAL(INDX(1:IX))
                        !    min_indz = MINVAL(INDZ(1:IZ))
                        !    max_indz = MAXVAL(INDZ(1:IZ))
                        !    IF (min_indx < 1 .OR. max_indx > npt_gfx .OR. min_indz < 1 .OR. max_indz > npt_gfx) THEN
                        !       WRITE (*, *) 'IND_FRECHET_BLOCK: index range out of bounds:', &
                        !          ' INDX=', min_indx, ' to ', max_indx, &
                        !          ' INDZ=', min_indz, ' to ', max_indz, &
                        !          ' valid=1..', npt_gfx, ' ID=', ID, ' IK=', IK, ' IX=', IX, ' IZ=', IZ
                        !       CALL MPI_Abort(MPI_COMM_WORLD, 905, ierr_abort)
                        !    END IF
                        ! END IF

                        ! Green's function at this GQ point for shot and receiver
                        GS(1) = GFX(IS, II, IK, IP)
                        GS(2) = GFY(IS, II, IK, IP)
                        GS(3) = GFZ(IS, II, IK, IP)

                        GR(1) = GFX(IR, JJ, IK, IP)
                        GR(2) = GFY(IR, JJ, IK, IP)
                        GR(3) = GFZ(IR, JJ, IK, IP)

!  for d_Gs/dx, d_Gg/dx
                        ! t_tmp = MPI_Wtime()
                        DXGS(1:3) = (0.D0, 0.D0)
                        DXGR(1:3) = (0.D0, 0.D0)
                        DZGS(1:3) = (0.D0, 0.D0)
                        DZGR(1:3) = (0.D0, 0.D0)
                        ! Set up the Green's functions matrix
                  DO LL = 1, IX !  computing spatial derivatives using the pre-computed DNX and IND
                     idx = INDX(LL)
                     tmp_cdnx = cDNX(LL)
                     IF (idx < 1 .OR. idx > NPT) THEN
                        WRITE (*, *) 'IND_FRECHET_BLOCK: index range out of bounds:', &
                           ' INDX=', idx, 'NBLOCK=', NBLOCK, 'ID', ID
                        CALL MPI_Abort(MPI_COMM_WORLD, 905, ierr_abort)
                     end if

                     DXGS(1) = DXGS(1) + tmp_cdnx*GFX(IS, II, IK, idx)
                     DXGS(2) = DXGS(2) + tmp_cdnx*GFY(IS, II, IK, idx)
                     DXGS(3) = DXGS(3) + tmp_cdnx*GFZ(IS, II, IK, idx)

                     DXGR(1) = DXGR(1) + tmp_cdnx*GFX(IR, JJ, IK, idx)
                     DXGR(2) = DXGR(2) + tmp_cdnx*GFY(IR, JJ, IK, idx)
                     DXGR(3) = DXGR(3) + tmp_cdnx*GFZ(IR, JJ, IK, idx)
                  END DO

                  ! !for d_Gs/dz, d_Gg/dz

                  DO LL = 1, IZ
                     idz = INDZ(LL)
                     tmp_cdnz = cDNZ(LL)
                     IF (idz < 1 .OR. idz > NPT) THEN
                        WRITE (*, *) 'IND_FRECHET_BLOCK: index range out of bounds:', &
                           ' INDZ=', idz, 'NBLOCK=', NBLOCK,'II', II, 'ID', ID,'JJ',JJ
                        CALL MPI_Abort(MPI_COMM_WORLD, 905, ierr_abort)
                     end if
                     DZGS(1) = DZGS(1) + tmp_cdnz*GFX(IS, II, IK, idz)
                     DZGS(2) = DZGS(2) + tmp_cdnz*GFY(IS, II, IK, idz)
                     DZGS(3) = DZGS(3) + tmp_cdnz*GFZ(IS, II, IK, idz)

                     DZGR(1) = DZGR(1) + tmp_cdnz*GFX(IR, JJ, IK, idz)
                     DZGR(2) = DZGR(2) + tmp_cdnz*GFY(IR, JJ, IK, idz)
                     DZGR(3) = DZGR(3) + tmp_cdnz*GFZ(IR, JJ, IK, idz)
                  END DO

                        DGRHO = DRHO0(GS, GR)   !d_Gsg/d_rho
                     END DO

                     ! General anisotropic media
                     ! t_fft = t_fft + (MPI_Wtime() - t_tmp)
                  END DO
               END DO
            END SELECT
            !---- constant-block's integration -----
            !FACT=DCMPLX(1.D+10,0.D0)  WIJ=DCMPLX(0.25D0*DX*DZ*WT(K)*WT(L),0.D0)
            !apply integration weights and accumulate into global FRECHET matrix
            ! t_tmp = MPI_Wtime()

            J0 = 0
            IF (NPAR >= 1 .AND. INVP(1) .EQ. 1) THEN
               J0 = J0 + 1; I1 = IMP(J0)! IMP is the index of block.
               FRECHET(ID, I1) = FRECHET(ID, I1) + FACT*OMIG2*WIJ*FRG(S, G, DG_RHO)
            END IF

            IF (NPAR >= 2 .AND. INVP(2) .EQ. 1) THEN
               J0 = J0 + 1; I2 = IMP(J0)
               FRECHET(ID, I2) = FRECHET(ID, I2) + FACT*WIJ*FRG(S, G, DG_ALP)!(c11)
            END IF

            IF (NPAR >= 3 .AND. INVP(3) .EQ. 1) THEN
               J0 = J0 + 1; I3 = IMP(J0)
               FRECHET(ID, I3) = FRECHET(ID, I3) + FACT*WIJ*FRG(S, G, DG_BET)!(c13)
            END IF

            IF (IANISO .EQ. 3) THEN
               IF (IVISCO .EQ. 1) THEN
                  IF (NPAR >= 4 .AND. INVP(4) .EQ. 1) THEN
                     J0 = J0 + 1; I4 = IMP(J0)
                     FRECHET(ID, I4) = FRECHET(ID, I4) + FACT*WIJ*FRG(S, G, DG_QVP)
                  END IF
                  IF (NPAR >= 5 .AND. INVP(5) .EQ. 1) THEN
                     J0 = J0 + 1; I5 = IMP(J0)
                     FRECHET(ID, I5) = FRECHET(ID, I5) + FACT*WIJ*FRG(S, G, DG_QVS)
                  END IF
               END IF
               GO TO 16
            END IF

            IF (NPAR >= 4 .AND. INVP(4) .EQ. 1) THEN
               J0 = J0 + 1; I4 = IMP(J0)
               FRECHET(ID, I4) = FRECHET(ID, I4) + FACT*WIJ*FRG(S, G, DG_EPS)!(c33)
            END IF

            IF (NPAR >= 5 .AND. INVP(5) .EQ. 1) THEN
               J0 = J0 + 1; I5 = IMP(J0)
               FRECHET(ID, I5) = FRECHET(ID, I5) + FACT*WIJ*FRG(S, G, DG_DEL)!(c44)
            END IF

            IF (NPAR >= 6 .AND. INVP(6) .EQ. 1) THEN
               J0 = J0 + 1; I6 = IMP(J0)
               FRECHET(ID, I6) = FRECHET(ID, I6) + FACT*WIJ*FRG(S, G, DG_GAM)!(c66)
            END IF

            ! For IANISO=7, slot 7 is theta; for IANISO=6, slot 7 is Q11.
            IF (IANISO .EQ. 7) THEN
               IF (NPAR >= 7 .AND. INVP(7) .EQ. 1) THEN
                  J0 = J0 + 1; I7 = IMP(J0)
                  FRECHET(ID, I7) = FRECHET(ID, I7) + FACT*WIJ*FRG(S, G, DG_ZTA)
               END IF
            END IF

            IF (IVISCO .EQ. 1) THEN
               IF (IANISO .EQ. 6) THEN
                  IF (NPAR >= 7 .AND. INVP(7) .EQ. 1) THEN
                     J0 = J0 + 1; I8 = IMP(J0)
                     FRECHET(ID, I8) = FRECHET(ID, I8) + FACT*WIJ*FRG(S, G, DG_Q11)!(Q11)
                  END IF
                  IF (NPAR >= 8 .AND. INVP(8) .EQ. 1) THEN
                     J0 = J0 + 1; I9 = IMP(J0)
                     FRECHET(ID, I9) = FRECHET(ID, I9) + FACT*WIJ*FRG(S, G, DG_Q13)!(Q13)
                  END IF
                  IF (NPAR >= 9 .AND. INVP(9) .EQ. 1) THEN
                     J0 = J0 + 1; I10 = IMP(J0)
                     FRECHET(ID, I10) = FRECHET(ID, I10) + FACT*WIJ*FRG(S, G, DG_Q33)!(Q33)
                  END IF
                  IF (NPAR >= 10 .AND. INVP(10) .EQ. 1) THEN
                     J0 = J0 + 1; I11 = IMP(J0)
                     FRECHET(ID, I11) = FRECHET(ID, I11) + FACT*WIJ*FRG(S, G, DG_Q44)!(Q44)
                  END IF
                  IF (NPAR >= 11 .AND. INVP(11) .EQ. 1) THEN
                     J0 = J0 + 1; I12 = IMP(J0)
                     FRECHET(ID, I12) = FRECHET(ID, I12) + FACT*WIJ*FRG(S, G, DG_Q66)!(Q66)
                  END IF
               ELSE IF (IANISO .EQ. 7) THEN
                  IF (NPAR >= 8 .AND. INVP(8) .EQ. 1) THEN
                     J0 = J0 + 1; I8 = IMP(J0)
                     FRECHET(ID, I8) = FRECHET(ID, I8) + FACT*WIJ*FRG(S, G, DG_Q11)!(Q11)
                  END IF
                  IF (NPAR >= 9 .AND. INVP(9) .EQ. 1) THEN
                     J0 = J0 + 1; I9 = IMP(J0)
                     FRECHET(ID, I9) = FRECHET(ID, I9) + FACT*WIJ*FRG(S, G, DG_Q13)!(Q13)
                  END IF
                  IF (NPAR >= 10 .AND. INVP(10) .EQ. 1) THEN
                     J0 = J0 + 1; I10 = IMP(J0)
                     FRECHET(ID, I10) = FRECHET(ID, I10) + FACT*WIJ*FRG(S, G, DG_Q33)!(Q33)
                  END IF
                  IF (NPAR >= 11 .AND. INVP(11) .EQ. 1) THEN
                     J0 = J0 + 1; I11 = IMP(J0)
                     FRECHET(ID, I11) = FRECHET(ID, I11) + FACT*WIJ*FRG(S, G, DG_Q44)!(Q44)
                  END IF
                  IF (NPAR >= 12 .AND. INVP(12) .EQ. 1) THEN
                     J0 = J0 + 1; I12 = IMP(J0)
                     FRECHET(ID, I12) = FRECHET(ID, I12) + FACT*WIJ*FRG(S, G, DG_Q66)!(Q66)
                  END IF
               END IF
            END IF
            ! t_accum = t_accum + (MPI_Wtime() - t_tmp)

17          CONTINUE !end of GQG-points
16          CONTINUE !end of GQG-points

            RETURN
            END SUBROUTINE IND_FRECHET_BLOCK
            ! !                     ! !                     ! !

            ! !----------------------------------------------------------------------C
            ! !                                                                      C
            ! !     Equations (41) for the Frechet derivative: dGij/dm               C
            ! !                                                  FRG=ST⋅GF⋅G              C
            ! !----------------------------------------------------------------------C
            FUNCTION FRG(S, G, GF)
               IMPLICIT NONE
               REAL(dp), INTENT(IN) :: S(3), G(3)
               COMPLEX(dp), INTENT(IN) :: GF(3, 3)
               COMPLEX(dp) :: FRG
               FRG = CMPLX(S(1), 0.0_dp, dp)*(GF(1, 1)*CMPLX(G(1), 0.0_dp, dp) &
                                              + GF(1, 2)*CMPLX(G(2), 0.0_dp, dp) &
                                              + GF(1, 3)*CMPLX(G(3), 0.0_dp, dp)) &
                     + CMPLX(S(2), 0.0_dp, dp)*(GF(2, 1)*CMPLX(G(1), 0.0_dp, dp) &
                                                + GF(2, 2)*CMPLX(G(2), 0.0_dp, dp) &
                                                + GF(2, 3)*CMPLX(G(3), 0.0_dp, dp)) &
                     + CMPLX(S(3), 0.0_dp, dp)*(GF(3, 1)*CMPLX(G(1), 0.0_dp, dp) &
                                                + GF(3, 2)*CMPLX(G(2), 0.0_dp, dp) &
                                                + GF(3, 3)*CMPLX(G(3), 0.0_dp, dp))
               RETURN
            END

            end module Frechet_mod!  ! Calculate expected memory for the full Hessian (each element is 16 bytes)
