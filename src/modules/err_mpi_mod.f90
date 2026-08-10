MODULE err_mpi_mod
  USE mpi
  IMPLICIT NONE
  PRIVATE

  INTEGER :: m_comm   = MPI_COMM_WORLD
  INTEGER :: m_rank   = -1
  INTEGER :: m_mpierr = 0

  PUBLIC :: err_init, fail, abort_rank

CONTAINS

  ! Initialize communicator and rank (same signature as before)
  SUBROUTINE err_init(comm, rank)
    INTEGER, INTENT(IN) :: comm, rank
    m_comm = comm
    m_rank = rank
  END SUBROUTINE err_init

  ! Abort whole job with default exit code 999 (same call as before)
  SUBROUTINE fail(msg)
    CHARACTER(*), INTENT(IN) :: msg
    WRITE (*,'(A,I0,A,1X,A)') 'Rank', m_rank, ':', TRIM(msg)
    CALL FLUSH(6)
    CALL MPI_Abort(m_comm, 999, m_mpierr)
  END SUBROUTINE fail

  ! Abort whole job with a specific exit code (same call as before)
  SUBROUTINE abort_rank(msg, code)
    CHARACTER(*), INTENT(IN) :: msg
    INTEGER,      INTENT(IN) :: code
    WRITE (*,'(A,I0,A,1X,A)') 'Rank', m_rank, ':', TRIM(msg)
    CALL FLUSH(6)
    CALL MPI_Abort(m_comm, code, m_mpierr)
  END SUBROUTINE abort_rank

END MODULE err_mpi_mod
