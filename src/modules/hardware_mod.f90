module hardware_mod
    implicit none
! In e.g. hardware_mod
character(len=*), parameter :: BUILD_MODE = &
#ifdef SEIS_BUILD
    SEIS_BUILD
#else
    'UNKNOWN'
#endif

character(len=*), parameter :: COMPILE_DATE = &
#ifdef __DATE__
    __DATE__
#else
    'N/A'
#endif

character(len=*), parameter :: COMPILE_TIME = &
#ifdef __TIME__
    __TIME__
#else
    'N/A'
#endif

    ! ---- compile-time defaults ----
    integer, parameter :: NTHREAD_MAX_DEFAULT = 52
    integer, parameter :: MEM_MAX_DEFAULT     = 312   ! whatever units you’re using

    ! ---- runtime values (overridable via env) ----
    integer :: NTHREAD_MAX = NTHREAD_MAX_DEFAULT
    integer :: MEM_MAX     = MEM_MAX_DEFAULT

contains

   subroutine initialize_hardware(comm)
    use mpi
    use omp_lib
    implicit none

    ! Arguments
    integer, intent(in) :: comm

    ! Local variables
    character(len=32) :: buf
    integer :: istat, iostat_env, tmp, ierr
    integer :: omp_cores, my_rank

    ! Initialize defaults
    NTHREAD_MAX = NTHREAD_MAX_DEFAULT
    MEM_MAX     = MEM_MAX_DEFAULT

    ! Get MPI rank (for optional logging)
    call MPI_Comm_rank(comm, my_rank, ierr)

    ! -----------------------------
    ! Check environment override
    ! -----------------------------
    call get_environment_variable('NTHREAD_MAX', buf, status=istat)
    if (istat == 0 .and. len_trim(buf) > 0) then
        read(buf, *, iostat=iostat_env) tmp
        if (iostat_env == 0 .and. tmp > 0) then
            NTHREAD_MAX = tmp
            if (my_rank == 0) then
                write(*,*) '[hardware_mod] NTHREAD_MAX overridden by environment:', NTHREAD_MAX
            end if
            return
        end if
    end if

    ! -----------------------------
    ! Get available OpenMP cores per rank
    ! -----------------------------
    omp_cores = omp_get_num_procs()
    NTHREAD_MAX = MAX(1, omp_cores)

    ! -----------------------------
    ! Optional logging
    ! -----------------------------
    if (my_rank == 0) then
        write(*,*) '[hardware_mod] Detected OMP cores per rank: ', omp_cores
        write(*,*) '[hardware_mod] NTHREAD_MAX set to: ', NTHREAD_MAX
    end if
end subroutine initialize_hardware



end module hardware_mod
