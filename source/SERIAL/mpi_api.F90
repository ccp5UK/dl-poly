Module mpi_api

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for faking MPI calls needed for serial compilation
  !
  ! NOTE: The MPI-I/O interfaces are only valid for a single
  !       datarep='native' ASCII stream and thus can only deal with
  !       one type per one file at a time (a serious restriction) !!!
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov june 2010
  ! contrib   - a.m.elena march 2016
  ! contrib   - m.a.seaton june 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use asserts,  Only: assert
  Use kinds,    Only: li,&
                      si,&
                      wp
  Use particle, Only: corePart

  Implicit None

  Include 'mpif.h' !  Include 'mpiof.h' ! Needed for ScaliMPI

  ! MPI thread level support
  Integer, Parameter   :: MPI_THREAD_FUNNELED = 0

  ! MPI address kind

  Integer, Parameter   :: MPI_ADDRESS_KIND = li

  ! MPI-I/O SPECIFICS
  ! MPI-I/O kind = high-precision integer 'kinds, only : wp' module

  Integer, Parameter   :: MPI_OFFSET_KIND = li

  ! MPI-I/O info

  Integer, Save        :: MPI_INFO_NULL = 0

  ! MPI-I/O bookkeeping

  Integer, Parameter                                               :: mpi_io_max = 64
  Integer, Save :: mpi_io_cnt = 1
  Integer, Dimension(1:mpi_io_max), Save :: mpi_io_rec_len = 0
  Integer, Dimension(1:mpi_io_max), Save :: mpi_io_etype = 0
  Integer, Dimension(1:mpi_io_max), Save :: mpi_io_ftype = 0
  Integer, Dimension(1:mpi_io_max), Save :: mpi_io_comm = 0
  Character(Len=40), Dimension(1:mpi_io_max), Save :: mpi_io_fname = ' '
  Integer, Dimension(1:mpi_io_max), Save :: mpi_io_mode = 0
  Integer, Dimension(1:mpi_io_max), Save :: mpi_io_fh = 0
  Integer(Kind=MPI_OFFSET_KIND), Dimension(1:mpi_io_max), Save :: mpi_io_disp = 0_MPI_OFFSET_KIND
  Character(Len=6), Dimension(1:mpi_io_max), Save :: mpi_io_datarep = ' '
  Public :: MPI_INIT, MPI_INIT_THREAD, MPI_FINALIZE, &
            MPI_ABORT, MPI_COMM_RANK, &
            MPI_COMM_SIZE, MPI_COMM_DUP, MPI_COMM_SPLIT, &
            MPI_COMM_FREE, MPI_BARRIER, MPI_WAIT, MPI_WAITALL, &
            MPI_WTIME, MPI_BCAST, MPI_ALLREDUCE, MPI_GATHERV, &
            MPI_ALLGATHER, MPI_ALLTOALL, MPI_ALLTOALLV, &
            MPI_SEND, MPI_ISEND, MPI_ISSEND, MPI_RECV, MPI_IRECV, &
            MPI_SCATTER, MPI_SCATTERV, &
            MPI_OFFSET_KIND, MPI_INFO_NULL, &
            MPI_TYPE_GET_EXTENT, MPI_TYPE_EXTENT, &
            MPI_TYPE_CREATE_HVECTOR, MPI_TYPE_HVECTOR, &
            MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_FREE, &
            MPI_FILE_DELETE, MPI_FILE_OPEN, MPI_FILE_CLOSE, &
            MPI_FILE_SET_VIEW, MPI_FILE_GET_VIEW, MPI_GET_COUNT, &
            MPI_FILE_WRITE_AT, MPI_FILE_READ_AT, &
            MPI_GET_PROCESSOR_NAME, MPI_GET_LIBRARY_VERSION, &
            MPI_GET_VERSION, MPI_TYPE_CREATE_STRUCT, &
            MPI_GET_ADDRESS, MPI_TYPE_CREATE_RESIZED

  Private :: mpi_io_max, mpi_io_cnt, mpi_io_rec_len, &
             mpi_io_etype, mpi_io_ftype, mpi_io_comm, &
             mpi_io_fname, mpi_io_mode, mpi_io_fh, &
             mpi_io_disp, mpi_io_datarep

  Interface MPI_ABORT
    Module Procedure MPI_ABORT_s
    Module Procedure MPI_ABORT_l
  End Interface !MPI_ABORT

  Interface MPI_BCAST
    Module Procedure MPI_BCAST_log_s
    Module Procedure MPI_BCAST_log_v
    Module Procedure MPI_BCAST_chr_s
    Module Procedure MPI_BCAST_chr_v
    Module Procedure MPI_BCAST_int_s
    Module Procedure MPI_BCAST_int_v
    Module Procedure MPI_BCAST_int_s_16
    Module Procedure MPI_BCAST_rwp_s
    Module Procedure MPI_BCAST_rwp_v
    Module Procedure MPI_BCAST_rwp_m
  End Interface !MPI_BCAST

  Interface MPI_ALLREDUCE
    Module Procedure MPI_ALLREDUCE_log_s
    Module Procedure MPI_ALLREDUCE_log_v
    Module Procedure MPI_ALLREDUCE_int_s
    Module Procedure MPI_ALLREDUCE_int_v
    Module Procedure MPI_ALLREDUCE_rwp_s
    Module Procedure MPI_ALLREDUCE_rwp_v
    Module Procedure MPI_ALLREDUCE_rwp_m
    Module Procedure MPI_ALLREDUCE_cwp_v
  End Interface !MPI_ALLREDUCE

  Interface MPI_GATHERV
    Module Procedure MPI_GATHERV_log11
    Module Procedure MPI_GATHERV_log22
    Module Procedure MPI_GATHERV_chr11
    Module Procedure MPI_GATHERV_chr22
    Module Procedure MPI_GATHERV_int11
    Module Procedure MPI_GATHERV_int22
    Module Procedure MPI_GATHERV_rwp11
    Module Procedure MPI_GATHERV_rwp22
    Module Procedure MPI_GATHERV_cwp11
    Module Procedure MPI_GATHERV_cwp22
  End Interface !MPI_GATHERV

  Interface MPI_ALLGATHER
    Module Procedure MPI_ALLGATHER_log_s
    Module Procedure MPI_ALLGATHER_log_v
    Module Procedure MPI_ALLGATHER_chr_s
    Module Procedure MPI_ALLGATHER_chr_v
    Module Procedure MPI_ALLGATHER_int_s
    Module Procedure MPI_ALLGATHER_int_v
    Module Procedure MPI_ALLGATHER_rwp_s
    Module Procedure MPI_ALLGATHER_rwp_v
    Module Procedure MPI_ALLGATHER_cwp_s
    Module Procedure MPI_ALLGATHER_cwp_v
  End Interface !MPI_ALLGATHER

  Interface MPI_ALLGATHERV
    Module Procedure MPI_ALLGATHER_chr_vv
    Module Procedure MPI_ALLGATHER_rwp_vv
    Module Procedure MPI_ALLGATHER_rwp_mm
  End Interface MPI_ALLGATHERV !MPI_ALLGATHERV

  Interface MPI_ALLTOALL
    Module Procedure MPI_ALLTOALL_log11
    Module Procedure MPI_ALLTOALL_log22
    Module Procedure MPI_ALLTOALL_chr11
    Module Procedure MPI_ALLTOALL_chr22
    Module Procedure MPI_ALLTOALL_int11
    Module Procedure MPI_ALLTOALL_int22
    Module Procedure MPI_ALLTOALL_rwp11
    Module Procedure MPI_ALLTOALL_rwp22
    Module Procedure MPI_ALLTOALL_cwp11
    Module Procedure MPI_ALLTOALL_cwp22
  End Interface !MPI_ALLTOALL

  Interface MPI_ALLTOALLV
    Module Procedure MPI_ALLTOALLV_log11
    Module Procedure MPI_ALLTOALLV_log22
    Module Procedure MPI_ALLTOALLV_chr11
    Module Procedure MPI_ALLTOALLV_chr22
    Module Procedure MPI_ALLTOALLV_int11
    Module Procedure MPI_ALLTOALLV_int22
    Module Procedure MPI_ALLTOALLV_rwp11
    Module Procedure MPI_ALLTOALLV_rwp22
    Module Procedure MPI_ALLTOALLV_cwp11
    Module Procedure MPI_ALLTOALLV_cwp22
  End Interface !MPI_ALLTOALLV

  Interface MPI_SEND
    Module Procedure MPI_SEND_log_s
    Module Procedure MPI_SEND_log_v
    Module Procedure MPI_SEND_chr_s
    Module Procedure MPI_SEND_chr_v
    Module Procedure MPI_SEND_int_s
    Module Procedure MPI_SEND_int_v
    Module Procedure MPI_SEND_rwp_s
    Module Procedure MPI_SEND_rwp_v
    Module Procedure MPI_SEND_rwp_m
    Module Procedure MPI_SEND_rwp_c
    Module Procedure MPI_SEND_cwp_s
    Module Procedure MPI_SEND_cwp_v
    Module Procedure MPI_SEND_cwp_m
    Module Procedure MPI_SEND_cwp_c
    Module Procedure MPI_SEND_part_s
    Module Procedure MPI_SEND_part_v
  End Interface !MPI_SEND

  Interface MPI_ISEND
    Module Procedure MPI_ISEND_log_s
    Module Procedure MPI_ISEND_log_v
    Module Procedure MPI_ISEND_chr_s
    Module Procedure MPI_ISEND_chr_v
    Module Procedure MPI_ISEND_int_s
    Module Procedure MPI_ISEND_int_v
    Module Procedure MPI_ISEND_rwp_s
    Module Procedure MPI_ISEND_rwp_v
    Module Procedure MPI_ISEND_rwp_m
    Module Procedure MPI_ISEND_rwp_c
    Module Procedure MPI_ISEND_rwp_f
    Module Procedure MPI_ISEND_cwp_s
    Module Procedure MPI_ISEND_cwp_v
    Module Procedure MPI_ISEND_cwp_m
    Module Procedure MPI_ISEND_cwp_c
  End Interface !MPI_ISEND

  Interface MPI_ISSEND
    Module Procedure MPI_ISSEND_log_s
    Module Procedure MPI_ISSEND_log_v
    Module Procedure MPI_ISSEND_chr_s
    Module Procedure MPI_ISSEND_chr_v
    Module Procedure MPI_ISSEND_int_s
    Module Procedure MPI_ISSEND_int_v
    Module Procedure MPI_ISSEND_rwp_s
    Module Procedure MPI_ISSEND_rwp_v
    Module Procedure MPI_ISSEND_rwp_m
    Module Procedure MPI_ISSEND_rwp_c
    Module Procedure MPI_ISSEND_cwp_s
    Module Procedure MPI_ISSEND_cwp_v
    Module Procedure MPI_ISSEND_cwp_m
    Module Procedure MPI_ISSEND_cwp_c
  End Interface !MPI_ISSEND

  Interface MPI_RECV
    Module Procedure MPI_RECV_log_s
    Module Procedure MPI_RECV_log_v
    Module Procedure MPI_RECV_chr_s
    Module Procedure MPI_RECV_chr_v
    Module Procedure MPI_RECV_int_s
    Module Procedure MPI_RECV_int_v
    Module Procedure MPI_RECV_rwp_s
    Module Procedure MPI_RECV_rwp_v
    Module Procedure MPI_RECV_rwp_m
    Module Procedure MPI_RECV_rwp_c
    Module Procedure MPI_RECV_cwp_s
    Module Procedure MPI_RECV_cwp_v
    Module Procedure MPI_RECV_cwp_m
    Module Procedure MPI_RECV_cwp_c
    Module Procedure MPI_RECV_part_s
    Module Procedure MPI_RECV_part_v
  End Interface !MPI_RECV

  Interface MPI_IRECV
    Module Procedure MPI_IRECV_log_s
    Module Procedure MPI_IRECV_log_v
    Module Procedure MPI_IRECV_chr_s
    Module Procedure MPI_IRECV_chr_v
    Module Procedure MPI_IRECV_int_s
    Module Procedure MPI_IRECV_int_v
    Module Procedure MPI_IRECV_rwp_s
    Module Procedure MPI_IRECV_rwp_v
    Module Procedure MPI_IRECV_rwp_m
    Module Procedure MPI_IRECV_rwp_c
    Module Procedure MPI_IRECV_rwp_f
    Module Procedure MPI_IRECV_cwp_s
    Module Procedure MPI_IRECV_cwp_v
    Module Procedure MPI_IRECV_cwp_m
    Module Procedure MPI_IRECV_cwp_c
  End Interface !MPI_IRECV

  Interface MPI_SCATTER
    Module Procedure MPI_SCATTER_log_ss
    Module Procedure MPI_SCATTER_log_sv
    Module Procedure MPI_SCATTER_log_vs
    Module Procedure MPI_SCATTER_log_vv
    Module Procedure MPI_SCATTER_chr_ss
    Module Procedure MPI_SCATTER_chr_sv
    Module Procedure MPI_SCATTER_chr_vs
    Module Procedure MPI_SCATTER_chr_vv
    Module Procedure MPI_SCATTER_int_ss
    Module Procedure MPI_SCATTER_int_sv
    Module Procedure MPI_SCATTER_int_vs
    Module Procedure MPI_SCATTER_int_vv
    Module Procedure MPI_SCATTER_rwp_ss
    Module Procedure MPI_SCATTER_rwp_sv
    Module Procedure MPI_SCATTER_rwp_vs
    Module Procedure MPI_SCATTER_rwp_vv
    Module Procedure MPI_SCATTER_cwp_ss
    Module Procedure MPI_SCATTER_cwp_sv
    Module Procedure MPI_SCATTER_cwp_vs
    Module Procedure MPI_SCATTER_cwp_vv
  End Interface !MPI_SCATTER

  Interface MPI_SCATTERV
    Module Procedure MPI_SCATTERV_log_vv
    Module Procedure MPI_SCATTERV_chr_vv
    Module Procedure MPI_SCATTERV_int_vv
    Module Procedure MPI_SCATTERV_rwp_vv
    Module Procedure MPI_SCATTERV_rwp_mm
    Module Procedure MPI_SCATTERV_cwp_vv
    Module Procedure MPI_SCATTERV_cwp_mm
  End Interface !MPI_SCATTERV

  Interface MPI_FILE_WRITE_AT
    Module Procedure MPI_FILE_WRITE_AT_chr_s
    Module Procedure MPI_FILE_WRITE_AT_chr_v_1
    Module Procedure MPI_FILE_WRITE_AT_chr_m_1
  End Interface !MPI_FILE_WRITE_AT

  Interface MPI_FILE_READ_AT
    Module Procedure MPI_FILE_READ_AT_chr_s
    Module Procedure MPI_FILE_READ_AT_chr_v_1
    Module Procedure MPI_FILE_READ_AT_chr_m_1
  End Interface !MPI_FILE_READ_AT

  Interface MPI_GET_ADDRESS
    Module Procedure MPI_GET_ADDRESS_class
  End Interface MPI_GET_ADDRESS

Contains

  Subroutine MPI_INIT(ierr)

    Integer, Intent(  Out) :: ierr

    ierr = 0

  End Subroutine MPI_INIT

  Subroutine MPI_INIT_THREAD(required_threading, provided_threading, ierr)

    Integer, Intent(In   ) :: required_threading
    Integer, Intent(  Out) :: provided_threading, ierr

    provided_threading = 0
    ierr = 0

  End Subroutine MPI_INIT_THREAD

  Subroutine MPI_FINALIZE(ierr)

    Integer, Intent(  Out) :: ierr

    ierr = 0

  End Subroutine MPI_FINALIZE

  Subroutine MPI_ABORT_s(ierr)

    Integer, Intent(  Out) :: ierr

    ierr = 99
    Stop

  End Subroutine MPI_ABORT_s
  Subroutine MPI_ABORT_l(MPI_COMM_WORLD, sh_error, ierr)

    Integer, Intent(In   ) :: MPI_COMM_WORLD, sh_error
    Integer, Intent(  Out) :: ierr

    ierr = 99
    Stop

  End Subroutine MPI_ABORT_l

  Subroutine MPI_COMM_RANK(MPI_COMM_WORLD, idnode, ierr)

    Integer, Intent(In   ) :: MPI_COMM_WORLD
    Integer, Intent(  Out) :: idnode, ierr

    idnode = 0
    ierr = 0

  End Subroutine MPI_COMM_RANK

  Subroutine MPI_COMM_SIZE(MPI_COMM_WORLD, mxnode, ierr)

    Integer, Intent(In   ) :: MPI_COMM_WORLD
    Integer, Intent(  Out) :: mxnode, ierr

    mxnode = 1
    ierr = 0

  End Subroutine MPI_COMM_SIZE

  Subroutine MPI_COMM_DUP(COMM, NEW_COMM, ierr)

    Integer, Intent(In   ) :: COMM
    Integer, Intent(  Out) :: NEW_COMM, ierr

    NEW_COMM = COMM + 1
    ierr = 0

  End Subroutine MPI_COMM_DUP

  Subroutine MPI_COMM_SPLIT(COMM, COLOR, KEY, NEW_COMM, ierr)

    Integer, Intent(In   ) :: COMM, COLOR, KEY
    Integer, Intent(  Out) :: NEW_COMM, ierr

    NEW_COMM = COMM + 1
    ierr = 0

  End Subroutine MPI_COMM_SPLIT

  Subroutine MPI_COMM_FREE(COMM, ierr)

    Integer, Intent(InOut) :: COMM
    Integer, Intent(  Out) :: ierr

    ierr = 0

  End Subroutine MPI_COMM_FREE

  Subroutine MPI_BARRIER(MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0

  End Subroutine MPI_BARRIER

  Subroutine MPI_WAIT(request, status, ierr)

    Integer, Intent(In   ) :: request
    Integer, Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0

  End Subroutine MPI_WAIT

  Subroutine MPI_WAITALL(n, request, status, ierr)

    Integer, Intent(In   ) :: n, request(:)
    Integer, Intent(  Out) :: status(MPI_STATUS_SIZE, *), ierr

    ierr = 0
    status(MPI_STATUS_SIZE, 1:n) = 0

  End Subroutine MPI_WAITALL

  Function MPI_WTIME()

    Real(Kind=wp) :: MPI_WTIME

    Character(Len=10)      :: time
    Character(Len=5)       :: zone
    Character(Len=8)       :: date1
    Character(Len=8), Save :: date = ' '
    Integer                :: value(1:8)
    Integer, Save          :: days = 0
    Logical, Save          :: newjob = .true.

    If (newjob) Then
      newjob = .false.

      Call Date_and_time(date, time, zone, value)

      MPI_WTIME = Real(value(5), wp) * 3600.0_wp + Real(value(6), wp) * 60.0_wp + &
                  Real(value(7), wp) + Real(value(8), wp) / 1000.0_wp
    Else
      Call Date_and_time(date1, time, zone, value)

      ! time-per-timestep & start-up and close-down times
      ! are assumed to be shorter than 24h

      If (date /= date1) Then
        date = date1
        days = days + 1
      End If

      MPI_WTIME = Real(value(5), wp) * 3600.0_wp + Real(value(6), wp) * 60.0_wp + &
                  Real(value(7), wp) + Real(value(8), wp) / 1000.0_wp + &
                  Real(days, wp) * 86400.0_wp

    End If

  End Function MPI_WTIME

  Subroutine MPI_BCAST_log_s(aaa, n, MPI_LOGICAL, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_log_s
  Subroutine MPI_BCAST_log_v(aaa, n, MPI_LOGICAL, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n > Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_log_v
  Subroutine MPI_BCAST_chr_s(aaa, n, MPI_CHARACTER, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n > Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_chr_s
  Subroutine MPI_BCAST_chr_v(aaa, n, MPI_CHARACTER, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n > Size(aaa) * Len(aaa) .or. n < 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_chr_v
  Subroutine MPI_BCAST_int_s(aaa, n, MPI_INTEGER, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa, n, MPI_INTEGER, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_int_s
  Subroutine MPI_BCAST_int_s_16(aaa, n, MPI_INTEGER, idnode, MPI_COMM_WORLD, ierr)

    Integer(Kind=si), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_INTEGER, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_int_s_16
  Subroutine MPI_BCAST_int_v(aaa, n, MPI_INTEGER, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), n, MPI_INTEGER, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n > Size(aaa) .or. n < 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_int_v
  Subroutine MPI_BCAST_rwp_s(aaa, n, MPI_WP, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_rwp_s

  Subroutine MPI_BCAST_rwp_v(aaa, n, MPI_WP, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n > Size(aaa) .or. n < 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_rwp_v

  Subroutine MPI_BCAST_rwp_m(aaa, n, MPI_WP, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n > Size(aaa) .or. n < 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_BCAST_rwp_m

  Subroutine MPI_ALLREDUCE_log_s(aaa, bbb, n, MPI_LOGICAL, MPI, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa
    Logical, Intent(InOut) :: bbb
    Integer, Intent(In   ) :: n, MPI_LOGICAL, MPI, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (n /= 1) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLREDUCE_log_s
  Subroutine MPI_ALLREDUCE_log_v(aaa, bbb, n, MPI_LOGICAL, MPI, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Logical, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, MPI, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (n > Size(aaa) .or. n > Size(bbb) .or. n < 1) Then
      ierr = 1
      Stop
    End If
    bbb(1:n) = aaa(1:n)

  End Subroutine MPI_ALLREDUCE_log_v
  Subroutine MPI_ALLREDUCE_int_s(aaa, bbb, n, MPI_INTEGER, MPI, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa
    Integer, Intent(InOut) :: bbb
    Integer, Intent(In   ) :: n, MPI_INTEGER, MPI, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (n /= 1) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLREDUCE_int_s
  Subroutine MPI_ALLREDUCE_int_v(aaa, bbb, n, MPI_INTEGER, MPI, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:)
    Integer, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: n, MPI_INTEGER, MPI, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (n > Size(aaa) .or. n > Size(bbb) .or. n < 1) Then
      ierr = 1
      Stop
    End If
    bbb(1:n) = aaa(1:n)

  End Subroutine MPI_ALLREDUCE_int_v
  Subroutine MPI_ALLREDUCE_rwp_s(aaa, bbb, n, MPI_WP, MPI, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Real(Kind=wp), Intent(InOut) :: bbb
    Integer,       Intent(In   ) :: n, MPI_WP, MPI, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (n /= 1) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLREDUCE_rwp_s

  Subroutine MPI_ALLREDUCE_rwp_v(aaa, bbb, n, MPI_WP, MPI, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Real(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,       Intent(In   ) :: n, MPI_WP, MPI, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (n > Size(aaa) .or. n > Size(bbb) .or. n < 1) Then
      ierr = 1
      Stop
    End If
    bbb(1:n) = aaa(1:n)

  End Subroutine MPI_ALLREDUCE_rwp_v

  Subroutine MPI_ALLREDUCE_cwp_v(aaa, bbb, n, MPI_WP, MPI, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Complex(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: n, MPI_WP, MPI, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (n > Size(aaa) .or. n > Size(bbb) .or. n < 1) Then
      ierr = 1
      Stop
    End If
    bbb(1:n) = aaa(1:n)

  End Subroutine MPI_ALLREDUCE_cwp_v

  Subroutine MPI_ALLREDUCE_rwp_m(aaa, bbb, n, MPI_WP, MPI, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Real(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, MPI, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (n /= Size(aaa) .or. n /= Size(bbb) .or. n < 1) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLREDUCE_rwp_m

  Subroutine MPI_GATHERV_log11(aaa, s_a, MPI_LOGICALa, bbb, r_b, disp, MPI_LOGICALb, root, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: s_a, MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: r_b(:), disp(:), MPI_LOGICALb, root, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: i, j, k

    i = Ubound(r_b, Dim=1)
    j = Ubound(disp, Dim=1)

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_a /= r_b(i) .or. i < 1 .or. j < 1 .or. root < 0) Then
      ierr = 1
      Stop
    End If

    k = disp(j)
    bbb(k + 1:k + s_a) = aaa(1:s_a)

  End Subroutine MPI_GATHERV_log11
  Subroutine MPI_GATHERV_log22(aaa, s_a, MPI_LOGICALa, bbb, r_b, disp, MPI_LOGICALb, root, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:, :)
    Integer, Intent(In   ) :: s_a, MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:, :)
    Integer, Intent(In   ) :: r_b(:), disp(:), MPI_LOGICALb, root, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer              :: i, j, k, l, m, n
    Logical, Allocatable :: aa1(:), bb1(:)

    ierr = 0

    i = Size(aaa, Dim=1)
    j = Size(aaa, Dim=2)
    k = i * j

    If (k < s_a) Then
      ierr = 1
      Stop
    End If

    Allocate (aa1(1:k), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 2
      Stop
    End If

    aa1(1:k) = Reshape(aaa, (/k/))

    l = Size(bbb, Dim=1)
    m = Size(bbb, Dim=2)
    n = l * m

    Allocate (bb1(1:n), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 3
      Stop
    End If

    bb1(1:n) = Reshape(bbb, (/n/))

    Call MPI_GATHERV_log11(aa1, s_a, MPI_LOGICALa, bb1, r_b, disp, MPI_LOGICALb, root, MPI_COMM_WORLD, ierr)

    bbb = Reshape(bb1(1:n), (/l, m/))

    Deallocate (aa1, bb1, Stat=ierr)
    If (ierr /= 0) Then
      ierr = 4
      Stop
    End If

  End Subroutine MPI_GATHERV_log22
  Subroutine MPI_GATHERV_chr11(aaa, s_a, MPI_CHARACTERa, bbb, r_b, disp, MPI_CHARACTERb, root, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: s_a, MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: r_b(:), disp(:), MPI_CHARACTERb, root, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: i, j, k

    i = Ubound(r_b, Dim=1)
    j = Ubound(disp, Dim=1)

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_a /= r_b(i) .or. Len(aaa) /= Len(bbb) .or. i < 1 .or. j < 1 .or. root < 0) Then
      ierr = 1
      Stop
    End If

    k = disp(j)
    bbb(k + 1:k + s_a) = aaa(1:s_a)

  End Subroutine MPI_GATHERV_chr11
  Subroutine MPI_GATHERV_chr22(aaa, s_a, MPI_CHARACTERa, bbb, r_b, disp, MPI_CHARACTERb, root, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: s_a, MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:, :)
    Integer,          Intent(In   ) :: r_b(:), disp(:), MPI_CHARACTERb, root, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Character(Len=Len(aaa)), Allocatable :: aa1(:), bb1(:)
    Integer                              :: i, j, k, l, m, n

    ierr = 0

    i = Size(aaa, Dim=1)
    j = Size(aaa, Dim=2)
    k = i * j

    If (k < s_a) Then
      ierr = 1
      Stop
    End If

    Allocate (aa1(1:k), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 2
      Stop
    End If

    aa1(1:k) = Reshape(aaa, (/k/))

    l = Size(bbb, Dim=1)
    m = Size(bbb, Dim=2)
    n = l * m

    Allocate (bb1(1:n), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 3
      Stop
    End If

    bb1(1:n) = Reshape(bbb, (/n/))

    Call MPI_GATHERV_chr11(aa1, s_a, MPI_CHARACTERa, bb1, r_b, disp, MPI_CHARACTERb, root, MPI_COMM_WORLD, ierr)

    bbb = Reshape(bb1(1:n), (/l, m/))

    Deallocate (aa1, bb1, Stat=ierr)
    If (ierr /= 0) Then
      ierr = 4
      Stop
    End If

  End Subroutine MPI_GATHERV_chr22
  Subroutine MPI_GATHERV_int11(aaa, s_a, MPI_INTEGERa, bbb, r_b, disp, MPI_INTEGERb, root, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), s_a, MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: r_b(:), disp(:), MPI_INTEGERb, root, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: i, j, k

    i = Ubound(r_b, Dim=1)
    j = Ubound(disp, Dim=1)

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_a /= r_b(i) .or. i < 1 .or. j < 1 .or. root < 0) Then
      ierr = 1
      Stop
    End If

    k = disp(j)
    bbb(k + 1:k + s_a) = aaa(1:s_a)

  End Subroutine MPI_GATHERV_int11
  Subroutine MPI_GATHERV_int22(aaa, s_a, MPI_INTEGERa, bbb, r_b, disp, MPI_INTEGERb, root, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:, :), s_a, MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:, :)
    Integer, Intent(In   ) :: r_b(:), disp(:), MPI_INTEGERb, root, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer              :: i, j, k, l, m, n
    Integer, Allocatable :: aa1(:), bb1(:)

    ierr = 0

    i = Size(aaa, Dim=1)
    j = Size(aaa, Dim=2)
    k = i * j

    If (k < s_a) Then
      ierr = 1
      Stop
    End If

    Allocate (aa1(1:k), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 2
      Stop
    End If

    aa1(1:k) = Reshape(aaa, (/k/))

    l = Size(bbb, Dim=1)
    m = Size(bbb, Dim=2)
    n = l * m

    Allocate (bb1(1:n), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 3
      Stop
    End If

    bb1(1:n) = Reshape(bbb, (/n/))

    Call MPI_GATHERV_int11(aa1, s_a, MPI_INTEGERa, bb1, r_b, disp, MPI_INTEGERb, root, MPI_COMM_WORLD, ierr)

    bbb = Reshape(bb1(1:n), (/l, m/))

    Deallocate (aa1, bb1, Stat=ierr)
    If (ierr /= 0) Then
      ierr = 4
      Stop
    End If

  End Subroutine MPI_GATHERV_int22
  Subroutine MPI_GATHERV_rwp11(aaa, s_a, MPI_WPa, bbb, r_b, disp, MPI_WPb, root, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: s_a, MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,       Intent(In   ) :: r_b(:), disp(:), MPI_WPb, root, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: i, j, k

    i = Ubound(r_b, Dim=1)
    j = Ubound(disp, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_a /= r_b(i) .or. i < 1 .or. j < 1 .or. root < 0) Then
      ierr = 1
      Stop
    End If

    k = disp(j)
    bbb(k + 1:k + s_a) = aaa(1:s_a)

  End Subroutine MPI_GATHERV_rwp11
  Subroutine MPI_GATHERV_rwp22(aaa, s_a, MPI_WPa, bbb, r_b, disp, MPI_WPb, root, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: s_a, MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,       Intent(In   ) :: r_b(:), disp(:), MPI_WPb, root, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer                    :: i, j, k, l, m, n
    Real(Kind=wp), Allocatable :: aa1(:), bb1(:)

    ierr = 0

    i = Size(aaa, Dim=1)
    j = Size(aaa, Dim=2)
    k = i * j

    If (k < s_a) Then
      ierr = 1
      Stop
    End If

    Allocate (aa1(1:k), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 2
      Stop
    End If

    aa1(1:k) = Reshape(aaa, (/k/))

    l = Size(bbb, Dim=1)
    m = Size(bbb, Dim=2)
    n = l * m

    Allocate (bb1(1:n), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 3
      Stop
    End If

    bb1(1:n) = Reshape(bbb, (/n/))

    Call MPI_GATHERV_rwp11(aa1, s_a, MPI_WPa, bb1, r_b, disp, MPI_WPb, root, MPI_COMM_WORLD, ierr)

    bbb = Reshape(bb1(1:n), (/l, m/))

    Deallocate (aa1, bb1, Stat=ierr)
    If (ierr /= 0) Then
      ierr = 4
      Stop
    End If

  End Subroutine MPI_GATHERV_rwp22
  Subroutine MPI_GATHERV_cwp11(aaa, s_a, MPI_WPa, bbb, r_b, disp, MPI_WPb, root, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: s_a, MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: r_b(:), disp(:), MPI_WPb, root, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: i, j, k

    i = Ubound(r_b, Dim=1)
    j = Ubound(disp, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_a /= r_b(i) .or. i < 1 .or. j < 1 .or. root < 0) Then
      ierr = 1
      Stop
    End If

    k = disp(j)
    bbb(k + 1:k + s_a) = aaa(1:s_a)

  End Subroutine MPI_GATHERV_cwp11
  Subroutine MPI_GATHERV_cwp22(aaa, s_a, MPI_WPa, bbb, r_b, disp, MPI_WPb, root, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: s_a, MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,          Intent(In   ) :: r_b(:), disp(:), MPI_WPb, root, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Complex(Kind=wp), Allocatable :: aa1(:), bb1(:)
    Integer                       :: i, j, k, l, m, n

    ierr = 0

    i = Size(aaa, Dim=1)
    j = Size(aaa, Dim=2)
    k = i * j

    If (k < s_a) Then
      ierr = 1
      Stop
    End If

    Allocate (aa1(1:k), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 2
      Stop
    End If

    aa1(1:k) = Reshape(aaa, (/k/))

    l = Size(bbb, Dim=1)
    m = Size(bbb, Dim=2)
    n = l * m

    Allocate (bb1(1:n), Stat=ierr)
    If (ierr /= 0) Then
      ierr = 3
      Stop
    End If

    bb1(1:n) = Reshape(bbb, (/n/))

    Call MPI_GATHERV_cwp11(aa1, s_a, MPI_WPa, bb1, r_b, disp, MPI_WPb, root, MPI_COMM_WORLD, ierr)

    bbb = Reshape(bb1(1:n), (/l, m/))

    Deallocate (aa1, bb1, Stat=ierr)
    If (ierr /= 0) Then
      ierr = 4
      Stop
    End If

  End Subroutine MPI_GATHERV_cwp22

  Subroutine MPI_ALLGATHER_log_s(aaa, s_a, MPI_LOGICALa, bbb, s_b, MPI_LOGICALb, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: s_a, MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: s_b, MPI_LOGICALb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. s_a /= 1 .or. s_b < 1 .or. s_b > Size(bbb)) Then
      ierr = 1
      Stop
    End If
    bbb(1) = aaa

  End Subroutine MPI_ALLGATHER_log_s
  Subroutine MPI_ALLGATHER_log_v(aaa, s_a, MPI_LOGICALa, bbb, s_b, MPI_LOGICALb, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: s_a, MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: s_b, MPI_LOGICALb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a < s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_b) = aaa(1:s_b)

  End Subroutine MPI_ALLGATHER_log_v
  Subroutine MPI_ALLGATHER_chr_s(aaa, s_a, MPI_CHARACTERa, bbb, s_b, MPI_CHARACTERb, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: s_a, MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: s_b, MPI_CHARACTERb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: i

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. s_a > Len(aaa) .or. s_a < Len(bbb) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a < s_b) Then
      ierr = 1
      Stop
    End If

    Do i = 1, s_a
      bbb(1) (i:i) = aaa(i:i)
    End Do

  End Subroutine MPI_ALLGATHER_chr_s
  Subroutine MPI_ALLGATHER_chr_v(aaa, s_a, MPI_CHARACTERa, bbb, s_b, MPI_CHARACTERb, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: s_a, MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: s_b, MPI_CHARACTERb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. Len(aaa) /= Len(bbb) .or. &
        s_a < 1 .or. s_a > Size(aaa) .or. s_b < 1 .or. s_b > Size(bbb) .or. s_a < s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_b) = aaa(1:s_b)

  End Subroutine MPI_ALLGATHER_chr_v
  Subroutine MPI_ALLGATHER_int_s(aaa, s_a, MPI_INTEGERa, bbb, s_b, MPI_INTEGERb, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa, s_a, MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: s_b, MPI_INTEGERb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. s_a /= 1 .or. s_b < 1 .or. s_b > Size(bbb)) Then
      ierr = 1
      Stop
    End If
    bbb(1) = aaa

  End Subroutine MPI_ALLGATHER_int_s
  Subroutine MPI_ALLGATHER_int_v(aaa, s_a, MPI_INTEGERa, bbb, s_b, MPI_INTEGERb, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), s_a, MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: s_b, MPI_INTEGERb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a < s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_b) = aaa(1:s_b)

  End Subroutine MPI_ALLGATHER_int_v
  Subroutine MPI_ALLGATHER_rwp_s(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: s_a, MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,       Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a /= 1 .or. s_b < 1 .or. s_b > Size(bbb)) Then
      ierr = 1
      Stop
    End If
    bbb(1) = aaa

  End Subroutine MPI_ALLGATHER_rwp_s
  Subroutine MPI_ALLGATHER_rwp_v(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: s_a, MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,       Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a < s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_b) = aaa(1:s_b)

  End Subroutine MPI_ALLGATHER_rwp_v
  Subroutine MPI_ALLGATHER_cwp_s(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: s_a, MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a /= 1 .or. s_b < 1 .or. s_b > Size(bbb)) Then
      ierr = 1
      Stop
    End If
    bbb(1) = aaa

  End Subroutine MPI_ALLGATHER_cwp_s
  Subroutine MPI_ALLGATHER_cwp_v(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: s_a, MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a < s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_b) = aaa(1:s_b)

  End Subroutine MPI_ALLGATHER_cwp_v

  Subroutine MPI_ALLGATHER_chr_vv(send_buf, send_size, MPI_CHARACTER_send, &
                                  recv_buf, recv_counts, disps, MPI_CHARACTER_recv, &
                                  MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: send_buf(:)
    Integer,          Intent(In   ) :: send_size, MPI_CHARACTER_send
    Character(Len=*), Intent(InOut) :: recv_buf(:)
    Integer,          Intent(In   ) :: recv_counts(:), disps(:), MPI_CHARACTER_recv, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Call assert(Sum(recv_counts) == Size(recv_buf), &
                message="Sum of send buffer sizes in receive_counts does not equal total size of receive buffer")
    Call assert(MPI_CHARACTER_send == MPI_CHARACTER_recv, &
                message="MPI_CHARACTER_send /= MPI_CHARACTER_recv")
    Call assert(Len(send_buf(1)) == Len(recv_buf(1)), &
                message='Character length of send and receive buffers differ')
    Call assert(Size(send_buf) == Size(recv_buf), message="Send and receive sizes differ")

    ierr = 0
    recv_buf = send_buf

  End Subroutine MPI_ALLGATHER_chr_vv

  Subroutine MPI_ALLGATHER_rwp_vv(send_buf, send_size, MPI_REAL_WP_send, &
                                  recv_buf, recv_counts, disps, MPI_REAL_WP_recv, &
                                  MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: send_buf(:)
    Integer,       Intent(In   ) :: send_size, MPI_REAL_WP_send
    Real(Kind=wp), Intent(InOut) :: recv_buf(:)
    Integer,       Intent(In   ) :: recv_counts(:), disps(:), MPI_REAL_WP_recv, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Call assert(Sum(recv_counts) == Size(recv_buf), &
                message="Sum of send buffer sizes in receiv _counts does not equal total size of receive buffer")
    Call assert(MPI_REAL_WP_send == MPI_REAL_WP_recv, message="MPI_REAL_WP_send /= MPI_REAL_WP_recv")
    Call assert(Size(send_buf) == Size(recv_buf), message="Send and receive sizes differ")

    ierr = 0
    recv_buf = send_buf

  End Subroutine MPI_ALLGATHER_rwp_vv

  Subroutine mpi_allgather_rwp_mm(send_buf, send_size, MPI_REAL_WP_send, &
                                  recv_buf, recv_counts, disps, MPI_REAL_WP_recv, &
                                  MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: send_buf(:, :)
    Integer,       Intent(In   ) :: send_size, MPI_REAL_WP_send
    Real(Kind=wp), Intent(InOut) :: recv_buf(:, :)
    Integer,       Intent(In   ) :: recv_counts(:), disps(:), MPI_REAL_WP_recv, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Call assert(Sum(recv_counts) == Size(recv_buf), &
                message="Sum of send buffer sizes in receive counts does not equal total size of receive buffer")
    Call assert(MPI_REAL_WP_send == MPI_REAL_WP_recv, &
                message="MPI_REAL_WP_send /= MPI_REAL_WP_recv")
    Call assert(Size(send_buf, Dim=1) == Size(recv_buf, Dim=1), &
                message="Dimension 1 of send and receive arrays differ")
    Call assert(Size(send_buf, Dim=2) == Size(recv_buf, Dim=2), &
                message="Dimension 2 of send and receive arrays differ")

    ierr = 0
    recv_buf = send_buf

  End Subroutine mpi_allgather_rwp_mm

  Subroutine MPI_ALLTOALL_log11(aaa, s_a, MPI_LOGICALa, bbb, s_b, MPI_LOGICALb, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: s_a, MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: s_b, MPI_LOGICALb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_a) = aaa(1:s_a)

  End Subroutine MPI_ALLTOALL_log11
  Subroutine MPI_ALLTOALL_log22(aaa, s_a, MPI_LOGICALa, bbb, s_b, MPI_LOGICALb, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:, :)
    Integer, Intent(In   ) :: s_a, MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:, :)
    Integer, Intent(In   ) :: s_b, MPI_LOGICALb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: i, j, k

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If

    k = 0
    X: Do i = 1, Ubound(aaa, Dim=2)
      Do j = 1, Ubound(bbb, Dim=2)
        bbb(j, i) = aaa(j, i)
        k = k + 1
        If (k == s_a) Then
          Exit X
        End If
      End Do
    End Do X

  End Subroutine MPI_ALLTOALL_log22
  Subroutine MPI_ALLTOALL_chr11(aaa, s_a, MPI_CHARACTERa, bbb, s_b, MPI_CHARACTERb, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: s_a, MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: s_b, MPI_CHARACTERb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b .or. Len(aaa) /= Len(bbb)) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_a) = aaa(1:s_a)

  End Subroutine MPI_ALLTOALL_chr11
  Subroutine MPI_ALLTOALL_chr22(aaa, s_a, MPI_CHARACTERa, bbb, s_b, MPI_CHARACTERb, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: s_a, MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:, :)
    Integer,          Intent(In   ) :: s_b, MPI_CHARACTERb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: i, j, k

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b .or. Len(aaa) /= Len(bbb)) Then
      ierr = 1
      Stop
    End If

    k = 0
    X: Do i = 1, Ubound(aaa, Dim=2)
      Do j = 1, Ubound(bbb, Dim=2)
        bbb(j, i) = aaa(j, i)
        k = k + 1
        If (k == s_a) Then
          Exit X
        End If
      End Do
    End Do X

  End Subroutine MPI_ALLTOALL_chr22
  Subroutine MPI_ALLTOALL_int11(aaa, s_a, MPI_INTEGERa, bbb, s_b, MPI_INTEGERb, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), s_a, MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: s_b, MPI_INTEGERb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_a) = aaa(1:s_a)

  End Subroutine MPI_ALLTOALL_int11
  Subroutine MPI_ALLTOALL_int22(aaa, s_a, MPI_INTEGERa, bbb, s_b, MPI_INTEGERb, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:, :), s_a, MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:, :)
    Integer, Intent(In   ) :: s_b, MPI_INTEGERb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: i, j, k

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If

    k = 0
    X: Do i = 1, Ubound(aaa, Dim=2)
      Do j = 1, Ubound(bbb, Dim=2)
        bbb(j, i) = aaa(j, i)
        k = k + 1
        If (k == s_a) Then
          Exit X
        End If
      End Do
    End Do X

  End Subroutine MPI_ALLTOALL_int22
  Subroutine MPI_ALLTOALL_rwp11(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: s_a, MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,       Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_a) = aaa(1:s_a)

  End Subroutine MPI_ALLTOALL_rwp11
  Subroutine MPI_ALLTOALL_rwp22(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: s_a, MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,       Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: i, j, k

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If

    k = 0
    X: Do i = 1, Ubound(aaa, Dim=2)
      Do j = 1, Ubound(bbb, Dim=2)
        bbb(j, i) = aaa(j, i)
        k = k + 1
        If (k == s_a) Then
          Exit X
        End If
      End Do
    End Do X

  End Subroutine MPI_ALLTOALL_rwp22
  Subroutine MPI_ALLTOALL_cwp11(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: s_a, MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If
    bbb(1:s_a) = aaa(1:s_a)

  End Subroutine MPI_ALLTOALL_cwp11
  Subroutine MPI_ALLTOALL_cwp22(aaa, s_a, MPI_WPa, bbb, s_b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: s_a, MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,          Intent(In   ) :: s_b, MPI_WPb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: i, j, k

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s_a < 1 .or. s_a > Size(aaa) .or. &
        s_b < 1 .or. s_b > Size(bbb) .or. s_a /= s_b) Then
      ierr = 1
      Stop
    End If

    k = 0
    X: Do i = 1, Ubound(aaa, Dim=2)
      Do j = 1, Ubound(bbb, Dim=2)
        bbb(j, i) = aaa(j, i)
        k = k + 1
        If (k == s_a) Then
          Exit X
        End If
      End Do
    End Do X

  End Subroutine MPI_ALLTOALL_cwp22

  Subroutine MPI_ALLTOALLV_log11(aaa, aa, a, MPI_LOGICALa, bbb, bb, b, MPI_LOGICALb, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: aa(:), a(:), MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: bb(:), b(:), MPI_LOGICALb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_log11
  Subroutine MPI_ALLTOALLV_log22(aaa, aa, a, MPI_LOGICALa, bbb, bb, b, MPI_LOGICALb, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:, :)
    Integer, Intent(In   ) :: aa(:), a(:), MPI_LOGICALa
    Logical, Intent(InOut) :: bbb(:, :)
    Integer, Intent(In   ) :: bb(:), b(:), MPI_LOGICALb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_log22
  Subroutine MPI_ALLTOALLV_chr11(aaa, aa, a, MPI_CHARACTERa, bbb, bb, b, MPI_CHARACTERb, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: aa(:), a(:), MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: bb(:), b(:), MPI_CHARACTERb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_chr11
  Subroutine MPI_ALLTOALLV_chr22(aaa, aa, a, MPI_CHARACTERa, bbb, bb, b, MPI_CHARACTERb, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: aa(:), a(:), MPI_CHARACTERa
    Character(Len=*), Intent(InOut) :: bbb(:, :)
    Integer,          Intent(In   ) :: bb(:), b(:), MPI_CHARACTERb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_chr22
  Subroutine MPI_ALLTOALLV_int11(aaa, aa, a, MPI_INTEGERa, bbb, bb, b, MPI_INTEGERb, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), aa(:), a(:), MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:)
    Integer, Intent(In   ) :: bb(:), b(:), MPI_INTEGERb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_int11
  Subroutine MPI_ALLTOALLV_int22(aaa, aa, a, MPI_INTEGERa, bbb, bb, b, MPI_INTEGERb, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:, :), aa(:), a(:), MPI_INTEGERa
    Integer, Intent(InOut) :: bbb(:, :)
    Integer, Intent(In   ) :: bb(:), b(:), MPI_INTEGERb, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_int22
  Subroutine MPI_ALLTOALLV_rwp11(aaa, aa, a, MPI_WPa, bbb, bb, b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: aa(:), a(:), MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,       Intent(In   ) :: bb(:), b(:), MPI_WPb, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_rwp11
  Subroutine MPI_ALLTOALLV_rwp22(aaa, aa, a, MPI_WPa, bbb, bb, b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: aa(:), a(:), MPI_WPa
    Real(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,       Intent(In   ) :: bb(:), b(:), MPI_WPb, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_rwp22
  Subroutine MPI_ALLTOALLV_cwp11(aaa, aa, a, MPI_WPa, bbb, bb, b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: aa(:), a(:), MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:)
    Integer,          Intent(In   ) :: bb(:), b(:), MPI_WPb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_cwp11
  Subroutine MPI_ALLTOALLV_cwp22(aaa, aa, a, MPI_WPa, bbb, bb, b, MPI_WPb, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: aa(:), a(:), MPI_WPa
    Complex(Kind=wp), Intent(InOut) :: bbb(:, :)
    Integer,          Intent(In   ) :: bb(:), b(:), MPI_WPb, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. Size(aaa) /= Size(bbb) .or. &
        Size(aa) /= Size(bb) .or. Size(a) /= Size(b)) Then
      ierr = 1
      Stop
    End If
    bbb = aaa

  End Subroutine MPI_ALLTOALLV_cwp22

  Subroutine MPI_SEND_log_s(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_log_s
  Subroutine MPI_SEND_log_v(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_log_v
  Subroutine MPI_SEND_chr_s(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_chr_s
  Subroutine MPI_SEND_chr_v(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Size(aaa) * Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_chr_v
  Subroutine MPI_SEND_int_s(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_int_s
  Subroutine MPI_SEND_int_v(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_int_v
  Subroutine MPI_SEND_rwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_rwp_s
  Subroutine MPI_SEND_rwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_rwp_v
  Subroutine MPI_SEND_rwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_rwp_m
  Subroutine MPI_SEND_rwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_rwp_c
  Subroutine MPI_SEND_cwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 2) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_cwp_s
  Subroutine MPI_SEND_cwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_cwp_v
  Subroutine MPI_SEND_cwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_cwp_m
  Subroutine MPI_SEND_cwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_SEND_cwp_c
  Subroutine MPI_SEND_part_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Type(corePart), Intent(In   ) :: aaa
    Integer,        Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,        Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0) Then
      ierr = 1
      Stop
    End If
  End Subroutine MPI_SEND_part_s
  Subroutine MPI_SEND_part_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, ierr)

    Type(corePart), Intent(In   ) :: aaa(:)
    Integer,        Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,        Intent(  Out) :: ierr

    ierr = 0
    If (idnode /= 0) Then
      ierr = 1
      Stop
    End If
  End Subroutine MPI_SEND_part_v

  Subroutine MPI_ISEND_log_s(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_log_s
  Subroutine MPI_ISEND_log_v(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_log_v
  Subroutine MPI_ISEND_chr_s(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_chr_s
  Subroutine MPI_ISEND_chr_v(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa) * Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_chr_v
  Subroutine MPI_ISEND_int_s(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Integer, Intent(In   ) :: aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_int_s
  Subroutine MPI_ISEND_int_v(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Integer, Intent(In   ) :: aaa(:), n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_int_v
  Subroutine MPI_ISEND_rwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_rwp_s
  Subroutine MPI_ISEND_rwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_rwp_v
  Subroutine MPI_ISEND_rwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_rwp_m
  Subroutine MPI_ISEND_rwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_rwp_c
  Subroutine MPI_ISEND_rwp_f(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_rwp_f
  Subroutine MPI_ISEND_cwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_cwp_s
  Subroutine MPI_ISEND_cwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_cwp_v
  Subroutine MPI_ISEND_cwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_cwp_m
  Subroutine MPI_ISEND_cwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISEND_cwp_c

  Subroutine MPI_ISSEND_log_s(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_log_s
  Subroutine MPI_ISSEND_log_v(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_log_v
  Subroutine MPI_ISSEND_chr_s(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_chr_s
  Subroutine MPI_ISSEND_chr_v(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa) * Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_chr_v
  Subroutine MPI_ISSEND_int_s(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Integer, Intent(In   ) :: aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_int_s
  Subroutine MPI_ISSEND_int_v(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Integer, Intent(In   ) :: aaa(:), n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_int_v
  Subroutine MPI_ISSEND_rwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_rwp_s
  Subroutine MPI_ISSEND_rwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_rwp_v
  Subroutine MPI_ISSEND_rwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_rwp_m
  Subroutine MPI_ISSEND_rwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_rwp_c
  Subroutine MPI_ISSEND_cwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_cwp_s
  Subroutine MPI_ISSEND_cwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_cwp_v
  Subroutine MPI_ISSEND_cwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_cwp_m
  Subroutine MPI_ISSEND_cwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:, :, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_ISSEND_cwp_c

  Subroutine MPI_RECV_log_s(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Logical, Intent(InOut) :: aaa
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_log_s
  Subroutine MPI_RECV_log_v(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Logical, Intent(InOut) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_log_v
  Subroutine MPI_RECV_chr_s(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Character(Len=*), Intent(InOut) :: aaa
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_chr_s
  Subroutine MPI_RECV_chr_v(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Character(Len=*), Intent(InOut) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Size(aaa) * Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_chr_v
  Subroutine MPI_RECV_int_s(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Integer, Intent(InOut) :: aaa
    Integer, Intent(In   ) :: n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_int_s
  Subroutine MPI_RECV_int_v(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Integer, Intent(InOut) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_int_v
  Subroutine MPI_RECV_rwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_rwp_s
  Subroutine MPI_RECV_rwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_rwp_v
  Subroutine MPI_RECV_rwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_rwp_m
  Subroutine MPI_RECV_rwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_rwp_c
  Subroutine MPI_RECV_cwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 2) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_cwp_s
  Subroutine MPI_RECV_cwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_cwp_v
  Subroutine MPI_RECV_cwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa(:, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_cwp_m
  Subroutine MPI_RECV_cwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, status, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa(:, :, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_RECV_cwp_c
  Subroutine MPI_RECV_part_s(aaa, n, part, idnode, tag, MPI_COMM_WORLD, status, ierr)
    Type(corePart), Intent(InOut) :: aaa
    Integer,        Intent(In   ) :: n, part, idnode, tag, MPI_COMM_WORLD
    Integer,        Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0) Then
      ierr = 1
      Stop
    End If
  End Subroutine MPI_RECV_part_s
  Subroutine MPI_RECV_part_v(aaa, n, part, idnode, tag, MPI_COMM_WORLD, status, ierr)
    Type(corePart), Intent(InOut) :: aaa(:)
    Integer,        Intent(In   ) :: n, part, idnode, tag, MPI_COMM_WORLD
    Integer,        Intent(  Out) :: status(:), ierr

    ierr = 0
    status = 0
    If (idnode /= 0) Then
      ierr = 1
      Stop
    End If
  End Subroutine MPI_RECV_part_v

  Subroutine MPI_IRECV_log_s(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Logical, Intent(InOut) :: aaa
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_log_s
  Subroutine MPI_IRECV_log_v(aaa, n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Logical, Intent(InOut) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_LOGICAL, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_log_v
  Subroutine MPI_IRECV_chr_s(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Character(Len=*), Intent(InOut) :: aaa
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_chr_s
  Subroutine MPI_IRECV_chr_v(aaa, n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Character(Len=*), Intent(InOut) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_CHARACTER, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa) * Len(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_chr_v
  Subroutine MPI_IRECV_int_s(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Integer, Intent(InOut) :: aaa
    Integer, Intent(In   ) :: n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_int_s
  Subroutine MPI_IRECV_int_v(aaa, n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Integer, Intent(InOut) :: aaa(:)
    Integer, Intent(In   ) :: n, MPI_INTEGER, idnode, tag, MPI_COMM_WORLD
    Integer, Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_int_v
  Subroutine MPI_IRECV_rwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 1) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_rwp_s
  Subroutine MPI_IRECV_rwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_rwp_v
  Subroutine MPI_IRECV_rwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_rwp_m
  Subroutine MPI_IRECV_rwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_rwp_c
  Subroutine MPI_IRECV_rwp_f(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Real(Kind=wp), Intent(InOut) :: aaa(:, :, :, :)
    Integer,       Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_rwp_f
  Subroutine MPI_IRECV_cwp_s(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_cwp_s
  Subroutine MPI_IRECV_cwp_v(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa(:)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_cwp_v
  Subroutine MPI_IRECV_cwp_m(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa(:, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_cwp_m
  Subroutine MPI_IRECV_cwp_c(aaa, n, MPI_WP, idnode, tag, MPI_COMM_WORLD, request, ierr)

    Complex(Kind=wp), Intent(InOut) :: aaa(:, :, :)
    Integer,          Intent(In   ) :: n, MPI_WP, idnode, tag, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: request, ierr

    ierr = 0
    request = 0
    If (idnode /= 0 .or. n /= 2 * Size(aaa)) Then
      ierr = 1
      Stop
    End If

  End Subroutine MPI_IRECV_cwp_c

  Subroutine MPI_SCATTER_log_ss(aaa, na, MPI_LOGICALa, bbb, nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: na, MPI_LOGICALa
    Logical, Intent(  Out) :: bbb
    Integer, Intent(In   ) :: nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. na /= 1 .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa

  End Subroutine MPI_SCATTER_log_ss
  Subroutine MPI_SCATTER_log_sv(aaa, na, MPI_LOGICALa, bbb, nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa
    Integer, Intent(In   ) :: na, MPI_LOGICALa
    Logical, Intent(  Out) :: bbb(:)
    Integer, Intent(In   ) :: nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: bl, bs

    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. na /= 1 .or. na /= nb .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl) = aaa

  End Subroutine MPI_SCATTER_log_sv
  Subroutine MPI_SCATTER_log_vs(aaa, na, MPI_LOGICALa, bbb, nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: na, MPI_LOGICALa
    Logical, Intent(  Out) :: bbb
    Integer, Intent(In   ) :: nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: al, as

    as = Size(aaa); al = Lbound(aaa, Dim=1)

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. na /= 1 .or. na /= nb .or. na > as) Then
      ierr = 1
      Stop
    End If

    bbb = aaa(al)

  End Subroutine MPI_SCATTER_log_vs
  Subroutine MPI_SCATTER_log_vv(aaa, na, MPI_LOGICALa, bbb, nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: aaa(:)
    Integer, Intent(In   ) :: na, MPI_LOGICALa
    Logical, Intent(  Out) :: bbb(:)
    Integer, Intent(In   ) :: nb, MPI_LOGICALb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: al, as, bl, bs

    as = Size(aaa); al = Lbound(aaa, Dim=1)
    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. na /= nb .or. na > as .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl:bl + na - 1) = aaa(al:al + na - 1)

  End Subroutine MPI_SCATTER_log_vv
  Subroutine MPI_SCATTER_chr_ss(aaa, na, MPI_CHARACTERa, bbb, nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: na, MPI_CHARACTERa
    Character(Len=*), Intent(  Out) :: bbb
    Integer,          Intent(In   ) :: nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: la, lb

    la = Len(aaa); lb = Len(bbb)

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. na /= nb .or. la /= lb .or. na /= la) Then
      ierr = 1
      Stop
    End If

    bbb = aaa

  End Subroutine MPI_SCATTER_chr_ss
  Subroutine MPI_SCATTER_chr_sv(aaa, na, MPI_CHARACTERa, bbb, nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: na, MPI_CHARACTERa
    Character(Len=*), Intent(  Out) :: bbb(:)
    Integer,          Intent(In   ) :: nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: bl, bs, la, lb

    la = Len(aaa); lb = Len(bbb)
    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. na /= nb .or. la /= lb .or. na /= la .or. nb > bs * lb) Then
      ierr = 1
      Stop
    End If

    bbb(bl:bl) = aaa

  End Subroutine MPI_SCATTER_chr_sv
  Subroutine MPI_SCATTER_chr_vs(aaa, na, MPI_CHARACTERa, bbb, nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: na, MPI_CHARACTERa
    Character(Len=*), Intent(  Out) :: bbb
    Integer,          Intent(In   ) :: nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: al, as, la, lb

    la = Len(aaa); lb = Len(bbb)
    as = Size(aaa); al = Lbound(aaa, Dim=1)

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. na /= nb .or. na > as * la .or. la /= lb .or. na /= la) Then
      ierr = 1
      Stop
    End If

    bbb = aaa(al)

  End Subroutine MPI_SCATTER_chr_vs
  Subroutine MPI_SCATTER_chr_vv(aaa, na, MPI_CHARACTERa, bbb, nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: na, MPI_CHARACTERa
    Character(Len=*), Intent(  Out) :: bbb(:)
    Integer,          Intent(In   ) :: nb, MPI_CHARACTERb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: al, as, bl, bs, la, lb

    la = Len(aaa); lb = Len(bbb)
    as = Size(aaa); al = Lbound(aaa, Dim=1)
    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. na /= nb .or. na > as * la .or. la /= lb .or. na /= la .or. nb > bs * lb) Then
      ierr = 1
      Stop
    End If

    bbb(bl) = aaa(al)

  End Subroutine MPI_SCATTER_chr_vv
  Subroutine MPI_SCATTER_int_ss(aaa, na, MPI_INTEGERa, bbb, nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa, na, MPI_INTEGERa
    Integer, Intent(  Out) :: bbb
    Integer, Intent(In   ) :: nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. na /= 1 .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa

  End Subroutine MPI_SCATTER_int_ss
  Subroutine MPI_SCATTER_int_sv(aaa, na, MPI_INTEGERa, bbb, nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa, na, MPI_INTEGERa
    Integer, Intent(  Out) :: bbb(:)
    Integer, Intent(In   ) :: nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: bl, bs

    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. na /= 1 .or. na /= nb .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl) = aaa

  End Subroutine MPI_SCATTER_int_sv
  Subroutine MPI_SCATTER_int_vs(aaa, na, MPI_INTEGERa, bbb, nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), na, MPI_INTEGERa
    Integer, Intent(  Out) :: bbb
    Integer, Intent(In   ) :: nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: al, as

    as = Size(aaa); al = Lbound(aaa, Dim=1)

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. na /= 1 .or. na > as .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa(al)

  End Subroutine MPI_SCATTER_int_vs
  Subroutine MPI_SCATTER_int_vv(aaa, na, MPI_INTEGERa, bbb, nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: aaa(:), na, MPI_INTEGERa
    Integer, Intent(  Out) :: bbb(:)
    Integer, Intent(In   ) :: nb, MPI_INTEGERb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: al, as, bl, bs

    as = Size(aaa); al = Lbound(aaa, Dim=1)
    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. na /= nb .or. na > as .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl:bl + na - 1) = aaa(al:al + na - 1)

  End Subroutine MPI_SCATTER_int_vv
  Subroutine MPI_SCATTER_rwp_ss(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: na, MPI_WPa
    Real(Kind=wp), Intent(  Out) :: bbb
    Integer,       Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= 1 .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa

  End Subroutine MPI_SCATTER_rwp_ss
  Subroutine MPI_SCATTER_rwp_sv(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa
    Integer,       Intent(In   ) :: na, MPI_WPa
    Real(Kind=wp), Intent(  Out) :: bbb(:)
    Integer,       Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: bl, bs

    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= 1 .or. na /= nb .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl) = aaa

  End Subroutine MPI_SCATTER_rwp_sv
  Subroutine MPI_SCATTER_rwp_vs(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: na, MPI_WPa
    Real(Kind=wp), Intent(  Out) :: bbb
    Integer,       Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: al, as

    as = Size(aaa); al = Lbound(aaa, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= 1 .or. na > as .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa(al)

  End Subroutine MPI_SCATTER_rwp_vs
  Subroutine MPI_SCATTER_rwp_vv(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,       Intent(In   ) :: na, MPI_WPa
    Real(Kind=wp), Intent(  Out) :: bbb(:)
    Integer,       Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: al, as, bl, bs

    as = Size(aaa); al = Lbound(aaa, Dim=1)
    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= nb .or. na > as .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl:bl + na - 1) = aaa(al:al + na - 1)

  End Subroutine MPI_SCATTER_rwp_vv
  Subroutine MPI_SCATTER_cwp_ss(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: na, MPI_WPa
    Complex(Kind=wp), Intent(  Out) :: bbb
    Integer,          Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= 1 .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa

  End Subroutine MPI_SCATTER_cwp_ss
  Subroutine MPI_SCATTER_cwp_sv(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa
    Integer,          Intent(In   ) :: na, MPI_WPa
    Complex(Kind=wp), Intent(  Out) :: bbb(:)
    Integer,          Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: bl, bs

    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= 1 .or. na /= nb .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl) = aaa

  End Subroutine MPI_SCATTER_cwp_sv
  Subroutine MPI_SCATTER_cwp_vs(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: na, MPI_WPa
    Complex(Kind=wp), Intent(  Out) :: bbb
    Integer,          Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: al, as

    as = Size(aaa); al = Lbound(aaa, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= 1 .or. na > as .or. na /= nb) Then
      ierr = 1
      Stop
    End If

    bbb = aaa(al)

  End Subroutine MPI_SCATTER_cwp_vs
  Subroutine MPI_SCATTER_cwp_vv(aaa, na, MPI_WPa, bbb, nb, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: aaa(:)
    Integer,          Intent(In   ) :: na, MPI_WPa
    Complex(Kind=wp), Intent(  Out) :: bbb(:)
    Integer,          Intent(In   ) :: nb, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: al, as, bl, bs

    as = Size(aaa); al = Lbound(aaa, Dim=1)
    bs = Size(bbb); bl = Lbound(bbb, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. na /= nb .or. na > as .or. nb > bs) Then
      ierr = 1
      Stop
    End If

    bbb(bl:bl + na - 1) = aaa(al:al + na - 1)

  End Subroutine MPI_SCATTER_cwp_vv

  Subroutine MPI_SCATTERV_log_vv(isend, iscnt, idisp, MPI_LOGICALa, irecv, ircnt, MPI_LOGICALb, idnode, MPI_COMM_WORLD, ierr)

    Logical, Intent(In   ) :: isend(:)
    Integer, Intent(In   ) :: iscnt(:), idisp(:), MPI_LOGICALa
    Logical, Intent(  Out) :: irecv(:)
    Integer, Intent(In   ) :: ircnt, MPI_LOGICALb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: l_idisp, l_irecv, l_iscnt, l_isend, s_idisp, s_irecv, s_iscnt, s_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); l_isend = Lbound(isend, Dim=1)
    s_irecv = Size(irecv); l_irecv = Lbound(irecv, Dim=1)

    ierr = 0
    If (MPI_LOGICALa /= MPI_LOGICALb .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv < ircnt) Then
      ierr = 1
      Stop
    End If

    irecv(l_irecv:l_irecv + ircnt - 1) = isend(l_isend:l_isend + ircnt - 1)

  End Subroutine MPI_SCATTERV_log_vv
  Subroutine MPI_SCATTERV_chr_vv(isend, iscnt, idisp, MPI_CHARACTERa, irecv, ircnt, MPI_CHARACTERb, idnode, MPI_COMM_WORLD, ierr)

    Character(Len=*), Intent(In   ) :: isend(:)
    Integer,          Intent(In   ) :: iscnt(:), idisp(:), MPI_CHARACTERa
    Character(Len=*), Intent(  Out) :: irecv(:)
    Integer,          Intent(In   ) :: ircnt, MPI_CHARACTERb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: k, l_idisp, l_irecv, l_iscnt, l_isend, len_ir, len_is, s_idisp, s_irecv, s_iscnt, &
               s_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); l_isend = Lbound(isend, Dim=1); len_is = Len(isend)
    s_irecv = Size(irecv); l_irecv = Lbound(irecv, Dim=1); len_ir = Len(irecv)

    ierr = 0
    If (MPI_CHARACTERa /= MPI_CHARACTERb .or. len_is /= len_ir .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend * len_is .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv * len_ir < ircnt) Then
      ierr = 1
      Stop
    End If

    k = ircnt / len_ir - 1
    irecv(l_irecv:l_irecv + k) = isend(l_isend:l_isend + k)

  End Subroutine MPI_SCATTERV_chr_vv
  Subroutine MPI_SCATTERV_int_vv(isend, iscnt, idisp, MPI_INTEGERa, irecv, ircnt, MPI_INTEGERb, idnode, MPI_COMM_WORLD, ierr)

    Integer, Intent(In   ) :: isend(:), iscnt(:), idisp(:), MPI_INTEGERa
    Integer, Intent(  Out) :: irecv(:)
    Integer, Intent(In   ) :: ircnt, MPI_INTEGERb, idnode, MPI_COMM_WORLD
    Integer, Intent(  Out) :: ierr

    Integer :: l_idisp, l_irecv, l_iscnt, l_isend, s_idisp, s_irecv, s_iscnt, s_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); l_isend = Lbound(isend, Dim=1)
    s_irecv = Size(irecv); l_irecv = Lbound(irecv, Dim=1)

    ierr = 0
    If (MPI_INTEGERa /= MPI_INTEGERb .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv < ircnt) Then
      ierr = 1
      Stop
    End If

    irecv(l_irecv:l_irecv + ircnt - 1) = isend(l_isend:l_isend + ircnt - 1)

  End Subroutine MPI_SCATTERV_int_vv
  Subroutine MPI_SCATTERV_rwp_vv(isend, iscnt, idisp, MPI_WPa, irecv, ircnt, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: isend(:)
    Integer,       Intent(In   ) :: iscnt(:), idisp(:), MPI_WPa
    Real(Kind=wp), Intent(  Out) :: irecv(:)
    Integer,       Intent(In   ) :: ircnt, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: l_idisp, l_irecv, l_iscnt, l_isend, s_idisp, s_irecv, s_iscnt, s_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); l_isend = Lbound(isend, Dim=1)
    s_irecv = Size(irecv); l_irecv = Lbound(irecv, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv < ircnt) Then
      ierr = 1
      Stop
    End If

    irecv(l_irecv:l_irecv + ircnt - 1) = isend(l_isend:l_isend + ircnt - 1)

  End Subroutine MPI_SCATTERV_rwp_vv
  Subroutine MPI_SCATTERV_rwp_mm(isend, iscnt, idisp, MPI_WPa, irecv, ircnt, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Real(Kind=wp), Intent(In   ) :: isend(:, :)
    Integer,       Intent(In   ) :: iscnt(:), idisp(:), MPI_WPa
    Real(Kind=wp), Intent(  Out) :: irecv(:, :)
    Integer,       Intent(In   ) :: ircnt, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,       Intent(  Out) :: ierr

    Integer :: k, l1_irecv, l1_isend, l2_irecv, l2_isend, l_idisp, l_iscnt, s1_irecv, s1_isend, &
               s_idisp, s_irecv, s_iscnt, s_isend, u1_irecv, u1_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); s1_isend = Size(isend, Dim=1)
    l1_isend = Lbound(isend, Dim=1); u1_isend = Ubound(isend, Dim=1); l2_isend = Lbound(isend, Dim=2)
    s_irecv = Size(irecv); s1_irecv = Size(irecv, Dim=1)
    l1_irecv = Lbound(irecv, Dim=1); u1_irecv = Ubound(irecv, Dim=1); l2_irecv = Lbound(irecv, Dim=2)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s1_irecv /= s1_isend .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv < ircnt) Then
      ierr = 1
      Stop
    End If

    k = ircnt / s1_irecv
    irecv(l1_irecv:u1_irecv, l2_irecv:l2_irecv + k - 1) = isend(l1_isend:u1_isend, l2_isend:l2_isend + k - 1)

  End Subroutine MPI_SCATTERV_rwp_mm
  Subroutine MPI_SCATTERV_cwp_vv(isend, iscnt, idisp, MPI_WPa, irecv, ircnt, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: isend(:)
    Integer,          Intent(In   ) :: iscnt(:), idisp(:), MPI_WPa
    Complex(Kind=wp), Intent(  Out) :: irecv(:)
    Integer,          Intent(In   ) :: ircnt, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: l_idisp, l_irecv, l_iscnt, l_isend, s_idisp, s_irecv, s_iscnt, s_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); l_isend = Lbound(isend, Dim=1)
    s_irecv = Size(irecv); l_irecv = Lbound(irecv, Dim=1)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv < ircnt) Then
      ierr = 1
      Stop
    End If

    irecv(l_irecv:l_irecv + ircnt - 1) = isend(l_isend:l_isend + ircnt - 1)

  End Subroutine MPI_SCATTERV_cwp_vv

  Subroutine MPI_SCATTERV_cwp_mm(isend, iscnt, idisp, MPI_WPa, irecv, ircnt, MPI_WPb, idnode, MPI_COMM_WORLD, ierr)

    Complex(Kind=wp), Intent(In   ) :: isend(:, :)
    Integer,          Intent(In   ) :: iscnt(:), idisp(:), MPI_WPa
    Complex(Kind=wp), Intent(  Out) :: irecv(:, :)
    Integer,          Intent(In   ) :: ircnt, MPI_WPb, idnode, MPI_COMM_WORLD
    Integer,          Intent(  Out) :: ierr

    Integer :: k, l1_irecv, l1_isend, l2_irecv, l2_isend, l_idisp, l_iscnt, s1_irecv, s1_isend, &
               s_idisp, s_irecv, s_iscnt, s_isend, u1_irecv, u1_isend

    s_iscnt = Size(iscnt); l_iscnt = Lbound(iscnt, Dim=1)
    s_idisp = Size(idisp); l_idisp = Lbound(idisp, Dim=1)
    s_isend = Size(isend); s1_isend = Size(isend, Dim=1)
    l1_isend = Lbound(isend, Dim=1); u1_isend = Ubound(isend, Dim=1); l2_isend = Lbound(isend, Dim=2)
    s_irecv = Size(irecv); s1_irecv = Size(irecv, Dim=1)
    l1_irecv = Lbound(irecv, Dim=1); u1_irecv = Ubound(irecv, Dim=1); l2_irecv = Lbound(irecv, Dim=2)

    ierr = 0
    If (MPI_WPa /= MPI_WPb .or. s1_irecv /= s1_isend .or. &
        s_idisp /= 1 .or. idisp(l_idisp) /= 0 .or. &
        s_iscnt /= 1 .or. iscnt(l_iscnt) > s_isend .or. &
        iscnt(l_iscnt) /= ircnt .or. s_irecv < ircnt) Then
      ierr = 1
      Stop
    End If

    k = ircnt / s1_irecv
    irecv(l1_irecv:u1_irecv, l2_irecv:l2_irecv + k - 1) = isend(l1_isend:u1_isend, l2_isend:l2_isend + k - 1)

  End Subroutine MPI_SCATTERV_cwp_mm

  Subroutine MPI_TYPE_GET_EXTENT(size_in, mpi_type_out1, mpi_type_out2, ierr)

    Integer,                        Intent(In   ) :: size_in
    Integer(Kind=MPI_ADDRESS_KIND), Intent(  Out) :: mpi_type_out1, mpi_type_out2
    Integer,                        Intent(  Out) :: ierr

    ierr = 0
    mpi_type_out1 = size_in
    mpi_type_out2 = size_in

  End Subroutine MPI_TYPE_GET_EXTENT

  Subroutine MPI_TYPE_EXTENT(size_in, mpi_type_out, ierr)

    Integer,                        Intent(In   ) :: size_in
    Integer(Kind=MPI_ADDRESS_KIND), Intent(  Out) :: mpi_type_out
    Integer,                        Intent(  Out) :: ierr

    ierr = 0
    mpi_type_out = size_in

  End Subroutine MPI_TYPE_EXTENT

  Subroutine MPI_TYPE_CONTIGUOUS(size_in, mpi_type_in, mpi_type_out, ierr)

    Integer, Intent(In   ) :: size_in, mpi_type_in
    Integer, Intent(  Out) :: mpi_type_out, ierr

    ierr = 0
    mpi_type_out = size_in * mpi_type_in

  End Subroutine MPI_TYPE_CONTIGUOUS

  Subroutine MPI_TYPE_COMMIT(mpi_type, ierr)

    Integer, Intent(In   ) :: mpi_type
    Integer, Intent(  Out) :: ierr

    ierr = 0

    ! Bookkeeping

    mpi_io_ftype(mpi_io_cnt) = mpi_type

  End Subroutine MPI_TYPE_COMMIT

  Subroutine MPI_TYPE_VECTOR(nnn, block_length, stride, mpi_type_in, mpi_type_out, ierr)

    Integer, Intent(In   ) :: nnn, block_length, stride, mpi_type_in
    Integer, Intent(  Out) :: mpi_type_out, ierr

    ierr = 0
    mpi_type_out = mpi_type_in * nnn

  End Subroutine MPI_TYPE_VECTOR

  Subroutine MPI_TYPE_CREATE_HVECTOR(nnn, block_length, stride, mpi_type_in, mpi_type_out, ierr)

    Integer,                        Intent(In   ) :: nnn, block_length
    Integer(Kind=MPI_ADDRESS_KIND), Intent(In   ) :: stride
    Integer,                        Intent(In   ) :: mpi_type_in
    Integer,                        Intent(  Out) :: mpi_type_out, ierr

    ierr = 0
    mpi_type_out = mpi_type_in * nnn

  End Subroutine MPI_TYPE_CREATE_HVECTOR

  Subroutine MPI_TYPE_HVECTOR(nnn, block_length, stride, mpi_type_in, mpi_type_out, ierr)

    Integer,                        Intent(In   ) :: nnn, block_length
    Integer(Kind=MPI_ADDRESS_KIND), Intent(In   ) :: stride
    Integer,                        Intent(In   ) :: mpi_type_in
    Integer,                        Intent(  Out) :: mpi_type_out, ierr

    ierr = 0
    mpi_type_out = mpi_type_in * nnn

  End Subroutine MPI_TYPE_HVECTOR

  Subroutine MPI_TYPE_FREE(mpi_type, ierr)

    Integer, Intent(In   ) :: mpi_type
    Integer, Intent(  Out) :: ierr

    ierr = 0

  End Subroutine MPI_TYPE_FREE

  Subroutine MPI_FILE_DELETE(file_name, mpi_info, ierr)

    Character(Len=*), Intent(In   ) :: file_name
    Integer,          Intent(InOut) :: mpi_info
    Integer,          Intent(  Out) :: ierr

    Integer :: file_unit = 30
    Logical :: lexist

    ierr = 0
    Inquire (File=file_name, Exist=lexist)
    If (.not. lexist) Return

    10 Continue
    Inquire (Unit=file_unit, Opened=lexist)
    If (lexist) Then
      file_unit = file_unit + 1
      Go To 10
    End If
    Open (Unit=file_unit, File=file_name, Status="replace")
    Close (Unit=file_unit)

  End Subroutine MPI_FILE_DELETE

  Subroutine MPI_FILE_OPEN(MPI_COMM_WORLD, file_name, mpi_mode, mpi_info, fh, ierr)

    Integer,          Intent(In   ) :: MPI_COMM_WORLD
    Character(Len=*), Intent(In   ) :: file_name
    Integer,          Intent(In   ) :: mpi_mode
    Integer,          Intent(InOut) :: mpi_info
    Integer,          Intent(  Out) :: fh, ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize) :: a_ccess, a_ction, f_orm, p_osition, s_tatus, word
    Integer              :: i, wlen
    Logical              :: lexist

    ! Defaults

    fh = 30
    ierr = 0

    Do i = 1, mpi_io_cnt - 1
      If (mpi_io_fh(i) == fh) fh = fh + 1
    End Do

    ! Bookkeeping

    mpi_io_comm(mpi_io_cnt) = MPI_COMM_WORLD
    mpi_io_fname(mpi_io_cnt) = file_name
    mpi_io_mode(mpi_io_cnt) = mpi_mode
    mpi_io_fh(mpi_io_cnt) = fh

    mpi_io_cnt = mpi_io_cnt + 1

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_mode
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1' .or. word(wlen:wlen) == '3') Then
      Inquire (File=file_name, Exist=lexist)
      If (lexist) Then
        s_tatus = 'old'
        If (word(wlen:wlen) == '3') Then
          ierr = 1
          Stop
        End If
      Else
        s_tatus = 'new'
      End If
    Else If (word(wlen:wlen) /= '0' .and. word(wlen:wlen) /= '2') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 3
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 6
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 8
      Stop
    End If

    ! Check assumption

    Inquire (Unit=fh, Opened=lexist)
    If (lexist) Then
      ierr = 9
      Stop
    End If

  End Subroutine MPI_FILE_OPEN

  Subroutine MPI_FILE_SET_VIEW(fh, disp, etype, filetype, datarep, mpi_info, ierr)

    Integer,                       Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND), Intent(In   ) :: disp
    Integer,                       Intent(In   ) :: etype, filetype
    Character(Len=*),              Intent(In   ) :: datarep
    Integer,                       Intent(InOut) :: mpi_info
    Integer,                       Intent(  Out) :: ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize) :: a_ccess, a_ction, f_orm, p_osition, s_tatus, word
    Integer              :: cnt, i, wlen
    Logical              :: lexist

    ! Defaults

    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    If (mpi_io_ftype(cnt) /= filetype) Then
      ierr = 2
      Stop
    End If

    mpi_io_etype(cnt) = etype
    mpi_io_disp(cnt) = disp
    mpi_io_rec_len(cnt) = (filetype / etype) * (etype / 1) !ALL MPI_PRIMITIVE_TYPES=1 in MPIF.H
    mpi_io_datarep(cnt) = datarep

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 3
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    Inquire (Unit=fh, Opened=lexist)
    If (lexist) Close (Unit=fh)

    Open (Unit=fh, File=mpi_io_fname(cnt), Status=s_tatus, Action=a_ction, Access=a_ccess, Form=f_orm, Recl=mpi_io_rec_len(cnt))

  End Subroutine MPI_FILE_SET_VIEW

  Subroutine MPI_FILE_GET_VIEW(fh, disp, etype, filetype, datarep, ierr)

    Integer,                       Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND), Intent(  Out) :: disp
    Integer,                       Intent(  Out) :: etype, filetype
    Character(Len=*),              Intent(  Out) :: datarep
    Integer,                       Intent(  Out) :: ierr

    Integer, Parameter :: wsize = 12

    Integer :: cnt, i

    ! Defaults

    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    disp = mpi_io_disp(cnt)
    etype = mpi_io_etype(cnt)
    filetype = mpi_io_ftype(cnt)
    datarep = mpi_io_datarep(cnt)

  End Subroutine MPI_FILE_GET_VIEW

  Subroutine MPI_GET_COUNT(status, etype, count, ierr)

    Integer, Intent(In   ) :: status(:), etype
    Integer, Intent(  Out) :: count, ierr

    ! Defaults

    ierr = 0
    count = status(1) / etype

  End Subroutine MPI_GET_COUNT

  Subroutine MPI_FILE_CLOSE(fh, ierr)

    Integer, Intent(In   ) :: fh
    Integer, Intent(  Out) :: ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize) :: a_ccess, a_ction, f_orm, p_osition, s_tatus, word
    Integer              :: cnt, i, wlen
    Logical              :: ldelete = .false., lexist

    ! Defaults

    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 3
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
      ldelete = .true.
      s_tatus = 'replace'
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 6
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 8
      Stop
    End If

    Inquire (Unit=fh, Opened=lexist)
    If (.not. lexist) Then
      ierr = 9
      Stop
    End If

    Close (Unit=fh)

    If (ldelete) Then
      Open (Unit=fh, File=mpi_io_fname(cnt), Status=s_tatus)
      Close (Unit=fh)
    End If

    ! Initialise bookkeeping

    If (cnt == mpi_io_cnt - 1) Then
      mpi_io_rec_len(cnt) = 0
      mpi_io_etype(cnt) = 0
      mpi_io_ftype(cnt) = 0
      mpi_io_comm(cnt) = 0
      mpi_io_fname(cnt) = ' '
      mpi_io_mode(cnt) = 0
      mpi_io_fh(cnt) = 0
      mpi_io_disp(cnt) = 0_MPI_OFFSET_KIND
      mpi_io_datarep(cnt) = ' '
    Else
      mpi_io_rec_len(cnt) = mpi_io_rec_len(mpi_io_cnt - 1)
      mpi_io_etype(cnt) = mpi_io_etype(mpi_io_cnt - 1)
      mpi_io_ftype(cnt) = mpi_io_ftype(mpi_io_cnt - 1)
      mpi_io_comm(cnt) = mpi_io_comm(mpi_io_cnt - 1)
      mpi_io_fname(cnt) = mpi_io_fname(mpi_io_cnt - 1)
      mpi_io_mode(cnt) = mpi_io_mode(mpi_io_cnt - 1)
      mpi_io_fh(cnt) = mpi_io_fh(mpi_io_cnt - 1)
      mpi_io_disp(cnt) = mpi_io_disp(mpi_io_cnt - 1)
      mpi_io_datarep(cnt) = mpi_io_datarep(mpi_io_cnt - 1)

      mpi_io_rec_len(mpi_io_cnt - 1) = 0
      mpi_io_etype(mpi_io_cnt - 1) = 0
      mpi_io_ftype(mpi_io_cnt - 1) = 0
      mpi_io_comm(mpi_io_cnt - 1) = 0
      mpi_io_fname(mpi_io_cnt - 1) = ' '
      mpi_io_mode(mpi_io_cnt - 1) = 0
      mpi_io_fh(mpi_io_cnt - 1) = 0
      mpi_io_disp(mpi_io_cnt - 1) = 0_MPI_OFFSET_KIND
      mpi_io_datarep(mpi_io_cnt - 1) = ' '
    End If
    mpi_io_cnt = mpi_io_cnt - 1

  End Subroutine MPI_FILE_CLOSE

  Subroutine MPI_FILE_WRITE_AT_chr_s(fh, offset, record, count, datatype, status, ierr)

    Integer,                       Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND), Intent(In   ) :: offset
    Character(Len=*),              Intent(In   ) :: record
    Integer,                       Intent(In   ) :: count, datatype
    Integer,                       Intent(  Out) :: status(:), ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize)          :: a_ccess, a_ction, f_orm, forma, p_osition, s_tatus, word
    Integer                       :: cnt, i, reclen, wlen
    Integer(Kind=MPI_OFFSET_KIND) :: line

    ! Defaults

    status = 0
    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
      ierr = 3
      Stop
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    ! record size and format type

    reclen = mpi_io_rec_len(cnt)
    forma = ' '; Write (forma, 20) reclen

    line = mpi_io_disp(cnt) + offset + Int(1, MPI_OFFSET_KIND)
    Write (Unit=fh, Fmt=forma, Rec=line, Err=30) record(1:reclen)

    status(1) = count * mpi_io_etype(cnt)

    Return

    20 Format('(a', i0, ')')

    30 Continue
    ierr = 12
    status(1) = -1

  End Subroutine MPI_FILE_WRITE_AT_chr_s
  Subroutine MPI_FILE_WRITE_AT_chr_v_1(fh, offset, record, count, datatype, status, ierr)

    Integer,                        Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND),  Intent(In   ) :: offset
    Character(Len=1), Dimension(:), Intent(In   ) :: record
    Integer,                        Intent(In   ) :: count, datatype
    Integer,                        Intent(  Out) :: status(:), ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize)          :: a_ccess, a_ction, f_orm, forma, p_osition, s_tatus, word
    Integer                       :: cnt, i, reclen, wlen
    Integer(Kind=MPI_OFFSET_KIND) :: line

    ! Defaults

    status = 0
    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
      ierr = 3
      Stop
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    ! record size and format type

    reclen = mpi_io_rec_len(cnt)
    forma = ' '; Write (forma, 20) reclen

    line = mpi_io_disp(cnt) + offset + Int(1, MPI_OFFSET_KIND)
    Write (Unit=fh, Fmt=forma, Rec=line, Err=30) record(1:reclen)

    status(1) = count * mpi_io_etype(cnt)

    Return

    20 Format('(', i0, 'a)')

    30 Continue
    ierr = 12
    status(1) = -1

  End Subroutine MPI_FILE_WRITE_AT_chr_v_1
  Subroutine MPI_FILE_WRITE_AT_chr_m_1(fh, offset, record, count, datatype, status, ierr)

    Integer,                           Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND),     Intent(In   ) :: offset
    Character(Len=1), Dimension(:, :), Intent(In   ) :: record
    Integer,                           Intent(In   ) :: count, datatype
    Integer,                           Intent(  Out) :: status(:), ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize)          :: a_ccess, a_ction, f_orm, forma, p_osition, s_tatus, word
    Integer                       :: cnt, i, reclen, wlen
    Integer(Kind=MPI_OFFSET_KIND) :: line

    ! Defaults

    status = 0
    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
      ierr = 3
      Stop
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    ! record size and format type

    reclen = mpi_io_rec_len(cnt)
    If (Ubound(record, Dim=1) /= reclen) Then
      ierr = 12
      Stop
    End If
    forma = ' '; Write (forma, 20) reclen

    line = mpi_io_disp(cnt) + offset + Int(1, MPI_OFFSET_KIND)
    Write (Unit=fh, Fmt=forma, Rec=line, Err=30) record(1:reclen, 1:count)

    status(1) = count * mpi_io_etype(cnt)

    Return

    20 Format('(', i0, 'a)')

    30 Continue
    ierr = 13
    status(1) = -1

  End Subroutine MPI_FILE_WRITE_AT_chr_m_1

  Subroutine MPI_FILE_READ_AT_chr_s(fh, offset, record, count, datatype, status, ierr)

    Integer,                       Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND), Intent(In   ) :: offset
    Character(Len=*),              Intent(  Out) :: record
    Integer,                       Intent(In   ) :: count, datatype
    Integer,                       Intent(  Out) :: status(:), ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize)          :: a_ccess, a_ction, f_orm, forma, p_osition, s_tatus, word
    Integer                       :: cnt, i, reclen, wlen
    Integer(Kind=MPI_OFFSET_KIND) :: line

    ! Defaults

    status = 0
    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
      ierr = 3
      Stop
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    ! record size and format type

    reclen = mpi_io_rec_len(cnt)
    forma = ' '; Write (forma, 20) reclen

    line = mpi_io_disp(cnt) + offset + Int(1, MPI_OFFSET_KIND)
    Read (Unit=fh, Fmt=forma, Rec=line, Err=30) record(1:reclen)

    status(1) = count * mpi_io_etype(cnt)

    Return

    20 Format('(a', i0, ')')

    30 Continue
    ierr = 12
    status(1) = -1

  End Subroutine MPI_FILE_READ_AT_chr_s
  Subroutine MPI_FILE_READ_AT_chr_v_1(fh, offset, record, count, datatype, status, ierr)

    Integer,                        Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND),  Intent(In   ) :: offset
    Character(Len=1), Dimension(:), Intent(  Out) :: record
    Integer,                        Intent(In   ) :: count, datatype
    Integer,                        Intent(  Out) :: status(:), ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize)          :: a_ccess, a_ction, f_orm, forma, p_osition, s_tatus, word
    Integer                       :: cnt, i, reclen, wlen
    Integer(Kind=MPI_OFFSET_KIND) :: line

    ! Defaults

    status = 0
    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
      ierr = 3
      Stop
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    ! record size and format type

    reclen = mpi_io_rec_len(cnt)
    forma = ' '; Write (forma, 20) reclen

    line = mpi_io_disp(cnt) + offset + Int(1, MPI_OFFSET_KIND)
    Read (Unit=fh, Fmt=forma, Rec=line, Err=30) record(1:reclen)

    status(1) = count * mpi_io_etype(cnt)

    Return

    20 Format('(', i0, 'a)')

    30 Continue
    ierr = 12
    status(1) = -1

  End Subroutine MPI_FILE_READ_AT_chr_v_1
  Subroutine MPI_FILE_READ_AT_chr_m_1(fh, offset, record, count, datatype, status, ierr)

    Integer,                           Intent(In   ) :: fh
    Integer(Kind=MPI_OFFSET_KIND),     Intent(In   ) :: offset
    Character(Len=1), Dimension(:, :), Intent(  Out) :: record
    Integer,                           Intent(In   ) :: count, datatype
    Integer,                           Intent(  Out) :: status(:), ierr

    Integer, Parameter :: wsize = 12

    Character(Len=wsize)          :: a_ccess, a_ction, f_orm, forma, p_osition, s_tatus, word
    Integer                       :: cnt, i, reclen, wlen
    Integer(Kind=MPI_OFFSET_KIND) :: line

    ! Defaults

    status = 0
    ierr = 0

    ! Bookkeeping

    cnt = 0
    Do i = 1, mpi_io_cnt - 1
      If (fh == mpi_io_fh(i)) Then
        cnt = i
        Exit
      End If
    End Do

    If (cnt == 0) Then
      ierr = 1
      Stop
    End If

    a_ction = 'readwrite'
    s_tatus = 'unknown'
    a_ccess = 'direct'
    f_orm = 'formatted'
    p_osition = ' '

    word = ' '
    Write (word, Fmt='(i6.6)') mpi_io_mode(cnt)
    wlen = Len_trim(word)

    ! MPI_MODES MUST BE DEFINED LIKE THIS IN MPIF.H
    !
    ! MPI_MODE_CREATE          =      1
    ! MPI_MODE_EXCL            =      2
    ! MPI_MODE_RDONLY          =     10
    ! MPI_MODE_WRONLY          =     20
    ! MPI_MODE_RDWR            =     30
    ! MPI_MODE_DELETE_ON_CLOSE =    100
    ! MPI_MODE_SEQUENTIAL      =   1000
    ! MPI_MODE_APPEND          =  10000
    ! MPI_MODE_UNIQUE_OPEN     = 100000

    If (word(wlen:wlen) == '1') Then
      s_tatus = 'old'
    Else If (word(wlen:wlen) /= '0') Then
      ierr = 2
      Stop
    End If

    If (word(wlen - 1:wlen - 1) == '1') Then
      a_ction = 'read'
    Else If (word(wlen - 1:wlen - 1) == '2') Then
      a_ction = 'write'
      ierr = 3
      Stop
    Else If (word(wlen - 1:wlen - 1) == '3') Then
      a_ction = 'readwrite'
    Else If (word(wlen - 1:wlen - 1) /= '0') Then
      ierr = 4
      Stop
    End If

    If (word(wlen - 2:wlen - 2) == '1') Then
    Else If (word(wlen - 2:wlen - 2) /= '0') Then
      ierr = 5
      Stop
    End If

    If (word(wlen - 3:wlen - 3) == '1') Then
      a_ccess = 'sequential'
      p_osition = 'asis'
      ierr = 6
      Stop
    Else If (word(wlen - 3:wlen - 3) /= '0') Then
      ierr = 7
      Stop
    End If

    If (word(wlen - 4:wlen - 4) == '1') Then
      p_osition = 'append'
      ierr = 8
      Stop
      If (word(wlen - 3:wlen - 3) /= '1') Then
        ierr = 9
        Stop
      End If
    Else If (word(wlen - 4:wlen - 4) /= '0') Then
      ierr = 10
      Stop
    End If

    If (word(wlen - 5:wlen - 5) == '1') Then
    Else If (word(wlen - 5:wlen - 5) /= '0') Then
      ierr = 11
      Stop
    End If

    ! record size and format type

    reclen = mpi_io_rec_len(cnt)
    If (Ubound(record, Dim=1) /= reclen) Then
      ierr = 12
      Stop
    End If
    forma = ' '; Write (forma, 20) reclen

    line = mpi_io_disp(cnt) + offset + Int(1, MPI_OFFSET_KIND)
    Read (Unit=fh, Fmt=forma, Rec=line, Err=30) record(1:reclen, 1:count)

    status(1) = count * mpi_io_etype(cnt)

    Return

    20 Format('(', i0, 'a)')

    30 Continue
    ierr = 13
    status(1) = -1

  End Subroutine MPI_FILE_READ_AT_chr_m_1

  Subroutine MPI_GET_VERSION(mpi_ver, mpi_subver, ierr)
    Integer, Intent(  Out) :: mpi_ver, mpi_subver, ierr

    mpi_ver = 0
    mpi_subver = 0
    ierr = 0
  End Subroutine MPI_GET_VERSION

  Subroutine MPI_GET_PROCESSOR_NAME(proc_name, lname, ierr)
    Character(Len=MPI_MAX_PROCESSOR_NAME), Intent(  Out) :: proc_name
    Integer,                               Intent(  Out) :: lname, ierr

    proc_name = "*"
    lname = 1
    ierr = 0
  End Subroutine MPI_GET_PROCESSOR_NAME

  Subroutine MPI_GET_LIBRARY_VERSION(lib_version, lversion, ierr)
    Character(Len=MPI_MAX_LIBRARY_VERSION_STRING), Intent(  Out) :: lib_version
    Integer,                                       Intent(  Out) :: lversion, ierr

    lib_version = "*"
    lversion = 1
    ierr = 0
  End Subroutine MPI_GET_LIBRARY_VERSION

  Subroutine MPI_TYPE_CREATE_STRUCT(counts, array_of_blocklengths, array_of_displacements, &
                                    array_of_types, newtype, ierr)

    Integer,                        Intent(In   ) :: counts, array_of_blocklengths(:)
    Integer(kind=MPI_ADDRESS_KIND), Intent(In   ) :: array_of_displacements(:)
    Integer,                        Intent(In   ) :: array_of_types(:)
    Integer,                        Intent(  Out) :: newtype, ierr

    newtype = 1
    ierr = 1

  End Subroutine MPI_TYPE_CREATE_STRUCT

  Subroutine MPI_GET_ADDRESS_int(location, address, ierr)
    Integer,                   Intent(In   ) :: location
    Integer(MPI_ADDRESS_KIND), Intent(  Out) :: address
    Integer,                   Intent(  Out) :: ierr

    address = 1
    ierr = 1
  End Subroutine MPI_GET_ADDRESS_int
  Subroutine MPI_GET_ADDRESS_real(location, address, ierr)
    Real(Kind=wp),             Intent(In   ) :: location
    Integer(MPI_ADDRESS_KIND), Intent(  Out) :: address
    Integer,                   Intent(  Out) :: ierr

    address = 1
    ierr = 1
  End Subroutine MPI_GET_ADDRESS_real
  Subroutine MPI_GET_ADDRESS_class(location, address, ierr)
    Class(*),                  Intent(In   ) :: location
    Integer(MPI_ADDRESS_KIND), Intent(  Out) :: address
    Integer,                   Intent(  Out) :: ierr

    address = 1
    ierr = 1
  End Subroutine MPI_GET_ADDRESS_class

  Subroutine MPI_TYPE_CREATE_RESIZED(oldtype, lb, extent, newtype, ierr)
    Integer,                   Intent(In   ) :: oldtype
    Integer(MPI_ADDRESS_KIND), Intent(In   ) :: lb, extent
    Integer,                   Intent(  Out) :: newtype, ierr

    newtype = 1
    ierr = 1
  End Subroutine MPI_TYPE_CREATE_RESIZED

End Module mpi_api
