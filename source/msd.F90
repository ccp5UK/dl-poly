Module msd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring msd routines variables
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2008
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use comms,           Only: MsdWrite_tag,&
                             comm_self,&
                             comms_type,&
                             gbcast,&
                             gcheck,&
                             grecv,&
                             gsend,&
                             gsum,&
                             gsync,&
                             mode_wronly,&
                             offset_kind
  Use configuration,   Only: configuration_type
  Use constants,       Only: boltz,&
                             zero_plus
  Use errors_warnings, Only: error
  Use filename,        Only: FILE_MSD,&
                             file_type
  Use flow_control,    Only: RESTART_KEY_OLD
  Use io,              Only: &
                             IO_ALLOCATION_ERROR, IO_BASE_COMM_NOT_SET, IO_MSDTMP, &
                             IO_UNKNOWN_WRITE_LEVEL, IO_UNKNOWN_WRITE_OPTION, &
                             IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
                             IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT, &
                             IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, io_close, &
                             io_finalize, io_get_parameters, io_init, io_open, io_set_parameters, &
                             io_type, io_write_batch, io_write_record, io_write_sorted_file
  Use kinds,           Only: li,&
                             wi,&
                             wp
  Use parse,           Only: get_word,&
                             tabs_2_blanks,&
                             word_2_real

  Implicit None

  Private

  Type, Public :: msd_type
    Private

    !> MSD recording switch
    Logical, Public          :: l_msd = .false.
    !> Step to begin recording MSD
    Integer(Kind=wi), Public :: start
    !> Frequency to record MSD (steps)
    Integer(Kind=wi), Public :: freq
    Logical                  :: newjob = .true., &
                                fast = .true.
    Character(Len=1024)        :: fname
    Integer(Kind=li)         :: rec = 0_li, &
                                frm = 0_li
  End Type msd_type

  Public :: msd_write
Contains

  Subroutine msd_write(config, keyres, megatm, nstep, tstep, time, stpval, dof_site, io, msd_data, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for writing MSDTMP file at selected intervals
    ! in simulation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! contrib   - i.j.bush
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type),    Intent(InOut) :: config
    Integer,                     Intent(In   ) :: keyres, megatm, nstep
    Real(Kind=wp),               Intent(In   ) :: tstep, time
    Real(Kind=wp),               Intent(InOut) :: stpval(0:)
    Real(Kind=wp), Dimension(:), Intent(In   ) :: dof_site
    Type(io_type),               Intent(InOut) :: io
    Type(msd_type),              Intent(Inout) :: msd_data
    Type(file_type),             Intent(InOut) :: files(:)
    Type(comms_type),            Intent(InOut) :: comm

    Integer, Parameter :: recsz = 53

    Character                                      :: lf
    Character(Len=1), Allocatable, Dimension(:, :) :: chbat
    Character(Len=256)                             :: message
    Character(Len=40)                              :: word
    Character(Len=8), Allocatable, Dimension(:)    :: chbuf
    Character(Len=recsz)                           :: record
    Integer                                        :: batsz, fail(1:2), fh, i, ierr, io_write, &
                                                      jatms, jdnode, jj, k
    Integer(Kind=offset_kind)                      :: rec_mpi_io
    Integer, Allocatable, Dimension(:)             :: iwrk, n_atm
    Logical                                        :: lexist, ready, safe
    Real(Kind=wp)                                  :: buffer(1:2), tmp
    Real(Kind=wp), Allocatable, Dimension(:)       :: ddd, eee

! default record size
! Some parameters and variables needed by io interfaces

    If (.not. (nstep >= msd_data%start .and. Mod(nstep - msd_data%start, msd_data%freq) == 0)) Return

    ! Get write buffer size and line feed character

    Call io_get_parameters(io, user_method_write=io_write)
    Call io_get_parameters(io, user_buffer_size_write=batsz)
    Call io_get_parameters(io, user_line_feed=lf)

    If (msd_data%newjob) Then
      msd_data%newjob = .false.

      ! name convention

      msd_data%fname = files(FILE_MSD)%filename

      ! If keyres = RESTART_KEY_OLD, is MSDTMP old (does it exist) and
      ! how many frames and records are in there

      lexist = .true.
      If (keyres == RESTART_KEY_OLD) Then
        If (comm%idnode == 0) Inquire (File=msd_data%fname, Exist=lexist)
        Call gcheck(comm, lexist, "enforce")
      Else
        lexist = .false.
      End If

      ! Generate file is non-existent

      10 Continue
      If (.not. lexist) Then

        If (comm%idnode == 0) Then
          Open (Newunit=files(FILE_MSD)%unit_no, File=msd_data%fname, &
                Form='formatted', Access='direct', Status='replace', Recl=recsz)
          Write (Unit=files(FILE_MSD)%unit_no, Fmt='(a52,a1)', Rec=Int(1, li)) config%cfgname(1:52), lf
          Write (Unit=files(FILE_MSD)%unit_no, Fmt='(i10,2i21,a1)', Rec=Int(2, li)) megatm, msd_data%frm, msd_data%rec, lf
          Call files(FILE_MSD)%close ()
        End If
        msd_data%rec = Int(2, li)
        msd_data%frm = Int(0, li)

        ! Get some sense of it

      Else

        safe = .true.
        If (comm%idnode == 0) Then

          Open (Newunit=files(FILE_MSD)%unit_no, File=msd_data%fname, Form='formatted')

          Do

            record = ' '

            ! Assume new style of MSDTMP with bookkeeping.

            If (msd_data%fast) Then

              Read (Unit=files(FILE_MSD)%unit_no, Fmt=*, End=10) ! title record
              msd_data%rec = msd_data%rec + Int(1, li)
              record = ' '
              Read (Unit=files(FILE_MSD)%unit_no, Fmt='(a)', End=10) record ! bookkeeping record
              Call tabs_2_blanks(record)
              msd_data%rec = msd_data%rec + Int(1, li)

              Call get_word(record, word)
              If (word(1:Len_trim(word)) /= 'timestep') Then
                Call get_word(record, word); Call get_word(record, word)
                Call get_word(record, word); msd_data%frm = Nint(word_2_real(word, 0.0_wp), li)
                Call get_word(record, word); msd_data%rec = Nint(word_2_real(word, 0.0_wp), li)
                If (msd_data%frm /= Int(0, li) .and. msd_data%rec > Int(2, li)) Then
                  Go To 20 ! New style
                Else
                  msd_data%fast = .false. ! TOUGH, old style
                  msd_data%rec = Int(2, li)
                  msd_data%frm = Int(0, li)
                End If
              Else
                safe = .false. ! Overwrite the file, it's junk to me
                Go To 20
              End If

              ! TOUGH, it needs scanning through

            Else

              Read (Unit=files(FILE_MSD)%unit_no, Fmt='(a)', End=20) record ! timestep record
              Call tabs_2_blanks(record)
              msd_data%rec = msd_data%rec + Int(1, li)

              Call get_word(record, word)
              Call get_word(record, word); jj = Nint(word_2_real(word)) ! total number of lines to read

              word = ' '
              Write (word, '( "(", i0, "( / ) )" )') jj - 1
              Read (Unit=files(FILE_MSD)%unit_no, Fmt=word, End=20)
              msd_data%rec = msd_data%rec + Int(jj, li)
              msd_data%frm = msd_data%frm + Int(1, li)

            End If

          End Do

          20 Continue
          Call files(FILE_MSD)%close ()

        End If

        Call gcheck(comm, safe, "enforce")
        If (.not. safe) Then
          lexist = .false.

          msd_data%rec = Int(0, li)
          msd_data%frm = Int(0, li)

          Go To 10
        Else If (comm%mxnode > 1) Then
          buffer(1) = Real(msd_data%frm, wp)
          buffer(2) = Real(msd_data%rec, wp)

          Call gbcast(comm, buffer(1:2), 9)

          msd_data%frm = Nint(buffer(1), li)
          msd_data%rec = Nint(buffer(2), li)
        End If

      End If
    End If

    ! Get offsets and define batch

    fail = 0
    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT .or. &
        io_write == IO_WRITE_UNSORTED_MASTER) Then
      Allocate (n_atm(0:comm%mxnode), Stat=fail(1))
      Allocate (chbat(1:recsz, 1:batsz), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a,i0)') 'msd_write allocation failure 0'
        Call error(0, message)
      End If

      chbat = ' '
      n_atm = 0; n_atm(comm%idnode + 1) = config%natms
      Call gsum(comm, n_atm)
      n_atm(0) = Sum(n_atm(0:comm%idnode))
    End If

    ! Notes:
    ! the MPI-I/O records are numbered from 0 (not 1)
    ! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
    !   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
    !   both declared as: Integer(Kind = offset_kind)

    ! Update frame

    msd_data%frm = msd_data%frm + Int(1, li)

    ! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT) Then

      ! Write header information, where just one node is needed
      ! Start of file

      rec_mpi_io = Int(msd_data%rec, offset_kind)
      If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm=comm_self)
        Call io_init(io, recsz)
        Call io_open(io, io_write, comm_self, msd_data%fname, mode_wronly, fh)

        Write (record, Fmt='(a8,2i10,2f12.6,a1)') 'timestep', nstep, megatm, tstep, time, lf

        ! Dump header information

        Call io_write_record(io, fh, rec_mpi_io, record)

        Call io_close(io, fh)
        Call io_finalize(io)

      End If

      ! Start of file

      msd_data%rec = msd_data%rec + Int(1, li)
      rec_mpi_io = Int(msd_data%rec, offset_kind) + Int(n_atm(0), offset_kind)
      jj = 0

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_init(io, recsz)
      Call io_open(io, io_write, comm%comm, msd_data%fname, mode_wronly, fh)

      Do i = 1, config%natms
        k = 2 * i

        If (Abs(dof_site(config%lsite(i))) > zero_plus) tmp = config%weight(i) * stpval(36 + k) / (boltz * 3.0_wp)

        Write (record, Fmt='(a8,i10,1p,2e13.4,a8,a1)') config%atmnam(i), &
          config%ltg(i), Sqrt(stpval(36 + k - 1)), tmp, Repeat(' ', 8), lf
        jj = jj + 1
        Do k = 1, recsz
          chbat(k, jj) = record(k:k)
        End Do

        ! Dump batch and update start of file

        If (jj >= batsz .or. i == config%natms) Then
          Call io_write_batch(io, fh, rec_mpi_io, jj, chbat)
          rec_mpi_io = rec_mpi_io + Int(jj, offset_kind)
          jj = 0
        End If
      End Do

      ! Update and save offset pointer

      msd_data%rec = msd_data%rec + Int(megatm, li)
      If (comm%idnode == 0) Then
        Write (record, Fmt='(i10,2i21,a1)') megatm, msd_data%frm, msd_data%rec, lf
        Call io_write_record(io, fh, Int(1, offset_kind), record)
      End If

      Call io_close(io, fh)
      Call io_finalize(io)

      ! UNSORTED Serial Direct Access FORTRAN

    Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

      Allocate (chbuf(1:config%mxatms), iwrk(1:config%mxatms), ddd(1:config%mxatms), eee(1:config%mxatms), Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'msd_write allocation failure'
        Call error(0, message)
      End If

      ! node 0 handles I/O

      If (comm%idnode == 0) Then

        Open (Newunit=files(FILE_MSD)%unit_no, File=msd_data%fname, Form='formatted', Access='direct', Recl=recsz)

        msd_data%rec = msd_data%rec + Int(1, li)
        Write (Unit=files(FILE_MSD)%unit_no, Fmt='(a8,2i10,2f12.6,a1)', Rec=msd_data%rec) 'timestep', nstep, &
          megatm, tstep, time, lf

        Do i = 1, config%natms
          k = 2 * i

          iwrk(i) = config%ltg(i)
          chbuf(i) = config%atmnam(i)

          ddd(i) = Sqrt(stpval(36 + k - 1))

          If (Abs(dof_site(config%lsite(i))) > zero_plus) tmp = config%weight(i) * stpval(36 + k) / (boltz * 3.0_wp)
          eee(i) = tmp
        End Do

        jatms = config%natms
        ready = .true.
        Do jdnode = 0, comm%mxnode - 1
          If (jdnode > 0) Then
            Call gsend(comm, ready, jdnode, MsdWrite_tag)

            Call grecv(comm, jatms, jdnode, MsdWrite_tag)
            If (jatms > 0) Then
              Call grecv(comm, chbuf(1:jatms), jdnode, MsdWrite_tag)
              Call grecv(comm, iwrk(1:jatms), jdnode, MsdWrite_tag)

              Call grecv(comm, ddd(1:jatms), jdnode, MsdWrite_tag)
              Call grecv(comm, eee(1:jatms), jdnode, MsdWrite_tag)
            End If
          End If

          jj = 0
          Do i = 1, jatms
            Write (record, Fmt='(a8,i10,1p,2e13.4,a8,a1)') chbuf(i), iwrk(i), ddd(i), eee(i), Repeat(' ', 8), lf
            jj = jj + 1
            Do k = 1, recsz
              chbat(k, jj) = record(k:k)
            End Do

            ! Dump batch and update start of file

            If (jj >= batsz .or. i == jatms) Then
              Write (Unit=files(FILE_MSD)%unit_no, Fmt='(53a)', Rec=msd_data%rec + Int(1, li)) (chbat(:, k), k=1, jj)
              msd_data%rec = msd_data%rec + Int(jj, li)
              jj = 0
            End If
          End Do
        End Do

        ! Update main header

        Write (Unit=files(FILE_MSD)%unit_no, Fmt='(i10,2i21,a1)', Rec=Int(2, li)) megatm, msd_data%frm, msd_data%rec, lf

        Call files(FILE_MSD)%close ()

      Else

        Do i = 1, config%natms
          k = 2 * i

          ddd(i) = Sqrt(stpval(36 + k - 1))

          If (Abs(dof_site(config%lsite(i))) > zero_plus) tmp = config%weight(i) * stpval(36 + k) / (boltz * 3.0_wp)
          eee(i) = tmp
        End Do

        Call grecv(comm, ready, 0, MsdWrite_tag)

        Call gsend(comm, config%natms, 0, MsdWrite_tag)
        If (config%natms > 0) Then
          Call gsend(comm, config%atmnam(1:config%natms), 0, MsdWrite_tag)
          Call gsend(comm, config%ltg(1:config%natms), 0, MsdWrite_tag)

          Call gsend(comm, ddd(1:config%natms), 0, MsdWrite_tag)
          Call gsend(comm, eee(1:config%natms), 0, MsdWrite_tag)
        End If

        ! Save offset pointer

        msd_data%rec = msd_data%rec + Int(megatm + 1, li)

      End If

      Deallocate (chbuf, iwrk, ddd, eee, Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'msd_write deallocation failure'
        Call error(0, message)
      End If

      ! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

    Else If (io_write == IO_WRITE_SORTED_MPIIO .or. &
             io_write == IO_WRITE_SORTED_DIRECT) Then

      rec_mpi_io = Int(msd_data%rec, offset_kind)
      jj = 0
      If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm=comm_self)
        Call io_init(io, recsz)
        Call io_open(io, io_write, comm_self, msd_data%fname, mode_wronly, fh)

        ! Write header information

        Write (record, Fmt='(a8,2i10,2f12.6,a1)') 'timestep', nstep, megatm, tstep, time, lf
        Call io_write_record(io, fh, rec_mpi_io, record)

        Call io_close(io, fh)
        Call io_finalize(io)

      End If
      Call gsync(comm)

      ! Start of file

      rec_mpi_io = Int(msd_data%rec, offset_kind) + Int(1, li)

      ! Write the rest

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_init(io, recsz)
      Call io_open(io, io_write, comm%comm, msd_data%fname, mode_wronly, fh)

      Allocate (ddd(1:config%mxatms), eee(1:config%mxatms), Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'msd_write allocation failure'
        Call error(0, message)
      End If

      Do i = 1, config%natms
        k = 2 * i

        ddd(i) = Sqrt(stpval(36 + k - 1))

        If (Abs(dof_site(config%lsite(i))) > zero_plus) tmp = config%weight(i) * stpval(36 + k) / (boltz * 3.0_wp)
        eee(i) = tmp
      End Do

      Call io_write_sorted_file(io, fh, 0, IO_MSDTMP, rec_mpi_io, config%natms, &
                                config%ltg, config%atmnam, ddd, eee, (/0.0_wp/), &
                                (/0.0_wp/), (/0.0_wp/), (/0.0_wp/), &
                                (/0.0_wp/), (/0.0_wp/), (/0.0_wp/), &
                                (/0.0_wp/), (/0.0_wp/), (/0.0_wp/), ierr)

      If (ierr /= 0) Then
        Select Case (ierr)
        Case (IO_BASE_COMM_NOT_SET)
          Call error(1050)
        Case (IO_ALLOCATION_ERROR)
          Call error(1053)
        Case (IO_UNKNOWN_WRITE_OPTION)
          Call error(1056)
        Case (IO_UNKNOWN_WRITE_LEVEL)
          Call error(1059)
        End Select
      End If

      ! Update and save offset pointer

      msd_data%rec = msd_data%rec + Int(1, li) + Int(megatm, li)
      If (comm%idnode == 0) Then
        Write (record, Fmt='(i10,2i21,a1)') megatm, msd_data%frm, msd_data%rec, lf
        Call io_write_record(io, fh, Int(1, offset_kind), record)
      End If

      Call io_close(io, fh)
      Call io_finalize(io)

      Deallocate (ddd, eee, Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'msd_write deallocation failure'
        Call error(0, message)
      End If

      ! SORTED Serial Direct Access FORTRAN

    Else If (io_write == IO_WRITE_SORTED_MASTER) Then

      Allocate (chbuf(1:config%mxatms), iwrk(1:config%mxatms), ddd(1:config%mxatms), eee(1:config%mxatms), Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'msd_write allocation failure'
        Call error(0, message)
      End If

      ! node 0 handles I/O

      If (comm%idnode == 0) Then

        Open (Newunit=files(FILE_MSD)%unit_no, File='fname', Form='formatted', Access='direct', Recl=recsz)

        msd_data%rec = msd_data%rec + Int(1, li)
        Write (Unit=files(FILE_MSD)%unit_no, Fmt='(a8,2i10,2f12.6,a1)', Rec=msd_data%rec) 'timestep', nstep, &
          megatm, tstep, time, lf

        Do i = 1, config%natms
          k = 2 * i

          iwrk(i) = config%ltg(i)
          chbuf(i) = config%atmnam(i)

          ddd(i) = Sqrt(stpval(36 + k - 1))

          If (Abs(dof_site(config%lsite(i))) > zero_plus) tmp = config%weight(i) * stpval(36 + k) / (boltz * 3.0_wp)
          eee(i) = tmp
        End Do

        jatms = config%natms
        ready = .true.
        Do jdnode = 0, comm%mxnode - 1
          If (jdnode > 0) Then
            Call gsend(comm, ready, jdnode, MsdWrite_tag)

            Call grecv(comm, jatms, jdnode, MsdWrite_tag)
            If (jatms > 0) Then
              Call grecv(comm, chbuf(1:jatms), jdnode, MsdWrite_tag)
              Call grecv(comm, iwrk(1:jatms), jdnode, MsdWrite_tag)

              Call grecv(comm, ddd(1:jatms), jdnode, MsdWrite_tag)
              Call grecv(comm, eee(1:jatms), jdnode, MsdWrite_tag)
            End If
          End If

          Do i = 1, jatms
            Write (Unit=files(FILE_MSD)%unit_no, Fmt='(a8,i10,1p,2e13.4,a8,a1)', Rec=msd_data%rec + Int(iwrk(i), li)) &
              chbuf(i), iwrk(i), ddd(i), eee(i), Repeat(' ', 8), lf
          End Do
        End Do

        ! Update and save offset pointer

        msd_data%rec = msd_data%rec + Int(megatm, li)
        Write (Unit=files(FILE_MSD)%unit_no, Fmt='(i10,2i21,a1)', Rec=Int(2, li)) megatm, msd_data%frm, msd_data%rec, lf

        Call files(FILE_MSD)%close ()

      Else

        Do i = 1, config%natms
          k = 2 * i

          ddd(i) = Sqrt(stpval(36 + k - 1))

          If (Abs(dof_site(config%lsite(i))) > zero_plus) tmp = config%weight(i) * stpval(36 + k) / (boltz * 3.0_wp)
          eee(i) = tmp
        End Do

        Call grecv(comm, ready, 0, MsdWrite_tag)

        Call gsend(comm, config%natms, 0, MsdWrite_tag)
        If (config%natms > 0) Then
          Call gsend(comm, config%atmnam(1:config%natms), 0, MsdWrite_tag)
          Call gsend(comm, config%ltg(1:config%natms), 0, MsdWrite_tag)

          Call gsend(comm, ddd(1:config%natms), 0, MsdWrite_tag)
          Call gsend(comm, eee(1:config%natms), 0, MsdWrite_tag)
        End If

        msd_data%rec = msd_data%rec + Int(megatm + 1, li)

      End If

      Deallocate (chbuf, iwrk, ddd, eee, Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'msd_write deallocation failure'
        Call error(0, message)
      End If

    End If

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT .or. &
        io_write == IO_WRITE_UNSORTED_MASTER) Then
      Deallocate (n_atm, Stat=fail(1))
      Deallocate (chbat, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'msd_write deallocation failure 0'
        Call error(0, message)
      End If
    End If

    Call gsync(comm)

  End Subroutine msd_write
End Module msd
