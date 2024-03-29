Module rsds
  Use comms,           Only: RsdWrite_tag,&
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
  Use constants,       Only: nrsddt
  Use core_shell,      Only: core_shell_type
  ! this is assymetric with respect to the rest. will need probably rsd defined in this type
  Use errors_warnings, Only: error
  Use flow_control,    Only: RESTART_KEY_OLD
  Use io,              Only: &
                             IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
                             IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT, &
                             IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, io_close, &
                             io_finalize, io_get_parameters, io_init, io_open, io_set_parameters, &
                             io_type, io_write_batch, io_write_record
  Use kinds,           Only: li,STR_LEN,&
                             wi,&
                             wp
  Use parse,           Only: get_word,&
                             tabs_2_blanks,&
                             word_2_real

  Implicit None

  Private

  Type, Public :: rsd_type
    Private

    !> Active
    Logical, Public          :: lrsd = .false.
    !> Initialised
    Logical                  :: newjob = .true.
    Integer(Kind=li)         :: rec = 0_li, &
                                frm = 0_li
    !> Step to start at
    Integer(Kind=wi), Public :: nsrsd = 0
    !> Frequency of collection
    Integer(Kind=wi), Public :: isrsd = 1
    !> Cutoff for accepting an atom has moved
    Real(Kind=wp), Public    :: rrsd = 0.15_wp

  End Type
  Public :: rsd_write

Contains

  Subroutine rsd_write(keyres, nstep, tstep, io, rsdc, time, cshell, rsd, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for writing RSDDAT file at selected intervals
    ! in simulation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: keyres, nstep
    Real(Kind=wp),            Intent(In   ) :: tstep
    Type(io_type),            Intent(InOut) :: io
    Type(rsd_type),           Intent(InOut) :: rsdc
    Real(Kind=wp),            Intent(In   ) :: time
    Type(core_shell_type),    Intent(InOut) :: cshell
    Real(Kind=wp),            Intent(InOut) :: rsd(:)
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Integer, Parameter :: recsz = 73

    Character                                      :: lf
    Character(Len=1), Allocatable, Dimension(:, :) :: chbat
    Character(Len=STR_LEN)                             :: message
    Character(Len=40)                              :: word
    Character(Len=8), Allocatable, Dimension(:)    :: chbuf, nam
    Character(Len=recsz)                           :: record
    Integer                                        :: batsz, fail(1:2), fh, i, io_write, j, jatms, &
                                                      jdnode, k, megn, n
    Integer(Kind=offset_kind)                      :: rec_mpi_io
    Integer, Allocatable, Dimension(:)             :: ind, iwrk, n_n
    Logical                                        :: l_tmp, lexist, ready, safe
    Real(Kind=wp)                                  :: buffer(1:2)
    Real(Kind=wp), Allocatable, Dimension(:)       :: axx, ayy, azz, bxx, byy, bzz, dr

! default record size
! Some parameters and variables needed by io interfaces

    If (.not. (nstep >= rsdc%nsrsd .and. Mod(nstep - rsdc%nsrsd, rsdc%isrsd) == 0)) Return

    ! Get write buffer size and line feed character

    Call io_get_parameters(io, user_method_write=io_write)
    Call io_get_parameters(io, user_buffer_size_write=batsz)
    Call io_get_parameters(io, user_line_feed=lf)

    If (rsdc%newjob) Then
      rsdc%newjob = .false.

      ! If the keyres = RESTART_KEY_OLD, is RSDDAT old (does it exist) and
      ! how many frames and records are in there

      lexist = .true.
      If (keyres == RESTART_KEY_OLD) Then
        If (comm%idnode == 0) Inquire (File='RSDDAT', Exist=lexist)
        Call gcheck(comm, lexist, "enforce")
      Else
        lexist = .false.
      End If

      ! Generate file is non-existent

      10 Continue
      If (.not. lexist) Then

        If (comm%idnode == 0) Then
          Open (Unit=nrsddt, File='RSDDAT', Form='formatted', Access='direct', Status='replace', Recl=recsz)
          Write (Unit=nrsddt, Fmt='(a72,a1)', Rec=1) config%cfgname(1:72), lf
          Write (Unit=nrsddt, Fmt='(f11.3,a19,2i21,a1)', Rec=2) rsdc%rrsd, Repeat(' ', 19), rsdc%frm, rsdc%rec, lf
          Close (Unit=nrsddt)
        End If
        rsdc%rec = Int(2, li)
        rsdc%frm = Int(0, li)

        ! Get some sense of it

      Else

        safe = .true.
        If (comm%idnode == 0) Then

          Open (Unit=nrsddt, File='RSDDAT', Form='formatted')

          l_tmp = .true.
          Do

            record = ' '
            If (l_tmp) Then

              Read (Unit=nrsddt, Fmt=*, End=20) ! title record
              rsdc%rec = rsdc%rec + Int(1, li)
              Read (Unit=nrsddt, Fmt='(a)', End=20) record ! bookkeeping record
              rsdc%rec = rsdc%rec + Int(1, li)

              Call tabs_2_blanks(record); Call get_word(record, word)
              If (word(1:Len_trim(word)) /= 'timestep') Then
                Call get_word(record, word); Call get_word(record, word)
                Call get_word(record, word); rsdc%frm = Nint(word_2_real(word, 0.0_wp), li)
                Call get_word(record, word); rsdc%rec = Nint(word_2_real(word, 0.0_wp), li)
                If (rsdc%frm /= Int(0, li) .and. rsdc%rec > Int(2, li)) Then
                  Go To 20 ! New style
                Else
                  l_tmp = .false. ! TOUGH, old style
                  rsdc%rec = Int(2, li)
                  rsdc%frm = Int(0, li)
                End If
              Else
                safe = .false. ! Overwrite the file, it's junk to me
                Go To 20
              End If

            Else

              Read (Unit=nrsddt, Fmt=*, End=20) ! timestep record
              rsdc%rec = rsdc%rec + Int(1, li)

              Read (Unit=nrsddt, Fmt='(a)', End=20) record ! displacments record
              rsdc%rec = rsdc%rec + Int(1, li)

              Call tabs_2_blanks(record); Call get_word(record, word)
              Call get_word(record, word); j = Nint(word_2_real(word))

              Do i = 1, 3 + 2 * j ! 3 lines for config%cell parameters and 2*j entries for displacements
                Read (Unit=nrsddt, Fmt=*, End=20)
                rsdc%rec = rsdc%rec + Int(1, li)
              End Do
              rsdc%frm = rsdc%frm + Int(1, li)

            End If

          End Do

          20 Continue
          Close (Unit=nrsddt)

        End If

        Call gcheck(comm, safe, "enforce")
        If (.not. safe) Then
          lexist = .false.

          rsdc%rec = Int(0, li)
          rsdc%frm = Int(0, li)

          Go To 10
        Else If (comm%mxnode > 1) Then
          buffer(1) = Real(rsdc%frm, wp)
          buffer(2) = Real(rsdc%rec, wp)

          Call gbcast(comm, buffer, 0)

          rsdc%frm = Nint(buffer(1), li)
          rsdc%rec = Nint(buffer(2), li)
        End If

      End If
    End If

    fail = 0
    Allocate (nam(1:config%mxatms), ind(1:config%mxatms), dr(1:config%mxatms), Stat=fail(1))
    Allocate (axx(1:config%mxatms), ayy(1:config%mxatms), azz(1:config%mxatms), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'rsd_write allocation failure'
      Call error(0, message)
    End If

    n = 0
    Do i = 1, config%natms
      If (rsd(i) > rsdc%rrsd .and. cshell%legshl(0, i) >= 0) Then
        n = n + 1
        nam(n) = config%atmnam(i)
        ind(n) = config%ltg(i)
        dr(n) = rsd(i)

        axx(n) = config%parts(i)%xxx
        ayy(n) = config%parts(i)%yyy
        azz(n) = config%parts(i)%zzz
      End If
    End Do

    ! Sum global displacements values

    megn = n
    Call gsum(comm, megn)

    ! Get relative offsets for parallel printing

    Allocate (n_n(0:comm%mxnode), Stat=fail(1))
    Allocate (chbat(1:recsz, 1:batsz), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'rsd_write allocation failure 2'
      Call error(0, message)
    End If

    n_n = 0; n_n(comm%idnode + 1) = n
    Call gsum(comm, n_n)
    n_n(0) = Sum(n_n(0:comm%idnode))

    chbat = ' '

    ! Notes:
    ! the MPI-I/O records are numbersdc%rsdc%rec from 0 (not 1)
    ! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
    !   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
    !   both declarsdc%rsdc%rec as: Integer(Kind = offset_kind)

    ! Update frame

    rsdc%frm = rsdc%frm + Int(1, li)

    If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
        io_write == IO_WRITE_UNSORTED_DIRECT .or. &
        io_write == IO_WRITE_SORTED_MPIIO .or. &
        io_write == IO_WRITE_SORTED_DIRECT) Then

      ! Write header and cell information, where just one node is needed
      ! Start of file

      rec_mpi_io = Int(rsdc%rec, offset_kind)
      j = 0
      If (comm%idnode == 0) Then

        Call io_set_parameters(io, user_comm=comm_self)
        Call io_init(io, recsz)
        Call io_open(io, io_write, comm_self, 'RSDDAT', mode_wronly, fh)

        Write (record, Fmt='(a8,i10,2f20.6,i3,f11.3,a1)') &
          'timestep', nstep, tstep, time, config%imcon, rsdc%rrsd, lf
        j = j + 1
        Do k = 1, recsz
          chbat(k, j) = record(k:k)
        End Do

        Write (record, Fmt='(a14,i10,a48,a1)') 'displacements ', megn, Repeat(' ', 48), lf
        j = j + 1
        Do k = 1, recsz
          chbat(k, j) = record(k:k)
        End Do

        Do i = 0, 2
          Write (record, '( 3f20.10, a12, a1 )') &
            config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3), Repeat(' ', 12), lf
          j = j + 1
          Do k = 1, recsz
            chbat(k, j) = record(k:k)
          End Do
        End Do

        ! Dump header and cell information

        Call io_write_batch(io, fh, rec_mpi_io, j, chbat)

        Call io_close(io, fh)
        Call io_finalize(io)

      Else

        j = j + 5

      End If
      Call gsync(comm)

      ! Start of file

      rec_mpi_io = Int(rsdc%rec, offset_kind) + Int(j, offset_kind) + Int(2, offset_kind) * Int(n_n(0), offset_kind)
      j = 0

      Call io_set_parameters(io, user_comm=comm%comm)
      Call io_init(io, recsz)
      Call io_open(io, io_write, comm%comm, 'RSDDAT', mode_wronly, fh)

      Do i = 1, n
        Write (record, Fmt='(a8,i10,f11.3,a43,a1)') nam(i), ind(i), dr(i), Repeat(' ', 43), lf
        j = j + 1
        Do k = 1, recsz
          chbat(k, j) = record(k:k)
        End Do

        Write (record, Fmt='(3g20.10,a12,a1)') axx(i), ayy(i), azz(i), Repeat(' ', 12), lf
        j = j + 1
        Do k = 1, recsz
          chbat(k, j) = record(k:k)
        End Do

        ! Dump batch and update start of file

        If (j + 2 >= batsz .or. i == n) Then
          Call io_write_batch(io, fh, rec_mpi_io, j, chbat)
          rec_mpi_io = rec_mpi_io + Int(j, offset_kind)
          j = 0
        End If
      End Do

      ! Update and save offset pointer

      rsdc%rec = rsdc%rec + Int(5, li) + Int(2, li) * Int(megn, li)
      If (comm%idnode == 0) Then
        Write (record, Fmt='(f11.3,a19,2i21,a1)') rsdc%rrsd, Repeat(' ', 19), rsdc%frm, rsdc%rec, lf
        Call io_write_record(io, fh, Int(1, offset_kind), record)
      End If

      Call io_close(io, fh)
      Call io_finalize(io)

    Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
             io_write == IO_WRITE_SORTED_MASTER) Then

      Allocate (chbuf(1:config%mxatms), iwrk(1:config%mxatms), Stat=fail(1))
      Allocate (bxx(1:config%mxatms), byy(1:config%mxatms), bzz(1:config%mxatms), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'rsd_write allocation failure 3'
        Call error(0, message)
      End If

      ! node 0 handles I/O
      ! Start of file

      j = 0
      If (comm%idnode == 0) Then
        Open (Unit=nrsddt, File='RSDDAT', Form='formatted', Access='direct', Recl=recsz)

        ! Accumulate header

        Write (record, Fmt='(a8,i10,2f20.6,i3,f11.3,a1)') &
          'timestep', nstep, tstep, time, config%imcon, rsdc%rrsd, lf
        j = j + 1
        Do k = 1, recsz
          chbat(k, j) = record(k:k)
        End Do

        Write (record, Fmt='(a14,i10,a48,a1)') 'displacements ', megn, Repeat(' ', 48), lf
        j = j + 1
        Do k = 1, recsz
          chbat(k, j) = record(k:k)
        End Do

        Do i = 0, 2
          Write (record, Fmt='(3f20.10,a12,a1)') &
            config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3), Repeat(' ', 12), lf
          j = j + 1
          Do k = 1, recsz
            chbat(k, j) = record(k:k)
          End Do
        End Do

        ! Dump header and update start of file

        Write (Unit=nrsddt, Fmt='(73a)', Rec=rsdc%rec + Int(1, li)) (chbat(:, k), k=1, j)
        rsdc%rec = rsdc%rec + Int(j, li)
        j = 0

        Do i = 1, n
          chbuf(i) = nam(i)
          iwrk(i) = ind(i)
          dr(i) = rsd(i)

          bxx(i) = axx(i)
          byy(i) = ayy(i)
          bzz(i) = azz(i)
        End Do

        jatms = n
        ready = .true.
        Do jdnode = 0, comm%mxnode - 1
          If (jdnode > 0) Then
            Call gsend(comm, ready, jdnode, RsdWrite_tag)

            Call grecv(comm, jatms, jdnode, RsdWrite_tag)
            If (jatms > 0) Then
              Call grecv(comm, chbuf(1:jatms), jdnode, RsdWrite_tag)
              Call grecv(comm, iwrk(1:jatms), jdnode, RsdWrite_tag)
              Call grecv(comm, dr(1:jatms), jdnode, RsdWrite_tag)

              Call grecv(comm, bxx(1:jatms), jdnode, RsdWrite_tag)
              Call grecv(comm, byy(1:jatms), jdnode, RsdWrite_tag)
              Call grecv(comm, bzz(1:jatms), jdnode, RsdWrite_tag)
            End If
          End If

          Do i = 1, jatms
            Write (record, Fmt='(a8,i10,f11.3,a43,a1)') chbuf(i), iwrk(i), dr(i), Repeat(' ', 43), lf
            j = j + 1
            Do k = 1, recsz
              chbat(k, j) = record(k:k)
            End Do

            Write (record, Fmt='(3g20.10,a12,a1)') bxx(i), byy(i), bzz(i), Repeat(' ', 12), lf
            j = j + 1
            Do k = 1, recsz
              chbat(k, j) = record(k:k)
            End Do

            ! Dump batch and update start of file

            If (j + 2 >= batsz .or. i == jatms) Then
              Write (Unit=nrsddt, Fmt='(73a)', Rec=rsdc%rec + Int(1, li)) (chbat(:, k), k=1, j)
              rsdc%rec = rsdc%rec + Int(j, li)
              j = 0
            End If
          End Do
        End Do

        ! Update main header

        Write (Unit=nrsddt, Fmt='(f11.3,a19,2i21,a1)', Rec=2) rsdc%rrsd, Repeat(' ', 19), rsdc%frm, rsdc%rec, lf

        Close (Unit=nrsddt)

      Else

        Call grecv(comm, ready, 0, RsdWrite_tag)

        Call gsend(comm, n, 0, RsdWrite_tag)
        If (n > 0) Then
          Call gsend(comm, nam(1:n), 0, RsdWrite_tag)
          Call gsend(comm, ind(1:n), 0, RsdWrite_tag)
          Call gsend(comm, dr(1:n), 0, RsdWrite_tag)

          Call gsend(comm, axx(1:n), 0, RsdWrite_tag)
          Call gsend(comm, ayy(1:n), 0, RsdWrite_tag)
          Call gsend(comm, azz(1:n), 0, RsdWrite_tag)
        End If

        ! Save offset pointer

        rsdc%rec = rsdc%rec + Int(5, li) + Int(2, li) * Int(megn, li)

      End If

      Deallocate (chbuf, iwrk, Stat=fail(1))
      Deallocate (bxx, byy, bzz, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'rsd_write deallocation failure 3'
        Call error(0, message)
      End If

    End If

    Call gsync(comm)

    Deallocate (n_n, Stat=fail(1))
    Deallocate (chbat, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a,i0)') 'rsd_write deallocation failure 2'
      Call error(0, message)
    End If

    Deallocate (nam, ind, dr, Stat=fail(1))
    Deallocate (axx, ayy, azz, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'rsd_write deallocation failure'
      Call error(0, message)
    End If

  End Subroutine rsd_write
End Module rsds
