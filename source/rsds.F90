Module rsds
  Use kinds, Only : wp, li
  Use comms, Only : comms_type,gbcast,RsdWrite_tag,gsum,wp_mpi,gsync,gcheck, &
                    gsend,grecv
  Use setup
  Use configuration,     Only : cfgname,imcon,cell,natms, &
                                atmnam,ltg,xxx,yyy,zzz
  Use core_shell, Only : legshl
  Use statistics, Only : rsd

  Use parse,      Only : tabs_2_blanks, get_word, word_2_real
  Use io,         Only : io_set_parameters,        &
                                io_get_parameters,        &
                                io_init, io_open,         &
                                io_write_record,          &
                                io_write_batch,           &
                                io_close, io_finalize,    &
                                IO_WRITE_UNSORTED_MPIIO,  &
                                IO_WRITE_UNSORTED_DIRECT, &
                                IO_WRITE_UNSORTED_MASTER, &
                                IO_WRITE_SORTED_MPIIO,    &
                                IO_WRITE_SORTED_DIRECT,   &
                                IO_WRITE_SORTED_NETCDF,   &
                                IO_WRITE_SORTED_MASTER

  Use errors_warnings, Only : error
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif
  Implicit None

Contains


Subroutine rsd_write(keyres,nsrsd,isrsd,rrsd,nstep,tstep,time,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing RSDDAT file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: keyres,nsrsd,isrsd,nstep
  Real( Kind = wp ), Intent( In    ) :: rrsd,tstep,time
  Type( comms_type ), Intent( InOut ) :: comm

  Integer, Parameter :: recsz = 73 ! default record size

  Logical,           Save :: newjob = .true.
  Integer(Kind=li),  Save :: rec    = 0_li , &
                             frm    = 0_li

  Logical                 :: safe,lexist,l_tmp,ready
  Character( Len = 40 )   :: word
  Integer                 :: fail(1:2),i,j,k,n,megn,jdnode,jatms
  Real( Kind = wp )       :: buffer(1:2)

! Some parameters and variables needed by io interfaces

  Integer                           :: fh, io_write, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

  Character( Len = 8 ), Dimension( : ),    Allocatable :: nam
  Integer,              Dimension( : ),    Allocatable :: ind
  Real( Kind = wp )   , Dimension( : ),    Allocatable :: dr

  Integer,              Dimension( : ),    Allocatable :: n_n

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz

  Character( Len = 256 ) :: message
  If (.not.(nstep >= nsrsd .and. Mod(nstep-nsrsd,isrsd) == 0)) Return

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! netCDF not implemented for RSDDAT.  Switch to DEFAULT temporarily.

  If (io_write == IO_WRITE_SORTED_NETCDF) io_write = IO_WRITE_SORTED_MPIIO

  If (newjob) Then
     newjob = .false.

! If the keyres=1, is RSDDAT old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (comm%idnode == 0) Inquire(File='RSDDAT', Exist=lexist)
        Call gcheck(comm,lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (comm%idnode == 0) Then
           Open(Unit=nrsddt, File='RSDDAT', Form='formatted', Access='direct', Status='replace', Recl=recsz)
           Write(Unit=nrsddt, Fmt='(a72,a1)',           Rec=1) cfgname(1:72),lf
           Write(Unit=nrsddt, Fmt='(f11.3,a19,2i21,a1)', Rec=2) rrsd,Repeat(' ',19),frm,rec,lf
           Close(Unit=nrsddt)
        End If
        rec=Int(2,li)
        frm=Int(0,li)

! Get some sense of it

     Else

        safe=.true.
        If (comm%idnode == 0) Then

           Open(Unit=nrsddt, File='RSDDAT', Form='formatted')

           l_tmp =.true.
           Do

              record=' '
              If (l_tmp) Then

                 Read(Unit=nrsddt, Fmt=*, End=20)            ! title record
                 rec=rec+Int(1,li)
                 Read(Unit=nrsddt, Fmt='(a)', End=20) record ! bookkeeping record
                 rec=rec+Int(1,li)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 If (word(1:Len_Trim(word)) /= 'timestep') Then
                    Call get_word(record,word) ; Call get_word(record,word)
                    Call get_word(record,word) ; frm=Nint(word_2_real(word,0.0_wp),li)
                    Call get_word(record,word) ; rec=Nint(word_2_real(word,0.0_wp),li)
                    If (frm /= Int(0,li) .and. rec > Int(2,li)) Then
                       Go To 20 ! New style
                    Else
                       l_tmp=.false. ! TOUGH, old style
                       rec=Int(2,li)
                       frm=Int(0,li)
                    End If
                 Else
                    safe=.false. ! Overwrite the file, it's junk to me
                    Go To 20
                 End If

              Else

                 Read(Unit=nrsddt, Fmt=*, End=20)            ! timestep record
                 rec=rec+Int(1,li)

                 Read(Unit=nrsddt, Fmt='(a)', End=20) record ! displacments record
                 rec=rec+Int(1,li)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 Call get_word(record,word) ; j=Nint(word_2_real(word))

                 Do i=1,3+2*j ! 3 lines for cell parameters and 2*j entries for displacements
                    Read(Unit=nrsddt, Fmt=*, End=20)
                    rec=rec+Int(1,li)
                 End Do
                 frm=frm+Int(1,li)

              End If

           End Do

20         Continue
           Close(Unit=nrsddt)

        End If

        Call gcheck(comm,safe,"enforce")
        If (.not.safe) Then
           lexist=.false.

           rec=Int(0,li)
           frm=Int(0,li)

           Go To 10
        Else If (comm%mxnode > 1) Then
           buffer(1)=Real(frm,wp)
           buffer(2)=Real(rec,wp)

           Call gbcast(comm,buffer,0)

           frm=Nint(buffer(1),li)
           rec=Nint(buffer(2),li)
        End If

     End If
  End If

  fail=0
  Allocate (nam(1:mxatms),ind(1:mxatms),dr(1:mxatms),  Stat=fail(1))
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'rsd_write allocation failure'
     Call error(0,message)
  End If

  n=0
  Do i=1,natms
     If (rsd(i) > rrsd .and. legshl(0,i) >= 0) Then
        n=n+1
        nam(n)=atmnam(i)
        ind(n)=ltg(i)
        dr(n) =rsd(i)

        axx(n)=xxx(i)
        ayy(n)=yyy(i)
        azz(n)=zzz(i)
     End If
  End Do

! Sum global displacements values

  megn = n
  Call gsum(comm,megn)

! Get relative offsets for parallel printing

  Allocate (n_n(0:comm%mxnode),          Stat=fail(1))
  Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'rsd_write allocation failure 2'
     Call error(0,message)
  End If

  n_n=0 ; n_n(comm%idnode+1)=n
  Call gsum(comm,n_n)
  n_n(0)=Sum(n_n(0:comm%idnode))

  chbat=' '

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = MPI_OFFSET_KIND)

! Update frame

  frm=frm+Int(1,li)

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     j=0
     If (comm%idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_open( io_write, MPI_COMM_SELF, 'RSDDAT', MPI_MODE_WRONLY, fh )

        Write(record, Fmt='(a8,i10,2f20.6,i3,f11.3,a1)') &
           'timestep',nstep,tstep,time,imcon,rrsd,lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(a14,i10,a48,a1)') 'displacements ',megn,Repeat(' ',48),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Do i = 0, 2
           Write( record, '( 3f20.10, a12, a1 )' ) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           j=j+1
           Do k=1,recsz
              chbat(k,j) = record(k:k)
           End Do
        End Do

! Dump header and cell information

        Call io_write_batch( fh, rec_mpi_io, j, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        j=j+5

     End If
     Call gsync(comm)

! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(j,MPI_OFFSET_KIND)+Int(2,MPI_OFFSET_KIND)*Int(n_n(0),MPI_OFFSET_KIND)
     j=0

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, 'RSDDAT', MPI_MODE_WRONLY, fh )

     Do i=1,n
        Write(record, Fmt='(a8,i10,f11.3,a43,a1)') nam(i),ind(i),dr(i),Repeat(' ',43),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

! Dump batch and update start of file

        If (j + 2 >= batsz .or. i == n) Then
           Call io_write_batch( fh, rec_mpi_io, j, chbat )
           rec_mpi_io=rec_mpi_io+Int(j,MPI_OFFSET_KIND)
           j=0
        End If
     End Do

! Update and save offset pointer

     rec=rec+Int(5,li)+Int(2,li)*Int(megn,li)
     If (comm%idnode == 0) Then
        Write(record, Fmt='(f11.3,a19,2i21,a1)') rrsd,Repeat(' ',19),frm,rec,lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record )
     End If

     Call io_close( fh )
     Call io_finalize

  Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
           io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'rsd_write allocation failure 3'
        Call error(0,message)
     End If

! node 0 handles I/O
! Start of file

     j=0
     If (comm%idnode == 0) Then
        Open(Unit=nrsddt, File='RSDDAT', Form='formatted', Access='direct', Recl=recsz)

! Accumulate header

        Write(record, Fmt='(a8,i10,2f20.6,i3,f11.3,a1)') &
           'timestep',nstep,tstep,time,imcon,rrsd,lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Write(record, Fmt='(a14,i10,a48,a1)') 'displacements ',megn,Repeat(' ',48),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record, Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           j=j+1
           Do k=1,recsz
              chbat(k,j) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=nrsddt, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,j)
        rec=rec+Int(j,li)
        j=0

        Do i=1,n
           chbuf(i)=nam(i)
           iwrk(i)=ind(i)
           dr(i)=rsd(i)

           bxx(i)=axx(i)
           byy(i)=ayy(i)
           bzz(i)=azz(i)
        End Do

        jatms=n
        ready=.true.
        Do jdnode=0,comm%mxnode-1
           If (jdnode > 0) Then
              Call gsend(comm,ready,jdnode,RsdWrite_tag)

              Call grecv(comm,jatms,jdnode,RsdWrite_tag)
              If (jatms > 0) Then
                 Call grecv(comm,chbuf(1:jatms),jdnode,RsdWrite_tag)
                 Call grecv(comm,iwrk(1:jatms),jdnode,RsdWrite_tag)
                 Call grecv(comm,dr(1:jatms),jdnode,RsdWrite_tag)

                 Call grecv(comm,bxx(1:jatms),jdnode,RsdWrite_tag)
                 Call grecv(comm,byy(1:jatms),jdnode,RsdWrite_tag)
                 Call grecv(comm,bzz(1:jatms),jdnode,RsdWrite_tag)
              End If
           End If

           Do i=1,jatms
              Write(record, Fmt='(a8,i10,f11.3,a43,a1)') chbuf(i),iwrk(i),dr(i),Repeat(' ',43),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

! Dump batch and update start of file

              If (j + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nrsddt, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,j)
                 rec=rec+Int(j,li)
                 j=0
              End If
           End Do
        End Do

! Update main header

        Write(Unit=nrsddt, Fmt='(f11.3,a19,2i21,a1)', Rec=2) rrsd,Repeat(' ',19),frm,rec,lf

        Close(Unit=nrsddt)

     Else

        Call grecv(comm,ready,0,RsdWrite_tag)

        Call gsend(comm,n,0,RsdWrite_tag)
        If (n > 0) Then
           Call gsend(comm,nam(1:n),0,RsdWrite_tag)
           Call gsend(comm,ind(1:n),0,RsdWrite_tag)
           Call gsend(comm,dr(1:n),0,RsdWrite_tag)

           Call gsend(comm,axx(1:n),0,RsdWrite_tag)
           Call gsend(comm,ayy(1:n),0,RsdWrite_tag)
           Call gsend(comm,azz(1:n),0,RsdWrite_tag)
        End If

! Save offset pointer

        rec=rec+Int(5,li)+Int(2,li)*Int(megn,li)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (bxx,byy,bzz, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'rsd_write deallocation failure 3'
        Call error(0,message)
     End If

  End If

  Call gsync(comm)

  Deallocate (n_n,   Stat=fail(1))
  Deallocate (chbat, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a,i0)') 'rsd_write deallocation failure 2'
     Call error(0,message)
  End If

  Deallocate (nam,ind,dr,  Stat=fail(1))
  Deallocate (axx,ayy,azz, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'rsd_write deallocation failure'
     Call error(0,message)
  End If

End Subroutine rsd_write
End Module rsds
