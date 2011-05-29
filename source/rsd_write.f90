Subroutine rsd_write(imcon,keyres,nsrsd,isrsd,rrsd,nstep,tstep,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing RSDDAT file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use site_module,       Only : ntpshl,unqshl
  Use config_module,     Only : cfgname,cell,natms, &
                                atmnam,ltg,xxx,yyy,zzz
  Use statistics_module, Only : rsd

  Use parse_module,      Only : tabs_2_blanks, get_word, word_2_real
  Use io_module,         Only : io_set_parameters,        &
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

  Implicit None

  Integer,           Intent( In    ) :: imcon,keyres, &
                                        nsrsd,isrsd,nstep
  Real( Kind = wp ), Intent( In    ) :: rrsd,tstep,time

  Integer, Parameter :: recsz = 73 ! default record size

  Logical,           Save :: newjob = .true.
  Integer(Kind=ip),  Save :: rec    = 0_ip , &
                             frm    = 0_ip

  Logical                 :: safe,lexist,l_tmp,ready
  Character( Len = 40 )   :: word
  Integer                 :: fail(1:2),i,j,k,n,megn,jdnode,jatms
  Real( Kind = wp )       :: buffer(1:2)

! Some parameters and variables needed by io_module interfaces

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
        If (idnode == 0) Inquire(File='RSDDAT', Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist)
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (idnode == 0) Then
           Open(Unit=nrsddt, File='RSDDAT', Form='formatted', Access='direct', Status='replace', Recl=recsz)
           Write(Unit=nrsddt, Fmt='(a72,a1)',           Rec=1) cfgname(1:72),lf
           Write(Unit=nrsddt, Fmt='(f7.3,a23,2i21,a1)', Rec=2) rrsd,Repeat(' ',23),frm,rec,lf
           Close(Unit=nrsddt)
        End If
        rec=Int(2,ip)
        frm=Int(0,ip)

! Get some sense of it

     Else

        safe=.true.
        If (idnode == 0) Then

           Open(Unit=nrsddt, File='RSDDAT', Form='formatted')

           Do While (.true.)

              record=' '
              If (l_tmp) Then

                 Read(Unit=nrsddt, Fmt=*, End=20)            ! title record
                 rec=rec+Int(1,ip)
                 Read(Unit=nrsddt, Fmt='(a)', End=20) record ! bookkeeping record
                 rec=rec+Int(1,ip)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 If (word(1:Len_Trim(word)) /= 'timestep') Then
                    Call get_word(record,word) ; Call get_word(record,word)
                    Call get_word(record,word) ; frm=Nint(word_2_real(word,0.0_wp),ip)
                    Call get_word(record,word) ; rec=Nint(word_2_real(word,0.0_wp),ip)
                    If (frm /= Int(0,ip) .and. rec > Int(2,ip)) Then
                       Go To 20 ! New style
                    Else
                       l_tmp=.false. ! TOUGH, old style
                       rec=Int(2,ip)
                       frm=Int(0,ip)
                    End If
                 Else
                    safe=.false. ! Overwrite the file, it's junk to me
                    Go To 20
                 End If

              Else

                 Read(Unit=nrsddt, Fmt=*, End=20)            ! timestep record
                 rec=rec+Int(1,ip)

                 Read(Unit=nrsddt, Fmt='(a)', End=20) record ! displacments record
                 rec=rec+Int(1,ip)

                 Call tabs_2_blanks(record) ; Call get_word(record,word)
                 Call get_word(record,word) ; j=Nint(word_2_real(word))

                 Do i=1,3+2*j ! 3 lines for cell parameters and 2*j entries for displacements
                    Read(Unit=nrsddt, Fmt=*, End=20)
                    rec=rec+Int(1,ip)
                 End Do
                 frm=frm+Int(1,ip)

              End If

           End Do

20         Continue
           Close(Unit=nrsddt)

        End If

        If (mxnode > 1) Call gcheck(safe)
        If (.not.safe) Then
           lexist=.false.

           rec=Int(0,ip)
           frm=Int(0,ip)

           Go To 10
        Else If (mxnode > 1) Then
           buffer(1)=Real(frm,wp)
           buffer(2)=Real(rec,wp)

           Call MPI_BCAST(buffer(1:2), 2, wp_mpi, 0, dlp_comm_world, ierr)

           frm=Nint(buffer(1),ip)
           rec=Nint(buffer(2),ip)
        End If

     End If
  End If

  fail=0
  Allocate (nam(1:mxatms),ind(1:mxatms),dr(1:mxatms),  Stat=fail(1))
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'rsd_write allocation failure, node: ', idnode
     Call error(0)
  End If

  n=0
  Do i=1,natms
     If (rsd(i) > rrsd .and. (.not.Any(unqshl(1:ntpshl) == atmnam(i)))) Then
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
  If (mxnode > 1) Call gsum(megn)

! Get relative offsets for parallel printing

  Allocate (n_n(0:mxnode),          Stat=fail(1))
  Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'rsd_write allocation failure 2, node: ', idnode
     Call error(0)
  End If

  n_n=0 ; n_n(idnode+1)=n
  If (mxnode > 1) Call gsum(n_n)
  n_n(0)=Sum(n_n(0:idnode))

  chbat=' '

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(kind = MPI_OFFSET_KIND)

! Update frame

  frm=frm+Int(1,ip)

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     j=0
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_open( io_write, MPI_COMM_SELF, 'RSDDAT', MPI_MODE_WRONLY, fh )

        Write(record, Fmt='(a8,i10,2f12.6,i5,f7.3,a18,a1)') &
           'timestep',nstep,tstep,time,imcon,rrsd,Repeat(' ',18),lf
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
     Call gsync()

! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(j,MPI_OFFSET_KIND)+Int(2,MPI_OFFSET_KIND)*Int(n_n(0),MPI_OFFSET_KIND)
     j=0

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, 'RSDDAT', MPI_MODE_WRONLY, fh )

     Do i=1,n
        Write(record, Fmt='(a8,i10,f7.3,a47,a1)') nam(i),ind(i),dr(i),Repeat(' ',47),lf
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

     rec=rec+Int(5,ip)+Int(2,ip)*Int(megn,ip)
     If (idnode == 0) Then
        Write(record, Fmt='(f7.3,a23,2i21,a1)') rrsd,Repeat(' ',23),frm,rec,lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record )
     End If

     Call io_close( fh )
     Call io_finalize

  Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
           io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(2))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'rsd_write allocation failure 3, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O
! Start of file

     j=0
     If (idnode == 0) Then
        Open(Unit=nrsddt, File='RSDDAT', Form='formatted', Access='direct', Recl=recsz)

! Accumulate header

        Write(record, Fmt='(a8,i10,2f12.6,i5,f7.3,a18,a1)') &
           'timestep',nstep,tstep,time,imcon,rrsd,Repeat(' ',18),lf
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

        Write(Unit=nrsddt, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,j)
        rec=rec+Int(j,ip)
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
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,RsdWrite_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(dr,jatms,wp_mpi,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(byy,jatms,wp_mpi,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,RsdWrite_tag,dlp_comm_world,status,ierr)
              End If
           End If

           Do i=1,jatms
              Write(record, Fmt='(a8,i10,f7.3,a47,a1)') chbuf(i),iwrk(i),dr(i),Repeat(' ',47),lf
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
                 Write(Unit=nrsddt, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,j)
                 rec=rec+Int(j,ip)
                 j=0
              End If
           End Do
        End Do

! Update main header

        Write(Unit=nrsddt, Fmt='(f7.3,a23,2i21,a1)', Rec=2) rrsd,Repeat(' ',23),frm,rec,lf

        Close(Unit=nrsddt)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,RsdWrite_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(n,1,MPI_INTEGER,0,RsdWrite_tag,dlp_comm_world,ierr)
        If (n > 0) Then
           Call MPI_SEND(nam,8*n,MPI_CHARACTER,0,RsdWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ind,n,MPI_INTEGER,0,RsdWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(dr,n,wp_mpi,0,RsdWrite_tag,dlp_comm_world,ierr)

           Call MPI_SEND(axx,n,wp_mpi,0,RsdWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ayy,n,wp_mpi,0,RsdWrite_tag,dlp_comm_world,ierr)
           Call MPI_SEND(azz,n,wp_mpi,0,RsdWrite_tag,dlp_comm_world,ierr)
        End If

! Save offset pointer

        rec=rec+Int(5,ip)+Int(2,ip)*Int(megn,ip)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (bxx,byy,bzz, Stat=fail(2))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'rsd_write deallocation failure 3, node: ', idnode
        Call error(0)
     End If

  End If

  If (mxnode > 1) Call gsync()

  Deallocate (n_n,   Stat=fail(1))
  Deallocate (chbat, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'rsd_write deallocation failure 2, node: ', idnode
     Call error(0)
  End If

  Deallocate (nam,ind,dr,  Stat=fail(1))
  Deallocate (axx,ayy,azz, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'rsd_write deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine rsd_write
