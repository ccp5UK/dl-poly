Subroutine msd_write(keyres,nstmsd,istmsd,megatm,nstep,tstep,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing MSDTMP file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2010
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use site_module,       Only : dofsit
  Use config_module,     Only : cfgname,natms,atmnam,lsite,ltg, &
                                weight

  Use statistics_module, Only : ravval
  Use parse_module,      Only : tabs_2_blanks, get_word, word_2_real
  Use io_module,         Only : io_set_parameters,             &
                                io_get_parameters,             &
                                io_init, io_nc_create,         &
                                io_open, io_write_record,      &
                                io_write_batch,                &
                                io_nc_get_dim,                 &
                                io_nc_put_var,                 &
                                io_write_sorted_file,          &
                                io_close, io_finalize,         &
                                io_nc_get_real_precision,      &
                                io_nc_get_file_real_precision, &
                                IO_MSDTMP,                     &
                                IO_BASE_COMM_NOT_SET,          &
                                IO_ALLOCATION_ERROR,           &
                                IO_UNKNOWN_WRITE_OPTION,       &
                                IO_UNKNOWN_WRITE_LEVEL,        &
                                IO_WRITE_UNSORTED_MPIIO,       &
                                IO_WRITE_UNSORTED_DIRECT,      &
                                IO_WRITE_UNSORTED_MASTER,      &
                                IO_WRITE_SORTED_MPIIO,         &
                                IO_WRITE_SORTED_DIRECT,        &
                                IO_WRITE_SORTED_NETCDF,        &
                                IO_WRITE_SORTED_MASTER

  Implicit None

  Integer,           Intent( In    ) :: keyres,nstmsd,istmsd, &
                                        megatm,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep,time

  Integer, Parameter :: recsz = 53 ! default record size

  Logical,               Save :: newjob = .true. , &
                                 fast   = .true.
  Integer(Kind=ip),      Save :: rec = 0_ip , &
                                 frm = 0_ip

  Logical                :: lexist,safe,ready
  Character( Len = 40 )  :: word,fname
  Integer                :: fail(1:2),i,j,k,jdnode,jatms
  Real( Kind = wp )      :: buffer(1:2),tmp

! Some parameters and variables needed by io_module interfaces

  Integer                           :: fh, io_write, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: ddd,eee

! netCDF check

  Integer :: file_p, file_r
  Integer :: io_p, io_r


  If (.not.(nstep >= nstmsd .and. Mod(nstep-nstmsd,istmsd) == 0)) Return

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

  If (newjob) Then
     newjob = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = 'MSDTMP'
     Else
        fname = 'MSDTMP.nc'
     End If

! If keyres=1, is MSDTMP old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (idnode == 0) Inquire(File=fname, Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist)
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (idnode == 0) Then
              Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Status='replace', Recl=recsz)
              Write(Unit=nhist, Fmt='(a52,a1)',      Rec=Int(1,ip)) cfgname(1:52),lf
              Write(Unit=nhist, Fmt='(i10,2i21,a1)', Rec=Int(2,ip)) megatm,frm,rec,lf
              Close(Unit=nhist)
           End If
           rec=Int(2,ip)
           frm=Int(0,ip)
        Else
           Call io_set_parameters( user_comm = dlp_comm_world )
           Call io_nc_create( dlp_comm_world, fname, cfgname, megatm )
        End If

! Get some sense of it

     Else

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           safe=.true.
           If (idnode == 0) Then

              Open(Unit=nhist, File=fname, Form='formatted')

              Do While (.true.)

                 record(1:recsz)=' '

! Assume new style of MSDTMP with bookkeeping.

                 If (fast) Then

                    Read(Unit=nhist, Fmt=*, End=10)            ! title record
                    rec=rec+Int(1,ip)
                    record=' '
                    Read(Unit=nhist, Fmt='(a)', End=10) record ! bookkeeping record
                    Call tabs_2_blanks(record)
                    rec=rec+Int(1,ip)

                    Call get_word(record,word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record,word) ; Call get_word(record,word)
                       Call get_word(record,word) ; frm=Nint(word_2_real(word,0.0_wp),ip)
                       Call get_word(record,word) ; rec=Nint(word_2_real(word,0.0_wp),ip)
                       If (frm /= Int(0,ip) .and. rec > Int(2,ip)) Then
                          Go To 20 ! New style
                       Else
                          fast=.false. ! TOUGH, old style
                          rec=Int(2,ip)
                          frm=Int(0,ip)
                       End If
                    Else
                       safe=.false. ! Overwrite the file, it's junk to me
                       Go To 20
                    End If

! TOUGH, it needs scanning through

                 Else

                    Read(Unit=nhist, Fmt='(a)', End=20) record(1:recsz) ! timestep record
                    Call tabs_2_blanks(record)
                    rec=rec+Int(1,ip)

                    Call get_word(record,word)
                    Call get_word(record,word) ; j=Nint(word_2_real(word)) ! total number of lines to read

                    word=' '
                    Write(word,'( "(", i0, "( / ) )" )') j-1
                    Read(Unit=nhist, Fmt=word, End=20)
                    rec=rec+Int(i,ip)
                    frm=frm+Int(1,ip)

                 End If

              End Do

20            Continue
              Close(Unit=nhist)

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

              Call gsum(buffer(1:2))

              frm=Nint(buffer(1),ip)
              rec=Nint(buffer(2),ip)
           End If

        Else ! netCDF read

           Call io_set_parameters( user_comm = dlp_comm_world )
           Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_RDONLY, fh )

! Get the precision that the msdtmp file was written in
! and check it matches the requested precision

           Call io_nc_get_file_real_precision( fh, file_p, file_r, ierr )
           safe = (ierr == 0)
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Write(nrite,'(/,1x,a)') &
  "Can not determine precision in an exisiting MSDTMP.nc file in msd_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync
              Call error(0)
           End If

           Call io_nc_get_real_precision( io_p, io_r, ierr )
           safe = (ierr == 0)
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Write(nrite,'(/,1x,a)') &
  "Can not determine the desired writing precision in msd_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync
              Call error(0)
           End If

           safe = (io_p == file_p .and. io_r == file_r)
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Then
                 Write(nrite,'(/,1x,a)') &
  "Requested writing precision inconsistent with that in an existing MSDTMP.nc"

                 Write(nrite, Fmt='(1x,a)', Advance='No') "Precision requested:"
                 Select Case( Selected_real_kind( io_p, io_r ) )
                 Case( Kind( 1.0 ) )
                    Write(nrite,'(1x,a)') "Single"
                 Case( Kind( 1.0d0 ) )
                    Write(nrite,'(1x,a)') "Double"
                 End Select
                 Write(nrite, Fmt='(1x,a)', Advance='No') "Precision in file  :"
                 Select Case( Selected_real_kind( file_p, file_r ) )
                 Case( Kind( 1.0 ) )
                    Write(nrite,'(1x,a)') "Single"
                 Case( Kind( 1.0d0 ) )
                    Write(nrite,'(1x,a)') "Double"
                 End Select
              End If

              Call gsync
              Call error(0)
           End If

! Get the frame number to check
! For netcdf this is the "frame number" which is not a long integer!

           Call io_nc_get_dim( 'frame', fh, i )

           If (i > 0) Then
              frm=Int(i,ip)

              Call io_close( fh )
           Else ! Overwrite the file, it's junk to me
              lexist=.false.

              rec=Int(0,ip)
              frm=Int(0,ip)

              Go To 10
           End If
        End If

     End If
  End If

! Get offsets and define batch

  fail=0
  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Allocate (n_atm(0:mxnode),        Stat=fail(1))
     Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write allocation failure 0, node: ', idnode
        Call error(0)
     End If

     chbat=' '
     n_atm=0 ; n_atm(idnode+1)=natms
     If (mxnode > 1) Call gsum(n_atm)
     n_atm(0)=Sum(n_atm(0:idnode))
  End If

! Update frame

  frm=frm+Int(1,ip)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_WRONLY, fh )

! Write header and update start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     If (idnode == 0) Then
        Write(record, Fmt='(a8,2i10,2f12.6,a1)') 'timestep',nstep,megatm,tstep,time,lf
        Call io_write_record( fh, rec_mpi_io, record )
     End If
     rec=rec+Int(1,ip)
     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(n_atm(0),MPI_OFFSET_KIND)

     j=0
     Do i=1,natms
        k=2*i

        If (Abs(dofsit(lsite(i))) > zero_plus) tmp=weight(i)*ravval(27+k)/(boltz*3.0_wp)

        Write(record, Fmt='(a8,i10,1p,2e13.4,a8,a1)') atmnam(i),ltg(i),Sqrt(ravval(27+k-1)),tmp,Repeat(' ',8),lf
        j=j+1
        Do k=1,recsz
           chbat(k,j) = record(k:k)
        End Do

! Dump batch and update start of file

        If (j >= batsz .or. i == natms) Then
           Call io_write_batch( fh, rec_mpi_io, j, chbat )
           rec_mpi_io=rec_mpi_io+Int(j,MPI_OFFSET_KIND)
           j=0
        End If
     End Do

! Update and save offset pointer

     rec=rec+Int(megatm,ip)
     If (idnode == 0) Then
        Write(record, Fmt='(i10,2i21,a1)') megatm,frm,rec,lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record )
     End If

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),ddd(1:mxatms),eee(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O

     If (idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,ip)
        Write(Unit=nhist, Fmt='(a8,2i10,2f12.6,a1)', Rec=rec) 'timestep',nstep,megatm,tstep,time,lf

        Do i=1,natms
           k=2*i

           iwrk(i)=ltg(i)
           chbuf(i)=atmnam(i)

           ddd(i)=Sqrt(ravval(27+k-1))

           If (Abs(dofsit(lsite(i))) > zero_plus) tmp=weight(i)*ravval(27+k)/(boltz*3.0_wp)
           eee(i)=tmp
        End Do

        jatms=natms

        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(ddd,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(eee,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
           End If

           j=0
           Do i=1,jatms
              Write(record, Fmt='(a8,i10,1p,2e13.4,a8,a1)') chbuf(i),iwrk(i),ddd(i),eee(i),Repeat(' ',8),lf
              j=j+1
              Do k=1,recsz
                 chbat(k,j) = record(k:k)
              End Do

! Dump batch and update start of file

              If (j >= batsz .or. i == jatms) Then
                 Write(Unit=nhist, Fmt='(53a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,j)
                 rec=rec+Int(j,ip)
                 j=0
              End If
           End Do
        End Do

! Update main header

        Write(Unit=nhist, Fmt='(i10,2i21,a1)', Rec=Int(2,ip)) megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Do i=1,natms
           k=2*i

           ddd(i)=Sqrt(ravval(27+k-1))

           If (Abs(dofsit(lsite(i))) > zero_plus) tmp=weight(i)*ravval(27+k)/(boltz*3.0_wp)
           eee(i)=tmp
        End Do

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(ddd,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(eee,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

        rec=rec+Int(megatm+1,ip)

     End If

     Deallocate (chbuf,iwrk,ddd,eee, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write deallocation failure, node: ', idnode
        Call error(0)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname,  MPI_MODE_WRONLY, fh )

     Allocate (ddd(1:mxatms),eee(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write allocation failure, node: ', idnode
        Call error(0)
     End If

     Do i=1,natms
        k=2*i

        ddd(i)=Sqrt(ravval(27+k-1))

        If (Abs(dofsit(lsite(i))) > zero_plus) tmp=weight(i)*ravval(27+k)/(boltz*3.0_wp)
        eee(i)=tmp
     End Do

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and update start of file

        rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
        If (idnode == 0) Then
           Write(record, Fmt='(a8,2i10,2f12.6,a1)') 'timestep',nstep,megatm,tstep,time,lf
           Call io_write_record( fh, rec_mpi_io, record )
        End If
        rec=rec+Int(1,ip)
        rec_mpi_io=Int(rec,MPI_OFFSET_KIND)

     Else ! netCDF write

! Get the current and new frame numbers

        Call io_nc_get_dim( 'frame', fh, i )
        i = i + 1
        rec_mpi_io = Int(i,MPI_OFFSET_KIND)

        Call io_nc_put_var( 'time'           , fh,   time, i, 1 )
        Call io_nc_put_var( 'step'           , fh,  nstep, i, 1 )
        Call io_nc_put_var( 'timestep '      , fh,  tstep, i, 1 )

     End If

     Call io_write_sorted_file( fh, 0, IO_MSDTMP, rec_mpi_io, natms, &
          ltg, atmnam, ddd, eee, (/ 0.0_wp /),                       &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),  ierr )

     If ( ierr /= 0 ) Then
        Select Case( ierr )
        Case( IO_BASE_COMM_NOT_SET )
           Call error( 1050 )
        Case( IO_ALLOCATION_ERROR )
           Call error( 1053 )
        Case( IO_UNKNOWN_WRITE_OPTION )
           Call error( 1056 )
        Case( IO_UNKNOWN_WRITE_LEVEL )
           Call error( 1059 )
        End Select
     End If

! Update and save offset pointer

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec=rec+Int(megatm,ip)
        If (idnode == 0) Then
           Write(record, Fmt='(i10,2i21,a1)') megatm,frm,rec,lf
           Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record )
        End If
     End If

     Call io_close( fh )
     Call io_finalize

     Deallocate (ddd,eee, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write deallocation failure, node: ', idnode
        Call error(0)
     End If

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),ddd(1:mxatms),eee(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O

     If (idnode == 0) Then

        Open(Unit=nhist, File='fname', Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,ip)
        Write(Unit=nhist, Fmt='(a8,2i10,2f12.6,a1)', Rec=rec) 'timestep',nstep,megatm,tstep,time,lf

        Do i=1,natms
           k=2*i

           iwrk(i)=ltg(i)
           chbuf(i)=atmnam(i)

           ddd(i)=Sqrt(ravval(27+k-1))

           If (Abs(dofsit(lsite(i))) > zero_plus) tmp=weight(i)*ravval(27+k)/(boltz*3.0_wp)
           eee(i)=tmp
        End Do

        jatms=natms

        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(ddd,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(eee,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
           End If

           Do i=1,jatms
              Write(Unit=nhist, Fmt='(a8,i10,1p,2e13.4,a8,a1)', Rec=rec+Int(iwrk(i),ip)) &
                   chbuf(i),iwrk(i),ddd(i),eee(i),Repeat(' ',8),lf
           End Do
        End Do

! Update and save offset pointer

        rec=rec+Int(megatm,ip)
        Write(Unit=nhist, Fmt='(i10,2i21,a1)', Rec=Int(2,ip)) megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Do i=1,natms
           k=2*i

           ddd(i)=Sqrt(ravval(27+k-1))

           If (Abs(dofsit(lsite(i))) > zero_plus) tmp=weight(i)*ravval(27+k)/(boltz*3.0_wp)
           eee(i)=tmp
        End Do

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(ddd,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(eee,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

        rec=rec+Int(megatm+1,ip)

     End If

     Deallocate (chbuf,iwrk,ddd,eee, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write deallocation failure, node: ', idnode
        Call error(0)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'msd_write deallocation failure 0, node: ', idnode
        Call error(0)
     End If
  End If

  If (mxnode > 1) Call gsync()

End Subroutine msd_write
