Subroutine trajectory_write &
           (imcon,keyres,nstraj,istraj,keytrj,megatm,nstep,tstep,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing history file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2013
! contrib   - w.smith, i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use config_module,     Only : cfgname,cell,natms,     &
                                ltg,atmnam,chge,weight, &
                                xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use statistics_module, Only : rsd
  Use parse_module,      Only : tabs_2_blanks, get_word, &
                                strip_blanks, word_2_real
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
                                IO_HISTORY,IO_HISTORD,         &
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

  Integer,           Intent( In    ) :: imcon,keyres,         &
                                        nstraj,istraj,keytrj, &
                                        megatm,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep,time


  Logical,               Save :: newjob = .true. , &
                                 fast   = .true.
  Character( Len = 40 ), Save :: fname
  Integer,               Save :: recsz  = 73 ! default record size
  Integer(Kind=ip),      Save :: rec    = 0_ip , &
                                 frm    = 0_ip

  Logical                :: lexist,safe,ready
  Character( Len = 40 )  :: word
  Integer                :: fail(1:5),i,jj,k,jdnode,jatms
  Integer(Kind=ip)       :: rec1
  Real( Kind = wp )      :: buffer(1:2)
  Real( Kind = wp )      :: celprp(1:10),cell_vecs(1:3,1:3)
  Real( Kind = wp )      :: lengths(1:3), angles(1:3)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: fh, io_write, batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io,jj_io
  Character                         :: lf
  Character( Len = 200 )            :: record

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: cxx,cyy,czz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: ddd,eee,fff

! netCDF check

  Integer :: file_p, file_r
  Integer :: io_p, io_r


  If (.not.(nstep >= nstraj .and. Mod(nstep-nstraj,istraj) == 0)) Return

! Get write method, buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

  If (keytrj == 3) Then
     recsz = 35 ! record size for compressed HISTORY file
     Go To 100
  End If

  If (newjob) Then
     newjob = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = 'HISTORY'
     Else
        fname = 'HISTORY.nc'
     End If

! If keyres=1, is HISTORY old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (idnode == 0) Inquire(File=fname, Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

10   Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (idnode == 0) Then
              Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Status='replace', Recl=recsz)
              Write(Unit=nhist, Fmt='(a72,a1)',       Rec=Int(1,ip)) cfgname(1:72),lf
              Write(Unit=nhist, Fmt='(3i10,2i21,a1)', Rec=Int(2,ip)) keytrj,imcon,megatm,frm,rec,lf
              Close(Unit=nhist)
           End If
           rec=Int(2,ip)
           frm=Int(0,ip)
        Else
           If (idnode == 0) Then
              Call io_set_parameters( user_comm = MPI_COMM_SELF )
              Call io_nc_create( MPI_COMM_SELF, fname, cfgname, megatm )
           End If
        End If

! Get some sense out of it

     Else

        safe=.true.
        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           If (idnode == 0) Then

              Open(Unit=nhist, File=fname, Form='formatted')

              Do

                 record(1:recsz)=' '

! Assume new style of HISTORY with bookkeeping.

                 If (fast) Then

                    Read(Unit=nhist, Fmt=*, End=20)                     ! title record
                    rec=rec+Int(1,ip)
                    Read(Unit=nhist, Fmt='(a)', End=20) record(1:recsz) ! bookkeeping record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    rec=rec+Int(1,ip)

                    Call get_word(record(1:recsz),word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record(1:recsz),word) ; Call get_word(record(1:recsz),word)
                       Call get_word(record(1:recsz),word) ; frm=Nint(word_2_real(word,0.0_wp),ip)
                       Call get_word(record(1:recsz),word) ; rec=Nint(word_2_real(word,0.0_wp),ip)
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
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    rec=rec+Int(1,ip)

                    Call get_word(record(1:recsz),word) ; Call get_word(record(1:recsz),word)
                    Call get_word(record(1:recsz),word) ; jj=Nint(word_2_real(word))
                    Call get_word(record(1:recsz),word) ; k=Nint(word_2_real(word))

                    word=' '
                    i = 3 + (2+k)*jj ! total number of lines to read
                    Write(word,'( "(", i0, "( / ) )" )') i-1
                    Read(Unit=nhist, Fmt=word, End=20)
                    rec=rec+Int(i,ip)
                    frm=frm+Int(1,ip)

                 End If

              End Do

20            Continue
              Close(Unit=nhist)

           End If

           If (mxnode > 1) Call gcheck(safe,"enforce")
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

        Else ! netCDF read

           If (idnode == 0) Then
              Call io_set_parameters( user_comm = MPI_COMM_SELF )
              Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_RDONLY, fh )

! Get the precision that the history file was written in
! and check it matches the requested precision

              Call io_nc_get_file_real_precision( fh, file_p, file_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Write(nrite,'(/,1x,a)') &
  "Can not determine precision in an exisiting HISTORY.nc file in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync()
              Call error(0)
           End If

           If (idnode == 0) Then
              Call io_nc_get_real_precision( io_p, io_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Write(nrite,'(/,1x,a)') &
  "Can not determine the desired writing precision in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync()
              Call error(0)
           End If

           If (idnode == 0) safe = (io_p == file_p .and. io_r == file_r)
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Then
                 Write(nrite,'(/,1x,a)') &
  "Requested writing precision inconsistent with that in an existing HISTORY.nc"
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

              Call gsync()
              Call error(0)
           End If

! Get the frame number to check
! For netCDF this is the "frame number" which is not a long integer!

           If (idnode == 0) Call io_nc_get_dim( 'frame', fh, jj )
           Call MPI_BCAST(jj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

           If (jj > 0) Then
              frm=Int(jj,ip)

              If (idnode == 0) Call io_close( fh )
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
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write allocation failure 0, node: ', idnode
        Call error(0)
     End If

     chbat=' '
     n_atm=0 ; n_atm(idnode+1)=natms
     If (mxnode > 1) Call gsum(n_atm)
     n_atm(0)=Sum(n_atm(0:idnode))
  End If

! Notes:
! the MPI-I/O records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(Kind = MPI_OFFSET_KIND)

! Update frame

  frm=frm+Int(1,ip)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header and cell information, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     jj=0
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_WRONLY, fh )

        Write(record(1:recsz), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record(1:recsz), Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and cell information

        Call io_write_batch( fh, rec_mpi_io, jj, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        jj=jj+4

     End If
     Call gsync()

! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(jj,MPI_OFFSET_KIND)+Int(n_atm(0),MPI_OFFSET_KIND)*Int(keytrj+2,MPI_OFFSET_KIND)
     jj=0

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_WRONLY, fh )

     Do i=1,natms
        Write(record(1:recsz), Fmt='(a8,i10,3f12.6,a18,a1)') atmnam(i),ltg(i),weight(i),chge(i),rsd(i),Repeat(' ',18),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') xxx(i),yyy(i),zzz(i),Repeat(' ',12),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        If (keytrj > 0) Then
           Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') vxx(i),vyy(i),vzz(i),Repeat(' ',12),lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End If

        If (keytrj > 1) Then
           Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') fxx(i),fyy(i),fzz(i),Repeat(' ',12),lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End If

! Dump batch and update start of file

        If (jj + keytrj + 2 >= batsz .or. i == natms) Then
           Call io_write_batch( fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,MPI_OFFSET_KIND)
           jj=0
        End If
     End Do

! Update and save offset pointer

     rec=rec+Int(4,ip)+Int(megatm,ip)*Int(keytrj+2,ip)
     If (idnode == 0) Then
        Write(record(1:recsz), Fmt='(3i10,2i21,a1)') keytrj,imcon,megatm,frm,rec,lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record(1:recsz) )
     End If

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     Allocate (ddd(1:mxatms),eee(1:mxatms),fff(1:mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O
! Start of file

     jj=0
     If (idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

! Accumulate header

        Write(record(1:recsz), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Do i = 0, 2
           Write(record(1:recsz), Fmt='(3f20.10,a12,a1)') &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do
        End Do

! Dump header and update start of file

        Write(Unit=nhist, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,jj)
        rec=rec+Int(jj,ip)
        jj=0

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           chbuf(i)=atmnam(i)
           ddd(i)=chge(i)
           eee(i)=weight(i)
           fff(i)=rsd(i)
        End Do

        If (keytrj >= 1) Then
           Do i=1,natms
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)
           End Do
        End If

        If (keytrj >= 2) Then
           Do i=1,natms
              cxx(i)=fxx(i)
              cyy(i)=fyy(i)
              czz(i)=fzz(i)
           End Do
        End If

        jatms=natms
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(ddd,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(eee,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(fff,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 If (keytrj > 0) Then
                    Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(byy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 End If

                 If (keytrj > 1) Then
                    Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(czz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 End If
              End If
           End If

           Do i=1,jatms
              Write(record(1:recsz), Fmt='(a8,i10,3f12.6,a18,a1)') chbuf(i),iwrk(i),eee(i),ddd(i),fff(i),Repeat(' ',18),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              If (keytrj >= 1) Then
                 Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,recsz
                    chbat(k,jj) = record(k:k)
                 End Do
              End If

              If (keytrj >= 2) Then
                 Write(record(1:recsz), Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,recsz
                    chbat(k,jj) = record(k:k)
                 End Do
              End If

! Dump batch and update start of file

              If (jj + keytrj + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nhist, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,jj)
                 rec=rec+Int(jj,ip)
                 jj=0
              End If
           End Do
        End Do

! Update main header

        Write(Unit=nhist, Fmt='(3i10,2i21,a1)', Rec=Int(2,ip)) keytrj,imcon,megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

           Call MPI_SEND(chge,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(weight,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(rsd,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

           If (keytrj > 0) Then
              Call MPI_SEND(vxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           End If

           If (keytrj > 1) Then
              Call MPI_SEND(fxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(fyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(fzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           End If
        End If

! Save offset pointer

        rec=rec+Int(4,ip)+Int(megatm,ip)*Int(keytrj+2,ip)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     Deallocate (ddd,eee,fff, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write deallocation failure, node: ', idnode
        Call error(0)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     jj_io=rec_mpi_io
     jj=0 ! netCDF current frame
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_WRONLY, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and cell information

           Write(record(1:recsz), Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf
           Call io_write_record( fh, jj_io, record(1:recsz) )
           jj_io=jj_io + Int(1,MPI_OFFSET_KIND)

           Do i = 0, 2
              Write(record(1:recsz), Fmt='(3f20.10,a12,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              Call io_write_record( fh, jj_io, record(1:recsz) )
              jj_io=jj_io+Int(1,MPI_OFFSET_KIND)
           End Do

        Else ! netCDF write

! Get the current and new frame numbers

           Call io_nc_get_dim( 'frame', fh, jj )
           jj = jj + 1

           Call io_nc_put_var( 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var( 'datalevel'      , fh, keytrj, jj, 1 )
           Call io_nc_put_var( 'imageconvention', fh,  imcon, jj, 1 )
           Call io_nc_put_var( 'timestep '      , fh,  tstep, jj, 1 )

           Call dcell(cell,celprp) ! get cell properties

           cell_vecs = Reshape( cell, (/ 3, 3 /) )

           lengths( 1 ) = celprp( 1 )
           lengths( 2 ) = celprp( 2 )
           lengths( 3 ) = celprp( 3 )

           angles ( 1 ) = Acos( celprp( 5 ) )
           angles ( 2 ) = Acos( celprp( 6 ) )
           angles ( 3 ) = Acos( celprp( 4 ) )
           angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

! Print

           Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, jj /), (/ 3, 3, 1 /) )
           Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, jj /), (/    3, 1 /) )
           Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, jj /), (/    3, 1 /) )

        End If

        Call io_close( fh )
        Call io_finalize

     End If
     Call gsync()

! Start of file

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(4,MPI_OFFSET_KIND)
     Else ! netCDF write
        Call MPI_BCAST(jj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        rec_mpi_io = Int(jj,MPI_OFFSET_KIND)
     End If

! Write the rest

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_WRONLY, fh )

     Call io_write_sorted_file( fh, keytrj, IO_HISTORY, rec_mpi_io, natms, &
          ltg, atmnam, weight, chge, rsd, xxx, yyy, zzz,                   &
          vxx, vyy, vzz, fxx, fyy, fzz, ierr )

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
        rec=rec+Int(4,ip)+Int(megatm,ip)*Int(keytrj+2,ip)
        If (idnode == 0) Then
           Write(record(1:recsz), Fmt='(3i10,2i21,a1)') keytrj,imcon,megatm,frm,rec,lf
           Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record(1:recsz) )
        End If
     End If

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     Allocate (ddd(1:mxatms),eee(1:mxatms),fff(1:mxatms), Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O

     If (idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,ip)
        Write(Unit=nhist, Fmt='(a8,2i10,2i2,2f20.6,a1)') 'timestep',nstep,megatm,keytrj,imcon,tstep,time,lf

        Do i = 0, 2
           rec=rec+Int(1,ip)
           Write(Unit=nhist, Fmt='(3f20.10,a12,a1)', Rec=rec) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
        End Do

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           chbuf(i)=atmnam(i)
           ddd(i)=chge(i)
           eee(i)=weight(i)
           fff(i)=rsd(i)
        End Do

        If (keytrj >= 1) Then
           Do i=1,natms
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)
           End Do
        End If

        If (keytrj >= 2) Then
           Do i=1,natms
              cxx(i)=fxx(i)
              cyy(i)=fyy(i)
              czz(i)=fzz(i)
           End Do
        End If

        jatms=natms
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(ddd,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(eee,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(fff,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 If (keytrj > 0) Then
                    Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(byy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 End If

                 If (keytrj > 1) Then
                    Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(czz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 End If
              End If
           End If

           Do i=1,jatms
              rec1=rec+Int(iwrk(i)-1,ip)*Int(keytrj+2,ip)+Int(1,ip)
              Write(Unit=nhist, Fmt='(a8,i10,3f12.6,a18,a1)', Rec=rec1) chbuf(i),iwrk(i),eee(i),ddd(i),fff(i),Repeat(' ',18),lf

              rec1=rec1+Int(1,ip)
              Write(Unit=nhist, Fmt='(3g20.10,a12,a1)', Rec=rec1) axx(i),ayy(i),azz(i),Repeat(' ',12),lf

              If (keytrj >= 1) Then
                 rec1=rec1+Int(1,ip)
                 Write(Unit=nhist, Fmt='(3g20.10,a12,a1)', Rec=rec1) bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
              End If

              If (keytrj >= 2) Then
                 rec1=rec1+Int(1,ip)
                 Write(Unit=nhist, Fmt='(3g20.10,a12,a1)', Rec=rec1) cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
              End If
           End Do
        End Do

! Update main header

        rec=rec+Int(megatm,ip)*Int(keytrj+2,ip)
        Write(Unit=nhist, Fmt='(3i10,2i21,a1)', Rec=Int(2,ip)) keytrj,imcon,megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

           Call MPI_SEND(chge,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(weight,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(rsd,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

           If (keytrj > 0) Then
              Call MPI_SEND(vxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           End If

           If (keytrj > 1) Then
              Call MPI_SEND(fxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(fyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
              Call MPI_SEND(fzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           End If
        End If

! Save offset pointer

        rec=rec+Int(4,ip)+Int(megatm,ip)*Int(keytrj+2,ip)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     Deallocate (ddd,eee,fff, Stat=fail(5))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write deallocation failure, node: ', idnode
        Call error(0)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write deallocation failure 0, node: ', idnode
        Call error(0)
     End If
  End If

  If (mxnode > 1) Call gsync()

  Return

100 Continue

  If (newjob) Then
     newjob = .false.

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = 'HISTORY'
     Else
        fname = 'HISTORY.nc'
     End If

! If keyres=1, is HISTORY old (does it exist) and
! how many frames and records are in there

     lexist=.true.
     If (keyres == 1) Then
        If (idnode == 0) Inquire(File=fname, Exist=lexist)
        If (mxnode > 1) Call gcheck(lexist,"enforce")
     Else
        lexist=.false.
     End If

! Generate file is non-existent

110  Continue
     If (.not.lexist) Then

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           If (idnode == 0) Then
              Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Status='replace', Recl=recsz)
              Write(Unit=nhist, Fmt='(a34,a1)',      Rec=Int(1,ip)) cfgname(1:34),lf
              Write(Unit=nhist, Fmt='(2i2,3i10,a1)', Rec=Int(2,ip)) keytrj,imcon,megatm,frm,rec,lf
              Close(Unit=nhist)
           End If
           rec=Int(2,ip)
           frm=Int(0,ip)
        Else
           If (idnode == 0) Then
              Call io_set_parameters( user_comm = MPI_COMM_SELF )
              Call io_nc_create( MPI_COMM_SELF, fname, cfgname, megatm )
           End If
        End If

! Get some sense of it

     Else

        safe=.true.
        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

           If (idnode == 0) Then

              Open(Unit=nhist, File=fname, Form='formatted')

              Do

                 record(1:recsz)=' '

! Assume new style of HISTORY with bookkeeping.

                 If (fast) Then

                    Read(Unit=nhist, Fmt=*, End=120)                     ! title record
                    rec=rec+Int(1,ip)
                    Read(Unit=nhist, Fmt='(a)', End=120) record(1:recsz) ! bookkeeping record
                    Call tabs_2_blanks(record) ; Call strip_blanks(record)
                    rec=rec+Int(1,ip)

                    Call get_word(record(1:recsz),word)
                    If (word(1:Len_Trim(word)) /= 'timestep') Then
                       Call get_word(record(1:recsz),word) ; Call get_word(record(1:recsz),word)
                       Call get_word(record(1:recsz),word) ; frm=Nint(word_2_real(word,0.0_wp),ip)
                       Call get_word(record(1:recsz),word) ; rec=Nint(word_2_real(word,0.0_wp),ip)
                       If (frm /= Int(0,ip) .and. rec > Int(2,ip)) Then
                          Go To 120 ! New style
                       Else
                          fast=.false. ! TOUGH, old style
                          rec=Int(2,ip)
                          frm=Int(0,ip)
                       End If
                    Else
                       safe=.false. ! Overwrite the file, it's junk to me
                       Go To 120
                    End If

! TOUGH, it needs scanning through

                 Else

                    Read(Unit=nhist, Fmt=*, End=120)                 ! timestep record
                    rec=rec+Int(1,ip)

                    word=' '
                    i = 3 + megatm ! total number of lines to read
                    Write(word,'( "(", i0, "( / ) )" )') i-1
                    Read(Unit=nhist, Fmt=word, End=120)
                    rec=rec+Int(i,ip)
                    frm=frm+Int(1,ip)

                 End If

              End Do

120         Continue
            Close(Unit=nhist)

           End If

           If (mxnode > 1) Call gcheck(safe,"enforce")
           If (.not.safe) Then
              lexist=.false.

              rec=Int(0,ip)
              frm=Int(0,ip)

              Go To 110
           Else If (mxnode > 1) Then
              buffer(1)=Real(frm,wp)
              buffer(2)=Real(rec,wp)

              Call MPI_BCAST(buffer(1:2), 2, wp_mpi, 0, dlp_comm_world, ierr)

              frm=Nint(buffer(1),ip)
              rec=Nint(buffer(2),ip)
           End If

        Else ! netCDF read

           If (idnode == 0) Then
              Call io_set_parameters( user_comm = MPI_COMM_SELF )
              Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_RDONLY, fh )

! Get the precision that the history file was written in
! and check it matches the requested precision

              Call io_nc_get_file_real_precision( fh, file_p, file_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Write(nrite,'(/,1x,a)') &
  "Can not determine precision in an exisiting HISTORY.nc file in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync()
              Call error(0)
           End If

           If (idnode == 0) Then
              Call io_nc_get_real_precision( io_p, io_r, ierr )
              safe = (ierr == 0)
           End If
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Write(nrite,'(/,1x,a)') &
  "Can not determine the desired writing precision in trajectory_write"

! Sync before killing for the error in the hope that something sensible happens

              Call gsync()
              Call error(0)
           End If

           If (idnode == 0) safe = (io_p == file_p .and. io_r == file_r)
           Call gcheck(safe)
           If (.not.safe) Then
              If (idnode == 0) Then
                 Write(nrite,'(/,1x,a)') &
  "Requested writing precision inconsistent with that in an existing HISTORY.nc"
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

              Call gsync()
              Call error(0)
           End If

! Get the frame number to check
! For netCDF this is the "frame number" which is not a long integer!

           jj=0
           If (idnode == 0) Call io_nc_get_dim( 'frame', fh, jj )
           Call MPI_BCAST(jj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

           If (jj > 0) Then
              frm=Int(jj,ip)

           If (idnode == 0) Call io_close( fh )
           Else ! Overwrite the file, it's junk to me
              lexist=.false.

              rec=Int(0,ip)
              frm=Int(0,ip)

              Go To 110
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
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write allocation failure 0, node: ', idnode
        Call error(0)
     End If

     chbat=' '
     n_atm=0 ; n_atm(idnode+1)=natms
     If (mxnode > 1) Call gsum(n_atm)
     n_atm(0)=Sum(n_atm(0:idnode))
  End If

! Update frame

  frm=frm+Int(1,ip)

! NO UNSORTED WRITING AS NO INDICES ARE WRITTEN
! In other words if the user asked for unsorted I/O they get it sorted.

! SORTED MPI-I/O or Parallel Direct Access FORTRAN

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_SORTED_MPIIO    .or. &
      io_write == IO_WRITE_SORTED_DIRECT   .or. &
      io_write == IO_WRITE_SORTED_NETCDF) Then

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_WRONLY, fh )

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(rec,MPI_OFFSET_KIND)
     jj_io=rec_mpi_io
     jj=0 ! netCDF current frame
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_WRONLY, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header and cell information

           Write(record(1:recsz), Fmt='(a8,i6,f8.5,f12.5,a1)') 'timestep',nstep,tstep,time,lf
           Call io_write_record( fh, jj_io, record(1:recsz) )
           jj_io=jj_io + Int(1,MPI_OFFSET_KIND)

           Do i = 0, 2
              Write(record(1:recsz), Fmt='(3f10.3,a4,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 4 ), lf
              Call io_write_record( fh, jj_io, record(1:recsz) )
              jj_io=jj_io+Int(1,MPI_OFFSET_KIND)
           End Do

        Else ! netCDF write

! Get the current and new frame numbers

           Call io_nc_get_dim( 'frame', fh, jj )
           jj = jj + 1

           Call io_nc_put_var( 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var( 'timestep '      , fh,  tstep, jj, 1 )

           Call dcell(cell,celprp) ! get cell properties

           cell_vecs = Reshape( cell, (/ 3, 3 /) )

           lengths( 1 ) = celprp( 1 )
           lengths( 2 ) = celprp( 2 )
           lengths( 3 ) = celprp( 3 )

           angles ( 1 ) = Acos( celprp( 5 ) )
           angles ( 2 ) = Acos( celprp( 6 ) )
           angles ( 3 ) = Acos( celprp( 4 ) )
           angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

! Print

           Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, jj /), (/ 3, 3, 1 /) )
           Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, jj /), (/    3, 1 /) )
           Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, jj /), (/    3, 1 /) )

        End If

        Call io_close( fh )
        Call io_finalize

     End If
     Call gsync()

! Start of file

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        rec_mpi_io=Int(rec,MPI_OFFSET_KIND)+Int(4,MPI_OFFSET_KIND)
     Else ! netCDF write
        Call MPI_BCAST(jj, 1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        rec_mpi_io = Int(jj,MPI_OFFSET_KIND)
     End If

! Write the rest

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_WRONLY, fh )

     Call io_write_sorted_file( fh, 0*keytrj, IO_HISTORD, rec_mpi_io, natms, &
          ltg, atmnam,  (/ 0.0_wp /),  (/ 0.0_wp /), rsd, xxx, yyy, zzz,     &
          (/ 0.0_wp /),  (/ 0.0_wp /),  (/ 0.0_wp /),                        &
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
        rec=rec+Int(4,ip)+Int(megatm,ip)
        If (idnode == 0) Then
           Write(record(1:recsz), Fmt='(2i2,3i10,a1)') keytrj,imcon,megatm,frm,rec,lf
           Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record(1:recsz) )
        End If
     End If

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (fff(1:mxatms),                             Stat=fail(3))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O

     If (idnode == 0) Then

        Open(Unit=nhist, File=fname, Form='formatted', Access='direct', Recl=recsz)

        rec=rec+Int(1,ip)
        Write(Unit=nhist, Fmt='(a8,i6,f8.5,f12.5,a1)', Rec=rec) 'timestep',nstep,tstep,time,lf

        Do i = 0, 2
           rec=rec+Int(1,ip)
           Write(Unit=nhist, Fmt='(3f10.3,a4,a1)', Rec=rec) &
                cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 4 ), lf
        End Do

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           chbuf(i)=atmnam(i)
           fff(i)=rsd(i)
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(fff,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              End If
           End If

           Do i=1,jatms
              rec1=rec+Int(iwrk(i),ip)
              Write(Unit=nhist, Fmt='(a6,4f7.1,a1)', Rec=rec1) chbuf(i),axx(i),ayy(i),azz(i),fff(i),lf
           End Do
        End Do

! Update main header

        rec=rec+Int(megatm,ip)
        Write(Unit=nhist, Fmt='(2i2,3i10,a1)', Rec=Int(2,ip)) keytrj,imcon,megatm,frm,rec,lf

        Close(Unit=nhist)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

           Call MPI_SEND(rsd,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        End If

! Save offset pointer

        rec=rec+Int(4,ip)+Int(megatm,ip)

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (fff,         Stat=fail(3))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write deallocation failure, node: ', idnode
        Call error(0)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'trajectory_write deallocation failure 0, node: ', idnode
        Call error(0)
     End If
  End If

  If (mxnode > 1) Call gsync()

End Subroutine trajectory_write
