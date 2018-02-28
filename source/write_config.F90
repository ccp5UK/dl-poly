Subroutine write_config(name,levcfg,megatm,nstep,tstep,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing configuration file
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, li
  Use comms_module
  Use setup_module
  Use config_module, Only : cfgname,imcon,cell,natms,ltg,atmnam, &
                            xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use io_module,     Only : io_set_parameters,         &
                            io_get_parameters,         &
                            io_init, io_nc_create,     &
                            io_open, io_write_record,  &
                            io_write_batch,            &
                            io_nc_put_var,             &
                            io_write_sorted_file,      &
                            io_delete,                 &
                            io_close, io_finalize,     &
                            IO_RESTART,                &
                            IO_BASE_COMM_NOT_SET,      &
                            IO_ALLOCATION_ERROR,       &
                            IO_UNKNOWN_WRITE_OPTION,   &
                            IO_UNKNOWN_WRITE_LEVEL,    &
                            IO_WRITE_UNSORTED_MPIIO,   &
                            IO_WRITE_UNSORTED_DIRECT,  &
                            IO_WRITE_UNSORTED_MASTER,  &
                            IO_WRITE_SORTED_MPIIO,     &
                            IO_WRITE_SORTED_DIRECT,    &
                            IO_WRITE_SORTED_NETCDF,    &
                            IO_WRITE_SORTED_MASTER

  Implicit None

  Character( Len = * ), Intent( In    ) :: name
  Integer,              Intent( In    ) :: levcfg,megatm,nstep
  Real( Kind = wp ),    Intent( In    ) :: tstep,time

  Integer, Parameter :: recsz    =   73 ! default record size

  Logical               :: ready
  Character( Len = 40 ) :: fname
  Integer(Kind=li)      :: rec,rec1     ! record line

  Integer               :: fail(1:4),i,k,jj,jdnode,jatms

  Real( Kind = wp )     :: celprp(1:10),cell_vecs(1:3,1:3)
  Real( Kind = wp )     :: lengths(1:3), angles(1:3)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: fh
  Integer                           :: io_write,batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: bxx,byy,bzz
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: cxx,cyy,czz


! Get write method buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! Get offsets and define batch

  fail=0
  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Allocate (n_atm(0:mxnode),        Stat=fail(1))
     Allocate (chbat(1:recsz,1:batsz), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'write_config allocation failure 0, node: ', idnode
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

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(1-1,MPI_OFFSET_KIND)
     jj=0
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_delete( name ) ! Sort existence issues
        Call io_open( io_write, MPI_COMM_SELF, name, MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(4i10,1p,2e16.7,a1)') levcfg,imcon,megatm,nstep,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

! Accumulate header - optional cell information (if present)

        If (imcon > 0) Then
           Do i = 0, 2
              Write(record, Fmt='(3f20.10,a12,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do
           End Do
        End If

! Dump header

        Call io_write_batch( fh, rec_mpi_io, jj, chbat )

        Call io_close( fh )
        Call io_finalize

     Else

        jj=jj+2
        If (imcon > 0) jj=jj+3

     End If
     Call gsync()

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_delete( name )
     Call io_open( io_write, dlp_comm_world, name, MPI_MODE_WRONLY, fh )

! Start of file (updated)

     rec_mpi_io=Int(jj,MPI_OFFSET_KIND)+Int(n_atm(0),MPI_OFFSET_KIND)*Int(levcfg+2,MPI_OFFSET_KIND)
     jj=0
     Do i=1,natms
        Write(record, Fmt='(a8,i10,a54,a1)') atmnam(i),ltg(i),Repeat(' ',54),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') xxx(i),yyy(i),zzz(i),Repeat(' ',12),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        If (levcfg > 0) Then
           Write(record, Fmt='(3g20.10,a12,a1)') vxx(i),vyy(i),vzz(i),Repeat(' ',12),lf
           jj=jj+1
           Do k=1,recsz
              chbat(k,jj) = record(k:k)
           End Do

           If (levcfg > 1) Then
              Write(record, Fmt='(3g20.10,a12,a1)') fxx(i),fyy(i),fzz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do
           End If
        End If

! Dump batch and update start of file

        If (jj + levcfg + 2 >= batsz .or. i == natms) Then
           Call io_write_batch( fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,MPI_OFFSET_KIND)
           jj=0
        End If
     End Do

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'write_config allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O
! Start of file

     rec=Int(0,li)
     jj=0
     If (idnode == 0) Then

! Write configuration data to new configuration file

        Open(Unit=nconf, File=name, Form='formatted', Access='direct', Recl=recsz, Status='replace')

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(4i10,1p,2e16.7,a1)') levcfg,imcon,megatm,nstep,tstep,time,lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

! Accumulate header - optional cell information (if present)

        If (imcon > 0) Then
           Do i = 0, 2
              Write(record, Fmt='(3f20.10,a12,a1)') &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do
           End Do
        End If

! Dump header and update start of file

        Write(Unit=nconf, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
        rec=Int(jj,li)
        jj=0

        Do i=1,natms
           iwrk(i)=ltg(i)
           chbuf(i)=atmnam(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           If (levcfg > 0) Then
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)

              If (levcfg > 1) Then
                 cxx(i)=fxx(i)
                 cyy(i)=fyy(i)
                 czz(i)=fzz(i)
              End If
           End If
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,WriteConf_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)

                 If (levcfg > 0) Then
                    Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(byy,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)

                    If (levcfg > 1) Then
                       Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                       Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                       Call MPI_RECV(czz,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                    End If
                 End If
              End If
           End If

           jj=0
           Do i=1,jatms
              Write(record, Fmt='(a8,i10,a54,a1)') atmnam(i),iwrk(i),Repeat(' ',54),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              If (levcfg > 0) Then
                 Write(record, Fmt='(3g20.10,a12,a1)') bxx(i),byy(i),bzz(i),Repeat(' ',12),lf
                 jj=jj+1
                 Do k=1,recsz
                    chbat(k,jj) = record(k:k)
                 End Do

                 If (levcfg > 1) Then
                    Write(record, Fmt='(3g20.10,a12,a1)') cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
                    jj=jj+1
                    Do k=1,recsz
                       chbat(k,jj) = record(k:k)
                    End Do
                 End If
              End If

! Dump batch and update start of file

              If (jj + levcfg + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nconf, Fmt='(73a)', Rec=rec+Int(1,li)) (chbat(:,k), k=1,jj)
                 rec=rec+Int(jj,li)
                 jj=0
              End If
           End Do
        End Do

        Close(Unit=nconf)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,WriteConf_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,WriteConf_tag,dlp_comm_world,ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,WriteConf_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,WriteConf_tag,dlp_comm_world,ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)

           If (levcfg > 0) Then
              Call MPI_SEND(vxx,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vyy,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vzz,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)

              If (levcfg > 1) Then
                 Call MPI_SEND(fxx,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
                 Call MPI_SEND(fyy,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
                 Call MPI_SEND(fzz,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
              End If
           End If
        End If

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'write_config deallocation failure, node: ', idnode
        Call error(0)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! name convention

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        fname = name
     Else
        fname = name(1:Len_Trim(name)) // '.nc'
     End If

! Write header only at start, where just one node is needed
! Start of file

     rec_mpi_io=Int(1-1,MPI_OFFSET_KIND)
     jj=0
     If (idnode == 0) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_delete( fname ) ! Sort existence issues
        If (io_write == IO_WRITE_SORTED_NETCDF) Call io_nc_create( MPI_COMM_SELF, fname, cfgname, megatm )
        Call io_open( io_write, MPI_COMM_SELF, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

! Non netCDF

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Write header

           Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
           Call io_write_record( fh, Int(jj,MPI_OFFSET_KIND), record )
           jj=jj+1

           Write(record, Fmt='(4i10,1p,2e16.7,a1)') levcfg,imcon,megatm,nstep,tstep,time,lf
           Call io_write_record( fh, Int(jj,MPI_OFFSET_KIND), record )
           jj=jj+1

! Write optional cell information (if present)

           If (imcon > 0) Then
              Do i = 0, 2
                 Write( record, '( 3f20.10, a12, a1 )' ) &
                      cell( 1 + i * 3: 3 + i * 3 ), Repeat( ' ', 12 ), lf
                 Call io_write_record( fh, Int(jj,MPI_OFFSET_KIND), record )
                 jj=jj+1
              End Do
           End If

        Else ! netCDF write

           jj=1 ! For config there is only one frame

           Call io_nc_put_var( 'time'           , fh,   time, jj, 1 )
           Call io_nc_put_var( 'step'           , fh,  nstep, jj, 1 )
           Call io_nc_put_var( 'datalevel'      , fh, levcfg, jj, 1 )
           Call io_nc_put_var( 'imageconvention', fh,  imcon, jj, 1 )
           Call io_nc_put_var( 'timestep'       , fh,  tstep, jj, 1 )

           If (imcon > 0) Then
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

        End If

        Call io_close( fh )
        Call io_finalize

     Else

        If (io_write /= IO_WRITE_SORTED_NETCDF) Then
           jj=jj+2
           If (imcon > 0) jj=jj+3
        Else
           jj=1
        End If

     End If
     Call gsync()

! Write the rest

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fname, MPI_MODE_WRONLY, fh )

     rec_mpi_io=rec_mpi_io+Int(jj,MPI_OFFSET_KIND)
     Call io_write_sorted_file( fh, levcfg, IO_RESTART, rec_mpi_io, natms,      &
          ltg, atmnam, (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), xxx, yyy, zzz, &
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

     Call io_close( fh )
     Call io_finalize

! SORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_SORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms),            Stat=fail(1))
     Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
     Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
     Allocate (cxx(1:mxatms),cyy(1:mxatms),czz(1:mxatms), Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'write_config allocation failure, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O
! Start of file

     rec=Int(0,li)
     If (idnode == 0) Then

! Write configuration data to new configuration file

        Open(Unit=nconf, File=name, Form='formatted', Access='direct', Recl=recsz, Status='replace')

! Write header

        rec=rec+Int(1,li)
        Write(Unit=nconf, Fmt='(a72,a1)',            Rec=rec) cfgname(1:72),lf
        rec=rec+Int(1,li)
        Write(Unit=nconf, Fmt='(4i10,1p,2e16.7,a1)', Rec=rec) levcfg,imcon,megatm,nstep,tstep,time,lf

! Write optional cell information (if present)

        If (imcon > 0) Then
           Do i = 0, 2
              rec=rec+Int(1,li)
              Write(Unit=nconf, Fmt='(3f20.10,a12,a1)', Rec=rec) &
                   cell( 1 + i * 3 ), cell( 2 + i * 3 ), cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
           End Do
        End If

        Do i=1,natms
           iwrk(i)=ltg(i)
           chbuf(i)=atmnam(i)

           axx(i)=xxx(i)
           ayy(i)=yyy(i)
           azz(i)=zzz(i)

           If (levcfg > 0) Then
              bxx(i)=vxx(i)
              byy(i)=vyy(i)
              bzz(i)=vzz(i)

              If (levcfg > 1) Then
                 cxx(i)=fxx(i)
                 cyy(i)=fyy(i)
                 czz(i)=fzz(i)
              End If
           End If
        End Do

        jatms=natms
        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,WriteConf_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
              If (jatms > 0) Then
                 Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)

                 Call MPI_RECV(axx,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                 Call MPI_RECV(azz,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)

                 If (levcfg > 0) Then
                    Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(byy,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                    Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)

                    If (levcfg > 1) Then
                       Call MPI_RECV(cxx,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                       Call MPI_RECV(cyy,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                       Call MPI_RECV(czz,jatms,wp_mpi,jdnode,WriteConf_tag,dlp_comm_world,status,ierr)
                    End If
                 End If
              End If
           End If

           Do i=1,jatms
              rec1=rec+Int(iwrk(i)-1,li)*Int(levcfg+2)+Int(1,li)
              Write(Unit=nconf, Fmt='(a8,i10,a54,a1)',     Rec=rec1) chbuf(i),iwrk(i),Repeat(' ',54),lf
              rec1=rec1+Int(1,li)
              Write(Unit=nconf, Fmt='(3g20.10,a12,a1)',    Rec=rec1) axx(i),ayy(i),azz(i),Repeat(' ',12),lf

              If (levcfg > 0) Then
                 rec1=rec1+Int(1,li)
                 Write(Unit=nconf, Fmt='(3g20.10,a12,a1)', Rec=rec1) bxx(i),byy(i),bzz(i),Repeat(' ',12),lf

                 If (levcfg > 1) Then
                    rec1=rec1+Int(1,li)
                    Write(Unit=nconf, Fmt='(3g20.10,a12,a1)', Rec=rec1) cxx(i),cyy(i),czz(i),Repeat(' ',12),lf
                 End If
              End If
           End Do
        End Do

        Close(Unit=nconf)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,WriteConf_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,WriteConf_tag,dlp_comm_world,ierr)
        If (natms > 0) Then
           Call MPI_SEND(atmnam,8*natms,MPI_CHARACTER,0,WriteConf_tag,dlp_comm_world,ierr)
           Call MPI_SEND(ltg,natms,MPI_INTEGER,0,WriteConf_tag,dlp_comm_world,ierr)

           Call MPI_SEND(xxx,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
           Call MPI_SEND(yyy,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
           Call MPI_SEND(zzz,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)

           If (levcfg > 0) Then
              Call MPI_SEND(vxx,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vyy,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
              Call MPI_SEND(vzz,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)

              If (levcfg > 1) Then
                 Call MPI_SEND(fxx,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
                 Call MPI_SEND(fyy,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
                 Call MPI_SEND(fzz,natms,wp_mpi,0,WriteConf_tag,dlp_comm_world,ierr)
              End If
           End If
        End If

     End If

     Deallocate (chbuf,iwrk,  Stat=fail(1))
     Deallocate (axx,ayy,azz, Stat=fail(2))
     Deallocate (bxx,byy,bzz, Stat=fail(3))
     Deallocate (cxx,cyy,czz, Stat=fail(4))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'write_config deallocation failure, node: ', idnode
        Call error(0)
     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
      io_write == IO_WRITE_UNSORTED_DIRECT .or. &
      io_write == IO_WRITE_UNSORTED_MASTER) Then
     Deallocate (n_atm, Stat=fail(1))
     Deallocate (chbat, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'write_config deallocation failure 0, node: ', idnode
        Call error(0)
     End If
  End If

  If (mxnode > 1) Call gsync()

End Subroutine write_config
