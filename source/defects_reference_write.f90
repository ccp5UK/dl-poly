Subroutine defects_reference_write(name,imcon,megref,nrefs,namr,indr,xr,yr,zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for re-writing the reference file needed for
!           defects detection
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2010
! contrib   - i.j.bush
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use config_module, Only : cfgname,cell
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
                            IO_WRITE_SORTED_MASTER,    &
                            IO_READ_NETCDF

  Implicit None

  Character( Len = * ), Intent( In    ) :: name
  Integer,              Intent( In    ) :: imcon,megref,nrefs
  Character( Len = 8 ), Intent( In    ) :: namr(1:mxatms)
  Integer,              Intent( In    ) :: indr(1:mxatms)
  Real( Kind = wp ),    Intent( In    ) :: xr(1:mxatms),yr(1:mxatms),zr(1:mxatms)

  Integer, Parameter :: recsz = 73 ! default record size

  Logical               :: ready
  Integer(Kind=ip)      :: rec          ! record line

  Integer               :: fail(1:4),i,k,jj,jdnode,jatms

  Real( Kind = wp )     :: celprp(1:10),x,y,z,cell_vecs(1:3,1:3)
  Real( Kind = wp )     :: lengths(1:3), angles(1:3)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: fh
  Integer                           :: io_write,io_read,batsz
  Integer( Kind = MPI_OFFSET_KIND ) :: rec_mpi_io
  Character( Len = recsz )          :: record
  Character                         :: lf

  Character( Len = 1 ), Dimension( :, : ), Allocatable :: chbat
  Character( Len = 8 ), Dimension( : ),    Allocatable :: chbuf
  Integer,              Dimension( : ),    Allocatable :: iwrk,n_atm
  Real( Kind = wp ),    Dimension( : ),    Allocatable :: axx,ayy,azz


  fail=0
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_write allocation failure, node: ', idnode
     Call error(0)
  End If

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_method_read       = io_read  )
  Call io_get_parameters( user_buffer_size_write = batsz    )
  Call io_get_parameters( user_line_feed         = lf       )

! Nasty fix here as this is a REWRITE!!!

  If (io_read == IO_READ_NETCDF) io_write = IO_WRITE_SORTED_NETCDF

! Get offsets and define batch

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
     n_atm=0 ; n_atm(idnode+1)=nrefs
     If (mxnode > 1) Call gsum(n_atm)
     n_atm(0)=Sum(n_atm(0:idnode))
  End If

! Notes:
! the MPI_IO records are numbered from 0 (not 1)
! - the displacement (disp_mpi_io) in the MPI_FILE_SET_VIEW call, and
!   the record number (rec_mpi_io) in the MPI_WRITE_FILE_AT calls are
!   both declared as: Integer(kind = MPI_OFFSET_KIND)

! UNSORTED MPI-I/O or Parallel Direct Access FORTRAN

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_UNSORTED_DIRECT) Then

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_delete( name )
     Call io_open( io_write, dlp_comm_world, name, MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

! Start of file

     rec_mpi_io=Int(1-1,MPI_OFFSET_KIND)

! Accumulate header

     If (idnode == 0) Then
        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        Do k=1,recsz
           chbat(k,1) = record(k:k)
        End Do

        Write(record, Fmt='(3i10,a42,a1)') 0,imcon,megref,Repeat(' ',42),lf
        Do k=1,recsz
           chbat(k,2) = record(k:k)
        End Do

        Write(record, Fmt='(3f20.10,a12,a1)') cell(1),cell(2),cell(3),Repeat(' ',12),lf
        Do k=1,recsz
           chbat(k,3) = record(k:k)
        End Do

        Write(record, Fmt='(3f20.10,a12,a1)') cell(4),cell(5),cell(6),Repeat(' ',12),lf
        Do k=1,recsz
           chbat(k,4) = record(k:k)
        End Do

        Write(record, Fmt='(3f20.10,a12,a1)') cell(7),cell(8),cell(9),Repeat(' ',12),lf
        Do k=1,recsz
           chbat(k,5) = record(k:k)
        End Do
     End If
     jj=5

! Dump header and update start of file

     If (idnode == 0) Call io_write_batch( fh, rec_mpi_io, jj, chbat )
     rec_mpi_io=Int(jj,MPI_OFFSET_KIND)+Int(n_atm(0),MPI_OFFSET_KIND)*Int(2,MPI_OFFSET_KIND)

! Get to real space

     Do i=1,nrefs
        axx(i)=cell(1)*xr(i)+cell(4)*yr(i)+cell(7)*zr(i)
        ayy(i)=cell(2)*xr(i)+cell(5)*yr(i)+cell(8)*zr(i)
        azz(i)=cell(3)*xr(i)+cell(6)*yr(i)+cell(9)*zr(i)
     End Do

! DD bound

     Call pbcshift(imcon,cell,nrefs,axx,ayy,azz)

     jj=0
     Do i=1,nrefs
        Write(record, Fmt='(a8,i10,a54,a1)') namr(i),indr(i),Repeat(' ',54),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

        Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
        jj=jj+1
        Do k=1,recsz
           chbat(k,jj) = record(k:k)
        End Do

! Dump batch and update start of file

        If (jj + 2 >= batsz .or. i == nrefs) Then
           Call io_write_batch( fh, rec_mpi_io, jj, chbat )
           rec_mpi_io=rec_mpi_io+Int(jj,MPI_OFFSET_KIND)
           jj=0
        End If
     End Do

     Call io_close( fh )
     Call io_finalize

! UNSORTED Serial Direct Access FORTRAN

  Else If (io_write == IO_WRITE_UNSORTED_MASTER) Then

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_write allocation failure 1, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O

     If (idnode == 0) Then

! Obliterate old REFERENCE if not needed and print header

        Open(Unit=nrefdt, File=name, Form='formatted', Access='direct', Recl=recsz, Status='replace')

! Accumulate header

        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        Do k=1,recsz
           chbat(k,1) = record(k:k)
        End Do

        Write(record, Fmt='(3i10,a42,a1)') 0,imcon,megref,Repeat(' ',42),lf
        Do k=1,recsz
           chbat(k,2) = record(k:k)
        End Do

        Write(record, Fmt='(3f20.10,a12,a1)') cell(1),cell(2),cell(3),Repeat(' ',12),lf
        Do k=1,recsz
           chbat(k,3) = record(k:k)
        End Do

        Write(record, Fmt='(3f20.10,a12,a1)') cell(4),cell(5),cell(6),Repeat(' ',12),lf
        Do k=1,recsz
           chbat(k,4) = record(k:k)
        End Do

        Write(record, Fmt='(3f20.10,a12,a1)') cell(7),cell(8),cell(9),Repeat(' ',12),lf
        Do k=1,recsz
           chbat(k,5) = record(k:k)
        End Do

        jj = 5

! Dump header and update start of file

        Write(Unit=nrefdt, Fmt='(73a)', Rec=Int(1,ip)) (chbat(:,k), k=1,jj)
        rec=Int(5,ip)

! Get to real space

        Do i=1,nrefs
           chbuf(i)=namr(i)
           iwrk(i)=indr(i)

           axx(i)=cell(1)*xr(i)+cell(4)*yr(i)+cell(7)*zr(i)
           ayy(i)=cell(2)*xr(i)+cell(5)*yr(i)+cell(8)*zr(i)
           azz(i)=cell(3)*xr(i)+cell(6)*yr(i)+cell(9)*zr(i)
        End Do

        jatms=nrefs

        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

! Get to real space

              Do i=1,jatms
                 x=axx(i) ; y=ayy(i) ; z=azz(i)
                 axx(i)=cell(1)*x+cell(4)*y+cell(7)*z
                 ayy(i)=cell(2)*x+cell(5)*y+cell(8)*z
                 azz(i)=cell(3)*x+cell(6)*y+cell(9)*z
              End Do
           End If

! DD bound

           Call pbcshift(imcon,cell,jatms,axx,ayy,azz)

           jj=0
           Do i=1,jatms
              Write(record, Fmt='(a8,i10,a54,a1)') chbuf(i),iwrk(i),Repeat(' ',54),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

              Write(record, Fmt='(3g20.10,a12,a1)') axx(i),ayy(i),azz(i),Repeat(' ',12),lf
              jj=jj+1
              Do k=1,recsz
                 chbat(k,jj) = record(k:k)
              End Do

! Dump batch and update start of file

              If (jj + 2 >= batsz .or. i == jatms) Then
                 Write(Unit=nrefdt, Fmt='(73a)', Rec=rec+Int(1,ip)) (chbat(:,k), k=1,jj)
                 rec=rec+Int(jj,ip)
                 jj=0
              End If
           End Do
        End Do

! Close REFERENCE

        Close(Unit=nrefdt)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(nrefs,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(namr,8*nrefs,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(indr,nrefs,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(xr,nrefs,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(yr,nrefs,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(zr,nrefs,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

     End If

     Deallocate (chbuf,iwrk, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_write deallocation failure 1, node: ', idnode
        Call error(0)
     End If

! SORTED MPI-I/O or Parallel Direct Access FORTRAN or netCDF

  Else If (io_write == IO_WRITE_SORTED_MPIIO  .or. &
           io_write == IO_WRITE_SORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

! Write header only at start, where just one node is needed

     Call io_set_parameters( user_comm = MPI_COMM_SELF )
     Call io_init( recsz )

! Sort existence issues
     Call io_delete( name )
     If (io_write == IO_WRITE_SORTED_NETCDF) &
        Call io_nc_create( MPI_COMM_SELF, name, cfgname, megref )

     Call io_open( io_write, MPI_COMM_SELF, name, MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

! Non netCDF

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then

! Start of file

        rec_mpi_io=Int(1-1,MPI_OFFSET_KIND)

! Write header

        record=' '
        Write(record, Fmt='(a72,a1)') cfgname(1:72),lf
        Call io_write_record( fh, rec_mpi_io, record )
        rec_mpi_io=rec_mpi_io+Int(1,MPI_OFFSET_KIND)

        Write(record, Fmt='(3i10,a42,a1)') 0,imcon,megref,Repeat(' ',42),lf
        Call io_write_record( fh, rec_mpi_io, record )
        rec_mpi_io=rec_mpi_io+Int(1,MPI_OFFSET_KIND)

        Do i = 0, 2
           Write( record, '( 3f20.10, a12, a1 )' ) cell( 1 + i * 3: 3 + i * 3 ), Repeat( ' ', 12 ), lf
           Call io_write_record( fh, rec_mpi_io, record )
           rec_mpi_io=rec_mpi_io+Int(1,MPI_OFFSET_KIND)
        End Do

     Else ! netCDF write

        i=1 ! For config there is only one frame
        rec_mpi_io = Int(i,MPI_OFFSET_KIND) ! This is the "frame number" which is not a long integer!

        Call io_nc_put_var( 'time'           , fh, 0.0_wp, i, 1 )
        Call io_nc_put_var( 'step'           , fh,      0, i, 1 )
        Call io_nc_put_var( 'datalevel'      , fh,      0, i, 1 )
        Call io_nc_put_var( 'imageconvention', fh,  imcon, i, 1 )
        Call io_nc_put_var( 'timestep'       , fh, 0.0_wp, i, 1 )

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

        Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
        Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, i /), (/    3, 1 /) )
        Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, i /), (/    3, 1 /) )

     End If

     Call io_close( fh )
     Call io_finalize

     Call gsync()

! Get to real space

     Do i=1,nrefs
        axx(i)=cell(1)*xr(i)+cell(4)*yr(i)+cell(7)*zr(i)
        ayy(i)=cell(2)*xr(i)+cell(5)*yr(i)+cell(8)*zr(i)
        azz(i)=cell(3)*xr(i)+cell(6)*yr(i)+cell(9)*zr(i)
     End Do

! DD bound

     Call pbcshift(imcon,cell,nrefs,axx,ayy,azz)

! Write the rest

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, name, MPI_MODE_WRONLY, fh )

     Call io_write_sorted_file( fh, 0, IO_RESTART, rec_mpi_io, nrefs,          &
          indr, namr, (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), axx, ayy, azz, &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /),                            &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), ierr )

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

     Allocate (chbuf(1:mxatms),iwrk(1:mxatms), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_write allocation failure 1, node: ', idnode
        Call error(0)
     End If

! node 0 handles I/O

     If (idnode == 0) Then

! Obliterate old REFERENCE if not needed and print header

        Open(Unit=nrefdt, File=name, Form='formatted', Access='direct', Recl=73, Status='replace')
        Write(Unit=nrefdt, Fmt='(a72,a1)',         Rec=Int(1,ip)) cfgname(1:72),lf
        Write(Unit=nrefdt, Fmt='(3i10,a42,a1)',    Rec=Int(2,ip)) 0,imcon,megref,Repeat(' ',42),lf
        Write(Unit=nrefdt, Fmt='(3f20.10,a12,a1)', Rec=Int(3,ip)) cell(1),cell(2),cell(3),Repeat(' ',12),lf
        Write(Unit=nrefdt, Fmt='(3f20.10,a12,a1)', Rec=Int(4,ip)) cell(4),cell(5),cell(6),Repeat(' ',12),lf
        Write(Unit=nrefdt, Fmt='(3f20.10,a12,a1)', Rec=Int(5,ip)) cell(7),cell(8),cell(9),Repeat(' ',12),lf
        jj=5

! Get to real space

        Do i=1,nrefs
           chbuf(i)=namr(i)
           iwrk(i)=indr(i)

           axx(i)=cell(1)*xr(i)+cell(4)*yr(i)+cell(7)*zr(i)
           ayy(i)=cell(2)*xr(i)+cell(5)*yr(i)+cell(8)*zr(i)
           azz(i)=cell(3)*xr(i)+cell(6)*yr(i)+cell(9)*zr(i)
        End Do

        jatms=nrefs

        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Traject_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(chbuf,8*jatms,MPI_CHARACTER,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Traject_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Traject_tag,dlp_comm_world,status,ierr)

! Get to real space

              Do i=1,jatms
                 x=axx(i) ; y=ayy(i) ; z=azz(i)
                 axx(i)=cell(1)*x+cell(4)*y+cell(7)*z
                 ayy(i)=cell(2)*x+cell(5)*y+cell(8)*z
                 azz(i)=cell(3)*x+cell(6)*y+cell(9)*z
              End Do
           End If

! DD bound

           Call pbcshift(imcon,cell,jatms,axx,ayy,azz)

           Do i=1,jatms
              rec=Int(jj,ip)+Int(iwrk(i)-1,ip)*Int(2)+Int(1,ip)
              Write(Unit=nrefdt, Fmt='(a8,i10,a54,a1)', Rec=rec) chbuf(i),iwrk(i),Repeat(' ',54),lf

              rec=rec+Int(1,ip)
              Write(Unit=nrefdt, Fmt='(3g20.10,a12,a1)', Rec=rec) axx(i),ayy(i),azz(i),Repeat(' ',12),lf
           End Do
        End Do

! Close REFERENCE

        Close(Unit=nrefdt)

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Traject_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(nrefs,1,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(namr,8*nrefs,MPI_CHARACTER,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(indr,nrefs,MPI_INTEGER,0,Traject_tag,dlp_comm_world,ierr)

        Call MPI_SEND(xr,nrefs,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(yr,nrefs,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)
        Call MPI_SEND(zr,nrefs,wp_mpi,0,Traject_tag,dlp_comm_world,ierr)

     End If

     Deallocate (chbuf,iwrk, Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'defects_reference_write deallocation failure 1, node: ', idnode
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

  Deallocate (axx,ayy,azz, Stat=fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_write deallocation failure, node: ', idnode
     Call error(0)
  End If

  If (mxnode > 1) Call gsync()

End Subroutine defects_reference_write
