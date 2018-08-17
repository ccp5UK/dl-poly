Subroutine system_expand(l_str,rcut,nx,ny,nz,megatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 utility to expand the MD system by a nx*ny*nz volumetric
! replication of its contents along the MD cell lattice vectors,
! creating a new matching pair of topology-interaction (FIELD) and
! crystallographic (CONFIG) files, preserving FIELD's template intact
!
! supported image conditions: 1,2,3, 6(nz==1)
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
! contrib   - w.smith, i.j.bush
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module,        Only : nconf,nfield,nmpldt,nrite,config,field, &
                                  mxatms
  Use site_module
  Use config_module,       Only : cfgname,imcon,cell,natms,lsi,lsa, &
                                  atmnam,xxx,yyy,zzz
  Use core_shell_module,   Only : numshl,lstshl
  Use constraints_module,  Only : numcon,lstcon
  Use rigid_bodies_module, Only : numrgd,lstrgd
  Use bonds_module,        Only : numbonds,lstbnd,keybnd
  Use angles_module,       Only : numang,lstang,keyang
  Use dihedrals_module,    Only : numdih,lstdih
  Use inversions_module,   Only : numinv,lstinv
  Use parse_module,        Only : tabs_2_blanks, get_word, strip_blanks, &
                                  lower_case, word_2_real
  Use io_module,           Only : io_set_parameters,        &
                                  io_get_parameters,        &
                                  io_init, io_nc_create,    &
                                  io_open, io_write_record, &
                                  io_nc_put_var,            &
                                  io_write_sorted_file,     &
                                  io_close, io_finalize,    &
                                  io_delete,                &
                                  io_close, io_finalize,    &
                                  IO_RESTART,               &
                                  IO_BASE_COMM_NOT_SET,     &
                                  IO_ALLOCATION_ERROR,      &
                                  IO_UNKNOWN_WRITE_OPTION,  &
                                  IO_UNKNOWN_WRITE_LEVEL,   &
                                  IO_WRITE_UNSORTED_MPIIO,  &
                                  IO_WRITE_UNSORTED_DIRECT, &
                                  IO_WRITE_UNSORTED_MASTER, &
                                  IO_WRITE_SORTED_MPIIO,    &
                                  IO_WRITE_SORTED_DIRECT,   &
                                  IO_WRITE_SORTED_NETCDF,   &
                                  IO_WRITE_SORTED_MASTER

  Implicit None

  Logical,           Intent( In    ) :: l_str
  Integer,           Intent( In    ) :: nx,ny,megatm
  Real( Kind = wp ), Intent( In    ) :: rcut
  Integer,           Intent( InOut ) :: nz

  Integer, Parameter     :: recsz = 73 ! default record size

  Logical                :: safex,safey,safez,safer,safel,safem,safe,safeg, &
                            lmpldt=.false.

  Character( Len = 200 ) :: record,record1
  Character( Len = 40  ) :: word,fcfg,ffld,fmpl
  Integer                :: fail(1:5),nall, i,j,ix,iy,iz,m, &
                            itmols,setspc,imols,          &
                            indatm,indatm1,nattot,mxiter, &
                            sapmpt,sapmtt,iatm,jatm,      &
                            ishls,nshels,icnst,nconst,    &
                            nrigid,irgd,lrgd,             &
                            ibond,nbonds,iang,nangle,     &
                            idih,ndihed,iinv,ninver,      &
                            idm,loc_ind,index,at_scaled
  Integer(Kind=ip)       :: offset,rec
  Real( Kind = wp )      :: fx,fy,fz, x,y,z, t,c1,c2,c3,c4, &
                            celprp(1:10), hwx,hwy,hwz, r,   &
                            x1(1:1),y1(1:1),z1(1:1),        &
                            cell_vecs(1:3,1:3), lengths(1:3), angles(1:3)

! Some parameters and variables needed by io_module interfaces

  Integer                           :: fh, io_write
  Character( Len = recsz )          :: record2, record3
  Character                         :: lf

  Integer( Kind = MPI_OFFSET_KIND ) :: top_skip

  Real( Kind = wp ), Dimension( : ),     Allocatable :: f1,f2,f3, &
                                                        f4,f5,f6, &
                                                        f7,f8,f9, &
                                                        x_scaled, &
                                                        y_scaled, &
                                                        z_scaled, &
                                                        xm,ym,zm

  Integer,           Dimension( :,:,: ), Allocatable :: i_xyz

  Integer,           Dimension( : ),     Allocatable :: ltg_scaled

  Character( Len = Len( atmnam ) ), Dimension( : ), Allocatable :: atmnam_scaled


  fail=0
  Allocate (f1(1:nx),f2(1:nx),f3(1:nx),                      Stat=fail(1))
  Allocate (f4(1:ny),f5(1:ny),f6(1:ny),                      Stat=fail(2))
  Allocate (f7(1:nz),f8(1:nz),f9(1:nz),                      Stat=fail(3))
  Allocate (i_xyz(1:nx,1:ny,1:nz),                           Stat=fail(4))
  Allocate (xm(1:10*mxatms),ym(1:10*mxatms),zm(1:10*mxatms), Stat=fail(5))

  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'system_expand allocation failure, node: ', idnode
     Call error(0)
  End If

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_line_feed         = lf       )

! Print elapsed time and option header

  Call gtime(t)
  If (idnode == 0) Then
     Write(nrite,'(/,/,/,1x, "time elapsed since job start: ", f12.3, " sec",/)') t
     Write(nrite,'(4(1x,a,/))')                                                     &
     "*** Expanding the MD system by a nx*ny*nz volumetric replication        ***", &
     "*** of its contents along the MD cell lattice vectors, creating         ***", &
     "*** a new matching pair of topology-interaction (FIELD) and             ***", &
     "*** crystallographic (CONFIG) files, preserving FIELD's template intact ***"
     Write(nrite,'(1x,a,3i5,/)') '*** Replication dimensions (nx,ny,nz):', nx,ny,nz
  End If

! Holt or change execution if imcon is unsupported

  If (imcon == 0) Call error(570)
  If (imcon == 6 .and. nz > 1) Then
     nz=1
     Call warning(350,0.0_wp,0.0_wp,0.0_wp)
     Write(nrite,'(1x,a,3i5,/)') '*** Replication dimensions (nx,ny,nz):', nx,ny,nz
  End If

! Create names for the expanded CONFIG and FIELD

  record= ' ' ; Write(record,'(3(a1,i0))') '_',nx,'_',ny,'_',nz
  fcfg=' '
  fcfg=Trim(config) // record(1:Len_Trim(record))
  ffld=' '
  ffld=Trim(field) // record(1:Len_Trim(record))
  fmpl=' '
  fmpl="MPOLES" // record(1:Len_Trim(record))

! netCDF CONFIG name convention

  If (io_write == IO_WRITE_SORTED_NETCDF) fcfg=fcfg(1:Len_Trim(fcfg)) // '.nc'

  fx=Real(nx,wp)
  fy=Real(ny,wp)
  fz=Real(nz,wp)

  nall=nx*ny*nz

! Define cell vector displacement in z direction

  Do iz=1,nz
     z=Real(2*iz-nz-1,wp)
     f7(iz)=cell(7)*z/2.0_wp
     f8(iz)=cell(8)*z/2.0_wp
     f9(iz)=cell(9)*z/2.0_wp
  End Do

! Define cell vector displacement in y direction

  Do iy=1,ny
     y=Real(2*iy-ny-1,wp)
     f4(iy)=cell(4)*y/2.0_wp
     f5(iy)=cell(5)*y/2.0_wp
     f6(iy)=cell(6)*y/2.0_wp
  End Do

! Define cell vector displacement in x direction

  Do ix=1,nx
     x=Real(2*ix-nx-1,wp)
     f1(ix)=cell(1)*x/2.0_wp
     f2(ix)=cell(2)*x/2.0_wp
     f3(ix)=cell(3)*x/2.0_wp
  End Do

! Define hypercube counter

  Do iz=1,nz
     Do iy=1,ny
        Do ix=1,nx
           i_xyz(ix,iy,iz)=(ix-1)+nx*((iy-1)+ny*(iz-1))
        End Do
     End Do
  End Do

  Call dcell(cell,celprp) ! get cell properties

! define half cell widths and bond-length limit

  hwx = celprp(7)/2.0_wp
  hwy = celprp(8)/2.0_wp
  hwz = celprp(9)/2.0_wp
  c1  = Min(rcut/2.0_wp , 1.75_wp)
  c2  = c1*4.0_wp/3.0_wp
  c3  = c2*4.0_wp/3.0_wp
  c4  = c3*4.0_wp/3.0_wp

  If (idnode == 0) Then

! Make sure CONFIG(new) is empty and open it

     If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
         io_write == IO_WRITE_UNSORTED_DIRECT .or. &
         io_write == IO_WRITE_SORTED_MPIIO    .or. &
         io_write == IO_WRITE_SORTED_DIRECT   .or. &
         io_write == IO_WRITE_SORTED_NETCDF) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_delete( fcfg(1:Len_Trim(fcfg) ) )
        If (io_write == IO_WRITE_SORTED_NETCDF) Call io_nc_create( MPI_COMM_SELF, fcfg(1:Len_Trim(fcfg)), cfgname, megatm*nall )
        Call io_open( io_write, MPI_COMM_SELF, fcfg(1:Len_Trim(fcfg)), MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

     Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
              io_write == IO_WRITE_SORTED_MASTER ) Then

        Open(Unit=nconf, File=fcfg(1:Len_Trim(fcfg)), Status='replace')
        Close(Unit=nconf)
        Open(Unit=nconf, File=fcfg(1:Len_Trim(fcfg)), Form='formatted', Access='direct', Recl=recsz)
     End If

! Write configuration file headers

     Write(nrite,'(1x,2a)') '*** Expanding CONFIG in file ',fcfg(1:Len_Trim(fcfg))
     Write(nrite,'(1x,a)') '***'
     ! check if we expand a cube since not all time we end up with a cube back
     If ((imcon == 1) .and. ((nx /= ny) .or. (nx /= nz) .or. (ny /= nz ))) Then
        imcon = 3
     End If

     If      (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
              io_write == IO_WRITE_UNSORTED_DIRECT .or. &
              io_write == IO_WRITE_SORTED_MPIIO    .or. &
              io_write == IO_WRITE_SORTED_DIRECT) Then

        Write(record2, Fmt='(a72,a1)') cfgname(1:72),lf
        Call io_write_record( fh, Int(0,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3i10,a42,a1)') 0,imcon,nall*megatm,Repeat(' ',42),lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3f20.10,a12,a1)') fx*cell(1),fx*cell(2),fx*cell(3),Repeat(' ',12),lf
        Call io_write_record( fh, Int(2,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3f20.10,a12,a1)') fy*cell(4),fy*cell(5),fy*cell(6),Repeat(' ',12),lf
        Call io_write_record( fh, Int(3,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3f20.10,a12,a1)') fz*cell(7),fz*cell(8),fz*cell(9),Repeat(' ',12),lf
        Call io_write_record( fh, Int(4,MPI_OFFSET_KIND), record2 )

     Else If (io_write == IO_WRITE_SORTED_NETCDF) Then

        i=1 ! For config there is only one frame

        Call io_nc_put_var( 'time'           , fh, 0.0_wp, i, 1 )
        Call io_nc_put_var( 'step'           , fh,      0, i, 1 )
        Call io_nc_put_var( 'datalevel'      , fh,      0, i, 1 )
        Call io_nc_put_var( 'imageconvention', fh,  imcon, i, 1 )

        cell_vecs( :, 1 ) = fx * cell( 1:3 )
        cell_vecs( :, 2 ) = fy * cell( 4:6 )
        cell_vecs( :, 3 ) = fz * cell( 7:9 )

        lengths( 1 ) = fx * celprp( 1 )
        lengths( 2 ) = fy * celprp( 2 )
        lengths( 3 ) = fz * celprp( 3 )

        angles ( 1 ) = Acos( celprp( 5 ) )
        angles ( 2 ) = Acos( celprp( 6 ) )
        angles ( 3 ) = Acos( celprp( 4 ) )
        angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

        Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
        Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, i /), (/    3, 1 /) )
        Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, i /), (/    3, 1 /) )

     Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
              io_write == IO_WRITE_SORTED_MASTER ) Then

        Write(Unit=nconf, Fmt='(a72,a1)',         Rec=Int(1,ip)) cfgname(1:72),lf
        Write(Unit=nconf, Fmt='(3i10,a42,a1)',    Rec=Int(2,ip)) 0,imcon,nall*megatm,Repeat(' ',42),lf
        Write(Unit=nconf, Fmt='(3f20.12,a12,a1)', Rec=Int(3,ip)) fx*cell(1),fx*cell(2),fx*cell(3),Repeat(' ',12),lf
        Write(Unit=nconf, Fmt='(3f20.12,a12,a1)', Rec=Int(4,ip)) fy*cell(4),fy*cell(5),fy*cell(6),Repeat(' ',12),lf
        Write(Unit=nconf, Fmt='(3f20.12,a12,a1)', Rec=Int(5,ip)) fz*cell(7),fz*cell(8),fz*cell(9),Repeat(' ',12),lf

     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
      io_write == IO_WRITE_SORTED_MPIIO   .or. &
      io_write == IO_WRITE_SORTED_NETCDF) Then

     If (idnode == 0) Then
        Call io_close( fh )
        Call io_finalize
     End If

     Allocate ( atmnam_scaled( 1:natms * nall ), ltg_scaled( 1:natms * nall ), Stat = fail(1) )
     Allocate ( x_scaled( 1:natms * nall ), y_scaled( 1:natms * nall ), z_scaled( 1:natms * nall ), &
               Stat = fail(2) )
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'system_expand allocation failure 0, node: ', idnode
        Call error(0)
     End If

  End If

! Line counter in CONFIG(new):levcfg=0

  offset=Int(5,ip)

! Global atom counter

  nattot=0

! Local atom counter

  indatm=1

! local counter in scaled system

  at_scaled = 0

! running site and topology related intra-indices

  nshels=0
  nconst=0
  nrigid=0
  nbonds=0
  nangle=0
  ndihed=0
  ninver=0

  If (idnode == 0) Write(nrite,Fmt='(/,1x,a)') 'Checking topological contiguity of molecules...'

  safeg=.true. ! topology presumed safe

  sapmpt=0
  Do itmols=1,ntpmls
     setspc=nummols(itmols)*numsit(itmols)

     sapmtt=0
     Do imols=1,nummols(itmols)
        If (numsit(itmols) > 10*mxatms) Call error(0)

! Grab the coordinates of the atoms constituting this molecule

        indatm1=indatm
        Do m=1,numsit(itmols)
           nattot=nattot+1 ! Increase global atom counter in CONFIG(old)

           If (lsa(indatm1) == nattot) Then  ! If a local atom has a global index nattot
              loc_ind=lsi(indatm1)
              xm(m)=xxx(loc_ind)
              ym(m)=yyy(loc_ind)
              zm(m)=zzz(loc_ind)
              indatm1=indatm1+1 ! Increase local atom counter
           Else
              xm(m)=0.0_wp
              ym(m)=0.0_wp
              zm(m)=0.0_wp
           End If
        End Do
        nattot=nattot-numsit(itmols)

        If (mxnode > 1) Then
           Call gsum(xm(1:numsit(itmols)))
           Call gsum(ym(1:numsit(itmols)))
           Call gsum(zm(1:numsit(itmols)))
        End If

! Start unwrapping - not safe at start for each molecule

        indatm1=nattot-sapmpt-sapmtt
        safe=.false. ; mxiter=0
        Do While ((.not.safe) .and. mxiter < 42) ! meaning of LUEE is the limit
           If (.not.safe) mxiter=mxiter+1

           If ((mxiter == 42 .and. (.not.safe)) .and. (l_str .and. idnode == 0)) &
              Write(nrite,Fmt='(/,1x,2(a,i10),/)') 'MOLECULAR TYPE #: ',itmols, ' MOLECULE #: ',imols

           safe=.true.

           safel=.true.
           Do ishls=1,numshl(itmols)
              nshels=nshels+1

              iatm=lstshl(1,nshels)-indatm1
              jatm=lstshl(2,nshels)-indatm1

              safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
              safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
              safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
              safer=(safex .and. safey .and. safez)
              If (.not.safer) Then
                 x1(1)=xm(jatm)-xm(iatm)
                 y1(1)=ym(jatm)-ym(iatm)
                 z1(1)=zm(jatm)-zm(iatm)
                 Call images(imcon,cell,1,x1,y1,z1)
                 xm(jatm)=x1(1)+xm(iatm)
                 ym(jatm)=y1(1)+ym(iatm)
                 zm(jatm)=z1(1)+zm(iatm)
                 safer=.true.
              End If
              x=Abs(xm(jatm)-xm(iatm))
              y=Abs(ym(jatm)-ym(iatm))
              z=Abs(zm(jatm)-zm(iatm))
              t=c1
              safex=(x < t)
              safey=(y < t)
              safez=(z < t)
              r=Sqrt(x**2+y**2+z**2)
              If (safex .and. safey .and. safez) Then
                 safer=(r < t)
              Else
                 safer=.false.
              End If
              If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', r, ' > ', t, ' Angstroms'

                 t=c2
                 If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'

                 Write(nrite,Fmt='(1x,a,3i10)') 'CORE_SHELL UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',ishls,itmols,imols
                 Write(nrite,Fmt='(1x,a,3(1x,l1))') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z',safex,safey,safez
                 Write(nrite,Fmt='(1x,a,i10,3f10.1)')   'CORE  ',nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                 Write(nrite,Fmt='(1x,a,i10,3f10.1,/)') 'SHELL ',nattot+jatm,xm(jatm),ym(jatm),zm(jatm)
              End If
              safel=(safel .and. safer)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do icnst=1,numcon(itmols)
              nconst=nconst+1

              iatm=lstcon(1,nconst)-indatm1
              jatm=lstcon(2,nconst)-indatm1

              safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
              safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
              safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
              safer=(safex .and. safey .and. safez)
              If (.not.safer) Then
                 x1(1)=xm(jatm)-xm(iatm)
                 y1(1)=ym(jatm)-ym(iatm)
                 z1(1)=zm(jatm)-zm(iatm)
                 Call images(imcon,cell,1,x1,y1,z1)
                 xm(jatm)=x1(1)+xm(iatm)
                 ym(jatm)=y1(1)+ym(iatm)
                 zm(jatm)=z1(1)+zm(iatm)
                 safer=.true.
              End If
              x=Abs(xm(jatm)-xm(iatm))
              y=Abs(ym(jatm)-ym(iatm))
              z=Abs(zm(jatm)-zm(iatm))
              t=c2
              safex=(x < t)
              safey=(y < t)
              safez=(z < t)
              r=Sqrt(x**2+y**2+z**2)
              If (safex .and. safey .and. safez) Then
                 safer=(r < t)
              Else
                 safer=.false.
              End If
              If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', r, ' > ', t, ' Angstroms'

                 t=c3
                 If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'

                 Write(nrite,Fmt='(1x,a,3i10)') 'CONSTRAINT UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',icnst,itmols,imols
                 Write(nrite,Fmt='(1x,a,3(1x,l1))') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z',safex,safey,safez
                 Write(nrite,Fmt='(2i10,3f10.1)')   1,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                 Write(nrite,Fmt='(2i10,3f10.1,/)') 2,nattot+jatm,xm(jatm),ym(jatm),zm(jatm)
              End If
              safel=(safel .and. safer)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do irgd=1,numrgd(itmols)
              nrigid=nrigid+1

              safem=.true.
              lrgd=lstrgd(0,nrigid)
              Do i=1,lrgd-1
                 iatm=lstrgd(i,nrigid)-indatm1
                 Do j=i+1,lrgd
                    jatm=lstrgd(j,nrigid)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                          xm(jatm)=x1(1)+xm(iatm)
                          ym(jatm)=y1(1)+ym(iatm)
                          zm(jatm)=z1(1)+zm(iatm)
                          safer=.true.
                    End If
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    t=rcut
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) &
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') &
                            '*** WARNING **** WARNING *** DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'
                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'RIGID BODY UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',irgd,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,lrgd
                    iatm=lstrgd(i,nrigid)-indatm1
                    If (i < lrgd) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do ibond=1,numbonds(itmols)
              nbonds=nbonds+1

              iatm=lstbnd(1,nbonds)-indatm1
              jatm=lstbnd(2,nbonds)-indatm1

              safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
              safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
              safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
              safer=(safex .and. safey .and. safez)
              If (.not.safer) Then
                 x1(1)=xm(jatm)-xm(iatm)
                 y1(1)=ym(jatm)-ym(iatm)
                 z1(1)=zm(jatm)-zm(iatm)
                 Call images(imcon,cell,1,x1,y1,z1)
                 xm(jatm)=x1(1)+xm(iatm)
                 ym(jatm)=y1(1)+ym(iatm)
                 zm(jatm)=z1(1)+zm(iatm)
                 safer=.true.
              End If
              x=Abs(xm(jatm)-xm(iatm))
              y=Abs(ym(jatm)-ym(iatm))
              z=Abs(zm(jatm)-zm(iatm))
              If (keybnd(nbonds) > 0) Then
                 t=c2
              Else
                 t=3.0_wp*c1
              End If
              safex=(x < t)
              safey=(y < t)
              safez=(z < t)
              r=Sqrt(x**2+y**2+z**2)
              If (safex .and. safey .and. safez) Then
                 safer=(r < t)
              Else
                 safer=.false.
              End If
              If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', r, ' > ', t, ' Angstroms'

                 If (keybnd(nbonds) > 0) Then
                    t=c3
                 Else
                    t=3.0_wp*c2
                 End If
                 If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'

                 Write(nrite,Fmt='(1x,a,3i10)') 'BOND UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',ibond,itmols,imols
                 Write(nrite,Fmt='(1x,a,3(1x,l1))') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z',safex,safey,safez
                 Write(nrite,Fmt='(2i10,3f10.1)')   1,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                 Write(nrite,Fmt='(2i10,3f10.1,/)') 2,nattot+jatm,xm(jatm),ym(jatm),zm(jatm)
              End If
              safel=(safel .and. safer)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do iang=1,numang(itmols)
              nangle=nangle+1

              safem=.true.
              Do i=1,2
                 iatm=lstang(i,nangle)-indatm1
                 Do j=i+1,3
                    jatm=lstang(j,nangle)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                       xm(jatm)=x1(1)+xm(iatm)
                       ym(jatm)=y1(1)+ym(iatm)
                       zm(jatm)=z1(1)+zm(iatm)
                       safer=.true.
                    End If
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    If (keyang(nangle) > 0) Then
                       t=c1*Real(j-i+1,wp)
                    Else
                       t=c3*Real(j-i+1,wp)
                    End If
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) Then
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'

                       If (keyang(nangle) > 0) Then
                          t=c2*Real(j-i+1,wp)
                       Else
                          t=c4*Real(j-i+1,wp)
                       End If
                       If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'
                    End If
                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'ANGLE UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',iang,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,3
                    iatm=lstang(i,nangle)-indatm1
                    If (i < 3) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do idih=1,numdih(itmols)
              ndihed=ndihed+1

              safem=.true.
              Do i=1,3
                 iatm=lstdih(i,ndihed)-indatm1
                 Do j=i+1,4
                    jatm=lstdih(j,ndihed)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                       xm(jatm)=x1(1)+xm(iatm)
                       ym(jatm)=y1(1)+ym(iatm)
                       zm(jatm)=z1(1)+zm(iatm)
                       safer=.true.
                    End If
                    t=(c1+c2)*Real(j-i+1,wp)/2.0_wp
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) Then
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'

                       t=(c2+c3)*Real(j-i+1,wp)/2.0_wp
                       If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'
                    End If
                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'DIHEDRAL UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',idih,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,4
                    iatm=lstdih(i,ndihed)-indatm1
                    If (i < 4) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do iinv=1,numinv(itmols)
              ninver=ninver+1

              safem=.true.
              Do i=1,3
                 iatm=lstinv(i,ninver)-indatm1
                 Do j=i+1,4
                    jatm=lstinv(j,ninver)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                       xm(jatm)=x1(1)+xm(iatm)
                       ym(jatm)=y1(1)+ym(iatm)
                       zm(jatm)=z1(1)+zm(iatm)
                       safer=.true.
                    End If
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    t=c2
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. idnode == 0)) Then
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'

                       t=c3
                       If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'
                    End If

                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'INVERSION UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',iinv,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,4
                    iatm=lstinv(i,ninver)-indatm1
                    If (i < 4) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           If ( ((.not.safe) .and. imols <= nummols(itmols) .and. mxiter < 42) .or. &
                imols < nummols(itmols) ) Then
              nshels=nshels-numshl(itmols)
              nconst=nconst-numcon(itmols)
              nrigid=nrigid-numrgd(itmols)
              nbonds=nbonds-numbonds(itmols)
              nangle=nangle-numang(itmols)
              ndihed=ndihed-numdih(itmols)
              ninver=ninver-numinv(itmols)
           End If
        End Do

        safeg=(safeg.and.safe)

        Do m=1,numsit(itmols)
           nattot=nattot+1 ! Increase global atom counter in CONFIG(old)

           If (lsa(indatm) == nattot) Then ! If a local atom has a global index nattot

! Determine sending node for UN/SORTED MASTER

              If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                  io_write == IO_WRITE_SORTED_MASTER) Then
                 idm=idnode
                 If (mxnode > 1) Call gsum(idm)
              End If

! Get the local index of the particle

              loc_ind=lsi(indatm)

! Do particle replication by vector displacements in cyclic (z,y,x) directions

              Do iz=1,nz
                 Do iy=1,ny
                    Do ix=1,nx

                       x=xm(m)+f1(ix)+f4(iy)+f7(iz)
                       y=ym(m)+f2(ix)+f5(iy)+f8(iz)
                       z=zm(m)+f3(ix)+f6(iy)+f9(iz)

! Write 2 records @ line 'rec+1' and 'rec+2' for particle 'index' in CONFIG(new)

                       index = i_xyz(ix,iy,iz)*setspc + m
                       rec   = offset + Int(2,ip)*Int(index,ip) - Int(2,ip)
                       index = index + Int((offset - Int(5,ip))/Int(2,ip))

                       If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
                           io_write == IO_WRITE_SORTED_MPIIO   .or. &
                           io_write == IO_WRITE_SORTED_NETCDF) Then

                          at_scaled = at_scaled + 1

                          atmnam_scaled( at_scaled ) = atmnam( loc_ind )

                          ltg_scaled( at_scaled ) = index

                          x_scaled( at_scaled ) = x
                          y_scaled( at_scaled ) = y
                          z_scaled( at_scaled ) = z

                       Else

                          Write(record2, Fmt='(a8,i10,a54,a1)') atmnam(loc_ind),index,Repeat(' ',54),lf
                          Write(record3, Fmt='(3g20.12,a12,a1)') x,y,z,Repeat(' ',12),lf

                          If (io_write == IO_WRITE_UNSORTED_DIRECT .or. &
                              io_write == IO_WRITE_SORTED_DIRECT) Then

                             Call io_write_record( fh, Int(rec,MPI_OFFSET_KIND), record2 )
                             rec=rec+Int(1,ip)
                             Call io_write_record( fh, Int(rec,MPI_OFFSET_KIND), record3 )

                          Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                                   io_write == IO_WRITE_SORTED_MASTER) Then

                             If (idnode == 0) Then
                                rec=rec+Int(1,ip)
                                Write(Unit=nconf, Fmt='(a73)', Rec=rec) record2
                                rec=rec+Int(1,ip)
                                Write(Unit=nconf, Fmt='(a73)', Rec=rec) record3
                             Else
                                Call MPI_SEND(record2,recsz,MPI_CHARACTER,0,SysExpand_tag,dlp_comm_world,ierr)
                                Call MPI_SEND(record3,recsz,MPI_CHARACTER,0,SysExpand_tag,dlp_comm_world,ierr)
                             End If

                          End If

                       End If

                    End Do
                 End Do
              End Do

! Increase local atom counter

              indatm=indatm+1

           Else

! Determine sending node for UN/SORTED MASTER

              If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                  io_write == IO_WRITE_SORTED_MASTER) Then
                 idm=0 ! Initialise node number
                 If (mxnode > 1) Call gsum(idm)

                 Do iz=1,nz
                    Do iy=1,ny
                       Do ix=1,nx
                          rec   = offset + Int(2,ip)*(Int(i_xyz(ix,iy,iz),ip)*Int(setspc,ip) + Int(m,ip)) - Int(2,ip)

                          If (idnode == 0) Then
                             Call MPI_RECV(record2,recsz,MPI_CHARACTER,idm,SysExpand_tag,dlp_comm_world,status,ierr)
                             Call MPI_RECV(record3,recsz,MPI_CHARACTER,idm,SysExpand_tag,dlp_comm_world,status,ierr)

                             rec=rec+Int(1,ip)
                             Write(Unit=nconf, Fmt='(a73)', Rec=rec) record2
                             rec=rec+Int(1,ip)
                             Write(Unit=nconf, Fmt='(a73)', Rec=rec) record3
                          End If
                       End Do
                    End Do
                 End Do
              End If

           End If
        End Do

        sapmtt = sapmtt + numsit(itmols)
        offset = offset + Int(2,ip)*Int(numsit(itmols),ip)
     End Do

     sapmpt = sapmpt + sapmtt
     offset = offset + Int(2,ip)*Int(nall-1,ip)*Int(setspc,ip)
  End Do

  If ((.not.safeg) .and. idnode == 0) Write(nrite,Fmt='(/,1x,a)') &
     '*** warning - possible topological contiguity failures occurred !!! ***'

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_SORTED_MPIIO   .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

     Call io_set_parameters( user_comm = dlp_comm_world )
     Call io_init( recsz )
     Call io_open( io_write, dlp_comm_world, fcfg(1:Len_Trim(fcfg)), MPI_MODE_WRONLY, fh )

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        top_skip = Int(5,MPI_OFFSET_KIND)
     Else
        top_skip = Int(1,MPI_OFFSET_KIND) ! netCDF frame
     End If

     Call io_write_sorted_file( fh, 0, IO_RESTART, top_skip, at_scaled, &
          ltg_scaled, atmnam_scaled,                                    &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /),                     &
          x_scaled, y_scaled, z_scaled,                                 &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /),                     &
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

     Deallocate ( atmnam_scaled, ltg_scaled,    Stat = fail(1) )
     Deallocate ( x_scaled, y_scaled, z_scaled, Stat = fail(2) )
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'system_expand deallocation failure 0, node: ', idnode
        Call error(0)
     End If

  Else If (io_write == IO_WRITE_UNSORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_DIRECT) Then

     Call io_close( fh )
     Call io_finalize

  Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
           io_write == IO_WRITE_SORTED_MASTER) Then

     Close(Unit=nconf)

  End If

  Call gtime(t)

! Write summary data and proceed with FIELD

  x=0.5_wp*Min(fx*celprp(7),fy*celprp(8),fz*celprp(9))

  If (idnode == 0) Then
     Write(nrite,'(/,1x,3a)') '*** ', fcfg(1:Len_Trim(fcfg)), ' expansion completed !'
     Write(nrite,'(1x,a,i10,a)') '*** Size: ', nall*megatm, ' particles'
     Write(nrite,'(1x,a,f10.2,a)') '*** Maximum radius of cutoff: ', x, ' Angstroms'
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') t
     Open(Unit=nfield, File='FIELD', Status='old')
     Open(Unit=nconf, File=ffld(1:Len_Trim(ffld)), Status='replace')
     Write(nrite,'(/,1x,2a)')'*** Expanding FIELD in file ', ffld(1:Len_Trim(ffld))
     Write(nrite,'(1x,a)') '***'

! omit first line

     record=' '
     Read(Unit=nfield, Fmt='(a)', End=10) record
     Call tabs_2_blanks(record) ; Call strip_blanks(record)
     Write(nconf,'(a)') record(1:Len_Trim(record))

! read and process directives from field file

     Do
        record=' '
        Read(Unit=nfield, Fmt='(a)', End=10) record
        Call tabs_2_blanks(record) ; Call strip_blanks(record)
        record1=record
        Call get_word(record,word)
        Call lower_case(word)

        If      (word(1:5) == 'multi') Then

! MPOLES should exist

           lmpldt=.true.

! number of molecules of this type

        Else If (word(1:6) == 'nummol') Then

           Call get_word(record,word)
           index=Nint(word_2_real(word))
           Call get_word(record1,word)
           Write(Unit=nconf,Fmt='(a,i10)') word(1:Len_Trim(word)), nall*index

! close force field file

        Else If (word(1:5) == 'close') Then

           Call gtime(t)

           Write(Unit=nconf,Fmt='(a)') record1(1:Len_Trim(record1))
           Close(Unit=nfield)
           Close(Unit=nconf)
           Write(nrite,'(1x,3a)') '*** ', ffld(1:Len_Trim(ffld)), ' expansion done !'
           Exit

! just paste the copy

        Else

           Write(nconf,'(a)') record1(1:Len_Trim(record1))

        End If
     End Do

10   Continue

     If (lmpldt) Inquire(File='MPOLES', Exist=lmpldt)
     If (lmpldt) Then
        Open(Unit=nmpldt, File='MPOLES', Status='old')
        Open(Unit=nconf, File=fmpl(1:Len_Trim(fmpl)), Status='replace')
        Write(nrite,'(/,1x,2a)')'*** Expanding MPOLES in file ', fmpl(1:Len_Trim(fmpl))
        Write(nrite,'(1x,a)') '***'

! omit first line

        record=' '
        Read(Unit=nmpldt, Fmt='(a)', End=10) record
        Call tabs_2_blanks(record) ; Call strip_blanks(record)
        Write(nconf,'(a)') record(1:Len_Trim(record))

! read and process directives from mpoles file

        Do
           record=' '
           Read(Unit=nmpldt, Fmt='(a)', End=20) record
           Call tabs_2_blanks(record) ; Call strip_blanks(record)
           record1=record
           Call get_word(record,word)
           Call lower_case(word)

           If      (word(1:6) == 'nummol') Then

              Call get_word(record,word)
              index=Nint(word_2_real(word))
              Call get_word(record1,word)
              Write(Unit=nconf,Fmt='(a,i10)') word(1:Len_Trim(word)), nall*index

! close mpoles file

           Else If (word(1:5) == 'close') Then

              Call gtime(t)

              Write(Unit=nconf,Fmt='(a)') record1(1:Len_Trim(record1))
              Close(Unit=nmpldt)
              Close(Unit=nconf)
              Write(nrite,'(1x,3a)') '*** ', fmpl(1:Len_Trim(fmpl)), ' expansion done !'
              Exit

! just paste the copy

           Else

              Write(nconf,'(a)') record1(1:Len_Trim(record1))

           End If
        End Do
     End If

20   Continue

     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec",/)') t
     Write(nrite,'(1x,a)') '*** Simulation continues as scheduled...'
  End If
  If (mxnode > 1) Call gsync()

  Deallocate (f1,f2,f3, Stat=fail(1))
  Deallocate (f4,f5,f6, Stat=fail(2))
  Deallocate (f7,f8,f9, Stat=fail(3))
  Deallocate (i_xyz,    Stat=fail(4))
  Deallocate (xm,ym,zm, Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'system_expand dellocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine system_expand
