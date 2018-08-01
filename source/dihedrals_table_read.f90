Subroutine dihedrals_table_read(dihd_name)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABDIH file (for dihedral potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov january 2017
! amended   - a.v.brukhno & i.t.todorov december 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : pi,delth_max,ntable,nrite, &
                           mxtdih,mxgdih,zero_plus,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use dihedrals_module
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Character( Len = 32 ), Intent( In    ) :: dihd_name(1:mxtdih)

  Logical                :: safe,remake,zero
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 32 )  :: iddihd
  Character( Len = 8   ) :: atom1,atom2,atom3,atom4

  Integer                :: fail(1:2),ngrid,rtdih,itdih,jtdih,katom1,katom2,katom3,katom4,jtpatm,i,l
  Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,rrr0, &
                            ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

  Integer,           Allocatable :: read_type(:)
  Real( Kind = wp ), Allocatable :: bufpot(:),bufvir(:)


  If (idnode == 0) Open(Unit=ntable, File='TABDIH')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

! read mesh resolution not needed for dihedral angle dependent
! potentials/forces as delpot=360/ngrid running from -180 to 180

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABDIH if it's in .xvg format

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  delpot = 360.0_wp/Real(ngrid,wp)

  dlrpot = 360.0_wp/Real(mxgdih-4,wp)

! check grid spacing

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > delth_max .and. (.not.safe)) Then
     If (idnode == 0) Then
        Write(nrite,"(/,                                              &
             & ' expected (maximum) angular increment : ',1p,e15.7,/, &
             & ' TABDIH file actual angular increment : ',1p,e15.7)") &
             delth_max, delpot
        Write(nrite,"(/,                                                &
             & ' expected (minimum) number of grid points : ',0p,i10,/, &
             & ' TABDIH file actual number of grid points : ',0p,i10)") &
             mxgdih-4, ngrid
     End If
     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     If (idnode == 0) Write(nrite,"(/,' TABDIH arrays resized for mxgrid = ',i10)") mxgdih-4
  End If

! compare grids dimensions

  If (ngrid < mxgdih-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgdih-4,wp),0.0_wp)
     Call error(48)
  End If

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  fail=0
  Allocate (read_type(1:ltpdih(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - dihedrals_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  Call allocate_dihd_pot_arrays()

  read_type=0 ! initialise read_type
  Do rtdih=1,ltpdih(0)
     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABDIH if it's in .xvg format

     Call get_word(record,atom1)
     Call get_word(record,atom2)
     Call get_word(record,atom3)
     Call get_word(record,atom4)

     katom1=0
     katom2=0
     katom3=0
     katom4=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
        If (atom3 == unqatm(jtpatm)) katom3=jtpatm
        If (atom4 == unqatm(jtpatm)) katom4=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0 .or. katom3 == 0 .or. katom4 == 0) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABDIH'
        Call error(91)
     End If

! Construct unique name for the tabulated dihedral

     If (katom1 <= katom4) Then
        iddihd = atom1//atom2//atom3//atom4
     Else
        iddihd = atom4//atom3//atom2//atom1
     End If

! read potential arrays if potential is defined

     itdih=0
     Do jtdih=1,ltpdih(0)
        If (dihd_name(jtdih) == iddihd) Then
           Do itdih=1,mxtdih
              If (ltpdih(itdih) == jtdih) Exit
           End Do
           Exit
        End If
     End Do

     If (itdih == 0) Then ! All(dihd_name /= iddihd)
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABDIH'
        Call error(89)
     End If
     If (Any(read_type == jtdih)) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABDIH'
        Call error(172)
     Else
        read_type(jtdih)=jtdih
     End If

! read in potential & force arrays

     Do i=0,2
        bufpot(i) = 0.0_wp
        bufvir(i) = 0.0_wp
     End Do

! read in the zero and/or first & second data elements (potential & virial)

     zero=.false.
     If (idnode == 0) Then
        rrr=0.0_wp
        Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

        If (rrr+180.0_wp > zero_plus) Then ! no zero element data => extrapolate to zero
           If (Abs((rrr+180.0_wp-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABDIH stated  angular increment : ',1p,e15.7,/, &
                 & ' TABDIH read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0
           rrr0      = rrr

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABDIH stated  angular increment : ',1p,e15.7,/, &
                 & ' TABDIH read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr-rrr0
           End If

           bufpot(2) = bufp0
           bufvir(2) = bufv0
        Else ! zero element data found => read in the first element for checking delr
           zero=.true.
           bufpot(0) = bufp0
           bufvir(0) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABDIH stated  angular increment : ',1p,e15.7,/, &
                 & ' TABDIH read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           bufpot(2) = bufp0
           bufvir(2) = bufv0
        End If
     End If

! read in potential & force arrays

     Do i=3,ngrid
        If (idnode == 0) Then
           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufpot(i),bufvir(i)
        Else
           bufpot(i) = 0.0_wp
           bufvir(i) = 0.0_wp
        End If
     End Do

     If (mxnode > 1) Then
        Call MPI_BCAST(bufpot(0:ngrid), ngrid, wp_mpi, 0, dlp_comm_world, ierr)
        Call MPI_BCAST(bufvir(0:ngrid), ngrid, wp_mpi, 0, dlp_comm_world, ierr)
     End If

! linear extrapolation for grid point 0 at distances close to 0 -
! midpoint for a periodic function

     If (.not.zero) Then
        bufpot(0) = 0.5_wp*(bufpot(1)-bufpot(ngrid))
        bufvir(0) = 0.5_wp*(bufvir(1)-bufvir(ngrid))
     End If

! reconstruct arrays using 3pt interpolation

     If (remake) Then
        Do i=0,mxgdih-4
           rrr = Real(i,wp)*delth_max
           l   = Int(rrr*rdr)
           ppp = rrr*rdr-Real(l,wp)

! cyclic grid

           If (l  <= ngrid) Then
              vk  = bufpot(l)
           Else
              vk  = bufpot(Mod(l,ngrid+1))
           End If
           If (l+1 <= ngrid) Then
              vk1 = bufpot(l+1)
           Else
              vk1 = bufpot(Mod(l+1,ngrid+1))
           End If
           If (l+2 <= ngrid) Then
              vk2 = bufpot(l+2)
           Else
              vk2 = bufpot(Mod(l+2,ngrid+1))
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           vdih(i,jtdih) = t1 + (t2-t1)*ppp*0.5_wp
           vdih(i,jtdih) = vdih(i,jtdih)*engunit ! convert to internal units

! cyclic grid

           If (l  <= ngrid) Then
              vk  = bufvir(l)
           Else
              vk  = bufvir(Mod(l,ngrid+1))
           End If
           If (l+1 <= ngrid) Then
              vk1 = bufvir(l+1)
           Else
              vk1 = bufvir(Mod(l+1,ngrid+1))
           End If
           If (l+2 <= ngrid) Then
              vk2 = bufvir(l+2)
           Else
              vk2 = bufvir(Mod(l+2,ngrid+1))
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           gdih(i,jtdih) = t1 + (t2-t1)*ppp*0.5_wp
           gdih(i,jtdih) = gdih(i,jtdih)*engunit*rad2dgr ! convert to internal units
        End Do

        gdih(-1,jtdih) = rad2dgr/dlrpot
     Else
        Do i=0,mxgdih-4
           vdih(i,jtdih) = bufpot(i)*engunit ! convert to internal units
           gdih(i,jtdih) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do

! cyclic grid

        vdih(mxgdih-3,jtdih) = vdih(0,jtdih)
        vdih(mxgdih-2,jtdih) = vdih(1,jtdih)
        gdih(mxgdih-3,jtdih) = gdih(0,jtdih)
        gdih(mxgdih-2,jtdih) = gdih(1,jtdih)

        gdih(-1,jtdih) = rad2dgr/delpot
     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABDIH file'
  End If

! Break if not safe

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - dihedrals_table_read deallocation failure, node: ', idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine dihedrals_table_read
