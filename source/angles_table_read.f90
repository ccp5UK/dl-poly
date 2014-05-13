Subroutine angles_table_read(angl_name)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABANG file (for angle potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov may 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : pi,delth_max,ntable,nrite, &
                           mxtang,mxgang,zero_plus,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use angles_module
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Character( Len = 24 ), Intent( In    ) :: angl_name(1:mxtang)

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 24 )  :: idangl
  Character( Len = 8   ) :: atom1,atom2,atom3

  Integer                :: fail(1:2),ngrid,rtang,itang,jtang,katom1,katom2,katom3,jtpatm,i,l
  Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

  Integer,                           Allocatable :: read_type(:)
  Real( Kind = wp ), Dimension( : ), Allocatable :: bufpot(:),bufvir(:)


  If (idnode == 0) Open(Unit=ntable, File='TABANG')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

! read mesh resolution not needed for angle dependent
! potentials/forces as delpot=180/ngrid running from 0 to 180

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABANG if it's in .xvg format

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  delpot = 180.0_wp/Real(ngrid,wp)

  dlrpot = 180.0_wp/Real(mxgang-4,wp)

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
             & ' TABANG file actual angular increment : ',1p,e15.7)") &
             delth_max, delpot
        Write(nrite,"(/,                                                &
             & ' expected (minimum) number of grid points : ',0p,i10,/, &
             & ' TABANG file actual number of grid points : ',0p,i10)") &
             mxgang-4, ngrid
     End If

     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     If (idnode == 0) Write(nrite,"(/,' TABANG arrays resized for mxgrid = ',i10)") mxgang-4
  End If

! compare grids dimensions

  If (ngrid < mxgang-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgang-4,wp),0.0_wp)
     Call error(48)
  End If

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  fail=0
  Allocate (read_type(1:ltpang(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - angles_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  Call allocate_angl_pot_arrays()

  read_type=0 ! initialise read_type
  Do rtang=1,ltpang(0)
     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABANG if it's in .xvg format

     Call get_word(record,atom1)
     Call get_word(record,atom2)
     Call get_word(record,atom3)

     katom1=0
     katom2=0
     katom3=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
        If (atom3 == unqatm(jtpatm)) katom3=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0 .or. katom3 == 0) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call error(84)
     End If

! Construct unique name for the tabulated angle

     If (katom1 <= katom3) Then
        idangl = atom1//atom2//atom3
     Else
        idangl = atom3//atom2//atom1
     End If

! read potential arrays if potential is defined

     itang=0
     Do jtang=1,ltpang(0)
        If (angl_name(jtang) == idangl) Then
           Do itang=1,mxtang
              If (ltpang(itang) == jtang) Exit
           End Do
           Exit
        End If
     End Do

     If (itang == 0) Then ! All(angl_name /= idangl)
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call error(83)
     End If
     If (Any(read_type == jtang)) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'**** entry in TABANG'
        Call error(172)
     Else
        read_type(jtang)=jtang
     End If

! read in potential & force arrays

     Do i=0,2
        bufpot(0) = 0.0_wp
        bufvir(0) = 0.0_wp
     End Do

! read in the zero and/or first & second data elements (potential & virial)

     If (idnode == 0) Then
        Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

        If (rrr > zero_plus) Then ! no zero element data => extrapolate to zero
           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABANG stated  angular increment : ',1p,e15.7,/, &
                 & ' TABANG read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           bufpot(2) = bufp0
           bufvir(2) = bufv0

! linear extrapolation for distance close to 0

           bufpot(0) = 2.0_wp*bufpot(1)-bufpot(2)
           bufvir(0) = (2.0_wp*bufvir(1)-0.5_wp*bufvir(2))*rdr
        Else ! zero element data found => read in the first element for checking delr
           bufpot(0) = bufp0
           bufvir(0) = bufv0 ! virial/angle - finite force!

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABANG stated  angular increment : ',1p,e15.7,/, &
                 & ' TABANG read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           bufpot(2) = bufp0
           bufvir(2) = bufv0
        End If
     End If

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

! reconstruct arrays using 3pt interpolation

     If (remake) Then
        Do i=1,mxgang-3
           rrr = Real(i,wp)*dlrpot
           l   = Int(rrr*rdr)
           ppp = rrr*rdr-Real(l,wp)

           vk  = bufpot(l)

! linear extrapolation for the grid points just beyond the cutoff

           If (l+2 > ngrid) Then
              If (l+1 > ngrid) Then
                 vk1 = 2.0_wp*bufpot(l)-bufpot(l-1)
                 vk2 = 2.0_wp*vk1-bufpot(l)
              Else
                 vk1 = bufpot(l+1)
                 vk2 = 2.0_wp*bufpot(l+1)-bufpot(l)
              End If
           Else
              vk1 = bufpot(l+1)
              vk2 = bufpot(l+2)
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           vang(i,jtang) = t1 + (t2-t1)*ppp*0.5_wp
           vang(i,jtang) = vang(i,jtang)*engunit ! convert to internal units

           vk  = bufvir(l)

! linear extrapolation for the grid points just beyond the cutoff

           If (l+2 > ngrid) Then
              If (l+1 > ngrid) Then
                 vk1 = 2.0_wp*bufvir(l)-bufvir(l-1)
                 vk2 = 2.0_wp*vk1-bufvir(l)
              Else
                 vk1 = bufvir(l+1)
                 vk2 = 2.0_wp*bufvir(l+1)-bufvir(l)
              End If
           Else
              vk1 = bufvir(l+1)
              vk2 = bufvir(l+2)
           End If

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           gang(i,jtang) = t1 + (t2-t1)*ppp*0.5_wp
           gang(i,jtang) = gang(i,jtang)*engunit*rad2dgr ! convert to internal units
        End Do

        gang(-1,jtang) = rad2dgr/dlrpot
     Else
        Do i=1,mxgang-4
           vang(i,jtang) = bufpot(i)*engunit         ! convert to internal units
           gang(i,jtang) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do

! linear extrapolation for the grid point just beyond the cutoff

        vang(mxgang-3,jtang) = 2.0_wp*vang(mxgang-4,jtang) - vang(mxgang-5,jtang)
        gang(mxgang-3,jtang) = 2.0_wp*gang(mxgang-4,jtang) - gang(mxgang-5,jtang)

        gang(-1,jtang) = rad2dgr/delpot
     End If

! grid point at 0 and linear extrapolation for the grid point at mxgang-2

     vang(0,jtang) = bufpot(0)
     gang(0,jtang) = bufvir(0)

     vang(mxgang-2,jtang) = 2.0_wp*vang(mxgang-3,jtang) - vang(mxgang-4,jtang)
     gang(mxgang-2,jtang) = 2.0_wp*gang(mxgang-3,jtang) - gang(mxgang-4,jtang)
  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABANG file'
  End If

! Break if not safe

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - angles_table_read deallocation failure, node: ', idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine angles_table_read
