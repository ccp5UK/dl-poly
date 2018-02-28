Subroutine inversions_table_read(invr_name)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABINV file (for inversion potentials & forces only)
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module
  Use setup_module, Only : pi,delth_max,ntable,nrite, &
                           mxtinv,mxginv,zero_plus,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use inversions_module
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Character( Len = 32 ), Intent( In    ) :: invr_name(1:mxtinv)

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 32 )  :: idinvr
  Character( Len = 8   ) :: atom1,atom2,atom3,atom4

  Integer                :: fail(1:2),ngrid,rtinv,itinv,jtinv,katom1,katom2,katom3,katom4,jtpatm,i,l
  Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,rrr0, &
                            ppp,vk,vk1,vk2,t1,t2,bufp0,bufv0

  Integer,           Allocatable :: read_type(:)
  Real( Kind = wp ), Allocatable :: bufpot(:),bufvir(:)


  If (idnode == 0) Open(Unit=ntable, File='TABINV')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

! read mesh resolution not needed for inversion angle dependent
! potentials/forces as delpot=180/ngrid running from 0 to 180

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

  delpot = 180.0_wp/Real(ngrid,wp)

  dlrpot = 180.0_wp/Real(mxginv-4,wp)

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
             & ' TABINV file actual angular increment : ',1p,e15.7)") &
             delth_max, delpot
        Write(nrite,"(/,                                                &
             & ' expected (minimum) number of grid points : ',0p,i10,/, &
             & ' TABINV file actual number of grid points : ',0p,i10)") &
             mxginv-4, ngrid
     End If
     Call error(22)
  End If
  safe=.true.

  remake=.false.
  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     rdr=1.0_wp/delpot
     If (idnode == 0) Write(nrite,"(/,' TABINV arrays resized for mxgrid = ',i10)") mxginv-4
  End If

! compare grids dimensions

  If (ngrid < mxginv-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxginv-4,wp),0.0_wp)
     Call error(48)
  End If

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  fail=0
  Allocate (read_type(1:ltpinv(0)),          Stat=fail(1))
  Allocate (bufpot(0:ngrid),bufvir(0:ngrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - inversions_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  Call allocate_invr_pot_arrays()

  read_type=0 ! initialise read_type
  Do rtinv=1,ltpinv(0)
     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

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
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABINV'
        Call error(91)
     End If

! Construct unique name for the tabulated inversion

     If      (Min(katom2,katom3,katom4) == katom2) Then
        If (katom3 <= katom4) Then
           idinvr = atom1//atom2//atom3//atom4
        Else
           idinvr = atom1//atom2//atom4//atom3
        End If
     Else If (Min(katom2,katom3,katom4) == katom3) Then
        If (katom2 <= katom4) Then
           idinvr = atom1//atom3//atom2//atom4
        Else
           idinvr = atom1//atom3//atom4//atom2
        End If
     Else
        If (katom2 <= katom3) Then
           idinvr = atom1//atom4//atom2//atom3
        Else
           idinvr = atom1//atom4//atom3//atom2
        End If
     End If

! read potential arrays if potential is defined

     itinv=0
     Do jtinv=1,ltpinv(0)
        If (invr_name(jtinv) == idinvr) Then
           Do itinv=1,mxtinv
              If (ltpinv(itinv) == jtinv) Exit
           End Do
           Exit
        End If
     End Do

     If (itinv == 0) Then ! All(invr_name /= idinvr)
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABINV'
        Call error(89)
     End If
     If (Any(read_type == jtinv)) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'***',atom3,'***',atom4,'**** entry in TABINV'
        Call error(172)
     Else
        read_type(jtinv)=jtinv
     End If

! read in potential & force arrays

     Do i=0,2
        bufpot(i) = 0.0_wp
        bufvir(i) = 0.0_wp
     End Do

! read in the zero and/or first & second data elements (potential & virial)

     If (idnode == 0) Then
        rrr=0.0_wp
        Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

        If (rrr > zero_plus) Then ! no zero element data => extrapolate to zero
           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABINV stated  angular increment : ',1p,e15.7,/, &
                 & ' TABINV read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr
           End If

           bufpot(1) = bufp0
           bufvir(1) = bufv0
           rrr0      = rrr

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-rrr0-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABINV stated  angular increment : ',1p,e15.7,/, &
                 & ' TABINV read-in angular increment : ',1p,e15.7)") &
                 delpot,rrr-rrr0
           End If

           bufpot(2) = bufp0
           bufvir(2) = bufv0

! linear extrapolation for grid point 0 at distances close to 0

           bufpot(0) = 2.0_wp*bufpot(1)-bufpot(2)
           bufvir(0) = (2.0_wp*bufvir(1)-0.5_wp*bufvir(2))/dlrpot
        Else ! zero element data found => read in the first element for checking delr
           bufpot(0) = bufp0
           bufvir(0) = bufv0

           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufp0,bufv0

           If (Abs((rrr-delpot)/delpot) > 1.0e-8_wp) Then
              safe=.false.
              If (idnode == 0) Write(nrite,"(/,                       &
                 & ' TABINV stated  angular increment : ',1p,e15.7,/, &
                 & ' TABINV read-in angular increment : ',1p,e15.7)") &
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

! reconstruct arrays using 3pt interpolation

     If (remake) Then
        Do i=1,mxginv-4
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
           vinv(i,jtinv) = t1 + (t2-t1)*ppp*0.5_wp
           vinv(i,jtinv) = vinv(i,jtinv)*engunit ! convert to internal units

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
           ginv(i,jtinv) = t1 + (t2-t1)*ppp*0.5_wp
           ginv(i,jtinv) = ginv(i,jtinv)*engunit*rad2dgr ! convert to internal units
        End Do

        ginv(-1,jtinv) = rad2dgr/dlrpot
     Else
        Do i=1,mxginv-4
           vinv(i,jtinv) = bufpot(i)*engunit ! convert to internal units
           ginv(i,jtinv) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do

! linear extrapolation for the grid point just beyond the cutoff

        vinv(mxginv-3,jtinv) = 2.0_wp*vinv(mxginv-4,jtinv) - vinv(mxginv-5,jtinv)
        ginv(mxginv-3,jtinv) = 2.0_wp*ginv(mxginv-4,jtinv) - ginv(mxginv-5,jtinv)

        ginv(-1,jtinv) = rad2dgr/delpot
     End If

! grid point at 0 and linear extrapolation for the grid point at mxginv-2


     vinv(0,jtinv) = bufpot(0)
     ginv(0,jtinv) = bufvir(0)

     vinv(mxginv-2,jtinv) = 2.0_wp*vinv(mxginv-3,jtinv) - vinv(mxginv-4,jtinv)
     ginv(mxginv-2,jtinv) = 2.0_wp*ginv(mxginv-3,jtinv) - ginv(mxginv-4,jtinv)
  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABINV file'
  End If

! Break if not safe

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(22)

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - inversions_table_read deallocation failure, node: ', idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine inversions_table_read
