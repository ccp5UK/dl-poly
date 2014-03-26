Subroutine inversions_table_read(invr_name)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABINV file (for inversion potentials & forces only)
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : pi,ntable,nrite,mxtinv,mxgrid,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use inversions_module
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Character( Len = 32 ), Intent( In    ) :: invr_name(:)

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 32 )  :: idinvr
  Character( Len = 8   ) :: atom1,atom2,atom3,atom4

  Integer                :: fail(1:2),ngrid,rtinv,itinv,jtinv,katom1,katom2,katom3,katom4,jtpatm,i,l
  Real( Kind = wp )      :: delpot,dlrpot,rad2dgr,dgr2rad,rdr,rrr,ppp,vk,vk1,vk2,t1,t2

  Integer,                           Allocatable :: read_type(:)
  Real( Kind = wp ), Dimension( : ), Allocatable :: bufpot(:),bufvir(:)

  fail=0
  Allocate (read_type(1:ltpinv(0)),            Stat=fail(1))
  Allocate (bufpot(0:mxgrid),bufvir(0:mxgrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - inversions_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  Call allocate_invr_pot_arrays()

  remake=.false.

  If (idnode == 0) Open(Unit=ntable, File='TABINV')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABINV if it's in .xvg format

! read mesh resolution not needed for inversion angle dependent
! potentials/forces as delpot=180/ngrid running from 0 to 180

  rad2dgr= 180.0_wp/pi
  dgr2rad= pi/180.0_wp

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

! check array dimensions

  If (ngrid > mxgrid-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgrid-4,wp),0.0_wp)
     Call error(48)
  End If

  delpot = 180.0_wp/Real(ngrid,wp)
  dlrpot = 180.0_wp/Real(mxgrid-4,wp)

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > dlrpot .and. (.not.safe)) Then
     If (idnode == 0) Then
        Write(nrite,"(/,                                          &
           & 'expected (minimum) angular increment : ',1p,e15.7,/, &
           & 'TABINV file        angular increment : ',1p,e15.7)") &
           dlrpot, delpot
        Write(nrite,"(/,                                             &
           & 'expected (minimum) number of grid points : ',0p,i10,/, &
           & 'TABINV file actual number of grid points : ',0p,i10)") &
           mxgrid-4, ngrid
     End If
     Call error(22)
  End If
  safe=.true.

  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     If (idnode == 0) Write(nrite,"(/,' TABINV arrays resized for mxgrid = ',i10)") mxgrid-4
  End If

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

     Do i=1,ngrid
        If (idnode == 0) Then
           Read(Unit=ntable, Fmt=*, End=100, Err=100) rrr,bufpot(i),bufvir(i)
        Else
           bufpot(i) = 0.0_wp
           bufvir(i) = 0.0_wp
        End If
     End Do

! just in case, linear interpolation for for angle -180 (missing in the TABINV)

     bufpot(0) = 2.0_wp*bufpot(1)-bufpot(2)
     bufvir(0) = 2.0_wp*bufvir(1)-bufvir(2)

     If (mxnode > 1) Then
        Call MPI_BCAST(bufpot(0:ngrid), ngrid, wp_mpi, 0, dlp_comm_world, ierr)
        Call MPI_BCAST(bufvir(0:ngrid), ngrid, wp_mpi, 0, dlp_comm_world, ierr)
     Endif

! reconstruct arrays using 3pt interpolation

     If (remake) Then
        rdr=1.0_wp/delpot
        Do i=1,mxgrid-4
           rrr = Real(i,wp)*dlrpot
           l   = Int(rrr*rdr)
           ppp = rrr*rdr-Real(l,wp)

           vk  = bufpot(l)
           vk1 = bufpot(l+1)
           vk2 = bufpot(l+2)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           vinv(i,jtinv) = t1 + (t2-t1)*ppp*0.5_wp
           vinv(i,jtinv) = vinv(i,jtinv)*engunit ! convert to internal units

           vk  = bufvir(l)
           vk1 = bufvir(l+1)
           vk2 = bufvir(l+2)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           ginv(i,jtinv) = t1 + (t2-t1)*ppp*0.5_wp
           ginv(i,jtinv) = ginv(i,jtinv)*engunit*rad2dgr ! convert to internal units
        End Do
        ginv(0,jtinv) = rad2dgr/dlrpot
     Else
        Do i=1,ngrid
           vinv(i,jtinv) = bufpot(i)*engunit ! convert to internal units
           ginv(i,jtinv) = bufvir(i)*engunit*rad2dgr ! convert to internal units
        End Do
        ginv(0,jtinv) = rad2dgr/delpot
     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABINV file'
  End If

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
