Subroutine bonds_table_read(bond_name,rcut)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABBND file (for bond potentials & forces only)
!
! copyright - daresbury laboratory
! author    - a.v.brukhno & i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module, Only : ntable,nrite,mxtbnd,mxgrid,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use bonds_module
  Use parse_module, Only : get_line,get_word,word_2_real

  Implicit None

  Character( Len = 16 ), Intent( In    ) :: bond_name(:)
  Real( Kind = wp ),     Intent( In    ) :: rcut

  Logical                :: safe,remake
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 16 )  :: idbond
  Character( Len = 8   ) :: atom1,atom2

  Integer                :: fail(1:2),ngrid,rtbnd,itbnd,jtbnd,katom1,katom2,jtpatm,i,l
  Real( Kind = wp )      :: cutpot,delpot,dlrpot,rdr,rrr,ppp,vk,vk1,vk2,t1,t2

  Integer,                           Allocatable :: read_type(:)
  Real( Kind = wp ), Dimension( : ), Allocatable :: bufpot(:),bufvir(:)

  fail=0
  Allocate (read_type(1:ltpbnd(0)),            Stat=fail(1))
  Allocate (bufpot(0:mxgrid),bufvir(0:mxgrid), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - bonds_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  Call allocate_bond_pot_arrays()

  remake=.false.

  If (idnode == 0) Open(Unit=ntable, File='TABBND')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

! read mesh resolution

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

  i = Index(record,'#')      ! replace hash as it may occur in
  If (i > 0) record(i:i)=' ' ! TABBND if it's in .xvg format

  Call get_word(record,word)
  cutpot = word_2_real(word)

  If(2.0_wp*cutpot > rcut) Call error(504)

  Call get_word(record,word)
  ngrid = Nint(word_2_real(word))

! check array dimensions

  If (ngrid > mxgrid-4) Then
     Call warning(270,Real(ngrid,wp),Real(mxgrid-4,wp),0.0_wp)
     Call error(48)
  End If

  delpot = cutpot/Real(ngrid,wp)
  dlrpot = cutpot/Real(mxgrid-4,wp)

  safe = .false.
  If( Abs(delpot-dlrpot) < 1.0e-8_wp ) Then
     safe   = .true.
     delpot = dlrpot
  End If
  If (delpot > dlrpot .and. (.not.safe)) Then
     If (idnode == 0) Then
        Write(nrite,"(/,                                          &
           & 'expected (minimum) radial increment : ',1p,e15.7,/, &
           & 'TABBND file        radial increment : ',1p,e15.7)") &
           dlrpot, delpot
        Write(nrite,"(/,                                             &
           & 'expected (minimum) number of grid points : ',0p,i10,/, &
           & 'TABBND file actual number of grid points : ',0p,i10)") &
           mxgrid-4, ngrid
     End If
     Call error(22)
  End If
  safe=.true.

  If (Abs(1.0_wp-(delpot/dlrpot)) > 1.0e-8_wp) Then
     remake=.true.
     If (idnode == 0) Write(nrite,"(/,' TABBND arrays resized for mxgrid = ',i10)") mxgrid-4
  End If

  read_type=0 ! initialise read_type
  Do rtbnd=1,ltpbnd(0)
     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

     i = Index(record,'#')      ! replace hash as it may occur in
     If (i > 0) record(i:i)=' ' ! TABBND if it's in .xvg format

     Call get_word(record,atom1)
     Call get_word(record,atom2)

     katom1=0
     katom2=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'**** entry in TABBND'
        Call error(81)
     End If

! Construct unique name for the tabulated bond

     If (katom1 <= katom2) Then
        idbond = atom1//atom2
     Else
        idbond = atom2//atom1
     End If

! read potential arrays if potential is defined

     itbnd=0
     Do jtbnd=1,ltpbnd(0)
        If (bond_name(jtbnd) == idbond) Then
           Do itbnd=1,mxtbnd
              If (ltpbnd(itbnd) == jtbnd) Exit
           End Do
           Exit
        End If
     End Do

     If (itbnd == 0) Then ! All(bond_name /= idbond)
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'**** entry in TABBND'
        Call error(80)
     End If
     If (Any(read_type == jtbnd)) Then
        If (idnode == 0) Write(nrite,'(a)') '****',atom1,'***',atom2,'**** entry in TABBND'
        Call error(172)
     Else
        read_type(jtbnd)=jtbnd
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

! just in case, linear interpolation for distance ~> 0 (missing in the TABLE)

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
           vbnd(i,jtbnd) = t1 + (t2-t1)*ppp*0.5_wp
           vbnd(i,jtbnd) = vbnd(i,jtbnd)*engunit ! convert to internal units

           vk  = bufvir(l)
           vk1 = bufvir(l+1)
           vk2 = bufvir(l+2)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)
           gbnd(i,jtbnd) = t1 + (t2-t1)*ppp*0.5_wp
           gbnd(i,jtbnd) = gbnd(i,jtbnd)*engunit ! convert to internal units
        End Do
        vbnd(0,jtbnd) = cutpot
        gbnd(0,jtbnd) = 1.0_wp/dlrpot
     Else
        Do i=1,ngrid
           vbnd(i,jtbnd) = bufpot(i)*engunit ! convert to internal units
           gbnd(i,jtbnd) = bufvir(i)*engunit ! convert to internal units
        End Do
        vbnd(0,jtbnd) = cutpot
        gbnd(0,jtbnd) = 1.0_wp/delpot
     End If
  End Do

  If (idnode == 0) Then
     Close(Unit=ntable)
     Write(nrite,'(/,1x,a)') 'potential tables read from TABBND file'
  End If

  Deallocate (read_type,     Stat=fail(1))
  Deallocate (bufpot,bufvir, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'error - bonds_table_read deallocation failure, node: ', idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine bonds_table_read
