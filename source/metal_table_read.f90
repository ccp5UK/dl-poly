Subroutine metal_table_read(l_top)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABEAM file (for metal EAM forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 2006
! amended   - i.t.todorov july 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode,mxnode,gsum
  Use setup_module, Only : ntable,nrite,mxgrid,mxbuff,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use metal_module, Only : ntpmet,lstmet,vmet,dmet,fmet
  Use parse_module, Only : get_line,get_word,lower_case,word_2_real

  Implicit None

  Logical, Intent( In    ) :: l_top

  Logical                :: safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 4   ) :: keyword
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail(1:2),i,j,k,ipot,numpot,ktype,ngrid, &
                            cp,cd,ce,katom1,katom2,keymet,k0,jtpatm
  Real( Kind = wp )      :: start,finish

  Integer,           Dimension( : ), Allocatable :: cpair,cdens,cembed
  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet),cembed(1:ntpmet), Stat=fail(1))
  Allocate (buffer(1:mxbuff),                                                Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  cpair=0 ; cp=0
  cdens=0 ; cd=0
  cembed=0 ; ce=0

  If (idnode == 0) Open(Unit=ntable, File='TABEAM')

! skip header record

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100

! read number of potential functions in file

  Call get_line(safe,ntable,record)
  If (.not.safe) Go To 100
  Call get_word(record,word)
  numpot = Nint(word_2_real(word))

  Do ipot=1,numpot

! read data type, atom labels, number of points, start and end

     Call get_line(safe,ntable,record)
     If (.not.safe) Go To 100

! indentify data type

     Call get_word(record,keyword)
     Call lower_case(keyword)
     If      (keyword == 'pair') Then
          ktype = 1
     Else If (keyword == 'dens') Then
          ktype = 2
     Else If (keyword == 'embe') Then
          ktype = 3
     Else
          Call error(151)
     End If

! identify atom types

     Call get_word(record,atom1)
     If (ktype == 1) Then
        Call get_word(record,atom2)
     Else
        atom2 = atom1
     End If

! data specifiers

     Call get_word(record,word)
     ngrid = Nint(word_2_real(word))
     Call get_word(record,word)
     start  = word_2_real(word)
     Call get_word(record,word)
     finish = word_2_real(word)

! check atom identities

     katom1=0
     katom2=0

     Do jtpatm=1,ntpatm
        If (atom1 == unqatm(jtpatm)) katom1=jtpatm
        If (atom2 == unqatm(jtpatm)) katom2=jtpatm
     End Do

     If (katom1 == 0 .or. katom2 == 0) Then
        If (idnode == 0 .and. l_top) &
           Write(nrite,'(a)') '****',atom1,'***',atom2,'****'
        Call  error(81)
     End If

! store working parameters

     buffer(1)=Real(ngrid+4,wp) ! as if there are 4 extra elements after finish
     buffer(2)=start
     buffer(3)=finish
     buffer(4)=(finish-start)/Real(ngrid-1,wp)

     If (idnode == 0 .and. l_top) &
        Write(nrite,"(1x,i10,4x,2a8,3x,2a4,2x,i6,1p,3e15.6)") &
        ipot,atom1,atom2,'EAM-',keyword,ngrid,start,finish,buffer(4)

! limits shifted for DL_POLY interpolation

     buffer(2)=buffer(2)-buffer(4) ! l_int(min) >= 1
     buffer(3)=buffer(3)-buffer(4) ! l_int(max) < ngrid

! check array dimensions

     If (ngrid+4 > mxbuff) Then
        Call warning(270,Real(ngrid+4,wp),Real(mxbuff,wp),0.0_wp)
        Call error(48)
     End If

     keymet=(Max(katom1,katom2)*(Max(katom1,katom2)-1))/2 + Min(katom1,katom2)
     k0=lstmet(keymet)

! check for undefined potential

     If (k0 == 0) Call error(508)

! read in potential arrays

     Do i=1,(ngrid+3)/4
        k=Min(4,ngrid-(i-1)*4)
        If (idnode == 0) Then
           Read(Unit=ntable, Fmt=*, End=100) (buffer(4*i+j),j=1,k)
        Else
           buffer(4*i+1:4*i+k)=0.0_wp
        End If
     End Do

     If (mxnode > 1) Call gsum(buffer(5:ngrid+4))

! copy data to internal arrays

     If       (ktype == 1) Then

! pair potential terms

        cp=cp+1
        If (Any(cpair(1:cp-1) == k0)) Then
           Call error(509)
        Else
           cpair(cp)=k0
        End If

        vmet(1,k0,1)=buffer(1)
        vmet(2,k0,1)=buffer(2)
        vmet(3,k0,1)=buffer(3)
        vmet(4,k0,1)=buffer(4)

        Do i=5,mxgrid
           If (i-4 > ngrid) Then
             vmet(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             vmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of pair potential function

        Call metal_table_derivatives(k0,buffer,vmet)

! adapt derivatives for use in interpolation

        Do i=5,ngrid+4
           vmet(i,k0,2)=-(Real(i-4,wp)*buffer(4)+buffer(2))*vmet(i,k0,2)
        End Do

     Else If (ktype == 2) Then

! density function terms

        cd=cd+1
        If (Any(cdens(1:cd-1) == katom1)) Then
           Call error(510)
        Else
           cdens(cd)=katom1
        End If

        dmet(1,katom1,1)=buffer(1)
        dmet(2,katom1,1)=buffer(2)
        dmet(3,katom1,1)=buffer(3)
        dmet(4,katom1,1)=buffer(4)

        Do i=5,mxgrid
           If (i-4 > ngrid) Then
             dmet(i,katom1,1)=0.0_wp
           Else
             dmet(i,katom1,1)=buffer(i)
           End If
        End Do

! calculate derivative of density function

        Call metal_table_derivatives(katom1,buffer,dmet)

! adapt derivatives for use in interpolation

        dmet(1,katom1,2)=0.0_wp
        dmet(2,katom1,2)=0.0_wp
        dmet(3,katom1,2)=0.0_wp
        dmet(4,katom1,2)=0.0_wp

        Do i=5,ngrid+4
           dmet(i,katom1,2)=-(Real(i-4,wp)*buffer(4)+buffer(2))*dmet(i,katom1,2)
        End Do

     Else If (ktype == 3) Then

! embedding function arrays

        ce=ce+1
        If (Any(cembed(1:ce-1) == katom1)) Then
           Call error(511)
        Else
           cembed(ce)=katom1
        End If

        fmet(1,katom1,1)=buffer(1)
        fmet(2,katom1,1)=buffer(2)
        fmet(3,katom1,1)=buffer(3)
        fmet(4,katom1,1)=buffer(4)

        Do i=5,mxgrid
           If (i-4 > ngrid) Then
             fmet(i,katom1,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             fmet(i,katom1,1)=buffer(i)
           End If
        End Do

! calculate derivative of embedding function

        Call metal_table_derivatives(katom1,buffer,fmet)

     End If

  End Do

  If (idnode == 0) Close(Unit=ntable)
  If (idnode == 0 .and. l_top) Write(nrite,'(/,1x,a)') 'potential tables read from TABEAM file'

  Deallocate (cpair,cdens,cembed, Stat=fail(1))
  Deallocate (buffer,             Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_table_read deallocation failure, node: ', idnode
     Call error(0)
  End If

  Return

! end of file error exit

100 Continue

  If (idnode == 0) Close(Unit=ntable)
  Call error(24)

End Subroutine metal_table_read
