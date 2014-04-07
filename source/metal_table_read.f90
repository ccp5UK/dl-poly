Subroutine metal_table_read(l_top)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading potential energy and force arrays
! from TABEAM file (for metal EAM & EEAM forces only)
!
! copyright - daresbury laboratory
! author    - w.smith march 2006
! amended   - i.t.todorov april 2014
! contrib   - r.davidchak (eeam) june 2012
! contrib   - b.palmer (2band) may 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode,mxnode,gsum
  Use setup_module, Only : ntable,nrite,mxgmet,engunit
  Use site_module,  Only : ntpatm,unqatm
  Use metal_module, Only : ntpmet,tabmet,lstmet,vmet,dmet,dmes,fmet,fmes
  Use parse_module, Only : get_line,get_word,lower_case,word_2_real

  Implicit None

  Logical, Intent( In    ) :: l_top

  Logical                :: safe
  Character( Len = 200 ) :: record
  Character( Len = 40  ) :: word
  Character( Len = 4   ) :: keyword
  Character( Len = 8   ) :: atom1,atom2
  Integer                :: fail(1:2),i,j,k,ipot,numpot,ktype,ngrid, &
                            cp,cd,cds,ce,ces,katom1,katom2,keymet,k0,jtpatm
  Real( Kind = wp )      :: start,finish

  Integer,           Dimension( : ), Allocatable :: cpair, cdens,cdnss, &
                                                    cembed,cembds
  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  fail=0
  If      (tabmet == 1) Then ! EAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet),                              &
                                              cembed(1:ntpmet),                  Stat=fail(1))
  Else If (tabmet == 2) Then ! EEAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet**2),                           &
                                              cembed(1:ntpmet),                  Stat=fail(1))
  Else If (tabmet == 3) Then ! 2BEAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet),cdnss(1:ntpmet*(ntpmet+1)/2), &
                                              cembed(1:ntpmet),cembds(1:ntpmet), Stat=fail(1))
  Else If (tabmet == 4) Then ! 2BEEAM
     Allocate (cpair(1:(ntpmet*(ntpmet+1))/2),cdens(1:ntpmet**2),cdnss(1:ntpmet**2), &
                                              cembed(1:ntpmet),cembds(1:ntpmet), Stat=fail(1))
  End If
  Allocate (buffer(1:mxgmet),                                                    Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_table_read allocation failure, node: ', idnode
     Call error(0)
  End If
  cpair=0 ; cp=0
  cdens=0 ; cd=0
  cembed=0 ; ce=0
  If (tabmet == 3 .or. tabmet == 4) Then
    cdnss=0 ; cds=0
    cembds=0 ; ces=0
  End If

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

! identify data type

     Call get_word(record,keyword)
     Call lower_case(keyword)
     If      (keyword == 'pair') Then
          ktype = 1
     Else If (keyword == 'dens' .or. keyword == 'dden') Then
          ktype = 2
     Else If (keyword == 'embe' .or. keyword == 'demb') Then
          ktype = 3
     Else If (keyword == 'sden') Then
          ktype = 4
     Else If (keyword == 'semb') Then
          ktype = 5
     Else
          Call error(151)
     End If

! identify atom types

     Call get_word(record,atom1)
     If (ktype == 1 .or.                                        & ! pair
         (ktype == 2 .and. (tabmet == 2 .or. tabmet == 4)) .or. & ! den for EEAM and dden for 2BEEAM
         (ktype == 4 .and. (tabmet == 3 .or. tabmet == 4))) Then  ! sden for 2B(EAM and EEAM)
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
           Write(nrite,'(a)') '****',atom1,'***',atom2,'**** entry in TABEAM'
        Call error(81)
     End If

! store working parameters

     buffer(1)=Real(ngrid+4,wp) ! as if there are 4 extra elements after finish
     buffer(4)=(finish-start)/Real(ngrid-1,wp)
     buffer(2)=start-5.0_wp*buffer(4)
     buffer(3)=finish

     If (idnode == 0 .and. l_top) &
        Write(nrite,"(1x,i10,4x,2a8,3x,2a4,2x,i6,1p,3e15.6)") &
        ipot,atom1,atom2,'EAM-',keyword,ngrid,start,finish,buffer(4)

! check array dimensions

     If (ngrid+4 > mxgmet) Then
        Call warning(270,Real(ngrid+4,wp),Real(mxgmet,wp),0.0_wp)
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

! Set indices

!        k0=lstmet(keymet)

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

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             vmet(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             vmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of pair potential function

        Call metal_table_derivatives(k0,buffer,Size(vmet,2),vmet)

! adapt derivatives for use in interpolation

        Do i=5,ngrid+4
           vmet(i,k0,2)=-(Real(i,wp)*buffer(4)+buffer(2))*vmet(i,k0,2)
        End Do

     Else If (ktype == 2) Then

! density function terms
! s-density density function terms for EAM & EEAM
! d-density density function terms for 2B(EAM & EEAM)

! Set indices

        If      (tabmet == 1 .or. tabmet == 3) Then ! EAM
           k0=katom1
        Else If (tabmet == 2 .or. tabmet == 4) Then ! EEAM
           k0=(katom1-1)*ntpatm+katom2
        End If

        cd=cd+1
        If (Any(cdens(1:cd-1) == k0)) Then
           Call error(510)
        Else
           cdens(cd)=k0
        End If

        dmet(1,k0,1)=buffer(1)
        dmet(2,k0,1)=buffer(2)
        dmet(3,k0,1)=buffer(3)
        dmet(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             dmet(i,k0,1)=0.0_wp
           Else
             dmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of density function

        Call metal_table_derivatives(k0,buffer,Size(dmet,2),dmet)

! adapt derivatives for use in interpolation

        dmet(1,k0,2)=0.0_wp
        dmet(2,k0,2)=0.0_wp
        dmet(3,k0,2)=0.0_wp
        dmet(4,k0,2)=0.0_wp

        Do i=5,ngrid+4
           dmet(i,k0,2)=-(Real(i,wp)*buffer(4)+buffer(2))*dmet(i,k0,2)
        End Do

     Else If (ktype == 3) Then

! embedding function terms
! s-density embedding function terms for EAM & EEAM
! d-density embedding function terms for 2B(EAM & EEAM)

! Set indices

        k0=katom1

        ce=ce+1
        If (Any(cembed(1:ce-1) == k0)) Then
           Call error(511)
        Else
           cembed(ce)=k0
        End If

        fmet(1,k0,1)=buffer(1)
        fmet(2,k0,1)=buffer(2)
        fmet(3,k0,1)=buffer(3)
        fmet(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             fmet(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             fmet(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of embedding function

        Call metal_table_derivatives(k0,buffer,Size(fmet,2),fmet)

     Else If (ktype == 4) Then

! s-density function terms

! The 2BM formalism for alloys allows for a mixed s-band density: rho_{atom1,atom2} /= 0
! (and in general for the EEAM it may be non-symmetric: rho_{atom1,atom2} may be /= rho_{atom2,atom2})
! Some 2BM models rho_{atom1,atom1}=rho_{atom2,atom2}==0 with rho_{atom1,atom2} /= 0
! whereas others choose not to have mixed s-band densities.

! Set indices

        If (tabmet == 3) Then ! 2BMEAM
!           k0=lstmet(keymet)
        Else If (tabmet == 4) Then ! 2BMEEAM
           k0=(katom1-1)*ntpatm+katom2
        End If

        cds=cds+1
        If (Any(cdnss(1:cds-1) == k0)) Then
           Call error(510)
        Else
           cdnss(cds)=k0
        End If

        dmes(1,k0,1)=buffer(1)
        dmes(2,k0,1)=buffer(2)
        dmes(3,k0,1)=buffer(3)
        dmes(4,k0,1)=buffer(4)

        If (buffer(1) > 5) Then

           Do i=5,mxgmet
              If (i-4 > ngrid) Then
                 dmes(i,k0,1)=0.0_wp
              Else
                 dmes(i,k0,1)=buffer(i)
              End If
           End Do
! calculate derivative of density function

           Call metal_table_derivatives(k0,buffer,Size(dmes,2),dmes)

! adapt derivatives for use in interpolation

           dmes(1,k0,2)=0.0_wp
           dmes(2,k0,2)=0.0_wp
           dmes(3,k0,2)=0.0_wp
           dmes(4,k0,2)=0.0_wp

           Do i=5,ngrid+4
              dmes(i,k0,2)=-(Real(i,wp)*buffer(4)+buffer(2))*dmes(i,k0,2)
           End Do

        End If

     Else If (ktype == 5) Then

! s-embedding function terms

! Set index

        k0=katom1

        ces=ces+1
        If (Any(cembds(1:ces-1) == k0)) Then
           Call error(511)
        Else
           cembds(ces)=k0
        End If

        fmes(1,k0,1)=buffer(1)
        fmes(2,k0,1)=buffer(2)
        fmes(3,k0,1)=buffer(3)
        fmes(4,k0,1)=buffer(4)

        Do i=5,mxgmet
           If (i-4 > ngrid) Then
             fmes(i,k0,1)=0.0_wp
           Else
             buffer(i)=buffer(i)*engunit
             fmes(i,k0,1)=buffer(i)
           End If
        End Do

! calculate derivative of embedding function

        Call metal_table_derivatives(k0,buffer,Size(fmes,2),fmes)

     End If

  End Do

  If (idnode == 0) Close(Unit=ntable)
  If (idnode == 0 .and. l_top) Write(nrite,'(/,1x,a)') 'potential tables read from TABEAM file'

  If      (tabmet == 1 .or. tabmet == 2) Then ! EAM & EEAM
     Deallocate (cpair,cdens,cembed,              Stat=fail(1))
  Else If (tabmet == 3 .or. tabmet == 4) Then ! 2B(EAM & EEAM)
     Deallocate (cpair,cdens,cdnss,cembed,cembds, Stat=fail(1))
  End If
  Deallocate (buffer,                             Stat=fail(2))
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
