Subroutine rigid_bodies_setup(megatm,megfrz,megrgd,degtra,degrot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing RBs' rotational inertia tesnors
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsum,gmin
  Use site_module
  Use setup_module
  Use config_module,      Only : cell,nlast,ltg,lsite,lfrzn,xxx,yyy,zzz
  Use rigid_bodies_module

  Implicit None

  Integer,           Intent( In    ) :: megatm
  Integer,           Intent( InOut ) :: megfrz,megrgd
  Integer(Kind=ip),  Intent( InOut ) :: degtra,degrot

  Logical           :: safe,pass1,pass2
  Integer           :: fail(1:2),imcon,irgd,jrgd,krgd,lrgd,rgdtyp, &
                       i,ill,i1,i2,i3, nsite,itmols,nrigid,frzrgd, &
                       ifrz, rotrgd,trargd,iatm1,isite1,ntmp
  Real( Kind = wp ) :: rcut,tmp,weight,                            &
                       rotinr(1:3,1:3),rot1(1:3,1:3),rot(1:9),     &
                       rotall,rotxyz,aa(1:9),bb(1:9),det,rsq,dettest

  Integer,           Allocatable :: allrgd(:),fstrgd(:),lstsit(:)
  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)
  Real( Kind = wp ), Allocatable :: buffer(:)

  fail = 0 ; ntmp = mxlrgd*Max(mxrgd,mxtrgd)
  Allocate (allrgd(1:mxtrgd),fstrgd(1:mxrgd),      Stat = fail(1))
  Allocate (gxx(1:ntmp),gyy(1:ntmp), gzz(1:ntmp),  Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_setup allocation failure, node: ', idnode
     Call error(0)
  End If


! Recover/localise imcon and rcut

  imcon=rgdimc
  rcut=rgdrct

! Initialise safety flag

  safe=.true.

! Tag RBs, find their COMs and check their widths to rcut (system cutoff)

  Call rigid_bodies_tags()
  Call rigid_bodies_coms(imcon,xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz)
  Call rigid_bodies_widths(imcon,rcut)

! Find as many as possible different groups of RB units on this domain
! and qualify a representative by the oldest copy of the very first one

  allrgd=0        ! Initialise presence counter (un-encountered yet)
  fstrgd=megatm+1 ! Initialise order of presence (outside particle range)
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     If (allrgd(rgdtyp) == 0) allrgd(rgdtyp)=1

     If (allrgd(rgdtyp) == 1) Then
        i1=indrgd(1,irgd) ! local index of first member
        iatm1=ltg(i1)     ! global index of first member
        If (iatm1 < fstrgd(rgdtyp)) fstrgd(rgdtyp) = iatm1
     End If
  End Do

! Loop over all local, only domain present representatives of unique RB types
! and get principal axis systems of these RB unit types

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     If (allrgd(rgdtyp) == 1 .and. fstrgd(rgdtyp) == ltg(indrgd(1,irgd))) Then

! Get in the local scope of the unit if not fully frozen

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then ! If not fully frozen
           rgdind(0,rgdtyp)=0             ! Not fully frozen (yet)

           Do jrgd=1,lrgd
              krgd=krgd+1

              i=indrgd(jrgd,irgd) ! local index of particle/site

              gxx(krgd)=xxx(i)-rgdxxx(irgd)
              gyy(krgd)=yyy(i)-rgdyyy(irgd)
              gzz(krgd)=zzz(i)-rgdzzz(irgd)
           End Do
        Else                              ! Fully frozen (as if frozen point particle)
           rgdind(0,rgdtyp)=5
           rgdind(1,rgdtyp)=1
           rgdind(2,rgdtyp)=1
           rgdind(3,rgdtyp)=1
        End If
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! Get RB members internal coordinates for these unique RB unit types
! that are not fully frozen

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     If (allrgd(rgdtyp) == 1 .and. fstrgd(rgdtyp) == ltg(indrgd(1,irgd))) Then
        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then ! If not fully frozen
           rotinr=0.0_wp

           Do jrgd=1,lrgd
              krgd=krgd+1

! A slight shortcut here as we care about the diagonaliser only
! M = D(diagonal(d,d,d)) +/- R(rest) = d*1 +/- R.  The diagonaliser
! C does not care about d*1 & the +/- as C^-1*M*C == d*1 +/- C^-1*R*C

              weight=Real(1-rgdfrz(jrgd,rgdtyp),wp)*rgdwgt(jrgd,rgdtyp)

              rotinr(1,1)=rotinr(1,1)+weight*gxx(krgd)**2
              rotinr(2,1)=rotinr(2,1)+weight*gxx(krgd)*gyy(krgd)
              rotinr(3,1)=rotinr(3,1)+weight*gxx(krgd)*gzz(krgd)
              rotinr(2,2)=rotinr(2,2)+weight*gyy(krgd)**2
              rotinr(3,2)=rotinr(3,2)+weight*gyy(krgd)*gzz(krgd)
              rotinr(3,3)=rotinr(3,3)+weight*gzz(krgd)**2
           End Do
           rotinr(1,2)=rotinr(2,1)
           rotinr(1,3)=rotinr(3,1)
           rotinr(2,3)=rotinr(3,2)

! Diagonalise to get eigen values and vectors

           Call jacobi(3,rotinr,rot1)

           rot(1)=rot1(1,1)
           rot(2)=rot1(1,2)
           rot(3)=rot1(1,3)
           rot(4)=rot1(2,1)
           rot(5)=rot1(2,2)
           rot(6)=rot1(2,3)
           rot(7)=rot1(3,1)
           rot(8)=rot1(3,2)
           rot(9)=rot1(3,3)

           krgd=krgd-lrgd
           Do jrgd=1,lrgd
              krgd=krgd+1

              rgdx(jrgd,rgdtyp)=rot(1)*gxx(krgd)+rot(4)*gyy(krgd)+rot(7)*gzz(krgd)
              rgdy(jrgd,rgdtyp)=rot(2)*gxx(krgd)+rot(5)*gyy(krgd)+rot(8)*gzz(krgd)
              rgdz(jrgd,rgdtyp)=rot(3)*gxx(krgd)+rot(6)*gyy(krgd)+rot(9)*gzz(krgd)
           End Do
        End If
     End If
  End Do

! OUT OF DOMAIN SCOPE & DIVE IN GLOBAL SCOPE

! Globalise internal coordinates and rgdind for all RB unit types
! by broadcasting from the lowest rank processor that holds a copy
! of the original representative via global summation.

  If (mxnode > 1) Then
     Allocate (buffer(1:mxtrgd*(4+3*mxlrgd)), Stat = fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_setup allocation failure 1, node: ', idnode
        Call error(0)
     End If

     krgd=0
     Do irgd=1,mxtrgd
        lrgd=lstrgd(0,irgd)

        iatm1=fstrgd(irgd)
        Call gmin(iatm1)

        If (allrgd(irgd) == 1 .and. fstrgd(irgd) == iatm1) Then
           ntmp=idnode
        Else
           ntmp=mxnode
        End If
        Call gmin(ntmp)

        If (idnode == ntmp) Then
           buffer(krgd+1)=Real(rgdind(0,irgd),wp)
           buffer(krgd+2)=Real(rgdind(1,irgd),wp)
           buffer(krgd+3)=Real(rgdind(2,irgd),wp)
           buffer(krgd+4)=Real(rgdind(3,irgd),wp)
           krgd=krgd+4

           Do jrgd=1,lrgd
              buffer(krgd+1)=rgdx(jrgd,irgd)
              buffer(krgd+2)=rgdy(jrgd,irgd)
              buffer(krgd+3)=rgdz(jrgd,irgd)
              krgd=krgd+3
           End Do
        Else
           buffer(krgd+1)=0.0_wp
           buffer(krgd+2)=0.0_wp
           buffer(krgd+3)=0.0_wp
           buffer(krgd+4)=0.0_wp
           krgd=krgd+4

           Do jrgd=1,lrgd
              buffer(krgd+1)=0.0_wp
              buffer(krgd+2)=0.0_wp
              buffer(krgd+3)=0.0_wp
              krgd=krgd+3
           End Do
        End If
     End Do

     Call gsum(buffer(1:krgd))

     krgd=0
     Do irgd=1,mxtrgd
        rgdind(0,irgd)=Nint(buffer(krgd+1))
        rgdind(1,irgd)=Nint(buffer(krgd+2))
        rgdind(2,irgd)=Nint(buffer(krgd+3))
        rgdind(3,irgd)=Nint(buffer(krgd+4))
        krgd=krgd+4

        lrgd=lstrgd(0,irgd)
        Do jrgd=1,lrgd
           rgdx(jrgd,irgd)=buffer(krgd+1)
           rgdy(jrgd,irgd)=buffer(krgd+2)
           rgdz(jrgd,irgd)=buffer(krgd+3)
           krgd=krgd+3
        End Do
     End Do

     Deallocate (buffer, Stat = fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_setup deallocation failure 1, node: ', idnode
        Call error(0)
     End If
  End If

  Deallocate (allrgd,fstrgd, Stat = fail(1))
  Deallocate (gxx,gyy,gzz,   Stat = fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_setup deallocation failure, node: ', idnode
     Call error(0)
  End If

  Do irgd=1,mxtrgd
     lrgd=lstrgd(0,irgd)
     If (rgdfrz(0,irgd) < lrgd) Then ! If not fully frozen
        Do jrgd=1,lrgd

! Impose rounding

           If (Abs(rgdx(jrgd,irgd)) < 1.0e-8_wp) rgdx(jrgd,irgd)=0.0_wp
           If (Abs(rgdy(jrgd,irgd)) < 1.0e-8_wp) rgdy(jrgd,irgd)=0.0_wp
           If (Abs(rgdz(jrgd,irgd)) < 1.0e-8_wp) rgdz(jrgd,irgd)=0.0_wp

! rotational inertia tensor of group type

           weight=Real(1-rgdfrz(jrgd,irgd),wp)*rgdwgt(jrgd,irgd)

           rgdrix(1,irgd) = rgdrix(1,irgd) + weight*(rgdy(jrgd,irgd)**2+rgdz(jrgd,irgd)**2)
           rgdriy(1,irgd) = rgdriy(1,irgd) + weight*(rgdz(jrgd,irgd)**2+rgdx(jrgd,irgd)**2)
           rgdriz(1,irgd) = rgdriz(1,irgd) + weight*(rgdx(jrgd,irgd)**2+rgdy(jrgd,irgd)**2)

        End Do

! set axis system such that: Ixx >= Iyy >= Izz

        rotxyz=Max(rgdrix(1,irgd),rgdriy(1,irgd),rgdriz(1,irgd))

        If (rotxyz >= rgdrix(1,irgd)) Then
           If (rotxyz <= rgdriy(1,irgd)) Then
              Do jrgd=1,lrgd
                 tmp=rgdx(jrgd,irgd)
                 rgdx(jrgd,irgd)=rgdy(jrgd,irgd)
                 rgdy(jrgd,irgd)=-tmp
              End Do
              rgdriy(1,irgd)=rgdrix(1,irgd)
              rgdrix(1,irgd)=rotxyz
           Else If (rotxyz <= rgdriz(1,irgd)) Then
              Do jrgd=1,lrgd
                 tmp=rgdx(jrgd,irgd)
                 rgdx(jrgd,irgd)=rgdz(jrgd,irgd)
                 rgdz(jrgd,irgd)=-tmp
              End Do
              rgdriz(1,irgd)=rgdrix(1,irgd)
              rgdrix(1,irgd)=rotxyz
           End If
        End If

        If (rgdriz(1,irgd) > rgdriy(1,irgd)) Then
           Do jrgd=1,lrgd
              tmp=rgdy(jrgd,irgd)
              rgdy(jrgd,irgd)=rgdz(jrgd,irgd)
              rgdz(jrgd,irgd)=-tmp
           End Do
           tmp=rgdriz(1,irgd)
           rgdriz(1,irgd)=rgdriy(1,irgd)
           rgdriy(1,irgd)=tmp
        End If

        rotall=rgdrix(1,irgd)+rgdriy(1,irgd)+rgdriz(1,irgd)
        If (rotall <= 1.0e-5_wp) rotall=1.0_wp

! test for type of unit (point/linear/bulk RB == ill=2/1/0)
! and get reciprocal of RI in RB unit internal frame of axis

        ill=0
        If (rgdrix(1,irgd)/rotall < 1.0e-5_wp) Then
           ill=ill+1
        Else
           rgdrix(2,irgd)=1.0_wp/rgdrix(1,irgd)
        End If
        If (rgdriy(1,irgd)/rotall < 1.0e-5_wp) Then
           ill=ill+1
        Else
           rgdriy(2,irgd)=1.0_wp/rgdriy(1,irgd)
        End If
        If (rgdriz(1,irgd)/rotall < 1.0e-5_wp) Then
           ill=ill+1
        Else
           rgdriz(2,irgd)=1.0_wp/rgdriz(1,irgd)
        End If

        rgdind(0,irgd)=ill

        If (ill > 1) Then

! point molecules and one particle RBs are not allowed by default!
! also, partly frozen RBs with only massless unfrozen particles are
! caught in read_field!!!

           safe=.false.
           Exit

        Else If (ill == 1) Then

           If      (rgdfrz(0,irgd) == 0) Then

! linear unfrozen molecule

              rgdind(1,irgd)=1
              rgdind(2,irgd)=2

           Else If (rgdfrz(0,irgd) >  1) Then

! RB with 2+ frozen sites in line (not possible for 1 frozen site only)

              i=0
              Do jrgd=1,lrgd
                 If (rgdfrz(jrgd,irgd) == 1) Then
                    i=i+1
                    rgdind(i,irgd)=jrgd
                    If (i == 3) Exit
                 End If
              End Do

           End If

           If (rgdind(3,irgd) == 0) rgdind(3,irgd)=rgdind(1,irgd)

           i1=rgdind(1,irgd)
           i2=rgdind(2,irgd)

           aa(1)=rgdx(i1,irgd)-rgdx(i2,irgd)
           aa(4)=rgdy(i1,irgd)-rgdy(i2,irgd)
           aa(7)=rgdz(i1,irgd)-rgdz(i2,irgd)
           rsq=Sqrt(aa(1)**2+aa(4)**2+aa(7)**2)

           If      (Abs(aa(7)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(4)**2+aa(7)**2)
              aa(2)= 0.0_wp
              aa(5)= aa(7)/rsq
              aa(8)=-aa(4)/rsq
           Else If (Abs(aa(4)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(4)**2+aa(1)**2)
              aa(2)=-aa(4)/rsq
              aa(5)= aa(1)/rsq
              aa(8)= 0.0_wp
           Else If (Abs(aa(1)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(1)**2+aa(7)**2)
              aa(2)=-aa(7)/rsq
              aa(5)= 0.0_wp
              aa(8)= aa(1)/rsq
           End If

           aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
           aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
           aa(9)=aa(1)*aa(5)-aa(4)*aa(2)

           Call invert(aa,bb,det)

! Check aa validity

           If (rgdfrz(0,irgd) /= 1 .and. Abs(det) < 1.0e-5_wp) Then
              safe=.false.
              Exit
           End If

! Store tensor

           Do i=1,9
             rgdaxs(i,irgd)=bb(i)
           End Do

        Else If (ill == 0) Then

! (1) non-linear molecule or (2) RB with 3+ frozen sites
! as at least 3 not in line (as if all were)

           If (rgdfrz(0,irgd) > 1) Then
              rgdind(0,irgd)=4
              rgdind(1,irgd)=1
              rgdind(2,irgd)=1
              rgdind(3,irgd)=1
           Else
              i1=1
              i2=1
              i3=1

              pass1=.true.
              dettest=1.0e-1_wp

              Do While (pass1 .and. i2 < lrgd-1)

                 i2=i2+1
                 i3=i2
                 pass2=.true.

                 Do While (pass2 .and. i3 < lrgd)

                    i3=i3+1

                    aa(1)=rgdx(i1,irgd)-rgdx(i2,irgd)
                    aa(4)=rgdy(i1,irgd)-rgdy(i2,irgd)
                    aa(7)=rgdz(i1,irgd)-rgdz(i2,irgd)

                    aa(2)=rgdx(i1,irgd)-rgdx(i3,irgd)
                    aa(5)=rgdy(i1,irgd)-rgdy(i3,irgd)
                    aa(8)=rgdz(i1,irgd)-rgdz(i3,irgd)

                    aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
                    aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
                    aa(9)=aa(1)*aa(5)-aa(4)*aa(2)

! invert matrix

                    Call invert(aa,bb,det)

! check on size of determinant - to see if the 3 sites are
! too close to being linear for safety.

                    pass2=(Abs(det) < dettest)

                 End Do

                 pass1=(Abs(det) < dettest)

              End Do

! Check aa validity

              If (Abs(det) < dettest) Then
                 safe=.false.
                 Exit
              End If

! store indices used

              rgdind(1,irgd)=i1
              rgdind(2,irgd)=i2
              rgdind(3,irgd)=i3

! store coefficients

              Do i=1,9
                rgdaxs(i,irgd)=bb(i)
              End Do
           End If

        End If
     End If
  End Do

! Report the naughty RB type, number of unit in molecule type and abort

  If (.not.safe) Then
     nrigid=0
  S: Do itmols=1,ntpmls
        Do i=1,numrgd(itmols)
           nrigid=nrigid+1
           If (nrigid == irgd) Exit S
        End Do
     End Do S

     Call warning(307,Real(irgd,wp),Real(itmols,wp),Real(ill,wp))
     Call error(650)
  End If

! Sort out degrees of freedom (redundancy & corrections)
! correct frzsit,dofsit,rgdwgt,rgdwg1 if needed

  Allocate (lstsit(0:mxlrgd*mxtrgd), Stat = fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_setup allocation failure 2, node: ', idnode
     Call error(0)
  End If
  lstsit=0

! initialise rotational and translational DoF

  degtra=Int(0,ip)
  degrot=Int(0,ip)

  nsite =0
  nrigid=0
  Do itmols=1,ntpmls
     ifrz=0

     frzrgd=0
     trargd=0
     rotrgd=0

     Do irgd=1,numrgd(itmols)
        nrigid=nrigid+1

        lrgd=lstrgd(0,nrigid)
        ill=rgdind(0,nrigid)

        If (ill == 1) Then

           If (rgdfrz(0,nrigid) == 0) Then

! linear molecules lose one axis of rotation but com moves about

              rotrgd=rotrgd+2
              trargd=trargd+3

           Else

! not fully frozen RB with 2+ frozen sites in line have
! no COM momentum and rotate only around 1 axis

              rotrgd=rotrgd+1

           End If

        Else If (ill == 0) Then

           If (rgdfrz(0,nrigid) == 0) Then

! proper unfrozen RB with 3 rot DoFs (rot axis)
! and 3 tra DoF (COM moves about)

              trargd=trargd+3
              rotrgd=rotrgd+3

            Else

! not fully frozen RB with 1 frozen site
! no COM momentum with 3 rot DoFs (rot axis)

              rotrgd=rotrgd+3

            End If

        Else If (ill == 4) Then

! As if fully frozen RBs - must get fully frozen then

           rgdind(0,nrigid)=5
           rgdind(1,rgdtyp)=1
           rgdind(2,rgdtyp)=1
           rgdind(3,rgdtyp)=1

           rgdfrz(0,nrigid)=lrgd
           rgdwgt(0,nrigid)=0.0_wp
           rgdwg1(0,nrigid)=Real(lrgd,wp)

           rgdx(:,nrigid)=0.0_wp
           rgdy(:,nrigid)=0.0_wp
           rgdz(:,nrigid)=0.0_wp

           rgdrix(:,nrigid)=0.0_wp
           rgdriy(:,nrigid)=0.0_wp
           rgdriz(:,nrigid)=0.0_wp

           Call warning(305,Real(irgd,wp),Real(itmols,wp),0.0_wp)

           frzrgd=frzrgd+1

           Do jrgd=1,lrgd
              iatm1=lstrgd(jrgd,nrigid)
              isite1=nsite+iatm1

              If (frzsit(isite1) == 0) Then
                 ifrz=ifrz+1

                 frzsit(isite1)=1
                 dofsit(isite1)=0.0_wp

                 lstsit(0)=lstsit(0)+1
                 lstsit(lstsit(0))=isite1
              End If
           End Do

        End If
     End Do

     megfrz=megfrz+ifrz*nummols(itmols)

     megrgd=megrgd-frzrgd*nummols(itmols)
     degtra=degtra+Int(nummols(itmols),ip)*Int(trargd,ip)
     degrot=degrot+Int(nummols(itmols),ip)*Int(rotrgd,ip)

     nsite=nsite+numsit(itmols)
  End Do

! In case of any refreezing changes refresh the local list of frozen atoms

  If (lstsit(0) > 0) Then
     Do i=1,nlast
        lfrzn(i)=frzsit(lsite(i))
     End Do
  End If

  Deallocate (lstsit, Stat = fail(1))
  If (fail(1) > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'rigid_bodies_setup deallocation failure 2, node: ', idnode
     Call error(0)
  End If

! summarise results

  If (idnode==0) Then
     Write(nrite,"(/,/,1x,'summary of rigid body set up')")

     nrigid=0
     Do itmols=1,ntpmls
        Write(nrite,"(/,2x,'in molecule',i6)") itmols

        If (numrgd(itmols) == 0) Write(nrite,"(/,2x,'no rigid bodies specified')")

        Do i=1,numrgd(itmols)
           If (Mod(nrigid,5) == 0) Then
 Write(nrite,"(/,1x,' type :: members :: frozen status :: unfrozen mass :: translational DoF :: rotational DoF')")
 Write(nrite,"(  1x,'               rotational inertia:        x                   y                   z      ')")
           End If
           nrigid=nrigid+1

           lrgd=lstrgd(0,nrigid)
           ifrz=rgdfrz(0,nrigid)
           ill =rgdind(0,nrigid)

           If      (ifrz < lrgd) Then
              If      (ill == 1) Then ! Linear RB
                 If      (ifrz == 0) Then
                    trargd=3
                    rotrgd=2
                 Else If (ifrz >  1) Then
                    trargd=0
                    rotrgd=1
                 End If
              Else If (ill == 0) Then ! Proper RB
                 If      (ifrz == 0) Then
                    trargd=3
                    rotrgd=3
                 Else If (ifrz == 1) Then
                    trargd=0
                    rotrgd=3
                 End If
              End If
           Else If (ifrz == lrgd) Then
              trargd=0
              rotrgd=0
           End If

 Write(nrite,"(/,i5,2x,i6,9x,i6,9x,f13.6,8x,i6,14x,i6)") nrigid,lrgd,ifrz,rgdwgt(0,nrigid),trargd,rotrgd
 Write(nrite,"(30x,3f20.10)") rgdrix(1,nrigid),rgdriy(1,nrigid),rgdriz(1,nrigid)
           If (lrgd > ifrz) Then
 Write(nrite,"(  1x,'         member  ::   coordinates:        x                   y                   z      ')")
              Do jrgd=1,lrgd
 Write(nrite,"(7x,i6,17x,3f20.10)") jrgd,rgdx(jrgd,nrigid),rgdy(jrgd,nrigid),rgdz(jrgd,nrigid)
              End Do
           End If
        End Do
     End Do
  End If
  If (idnode == 0) Write(nrite,Fmt=*)

! equalise sites DoF due to participating in a good RB (not fully frozen)

  nsite =0
  nrigid=0
  Do itmols=1,ntpmls
     Do irgd=1,numrgd(itmols)
        nrigid=nrigid+1

        lrgd=lstrgd(0,nrigid)
        ill=rgdind(0,nrigid)

! 6 = 3(rot) + 3(tra) DoF per RB

        If (ill == 1) Then

           ntmp=0
           Do jrgd=1,lrgd
              iatm1=lstrgd(jrgd,nrigid)
              isite1=nsite+iatm1

              If (dofsit(isite1) > zero_plus) ntmp=ntmp+1
           End Do

           If (rgdfrz(0,nrigid) == 0) Then

! linear molecule losing one axis of rotation - losing 1(rot) DoF

              Do jrgd=1,lrgd
                 iatm1=lstrgd(jrgd,nrigid)
                 isite1=nsite+iatm1

                 If (dofsit(isite1) > zero_plus) dofsit(isite1)=5.0_wp/Real(ntmp,wp)
              End Do

           Else

! RB with 2+ frozen sites in line with members restricted to
! a circular line around one axis - losing 2(rot) & 3(tra) DoF

              Do jrgd=1,lrgd
                 iatm1=lstrgd(jrgd,nrigid)
                 isite1=nsite+iatm1

                 If (dofsit(isite1) > zero_plus) dofsit(isite1)=1.0_wp/Real(ntmp,wp)
              End Do

           End If

        Else If (ill == 0) Then

           ntmp=0
           Do jrgd=1,lrgd
              iatm1=lstrgd(jrgd,nrigid)
              isite1=nsite+iatm1

              If (dofsit(isite1) > zero_plus) ntmp=ntmp+1
           End Do

           If (rgdfrz(0,nrigid) == 0) Then

! Proper RB

              Do jrgd=1,lrgd
                 iatm1=lstrgd(jrgd,nrigid)
                 isite1=nsite+iatm1

                 If (dofsit(isite1) > zero_plus) dofsit(isite1)=6.0_wp/Real(ntmp,wp)
              End Do

           Else If (rgdfrz(0,nrigid) == 1) Then

! RB with 1 frozen site with members restricted
! to a spherical surface - losing 3(tra) DoF

              Do jrgd=1,lrgd
                 iatm1=lstrgd(jrgd,nrigid)
                 isite1=nsite+iatm1

                 If (dofsit(isite1) > zero_plus) dofsit(isite1)=3.0_wp/Real(ntmp,wp)
              End Do

           Else

              Do jrgd=1,lrgd
                 iatm1=lstrgd(jrgd,nrigid)
                 isite1=nsite+iatm1

                 If (dofsit(isite1) > zero_plus) safe=.false.
              End Do

           End If

        End If
     End Do

! Report the naughty RB type, number of unit in molecule type and abort

     If (.not.safe) Then
        Call warning(307,Real(irgd,wp),Real(itmols,wp),Real(ill,wp))
        Call error(644)
     End If

     nsite=nsite+numsit(itmols)
  End Do

! OUT OF GLOBAL SCOPE & BACK IN DOMAIN SCOPE

! set-up quaternions

  Call q_setup()

End Subroutine rigid_bodies_setup
