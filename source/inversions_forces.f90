Subroutine inversions_forces(isw,imcon,enginv,virinv,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating inversion energy and force terms
!
! isw = 0 - collect statistics
! isw = 1 - calculate forces
! isw = 2 - do both
!
! copyright - daresbury laboratory
! author    - w.smith may 1996
! amended   - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsync,gsum,gcheck
  Use setup_module,      Only : nrite,mxinv,mxginv1,pi
  Use config_module,     Only : cell,natms,nlast,lsi,lsa,lfrzn, &
                                xxx,yyy,zzz,fxx,fyy,fzz
  Use inversions_module, Only : ntinv,keyinv,listinv,prminv, &
                                ltpinv,vinv,ginv,ncfinv,ldfinv,dstinv

  Implicit None

  Integer,                             Intent( In    ) :: isw,imcon
  Real( Kind = wp ),                   Intent(   Out ) :: enginv,virinv
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Logical           :: safe
  Integer           :: fail(1:4),i,j,l,ia,ib,ic,id,kk,keyi,local_index
  Real( Kind = wp ) :: xab,yab,zab,rab2,rrab, xac,yac,zac,rac2,rrac,   &
                       xad,yad,zad,rad2,rrad, rbc,rcd,rdb,             &
                       ubx,uby,ubz,ubn,rub, vbx,vby,vbz,vbn,rvb,wwb,   &
                       ucx,ucy,ucz,ucn,ruc, vcx,vcy,vcz,vcn,rvc,wwc,   &
                       udx,udy,udz,udn,rud, vdx,vdy,vdz,vdn,rvd,wwd,   &
                       cosb,cosc,cosd, rubc,rubd,rucd,rucb,rudb,rudc,  &
                       rvbc,rvbd,rvcd,rvcb,rvdb,rvdc,                  &
                       thb,thc,thd,  k,th0,cos0,a,b,m,                 &
                       uuu,uu2,uun,uux,uuy,uuz,                        &
                       pterm,vterm,gamma,gamb,gamc,gamd,               &
                       rdelth,rdr,ppp,vk,vk1,vk2,t1,t2,                &
                       fax,fay,faz, fbx,fby,fbz,                       &
                       fcx,fcy,fcz, fdx,fdy,fdz,                       &
                       strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Real( Kind = wp ), Allocatable :: xdac(:),ydac(:),zdac(:)
  Real( Kind = wp ), Allocatable :: xdad(:),ydad(:),zdad(:)

  fail=0
  Allocate (lunsafe(1:mxinv),lstopt(0:4,1:mxinv),      Stat=fail(1))
  Allocate (xdab(1:mxinv),ydab(1:mxinv),zdab(1:mxinv), Stat=fail(2))
  Allocate (xdac(1:mxinv),ydac(1:mxinv),zdac(1:mxinv), Stat=fail(3))
  Allocate (xdad(1:mxinv),ydad(1:mxinv),zdad(1:mxinv), Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'inversions_forces allocation failure, node: ', idnode
     Call error(0)
  End If


! calculate atom separation vectors

  Do i=1,ntinv
     lunsafe(i)=.false.

! indices of atoms involved

     ia=local_index(listinv(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listinv(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib
     ic=local_index(listinv(3,i),nlast,lsi,lsa) ; lstopt(3,i)=ic
     id=local_index(listinv(4,i),nlast,lsi,lsa) ; lstopt(4,i)=id

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0 .and. ic > 0 .and. id > 0) Then !Tag
        If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic)*lfrzn(id) == 0) Then
           If (ia <= natms .or. ib <= natms .or. ic <= natms .or. id <= natms) Then
              lstopt(0,i)=1
            End If
        End If
     Else                                                    ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
              (ib > 0 .and. ib <= natms) .or.   &
              (ic > 0 .and. ic <= natms) .or.   &
              (id > 0 .and. id <= natms)) .and. &
             (ia == 0 .or. ib == 0 .or. ic == 0 .or. id == 0) ) lunsafe(i)=.true.
     End If

! define components of bond vectors

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ib)-xxx(ia)
        ydab(i)=yyy(ib)-yyy(ia)
        zdab(i)=zzz(ib)-zzz(ia)

! select potential energy function type

        kk=listinv(0,i)
        keyi = Abs(keyinv(kk))

        If (keyi == 5) Then
           xdac(i)=xxx(ic)-xxx(ib)
           ydac(i)=yyy(ic)-yyy(ib)
           zdac(i)=zzz(ic)-zzz(ib)

           xdad(i)=xxx(id)-xxx(ib)
           ydad(i)=yyy(id)-yyy(ib)
           zdad(i)=zzz(id)-zzz(ib)
        Else
           xdac(i)=xxx(ic)-xxx(ia)
           ydac(i)=yyy(ic)-yyy(ia)
           zdac(i)=zzz(ic)-zzz(ia)

           xdad(i)=xxx(id)-xxx(ia)
           ydad(i)=yyy(id)-yyy(ia)
           zdad(i)=zzz(id)-zzz(ia)
        End If
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
!
!        xdac(i)=0.0_wp
!        ydac(i)=0.0_wp
!        zdac(i)=0.0_wp
!
!        xdad(i)=0.0_wp
!        ydad(i)=0.0_wp
!        zdad(i)=0.0_wp
     End If
  End Do

! Check for uncompressed units

  safe = .not. Any(lunsafe(1:ntinv))
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Then
     Do j=0,mxnode-1
        If (idnode == j) Then
           Do i=1,ntinv
              If (lunsafe(i)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                 '*** warning - global unit number', listinv(0,i), &
                 ' , with a head particle number', listinv(1,i),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(134)
  End If

! periodic boundary condition

  Call images(imcon,cell,ntinv,xdab,ydab,zdab)
  Call images(imcon,cell,ntinv,xdac,ydac,zdac)
  Call images(imcon,cell,ntinv,xdad,ydad,zdad)

  If (Mod(isw,3) > 0) Then

! Initialise safety flag

     safe=.true.

! zero inversion energy accumulator

     enginv=0.0_wp
     virinv=0.0_wp

! initialise stress tensor accumulators

     strs1=0.0_wp
     strs2=0.0_wp
     strs3=0.0_wp
     strs5=0.0_wp
     strs6=0.0_wp
     strs9=0.0_wp

  End If

! Recover bin size and increment counter

  If (Mod(isw,2) == 0) Then
     rdelth = Real(mxginv1,wp)/pi
     ncfinv = ncfinv + 1
  End If

! loop over all specified inversions

  Do i=1,ntinv
     If (lstopt(0,i) > 0) Then

! indices of atoms involved

        ia=lstopt(1,i)
        ib=lstopt(2,i)
        ic=lstopt(3,i)
        id=lstopt(4,i)

! define components of bond vectors

        xab=xdab(i)
        yab=ydab(i)
        zab=zdab(i)
        rab2=xab*xab+yab*yab+zab*zab
        rrab=1.0_wp/Sqrt(rab2)

        xac=xdac(i)
        yac=ydac(i)
        zac=zdac(i)
        rac2=xac*xac+yac*yac+zac*zac
        rrac=1.0_wp/Sqrt(rac2)

        xad=xdad(i)
        yad=ydad(i)
        zad=zdad(i)
        rad2=xad*xad+yad*yad+zad*zad
        rrad=1.0_wp/Sqrt(rad2)

! select potential energy function type

        kk=listinv(0,i)
        keyi=keyinv(kk)

        If (keyi == 5) Then

! calculate vector normal to plane

           uux=yac*zad-zac*yad
           uuy=zac*xad-xac*zad
           uuz=xac*yad-yac*xad
           uun=1.0_wp/Sqrt(uux**2+uuy**2+uuz**2)
           uux=uun*uux
           uuy=uun*uuy
           uuz=uun*uuz
           uuu=xab*uux+yab*uuy+zab*uuz

        Else

! scalar products of bond vectors

           rbc=xab*xac+yab*yac+zab*zac
           rcd=xac*xad+yac*yad+zac*zad
           rdb=xad*xab+yad*yab+zad*zab

! calculate bond-angle-plane vectors

           ubx=xac*rrac+xad*rrad
           uby=yac*rrac+yad*rrad
           ubz=zac*rrac+zad*rrad
           ubn=1.0_wp/Sqrt(ubx**2+uby**2+ubz**2)
           ubx=ubn*ubx
           uby=ubn*uby
           ubz=ubn*ubz
           rub=xab*ubx+yab*uby+zab*ubz

           vbx=xac*rrac-xad*rrad
           vby=yac*rrac-yad*rrad
           vbz=zac*rrac-zad*rrad
           vbn=1.0_wp/Sqrt(vbx**2+vby**2+vbz**2)
           vbx=vbn*vbx
           vby=vbn*vby
           vbz=vbn*vbz
           rvb=xab*vbx+yab*vby+zab*vbz
           wwb=Sqrt(rub**2+rvb**2)

           ucx=xad*rrad+xab*rrab
           ucy=yad*rrad+yab*rrab
           ucz=zad*rrad+zab*rrab
           ucn=1.0_wp/Sqrt(ucx**2+ucy**2+ucz**2)
           ucx=ucn*ucx
           ucy=ucn*ucy
           ucz=ucn*ucz
           ruc=xac*ucx+yac*ucy+zac*ucz

           vcx=xad*rrad-xab*rrab
           vcy=yad*rrad-yab*rrab
           vcz=zad*rrad-zab*rrab
           vcn=1.0_wp/Sqrt(vcx**2+vcy**2+vcz**2)
           vcx=vcn*vcx
           vcy=vcn*vcy
           vcz=vcn*vcz
           rvc=xac*vcx+yac*vcy+zac*vcz
           wwc=Sqrt(ruc**2+rvc**2)

           udx=xab*rrab+xac*rrac
           udy=yab*rrab+yac*rrac
           udz=zab*rrab+zac*rrac
           udn=1.0_wp/Sqrt(udx**2+udy**2+udz**2)
           udx=udn*udx
           udy=udn*udy
           udz=udn*udz
           rud=xad*udx+yad*udy+zad*udz

           vdx=xab*rrab-xac*rrac
           vdy=yab*rrab-yac*rrac
           vdz=zab*rrab-zac*rrac
           vdn=1.0_wp/Sqrt(vdx**2+vdy**2+vdz**2)
           vdx=vdn*vdx
           vdy=vdn*vdy
           vdz=vdn*vdz
           rvd=xad*vdx+yad*vdy+zad*vdz
           wwd=Sqrt(rud**2+rvd**2)

! calculate inversion angle cosines

           cosb=wwb*rrab ; If (Abs(cosb) > 1.0_wp) cosb=Sign(1.0_wp,cosb)
           cosc=wwc*rrac ; If (Abs(cosc) > 1.0_wp) cosc=Sign(1.0_wp,cosb)
           cosd=wwd*rrad ; If (Abs(cosd) > 1.0_wp) cosd=Sign(1.0_wp,cosb)

! accumulate the histogram (distribution)

           If (Mod(isw,2) == 0 .and. ib <= natms) Then
              j = ldfinv(kk)

              thb=Acos(cosb)
              l = Min(1+Int(thb*rdelth),mxginv1)
              dstinv(l,j) = dstinv(l,j) + 1.0_wp/3.0_wp

              thc=Acos(cosc)
              l = Min(1+Int(thc*rdelth),mxginv1)
              dstinv(l,j) = dstinv(l,j) + 1.0_wp/3.0_wp

              thd=Acos(cosd)
              l = Min(1+Int(thd*rdelth),mxginv1)
              dstinv(l,j) = dstinv(l,j) + 1.0_wp/3.0_wp
           End If

        End If
        If (isw == 0) Cycle

! calculate potential energy and scalar force term

        If      (keyi == 1) Then

! harmonic inversion potential

           k  =prminv(1,kk)/6.0_wp
           th0=prminv(2,kk)

           thb=Acos(cosb)
           thc=Acos(cosc)
           thd=Acos(cosd)

           pterm=k*((thb-th0)**2+(thc-th0)**2+(thd-th0)**2)
           vterm=0.0_wp
           gamma=0.0_wp
           gamb=0.0_wp
           If (Abs(thb) > 1.0e-10_wp) gamb=2.0_wp*k*(thb-th0)/Sin(thb)
           gamc=0.0_wp
           If (Abs(thc) > 1.0e-10_wp) gamc=2.0_wp*k*(thc-th0)/Sin(thc)
           gamd=0.0_wp
           If (Abs(thd) > 1.0e-10_wp) gamd=2.0_wp*k*(thd-th0)/Sin(thd)

        Else If (keyi == 2) Then

! harmonic cosine inversion potential

           k   =prminv(1,kk)/6.0_wp
           cos0=prminv(2,kk)

           pterm=k*((cosb-cos0)**2+(cosc-cos0)**2+(cosd-cos0)**2)
           vterm=0.0_wp
           gamma=0.0_wp
           gamb=-2.0_wp*k*(cosb-cos0)
           gamc=-2.0_wp*k*(cosc-cos0)
           gamd=-2.0_wp*k*(cosd-cos0)

        Else If (keyi == 3) Then

! planar inversion potentials

           a=prminv(1,kk)

           pterm=a*(1.0_wp-(cosb+cosc+cosd)/3.0_wp)
           vterm=0.0_wp
           gamma=0.0_wp
           gamb=a/3.0_wp
           gamc=a/3.0_wp
           gamd=a/3.0_wp

        Else If (keyi == 4) Then

! extended planar inversion potentials

           k  =prminv(1,kk)/6.0_wp
           th0=prminv(2,kk)
           m  =prminv(3,kk)

           thb=Acos(cosb)
           thc=Acos(cosc)
           thd=Acos(cosd)

           pterm=3.0_wp*k*(1.0_wp-(Cos(m*thb-th0)+Cos(m*thc-th0)+Cos(m*thd-th0))/3.0_wp)
           vterm=0.0_wp
           gamma=0.0_wp
           gamb=0.0_wp
           If (Abs(thb) > 1.0e-10_wp) gamb=k*Sin(m*thb-th0)/Sin(thb)
           gamc=0.0_wp
           If (Abs(thc) > 1.0e-10_wp) gamc=k*Sin(m*thc-th0)/Sin(thc)
           gamd=0.0_wp
           If (Abs(thd) > 1.0e-10_wp) gamd=k*Sin(m*thd-th0)/Sin(thd)

        Else If (keyi == 5) Then

! planar calcite potential

           a=prminv(1,kk)
           b=prminv(2,kk)

           uu2=uuu*uuu
           m=2.0_wp*a+4.0_wp*b*uu2

           pterm=uu2*(a+b*uu2)
           vterm=uu2*m
           gamma=-uuu*m
           gamb=0.0_wp
           gamc=0.0_wp
           gamd=0.0_wp

        Else If (keyi == 20) Then

! TABINV potential

           pterm=0.0_wp

           j = ltpinv(kk)
           rdr = ginv(-1,j) ! 1.0_wp/delpot (in rad^-1)

           thb=Acos(cosb)
           thc=Acos(cosc)
           thd=Acos(cosd)

           l   = Int(thb*rdr)
           ppp = thb*rdr - Real(l,wp)

           vk  = vinv(l,j)
           vk1 = vinv(l+1,j)
           vk2 = vinv(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           pterm = pterm + t1 + (t2-t1)*ppp*0.5_wp

           vk  = ginv(l,j) ; If (l == 0) vk = vk*thb
           vk1 = ginv(l+1,j)
           vk2 = ginv(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           gamb = -t1 + (t2-t1)*ppp*0.5_wp

           l   = Int(thc*rdr)
           ppp = thc*rdr - Real(l,wp)

           vk  = vinv(l,j)
           vk1 = vinv(l+1,j)
           vk2 = vinv(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           pterm = pterm + t1 + (t2-t1)*ppp*0.5_wp

           vk  = ginv(l,j) ; If (l == 0) vk = vk*thc
           vk1 = ginv(l+1,j)
           vk2 = ginv(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           gamc = -t1 + (t2-t1)*ppp*0.5_wp

           l   = Int(thd*rdr)
           ppp = thd*rdr - Real(l,wp)

           vk  = Merge(vinv(l,j), 0.0_wp, l > 0)
           vk1 = vinv(l+1,j)
           vk2 = vinv(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           pterm = pterm + t1 + (t2-t1)*ppp*0.5_wp

           vk  = ginv(l,j) ; If (l == 0) vk = vk*thd
           vk1 = ginv(l+1,j)
           vk2 = ginv(l+2,j)

           t1 = vk  + (vk1 - vk)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           gamd = -t1 + (t2-t1)*ppp*0.5_wp

           vterm=0.0_wp
           gamma=0.0_wp

        Else

! undefined potential

           safe=.false.
           pterm=0.0_wp
           vterm=0.0_wp
           gamb=0.0_wp
           gamc=0.0_wp
           gamd=0.0_wp

        End If

        If (keyinv(kk) == 5) Then

! calculate atomic forces

           fax=-gamma*uux
           fay=-gamma*uuy
           faz=-gamma*uuz

           fcx=gamma*uun*((yad*zab-zad*yab)-uuu*(yad*uuz-zad*uuy))
           fcy=gamma*uun*((zad*xab-xad*zab)-uuu*(zad*uux-xad*uuz))
           fcz=gamma*uun*((xad*yab-yad*xab)-uuu*(xad*uuy-yad*uux))

           fdx=gamma*uun*((yab*zac-zab*yac)-uuu*(zac*uuy-yac*uuz))
           fdy=gamma*uun*((zab*xac-xab*zac)-uuu*(xac*uuz-zac*uux))
           fdz=gamma*uun*((xab*yac-yab*xac)-uuu*(yac*uux-xac*uuy))

           fbx=-(fax+fcx+fdx)
           fby=-(fay+fcy+fdy)
           fbz=-(faz+fcz+fdz)

        Else

! calculate bond and u,v scalar products

           rubc=xab*ucx+yab*ucy+zab*ucz
           rubd=xab*udx+yab*udy+zab*udz
           rucd=xac*udx+yac*udy+zac*udz
           rucb=xac*ubx+yac*uby+zac*ubz
           rudb=xad*ubx+yad*uby+zad*ubz
           rudc=xad*ucx+yad*ucy+zad*ucz

           rvbc=xab*vcx+yab*vcy+zab*vcz
           rvbd=xab*vdx+yab*vdy+zab*vdz
           rvcd=xac*vdx+yac*vdy+zac*vdz
           rvcb=xac*vbx+yac*vby+zac*vbz
           rvdb=xad*vbx+yad*vby+zad*vbz
           rvdc=xad*vcx+yad*vcy+zad*vcz

! calculate atomic forces

           fbx = gamb*(-cosb*xab*rrab**2+rrab*(rub*ubx+rvb*vbx)/wwb) +       &
                 ( ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2) -   &
                   rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2) ) * &
                 gamc*rrac/wwc +                                             &
                 ( rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2) +   &
                   rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2) ) * &
                 gamd*rrad/wwd

           fby = gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb) +       &
                 ( ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2) -   &
                   rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2) ) * &
                 gamc*rrac/wwc +                                             &
                 ( rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2) +   &
                   rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2) ) * &
                 gamd*rrad/wwd

           fbz = gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb) +       &
                 ( ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2) -   &
                   rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2) ) * &
                 gamc*rrac/wwc +                                             &
                 ( rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2) +   &
                   rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2) ) * &
                 gamd*rrad/wwd

           fcx = gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc) +       &
                 ( rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2) -   &
                   rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2) ) * &
                 gamd*rrad/wwd +                                             &
                 ( rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2) +   &
                   rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2) ) * &
                 gamb*rrab/wwb

           fcy = gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc) +       &
                 ( rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2) -   &
                   rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2) ) * &
                 gamd*rrad/wwd +                                             &
                 ( rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2) +   &
                   rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2) ) * &
                 gamb*rrab/wwb

           fcz = gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc) +       &
                 ( rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2) -   &
                   rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2) ) * &
                 gamd*rrad/wwd +                                             &
                 ( rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2) +   &
                   rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2) ) * &
                 gamb*rrab/wwb

           fdx = gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd) +       &
                 ( rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2) -   &
                   rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2) ) * &
                 gamb*rrab/wwb +                                             &
                 ( ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2) +   &
                   rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2) ) * &
                 gamc*rrac/wwc

           fdy = gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd) +       &
                 ( rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2) -   &
                   rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2) ) * &
                 gamb*rrab/wwb +                                             &
                 ( ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2) +   &
                   rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2) ) * &
                 gamc*rrac/wwc

           fdz = gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd) +       &
                 ( rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2) -   &
                   rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2) ) * &
                 gamb*rrab/wwb +                                             &
                 ( ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2) +   &
                   rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2) ) * &
                 gamc*rrac/wwc

           fax = -(fbx+fcx+fdx)
           fay = -(fby+fcy+fdy)
           faz = -(fbz+fcz+fdz)

        End If

        If (ia <= natms) Then

! inversion energy and virial (associated to the head atom)

           enginv=enginv+pterm
           virinv=virinv+vterm

! stress tensor calculation for inversion terms

           If (keyinv(kk) == 5) Then
              strs1 = strs1 + uuu*gamma*uux*uux
              strs2 = strs2 + uuu*gamma*uux*uuy
              strs3 = strs3 + uuu*gamma*uux*uuz
              strs5 = strs5 + uuu*gamma*uuy*uuy
              strs6 = strs6 + uuu*gamma*uuy*uuz
              strs9 = strs9 + uuu*gamma*uuz*uuz
           Else
              strs1 = strs1 + xab*fbx + xac*fcx + xad*fdx
              strs2 = strs2 + yab*fbx + yac*fcx + yad*fdx
              strs3 = strs3 + zab*fbx + zac*fcx + zad*fdx
              strs5 = strs5 + yab*fby + yac*fcy + yad*fdy
              strs6 = strs6 + yab*fbz + yac*fcz + yad*fdz
              strs9 = strs9 + zab*fbz + zac*fcz + zad*fdz
           End If

           fxx(ia)=fxx(ia)+fax
           fyy(ia)=fyy(ia)+fay
           fzz(ia)=fzz(ia)+faz

        End If

        If (ib <= natms) Then

           fxx(ib)=fxx(ib)+fbx
           fyy(ib)=fyy(ib)+fby
           fzz(ib)=fzz(ib)+fbz

        End If

        If (ic <= natms) Then

           fxx(ic)=fxx(ic)+fcx
           fyy(ic)=fyy(ic)+fcy
           fzz(ic)=fzz(ic)+fcz

        End If

        If (id <= natms) Then

           fxx(id)=fxx(id)+fdx
           fyy(id)=fyy(id)+fdy
           fzz(id)=fzz(id)+fdz

        End If

     End If
  End Do

  If (Mod(isw,3) > 0) Then

! complete stress tensor

     stress(1) = stress(1) + strs1
     stress(2) = stress(2) + strs2
     stress(3) = stress(3) + strs3
     stress(4) = stress(4) + strs2
     stress(5) = stress(5) + strs5
     stress(6) = stress(6) + strs6
     stress(7) = stress(7) + strs3
     stress(8) = stress(8) + strs6
     stress(9) = stress(9) + strs9

! check for undefined potentials

     If (mxnode > 1) Call gcheck(safe)
     If (.not.safe) Call error(449)

! global sum of inversion potential and virial

     If (mxnode > 1) Then
        buffer(1)=enginv
        buffer(2)=virinv
        Call gsum(buffer(1:2))
        enginv=buffer(1)
        virinv=buffer(2)
     End If

  End If

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  Deallocate (xdac,ydac,zdac, Stat=fail(3))
  Deallocate (xdad,ydad,zdad, Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'inversions_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine inversions_forces
