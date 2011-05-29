Subroutine nst_h1_lfv                             &
           (lvar,mndis,mxdis,tstep,               &
           iso,degfre,sigma,taut,chit,cint,consv, &
           press,taup,chip,eta,                   &
           stress,strext,ten,elrc,virlrc,         &
           strkin,strknf,strknt,engke,engrot,     &
           imcon,mxshak,tolnce,mxquat,quattol,    &
           megcon,strcon,vircon,                  &
           megpmf,strpmf,virpmf,                  &
           strcom,vircom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian and rotational, singled
! RBs, equations of motion in molecular dynamics
! - leapfrog verlet with Hoover thermostat and barostat (anisotropic
! pressure control)
!
! Parrinello-Rahman type: changing cell shape
!
! Note: (1) this ensemble is modified from its original form as in
!           reference1 to that shown in reference2, and now there is
!           coupling between the thermostat and the barostat
!       (2) this ensemble is not correct when there is an external
!           field applied on the system
!
! reference1: Melchionna, Ciccotti and Holian, Mol Phys 1993, 78, p533
! reference2: Martyna, Tuckerman, Tobias, Klein
!             Mol. Phys., 1996, Vol. 87 (5), p. 1117
!
! iso=0 fully anisotropic barostat
! iso=1 semi-isotropic barostat to constant normal pressure & surface area
! iso=2 semi-isotropic barostat to constant normal pressure & surface tension
!
! reference: Mitsunori Ikeguchi, J Comp Chem 2004, 25, p529
!
! copyright - daresbury laboratory
! author    - w.smith march 2009
! amended   - i.t.todorov may 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax,gcheck
  Use setup_module
  Use domains_module,     Only : map
  Use site_module,        Only : ntpatm,dens,ntpshl,unqshl
  Use config_module,      Only : cell,volm,natms,nlast,nfree, &
                                 lfrzn,lstfre,atmnam,weight,  &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use kinetic_module,     Only : getcom,getvom,getknr, &
                                 kinstresf,kinstrest

  Implicit None

  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,sigma,taut,press,taup
  Integer,           Intent( In    ) :: iso
  Integer(Kind=ip),  Intent( In    ) :: degfre
  Real( Kind = wp ), Intent( InOut ) :: chit,cint,eta(1:9)
  Real( Kind = wp ), Intent(   Out ) :: consv,chip
  Real( Kind = wp ), Intent( InOut ) :: tstep
  Real( Kind = wp ), Intent( In    ) :: stress(1:9),strext(1:9),ten
  Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc
  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                        strknf(1:9),strknt(1:9),engrot

  Integer,           Intent( In    ) :: imcon,mxshak,mxquat
  Real( Kind = wp ), Intent( In    ) :: tolnce,quattol
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon, &
                                        strpmf(1:9),virpmf, &
                                        strcom(1:9),vircom


  Logical,           Save :: newjob = .true. , &
                             unsafe = .false.
  Logical                 :: safe,lv_up,lv_dn
  Integer,           Save :: mxiter,mxkit
  Integer                 :: fail(1:20),matms,iter,kit,i,j,i1,i2, &
                             irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ), Save :: volm0,elrc0,virlrc0,h_z
  Real( Kind = wp ), Save :: qmass,ceng,pmass
  Real( Kind = wp )       :: hstep,rstep
  Real( Kind = wp )       :: chit1,chit2,eta1(1:9),eta2(1:9),chip0
  Real( Kind = wp )       :: cell0(1:9),celprp(1:10)
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             vom(1:3),com(1:3),aaa(1:9),bbb(1:9)
  Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                             opx,opy,opz,fmx,fmy,fmz,       &
                             tqx,tqy,tqz,trx,try,trz,       &
                             oxp,oyp,ozp,oxq,oyq,ozq,       &
                             vpx,vpy,vpz

! uni is the diagonal unit matrix

  Real( Kind = wp ), Parameter :: &
  uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


  Logical,           Allocatable :: lstitr(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)
  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)
  Real( Kind = wp ), Allocatable :: uxt(:),uyt(:),uzt(:)

  Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
  Real( Kind = wp ), Allocatable :: q0t(:),q1t(:),q2t(:),q3t(:)
  Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
  Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
  Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)
  Real( Kind = wp ), Allocatable :: rgdfxx(:),rgdfyy(:),rgdfzz(:)
  Real( Kind = wp ), Allocatable :: rgdtxx(:),rgdtyy(:),rgdtzz(:)
  Real( Kind = wp ), Allocatable :: rgdxxo(:),rgdyyo(:),rgdzzo(:)
  Real( Kind = wp ), Allocatable :: rgdvxo(:),rgdvyo(:),rgdvzo(:)
  Real( Kind = wp ), Allocatable :: rgdoxo(:),rgdoyo(:),rgdozo(:)

  Real( Kind = wp ), Allocatable, Save :: dens0(:)

  fail=0
  If (megcon > 0 .or. megpmf > 0) Then
     Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
     If (megcon > 0) Then
        Allocate (lstopt(1:2,1:mxcons),listot(1:mxatms),          Stat=fail( 2))
        Allocate (dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons),      Stat=fail( 3))
     End If
     If (megpmf > 0) Then
        Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat=fail( 4))
        Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),         Stat=fail( 5))
     End If
  End If
  Allocate (ggx(1:mxlrgd*Max(mxrgd,mxtrgd)),ggy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
            ggz(1:mxlrgd*Max(mxrgd,mxtrgd)),                      Stat=fail( 6))
  Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),            Stat=fail( 7))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
  Allocate (uxt(1:mxatms),uyt(1:mxatms),uzt(1:mxatms),            Stat=fail(11))
  Allocate (q0t(1:mxrgd),q1t(1:mxrgd),q2t(1:mxrgd),q3t(1:mxrgd),  Stat=fail(12))
  Allocate (rgdxxt(1:mxrgd),rgdyyt(1:mxrgd),rgdzzt(1:mxrgd),      Stat=fail(13))
  Allocate (rgdvxt(1:mxrgd),rgdvyt(1:mxrgd),rgdvzt(1:mxrgd),      Stat=fail(14))
  Allocate (rgdoxt(1:mxrgd),rgdoyt(1:mxrgd),rgdozt(1:mxrgd),      Stat=fail(15))
  Allocate (rgdfxx(1:mxrgd),rgdfyy(1:mxrgd),rgdfzz(1:mxrgd),      Stat=fail(16))
  Allocate (rgdtxx(1:mxrgd),rgdtyy(1:mxrgd),rgdtzz(1:mxrgd),      Stat=fail(17))
  Allocate (rgdxxo(1:mxrgd),rgdyyo(1:mxrgd),rgdzzo(1:mxrgd),      Stat=fail(18))
  Allocate (rgdvxo(1:mxrgd),rgdvyo(1:mxrgd),rgdvzo(1:mxrgd),      Stat=fail(19))
  Allocate (rgdoxo(1:mxrgd),rgdoyo(1:mxrgd),rgdozo(1:mxrgd),      Stat=fail(20))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nst_h1 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! store initial values of volume, long range corrections and density

     volm0   = volm
     elrc0   = elrc
     virlrc0 = virlrc

     Allocate (dens0(1:mxatyp), Stat=fail(1))
     If (fail(1) > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'dens0 allocation failure, node: ', idnode
        Call error(0)
     End If
     Do i=1,ntpatm
        dens0(i) = dens(i)
     End Do

! Initialise and get h_z for iso=2

     h_z=0
     If (iso == 2) Then
        Call dcell(cell,celprp)
        h_z=celprp(9)
     End If

! inertia parameters for Nose-Hoover thermostat and barostat

     qmass = 2.0_wp*sigma*taut**2
     tmp   = 2.0_wp*sigma / (boltz*Real(degfre,wp))
     If      (iso == 0) Then
        ceng  = 2.0_wp*sigma + 3.0_wp**2*boltz*tmp
     Else If (iso == 1) Then
        ceng  = 2.0_wp*sigma + 1.0_wp*boltz*tmp
     Else If (iso == 2) Then
        ceng  = 2.0_wp*sigma + 3.0_wp*boltz*tmp
     End If
     pmass = ((2.0_wp*sigma + 3.0_wp*boltz*tmp)/3.0_wp)*taup**2

! trace[eta*transpose(eta)] = trace[eta*eta]: eta is symmetric

     chip0 = Sqrt( eta(1)**2 + 2*eta(2)**2 + 2*eta(3)**2 + eta(5)**2 + 2*eta(6)**2 + eta(9)**2 )

! set number of constraint+pmf shake iterations and general iteration cycles

     mxiter=7
     If (megcon > 0 .or.  megpmf > 0) Then
        mxkit=1
        mxiter=mxiter+1
     End If
     If (megcon > 0 .and. megpmf > 0) mxkit=mxshak

! unsafe positioning due to possibly locally shared RBs

     unsafe=(Any(map == idnode))
  End If

! set matms

  matms=nlast
  If (mxnode == 1) matms=natms

  If (megcon > 0 .or. megpmf > 0) Then
     lstitr(1:natms)=.false. ! initialise lstitr

! construct current bond vectors and listot array (shared
! constraint atoms) for iterative bond algorithms

     If (megcon > 0) Call constraints_tags(imcon,lstitr,lstopt,dxx,dyy,dzz,listot)

! construct current PMF constraint vectors and shared description
! for iterative PMF constraint algorithms

     If (megpmf > 0) Call pmf_tags(imcon,lstitr,indpmf,pxx,pyy,pzz)
  End If

! Get the RB particles vectors wrt the RB's COM

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

! COM distances

           ggx(krgd)=xxx(i)-rgdxxx(irgd)
           ggy(krgd)=yyy(i)-rgdyyy(irgd)
           ggz(krgd)=zzz(i)-rgdzzz(irgd)
        End Do
     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,ggx,ggy,ggz)

! timestep derivatives

  hstep = 0.5_wp*tstep
  rstep = 1.0_wp/tstep
  lv_up = .false.
  lv_dn = .false.

! Get RB COM stress and virial

  Call rigid_bodies_stress(strcom,ggx,ggy,ggz)
  vircom=-(strcom(1)+strcom(5)+strcom(9))

! leapfrog verlet algorithm (starts with velocities at half step)!!!

! store initial values

  Do i=1,matms
     xxt(i) = xxx(i)
     yyt(i) = yyy(i)
     zzt(i) = zzz(i)

     vxt(i) = vxx(i)
     vyt(i) = vyy(i)
     vzt(i) = vzz(i)

     fxt(i) = fxx(i)
     fyt(i) = fyy(i)
     fzt(i) = fzz(i)
  End Do

  Do irgd=1,ntrgd
     q0t(irgd)=q0(irgd)
     q1t(irgd)=q1(irgd)
     q2t(irgd)=q2(irgd)
     q3t(irgd)=q3(irgd)

     rgdxxt(irgd) = rgdxxx(irgd)
     rgdyyt(irgd) = rgdyyy(irgd)
     rgdzzt(irgd) = rgdzzz(irgd)

     rgdvxt(irgd) = rgdvxx(irgd)
     rgdvyt(irgd) = rgdvyy(irgd)
     rgdvzt(irgd) = rgdvzz(irgd)

     rgdoxt(irgd) = rgdoxx(irgd)
     rgdoyt(irgd) = rgdoyy(irgd)
     rgdozt(irgd) = rgdozz(irgd)
  End Do

! store temporary cell parameters

  cell0=cell

! Get RB force and torque

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then

        fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
        tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
        Do jrgd=1,lrgd
           krgd=krgd+1

           i=indrgd(jrgd,irgd) ! local index of particle/site

! If the RB has a frozen particle then no net force

           If (rgdfrz(0,rgdtyp) == 0) Then
              fmx=fmx+fxt(i)
              fmy=fmy+fyt(i)
              fmz=fmz+fzt(i)
           End If

           tqx=tqx+ggy(krgd)*fzt(i)-ggz(krgd)*fyt(i)
           tqy=tqy+ggz(krgd)*fxt(i)-ggx(krgd)*fzt(i)
           tqz=tqz+ggx(krgd)*fyt(i)-ggy(krgd)*fxt(i)
        End Do

! If the RB has 2+ frozen particles (ill=1) the net torque
! must align along the axis of rotation

        If (rgdfrz(0,rgdtyp) > 1) Then
           i1=indrgd(rgdind(1,rgdtyp),irgd)
           i2=indrgd(rgdind(2,rgdtyp),irgd)

           x(1)=xxt(i1)-xxt(i2)
           y(1)=yyt(i1)-yyt(i2)
           z(1)=zzt(i1)-zzt(i2)

           Call images(imcon,cell0,1,x,y,z)

           tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)/(x(1)**2+y(1)**2+z(1)**2)
           tqx=x(1)*tmp
           tqy=y(1)*tmp
           tqz=z(1)*tmp
        End If

! current rotation matrix

        Call getrotmat(q0t(irgd),q1t(irgd),q2t(irgd),q3t(irgd),rot)

! calculate torque in principal frame

        trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
        try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
        trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

! store COM force and torque

        rgdfxx(irgd)=fmx
        rgdfyy(irgd)=fmy
        rgdfzz(irgd)=fmz

        rgdtxx(irgd)=trx
        rgdtyy(irgd)=try
        rgdtzz(irgd)=trz

     Else

        rgdfxx(irgd)=0.0_wp
        rgdfyy(irgd)=0.0_wp
        rgdfzz(irgd)=0.0_wp

        rgdtxx(irgd)=0.0_wp
        rgdtyy(irgd)=0.0_wp
        rgdtzz(irgd)=0.0_wp

     End If
  End Do

! calculate system centre of mass

  Call getcom(natms,weight,xxx,yyy,zzz,com)

100 Continue

! constraint virial and stress tensor

  If (megcon > 0) Then
     vircon=0.0_wp
     strcon=0.0_wp
  End If

! PMF virial and stress tensor

  If (megpmf > 0) Then
     virpmf=0.0_wp
     strpmf=0.0_wp
  End If

! estimates for FPs

  Do j=1,nfree
     i=lstfre(j)

     If (weight(i) > 1.0e-6_wp) Then

! estimate velocity at full step

        tmp=hstep/weight(i)
        uxt(i)=vxt(i)+tmp*fxt(i)
        uyt(i)=vyt(i)+tmp*fyt(i)
        uzt(i)=vzt(i)+tmp*fzt(i)

! first estimate of new velocities - no thermostat

        tmp=tstep/weight(i)
        vxx(i)=vxt(i)+tmp*fxt(i)
        vyy(i)=vyt(i)+tmp*fyt(i)
        vzz(i)=vzt(i)+tmp*fzt(i)

! first estimate of new positions

        xxx(i)=xxt(i)+tstep*vxx(i)
        yyy(i)=yyt(i)+tstep*vyy(i)
        zzz(i)=zzt(i)+tstep*vzz(i)

! first estimate of position at half step

        oxt(i)=0.5_wp*(xxt(i)+xxx(i)) - com(1)
        oyt(i)=0.5_wp*(yyt(i)+yyy(i)) - com(2)
        ozt(i)=0.5_wp*(zzt(i)+zzz(i)) - com(3)

     Else

        uxt(i)=0.0_wp
        uyt(i)=0.0_wp
        uzt(i)=0.0_wp

     End If
  End Do

! estimate RB COM and angular velocities at full step

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then

! recover COM force and torque

        fmx=rgdfxx(irgd)
        fmy=rgdfyy(irgd)
        fmz=rgdfzz(irgd)

        trx=rgdtxx(irgd)
        try=rgdtyy(irgd)
        trz=rgdtzz(irgd)

! estimate COM velocity at time step n - no thermostat

        tmp=hstep/rgdwgt(0,rgdtyp)
        rgdvxo(irgd)=rgdvxt(irgd)+tmp*fmx
        rgdvyo(irgd)=rgdvyt(irgd)+tmp*fmy
        rgdvzo(irgd)=rgdvzt(irgd)+tmp*fmz

! first estimate new COM velocity - no thermostat

        tmp=tstep/rgdwgt(0,rgdtyp)
        rgdvxx(irgd)=rgdvxt(irgd)+tmp*fmx
        rgdvyy(irgd)=rgdvyt(irgd)+tmp*fmy
        rgdvzz(irgd)=rgdvzt(irgd)+tmp*fmz

! first estimate of new COM positions

        rgdxxx(irgd)=rgdxxt(irgd)+tstep*rgdvxx(irgd)
        rgdyyy(irgd)=rgdyyt(irgd)+tstep*rgdvyy(irgd)
        rgdzzz(irgd)=rgdzzt(irgd)+tstep*rgdvzz(irgd)

! first estimate of position at half step

        rgdxxo(irgd)=0.5_wp*(rgdxxt(irgd)+rgdxxx(irgd)) - com(1)
        rgdyyo(irgd)=0.5_wp*(rgdyyt(irgd)+rgdyyy(irgd)) - com(2)
        rgdzzo(irgd)=0.5_wp*(rgdzzt(irgd)+rgdzzz(irgd)) - com(3)

! calculate rigid body rotational motion
! iterate angular velocity for time step n

        oxp=rgdoxt(irgd)
        oyp=rgdoyt(irgd)
        ozp=rgdozt(irgd)
        Do i=1,5
           opx=rgdrix(2,rgdtyp)*(trx + (rgdriy(1,rgdtyp)-rgdriz(1,rgdtyp))*oyp*ozp)
           opy=rgdriy(2,rgdtyp)*(try + (rgdriz(1,rgdtyp)-rgdrix(1,rgdtyp))*ozp*oxp)
           opz=rgdriz(2,rgdtyp)*(trz + (rgdrix(1,rgdtyp)-rgdriy(1,rgdtyp))*oxp*oyp)

! improved angular velocity at time step n

           oxp=rgdoxt(irgd)+hstep*opx
           oyp=rgdoyt(irgd)+hstep*opy
           ozp=rgdozt(irgd)+hstep*opz
        End Do
        rgdoxo(irgd)=oxp
        rgdoyo(irgd)=oyp
        rgdozo(irgd)=ozp

     Else

        rgdxxo(irgd)=rgdxxt(irgd) - com(1)
        rgdyyo(irgd)=rgdyyt(irgd) - com(2)
        rgdzzo(irgd)=rgdzzt(irgd) - com(3)

        rgdvxo(irgd)=0.0_wp
        rgdvyo(irgd)=0.0_wp
        rgdvzo(irgd)=0.0_wp

        rgdoxo(irgd)=0.0_wp
        rgdoyo(irgd)=0.0_wp
        rgdozo(irgd)=0.0_wp

     End If
  End Do

! calculate kinetic energy and stress

  Call kinstresf(uxt,uyt,uzt,strknf)
  Call kinstrest(rgdvxo,rgdvyo,rgdvzo,strknt)

  strkin=strknf+strknt
  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! update rotational energy at full step

  engrot=getknr(rgdoxo,rgdoyo,rgdozo)

! propagate chit and eta sets and couple
! (strcon,strpmf,chit2,eta2 are zero!!!)

  chit2=0.0_wp
  eta2 =0.0_wp
  chip0=0.0_wp

! split anisotropic from semi-isotropic barostats (iso=0,1,2)

  chit1=chit+tstep*(-ceng)/qmass
  If (iso == 0) Then
     eta1=(eta+tstep*((stress+strkin+strcom) - (press*uni+strext)*volm)/pmass)
  Else
     If      (iso == 1) Then
        eta1(1:8)=0.0_wp
     Else If (iso == 2) Then
        eta1(1)=(eta(1)+tstep*((stress(1)+strkin(1)+strcom(1)) - (press+strext(1)-ten/h_z)*volm)/pmass)
        eta1(2:4)=0.0_wp
        eta1(5)=(eta(5)+tstep*((stress(5)+strkin(5)+strcom(5)) - (press+strext(5)-ten/h_z)*volm)/pmass)
        eta1(6:8)=0.0_wp
     End If
     eta1(9)=(eta(9)+tstep*((stress(9)+strkin(9)+strcom(9)) - (press+strext(9))*volm)/pmass)
  End If

  chit2 = 0.5_wp*(chit+chit1)
  eta2 = 0.5_wp*(eta+eta1)

! trace[eta*transpose(eta)] = trace[eta*eta]: eta is symmetric

  chip0 = Sqrt( eta2(1)**2 + 2*eta2(2)**2 + 2*eta2(3)**2 + eta2(5)**2 + 2*eta2(6)**2 + eta2(9)**2 )

! iterate forces, vircon, virpmf, chit and eta

  Do iter=1,mxiter

! update velocity and position using Nose-Hoover thermostating
! and barostating in leapfrog verlet scheme - FPs

     Do j=1,nfree
        i=lstfre(j)

        If (weight(i) > 1.0e-6_wp) Then
           tmp=1.0_wp/weight(i)

           vxx(i)=vxt(i) + tstep*( tmp*fxx(i) - (chit2*uxt(i) + eta2(1)*uxt(i)+eta2(2)*uyt(i)+eta2(3)*uzt(i)) )
           vyy(i)=vyt(i) + tstep*( tmp*fyy(i) - (chit2*uyt(i) + eta2(2)*uxt(i)+eta2(5)*uyt(i)+eta2(6)*uzt(i)) )
           vzz(i)=vzt(i) + tstep*( tmp*fzz(i) - (chit2*uzt(i) + eta2(3)*uxt(i)+eta2(6)*uyt(i)+eta2(9)*uzt(i)) )

           xxx(i)=xxt(i) + tstep*(vxx(i) + oxt(i)*eta1(1)+oyt(i)*eta1(2)+ozt(i)*eta1(3))
           yyy(i)=yyt(i) + tstep*(vyy(i) + oxt(i)*eta1(2)+oyt(i)*eta1(5)+ozt(i)*eta1(6))
           zzz(i)=zzt(i) + tstep*(vzz(i) + oxt(i)*eta1(3)+oyt(i)*eta1(6)+ozt(i)*eta1(9))
        End If
     End Do

! SHAKE procedures

     If (megcon > 0 .or. megpmf > 0) Then
        safe=.false.
        kit =0

! scale cell vectors: second order taylor expansion of Exp(tstep*eta1)

        aaa=tstep*eta1
        Call mat_mul(aaa,aaa,bbb)
        aaa=uni+aaa+0.5_wp*bbb
        Call mat_mul(aaa,cell0,cell)

! store integrated positions

        Do j=1,nfree
           i=lstfre(j)

           If (lstitr(i)) Then
              oxt(i)=xxx(i)
              oyt(i)=yyy(i)
              ozt(i)=zzz(i)
           End If
        End Do

        Do While ((.not.safe) .and. kit <= mxkit)
           kit=kit+1

           If (megcon > 0) Then

! apply constraint correction: vircon,strcon - constraint virial,stress

              Call constraints_shake_lfv &
           (imcon,mxshak,tolnce,tstep, &
           lstopt,dxx,dyy,dzz,listot,  &
           xxx,yyy,zzz,str,vir)

! constraint virial and stress tensor

              vircon=vircon+vir
              strcon=strcon+str

              safe=.true.
           End If

           If (megpmf > 0) Then

! apply PMF correction: virpmf,strpmf - PMF constraint virial,stress

              Call pmf_shake_lfv       &
           (imcon,mxshak,tolnce,tstep, &
           indpmf,pxx,pyy,pzz,         &
           xxx,yyy,zzz,str,vir)

! PMF virial and stress tensor

              virpmf=virpmf+vir
              strpmf=strpmf+str

              safe=(Abs(vir) <= zero_plus)
           End If
        End Do

        If (.not.safe) Call error(478)

! calculate velocity and force correction

        Do j=1,nfree
           i=lstfre(j)

           If (lstitr(i)) Then
              xt=(xxx(i)-oxt(i))*rstep
              yt=(yyy(i)-oyt(i))*rstep
              zt=(zzz(i)-ozt(i))*rstep

              vxx(i)=vxx(i)+xt
              vyy(i)=vyy(i)+yt
              vzz(i)=vzz(i)+zt

              tmp=weight(i)*rstep
              xt=xt*tmp
              yt=yt*tmp
              zt=zt*tmp

              fxx(i)=fxx(i)+xt
              fyy(i)=fyy(i)+yt
              fzz(i)=fzz(i)+zt
           End If
        End Do
     End If

! calculate velocity of FPs at full step

     Do j=1,nfree
        i=lstfre(j)

        uxt(i)=0.5_wp*(vxt(i)+vxx(i))
        uyt(i)=0.5_wp*(vyt(i)+vyy(i))
        uzt(i)=0.5_wp*(vzt(i)+vzz(i))
     End Do

! LF update RB COM and angular velocities, and COM positions

     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

! For all good RBs

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then

! recover COM force and torque

           fmx=rgdfxx(irgd)
           fmy=rgdfyy(irgd)
           fmz=rgdfzz(irgd)

           trx=rgdtxx(irgd)
           try=rgdtyy(irgd)
           trz=rgdtzz(irgd)

! calculate rigid body rotational motion
! iterate angular velocity for time step n

           oxp=rgdoxo(irgd)
           oyp=rgdoyo(irgd)
           ozp=rgdozo(irgd)
           Do i=1,5
              opx=rgdrix(2,rgdtyp)*(trx + (rgdriy(1,rgdtyp)-rgdriz(1,rgdtyp))*oyp*ozp)
              opy=rgdriy(2,rgdtyp)*(try + (rgdriz(1,rgdtyp)-rgdrix(1,rgdtyp))*ozp*oxp)
              opz=rgdriz(2,rgdtyp)*(trz + (rgdrix(1,rgdtyp)-rgdriy(1,rgdtyp))*oxp*oyp)

! improved angular velocity at time step n

              oxp=rgdoxt(irgd)+hstep*(opx-chit2*oxp)
              oyp=rgdoyt(irgd)+hstep*(opy-chit2*oyp)
              ozp=rgdozt(irgd)+hstep*(opz-chit2*ozp)
           End Do
           rgdoxo(irgd)=oxp
           rgdoyo(irgd)=oyp
           rgdozo(irgd)=ozp

! angular velocity at time step n+1/2

           rgdoxx(irgd)=rgdoxt(irgd)+tstep*(opx-chit2*oxp)
           rgdoyy(irgd)=rgdoyt(irgd)+tstep*(opy-chit2*oyp)
           rgdozz(irgd)=rgdozt(irgd)+tstep*(opz-chit2*ozp)

! advance RB COM velocity by leapfrog

           tmp=1.0_wp/rgdwgt(0,rgdtyp)
           rgdvxx(irgd)=rgdvxt(irgd) + tstep*( tmp*fmx - (chit2*rgdvxo(irgd) + &
                        eta2(1)*rgdvxo(irgd)+eta2(2)*rgdvyo(irgd)+eta2(3)*rgdvzo(irgd)) )
           rgdvyy(irgd)=rgdvyt(irgd) + tstep*( tmp*fmy - (chit2*rgdvyo(irgd) + &
                        eta2(2)*rgdvxo(irgd)+eta2(5)*rgdvyo(irgd)+eta2(6)*rgdvzo(irgd)) )
           rgdvzz(irgd)=rgdvzt(irgd) + tstep*( tmp*fmz - (chit2*rgdvzo(irgd) + &
                        eta2(3)*rgdvxo(irgd)+eta2(6)*rgdvyo(irgd)+eta2(9)*rgdvzo(irgd)) )

! Better estimate of full time step velocities

           rgdvxo(irgd) = 0.5_wp*(rgdvxx(irgd)+rgdvxt(irgd))
           rgdvyo(irgd) = 0.5_wp*(rgdvyy(irgd)+rgdvyt(irgd))
           rgdvzo(irgd) = 0.5_wp*(rgdvzz(irgd)+rgdvzt(irgd))

! advance RB position by leapfrog

           rgdxxx(irgd)=rgdxxt(irgd) + tstep*(rgdvxx(irgd) + rgdxxo(irgd)*eta1(1)+rgdyyo(irgd)*eta1(2)+rgdzzo(irgd)*eta1(3))
           rgdyyy(irgd)=rgdyyt(irgd) + tstep*(rgdvyy(irgd) + rgdxxo(irgd)*eta1(2)+rgdyyo(irgd)*eta1(5)+rgdzzo(irgd)*eta1(6))
           rgdzzz(irgd)=rgdzzt(irgd) + tstep*(rgdvzz(irgd) + rgdxxo(irgd)*eta1(3)+rgdyyo(irgd)*eta1(6)+rgdzzo(irgd)*eta1(9))

        Else

! advance RB position

           rgdxxx(irgd)=rgdxxt(irgd) + tstep*(rgdxxo(irgd)*eta1(1)+rgdyyo(irgd)*eta1(2)+rgdzzo(irgd)*eta1(3))
           rgdyyy(irgd)=rgdyyt(irgd) + tstep*(rgdxxo(irgd)*eta1(2)+rgdyyo(irgd)*eta1(5)+rgdzzo(irgd)*eta1(6))
           rgdzzz(irgd)=rgdzzt(irgd) + tstep*(rgdxxo(irgd)*eta1(3)+rgdyyo(irgd)*eta1(6)+rgdzzo(irgd)*eta1(9))

        End If
     End Do

! update rotational energy at full step

     engrot=getknr(rgdoxo,rgdoyo,rgdozo)

     If (iter < mxiter) Then

! calculate position of FPs at half step

        Do j=1,nfree
           i=lstfre(j)

           oxt(i)=0.5_wp*(xxt(i)+xxx(i)) - com(1)
           oyt(i)=0.5_wp*(yyt(i)+yyy(i)) - com(2)
           ozt(i)=0.5_wp*(zzt(i)+zzz(i)) - com(3)
        End Do

! calculate position of RBs at half step

        Do irgd=1,ntrgd
           rgdxxo(irgd)=0.5_wp*(rgdxxt(irgd)+rgdxxx(irgd)) - com(1)
           rgdyyo(irgd)=0.5_wp*(rgdyyt(irgd)+rgdyyy(irgd)) - com(2)
           rgdzzo(irgd)=0.5_wp*(rgdzzt(irgd)+rgdzzz(irgd)) - com(3)
        End Do

! calculate kinetic energy and stress

        Call kinstresf(uxt,uyt,uzt,strknf)
        Call kinstrest(rgdvxo,rgdvyo,rgdvzo,strknt)

        strkin=strknf+strknt
        engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! propagate chit and eta sets and couple
! (strcon,strpmf,chit2,eta2 are new here!!!)

        chit1=chit+tstep*(2.0_wp*(engke+engrot)+pmass*chip0**2-ceng)/qmass
        If (iso == 0) Then
           eta1=(eta+tstep*((stress+strkin+strcon+strpmf+strcom) - &
                            (press*uni+strext)*volm)/pmass)*Exp(-tstep*chit2)
        Else
           If      (iso == 1) Then
              eta1(1:8)=0.0_wp
           Else If (iso == 2) Then
              eta1(1)=(eta(1)+tstep*((stress(1)+strkin(1)+strcon(1)+strpmf(1)+strcom(1)) - &
                                       (press+strext(1)-ten/h_z)*volm)/pmass)*Exp(-tstep*chit2)
              eta1(2:4)=0.0_wp
              eta1(5)=(eta(5)+tstep*((stress(5)+strkin(5)+strcon(5)+strpmf(5)+strcom(5)) - &
                                       (press+strext(5)-ten/h_z)*volm)/pmass)*Exp(-tstep*chit2)
              eta1(6:8)=0.0_wp
           End If
           eta1(9)=(eta(9)+tstep*((stress(9)+strkin(9)+strcon(9)+strpmf(9)+strcom(9)) - &
                                 (press+strext(9))*volm)/pmass)*Exp(-tstep*chit2)
        End If

        chit2 = 0.5_wp*(chit+chit1)
        eta2 = 0.5_wp*(eta+eta1)

! trace[eta*transpose(eta)] = trace[eta*eta]: eta is symmetric

        chip0 = Sqrt( eta2(1)**2 + 2*eta2(2)**2 + 2*eta2(3)**2 + eta2(5)**2 + 2*eta2(6)**2 + eta2(9)**2 )

     End If

  End Do

! scale cell vectors: second order taylor expansion of Exp(tstep*eta1)

  If (megcon == 0 .and. megpmf == 0) Then
     aaa=tstep*eta1
     Call mat_mul(aaa,aaa,bbb)
     aaa=uni+aaa+0.5_wp*bbb
     Call mat_mul(aaa,cell0,cell)
  End If

! Update RB oreintation, and velocity and position of RB's contituents
! Initialise safety flag for quaternion convergence

  safe=.true.

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then

! recover angular velocity at timestep n

        oxp=rgdoxo(irgd)
        oyp=rgdoyo(irgd)
        ozp=rgdozo(irgd)

! angular velocity at time step n+1
! needed for quaternions algorithm

        oxq=2.0_wp*rgdoxx(irgd)-oxp
        oyq=2.0_wp*rgdoyy(irgd)-oyp
        ozq=2.0_wp*rgdozz(irgd)-ozp

! assign new quaternions and get new rotation matrix

        Call q_update                                     &
           (tstep,oxp,oyp,ozp,oxq,oyq,ozq,mxquat,quattol, &
           q0(irgd),q1(irgd),q2(irgd),q3(irgd),safe)
        Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! update RB members positions and halfstep velocities

        Do jrgd=1,lrgd
           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i <= natms) Then
              If (rgdfrz(jrgd,rgdtyp) == 0) Then
                 x(1)=rgdx(jrgd,rgdtyp)
                 y(1)=rgdy(jrgd,rgdtyp)
                 z(1)=rgdz(jrgd,rgdtyp)

! new atomic positions

                 xxx(i)=rot(1)*x(1)+rot(2)*y(1)+rot(3)*z(1) + rgdxxx(irgd)
                 yyy(i)=rot(4)*x(1)+rot(5)*y(1)+rot(6)*z(1) + rgdyyy(irgd)
                 zzz(i)=rot(7)*x(1)+rot(8)*y(1)+rot(9)*z(1) + rgdzzz(irgd)

! new atomic velocities in body frame

                 vpx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                 vpy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                 vpz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

! DD bound positions

                 If (unsafe) Then
                    xt=xxt(i)-com(1)
                    yt=yyt(i)-com(2)
                    zt=zzt(i)-com(3)
                    vxx(i)=xt*aaa(1)+yt*aaa(2)+zt*aaa(3)
                    vyy(i)=xt*aaa(2)+yt*aaa(5)+zt*aaa(6)
                    vzz(i)=xt*aaa(3)+yt*aaa(6)+zt*aaa(9)

                    x(1)=(xxx(i)-com(1))-vxx(i)
                    y(1)=(yyy(i)-com(2))-vyy(i)
                    z(1)=(zzz(i)-com(3))-vzz(i)
                    Call images(imcon,cell,1,x,y,z)
                    xxx(i)=(x(1)+com(1))+vxx(i)
                    yyy(i)=(y(1)+com(2))+vyy(i)
                    zzz(i)=(z(1)+com(3))+vzz(i)
                 End If

! new atomic velocities in lab frame

                 vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rgdvxx(irgd)
                 vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rgdvyy(irgd)
                 vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rgdvzz(irgd)
              Else
                 x(1)=rgdxxx(irgd)-rgdxxt(irgd)
                 y(1)=rgdyyy(irgd)-rgdyyt(irgd)
                 z(1)=rgdzzz(irgd)-rgdzzt(irgd)
                 If (unsafe) Call images(imcon,cell,1,x,y,z) ! DD bound positions
                 xxx(i)=xxt(i)+x(1)
                 yyy(i)=yyt(i)+y(1)
                 zzz(i)=zzt(i)+z(1)
              End If
           End If
        End Do

     Else

! new atomic positions

        Do jrgd=1,lrgd
           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i <= natms) Then
              x(1)=rgdxxx(irgd)-rgdxxt(irgd)
              y(1)=rgdyyy(irgd)-rgdyyt(irgd)
              z(1)=rgdzzz(irgd)-rgdzzt(irgd)
              If (unsafe) Call images(imcon,cell,1,x,y,z) ! DD bound positions
              xxx(i)=xxt(i)+x(1)
              yyy(i)=yyt(i)+y(1)
              zzz(i)=zzt(i)+z(1)
           End If
        End Do

     End If
  End Do

! Check safety for quaternion convergence

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(321)

! check timestep for variable timestep

  If (lvar) Then

! update maximum distance a particle has travelled

     mxdr = 0.0_wp
     Do i=1,natms
        If (.not.Any(unqshl(1:ntpshl) == atmnam(i))) &
           mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
     End Do
     mxdr=Sqrt(mxdr)
     If (mxnode > 1) Call gmax(mxdr)

     If (mxdr < mndis .or. mxdr > mxdis) Then

! scale tstep and derivatives

        If (mxdr > mxdis) Then
           lv_up = .true.
           If (lv_dn) Then
              tstep = 0.75_wp*tstep
              hstep = 0.50_wp*tstep
           Else
              tstep = hstep
              hstep = 0.50_wp*tstep
           End If
           If (idnode == 0) Write(nrite,"(/,1x, &
              & 'timestep decreased, new timestep is:',3x,1p,e12.4,/)") tstep
        End If
        If (mxdr < mndis) Then
           lv_dn = .true.
           If (lv_up) Then
              tstep = 1.50_wp*tstep
              hstep = 0.50_wp*tstep
           Else
              hstep = tstep
              tstep = 2.00_wp*tstep
           End If
           If (idnode == 0) Write(nrite,"(/,1x, &
              & 'timestep increased, new timestep is:',3x,1p,e12.4,/)") tstep
        End If
        rstep = 1.0_wp/tstep

! restore initial conditions

        Do i=1,matms
           fxx(i) = fxt(i)
           fyy(i) = fyt(i)
           fzz(i) = fzt(i)
        End Do

        Do irgd=1,ntrgd
           q0(irgd)=q0t(irgd)
           q1(irgd)=q1t(irgd)
           q2(irgd)=q2t(irgd)
           q3(irgd)=q3t(irgd)
        End Do

! restart lfv

        Go To 100
     End If
  End If

! update volume and construct a 'mock' chip=Tr[eta1]/3

  chip=eta1(1)+eta1(5)+eta1(9)
  volm=volm*Exp(tstep*chip)
  chip=chip/3.0_wp

! adjust long range corrections and number density

  tmp=(volm0/volm)
  elrc=elrc0*tmp
  virlrc=virlrc0*tmp
  Do i=1,ntpatm
     dens(i)=dens0(i)*tmp
  End Do

! get h_z for iso=2

  If (iso == 2) Then
     Call dcell(cell,celprp)
     h_z=celprp(9)
  End If

! update chit and cint (cint is at full timestep)

  chit=chit1
  cint=cint+tstep*chit2

! update eta

  eta=eta1

! conserved quantity less kinetic and potential energy terms
! (it is at full timestep)

  consv = 0.5_wp*qmass*chit2**2 + 0.5_wp*pmass*chip0**2 + ceng*cint + press*volm

! remove system centre of mass velocity

  Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz)

  Do j=1,nfree
     i=lstfre(j)

     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        vxx(i) = vxx(i) - vom(1)
        vyy(i) = vyy(i) - vom(2)
        vzz(i) = vzz(i) - vom(3)
     End If
  End Do

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     If (rgdfrz(0,rgdtyp) == 0) Then
        rgdvxx(irgd) = rgdvxx(irgd) - vom(1)
        rgdvyy(irgd) = rgdvyy(irgd) - vom(2)
        rgdvzz(irgd) = rgdvzz(irgd) - vom(3)

        lrgd=listrgd(-1,irgd)
        Do jrgd=1,lrgd
           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i <= natms) Then
              vxx(i) = vxx(i) - vom(1)
              vyy(i) = vyy(i) - vom(2)
              vzz(i) = vzz(i) - vom(3)
           End If
        End Do
     End If
  End Do

  Call getvom(vom,uxt,uyt,uzt,rgdvxo,rgdvyo,rgdvzo)

  Do j=1,nfree
     i=lstfre(j)

     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        uxt(i) = uxt(i) - vom(1)
        uyt(i) = uyt(i) - vom(2)
        uzt(i) = uzt(i) - vom(3)
     End If
  End Do

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

     If (rgdfrz(0,rgdtyp) == 0) Then
        rgdvxo(irgd) = rgdvxo(irgd) - vom(1)
        rgdvyo(irgd) = rgdvyo(irgd) - vom(2)
        rgdvzo(irgd) = rgdvzo(irgd) - vom(3)

        lrgd=listrgd(-1,irgd)
        Do jrgd=1,lrgd
           i=indrgd(jrgd,irgd) ! local index of particle/site

           If (i <= natms) Then
              uxt(i) = uxt(i) - vom(1)
              uyt(i) = uyt(i) - vom(2)
              uzt(i) = uzt(i) - vom(3)
           End If
        End Do
     End If
  End Do

! update kinetic energy and stress at full step

  Call kinstresf(uxt,uyt,uzt,strknf)
  Call kinstrest(rgdvxo,rgdvyo,rgdvzo,strknt)

  strkin=strknf+strknt
  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  If (megcon > 0 .or. megpmf > 0) Then
     Deallocate (lstitr,            Stat=fail( 1))
     If (megcon > 0) Then
        Deallocate (lstopt,listot,  Stat=fail( 2))
        Deallocate (dxx,dyy,dzz,    Stat=fail( 3))
     End If
     If (megpmf > 0) Then
        Deallocate (indpmf,         Stat=fail( 4))
        Deallocate (pxx,pyy,pzz,    Stat=fail( 5))
     End If
  End If
  Deallocate (ggx,ggy,ggz,          Stat=fail( 6))
  Deallocate (oxt,oyt,ozt,          Stat=fail( 7))
  Deallocate (xxt,yyt,zzt,          Stat=fail( 8))
  Deallocate (vxt,vyt,vzt,          Stat=fail( 9))
  Deallocate (fxt,fyt,fzt,          Stat=fail(10))
  Deallocate (uxt,uyt,uzt,          Stat=fail(11))
  Deallocate (q0t,q1t,q2t,q3t,      Stat=fail(12))
  Deallocate (rgdxxt,rgdyyt,rgdzzt, Stat=fail(13))
  Deallocate (rgdvxt,rgdvyt,rgdvzt, Stat=fail(14))
  Deallocate (rgdoxt,rgdoyt,rgdozt, Stat=fail(15))
  Deallocate (rgdfxx,rgdfyy,rgdfzz, Stat=fail(16))
  Deallocate (rgdtxx,rgdtyy,rgdtzz, Stat=fail(17))
  Deallocate (rgdxxo,rgdyyo,rgdzzo, Stat=fail(18))
  Deallocate (rgdvxo,rgdvyo,rgdvzo, Stat=fail(19))
  Deallocate (rgdoxo,rgdoyo,rgdozo, Stat=fail(20))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nst_h1 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nst_h1_lfv
