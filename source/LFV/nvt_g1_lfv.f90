Subroutine nvt_g1_lfv                          &
           (lvar,mndis,mxdis,mxstp,temp,tstep, &
           sigma,taut,gama,chit,cint,consv,    &
           strkin,strknf,strknt,engke,engrot,  &
           imcon,mxshak,tolnce,mxquat,quattol, &
           megcon,strcon,vircon,               &
           megpmf,strpmf,virpmf,               &
           strcom,vircom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian and rotational, singled
! RBs, equations of motion in molecular dynamics - leapfrog verlet with
! gentle ergodic thermostat - a Nose-Hoover thermostat chained to a
! Langevin thermostat
!
! reference1: B. Leimkuhler, E. Noorizadeh, F. Theil
!             J. Stat. Phys. (2009) 135: 261-277
!
! reference2: A. Samoletov, M.A.J. Chaplain, C.P. Dettmann
!             J. Stat. Phys. (2007) 128, 1321-1336
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax,gcheck,gsum
  Use setup_module
  Use domains_module,     Only : map
  Use site_module,        Only : ntpshl,unqshl
  Use config_module,      Only : cell,natms,nlast,nfree, &
                                 lstfre,atmnam,weight,   &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use langevin_module,    Only : r_0
  Use kinetic_module,     Only : getknf,getknt,getknr,kinstresf,kinstrest

  Implicit None


  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp, &
                                        temp,sigma,taut,gama
  Real( Kind = wp ), Intent( InOut ) :: chit,cint
  Real( Kind = wp ), Intent(   Out ) :: consv
  Real( Kind = wp ), Intent( InOut ) :: tstep
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
  Integer                 :: fail(1:18),matms,iter,kit,i,j,i1,i2, &
                             irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ), Save :: qmass,ceng
  Real( Kind = wp )       :: hstep,rstep,uni
  Real( Kind = wp )       :: chit1,chit2,engkf,engkt
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp,fex
  Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                             opx,opy,opz,fmx,fmy,fmz,       &
                             tqx,tqy,tqz,trx,try,trz,       &
                             oxp,oyp,ozp,oxq,oyq,ozq,       &
                             vpx,vpy,vpz


  Logical,           Allocatable :: lstitr(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)
  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

  Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
  Real( Kind = wp ), Allocatable :: q0t(:),q1t(:),q2t(:),q3t(:)
  Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
  Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
  Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)
  Real( Kind = wp ), Allocatable :: rgdfxx(:),rgdfyy(:),rgdfzz(:)
  Real( Kind = wp ), Allocatable :: rgdtxx(:),rgdtyy(:),rgdtzz(:)
  Real( Kind = wp ), Allocatable :: rgdvxo(:),rgdvyo(:),rgdvzo(:)
  Real( Kind = wp ), Allocatable :: rgdoxo(:),rgdoyo(:),rgdozo(:)

  fail=0
  If (megcon > 0 .or. megpmf > 0) Then
     Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
     If (megcon > 0) Then
        Allocate (lstopt(0:2,1:mxcons),listot(1:mxatms),          Stat=fail( 2))
        Allocate (dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons),      Stat=fail( 3))
     End If
     If (megpmf > 0) Then
        Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat=fail( 4))
        Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),         Stat=fail( 5))
     End If
  End If
  Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), &
                                                                  Stat=fail( 6))
  Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),            Stat=fail( 7))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
  Allocate (q0t(1:mxrgd),q1t(1:mxrgd),q2t(1:mxrgd),q3t(1:mxrgd),  Stat=fail(11))
  Allocate (rgdxxt(1:mxrgd),rgdyyt(1:mxrgd),rgdzzt(1:mxrgd),      Stat=fail(12))
  Allocate (rgdvxt(1:mxrgd),rgdvyt(1:mxrgd),rgdvzt(1:mxrgd),      Stat=fail(13))
  Allocate (rgdoxt(1:mxrgd),rgdoyt(1:mxrgd),rgdozt(1:mxrgd),      Stat=fail(14))
  Allocate (rgdfxx(1:mxrgd),rgdfyy(1:mxrgd),rgdfzz(1:mxrgd),      Stat=fail(15))
  Allocate (rgdtxx(1:mxrgd),rgdtyy(1:mxrgd),rgdtzz(1:mxrgd),      Stat=fail(16))
  Allocate (rgdvxo(1:mxrgd),rgdvyo(1:mxrgd),rgdvzo(1:mxrgd),      Stat=fail(17))
  Allocate (rgdoxo(1:mxrgd),rgdoyo(1:mxrgd),rgdozo(1:mxrgd),      Stat=fail(18))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_g1 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! inertia parameter for Nose-Hoover thermostat

     qmass = 2.0_wp*sigma*taut**2
     ceng  = 2.0_wp*sigma

! set number of constraint+pmf shake iterations and general iteration cycles

     mxiter=3
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

! generate a Gaussian random number for use in the
! Langevin process on the thermostat friction

  r_0=-6.0_wp
  Do i=1,12
     r_0=r_0+uni()
  End Do
  If (mxnode > 1) Then
     Call gsum(r_0)
     r_0=r_0/Sqrt(Real(mxnode,wp))
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

           Call images(imcon,cell,1,x,y,z)

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
        oxt(i)=vxt(i)+tmp*fxt(i)
        oyt(i)=vyt(i)+tmp*fyt(i)
        ozt(i)=vzt(i)+tmp*fzt(i)

! first estimate of new velocities - no thermostat

        tmp=tstep/weight(i)
        vxx(i)=vxt(i)+tmp*fxt(i)
        vyy(i)=vyt(i)+tmp*fyt(i)
        vzz(i)=vzt(i)+tmp*fzt(i)

     Else

        oxt(i)=0.0_wp
        oyt(i)=0.0_wp
        ozt(i)=0.0_wp

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

! estimate COM velocity at time step n

        tmp=hstep/rgdwgt(0,rgdtyp)
        rgdvxo(irgd)=rgdvxt(irgd)+tmp*fmx
        rgdvyo(irgd)=rgdvyt(irgd)+tmp*fmy
        rgdvzo(irgd)=rgdvzt(irgd)+tmp*fmz

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

        rgdvxo(irgd)=0.0_wp
        rgdvyo(irgd)=0.0_wp
        rgdvzo(irgd)=0.0_wp

        rgdoxo(irgd)=0.0_wp
        rgdoyo(irgd)=0.0_wp
        rgdozo(irgd)=0.0_wp

     End If
  End Do

! update rotational energy at full step

  engrot=getknr(rgdoxo,rgdoyo,rgdozo)

! calculate kinetic energy

  engkf=getknf(oxt,oyt,ozt)
  engkt=getknt(rgdvxo,rgdvyo,rgdvzo)

  engke=engkf+engkt

! propagate chit set

  fex=Exp(-gama*tstep)
  chit1=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
        tstep*(2.0_wp*(engke+engrot)-ceng)/qmass
  chit2 = 0.5_wp*(chit+chit1)

! iterate chit and forces

  Do iter=1,mxiter

! update velocity and position using Nose-Hoover thermostatting - FPs

     Do j=1,nfree
        i=lstfre(j)

        If (weight(i) > 1.0e-6_wp) Then
           tmp=1.0_wp/weight(i)

           vxx(i)=vxt(i)+tstep*(tmp*fxx(i) - chit2*oxt(i))
           vyy(i)=vyt(i)+tstep*(tmp*fyy(i) - chit2*oyt(i))
           vzz(i)=vzt(i)+tstep*(tmp*fzz(i) - chit2*ozt(i))

           xxx(i)=xxt(i)+tstep*vxx(i)
           yyy(i)=yyt(i)+tstep*vyy(i)
           zzz(i)=zzt(i)+tstep*vzz(i)
        End If
     End Do

! SHAKE procedures

     If (megcon > 0 .or. megpmf > 0) Then
        safe=.false.
        kit =0

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

        oxt(i)=0.5_wp*(vxt(i)+vxx(i))
        oyt(i)=0.5_wp*(vyt(i)+vyy(i))
        ozt(i)=0.5_wp*(vzt(i)+vzz(i))
     End Do

! LF update RB COM and angular velocities

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

              oxp=rgdoxt(irgd)+hstep*(opx - chit2*oxp)
              oyp=rgdoyt(irgd)+hstep*(opy - chit2*oyp)
              ozp=rgdozt(irgd)+hstep*(opz - chit2*ozp)
           End Do
           rgdoxo(irgd)=oxp
           rgdoyo(irgd)=oyp
           rgdozo(irgd)=ozp

! angular velocity at time step n+1/2

           rgdoxx(irgd)=rgdoxt(irgd)+tstep*(opx - chit2*oxp)
           rgdoyy(irgd)=rgdoyt(irgd)+tstep*(opy - chit2*oyp)
           rgdozz(irgd)=rgdozt(irgd)+tstep*(opz - chit2*ozp)

! advance RB COM velocity by leapfrog

           tmp=1.0_wp/rgdwgt(0,rgdtyp)
           rgdvxx(irgd)=rgdvxt(irgd)+tstep*(tmp*fmx - chit2*rgdvxo(irgd))
           rgdvyy(irgd)=rgdvyt(irgd)+tstep*(tmp*fmy - chit2*rgdvyo(irgd))
           rgdvzz(irgd)=rgdvzt(irgd)+tstep*(tmp*fmz - chit2*rgdvzo(irgd))

! Better estimate of full time step velocities

           rgdvxo(irgd) = 0.5_wp*(rgdvxx(irgd)+rgdvxt(irgd))
           rgdvyo(irgd) = 0.5_wp*(rgdvyy(irgd)+rgdvyt(irgd))
           rgdvzo(irgd) = 0.5_wp*(rgdvzz(irgd)+rgdvzt(irgd))

        End If
     End Do

! update rotational energy at full step

     engrot=getknr(rgdoxo,rgdoyo,rgdozo)

     If (iter < mxiter) Then

! update kinetic energy

        engkf=getknf(oxt,oyt,ozt)
        engkt=getknt(rgdvxo,rgdvyo,rgdvzo)

        engke=engkf+engkt

! propagate chit set

        chit1=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
              tstep*(2.0_wp*(engke+engrot)-ceng)/qmass
        chit2 = 0.5_wp*(chit+chit1)

     End If

  End Do

! Update RB orientation and COM position,
! velocity and position of RB's constituents
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

! advance RB position by leapfrog

        rgdxxx(irgd)=rgdxxt(irgd)+tstep*rgdvxx(irgd)
        rgdyyy(irgd)=rgdyyt(irgd)+tstep*rgdvyy(irgd)
        rgdzzz(irgd)=rgdzzt(irgd)+tstep*rgdvzz(irgd)

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
                    x(1)=xxx(i)-xxt(i)
                    y(1)=yyy(i)-yyt(i)
                    z(1)=zzz(i)-zzt(i)
                    Call images(imcon,cell,1,x,y,z)
                    xxx(i)=x(1)+xxt(i)
                    yyy(i)=y(1)+yyt(i)
                    zzz(i)=z(1)+zzt(i)
                 End If

! new atomic velocities in lab frame

                 vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rgdvxx(irgd)
                 vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rgdvyy(irgd)
                 vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rgdvzz(irgd)
              End If
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

     If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

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
           If (tstep > mxstp) Then
              tstep = mxstp
              hstep = 0.50_wp*tstep
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

! update chit and cint (cint is at full timestep)

  chit=chit1
  cint=cint+tstep*chit2

! conserved quantity less kinetic and potential energy terms
! (it is at full timestep)

  consv = 0.5_wp*qmass*chit2**2 + ceng*cint

! update kinetic energy and stress at full step

  Call kinstresf(oxt,oyt,ozt,strknf)
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
  Deallocate (q0t,q1t,q2t,q3t,      Stat=fail(11))
  Deallocate (rgdxxt,rgdyyt,rgdzzt, Stat=fail(12))
  Deallocate (rgdvxt,rgdvyt,rgdvzt, Stat=fail(13))
  Deallocate (rgdoxt,rgdoyt,rgdozt, Stat=fail(14))
  Deallocate (rgdfxx,rgdfyy,rgdfzz, Stat=fail(15))
  Deallocate (rgdtxx,rgdtyy,rgdtzz, Stat=fail(16))
  Deallocate (rgdvxo,rgdvyo,rgdvzo, Stat=fail(17))
  Deallocate (rgdoxo,rgdoyo,rgdozo, Stat=fail(18))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_g1 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_g1_lfv
