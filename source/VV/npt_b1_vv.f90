Subroutine npt_b1_vv                          &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           press,taup,chip,eta,               &
           virtot,elrc,virlrc,                &
           strkin,strknf,strknt,engke,engrot, &
           imcon,mxshak,tolnce,               &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian and rotational, singled
! RBs, equations of motion in molecular dynamics
! - velocity verlet with Berendsen thermostat and isotropic pressure
! control (not symplectic)
!
! isothermal compressibility (beta) set to that of liquid water
! = 0.007372 dl_poly units
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax
  Use setup_module
  Use domains_module,     Only : map
  Use site_module,        Only : ntpatm,dens,ntpshl,unqshl
  Use config_module,      Only : cell,volm,natms,nlast,nfree, &
                                 lfrzn,lstfre,atmnam,weight,  &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use kinetic_module,     Only : getvom
  Use constraints_module, Only : passcon
  Use pmf_module,         Only : passpmf

  Implicit None

  Integer,           Intent( In    ) :: isw
  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp, &
                                        sigma,taut,press,taup
  Real( Kind = wp ), Intent(   Out ) :: chit,chip,eta(1:9)
  Real( Kind = wp ), Intent( InOut ) :: tstep
  Real( Kind = wp ), Intent( In    ) :: virtot
  Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc
  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                        strknf(1:9),strknt(1:9),engrot

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon, &
                                        strpmf(1:9),virpmf, &
                                        strcom(1:9),vircom


  Logical,           Save :: newjob = .true. , &
                             unsafe = .false.
  Logical                 :: safe,lcol,lfst,lv_up,lv_dn
  Integer,           Save :: mxiter,mxkit,kit
  Integer                 :: fail(1:15),matms,iter,i,j,i1,i2, &
                             irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ), Save :: volm0,elrc0,virlrc0
  Real( Kind = wp )       :: hstep,rstep
  Real( Kind = wp )       :: czero(1:9)
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             scale,vom(1:3),pr
  Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                             opx,opy,opz,fmx,fmy,fmz,       &
                             tqx,tqy,tqz,trx,try,trz,       &
                             qt0,qt1,qt2,qt3,               &
                             vpx,vpy,vpz

  Real( Kind = wp ), Parameter :: beta = 7.3728e-3_wp


  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

  Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
  Real( Kind = wp ), Allocatable :: q0t(:),q1t(:),q2t(:),q3t(:)
  Real( Kind = wp ), Allocatable :: p0(:),p1(:),p2(:),p3(:)
  Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
  Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
  Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)

  Real( Kind = wp ), Allocatable, Save :: dens0(:)

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
     Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail( 6))
  End If
  Allocate (ggx(1:mxlrgd*mxrgd),ggy(1:mxlrgd*mxrgd),ggz(1:mxlrgd*mxrgd), &
                                                                  Stat=fail( 7))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
  Allocate (q0t(1:mxrgd),q1t(1:mxrgd),q2t(1:mxrgd),q3t(1:mxrgd),  Stat=fail(11))
  Allocate (p0(1:mxrgd),p1(1:mxrgd),p2(1:mxrgd),p3(1:mxrgd),      Stat=fail(12))
  Allocate (rgdxxt(1:mxrgd),rgdyyt(1:mxrgd),rgdzzt(1:mxrgd),      Stat=fail(13))
  Allocate (rgdvxt(1:mxrgd),rgdvyt(1:mxrgd),rgdvzt(1:mxrgd),      Stat=fail(14))
  Allocate (rgdoxt(1:mxrgd),rgdoyt(1:mxrgd),rgdozt(1:mxrgd),      Stat=fail(15))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'npt_b1 allocation failure, node: ', idnode
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

! set number of constraint+pmf shake iterations and general iteration cycles

     mxiter=1
     If (megcon > 0 .or.  megpmf > 0) Then
        mxkit=1
        mxiter=mxiter+12
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

! first pass of velocity verlet algorithm

  If (isw == 0) Then

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

     czero=cell

100  Continue

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

     krgd=0
     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

! For all good RBs

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then

! calculate COM force and torque

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

              Call images(imcon,czero,1,x,y,z)

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

! calculate quaternion torques

           qt0=2.0_wp*(-q1t(irgd)*trx-q2t(irgd)*try-q3t(irgd)*trz)
           qt1=2.0_wp*( q0t(irgd)*trx-q3t(irgd)*try+q2t(irgd)*trz)
           qt2=2.0_wp*( q3t(irgd)*trx+q0t(irgd)*try-q1t(irgd)*trz)
           qt3=2.0_wp*(-q2t(irgd)*trx+q1t(irgd)*try+q0t(irgd)*trz)

! recover quaternion momenta at start of time step

           opx=rgdoxt(irgd)*rgdrix(1,rgdtyp)
           opy=rgdoyt(irgd)*rgdriy(1,rgdtyp)
           opz=rgdozt(irgd)*rgdriz(1,rgdtyp)

           p0(irgd)=2.0_wp*(-q1t(irgd)*opx-q2t(irgd)*opy-q3t(irgd)*opz)
           p1(irgd)=2.0_wp*( q0t(irgd)*opx-q3t(irgd)*opy+q2t(irgd)*opz)
           p2(irgd)=2.0_wp*( q3t(irgd)*opx+q0t(irgd)*opy-q1t(irgd)*opz)
           p3(irgd)=2.0_wp*(-q2t(irgd)*opx+q1t(irgd)*opy+q0t(irgd)*opz)

! update quaternion momenta to half step

           p0(irgd)=p0(irgd)+hstep*qt0
           p1(irgd)=p1(irgd)+hstep*qt1
           p2(irgd)=p2(irgd)+hstep*qt2
           p3(irgd)=p3(irgd)+hstep*qt3

! update RB COM velocities to half step

           tmp=hstep/rgdwgt(0,rgdtyp)
           rgdvxx(irgd)=rgdvxt(irgd)+tmp*fmx
           rgdvyy(irgd)=rgdvyt(irgd)+tmp*fmy
           rgdvzz(irgd)=rgdvzt(irgd)+tmp*fmz
        End If
     End Do

! iterate forces, vircon and chip

     Do iter=1,mxiter

! Berendsen barostat and thermostat are not coupled
! calculate system pressure: iterate vircon and virpmf

        pr = (2.0_wp*engke-virtot-vircon-virpmf-vircom) / (3.0_wp*volm)

! pressure control variable

        chip = 1.0_wp + beta*tstep*(pr-press)/taup
        scale = chip**(1.0_wp/3.0_wp)

! update cell parameters: isotropic

        cell=scale*czero

! update velocity and position of FPs

        Do j=1,nfree
           i=lstfre(j)

           If (weight(i) > 1.0e-6_wp) Then
              tmp=hstep/weight(i)
              vxx(i)=vxt(i)+tmp*fxx(i)
              vyy(i)=vyt(i)+tmp*fyy(i)
              vzz(i)=vzt(i)+tmp*fzz(i)

              xxx(i)=scale*xxt(i)+tstep*vxx(i)
              yyy(i)=scale*yyt(i)+tstep*vyy(i)
              zzz(i)=scale*zzt(i)+tstep*vzz(i)
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

                 Call constraints_shake_vv &
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

                 Call pmf_shake_vv     &
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

! Collect per step passage statistics for bond and pmf constraints

           If (iter == mxiter) Then
              If (megcon > 0) Then
                 passcon(3,2,1)=passcon(2,2,1)*passcon(3,2,1)
                 passcon(2,2,1)=passcon(2,2,1)+1
                 passcon(3,2,1)=passcon(3,2,1)/passcon(2,2,1)+passcon(1,2,1)/passcon(2,2,1)
                 passcon(4,2,1)=Min(passcon(1,2,1),passcon(4,2,1))
                 passcon(5,2,1)=Max(passcon(1,2,1),passcon(5,2,1))
                 passcon(1,2,1)=0.0_wp ! Reset
              End If

              If (megpmf > 0) Then
                 passpmf(3,2,1)=passpmf(2,2,1)*passpmf(3,2,1)
                 passpmf(2,2,1)=passpmf(2,2,1)+1
                 passpmf(3,2,1)=passpmf(3,2,1)/passpmf(2,2,1)+passpmf(1,2,1)/passpmf(2,2,1)
                 passpmf(4,2,1)=Min(passpmf(1,2,1),passpmf(4,2,1))
                 passpmf(5,2,1)=Max(passpmf(1,2,1),passpmf(5,2,1))
                 passpmf(1,2,1)=0.0_wp ! Reset
              End If
           End If

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

                 tmp=weight(i)/hstep
                 xt=xt*tmp
                 yt=yt*tmp
                 zt=zt*tmp

                 fxx(i)=fxx(i)+xt
                 fyy(i)=fyy(i)+yt
                 fzz(i)=fzz(i)+zt
              End If
           End Do
        End If

     End Do

! update velocity and position of RBs

     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

! For all good RBs

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then

! rotate RB quaternions - update q to full timestep & amend p
! and get new rotation matrix

           Call no_squish                                             &
           (tstep,rgdrix(2,rgdtyp),rgdriy(2,rgdtyp),rgdriz(2,rgdtyp), &
           q0(irgd),q1(irgd),q2(irgd),q3(irgd),p0(irgd),p1(irgd),p2(irgd),p3(irgd))
           Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! update RB angular velocities to half step

           opx=0.5_wp*(-q1(irgd)*p0(irgd)+q0(irgd)*p1(irgd)+q3(irgd)*p2(irgd)-q2(irgd)*p3(irgd))
           opy=0.5_wp*(-q2(irgd)*p0(irgd)-q3(irgd)*p1(irgd)+q0(irgd)*p2(irgd)+q1(irgd)*p3(irgd))
           opz=0.5_wp*(-q3(irgd)*p0(irgd)+q2(irgd)*p1(irgd)-q1(irgd)*p2(irgd)+q0(irgd)*p3(irgd))

           rgdoxx(irgd)=opx*rgdrix(2,rgdtyp)
           rgdoyy(irgd)=opy*rgdriy(2,rgdtyp)
           rgdozz(irgd)=opz*rgdriz(2,rgdtyp)

! update RB COM to full step

           rgdxxx(irgd)=scale*rgdxxt(irgd)+tstep*rgdvxx(irgd)
           rgdyyy(irgd)=scale*rgdyyt(irgd)+tstep*rgdvyy(irgd)
           rgdzzz(irgd)=scale*rgdzzt(irgd)+tstep*rgdvzz(irgd)

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
                       vxx(i)=scale*xxt(i)
                       vyy(i)=scale*yyt(i)
                       vzz(i)=scale*zzt(i)

                       x(1)=xxx(i)-vxx(i)
                       y(1)=yyy(i)-vyy(i)
                       z(1)=zzz(i)-vzz(i)
                       Call images(imcon,cell,1,x,y,z)
                       xxx(i)=x(1)+vxx(i)
                       yyy(i)=y(1)+vyy(i)
                       zzz(i)=z(1)+vzz(i)
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

! update RB COM to full step

           rgdxxx(irgd)=scale*rgdxxt(irgd)
           rgdyyy(irgd)=scale*rgdyyt(irgd)
           rgdzzz(irgd)=scale*rgdzzt(irgd)

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

! restart vv1

           Go To 100
        End If
     End If

! update volume: chip=scale^3

     volm=chip*volm

! adjust long range corrections and number density

     tmp=(volm0/volm)
     elrc=elrc0*tmp
     virlrc=virlrc0*tmp
     Do i=1,ntpatm
        dens(i)=dens0(i)*tmp
     End Do

! construct a 'mock' scaling tensor for xscale

     Do i=2,8
        eta(i)=0.0_wp
     End Do
     eta(1)=scale
     eta(5)=scale
     eta(9)=scale

! second pass of velocity verlet algorithm

  Else

! update velocity of FPs

     Do j=1,nfree
        i=lstfre(j)

        If (weight(i) > 1.0e-6_wp) Then
           tmp=hstep/weight(i)
           vxx(i)=vxx(i)+tmp*fxx(i)
           vyy(i)=vyy(i)+tmp*fyy(i)
           vzz(i)=vzz(i)+tmp*fzz(i)
        End If
     End Do

! RATTLE procedures
! apply velocity corrections to bond and PMF constraints

     If (megcon > 0 .or. megpmf > 0) Then
        Do i=1,kit
           lfst = (i == 1)
           lcol = (i == kit)

           If (megcon > 0) Call constraints_rattle &
           (mxshak,tolnce,tstep,lfst,lcol, &
           lstopt,dxx,dyy,dzz,listot,      &
           vxx,vyy,vzz)

           If (megpmf > 0) Call pmf_rattle &
           (mxshak,tolnce,tstep,lfst,lcol, &
           indpmf,pxx,pyy,pzz,             &
           vxx,vyy,vzz)
        End Do
     End If

! Get RB COM stress and virial

     Call rigid_bodies_stress(strcom,ggx,ggy,ggz)
     vircom=-(strcom(1)+strcom(5)+strcom(9))

! update velocity of RBs

     krgd=0
     Do irgd=1,ntrgd
        rgdtyp=listrgd(0,irgd)

! For all good RBs

        lrgd=listrgd(-1,irgd)
        If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

! calculate COM force and torque

           fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
           tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
           Do jrgd=1,lrgd
              krgd=krgd+1

              i=indrgd(jrgd,irgd) ! local index of particle/site

! If the RB has a frozen particle then no net force

              If (rgdfrz(0,rgdtyp) == 0) Then
                 fmx=fmx+fxx(i)
                 fmy=fmy+fyy(i)
                 fmz=fmz+fzz(i)
              End If

              tqx=tqx+ggy(krgd)*fzz(i)-ggz(krgd)*fyy(i)
              tqy=tqy+ggz(krgd)*fxx(i)-ggx(krgd)*fzz(i)
              tqz=tqz+ggx(krgd)*fyy(i)-ggy(krgd)*fxx(i)
           End Do

! If the RB has 2+ frozen particles (ill=1) the net torque
! must align along the axis of rotation

           If (rgdfrz(0,rgdtyp) > 1) Then
              i1=indrgd(rgdind(1,rgdtyp),irgd)
              i2=indrgd(rgdind(2,rgdtyp),irgd)

              x(1)=xxx(i1)-xxx(i2)
              y(1)=yyy(i1)-yyy(i2)
              z(1)=zzz(i1)-zzz(i2)

              Call images(imcon,cell,1,x,y,z)

              tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)/(x(1)**2+y(1)**2+z(1)**2)
              tqx=x(1)*tmp
              tqy=y(1)*tmp
              tqz=z(1)*tmp
           End If

! current rotation matrix

           Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! calculate torque in principal frame

           trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
           try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
           trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

! calculate quaternion torques

           qt0=2.0_wp*(-q1(irgd)*trx-q2(irgd)*try-q3(irgd)*trz)
           qt1=2.0_wp*( q0(irgd)*trx-q3(irgd)*try+q2(irgd)*trz)
           qt2=2.0_wp*( q3(irgd)*trx+q0(irgd)*try-q1(irgd)*trz)
           qt3=2.0_wp*(-q2(irgd)*trx+q1(irgd)*try+q0(irgd)*trz)

! recover quaternion momenta at half time step

           opx=rgdoxx(irgd)*rgdrix(1,rgdtyp)
           opy=rgdoyy(irgd)*rgdriy(1,rgdtyp)
           opz=rgdozz(irgd)*rgdriz(1,rgdtyp)

           p0(irgd)=2.0_wp*(-q1(irgd)*opx-q2(irgd)*opy-q3(irgd)*opz)
           p1(irgd)=2.0_wp*( q0(irgd)*opx-q3(irgd)*opy+q2(irgd)*opz)
           p2(irgd)=2.0_wp*( q3(irgd)*opx+q0(irgd)*opy-q1(irgd)*opz)
           p3(irgd)=2.0_wp*(-q2(irgd)*opx+q1(irgd)*opy+q0(irgd)*opz)

! update quaternion momenta to full step

           p0(irgd)=p0(irgd)+hstep*qt0
           p1(irgd)=p1(irgd)+hstep*qt1
           p2(irgd)=p2(irgd)+hstep*qt2
           p3(irgd)=p3(irgd)+hstep*qt3

! update RB angular & COM velocities to full step

           opx=0.5_wp*(-q1(irgd)*p0(irgd)+q0(irgd)*p1(irgd)+q3(irgd)*p2(irgd)-q2(irgd)*p3(irgd))
           opy=0.5_wp*(-q2(irgd)*p0(irgd)-q3(irgd)*p1(irgd)+q0(irgd)*p2(irgd)+q1(irgd)*p3(irgd))
           opz=0.5_wp*(-q3(irgd)*p0(irgd)+q2(irgd)*p1(irgd)-q1(irgd)*p2(irgd)+q0(irgd)*p3(irgd))

           rgdoxx(irgd)=opx*rgdrix(2,rgdtyp)
           rgdoyy(irgd)=opy*rgdriy(2,rgdtyp)
           rgdozz(irgd)=opz*rgdriz(2,rgdtyp)

           tmp=hstep/rgdwgt(0,rgdtyp)
           rgdvxx(irgd)=rgdvxx(irgd)+tmp*fmx
           rgdvyy(irgd)=rgdvyy(irgd)+tmp*fmy
           rgdvzz(irgd)=rgdvzz(irgd)+tmp*fmz

! update RB members velocities

           Do jrgd=1,lrgd
              If (rgdfrz(jrgd,rgdtyp) == 0) Then
                 i=indrgd(jrgd,irgd) ! local index of particle/site

                 If (i <= natms) Then
                    x(1)=rgdx(jrgd,rgdtyp)
                    y(1)=rgdy(jrgd,rgdtyp)
                    z(1)=rgdz(jrgd,rgdtyp)

! new atomic velocities in body frame

                    vpx=rgdoyy(irgd)*z(1)-rgdozz(irgd)*y(1)
                    vpy=rgdozz(irgd)*x(1)-rgdoxx(irgd)*z(1)
                    vpz=rgdoxx(irgd)*y(1)-rgdoyy(irgd)*x(1)

! new atomic velocities in lab frame

                    vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rgdvxx(irgd)
                    vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rgdvyy(irgd)
                    vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rgdvzz(irgd)
                 End If
              End If
           End Do

        End If
     End Do

! integrate and apply nvt_b1_scl thermostat - full step

     Call nvt_b1_scl &
           (1,tstep,sigma,taut,vxx,vyy,vzz,           &
           rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz, &
           chit,strkin,strknf,strknt,engke,engrot)

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

! update kinetic energy and stress

     Call nvt_b1_scl &
           (0,tstep,sigma,taut,vxx,vyy,vzz,           &
           rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz, &
           chit,strkin,strknf,strknt,engke,engrot)

  End If

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
     Deallocate (oxt,oyt,ozt,       Stat=fail( 6))
  End If
  Deallocate (ggx,ggy,ggz,          Stat=fail( 7))
  Deallocate (xxt,yyt,zzt,          Stat=fail( 8))
  Deallocate (vxt,vyt,vzt,          Stat=fail( 9))
  Deallocate (fxt,fyt,fzt,          Stat=fail(10))
  Deallocate (q0t,q1t,q2t,q3t,      Stat=fail(11))
  Deallocate (p0,p1,p2,p3,          Stat=fail(12))
  Deallocate (rgdxxt,rgdyyt,rgdzzt, Stat=fail(13))
  Deallocate (rgdvxt,rgdvyt,rgdvzt, Stat=fail(14))
  Deallocate (rgdoxt,rgdoyt,rgdozt, Stat=fail(15))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'npt_b1 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine npt_b1_vv
