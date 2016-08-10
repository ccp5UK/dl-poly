Subroutine nvt_a1_vv                          &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,keyshl,taut,soft,       &
           strkin,strknf,strknt,engke,engrot, &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           strcom,vircom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian and rotational, singled
! RBs, equations of motion in molecular dynamics
! - velocity verlet with Andersen thermostat
! (standard brownian dynamics)
!
! Ref: Molecular dynamics at constant pressure and/or temperature,
!      H.C. Andersen. J. Chem. Phys., 72:2384-2393, 1980.
!
! (dynamics not symplectic due to the pseudo-gaussian, resampled
!  particles' momenta of a particle subset on each domain)
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsum,gmax
  Use setup_module
  Use domains_module,     Only : map
  Use site_module,        Only : dofsit
  Use config_module,      Only : imcon,cell,natms,nlast,nfree,lsite, &
                                 lsi,lsa,ltg,lfrzn,lfree,lstfre,     &
                                 weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use kinetic_module,     Only : getvom,getknr,kinstresf,kinstrest
  Use core_shell_module,  Only : ntshl,listshl,legshl,lshmv_shl,lishp_shl,lashp_shl
  Use constraints_module, Only : passcon
  Use pmf_module,         Only : passpmf

  Implicit None

  Integer,           Intent( In    ) :: isw

  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
  Real( Kind = wp ), Intent( InOut ) :: tstep

  Integer,           Intent( In    ) :: nstep,keyshl
  Real( Kind = wp ), Intent( In    ) :: temp,taut,soft

  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                        strknf(1:9),strknt(1:9),engrot

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon, &
                                        strpmf(1:9),virpmf

  Real( Kind = wp ), Intent( InOut ) :: strcom(1:9),vircom


  Logical,           Save :: newjob = .true. , &
                             unsafe = .false.
  Logical                 :: safe,lcol,lfst,lv_up,lv_dn
  Integer,           Save :: mxkit,kit
  Integer                 :: fail(1:17),i,j,k,ntp,  &
                             stp,i1,i2,local_index, &
                             matms,rtp,irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp )       :: hstep,rstep,sarurnd
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             scale,tkin,vom(1:3)
  Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                             opx,opy,opz,fmx,fmy,fmz,       &
                             tqx,tqy,tqz,trx,try,trz,       &
                             qt0,qt1,qt2,qt3,p0,p1,p2,p3,   &
                             vpx,vpy,vpz


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
  Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
  Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
  Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)

! q. index arrays and tp. sum arrays

  Integer,           Allocatable :: qn(:),tpn(:)
  Integer,           Allocatable :: qs(:,:),tps(:)
  Integer,           Allocatable :: qr(:),tpr(:)

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
  Allocate (rgdxxt(1:mxrgd),rgdyyt(1:mxrgd),rgdzzt(1:mxrgd),      Stat=fail(12))
  Allocate (rgdvxt(1:mxrgd),rgdvyt(1:mxrgd),rgdvzt(1:mxrgd),      Stat=fail(13))
  Allocate (rgdoxt(1:mxrgd),rgdoyt(1:mxrgd),rgdozt(1:mxrgd),      Stat=fail(14))
  Allocate (qn(1:mxatms),tpn(0:mxnode-1),                         Stat=fail(15))
  Allocate (qs(0:2,1:mxshl),tps(0:mxnode-1),                      Stat=fail(16))
  Allocate (qr(1:mxrgd),tpr(0:mxnode-1),                          Stat=fail(17))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_a1 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! set number of constraint+pmf shake iterations

     If (megcon > 0 .or.  megpmf > 0) mxkit=1
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

     If (megcon > 0) Call constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot)

! construct current PMF constraint vectors and shared description
! for iterative PMF constraint algorithms

     If (megpmf > 0) Call pmf_tags(lstitr,indpmf,pxx,pyy,pzz)
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

! update velocity and position of FPs

     Do j=1,nfree
        i=lstfre(j)

        If (weight(i) > 1.0e-6_wp) Then
           tmp=hstep/weight(i)
           vxx(i)=vxt(i)+tmp*fxt(i)
           vyy(i)=vyt(i)+tmp*fyt(i)
           vzz(i)=vzt(i)+tmp*fzt(i)

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

              Call constraints_shake_vv &
           (mxshak,tolnce,tstep,      &
           lstopt,dxx,dyy,dzz,listot, &
           xxx,yyy,zzz,str,vir)

! constraint virial and stress tensor

              vircon=vircon+vir
              strcon=strcon+str

              safe=.true.
           End If

           If (megpmf > 0) Then

! apply PMF correction: virpmf,strpmf - PMF constraint virial,stress

              Call pmf_shake_vv  &
           (mxshak,tolnce,tstep, &
           indpmf,pxx,pyy,pzz,   &
           xxx,yyy,zzz,str,vir)

! PMF virial and stress tensor

              virpmf=virpmf+vir
              strpmf=strpmf+str

              safe=(Abs(vir) <= zero_plus)
           End If
        End Do

        If (.not.safe) Call error(478)

! Collect per step passage statistics for bond and pmf constraints

        If (megcon > 0) Then
           passcon(3,2,1)=passcon(2,2,1)*passcon(3,2,1)
           passcon(2,2,1)=passcon(2,2,1)+1.0_wp
           passcon(3,2,1)=passcon(3,2,1)/passcon(2,2,1)+passcon(1,2,1)/passcon(2,2,1)
           passcon(4,2,1)=Min(passcon(1,2,1),passcon(4,2,1))
           passcon(5,2,1)=Max(passcon(1,2,1),passcon(5,2,1))
           passcon(1,2,1)=0.0_wp ! Reset
        End If

        If (megpmf > 0) Then
           passpmf(3,2,1)=passpmf(2,2,1)*passpmf(3,2,1)
           passpmf(2,2,1)=passpmf(2,2,1)+1.0_wp
           passpmf(3,2,1)=passpmf(3,2,1)/passpmf(2,2,1)+passpmf(1,2,1)/passpmf(2,2,1)
           passpmf(4,2,1)=Min(passpmf(1,2,1),passpmf(4,2,1))
           passpmf(5,2,1)=Max(passpmf(1,2,1),passpmf(5,2,1))
           passpmf(1,2,1)=0.0_wp ! Reset
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

! update velocity and position of RBs

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

! calculate quaternion torques

           qt0=2.0_wp*(-q1t(irgd)*trx-q2t(irgd)*try-q3t(irgd)*trz)
           qt1=2.0_wp*( q0t(irgd)*trx-q3t(irgd)*try+q2t(irgd)*trz)
           qt2=2.0_wp*( q3t(irgd)*trx+q0t(irgd)*try-q1t(irgd)*trz)
           qt3=2.0_wp*(-q2t(irgd)*trx+q1t(irgd)*try+q0t(irgd)*trz)

! recover quaternion momenta at start of time step

           opx=rgdoxt(irgd)*rgdrix(1,rgdtyp)
           opy=rgdoyt(irgd)*rgdriy(1,rgdtyp)
           opz=rgdozt(irgd)*rgdriz(1,rgdtyp)

           p0=2.0_wp*(-q1t(irgd)*opx-q2t(irgd)*opy-q3t(irgd)*opz)
           p1=2.0_wp*( q0t(irgd)*opx-q3t(irgd)*opy+q2t(irgd)*opz)
           p2=2.0_wp*( q3t(irgd)*opx+q0t(irgd)*opy-q1t(irgd)*opz)
           p3=2.0_wp*(-q2t(irgd)*opx+q1t(irgd)*opy+q0t(irgd)*opz)

! update quaternion momenta to half step

           p0=p0+hstep*qt0
           p1=p1+hstep*qt1
           p2=p2+hstep*qt2
           p3=p3+hstep*qt3

! rotate RB quaternions - update q to full timestep & amend p
! and get new rotation matrix

           Call no_squish                                             &
           (tstep,rgdrix(2,rgdtyp),rgdriy(2,rgdtyp),rgdriz(2,rgdtyp), &
           q0(irgd),q1(irgd),q2(irgd),q3(irgd),p0,p1,p2,p3)
           Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! update RB angular & COM velocities to half step

           opx=0.5_wp*(-q1(irgd)*p0+q0(irgd)*p1+q3(irgd)*p2-q2(irgd)*p3)
           opy=0.5_wp*(-q2(irgd)*p0-q3(irgd)*p1+q0(irgd)*p2+q1(irgd)*p3)
           opz=0.5_wp*(-q3(irgd)*p0+q2(irgd)*p1-q1(irgd)*p2+q0(irgd)*p3)

           rgdoxx(irgd)=opx*rgdrix(2,rgdtyp)
           rgdoyy(irgd)=opy*rgdriy(2,rgdtyp)
           rgdozz(irgd)=opz*rgdriz(2,rgdtyp)

           tmp=hstep/rgdwgt(0,rgdtyp)
           rgdvxx(irgd)=rgdvxt(irgd)+tmp*fmx
           rgdvyy(irgd)=rgdvyt(irgd)+tmp*fmy
           rgdvzz(irgd)=rgdvzt(irgd)+tmp*fmz

! update RB COM to full step

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

! check timestep for variable timestep

     If (lvar) Then

! update maximum distance a particle has travelled

        mxdr = 0.0_wp
        Do i=1,natms
           If (legshl(0,i) >= 0) &
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

! second stage of velocity verlet algorithm

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

           p0=2.0_wp*(-q1(irgd)*opx-q2(irgd)*opy-q3(irgd)*opz)
           p1=2.0_wp*( q0(irgd)*opx-q3(irgd)*opy+q2(irgd)*opz)
           p2=2.0_wp*( q3(irgd)*opx+q0(irgd)*opy-q1(irgd)*opz)
           p3=2.0_wp*(-q2(irgd)*opx+q1(irgd)*opy+q0(irgd)*opz)

! update quaternion momenta to full step

           p0=p0+hstep*qt0
           p1=p1+hstep*qt1
           p2=p2+hstep*qt2
           p3=p3+hstep*qt3

! update RB angular & COM velocities to full step

           opx=0.5_wp*(-q1(irgd)*p0+q0(irgd)*p1+q3(irgd)*p2-q2(irgd)*p3)
           opy=0.5_wp*(-q2(irgd)*p0-q3(irgd)*p1+q0(irgd)*p2+q1(irgd)*p3)
           opz=0.5_wp*(-q3(irgd)*p0+q2(irgd)*p1-q1(irgd)*p2+q0(irgd)*p3)

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

! Andersen Thermostat
!
! qualify non-shell, non-frozen free particles (n) for a random kick
! qualify RBs (r) by them and derive related shells (s)

! tpn(idnode) number of thermostatted particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen particles to thermostat

     qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
     qs(0:2,1:ntshl) = 0 ! unqualified core-shell unit with a local shell
     qr(1:ntrgd)     = 0 ! unqualified RB

     j = 0
     scale = tstep/taut
     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. legshl(0,i) >= 0) Then
           If (sarurnd(ltg(i),0,nstep) <= scale) Then
              j = j + 1
              qn(i) = 1
           End If
        End If
     End Do
     tpn(idnode) = j
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tpn(i) = 0
        End Do
        Call gsum(tpn)
     End If
     ntp = Sum(tpn)

     If (ntp == 0) Go To 200

     j = 0 ! no qualified good RB (one qualified RB is enough to trigger all)
     Do i=1,matms
        If (qn(i) == 1) Then
           If (lfree(i) == 1) j = j + 1
        End If
     End Do
     Call gsum(j)

! tpr(idnode) number of thermostatted RB units on this node (idnode)
! rtp - grand total of RB units to thermostat
! (can be larger than megrgd due to sharing)

     k = 0
     If (j > 0) Then
        If (lshmv_rgd) Then
           qn(natms+1:nlast) = 0 ! refresh the q array for shared RB units
           Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,qn)
        End If

        j = 0
        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) < lrgd) Then
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! local index of particle/site
                 If (qn(i) == 1) Then
                    If (qr(irgd) == 0) Then ! An overall hit is registered
                        qr(irgd) = 1
                        j = j + 1
                    End If

                    If (i <= natms) tpn(idnode) = tpn(idnode) - 1 ! Less free particles are hit
                 End If
              End Do

              If (qr(irgd) == 1) Then ! accounting for a random kick on the RB
                 i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
                 i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum
                 If (rgdfrz(0,rgdtyp) == 0) Then
                    If (i1 <= natms) k = k + 1
                 End If
                 If (i2 <= natms) k = k + 1
              End If
           End If
        End Do
     End If
! tpn(idnode) number of thermostatted free particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen free particles to thermostat
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tpn(i) = 0
        End Do
        Call gsum(tpn)
     End If
     ntp = Sum(tpn)
     tpr(idnode) = j
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tpr(i) = 0
        End Do
        Call gsum(tpr)
     End If
     rtp = Sum(tpr)

! Get gaussian distribution

     tkin = 0.0_wp
     mxdr = 0.0_wp
     Do i=1,natms
        If (qn(i) == 1 .and. lfree(i) == 0) Then
           If (dofsit(lsite(i)) > zero_plus) mxdr = mxdr + dofsit(lsite(i))

! Get gaussian distribution (unit variance)

           Call box_mueller_saru3(ltg(i),nstep,xxt(i),yyt(i),zzt(i))

! Get scaler to target variance/Sqrt(weight)

           tmp = 1.0_wp/Sqrt(weight(i))
           xxt(i) = xxt(i)*tmp
           yyt(i) = yyt(i)*tmp
           zzt(i) = zzt(i)*tmp

           tkin = tkin + weight(i)*(xxt(i)**2+yyt(i)**2+zzt(i)**2)
        End If
     End Do

     If (rtp > 0) Then
        Do irgd=1,ntrgd
           If (qr(irgd) == 1) Then
              rgdtyp=listrgd(0,irgd)

              lrgd=listrgd(-1,irgd)
              Do jrgd=1,lrgd
                 i=indrgd(jrgd,irgd) ! particle index
                 If (i <= natms) Then
                    If (dofsit(lsite(i)) > zero_plus) mxdr = mxdr + dofsit(lsite(i))
                 End If
              End Do

              i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
              i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

              If (rgdfrz(0,rgdtyp) == 0 .and. i1 <= natms) Then

! Get gaussian distribution (unit variance)

                 Call box_mueller_saru3(ltg(i1),nstep,xxt(i1),yyt(i1),zzt(i1))

! Get scaler to target variance/Sqrt(weight)

                 tmp = 1.0_wp/Sqrt(rgdwgt(0,rgdtyp))

                 xxt(i1) = xxt(i1)*tmp
                 yyt(i1) = yyt(i1)*tmp
                 zzt(i1) = zzt(i1)*tmp

                 tkin = tkin + rgdwgt(0,rgdtyp)*(xxt(i1)**2+yyt(i1)**2+zzt(i1)**2)

              End If

              If (i2 <= natms) Then

! Get gaussian distribution (unit variance)

                 Call box_mueller_saru3(ltg(i2),nstep,xxt(i2),yyt(i2),zzt(i2))

! Get scaler to target variance/Sqrt(weight) -
! 3 different reciprocal moments of inertia

                 xxt(i2) = xxt(i2)*Sqrt(rgdrix(2,rgdtyp))
                 yyt(i2) = yyt(i2)*Sqrt(rgdriy(2,rgdtyp))
                 zzt(i2) = zzt(i2)*Sqrt(rgdriz(2,rgdtyp))

                 tkin = tkin + (rgdrix(1,rgdtyp)*xxt(i2)**2+rgdriy(1,rgdtyp)*yyt(i2)**2+rgdriz(1,rgdtyp)*zzt(i2)**2)

              End If
           End If
        End Do
     End If

     If (mxnode > 1) Call gsum(tkin)
     If (tkin <= zero_plus) tkin = 1.0_wp
     If (mxnode > 1) Call gsum(mxdr)

! Scale to target temperature and apply thermostat

     scale = Sqrt(mxdr * boltz * temp / tkin)
     tmp = Sqrt(1.0_wp-soft**2)*scale

     j = 0
     Do i=1,natms
        If (qn(i) == 1 .and. lfree(i) == 0) Then
           If (soft <= zero_plus) Then ! New target velocity
              vxx(i) = xxt(i)*scale
              vyy(i) = yyt(i)*scale
              vzz(i) = zzt(i)*scale
           Else ! Softened velocity (mixture between old & new)
              vxx(i) = soft*vxx(i) + tmp*xxt(i)
              vyy(i) = soft*vyy(i) + tmp*yyt(i)
              vzz(i) = soft*vzz(i) + tmp*zzt(i)
           End If
        End If
     End Do

     If (rtp > 0) Then

! Update shared RBs' velocities

        If (lshmv_rgd) Call update_shared_units(natms,nlast,lsi,lsa,lishp_rgd,lashp_rgd,xxt,yyt,zzt)

! calculate new RBs' COM and angular velocities

        Do irgd=1,ntrgd
           If (qr(irgd) == 1) Then
              rgdtyp=listrgd(0,irgd)

              i1=indrgd(1,irgd) ! particle to bare the random RB COM momentum
              i2=indrgd(2,irgd) ! particle to bare the random RB angular momentum

              If (rgdfrz(0,rgdtyp) == 0) Then
                 If (soft <= zero_plus) Then ! New target velocity
                    rgdvxx(irgd) = xxt(i1)*scale
                    rgdvyy(irgd) = yyt(i1)*scale
                    rgdvzz(irgd) = zzt(i1)*scale
                 Else ! Softened velocity (mixture between old & new)
                    rgdvxx(irgd) = soft*rgdvxx(irgd) + tmp*xxt(i1)
                    rgdvyy(irgd) = soft*rgdvyy(irgd) + tmp*yyt(i1)
                    rgdvzz(irgd) = soft*rgdvzz(irgd) + tmp*zzt(i1)
                 End If
              End If

              If (soft <= zero_plus) Then ! New target velocity
                 rgdoxx(irgd) = xxt(i2)*scale
                 rgdoyy(irgd) = yyt(i2)*scale
                 rgdozz(irgd) = zzt(i2)*scale
              Else ! Softened velocity (mixture between old & new)
                 rgdoxx(irgd) = soft*rgdoxx(irgd) + tmp*xxt(i2)
                 rgdoyy(irgd) = soft*rgdoyy(irgd) + tmp*yyt(i2)
                 rgdozz(irgd) = soft*rgdozz(irgd) + tmp*zzt(i2)
              End If

! get new rotation matrix

              Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

! update RB members new velocities

              lrgd=listrgd(-1,irgd)
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
     End If

! tps(idnode) number of thermostatted core-shell units on this node (idnode)
! stp - grand total of core-shell units to thermostat

     j = 0
     If (keyshl == 1) Then
        If (lshmv_shl) Then ! refresh the q array for shared core-shell units
           qn(natms+1:nlast) = 0
           Call update_shared_units_int(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,qn)
        End If

        If (ntshl > 0) Then
           Do k=1,ntshl
              i1=local_index(listshl(1,k),matms,lsi,lsa)
              i2=local_index(listshl(2,k),matms,lsi,lsa)

              If (qn(i1) == 1 .and. i2 > 0 .and. i2 <= natms) Then
                 j = j + 1

                 qs(0,k)=1
                 qs(1,k)=i1
                 qs(2,k)=i2
              End If
           End Do
        End If
     End If
     tps(idnode) = j
     If (mxnode > 1) Then
        Do i=0,mxnode-1
           If (i /= idnode) tps(i) = 0
        End Do
        Call gsum(tps)
     End If
     stp = Sum(tps)

! Thermalise the shells on hit cores

     If (stp > 0) Then
        If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)

        If (tps(idnode) > 0) Then
           Do k=1,ntshl
              If (qs(0,k) == 1) Then
                 i1=qs(1,k)
                 i2=qs(2,k)

                 vxx(i2)=vxx(i1)
                 vyy(i2)=vyy(i1)
                 vzz(i2)=vzz(i1)
              End If
           End Do
        End If
     End If

! remove system centre of mass velocity (random momentum walk)

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

200 Continue

! update kinetic energy and stress

     Call kinstresf(vxx,vyy,vzz,strknf)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

     strkin=strknf+strknt
     engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! update rotational energy

     engrot=getknr(rgdoxx,rgdoyy,rgdozz)

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
  Deallocate (rgdxxt,rgdyyt,rgdzzt, Stat=fail(12))
  Deallocate (rgdvxt,rgdvyt,rgdvzt, Stat=fail(13))
  Deallocate (rgdoxt,rgdoyt,rgdozt, Stat=fail(14))
  Deallocate (qn,tpn,               Stat=fail(15))
  Deallocate (qs,tps,               Stat=fail(16))
  Deallocate (qr,tpr,               Stat=fail(17))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_a1 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_a1_vv
