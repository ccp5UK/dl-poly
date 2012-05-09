Subroutine nvt_g1_vv                               &
           (isw,lvar,mndis,mxdis,mxstp,temp,tstep, &
           sigma,taut,gama,chit,cint,consv,        &
           strkin,strknf,strknt,engke,engrot,      &
           imcon,mxshak,tolnce,                    &
           megcon,strcon,vircon,                   &
           megpmf,strpmf,virpmf,                   &
           strcom,vircom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian and rotational, singled
! RBs, equations of motion in molecular dynamics
! - velocity verlet with gentle ergodic thermostat -
! a Nose-Hoover thermostat chained to a Langevin thermostat -
!
! reference1: B. Leimkuhler, E. Noorizadeh, F. Theil
!             J. Stat. Phys. (2009) 135: 261–277
!
! reference2: A. Samoletov, M.A.J. Chaplain, C.P. Dettmann
!             J. Stat. Phys. (2007) 128, 1321–1336

! copyright - daresbury laboratory
! author    - i.t.todorov january 2012

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax,gsum
  Use setup_module
  Use domains_module,     Only : map
  Use site_module,        Only : ntpshl,unqshl
  Use config_module,      Only : cell,natms,nlast,nfree, &
                                 lstfre,atmnam,weight,   &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use rigid_bodies_module
  Use langevin_module,    Only : l_gst_s,r_0
  Use kinetic_module,     Only : kinstresf,kinstrest

  Implicit None

  Integer,           Intent( In    ) :: isw
  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp, &
                                        sigma,temp,taut,gama
  Real( Kind = wp ), Intent( InOut ) :: chit,cint
  Real( Kind = wp ), Intent(   Out ) :: consv
  Real( Kind = wp ), Intent( InOut ) :: tstep
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
  Logical                 :: safe,lv_up,lv_dn
  Integer,           Save :: mxkit,kit
  Integer                 :: fail(1:14),matms,i,j,i1,i2, &
                             irgd,jrgd,krgd,lrgd,rgdtyp
  Real( Kind = wp ), Save :: qmass,ceng
  Real( Kind = wp )       :: hstep,rstep,uni
  Real( Kind = wp )       :: chitdr,cintdr
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp
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
  Allocate (ggx(1:mxlrgd*Max(mxrgd,mxtrgd)),ggy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
            ggz(1:mxlrgd*Max(mxrgd,mxtrgd)),                      Stat=fail( 7))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
  Allocate (q0t(1:mxrgd),q1t(1:mxrgd),q2t(1:mxrgd),q3t(1:mxrgd),  Stat=fail(11))
  Allocate (rgdxxt(1:mxrgd),rgdyyt(1:mxrgd),rgdzzt(1:mxrgd),      Stat=fail(12))
  Allocate (rgdvxt(1:mxrgd),rgdvyt(1:mxrgd),rgdvzt(1:mxrgd),      Stat=fail(13))
  Allocate (rgdoxt(1:mxrgd),rgdoyt(1:mxrgd),rgdozt(1:mxrgd),      Stat=fail(14))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_g1 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! inertia parameter for Nose-Hoover thermostat

     qmass = 2.0_wp*sigma*taut**2
     ceng  = 2.0_wp*sigma

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

     If (megcon > 0) Call constraints_tags(imcon,lstitr,lstopt,dxx,dyy,dzz,listot)

! construct current PMF constraint vectors and shared description
! for iterative PMF constraint algorithms

     If (megpmf > 0) Call pmf_tags(imcon,lstitr,indpmf,pxx,pyy,pzz)

! generate a Gaussian random number for use in the
! Langevin process on the thermostat friction
! if not read from REVOLD

     If (l_gst_s) Then
        r_0=-6.0_wp
        Do i=1,12
           r_0=r_0+uni()
        End Do
        If (mxnode > 1) Then
           Call gsum(r_0)
           r_0=r_0/Sqrt(Real(mxnode,wp))
        End If
     End If
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

! store initial values of integration variables

     chitdr=chit
     cintdr=cint

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

! integrate and apply nvt_g1_scl thermostat - 1/2 step

     Call nvt_g1_scl &
           (hstep,ceng,qmass,temp,gama,r_0,0.0_wp,0.0_wp, &
           vxx,vyy,vzz,                                &
           rgdvxx,rgdvyy,rgdvzz,                       &
           rgdoxx,rgdoyy,rgdozz,                       &
           chit,cint,engke,engrot)

! update velocity and position of FPs

     Do j=1,nfree
        i=lstfre(j)

        If (weight(i) > 1.0e-6_wp) Then
           tmp=hstep/weight(i)
           vxx(i)=vxx(i)+tmp*fxt(i)
           vyy(i)=vyy(i)+tmp*fyt(i)
           vzz(i)=vzz(i)+tmp*fzt(i)

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

              Call pmf_shake_vv        &
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

! recover quaternion momenta

           opx=rgdoxx(irgd)*rgdrix(1,rgdtyp)
           opy=rgdoyy(irgd)*rgdriy(1,rgdtyp)
           opz=rgdozz(irgd)*rgdriz(1,rgdtyp)

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
           rgdvxx(irgd)=rgdvxx(irgd)+tmp*fmx
           rgdvyy(irgd)=rgdvyy(irgd)+tmp*fmy
           rgdvzz(irgd)=rgdvzz(irgd)+tmp*fmz

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

           chit=chitdr
           cint=cintdr

           Do i=1,matms
              vxx(i) = vxt(i)
              vyy(i) = vyt(i)
              vzz(i) = vzt(i)
           End Do

           Do irgd=1,ntrgd
              q0(irgd)=q0t(irgd)
              q1(irgd)=q1t(irgd)
              q2(irgd)=q2t(irgd)
              q3(irgd)=q3t(irgd)

              rgdvxx(irgd) = rgdvxt(irgd)
              rgdvyy(irgd) = rgdvyt(irgd)
              rgdvzz(irgd) = rgdvzt(irgd)

              rgdoxx(irgd) = rgdoxt(irgd)
              rgdoyy(irgd) = rgdoyt(irgd)
              rgdozz(irgd) = rgdozt(irgd)
           End Do

! restart vv1

           Go To 100
        End If
     End If

! second stage of velocity verlet algorithm

  Else

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
! apply velocity corrections to constraints

     If (megcon > 0 .or. megpmf > 0) Then
        Do i=1,kit
           If (megcon > 0) Call constraints_rattle &
           (mxshak,tolnce,tstep,      &
           lstopt,dxx,dyy,dzz,listot, &
           vxx,vyy,vzz)

! apply velocity corrections to PMFs

           If (megpmf > 0) Call pmf_rattle &
           (mxshak,tolnce,tstep, &
           indpmf,pxx,pyy,pzz,   &
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

! integrate and apply nvt_g1_scl thermostat - 1/2 step

     Call nvt_g1_scl &
           (hstep,ceng,qmass,temp,gama,r_0,0.0_wp,0.0_wp, &
           vxx,vyy,vzz,                                &
           rgdvxx,rgdvyy,rgdvzz,                       &
           rgdoxx,rgdoyy,rgdozz,                       &
           chit,cint,engke,engrot)

! conserved quantity less kinetic and potential energy terms

     consv = 0.5_wp*qmass*chit**2 + ceng*cint

! update kinetic energy and stress

     Call kinstresf(vxx,vyy,vzz,strknf)
     Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

     strkin=strknf+strknt

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
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_g1 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_g1_vv