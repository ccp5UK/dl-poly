Subroutine nvt_a0_vv                                 &
           (isw,lvar,mndis,mxdis,temp,tstep,         &
           keyshl,taut,soft,strkin,engke,            &
           imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - velocity verlet with Andersen thermostat
! (standard brownian dynamics)
!
! Ref: Molecular dynamics at constant pressure and/or temperature,
!      H.C. Andersen. J. Chem. Phys., 72:2384-2393, 1980.
!
! (dynamics not symplectic due to the pseudo-gaussian, resampled
!  particles' momenta of a particle subset on each domain)
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsum,gmax
  Use setup_module
  Use site_module,       Only : dofsit,ntpshl,unqshl
  Use config_module,     Only : natms,nlast,lsite,lsi,lsa,lfrzn,atmnam, &
                                weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use core_shell_module, Only : ntshl,listshl,lshmv_shl,lishp_shl,lashp_shl
  Use kinetic_module,    Only : getvom,kinstress

  Implicit None

  Integer,           Intent( In    ) :: isw,keyshl
  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,temp,taut,soft
  Real( Kind = wp ), Intent( InOut ) :: tstep
  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon,strpmf(1:9),virpmf


  Logical,           Save :: newjob = .true.
  Logical                 :: safe,lv_up,lv_dn
  Integer,           Save :: mxkit,kit
  Integer                 :: fail(1:11),i,j,k,ntp,  &
                             stp,i1,i2,local_index, &
                             matms
  Real( Kind = wp )       :: hstep,rstep,uni
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             scale,tkin,vom(1:3)


  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

! q. index arrays and tp. sum arrays

  Integer,           Allocatable :: qn(:),tpn(:)
  Integer,           Allocatable :: qs(:,:),tps(:)

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
     Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail( 6))
  End If
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 7))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 8))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail( 9))
  Allocate (qn(1:mxatms),tpn(0:mxnode-1),                         Stat=fail(10))
  Allocate (qs(0:2,1:mxshl),tps(0:mxnode-1),                      Stat=fail(11))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_a0 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! set number of constraint+pmf shake iterations

     If (megcon > 0 .or.  megpmf > 0) mxkit=1
     If (megcon > 0 .and. megpmf > 0) mxkit=mxshak
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

! timestep derivatives

  hstep = 0.5_wp*tstep
  rstep = 1.0_wp/tstep
  lv_up = .false.
  lv_dn = .false.

! first pass of velocity verlet algorithm

  If (isw == 0) Then

! store initial values

     Do i=1,natms
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

! update velocity and position

     Do i=1,natms
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

        Do i=1,natms
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

        Do i=1,natms
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

! restart vv1

           Go To 100
        End If
     End If

! second stage of velocity verlet algorithm

  Else

! update velocity

     Do i=1,natms
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

! Andersen Thermostat
!
! Here we become node-dependent (using of uni and gauss) - i.e.
! pseudo-randomness depends on the DD mapping which depends on
! the number of nodes and system size
!
! qualify non-shell, non-frozen particles (n) for a random kick
! derive related shells (s)

! tpn(idnode) number of thermostated particles on this node (idnode)
! ntp - grand total of non-shell, non-frozen particles to thermostat

     qn(1:natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
     qs(0:2,1:ntshl) = 0 ! unqualified core-shell unit with a local shell

     j = 0
     tmp = tstep/taut
     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. (.not.Any(unqshl(1:ntpshl) == atmnam(i)))) Then
           If (uni() <= tmp) Then
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

! tps(idnode) number of thermostated core-shell units on this node (idnode)
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

! Get gaussian distribution (unit variance)

     Call gauss(tpn(idnode),xxt,yyt,zzt)

     tkin = 0.0_wp
     mxdr = 0.0_wp
     j = 0
     Do i=1,natms
        If (qn(i) == 1) Then
           If (dofsit(lsite(i)) > zero_plus) mxdr = mxdr + dofsit(lsite(i))

           j = j + 1

! Get scaler to target variance/Sqrt(weight)

           tmp = 1.0_wp/Sqrt(weight(i))
           xxt(j) = xxt(j)*tmp
           yyt(j) = yyt(j)*tmp
           zzt(j) = zzt(j)*tmp

           tkin = tkin + weight(i)*(xxt(j)**2+yyt(j)**2+zzt(j)**2)
        End If
     End Do

     If (mxnode > 1) Call gsum(tkin)
     If (tkin <= zero_plus) tkin = 1.0_wp
     If (mxnode > 1) Call gsum(mxdr)

! Scale to target tempearture and apply thermostat

     scale = Sqrt(mxdr * boltz * temp / tkin)
     tmp = Sqrt(1.0_wp-soft**2)*scale

     j = 0
     Do i=1,natms
        If (qn(i) == 1) Then
           j = j + 1

           If (soft <= zero_plus) Then ! New target velocity
              vxx(i) = xxt(j)*scale
              vyy(i) = yyt(j)*scale
              vzz(i) = zzt(j)*scale
           Else ! Softened velocity (mixture between old & new)
              vxx(i) = soft*vxx(i) + tmp*xxt(j)
              vyy(i) = soft*vyy(i) + tmp*yyt(j)
              vzz(i) = soft*vzz(i) + tmp*zzt(j)
           End If
        End If
     End Do

! Thermalise the shells on hit cores

     If (stp > 0) Then
        If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,vxx,vyy,vzz)

        If (tps(idnode) > 0) Then
           j = 0
           Do k=1,ntshl
              If (qs(0,k) == 1) Then
                 j = j + 1

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

     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i)=vxx(i)-vom(1)
           vyy(i)=vyy(i)-vom(2)
           vzz(i)=vzz(i)-vom(3)
        End If
     End Do

200  Continue

! update kinetic energy and stress

     Call kinstress(vxx,vyy,vzz,strkin)
     engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  End If

  If (megcon > 0 .or. megpmf > 0) Then
     Deallocate (lstitr,           Stat=fail( 1))
     If (megcon > 0) Then
        Deallocate (lstopt,listot, Stat=fail( 2))
        Deallocate (dxx,dyy,dzz,   Stat=fail( 3))
     End If
     If (megpmf > 0) Then
        Deallocate (indpmf,        Stat=fail( 4))
        Deallocate (pxx,pyy,pzz,   Stat=fail( 5))
     End If
     Deallocate (oxt,oyt,ozt,      Stat=fail( 6))
  End If
  Deallocate (xxt,yyt,zzt,         Stat=fail( 7))
  Deallocate (vxt,vyt,vzt,         Stat=fail( 8))
  Deallocate (fxt,fyt,fzt,         Stat=fail( 9))
  Deallocate (qn,tpn,              Stat=fail(10))
  Deallocate (qs,tps,              Stat=fail(11))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_a0 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_a0_vv
