Subroutine nvt_l2_vv                          &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,chi_ep,chi_es,vel_es2,  &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - velocity verlet with inhomogeneous
! (two-temperature) Langevin thermostat, based on Langevin impulse 
! (LI) integration
!
! Ref: Jesus A. Izaguirre, `Langevin Stabilisation of Multiscale
!      Mollified Dynamics', editors: A. Brandt, K. Binder, J. Bernholc,
!      Multiscale Computational Methods in Chemistry and Physics,
!      January 2001, vol. 117 of NATO Science Series: Series III -
!      Computer and System Sciences, pp. 34-47, IOS Press, Amsterdam
!
! (brownian dynamics is not symplectic due to the random forces)
!
! copyright - daresbury laboratory
! authors   - i.t.todorov & m.a.seaton june 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,       Only : idnode,mxnode,gmax
  Use setup_module
  Use config_module,      Only : natms,lfrzn,weight, &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use kinetic_module,     Only : getvom,kinstress
  Use core_shell_module,  Only : legshl
  Use constraints_module, Only : passcon
  Use pmf_module,         Only : passpmf
  Use ttm_module
  Use ttm_utils,          Only : Gep,calcchies,eltemp_max

  Implicit None

  Integer,           Intent( In    ) :: isw

  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
  Real( Kind = wp ), Intent( InOut ) :: tstep

  Integer,           Intent( In    ) :: nstep
  Real( Kind = wp ), Intent( In    ) :: temp,chi_es,vel_es2

  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke,chi_ep

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon, &
                                        strpmf(1:9),virpmf


  Logical,           Save :: newjob = .true.
  Logical                 :: safe,lcol,lfst,lv_up,lv_dn,lrand,lvel
  Integer,           Save :: mxkit,kit
  Integer                 :: fail(1:9),i,ia,ja,ka,ijk
  Real( Kind = wp )       :: hstep,rstep,chi
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp,vom(1:3), &
                             t0,t1,t2,scr1,scl1,scv1,                 &
                             t0a,t1a,t2a,t0b,t1b,t2b,scr1a,scl1a,     &
                             scv1a,scr1b,scl1b,scv1b,velsq,eltempmax


  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)
  Real( Kind = wp ), Allocatable :: fxr(:),fyr(:),fzr(:)
  Real( Kind = wp ), Allocatable :: fxl(:),fyl(:),fzl(:)

  fail=0
  If (megcon > 0 .or. megpmf > 0) Then
     Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
     If (megcon > 0) Then
        Allocate (lstopt(0:2,1:mxcons),listot(1:mxatms),          Stat=fail(2))
        Allocate (dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons),      Stat=fail(3))
     End If
     If (megpmf > 0) Then
        Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat=fail(4))
        Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),         Stat=fail(5))
     End If
     Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
  End If
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail(7))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail(8))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(9))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_l2 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! set number of constraint+pmf shake iterations

     If (megcon > 0 .or.  megpmf > 0) mxkit=1
     If (megcon > 0 .and. megpmf > 0) mxkit=mxshak
  End If

  If (megcon > 0 .or. megpmf > 0) Then
     lstitr(1:natms)=.false. ! initialise lstitr

! construct current bond vectors and listot array (shared
! constraint atoms) for iterative bond algorithms

     If (megcon > 0) Call constraints_tags(lstitr,lstopt,dxx,dyy,dzz,listot)

! construct current PMF constraint vectors and shared description
! for iterative PMF constraint algorithms

     If (megpmf > 0) Call pmf_tags(lstitr,indpmf,pxx,pyy,pzz)
  End If

! timestep derivatives

  hstep = 0.5_wp*tstep
  rstep = 1.0_wp/tstep
  lv_up = .false.
  lv_dn = .false.

! Rescale chi to match average electronic temperature if
! using homogeneous electron-phonon coupling

  If (l_ttm .and. gvar==1) Call calcchies(chi_ep)

! check whether or not Langevin forces are needed: if electron-phonon
! friction coefficient is/will be greater than zero and coupling is 
! switched on after time offset

  lrand = ((chi_ep>zero_plus .or. gvar==2) .and. l_epcp)

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

! Set afresh Langevin random forces

     Allocate (fxr(1:mxatms),fyr(1:mxatms),fzr(1:mxatms), Stat=fail(1))
     Allocate (fxl(1:mxatms),fyl(1:mxatms),fzl(1:mxatms), Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'nvt_l2 allocation failure+, node: ', idnode
        Call error(0)
     End If

     If (lrand) Then
       Call langevin_forces(nstep,temp,tstep,chi_ep,fxr,fyr,fzr)
       Call langevin_forces(-nstep,temp,tstep,chi_ep,fxl,fyl,fzl)
     Else
       fxr = 0.0_wp; fyr = 0.0_wp; fzr = 0.0_wp
       fxl = 0.0_wp; fyl = 0.0_wp; fzl = 0.0_wp
     End If

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

! Create primitive scalers and adjust/increase timestep if need be
! when Cholesky factorisation is compromised

     Select Case (gvar)
     Case (0,1)
       chi = Max (chi_ep, chi_ep+chi_es)
     Case (2)
       Call eltemp_max (eltempmax)
       chi = Gep(eltempmax)
       chi = Max (chi, chi+chi_es)
     End Select
     t0 = Exp(-chi*tstep)
     t1 = (1.0_wp-t0   )/(  chi)
     t2 = (1.0_wp-t0**2)/(2*chi)

     safe=.true.
     Do
        tmp=t1**2/t2
        If (tstep-tmp >= zero_plus) Then
           If ((.not.safe) .and. idnode == 0) Write(nrite,"(/,1x, &
              & 'timestep increased due to impossibility of integration, new timestep is:',3x,1p,e16.8,/)") tstep
           Exit
        Else
           safe=.false.
           tstep=tmp+1.0e-10_wp
           t0 = Exp(-chi*tstep)
           t1 = (1.0_wp-t0   )/(  chi)
           t2 = (1.0_wp-t0**2)/(2*chi)
        End If
     End Do

! Create complex scalers: only for constant and homogeneous 
! electron-phonon coupling

     Select Case (gvar)
     Case (0,1)
       chi = Merge(chi_ep,0.0_wp,l_epcp)+chi_es
       t0a = Exp(-chi*tstep)
       If (chi>zero_plus) Then
         t1a = (1.0_wp-t0a   )/(  chi)
         t2a = (1.0_wp-t0a**2)/(2*chi)
         scr1a = (t1a-t2a)/Sqrt(t2a*tstep)/chi
         scl1a = Sqrt(1.0_wp-(t1a**2)/(t2a*tstep))/chi
         scv1a = Sqrt(t2a/tstep)
       Else
         t1a = tstep
         t2a = 0.0_wp
         scr1a = 0.0_wp
         scl1a = 0.0_wp
         scv1a = 0.0_wp
       End If
       If (lrand) Then
         t0b = Exp(-chi_ep*tstep)
         t1b = (1.0_wp-t0b)/chi_ep
         t2b = (1.0_wp-t0b**2)/(2*chi_ep)
         scr1b = (t1b-t2b)/Sqrt(t2b*tstep)/chi_ep
         scl1b = Sqrt(1.0_wp-(t1b**2)/(t2b*tstep))/chi_ep
         scv1b = Sqrt(t2b/tstep)
       Else
         t0b = 0.0_wp
         t1b = tstep
         t2b = 0.0_wp
         scr1b = 0.0_wp
         scl1b = 0.0_wp
         scv1b = 0.0_wp
       End If
     End Select

! update velocity and position

     If (l_ttm) Then

       If (oneway) Then
       ! one-way electron-phonon coupling
         Do i=1,natms
            velsq = vxx(i)*vxx(i)+vyy(i)*vyy(i)+vzz(i)*vzz(i)
            lvel = (velsq>vel_es2 .and. chi_es>zero_plus)
            If (weight(i) > 1.0e-6_wp) Then
               ! check for active cell and electronic temperature is
               ! higher than ionic tmeperature: if not, switch off thermostat
               ia = Floor((xxx(i)+zerocell(1))/delx) + 1
               ja = Floor((yyy(i)+zerocell(2))/dely) + 1
               ka = Floor((zzz(i)+zerocell(3))/delz) + 1
               ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
               If (act_ele_cell(ijk,0,0,0)>zero_plus .and. eltemp(ijk,0,0,0)>tempion(ijk)) Then
                 Select Case (gvar)
                 Case (0,1)
                   t0 = Merge(t0a,t0b,lvel)
                   t1 = Merge(t1a,t1b,lvel)
                   scr1 = Merge(scr1a,scr1b,lvel)
                   scl1 = Merge(scl1a,scl1b,lvel)
                   scv1 = Merge(scv1a,scv1b,lvel)
                 Case (2)
                   chi = Merge(Gep(eltemp(ijk,0,0,0)),0.0_wp,l_epcp) + Merge(chi_es,0.0_wp,lvel)
                   If (l_epcp) Then
                     t0 = Exp(-tstep*chi)
                     t1 = (1.0_wp-t0)/chi
                     t2 = (1.0_wp-t0**2)/(2*chi)
                     scr1 = (t1-t2)/Sqrt(t2*tstep)/chi
                     scl1 = Sqrt(1.0_wp-(t1**2)/(t2*tstep))/chi
                     scv1 = Sqrt(t2/tstep)
                   Else
                     t0 = 1.0_wp
                     t1 = tstep
                     scr1 = 0.0_wp
                     scl1 = 0.0_wp
                     scv1 = 0.0_wp
                   End If
                 End Select
               Else
                 t0 = 1.0_wp
                 t1 = tstep
                 scr1 = 0.0_wp
                 scl1 = 0.0_wp
                 scv1 = 0.0_wp
               End If

! Half-kick velocity

               tmp=hstep/weight(i)
               vxx(i)=vxt(i)+tmp*fxt(i)
               vyy(i)=vyt(i)+tmp*fyt(i)
               vzz(i)=vzt(i)+tmp*fzt(i)

               tmp=tstep/weight(i)

! Full time fluctuations on positions using half-kick velocity

               xxx(i)=xxt(i)+vxx(i)*t1+tmp*(fxr(i)*scr1+fxl(i)*scl1)
               yyy(i)=yyt(i)+vyy(i)*t1+tmp*(fyr(i)*scr1+fyl(i)*scl1)
               zzz(i)=zzt(i)+vzz(i)*t1+tmp*(fzr(i)*scr1+fzl(i)*scl1)

! Full time fluctuations on half-kick perculiar (thermal only) velocity

               vxx(i)=(vxx(i)-ttmvom(ijk,1))*t0+tmp*fxr(i)*scv1+ttmvom(ijk,1)
               vyy(i)=(vyy(i)-ttmvom(ijk,2))*t0+tmp*fyr(i)*scv1+ttmvom(ijk,2)
               vzz(i)=(vzz(i)-ttmvom(ijk,3))*t0+tmp*fzr(i)*scv1+ttmvom(ijk,3)

            End If
         End Do

       Else

         Do i=1,natms
            velsq = vxx(i)*vxx(i)+vyy(i)*vyy(i)+vzz(i)*vzz(i)
            lvel = (velsq>vel_es2 .and. chi_es>zero_plus)
            If (weight(i) > 1.0e-6_wp) Then
               ! check for active cell: if not, switch off thermostat
               ia = Floor((xxx(i)+zerocell(1))/delx) + 1
               ja = Floor((yyy(i)+zerocell(2))/dely) + 1
               ka = Floor((zzz(i)+zerocell(3))/delz) + 1
               ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
               If (act_ele_cell(ijk,0,0,0)>zero_plus) Then
                 Select Case (gvar)
                 Case (0,1)
                   t0 = Merge(t0a,t0b,lvel)
                   t1 = Merge(t1a,t1b,lvel)
                   scr1 = Merge(scr1a,scr1b,lvel)
                   scl1 = Merge(scl1a,scl1b,lvel)
                   scv1 = Merge(scv1a,scv1b,lvel)
                 Case (2)
                   chi = Merge(Gep(eltemp(ijk,0,0,0)),0.0_wp,l_epcp) + Merge(chi_es,0.0_wp,lvel)
                   If (l_epcp) Then
                     t0 = Exp(-tstep*chi)
                     t1 = (1.0_wp-t0)/chi
                     t2 = (1.0_wp-t0**2)/(2*chi)
                     scr1 = (t1-t2)/Sqrt(t2*tstep)/chi
                     scl1 = Sqrt(1.0_wp-(t1**2)/(t2*tstep))/chi
                     scv1 = Sqrt(t2/tstep)
                   Else
                     t0 = 1.0_wp
                     t1 = tstep
                     t2 = 0.0_wp
                     scr1 = 0.0_wp
                     scl1 = 0.0_wp
                     scv1 = 0.0_wp
                   End If
                 End Select
               Else
                 t0 = 1.0_wp
                 t1 = tstep
                 t2 = 0.0_wp
                 scr1 = 0.0_wp
                 scl1 = 0.0_wp
                 scv1 = 0.0_wp
               End If

! Half-kick velocity

               tmp=hstep/weight(i)
               vxx(i)=vxt(i)+tmp*fxt(i)
               vyy(i)=vyt(i)+tmp*fyt(i)
               vzz(i)=vzt(i)+tmp*fzt(i)

               tmp=tstep/weight(i)

! Full time fluctuations on positions using half-kick velocity

               xxx(i)=xxt(i)+vxx(i)*t1+tmp*(fxr(i)*scr1+fxl(i)*scl1)
               yyy(i)=yyt(i)+vyy(i)*t1+tmp*(fyr(i)*scr1+fyl(i)*scl1)
               zzz(i)=zzt(i)+vzz(i)*t1+tmp*(fzr(i)*scr1+fzl(i)*scl1)

! Full time fluctuations on half-kick peculiar (thermal only) velocity

               vxx(i)=(vxx(i)-ttmvom(ijk,1))*t0+tmp*fxr(i)*scv1+ttmvom(ijk,1)
               vyy(i)=(vyy(i)-ttmvom(ijk,2))*t0+tmp*fyr(i)*scv1+ttmvom(ijk,2)
               vzz(i)=(vzz(i)-ttmvom(ijk,3))*t0+tmp*fzr(i)*scv1+ttmvom(ijk,3)

            End If
         End Do
       End If

     Else

! no ttm option: just inhomogeneous Langevin thermostat

       Do i=1,natms
          velsq = vxx(i)*vxx(i)+vyy(i)*vyy(i)+vzz(i)*vzz(i)
          lvel = (velsq>vel_es2 .and. chi_es>zero_plus)
          If (weight(i) > 1.0e-6_wp) Then

! Half-kick velocity

             tmp=hstep/weight(i)
             vxx(i)=vxt(i)+tmp*fxt(i)
             vyy(i)=vyt(i)+tmp*fyt(i)
             vzz(i)=vzt(i)+tmp*fzt(i)

             t0 = Merge(t0a,t0b,lvel)
             t1 = Merge(t1a,t1b,lvel)
             scr1 = Merge(scr1a,scr1b,lvel)
             scl1 = Merge(scl1a,scl1b,lvel)
             scv1 = Merge(scv1a,scv1b,lvel)
             tmp=tstep/weight(i)*(1.0_wp+Merge(chi_es/chi_ep,0.0_wp,lvel))

! Full time fluctuations on positions using half-kick velocity
! (time multipler adjusted to increase random forces when
! electron stopping is required due to increase in chi)

             xxx(i)=xxt(i)+vxx(i)*t1+tmp*(fxr(i)*scr1+fxl(i)*scl1)
             yyy(i)=yyt(i)+vyy(i)*t1+tmp*(fyr(i)*scr1+fyl(i)*scl1)
             zzz(i)=zzt(i)+vzz(i)*t1+tmp*(fzr(i)*scr1+fzl(i)*scl1)

! Full time fluctuations on half-kick velocity

             vxx(i)=vxx(i)*t0+tmp*fxr(i)*scv1
             vyy(i)=vyy(i)*t0+tmp*fyr(i)*scv1
             vzz(i)=vzz(i)*t0+tmp*fzr(i)*scv1

          End If
       End Do

     End If

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
           If (legshl(0,i) >= 0) &
              mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
        End Do
        mxdr=Sqrt(mxdr)
        If (mxnode > 1) Call gmax(mxdr)

        If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

! scale tstep and derivatives, and get scaler for Langevin random forces

           If (mxdr > mxdis) Then
              lv_up = .true.
              If (lv_dn) Then
                 tmp = Sqrt(4.0_wp/3.0_wp)
                 tstep = 0.75_wp*tstep
                 hstep = 0.50_wp*tstep
              Else
                 tmp = Sqrt(2.0_wp)
                 tstep = hstep
                 hstep = 0.50_wp*tstep
              End If
              If (idnode == 0) Write(nrite,"(/,1x, &
                 & 'timestep decreased, new timestep is:',3x,1p,e12.4,/)") tstep
           End If
           If (mxdr < mndis) Then
              lv_dn = .true.
              If (lv_up) Then
                 tmp = Sqrt(2.0_wp/3.0_wp)
                 tstep = 1.50_wp*tstep
                 hstep = 0.50_wp*tstep
              Else
                 tmp = Sqrt(0.5_wp)
                 hstep = tstep
                 tstep = 2.00_wp*tstep
              End If
              If (tstep > mxstp) Then
                 tmp = tmp*Sqrt(tstep/mxstp)
                 tstep = mxstp
                 hstep = 0.50_wp*tstep
              End If
              If (idnode == 0) Write(nrite,"(/,1x, &
                 & 'timestep increased, new timestep is:',3x,1p,e12.4,/)") tstep
           End If
           rstep = 1.0_wp/tstep

! scale Langevin random forces

           Do i=1,natms
              fxr(i) = fxr(i)*tmp
              fyr(i) = fyr(i)*tmp
              fzr(i) = fzr(i)*tmp

              fxl(i) = fxl(i)*tmp
              fyl(i) = fyl(i)*tmp
              fzl(i) = fzl(i)*tmp
           End Do

! restart vv1

           Go To 100
        End If
     End If

     Deallocate (fxr,fyr,fzr, Stat=fail(1))
     Deallocate (fxl,fyl,fzl, Stat=fail(2))
     If (Any(fail > 0)) Then
        Write(nrite,'(/,1x,a,i0)') 'nvt_l2 deallocation failure+, node: ', idnode
        Call error(0)
     End If

! second stage of velocity verlet algorithm

  Else

! update velocity (another half-kick)

     Do i=1,natms
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

! remove system centre of mass velocity

     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i)=vxx(i)-vom(1)
           vyy(i)=vyy(i)-vom(2)
           vzz(i)=vzz(i)-vom(3)
        End If
     End Do

! update kinetic energy and stress

     Call kinstress(vxx,vyy,vzz,strkin)
     engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  End If

  If (megcon > 0 .or. megpmf > 0) Then
     Deallocate (lstitr,           Stat=fail(1))
     If (megcon > 0) Then
        Deallocate (lstopt,listot, Stat=fail(2))
        Deallocate (dxx,dyy,dzz,   Stat=fail(3))
     End If
     If (megpmf > 0) Then
        Deallocate (indpmf,        Stat=fail(4))
        Deallocate (pxx,pyy,pzz,   Stat=fail(5))
     End If
     Deallocate (oxt,oyt,ozt,      Stat=fail(6))
  End If
  Deallocate (xxt,yyt,zzt,         Stat=fail(7))
  Deallocate (vxt,vyt,vzt,         Stat=fail(8))
  Deallocate (fxt,fyt,fzt,         Stat=fail(9))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_l2 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_l2_vv
