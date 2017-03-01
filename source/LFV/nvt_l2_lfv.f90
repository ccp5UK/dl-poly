Subroutine nvt_l2_lfv                     &
           (lvar,mndis,mxdis,mxstp,tstep, &
           nstep,temp,chi_ep,chi_es,      &
           vel_es2,strkin,engke,          &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - leapfrog verlet with Langevin thermostat
! (standard brownian dynamics)
!
! copyright - daresbury laboratory
! authors   - i.t.todorov,s.daraszewicz,m.a.seaton & g.khara
!             february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax
  Use setup_module
  Use config_module,      Only : natms,lfrzn,weight, &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use langevin_module,    Only : fxl,fyl,fzl
  Use kinetic_module,     Only : getvom,kinstress
  Use core_shell_module,  Only : legshl
  Use constraints_module, Only : passcon
  Use pmf_module,         Only : passpmf
  Use ttm_module
  Use ttm_utils,          Only : Gep,calcchies

  Implicit None

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
  Logical                 :: safe,lv_up,lv_dn,lrand,lvel
  Integer,           Save :: mxkit
  Integer                 :: fail(1:9),kit,i,ia,ja,ka,ijk
  Real( Kind = wp )       :: rstep,chi
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             vom(1:3),sclv,sclf,gscale,      &
                             gscale1,velsq


  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

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
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 7))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 8))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail( 9))
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

  rstep = 1.0_wp/tstep
  lv_up = .false.
  lv_dn = .false.

! leapfrog verlet algorithm (starts with velocities at half step)!!!

! check whether or not Langevin forces are needed: if electron-phonon
! friction coefficient is/will be greater than zero and coupling is 
! switched on after time offset

  lrand = ((chi_ep>zero_plus .or. gvar==2) .and. l_epcp)

! Rescale chi to match average electronic temperature if
! using homogeneous electron-phonon coupling

  If (l_ttm .and. gvar==1) Call calcchies(chi_ep)

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

  If (lrand) Then
    Call langevin_forces(nstep-1,temp,tstep,chi_ep,fxl,fyl,fzl)
  Else
    fxl = 0.0_wp; fyl = 0.0_wp; fzl = 0.0_wp
  End If

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

! update velocity and position

  If (l_ttm) Then
    Do i=1,natms
       velsq = vxx(i)*vxx(i)+vyy(i)*vyy(i)+vzz(i)*vzz(i)
       lvel = (velsq>vel_es2 .and. chi_es>zero_plus)
       If (weight(i) > 1.0e-6_wp) Then
          Select Case (gvar)
          Case (0,1)
            chi = Merge(chi_ep,0.0_wp,l_epcp) + Merge(chi_es,0.0_wp,lvel)
          Case (2)
            ia = Floor((xxx(i)+zerocell(1))/delx) + 1
            ja = Floor((yyy(i)+zerocell(2))/dely) + 1
            ka = Floor((zzz(i)+zerocell(3))/delz) + 1
            ijk = 1 + ia + (ntcell(1)+2) * (ja + (ntcell(2)+2) * ka)
            chi = Merge(Gep(eltemp(ijk,0,0,0)),0.0_wp,l_epcp) + Merge(chi_es,0.0_wp,lvel)
          End Select
          tmp  = 1.0_wp/(1.0_wp+0.5_wp*chi*tstep)
          sclv = 2.0_wp*tmp-1.0_wp
          sclf = tstep*tmp
          tmp=sclf/weight(i)
          vxx(i)=vxt(i)*sclv+tmp*(fxt(i)+fxl(i))
          vyy(i)=vyt(i)*sclv+tmp*(fyt(i)+fyl(i))
          vzz(i)=vzt(i)*sclv+tmp*(fzt(i)+fzl(i))

          xxx(i)=xxt(i)+tstep*vxx(i)
          yyy(i)=yyt(i)+tstep*vyy(i)
          zzz(i)=zzt(i)+tstep*vzz(i)
       End If
    End Do
  Else
    gscale1 = Sqrt(1.0_wp+chi_es/chi_ep)
    Do i=1,natms
       velsq = vxx(i)*vxx(i)+vyy(i)*vyy(i)+vzz(i)*vzz(i)
       lvel = (velsq>vel_es2 .and. chi_es>zero_plus)
       If (weight(i) > 1.0e-6_wp) Then
          chi    = Merge(chi_ep,0.0_wp,l_epcp) + Merge(chi_es,0.0_wp,lvel)
          tmp    = 1.0_wp/(1.0_wp+0.5_wp*chi*tstep)
          sclv   = 2.0_wp*tmp-1.0_wp
          sclf   = tstep*tmp
          gscale = Merge(gscale1,1.0_wp,lvel)
          tmp=sclf/weight(i)
          vxx(i)=vxt(i)*sclv+tmp*(fxt(i)+gscale*fxl(i))
          vyy(i)=vyt(i)*sclv+tmp*(fyt(i)+gscale*fyl(i))
          vzz(i)=vzt(i)*sclv+tmp*(fzt(i)+gscale*fzl(i))

          xxx(i)=xxt(i)+tstep*vxx(i)
          yyy(i)=yyt(i)+tstep*vyy(i)
          zzz(i)=zzt(i)+tstep*vzz(i)
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

           Call constraints_shake_lfv &
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

           Call pmf_shake_lfv    &
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
           Else
              tmp = Sqrt(2.0_wp)
              tstep = 0.50_wp*tstep
           End If
           If (idnode == 0) Write(nrite,"(/,1x, &
              & 'timestep decreased, new timestep is:',3x,1p,e12.4,/)") tstep
        End If
        If (mxdr < mndis) Then
           lv_dn = .true.
           If (lv_up) Then
              tmp = Sqrt(2.0_wp/3.0_wp)
              tstep = 1.50_wp*tstep
           Else
              tmp = Sqrt(0.5_wp)
              tstep = 2.00_wp*tstep
           End If
           If (tstep > mxstp) Then
              tmp = tmp*Sqrt(tstep/mxstp)
              tstep = mxstp
           End If
           If (idnode == 0) Write(nrite,"(/,1x, &
              & 'timestep increased, new timestep is:',3x,1p,e12.4,/)") tstep
        End If
        rstep = 1.0_wp/tstep

! scale Langevin random forces

        Do i=1,natms
           fxl(i)=fxl(i)*tmp
           fyl(i)=fyl(i)*tmp
           fzl(i)=fzl(i)*tmp
        End Do

! restart lfv

        Go To 100
     End If
  End If

! full step velocity

  Do i=1,natms
     vxt(i) = 0.5_wp*(vxx(i)+vxt(i))
     vyt(i) = 0.5_wp*(vyy(i)+vyt(i))
     vzt(i) = 0.5_wp*(vzz(i)+vzt(i))
  End Do

! remove system centre of mass velocity

  Call getvom(vom,vxx,vyy,vzz)
  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        vxx(i)=vxx(i)-vom(1)
        vyy(i)=vyy(i)-vom(2)
        vzz(i)=vzz(i)-vom(3)
     End If
  End Do

  Call getvom(vom,vxt,vyt,vzt)
  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        vxt(i)=vxt(i)-vom(1)
        vyt(i)=vyt(i)-vom(2)
        vzt(i)=vzt(i)-vom(3)
     End If
  End Do

! update kinetic energy and stress at full step

  Call kinstress(vxt,vyt,vzt,strkin)
  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

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
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_l2 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_l2_lfv
