Subroutine npt_m0_lfv                     &
           (lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,cint,          &
           press,taup,chip,eta,           &
           degfre,virtot,                 &
           consv,                         &
           strkin,engke,                  &
           mxshak,tolnce,                 &
           megcon,strcon,vircon,          &
           megpmf,strpmf,virpmf,          &
           elrc,virlrc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - leapfrog verlet with Nose-Hoover thermostat and
! barostat (isotropic pressure control) and MTK coupling (symplectic)
!
! isotropic cell fluctuations
!
! reference: Martyna, Tuckerman, Tobias, Klein
!            Mol. Phys., 1996, Vol. 87 (5), p. 1117
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax
  Use setup_module
  Use site_module,        Only : ntpatm,dens,legshl
  Use config_module,      Only : cell,volm,natms,lfrzn,weight, &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use kinetic_module,     Only : getvom,getkin,kinstress
  Use constraints_module, Only : passcon
  Use pmf_module,         Only : passpmf

  Implicit None

  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
  Real( Kind = wp ), Intent( InOut ) :: tstep

  Real( Kind = wp ), Intent( In    ) :: sigma,taut
  Real( Kind = wp ), Intent( InOut ) :: chit,cint

  Real( Kind = wp ), Intent( In    ) :: press,taup
  Real( Kind = wp ), Intent( InOut ) :: chip
  Real( Kind = wp ), Intent(   Out ) :: eta(1:9)

  Integer(Kind=ip),  Intent( In    ) :: degfre
  Real( Kind = wp ), Intent( In    ) :: virtot

  Real( Kind = wp ), Intent(   Out ) :: consv

  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon, &
                                        strpmf(1:9),virpmf

  Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc


  Logical,           Save :: newjob = .true.
  Logical                 :: safe,lv_up,lv_dn
  Integer,           Save :: mxiter,mxkit
  Integer                 :: fail(1:10),iter,kit,i
  Real( Kind = wp ), Save :: cell0(1:9),volm0,elrc0,virlrc0
  Real( Kind = wp ), Save :: qmass,ceng,pmass,factor
  Real( Kind = wp )       :: hstep,rstep
  Real( Kind = wp )       :: chit1,chit2,chip1,chip2
  Real( Kind = wp )       :: vzero
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             scale,vom(1:3)


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
  End If
  Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),            Stat=fail( 6))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 7))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 8))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail( 9))
  Allocate (uxt(1:mxatms),uyt(1:mxatms),uzt(1:mxatms),            Stat=fail(10))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'npt_m0 allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! store initial values of volume, long range corrections and density

     cell0   = cell
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

! inertia parameters for Nose-Hoover thermostat and barostat

     qmass = 2.0_wp*sigma*taut**2
     tmp   = 2.0_wp*sigma / (boltz*Real(degfre,wp))
     ceng  = 2.0_wp*sigma + boltz*tmp
     pmass = (2.0_wp*sigma + 3.0_wp*boltz*tmp)*taup**2
     factor= 3.0_wp/Real(degfre,wp)

! set number of constraint+pmf shake iterations and general iteration cycles

     mxiter=7
     If (megcon > 0 .or.  megpmf > 0) Then
        mxkit=1
        mxiter=mxiter+1
     End If
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

! leapfrog verlet algorithm (starts with velocities at half step)!!!

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

! store temporary volume

  vzero=volm

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

  Do i=1,natms
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

        oxt(i)=0.5_wp*(xxt(i)+xxx(i))
        oyt(i)=0.5_wp*(yyt(i)+yyy(i))
        ozt(i)=0.5_wp*(zzt(i)+zzz(i))

     Else

        uxt(i)=0.0_wp
        uyt(i)=0.0_wp
        uzt(i)=0.0_wp

     End If
  End Do

! calculate kinetic energy

  engke=getkin(uxt,uyt,uzt)

! propagate chit and chip sets and couple
! (vircon,virpmf,chit2,chip2 are zero!!!)

  chit2=0.0_wp
  chip2=0.0_wp

  chit1 = chit + tstep*(-ceng)/qmass
  chip1 = chip + tstep*((2.0_wp*(1.0_wp+factor)*engke-virtot) - 3.0_wp*press*vzero)/pmass

  chit2 = 0.5_wp*(chit+chit1)
  chip2 = 0.5_wp*(chip+chip1)

! iterate forces, vircon, virpmf, chit and chip

  Do iter=1,mxiter

! update velocity and position using Nose-Hoover thermostatting
! and barostating in leapfrog verlet scheme

     Do i=1,natms
        If (weight(i) > 1.0e-6_wp) Then
           tmp=1.0_wp/weight(i)

           vxx(i)=vxt(i)+tstep*(tmp*fxx(i)-(chit2+(1.0_wp+factor)*chip2)*uxt(i))
           vyy(i)=vyt(i)+tstep*(tmp*fyy(i)-(chit2+(1.0_wp+factor)*chip2)*uyt(i))
           vzz(i)=vzt(i)+tstep*(tmp*fzz(i)-(chit2+(1.0_wp+factor)*chip2)*uzt(i))

           xxx(i)=xxt(i)+tstep*(vxx(i)+chip1*oxt(i))
           yyy(i)=yyt(i)+tstep*(vyy(i)+chip1*oyt(i))
           zzz(i)=zzt(i)+tstep*(vzz(i)+chip1*ozt(i))
        End If
     End Do

! SHAKE procedures

     If (megcon > 0 .or. megpmf > 0) Then
        safe=.false.
        kit =0

! update volume

        volm=vzero*Exp(3.0_wp*tstep*chip1)

! scale cell vectors - isotropic

        scale=(volm/volm0)**(1.0_wp/3.0_wp)
        cell=cell0*scale

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

              Call pmf_shake_lfv &
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

        If (iter == mxiter) Then
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

! calculate velocity at full step

     Do i=1,natms
        uxt(i)=0.5_wp*(vxt(i)+vxx(i))
        uyt(i)=0.5_wp*(vyt(i)+vyy(i))
        uzt(i)=0.5_wp*(vzt(i)+vzz(i))
     End Do

     If (iter < mxiter) Then

! calculate position at half step

        Do i=1,natms
           oxt(i)=0.5_wp*(xxt(i)+xxx(i))
           oyt(i)=0.5_wp*(yyt(i)+yyy(i))
           ozt(i)=0.5_wp*(zzt(i)+zzz(i))
        End Do

! calculate kinetic energy

        engke=getkin(uxt,uyt,uzt)

! propagate chit and chip sets and couple
! (volm=vzero, vircon,virpmf,chit2,chip2 are freshly new here!!!)

        chit1 = chit + tstep*(2.0_wp*engke+pmass*chip2**2-ceng)/qmass
        chip1 = chip*Exp(-tstep*chit2) +                                      &
                tstep*( (2.0_wp*(1.0_wp+factor)*engke-virtot-vircon-virpmf) - &
                       3.0_wp*press*vzero )/pmass

        chit2 = 0.5_wp*(chit+chit1)
        chip2 = 0.5_wp*(chip+chip1)

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

! scale tstep

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

! restore forces

        Do i=1,natms
           fxx(i) = fxt(i)
           fyy(i) = fyt(i)
           fzz(i) = fzt(i)
        End Do

! restart lfv

        Go To 100
     End If
  End If

! update volume and scale cell vectors - isotropic

  If (megcon == 0 .and. megpmf == 0) Then
     volm=vzero*Exp(3.0_wp*tstep*chip1)

     scale=(volm/volm0)**(1.0_wp/3.0_wp)
     cell=cell0*scale
  End If

! adjust long range corrections and number density

  tmp=(volm0/volm)
  elrc=elrc0*tmp
  virlrc=virlrc0*tmp
  Do i=1,ntpatm
     dens(i)=dens0(i)*tmp
  End Do

! update chit and cint (cint is at full timestep)

  chit=chit1
  cint=cint+tstep*chit2

! update chip

  chip=chip1

! construct a 'mock' scaling tensor

  Do i=2,8
     eta(i)=0.0_wp
  End Do
  eta(1)=chip
  eta(5)=chip
  eta(9)=chip

! conserved quantity less kinetic and potential energy terms
! (it is at full timestep)

  consv = 0.5_wp*qmass*chit2**2 + 0.5_wp*pmass*chip2**2 + ceng*cint + press*volm

! remove system centre of mass velocity

  Call getvom(vom,vxx,vyy,vzz)
  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        vxx(i)=vxx(i)-vom(1)
        vyy(i)=vyy(i)-vom(2)
        vzz(i)=vzz(i)-vom(3)
     End If
  End Do

  Call getvom(vom,uxt,uyt,uzt)
  Do i=1,natms
     If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
        uxt(i)=uxt(i)-vom(1)
        uyt(i)=uyt(i)-vom(2)
        uzt(i)=uzt(i)-vom(3)
     End If
  End Do

! update kinetic energy and stress at full step

  Call kinstress(uxt,uyt,uzt,strkin)
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
  End If
  Deallocate (oxt,oyt,ozt,         Stat=fail( 6))
  Deallocate (xxt,yyt,zzt,         Stat=fail( 7))
  Deallocate (vxt,vyt,vzt,         Stat=fail( 8))
  Deallocate (fxt,fyt,fzt,         Stat=fail( 9))
  Deallocate (uxt,uyt,uzt,         Stat=fail(10))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'npt_m0 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine npt_m0_lfv
