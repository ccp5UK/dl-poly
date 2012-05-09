Subroutine nvt_g0_lfv                                       &
           (lvar,mndis,mxdis,mxstp,temp,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon,        &
           megpmf,strpmf,virpmf,                            &
           sigma,taut,gama,chit,cint,consv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - leapfrog verlet with gentle ergodic thermostat -
! a Nose-Hoover thermostat chained to a Langevin thermostat
!
! reference1: B. Leimkuhler, E. Noorizadeh, F. Theil
!             J. Stat. Phys. (2009) 135: 261–277
!
! reference2: A. Samoletov, M.A.J. Chaplain, C.P. Dettmann
!             J. Stat. Phys. (2007) 128, 1321–1336
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,    Only : idnode,mxnode,gmax,gsum
  Use setup_module
  Use site_module,     Only : ntpshl,unqshl
  Use config_module,   Only : natms,atmnam,weight, &
                              xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use langevin_module, Only : r_0
  Use kinetic_module,  Only : getkin,kinstress

  Implicit None

  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp,temp
  Real( Kind = wp ), Intent( InOut ) :: tstep
  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon,strpmf(1:9),virpmf

  Real( Kind = wp ), Intent( In    ) :: sigma,taut,gama
  Real( Kind = wp ), Intent( InOut ) :: chit,cint
  Real( Kind = wp ), Intent(   Out ) :: consv


  Logical,           Save :: newjob = .true.
  Logical                 :: safe,lv_up,lv_dn
  Integer,           Save :: mxiter,mxkit
  Integer                 :: fail(1:9),iter,kit,i
  Real( Kind = wp ), Save :: qmass,ceng
  Real( Kind = wp )       :: hstep,rstep,uni
  Real( Kind = wp )       :: chit1,chit2
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp,fex


  Logical,           Allocatable :: lstitr(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)
  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

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
  End If
  Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),            Stat=fail(6))
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail(7))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail(8))
  Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(9))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_g0 allocation failure, node: ', idnode
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
  End If

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

! calculate kinetic energy

  engke=getkin(oxt,oyt,ozt)

! propagate chit set

  fex=Exp(-gama*tstep)
  chit1=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
        tstep*(2.0_wp*engke-ceng)/qmass
  chit2 = 0.5_wp*(chit+chit1)

! iterate chit and forces

  Do iter=1,mxiter

! update velocity and position using Nose-Hoover thermostating

     Do i=1,natms
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
        oxt(i)=0.5_wp*(vxt(i)+vxx(i))
        oyt(i)=0.5_wp*(vyt(i)+vyy(i))
        ozt(i)=0.5_wp*(vzt(i)+vzz(i))
     End Do

     If (iter < mxiter) Then

! calculate kinetic energy

        engke=getkin(oxt,oyt,ozt)

! propagate chit set

        chit1=fex*chit + Sqrt((1.0_wp-fex**2) * boltz*temp/qmass)*r_0 + &
              tstep*(2.0_wp*engke-ceng)/qmass
        chit2 = 0.5_wp*(chit+chit1)

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

! update chit and cint (cint is at full timestep)

  chit=chit1
  cint=cint+tstep*chit2

! conserved quantity less kinetic and potential energy terms
! (it is at full timestep)

  consv = 0.5_wp*qmass*chit2**2 + ceng*cint

! update kinetic energy and stress at full step

  Call kinstress(oxt,oyt,ozt,strkin)
  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

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
  End If
  Deallocate (oxt,oyt,ozt,         Stat=fail(6))
  Deallocate (xxt,yyt,zzt,         Stat=fail(7))
  Deallocate (vxt,vyt,vzt,         Stat=fail(8))
  Deallocate (fxt,fyt,fzt,         Stat=fail(9))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'nvt_g0 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nvt_g0_lfv