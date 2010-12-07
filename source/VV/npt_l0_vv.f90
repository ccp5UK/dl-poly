Subroutine npt_l0_vv                                 &
           (isw,lvar,mndis,mxdis,tstep,strkin,engke, &
           imcon,mxshak,tolnce,megcon,strcon,vircon, &
           megpmf,strpmf,virpmf,                     &
           degfre,sigma,chi,consv,                   &
           press,tai,chip,eta,virtot,                &
           elrc,virlrc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - velocity verlet with Langevin thermostat and
! barostat (isotropic pressure control) (symplectic)
!
! isotropic cell fluctuations
!
! reference: D. Quigley and M.I.J. Probert
!            J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsum,gmax
  Use setup_module
  Use site_module,        Only : ntpatm,dens
  Use config_module,      Only : cell,volm,natms,ltg,lfrzn,weight, &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use langevin_module,    Only : l_lan_s,fxl,fyl,fzl,fpl
  Use core_shell_module,  Only : ntshl,listshl
  Use kinetic_module,     Only : getvom,getkin,kinstress

  Implicit None

  Integer,           Intent( In    ) :: isw
  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis
  Real( Kind = wp ), Intent( InOut ) :: tstep
  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke

  Integer,           Intent( In    ) :: imcon,mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon,strpmf(1:9),virpmf

  Integer(Kind=ip),  Intent( In    ) :: degfre
  Real( Kind = wp ), Intent( In    ) :: sigma,chi
  Real( Kind = wp ), Intent(   Out ) :: consv
  Real( Kind = wp ), Intent( In    ) :: press,tai
  Real( Kind = wp ), Intent( InOut ) :: chip
  Real( Kind = wp ), Intent(   Out ) :: eta(1:9)
  Real( Kind = wp ), Intent( In    ) :: virtot
  Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc


  Logical,           Save :: newjob = .true.
  Logical                 :: safe,lv_up,lv_dn
  Integer,           Save :: mxiter,mxkit,kit
  Integer                 :: fail(1:9),iter,i
  Real( Kind = wp ), Save :: cell0(1:9),volm0,elrc0,virlrc0
  Real( Kind = wp ), Save :: temp,pmass
  Real( Kind = wp )       :: hstep,qstep,rstep,uni
  Real( Kind = wp )       :: chip0,engke0
  Real( Kind = wp )       :: vzero
  Real( Kind = wp )       :: xt,yt,zt,vir,vir1,str(1:9),mxdr,tmp, &
                             scale,vom(1:3)


  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

  Integer,           Allocatable :: lstopt(:,:),listot(:)
  Real( Kind = wp ), Allocatable :: dxx(:),dyy(:),dzz(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)

  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

  Real( Kind = wp ), Allocatable, Save :: dens0(:)

  fail=0
  If (megcon > 0 .or. megpmf > 0) Then
     Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
     If (megcon > 0) Then
        Allocate (lstopt(1:2,1:mxcons),listot(1:mxatms),          Stat=fail(2))
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
     Write(nrite,'(/,1x,a,i0)') 'npt_l0 allocation failure, node: ', idnode
     Call error(0)
  End If


! timestep derivatives

  hstep = 0.5_wp*tstep
  qstep = 0.5_wp*hstep
  rstep = 1.0_wp/tstep
  lv_up = .false.
  lv_dn = .false.

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

! inertia parameter for barostat

     temp  = 2.0_wp*sigma / (boltz*Real(degfre,wp))
     pmass = (2.0_wp*sigma + 3.0_wp*boltz*temp)*(2.0_wp*pi/tai)**2

! set number of constraint+pmf shake iterations and general iteration cycles

     mxiter=1
     If (megcon > 0 .or.  megpmf > 0) Then
        mxkit=1
        mxiter=mxiter+3
     End If
     If (megcon > 0 .and. megpmf > 0) mxkit=mxshak

! Generate Langevin forces for particles and
! Langevin pseudo-tensor force for barostat piston
! if not read from REVOLD

     If (l_lan_s) Then
        Call langevin_forces(temp,tstep,chi,fxl,fyl,fzl)

        fpl=0.0_wp
        tmp=-6.0_wp
        Do i=1,12
           tmp=tmp+uni()
        End Do
        If (mxnode > 1) Then
           Call gsum(tmp)
           tmp=tmp/Sqrt(Real(mxnode,wp))
        End If
        tmp=tmp*Sqrt(2.0_wp*tai*boltz*temp*pmass*rstep)/3.0_wp
        fpl(1)=tmp
        fpl(5)=tmp
        fpl(9)=tmp
     Else
        tmp=(fpl(1)+fpl(5)+fpl(9))/3.0_wp
        fpl=0.0_wp
     End If
     fpl(1)=tmp
     fpl(5)=tmp
     fpl(9)=tmp
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

! store current integration variables

     vzero=volm
     chip0=chip
     engke0=engke

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

     Do iter=1,mxiter

! integrate and apply Langevin thermostat - 1/4 step

        scale=Exp(-qstep*chi)
        Do i=1,natms
           vxx(i)=scale*vxx(i)
           vyy(i)=scale*vyy(i)
           vzz(i)=scale*vzz(i)
        End Do
        engke=engke*scale**2

! constraint+pmf virial and stress

        vir=vircon+virpmf
        str=strcon+strpmf

! integrate and apply npt_h0_scl barostat - 1/2 step
! augment vir to include the random force on the barostat

        vir1=vir-3.0_wp*fpl(1)
        Call npt_h0_scl &
           (1,hstep,degfre,pmass,tai,volm,press,vir1,virtot, &
           vxx,vyy,vzz,chip,engke)

! integrate and apply Langevin thermostat - 1/4 step

        scale=Exp(-qstep*chi)
        Do i=1,natms
           vxx(i)=scale*vxx(i)
           vyy(i)=scale*vyy(i)
           vzz(i)=scale*vzz(i)
        End Do
        engke=engke*scale**2

! update velocities

        Do i=1,natms
           If (weight(i) > 1.0e-6_wp) Then
              tmp=hstep/weight(i)
              vxx(i)=vxx(i)+tmp*(fxx(i)+fxl(i))
              vyy(i)=vyy(i)+tmp*(fyy(i)+fyl(i))
              vzz(i)=vzz(i)+tmp*(fzz(i)+fzl(i))
           End If
        End Do

! update volume

        volm=volm*Exp(3.0_wp*tstep*chip)

! scale cell vectors - isotropic

        scale=(volm/volm0)**(1.0_wp/3.0_wp)
        cell=cell0*scale

! update positions

        scale=Exp(tstep*chip)
        Do i=1,natms
           If (weight(i) > 1.0e-6_wp) Then
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

! restore original integration parameters as well as
! velocities if iter < mxiter
! in the next iteration vircon and virpmf are freshly new

        If (iter < mxiter) Then
           volm=vzero
           chip=chip0
           engke=engke0

           Do i=1,natms
              vxx(i) = vxt(i)
              vyy(i) = vyt(i)
              vzz(i) = vzt(i)
           End Do
        End If
     End Do

! check timestep for variable timestep

     If (lvar) Then

! update maximum distance a particle has travelled

        mxdr = 0.0_wp
        Do i=1,natms
           If (All(listshl(2,1:ntshl) /= ltg(i))) &
              mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
        End Do
        mxdr=Sqrt(mxdr)
        If (mxnode > 1) Call gmax(mxdr)

        If (mxdr < mndis .or. mxdr > mxdis) Then

! scale tstep and derivatives, and get scaler for Langevin random forces

           If (mxdr > mxdis) Then
              lv_up = .true.
              If (lv_dn) Then
                 tstep = 0.75_wp*tstep
                 hstep = 0.50_wp*tstep
                 qstep = 0.50_wp*hstep
                 tmp = Sqrt(4.0_wp/3.0_wp)
              Else
                 tstep = hstep
                 hstep = 0.50_wp*tstep
                 qstep = 0.50_wp*hstep
                 tmp = Sqrt(2.0_wp)
              End If
              If (idnode == 0) Write(nrite,"(/,1x, &
                 & 'timestep decreased, new timestep is:',3x,1p,e12.4,/)") tstep
           End If
           If (mxdr < mndis) Then
              lv_dn = .true.
              If (lv_up) Then
                 tstep = 1.50_wp*tstep
                 hstep = 0.50_wp*tstep
                 qstep = 0.50_wp*hstep
                 tmp = Sqrt(2.0_wp/3.0_wp)
              Else
                 qstep = hstep
                 hstep = tstep
                 tstep = 2.00_wp*tstep
                 tmp = Sqrt(0.5_wp)
              End If
              If (idnode == 0) Write(nrite,"(/,1x, &
                 & 'timestep increased, new timestep is:',3x,1p,e12.4,/)") tstep
           End If
           rstep = 1.0_wp/tstep

! scale Langevin random forces

           Do i=1,natms
              fxl(i) = fxl(i)*tmp
              fyl(i) = fyl(i)*tmp
              fzl(i) = fzl(i)*tmp
           End Do
           fpl = fpl*tmp

! restore initial conditions

           volm=vzero
           chip=chip0
           engke=engke0

           Do i=1,natms
              vxx(i) = vxt(i)
              vyy(i) = vyt(i)
              vzz(i) = vzt(i)

              fxx(i) = fxt(i)
              fyy(i) = fyt(i)
              fzz(i) = fzt(i)
           End Do

! restart vv1

           Go To 100
        End If
     End If

! adjust long range corrections and number density

     tmp=(volm0/volm)
     elrc=elrc0*tmp
     virlrc=virlrc0*tmp
     Do i=1,ntpatm
        dens(i)=dens0(i)*tmp
     End Do

! second stage of velocity verlet algorithm

  Else

! Generate Langevin forces for particles and
! Langevin pseudo-tensor force for barostat piston

     Call langevin_forces(temp,tstep,chi,fxl,fyl,fzl)

     fpl=0.0_wp
     tmp=-6.0_wp
     Do i=1,12
        tmp=tmp+uni()
     End Do
     If (mxnode > 1) Then
        Call gsum(tmp)
        tmp=tmp/Sqrt(Real(mxnode,wp))
     End If
     tmp=tmp*Sqrt(2.0_wp*tai*boltz*temp*pmass*rstep)/3.0_wp
     fpl(1)=tmp
     fpl(5)=tmp
     fpl(9)=tmp

! update velocity

     Do i=1,natms
        If (weight(i) > 1.0e-6_wp) Then
           tmp=hstep/weight(i)
           vxx(i)=vxx(i)+tmp*(fxx(i)+fxl(i))
           vyy(i)=vyy(i)+tmp*(fyy(i)+fyl(i))
           vzz(i)=vzz(i)+tmp*(fzz(i)+fzl(i))
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

! update kinetic energy

     engke=getkin(vxx,vyy,vzz)

! integrate and apply Langevin thermostat - 1/4 step

     scale=Exp(-qstep*chi)
     Do i=1,natms
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
     End Do
     engke=engke*scale**2

! constraint+pmf virial and stress

     vir=vircon+virpmf
     str=strcon+strpmf

! integrate and apply npt_h0_scl barostat - 1/2 step
! augment vir to include the random force on the barostat

     vir1=vir-3.0_wp*fpl(1)
     Call npt_h0_scl &
           (1,hstep,degfre,pmass,tai,volm,press,vir1,virtot, &
           vxx,vyy,vzz,chip,engke)

! integrate and apply Langevin thermostat - 1/4 step

     scale=Exp(-qstep*chi)
     Do i=1,natms
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
     End Do
     engke=engke*scale**2

! conserved quantity less kinetic and potential energy terms

     consv = 0.5_wp*pmass*chip**2 + press*volm

! remove system centre of mass velocity

     Call getvom(vom,vxx,vyy,vzz)

     Do i=1,natms
        If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
           vxx(i)=vxx(i)-vom(1)
           vyy(i)=vyy(i)-vom(2)
           vzz(i)=vzz(i)-vom(3)
        End If
     End Do

! construct a 'mock' scaling tensor

     Do i=2,8
        eta(i)=0.0_wp
     End Do
     eta(1)=chip
     eta(5)=chip
     eta(9)=chip

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
     Write(nrite,'(/,1x,a,i0)') 'npt_l0 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine npt_l0_vv
