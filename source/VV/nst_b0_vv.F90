Subroutine nst_b0_vv                          &
           (isw,lvar,mndis,mxdis,mxstp,tstep, &
           sigma,taut,chit,                   &
           press,strext,taup,chip,eta,        &
           iso,ten,stress,                    &
           strkin,engke,                      &
           mxshak,tolnce,                     &
           megcon,strcon,vircon,              &
           megpmf,strpmf,virpmf,              &
           elrc,virlrc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for integrating newtonian equations of motion in
! molecular dynamics - velocity verlet with Berendsen thermostat and
! anisotropic pressure control (not symplectic)
!
! isothermal compressibility (beta) set to that of liquid water
! = 0.007372 dl_poly units
!
! iso=0 fully anisotropic barostat
! iso=1 semi-isotropic barostat to constant normal pressure & surface area
! iso=2 semi-isotropic barostat to constant normal pressure & surface tension
!                               or with orthorhombic constraints (ten=0.0_wp)
! iso=3 semi-isotropic barostat with semi-orthorhombic constraints
!
! reference: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2016
! contrib   - a.m.elena december 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,       Only : idnode,mxnode,gmax
  Use setup_module
  Use site_module,        Only : ntpatm,dens
  Use config_module,      Only : cell,volm,natms,lfrzn,weight, &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use kinetic_module,     Only : getvom
  Use core_shell_module,  Only : legshl
  Use constraints_module, Only : passcon
  Use pmf_module,         Only : passpmf

  Implicit None

  Integer,           Intent( In    ) :: isw

  Logical,           Intent( In    ) :: lvar
  Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
  Real( Kind = wp ), Intent( InOut ) :: tstep

  Real( Kind = wp ), Intent( In    ) :: sigma,taut
  Real( Kind = wp ), Intent(   Out ) :: chit

  Real( Kind = wp ), Intent( In    ) :: press,strext(1:9),taup
  Real( Kind = wp ), Intent(   Out ) :: chip,eta(1:9)

  Integer,           Intent( In    ) :: iso
  Real( Kind = wp ), Intent( In    ) :: ten,stress(1:9)

  Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Integer,           Intent( In    ) :: megcon,megpmf
  Real( Kind = wp ), Intent( InOut ) :: strcon(1:9),vircon, &
                                        strpmf(1:9),virpmf

  Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc


  Logical,           Save :: newjob = .true.
  Logical                 :: safe,lcol,lfst,lv_up,lv_dn
  Integer,           Save :: mxiter,mxkit,kit
  Integer                 :: fail(1:9),iter,i
  Real( Kind = wp ), Save :: volm0,elrc0,virlrc0,h_z
  Real( Kind = wp )       :: hstep,rstep
  Real( Kind = wp )       :: czero(1:9),celprp(1:10)
  Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                             vom(1:3)

  Real( Kind = wp ), Parameter :: beta = 7.3728e-3_wp

! uni is the diagonal unit matrix

  Real( Kind = wp ), Parameter :: &
  uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


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
     Write(nrite,'(/,1x,a,i0)') 'nst_b0 allocation failure, node: ', idnode
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

! Sort eta for iso>=1
! Initialise and get h_z for iso>1

     h_z=0
     If      (iso == 1) Then
        eta(1) = 1.0_wp ; eta(2:4) = 0.0_wp
        eta(5) = 1.0_wp ; eta(6:8) = 0.0_wp
     Else If (iso >  1) Then
        eta(2:4) = 0.0_wp
        eta(6:8) = 0.0_wp

        Call dcell(cell,celprp)
        h_z=celprp(9)
     End If

! set number of constraint+pmf shake iterations and general iteration cycles

     mxiter=1
     If (megcon > 0 .or.  megpmf > 0) Then
        mxkit=1
        mxiter=mxiter+12
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

! iterate forces, strcon and eta

     Do iter=1,mxiter

! Berendsen barostat and thermostat are not coupled
! calculate Berendsen barostat: eta, iterate strcon and strpmf

! split anisotropic from semi-isotropic barostats (iso=0,1,2,3)

        If (iso == 0) Then
           eta=uni + tstep*beta*(strcon+strpmf+stress+strkin-(press*uni+strext)*volm)/(taup*volm)
        Else
           If      (iso == 2) Then
              eta(1)=1.0_wp + tstep*beta*(strcon(1)+strpmf(1)+stress(1)+strkin(1)-(press+strext(1)-ten/h_z)*volm)/(taup*volm)
              eta(5)=1.0_wp + tstep*beta*(strcon(5)+strpmf(5)+stress(5)+strkin(5)-(press+strext(5)-ten/h_z)*volm)/(taup*volm)
           Else If (iso == 3) Then
              eta(1)=1.0_wp + tstep*beta*( 0.5_wp*                                                       &
                     (strcon(1)+strpmf(1)+stress(1)+strkin(1)+strcon(5)+strpmf(5)+stress(5)+strkin(5)) - &
                     (press+0.5_wp*(strext(1)+strext(5))-ten/h_z)*volm ) / (taup*volm)
              eta(5)=eta(1)
           End If
           eta(9)=1.0_wp + tstep*beta*(strcon(9)+strpmf(9)+stress(9)+strkin(9)-(press+strext(9))*volm)/(taup*volm)
        End If

! update velocity and position

        Do i=1,natms
           If (weight(i) > 1.0e-6_wp) Then
              tmp=hstep/weight(i)
              vxx(i)=vxt(i)+tmp*fxx(i)
              vyy(i)=vyt(i)+tmp*fyy(i)
              vzz(i)=vzt(i)+tmp*fzz(i)

              xxx(i)=tstep*vxx(i)+xxt(i)*eta(1)+yyt(i)*eta(2)+zzt(i)*eta(3)
              yyy(i)=tstep*vyy(i)+xxt(i)*eta(2)+yyt(i)*eta(5)+zzt(i)*eta(6)
              zzz(i)=tstep*vzz(i)+xxt(i)*eta(3)+yyt(i)*eta(6)+zzt(i)*eta(9)
           End If
        End Do

! SHAKE procedures

        If (megcon > 0 .or. megpmf > 0) Then
           safe=.false.
           kit =0

! update cell parameters: anisotropic

           Call mat_mul(eta,czero,cell)

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
           (mxshak,tolnce,tstep, &
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
           (mxshak,tolnce,tstep, &
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

           Do i=1,natms
              fxx(i) = fxt(i)
              fyy(i) = fyt(i)
              fzz(i) = fzt(i)
           End Do

! restart vv1

           Go To 100
        End If
     End If

! update cell parameters: anisotropic

     If (megcon == 0 .and. megpmf == 0) Call mat_mul(eta,czero,cell)

! update volume and construct a 'mock' chip

     Call dcell(cell,celprp)
     chip=celprp(10)/volm
     volm=celprp(10)

! adjust long range corrections and number density

     tmp=(volm0/volm)
     elrc=elrc0*tmp
     virlrc=virlrc0*tmp
     Do i=1,ntpatm
        dens(i)=dens0(i)*tmp
     End Do

! get h_z for iso>1

     If (iso > 1) Then
        h_z=celprp(9)
     End If

! second pass of velocity verlet algorithm

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

! integrate and apply nvt_b0_scl thermostat - full step

     Call nvt_b0_scl(1,tstep,sigma,taut,vxx,vyy,vzz,chit,strkin,engke)

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

     Call nvt_b0_scl(0,tstep,sigma,taut,vxx,vyy,vzz,chit,strkin,engke)

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
     Write(nrite,'(/,1x,a,i0)') 'nst_b0 deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine nst_b0_vv
