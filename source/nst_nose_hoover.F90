Module nst_nose_hoover

  Use kinds, Only : wp, li
  Use comms,       Only : comms_type,gmax
  Use setup
  Use site, Only : site_type
  Use configuration,      Only : imcon,cell,volm,natms,nlast,nfree, &
                                 lfrzn,lstfre,weight,               &
                                 xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use domains,     Only : map
  Use kinetics,     Only : getcom,getvom,kinstress,kinstresf,kinstrest
  Use core_shell,  Only : legshl
  Use constraints, Only : constraints_tags,apply_rattle,&
                          apply_shake, constraints_type
  Use pmf,         Only : pmf_tags,pmf_type
  Use rigid_bodies
  Use numerics,        Only : dcell, mat_mul
  Use nvt_nose_hoover, Only : nvt_h0_scl, nvt_h1_scl
  Use errors_warnings, Only : error,info
  Use thermostat, Only : thermostat_type
  Use statistics, Only : stats_type
  Use timer, Only : timer_type
  Implicit None

  Private

  Public :: nst_h0_vv, nst_h1_vv, nst_h0_scl, nst_h1_scl

Contains

  Subroutine nst_h0_vv                          &
             (isw,lvar,mndis,mxdis,mxstp,tstep, &
             degfre,stress,             &
             consv,                             &
             strkin,engke,                      &
             elrc,virlrc,cons,pmf,stat,thermo,site,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian equations of motion in
  ! molecular dynamics - velocity verlet with Nose-Hoover thermostat and
  ! barostat (anisotropic pressure control) (symplectic)
  !
  ! Parrinello-Rahman type: changing cell shape
  !
  ! thermo%iso=0 fully anisotropic barostat
  ! thermo%iso=1 semi-isotropic barostat to constant normal pressure & surface area
  ! thermo%iso=2 semi-isotropic barostat to constant normal pressure & surface tension
  !                               or with orthorhombic constraints (thermo%tension=0.0_wp)
  ! thermo%iso=3 semi-isotropic barostat with semi-orthorhombic constraints
  !
  ! Note: (1) this ensemble is modified from its original form as in
  !           reference1 to that shown in reference2, and now there is
  !           coupling between the thermostat and the barostat
  !       (2) this ensemble is not correct when there is an external
  !           field applied on the system
  !
  ! reference1: Melchionna, Ciccotti and Holian, Mol Phys 1993, 78, p533
  ! reference2: Martyna, Tuckerman, Tobias, Klein
  !             Mol. Phys., 1996, Vol. 87 (5), p. 1117
  ! reference3: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: isw

    Logical,           Intent( In    ) :: lvar
    Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ), Intent( InOut ) :: tstep

    Integer(Kind=li),  Intent( In    ) :: degfre
    Real( Kind = wp ), Intent( In    ) :: stress(1:9)

    Real( Kind = wp ), Intent(   Out ) :: consv

    Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke


    Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc
    Type( stats_type), Intent( InOut ) :: stat
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: site
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut) :: comm


    Logical,           Save :: newjob = .true.
    Logical                 :: safe,lcol,lfst,lv_up,lv_dn
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:9),iter,i
    Real( Kind = wp ), Save :: volm0,elrc0,virlrc0,h_z
    Real( Kind = wp ), Save :: qmass,ceng,pmass,chip0
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: chit0,cint0,chpzr,eta0(1:9)
    Real( Kind = wp )       :: cell0(1:9),vzero,celprp(1:10)
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                               vom(1:3),com(1:3),aaa(1:9),bbb(1:9)

  ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
    uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len =256 )  ::  message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail(1))
       Call cons%allocate_work(mxatms)
       Call pmf%allocate_work()
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
    End If
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail(7))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail(8))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nst_h0 allocation failure'
       Call error(0,message)
    End If


    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       volm0   = volm
       elrc0   = elrc
       virlrc0 = virlrc

       Allocate (dens0(1:mxatyp), Stat=fail(1))
       If (fail(1) > 0) Then
          Write(message,'(a)') 'dens0 allocation failure'
          Call error(0,message)
       End If
       Do i=1,site%ntype_atom
          dens0(i) = site%dens(i)
       End Do

  ! Sort thermo%eta for thermo%iso>=1
  ! Initialise and get h_z for thermo%iso>1

       h_z=0
       If      (thermo%iso == 1) Then
          thermo%eta(1:8) = 0.0_wp
       Else If (thermo%iso >  1) Then
          thermo%eta(2:4) = 0.0_wp
          thermo%eta(6:8) = 0.0_wp

          Call dcell(cell,celprp)
          h_z=celprp(9)
       End If

  ! inertia parameters for Nose-Hoover thermostat and barostat

       qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
       tmp   = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
       If      (thermo%iso == 0) Then
          ceng  = 2.0_wp*thermo%sigma + 3.0_wp**2*boltz*tmp
       Else If (thermo%iso == 1) Then
          ceng  = 2.0_wp*thermo%sigma + 1.0_wp*boltz*tmp
       Else If (thermo%iso == 2) Then
          ceng  = 2.0_wp*thermo%sigma + 3.0_wp*boltz*tmp
       Else If (thermo%iso == 3) Then
          ceng  = 2.0_wp*thermo%sigma + 2.0_wp*boltz*tmp
       End If
       pmass = ((2.0_wp*thermo%sigma + 3.0_wp*boltz*tmp)/3.0_wp)*thermo%tau_p**2

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

       chip0 = Sqrt( &
         thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
         thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+3
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0)Then
         Call constraints_tags(lstitr,cons,comm)
       End if

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0)Then
         Call pmf_tags(lstitr,pmf,comm)
       End If
    End If

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    qstep = 0.5_wp*hstep
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

  ! store current integration variables

       cell0=cell
       vzero=volm
       chit0=thermo%chi_t
       cint0=thermo%cint
       eta0 =thermo%eta
       chpzr=chip0

  ! calculate system centre of mass

       Call getcom(xxx,yyy,zzz,com,comm)

  100  Continue

  ! constraint virial and stress tensor

       If (cons%megcon > 0) Then
          stat%vircon=0.0_wp
          stat%strcon=0.0_wp
       End If

  ! PMF virial and stress tensor

       If (pmf%megpmf > 0) Then
          stat%virpmf=0.0_wp
          stat%strpmf=0.0_wp
       End If

       Do iter=1,mxiter

  ! integrate and apply nvt_h0_scl thermostat - 1/4 step

          Call nvt_h0_scl &
            (qstep,ceng,qmass,pmass,chip0, &
            vxx,vyy,vzz,engke,thermo,comm)

  ! constraint+pmf virial and stress

          vir=stat%vircon+stat%virpmf
          str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h0_scl barostat - 1/2 step

          Call nst_h0_scl &
             (0,hstep,degfre,pmass,thermo%chi_t,volm, &
             h_z,str,stress,         &
             vxx,vyy,vzz,strkin,engke,thermo,comm)

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

          chip0 = Sqrt( &
            thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
            thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! integrate and apply nvt_h0_scl thermostat - 1/4 step

          Call nvt_h0_scl &
            (qstep,ceng,qmass,pmass,chip0, &
            vxx,vyy,vzz,engke,thermo,comm)

  ! update velocities

          Do i=1,natms
             If (weight(i) > 1.0e-6_wp) Then
                tmp=hstep/weight(i)
                vxx(i)=vxx(i)+tmp*fxx(i)
                vyy(i)=vyy(i)+tmp*fyy(i)
                vzz(i)=vzz(i)+tmp*fzz(i)
             End If
          End Do

  ! scale cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

          aaa=tstep*thermo%eta
          Call mat_mul(aaa,aaa,bbb)
          aaa=uni+aaa+0.5_wp*bbb
          Call mat_mul(aaa,cell0,cell)

  ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

          thermo%chi_p = thermo%eta(1)+thermo%eta(5)+thermo%eta(9)
          volm = volm*Exp(tstep*thermo%chi_p)
          thermo%chi_p = thermo%chi_p / 3.0_wp

  ! update positions: second order taylor expansion of Exp(tstep*thermo%eta)

          Do i=1,natms
             If (weight(i) > 1.0e-6_wp) Then
                xt=xxt(i)-com(1)
                yt=yyt(i)-com(2)
                zt=zzt(i)-com(3)
                xxx(i) = tstep*vxx(i) + com(1) + xt*aaa(1) + yt*aaa(2) + zt*aaa(3)
                yyy(i) = tstep*vyy(i) + com(2) + xt*aaa(2) + yt*aaa(5) + zt*aaa(6)
                zzz(i) = tstep*vzz(i) + com(3) + xt*aaa(3) + yt*aaa(6) + zt*aaa(9)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
           Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
          lstitr,&
          stat,pmf,cons,tmr,comm)
          End If

  ! restore original integration parameters as well as
  ! velocities if iter < mxiter
  ! in the next iteration stat%strcon and stat%strpmf are freshly new

          If (iter < mxiter) Then
             volm=vzero
             thermo%chi_t=chit0
             thermo%cint=cint0
             thermo%eta =eta0
             chip0=chpzr

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
             If (legshl(0,i) >= 0) &
                mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
          End Do
          mxdr=Sqrt(mxdr)
          Call gmax(comm,mxdr)

          If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

  ! scale tstep and derivatives

             If (mxdr > mxdis) Then
                lv_up = .true.
                If (lv_dn) Then
                   tstep = 0.75_wp*tstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                Else
                   tstep = hstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                End If
                Write(message,'(a,1p,e12.4)') &
                  'timestep decreased, new timestep is: ', tstep
                Call info(message,.true.)
             End If
             If (mxdr < mndis) Then
                lv_dn = .true.
                If (lv_up) Then
                   tstep = 1.50_wp*tstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                Else
                   qstep = hstep
                   hstep = tstep
                   tstep = 2.00_wp*tstep
                End If
                If (tstep > mxstp) Then
                   tstep = mxstp
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                End If
                Write(message,'(a,1p,e12.4)') &
                  'timestep increased, new timestep is: ', tstep
                Call info(message,.true.)
             End If
             rstep = 1.0_wp/tstep

  ! restore initial conditions

             volm=vzero
             thermo%chi_t=chit0
             thermo%cint=cint0
             thermo%eta =eta0
             chip0=chpzr

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
       Do i=1,site%ntype_atom
          site%dens(i)=dens0(i)*tmp
       End Do

  ! get h_z for thermo%iso>1

       If (thermo%iso > 1) Then
          Call dcell(cell,celprp)
          h_z=celprp(9)
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
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_rattle(tstep,kit,&
                          pmf,cons,stat,tmr,comm)
       End If

  ! integrate and apply nvt_h0_scl thermostat - 1/4 step

       Call nvt_h0_scl &
             (qstep,ceng,qmass,pmass,chip0, &
             vxx,vyy,vzz,engke,thermo,comm)

  ! constraint+pmf virial and stress

       vir=stat%vircon+stat%virpmf
       str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h0_scl barostat - 1/2 step

       Call nst_h0_scl &
             (0,hstep,degfre,pmass,thermo%chi_t,volm, &
             h_z,str,stress,         &
             vxx,vyy,vzz,strkin,engke,thermo,comm)

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

       chip0 = Sqrt( &
         thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
         thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! integrate and apply nvt_h0_scl thermostat - 1/4 step

       Call nvt_h0_scl &
             (qstep,ceng,qmass,pmass,chip0, &
             vxx,vyy,vzz,engke,thermo,comm)

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*qmass*thermo%chi_t**2 + 0.5_wp*pmass*chip0**2 + ceng*thermo%cint + thermo%press*volm

  ! remove system centre of mass velocity

       Call getvom(vom,vxx,vyy,vzz,comm)

       Do i=1,natms
          If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
             vxx(i)=vxx(i)-vom(1)
             vyy(i)=vyy(i)-vom(2)
             vzz(i)=vzz(i)-vom(3)
          End If
       End Do

  ! update kinetic energy and stress

       Call kinstress(vxx,vyy,vzz,strkin,comm)
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,           Stat=fail(1))
       Call cons%deallocate_work()
Call pmf%deallocate_work()
Deallocate (oxt,oyt,ozt,       Stat=fail( 6))
    End If
    Deallocate (xxt,yyt,zzt,         Stat=fail(7))
    Deallocate (vxt,vyt,vzt,         Stat=fail(8))
    Deallocate (fxt,fyt,fzt,         Stat=fail(9))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nst_h0 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nst_h0_vv

  Subroutine nst_h1_vv                          &
             (isw,lvar,mndis,mxdis,mxstp,tstep, &
             degfre,degrot,stress,      &
             consv,                             &
             strkin,strknf,strknt,engke,engrot, &
             strcom,vircom,                     &
             elrc,virlrc,cons,pmf,stat,thermo,site,tmr,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
  ! RBs, equations of motion in molecular dynamics
  ! - velocity verlet with Nose-Hoover thermostat and
  ! barostat (anisotropic pressure control) (symplectic)
  !
  ! Parrinello-Rahman type: changing cell shape
  !
  ! thermo%iso=0 fully anisotropic barostat
  ! thermo%iso=1 semi-isotropic barostat to constant normal pressure & surface area
  ! thermo%iso=2 semi-isotropic barostat to constant normal pressure & surface tension
  !                               or with orthorhombic constraints (thermo%tension=0.0_wp)
  ! thermo%iso=3 semi-isotropic barostat with semi-orthorhombic constraints
  !
  ! Note: (1) this ensemble is modified from its original form as in
  !           reference1 to that shown in reference2, and now there is
  !           coupling between the thermostat and the barostat
  !       (2) this ensemble is not correct when there is an external
  !           field applied on the system
  !
  ! reference1: Melchionna, Ciccotti and Holian, Mol Phys 1993, 78, p533
  ! reference2: Martyna, Tuckerman, Tobias, Klein
  !             Mol. Phys., 1996, Vol. 87 (5), p. 1117
  ! reference3: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: isw

    Logical,           Intent( In    ) :: lvar
    Real( Kind = wp ), Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ), Intent( InOut ) :: tstep

    Integer(Kind=li),  Intent( In    ) :: degfre,degrot
    Real( Kind = wp ), Intent( In    ) :: stress(1:9)

    Real( Kind = wp ), Intent(   Out ) :: consv

    Real( Kind = wp ), Intent( InOut ) :: strkin(1:9),engke, &
                                          strknf(1:9),strknt(1:9),engrot


    Real( Kind = wp ), Intent( InOut ) :: strcom(1:9),vircom

    Real( Kind = wp ), Intent( InOut ) :: elrc,virlrc
    Type( stats_type), Intent( InOut ) :: stat
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: site
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( comms_type ), Intent( InOut) :: comm


    Logical,           Save :: newjob = .true. , &
                               unsafe = .false.
    Logical                 :: safe,lcol,lfst,lv_up,lv_dn
    Integer,           Save :: mxiter,mxkit,kit
    Integer                 :: fail(1:14),matms,iter,i,j,i1,i2, &
                               irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp ), Save :: volm0,elrc0,virlrc0,h_z
    Real( Kind = wp ), Save :: qmass,ceng,pmass,chip0
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: chit0,cint0,chpzr,eta0(1:9)
    Real( Kind = wp )       :: cell0(1:9),vzero,celprp(1:10)
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                               vom(1:3),com(1:3),aaa(1:9),bbb(1:9)
    Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                               opx,opy,opz,fmx,fmy,fmz,       &
                               tqx,tqy,tqz,trx,try,trz,       &
                               qt0,qt1,qt2,qt3,p0,p1,p2,p3,   &
                               vpx,vpy,vpz

  ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
    uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

    Real( Kind = wp ), Allocatable :: ggx(:),ggy(:),ggz(:)
    Real( Kind = wp ), Allocatable :: q0t(:),q1t(:),q2t(:),q3t(:)
    Real( Kind = wp ), Allocatable :: rgdxxt(:),rgdyyt(:),rgdzzt(:)
    Real( Kind = wp ), Allocatable :: rgdvxt(:),rgdvyt(:),rgdvzt(:)
    Real( Kind = wp ), Allocatable :: rgdoxt(:),rgdoyt(:),rgdozt(:)

    Real( Kind = wp ), Allocatable, Save :: dens0(:)
    Character ( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       Call cons%allocate_work(mxatms)
       Call pmf%allocate_work()
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
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
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nst_h1 allocation failure'
       Call error(0,message)
    End If


    If (newjob) Then
       newjob = .false.

  ! store initial values of volume, long range corrections and density

       volm0   = volm
       elrc0   = elrc
       virlrc0 = virlrc

       Allocate (dens0(1:mxatyp), Stat=fail(1))
       If (fail(1) > 0) Then
          Write(message,'(a)') 'dens0 allocation failure' 
          Call error(0,message)
       End If
       Do i=1,site%ntype_atom
          dens0(i) = site%dens(i)
       End Do

  ! Sort thermo%eta for thermo%iso>=1
  ! Initialise and get h_z for thermo%iso>1

       h_z=0
       If      (thermo%iso == 1) Then
          thermo%eta(1:8) = 0.0_wp
       Else If (thermo%iso >  1) Then
          thermo%eta(2:4) = 0.0_wp
          thermo%eta(6:8) = 0.0_wp

          Call dcell(cell,celprp)
          h_z=celprp(9)
       End If

  ! inertia parameters for Nose-Hoover thermostat and barostat

       qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
       tmp   = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
       If      (thermo%iso == 0) Then
          ceng  = 2.0_wp*thermo%sigma + 3.0_wp**2*boltz*tmp
       Else If (thermo%iso == 1) Then
          ceng  = 2.0_wp*thermo%sigma + 1.0_wp*boltz*tmp
       Else If (thermo%iso == 2) Then
          ceng  = 2.0_wp*thermo%sigma + 3.0_wp*boltz*tmp
       Else If (thermo%iso == 3) Then
          ceng  = 2.0_wp*thermo%sigma + 2.0_wp*boltz*tmp
       End If
       pmass = ((Real(degfre-degrot,wp) + 3.0_wp)/3.0_wp)*boltz*tmp*thermo%tau_p**2

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

       chip0 = Sqrt( &
         thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
         thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! set number of constraint+pmf shake iterations and general iteration cycles

       mxiter=1
       If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
          mxkit=1
          mxiter=mxiter+3
       End If
       If (cons%megcon > 0 .and. pmf%megpmf > 0) mxkit=cons%max_iter_shake

  ! unsafe positioning due to possibly locally shared RBs

       unsafe=(Any(map == comm%idnode))
    End If

  ! set matms

    matms=nlast
    If (comm%mxnode == 1) matms=natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0)Then
         Call constraints_tags(lstitr,cons,comm)
       End If 

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0)Then
         Call pmf_tags(lstitr,pmf,comm)
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
    qstep = 0.5_wp*hstep
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

  ! store current integration variables

       cell0=cell
       vzero=volm
       chit0=thermo%chi_t
       cint0=thermo%cint
       eta0 =thermo%eta
       chpzr=chip0

  ! calculate system centre of mass

       Call getcom(xxx,yyy,zzz,com,comm)

  100  Continue

  ! constraint virial and stress tensor

       If (cons%megcon > 0) Then
          stat%vircon=0.0_wp
          stat%strcon=0.0_wp
       End If

  ! PMF virial and stress tensor

       If (pmf%megpmf > 0) Then
          stat%virpmf=0.0_wp
          stat%strpmf=0.0_wp
       End If

       Do iter=1,mxiter

  ! integrate and apply nvt_h1_scl thermostat - 1/4 step

          Call nvt_h1_scl &
            (qstep,ceng,qmass,pmass,chip0, &
            vxx,vyy,vzz,                  &
            rgdvxx,rgdvyy,rgdvzz,         &
            rgdoxx,rgdoyy,rgdozz,         &
            engke,engrot,thermo,comm)

  ! constraint+pmf virial and stress

          vir=stat%vircon+stat%virpmf
          str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h1_scl barostat - 1/2 step

          Call nst_h1_scl &
             (0,hstep,degfre,degrot,pmass,thermo%chi_t,volm, &
             h_z,str,stress,strcom,         &
             vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,strkin,strknf,strknt,engke,thermo,comm)

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

          chip0 = Sqrt( &
            thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
            thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! integrate and apply nvt_h1_scl thermostat - 1/4 step

          Call nvt_h1_scl &
            (qstep,ceng,qmass,pmass,chip0, &
            vxx,vyy,vzz,                  &
            rgdvxx,rgdvyy,rgdvzz,         &
            rgdoxx,rgdoyy,rgdozz,         &
            engke,engrot,thermo,comm)

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

  ! scale cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

          aaa=tstep*thermo%eta
          Call mat_mul(aaa,aaa,bbb)
          aaa=uni+aaa+0.5_wp*bbb
          Call mat_mul(aaa,cell0,cell)

  ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

          thermo%chi_p = thermo%eta(1)+thermo%eta(5)+thermo%eta(9)
          volm = volm*Exp(tstep*thermo%chi_p)
          thermo%chi_p = thermo%chi_p / 3.0_wp

  ! update position of FPs: second order taylor expansion of Exp(tstep*thermo%eta)

          Do j=1,nfree
             i=lstfre(j)

             If (weight(i) > 1.0e-6_wp) Then
                xt=xxt(i)-com(1)
                yt=yyt(i)-com(2)
                zt=zzt(i)-com(3)
                xxx(i) = tstep*vxx(i) + com(1) + xt*aaa(1) + yt*aaa(2) + zt*aaa(3)
                yyy(i) = tstep*vyy(i) + com(2) + xt*aaa(2) + yt*aaa(5) + zt*aaa(6)
                zzz(i) = tstep*vzz(i) + com(3) + xt*aaa(3) + yt*aaa(6) + zt*aaa(9)
             End If
          End Do

  ! SHAKE procedures

          If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
           Call apply_shake(tstep,mxkit,kit,oxt,oyt,ozt,&
          lstitr,&
          stat,pmf,cons,tmr,comm)
          End If

  ! restore original integration parameters as well as
  ! velocities if iter < mxiter
  ! in the next iteration stat%strcon and stat%strpmf are freshly new

          If (iter < mxiter) Then
             volm=vzero
             thermo%chi_t=chit0
             thermo%cint=cint0
             thermo%eta =eta0
             chip0=chpzr

             Do j=1,nfree
                i=lstfre(j)

                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)
             End Do

             Do irgd=1,ntrgd
                rgdvxx(irgd) = rgdvxt(irgd)
                rgdvyy(irgd) = rgdvyt(irgd)
                rgdvzz(irgd) = rgdvzt(irgd)

                rgdoxx(irgd) = rgdoxt(irgd)
                rgdoyy(irgd) = rgdoyt(irgd)
                rgdozz(irgd) = rgdozt(irgd)
             End Do
          End If
       End Do

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

                Call images(imcon,cell0,1,x,y,z)

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

             xt=rgdxxt(irgd)-com(1)
             yt=rgdyyt(irgd)-com(2)
             zt=rgdzzt(irgd)-com(3)
             rgdxxx(irgd) = tstep*rgdvxx(irgd) + com(1) + xt*aaa(1) + yt*aaa(2) + zt*aaa(3)
             rgdyyy(irgd) = tstep*rgdvyy(irgd) + com(2) + xt*aaa(2) + yt*aaa(5) + zt*aaa(6)
             rgdzzz(irgd) = tstep*rgdvzz(irgd) + com(3) + xt*aaa(3) + yt*aaa(6) + zt*aaa(9)

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
                         xt=xxt(i)-com(1)
                         yt=yyt(i)-com(2)
                         zt=zzt(i)-com(3)
                         vxx(i)=xt*aaa(1)+yt*aaa(2)+zt*aaa(3)
                         vyy(i)=xt*aaa(2)+yt*aaa(5)+zt*aaa(6)
                         vzz(i)=xt*aaa(3)+yt*aaa(6)+zt*aaa(9)

                         x(1)=(xxx(i)-com(1))-vxx(i)
                         y(1)=(yyy(i)-com(2))-vyy(i)
                         z(1)=(zzz(i)-com(3))-vzz(i)
                         Call images(imcon,cell,1,x,y,z)
                         xxx(i)=(x(1)+com(1))+vxx(i)
                         yyy(i)=(y(1)+com(2))+vyy(i)
                         zzz(i)=(z(1)+com(3))+vzz(i)
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

             xt=rgdxxt(irgd)-com(1)
             yt=rgdyyt(irgd)-com(2)
             zt=rgdzzt(irgd)-com(3)
             rgdxxx(irgd) = com(1) + xt*aaa(1)+yt*aaa(2)+zt*aaa(3)
             rgdyyy(irgd) = com(2) + xt*aaa(2)+yt*aaa(5)+zt*aaa(6)
             rgdzzz(irgd) = com(3) + xt*aaa(3)+yt*aaa(6)+zt*aaa(9)

  ! update RB members positions

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
             If (legshl(0,i) >= 0) &
                mxdr=Max(mxdr,(xxx(i)-xxt(i))**2 + (yyy(i)-yyt(i))**2 + (zzz(i)-zzt(i))**2)
          End Do
          mxdr=Sqrt(mxdr)
          Call gmax(comm,mxdr)

          If ((mxdr < mndis .or. mxdr > mxdis) .and. tstep < mxstp) Then

  ! scale tstep and derivatives

             If (mxdr > mxdis) Then
                lv_up = .true.
                If (lv_dn) Then
                   tstep = 0.75_wp*tstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                Else
                   tstep = hstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                End If
                Write(message,'(a,1p,e12.4)') &
                  'timestep decreased, new timestep is: ', tstep
                Call info(message,.true.)
             End If
             If (mxdr < mndis) Then
                lv_dn = .true.
                If (lv_up) Then
                   tstep = 1.50_wp*tstep
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                Else
                   qstep = hstep
                   hstep = tstep
                   tstep = 2.00_wp*tstep
                End If
                If (tstep > mxstp) Then
                   tstep = mxstp
                   hstep = 0.50_wp*tstep
                   qstep = 0.50_wp*hstep
                End If
                Write(message,'(a,1p,e12.4)') &
                  'timestep increased, new timestep is: ', tstep
                Call info(message,.true.)
             End If
             rstep = 1.0_wp/tstep

  ! restore initial conditions

             volm=vzero
             thermo%chi_t=chit0
             thermo%cint=cint0
             thermo%eta =eta0
             chip0=chpzr

             Do i=1,matms
                vxx(i) = vxt(i)
                vyy(i) = vyt(i)
                vzz(i) = vzt(i)

                fxx(i) = fxt(i)
                fyy(i) = fyt(i)
                fzz(i) = fzt(i)
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

  ! adjust long range corrections and number density

       tmp=(volm0/volm)
       elrc=elrc0*tmp
       virlrc=virlrc0*tmp
       Do i=1,site%ntype_atom
          site%dens(i)=dens0(i)*tmp
       End Do

  ! get h_z for thermo%iso>1

       If (thermo%iso > 1) Then
          Call dcell(cell,celprp)
          h_z=celprp(9)
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

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_rattle(tstep,kit,&
                          pmf,cons,stat,tmr,comm)
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stress(strcom,ggx,ggy,ggz,comm)
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

  ! integrate and apply nvt_h1_scl thermostat - 1/4 step

       Call nvt_h1_scl &
             (qstep,ceng,qmass,pmass,chip0, &
             vxx,vyy,vzz,                  &
             rgdvxx,rgdvyy,rgdvzz,         &
             rgdoxx,rgdoyy,rgdozz,         &
             engke,engrot,thermo,comm)

  ! constraint+pmf virial and stress

       vir=stat%vircon+stat%virpmf
       str=stat%strcon+stat%strpmf

  ! integrate and apply nst_h1_scl barostat - 1/2 step

       Call nst_h1_scl &
             (0,hstep,degfre,degrot,pmass,thermo%chi_t,volm, &
             h_z,str,stress,strcom,         &
             vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,strkin,strknf,strknt,engke,thermo,comm)

  ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

       chip0 = Sqrt( &
         thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
         thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

  ! integrate and apply nvt_h1_scl thermostat - 1/4 step

       Call nvt_h1_scl &
             (qstep,ceng,qmass,pmass,chip0, &
             vxx,vyy,vzz,                  &
             rgdvxx,rgdvyy,rgdvzz,         &
             rgdoxx,rgdoyy,rgdozz,         &
             engke,engrot,thermo,comm)

  ! conserved quantity less kinetic and potential energy terms

       consv = 0.5_wp*qmass*thermo%chi_t**2 + 0.5_wp*pmass*chip0**2 + &
         ceng*thermo%cint + thermo%press*volm

  ! remove system centre of mass velocity

       Call getvom(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,comm)

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

       Call kinstresf(vxx,vyy,vzz,strknf,comm)
       Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

       strkin=strknf+strknt
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,            Stat=fail( 1))
       Call cons%deallocate_work( )
       Call pmf%deallocate_work()
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
       Write(message,'(a)') 'nst_h1 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nst_h1_vv

  Subroutine nst_h0_scl &
             (sw,tstep,degfre,pmass,chit,volm, &
             h_z,strcon,stress,       &
             vxx,vyy,vzz,strkin,engke,thermo,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to integrate and apply NsT barostat
  !
  ! sw=1 coupling to NVT thermostat for nst_m ensemble and
  !                                     additional scaling factor
  !
  ! sw=0 coupling to NVT thermostat for nst_h ensemble and
  !                                     no additional scaling factor
  !
  ! thermo%iso=0 fully anisotropic barostat
  ! thermo%iso=1 semi-isotropic barostat to constant normal pressure & surface area
  ! thermo%iso=2 semi-isotropic barostat to constant normal pressure & surface tension
  !                               or with orthorhombic constraints (thermo%tension=0.0_wp)
  ! thermo%iso=3 semi-isotropic barostat with semi-orthorhombic constraints
  !
  ! reference: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov december 2012
  ! contrib.  - a.m.elena december 2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: sw
    Integer(Kind=li),  Intent( In    ) :: degfre

    Real( Kind = wp ), Intent( In    ) :: tstep,pmass,chit,volm,h_z
    Real( Kind = wp ), Intent( In    ) :: strcon(1:9),stress(1:9)
    Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
    Real( Kind = wp ), Intent(   Out ) :: strkin(1:9)
    Real( Kind = wp ), Intent(   Out ) :: engke
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( comms_type ), Intent( InOut) :: comm


    Logical,     Save :: newjob = .true.

    Integer           :: i

    Real( Kind = wp ) :: a1,a2,a3,a5,a6,a9,b1,b2,b3,b5,b6,b9, vxt,vyt,vzt

    Real( Kind = wp ), Save :: rf, factor

  ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
    uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

    Real( Kind = wp ) :: hstep,qstep


  ! Initialise factor and 1/Nf for Nose-Hoover ensembles

    If (newjob) Then
       newjob = .false.

       factor = 0.0_wp
       rf = 0.0_wp
       If (sw == 1) rf=1.0_wp/Real(degfre,wp)
    End If

  ! timestep derivatives

    hstep=0.5_wp*tstep
    qstep=0.5_wp*hstep

  ! thermostat thermo%eta to 1/4*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  ! calculate kinetic contribution to stress tensor

    Call kinstress(vxx,vyy,vzz,strkin,comm)

  ! kinetic energy

    engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  ! barostat thermo%eta to 1/2*tstep

    If (sw == 1) factor = 2.0_wp*engke*rf

  ! split anisotropic from semi-isotropic barostats (thermo%iso=0,1,2,3)

    If (thermo%iso == 0) Then
       thermo%eta=thermo%eta + hstep*(strcon+stress+strkin + factor*uni - &
         (thermo%press*uni+thermo%stress)*volm)/pmass
    Else
       If      (thermo%iso == 2) Then
          thermo%eta(1)=thermo%eta(1) + hstep*(strcon(1)+stress(1)+strkin(1) + &
            factor - (thermo%press+thermo%stress(1)-thermo%tension/h_z)*volm)/pmass
          thermo%eta(5)=thermo%eta(5) + hstep*(strcon(5)+stress(5)+strkin(5) + &
            factor - (thermo%press+thermo%stress(5)-thermo%tension/h_z)*volm)/pmass
       Else If (thermo%iso == 3) Then
          thermo%eta(1)=0.5_wp*(thermo%eta(1)+thermo%eta(5)) + hstep*( 0.5_wp*                        &
                 (strcon(1)+stress(1)+strkin(1)+strcon(5)+stress(5)+strkin(5)) + &
                 factor - (thermo%press+0.5_wp*(thermo%stress(1)+thermo%stress(5))-thermo%tension/h_z)*volm ) / pmass
          thermo%eta(5)=thermo%eta(1)
       End If
       thermo%eta(9)=thermo%eta(9) + hstep*(strcon(9)+stress(9)+strkin(9) + &
         factor - (thermo%press+thermo%stress(9))*volm)/pmass
    End If

  ! thermostat thermo%eta to 2/4*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  ! barostat the velocities to full 1*tstep
  ! second order taylor expansion of Exp[-tstep*(thermo%eta+factor*I)],
  ! where I is the unit tensor
  ! factor = Tr(thermo%eta)/Nf if sw=1, where Nf is degfre,
  ! else if sw=0 then factor=0, by default

    If (sw == 1) factor = (thermo%eta(1)+thermo%eta(5)+thermo%eta(9))*rf

    a1 = -tstep*(thermo%eta(1)+factor)
    a2 = -tstep*thermo%eta(2)
    a3 = -tstep*thermo%eta(3)
    a5 = -tstep*(thermo%eta(5)+factor)
    a6 = -tstep*thermo%eta(6)
    a9 = -tstep*(thermo%eta(9)+factor)

    b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
    b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
    b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
    b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
    b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
    b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

    Do i=1,natms
       vxt=vxx(i)
       vyt=vyy(i)
       vzt=vzz(i)

       vxx(i) = b1*vxt + b2*vyt + b3*vzt
       vyy(i) = b2*vxt + b5*vyt + b6*vzt
       vzz(i) = b3*vxt + b6*vyt + b9*vzt
    End Do

  ! thermostat thermo%eta to 2/4*tstep

  ! calculate kinetic contribution to stress tensor

    Call kinstress(vxx,vyy,vzz,strkin,comm)

  ! kinetic energy

    engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  ! thermostat thermo%eta to 3/4*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  ! barostat thermo%eta to full (2/2)*tstep

    If (sw == 1) factor = 2.0_wp*engke*rf

  ! split anisotropic from semi-isotropic barostats (thermo%iso=0,1,2,3)

    If (thermo%iso == 0) Then
       thermo%eta=thermo%eta + hstep*(strcon+stress+strkin + factor*uni - &
         (thermo%press*uni+thermo%stress)*volm)/pmass
    Else
       If      (thermo%iso == 2) Then
          thermo%eta(1)=thermo%eta(1) + hstep*(strcon(1)+stress(1)+strkin(1) + &
            factor - (thermo%press+thermo%stress(1)-thermo%tension/h_z)*volm)/pmass
          thermo%eta(5)=thermo%eta(5) + hstep*(strcon(5)+stress(5)+strkin(5) + &
            factor - (thermo%press+thermo%stress(5)-thermo%tension/h_z)*volm)/pmass
       Else If (thermo%iso == 3) Then
          thermo%eta(1)=0.5_wp*(thermo%eta(1)+thermo%eta(5)) + hstep*( 0.5_wp*                        &
                 (strcon(1)+stress(1)+strkin(1)+strcon(5)+stress(5)+strkin(5)) + &
                 factor - (thermo%press+0.5_wp*(thermo%stress(1)+thermo%stress(5))-thermo%tension/h_z)*volm ) / pmass
          thermo%eta(5)=thermo%eta(1)
       End If
       thermo%eta(9)=thermo%eta(9) + hstep*(strcon(9)+stress(9)+strkin(9) + &
         factor - (thermo%press+thermo%stress(9))*volm)/pmass
    End If

  ! thermostat thermo%eta to full (4/4)*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  End Subroutine nst_h0_scl

  Subroutine nst_h1_scl &
             (sw,tstep,degfre,degrot,pmass,chit,volm,  &
             h_z,strcon,stress,strcom,        &
             vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,strkin,strknf,strknt,engke,thermo,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 routine to integrate and apply NsT barostat
  ! when singled RBs are present
  !
  ! sw=1 coupling to NVT thermostat for nst_m ensemble and
  !                                     additional scaling factor
  !
  ! sw=0 coupling to NVT thermostat for nst_h ensemble and
  !                                     no additional scaling factor
  !
  ! thermo%iso=0 fully anisotropic barostat
  ! thermo%iso=1 semi-isotropic barostat to constant normal pressure & surface area
  ! thermo%iso=2 semi-isotropic barostat to constant normal pressure & surface tension
  !                               or with orthorhombic constraints (thermo%tension=0.0_wp)
  ! thermo%iso=3 semi-isotropic barostat with semi-orthorhombic constraints
  !
  ! reference: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov december 2012
  ! contrib   - a.m.elena december 2017
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,           Intent( In    ) :: sw
    Integer(Kind=li),  Intent( In    ) :: degfre,degrot

    Real( Kind = wp ), Intent( In    ) :: tstep,pmass,chit,volm,h_z
    Real( Kind = wp ), Intent( In    ) :: strcon(1:9),stress(1:9),strcom(1:9)
    Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
    Real( Kind = wp ), Intent( InOut ) :: rgdvxx(1:mxrgd),rgdvyy(1:mxrgd),rgdvzz(1:mxrgd)
    Real( Kind = wp ), Intent(   Out ) :: strkin(1:9),strknf(1:9),strknt(1:9)
    Real( Kind = wp ), Intent(   Out ) :: engke
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( comms_type ), Intent( InOut) :: comm


    Logical,     Save :: newjob = .true.

    Integer           :: i,j,irgd

    Real( Kind = wp ) :: a1,a2,a3,a5,a6,a9,b1,b2,b3,b5,b6,b9, vxt,vyt,vzt

  ! initialise factor for Nose-Hoover ensembles

    Real( Kind = wp ), Save :: rf, factor

  ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
    uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

    Real( Kind = wp ) :: hstep,qstep


  ! Initialise factor and 1/Nf for Nose-Hoover ensembles

    If (newjob) Then
       newjob = .false.

       factor = 0.0_wp
       rf = 0.0_wp
       If (sw == 1) rf=1.0_wp/Real(degfre-degrot,wp)
    End If

  ! timestep derivatives

    hstep=0.5_wp*tstep
    qstep=0.5_wp*hstep

  ! thermostat thermo%eta to 1/4*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  ! calculate kinetic contributions to stress tensor

    Call kinstresf(vxx,vyy,vzz,strknf,comm)
    Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

    strkin=strknf+strknt

  ! kinetic energy

    engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  ! barostat thermo%eta to 1/2*tstep

    If (sw == 1) factor = 2.0_wp*engke*rf

  ! split anisotropic from semi-isotropic barostats (thermo%iso=0,1,2,3)

    If (thermo%iso == 0) Then
       thermo%eta=thermo%eta + hstep*(strcom+strcon+stress+strkin + factor*uni - &
         (thermo%press*uni+thermo%stress)*volm)/pmass
    Else
       If      (thermo%iso == 2) Then
          thermo%eta(1)=thermo%eta(1) + hstep*(strcom(1)+strcon(1)+stress(1)+strkin(1) + &
            factor - (thermo%press+thermo%stress(1)-thermo%tension/h_z)*volm)/pmass
          thermo%eta(5)=thermo%eta(5) + hstep*(strcom(5)+strcon(5)+stress(5)+strkin(5) + &
            factor - (thermo%press+thermo%stress(5)-thermo%tension/h_z)*volm)/pmass
       Else If (thermo%iso == 3) Then
          thermo%eta(1)=0.5_wp*(thermo%eta(1)+thermo%eta(5)) + hstep*( 0.5_wp* &
                 (strcom(1)+strcon(1)+stress(1)+strkin(1)+strcom(5)+strcon(5)+stress(5)+strkin(5)) + &
                 factor - (thermo%press+0.5_wp*(thermo%stress(1)+thermo%stress(5))-thermo%tension/h_z)*volm ) / pmass
          thermo%eta(5)=thermo%eta(1)
       End If
       thermo%eta(9)=thermo%eta(9) + hstep*(strcom(9)+strcon(9)+stress(9)+strkin(9) + &
         factor - (thermo%press+thermo%stress(9))*volm)/pmass
    End If

  ! thermostat thermo%eta to 2/4*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  ! barostat the velocities to full 1*tstep
  ! second order taylor expansion of Exp[-tstep*(thermo%eta+factor*I)],
  ! where I is the unit tensor
  ! factor = Tr(thermo%eta)/Nf if sw=1, where Nf is degfre-degrot,
  ! else if sw=0 then factor=0, by default

    If (sw == 1) factor = (thermo%eta(1)+thermo%eta(5)+thermo%eta(9))*rf

    a1 = -tstep*(thermo%eta(1)+factor)
    a2 = -tstep*thermo%eta(2)
    a3 = -tstep*thermo%eta(3)
    a5 = -tstep*(thermo%eta(5)+factor)
    a6 = -tstep*thermo%eta(6)
    a9 = -tstep*(thermo%eta(9)+factor)

    b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
    b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
    b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
    b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
    b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
    b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

    Do j=1,nfree
       i=lstfre(j)

       vxt=vxx(i)
       vyt=vyy(i)
       vzt=vzz(i)

       vxx(i) = b1*vxt + b2*vyt + b3*vzt
       vyy(i) = b2*vxt + b5*vyt + b6*vzt
       vzz(i) = b3*vxt + b6*vyt + b9*vzt
    End Do

    Do irgd=1,ntrgd
       vxt=rgdvxx(irgd)
       vyt=rgdvyy(irgd)
       vzt=rgdvzz(irgd)

       rgdvxx(irgd) = b1*vxt + b2*vyt + b3*vzt
       rgdvyy(irgd) = b2*vxt + b5*vyt + b6*vzt
       rgdvzz(irgd) = b3*vxt + b6*vyt + b9*vzt
    End Do

  ! thermostat thermo%eta to 2/4*tstep

  ! calculate kinetic contributions to stress tensor

    Call kinstresf(vxx,vyy,vzz,strknf,comm)
    Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt,comm)

    strkin=strknf+strknt

  ! kinetic energy

    engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  ! thermostat thermo%eta to 3/4*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  ! barostat thermo%eta to full (2/2)*tstep

    If (sw == 1) factor = 2.0_wp*engke*rf

  ! split anisotropic from semi-isotropic barostats (thermo%iso=0,1,2,3)

    If (thermo%iso == 0) Then
       thermo%eta=thermo%eta + hstep*(strcom+strcon+stress+strkin + factor*uni - (thermo%press*uni+thermo%stress)*volm)/pmass
    Else
       If      (thermo%iso == 2) Then
          thermo%eta(1)=thermo%eta(1) + hstep*(strcom(1)+strcon(1)+stress(1)+strkin(1) + &
            factor - (thermo%press+thermo%stress(1)-thermo%tension/h_z)*volm)/pmass
          thermo%eta(5)=thermo%eta(5) + hstep*(strcom(5)+strcon(5)+stress(5)+strkin(5) + &
            factor - (thermo%press+thermo%stress(5)-thermo%tension/h_z)*volm)/pmass
       Else If (thermo%iso == 3) Then
          thermo%eta(1)=0.5_wp*(thermo%eta(1)+thermo%eta(5)) + hstep*( 0.5_wp* &
                 (strcom(1)+strcon(1)+stress(1)+strkin(1)+strcom(5)+strcon(5)+stress(5)+strkin(5)) + &
                 factor - (thermo%press+0.5_wp*(thermo%stress(1)+thermo%stress(5))-thermo%tension/h_z)*volm ) / pmass
          thermo%eta(5)=thermo%eta(1)
       End If
       thermo%eta(9)=thermo%eta(9) + hstep*(strcom(9)+strcon(9)+stress(9)+strkin(9) + &
         factor - (thermo%press+thermo%stress(9))*volm)/pmass
    End If

  ! thermostat thermo%eta to full (4/4)*tstep

    thermo%eta=thermo%eta*Exp(-qstep*chit)

  End Subroutine nst_h1_scl
End Module nst_nose_hoover
