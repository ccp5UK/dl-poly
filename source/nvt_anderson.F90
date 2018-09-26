Module nvt_anderson
  Use kinds,         Only : wp
  Use comms,         Only : comms_type,gsum,gmax
  Use domains,       Only : domains_type
  Use setup,         Only : boltz,zero_plus,mxatms
  Use site, Only : site_type
  Use configuration, Only : configuration_type
  Use particle,      Only : corePart
  Use kinetics,      Only : getvom,getknr,kinstress,kinstresf,kinstrest
  Use constraints,   Only : constraints_tags,apply_shake, &
                            apply_rattle,constraints_type
  Use pmf,           Only : pmf_tags,pmf_type
  Use rigid_bodies,  Only : rigid_bodies_type,getrotmat,no_squish,rigid_bodies_stress
  Use numerics, Only : seed_type,images,local_index,box_mueller_saru3,sarurnd
  Use shared_units, Only : update_shared_units,update_shared_units_int
  Use errors_warnings, Only : error,info
  Use core_shell, Only : core_shell_type,SHELL_ADIABATIC
Use statistics, Only : stats_type
  Use timer, Only : timer_type
  Use thermostat, Only : adjust_timestep,thermostat_type
Use core_shell, Only : core_shell_type
  Implicit None

  Private

  Public :: nvt_a0_vv, nvt_a1_vv

Contains

  Subroutine nvt_a0_vv(isw,lvar,mndis,mxdis,mxstp,tstep,nstep, &
      strkin,engke,cshell,cons,pmf,stat,thermo,sites,domain,   &
      tmr,config,seed,comm)

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
  ! author    - i.t.todorov august 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: isw
    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep
    Integer,            Intent( In    ) :: nstep
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( configuration_type ),   Intent( InOut ) :: config
    Type(seed_type), Intent(InOut) :: seed
    Type( comms_type ), Intent( InOut) :: comm


    Logical                 :: safe,lcol,lfst
    Integer                 :: fail(1:11),i,j,k,ntp,  &
                               stp,i1,i2, &
                               matms
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                               scale,tkin,vom(1:3)


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)

  ! q. index arrays and tp. sum arrays

    Integer,           Allocatable :: qn(:),tpn(:)
    Integer,           Allocatable :: qs(:,:),tps(:)
    Character( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       Call cons%allocate_work(mxatms)
       Call pmf%allocate_work()
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
    End If
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 7))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 8))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail( 9))
    Allocate (qn(1:mxatms),tpn(0:comm%mxnode-1),                    Stat=fail(10))
    Allocate (qs(0:2,1:cshell%mxshl),tps(0:comm%mxnode-1),                 Stat=fail(11))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_a0 allocation failure'
       Call error(0,message)
    End If


    If (thermo%newjob) Then
       thermo%newjob = .false.

  ! set number of constraint+pmf shake iterations

       If (cons%megcon > 0 .or.  pmf%megpmf > 0) thermo%mxkit=1
       If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit=cons%max_iter_shake
    End If

  ! set matms

    matms=config%nlast
    If (comm%mxnode == 1) matms=config%natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:config%natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms

       If (cons%megcon > 0) Then
         Call constraints_tags(lstitr,cons,config,comm)
       End If

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms

       If (pmf%megpmf > 0) Then
         Call pmf_tags(lstitr,pmf,config,comm)
       End if
    End If

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! store initial values

       Do i=1,config%natms
          xxt(i) = config%parts(i)%xxx
          yyt(i) = config%parts(i)%yyy
          zzt(i) = config%parts(i)%zzz

          vxt(i) = config%vxx(i)
          vyt(i) = config%vyy(i)
          vzt(i) = config%vzz(i)

          fxt(i) = config%parts(i)%fxx
          fyt(i) = config%parts(i)%fyy
          fzt(i) = config%parts(i)%fzz
       End Do

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

  ! update velocity and position

       Do i=1,config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
             tmp=hstep/config%weight(i)
             config%vxx(i)=vxt(i)+tmp*fxt(i)
             config%vyy(i)=vyt(i)+tmp*fyt(i)
             config%vzz(i)=vzt(i)+tmp*fzt(i)

             config%parts(i)%xxx=xxt(i)+tstep*config%vxx(i)
             config%parts(i)%yyy=yyt(i)+tstep*config%vyy(i)
             config%parts(i)%zzz=zzt(i)+tstep*config%vzz(i)
          End If
       End Do

  ! SHAKE procedures
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep,thermo%mxkit,thermo%kit,oxt,oyt,ozt,&
          lstitr,stat,pmf,cons,domain,tmr,config,comm)
      End If

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,config%natms,config%parts,&
 xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
            Call info(message,.true.)

  ! restart vv1

             Go To 100
          End If
       End If

  ! second stage of velocity verlet algorithm

    Else

  ! update velocity

       Do i=1,config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
             tmp=hstep/config%weight(i)
             config%vxx(i)=config%vxx(i)+tmp*config%parts(i)%fxx
             config%vyy(i)=config%vyy(i)+tmp*config%parts(i)%fyy
             config%vzz(i)=config%vzz(i)+tmp*config%parts(i)%fzz
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,thermo%kit,pmf,cons,stat,domain,tmr,config,comm)
       End If

  ! Andersen Thermostat
  !
  ! qualify non-shell, non-frozen particles (n) for a random kick
  ! derive related shells (s)

  ! tpn(comm%idnode) number of thermostatted particles on this node (comm%idnode)
  ! ntp - grand total of non-shell, non-frozen particles to thermostat

       qn(1:config%natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
       qs(0:2,1:cshell%ntshl) = 0 ! unqualified core-shell unit with a local shell

       j = 0
       tkin = 0.0_wp
       mxdr = 0.0_wp
       scale = tstep/thermo%tau_t
       Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0) Then
             If (sarurnd(seed,config%ltg(i),0,nstep) <= scale) Then
                j = j + 1
                qn(i) = 1

                If (sites%dof_site(config%lsite(i)) > zero_plus) mxdr = mxdr + sites%dof_site(config%lsite(i))

  ! Get gaussian distribution (unit variance)

                Call box_mueller_saru3(seed,config%ltg(i),nstep,xxt(i),yyt(i),zzt(i))

  ! Get scaler to target variance/Sqrt(config%weight)

                tmp = 1.0_wp/Sqrt(config%weight(i))
                xxt(i) = xxt(i)*tmp
                yyt(i) = yyt(i)*tmp
                zzt(i) = zzt(i)*tmp

                tkin = tkin + config%weight(i)*(xxt(i)**2+yyt(i)**2+zzt(i)**2)
             End If
          End If
       End Do
       tpn(comm%idnode) = j
       Do i=0,comm%mxnode-1
         If (i /= comm%idnode) tpn(i) = 0
       End Do
       Call gsum(comm,tpn)
       ntp = Sum(tpn)

       If (ntp == 0) Go To 200

       Call gsum(comm,tkin)
       If (tkin <= zero_plus) tkin = 1.0_wp
       Call gsum(comm,mxdr)

  ! Scale to target temperature and apply thermostat

       scale = Sqrt(mxdr * boltz * thermo%temp / tkin)
       tmp = Sqrt(1.0_wp-thermo%soft**2)*scale

       Do i=1,config%natms
          If (qn(i) == 1) Then
             If (thermo%soft <= zero_plus) Then ! New target velocity
                config%vxx(i) = xxt(i)*scale
                config%vyy(i) = yyt(i)*scale
                config%vzz(i) = zzt(i)*scale
             Else ! Softened velocity (mixture between old & new)
                config%vxx(i) = thermo%soft*config%vxx(i) + tmp*xxt(i)
                config%vyy(i) = thermo%soft*config%vyy(i) + tmp*yyt(i)
                config%vzz(i) = thermo%soft*config%vzz(i) + tmp*zzt(i)
             End If
          End If
       End Do

  ! tps(comm%idnode) number of thermostatted core-shell units on this node (comm%idnode)
  ! stp - grand total of core-shell units to thermostat

       j = 0
       If (cshell%keyshl == SHELL_ADIABATIC) Then
          If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
             qn(config%natms+1:config%nlast) = 0
             Call update_shared_units_int(config,cshell%lishp_shl, &
               cshell%lashp_shl,qn,domain,comm)
          End If

          If (cshell%ntshl > 0) Then
             Do k=1,cshell%ntshl
                i1=local_index(cshell%listshl(1,k),matms,config%lsi,config%lsa)
                i2=local_index(cshell%listshl(2,k),matms,config%lsi,config%lsa)

                If (qn(i1) == 1 .and. i2 > 0 .and. i2 <= config%natms) Then
                   j = j + 1

                   qs(0,k)=1
                   qs(1,k)=i1
                   qs(2,k)=i2
                End If
             End Do
          End If
       End If
       tps(comm%idnode) = j
       Do i=0,comm%mxnode-1
         If (i /= comm%idnode) tps(i) = 0
       End Do
       Call gsum(comm,tps)
       stp = Sum(tps)

  ! Thermalise the shells on hit cores

       If (stp > 0) Then
          If (cshell%lshmv_shl) Then
            Call update_shared_units(config,cshell%lishp_shl, &
              cshell%lashp_shl,config%vxx,config%vyy,config%vzz,domain,comm)
          End If

          If (tps(comm%idnode) > 0) Then
             Do k=1,cshell%ntshl
                If (qs(0,k) == 1) Then
                   i1=qs(1,k)
                   i2=qs(2,k)

                   config%vxx(i2)=config%vxx(i1)
                   config%vyy(i2)=config%vyy(i1)
                   config%vzz(i2)=config%vzz(i1)
                End If
             End Do
          End If
       End If

  ! remove system centre of mass velocity (random momentum walk)

       Call getvom(vom,config%vxx,config%vyy,config%vzz,config,comm)

       Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
             config%vxx(i)=config%vxx(i)-vom(1)
             config%vyy(i)=config%vyy(i)-vom(2)
             config%vzz(i)=config%vzz(i)-vom(3)
          End If
       End Do

  200  Continue

  ! update kinetic energy and stress

       Call kinstress(config%vxx,config%vyy,config%vzz,strkin,config,comm)
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,           Stat=fail( 1))
       Call cons%deallocate_work()
       Call pmf%deallocate_work()
       Deallocate (oxt,oyt,ozt,       Stat=fail( 6))
    End If
    Deallocate (xxt,yyt,zzt,         Stat=fail( 7))
    Deallocate (vxt,vyt,vzt,         Stat=fail( 8))
    Deallocate (fxt,fyt,fzt,         Stat=fail( 9))
    Deallocate (qn,tpn,              Stat=fail(10))
    Deallocate (qs,tps,              Stat=fail(11))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_a0 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nvt_a0_vv

  Subroutine nvt_a1_vv(isw,lvar,mndis,mxdis,mxstp,tstep,nstep, &
      strkin,strknf,strknt,engke,engrot,strcom,vircom, &
      cshell,cons,pmf,stat,thermo,sites,rigid,domain,  &
      tmr,config,seed,comm)

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

    Integer,            Intent( In    ) :: isw
    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep
    Integer,            Intent( In    ) :: nstep
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke, &
                                           strknf(1:9),strknt(1:9),engrot
    Real( Kind = wp ),  Intent( InOut ) :: strcom(1:9),vircom
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( configuration_type ),   Intent( InOut ) :: config
    Type(seed_type), Intent(InOut) :: seed
    Type( comms_type ), Intent( InOut ) :: comm


    Logical                 :: safe,lcol,lfst
    Integer                 :: fail(1:17),i,j,k,ntp,  &
                               stp,i1,i2, &
                               matms,rtp,irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp )       :: hstep,rstep
    Real( Kind = wp )       :: xt,yt,zt,vir,str(1:9),mxdr,tmp, &
                               scale,tkin,vom(1:3)
    Real( Kind = wp )       :: x(1:1),y(1:1),z(1:1),rot(1:9), &
                               opx,opy,opz,fmx,fmy,fmz,       &
                               tqx,tqy,tqz,trx,try,trz,       &
                               qt0,qt1,qt2,qt3,p0,p1,p2,p3,   &
                               vpx,vpy,vpz


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

  ! q. index arrays and tp. sum arrays

    Integer,           Allocatable :: qn(:),tpn(:)
    Integer,           Allocatable :: qs(:,:),tps(:)
    Integer,           Allocatable :: qr(:),tpr(:)
    Character( Len = 256 ) :: message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Allocate (lstitr(1:mxatms),                                  Stat=fail( 1))
       Call cons%allocate_work(mxatms)
       Call pmf%allocate_work()
       Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms),         Stat=fail(6))
    End If
    Allocate (ggx(1:rigid%max_list*rigid%max_rigid), &
      ggy(1:rigid%max_list*rigid%max_rigid), &
      ggz(1:rigid%max_list*rigid%max_rigid), Stat=fail( 7))
    Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms),            Stat=fail( 8))
    Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),            Stat=fail( 9))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms),            Stat=fail(10))
    Allocate (q0t(1:rigid%max_rigid), &
      q1t(1:rigid%max_rigid), &
      q2t(1:rigid%max_rigid), &
      q3t(1:rigid%max_rigid),  Stat=fail(11))
    Allocate (rgdxxt(1:rigid%max_rigid), &
      rgdyyt(1:rigid%max_rigid), &
      rgdzzt(1:rigid%max_rigid),      Stat=fail(12))
    Allocate (rgdvxt(1:rigid%max_rigid), &
      rgdvyt(1:rigid%max_rigid), &
      rgdvzt(1:rigid%max_rigid),      Stat=fail(13))
    Allocate (rgdoxt(1:rigid%max_rigid), &
      rgdoyt(1:rigid%max_rigid), &
      rgdozt(1:rigid%max_rigid),      Stat=fail(14))
    Allocate (qn(1:mxatms),tpn(0:comm%mxnode-1),                    Stat=fail(15))
    Allocate (qs(0:2,1:cshell%mxshl),tps(0:comm%mxnode-1),                 Stat=fail(16))
    Allocate (qr(1:rigid%max_rigid),tpr(0:comm%mxnode-1),                     Stat=fail(17))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_a1 allocation failure'
       Call error(0,message)
    End If


    If (thermo%newjob) Then
       thermo%newjob = .false.

  ! set number of constraint+pmf shake iterations

       If (cons%megcon > 0 .or.  pmf%megpmf > 0) thermo%mxkit=1
       If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit=cons%max_iter_shake

  ! thermo%unsafe positioning due to possibly locally shared RBs

       thermo%unsafe=(Any(domain%map == comm%idnode))
    End If

  ! set matms

    matms=config%nlast
    If (comm%mxnode == 1) matms=config%natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       lstitr(1:config%natms)=.false. ! initialise lstitr

  ! construct current bond vectors and listot array (shared
  ! constraint atoms) for iterative bond algorithms
       If (cons%megcon > 0) Then
         Call constraints_tags(lstitr,cons,config,comm)
       End If

  ! construct current PMF constraint vectors and shared description
  ! for iterative PMF constraint algorithms
       If (pmf%megpmf > 0) Then
         Call pmf_tags(lstitr,pmf,config,comm)
       End If
    End If

  ! Get the RB particles vectors wrt the RB's COM

    krgd=0
    Do irgd=1,rigid%n_types
       rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

       lrgd=rigid%list(-1,irgd)
       If (rigid%frozen(0,rgdtyp) < lrgd) Then
          Do jrgd=1,lrgd
             krgd=krgd+1

             i=rigid%index_local(jrgd,irgd) ! local index of particle/site

  ! COM distances

             ggx(krgd)=config%parts(i)%xxx-rigid%xxx(irgd)
             ggy(krgd)=config%parts(i)%yyy-rigid%yyy(irgd)
             ggz(krgd)=config%parts(i)%zzz-rigid%zzz(irgd)
          End Do
       End If
    End Do

  ! minimum image convention for bond vectors

    Call images(config%imcon,config%cell,krgd,ggx,ggy,ggz)

  ! timestep derivatives

    hstep = 0.5_wp*tstep
    rstep = 1.0_wp/tstep

  ! first pass of velocity verlet algorithm

    If (isw == 0) Then

  ! store initial values

       Do i=1,matms
          xxt(i) = config%parts(i)%xxx
          yyt(i) = config%parts(i)%yyy
          zzt(i) = config%parts(i)%zzz

          vxt(i) = config%vxx(i)
          vyt(i) = config%vyy(i)
          vzt(i) = config%vzz(i)

          fxt(i) = config%parts(i)%fxx
          fyt(i) = config%parts(i)%fyy
          fzt(i) = config%parts(i)%fzz
       End Do

       Do irgd=1,rigid%n_types
          q0t(irgd)=rigid%q0(irgd)
          q1t(irgd)=rigid%q1(irgd)
          q2t(irgd)=rigid%q2(irgd)
          q3t(irgd)=rigid%q3(irgd)

          rgdxxt(irgd) = rigid%xxx(irgd)
          rgdyyt(irgd) = rigid%yyy(irgd)
          rgdzzt(irgd) = rigid%zzz(irgd)

          rgdvxt(irgd) = rigid%vxx(irgd)
          rgdvyt(irgd) = rigid%vyy(irgd)
          rgdvzt(irgd) = rigid%vzz(irgd)

          rgdoxt(irgd) = rigid%oxx(irgd)
          rgdoyt(irgd) = rigid%oyy(irgd)
          rgdozt(irgd) = rigid%ozz(irgd)
       End Do

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

  ! update velocity and position of FPs

       Do j=1,config%nfree
          i=config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
             tmp=hstep/config%weight(i)
             config%vxx(i)=vxt(i)+tmp*fxt(i)
             config%vyy(i)=vyt(i)+tmp*fyt(i)
             config%vzz(i)=vzt(i)+tmp*fzt(i)

             config%parts(i)%xxx=xxt(i)+tstep*config%vxx(i)
             config%parts(i)%yyy=yyt(i)+tstep*config%vyy(i)
             config%parts(i)%zzz=zzt(i)+tstep*config%vzz(i)
          End If
       End Do

  ! SHAKE procedures
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep,thermo%mxkit,thermo%kit,oxt,oyt,ozt,&
          lstitr,stat,pmf,cons,domain,tmr,config,comm)
      End If
  ! update velocity and position of RBs

       krgd=0
       Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) < lrgd) Then

  ! calculate COM force and torque

             fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
             tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rigid%frozen(0,rgdtyp) == 0) Then
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

             If (rigid%frozen(0,rgdtyp) > 1) Then
                i1=rigid%index_local(rigid%index_global(1,rgdtyp),irgd)
                i2=rigid%index_local(rigid%index_global(2,rgdtyp),irgd)

                x(1)=xxt(i1)-xxt(i2)
                y(1)=yyt(i1)-yyt(i2)
                z(1)=zzt(i1)-zzt(i2)

                Call images(config%imcon,config%cell,1,x,y,z)

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

             opx=rgdoxt(irgd)*rigid%rix(1,rgdtyp)
             opy=rgdoyt(irgd)*rigid%riy(1,rgdtyp)
             opz=rgdozt(irgd)*rigid%riz(1,rgdtyp)

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
             (tstep,rigid%rix(2,rgdtyp),rigid%riy(2,rgdtyp),rigid%riz(2,rgdtyp), &
             rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),p0,p1,p2,p3)
             Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

  ! update RB angular & COM velocities to half step

             opx=0.5_wp*(-rigid%q1(irgd)*p0+rigid%q0(irgd)*p1+rigid%q3(irgd)*p2-rigid%q2(irgd)*p3)
             opy=0.5_wp*(-rigid%q2(irgd)*p0-rigid%q3(irgd)*p1+rigid%q0(irgd)*p2+rigid%q1(irgd)*p3)
             opz=0.5_wp*(-rigid%q3(irgd)*p0+rigid%q2(irgd)*p1-rigid%q1(irgd)*p2+rigid%q0(irgd)*p3)

             rigid%oxx(irgd)=opx*rigid%rix(2,rgdtyp)
             rigid%oyy(irgd)=opy*rigid%riy(2,rgdtyp)
             rigid%ozz(irgd)=opz*rigid%riz(2,rgdtyp)

             tmp=hstep/rigid%weight(0,rgdtyp)
             rigid%vxx(irgd)=rgdvxt(irgd)+tmp*fmx
             rigid%vyy(irgd)=rgdvyt(irgd)+tmp*fmy
             rigid%vzz(irgd)=rgdvzt(irgd)+tmp*fmz

  ! update RB COM to full step

             rigid%xxx(irgd)=rgdxxt(irgd)+tstep*rigid%vxx(irgd)
             rigid%yyy(irgd)=rgdyyt(irgd)+tstep*rigid%vyy(irgd)
             rigid%zzz(irgd)=rgdzzt(irgd)+tstep*rigid%vzz(irgd)

  ! update RB members positions and halfstep velocities

             Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= config%natms) Then
                   If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                      x(1)=rigid%x(jrgd,rgdtyp)
                      y(1)=rigid%y(jrgd,rgdtyp)
                      z(1)=rigid%z(jrgd,rgdtyp)

  ! new atomic positions

                      config%parts(i)%xxx=rot(1)*x(1)+rot(2)*y(1)+rot(3)*z(1) + rigid%xxx(irgd)
                      config%parts(i)%yyy=rot(4)*x(1)+rot(5)*y(1)+rot(6)*z(1) + rigid%yyy(irgd)
                      config%parts(i)%zzz=rot(7)*x(1)+rot(8)*y(1)+rot(9)*z(1) + rigid%zzz(irgd)

  ! new atomic velocities in body frame

                      vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                      vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                      vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

  ! DD bound positions

                      If (thermo%unsafe) Then
                         x(1)=config%parts(i)%xxx-xxt(i)
                         y(1)=config%parts(i)%yyy-yyt(i)
                         z(1)=config%parts(i)%zzz-zzt(i)
                         Call images(config%imcon,config%cell,1,x,y,z)
                         config%parts(i)%xxx=x(1)+xxt(i)
                         config%parts(i)%yyy=y(1)+yyt(i)
                         config%parts(i)%zzz=z(1)+zzt(i)
                      End If

  ! new atomic velocities in lab frame

                      config%vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                      config%vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                      config%vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                   End If
                End If
             End Do

          End If
       End Do

  ! check timestep for variable timestep

       If (lvar) Then
If ( adjust_timestep(tstep,hstep,rstep,mndis,mxdis,mxstp,config%natms,config%parts,&
 xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
            Call info(message,.true.)

  ! restore initial conditions

             Do irgd=1,rigid%n_types
                rigid%q0(irgd)=q0t(irgd)
                rigid%q1(irgd)=q1t(irgd)
                rigid%q2(irgd)=q2t(irgd)
                rigid%q3(irgd)=q3t(irgd)
             End Do

  ! restart vv1

             Go To 100
          End If
       End If

  ! second stage of velocity verlet algorithm

    Else

  ! update velocity of FPs

       Do j=1,config%nfree
          i=config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
             tmp=hstep/config%weight(i)
             config%vxx(i)=config%vxx(i)+tmp*config%parts(i)%fxx
             config%vyy(i)=config%vyy(i)+tmp*config%parts(i)%fyy
             config%vzz(i)=config%vzz(i)+tmp*config%parts(i)%fzz
          End If
       End Do

  ! RATTLE procedures
  ! apply velocity corrections to bond and PMF constraints

       If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
         Call apply_rattle(tstep,thermo%kit,pmf,cons,stat,domain,tmr,config,comm)
       End If

  ! Get RB COM stress and virial

       Call rigid_bodies_stress(strcom,ggx,ggy,ggz,rigid,config,comm)
       vircom=-(strcom(1)+strcom(5)+strcom(9))

  ! update velocity of RBs

       krgd=0
       Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) < lrgd) Then ! Not that it matters

  ! calculate COM force and torque

             fmx=0.0_wp ; fmy=0.0_wp ; fmz=0.0_wp
             tqx=0.0_wp ; tqy=0.0_wp ; tqz=0.0_wp
             Do jrgd=1,lrgd
                krgd=krgd+1

                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

  ! If the RB has a frozen particle then no net force

                If (rigid%frozen(0,rgdtyp) == 0) Then
                   fmx=fmx+config%parts(i)%fxx
                   fmy=fmy+config%parts(i)%fyy
                   fmz=fmz+config%parts(i)%fzz
                End If

                tqx=tqx+ggy(krgd)*config%parts(i)%fzz-ggz(krgd)*config%parts(i)%fyy
                tqy=tqy+ggz(krgd)*config%parts(i)%fxx-ggx(krgd)*config%parts(i)%fzz
                tqz=tqz+ggx(krgd)*config%parts(i)%fyy-ggy(krgd)*config%parts(i)%fxx
             End Do

  ! If the RB has 2+ frozen particles (ill=1) the net torque
  ! must align along the axis of rotation

             If (rigid%frozen(0,rgdtyp) > 1) Then
                i1=rigid%index_local(rigid%index_global(1,rgdtyp),irgd)
                i2=rigid%index_local(rigid%index_global(2,rgdtyp),irgd)

                x(1)=config%parts(i1)%xxx-config%parts(i2)%xxx
                y(1)=config%parts(i1)%yyy-config%parts(i2)%yyy
                z(1)=config%parts(i1)%zzz-config%parts(i2)%zzz

                Call images(config%imcon,config%cell,1,x,y,z)

                tmp=(x(1)*tqx+y(1)*tqy+z(1)*tqz)/(x(1)**2+y(1)**2+z(1)**2)
                tqx=x(1)*tmp
                tqy=y(1)*tmp
                tqz=z(1)*tmp
             End If

  ! current rotation matrix

             Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

  ! calculate torque in principal frame

             trx=tqx*rot(1)+tqy*rot(4)+tqz*rot(7)
             try=tqx*rot(2)+tqy*rot(5)+tqz*rot(8)
             trz=tqx*rot(3)+tqy*rot(6)+tqz*rot(9)

  ! calculate quaternion torques

             qt0=2.0_wp*(-rigid%q1(irgd)*trx-rigid%q2(irgd)*try-rigid%q3(irgd)*trz)
             qt1=2.0_wp*( rigid%q0(irgd)*trx-rigid%q3(irgd)*try+rigid%q2(irgd)*trz)
             qt2=2.0_wp*( rigid%q3(irgd)*trx+rigid%q0(irgd)*try-rigid%q1(irgd)*trz)
             qt3=2.0_wp*(-rigid%q2(irgd)*trx+rigid%q1(irgd)*try+rigid%q0(irgd)*trz)

  ! recover quaternion momenta at half time step

             opx=rigid%oxx(irgd)*rigid%rix(1,rgdtyp)
             opy=rigid%oyy(irgd)*rigid%riy(1,rgdtyp)
             opz=rigid%ozz(irgd)*rigid%riz(1,rgdtyp)

             p0=2.0_wp*(-rigid%q1(irgd)*opx-rigid%q2(irgd)*opy-rigid%q3(irgd)*opz)
             p1=2.0_wp*( rigid%q0(irgd)*opx-rigid%q3(irgd)*opy+rigid%q2(irgd)*opz)
             p2=2.0_wp*( rigid%q3(irgd)*opx+rigid%q0(irgd)*opy-rigid%q1(irgd)*opz)
             p3=2.0_wp*(-rigid%q2(irgd)*opx+rigid%q1(irgd)*opy+rigid%q0(irgd)*opz)

  ! update quaternion momenta to full step

             p0=p0+hstep*qt0
             p1=p1+hstep*qt1
             p2=p2+hstep*qt2
             p3=p3+hstep*qt3

  ! update RB angular & COM velocities to full step

             opx=0.5_wp*(-rigid%q1(irgd)*p0+rigid%q0(irgd)*p1+rigid%q3(irgd)*p2-rigid%q2(irgd)*p3)
             opy=0.5_wp*(-rigid%q2(irgd)*p0-rigid%q3(irgd)*p1+rigid%q0(irgd)*p2+rigid%q1(irgd)*p3)
             opz=0.5_wp*(-rigid%q3(irgd)*p0+rigid%q2(irgd)*p1-rigid%q1(irgd)*p2+rigid%q0(irgd)*p3)

             rigid%oxx(irgd)=opx*rigid%rix(2,rgdtyp)
             rigid%oyy(irgd)=opy*rigid%riy(2,rgdtyp)
             rigid%ozz(irgd)=opz*rigid%riz(2,rgdtyp)

             tmp=hstep/rigid%weight(0,rgdtyp)
             rigid%vxx(irgd)=rigid%vxx(irgd)+tmp*fmx
             rigid%vyy(irgd)=rigid%vyy(irgd)+tmp*fmy
             rigid%vzz(irgd)=rigid%vzz(irgd)+tmp*fmz

  ! update RB members velocities

             Do jrgd=1,lrgd
                If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                   i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                   If (i <= config%natms) Then
                      x(1)=rigid%x(jrgd,rgdtyp)
                      y(1)=rigid%y(jrgd,rgdtyp)
                      z(1)=rigid%z(jrgd,rgdtyp)

  ! new atomic velocities in body frame

                      vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                      vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                      vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

  ! new atomic velocities in lab frame

                      config%vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                      config%vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                      config%vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                   End If
                End If
             End Do

          End If
       End Do

  ! Andersen Thermostat
  !
  ! qualify non-shell, non-frozen free particles (n) for a random kick
  ! qualify RBs (r) by them and derive related shells (s)

  ! tpn(comm%idnode) number of thermostatted particles on this node (comm%idnode)
  ! ntp - grand total of non-shell, non-frozen particles to thermostat

       qn(1:config%natms)     = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
       qs(0:2,1:cshell%ntshl) = 0 ! unqualified core-shell unit with a local shell
       qr(1:rigid%n_types)     = 0 ! unqualified RB

       j = 0
       scale = tstep/thermo%tau_t
       Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. cshell%legshl(0,i) >= 0) Then
             If (sarurnd(seed,config%ltg(i),0,nstep) <= scale) Then
                j = j + 1
                qn(i) = 1
             End If
          End If
       End Do
       tpn(comm%idnode) = j
       Do i=0,comM%mxnode-1
         If (i /= comm%idnode) tpn(i) = 0
       End Do
       Call gsum(comm,tpn)
       ntp = Sum(tpn)

       If (ntp == 0) Go To 200

       j = 0 ! no qualified good RB (one qualified RB is enough to trigger all)
       Do i=1,matms
          If (qn(i) == 1) Then
             If (config%lfree(i) == 1) j = j + 1
          End If
       End Do
       Call gsum(comm,j)

  ! tpr(comm%idnode) number of thermostatted RB units on this node (comm%idnode)
  ! rtp - grand total of RB units to thermostat
  ! (can be larger than rigid%total due to sharing)

       k = 0
       If (j > 0) Then
          If (rigid%share) Then
             qn(config%natms+1:config%nlast) = 0 ! refresh the q array for shared RB units
             Call update_shared_units_int(config,rigid%list_shared, &
               rigid%map_shared,qn,domain,comm)
          End If

          j = 0
          Do irgd=1,rigid%n_types
             rgdtyp=rigid%list(0,irgd)

  ! For all good RBs

             lrgd=rigid%list(-1,irgd)
             If (rigid%frozen(0,rgdtyp) < lrgd) Then
                Do jrgd=1,lrgd
                   i=rigid%index_local(jrgd,irgd) ! local index of particle/site
                   If (qn(i) == 1) Then
                      If (qr(irgd) == 0) Then ! An overall hit is registered
                          qr(irgd) = 1
                          j = j + 1
                      End If

                      If (i <= config%natms) tpn(comm%idnode) = tpn(comm%idnode) - 1 ! Less free particles are hit
                   End If
                End Do

                If (qr(irgd) == 1) Then ! accounting for a random kick on the RB
                   i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                   i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum
                   If (rigid%frozen(0,rgdtyp) == 0) Then
                      If (i1 <= config%natms) k = k + 1
                   End If
                   If (i2 <= config%natms) k = k + 1
                End If
             End If
          End Do
       End If
  ! tpn(comm%idnode) number of thermostatted free particles on this node (comm%idnode)
  ! ntp - grand total of non-shell, non-frozen free particles to thermostat
       Do i=0,comm%mxnode-1
          If (i /= comm%idnode) tpn(i) = 0
       End Do
       Call gsum(comm,tpn)
       ntp = Sum(tpn)
       tpr(comm%idnode) = j
       Do i=0,comm%mxnode-1
         If (i /= comm%idnode) tpr(i) = 0
       End Do
       Call gsum(comm,tpr)
       rtp = Sum(tpr)

  ! Get gaussian distribution

       tkin = 0.0_wp
       mxdr = 0.0_wp
       Do i=1,config%natms
          If (qn(i) == 1 .and. config%lfree(i) == 0) Then
             If (sites%dof_site(config%lsite(i)) > zero_plus) mxdr = mxdr + sites%dof_site(config%lsite(i))

  ! Get gaussian distribution (unit variance)

             Call box_mueller_saru3(seed,config%ltg(i),nstep,xxt(i),yyt(i),zzt(i))

  ! Get scaler to target variance/Sqrt(config%weight)

             tmp = 1.0_wp/Sqrt(config%weight(i))
             xxt(i) = xxt(i)*tmp
             yyt(i) = yyt(i)*tmp
             zzt(i) = zzt(i)*tmp

             tkin = tkin + config%weight(i)*(xxt(i)**2+yyt(i)**2+zzt(i)**2)
          End If
       End Do

       If (rtp > 0) Then
          Do irgd=1,rigid%n_types
             If (qr(irgd) == 1) Then
                rgdtyp=rigid%list(0,irgd)

                lrgd=rigid%list(-1,irgd)
                Do jrgd=1,lrgd
                   i=rigid%index_local(jrgd,irgd) ! particle index
                   If (i <= config%natms) Then
                      If (sites%dof_site(config%lsite(i)) > zero_plus) mxdr = mxdr + sites%dof_site(config%lsite(i))
                   End If
                End Do

                i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                If (rigid%frozen(0,rgdtyp) == 0 .and. i1 <= config%natms) Then

  ! Get gaussian distribution (unit variance)

                   Call box_mueller_saru3(seed,config%ltg(i1),nstep,xxt(i1),yyt(i1),zzt(i1))

  ! Get scaler to target variance/Sqrt(config%weight)

                   tmp = 1.0_wp/Sqrt(rigid%weight(0,rgdtyp))

                   xxt(i1) = xxt(i1)*tmp
                   yyt(i1) = yyt(i1)*tmp
                   zzt(i1) = zzt(i1)*tmp

                   tkin = tkin + rigid%weight(0,rgdtyp)*(xxt(i1)**2+yyt(i1)**2+zzt(i1)**2)

                End If

                If (i2 <= config%natms) Then

  ! Get gaussian distribution (unit variance)

                   Call box_mueller_saru3(seed,config%ltg(i2),nstep,xxt(i2),yyt(i2),zzt(i2))

  ! Get scaler to target variance/Sqrt(config%weight) -
  ! 3 different reciprocal moments of inertia

                   xxt(i2) = xxt(i2)*Sqrt(rigid%rix(2,rgdtyp))
                   yyt(i2) = yyt(i2)*Sqrt(rigid%riy(2,rgdtyp))
                   zzt(i2) = zzt(i2)*Sqrt(rigid%riz(2,rgdtyp))

                   tkin = tkin + (rigid%rix(1,rgdtyp)*xxt(i2)**2 + &
                     rigid%riy(1,rgdtyp)*yyt(i2)**2 + &
                     rigid%riz(1,rgdtyp)*zzt(i2)**2)

                End If
             End If
          End Do
       End If

       Call gsum(comm,tkin)
       If (tkin <= zero_plus) tkin = 1.0_wp
       Call gsum(comm,mxdr)

  ! Scale to target temperature and apply thermostat

       scale = Sqrt(mxdr * boltz * thermo%temp / tkin)
       tmp = Sqrt(1.0_wp-thermo%soft**2)*scale

       j = 0
       Do i=1,config%natms
          If (qn(i) == 1 .and. config%lfree(i) == 0) Then
             If (thermo%soft <= zero_plus) Then ! New target velocity
                config%vxx(i) = xxt(i)*scale
                config%vyy(i) = yyt(i)*scale
                config%vzz(i) = zzt(i)*scale
             Else ! Softened velocity (mixture between old & new)
                config%vxx(i) = thermo%soft*config%vxx(i) + tmp*xxt(i)
                config%vyy(i) = thermo%soft*config%vyy(i) + tmp*yyt(i)
                config%vzz(i) = thermo%soft*config%vzz(i) + tmp*zzt(i)
             End If
          End If
       End Do

       If (rtp > 0) Then

  ! Update shared RBs' velocities

          If (rigid%share) Then
            Call update_shared_units(config,rigid%list_shared, &
              rigid%map_shared,xxt,yyt,zzt,domain,comm)
          End If

  ! calculate new RBs' COM and angular velocities

          Do irgd=1,rigid%n_types
             If (qr(irgd) == 1) Then
                rgdtyp=rigid%list(0,irgd)

                i1=rigid%index_local(1,irgd) ! particle to bare the random RB COM momentum
                i2=rigid%index_local(2,irgd) ! particle to bare the random RB angular momentum

                If (rigid%frozen(0,rgdtyp) == 0) Then
                   If (thermo%soft <= zero_plus) Then ! New target velocity
                      rigid%vxx(irgd) = xxt(i1)*scale
                      rigid%vyy(irgd) = yyt(i1)*scale
                      rigid%vzz(irgd) = zzt(i1)*scale
                   Else ! Softened velocity (mixture between old & new)
                      rigid%vxx(irgd) = thermo%soft*rigid%vxx(irgd) + tmp*xxt(i1)
                      rigid%vyy(irgd) = thermo%soft*rigid%vyy(irgd) + tmp*yyt(i1)
                      rigid%vzz(irgd) = thermo%soft*rigid%vzz(irgd) + tmp*zzt(i1)
                   End If
                End If

                If (thermo%soft <= zero_plus) Then ! New target velocity
                   rigid%oxx(irgd) = xxt(i2)*scale
                   rigid%oyy(irgd) = yyt(i2)*scale
                   rigid%ozz(irgd) = zzt(i2)*scale
                Else ! Softened velocity (mixture between old & new)
                   rigid%oxx(irgd) = thermo%soft*rigid%oxx(irgd) + tmp*xxt(i2)
                   rigid%oyy(irgd) = thermo%soft*rigid%oyy(irgd) + tmp*yyt(i2)
                   rigid%ozz(irgd) = thermo%soft*rigid%ozz(irgd) + tmp*zzt(i2)
                End If

  ! get new rotation matrix

                Call getrotmat(rigid%q0(irgd),rigid%q1(irgd),rigid%q2(irgd),rigid%q3(irgd),rot)

  ! update RB members new velocities

                lrgd=rigid%list(-1,irgd)
                Do jrgd=1,lrgd
                   If (rigid%frozen(jrgd,rgdtyp) == 0) Then
                      i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                      If (i <= config%natms) Then
                         x(1)=rigid%x(jrgd,rgdtyp)
                         y(1)=rigid%y(jrgd,rgdtyp)
                         z(1)=rigid%z(jrgd,rgdtyp)

  ! new atomic velocities in body frame

                         vpx=rigid%oyy(irgd)*z(1)-rigid%ozz(irgd)*y(1)
                         vpy=rigid%ozz(irgd)*x(1)-rigid%oxx(irgd)*z(1)
                         vpz=rigid%oxx(irgd)*y(1)-rigid%oyy(irgd)*x(1)

  ! new atomic velocities in lab frame

                         config%vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                         config%vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                         config%vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
                      End If
                   End If
                End Do
             End If
          End Do
       End If

  ! tps(comm%idnode) number of thermostatted core-shell units on this node (comm%idnode)
  ! stp - grand total of core-shell units to thermostat

       j = 0
       If (cshell%keyshl == SHELL_ADIABATIC) Then
          If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
             qn(config%natms+1:config%nlast) = 0
             Call update_shared_units_int(config,cshell%lishp_shl, &
               cshell%lashp_shl,qn,domain,comm)
          End If

          If (cshell%ntshl > 0) Then
             Do k=1,cshell%ntshl
                i1=local_index(cshell%listshl(1,k),matms,config%lsi,config%lsa)
                i2=local_index(cshell%listshl(2,k),matms,config%lsi,config%lsa)

                If (qn(i1) == 1 .and. i2 > 0 .and. i2 <= config%natms) Then
                   j = j + 1

                   qs(0,k)=1
                   qs(1,k)=i1
                   qs(2,k)=i2
                End If
             End Do
          End If
       End If
       tps(comm%idnode) = j
       Do i=0,comm%mxnode-1
          If (i /= comm%idnode) tps(i) = 0
       End Do
       Call gsum(comm,tps)
       stp = Sum(tps)

  ! Thermalise the shells on hit cores

       If (stp > 0) Then
          If (cshell%lshmv_shl) Then
            Call update_shared_units(config,cshell%lishp_shl, &
              cshell%lashp_shl,config%vxx,config%vyy,config%vzz,domain,comm)
          End If

          If (tps(comm%idnode) > 0) Then
             Do k=1,cshell%ntshl
                If (qs(0,k) == 1) Then
                   i1=qs(1,k)
                   i2=qs(2,k)

                   config%vxx(i2)=config%vxx(i1)
                   config%vyy(i2)=config%vyy(i1)
                   config%vzz(i2)=config%vzz(i1)
                End If
             End Do
          End If
       End If

  ! remove system centre of mass velocity (random momentum walk)

       Call getvom(vom,config%vxx,config%vyy,config%vzz,rigid,config,comm)

       Do j=1,config%nfree
          i=config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
             config%vxx(i) = config%vxx(i) - vom(1)
             config%vyy(i) = config%vyy(i) - vom(2)
             config%vzz(i) = config%vzz(i) - vom(3)
          End If
       End Do

       Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

          If (rigid%frozen(0,rgdtyp) == 0) Then
             rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
             rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
             rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

             lrgd=rigid%list(-1,irgd)
             Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= config%natms) Then
                   config%vxx(i) = config%vxx(i) - vom(1)
                   config%vyy(i) = config%vyy(i) - vom(2)
                   config%vzz(i) = config%vzz(i) - vom(3)
                End If
             End Do
          End If
       End Do

  200 Continue

  ! update kinetic energy and stress

       Call kinstresf(config%vxx,config%vyy,config%vzz,strknf,config,comm)
       Call kinstrest(rigid,strknt,comm)

       strkin=strknf+strknt
       engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

  ! update rotational energy

       engrot=getknr(rigid,comm)

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
       Deallocate (lstitr,            Stat=fail( 1))
       Call cons%deallocate_work()
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
    Deallocate (qn,tpn,               Stat=fail(15))
    Deallocate (qs,tps,               Stat=fail(16))
    Deallocate (qr,tpr,               Stat=fail(17))
    If (Any(fail > 0)) Then
       Write(message,'(a)') 'nvt_a1 deallocation failure'
       Call error(0,message)
    End If

  End Subroutine nvt_a1_vv
End Module nvt_anderson
