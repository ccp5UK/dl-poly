Module nst_mtk
  Use kinds,              Only : wp, li
  Use comms,              Only : comms_type
  Use constants,          Only : boltz  
  Use site, Only : site_type
  Use configuration,      Only : configuration_type
  Use domains,            Only : domains_type
  Use kinetics,           Only : getvom,kinstress,kinstresf,kinstrest
  Use constraints,        Only : apply_rattle,apply_shake,&
    constraints_tags, constraints_type
  Use pmf,                Only : pmf_tags,pmf_type
  Use rigid_bodies, Only : rigid_bodies_type,getrotmat,no_squish,rigid_bodies_stress
  Use nvt_nose_hoover,    Only : nvt_h0_scl,nvt_h1_scl
  Use nst_nose_hoover,    Only : nst_h0_scl,nst_h1_scl
  Use numerics,           Only : dcell, mat_mul,images
  Use core_shell, Only : core_shell_type
  Use statistics, Only : stats_type
  Use timer, Only : timer_type
  Use thermostat, Only : adjust_timestep,thermostat_type, &
    CONSTRAINT_NONE, CONSTRAINT_SURFACE_AREA, &
    CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC, &
    VV_FIRST_STAGE
  Use vdw, Only : vdw_type
  Use errors_warnings, Only : error,info
  Implicit None

  Private

  Public :: nst_m0_vv, nst_m1_vv

Contains

  Subroutine nst_m0_vv(stage,lvar,mndis,mxdis,mxstp,tstep, &
      degfre,stress,             &
      consv,                             &
      strkin,engke,                      &
      cshell,cons,pmf,stat,thermo,sites,&
      vdws,domain,tmr,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet with Nose-Hoover thermostat and
    ! barostat (anisotropic pressure control) and MTK coupling (symplectic)
    !
    ! Parrinello-Rahman type: changing config%cell shape
    !
    ! reference1: Martyna, Tuckerman, Tobias, Klein
    !             Mol. Phys., 1996, Vol. 87 (5), p. 1117
    ! reference2: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: stage
    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke
    Integer(Kind=li),   Intent( In    ) :: degfre
    Real( Kind = wp ),  Intent( In    ) :: stress(1:9)
    Real( Kind = wp ),  Intent(   Out ) :: consv
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( configuration_type ),   Intent( InOut ) :: config
    Type( comms_type ), Intent( InOut ) :: comm

    Integer                 :: fail(1:9),iter,i
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: chit0,cint0,chpzr,eta0(1:9)
    Real( Kind = wp )       :: cell0(1:9),vzero,celprp(1:10)
    Real( Kind = wp )       :: vir,str(1:9),mxdr,tmp, &
      vom(1:3),aaa(1:9),bbb(1:9)

    ! uni is the diagonal unit matrix

    Real( Kind = wp ), Parameter :: &
      uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)


    Logical,           Allocatable :: lstitr(:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)


    Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
    Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
    Real( Kind = wp ), Allocatable :: fxt(:),fyt(:),fzt(:)


    Character ( Len = 256 )  ::  message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms),                                  Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms),oyt(1:config%mxatms),ozt(1:config%mxatms),         Stat=fail(6))
    End If
    Allocate (xxt(1:config%mxatms),yyt(1:config%mxatms),zzt(1:config%mxatms),            Stat=fail(7))
    Allocate (vxt(1:config%mxatms),vyt(1:config%mxatms),vzt(1:config%mxatms),            Stat=fail(8))
    Allocate (fxt(1:config%mxatms),fyt(1:config%mxatms),fzt(1:config%mxatms),            Stat=fail(9))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'nst_m0 allocation failure'
      Call error(0,message)
    End If


    If (thermo%newjob) Then
      thermo%newjob = .false.

      ! store initial values of volume, long range corrections and density

      thermo%volm0   = config%volm
      thermo%elrc0   = vdws%elrc
      thermo%virlrc0 = vdws%vlrc

      Allocate (thermo%dens0(1:sites%mxatyp), Stat=fail(1))
      If (fail(1) > 0) Then
        Write(message,'(a)') 'dens0 allocation failure'
        Call error(0,message)
      End If
      Do i=1,sites%ntype_atom
        thermo%dens0(i) = sites%dens(i)
      End Do

      ! Sort thermo%eta for thermo%iso /= CONSTRAINT_NONE
      ! Initialise and get h_z for orthorhombic constraints

      thermo%h_z=0
      If      (thermo%iso == CONSTRAINT_SURFACE_AREA) Then
        thermo%eta(1:8) = 0.0_wp
      Else If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION,CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
        thermo%eta(2:4) = 0.0_wp
        thermo%eta(6:8) = 0.0_wp

        Call dcell(config%cell,celprp)
        thermo%h_z=celprp(9)
      End If

      ! inertia parameters for Nose-Hoover thermostat and barostat

      thermo%qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
      tmp   = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
      If      (thermo%iso == CONSTRAINT_NONE) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 3.0_wp**2*boltz*tmp
      Else If (thermo%iso == CONSTRAINT_SURFACE_AREA) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 1.0_wp*boltz*tmp
      Else If (thermo%iso == CONSTRAINT_SURFACE_TENSION) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 3.0_wp*boltz*tmp
      Else If (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 2.0_wp*boltz*tmp
      End If
      thermo%pmass = ((2.0_wp*thermo%sigma + 3.0_wp*boltz*tmp)/3.0_wp)*thermo%tau_p**2

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 &
        + thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

      ! set number of constraint+pmf shake iterations and general iteration cycles

      thermo%mxiter=1
      If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
        thermo%mxkit=1
        thermo%mxiter=thermo%mxiter+3
      End If
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit=cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms)=.false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0)Then
        Call constraints_tags(lstitr,cons,config,comm)
      End If 

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0)Then
        Call pmf_tags(lstitr,pmf,config,comm)
      End If
    End If

    ! timestep derivatives

    hstep = 0.5_wp*tstep
    qstep = 0.5_wp*hstep
    rstep = 1.0_wp/tstep

    ! first pass of velocity verlet algorithm

    If (stage == VV_FIRST_STAGE) Then

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

      ! store current integration variables

      cell0=config%cell
      vzero=config%volm
      chit0=thermo%chi_t
      cint0=thermo%cint
      eta0 =thermo%eta
      chpzr=thermo%chip0

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

      Do iter=1,thermo%mxiter

        ! integrate and apply nvt_h0_scl thermostat - 1/4 step

        Call nvt_h0_scl(qstep,thermo%pmass,thermo%chip0,engke,thermo,config,comm)

        ! constraint+pmf virial and stress

        vir=stat%vircon+stat%virpmf
        str=stat%strcon+stat%strpmf

        ! integrate and apply nst_h0_scl barostat - 1/2 step

        Call nst_h0_scl(1,hstep,degfre,thermo%chi_t,str,stress, &
          strkin,engke,thermo,config,comm)

        ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

        thermo%chip0 = Sqrt( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
          thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

        ! integrate and apply nvt_h0_scl thermostat - 1/4 step

        Call nvt_h0_scl(qstep,thermo%pmass,thermo%chip0,engke,thermo,config,comm)

        ! update velocities

        Do i=1,config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
            tmp=hstep/config%weight(i)
            config%vxx(i)=config%vxx(i)+tmp*config%parts(i)%fxx
            config%vyy(i)=config%vyy(i)+tmp*config%parts(i)%fyy
            config%vzz(i)=config%vzz(i)+tmp*config%parts(i)%fzz
          End If
        End Do

        ! scale config%cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

        aaa=tstep*thermo%eta
        Call mat_mul(aaa,aaa,bbb)
        aaa=uni+aaa+0.5_wp*bbb
        Call mat_mul(aaa,cell0,config%cell)

        ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

        thermo%chi_p = thermo%eta(1)+thermo%eta(5)+thermo%eta(9)
        config%volm = config%volm*Exp(tstep*thermo%chi_p)
        thermo%chi_p = thermo%chi_p / 3.0_wp

        ! update positions: second order taylor expansion of Exp(tstep*thermo%eta)

        Do i=1,config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx=tstep*config%vxx(i)+xxt(i)*aaa(1)+yyt(i)*aaa(2)+zzt(i)*aaa(3)
            config%parts(i)%yyy=tstep*config%vyy(i)+xxt(i)*aaa(2)+yyt(i)*aaa(5)+zzt(i)*aaa(6)
            config%parts(i)%zzz=tstep*config%vzz(i)+xxt(i)*aaa(3)+yyt(i)*aaa(6)+zzt(i)*aaa(9)
          End If
        End Do

        ! SHAKE procedures

        If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_shake(tstep,oxt,oyt,ozt,&
            lstitr,stat,pmf,cons,domain,tmr,config,thermo,comm)
        End If

        ! restore original integration parameters as well as
        ! velocities if iter < mxiter
        ! in the next iteration stat%strcon and stat%strpmf are freshly new

        If (iter < thermo%mxiter) Then
          config%volm=vzero
          thermo%chi_t=chit0
          thermo%cint=cint0
          thermo%eta =eta0
          thermo%chip0=chpzr

          Do i=1,config%natms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
          End Do
        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If ( adjust_timestep(tstep,hstep,rstep,qstep,mndis,mxdis,mxstp,config%natms,config%parts,&
          xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
          Call info(message,.true.)

          ! restore initial conditions

          config%volm=vzero
          thermo%chi_t=chit0
          thermo%cint=cint0
          thermo%eta =eta0
          thermo%chip0=chpzr

          Do i=1,config%natms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)

            config%parts(i)%fxx = fxt(i)
            config%parts(i)%fyy = fyt(i)
            config%parts(i)%fzz = fzt(i)
          End Do

          ! restart vv1

          Go To 100
        End If
      End If

      ! adjust long range corrections and number density

      tmp=(thermo%volm0/config%volm)
      vdws%elrc=thermo%elrc0*tmp
      vdws%vlrc=thermo%virlrc0*tmp
      Do i=1,sites%ntype_atom
        sites%dens(i)=thermo%dens0(i)*tmp
      End Do

      ! get h_z for orthorhombic constraints

      If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION,CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
        Call dcell(config%cell,celprp)
        thermo%h_z=celprp(9)
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

      ! integrate and apply nvt_h0_scl thermostat - 1/4 step

      Call nvt_h0_scl(qstep,thermo%pmass,thermo%chip0,engke,thermo,config,comm)

      ! constraint+pmf virial and stress

      vir=stat%vircon+stat%virpmf
      str=stat%strcon+stat%strpmf

      ! integrate and apply nst_h0_scl barostat - 1/2 step

      Call nst_h0_scl(1,hstep,degfre,thermo%chi_t,str,stress, &
        strkin,engke,thermo,config,comm)

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
        thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

      ! integrate and apply nvt_h0_scl thermostat - 1/4 step

      Call nvt_h0_scl(qstep,thermo%pmass,thermo%chip0,engke,thermo,config,comm)

      ! conserved quantity less kinetic and potential energy terms
      consv = 0.5_wp*thermo%qmass*thermo%chi_t**2 + 0.5_wp*thermo%pmass*thermo%chip0**2 + &
        thermo%ceng*thermo%cint + thermo%press*config%volm


      ! remove system centre of mass velocity

      Call getvom(vom,config%vxx,config%vyy,config%vzz,config,comm)

      Do i=1,config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i)=config%vxx(i)-vom(1)
          config%vyy(i)=config%vyy(i)-vom(2)
          config%vzz(i)=config%vzz(i)-vom(3)
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstress(config%vxx,config%vyy,config%vzz,strkin,config,comm)
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
      Write(message,'(a)') 'nst_m0 deallocation failure'
      Call error(0,message)
    End If

  End Subroutine nst_m0_vv

  Subroutine nst_m1_vv(stage,lvar,mndis,mxdis,mxstp,tstep, &
      degfre,degrot,stress,      &
      consv,                             &
      strkin,strknf,strknt,engke,engrot, &
      strcom,vircom,                     &
      cshell,cons,pmf,stat,thermo,sites, &
      vdws,rigid,domain,tmr,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, equations of motion in molecular dynamics
    ! - velocity verlet with Nose-Hoover thermostat and
    ! barostat (anisotropic pressure control) and MTK coupling (symplectic)
    !
    ! Parrinello-Rahman type: changing config%cell shape
    !
    ! reference1: Martyna, Tuckerman, Tobias, Klein
    !             Mol. Phys., 1996, Vol. 87 (5), p. 1117
    ! reference2: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: stage
    Logical,            Intent( In    ) :: lvar
    Real( Kind = wp ),  Intent( In    ) :: mndis,mxdis,mxstp
    Real( Kind = wp ),  Intent( InOut ) :: tstep
    Integer(Kind=li),   Intent( In    ) :: degfre,degrot
    Real( Kind = wp ),  Intent( In    ) :: stress(1:9)
    Real( Kind = wp ),  Intent(   Out ) :: consv
    Real( Kind = wp ),  Intent( InOut ) :: strkin(1:9),engke, &
      strknf(1:9),strknt(1:9),engrot
    Real( Kind = wp ),  Intent( InOut ) :: strcom(1:9),vircom
    Type( stats_type), Intent( InOut ) :: stat
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( pmf_type ), Intent( InOut ) :: pmf
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( site_type ), Intent( InOut ) :: sites
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( timer_type ), Intent( InOut ) :: tmr
    Type( configuration_type ),   Intent( InOut ) :: config
    Type( comms_type ), Intent( InOut ) :: comm

    Integer                 :: fail(1:14),matms,iter,i,j,i1,i2, &
      irgd,jrgd,krgd,lrgd,rgdtyp
    Real( Kind = wp )       :: hstep,qstep,rstep
    Real( Kind = wp )       :: chit0,cint0,chpzr,eta0(1:9)
    Real( Kind = wp )       :: cell0(1:9),vzero,celprp(1:10)
    Real( Kind = wp )       :: vir,str(1:9),mxdr,tmp, &
      vom(1:3),aaa(1:9),bbb(1:9)
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

    Character ( Len = 256 )        :: message

    fail=0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms),                                  Stat=fail( 1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms),oyt(1:config%mxatms),ozt(1:config%mxatms),         Stat=fail(6))
    End If
    Allocate (ggx(1:rigid%max_list*rigid%max_rigid), &
      ggy(1:rigid%max_list*rigid%max_rigid), &
      ggz(1:rigid%max_list*rigid%max_rigid), &
      Stat=fail( 7))
    Allocate (xxt(1:config%mxatms),yyt(1:config%mxatms),zzt(1:config%mxatms),            Stat=fail( 8))
    Allocate (vxt(1:config%mxatms),vyt(1:config%mxatms),vzt(1:config%mxatms),            Stat=fail( 9))
    Allocate (fxt(1:config%mxatms),fyt(1:config%mxatms),fzt(1:config%mxatms),            Stat=fail(10))
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
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'nst_m1 allocation failure'
      Call error(0,message)
    End If


    !If (newjob) Then
    !   newjob = .false.
    If (thermo%newjob) Then
      thermo%newjob = .false.

      ! store initial values of volume, long range corrections and density

      thermo%volm0   = config%volm
      thermo%elrc0   = vdws%elrc
      thermo%virlrc0 = vdws%vlrc

      Allocate (thermo%dens0(1:sites%mxatyp), Stat=fail(1))
      If (fail(1) > 0) Then
        Write(message,'(a)') 'dens0 allocation failure'
        Call error(0,message)
      End If
      Do i=1,sites%ntype_atom
        thermo%dens0(i) = sites%dens(i)
      End Do

      ! Sort thermo%eta for thermo%iso /= CONSTRAINT_NONE
      ! Initialise and get h_z for orthorhombic constraints

      thermo%h_z=0
      If      (thermo%iso == CONSTRAINT_SURFACE_AREA) Then
        thermo%eta(1:8) = 0.0_wp
      Else If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION,CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
        thermo%eta(2:4) = 0.0_wp
        thermo%eta(6:8) = 0.0_wp

        Call dcell(config%cell,celprp)
        thermo%h_z=celprp(9)
      End If

      ! inertia parameters for Nose-Hoover thermostat and barostat

      thermo%qmass = 2.0_wp*thermo%sigma*thermo%tau_t**2
      tmp   = 2.0_wp*thermo%sigma / (boltz*Real(degfre,wp))
      If      (thermo%iso == CONSTRAINT_NONE) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 3.0_wp**2*boltz*tmp
      Else If (thermo%iso == CONSTRAINT_SURFACE_AREA) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 1.0_wp*boltz*tmp
      Else If (thermo%iso == CONSTRAINT_SURFACE_TENSION) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 3.0_wp*boltz*tmp
      Else If (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) Then
        thermo%ceng  = 2.0_wp*thermo%sigma + 2.0_wp*boltz*tmp
      End If
      thermo%pmass = ((Real(degfre-degrot,wp) + 3.0_wp)/3.0_wp)*boltz*tmp*thermo%tau_p**2

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
        thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

      ! set number of constraint+pmf shake iterations and general iteration cycles

      thermo%mxiter=1
      If (cons%megcon > 0 .or.  pmf%megpmf > 0) Then
        thermo%mxkit=1
        thermo%mxiter=thermo%mxiter+3
      End If
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit=cons%max_iter_shake

      ! unsafe positioning due to possibly locally shared RBs

      thermo%unsafe=(Any(domain%map == comm%idnode))
    End If

    ! set matms

    matms=config%nlast
    If (comm%mxnode == 1) matms=config%natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms)=.false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0)Then
        Call constraints_tags(lstitr,cons,config,comm)
      End If

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0)Then
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
    qstep = 0.5_wp*hstep
    rstep = 1.0_wp/tstep

    ! first pass of velocity verlet algorithm

    If (stage == VV_FIRST_STAGE) Then

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

      ! store current integration variables

      cell0=config%cell
      vzero=config%volm
      chit0=thermo%chi_t
      cint0=thermo%cint
      eta0 =thermo%eta
      chpzr=thermo%chip0

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

      Do iter=1,thermo%mxiter

        ! integrate and apply nvt_h1_scl thermostat - 1/4 step

        Call nvt_h1_scl(qstep,thermo%pmass,thermo%chip0,engke,engrot,thermo,rigid,config,comm)

        ! constraint+pmf virial and stress

        vir=stat%vircon+stat%virpmf
        str=stat%strcon+stat%strpmf

        ! integrate and apply nst_h1_scl barostat - 1/2 step

        Call nst_h1_scl(1,hstep,degfre,degrot,thermo%chi_t,str,stress,strcom, &
          strkin,strknf,strknt,engke,thermo,rigid,config,comm)

        ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

        thermo%chip0 = Sqrt( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
          thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

        ! integrate and apply nvt_h1_scl thermostat - 1/4 step

        Call nvt_h1_scl(qstep,thermo%pmass,thermo%chip0,engke,engrot,thermo,rigid,config,comm)

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

        ! scale config%cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

        aaa=tstep*thermo%eta
        Call mat_mul(aaa,aaa,bbb)
        aaa=uni+aaa+0.5_wp*bbb
        Call mat_mul(aaa,cell0,config%cell)

        ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

        thermo%chi_p = thermo%eta(1)+thermo%eta(5)+thermo%eta(9)
        config%volm = config%volm*Exp(tstep*thermo%chi_p)
        thermo%chi_p = thermo%chi_p / 3.0_wp

        ! update position of FPs: second order taylor expansion of Exp(tstep*thermo%eta)

        Do j=1,config%nfree
          i=config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx=tstep*config%vxx(i)+xxt(i)*aaa(1)+yyt(i)*aaa(2)+zzt(i)*aaa(3)
            config%parts(i)%yyy=tstep*config%vyy(i)+xxt(i)*aaa(2)+yyt(i)*aaa(5)+zzt(i)*aaa(6)
            config%parts(i)%zzz=tstep*config%vzz(i)+xxt(i)*aaa(3)+yyt(i)*aaa(6)+zzt(i)*aaa(9)
          End If
        End Do

        ! SHAKE procedures

        If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_shake(tstep,oxt,oyt,ozt,&
            lstitr,stat,pmf,cons,domain,tmr,config,thermo,comm)
        End If

        ! restore original integration parameters as well as
        ! velocities if iter < mxiter
        ! in the next iteration stat%strcon and stat%strpmf are freshly new

        If (iter < thermo%mxiter) Then
          config%volm=vzero
          thermo%chi_t=chit0
          thermo%cint=cint0
          thermo%eta =eta0
          thermo%chip0=chpzr

          Do j=1,config%nfree
            i=config%lstfre(j)

            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
          End Do

          Do irgd=1,rigid%n_types
            rigid%vxx(irgd) = rgdvxt(irgd)
            rigid%vyy(irgd) = rgdvyt(irgd)
            rigid%vzz(irgd) = rgdvzt(irgd)

            rigid%oxx(irgd) = rgdoxt(irgd)
            rigid%oyy(irgd) = rgdoyt(irgd)
            rigid%ozz(irgd) = rgdozt(irgd)
          End Do
        End If
      End Do

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

            Call images(config%imcon,cell0,1,x,y,z)

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

          opx=rigid%oxx(irgd)*rigid%rix(1,rgdtyp)
          opy=rigid%oyy(irgd)*rigid%riy(1,rgdtyp)
          opz=rigid%ozz(irgd)*rigid%riz(1,rgdtyp)

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
          rigid%vxx(irgd)=rigid%vxx(irgd)+tmp*fmx
          rigid%vyy(irgd)=rigid%vyy(irgd)+tmp*fmy
          rigid%vzz(irgd)=rigid%vzz(irgd)+tmp*fmz

          ! update RB COM to full step

          rigid%xxx(irgd) = tstep*rigid%vxx(irgd) + &
            rgdxxt(irgd)*aaa(1) + rgdyyt(irgd)*aaa(2) + rgdzzt(irgd)*aaa(3)
          rigid%yyy(irgd) = tstep*rigid%vyy(irgd) + &
            rgdxxt(irgd)*aaa(2) + rgdyyt(irgd)*aaa(5) + rgdzzt(irgd)*aaa(6)
          rigid%zzz(irgd) = tstep*rigid%vzz(irgd) + &
            rgdxxt(irgd)*aaa(3) + rgdyyt(irgd)*aaa(6) + rgdzzt(irgd)*aaa(9)

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
                  config%vxx(i)=xxt(i)*aaa(1)+yyt(i)*aaa(2)+zzt(i)*aaa(3)
                  config%vyy(i)=xxt(i)*aaa(2)+yyt(i)*aaa(5)+zzt(i)*aaa(6)
                  config%vzz(i)=xxt(i)*aaa(3)+yyt(i)*aaa(6)+zzt(i)*aaa(9)

                  x(1)=config%parts(i)%xxx-config%vxx(i)
                  y(1)=config%parts(i)%yyy-config%vyy(i)
                  z(1)=config%parts(i)%zzz-config%vzz(i)
                  Call images(config%imcon,config%cell,1,x,y,z)
                  config%parts(i)%xxx=x(1)+config%vxx(i)
                  config%parts(i)%yyy=y(1)+config%vyy(i)
                  config%parts(i)%zzz=z(1)+config%vzz(i)
                End If

                ! new atomic velocities in lab frame

                config%vxx(i)=rot(1)*vpx+rot(2)*vpy+rot(3)*vpz+rigid%vxx(irgd)
                config%vyy(i)=rot(4)*vpx+rot(5)*vpy+rot(6)*vpz+rigid%vyy(irgd)
                config%vzz(i)=rot(7)*vpx+rot(8)*vpy+rot(9)*vpz+rigid%vzz(irgd)
              Else
                x(1)=rigid%xxx(irgd)-rgdxxt(irgd)
                y(1)=rigid%yyy(irgd)-rgdyyt(irgd)
                z(1)=rigid%zzz(irgd)-rgdzzt(irgd)
                If (thermo%unsafe) Call images(config%imcon,config%cell,1,x,y,z) ! DD bound positions
                config%parts(i)%xxx=xxt(i)+x(1)
                config%parts(i)%yyy=yyt(i)+y(1)
                config%parts(i)%zzz=zzt(i)+z(1)
              End If
            End If
          End Do

        Else

          ! update RB COM to full step

          rigid%xxx(irgd) = rgdxxt(irgd)*aaa(1) + rgdyyt(irgd)*aaa(2) + rgdzzt(irgd)*aaa(3)
          rigid%yyy(irgd) = rgdxxt(irgd)*aaa(2) + rgdyyt(irgd)*aaa(5) + rgdzzt(irgd)*aaa(6)
          rigid%zzz(irgd) = rgdxxt(irgd)*aaa(3) + rgdyyt(irgd)*aaa(6) + rgdzzt(irgd)*aaa(9)

          ! update RB members positions

          Do jrgd=1,lrgd
            i=rigid%index_local(jrgd,irgd) ! local index of particle/site

            If (i <= config%natms) Then
              x(1)=rigid%xxx(irgd)-rgdxxt(irgd)
              y(1)=rigid%yyy(irgd)-rgdyyt(irgd)
              z(1)=rigid%zzz(irgd)-rgdzzt(irgd)
              If (thermo%unsafe) Call images(config%imcon,config%cell,1,x,y,z) ! DD bound positions
              config%parts(i)%xxx=xxt(i)+x(1)
              config%parts(i)%yyy=yyt(i)+y(1)
              config%parts(i)%zzz=zzt(i)+z(1)
            End If
          End Do

        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If ( adjust_timestep(tstep,hstep,rstep,qstep,mndis,mxdis,mxstp,config%natms,config%parts,&
          xxt,yyt,zzt,cshell%legshl,message,mxdr,comm)) Then 
          Call info(message,.true.)

          ! restore initial conditions

          config%volm=vzero
          thermo%chi_t=chit0
          thermo%cint=cint0
          thermo%eta =eta0
          thermo%chip0=chpzr

          Do i=1,matms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)

            config%parts(i)%fxx = fxt(i)
            config%parts(i)%fyy = fyt(i)
            config%parts(i)%fzz = fzt(i)
          End Do

          Do irgd=1,rigid%n_types
            rigid%q0(irgd)=q0t(irgd)
            rigid%q1(irgd)=q1t(irgd)
            rigid%q2(irgd)=q2t(irgd)
            rigid%q3(irgd)=q3t(irgd)

            rigid%vxx(irgd) = rgdvxt(irgd)
            rigid%vyy(irgd) = rgdvyt(irgd)
            rigid%vzz(irgd) = rgdvzt(irgd)

            rigid%oxx(irgd) = rgdoxt(irgd)
            rigid%oyy(irgd) = rgdoyt(irgd)
            rigid%ozz(irgd) = rgdozt(irgd)
          End Do

          ! restart vv1

          Go To 100
        End If
      End If

      ! adjust long range corrections and number density

      tmp=(thermo%volm0/config%volm)
      vdws%elrc=thermo%elrc0*tmp
      vdws%vlrc=thermo%virlrc0*tmp
      Do i=1,sites%ntype_atom
        sites%dens(i)=thermo%dens0(i)*tmp
      End Do

      ! get h_z for orthorhombic constraints

      If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION,CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
        Call dcell(config%cell,celprp)
        thermo%h_z=celprp(9)
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

      ! update RB members positions

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

      ! integrate and apply nvt_h1_scl thermostat - 1/4 step

      Call nvt_h1_scl(qstep,thermo%pmass,thermo%chip0,engke,engrot,thermo,rigid,config,comm)

      ! constraint+pmf virial and stress

      vir=stat%vircon+stat%virpmf
      str=stat%strcon+stat%strpmf

      ! integrate and apply nst_h1_scl barostat - 1/2 step

      Call nst_h1_scl(1,hstep,degfre,degrot,thermo%chi_t,str,stress,strcom, &
        strkin,strknf,strknt,engke,thermo,rigid,config,comm)

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( thermo%eta(1)**2 + 2*thermo%eta(2)**2 + 2*thermo%eta(3)**2 + &
        thermo%eta(5)**2 + 2*thermo%eta(6)**2 + thermo%eta(9)**2 )

      ! integrate and apply nvt_h1_scl thermostat - 1/4 step

      Call nvt_h1_scl(qstep,thermo%pmass,thermo%chip0,engke,engrot,thermo,rigid,config,comm)

      ! conserved quantity less kinetic and potential energy terms

      consv = 0.5_wp*thermo%qmass*thermo%chi_t**2 + 0.5_wp*thermo%pmass*thermo%chip0**2 + &
        thermo%ceng*thermo%cint + thermo%press*config%volm

      ! remove system centre of mass velocity

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

      ! update kinetic energy and stress

      Call kinstresf(config%vxx,config%vyy,config%vzz,strknf,config,comm)
      Call kinstrest(rigid,strknt,comm)

      strkin=strknf+strknt
      engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

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
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'nst_m1 deallocation failure'
      Call error(0,message)
    End If

  End Subroutine nst_m1_vv
End Module nst_mtk
