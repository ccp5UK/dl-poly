Module nst_nose_hoover

  Use comms,           Only: comms_type
  Use configuration,   Only: configuration_type,&
                             getcom
  Use constants,       Only: boltz
  Use constraints,     Only: apply_rattle,&
                             apply_shake,&
                             constraints_tags,&
                             constraints_type
  Use core_shell,      Only: core_shell_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             info
  Use kinds,           Only: li,STR_LEN,&
                             wp
  Use kinetics,        Only: getvom,&
                             kinstresf,&
                             kinstress,&
                             kinstrest
  Use numerics,        Only: dcell,&
                             images,&
                             mat_mul
  Use nvt_nose_hoover, Only: nvt_h0_scl,&
                             nvt_h1_scl
  Use pmf,             Only: pmf_tags,&
                             pmf_type
  Use rigid_bodies,    Only: getrotmat,&
                             no_squish,&
                             rigid_bodies_stress,&
                             rigid_bodies_type
  Use site,            Only: site_type
  Use statistics,      Only: stats_type
  Use thermostat,      Only: CONSTRAINT_NONE,&
                             CONSTRAINT_SEMI_ORTHORHOMBIC,&
                             CONSTRAINT_SURFACE_AREA,&
                             CONSTRAINT_SURFACE_TENSION,&
                             VV_FIRST_STAGE,&
                             adjust_timestep,&
                             thermostat_type
  Use timer,           Only: timer_type
  Use vdw,             Only: vdw_type

  Implicit None

  Private

  Public :: nst_h0_vv, nst_h1_vv, nst_h0_scl, nst_h1_scl

Contains

  Subroutine nst_h0_vv(stage, lvar, mndis, mxdis, mxstp, tstep, &
                       degfre, stress, &
                       consv, &
                       strkin, engke, &
                       cshell, cons, pmf, stat, thermo, sites, &
                       vdws, domain, tmr, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet with Nose-Hoover thermostat and
    ! barostat (anisotropic pressure control) (symplectic)
    !
    ! Parrinello-Rahman type: changing config%cell shape
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
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: lvar
    Real(Kind=wp),            Intent(In   ) :: mndis, mxdis, mxstp
    Real(Kind=wp),            Intent(InOut) :: tstep
    Integer(Kind=li),         Intent(In   ) :: degfre
    Real(Kind=wp),            Intent(In   ) :: stress(1:9)
    Real(Kind=wp),            Intent(  Out) :: consv
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), engke
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(site_type),          Intent(InOut) :: sites
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Real(Kind=wp), Parameter :: uni(1:9) = (/1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp&
                                , 0.0_wp, 1.0_wp/)

    Character(Len=STR_LEN)         :: message
    Integer                    :: fail(1:9), i, iter
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: aaa(1:9), bbb(1:9), cell0(1:9), celprp(1:10), chit0, chpzr, &
                                  cint0, com(1:3), eta0(1:9), hstep, mxdr, qstep, rstep, str(1:9), &
                                  tmp, vir, vom(1:3), vzero, xt, yt, zt
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), oxt(:), oyt(:), ozt(:), vxt(:), vyt(:), &
                                  vzt(:), xxt(:), yyt(:), zzt(:)

! uni is the diagonal unit matrix

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms), Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms), oyt(1:config%mxatms), ozt(1:config%mxatms), Stat=fail(6))
    End If
    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail(7))
    Allocate (vxt(1:config%mxatms), vyt(1:config%mxatms), vzt(1:config%mxatms), Stat=fail(8))
    Allocate (fxt(1:config%mxatms), fyt(1:config%mxatms), fzt(1:config%mxatms), Stat=fail(9))
    If (Any(fail > 0)) Then
      Call error_alloc('pos/vel/force arrays', 'nst_h0_vv')
    End If

    If (thermo%newjob_0) Then
      thermo%newjob_0 = .false.

      ! store initial values of volume, long range corrections and density

      thermo%volm0 = config%volm
      thermo%elrc0 = vdws%elrc
      thermo%virlrc0 = vdws%vlrc

      Allocate (thermo%dens0(1:sites%mxatyp), Stat=fail(1))
      If (fail(1) > 0) Then
        Call error_alloc('thermo%dens0', 'nst_h0_vv')
      End If

      Do i = 1, sites%ntype_atom
        thermo%dens0(i) = sites%dens(i)
      End Do

      ! Sort thermo%eta for thermo%iso /= CONSTRAINT_NONE
      ! Initialise and get thermo%h_z for orthorhombic constraints

      thermo%h_z = 0
      Select Case (thermo%iso)
      Case (CONSTRAINT_SURFACE_AREA)
        thermo%eta(1:8) = 0.0_wp
      Case (CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC)
        thermo%eta(2:4) = 0.0_wp
        thermo%eta(6:8) = 0.0_wp

        Call dcell(config%cell, celprp)
        thermo%h_z = celprp(9)
      End Select

      ! inertia parameters for Nose-Hoover thermostat and barostat

      thermo%qmass = 2.0_wp * thermo%sigma * thermo%tau_t**2
      tmp = 2.0_wp * thermo%sigma / (boltz * Real(degfre, wp))

      Select Case (thermo%iso)
      Case (CONSTRAINT_NONE)
        thermo%ceng = 2.0_wp * thermo%sigma + 3.0_wp**2 * boltz * tmp
      Case (CONSTRAINT_SURFACE_AREA)
        thermo%ceng = 2.0_wp * thermo%sigma + 1.0_wp * boltz * tmp
      Case (CONSTRAINT_SURFACE_TENSION)
        thermo%ceng = 2.0_wp * thermo%sigma + 3.0_wp * boltz * tmp
      Case (CONSTRAINT_SEMI_ORTHORHOMBIC)
        thermo%ceng = 2.0_wp * thermo%sigma + 2.0_wp * boltz * tmp
      End Select

      thermo%pmass = ((2.0_wp * thermo%sigma + 3.0_wp * boltz * tmp) / 3.0_wp) * thermo%tau_p**2

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( &
                     thermo%eta(1)**2 + 2 * thermo%eta(2)**2 + 2 * thermo%eta(3)**2 + &
                     thermo%eta(5)**2 + 2 * thermo%eta(6)**2 + thermo%eta(9)**2)

      ! set number of constraint+pmf shake iterations and general iteration cycles

      thermo%mxiter = 1
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        thermo%mxkit = 1
        thermo%mxiter = thermo%mxiter + 3
      End If
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0) Then
        Call constraints_tags(lstitr, cons, config, comm)
      End If

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr, pmf, config, comm)
      End If
    End If

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    qstep = 0.5_wp * hstep
    rstep = 1.0_wp / tstep

    ! first pass of velocity verlet algorithm

    If (stage == VV_FIRST_STAGE) Then

      ! store initial values

      Do i = 1, config%natms
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

      cell0 = config%cell
      vzero = config%volm
      chit0 = thermo%chi_t
      cint0 = thermo%cint
      eta0 = thermo%eta
      chpzr = thermo%chip0

      ! calculate system centre of mass

      Call getcom(config, com, comm)

      100 Continue

      ! constraint virial and stress tensor

      If (cons%megcon > 0) Then
        stat%vircon = 0.0_wp
        stat%strcon = 0.0_wp
      End If

      ! PMF virial and stress tensor

      If (pmf%megpmf > 0) Then
        stat%virpmf = 0.0_wp
        stat%strpmf = 0.0_wp
      End If

      Do iter = 1, thermo%mxiter

        ! integrate and apply nvt_h0_scl thermostat - 1/4 step

        Call nvt_h0_scl(qstep, thermo%pmass, thermo%chip0, engke, thermo, config, comm)

        ! constraint+pmf virial and stress

        vir = stat%vircon + stat%virpmf
        str = stat%strcon + stat%strpmf

        ! integrate and apply nst_h0_scl barostat - 1/2 step

        Call nst_h0_scl(0, hstep, degfre, thermo%chi_t, str, stress, &
                        strkin, engke, thermo, config, comm)

        ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

        thermo%chip0 = Sqrt( &
                       thermo%eta(1)**2 + 2 * thermo%eta(2)**2 + 2 * thermo%eta(3)**2 + &
                       thermo%eta(5)**2 + 2 * thermo%eta(6)**2 + thermo%eta(9)**2)

        ! integrate and apply nvt_h0_scl thermostat - 1/4 step

        Call nvt_h0_scl(qstep, thermo%pmass, thermo%chip0, engke, thermo, config, comm)

        ! update velocities

        Do i = 1, config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
            tmp = hstep / config%weight(i)
            config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
            config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
            config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
          End If
        End Do

        ! scale config%cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

        aaa = tstep * thermo%eta
        Call mat_mul(aaa, aaa, bbb)
        aaa = uni + aaa + 0.5_wp * bbb
        Call mat_mul(aaa, cell0, config%cell)

        ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

        thermo%chi_p = thermo%eta(1) + thermo%eta(5) + thermo%eta(9)
        config%volm = config%volm * Exp(tstep * thermo%chi_p)
        thermo%chi_p = thermo%chi_p / 3.0_wp

        ! update positions: second order taylor expansion of Exp(tstep*thermo%eta)

        Do i = 1, config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
            xt = xxt(i) - com(1)
            yt = yyt(i) - com(2)
            zt = zzt(i) - com(3)
            config%parts(i)%xxx = tstep * config%vxx(i) + com(1) + xt * aaa(1) + yt * aaa(2) + zt * aaa(3)
            config%parts(i)%yyy = tstep * config%vyy(i) + com(2) + xt * aaa(2) + yt * aaa(5) + zt * aaa(6)
            config%parts(i)%zzz = tstep * config%vzz(i) + com(3) + xt * aaa(3) + yt * aaa(6) + zt * aaa(9)
          End If
        End Do

        ! SHAKE procedures

        If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_shake(tstep, oxt, oyt, ozt, &
                           lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
        End If

        ! restore original integration parameters as well as
        ! velocities if iter < thermo%mxiter
        ! in the next iteration stat%strcon and stat%strpmf are freshly new

        If (iter < thermo%mxiter) Then
          config%volm = vzero
          thermo%chi_t = chit0
          thermo%cint = cint0
          thermo%eta = eta0
          thermo%chip0 = chpzr

          Do i = 1, config%natms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
          End Do
        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, qstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, mxdr, comm)) Then
          Call info(message, .true.)

          ! restore initial conditions

          config%volm = vzero
          thermo%chi_t = chit0
          thermo%cint = cint0
          thermo%eta = eta0
          thermo%chip0 = chpzr

          Do i = 1, config%natms
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

      tmp = (thermo%volm0 / config%volm)
      vdws%elrc = thermo%elrc0 * tmp
      vdws%vlrc = thermo%virlrc0 * tmp
      Do i = 1, sites%ntype_atom
        sites%dens(i) = thermo%dens0(i) * tmp
      End Do

      ! get thermo%h_z for orthorhombic constraints

      If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
        Call dcell(config%cell, celprp)
        thermo%h_z = celprp(9)
      End If

      ! second stage of velocity verlet algorithm

    Else

      ! update velocity

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
          config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
          config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! integrate and apply nvt_h0_scl thermostat - 1/4 step

      Call nvt_h0_scl(qstep, thermo%pmass, thermo%chip0, engke, thermo, config, comm)

      ! constraint+pmf virial and stress

      vir = stat%vircon + stat%virpmf
      str = stat%strcon + stat%strpmf

      ! integrate and apply nst_h0_scl barostat - 1/2 step

      Call nst_h0_scl(0, hstep, degfre, thermo%chi_t, str, stress, &
                      strkin, engke, thermo, config, comm)

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( &
                     thermo%eta(1)**2 + 2 * thermo%eta(2)**2 + 2 * thermo%eta(3)**2 + &
                     thermo%eta(5)**2 + 2 * thermo%eta(6)**2 + thermo%eta(9)**2)

      ! integrate and apply nvt_h0_scl thermostat - 1/4 step

      Call nvt_h0_scl(qstep, thermo%pmass, thermo%chip0, engke, thermo, config, comm)

      ! conserved quantity less kinetic and potential energy terms
      consv = 0.5_wp * thermo%qmass * thermo%chi_t**2 + 0.5_wp * thermo%pmass * thermo%chip0**2 + &
              thermo%ceng * thermo%cint + (thermo%press + sum(thermo%stress(1:9:4))/3.0_wp) * config%volm

      ! remove system centre of mass velocity

      Call getvom(vom, config, comm)

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstress(strkin, config, comm)
      engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Deallocate (lstitr, Stat=fail(1))
      Call cons%deallocate_work()
      Call pmf%deallocate_work()
      Deallocate (oxt, oyt, ozt, Stat=fail(6))
    End If
    Deallocate (xxt, yyt, zzt, Stat=fail(7))
    Deallocate (vxt, vyt, vzt, Stat=fail(8))
    Deallocate (fxt, fyt, fzt, Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nst_h0 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nst_h0_vv

  Subroutine nst_h1_vv(stage, lvar, mndis, mxdis, mxstp, tstep, &
                       degfre, degrot, stress, &
                       consv, &
                       strkin, strknf, strknt, engke, engrot, &
                       strcom, vircom, &
                       cshell, cons, pmf, stat, thermo, sites, &
                       vdws, rigid, domain, tmr, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, equations of motion in molecular dynamics
    ! - velocity verlet with Nose-Hoover thermostat and
    ! barostat (anisotropic pressure control) (symplectic)
    !
    ! Parrinello-Rahman type: changing config%cell shape
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
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    ! amended   - i.t.todorov november 2019 (RBs unsafe haloing)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: lvar
    Real(Kind=wp),            Intent(In   ) :: mndis, mxdis, mxstp
    Real(Kind=wp),            Intent(InOut) :: tstep
    Integer(Kind=li),         Intent(In   ) :: degfre, degrot
    Real(Kind=wp),            Intent(In   ) :: stress(1:9)
    Real(Kind=wp),            Intent(  Out) :: consv
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), strknf(1:9), strknt(1:9), engke, &
                                               engrot, strcom(1:9), vircom
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(site_type),          Intent(InOut) :: sites
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Real(Kind=wp), Parameter :: uni(1:9) = (/1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp&
                                , 0.0_wp, 1.0_wp/)

    Character(Len=STR_LEN)     :: message
    Integer                    :: fail(1:14), i, i1, i2, irgd, iter, j, jrgd, krgd, lrgd, matms, &
                                  rgdtyp
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: aaa(1:9), bbb(1:9), cell0(1:9), celprp(1:10), chit0, chpzr, &
                                  cint0, com(1:3), eta0(1:9), fmx, fmy, fmz, hstep, mxdr, opx, &
                                  opy, opz, p0, p1, p2, p3, qstep, qt0, qt1, qt2, qt3, rot(1:9), &
                                  rstep, str(1:9), tmp, tqx, tqy, tqz, trx, try, trz, vir, &
                                  vom(1:3), vpx, vpy, vpz, vzero, x(1:1), xt, y(1:1), yt, z(1:1), &
                                  zt
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), ggx(:), ggy(:), ggz(:), oxt(:), oyt(:), &
                                  ozt(:), q0t(:), q1t(:), q2t(:), q3t(:), rgdoxt(:), rgdoyt(:), &
                                  rgdozt(:), rgdvxt(:), rgdvyt(:), rgdvzt(:), rgdxxt(:), &
                                  rgdyyt(:), rgdzzt(:), vxt(:), vyt(:), vzt(:), xxt(:), yyt(:), &
                                  zzt(:)

! uni is the diagonal unit matrix

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%mxatms), Stat=fail(1))
      Call cons%allocate_work(config%mxatms)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%mxatms), oyt(1:config%mxatms), ozt(1:config%mxatms), Stat=fail(6))
    End If

    Allocate (ggx(1:rigid%max_list * rigid%max_rigid), &
              ggy(1:rigid%max_list * rigid%max_rigid), &
              ggz(1:rigid%max_list * rigid%max_rigid), &
              Stat=fail(7))
    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail(8))
    Allocate (vxt(1:config%mxatms), vyt(1:config%mxatms), vzt(1:config%mxatms), Stat=fail(9))
    Allocate (fxt(1:config%mxatms), fyt(1:config%mxatms), fzt(1:config%mxatms), Stat=fail(10))
    Allocate (q0t(1:rigid%max_rigid), &
              q1t(1:rigid%max_rigid), &
              q2t(1:rigid%max_rigid), &
              q3t(1:rigid%max_rigid), Stat=fail(11))
    Allocate (rgdxxt(1:rigid%max_rigid), &
              rgdyyt(1:rigid%max_rigid), &
              rgdzzt(1:rigid%max_rigid), Stat=fail(12))
    Allocate (rgdvxt(1:rigid%max_rigid), &
              rgdvyt(1:rigid%max_rigid), &
              rgdvzt(1:rigid%max_rigid), Stat=fail(13))
    Allocate (rgdoxt(1:rigid%max_rigid), &
              rgdoyt(1:rigid%max_rigid), &
              rgdozt(1:rigid%max_rigid), Stat=fail(14))
    If (Any(fail > 0)) Then
      Call error_alloc('work arrays', 'nst_h1_vv')
    End If

    If (thermo%newjob_1) Then
      thermo%newjob_1 = .false.

      ! store initial values of volume, long range corrections and density

      thermo%volm0 = config%volm
      thermo%elrc0 = vdws%elrc
      thermo%virlrc0 = vdws%vlrc

      Allocate (thermo%dens0(1:sites%mxatyp), Stat=fail(1))
      If (fail(1) > 0) Then
        Call error_alloc('thermo%dens0', 'nst_h1_vv')
      End If

      Do i = 1, sites%ntype_atom
        thermo%dens0(i) = sites%dens(i)
      End Do

      ! Sort thermo%eta for thermo%iso /= CONSTRAINT_NONE
      ! Initialise and get thermo%h_z for orthorhombic constraints

      thermo%h_z = 0
      Select Case (thermo%iso)
      Case (CONSTRAINT_SURFACE_AREA)
        thermo%eta(1:8) = 0.0_wp
      Case (CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC)
        thermo%eta(2:4) = 0.0_wp
        thermo%eta(6:8) = 0.0_wp

        Call dcell(config%cell, celprp)
        thermo%h_z = celprp(9)
      End Select

      ! inertia parameters for Nose-Hoover thermostat and barostat

      thermo%qmass = 2.0_wp * thermo%sigma * thermo%tau_t**2
      tmp = 2.0_wp * thermo%sigma / (boltz * Real(degfre, wp))
      Select Case (thermo%iso)
      Case (CONSTRAINT_NONE)
        thermo%ceng = 2.0_wp * thermo%sigma + 3.0_wp**2 * boltz * tmp
      Case (CONSTRAINT_SURFACE_AREA)
        thermo%ceng = 2.0_wp * thermo%sigma + 1.0_wp * boltz * tmp
      Case (CONSTRAINT_SURFACE_TENSION)
        thermo%ceng = 2.0_wp * thermo%sigma + 3.0_wp * boltz * tmp
      Case (CONSTRAINT_SEMI_ORTHORHOMBIC)
        thermo%ceng = 2.0_wp * thermo%sigma + 2.0_wp * boltz * tmp
      End Select
      thermo%pmass = ((Real(degfre - degrot, wp) + 3.0_wp) / 3.0_wp) * boltz * tmp * thermo%tau_p**2

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( &
                     thermo%eta(1)**2 + 2 * thermo%eta(2)**2 + 2 * thermo%eta(3)**2 + &
                     thermo%eta(5)**2 + 2 * thermo%eta(6)**2 + thermo%eta(9)**2)

      ! set number of constraint+pmf shake iterations and general iteration cycles

      thermo%mxiter = 1
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        thermo%mxkit = 1
        thermo%mxiter = thermo%mxiter + 3
      End If
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake

      ! thermo%unsafe positioning due to possibly locally shared RBs

      thermo%unsafe = Any(domain%map_unique > 0)
    End If

    ! set matms

    matms = config%nlast
    If (comm%mxnode == 1) matms = config%natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0) Then
        Call constraints_tags(lstitr, cons, config, comm)
      End If

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0) Then
        Call pmf_tags(lstitr, pmf, config, comm)
      End If

    End If

    ! Get the RB particles vectors wrt the RB's COM

    krgd = 0
    Do irgd = 1, rigid%n_types
      rgdtyp = rigid%list(0, irgd)

      ! For all good RBs

      lrgd = rigid%list(-1, irgd)
      If (rigid%frozen(0, rgdtyp) < lrgd) Then
        Do jrgd = 1, lrgd
          krgd = krgd + 1

          i = rigid%index_local(jrgd, irgd) ! local index of particle/site

          ! COM distances

          ggx(krgd) = config%parts(i)%xxx - rigid%xxx(irgd)
          ggy(krgd) = config%parts(i)%yyy - rigid%yyy(irgd)
          ggz(krgd) = config%parts(i)%zzz - rigid%zzz(irgd)
        End Do
      End If
    End Do

    ! minimum image convention for bond vectors

    Call images(config%imcon, config%cell, krgd, ggx, ggy, ggz)

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    qstep = 0.5_wp * hstep
    rstep = 1.0_wp / tstep

    ! first pass of velocity verlet algorithm

    If (stage == VV_FIRST_STAGE) Then

      ! store initial values

      Do i = 1, matms
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

      Do irgd = 1, rigid%n_types
        q0t(irgd) = rigid%q0(irgd)
        q1t(irgd) = rigid%q1(irgd)
        q2t(irgd) = rigid%q2(irgd)
        q3t(irgd) = rigid%q3(irgd)

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

      cell0 = config%cell
      vzero = config%volm
      chit0 = thermo%chi_t
      cint0 = thermo%cint
      eta0 = thermo%eta
      chpzr = thermo%chip0

      ! calculate system centre of mass

      Call getcom(config, com, comm)

      100 Continue

      ! constraint virial and stress tensor

      If (cons%megcon > 0) Then
        stat%vircon = 0.0_wp
        stat%strcon = 0.0_wp
      End If

      ! PMF virial and stress tensor

      If (pmf%megpmf > 0) Then
        stat%virpmf = 0.0_wp
        stat%strpmf = 0.0_wp
      End If

      Do iter = 1, thermo%mxiter

        ! integrate and apply nvt_h1_scl thermostat - 1/4 step

        Call nvt_h1_scl(qstep, thermo%pmass, thermo%chip0, engke, engrot, thermo, rigid, config, comm)

        ! constraint+pmf virial and stress

        vir = stat%vircon + stat%virpmf
        str = stat%strcon + stat%strpmf

        ! integrate and apply nst_h1_scl barostat - 1/2 step

        Call nst_h1_scl(0, hstep, degfre, degrot, thermo%chi_t, str, stress, strcom, &
                        strkin, strknf, strknt, engke, thermo, rigid, config, comm)

        ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

        thermo%chip0 = Sqrt( &
                       thermo%eta(1)**2 + 2 * thermo%eta(2)**2 + 2 * thermo%eta(3)**2 + &
                       thermo%eta(5)**2 + 2 * thermo%eta(6)**2 + thermo%eta(9)**2)

        ! integrate and apply nvt_h1_scl thermostat - 1/4 step

        Call nvt_h1_scl(qstep, thermo%pmass, thermo%chip0, engke, engrot, thermo, rigid, config, comm)

        ! update velocity of FPs

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
            tmp = hstep / config%weight(i)
            config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
            config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
            config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
          End If
        End Do

        ! scale config%cell vectors: second order taylor expansion of Exp(tstep*thermo%eta)

        aaa = tstep * thermo%eta
        Call mat_mul(aaa, aaa, bbb)
        aaa = uni + aaa + 0.5_wp * bbb
        Call mat_mul(aaa, cell0, config%cell)

        ! update volume and construct a 'mock' thermo%chi_p=Tr[thermo%eta]/3

        thermo%chi_p = thermo%eta(1) + thermo%eta(5) + thermo%eta(9)
        config%volm = config%volm * Exp(tstep * thermo%chi_p)
        thermo%chi_p = thermo%chi_p / 3.0_wp

        ! update position of FPs: second order taylor expansion of Exp(tstep*thermo%eta)

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
            xt = xxt(i) - com(1)
            yt = yyt(i) - com(2)
            zt = zzt(i) - com(3)
            config%parts(i)%xxx = tstep * config%vxx(i) + com(1) + xt * aaa(1) + yt * aaa(2) + zt * aaa(3)
            config%parts(i)%yyy = tstep * config%vyy(i) + com(2) + xt * aaa(2) + yt * aaa(5) + zt * aaa(6)
            config%parts(i)%zzz = tstep * config%vzz(i) + com(3) + xt * aaa(3) + yt * aaa(6) + zt * aaa(9)
          End If
        End Do

        ! SHAKE procedures

        If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_shake(tstep, oxt, oyt, ozt, &
                           lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
        End If

        ! restore original integration parameters as well as
        ! velocities if iter < thermo%mxiter
        ! in the next iteration stat%strcon and stat%strpmf are freshly new

        If (iter < thermo%mxiter) Then
          config%volm = vzero
          thermo%chi_t = chit0
          thermo%cint = cint0
          thermo%eta = eta0
          thermo%chip0 = chpzr

          Do j = 1, config%nfree
            i = config%lstfre(j)

            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
          End Do

          Do irgd = 1, rigid%n_types
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

      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        ! For all good RBs

        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then

          ! calculate COM force and torque

          fmx = 0.0_wp; fmy = 0.0_wp; fmz = 0.0_wp
          tqx = 0.0_wp; tqy = 0.0_wp; tqz = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force

            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + fxt(i)
              fmy = fmy + fyt(i)
              fmz = fmz + fzt(i)
            End If

            tqx = tqx + ggy(krgd) * fzt(i) - ggz(krgd) * fyt(i)
            tqy = tqy + ggz(krgd) * fxt(i) - ggx(krgd) * fzt(i)
            tqz = tqz + ggx(krgd) * fyt(i) - ggy(krgd) * fxt(i)
          End Do

          ! If the RB has 2+ frozen particles (ill=1) the net torque
          ! must align along the axis of rotation

          If (rigid%frozen(0, rgdtyp) > 1) Then
            i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
            i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

            x(1) = xxt(i1) - xxt(i2)
            y(1) = yyt(i1) - yyt(i2)
            z(1) = zzt(i1) - zzt(i2)

            Call images(config%imcon, cell0, 1, x, y, z)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) / (x(1)**2 + y(1)**2 + z(1)**2)
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp
          End If

          ! current rotation matrix

          Call getrotmat(q0t(irgd), q1t(irgd), q2t(irgd), q3t(irgd), rot)

          ! calculate torque in principal frame

          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          ! calculate quaternion torques

          qt0 = 2.0_wp * (-q1t(irgd) * trx - q2t(irgd) * try - q3t(irgd) * trz)
          qt1 = 2.0_wp * (q0t(irgd) * trx - q3t(irgd) * try + q2t(irgd) * trz)
          qt2 = 2.0_wp * (q3t(irgd) * trx + q0t(irgd) * try - q1t(irgd) * trz)
          qt3 = 2.0_wp * (-q2t(irgd) * trx + q1t(irgd) * try + q0t(irgd) * trz)

          ! recover quaternion momenta

          opx = rigid%oxx(irgd) * rigid%rix(1, rgdtyp)
          opy = rigid%oyy(irgd) * rigid%riy(1, rgdtyp)
          opz = rigid%ozz(irgd) * rigid%riz(1, rgdtyp)

          p0 = 2.0_wp * (-q1t(irgd) * opx - q2t(irgd) * opy - q3t(irgd) * opz)
          p1 = 2.0_wp * (q0t(irgd) * opx - q3t(irgd) * opy + q2t(irgd) * opz)
          p2 = 2.0_wp * (q3t(irgd) * opx + q0t(irgd) * opy - q1t(irgd) * opz)
          p3 = 2.0_wp * (-q2t(irgd) * opx + q1t(irgd) * opy + q0t(irgd) * opz)

          ! update quaternion momenta to half step

          p0 = p0 + hstep * qt0
          p1 = p1 + hstep * qt1
          p2 = p2 + hstep * qt2
          p3 = p3 + hstep * qt3

          ! rotate RB quaternions - update q to full timestep & amend p
          ! and get new rotation matrix

          Call no_squish &
            (tstep, rigid%rix(2, rgdtyp), rigid%riy(2, rgdtyp), rigid%riz(2, rgdtyp), &
             rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), p0, p1, p2, p3)
          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

          ! update RB angular & COM velocities to half step

          opx = 0.5_wp * (-rigid%q1(irgd) * p0 + rigid%q0(irgd) * p1 + rigid%q3(irgd) * p2 - rigid%q2(irgd) * p3)
          opy = 0.5_wp * (-rigid%q2(irgd) * p0 - rigid%q3(irgd) * p1 + rigid%q0(irgd) * p2 + rigid%q1(irgd) * p3)
          opz = 0.5_wp * (-rigid%q3(irgd) * p0 + rigid%q2(irgd) * p1 - rigid%q1(irgd) * p2 + rigid%q0(irgd) * p3)

          rigid%oxx(irgd) = opx * rigid%rix(2, rgdtyp)
          rigid%oyy(irgd) = opy * rigid%riy(2, rgdtyp)
          rigid%ozz(irgd) = opz * rigid%riz(2, rgdtyp)

          tmp = hstep / rigid%weight(0, rgdtyp)
          rigid%vxx(irgd) = rigid%vxx(irgd) + tmp * fmx
          rigid%vyy(irgd) = rigid%vyy(irgd) + tmp * fmy
          rigid%vzz(irgd) = rigid%vzz(irgd) + tmp * fmz

          ! update RB COM to full step

          xt = rgdxxt(irgd) - com(1)
          yt = rgdyyt(irgd) - com(2)
          zt = rgdzzt(irgd) - com(3)
          rigid%xxx(irgd) = tstep * rigid%vxx(irgd) + com(1) + xt * aaa(1) + yt * aaa(2) + zt * aaa(3)
          rigid%yyy(irgd) = tstep * rigid%vyy(irgd) + com(2) + xt * aaa(2) + yt * aaa(5) + zt * aaa(6)
          rigid%zzz(irgd) = tstep * rigid%vzz(irgd) + com(3) + xt * aaa(3) + yt * aaa(6) + zt * aaa(9)

          ! update RB members positions and halfstep velocities

          Do jrgd = 1, lrgd
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              If (rigid%frozen(jrgd, rgdtyp) == 0) Then
                x(1) = rigid%x(jrgd, rgdtyp)
                y(1) = rigid%y(jrgd, rgdtyp)
                z(1) = rigid%z(jrgd, rgdtyp)

                ! new atomic positions

                config%parts(i)%xxx = rot(1) * x(1) + rot(2) * y(1) + rot(3) * z(1) + rigid%xxx(irgd)
                config%parts(i)%yyy = rot(4) * x(1) + rot(5) * y(1) + rot(6) * z(1) + rigid%yyy(irgd)
                config%parts(i)%zzz = rot(7) * x(1) + rot(8) * y(1) + rot(9) * z(1) + rigid%zzz(irgd)

                ! new atomic velocities in body frame

                vpx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
                vpy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
                vpz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

                ! DD bound positions

                If (thermo%unsafe) Then
                  xt = xxt(i) - com(1)
                  yt = yyt(i) - com(2)
                  zt = zzt(i) - com(3)
                  config%vxx(i) = xt * aaa(1) + yt * aaa(2) + zt * aaa(3)
                  config%vyy(i) = xt * aaa(2) + yt * aaa(5) + zt * aaa(6)
                  config%vzz(i) = xt * aaa(3) + yt * aaa(6) + zt * aaa(9)

                  x(1) = (config%parts(i)%xxx - com(1)) - config%vxx(i)
                  y(1) = (config%parts(i)%yyy - com(2)) - config%vyy(i)
                  z(1) = (config%parts(i)%zzz - com(3)) - config%vzz(i)
                  Call images(config%imcon, config%cell, 1, x, y, z)
                  config%parts(i)%xxx = (x(1) + com(1)) + config%vxx(i)
                  config%parts(i)%yyy = (y(1) + com(2)) + config%vyy(i)
                  config%parts(i)%zzz = (z(1) + com(3)) + config%vzz(i)
                End If

                ! new atomic velocities in lab frame

                config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
                config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
                config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
              Else
                x(1) = rigid%xxx(irgd) - rgdxxt(irgd)
                y(1) = rigid%yyy(irgd) - rgdyyt(irgd)
                z(1) = rigid%zzz(irgd) - rgdzzt(irgd)
                If (thermo%unsafe) Call images(config%imcon, config%cell, 1, x, y, z) ! DD bound positions
                config%parts(i)%xxx = xxt(i) + x(1)
                config%parts(i)%yyy = yyt(i) + y(1)
                config%parts(i)%zzz = zzt(i) + z(1)
              End If
            End If
          End Do

        Else

          ! update RB COM to full step

          xt = rgdxxt(irgd) - com(1)
          yt = rgdyyt(irgd) - com(2)
          zt = rgdzzt(irgd) - com(3)
          rigid%xxx(irgd) = com(1) + xt * aaa(1) + yt * aaa(2) + zt * aaa(3)
          rigid%yyy(irgd) = com(2) + xt * aaa(2) + yt * aaa(5) + zt * aaa(6)
          rigid%zzz(irgd) = com(3) + xt * aaa(3) + yt * aaa(6) + zt * aaa(9)

          ! update RB members positions

          Do jrgd = 1, lrgd
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              x(1) = rigid%xxx(irgd) - rgdxxt(irgd)
              y(1) = rigid%yyy(irgd) - rgdyyt(irgd)
              z(1) = rigid%zzz(irgd) - rgdzzt(irgd)
              If (thermo%unsafe) Call images(config%imcon, config%cell, 1, x, y, z) ! DD bound positions
              config%parts(i)%xxx = xxt(i) + x(1)
              config%parts(i)%yyy = yyt(i) + y(1)
              config%parts(i)%zzz = zzt(i) + z(1)
            End If
          End Do

        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, qstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, mxdr, comm)) Then
          Call info(message, .true.)

          ! restore initial conditions

          config%volm = vzero
          thermo%chi_t = chit0
          thermo%cint = cint0
          thermo%eta = eta0
          thermo%chip0 = chpzr

          Do i = 1, matms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)

            config%parts(i)%fxx = fxt(i)
            config%parts(i)%fyy = fyt(i)
            config%parts(i)%fzz = fzt(i)
          End Do

          Do irgd = 1, rigid%n_types
            rigid%q0(irgd) = q0t(irgd)
            rigid%q1(irgd) = q1t(irgd)
            rigid%q2(irgd) = q2t(irgd)
            rigid%q3(irgd) = q3t(irgd)

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

      tmp = (thermo%volm0 / config%volm)
      vdws%elrc = thermo%elrc0 * tmp
      vdws%vlrc = thermo%virlrc0 * tmp
      Do i = 1, sites%ntype_atom
        sites%dens(i) = thermo%dens0(i) * tmp
      End Do

      ! get thermo%h_z for orthorhombic constraints

      If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
        Call dcell(config%cell, celprp)
        thermo%h_z = celprp(9)
      End If

      ! second stage of velocity verlet algorithm

    Else

      ! update velocity of FPs

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * config%parts(i)%fxx
          config%vyy(i) = config%vyy(i) + tmp * config%parts(i)%fyy
          config%vzz(i) = config%vzz(i) + tmp * config%parts(i)%fzz
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! Get RB COM stress and virial

      Call rigid_bodies_stress(strcom, ggx, ggy, ggz, rigid, config, comm)
      vircom = -(strcom(1) + strcom(5) + strcom(9))

      ! update velocity of RBs

      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        ! For all good RBs

        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then ! Not that it matters

          ! calculate COM force and torque

          fmx = 0.0_wp; fmy = 0.0_wp; fmz = 0.0_wp
          tqx = 0.0_wp; tqy = 0.0_wp; tqz = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force

            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + config%parts(i)%fxx
              fmy = fmy + config%parts(i)%fyy
              fmz = fmz + config%parts(i)%fzz
            End If

            tqx = tqx + ggy(krgd) * config%parts(i)%fzz - ggz(krgd) * config%parts(i)%fyy
            tqy = tqy + ggz(krgd) * config%parts(i)%fxx - ggx(krgd) * config%parts(i)%fzz
            tqz = tqz + ggx(krgd) * config%parts(i)%fyy - ggy(krgd) * config%parts(i)%fxx
          End Do

          ! If the RB has 2+ frozen particles (ill=1) the net torque
          ! must align along the axis of rotation

          If (rigid%frozen(0, rgdtyp) > 1) Then
            i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
            i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

            x(1) = config%parts(i1)%xxx - config%parts(i2)%xxx
            y(1) = config%parts(i1)%yyy - config%parts(i2)%yyy
            z(1) = config%parts(i1)%zzz - config%parts(i2)%zzz

            Call images(config%imcon, config%cell, 1, x, y, z)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) / (x(1)**2 + y(1)**2 + z(1)**2)
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp
          End If

          ! current rotation matrix

          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

          ! calculate torque in principal frame

          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          ! calculate quaternion torques

          qt0 = 2.0_wp * (-rigid%q1(irgd) * trx - rigid%q2(irgd) * try - rigid%q3(irgd) * trz)
          qt1 = 2.0_wp * (rigid%q0(irgd) * trx - rigid%q3(irgd) * try + rigid%q2(irgd) * trz)
          qt2 = 2.0_wp * (rigid%q3(irgd) * trx + rigid%q0(irgd) * try - rigid%q1(irgd) * trz)
          qt3 = 2.0_wp * (-rigid%q2(irgd) * trx + rigid%q1(irgd) * try + rigid%q0(irgd) * trz)

          ! recover quaternion momenta at half time step

          opx = rigid%oxx(irgd) * rigid%rix(1, rgdtyp)
          opy = rigid%oyy(irgd) * rigid%riy(1, rgdtyp)
          opz = rigid%ozz(irgd) * rigid%riz(1, rgdtyp)

          p0 = 2.0_wp * (-rigid%q1(irgd) * opx - rigid%q2(irgd) * opy - rigid%q3(irgd) * opz)
          p1 = 2.0_wp * (rigid%q0(irgd) * opx - rigid%q3(irgd) * opy + rigid%q2(irgd) * opz)
          p2 = 2.0_wp * (rigid%q3(irgd) * opx + rigid%q0(irgd) * opy - rigid%q1(irgd) * opz)
          p3 = 2.0_wp * (-rigid%q2(irgd) * opx + rigid%q1(irgd) * opy + rigid%q0(irgd) * opz)

          ! update quaternion momenta to full step

          p0 = p0 + hstep * qt0
          p1 = p1 + hstep * qt1
          p2 = p2 + hstep * qt2
          p3 = p3 + hstep * qt3

          ! update RB angular & COM velocities to full step

          opx = 0.5_wp * (-rigid%q1(irgd) * p0 + rigid%q0(irgd) * p1 + rigid%q3(irgd) * p2 - rigid%q2(irgd) * p3)
          opy = 0.5_wp * (-rigid%q2(irgd) * p0 - rigid%q3(irgd) * p1 + rigid%q0(irgd) * p2 + rigid%q1(irgd) * p3)
          opz = 0.5_wp * (-rigid%q3(irgd) * p0 + rigid%q2(irgd) * p1 - rigid%q1(irgd) * p2 + rigid%q0(irgd) * p3)

          rigid%oxx(irgd) = opx * rigid%rix(2, rgdtyp)
          rigid%oyy(irgd) = opy * rigid%riy(2, rgdtyp)
          rigid%ozz(irgd) = opz * rigid%riz(2, rgdtyp)

          tmp = hstep / rigid%weight(0, rgdtyp)
          rigid%vxx(irgd) = rigid%vxx(irgd) + tmp * fmx
          rigid%vyy(irgd) = rigid%vyy(irgd) + tmp * fmy
          rigid%vzz(irgd) = rigid%vzz(irgd) + tmp * fmz

          ! update RB members velocities

          Do jrgd = 1, lrgd
            If (rigid%frozen(jrgd, rgdtyp) == 0) Then
              i = rigid%index_local(jrgd, irgd) ! local index of particle/site

              If (i <= config%natms) Then
                x(1) = rigid%x(jrgd, rgdtyp)
                y(1) = rigid%y(jrgd, rgdtyp)
                z(1) = rigid%z(jrgd, rgdtyp)

                ! new atomic velocities in body frame

                vpx = rigid%oyy(irgd) * z(1) - rigid%ozz(irgd) * y(1)
                vpy = rigid%ozz(irgd) * x(1) - rigid%oxx(irgd) * z(1)
                vpz = rigid%oxx(irgd) * y(1) - rigid%oyy(irgd) * x(1)

                ! new atomic velocities in lab frame

                config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
                config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
                config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
              End If
            End If
          End Do

        End If
      End Do

      ! integrate and apply nvt_h1_scl thermostat - 1/4 step

      Call nvt_h1_scl(qstep, thermo%pmass, thermo%chip0, engke, engrot, thermo, rigid, config, comm)

      ! constraint+pmf virial and stress

      vir = stat%vircon + stat%virpmf
      str = stat%strcon + stat%strpmf

      ! integrate and apply nst_h1_scl barostat - 1/2 step

      Call nst_h1_scl(0, hstep, degfre, degrot, thermo%chi_t, str, stress, strcom, &
                      strkin, strknf, strknt, engke, thermo, rigid, config, comm)

      ! trace[thermo%eta*transpose(thermo%eta)] = trace[thermo%eta*thermo%eta]: thermo%eta is symmetric

      thermo%chip0 = Sqrt( &
                     thermo%eta(1)**2 + 2 * thermo%eta(2)**2 + 2 * thermo%eta(3)**2 + &
                     thermo%eta(5)**2 + 2 * thermo%eta(6)**2 + thermo%eta(9)**2)

      ! integrate and apply nvt_h1_scl thermostat - 1/4 step

      Call nvt_h1_scl(qstep, thermo%pmass, thermo%chip0, engke, engrot, thermo, rigid, config, comm)

      ! conserved quantity less kinetic and potential energy terms

      consv = 0.5_wp * thermo%qmass * thermo%chi_t**2 + 0.5_wp * thermo%pmass * thermo%chip0**2 + &
              thermo%ceng * thermo%cint + (thermo%press + sum(thermo%stress(1:9:4))/3.0_wp) * config%volm

      ! remove system centre of mass velocity

      Call getvom(vom, rigid, config, comm)

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        If (rigid%frozen(0, rgdtyp) == 0) Then
          rigid%vxx(irgd) = rigid%vxx(irgd) - vom(1)
          rigid%vyy(irgd) = rigid%vyy(irgd) - vom(2)
          rigid%vzz(irgd) = rigid%vzz(irgd) - vom(3)

          lrgd = rigid%list(-1, irgd)
          Do jrgd = 1, lrgd
            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            If (i <= config%natms) Then
              config%vxx(i) = config%vxx(i) - vom(1)
              config%vyy(i) = config%vyy(i) - vom(2)
              config%vzz(i) = config%vzz(i) - vom(3)
            End If
          End Do
        End If
      End Do

      ! update kinetic energy and stress

      Call kinstresf(strknf, config, comm)
      Call kinstrest(rigid, strknt, comm)

      strkin = strknf + strknt
      engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Deallocate (lstitr, Stat=fail(1))
      Call cons%deallocate_work()
      Call pmf%deallocate_work()
      Deallocate (oxt, oyt, ozt, Stat=fail(6))
    End If
    Deallocate (ggx, ggy, ggz, Stat=fail(7))
    Deallocate (xxt, yyt, zzt, Stat=fail(8))
    Deallocate (vxt, vyt, vzt, Stat=fail(9))
    Deallocate (fxt, fyt, fzt, Stat=fail(10))
    Deallocate (q0t, q1t, q2t, q3t, Stat=fail(11))
    Deallocate (rgdxxt, rgdyyt, rgdzzt, Stat=fail(12))
    Deallocate (rgdvxt, rgdvyt, rgdvzt, Stat=fail(13))
    Deallocate (rgdoxt, rgdoyt, rgdozt, Stat=fail(14))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nst_h1 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nst_h1_vv

  Subroutine nst_h0_scl(sw, tstep, degfre, chit, strcon, stress, strkin, engke, &
                        thermo, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to integrate and apply NsT barostat
    !
    ! sw=1 coupling to NVT thermostat for nst_m ensemble and
    !                                     additional scaling thermo%factor
    !
    ! sw=0 coupling to NVT thermostat for nst_h ensemble and
    !                                     no additional scaling thermo%factor
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
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: sw
    Real(Kind=wp),            Intent(In   ) :: tstep
    Integer(Kind=li),         Intent(In   ) :: degfre
    Real(Kind=wp),            Intent(In   ) :: chit, strcon(1:9), stress(1:9)
    Real(Kind=wp),            Intent(  Out) :: strkin(1:9), engke
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Real(Kind=wp), Parameter :: uni(1:9) = (/1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp&
                                , 0.0_wp, 1.0_wp/)

    Integer       :: i
    Real(Kind=wp) :: a1, a2, a3, a5, a6, a9, b1, b2, b3, b5, b6, b9, hstep, qstep, vxt, vyt, vzt

! uni is the diagonal unit matrix

    ! Initialise thermo%factor and 1/Nf for Nose-Hoover ensembles

    If (thermo%newjob_nst_scl_0) Then
      thermo%newjob_nst_scl_0 = .false.

      thermo%factor = 0.0_wp
      thermo%rf = 0.0_wp
      If (sw == 1) thermo%rf = 1.0_wp / Real(degfre, wp)
    End If

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    qstep = 0.5_wp * hstep

    ! thermostat thermo%eta to 1/4*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

    ! calculate kinetic contribution to stress tensor

    Call kinstress(strkin, config, comm)

    ! kinetic energy

    engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    ! barostat thermo%eta to 1/2*tstep

    If (sw == 1) thermo%factor = 2.0_wp * engke * thermo%rf

    ! split anisotropic from semi-isotropic barostats

    If (thermo%iso == CONSTRAINT_NONE) Then
      thermo%eta = thermo%eta + hstep * (strcon + stress + strkin + thermo%factor * uni - &
                                         (thermo%press * uni + thermo%stress) * config%volm) / thermo%pmass
    Else
      If (thermo%iso == CONSTRAINT_SURFACE_TENSION) Then
        thermo%eta(1) = thermo%eta(1) + hstep * (strcon(1) + stress(1) + strkin(1) + &
                                                 thermo%factor - (thermo%press + thermo%stress(1) - thermo%tension / thermo%h_z) * &
                                                 config%volm) / thermo%pmass
        thermo%eta(5) = thermo%eta(5) + hstep * (strcon(5) + stress(5) + strkin(5) + &
                                                 thermo%factor - (thermo%press + thermo%stress(5) - thermo%tension / thermo%h_z) * &
                                                 config%volm) / thermo%pmass
      Else If (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) Then
        thermo%eta(1) = 0.5_wp * (thermo%eta(1) + thermo%eta(5)) + hstep * (0.5_wp * &
                                         (strcon(1) + stress(1) + strkin(1) + strcon(5) + stress(5) + strkin(5)) + thermo%factor - &
                                   (thermo%press + 0.5_wp * (thermo%stress(1) + thermo%stress(5)) - thermo%tension / thermo%h_z) * &
                                                                            config%volm) / thermo%pmass
        thermo%eta(5) = thermo%eta(1)
      End If
      thermo%eta(9) = thermo%eta(9) + hstep * (strcon(9) + stress(9) + strkin(9) + &
                                               thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass
    End If

    ! thermostat thermo%eta to 2/4*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

    ! barostat the velocities to full 1*tstep
    ! second order taylor expansion of Exp[-tstep*(thermo%eta+thermo%factor*I)],
    ! where I is the unit tensor
    ! thermo%factor = Tr(thermo%eta)/Nf if sw=1, where Nf is degfre,
    ! else if sw=0 then thermo%factor=0, by default

    If (sw == 1) thermo%factor = (thermo%eta(1) + thermo%eta(5) + thermo%eta(9)) * thermo%rf

    a1 = -tstep * (thermo%eta(1) + thermo%factor)
    a2 = -tstep * thermo%eta(2)
    a3 = -tstep * thermo%eta(3)
    a5 = -tstep * (thermo%eta(5) + thermo%factor)
    a6 = -tstep * thermo%eta(6)
    a9 = -tstep * (thermo%eta(9) + thermo%factor)

    b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
    b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
    b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
    b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
    b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
    b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

    Do i = 1, config%natms
      vxt = config%vxx(i)
      vyt = config%vyy(i)
      vzt = config%vzz(i)

      config%vxx(i) = b1 * vxt + b2 * vyt + b3 * vzt
      config%vyy(i) = b2 * vxt + b5 * vyt + b6 * vzt
      config%vzz(i) = b3 * vxt + b6 * vyt + b9 * vzt
    End Do

    ! thermostat thermo%eta to 2/4*tstep

    ! calculate kinetic contribution to stress tensor

    Call kinstress(strkin, config, comm)

    ! kinetic energy

    engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    ! thermostat thermo%eta to 3/4*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

    ! barostat thermo%eta to full (2/2)*tstep

    If (sw == 1) thermo%factor = 2.0_wp * engke * thermo%rf

    ! split anisotropic from semi-isotropic barostats

    If (thermo%iso == CONSTRAINT_NONE) Then
      thermo%eta = thermo%eta + hstep * (strcon + stress + strkin + thermo%factor * uni - &
                                         (thermo%press * uni + thermo%stress) * config%volm) / thermo%pmass
    Else
      If (thermo%iso == CONSTRAINT_SURFACE_TENSION) Then
        thermo%eta(1) = thermo%eta(1) + hstep * (strcon(1) + stress(1) + strkin(1) + &
                       thermo%factor - (thermo%press + thermo%stress(1) - thermo%tension / thermo%h_z) * config%volm) / thermo%pmass
        thermo%eta(5) = thermo%eta(5) + hstep * (strcon(5) + stress(5) + strkin(5) + &
                       thermo%factor - (thermo%press + thermo%stress(5) - thermo%tension / thermo%h_z) * config%volm) / thermo%pmass
      Else If (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) Then
        thermo%eta(1) = 0.5_wp * (thermo%eta(1) + thermo%eta(5)) + hstep * (0.5_wp * &
                                         (strcon(1) + stress(1) + strkin(1) + strcon(5) + stress(5) + strkin(5)) + thermo%factor - &
                                   (thermo%press + 0.5_wp * (thermo%stress(1) + thermo%stress(5)) - thermo%tension / thermo%h_z) * &
                                                                            config%volm) / thermo%pmass
        thermo%eta(5) = thermo%eta(1)
      End If
      thermo%eta(9) = thermo%eta(9) + hstep * (strcon(9) + stress(9) + strkin(9) + &
                                               thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass
    End If

    ! thermostat thermo%eta to full (4/4)*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

  End Subroutine nst_h0_scl

  Subroutine nst_h1_scl(sw, tstep, degfre, degrot, chit, strcon, stress, strcom, &
                        strkin, strknf, strknt, engke, thermo, rigid, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to integrate and apply NsT barostat
    ! when singled RBs are present
    !
    ! sw=1 coupling to NVT thermostat for nst_m ensemble and
    !                                     additional scaling thermo%factor
    !
    ! sw=0 coupling to NVT thermostat for nst_h ensemble and
    !                                     no additional scaling thermo%factor
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
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: sw
    Real(Kind=wp),            Intent(In   ) :: tstep
    Integer(Kind=li),         Intent(In   ) :: degfre, degrot
    Real(Kind=wp),            Intent(In   ) :: chit, strcon(1:9), stress(1:9), strcom(1:9)
    Real(Kind=wp),            Intent(  Out) :: strkin(1:9), strknf(1:9), strknt(1:9), engke
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Real(Kind=wp), Parameter :: uni(1:9) = (/1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp&
                                , 0.0_wp, 1.0_wp/)

    Integer       :: i, irgd, j
    Real(Kind=wp) :: a1, a2, a3, a5, a6, a9, b1, b2, b3, b5, b6, b9, hstep, qstep, vxt, vyt, vzt

! initialise thermo%factor for Nose-Hoover ensembles
! uni is the diagonal unit matrix

    ! Initialise thermo%factor and 1/Nf for Nose-Hoover ensembles

    If (thermo%newjob_nst_scl_1) Then
      thermo%newjob_nst_scl_1 = .false.

      thermo%factor = 0.0_wp
      thermo%rf = 0.0_wp
      If (sw == 1) thermo%rf = 1.0_wp / Real(degfre - degrot, wp)
    End If

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    qstep = 0.5_wp * hstep

    ! thermostat thermo%eta to 1/4*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

    ! calculate kinetic contributions to stress tensor

    Call kinstresf(strknf, config, comm)
    Call kinstrest(rigid, strknt, comm)

    strkin = strknf + strknt

    ! kinetic energy

    engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    ! barostat thermo%eta to 1/2*tstep

    If (sw == 1) thermo%factor = 2.0_wp * engke * thermo%rf

    ! split anisotropic from semi-isotropic barostats

    Select Case (thermo%iso)
    Case (CONSTRAINT_NONE)
      thermo%eta = thermo%eta + hstep * (strcom + strcon + stress + strkin + thermo%factor * uni - &
        (thermo%press * uni + thermo%stress) * config%volm) / thermo%pmass

    Case (CONSTRAINT_SURFACE_TENSION)
      thermo%eta(1) = thermo%eta(1) + hstep * (strcom(1) + strcon(1) + stress(1) + strkin(1) + &
        thermo%factor - (thermo%press + thermo%stress(1) - thermo%tension / thermo%h_z) * &
        config%volm) / thermo%pmass
      thermo%eta(5) = thermo%eta(5) + hstep * (strcom(5) + strcon(5) + stress(5) + strkin(5) + &
        thermo%factor - (thermo%press + thermo%stress(5) - thermo%tension / thermo%h_z) * &
        config%volm) / thermo%pmass
      thermo%eta(9) = thermo%eta(9) + hstep * (strcom(9) + strcon(9) + stress(9) + strkin(9) + &
        thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass

    Case (CONSTRAINT_SEMI_ORTHORHOMBIC)
      thermo%eta(1) = 0.5_wp * (thermo%eta(1) + thermo%eta(5)) + hstep * (0.5_wp * &
        (strcom(1) + strcon(1) + stress(1) + strkin(1) + strcom(5) + strcon(5) + stress(5) + strkin(5)) + &
        thermo%factor - (thermo%press + 0.5_wp * (thermo%stress(1) + thermo%stress(5)) - &
        thermo%tension / thermo%h_z) * config%volm) / thermo%pmass
      thermo%eta(5) = thermo%eta(1)
      thermo%eta(9) = thermo%eta(9) + hstep * (strcom(9) + strcon(9) + stress(9) + strkin(9) + &
        thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass

    Case Default
      thermo%eta(9) = thermo%eta(9) + hstep * (strcom(9) + strcon(9) + stress(9) + strkin(9) + &
           thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass

    End Select

    ! thermostat thermo%eta to 2/4*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

    ! barostat the velocities to full 1*tstep
    ! second order taylor expansion of Exp[-tstep*(thermo%eta+thermo%factor*I)],
    ! where I is the unit tensor
    ! thermo%factor = Tr(thermo%eta)/Nf if sw=1, where Nf is degfre-degrot,
    ! else if sw=0 then thermo%factor=0, by default

    If (sw == 1) thermo%factor = (thermo%eta(1) + thermo%eta(5) + thermo%eta(9)) * thermo%rf

    a1 = -tstep * (thermo%eta(1) + thermo%factor)
    a2 = -tstep * thermo%eta(2)
    a3 = -tstep * thermo%eta(3)
    a5 = -tstep * (thermo%eta(5) + thermo%factor)
    a6 = -tstep * thermo%eta(6)
    a9 = -tstep * (thermo%eta(9) + thermo%factor)

    b1 = (a1 * a1 + a2 * a2 + a3 * a3) * 0.5_wp + a1 + 1.0_wp
    b2 = (a1 * a2 + a2 * a5 + a3 * a6) * 0.5_wp + a2
    b3 = (a1 * a3 + a2 * a6 + a3 * a9) * 0.5_wp + a3
    b5 = (a2 * a2 + a5 * a5 + a6 * a6) * 0.5_wp + a5 + 1.0_wp
    b6 = (a2 * a3 + a5 * a6 + a6 * a9) * 0.5_wp + a6
    b9 = (a3 * a3 + a6 * a6 + a9 * a9) * 0.5_wp + a9 + 1.0_wp

    Do j = 1, config%nfree
      i = config%lstfre(j)

      vxt = config%vxx(i)
      vyt = config%vyy(i)
      vzt = config%vzz(i)

      config%vxx(i) = b1 * vxt + b2 * vyt + b3 * vzt
      config%vyy(i) = b2 * vxt + b5 * vyt + b6 * vzt
      config%vzz(i) = b3 * vxt + b6 * vyt + b9 * vzt
    End Do

    Do irgd = 1, rigid%n_types
      vxt = rigid%vxx(irgd)
      vyt = rigid%vyy(irgd)
      vzt = rigid%vzz(irgd)

      rigid%vxx(irgd) = b1 * vxt + b2 * vyt + b3 * vzt
      rigid%vyy(irgd) = b2 * vxt + b5 * vyt + b6 * vzt
      rigid%vzz(irgd) = b3 * vxt + b6 * vyt + b9 * vzt
    End Do

    ! thermostat thermo%eta to 2/4*tstep

    ! calculate kinetic contributions to stress tensor

    Call kinstresf(strknf, config, comm)
    Call kinstrest(rigid, strknt, comm)

    strkin = strknf + strknt

    ! kinetic energy

    engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

    ! thermostat thermo%eta to 3/4*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

    ! barostat thermo%eta to full (2/2)*tstep

    If (sw == 1) thermo%factor = 2.0_wp * engke * thermo%rf

    ! split anisotropic from semi-isotropic barostats

    Select Case (thermo%iso)

    Case (CONSTRAINT_NONE)
      thermo%eta = thermo%eta + hstep * &
        (strcom + strcon + stress + strkin + thermo%factor * uni - &
        (thermo%press * uni + thermo%stress) * config%volm) / thermo%pmass

    Case (CONSTRAINT_SURFACE_TENSION)
      thermo%eta(1) = thermo%eta(1) + hstep * (strcom(1) + strcon(1) + stress(1) + strkin(1) + &
        thermo%factor - (thermo%press + thermo%stress(1) - thermo%tension / thermo%h_z) * &
        config%volm) / thermo%pmass
      thermo%eta(5) = thermo%eta(5) + hstep * (strcom(5) + strcon(5) + stress(5) + strkin(5) + &
        thermo%factor - (thermo%press + thermo%stress(5) - thermo%tension / thermo%h_z) * &
        config%volm) / thermo%pmass
      thermo%eta(9) = thermo%eta(9) + hstep * (strcom(9) + strcon(9) + stress(9) + strkin(9) + &
        thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass

    Case (CONSTRAINT_SEMI_ORTHORHOMBIC)
      thermo%eta(1) = 0.5_wp * (thermo%eta(1) + thermo%eta(5)) + hstep * (0.5_wp * &
        (strcom(1) + strcon(1) + stress(1) + strkin(1) + strcom(5) + strcon(5) + stress(5) + strkin(5)) + &
        thermo%factor - (thermo%press + 0.5_wp * (thermo%stress(1) + thermo%stress(5)) - &
        thermo%tension / thermo%h_z) * config%volm) / thermo%pmass
      thermo%eta(5) = thermo%eta(1)
      thermo%eta(9) = thermo%eta(9) + hstep * (strcom(9) + strcon(9) + stress(9) + strkin(9) + &
        thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass

    Case Default
      thermo%eta(9) = thermo%eta(9) + hstep * (strcom(9) + strcon(9) + stress(9) + strkin(9) + &
           thermo%factor - (thermo%press + thermo%stress(9)) * config%volm) / thermo%pmass

    End Select

    ! thermostat thermo%eta to full (4/4)*tstep

    thermo%eta = thermo%eta * Exp(-qstep * chit)

  End Subroutine nst_h1_scl
End Module nst_nose_hoover
