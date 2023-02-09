Module nvt_ekin
  Use comms,           Only: comms_type,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constraints,     Only: apply_rattle,&
                             apply_shake,&
                             constraints_tags,&
                             constraints_type
  Use core_shell,      Only: core_shell_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info
  Use kinds,           Only: wp,STR_LEN
  Use kinetics,        Only: kinstresf,&
                             kinstress,&
                             kinstrest
  Use numerics,        Only: images
  Use pmf,             Only: pmf_tags,&
                             pmf_type
  Use rigid_bodies,    Only: getrotmat,&
                             no_squish,&
                             rigid_bodies_stress,&
                             rigid_bodies_type
  Use statistics,      Only: stats_type
  Use thermostat,      Only: VV_FIRST_STAGE,&
                             adjust_timestep,&
                             thermostat_type
  Use timer,           Only: timer_type

  Implicit None

  Private

  Public :: nvt_e0_vv, nvt_e1_vv, nvt_e0_scl, nvt_e1_scl

  Interface nvt_e0_scl
    Module Procedure nvt_e0_scl_config
    Module Procedure nvt_e0_scl_arrays
  End Interface nvt_e0_scl

  Interface nvt_e1_scl
    Module Procedure nvt_e1_scl_config
    Module Procedure nvt_e1_scl_arrays
  End Interface nvt_e1_scl

Contains

  Subroutine nvt_e0_vv(stage, lvar, mndis, mxdis, mxstp, tstep, chit, strkin, engke, &
                       thermo, cshell, cons, pmf, stat, domain, tmr, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet with Gaussian temperature
    ! constraints (Ekin conservation, symplectic)
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
    Real(Kind=wp),            Intent(  Out) :: chit
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), engke
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)         :: message
    Integer                    :: fail(1:9), i
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: hstep, mxdr, rstep, tmp
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), oxt(:), oyt(:), ozt(:), vxt(:), vyt(:), &
                                  vzt(:), xxt(:), yyt(:), zzt(:)

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%nlast), Stat=fail(1))
      Call cons%allocate_work(config%nlast)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%nlast), oyt(1:config%nlast), ozt(1:config%nlast), Stat=fail(6))
    End If
    Allocate (xxt(1:config%nlast), yyt(1:config%nlast), zzt(1:config%nlast), Stat=fail(7))
    Allocate (vxt(1:config%nlast), vyt(1:config%nlast), vzt(1:config%nlast), Stat=fail(8))
    Allocate (fxt(1:config%nlast), fyt(1:config%nlast), fzt(1:config%nlast), Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_e0 allocation failure'
      Call error(0, message)
    End If

    If (thermo%newjob_0) Then
      thermo%newjob_0 = .false.

      ! set number of constraint+pmf shake iterations

      If (cons%megcon > 0 .or. pmf%megpmf > 0) thermo%mxkit = 1
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

      ! integrate and apply nvt_e0_scl thermostat - half a step

      Call nvt_e0_scl(1, hstep, fxt, fyt, fzt, config%vxx, config%vyy, config%vzz, chit, engke, config, comm)

      ! update velocity and position

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * fxt(i)
          config%vyy(i) = config%vyy(i) + tmp * fyt(i)
          config%vzz(i) = config%vzz(i) + tmp * fzt(i)

          config%parts(i)%xxx = xxt(i) + tstep * config%vxx(i)
          config%parts(i)%yyy = yyt(i) + tstep * config%vyy(i)
          config%parts(i)%zzz = zzt(i) + tstep * config%vzz(i)
        End If
      End Do

      ! SHAKE procedures

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep, oxt, oyt, ozt, &
                         lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
      End If

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, mxdr, comm)) Then
          Call info(message, .true.)

          ! restore initial conditions

          Do i = 1, config%natms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
          End Do

          ! restart vv1

          Go To 100
        End If
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

      ! integrate and apply nvt_e0_scl thermostat - half a step

      Call nvt_e0_scl(1, hstep, config, config%vxx, config%vyy, config%vzz, chit, engke, comm)

      ! kinetic contribution to stress tensor

      Call kinstress(strkin, config, comm)

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
      Write (message, '(a)') 'nvt_e0 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nvt_e0_vv

  Subroutine nvt_e1_vv(stage, lvar, mndis, mxdis, mxstp, tstep, chit, strkin, strknf, &
                       strknt, engke, engrot, strcom, vircom, &
                       thermo, cshell, cons, pmf, stat, rigid, domain, tmr, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, equations of motion in molecular dynamics
    ! - velocity verlet with Gaussian temperature constraints
    ! (Ekin conservation, symplectic)
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
    Real(Kind=wp),            Intent(  Out) :: chit
    Real(Kind=wp),            Intent(InOut) :: strkin(1:9), strknf(1:9), strknt(1:9), engke, &
                                               engrot, strcom(1:9), vircom
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stat
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)         :: message
    Integer                    :: fail(1:16), i, i1, i2, irgd, j, jrgd, krgd, lrgd, matms, rgdtyp
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: fmx, fmy, fmz, hstep, mxdr, opx, opy, opz, p0, p1, p2, p3, qt0, &
                                  qt1, qt2, qt3, rot(1:9), rstep, tmp, tqx, tqy, tqz, trx, try, &
                                  trz, vpx, vpy, vpz, x(1:1), y(1:1), z(1:1)
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), ggx(:), ggy(:), ggz(:), oxt(:), oyt(:), &
                                  ozt(:), q0t(:), q1t(:), q2t(:), q3t(:), rgdfxx(:), rgdfyy(:), &
                                  rgdfzz(:), rgdoxt(:), rgdoyt(:), rgdozt(:), rgdtxx(:), &
                                  rgdtyy(:), rgdtzz(:), rgdvxt(:), rgdvyt(:), rgdvzt(:), &
                                  rgdxxt(:), rgdyyt(:), rgdzzt(:), vxt(:), vyt(:), vzt(:), xxt(:), &
                                  yyt(:), zzt(:)

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%nlast), Stat=fail(1))
      Call cons%allocate_work(config%nlast)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%nlast), oyt(1:config%nlast), ozt(1:config%nlast), Stat=fail(6))
    End If
    Allocate (ggx(1:rigid%max_list * rigid%max_rigid), &
              ggy(1:rigid%max_list * rigid%max_rigid), &
              ggz(1:rigid%max_list * rigid%max_rigid), Stat=fail(7))
    Allocate (xxt(1:config%nlast), yyt(1:config%nlast), zzt(1:config%nlast), Stat=fail(8))
    Allocate (vxt(1:config%nlast), vyt(1:config%nlast), vzt(1:config%nlast), Stat=fail(9))
    Allocate (fxt(1:config%nlast), fyt(1:config%nlast), fzt(1:config%nlast), Stat=fail(10))
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
    Allocate (rgdfxx(1:rigid%max_rigid), &
              rgdfyy(1:rigid%max_rigid), &
              rgdfzz(1:rigid%max_rigid), Stat=fail(15))
    Allocate (rgdtxx(1:rigid%max_rigid), &
              rgdtyy(1:rigid%max_rigid), &
              rgdtzz(1:rigid%max_rigid), Stat=fail(16))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_e1 allocation failure'
      Call error(0, message)
    End If

    If (thermo%newjob_1) Then
      thermo%newjob_1 = .false.

      ! set number of constraint+pmf shake iterations

      If (cons%megcon > 0 .or. pmf%megpmf > 0) thermo%mxkit = 1
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
    rstep = 1.0_wp / tstep

    ! Get RB force and torque

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

        ! store COM force and torque

        rgdfxx(irgd) = fmx
        rgdfyy(irgd) = fmy
        rgdfzz(irgd) = fmz

        rgdtxx(irgd) = trx
        rgdtyy(irgd) = try
        rgdtzz(irgd) = trz

      Else

        rgdfxx(irgd) = 0.0_wp
        rgdfyy(irgd) = 0.0_wp
        rgdfzz(irgd) = 0.0_wp

        rgdtxx(irgd) = 0.0_wp
        rgdtyy(irgd) = 0.0_wp
        rgdtzz(irgd) = 0.0_wp

      End If
    End Do

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

      ! integrate and apply nvt_e1_scl thermostat - half a step

      Call nvt_e1_scl &
        (1, hstep, fxt, fyt, fzt, config%vxx, config%vyy, config%vzz, &
         rgdfxx, rgdfyy, rgdfzz, rgdtxx, rgdtyy, rgdtzz, &
         chit, engke, engrot, rigid, config, comm)

      ! update velocity and position of FPs

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * fxt(i)
          config%vyy(i) = config%vyy(i) + tmp * fyt(i)
          config%vzz(i) = config%vzz(i) + tmp * fzt(i)

          config%parts(i)%xxx = xxt(i) + tstep * config%vxx(i)
          config%parts(i)%yyy = yyt(i) + tstep * config%vyy(i)
          config%parts(i)%zzz = zzt(i) + tstep * config%vzz(i)
        End If
      End Do

      ! SHAKE procedures

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_shake(tstep, oxt, oyt, ozt, &
                         lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
      End If

      ! update velocity and position of RBs

      krgd = 0
      Do irgd = 1, rigid%n_types
        rgdtyp = rigid%list(0, irgd)

        ! For all good RBs

        lrgd = rigid%list(-1, irgd)
        If (rigid%frozen(0, rgdtyp) < lrgd) Then

          ! recover COM force and torque in principal frame

          fmx = rgdfxx(irgd)
          fmy = rgdfyy(irgd)
          fmz = rgdfzz(irgd)

          trx = rgdtxx(irgd)
          try = rgdtyy(irgd)
          trz = rgdtzz(irgd)

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

          rigid%xxx(irgd) = rgdxxt(irgd) + tstep * rigid%vxx(irgd)
          rigid%yyy(irgd) = rgdyyt(irgd) + tstep * rigid%vyy(irgd)
          rigid%zzz(irgd) = rgdzzt(irgd) + tstep * rigid%vzz(irgd)

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
                  x(1) = config%parts(i)%xxx - xxt(i)
                  y(1) = config%parts(i)%yyy - yyt(i)
                  z(1) = config%parts(i)%zzz - zzt(i)
                  Call images(config%imcon, config%cell, 1, x, y, z)
                  config%parts(i)%xxx = x(1) + xxt(i)
                  config%parts(i)%yyy = y(1) + yyt(i)
                  config%parts(i)%zzz = z(1) + zzt(i)
                End If

                ! new atomic velocities in lab frame

                config%vxx(i) = rot(1) * vpx + rot(2) * vpy + rot(3) * vpz + rigid%vxx(irgd)
                config%vyy(i) = rot(4) * vpx + rot(5) * vpy + rot(6) * vpz + rigid%vyy(irgd)
                config%vzz(i) = rot(7) * vpx + rot(8) * vpy + rot(9) * vpz + rigid%vzz(irgd)
              End If
            End If
          End Do

        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, mxdr, comm)) Then
          Call info(message, .true.)

          ! restore initial conditions

          Do i = 1, matms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
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

          ! recover COM force and torque in principal frame

          fmx = rgdfxx(irgd)
          fmy = rgdfyy(irgd)
          fmz = rgdfzz(irgd)

          trx = rgdtxx(irgd)
          try = rgdtyy(irgd)
          trz = rgdtzz(irgd)

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

          ! current rotation matrix

          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

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

      ! integrate and apply nvt_e1_scl thermostat - half a step

      Call nvt_e1_scl &
        (1, hstep, config, config%vxx, config%vyy, config%vzz, &
         rgdfxx, rgdfyy, rgdfzz, rgdtxx, rgdtyy, rgdtzz, &
         chit, engke, engrot, rigid, comm)

      ! kinetic contributions to stress tensor

      Call kinstresf(strknf, config, comm)
      Call kinstrest(rigid, strknt, comm)

      strkin = strknf + strknt

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
    Deallocate (rgdfxx, rgdfyy, rgdfzz, Stat=fail(15))
    Deallocate (rgdtxx, rgdtyy, rgdtzz, Stat=fail(16))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nvt_e1 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nvt_e1_vv

  Subroutine nvt_e0_scl_arrays(stage, tstep, fxx, fyy, fzz, vxx, vyy, vzz, chit, engke, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to integrate and apply NVT Evans thermostat
    ! when no singled rigid bodies are present
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                     Intent(In   ) :: stage
    Real(Kind=wp),               Intent(In   ) :: tstep
    Real(Kind=wp), Dimension(:), Intent(In   ) :: fxx, fyy, fzz
    Real(Kind=wp), Dimension(:), Intent(InOut) :: vxx, vyy, vzz
    Real(Kind=wp),               Intent(  Out) :: chit, engke
    Type(configuration_type),    Intent(InOut) :: config
    Type(comms_type),            Intent(InOut) :: comm

    Integer       :: i
    Real(Kind=wp) :: buffer(1:2), scale, vdotf

    engke = 0.0_wp
    vdotf = 0.0_wp
    Do i = 1, config%natms
      engke = engke + config%weight(i) * (vxx(i)**2 + vyy(i)**2 + vzz(i)**2)
      vdotf = vdotf + vxx(i) * fxx(i) + vyy(i) * fyy(i) + vzz(i) * fzz(i)
    End Do

    buffer(1) = engke
    buffer(2) = vdotf
    Call gsum(comm, buffer)
    engke = buffer(1)
    vdotf = buffer(2)

    ! velocity friction and temperature scaling coefficient
    ! for Evans thermostat at tstep

    chit = vdotf / engke

    ! get corrected energy

    engke = 0.5_wp * engke

    If (stage == VV_FIRST_STAGE) Return

    ! thermostat velocities

    scale = Exp(-chit * tstep)
    Do i = 1, config%natms
      vxx(i) = vxx(i) * scale
      vyy(i) = vyy(i) * scale
      vzz(i) = vzz(i) * scale
    End Do

    ! thermostat kinetic energy

    engke = engke * scale**2

  End Subroutine nvt_e0_scl_arrays

  Subroutine nvt_e0_scl_config(stage, tstep, config, vxx, vyy, vzz, chit, engke, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to integrate and apply NVT Evans thermostat
    ! when no singled rigid bodies are present
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                   Intent(In   ) :: stage
    Real(Kind=wp),                             Intent(In   ) :: tstep
    Type(configuration_type),                  Intent(In   ) :: config
    Real(Kind=wp), Dimension(1:config%nlast), Intent(InOut) :: vxx, vyy, vzz
    Real(Kind=wp),                             Intent(  Out) :: chit, engke
    Type(comms_type),                          Intent(InOut) :: comm

    Integer       :: i
    Real(Kind=wp) :: buffer(1:2), scale, vdotf

    engke = 0.0_wp
    vdotf = 0.0_wp
    Do i = 1, config%natms
      engke = engke + config%weight(i) * (vxx(i)**2 + vyy(i)**2 + vzz(i)**2)
      vdotf = vdotf + vxx(i) * config%parts(i)%fxx + vyy(i) * config%parts(i)%fyy + vzz(i) * config%parts(i)%fzz
    End Do

    buffer(1) = engke
    buffer(2) = vdotf
    Call gsum(comm, buffer)
    engke = buffer(1)
    vdotf = buffer(2)

    ! velocity friction and temperature scaling coefficient
    ! for Evans thermostat at tstep

    chit = vdotf / engke

    ! get corrected energy

    engke = 0.5_wp * engke

    If (stage == VV_FIRST_STAGE) Return

    ! thermostat velocities

    scale = Exp(-chit * tstep)
    Do i = 1, config%natms
      vxx(i) = vxx(i) * scale
      vyy(i) = vyy(i) * scale
      vzz(i) = vzz(i) * scale
    End Do

    ! thermostat kinetic energy

    engke = engke * scale**2

  End Subroutine nvt_e0_scl_config

  Subroutine nvt_e1_scl_arrays &
    (stage, tstep, fxx, fyy, fzz, vxx, vyy, vzz, &
     rgdfxx, rgdfyy, rgdfzz, rgdtxx, rgdtyy, rgdtzz, &
     chit, engke, engrot, rigid, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to integrate and apply NVT Evans thermostat
    ! when singled rigid bodies are present
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                     Intent(In   ) :: stage
    Real(Kind=wp),                               Intent(In   ) :: tstep
    Real(Kind=wp), Dimension(:),                 Intent(In   ) :: fxx, fyy, fzz
    Real(Kind=wp), Dimension(:),                 Intent(InOut) :: vxx, vyy, vzz
    Type(rigid_bodies_type),                     Intent(InOut) :: rigid
    Real(Kind=wp),                               Intent(  Out) :: engrot, engke, chit
    Real(Kind=wp), Dimension(1:rigid%max_rigid), Intent(In   ) :: rgdtzz, rgdtyy, rgdtxx, rgdfzz, &
                                                                  rgdfyy, rgdfxx
    Type(configuration_type),                    Intent(InOut) :: config
    Type(comms_type),                            Intent(InOut) :: comm

    Integer       :: i, irgd, j, lrgd, rgdtyp
    Real(Kind=wp) :: buffer(1:4), odott, scale, tmp, vdotf

    engke = 0.0_wp
    vdotf = 0.0_wp
    engrot = 0.0_wp
    odott = 0.0_wp

    ! Free particles

    Do j = 1, config%nfree
      i = config%lstfre(j)

      engke = engke + config%weight(i) * (vxx(i)**2 + vyy(i)**2 + vzz(i)**2)
      vdotf = vdotf + vxx(i) * fxx(i) + vyy(i) * fyy(i) + vzz(i) * fzz(i)
    End Do

    ! RBs

    Do irgd = 1, rigid%n_types
      rgdtyp = rigid%list(0, irgd)

      ! For all good RBs

      lrgd = rigid%list(-1, irgd)
      If (rigid%frozen(0, rgdtyp) < lrgd) Then
        tmp = Real(rigid%index_local(0, irgd), wp) / Real(lrgd, wp)

        If (rigid%frozen(0, rgdtyp) == 0) Then
          engke = engke + &
                  tmp * rigid%weight(0, rgdtyp) * (rigid%vxx(irgd)**2 + &
                                                   rigid%vyy(irgd)**2 + &
                                                   rigid%vzz(irgd)**2)
          vdotf = vdotf + &
                  tmp * (rigid%vxx(irgd) * rgdfxx(irgd) + &
                         rigid%vyy(irgd) * rgdfyy(irgd) + &
                         rigid%vzz(irgd) * rgdfzz(irgd))
        End If

        engrot = engrot + &
                 tmp * (rigid%rix(1, rgdtyp) * rigid%oxx(irgd)**2 + &
                        rigid%riy(1, rgdtyp) * rigid%oyy(irgd)**2 + &
                        rigid%riz(1, rgdtyp) * rigid%ozz(irgd)**2)
        odott = odott + &
                tmp * (rigid%oxx(irgd) * rgdtxx(irgd) + &
                       rigid%oyy(irgd) * rgdtyy(irgd) + &
                       rigid%ozz(irgd) * rgdtzz(irgd))
      End If
    End Do

    buffer(1) = engke
    buffer(2) = vdotf
    buffer(3) = engrot
    buffer(4) = odott
    Call gsum(comm, buffer)
    engke = buffer(1)
    vdotf = buffer(2)
    engrot = buffer(3)
    odott = buffer(4)

    ! velocity friction and temperature scaling coefficient
    ! for Evans thermostat at tstep

    chit = (vdotf + odott) / (engke + engrot)

    ! get corrected energies

    engke = 0.5_wp * engke
    engrot = 0.5_wp * engrot

    If (stage == VV_FIRST_STAGE) Return

    ! thermostat velocities

    scale = Exp(-chit * tstep)

    Do j = 1, config%nfree
      i = config%lstfre(j)

      vxx(i) = vxx(i) * scale
      vyy(i) = vyy(i) * scale
      vzz(i) = vzz(i) * scale
    End Do

    Do irgd = 1, rigid%n_types
      rigid%vxx(irgd) = rigid%vxx(irgd) * scale
      rigid%vyy(irgd) = rigid%vyy(irgd) * scale
      rigid%vzz(irgd) = rigid%vzz(irgd) * scale

      rigid%oxx(irgd) = rigid%oxx(irgd) * scale
      rigid%oyy(irgd) = rigid%oyy(irgd) * scale
      rigid%ozz(irgd) = rigid%ozz(irgd) * scale
    End Do

    ! thermostat kinetic and rotational energy

    tmp = scale**2
    engke = engke * tmp
    engrot = engrot * tmp

  End Subroutine nvt_e1_scl_arrays
  Subroutine nvt_e1_scl_config &
    (stage, tstep, config, vxx, vyy, vzz, &
     rgdfxx, rgdfyy, rgdfzz, rgdtxx, rgdtyy, rgdtzz, &
     chit, engke, engrot, rigid, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to integrate and apply NVT Evans thermostat
    ! when singled rigid bodies are present
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                     Intent(In   ) :: stage
    Real(Kind=wp),                               Intent(In   ) :: tstep
    Type(configuration_type),                    Intent(In   ) :: config
    Real(Kind=wp), Dimension(1:config%nlast),   Intent(InOut) :: vxx, vyy, vzz
    Type(rigid_bodies_type),                     Intent(InOut) :: rigid
    Real(Kind=wp),                               Intent(  Out) :: engrot, engke, chit
    Real(Kind=wp), Dimension(1:rigid%max_rigid), Intent(In   ) :: rgdtzz, rgdtyy, rgdtxx, rgdfzz, &
                                                                  rgdfyy, rgdfxx
    Type(comms_type),                            Intent(InOut) :: comm

    Integer       :: i, irgd, j, lrgd, rgdtyp
    Real(Kind=wp) :: buffer(1:4), odott, scale, tmp, vdotf

    engke = 0.0_wp
    vdotf = 0.0_wp
    engrot = 0.0_wp
    odott = 0.0_wp

    ! Free particles

    Do j = 1, config%nfree
      i = config%lstfre(j)

      engke = engke + config%weight(i) * (vxx(i)**2 + vyy(i)**2 + vzz(i)**2)
      vdotf = vdotf + vxx(i) * config%parts(i)%fxx + vyy(i) * config%parts(i)%fyy + vzz(i) * config%parts(i)%fzz
    End Do

    ! RBs

    Do irgd = 1, rigid%n_types
      rgdtyp = rigid%list(0, irgd)

      ! For all good RBs

      lrgd = rigid%list(-1, irgd)
      If (rigid%frozen(0, rgdtyp) < lrgd) Then
        tmp = Real(rigid%index_local(0, irgd), wp) / Real(lrgd, wp)

        If (rigid%frozen(0, rgdtyp) == 0) Then
          engke = engke + &
                  tmp * rigid%weight(0, rgdtyp) * (rigid%vxx(irgd)**2 + &
                                                   rigid%vyy(irgd)**2 + &
                                                   rigid%vzz(irgd)**2)
          vdotf = vdotf + &
                  tmp * (rigid%vxx(irgd) * rgdfxx(irgd) + &
                         rigid%vyy(irgd) * rgdfyy(irgd) + &
                         rigid%vzz(irgd) * rgdfzz(irgd))
        End If

        engrot = engrot + &
                 tmp * (rigid%rix(1, rgdtyp) * rigid%oxx(irgd)**2 + &
                        rigid%riy(1, rgdtyp) * rigid%oyy(irgd)**2 + &
                        rigid%riz(1, rgdtyp) * rigid%ozz(irgd)**2)
        odott = odott + &
                tmp * (rigid%oxx(irgd) * rgdtxx(irgd) + &
                       rigid%oyy(irgd) * rgdtyy(irgd) + &
                       rigid%ozz(irgd) * rgdtzz(irgd))
      End If
    End Do

    buffer(1) = engke
    buffer(2) = vdotf
    buffer(3) = engrot
    buffer(4) = odott
    Call gsum(comm, buffer)
    engke = buffer(1)
    vdotf = buffer(2)
    engrot = buffer(3)
    odott = buffer(4)

    ! velocity friction and temperature scaling coefficient
    ! for Evans thermostat at tstep

    chit = (vdotf + odott) / (engke + engrot)

    ! get corrected energies

    engke = 0.5_wp * engke
    engrot = 0.5_wp * engrot

    If (stage == VV_FIRST_STAGE) Return

    ! thermostat velocities

    scale = Exp(-chit * tstep)

    Do j = 1, config%nfree
      i = config%lstfre(j)

      vxx(i) = vxx(i) * scale
      vyy(i) = vyy(i) * scale
      vzz(i) = vzz(i) * scale
    End Do

    Do irgd = 1, rigid%n_types
      rigid%vxx(irgd) = rigid%vxx(irgd) * scale
      rigid%vyy(irgd) = rigid%vyy(irgd) * scale
      rigid%vzz(irgd) = rigid%vzz(irgd) * scale

      rigid%oxx(irgd) = rigid%oxx(irgd) * scale
      rigid%oyy(irgd) = rigid%oyy(irgd) * scale
      rigid%ozz(irgd) = rigid%ozz(irgd) * scale
    End Do

    ! thermostat kinetic and rotational energy

    tmp = scale**2
    engke = engke * tmp
    engrot = engrot * tmp

  End Subroutine nvt_e1_scl_config
End Module nvt_ekin
