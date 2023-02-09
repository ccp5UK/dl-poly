Module nve
  Use comms,           Only: comms_type
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
  Use kinetics,        Only: getknr,&
                             kinstresf,&
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

  Public :: nve_0_vv, nve_1_vv

Contains

  Subroutine nve_0_vv(stage, lvar, mndis, mxdis, mxstp, tstep, &
                      strkin, engke, thermo, &
                      cshell, cons, pmf, stat, domain, tmr, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet (symplectic)
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
    Real(Kind=wp),            Intent(InOut) :: tstep, strkin(1:9), engke
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
      Allocate (lstitr(1:config%natms), Stat=fail(1))
      Call cons%allocate_work(config%nlast)
      Call pmf%allocate_work()
      Allocate (oxt(1:config%nlast), oyt(1:config%nlast), ozt(1:config%nlast), Stat=fail(6))
    End If
    Allocate (xxt(1:config%nlast), yyt(1:config%nlast), zzt(1:config%nlast), Stat=fail(7))
    Allocate (vxt(1:config%nlast), vyt(1:config%nlast), vzt(1:config%nlast), Stat=fail(8))
    Allocate (fxt(1:config%nlast), fyt(1:config%nlast), fzt(1:config%nlast), Stat=fail(9))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nve_0 allocation failure'
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

      ! update velocity and position

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = vxt(i) + tmp * fxt(i)
          config%vyy(i) = vyt(i) + tmp * fyt(i)
          config%vzz(i) = vzt(i) + tmp * fzt(i)

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
      Write (message, '(a)') 'nve_0 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nve_0_vv

  Subroutine nve_1_vv(stage, lvar, mndis, mxdis, mxstp, tstep, &
                      strkin, strknf, strknt, engke, engrot, &
                      strcom, vircom, thermo, cshell, cons, pmf, stat, &
                      rigid, domain, tmr, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, equations of motion in molecular dynamics
    ! - velocity verlet (symplectic)
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
    Real(Kind=wp),            Intent(InOut) :: tstep, strkin(1:9), strknf(1:9), strknt(1:9), &
                                               engke, engrot, strcom(1:9), vircom
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
    Integer                    :: fail(1:14), i, i1, i2, irgd, j, jrgd, krgd, lrgd, matms, rgdtyp
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: fmx, fmy, fmz, hstep, mxdr, opx, opy, opz, p0, p1, p2, p3, qt0, &
                                  qt1, qt2, qt3, rot(1:9), rstep, tmp, tqx, tqy, tqz, trx, try, &
                                  trz, vpx, vpy, vpz, x(1:1), y(1:1), z(1:1)
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), ggx(:), ggy(:), ggz(:), oxt(:), oyt(:), &
                                  ozt(:), q0t(:), q1t(:), q2t(:), q3t(:), rgdoxt(:), rgdoyt(:), &
                                  rgdozt(:), rgdvxt(:), rgdvyt(:), rgdvzt(:), rgdxxt(:), &
                                  rgdyyt(:), rgdzzt(:), vxt(:), vyt(:), vzt(:), xxt(:), yyt(:), &
                                  zzt(:)

    fail = 0
    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      Allocate (lstitr(1:config%natms), Stat=fail(1))
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
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'nve_1 allocation failure'
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
    If (comm%mxnode == 1) Then
      matms = config%natms
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

    ! update velocity and position of FPs

    Do j = 1, config%nfree
      i = config%lstfre(j)

      If (config%weight(i) > 1.0e-6_wp) Then
        tmp = hstep / config%weight(i)
        config%vxx(i) = vxt(i) + tmp * fxt(i)
        config%vyy(i) = vyt(i) + tmp * fyt(i)
        config%vzz(i) = vzt(i) + tmp * fzt(i)

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

            Call images(config%imcon, config%cell, 1, x, y, z)

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

          ! recover quaternion momenta at start of time step

          opx = rgdoxt(irgd) * rigid%rix(1, rgdtyp)
          opy = rgdoyt(irgd) * rigid%riy(1, rgdtyp)
          opz = rgdozt(irgd) * rigid%riz(1, rgdtyp)

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
          rigid%vxx(irgd) = rgdvxt(irgd) + tmp * fmx
          rigid%vyy(irgd) = rgdvyt(irgd) + tmp * fmy
          rigid%vzz(irgd) = rgdvzt(irgd) + tmp * fmz

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

          Do irgd = 1, rigid%n_types
            rigid%q0(irgd) = q0t(irgd)
            rigid%q1(irgd) = q1t(irgd)
            rigid%q2(irgd) = q2t(irgd)
            rigid%q3(irgd) = q3t(irgd)
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

      ! update kinetic energy and stress

      Call kinstresf(strknf, config, comm)
      Call kinstrest(rigid, strknt, comm)

      strkin = strknf + strknt
      engke = 0.5_wp * (strkin(1) + strkin(5) + strkin(9))

      ! update rotational energy

      engrot = getknr(rigid, comm)

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
      Write (message, '(a)') 'nve_1 deallocation failure'
      Call error(0, message)
    End If

  End Subroutine nve_1_vv
End Module nve
