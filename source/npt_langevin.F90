Module npt_langevin
  Use comms,           Only: comms_type
  Use configuration,   Only: configuration_type
  Use constants,       Only: boltz,&
                             pi
  Use constraints,     Only: apply_rattle,&
                             apply_shake,&
                             constraints_tags,&
                             constraints_type
  Use core_shell,      Only: core_shell_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info
  Use kinds,           Only: li,&
                             wp
  Use kinetics,        Only: getkin,&
                             getknf,&
                             getknr,&
                             getknt,&
                             getvom,&
                             kinstresf,&
                             kinstress,&
                             kinstrest
  Use langevin,        Only: langevin_forces
  Use npt_nose_hoover, Only: npt_h0_scl,&
                             npt_h1_scl
  Use numerics,        Only: box_mueller_saru1,&
                             images,&
                             seed_type
  Use pmf,             Only: pmf_tags,&
                             pmf_type
  Use rigid_bodies,    Only: getrotmat,&
                             no_squish,&
                             rigid_bodies_stress,&
                             rigid_bodies_type
  Use shared_units,    Only: update_shared_units
  Use site,            Only: site_type
  Use statistics,      Only: stats_type
  Use thermostat,      Only: VV_FIRST_STAGE,&
                             adjust_timestep,&
                             thermostat_type
  Use timer,           Only: timer_type
  Use vdw,             Only: vdw_type

  Implicit None

  Private

  Public :: npt_l0_vv, npt_l1_vv

Contains

  Subroutine npt_l0_vv(stage, lvar, mndis, mxdis, mxstp, tstep, &
                       nstep, &
                       degfre, virtot, &
                       consv, &
                       strkin, engke, &
                       cshell, cons, pmf, stat, thermo, sites, &
                       vdws, domain, tmr, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian equations of motion in
    ! molecular dynamics - velocity verlet with Langevin thermostat and
    ! barostat (isotropic pressure control) (symplectic)
    !
    ! isotropic config%cell fluctuations
    !
    ! reference: D. Quigley and M.I.J. Probert
    !            J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
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
    Integer,                  Intent(In   ) :: nstep
    Integer(Kind=li),         Intent(In   ) :: degfre
    Real(Kind=wp),            Intent(In   ) :: virtot
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
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail(1:9), i, iter
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: chip0, engke0, hstep, qstep, rstep, scale, str(1:9), tmp, vir, &
                                  vir1, vom(1:3), vzero
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), oxt(:), oyt(:), ozt(:), vxt(:), vyt(:), &
                                  vzt(:), xxt(:), yyt(:), zzt(:)

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
      Write (message, '(a)') 'npt_l  Use statistics, Only : stats_type0 allocation failure'
      Call error(0, message)
    End If

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    qstep = 0.5_wp * hstep
    rstep = 1.0_wp / tstep

    If (thermo%newjob_0) Then
      thermo%newjob_0 = .false.

      ! store initial values of volume, long range corrections and density

      thermo%cell0 = config%cell
      thermo%volm0 = config%volm
      thermo%elrc0 = vdws%elrc
      thermo%virlrc0 = vdws%vlrc

      Allocate (thermo%dens0(1:sites%mxatyp), Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'thermo%dens0 allocation failure'
        Call error(0, message)
      End If
      Do i = 1, sites%ntype_atom
        thermo%dens0(i) = sites%dens(i)
      End Do

      ! inertia parameter for barostat

      thermo%temp_lang = 2.0_wp * thermo%sigma / (boltz * Real(degfre, wp))
      thermo%pmass = (2.0_wp * thermo%sigma + 3.0_wp * boltz * thermo%temp_lang) / (2.0_wp * pi * thermo%tai)**2

      ! set number of constraint+pmf shake iterations and general iteration cycles

      thermo%mxiter = 1
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        thermo%mxkit = 1
        thermo%mxiter = thermo%mxiter + 3
      End If
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake

      ! Langevin forces for particles are now generated in w_calculate_forces
      ! Generate Langevin pseudo-tensor force for barostat piston

      thermo%fpl = 0.0_wp
      Call box_mueller_saru1(seed, Int(degfre / 3_li), nstep - 1, tmp)
      tmp = tmp * Sqrt(2.0_wp * thermo%tai * boltz * thermo%temp_lang * thermo%pmass * rstep) / 3.0_wp
      thermo%fpl(1) = tmp
      thermo%fpl(5) = tmp
      thermo%fpl(9) = tmp
    End If

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0) Call constraints_tags(lstitr, cons, config, comm)

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0) Call pmf_tags(lstitr, pmf, config, comm)
    End If

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

      vzero = config%volm
      chip0 = thermo%chi_p
      engke0 = engke

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

        ! integrate and apply Langevin thermostat - 1/4 step

        scale = Exp(-qstep * thermo%chi)
        Do i = 1, config%natms
          config%vxx(i) = scale * config%vxx(i)
          config%vyy(i) = scale * config%vyy(i)
          config%vzz(i) = scale * config%vzz(i)
        End Do
        engke = engke * scale**2

        ! constraint+pmf virial and stress

        vir = stat%vircon + stat%virpmf
        str = stat%strcon + stat%strpmf

        ! integrate and apply npt_h0_scl barostat - 1/2 step
        ! augment vir to include the random force on the barostat

        vir1 = vir - 3.0_wp * thermo%fpl(1)
        Call npt_h0_scl(1, hstep, degfre, thermo%tai, vir1, virtot, &
                        engke, config, thermo)

        ! integrate and apply Langevin thermostat - 1/4 step

        scale = Exp(-qstep * thermo%chi)
        Do i = 1, config%natms
          config%vxx(i) = scale * config%vxx(i)
          config%vyy(i) = scale * config%vyy(i)
          config%vzz(i) = scale * config%vzz(i)
        End Do
        engke = engke * scale**2

        ! update velocities

        Do i = 1, config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
            tmp = hstep / config%weight(i)
            config%vxx(i) = config%vxx(i) + tmp * (config%parts(i)%fxx + thermo%fxl(i))
            config%vyy(i) = config%vyy(i) + tmp * (config%parts(i)%fyy + thermo%fyl(i))
            config%vzz(i) = config%vzz(i) + tmp * (config%parts(i)%fzz + thermo%fzl(i))
          End If
        End Do

        ! update volume

        config%volm = config%volm * Exp(3.0_wp * tstep * thermo%chi_p)

        ! scale config%cell vectors - isotropic

        scale = (config%volm / thermo%volm0)**(1.0_wp / 3.0_wp)
        config%cell = thermo%cell0 * scale

        ! update positions

        scale = Exp(tstep * thermo%chi_p)
        Do i = 1, config%natms
          If (config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = scale * xxt(i) + tstep * config%vxx(i)
            config%parts(i)%yyy = scale * yyt(i) + tstep * config%vyy(i)
            config%parts(i)%zzz = scale * zzt(i) + tstep * config%vzz(i)
          End If
        End Do

        ! SHAKE procedures

        If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_shake(tstep, oxt, oyt, ozt, &
                           lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
        End If

        ! restore original integration parameters as well as
        ! velocities if iter < thermo%mxiter
        ! in the next iteration stat%vircon and stat%virpmf are freshly new

        If (iter < thermo%mxiter) Then
          config%volm = vzero
          thermo%chi_p = chip0
          engke = engke0

          Do i = 1, config%natms
            config%vxx(i) = vxt(i)
            config%vyy(i) = vyt(i)
            config%vzz(i) = vzt(i)
          End Do
        End If
      End Do

      ! check timestep for variable timestep

      If (lvar) Then
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, tmp, comm)) Then
          Call info(message, .true.)

          ! scale Langevin random forces

          Do i = 1, config%natms
            thermo%fxl(i) = thermo%fxl(i) * tmp
            thermo%fyl(i) = thermo%fyl(i) * tmp
            thermo%fzl(i) = thermo%fzl(i) * tmp
          End Do
          thermo%fpl = thermo%fpl * tmp

          ! restore initial conditions

          config%volm = vzero
          thermo%chi_p = chip0
          engke = engke0

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

      ! second stage of velocity verlet algorithm

    Else

      ! Generate Langevin forces for particles and
      ! Langevin pseudo-tensor force for barostat piston

      Call langevin_forces(nstep, thermo%temp_lang, tstep, thermo%chi, thermo%fxl, thermo%fyl, thermo%fzl, cshell, config, seed)

      thermo%fpl = 0.0_wp
      Call box_mueller_saru1(seed, Int(degfre / 3_li), nstep, tmp)
      tmp = tmp * Sqrt(2.0_wp * thermo%tai * boltz * thermo%temp_lang * thermo%pmass * rstep) / 3.0_wp
      thermo%fpl(1) = tmp
      thermo%fpl(5) = tmp
      thermo%fpl(9) = tmp

      ! update velocity

      Do i = 1, config%natms
        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * (config%parts(i)%fxx + thermo%fxl(i))
          config%vyy(i) = config%vyy(i) + tmp * (config%parts(i)%fyy + thermo%fyl(i))
          config%vzz(i) = config%vzz(i) + tmp * (config%parts(i)%fzz + thermo%fzl(i))
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! update kinetic energy

      engke = getkin(config, config%vxx, config%vyy, config%vzz, comm)

      ! integrate and apply Langevin thermostat - 1/4 step

      scale = Exp(-qstep * thermo%chi)
      Do i = 1, config%natms
        config%vxx(i) = scale * config%vxx(i)
        config%vyy(i) = scale * config%vyy(i)
        config%vzz(i) = scale * config%vzz(i)
      End Do
      engke = engke * scale**2

      ! constraint+pmf virial and stress

      vir = stat%vircon + stat%virpmf
      str = stat%strcon + stat%strpmf

      ! integrate and apply npt_h0_scl barostat - 1/2 step
      ! augment vir to include the random force on the barostat

      vir1 = vir - 3.0_wp * thermo%fpl(1)
      Call npt_h0_scl(1, hstep, degfre, thermo%tai, vir1, virtot, &
                      engke, config, thermo)

      ! integrate and apply Langevin thermostat - 1/4 step

      scale = Exp(-qstep * thermo%chi)
      Do i = 1, config%natms
        config%vxx(i) = scale * config%vxx(i)
        config%vyy(i) = scale * config%vyy(i)
        config%vzz(i) = scale * config%vzz(i)
      End Do
      engke = engke * scale**2

      ! conserved quantity less kinetic and potential energy terms

      consv = 0.5_wp * thermo%pmass * thermo%chi_p**2 + thermo%press * config%volm

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

    ! construct a 'mock' scaling tensor for xscale

    Do i = 2, 8
      thermo%eta(i) = 0.0_wp
    End Do
    thermo%eta(1) = thermo%chi_p
    thermo%eta(5) = thermo%chi_p
    thermo%eta(9) = thermo%chi_p

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
      Write (message, '(a)') 'npt_l0 deallocation failure, node'
      Call error(0, message)
    End If
  End Subroutine npt_l0_vv

  Subroutine npt_l1_vv(stage, lvar, mndis, mxdis, mxstp, tstep, &
                       nstep, &
                       degfre, degrot, virtot, &
                       consv, &
                       strkin, strknf, strknt, engke, engrot, &
                       strcom, vircom, &
                       cshell, cons, pmf, stat, thermo, sites, &
                       vdws, rigid, domain, tmr, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for integrating newtonian and rotational, singled
    ! RBs, equations of motion in molecular dynamics
    ! - velocity verlet with Langevin thermostat and
    ! barostat (isotropic pressure control) (symplectic)
    !
    ! isotropic config%cell fluctuations
    !
    ! reference: D. Quigley and M.I.J. Probert
    !            J. Chem. Phys., 2004, Vol. 120 (24), p. 11432
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    ! amended   - i.t.todorov november 2019 (RBs unsafe haloing)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: stage
    Logical,                  Intent(In   ) :: lvar
    Real(Kind=wp),            Intent(In   ) :: mndis, mxdis, mxstp
    Real(Kind=wp),            Intent(InOut) :: tstep
    Integer,                  Intent(In   ) :: nstep
    Integer(Kind=li),         Intent(In   ) :: degfre, degrot
    Real(Kind=wp),            Intent(In   ) :: virtot
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
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)         :: message
    Integer                    :: fail(1:14), i, i1, i2, irgd, iter, j, jrgd, krgd, lrgd, matms, &
                                  rgdtyp
    Logical, Allocatable       :: lstitr(:)
    Real(Kind=wp)              :: chip0, czero(1:9), engke0, engknf, engknt, engrot0, fmx, fmxl, &
                                  fmy, fmyl, fmz, fmzl, hstep, mxdr, opx, opy, opz, p0, p1, p2, &
                                  p3, qstep, qt0, qt0l, qt1, qt1l, qt2, qt2l, qt3, qt3l, rot(1:9), &
                                  rstep, scale, str(1:9), tmp, tqx, tqxl, tqy, tqyl, tqz, tqzl, &
                                  trx, trxl, try, tryl, trz, trzl, vir, vir1, vom(1:3), vpx, vpy, &
                                  vpz, vzero, x(1:1), y(1:1), z(1:1)
    Real(Kind=wp), Allocatable :: fxt(:), fyt(:), fzt(:), ggx(:), ggy(:), ggz(:), oxt(:), oyt(:), &
                                  ozt(:), q0t(:), q1t(:), q2t(:), q3t(:), rgdoxt(:), rgdoyt(:), &
                                  rgdozt(:), rgdvxt(:), rgdvyt(:), rgdvzt(:), rgdxxt(:), &
                                  rgdyyt(:), rgdzzt(:), vxt(:), vyt(:), vzt(:), xxt(:), yyt(:), &
                                  zzt(:)

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
      Write (message, '(a)') 'npt_l1 allocation failure, node'
      Call error(0, message)
    End If

    ! timestep derivatives

    hstep = 0.5_wp * tstep
    qstep = 0.5_wp * hstep
    rstep = 1.0_wp / tstep

    If (thermo%newjob_1) Then
      thermo%newjob_1 = .false.

      ! store initial values of volume, long range corrections and density

      thermo%cell0 = config%cell
      thermo%volm0 = config%volm
      thermo%elrc0 = vdws%elrc
      thermo%virlrc0 = vdws%vlrc

      Allocate (thermo%dens0(1:sites%mxatyp), Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message, '(a)') 'thermo%dens0 allocation failure'
        Call error(0, message)
      End If
      Do i = 1, sites%ntype_atom
        thermo%dens0(i) = sites%dens(i)
      End Do

      ! inertia parameter for barostat

      thermo%temp_lang = 2.0_wp * thermo%sigma / (boltz * Real(degfre, wp))
      thermo%pmass = (Real(degfre - degrot, wp) + 3.0_wp) * boltz * thermo%temp_lang / (2.0_wp * pi * thermo%tai)**2

      ! set number of constraint+pmf shake iterations and general iteration cycles

      thermo%mxiter = 1
      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        thermo%mxkit = 1
        thermo%mxiter = thermo%mxiter + 3
      End If
      If (cons%megcon > 0 .and. pmf%megpmf > 0) thermo%mxkit = cons%max_iter_shake

      ! thermo%unsafe positioning due to possibly locally shared RBs

      thermo%unsafe = Any(domain%map_unique > 0)

      ! Langevin forces for particles are now generated in w_calculate_forces
      ! Generate Langevin pseudo-tensor force for barostat piston

      thermo%fpl = 0.0_wp
      Call box_mueller_saru1(seed, Int(degfre / 3_li), nstep - 1, tmp)
      tmp = tmp * Sqrt(2.0_wp * thermo%tai * boltz * thermo%temp_lang * thermo%pmass * rstep) / 3.0_wp
      thermo%fpl(1) = tmp
      thermo%fpl(5) = tmp
      thermo%fpl(9) = tmp
    End If

    ! set matms

    matms = config%nlast
    If (comm%mxnode == 1) matms = config%natms

    If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
      lstitr(1:config%natms) = .false. ! initialise lstitr

      ! construct current bond vectors and listot array (shared
      ! constraint atoms) for iterative bond algorithms

      If (cons%megcon > 0) Call constraints_tags(lstitr, cons, config, comm)

      ! construct current PMF constraint vectors and shared description
      ! for iterative PMF constraint algorithms

      If (pmf%megpmf > 0) Call pmf_tags(lstitr, pmf, config, comm)
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

    ! first pass of velocity verlet algorithm

    If (stage == VV_FIRST_STAGE) Then

      ! Globalise Langevin random forces for shared RBs

      If (rigid%share) Then
        Call update_shared_units(config, rigid%list_shared, &
                                 rigid%map_shared, thermo%fxl, thermo%fyl, thermo%fzl, domain, comm)
      End If

      ! Get strcom & vircom when starting afresh now done in w_calculate_forces
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

      czero = config%cell
      vzero = config%volm
      chip0 = thermo%chi_p
      engke0 = engke
      engrot0 = engrot

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

        ! integrate and apply Langevin thermostat - 1/4 step

        scale = Exp(-qstep * thermo%chi)
        Do j = 1, config%nfree
          i = config%lstfre(j)

          config%vxx(i) = scale * config%vxx(i)
          config%vyy(i) = scale * config%vyy(i)
          config%vzz(i) = scale * config%vzz(i)
        End Do

        Do irgd = 1, rigid%n_types
          rigid%vxx(irgd) = scale * rigid%vxx(irgd)
          rigid%vyy(irgd) = scale * rigid%vyy(irgd)
          rigid%vzz(irgd) = scale * rigid%vzz(irgd)

          rigid%oxx(irgd) = scale * rigid%oxx(irgd)
          rigid%oyy(irgd) = scale * rigid%oyy(irgd)
          rigid%ozz(irgd) = scale * rigid%ozz(irgd)
        End Do
        engke = engke * scale**2
        engrot = engrot * scale**2

        ! constraint+pmf virial and stress

        vir = stat%vircon + stat%virpmf
        str = stat%strcon + stat%strpmf

        ! integrate and apply npt_h1_scl barostat - 1/2 step
        ! augment vir to include the random force on the barostat

        vir1 = vir - 3.0_wp * thermo%fpl(1)
        Call npt_h1_scl(1, hstep, degfre, degrot, thermo%tai, vir1, virtot, vircom, &
                        engke, rigid, config, thermo)

        ! integrate and apply Langevin thermostat - 1/4 step

        scale = Exp(-qstep * thermo%chi)
        Do j = 1, config%nfree
          i = config%lstfre(j)

          config%vxx(i) = scale * config%vxx(i)
          config%vyy(i) = scale * config%vyy(i)
          config%vzz(i) = scale * config%vzz(i)
        End Do

        Do irgd = 1, rigid%n_types
          rigid%vxx(irgd) = scale * rigid%vxx(irgd)
          rigid%vyy(irgd) = scale * rigid%vyy(irgd)
          rigid%vzz(irgd) = scale * rigid%vzz(irgd)

          rigid%oxx(irgd) = scale * rigid%oxx(irgd)
          rigid%oyy(irgd) = scale * rigid%oyy(irgd)
          rigid%ozz(irgd) = scale * rigid%ozz(irgd)
        End Do
        engke = engke * scale**2
        engrot = engrot * scale**2

        ! update velocity of FPs

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
            tmp = hstep / config%weight(i)
            config%vxx(i) = config%vxx(i) + tmp * (config%parts(i)%fxx + thermo%fxl(i))
            config%vyy(i) = config%vyy(i) + tmp * (config%parts(i)%fyy + thermo%fyl(i))
            config%vzz(i) = config%vzz(i) + tmp * (config%parts(i)%fzz + thermo%fzl(i))
          End If
        End Do

        ! update volume

        config%volm = config%volm * Exp(3.0_wp * tstep * thermo%chi_p)

        ! scale config%cell vectors - isotropic

        scale = (config%volm / thermo%volm0)**(1.0_wp / 3.0_wp)
        config%cell = thermo%cell0 * scale

        ! update position of FPs

        scale = Exp(tstep * thermo%chi_p)
        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = scale * xxt(i) + tstep * config%vxx(i)
            config%parts(i)%yyy = scale * yyt(i) + tstep * config%vyy(i)
            config%parts(i)%zzz = scale * zzt(i) + tstep * config%vzz(i)
          End If
        End Do

        ! SHAKE procedures

        If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
          Call apply_shake(tstep, oxt, oyt, ozt, &
                           lstitr, stat, pmf, cons, domain, tmr, config, thermo, comm)
        End If

        ! restore original integration parameters as well as
        ! velocities if iter < thermo%mxiter
        ! in the next iteration stat%vircon and stat%virpmf are freshly new

        If (iter < thermo%mxiter) Then
          config%volm = vzero
          thermo%chi_p = chip0
          engke = engke0
          engrot = engrot0

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

          fmxl = 0.0_wp; fmyl = 0.0_wp; fmzl = 0.0_wp
          tqxl = 0.0_wp; tqyl = 0.0_wp; tqzl = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force

            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + fxt(i)
              fmy = fmy + fyt(i)
              fmz = fmz + fzt(i)

              fmxl = fmxl + thermo%fxl(i)
              fmyl = fmyl + thermo%fyl(i)
              fmzl = fmzl + thermo%fzl(i)
            End If

            tqx = tqx + ggy(krgd) * fzt(i) - ggz(krgd) * fyt(i)
            tqy = tqy + ggz(krgd) * fxt(i) - ggx(krgd) * fzt(i)
            tqz = tqz + ggx(krgd) * fyt(i) - ggy(krgd) * fxt(i)

            tqxl = tqxl + ggy(krgd) * thermo%fzl(i) - ggz(krgd) * thermo%fyl(i)
            tqyl = tqyl + ggz(krgd) * thermo%fxl(i) - ggx(krgd) * thermo%fzl(i)
            tqzl = tqzl + ggx(krgd) * thermo%fyl(i) - ggy(krgd) * thermo%fxl(i)
          End Do

          ! If the RB has 2+ frozen particles (ill=1) the net torque
          ! must align along the axis of rotation

          If (rigid%frozen(0, rgdtyp) > 1) Then
            i1 = rigid%index_local(rigid%index_global(1, rgdtyp), irgd)
            i2 = rigid%index_local(rigid%index_global(2, rgdtyp), irgd)

            x(1) = xxt(i1) - xxt(i2)
            y(1) = yyt(i1) - yyt(i2)
            z(1) = zzt(i1) - zzt(i2)

            Call images(config%imcon, czero, 1, x, y, z)

            mxdr = 1.0_wp / (x(1)**2 + y(1)**2 + z(1)**2)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) * mxdr
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp

            tmp = (x(1) * tqxl + y(1) * tqyl + z(1) * tqzl) * mxdr
            tqxl = x(1) * tmp
            tqyl = y(1) * tmp
            tqzl = z(1) * tmp
          End If

          ! current rotation matrix

          Call getrotmat(q0t(irgd), q1t(irgd), q2t(irgd), q3t(irgd), rot)

          ! calculate torque in principal frame

          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          trxl = tqxl * rot(1) + tqyl * rot(4) + tqzl * rot(7)
          tryl = tqxl * rot(2) + tqyl * rot(5) + tqzl * rot(8)
          trzl = tqxl * rot(3) + tqyl * rot(6) + tqzl * rot(9)

          ! calculate quaternion torques

          qt0 = 2.0_wp * (-q1t(irgd) * trx - q2t(irgd) * try - q3t(irgd) * trz)
          qt1 = 2.0_wp * (q0t(irgd) * trx - q3t(irgd) * try + q2t(irgd) * trz)
          qt2 = 2.0_wp * (q3t(irgd) * trx + q0t(irgd) * try - q1t(irgd) * trz)
          qt3 = 2.0_wp * (-q2t(irgd) * trx + q1t(irgd) * try + q0t(irgd) * trz)

          qt0l = 2.0_wp * (-q1t(irgd) * trxl - q2t(irgd) * tryl - q3t(irgd) * trzl)
          qt1l = 2.0_wp * (q0t(irgd) * trxl - q3t(irgd) * tryl + q2t(irgd) * trzl)
          qt2l = 2.0_wp * (q3t(irgd) * trxl + q0t(irgd) * tryl - q1t(irgd) * trzl)
          qt3l = 2.0_wp * (-q2t(irgd) * trxl + q1t(irgd) * tryl + q0t(irgd) * trzl)

          ! recover quaternion momenta

          opx = rigid%oxx(irgd) * rigid%rix(1, rgdtyp)
          opy = rigid%oyy(irgd) * rigid%riy(1, rgdtyp)
          opz = rigid%ozz(irgd) * rigid%riz(1, rgdtyp)

          p0 = 2.0_wp * (-q1t(irgd) * opx - q2t(irgd) * opy - q3t(irgd) * opz)
          p1 = 2.0_wp * (q0t(irgd) * opx - q3t(irgd) * opy + q2t(irgd) * opz)
          p2 = 2.0_wp * (q3t(irgd) * opx + q0t(irgd) * opy - q1t(irgd) * opz)
          p3 = 2.0_wp * (-q2t(irgd) * opx + q1t(irgd) * opy + q0t(irgd) * opz)

          ! update quaternion momenta to half step

          p0 = p0 + hstep * (qt0 + qt0l)
          p1 = p1 + hstep * (qt1 + qt1l)
          p2 = p2 + hstep * (qt2 + qt2l)
          p3 = p3 + hstep * (qt3 + qt3l)

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
          rigid%vxx(irgd) = rigid%vxx(irgd) + tmp * (fmx + fmxl)
          rigid%vyy(irgd) = rigid%vyy(irgd) + tmp * (fmy + fmyl)
          rigid%vzz(irgd) = rigid%vzz(irgd) + tmp * (fmz + fmzl)

          ! update RB COM to full step

          rigid%xxx(irgd) = scale * rgdxxt(irgd) + tstep * rigid%vxx(irgd)
          rigid%yyy(irgd) = scale * rgdyyt(irgd) + tstep * rigid%vyy(irgd)
          rigid%zzz(irgd) = scale * rgdzzt(irgd) + tstep * rigid%vzz(irgd)

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
                  config%vxx(i) = scale * xxt(i)
                  config%vyy(i) = scale * yyt(i)
                  config%vzz(i) = scale * zzt(i)

                  x(1) = config%parts(i)%xxx - config%vxx(i)
                  y(1) = config%parts(i)%yyy - config%vyy(i)
                  z(1) = config%parts(i)%zzz - config%vzz(i)
                  Call images(config%imcon, config%cell, 1, x, y, z)
                  config%parts(i)%xxx = x(1) + config%vxx(i)
                  config%parts(i)%yyy = y(1) + config%vyy(i)
                  config%parts(i)%zzz = z(1) + config%vzz(i)
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

          rigid%xxx(irgd) = scale * rgdxxt(irgd)
          rigid%yyy(irgd) = scale * rgdyyt(irgd)
          rigid%zzz(irgd) = scale * rgdzzt(irgd)

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
        If (adjust_timestep(tstep, hstep, rstep, mndis, mxdis, mxstp, config%natms, config%parts, &
                            xxt, yyt, zzt, cshell%legshl, message, tmp, comm)) Then
          Call info(message, .true.)

          ! scale Langevin random forces

          Do i = 1, matms
            thermo%fxl(i) = thermo%fxl(i) * tmp
            thermo%fyl(i) = thermo%fyl(i) * tmp
            thermo%fzl(i) = thermo%fzl(i) * tmp
          End Do
          thermo%fpl = thermo%fpl * tmp

          ! restore initial conditions

          config%volm = vzero
          thermo%chi_p = chip0
          engke = engke0
          engrot = engrot0

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

      ! second stage of velocity verlet algorithm

    Else

      ! Generate Langevin forces for particles and
      ! Langevin pseudo-tensor force for barostat piston

      Call langevin_forces(nstep, thermo%temp_lang, tstep, thermo%chi, thermo%fxl, thermo%fyl, thermo%fzl, &
                           cshell, config, seed)
      If (rigid%share) Then
        Call update_shared_units(config, rigid%list_shared, &
                                 rigid%map_shared, thermo%fxl, thermo%fyl, thermo%fzl, domain, comm)
      Endif

      thermo%fpl = 0.0_wp
      Call box_mueller_saru1(seed, Int(degfre / 3_li), nstep, tmp)
      tmp = tmp * Sqrt(2.0_wp * thermo%tai * boltz * thermo%temp_lang * thermo%pmass * rstep) / 3.0_wp
      thermo%fpl(1) = tmp
      thermo%fpl(5) = tmp
      thermo%fpl(9) = tmp

      ! update velocity of FPs

      Do j = 1, config%nfree
        i = config%lstfre(j)

        If (config%weight(i) > 1.0e-6_wp) Then
          tmp = hstep / config%weight(i)
          config%vxx(i) = config%vxx(i) + tmp * (config%parts(i)%fxx + thermo%fxl(i))
          config%vyy(i) = config%vyy(i) + tmp * (config%parts(i)%fyy + thermo%fyl(i))
          config%vzz(i) = config%vzz(i) + tmp * (config%parts(i)%fzz + thermo%fzl(i))
        End If
      End Do

      ! RATTLE procedures
      ! apply velocity corrections to bond and PMF constraints

      If (cons%megcon > 0 .or. pmf%megpmf > 0) Then
        Call apply_rattle(tstep, thermo%kit, pmf, cons, stat, domain, tmr, config, comm)
      End If

      ! Get RB COM stress and virial

      Call rigid_bodies_stress(strcom, ggx, ggy, ggz, config, rigid, comm, thermo%fxl, thermo%fyl, thermo%fzl)
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

          fmxl = 0.0_wp; fmyl = 0.0_wp; fmzl = 0.0_wp
          tqxl = 0.0_wp; tqyl = 0.0_wp; tqzl = 0.0_wp
          Do jrgd = 1, lrgd
            krgd = krgd + 1

            i = rigid%index_local(jrgd, irgd) ! local index of particle/site

            ! If the RB has a frozen particle then no net force

            If (rigid%frozen(0, rgdtyp) == 0) Then
              fmx = fmx + config%parts(i)%fxx
              fmy = fmy + config%parts(i)%fyy
              fmz = fmz + config%parts(i)%fzz

              fmxl = fmxl + thermo%fxl(i)
              fmyl = fmyl + thermo%fyl(i)
              fmzl = fmzl + thermo%fzl(i)
            End If

            tqx = tqx + ggy(krgd) * config%parts(i)%fzz - ggz(krgd) * config%parts(i)%fyy
            tqy = tqy + ggz(krgd) * config%parts(i)%fxx - ggx(krgd) * config%parts(i)%fzz
            tqz = tqz + ggx(krgd) * config%parts(i)%fyy - ggy(krgd) * config%parts(i)%fxx

            tqxl = tqxl + ggy(krgd) * thermo%fzl(i) - ggz(krgd) * thermo%fyl(i)
            tqyl = tqyl + ggz(krgd) * thermo%fxl(i) - ggx(krgd) * thermo%fzl(i)
            tqzl = tqzl + ggx(krgd) * thermo%fyl(i) - ggy(krgd) * thermo%fxl(i)
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

            mxdr = 1.0_wp / (x(1)**2 + y(1)**2 + z(1)**2)

            tmp = (x(1) * tqx + y(1) * tqy + z(1) * tqz) * mxdr
            tqx = x(1) * tmp
            tqy = y(1) * tmp
            tqz = z(1) * tmp

            tmp = (x(1) * tqxl + y(1) * tqyl + z(1) * tqzl) * mxdr
            tqxl = x(1) * tmp
            tqyl = y(1) * tmp
            tqzl = z(1) * tmp
          End If

          ! current rotation matrix

          Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

          ! calculate torque in principal frame

          trx = tqx * rot(1) + tqy * rot(4) + tqz * rot(7)
          try = tqx * rot(2) + tqy * rot(5) + tqz * rot(8)
          trz = tqx * rot(3) + tqy * rot(6) + tqz * rot(9)

          trxl = tqxl * rot(1) + tqyl * rot(4) + tqzl * rot(7)
          tryl = tqxl * rot(2) + tqyl * rot(5) + tqzl * rot(8)
          trzl = tqxl * rot(3) + tqyl * rot(6) + tqzl * rot(9)

          ! calculate quaternion torques

          qt0 = 2.0_wp * (-rigid%q1(irgd) * trx - rigid%q2(irgd) * try - rigid%q3(irgd) * trz)
          qt1 = 2.0_wp * (rigid%q0(irgd) * trx - rigid%q3(irgd) * try + rigid%q2(irgd) * trz)
          qt2 = 2.0_wp * (rigid%q3(irgd) * trx + rigid%q0(irgd) * try - rigid%q1(irgd) * trz)
          qt3 = 2.0_wp * (-rigid%q2(irgd) * trx + rigid%q1(irgd) * try + rigid%q0(irgd) * trz)

          qt0l = 2.0_wp * (-rigid%q1(irgd) * trxl - rigid%q2(irgd) * tryl - rigid%q3(irgd) * trzl)
          qt1l = 2.0_wp * (rigid%q0(irgd) * trxl - rigid%q3(irgd) * tryl + rigid%q2(irgd) * trzl)
          qt2l = 2.0_wp * (rigid%q3(irgd) * trxl + rigid%q0(irgd) * tryl - rigid%q1(irgd) * trzl)
          qt3l = 2.0_wp * (-rigid%q2(irgd) * trxl + rigid%q1(irgd) * tryl + rigid%q0(irgd) * trzl)

          ! recover quaternion momenta at half time step

          opx = rigid%oxx(irgd) * rigid%rix(1, rgdtyp)
          opy = rigid%oyy(irgd) * rigid%riy(1, rgdtyp)
          opz = rigid%ozz(irgd) * rigid%riz(1, rgdtyp)

          p0 = 2.0_wp * (-rigid%q1(irgd) * opx - rigid%q2(irgd) * opy - rigid%q3(irgd) * opz)
          p1 = 2.0_wp * (rigid%q0(irgd) * opx - rigid%q3(irgd) * opy + rigid%q2(irgd) * opz)
          p2 = 2.0_wp * (rigid%q3(irgd) * opx + rigid%q0(irgd) * opy - rigid%q1(irgd) * opz)
          p3 = 2.0_wp * (-rigid%q2(irgd) * opx + rigid%q1(irgd) * opy + rigid%q0(irgd) * opz)

          ! update quaternion momenta to full step

          p0 = p0 + hstep * (qt0 + qt0l)
          p1 = p1 + hstep * (qt1 + qt1l)
          p2 = p2 + hstep * (qt2 + qt2l)
          p3 = p3 + hstep * (qt3 + qt3l)

          ! update RB angular & COM velocities to full step

          opx = 0.5_wp * (-rigid%q1(irgd) * p0 + rigid%q0(irgd) * p1 + rigid%q3(irgd) * p2 - rigid%q2(irgd) * p3)
          opy = 0.5_wp * (-rigid%q2(irgd) * p0 - rigid%q3(irgd) * p1 + rigid%q0(irgd) * p2 + rigid%q1(irgd) * p3)
          opz = 0.5_wp * (-rigid%q3(irgd) * p0 + rigid%q2(irgd) * p1 - rigid%q1(irgd) * p2 + rigid%q0(irgd) * p3)

          rigid%oxx(irgd) = opx * rigid%rix(2, rgdtyp)
          rigid%oyy(irgd) = opy * rigid%riy(2, rgdtyp)
          rigid%ozz(irgd) = opz * rigid%riz(2, rgdtyp)

          tmp = hstep / rigid%weight(0, rgdtyp)
          rigid%vxx(irgd) = rigid%vxx(irgd) + tmp * (fmx + fmxl)
          rigid%vyy(irgd) = rigid%vyy(irgd) + tmp * (fmy + fmyl)
          rigid%vzz(irgd) = rigid%vzz(irgd) + tmp * (fmz + fmzl)

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

      ! update kinetic energy

      engknf = getknf(config%vxx, config%vyy, config%vzz, config, comm)
      engknt = getknt(rigid, comm)

      engke = engknf + engknt

      ! update rotational energy

      engrot = getknr(rigid, comm)

      ! integrate and apply Langevin thermostat - 1/4 step

      scale = Exp(-qstep * thermo%chi)
      Do j = 1, config%nfree
        i = config%lstfre(j)

        config%vxx(i) = scale * config%vxx(i)
        config%vyy(i) = scale * config%vyy(i)
        config%vzz(i) = scale * config%vzz(i)
      End Do

      Do irgd = 1, rigid%n_types
        rigid%vxx(irgd) = scale * rigid%vxx(irgd)
        rigid%vyy(irgd) = scale * rigid%vyy(irgd)
        rigid%vzz(irgd) = scale * rigid%vzz(irgd)

        rigid%oxx(irgd) = scale * rigid%oxx(irgd)
        rigid%oyy(irgd) = scale * rigid%oyy(irgd)
        rigid%ozz(irgd) = scale * rigid%ozz(irgd)
      End Do
      engke = engke * scale**2
      engrot = engrot * scale**2

      ! constraint+pmf virial and stress

      vir = stat%vircon + stat%virpmf
      str = stat%strcon + stat%strpmf

      ! integrate and apply npt_h1_scl barostat - 1/2 step
      ! augment vir to include the random force on the barostat

      vir1 = vir - 3.0_wp * thermo%fpl(1)
      Call npt_h1_scl(1, hstep, degfre, degrot, thermo%tai, vir1, virtot, vircom, &
                      engke, rigid, config, thermo)

      ! integrate and apply Langevin thermostat - 1/4 step

      scale = Exp(-qstep * thermo%chi)
      Do j = 1, config%nfree
        i = config%lstfre(j)

        config%vxx(i) = scale * config%vxx(i)
        config%vyy(i) = scale * config%vyy(i)
        config%vzz(i) = scale * config%vzz(i)
      End Do

      Do irgd = 1, rigid%n_types
        rigid%vxx(irgd) = scale * rigid%vxx(irgd)
        rigid%vyy(irgd) = scale * rigid%vyy(irgd)
        rigid%vzz(irgd) = scale * rigid%vzz(irgd)

        rigid%oxx(irgd) = scale * rigid%oxx(irgd)
        rigid%oyy(irgd) = scale * rigid%oyy(irgd)
        rigid%ozz(irgd) = scale * rigid%ozz(irgd)
      End Do
      engke = engke * scale**2
      engrot = engrot * scale**2

      ! conserved quantity less kinetic and potential energy terms

      consv = 0.5_wp * thermo%pmass * thermo%chi_p**2 + thermo%press * config%volm

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

    ! construct a 'mock' scaling tensor for xscale

    Do i = 2, 8
      thermo%eta(i) = 0.0_wp
    End Do
    thermo%eta(1) = thermo%chi_p
    thermo%eta(5) = thermo%chi_p
    thermo%eta(9) = thermo%chi_p

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
      Write (message, '(a)') 'npt_l1 deallocation failure'
      Call error(0, message)
    End If
  End Subroutine npt_l1_vv
End Module npt_langevin
