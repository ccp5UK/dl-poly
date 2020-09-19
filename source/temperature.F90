Module temperature
  Use comms,           Only: comms_type,&
                             gsum
  Use configuration,   Only: configuration_type,&
                             getcom,&
                             IMCON_NOPBC,&
                             IMCON_CUBIC,&
                             IMCON_ORTHORHOMBIC,&
                             IMCON_PARALLELOPIPED,&
                             IMCON_SLAB,&
                             IMCON_TRUNC_OCTO,&
                             IMCON_RHOMBIC_DODEC,&
                             IMCON_HEXAGONAL
  Use constants,       Only: boltz,&
                             zero_plus
  Use constraints,     Only: constraints_quench,&
                             constraints_type
  Use core_shell,      Only: SHELL_ADIABATIC,&
                             core_shell_quench,&
                             core_shell_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use flow_control,    Only: RESTART_KEY_CLEAN,&
                             RESTART_KEY_SCALE
  Use kinds,           Only: li,&
                             wp
  Use kinetics,        Only: getkin,&
                             getknf,&
                             getknr,&
                             getknt,&
                             getvom
  Use minimise,        Only: minimise_type
  Use numerics,        Only: box_mueller_saru3,&
                             invert,&
                             local_index,&
                             seed_type,&
                             uni
  Use pmf,             Only: pmf_quench,&
                             pmf_type
  Use rigid_bodies,    Only: getrotmat,&
                             rigid_bodies_quench,&
                             rigid_bodies_type
  Use shared_units,    Only: update_shared_units,&
                             update_shared_units_int
  Use statistics,      Only: stats_type
  Use thermostat,      Only: DPD_NULL,&
                             thermostat_type

  Implicit None

  Private

  Public :: set_temperature, regauss_temperature, scale_temperature

Contains

  Subroutine set_temperature(keyres, nstep, nstrun, &
                             engrot, dof_site, cshell, stat, cons, pmf, thermo, minim, &
                             rigid, domain, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for setting the initial system temperature
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                     Intent(InOut) :: keyres
    Integer,                     Intent(In   ) :: nstep, nstrun
    Real(Kind=wp),               Intent(  Out) :: engrot
    Real(Kind=wp), Dimension(:), Intent(In   ) :: dof_site
    Type(core_shell_type),       Intent(InOut) :: cshell
    Type(stats_type),            Intent(InOut) :: stat
    Type(constraints_type),      Intent(InOut) :: cons
    Type(pmf_type),              Intent(InOut) :: pmf
    Type(thermostat_type),       Intent(InOut) :: thermo
    Type(minimise_type),         Intent(In   ) :: minim
    Type(rigid_bodies_type),     Intent(InOut) :: rigid
    Type(domains_type),          Intent(In   ) :: domain
    Type(configuration_type),    Intent(InOut) :: config
    Type(seed_type),             Intent(InOut) :: seed
    Type(comms_type),            Intent(InOut) :: comm

    Character(Len=256)   :: message, messages(10)
    Integer              :: fail(1:2), i, i1, i2, irgd, j, jrgd, k, lrgd, ntp, rgdtyp, stp
    Integer(Kind=li)     :: com, con, frz, meg, non
    Integer, Allocatable :: qn(:), qs(:, :), tpn(:), tps(:)
    Logical              :: no_min_0, safe
    Real(Kind=wp)        :: engf, engk, engr, engt, rot(1:9), tmp, vom(1:3), vpx, vpy, vpz, &
                            x(1:1), y(1:1), z(1:1)

! q. index arrays and tp. sum arrays

    ! initialise rotational and translational DoF if no RB are present
    ! or re-initialise if all are frozen (does no harm)

    If (rigid%total == 0) Then
      config%degtra = Int(0, li)
      config%degrot = Int(0, li)
    End If

    ! Degrees Of Freedom (DoF) - all free particles

    meg = Int(3, li) * Int(config%atmfre, li)

    ! lost to frozen free atoms

    frz = Int(3, li) * Int(config%atmfrz, li)

    ! 3 lost for fixing COM translation

    If (config%l_vom .and. thermo%key_dpd == DPD_NULL) Then
      com = Int(3, li)
    Else
      com = Int(0, li)
    End If

    ! 3 lost for fixing angular momentum about origin
    ! (non-periodic systems only)

    non = Int(0, li)
    If (config%imcon == IMCON_NOPBC) non = Int(3, li)

    ! lost to shells

    config%degshl = Int(3, li) * Int(cshell%megshl, li)

    ! lost to constrained atoms and PMF constraints

    con = Int(cons%megcon) + Int(pmf%megpmf)

    ! TOTAL DoF

    config%degfre = meg - com - non - frz - config%degshl - con + config%degrot + config%degtra

    ! Report DoF

    Write (messages(1), '(a)') 'degrees of freedom break-down neigh%list:'
    Write (messages(2), '(2x,a,i12)') 'free particles        ', meg
    Write (messages(3), '(2x,a,i12)') 'centre of mass        ', -com
    Write (messages(4), '(2x,a,i12)') 'non-periodicity       ', -non
    Write (messages(5), '(2x,a,i12)') 'frozen free particles ', -frz
    Write (messages(6), '(2x,a,i12)') 'shell-pseudo          ', -config%degshl
    Write (messages(7), '(2x,a,i12)') 'constrained           ', -con
    Write (messages(8), '(2x,a,i12)') 'RB translational      ', config%degtra
    Write (messages(9), '(2x,a,i12)') 'RB rotational         ', config%degrot
    Write (messages(10), '(2x,a,i12)') 'total (real)          ', config%degfre
    Call info(messages, 10, .true.)

    ! Check DoF distribution

    tmp = 0.0_wp
    Do i = 1, config%natms
      If (dof_site(config%lsite(i)) > zero_plus) & ! Omit shells' negative DoFs
        tmp = tmp + dof_site(config%lsite(i))
    End Do
    Call gsum(comm, tmp)
    If (Nint(tmp, li) - non - com /= config%degfre) Call error(360)

    ! warn for much restrain on the system [ config%degfre <= (com+non+frz+con) ]
    ! and catch an over-restrained system

    If (config%degfre <= Int(com + non + frz + con, li)) &
      Call warning(210, Real(config%degfre, wp), Real((com + non + frz + con), wp), 0.0_wp)
    If (config%degfre < Int(1, li)) Call error(350)

    ! desired kinetic energy

    thermo%sigma = 0.5_wp * Real(config%degfre, wp) * boltz * thermo%temp

    ! avoid user defined 0K field to break up anything

    If (rigid%total > 0) Then
      engf = getknf(config%vxx, config%vyy, config%vzz, config, comm) / Real(Max(1_li, config%degfre), wp)
      engt = getknt(rigid, comm) / Real(Max(1_li, config%degtra), wp)

      engr = getknr(rigid, comm) / Real(Max(1_li, config%degrot), wp)
      engk = engf + engt
    Else
      engk = getkin(config, config%vxx, config%vyy, config%vzz, comm) / Real(Max(1_li, config%degfre), wp)
    End If
    If (thermo%sigma > 1.0e-6_wp .and. engk < 1.0e-6_wp .and. (keyres /= RESTART_KEY_CLEAN .and. nstrun /= 0)) Then
      Call warning('0K velocity field detected in CONFIG with a restart at non 0K temperature in CONTROL', .true.)
      Call info('*** clean start enforced ***', .true.)

      keyres = RESTART_KEY_CLEAN
    End If

    If (keyres == RESTART_KEY_CLEAN) Then

      Allocate (qn(1:config%mxatms), tpn(0:comm%mxnode - 1), Stat=fail(1))
      Allocate (qs(0:2, 1:cshell%mxshl), tps(0:comm%mxnode - 1), Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'set_temperature allocation failure'
        Call error(0, message)
      End If

      ! tpn(idnode) number of particles on this node (idnode)
      ! ntp - grand total of non-shell, non-frozen particles

      qn(1:config%natms) = 0 ! unqualified particle (non-massless, non-shells, non-frozen)
      qs(0:2, 1:cshell%ntshl) = 0 ! unqualified core-shell unit with a local shell

      j = 0
      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. cshell%legshl(0, i) >= 0) Then
          j = j + 1
          qn(i) = 1
        End If
      End Do
      tpn(comm%idnode) = j
      Do i = 0, comm%mxnode - 1
        If (i /= comm%idnode) tpn(i) = 0
      End Do
      Call gsum(comm, tpn)
      ntp = Sum(tpn)

      ! Save core-shell internal energy vectors due to possible hits
      ! tps(idnode) number of thermostatted core-shell units on this node (idnode)
      ! stp - grand total of core-shell units to thermostat

      j = 0
      If (cshell%keyshl == SHELL_ADIABATIC) Then ! just for the adiabatic shell model
        If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
          qn(config%natms + 1:config%nlast) = 0
          Call update_shared_units_int(config, cshell%lishp_shl, &
                                       cshell%lashp_shl, qn, domain, comm)
          Call update_shared_units(config, cshell%lishp_shl, &
                                   cshell%lashp_shl, config%vxx, config%vyy, config%vzz, domain, comm)
        End If

        If (cshell%ntshl > 0) Then
          Do k = 1, cshell%ntshl
            i1 = local_index(cshell%listshl(1, k), config%nlast, config%lsi, config%lsa)
            i2 = local_index(cshell%listshl(2, k), config%nlast, config%lsi, config%lsa)

            If (qn(i1) == 1 .and. i2 > 0 .and. i2 <= config%natms) Then
              j = j + 1

              qs(0, k) = 1
              qs(1, k) = i1
              qs(2, k) = i2
            End If
          End Do
        End If
      End If
      tps(comm%idnode) = j
      Do i = 0, comm%mxnode - 1
        If (i /= comm%idnode) tps(i) = 0
      End Do
      Call gsum(comm, tps)
      stp = Sum(tps)

      If (rigid%total > 0) Then

        k = 0
        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          ! For all good RBs

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd) ! local index of particle/site
              If (i <= config%natms) tpn(comm%idnode) = tpn(comm%idnode) - 1 ! Less free particles are hit
            End Do

            i1 = rigid%index_local(1, irgd) ! particle to bare the random RB COM momentum
            i2 = rigid%index_local(2, irgd) ! particle to bare the random RB angular momentum
            If (rigid%frozen(0, rgdtyp) == 0) Then
              If (i1 <= config%natms) k = k + 1
            End If
            If (i2 <= config%natms) k = k + 1
          End If
        End Do
        ! tpn(idnode) number of thermostatted free particles on this node (idnode)
        ! ntp - grand total of non-shell, non-frozen free particles to thermostat
        Do i = 0, comm%mxnode - 1
          If (i /= comm%idnode) tpn(i) = 0
        End Do
        Call gsum(comm, tpn)
        ntp = Sum(tpn)

        ! generate starting velocities

        Do i = 1, config%natms

          ! frozen and massless atoms are either motionless
          ! (with no actual DoF - relaxed shells, frozen sites)
          ! or despite their DoF, their motion is defined by
          ! the particles with masses in a RB.  The rest are
          ! to have gaussian distribution of their Ekin
          ! (adiabatic shells, constraints, PMFs and RBs are
          ! to be sorted out later by quenching)

          If (qn(i) == 1 .and. config%lfree(i) == 0) Then
            Call box_mueller_saru3(seed, config%ltg(i), 0, config%vxx(i), config%vyy(i), config%vzz(i))

            ! Get scaler to target variance/Sqrt(config%weight)

            tmp = 1.0_wp / Sqrt(config%weight(i))
            config%vxx(i) = config%vxx(i) * tmp
            config%vyy(i) = config%vyy(i) * tmp
            config%vzz(i) = config%vzz(i) * tmp
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          ! For all good RBs

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then
            i1 = rigid%index_local(1, irgd) ! particle to bare the random RB COM momentum
            i2 = rigid%index_local(2, irgd) ! particle to bare the random RB angular momentum

            If (rigid%frozen(0, rgdtyp) == 0 .and. i1 <= config%natms) Then
              Call box_mueller_saru3(seed, config%ltg(i1), 0, config%vxx(i1), config%vyy(i1), config%vzz(i1))

              ! Get scaler to target variance/Sqrt(config%weight)

              tmp = 1.0_wp / Sqrt(rigid%weight(0, rgdtyp))
              config%vxx(i1) = config%vxx(i1) * tmp
              config%vyy(i1) = config%vyy(i1) * tmp
              config%vzz(i1) = config%vzz(i1) * tmp
            End If

            If (i2 <= config%natms) Then
              Call box_mueller_saru3(seed, config%ltg(i2), 0, config%vxx(i2), config%vyy(i2), config%vzz(i2))

              ! Get scaler to target variance/Sqrt(config%weight) -
              ! 3 different reciprocal moments of inertia

              config%vxx(i2) = config%vxx(i2) * Sqrt(rigid%rix(2, rgdtyp))
              config%vyy(i2) = config%vyy(i2) * Sqrt(rigid%riy(2, rgdtyp))
              config%vzz(i2) = config%vzz(i2) * Sqrt(rigid%riz(2, rgdtyp))
            End If
          End If
        End Do

        ! Update shared RBs' velocities

        If (rigid%share) Then
          Call update_shared_units(config, rigid%list_shared, &
                                   rigid%map_shared, config%vxx, config%vyy, config%vzz, domain, comm)
        End If

        ! calculate new RBs' COM and angular velocities

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          ! For all good RBs

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then
            i1 = rigid%index_local(1, irgd) ! particle to bare the random RB COM momentum
            i2 = rigid%index_local(2, irgd) ! particle to bare the random RB angular momentum

            If (rigid%frozen(0, rgdtyp) == 0) Then
              rigid%vxx(irgd) = config%vxx(i1)
              rigid%vyy(irgd) = config%vyy(i1)
              rigid%vzz(irgd) = config%vzz(i1)
            End If

            rigid%oxx(irgd) = config%vxx(i2)
            rigid%oyy(irgd) = config%vyy(i2)
            rigid%ozz(irgd) = config%vzz(i2)
            If (i2 <= config%natms) Then
              If (config%lfrzn(i2) > 0 .or. config%weight(i) < 1.0e-6_wp) Then
                config%vxx(i2) = 0.0_wp
                config%vyy(i2) = 0.0_wp
                config%vzz(i2) = 0.0_wp
              End If
            End If

            ! get new rotation matrix

            Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

            ! update RB members new velocities

            lrgd = rigid%list(-1, irgd)
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

      Else ! no RBs present in the system

        ! generate starting velocities

        Do i = 1, config%natms

          ! frozen and massless atoms are either motionless
          ! (with no actual DoF - relaxed shells, frozen sites)
          ! The rest are to have gaussian distribution of their Ekin
          ! (adiabatic shells, constraints, PMFs are
          ! to be sorted out later by quenching)

          If (qn(i) == 1) Then
            Call box_mueller_saru3(seed, config%ltg(i), 0, config%vxx(i), config%vyy(i), config%vzz(i))

            ! Get scaler to target variance/Sqrt(config%weight)

            tmp = 1.0_wp / Sqrt(config%weight(i))
            config%vxx(i) = config%vxx(i) * tmp
            config%vyy(i) = config%vyy(i) * tmp
            config%vzz(i) = config%vzz(i) * tmp
          End If
        End Do

      End If

      If (stp > 0) Then
        If (cshell%lshmv_shl) Then ! refresh the q array for shared core-shell units
          qn(config%natms + 1:config%nlast) = 0
          Call update_shared_units_int(config, cshell%lishp_shl, &
                                       cshell%lashp_shl, qn, domain, comm)
          Call update_shared_units(config, cshell%lishp_shl, &
                                   cshell%lashp_shl, config%vxx, config%vyy, config%vzz, domain, comm)
        End If

        If (tps(comm%idnode) > 0) Then
          Do k = 1, cshell%ntshl
            If (qs(0, k) == 1) Then
              i1 = qs(1, k)
              i2 = qs(2, k)

              config%vxx(i2) = config%vxx(i1)
              config%vyy(i2) = config%vyy(i1)
              config%vzz(i2) = config%vzz(i1)
            End If
          End Do
        End If
      End If

      Deallocate (qn, tpn, Stat=fail(1))
      Deallocate (qs, tps, Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message, '(a)') 'set_temperature deallocation failure'
        Call error(0, message)
      End If

      ! remove centre of mass motion

      If (rigid%total > 0) Then
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
      Else
        Call getvom(vom, config, comm)

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) - vom(1)
            config%vyy(i) = config%vyy(i) - vom(2)
            config%vzz(i) = config%vzz(i) - vom(3)
          End If
        End Do
      End If

    Else

      ! quench RBs

      If (rigid%total > 0) Then
        Call rigid_bodies_quench(rigid, domain, config, comm)
      End If

    End If

    ! config%levcfg must be equalised to ONE (velocities now exist)

    config%levcfg = 1

    ! scale velocities for keyres = RESTART_KEY_CLEAN and keyres = RESTART_KEY_SCALE

    If (keyres == RESTART_KEY_CLEAN .or. keyres == RESTART_KEY_SCALE) Then

      ! Detect pure molecular statics only == CGM minimisation at zero timestep
      ! and no rescaling

      no_min_0 = .not. (minim%minimise .and. minim%freq == 0 .and. nstep == 0 .and. nstrun == 0 .and. keyres == RESTART_KEY_SCALE)

      ! quench constraints & PMFs

      If (no_min_0) Then
        If (cons%megcon > 0) Then
          Call constraints_quench(cons, stat, domain, config, comm)
        End If
        If (pmf%megpmf > 0) Call pmf_quench(cons%max_iter_shake, cons%tolerance, stat, pmf, config, comm)
      End If

      ! quench core-shell units in adiabatic model

      If (cshell%megshl > 0 .and. cshell%keyshl == SHELL_ADIABATIC .and. no_min_0) Then
        Do
          Call scale_temperature(thermo%sigma, config%degtra, config%degrot, config%degfre, rigid, config, comm)
          Call core_shell_quench(config, safe, thermo%temp, cshell, domain, comm)
          If (cons%megcon > 0) Then
            Call constraints_quench(cons, stat, domain, config, comm)
          End If
          If (pmf%megpmf > 0) Then
            Call pmf_quench(cons%max_iter_shake, cons%tolerance, stat, pmf, config, comm)
          End If
          If (rigid%total > 0) Then
            Call rigid_bodies_quench(rigid, domain, config, comm)
          End If
          If (safe) Exit
        End Do
      Else
        Call scale_temperature(thermo%sigma, config%degtra, config%degrot, config%degfre, rigid, config, comm)
      End If

    End If

    If (config%l_vom) Then

      ! remove centre of mass motion

      If (rigid%total > 0) Then
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
      Else
        Call getvom(vom, config, comm)

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) - vom(1)
            config%vyy(i) = config%vyy(i) - vom(2)
            config%vzz(i) = config%vzz(i) - vom(3)
          End If
        End Do
      End If
    Else ! make getvom always return 0, which is not good for
      Call config%chvom() ! standard MD as the flying ice-cub effect may happen
    End If ! and/or T(MD) is loses its microscopic meaning!

    ! Initialise engrot and if RBs exist calculate it

    engrot = 0.0_wp
    If (rigid%total > 0) engrot = getknr(rigid, comm)

  End Subroutine set_temperature

  Subroutine regauss_temperature(rigid, domain, config, seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to regauss the instantaneous system temperature
    ! by random pairwise swaps of the energy scaled momenta of dynamically
    ! active particles (no massless shells or massless RB members)
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(seed_type),          Intent(InOut) :: seed
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)   :: message
    Integer              :: fail, i, irgd, is, j, jrgd, k, l, lrgd, rgdtyp
    Integer, Allocatable :: ind(:), pair(:, :)
    Real(Kind=wp)        :: tmp, vom(1:3)

    fail = 0
    Allocate (ind(1:config%natms), pair(1:2, config%natms / 2), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'regauss_temperature allocation failure'
      Call error(0, message)
    End If

    ! Create and index array containing the indices of the
    ! dynamically active particles and zeros for the inactive

    Do i = 1, config%natms
      If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
        ind(i) = i
      Else
        ind(i) = 0
      End If
    End Do

    ! Compress the index array

    i = 1
    j = config%natms
    Do While (i < j)
      Do While (ind(j) == 0 .and. j > i)
        j = j - 1
      End Do

      If (i < j) Then
        If (ind(i) == 0) Then
          ind(i) = ind(j)
          ind(j) = 0
          j = j - 1
        Else
          i = i + 1
        End If
      End If
    End Do
    k = j

    ! Create a non-overlapping array of random pairs
    ! by exhausting the index array

    Do i = 1, k / 2
      Do l = 1, 2
        is = 1 + Int(Real(j, wp) * uni(seed, comm))
        pair(l, i) = ind(is)
        ind(is) = ind(j)
        ind(j) = 0
        j = j - 1
      End Do
    End Do

    ! Swap particles energies in the pair array

    Do i = 1, k / 2
      j = pair(1, i)
      l = pair(2, i)

      vom(1) = config%vxx(j)
      vom(2) = config%vyy(j)
      vom(3) = config%vzz(j)

      tmp = Sqrt(config%weight(l) / config%weight(j))
      config%vxx(j) = config%vxx(l) * tmp
      config%vyy(j) = config%vyy(l) * tmp
      config%vzz(j) = config%vzz(l) * tmp

      tmp = 1.0_wp / tmp
      config%vxx(l) = vom(1) * tmp
      config%vyy(l) = vom(2) * tmp
      config%vzz(l) = vom(3) * tmp
    End Do

    If (rigid%total > 0) Then

      ! quench RBs

      Call rigid_bodies_quench(rigid, domain, config, comm)

      ! remove centre of mass motion

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

    Else

      ! remove centre of mass motion

      Call getvom(vom, config, comm)

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do

    End If

    Deallocate (ind, pair, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'regauss_temperature deallocation failure'
      Call error(0, message)
    End If

  End Subroutine regauss_temperature

  Subroutine scale_temperature(sigma, degtra, degrot, degfre, rigid, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to scale the instantaneous system temperature
    ! to the target temperature
    !
    ! Note: zeroes angular momentum in non-periodic systems, frozen
    ! particles are not considered
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),            Intent(In   ) :: sigma
    Integer(Kind=li),         Intent(In   ) :: degtra, degrot, degfre
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Integer                                  :: fail, i, irgd, j, jrgd, lrgd, rgdtyp
    Real(Kind=wp)                            :: amx, amy, amz, com(1:3), engke, engkf, engkt, &
                                                engrot, rot(1:9), rotinv(1:9), tmp, tmp1, &
                                                vom(1:3), wxx, wyy, wzz, x, y, z
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer

    ! remove centre of mass motion

    If (rigid%total > 0) Then
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
    Else
      Call getvom(vom, config, comm)

      Do i = 1, config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
          config%vxx(i) = config%vxx(i) - vom(1)
          config%vyy(i) = config%vyy(i) - vom(2)
          config%vzz(i) = config%vzz(i) - vom(3)
        End If
      End Do
    End If

    ! zero angular momentum about centre of mass - non-periodic system

    If (config%imcon == IMCON_NOPBC) Then
      fail = 0
      Allocate (buffer(1:12), Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'scale_temperature allocation failure'
        Call error(0, message)
      End If

      ! calculate centre of mass position

      Call getcom(config, com, comm)

      If (rigid%total > 0) Then

        ! move to centre of mass origin

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx - com(1)
            config%parts(i)%yyy = config%parts(i)%yyy - com(2)
            config%parts(i)%zzz = config%parts(i)%zzz - com(3)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            rigid%xxx(irgd) = rigid%xxx(irgd) - com(1)
            rigid%yyy(irgd) = rigid%yyy(irgd) - com(2)
            rigid%zzz(irgd) = rigid%zzz(irgd) - com(3)
          End If
        End Do

        ! angular momentum accumulators

        amx = 0.0_wp
        amy = 0.0_wp
        amz = 0.0_wp

        ! rotational inertia accumulators

        rot = 0.0_wp

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            amx = amx + config%weight(i) * (config%parts(i)%yyy * config%vzz(i) - config%parts(i)%zzz * config%vyy(i))
            amy = amy + config%weight(i) * (config%parts(i)%zzz * config%vxx(i) - config%parts(i)%xxx * config%vzz(i))
            amz = amz + config%weight(i) * (config%parts(i)%xxx * config%vyy(i) - config%parts(i)%yyy * config%vxx(i))

            tmp = config%parts(i)%xxx**2 + config%parts(i)%yyy**2 + config%parts(i)%zzz**2
            rot(1) = rot(1) + config%weight(i) * (config%parts(i)%xxx * config%parts(i)%xxx - tmp)
            rot(2) = rot(2) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%yyy
            rot(3) = rot(3) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%zzz
            rot(5) = rot(5) + config%weight(i) * (config%parts(i)%yyy * config%parts(i)%yyy - tmp)
            rot(6) = rot(6) + config%weight(i) * config%parts(i)%yyy * config%parts(i)%zzz
            rot(9) = rot(9) + config%weight(i) * (config%parts(i)%zzz * config%parts(i)%zzz - tmp)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            lrgd = rigid%list(-1, irgd)

            tmp1 = rigid%weight(0, rgdtyp) * Real(rigid%index_local(0, irgd), wp) / Real(lrgd, wp)

            amx = amx + tmp1 * (rigid%yyy(irgd) * rigid%vzz(irgd) - rigid%zzz(irgd) * rigid%vyy(irgd))
            amy = amy + tmp1 * (rigid%zzz(irgd) * rigid%vxx(irgd) - rigid%xxx(irgd) * rigid%vzz(irgd))
            amz = amz + tmp1 * (rigid%xxx(irgd) * rigid%vyy(irgd) - rigid%yyy(irgd) * rigid%vxx(irgd))

            tmp = rigid%xxx(irgd)**2 + rigid%yyy(irgd)**2 + rigid%zzz(irgd)**2

            rot(1) = rot(1) + tmp1 * (rigid%xxx(irgd) * rigid%xxx(irgd) - tmp)
            rot(2) = rot(2) + tmp1 * rigid%xxx(irgd) * rigid%yyy(irgd)
            rot(3) = rot(3) + tmp1 * rigid%xxx(irgd) * rigid%zzz(irgd)
            rot(5) = rot(5) + tmp1 * (rigid%yyy(irgd) * rigid%yyy(irgd) - tmp)
            rot(6) = rot(6) + tmp1 * rigid%yyy(irgd) * rigid%zzz(irgd)
            rot(9) = rot(9) + tmp1 * (rigid%zzz(irgd) * rigid%zzz(irgd) - tmp)
          End If
        End Do

        ! complete rotational inertia matrix

        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)

        ! global sum of rotation

        buffer(1) = amx
        buffer(2) = amy
        buffer(3) = amz
        Do i = 1, 9
          buffer(i + 3) = rot(i)
        End Do

        Call gsum(comm, buffer)

        amx = buffer(1)
        amy = buffer(2)
        amz = buffer(3)
        Do i = 1, 9
          rot(i) = buffer(i + 3)
        End Do

        ! invert rotational inertia matrix

        Call invert(rot, rotinv, tmp)

        ! correction to angular velocity

        wxx = rotinv(1) * amx + rotinv(2) * amy + rotinv(3) * amz
        wyy = rotinv(4) * amx + rotinv(5) * amy + rotinv(6) * amz
        wzz = rotinv(7) * amx + rotinv(8) * amy + rotinv(9) * amz

        ! correction to linear velocity

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) + (wyy * config%parts(i)%zzz - wzz * config%parts(i)%yyy)
            config%vyy(i) = config%vyy(i) + (wzz * config%parts(i)%xxx - wxx * config%parts(i)%zzz)
            config%vzz(i) = config%vzz(i) + (wxx * config%parts(i)%yyy - wyy * config%parts(i)%xxx)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            x = (wyy * rigid%zzz(irgd) - wzz * rigid%yyy(irgd))
            y = (wzz * rigid%xxx(irgd) - wxx * rigid%zzz(irgd))
            z = (wxx * rigid%yyy(irgd) - wyy * rigid%xxx(irgd))

            rigid%vxx(irgd) = rigid%vxx(irgd) + x
            rigid%vyy(irgd) = rigid%vyy(irgd) + y
            rigid%vzz(irgd) = rigid%vzz(irgd) + z

            lrgd = rigid%list(-1, irgd)
            Do jrgd = 1, lrgd
              i = rigid%index_local(jrgd, irgd) ! local index of particle/site

              If (i <= config%natms) Then
                config%vxx(i) = config%vxx(i) + x
                config%vyy(i) = config%vyy(i) + y
                config%vzz(i) = config%vzz(i) + z
              End If
            End Do
          End If
        End Do

        ! reset positions to original reference frame

        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx + com(1)
            config%parts(i)%yyy = config%parts(i)%yyy + com(2)
            config%parts(i)%zzz = config%parts(i)%zzz + com(3)
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          If (rigid%frozen(0, rgdtyp) == 0) Then
            rigid%xxx(irgd) = rigid%xxx(irgd) + com(1)
            rigid%yyy(irgd) = rigid%yyy(irgd) + com(2)
            rigid%zzz(irgd) = rigid%zzz(irgd) + com(3)
          End If
        End Do

      Else

        ! move to centre of mass origin

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx - com(1)
            config%parts(i)%yyy = config%parts(i)%yyy - com(2)
            config%parts(i)%zzz = config%parts(i)%zzz - com(3)
          End If
        End Do

        ! angular momentum accumulators

        amx = 0.0_wp
        amy = 0.0_wp
        amz = 0.0_wp

        ! rotational inertia accumulators

        rot = 0.0_wp

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            amx = amx + config%weight(i) * (config%parts(i)%yyy * config%vzz(i) - config%parts(i)%zzz * config%vyy(i))
            amy = amy + config%weight(i) * (config%parts(i)%zzz * config%vxx(i) - config%parts(i)%xxx * config%vzz(i))
            amz = amz + config%weight(i) * (config%parts(i)%xxx * config%vyy(i) - config%parts(i)%yyy * config%vxx(i))

            tmp = config%parts(i)%xxx**2 + config%parts(i)%yyy**2 + config%parts(i)%zzz**2
            rot(1) = rot(1) + config%weight(i) * (config%parts(i)%xxx * config%parts(i)%xxx - tmp)
            rot(2) = rot(2) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%yyy
            rot(3) = rot(3) + config%weight(i) * config%parts(i)%xxx * config%parts(i)%zzz
            rot(5) = rot(5) + config%weight(i) * (config%parts(i)%yyy * config%parts(i)%yyy - tmp)
            rot(6) = rot(6) + config%weight(i) * config%parts(i)%yyy * config%parts(i)%zzz
            rot(9) = rot(9) + config%weight(i) * (config%parts(i)%zzz * config%parts(i)%zzz - tmp)
          End If
        End Do

        ! complete rotational inertia matrix

        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)

        ! global sum of rotation

        buffer(1) = amx
        buffer(2) = amy
        buffer(3) = amz
        Do i = 1, 9
          buffer(i + 3) = rot(i)
        End Do

        Call gsum(comm, buffer)

        amx = buffer(1)
        amy = buffer(2)
        amz = buffer(3)
        Do i = 1, 9
          rot(i) = buffer(i + 3)
        End Do

        ! invert rotational inertia matrix

        Call invert(rot, rotinv, tmp)

        ! correction to angular velocity

        wxx = rotinv(1) * amx + rotinv(2) * amy + rotinv(3) * amz
        wyy = rotinv(4) * amx + rotinv(5) * amy + rotinv(6) * amz
        wzz = rotinv(7) * amx + rotinv(8) * amy + rotinv(9) * amz

        ! correction to linear velocity

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) + (wyy * config%parts(i)%zzz - wzz * config%parts(i)%yyy)
            config%vyy(i) = config%vyy(i) + (wzz * config%parts(i)%xxx - wxx * config%parts(i)%zzz)
            config%vzz(i) = config%vzz(i) + (wxx * config%parts(i)%yyy - wyy * config%parts(i)%xxx)
          End If
        End Do

        ! reset positions to original reference frame

        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%parts(i)%xxx = config%parts(i)%xxx + com(1)
            config%parts(i)%yyy = config%parts(i)%yyy + com(2)
            config%parts(i)%zzz = config%parts(i)%zzz + com(3)
          End If
        End Do

      End If

      Deallocate (buffer, Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'scale_temperature deallocation failure'
        Call error(0, message)
      End If
    End If

    ! ensure equipartitioning - all degrees of freedom are equal
    ! calculate energy: free particles and RB translational and rotational

    engrot = 0.0_wp
    If (rigid%total > 0) Then
      engkf = getknf(config%vxx, config%vyy, config%vzz, config, comm)
      engkt = getknt(rigid, comm)

      engrot = getknr(rigid, comm)

      ! temporary replacement for small engrot

      tmp1 = Max(engrot, 1.0e-6_wp)

      ! Scale rotational energy to translational energy
      ! according to their respective DoF

      If (degtra > Int(0, li)) Then ! engkt > 0 (degrot > 0 and tmp1(engrot) > 0)
        tmp = Sqrt((engkt * Real(degrot, wp)) / (tmp1 * Real(degtra, wp)))

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then
            rigid%oxx(irgd) = rigid%oxx(irgd) * tmp
            rigid%oyy(irgd) = rigid%oyy(irgd) * tmp
            rigid%ozz(irgd) = rigid%ozz(irgd) * tmp
          End If
        End Do
        engrot = engrot * tmp**2
        tmp1 = Max(engrot, 1.0e-6_wp)
      End If

      ! Scale the energy per DoF of the RBs to that of a DoF of the free particles

      If (degfre - degtra - degrot > Int(0, li)) Then ! engkf > 0
        tmp = Sqrt((engkf * Real(degrot, wp)) / (tmp1 * Real(degfre - degtra - degrot, wp)))

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then
            rigid%oxx(irgd) = rigid%oxx(irgd) * tmp
            rigid%oyy(irgd) = rigid%oyy(irgd) * tmp
            rigid%ozz(irgd) = rigid%ozz(irgd) * tmp

            If (rigid%frozen(0, rgdtyp) == 0) Then
              rigid%vxx(irgd) = rigid%vxx(irgd) * tmp
              rigid%vyy(irgd) = rigid%vyy(irgd) * tmp
              rigid%vzz(irgd) = rigid%vzz(irgd) * tmp
            End If
          End If
        End Do
        engrot = engrot * tmp**2
        If (degtra > Int(0, li)) engkt = engkt * tmp**2
      End If

      engke = engkf + engkt
    Else
      engke = getkin(config, config%vxx, config%vyy, config%vzz, comm)
    End If

    ! apply temperature scaling

    If (engke + engrot > 1.0e-6_wp .and. sigma > zero_plus) Then
      tmp = Sqrt(sigma / (engke + engrot))

      If (rigid%total > 0) Then
        Do j = 1, config%nfree
          i = config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) * tmp
            config%vyy(i) = config%vyy(i) * tmp
            config%vzz(i) = config%vzz(i) * tmp
          End If
        End Do

        Do irgd = 1, rigid%n_types
          rgdtyp = rigid%list(0, irgd)

          lrgd = rigid%list(-1, irgd)
          If (rigid%frozen(0, rgdtyp) < lrgd) Then

            ! new angular velocity

            rigid%oxx(irgd) = rigid%oxx(irgd) * tmp
            rigid%oyy(irgd) = rigid%oyy(irgd) * tmp
            rigid%ozz(irgd) = rigid%ozz(irgd) * tmp

            ! new translational velocity

            If (rigid%frozen(0, rgdtyp) == 0) Then
              rigid%vxx(irgd) = rigid%vxx(irgd) * tmp
              rigid%vyy(irgd) = rigid%vyy(irgd) * tmp
              rigid%vzz(irgd) = rigid%vzz(irgd) * tmp
            End If

            ! new rotational matrix

            Call getrotmat(rigid%q0(irgd), rigid%q1(irgd), rigid%q2(irgd), rigid%q3(irgd), rot)

            Do jrgd = 1, lrgd
              If (rigid%frozen(jrgd, rgdtyp) == 0) Then ! Apply restrictions
                i = rigid%index_local(jrgd, irgd) ! local index of particle/site

                If (i <= config%natms) Then
                  x = rigid%x(jrgd, rgdtyp)
                  y = rigid%y(jrgd, rgdtyp)
                  z = rigid%z(jrgd, rgdtyp)

                  ! site velocity in body frame

                  wxx = rigid%oyy(irgd) * z - rigid%ozz(irgd) * y
                  wyy = rigid%ozz(irgd) * x - rigid%oxx(irgd) * z
                  wzz = rigid%oxx(irgd) * y - rigid%oyy(irgd) * x

                  ! new atomic velocities in lab frame

                  config%vxx(i) = rot(1) * wxx + rot(2) * wyy + rot(3) * wzz + rigid%vxx(irgd)
                  config%vyy(i) = rot(4) * wxx + rot(5) * wyy + rot(6) * wzz + rigid%vyy(irgd)
                  config%vzz(i) = rot(7) * wxx + rot(8) * wyy + rot(9) * wzz + rigid%vzz(irgd)
                End If
              End If
            End Do
          End If
        End Do
      Else
        Do i = 1, config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            config%vxx(i) = config%vxx(i) * tmp
            config%vyy(i) = config%vyy(i) * tmp
            config%vzz(i) = config%vzz(i) * tmp
          End If
        End Do
      End If
    Else ! sigma must be zero
      Do i = 1, config%natms
        config%vxx(i) = 0.0_wp; config%vyy(i) = 0.0_wp; config%vzz(i) = 0.0_wp
      End Do

      If (rigid%total > 0) Then
        Do irgd = 1, rigid%n_types
          rigid%vxx(irgd) = 0.0_wp; rigid%vyy(irgd) = 0.0_wp; rigid%vzz(irgd) = 0.0_wp
          rigid%oxx(irgd) = 0.0_wp; rigid%oyy(irgd) = 0.0_wp; rigid%ozz(irgd) = 0.0_wp
        End Do
      End If
    End If
  End Subroutine scale_temperature
End Module temperature
