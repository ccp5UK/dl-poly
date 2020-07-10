Module two_body
  Use comms,           Only: comms_type,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: pi,&
                             r4pie0
  Use coul_mpole,      Only: coul_chrm_forces,&
                             coul_cp_mforces,&
                             coul_dddp_mforces,&
                             coul_fscp_mforces,&
                             coul_rfp_mforces,&
                             d_ene_trq_mpoles
  Use coul_spole,      Only: coul_cp_forces,&
                             coul_dddp_forces,&
                             coul_fscp_forces,&
                             coul_rfp_forces
  Use domains,         Only: domains_type
  Use electrostatic,   Only: ELECTROSTATIC_COULOMB,&
                             ELECTROSTATIC_COULOMB_FORCE_SHIFT,&
                             ELECTROSTATIC_COULOMB_REACTION_FIELD,&
                             ELECTROSTATIC_DDDP,&
                             ELECTROSTATIC_EWALD,&
                             ELECTROSTATIC_NULL,&
                             ELECTROSTATIC_POISSON,&
                             electrostatic_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc
  Use ewald,           Only: ewald_type,&
                             ewald_vdw_coeffs,&
                             ewald_vdw_count,&
                             ewald_vdw_init
  Use ewald_spole,     Only: ewald_excl_forces,&
                             ewald_real_forces_coul,&
                             ewald_spme_forces_coul
  Use ewald_general,   Only: ewald_real_forces_gen,&
                             ewald_spme_forces_gen
  Use kim,             Only: kim_energy_and_forces,&
                             kim_type
  Use kinds,           Only: wp
  Use metal,           Only: metal_forces,&
                             metal_ld_compute,&
                             metal_lrc,&
                             metal_type
  Use mpole,           Only: POLARISATION_CHARMM,&
                             mpole_type
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: calc_erfc,&
                             calc_erfc_deriv
  Use poisson,         Only: poisson_forces,&
                             poisson_type
  Use rdfs,            Only: rdf_collect,&
                             rdf_excl_collect,&
                             rdf_frzn_collect,&
                             rdf_increase_block_number,&
                             rdf_type
  Use site,            Only: site_type
  Use spme,            Only: init_spme_data,&
                             spme_self_interaction
  Use statistics,      Only: stats_type
  Use timer,           Only: start_timer,&
                             stop_timer,&
                             timer_type
  Use vdw,             Only: vdw_forces_tab,&
                             vdw_forces_direct,&
                             vdw_type

  Implicit None

  Private

  Public :: two_body_forces
Contains

  Subroutine two_body_forces(ensemble, lbook, megfrz, leql, nsteql, nstep, &
                             stats, ewld, met, pois, neigh, sites, vdws, rdf, mpoles, electro, &
                             domain, tmr, kim_data, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating interatomic forces and rdf%rdf
    ! using the verlet neighbour list
    !
    ! vdws%n_vdw > 0 ------ switch for vdw potentials calculation
    ! met%n_potentials > 0 ------ switch for metal local density and potentials
    !                   calculations
    !
    ! nstfce - the rate at which the k-space contributions of SPME are
    !          refreshed.  Once every 1 <= nstfce <= 7 steps.
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2017
    ! contrib   - h.a.boateng february 2016
    ! contrib   - p.s.petkov february 2015
    ! contrib   - a.b.g.chalk january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent(In) :: ensemble
    Logical, Intent(In) :: lbook
    Integer, Intent(In) :: megfrz
    Logical, Intent(In) :: leql
    Integer, Intent(In) :: nsteql, nstep
    Type(stats_type), Intent(InOut) :: stats
    Type(ewald_type), Intent(InOut) :: ewld
    Type(metal_type), Intent(InOut) :: met
    Type(poisson_type), Intent(InOut) :: pois
    Type(neighbours_type), Intent(InOut) :: neigh
    Type(site_type), Intent(In) :: sites
    Type(vdw_type), Intent(InOut) :: vdws
    Type(rdf_type), Intent(InOut) :: rdf
    Type(mpole_type), Intent(InOut) :: mpoles
    Type(electrostatic_type), Intent(InOut) :: electro
    Type(domains_type), Intent(In) :: domain
    Type(timer_type), Intent(InOut) :: tmr
    Type(kim_type), Intent(InOut) :: kim_data
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type), Intent(InOut) :: comm

    Integer                                     :: fail, i, ipot, j, k, limit
    Logical                                     :: l_do_outer_loop, l_do_rdf, safe
    Logical, Save                               :: newjob = .true.
    Real(Kind=wp)                               :: buffer(0:19), engacc, engcpe_ch, engcpe_ex, &
                                                   engcpe_fr, engcpe_nz, engcpe_rc, engcpe_rl, &
                                                   engden, engkim, engmet, engvdw, engvdw_rc, &
                                                   engvdw_rl, factor_nz, tmp, viracc, vircpe_ch, &
                                                   vircpe_dt, vircpe_ex, vircpe_fr, vircpe_nz, &
                                                   vircpe_rc, vircpe_rl, virden, virkim, virmet, &
                                                   virvdw, virvdw_rc, virvdw_rl
    Real(kind=wp), Allocatable, Dimension(:)    :: coul_coeffs, rrt, xxt, yyt, zzt
    Real(kind=wp), Allocatable, Dimension(:, :) :: vdw_coeffs

! Array to remap charge/disp to 2D array

#ifdef CHRONO
    Call start_timer(tmr, 'Two-Body Init')
#endif

    safe = .true.
    fail = 0

    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    If (fail > 0) Call error_alloc('distance arrays', 'two_body_forces')

    l_do_rdf = (rdf%l_collect .and. ((.not. leql) .or. nstep >= nsteql) .and. Mod(nstep, rdf%freq) == 0)

    ! If k-space SPME is evaluated infrequently check whether
    ! at this timestep to evaluate or "refresh" with old values.
    ! At restart allocate the "refresh" arrays and force a fresh
    ! evaluation.  Repeat the same but only for the SPME k-space
    ! frozen-frozen evaluations in constant volume ensembles only.

    ! Coulomb
    If (ewld%vdw .and. .not. ewld%active) Call error(0, 'Ewald VdW requested but ewald not enabled')

    If (ewld%active) Then
      Allocate (coul_coeffs(config%mxatms), stat=fail)
      If (fail > 0) Call error_alloc('coul_coeffs', 'two_body_forces')
      coul_coeffs = config%parts(:)%chge


      If (newjob) Then

        If (.not. ewld%direct) Then
          Call electro%erfcgen(neigh%cutoff, ewld%alpha)
        End If

        ! Assume no VDW
        ewld%num_pots = 0

        ! Set number of pots
        If (ewld%vdw) Call ewald_vdw_count(ewld, vdws)

        If (electro%key == ELECTROSTATIC_EWALD) Then
          Allocate (ewld%spme_data(0:ewld%num_pots), stat=fail)
          If (fail > 0) Call error_alloc('ewld%spme_data', 'two_body_forces')

          ewld%spme_data(0)%scaling = r4pie0 / electro%eps
          Call init_spme_data(ewld%spme_data(0), 1)

          If (electro%multipolar) Then
            Call spme_self_interaction(ewld%alpha, config%natms, coul_coeffs, comm, ewld%spme_data(0), electro%mpoles)
          Else
            Call spme_self_interaction(ewld%alpha, config%natms, coul_coeffs, comm, ewld%spme_data(0))
          End If

        Else
          Allocate (ewld%spme_data(1:ewld%num_pots), stat=fail)
          If (fail > 0) Call error_alloc('ewld%spme_data', 'two_body_forces')
        End If

        If (ewld%vdw) Then
          Call ewald_vdw_init(ewld, vdws)
          Call ewald_vdw_coeffs(config, vdws, ewld, vdw_coeffs)

          Do ipot = 1, ewld%num_pots
            Call spme_self_interaction(ewld%alpha, config%natms, vdw_coeffs(:, ipot), comm, ewld%spme_data(ipot))
          End Do

          ! Disable long-range corrections (calculating long range explicitly)
          vdws%elrc = 0.0_wp
          vdws%vlrc = 0.0_wp

        End If

        newjob = .false.

      End If

      If (ewld%vdw) Call ewald_vdw_coeffs(config, vdws, ewld, vdw_coeffs)

    End If

    ! initialise energy and virial accumulators

    engkim = 0.0_wp
    virkim = 0.0_wp

    engden = 0.0_wp
    virden = 0.0_wp

    engmet = 0.0_wp
    virmet = 0.0_wp

    engvdw = 0.0_wp
    virvdw = 0.0_wp
    engvdw_rc = 0.0_wp
    virvdw_rc = 0.0_wp
    engvdw_rl = 0.0_wp
    virvdw_rl = 0.0_wp

    engvdw = 0.0_wp
    virvdw = 0.0_wp

    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    engcpe_rl = 0.0_wp
    vircpe_rl = 0.0_wp

    engcpe_ch = 0.0_wp
    vircpe_ch = 0.0_wp

    engcpe_ex = 0.0_wp
    vircpe_ex = 0.0_wp

    engcpe_fr = 0.0_wp
    vircpe_fr = 0.0_wp

    engcpe_nz = 0.0_wp
    vircpe_nz = 0.0_wp

    vircpe_dt = 0.0_wp

#ifdef CHRONO
    Call stop_timer(tmr, 'Two-Body Init')
#endif

    ! Calculate all contributions from KIM
    If (kim_data%active) Then
#ifdef CHRONO
      Call start_timer(tmr, 'KIM')
#endif
      Call kim_energy_and_forces(kim_data, config%natms, config%nlast, config%parts, &
                                 neigh%list, domain%map, config%lsite, config%ltype, config%lsi, config%lsa, config%ltg, &
                                 sites%site_name, engkim, virkim, stats%stress, comm)
#ifdef CHRONO
      Call stop_timer(tmr, 'KIM')
#endif
    End If

    If (met%n_potentials > 0) Then

      ! Reset metal long-range corrections (constant pressure/stress only)

      If (ensemble >= 20) Call metal_lrc(met, sites, config, comm)

      ! calculate local density in metals

      Call metal_ld_compute(engden, virden, stats%stress, sites%ntype_atom, met, neigh, &
                            domain, config, comm)
    End If

    ! calculate coulombic forces, Ewald sum - fourier contribution
#ifdef CHRONO
    Call start_timer(tmr, 'Long Range')
#endif

    If (electro%key == ELECTROSTATIC_EWALD) Then

      Call ewald_spme_forces_coul(ewld, ewld%spme_data(0), electro, domain, config, comm, &
        & coul_coeffs, stats, engcpe_rc, vircpe_rc, tmr)

    End If

    If (ewld%vdw) Then
#ifdef CHRONO
      Call start_timer(tmr, 'SPME Order-n')
#endif
      Do ipot = 1, ewld%num_pots

        Call ewald_spme_forces_gen(ewld, ewld%spme_data(ipot), electro, domain, config, comm, &
          & vdw_coeffs(:, ipot), stats, engacc, viracc, tmr)

        engvdw_rc = engvdw_rc + engacc
        virvdw_rc = virvdw_rc + viracc

      End Do

#ifdef CHRONO
      Call stop_timer(tmr, 'SPME Order-n')
#endif
    End If

#ifdef CHRONO
    Call stop_timer(tmr, 'Long Range')
    Call start_timer(tmr, 'Short Range')
#endif

#ifdef KIM
    l_do_outer_loop = ((met%n_potentials > 0) .or. &
                       (vdws%n_vdw > 0) .or. &
                       (mpoles%max_mpoles > 0) .or. &
                       ((electro%key /= ELECTROSTATIC_NULL) .and. &
                        (electro%key /= ELECTROSTATIC_POISSON)) .or. &
                       l_do_rdf)

    If (l_do_outer_loop) Then
      ! outer loop over atoms
#endif
      Do i = 1, config%natms

        ! Get neigh%list limit
        limit = neigh%list(0, i)

        ! calculate interatomic distances

        Do k = 1, limit
          j = neigh%list(k, i)
          xxt(k) = config%parts(i)%xxx - config%parts(j)%xxx
          yyt(k) = config%parts(i)%yyy - config%parts(j)%yyy
          zzt(k) = config%parts(i)%zzz - config%parts(j)%zzz
        End Do

        ! periodic boundary conditions not needed by LC construction
        !
        !     Call images(imcon,cell,limit,xxt,yyt,zzt)

        ! distances, thanks to Alin Elena (one too many changes)

        Do k = 1, limit
          rrt(k) = Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
        End Do

        ! calculate metal forces and potential

        If (met%n_potentials > 0) Then
          Call metal_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, safe, sites%ntype_atom, met, neigh, config)

          engmet = engmet + engacc
          virmet = virmet + viracc
        End If

        ! calculate short-range force and potential terms

        If (vdws%n_vdw > 0) Then
          If (ewld%vdw) Then
            Do ipot = 1, ewld%num_pots

              Call ewald_real_forces_gen(ewld%alpha, ewld%spme_data(ipot), neigh, config, stats, &
                   & vdw_coeffs(:, ipot), i, xxt, yyt, zzt, rrt, engacc, viracc)

              engvdw_rl = engvdw_rl + engacc
              virvdw_rl = virvdw_rl + viracc

            End Do

          Else If (vdws%l_direct) Then ! direct calculation

            Call vdw_forces_direct(i, xxt, yyt, zzt, rrt, engacc, viracc, stats, neigh, vdws, config)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc

          else

            Call vdw_forces_tab(i, xxt, yyt, zzt, rrt, engacc, viracc, stats, neigh, vdws, config)
            engvdw = engvdw + engacc
            virvdw = virvdw + viracc

          end If

        End If

        !-------------------
        ! COULOMBIC CONTRIBUTIONS
        !-------------------
        If (.not. electro%no_elec) then
          If (mpoles%max_mpoles > 0) Then

            Select Case (electro%key)
            Case (ELECTROSTATIC_EWALD)

              If (ewld%direct) Then
                Call ewald_real_forces_gen(ewld%alpha, ewld%spme_data(0), neigh, config, stats, &
                     & coul_coeffs, i, xxt, yyt, zzt, rrt, engacc, viracc)
              Else
                Call ewald_real_forces_coul(electro, ewld%alpha, ewld%spme_data(0), neigh, config, stats, &
                     & i, xxt, yyt, zzt, rrt, engacc, viracc)
              End If

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_DDDP)

              ! distance dependant dielectric potential

              Call coul_dddp_mforces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_COULOMB)

              ! coulombic 1/r potential with no truncation or damping

              Call coul_cp_mforces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case(ELECTROSTATIC_COULOMB_FORCE_SHIFT)

              ! force-shifted coulomb potentials

              Call coul_fscp_mforces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, electro, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_COULOMB_REACTION_FIELD)

              ! reaction field potential

              Call coul_rfp_mforces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, electro, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            End Select

          Else

            Select Case (electro%key)
            Case (ELECTROSTATIC_EWALD)

              ! calculate coulombic forces, Ewald sum - real space contribution

              If (ewld%direct) Then
                Call ewald_real_forces_gen(ewld%alpha, ewld%spme_data(0), neigh, config, stats, &
                     & coul_coeffs, i, xxt, yyt, zzt, rrt, engacc, viracc)
              Else
                Call ewald_real_forces_coul(electro, ewld%alpha, ewld%spme_data(0), neigh, config, stats, &
                     & i, xxt, yyt, zzt, rrt, engacc, viracc)
              End If

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_DDDP)

              ! distance dependant dielectric potential

              Call coul_dddp_forces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats, neigh, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_COULOMB)

              ! coulombic 1/r potential with no truncation or damping

              Call coul_cp_forces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats, neigh, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_COULOMB_FORCE_SHIFT)

              ! force-shifted coulomb potentials

              Call coul_fscp_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats, neigh, electro, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            Case (ELECTROSTATIC_COULOMB_REACTION_FIELD)

              ! reaction field potential

              Call coul_rfp_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats, neigh, electro, config)

              engcpe_rl = engcpe_rl + engacc
              vircpe_rl = vircpe_rl + viracc

            End Select

          End If
        End If

        ! accumulate radial distribution functions

        If (l_do_rdf) Call rdf_collect(i, rrt, neigh, config, rdf)

      End Do

#ifdef KIM
    End If
#endif
    ! Poisson solver alternative to Ewald

    If (electro%key == ELECTROSTATIC_POISSON) Then
      Call poisson_forces(engacc, viracc, stats%stress, pois, electro, domain, config, comm)
      engcpe_rl = engcpe_rl + engacc
      vircpe_rl = vircpe_rl + viracc
    End If

    ! metal potential safety

    If (safe) Then
      tmp = 0.0_wp
    Else
      tmp = 1.0_wp
    End If

    ! in the case of bonded interactions 3 possible subcases
    ! cases for excluded interactions:
    ! (1) RDF accumulate further the short-range exclusions
    ! (2) Ewald corrections due to short-range exclusions
    ! (3) CHARMM core-shell self-induction additions

    If (lbook .and. &
        (l_do_rdf .or. (Any([ELECTROSTATIC_EWALD, ELECTROSTATIC_POISSON] == electro%key)) &
         .or. mpoles%key == POLARISATION_CHARMM)) Then
      Do i = 1, config%natms ! outer loop over atoms
        limit = neigh%list(-1, i) - neigh%list(0, i) ! Get neigh%list limit
        If (limit > 0) Then

          ! calculate interatomic distances

          Do k = 1, limit
            j = neigh%list(neigh%list(0, i) + k, i)

            xxt(k) = config%parts(i)%xxx - config%parts(j)%xxx
            yyt(k) = config%parts(i)%yyy - config%parts(j)%yyy
            zzt(k) = config%parts(i)%zzz - config%parts(j)%zzz
          End Do

          ! periodic boundary conditions not needed by LC construction
          !
          !           Call images(imcon,cell,limit,xxt,yyt,zzt)

          ! square of distances

          Do k = 1, limit
            rrt(k) = Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
          End Do

          ! accumulate radial distribution functions

          If (l_do_rdf) Call rdf_excl_collect(i, rrt, neigh, config, rdf)

          If (electro%key == ELECTROSTATIC_EWALD) Then ! Ewald corrections

            Call ewald_excl_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, ewld, ewld%spme_data(0), config)
            ! Call ewald_excl_forces(ewld, ewld%spme_data(0), neigh, electro, config, coul_coeffs, &
            !                        i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress)

            engcpe_ex = engcpe_ex + engacc
            vircpe_ex = vircpe_ex + viracc
          End If

          ! get CHARMM core-shell self-induction contributions

          If (mpoles%key == POLARISATION_CHARMM) Then
            If (neigh%list(-3, i) - neigh%list(0, i) > 0) Then
              Call coul_chrm_forces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, config)

              engcpe_ch = engcpe_ch + engacc
              vircpe_ch = vircpe_ch + viracc
            End If
          End If

        End If
      End Do
    End If

    ! counter for rdf%rdf statistics outside loop structures
    ! and frozen-frozen rdf%rdf completeness

    If (l_do_rdf) Then
      If (megfrz /= 0) Then

        ! outer loop over atoms

        Do i = 1, config%natms

          ! Get neigh%list limit

          limit = neigh%list(-2, i) - neigh%list(-1, i)
          If (limit > 0) Then

            ! calculate interatomic distances

            Do k = 1, limit
              j = neigh%list(neigh%list(-1, i) + k, i)

              xxt(k) = config%parts(i)%xxx - config%parts(j)%xxx
              yyt(k) = config%parts(i)%yyy - config%parts(j)%yyy
              zzt(k) = config%parts(i)%zzz - config%parts(j)%zzz
            End Do

            ! periodic boundary conditions not needed by LC construction
            !
            !              Call images(imcon,cell,limit,xxt,yyt,zzt)

            ! square of distances

            Do k = 1, limit
              rrt(k) = Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
            End Do

            ! accumulate radial distribution functions

            Call rdf_frzn_collect(i, rrt, neigh, config, rdf)
          End If

        End Do

      End If

      rdf%n_configs = rdf%n_configs + 1
    End If

#ifdef CHRONO
    Call stop_timer(tmr, 'Short Range')
    Call start_timer(tmr, 'Two-Body Final')
#endif

    Deallocate (xxt, yyt, zzt, rrt, Stat=fail)
    If (fail > 0) Call error_dealloc('distance arrays', 'two_body_forces')

    !Increase rdf%block_number when required
    If (l_do_rdf) Then
      Call rdf_increase_block_number(rdf, nstep)
    End If

    ! Further Ewald/Poisson Solver corrections or an infrequent refresh

    If (ewld%vdw) Then
      engvdw = engvdw_rc + engvdw_rl
      virvdw = virvdw_rc + virvdw_rl
    End If

    If (Any([ELECTROSTATIC_EWALD, ELECTROSTATIC_POISSON] == electro%key)) Then
      If (ELECTROSTATIC_EWALD == electro%key) Then
        Deallocate (coul_coeffs, stat=fail)
        If (fail > 0) Call error_dealloc('coul_coeffs', 'two_body_forces')
        If (ewld%vdw) Then
          Deallocate (vdw_coeffs, stat=fail)
          If (fail > 0) Call error_dealloc('vdw_coeffs', 'two_body_forces')
        End If
      End If

      ! non-zero total system charge correction (for the whole system)
      ! ( Fuchs, Proc. R. Soc., A, 151, (585),1935 )
      If (Abs(config%sumchg) > 1.0e-6_wp) Then
        factor_nz = -0.5_wp * (pi * r4pie0 / electro%eps) * (config%sumchg / ewld%alpha)**2

        engcpe_nz = factor_nz / config%volm
        vircpe_nz = -3.0_wp * engcpe_nz
      End If
    End If

    ! Find the change of energy produced by the torques on multipoles
    ! under infinitesimal rotations & convert to Cartesian coordinates

    If (mpoles%max_mpoles > 0) Then
      Call d_ene_trq_mpoles(vircpe_dt, stats%stress, mpoles, config)
    End If

    ! sum up contributions to domain,potentials

    buffer(0) = tmp
    buffer(1) = engkim
    buffer(2) = virkim
    buffer(3) = engden
    buffer(4) = virden
    buffer(5) = engmet
    buffer(6) = virmet
    buffer(7) = engvdw
    buffer(8) = virvdw
    buffer(9) = engcpe_rc
    buffer(10) = vircpe_rc
    buffer(11) = engcpe_rl
    buffer(12) = vircpe_rl
    buffer(13) = engcpe_ch
    buffer(14) = vircpe_ch
    buffer(15) = engcpe_ex
    buffer(16) = vircpe_ex
    buffer(17) = engcpe_fr
    buffer(18) = vircpe_fr
    buffer(19) = vircpe_dt

    Call gsum(comm, buffer(0:19))

    tmp = buffer(0)
    engkim = buffer(1)
    virkim = buffer(2)
    engden = buffer(3)
    virden = buffer(4)
    engmet = buffer(5)
    virmet = buffer(6)
    engvdw = buffer(7)
    virvdw = buffer(8)
    engcpe_rc = buffer(9)
    vircpe_rc = buffer(10)
    engcpe_rl = buffer(11)
    vircpe_rl = buffer(12)
    engcpe_ch = buffer(13)
    vircpe_ch = buffer(14)
    engcpe_ex = buffer(15)
    vircpe_ex = buffer(16)
    engcpe_fr = buffer(17)
    vircpe_fr = buffer(18)
    vircpe_dt = buffer(19)

    safe = (tmp < 0.5_wp)
    If (.not. safe) Call error(505)

    ! Self-interaction is constant for the default charges only SPME

    ! If (electro%key == ELECTROSTATIC_EWALD) Then ! Sum it up for multipolar SPME
    !    If (mpoles%max_mpoles > 0 .and. mpoles%max_order <= 2) Call gsum(comm,ewld%spme_data(0)%self_interaction)
    !Write(message,'(a,1p,e18.10)') 'Self-interaction term: ',engsic
    !Call info(message,.true.)
    ! End If

    ! Globalise coulombic contributions: cpe

    stats%engcpe = stats%engcpe + engcpe_rc + engcpe_rl + engcpe_ch + engcpe_ex + engcpe_fr + engcpe_nz
    stats%vircpe = stats%vircpe + vircpe_rc + vircpe_rl + vircpe_ch + vircpe_ex + vircpe_fr + vircpe_nz + vircpe_dt

    ! Add non-zero total system charge correction to
    ! diagonal terms of stress tensor (per node)

    tmp = -vircpe_nz / (3.0_wp * Real(comm%mxnode, wp))
    stats%stress(1) = stats%stress(1) + tmp
    stats%stress(5) = stats%stress(5) + tmp
    stats%stress(9) = stats%stress(9) + tmp

    ! Globalise short-range, KIM and metal interactions with
    ! their long-range corrections contributions: srp

    stats%engsrp = stats%engsrp + engkim + (engden + engmet + met%elrc(0)) + (engvdw + vdws%elrc)
    stats%virsrp = stats%virsrp + virkim + (virden + virmet + met%vlrc(0)) + (virvdw + vdws%vlrc)

    If (stats%collect_pp) stats%pp_energy = stats%pp_energy + vdws%elrc / config%megatm

    ! Add long-range corrections to diagonal terms of stress tensor (per node)

    tmp = -(vdws%vlrc + met%vlrc(0)) / (3.0_wp * Real(comm%mxnode, wp))
    stats%stress(1) = stats%stress(1) + tmp
    stats%stress(5) = stats%stress(5) + tmp
    stats%stress(9) = stats%stress(9) + tmp

#ifdef CHRONO
    Call stop_timer(tmr, 'Two-Body Final')
#endif

  End Subroutine two_body_forces
End Module two_body
