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
  Use errors_warnings, Only: error
  Use ewald,           Only: ewald_type
  Use ewald_mpole,     Only: ewald_excl_mforces,&
                             ewald_excl_mforces_d,&
                             ewald_frzn_mforces,&
                             ewald_real_mforces,&
                             ewald_real_mforces_d,&
                             ewald_spme_mforces,&
                             ewald_spme_mforces_d
  Use ewald_spole,     Only: ewald_excl_forces,&
                             ewald_frzn_forces,&
                             ewald_real_forces,&
                             ewald_spme_forces
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
  Use poisson,         Only: poisson_forces,&
                             poisson_frzn_forces,&
                             poisson_type
  Use rdfs,            Only: rdf_collect,&
                             rdf_excl_collect,&
                             rdf_frzn_collect,&
                             rdf_increase_block_number,&
                             rdf_type
  Use site,            Only: site_type
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

    Integer,                  Intent(In   ) :: ensemble
    Logical,                  Intent(In   ) :: lbook
    Integer,                  Intent(In   ) :: megfrz
    Logical,                  Intent(In   ) :: leql
    Integer,                  Intent(In   ) :: nsteql, nstep
    Type(stats_type),         Intent(InOut) :: stats
    Type(ewald_type),         Intent(InOut) :: ewld
    Type(metal_type),         Intent(InOut) :: met
    Type(poisson_type),       Intent(InOut) :: pois
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(site_type),          Intent(In   ) :: sites
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(mpole_type),         Intent(InOut) :: mpoles
    Type(electrostatic_type), Intent(InOut) :: electro
    Type(domains_type),       Intent(In   ) :: domain
    Type(timer_type),         Intent(InOut) :: tmr
    Type(kim_type),           Intent(InOut) :: kim_data
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Integer                                  :: fail, i, j, k, limit
    Logical                                  :: l_do_outer_loop, l_do_rdf, safe
    Real(Kind=wp)                            :: buffer(0:19), engacc, engcpe_ch, engcpe_ex, &
                                                engcpe_fr, engcpe_nz, engcpe_rc, engcpe_rl, &
                                                engden, engkim, engmet, engvdw, factor_nz, tmp, &
                                                viracc, vircpe_ch, vircpe_dt, vircpe_ex, &
                                                vircpe_fr, vircpe_nz, vircpe_rc, vircpe_rl, &
                                                virden, virkim, virmet, virvdw
    Real(Kind=wp), Allocatable, Dimension(:) :: rrt, xxt, yyt, zzt

    safe = .true.
    fail = 0
    Allocate (xxt(1:neigh%max_list), yyt(1:neigh%max_list), zzt(1:neigh%max_list), rrt(1:neigh%max_list), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'two_body_forces allocation failure'
      Call error(0, message)
    End If

    l_do_rdf = (rdf%l_collect .and. ((.not. leql) .or. nstep >= nsteql) .and. Mod(nstep, rdf%freq) == 0)

    ! If k-space SPME is evaluated infrequently check whether
    ! at this timestep to evaluate or "refresh" with old values.
    ! At restart allocate the "refresh" arrays and force a fresh
    ! evaluation.  Repeat the same but only for the SPME k-space
    ! frozen-frozen evaluations in constant volume ensembles only.

    If (Any([ELECTROSTATIC_EWALD, ELECTROSTATIC_POISSON] == electro%key)) Then
      Call ewld%check(ensemble, megfrz, nsteql, electro%nstfce, nstep, config%mxatms)
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

    ! Calculate all contributions from KIM
    If (kim_data%active) Then
      Call kim_energy_and_forces(kim_data, config%natms, config%nlast, config%parts, &
                                 neigh%list, domain%map, config%lsite, config%ltype, config%lsi, config%lsa, config%ltg, &
                                 sites%site_name, engkim, virkim, stats%stress, comm)
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

    If (electro%key == ELECTROSTATIC_EWALD .and. ewld%l_fce) Then
      If (mpoles%max_mpoles > 0) Then
        If (mpoles%max_order <= 2) Then
          Call ewald_spme_mforces_d(engcpe_rc, vircpe_rc, stats%stress, ewld, mpoles, &
                                    electro, domain, config, comm)
        Else
          Call ewald_spme_mforces(engcpe_rc, vircpe_rc, stats%stress, ewld, mpoles, &
                                  electro, domain, config, comm)
        End If
      Else
        Call ewald_spme_forces(engcpe_rc, vircpe_rc, stats%stress, ewld, electro, &
                               domain, config, comm)
      End If
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
          If (vdws%l_direct) Then ! direct calculation

            Call vdw_forces_direct(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, vdws, config)

          else

            Call vdw_forces_tab(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, vdws, config)

          end If


          engvdw = engvdw + engacc
          virvdw = virvdw + viracc
        End If

        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! COULOMBIC CONTRIBUTIONS
        !!!!!!!!!!!!!!!!!!!!!!!!1

        If (mpoles%max_mpoles > 0) Then

          !!! MULTIPOLAR ATOMIC SITES

          If (electro%key == ELECTROSTATIC_EWALD) Then

            ! calculate coulombic forces, Ewald sum - real space contribution

            If (mpoles%max_order <= 2) Then
              Call ewald_real_mforces_d(i, xxt, yyt, zzt, rrt, engacc, &
                                        viracc, stats%stress, ewld, neigh, mpoles, electro, config)
            Else
              Call ewald_real_mforces(i, xxt, yyt, zzt, rrt, engacc, &
                                      viracc, stats%stress, neigh, mpoles, electro, config)
            End If

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_DDDP) Then

            ! distance dependant dielectric potential

            Call coul_dddp_mforces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_COULOMB) Then

            ! coulombic 1/r potential with no truncation or damping

            Call coul_cp_mforces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then

            ! force-shifted coulomb potentials

            Call coul_fscp_mforces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, electro, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then

            ! reaction field potential

            Call coul_rfp_mforces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, mpoles, electro, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          End If

        Else

          If (electro%key == ELECTROSTATIC_EWALD) Then

            ! calculate coulombic forces, Ewald sum - real space contribution

            Call ewald_real_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, electro, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_DDDP) Then

            ! distance dependant dielectric potential

            Call coul_dddp_forces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_COULOMB) Then

            ! coulombic 1/r potential with no truncation or damping

            Call coul_cp_forces(i, electro%eps, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_COULOMB_FORCE_SHIFT) Then

            ! force-shifted coulomb potentials

            Call coul_fscp_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, electro, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

          Else If (electro%key == ELECTROSTATIC_COULOMB_REACTION_FIELD) Then

            ! reaction field potential

            Call coul_rfp_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, electro, config)

            engcpe_rl = engcpe_rl + engacc
            vircpe_rl = vircpe_rl + viracc

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
      Call poisson_forces(engacc, viracc, stats%stress, pois, electro, domain, config, ewld, comm)
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
            If (mpoles%max_mpoles > 0) Then
              If (mpoles%max_order <= 2) Then
                Call ewald_excl_mforces_d(i, xxt, yyt, zzt, rrt, engacc, viracc, &
                                          stats%stress, neigh, mpoles, electro, config)
              Else
                Call ewald_excl_mforces(i, xxt, yyt, zzt, rrt, engacc, viracc, &
                                        stats%stress, neigh, mpoles, electro, config)
              End If
            Else
              Call ewald_excl_forces(i, xxt, yyt, zzt, rrt, engacc, viracc, stats%stress, neigh, electro, config)
            End If

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

    Deallocate (xxt, yyt, zzt, rrt, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'two_body_forces deallocation failure'
      Call error(0, message)
    End If

    !Increase rdf%block_number when required
    If (l_do_rdf) Then
      Call rdf_increase_block_number(rdf, nstep)
    End If

    ! Further Ewald/Poisson Solver corrections or an infrequent refresh

    If (Any([ELECTROSTATIC_EWALD, ELECTROSTATIC_POISSON] == electro%key)) Then
      If (ewld%l_fce) Then

        ! frozen pairs corrections to coulombic forces

        If (megfrz /= 0) Then
          If (electro%key == ELECTROSTATIC_EWALD) Then ! Ewald
            If (mpoles%max_mpoles > 0) Then
              Call ewald_frzn_mforces(engcpe_fr, vircpe_fr, stats%stress, ewld, &
                                      neigh, mpoles, electro, config, comm)
            Else
              Call ewald_frzn_forces(engcpe_fr, vircpe_fr, stats%stress, ewld, neigh, electro, config, comm)
            End If
          Else !If (electro%key == ELECTROSTATIC_POISSON) Then ! Poisson Solver
            Call poisson_frzn_forces(electro%eps, engcpe_fr, vircpe_fr, stats%stress, ewld, neigh, config, comm)
          End If
        End If

      Else

        ! Refresh all Ewald k-space contributions

        Call ewld%refresh(engcpe_rc, vircpe_rc, engcpe_fr, vircpe_fr, stats%stress, config)

      End If

      ! non-zero total system charge correction (for the whole system)
      ! ( Fuchs, Proc. R. Soc., A, 151, (585),1935 )

      If (Abs(config%sumchg) > 1.0e-6_wp) Then
        factor_nz = -0.5_wp * (pi * r4pie0 / electro%eps) * (config%sumchg / electro%alpha)**2

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

    If (electro%key == ELECTROSTATIC_EWALD) Then ! Sum it up for multipolar SPME
      If (mpoles%max_mpoles > 0 .and. mpoles%max_order <= 2) Call gsum(comm, ewld%engsic)
      !Write(message,'(a,1p,e18.10)') 'Self-interaction term: ',engsic
      !Call info(message,.true.)
    End If

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

    ! Add long-range corrections to diagonal terms of stress tensor (per node)

    tmp = -(vdws%vlrc + met%vlrc(0)) / (3.0_wp * Real(comm%mxnode, wp))
    stats%stress(1) = stats%stress(1) + tmp
    stats%stress(5) = stats%stress(5) + tmp
    stats%stress(9) = stats%stress(9) + tmp

#ifdef CHRONO
    Call stop_timer(tmr, 'Short Range')
#endif

  End Subroutine two_body_forces
End Module two_body
