Module ewald_spole
  !!----------------------------------------------------------------------!
  !!
  !! dl_poly_4 module for calculating coulombic energy and force terms
  !! in a periodic system using the smooth particle mesh ewald method
  !! by essmann et al. j. chem. phys. 103 (1995) 8577
  !! https://doi.org/10.1063/1.470117
  !!
  !! It is important to note that the convention of Fourier transform follows
  !! that of the Essmann, Darden, et al. throughout, as does the form
  !! of the reciprocal lattice vectors.
  !!
  !! copyright - daresbury laboratory
  !! authors   - i.t.todorov & w.smith & i.j.bush & j.s.wilkins august 2018
  !!
  !!----------------------------------------------------------------------!
  Use comms,           Only: comms_type
  Use configuration,   Only: configuration_type
  Use constants,       Only: r4pie0,&
                             sqrpi,&
                             twopi,&
                             zero_plus
  Use domains,         Only: domains_type
  Use electrostatic,   Only: electrostatic_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc
  Use ewald,           Only: ewald_type
  Use kinds,           Only: wp
  Use kspace,          Only: kspace_type
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: calc_erf,&
                             calc_erf_deriv,&
                             calc_erfc,&
                             dcell,&
                             invert
  Use parallel_fft,    Only: initialize_fft,&
                             pfft,&
                             pfft_indices
  Use spme,            Only: spme_component
  Use statistics,      Only: calculate_stress,&
                             stats_type
  Use timer,           Only: start_timer,&
                             stop_timer,&
                             timer_type

  Implicit None

  Private

  Public ::  ewald_real_forces_coul_tab, ewald_real_forces_coul, ewald_real_forces_gen
  Public ::  ewald_spme_forces
  Public ::  ewald_excl_forces, ewald_frzn_forces

  !> Temp X,Y,Z Scaled Coords (U/mu)
  Real(kind=wp), Dimension(:, :), Allocatable       :: recip_coords
  !> Indices to avoid type conversion
  Integer, Dimension(:, :), Allocatable       :: recip_indices
  !> temporary workspace for parallel fft
  Complex(Kind=wp), Dimension(:, :, :), Allocatable :: pfft_work

Contains

  Subroutine ewald_real_forces_coul_tab(electro, alpha, spme_datum, neigh, config, stats, iatm, x_pos, y_pos, z_pos, mod_dr_ij, &
    & engcpe_rl, vircpe_rl)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using tabulated erfs for ewald's method
    !!
    !! note: Real space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith august 1998
    !! amended   - i.t.todorov april 2015
    !! amended   - j. wilkins september 2018
    !!
    !!-----------------------------------------------------------------------
    Type(electrostatic_type),                   Intent(In   ) :: electro
    Real(Kind=wp),                              Intent(In   ) :: alpha
    Type(spme_component),                       Intent(In   ) :: spme_datum
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(configuration_type),                   Intent(InOut) :: config
    Type(stats_type),                           Intent(InOut) :: stats
    Integer,                                    Intent(In   ) :: iatm
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: x_pos, y_pos, z_pos, mod_dr_ij
    Real(Kind=wp),                              Intent(  Out) :: engcpe_rl, vircpe_rl

    Integer                     :: global_id_i, global_id_j, jatm, m
    Real(Kind=wp)               :: alpha_r, atom_coeffs_i, e_comp, erf_gamma, &
                                   mod_r_ij, prefac
    Real(Kind=wp), Dimension(9) :: stress_temp, stress_temp_comp
    Real(Kind=wp), Dimension(3) :: force_temp, force_temp_comp, pos_j
    Integer                     :: nearest_sample_index
    Real(Kind=wp)               :: difference
    Real(Kind=wp), Dimension(2) :: temp
    Real(Kind=wp), Dimension(3) :: points

!! Current atom
!! Atoms positions (neighbours, not global) and inter-particle separations
!! Energy and virial for the Real component
!! Type containing stress and per-particle info
!! Atom storage of coeffs !!Dimension(size(coeffs,1))
!! Tempeorary stress tensor
!! Temporary force vectors
!! Position of ion j
!! g_p & d/dr[g_p]
!! Q*g_p
!! Coeffs*inv_mod_r_ij**n
!! Inter-particle distances

    ! initialise accumulators

    engcpe_rl = 0.0_wp
    vircpe_rl = 0.0_wp
    stress_temp = 0.0_wp
    force_temp = 0.0_wp

    ! global identity of iatm

    global_id_i = config%ltg(iatm)

    atom_coeffs_i = config%parts(iatm)%chge * spme_datum%scaling

    ! ignore interaction if the coeffs or scaling are zero
    If (Abs(atom_coeffs_i) < zero_plus) Return

    ! start of primary loop for forces evaluation

    Do m = 1, neigh%list(0, iatm)

      ! atomic index and charge

      jatm = neigh%list(m, iatm)
      global_id_j = config%ltg(jatm)

      ! interatomic distance
      mod_r_ij = mod_dr_ij(m)
      prefac = config%parts(jatm)%chge

      ! interaction validity and truncation of potential
      If (Abs(prefac) > zero_plus .and. mod_r_ij < neigh%cutoff) Then

        pos_j = [x_pos(m), y_pos(m), z_pos(m)]

        ! Complete prefactor
        prefac = atom_coeffs_i * prefac

        nearest_sample_index = Int(mod_r_ij * electro%erfc_deriv%recip_spacing)
        difference = mod_r_ij * electro%erfc_deriv%recip_spacing - Real(nearest_sample_index, wp)
        points = electro%erfc_deriv%table(nearest_sample_index:nearest_sample_index + 2)
        if (nearest_sample_index == 0) points(1) = points(1) * mod_r_ij
        temp(1) = points(1) + (points(2) - points(1)) * difference
        temp(2) = points(2) + (points(3) - points(2)) * (difference - 1.0_wp)
        erf_gamma = prefac * (temp(1) + (temp(2) - temp(1)) * difference * 0.5_wp)
        ! erf_gamma = prefac * electro%erfc_deric%calc(mod_r_ij)

        ! calculate forces ( dU * r/||r|| )

        force_temp_comp = erf_gamma * pos_j
        force_temp = force_temp + force_temp_comp

        If (jatm <= config%natms .or. global_id_i < global_id_j) Then
          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - force_temp_comp(1)
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - force_temp_comp(2)
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - force_temp_comp(3)

          End If

          ! calculate components of G
          nearest_sample_index = Int(mod_r_ij * electro%erfc%recip_spacing)
          difference = mod_r_ij * electro%erfc%recip_spacing - Real(nearest_sample_index, wp)
          points = electro%erfc%table(nearest_sample_index:nearest_sample_index + 2)
          if (nearest_sample_index == 0) points(1) = points(1) * mod_r_ij
          temp(1) = points(1) + (points(2) - points(1)) * difference
          temp(2) = points(2) + (points(3) - points(2)) * (difference - 1.0_wp)
          e_comp = prefac * (temp(1) + (temp(2) - temp(1)) * difference * 0.5_wp)

          !e_comp = prefac * electro%erfc%calc(mod_r_ij)

          ! calculate interaction energy
          engcpe_rl = engcpe_rl + erf_gamma
          ! calculate virial ( F.r )

          vircpe_rl = vircpe_rl - erf_gamma * mod_r_ij**2

          ! calculate stress tensor
          stress_temp(1) = stress_temp(1) + pos_j(1) * force_temp_comp(1)
          stress_temp(2) = stress_temp(2) + pos_j(1) * force_temp_comp(2)
          stress_temp(3) = stress_temp(3) + pos_j(1) * force_temp_comp(3)
          stress_temp(4) = stress_temp(4) + pos_j(2) * force_temp_comp(2)
          stress_temp(5) = stress_temp(5) + pos_j(2) * force_temp_comp(3)
          stress_temp(6) = stress_temp(6) + pos_j(3) * force_temp_comp(3)
          ! stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          ! stress_temp = stress_temp + stress_temp_comp

        End If

        If (stats%collect_pp) Then
          stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + e_comp * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
          If (jatm <= config%natms) Then
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + e_comp * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
          End If
        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = config%parts(iatm)%fxx + force_temp(1)
    config%parts(iatm)%fyy = config%parts(iatm)%fyy + force_temp(2)
    config%parts(iatm)%fzz = config%parts(iatm)%fzz + force_temp(3)

    ! complete stress tensor

    stats%stress(1) = stats%stress(1) + stress_temp(1)
    stats%stress(2) = stats%stress(2) + stress_temp(2)
    stats%stress(3) = stats%stress(3) + stress_temp(3)
    stats%stress(4) = stats%stress(4) + stress_temp(2)
    stats%stress(5) = stats%stress(5) + stress_temp(4)
    stats%stress(6) = stats%stress(6) + stress_temp(5)
    stats%stress(7) = stats%stress(7) + stress_temp(3)
    stats%stress(8) = stats%stress(8) + stress_temp(5)
    stats%stress(9) = stats%stress(9) + stress_temp(6)
    ! stats%stress = stats%stress + stress_temp

  End Subroutine ewald_real_forces_coul_tab


  Subroutine ewald_real_forces_coul(electro, alpha, spme_datum, neigh, config, stats, iatm, x_pos, y_pos, z_pos, mod_dr_ij, &
    & engcpe_rl, vircpe_rl)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using ewald's method
    !!
    !! note: Real space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith august 1998
    !! amended   - i.t.todorov april 2015
    !! amended   - j. wilkins september 2018
    !!
    !!-----------------------------------------------------------------------
    Use constants, Only: rsqrpi
    Type(electrostatic_type),                   Intent(In   ) :: electro
    Real(Kind=wp),                              Intent(In   ) :: alpha
    Type(spme_component),                       Intent(In   ) :: spme_datum
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(configuration_type),                   Intent(InOut) :: config
    Type(stats_type),                           Intent(InOut) :: stats
    Integer,                                    Intent(In   ) :: iatm
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: x_pos, y_pos, z_pos, mod_dr_ij
    Real(Kind=wp),                              Intent(  Out) :: engcpe_rl, vircpe_rl

    Integer                     :: global_id_i, global_id_j, jatm, m
    Real(Kind=wp)               :: alpha_r, atom_coeffs_i, e_comp, erf_gamma, inv_mod_r_ij, &
                                   mod_r_ij, prefac
    Real(Kind=wp), Dimension(9) :: stress_temp, stress_temp_comp
    Real(Kind=wp), Dimension(3) :: force_temp, force_temp_comp, pos_j
    Integer                     :: nearest_sample_index
    Real(Kind=wp)               :: difference
    Real(Kind=wp), Dimension(2) :: temp
    Real(Kind=wp), Dimension(3) :: points

!! Current atom
!! Atoms positions (neighbours, not global) and inter-particle separations
!! Energy and virial for the Real component
!! Type containing stress and per-particle info
!! Atom storage of coeffs !!Dimension(size(coeffs,1))
!! Tempeorary stress tensor
!! Temporary force vectors
!! Position of ion j
!! g_p & d/dr[g_p]
!! Q*g_p
!! Coeffs*inv_mod_r_ij**n
!! Inter-particle distances

    ! initialise accumulators

    engcpe_rl = 0.0_wp
    vircpe_rl = 0.0_wp
    stress_temp = 0.0_wp
    force_temp = 0.0_wp

    ! global identity of iatm

    global_id_i = config%ltg(iatm)

    atom_coeffs_i = config%parts(iatm)%chge * spme_datum%scaling

    ! ignore interaction if the coeffs or scaling are zero
    If (Abs(atom_coeffs_i) < zero_plus) Return

    ! start of primary loop for forces evaluation

    Do m = 1, neigh%list(0, iatm)

      ! atomic index and charge

      jatm = neigh%list(m, iatm)
      global_id_j = config%ltg(jatm)

      ! interatomic distance
      mod_r_ij = mod_dr_ij(m)
      prefac = config%parts(jatm)%chge

      ! interaction validity and truncation of potential
      If (Abs(prefac) > zero_plus .and. mod_r_ij < neigh%cutoff) Then

        pos_j = [x_pos(m), y_pos(m), z_pos(m)]
        alpha_r = mod_r_ij * alpha
        inv_mod_r_ij = 1.0_wp / mod_r_ij

        ! Complete prefactor
        prefac = atom_coeffs_i * prefac * inv_mod_r_ij

        ! calculate components of G
        nearest_sample_index = Int(alpha_r * electro%erfc%recip_spacing)
        difference = alpha_r * electro%erfc%recip_spacing - Real(nearest_sample_index, wp)
        points = electro%erfc%table(nearest_sample_index:nearest_sample_index + 2)
        if (nearest_sample_index == 0) points(1) = points(1) * alpha_r

        temp(1) = points(1) + (points(2) - points(1)) * difference
        temp(2) = points(2) + (points(3) - points(2)) * (difference - 1.0_wp)
        e_comp = prefac * (temp(1) + (temp(2) - temp(1)) * difference * 0.5_wp)
        ! e_comp = prefac  * electro%erfc%calc(alpha_r)

        ! Because function is g_p(ar)/(r^n)
        ! => -n*(g/r^(n + 1) + a(dg/dr))
        ! calculate components of G
        nearest_sample_index = Int(alpha_r * electro%erfc_deriv%recip_spacing)
        difference = alpha_r * electro%erfc_deriv%recip_spacing - Real(nearest_sample_index, wp)
        points = electro%erfc_deriv%table(nearest_sample_index:nearest_sample_index + 2)
        if (nearest_sample_index == 0) points(1) = points(1) * alpha_r
        temp(1) = points(1) + (points(2) - points(1)) * difference
        temp(2) = points(2) + (points(3) - points(2)) * (difference - 1.0_wp)
        erf_gamma = (e_comp * inv_mod_r_ij + prefac * alpha * (temp(1) + (temp(2) - temp(1)) * difference * 0.5_wp))
        !erf_gamma = (e_comp * inv_mod_r_ij + prefac * alpha * electro%erfc_deriv%calc(alpha_r)) !2.0_wp*rsqrpi*alpha*exp(-(alpha_r**2)) )

        ! calculate forces ( dU * r/||r|| )

        force_temp_comp = erf_gamma * pos_j * inv_mod_r_ij
        force_temp = force_temp + force_temp_comp

        If (jatm <= config%natms .or. global_id_i < global_id_j) Then
          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - force_temp_comp(1)
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - force_temp_comp(2)
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - force_temp_comp(3)

          End If

          ! calculate interaction energy
          engcpe_rl = engcpe_rl + e_comp

          ! calculate virial ( F.r )

          vircpe_rl = vircpe_rl - erf_gamma * mod_r_ij

          ! calculate stress tensor
          stress_temp(1) = stress_temp(1) + pos_j(1) * force_temp_comp(1)
          stress_temp(2) = stress_temp(2) + pos_j(1) * force_temp_comp(2)
          stress_temp(3) = stress_temp(3) + pos_j(1) * force_temp_comp(3)
          stress_temp(4) = stress_temp(4) + pos_j(2) * force_temp_comp(2)
          stress_temp(5) = stress_temp(5) + pos_j(2) * force_temp_comp(3)
          stress_temp(6) = stress_temp(6) + pos_j(3) * force_temp_comp(3)
          ! stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          ! stress_temp = stress_temp + stress_temp_comp

        End If

        If (stats%collect_pp) Then
          stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + e_comp * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
          If (jatm <= config%natms) Then
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + e_comp * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
          End If
        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = config%parts(iatm)%fxx + force_temp(1)
    config%parts(iatm)%fyy = config%parts(iatm)%fyy + force_temp(2)
    config%parts(iatm)%fzz = config%parts(iatm)%fzz + force_temp(3)

    ! complete stress tensor
    stats%stress(1) = stats%stress(1) + stress_temp(1)
    stats%stress(2) = stats%stress(2) + stress_temp(2)
    stats%stress(3) = stats%stress(3) + stress_temp(3)
    stats%stress(4) = stats%stress(4) + stress_temp(2)
    stats%stress(5) = stats%stress(5) + stress_temp(4)
    stats%stress(6) = stats%stress(6) + stress_temp(5)
    stats%stress(7) = stats%stress(7) + stress_temp(3)
    stats%stress(8) = stats%stress(8) + stress_temp(5)
    stats%stress(9) = stats%stress(9) + stress_temp(6)
    ! stats%stress = stats%stress + stress_temp

  End Subroutine ewald_real_forces_coul

  Subroutine ewald_real_forces_gen(alpha, spme_datum, neigh, config, stats, coeffs, iatm, x_pos, y_pos, z_pos, mod_dr_ij, &
    & engcpe_rl, vircpe_rl)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using ewald's method
    !!
    !! note: Real space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith august 1998
    !! amended   - i.t.todorov april 2015
    !! amended   - j. wilkins september 2018
    !!
    !!-----------------------------------------------------------------------
    ! use spme, only : g_p, g_p_d

    Real(Kind=wp),                              Intent(In   ) :: alpha
    Type(spme_component),                       Intent(In   ) :: spme_datum
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(configuration_type),                   Intent(InOut) :: config
    Type(stats_type),                           Intent(InOut) :: stats
    Real(Kind=wp), Dimension(:),                Intent(In   ) :: coeffs
    Integer,                                    Intent(In   ) :: iatm
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: x_pos, y_pos, z_pos, mod_dr_ij
    Real(Kind=wp),                              Intent(  Out) :: engcpe_rl, vircpe_rl

    Integer                     :: global_id_i, global_id_j, jatm, m
    Real(Kind=wp)               :: alpha_r, atom_coeffs_i, e_comp, erf_gamma, g_fac, inv_mod_r_ij, &
                                   mod_r_ij, prefac
    Real(Kind=wp), Dimension(9) :: stress_temp, stress_temp_comp
    Real(Kind=wp), Dimension(3) :: force_temp, force_temp_comp, pos_j

!! Current atom
!! Coulomb charges/ multipole coeffs, etc.
!! Atoms positions (neighbours, not global) and inter-particle separations
!! Energy and virial for the Real component
!! Type containing stress and per-particle info
!! Atom storage of coeffs !!Dimension(size(coeffs,1))
!! Tempeorary stress tensor
!! Temporary force vectors
!! Position of ion j
!! g_p & d/dr[g_p]
!! Q*g_p
!! Coeffs*inv_mod_r_ij**n
!! Inter-particle distances

    ! initialise accumulators

    engcpe_rl = 0.0_wp
    vircpe_rl = 0.0_wp
    stress_temp = 0.0_wp
    force_temp = 0.0_wp

    ! global identity of iatm

    global_id_i = config%ltg(iatm)

    atom_coeffs_i = coeffs(iatm) * spme_datum%scaling

    ! ignore interaction if the coeffs or scaling are zero
    If (Abs(atom_coeffs_i) < zero_plus) Return

    ! start of primary loop for forces evaluation

    Do m = 1, neigh%list(0, iatm)

      ! atomic index and charge

      jatm = neigh%list(m, iatm)
      global_id_j = config%ltg(jatm)

      ! interatomic distance
      mod_r_ij = mod_dr_ij(m)
      prefac = coeffs(jatm)

      ! interaction validity and truncation of potential
      If (Abs(prefac) > zero_plus .and. mod_r_ij < neigh%cutoff) Then

        pos_j = [x_pos(m), y_pos(m), z_pos(m)]
        alpha_r = mod_r_ij * alpha
        inv_mod_r_ij = 1.0_wp / mod_r_ij

        ! Complete prefactor
        prefac = atom_coeffs_i * prefac * inv_mod_r_ij**spme_datum%pot_order

        ! calculate components of G
        g_fac = g_p(alpha_r, spme_datum%pot_order)

        e_comp = prefac * g_fac

        ! Because function is g_p(ar)/(r^n)
        ! => -n*(g/r^(n + 1) + a(dg/dr))
        erf_gamma = prefac * (g_p_d(alpha_r, spme_datum%pot_order) * alpha + &
          & spme_datum%pot_order * g_fac * inv_mod_r_ij)

        ! calculate forces ( dU * r/||r|| )

        force_temp_comp = erf_gamma * pos_j * inv_mod_r_ij
        force_temp = force_temp + force_temp_comp

        If (jatm <= config%natms .or. global_id_i < config%ltg(jatm)) Then
          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - force_temp_comp(1)
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - force_temp_comp(2)
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - force_temp_comp(3)

          End If

          ! calculate interaction energy
          engcpe_rl = engcpe_rl + e_comp

          ! calculate virial ( F.r )

          vircpe_rl = vircpe_rl - erf_gamma * mod_r_ij

          ! calculate stress tensor
          stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          stress_temp = stress_temp + stress_temp_comp

        End If

        If (stats%collect_pp) Then
          stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + e_comp * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
          If (jatm <= config%natms) Then
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + e_comp * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
          End If
        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = config%parts(iatm)%fxx + force_temp(1)
    config%parts(iatm)%fyy = config%parts(iatm)%fyy + force_temp(2)
    config%parts(iatm)%fzz = config%parts(iatm)%fzz + force_temp(3)

    ! complete stress tensor

    stats%stress = stats%stress + stress_temp

  End Subroutine ewald_real_forces_gen

  Subroutine ewald_spme_forces(ewld, spme_datum, electro, domain, config, comm, coeffs, stats, &
    & engcpe_rc, vircpe_rc, tmr)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using the smooth particle mesh ewald method
    !! by essmann et al. j. chem. phys. 103 (1995) 8577
    !!
    !! note: (fourier) reciprocal space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - i.t.todorov & w.smith & i.j.bush february 2016
    !! re-written in per-particle formulation - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Use comms, Only: gsum, gcheck, gsync
    Use constants, Only: twopi, pi, sqrpi, zero_plus
    Use bspline, Only: bspline_splines_gen
    Type(ewald_type),            Intent(inout) :: ewld
    Type(spme_component),        Intent(inout) :: spme_datum
    Type(electrostatic_type),    Intent(In   ) :: electro
    Type(domains_type),          Intent(In   ) :: domain
    Type(configuration_type),    Intent(inout) :: config
    Type(comms_type),            Intent(inout) :: comm
    Real(kind=wp), Dimension(:), Intent(In   ) :: coeffs
    Type(stats_type),            Intent(inout) :: stats
    Real(kind=wp),               Intent(  Out) :: engcpe_rc, vircpe_rc
    Type(timer_type),            Intent(InOut) :: tmr

    Complex(kind=wp), Allocatable, Dimension(:, :, :), Save :: potential_grid, stress_grid
    Integer                                                 :: dim, i
    Integer, Allocatable, Dimension(:)                      :: to_calc
    Integer, Dimension(4)                                   :: fail
    Logical                                                 :: llspl
    Logical, Save                                           :: newjob = .true.
    Real(kind=wp)                                           :: det, eng, rvolm, scale
    Real(kind=wp), Allocatable, Dimension(:)                :: Q_abc
    Real(kind=wp), Allocatable, Dimension(:, :)             :: F_abc, S_abc
    Real(kind=wp), Allocatable, Dimension(:, :, :), Save    :: charge_grid
    Real(kind=wp), Dimension(9)                             :: rcell, stress_temp

! Inputs and Outputs
!! Coefficients such as charges or pot params
!! Type containing details of stress and per-particle
!! Number of steps taken since calculation start
!! Energy and virial of Coulomb interaction
!! List of points to calculate
!! Unknown? Does this want to be saved?
!! Loop counters
!! Reciprocal lattice vectors
! Data constants and intermediate variables
!! Temporary stress tensor
!! Determinant of inverse matrix
!! Reciprocal volume
!! Coulomb factor / epsq?
!! Energy contribution
!! Energies
!! Forces
!! Stress
!! Ierr

    Call start_timer(tmr, 'Setup')

    If (Any(Abs(coeffs) > zero_plus)) Then
      Continue
    Else
      Return
    End If

    If (newjob) Then
      Call ewald_spme_init(domain, config%mxatms, comm, ewld%kspace, &
        & ewld%bspline, charge_grid, potential_grid, stress_grid)
      newjob = .false.
    End If

    fail = 0
    Allocate (recip_coords(3, config%nlast), stat=fail(1))
    Allocate (recip_indices(3, config%nlast), stat=fail(2))
    Allocate (to_calc(0:config%nlast), stat=fail(3))

    ! If not per-particle only need global sum, else need everything
    If (.not. stats%collect_pp) Then
      Allocate (Q_abc(0:0), F_abc(3, config%natms), S_abc(9, 0:0), stat=fail(4))
    Else
      Allocate (Q_abc(0:config%natms), F_abc(3, config%natms), S_abc(9, 0:config%natms), stat=fail(4))
    End If
    If (Any(fail > 0)) Call error_alloc('output_arrays', 'ewald_spme_forces')

    Q_abc = 0.0_wp
    F_abc = 0.0_wp
    S_abc = 0.0_wp

    ! Initialise accumulator

    to_calc(0) = 0

    ! initialise SPME potential energy and virial

    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    ! set working parameters

    rvolm = 1.0_wp / config%volm

    ! set scaling constant

    scale = pi * sqrpi * ewld%alpha**(spme_datum%pot_order - 3) * (0.5_wp * rvolm) * spme_datum%scaling

    ! calculate reciprocal cell

    Call invert(config%cell, rcell, det)
    If (Abs(det) < 1.0e-6_wp) Call error(120)
    Call stop_timer(tmr, 'Setup')

    ! convert cell coordinates to fractional coordinates intervalled [0,1)
    ! (bottom left corner of md cell) and stretch over kmaxs in different
    ! directions.  only the halo (natms,nlast] has fractional coordinates
    ! outside the [0,1) interval.  in the worst case scenario of one
    ! "effective" link-cell per domain and one domain in the md cell only,
    ! the halo will have fractional coordinates intervalled as
    ! [n,0)u[1,2), where -1 <= n < 0.  only the positive halo is needed by
    ! the b-splines since they distribute/spread charge density in
    ! negative direction with length the length of the spline.
    !
    ! the story has become more complicated with cutoff padding and the
    ! conditional updates of the vnl and thus the halo as now a domain
    ! (1:natms) particle can enter the halo and vice versa.  so dd
    ! bounding is unsafe!!!

    Call start_timer(tmr, 'Recip')
    llspl = .true.
    Do i = 1, config%nlast
      Do dim = 1, 3
        recip_coords(dim, i) = ewld%kspace%k_vec_dim_real(dim) * ( &
          & rcell(dim) * config%parts(i)%xxx + &
          & rcell(dim + 3) * config%parts(i)%yyy + &
          & rcell(dim + 6) * config%parts(i)%zzz + 0.5_wp)
      End Do

      ! if not dd bound in kmax grid space when .not.llvnl = (ewld%bspline%num_spline_pad == ewld%bspline)

      If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines .and. i <= config%natms) &
        & llspl = llspl .and. ( &
        & All(recip_coords(:, i) > ewld%kspace%domain_bounds(:, 1)) .and. &
        & All(recip_coords(:, i) < ewld%kspace%domain_bounds(:, 2)))

      ! detect if a particle is charged and in the md cell or in its positive halo
      ! (coords(i) >= -zero_plus) as the b-splines are negative directionally by propagation
      If (All(recip_coords(:, i) >= -zero_plus) .and. Abs(coeffs(i)) > zero_plus) Then
        to_calc(0) = to_calc(0) + 1
        to_calc(to_calc(0)) = i
      End If
    End Do

    recip_indices = Int(recip_coords)

    Call stop_timer(tmr, 'Recip')

    Call start_timer(tmr, 'BSpline')
    ! check for breakage of llspl when .not.llvnl = (ewld%bspline%num_spline_pad == ewld%bspline)

    ewld%bspline%num_spline_padded = ewld%bspline%num_spline_pad

    If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines) Then
      Call gcheck(comm, llspl)
      If (.not. llspl) ewld%bspline%num_spline_padded = ewld%bspline%num_splines + 1
    End If

    ! construct b-splines for atoms

    Call bspline_splines_gen(config%nlast, recip_coords, ewld%bspline)

    Deallocate (recip_coords, stat=fail(1))
    If (fail(1) > 0) Call error_dealloc('recip_coords', 'ewald_spme_forces')
    Call stop_timer(tmr, 'BSpline')

    Call start_timer(tmr, 'Charge')
    Call spme_construct_charge_array(to_calc(0), ewld, to_calc(1:), recip_indices, electro, coeffs, charge_grid)
    Call stop_timer(tmr, 'Charge')

    If (.not. stats%collect_pp .or. spme_datum%pot_order /= 1) Then

      ! If we don't need per-particle data, we can use the old method of getting the stress (cheaper)
      Call start_timer(tmr, 'Potential')
      Call spme_construct_potential_grid(ewld, rcell, charge_grid, spme_datum, &
        & potential_kernel, potential_grid, s_abc(:, 0))
      Call stop_timer(tmr, 'Potential')
      Call start_timer(tmr, 'ForceEnergy')
      Call spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, potential_grid, stats%collect_pp, q_abc, f_abc)
      Call stop_timer(tmr, 'ForceEnergy')

    Else

      Call spme_construct_potential_grid(ewld, rcell, charge_grid, spme_datum, &
        & potential_kernel, potential_grid)
      Call spme_construct_potential_grid(ewld, rcell, charge_grid, spme_datum, &
        & stress_kernel, stress_grid)

      Call spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, potential_grid, stats%collect_pp, q_abc, f_abc)
      Call spme_calc_stress(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, stress_grid, s_abc)

    End If

    Call start_timer(tmr, 'Output')

    ! Rescale to real space
    q_abc = q_abc * scale
    f_abc = f_abc * scale * 2.0_wp
    s_abc = s_abc * scale

    Q_abc(0) = Q_abc(0) / Real(comm%mxnode)
    S_abc(:, 0) = S_abc(:, 0) / Real(comm%mxnode)

    eng = Q_abc(0)

    ! as only looped over local stuff, we need to gsum the eng

    Call gsum(comm, eng)

    engcpe_rc = eng + spme_datum%self_interaction

    ! add up forces

    config%parts(1:config%natms)%fxx = config%parts(1:config%natms)%fxx + f_abc(1, :)
    config%parts(1:config%natms)%fyy = config%parts(1:config%natms)%fyy + f_abc(2, :)
    config%parts(1:config%natms)%fzz = config%parts(1:config%natms)%fzz + f_abc(3, :)

    ! Put accumulated stress tensor into stress temp to save on cache problems and translate to linear regime

    stress_temp = s_abc(:, 0)

    ! as only looped over local stuff, we need to gsum stress

    Call gsum(comm, stress_temp)

    ! scale strs and distribute per node

    stress_temp(1:9:4) = stress_temp(1:9:4) + eng
    stats%stress = stats%stress + stress_temp
    vircpe_rc = -Sum(stress_temp(1:9:4))

    If (stats%collect_pp) Then
      stats%pp_energy = stats%pp_energy + Q_abc(1:) + comm%mxnode * spme_datum%self_interaction / config%megatm
      stats%pp_stress = stats%pp_stress + S_abc(:, 1:)
    End If

    Deallocate (recip_indices, stat=fail(1))
    ! deallocate (ewld%bspline%derivs, stat=fail(2))
    Deallocate (to_calc, stat=fail(3))
    Deallocate (Q_abc, F_abc, S_abc, stat=fail(4))
    If (Any(fail > 0)) Call error_dealloc('output_arrays', 'ewald_spme_forces')
    Call stop_timer(tmr, 'Output')

  End Subroutine ewald_spme_forces

  Subroutine ewald_excl_forces(iatm,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex,stress, &
      neigh,ewld,spme_datum,config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating coulombic energy and force terms
    ! in a periodic system using ewald's method
    !
    ! Note: exclusion correction terms
    !       frozen pairs are ignored by default, they are not dealt with here
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

    Integer,                                  Intent( In    ) :: iatm
    Type( neighbours_type ), Intent( In    ) :: neigh
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( spme_component ),                   Intent( In    ) :: spme_datum
    Type( ewald_type ),                       Intent( In    ) :: ewld
    Type( configuration_type ),               Intent( InOut ) :: config

    Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp =  0.3275911_wp
    Real( Kind = wp ), Parameter :: rr3  = 1.0_wp/3.0_wp
    Real( Kind = wp ), Parameter :: r10  = 0.1_wp
    Real( Kind = wp ), Parameter :: r42  = 1.0_wp/42.0_wp
    Real( Kind = wp ), Parameter :: r216 = 1.0_wp/216.0_wp

    Integer           :: limit,idi,jatm,m
    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,alpr,alpr2, &
      erfr,egamma,exp1,tt,             &
      fix,fiy,fiz,fx,fy,fz,            &
      strs1,strs2,strs3,strs5,strs6,strs9

    ! initialise potential energy and virial

    engcpe_ex=0.0_wp
    vircpe_ex=0.0_wp

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! global identity of iatm

    idi=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea*spme_datum%scaling

      ! load forces

      fix=config%parts(iatm)%fxx
      fiy=config%parts(iatm)%fyy
      fiz=config%parts(iatm)%fzz

      ! Get neigh%list limit

      limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

      ! start of primary loop for forces evaluation

      Do m=1,limit

        ! atomic index and charge

        jatm=neigh%list(neigh%list(0,iatm)+m,iatm)
        chgprd=config%parts(jatm)%chge

        ! interatomic distance

        rrr=rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd=chgprd*chgea

          ! Squared distance

          rsq=rrr**2

          ! calculate forces

          alpr =rrr*ewld%alpha
          alpr2=alpr*alpr

          ! calculate error function and derivative

          If (alpr < 1.0e-2_wp) Then

            ! close particles (core-shell units) - small distances limit

            erfr=2.0_wp*chgprd*(ewld%alpha/sqrpi) * &
              (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

            egamma=-4.0_wp*chgprd*(ewld%alpha**3/sqrpi) * &
              (rr3+alpr2*(-2.0_wp*r10+alpr2*(3.0_wp*r42-4.0_wp*alpr2*r216)))

          Else

            ! distant particles - traditional

            exp1=Exp(-(ewld%alpha*rrr)**2)
            tt  =1.0_wp/(1.0_wp+pp*ewld%alpha*rrr)

            erfr=chgprd * &
              (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

            egamma=-(erfr-2.0_wp*chgprd*(ewld%alpha/sqrpi)*exp1)/rsq

          End If

          ! calculate forces

          fx = egamma*xxt(m)
          fy = egamma*yyt(m)
          fz = egamma*zzt(m)

          fix=fix+fx
          fiy=fiy+fy
          fiz=fiz+fz

          If (jatm <= config%natms) Then

            config%parts(jatm)%fxx=config%parts(jatm)%fxx-fx
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-fy
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-fz

          End If

          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then

            ! add potential energy and virial

            engcpe_ex = engcpe_ex - erfr
            vircpe_ex = vircpe_ex - egamma*rsq

            ! add stress tensor

            strs1 = strs1 + xxt(m)*fx
            strs2 = strs2 + xxt(m)*fy
            strs3 = strs3 + xxt(m)*fz
            strs5 = strs5 + yyt(m)*fy
            strs6 = strs6 + yyt(m)*fz
            strs9 = strs9 + zzt(m)*fz

          End If

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx=fix
      config%parts(iatm)%fyy=fiy
      config%parts(iatm)%fzz=fiz

      ! complete stress tensor

      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9

    End If

  End Subroutine ewald_excl_forces

  ! Subroutine ewald_excl_forces(ewld, spme_datum, neigh, electro, config, coeffs, iatm, x_pos, y_pos, z_pos, dr_j, &
  !   & engcpe_ex, vircpe_ex, stress)

  !   !!-----------------------------------------------------------------------
  !   !!
  !   !! dl_poly_4 subroutine for calculating coulombic energy and force terms
  !   !! in a periodic system using ewald's method
  !   !!
  !   !! Note: exclusion correction terms
  !   !!       frozen pairs are ignored by default, they are not dealt with here
  !   !!
  !   !! copyright - daresbury laboratory
  !   !! author    - i.t.todorov february 2015
  !   !!
  !   !!-----------------------------------------------------------------------

  !   Type(ewald_type),                           Intent(In   ) :: ewld
  !   Type(spme_component),                       Intent(In   ) :: spme_datum
  !   Type(neighbours_type),                      Intent(In   ) :: neigh
  !   Type(electrostatic_type),                   Intent(In   ) :: electro
  !   Type(configuration_type),                   Intent(InOut) :: config
  !   Real(Kind=wp), Dimension(:),                Intent(In   ) :: coeffs
  !   Integer,                                    Intent(In   ) :: iatm
  !   Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: x_pos, y_pos, z_pos, dr_j
  !   Real(Kind=wp),                              Intent(  Out) :: engcpe_ex, vircpe_ex
  !   Real(Kind=wp), Dimension(1:9),              Intent(InOut) :: stress

  !   Real(Kind=wp), Parameter :: r10 = 0.1_wp, r216 = 1.0_wp / 216.0_wp, r42 = 1.0_wp / 42.0_wp, &
  !        rr3 = 1.0_wp / 3.0_wp! , &
  !        ! a1 = 0.254829592_wp, a2 = -0.284496736_wp, a3 = 1.421413741_wp, &
  !        ! a4 = -1.453152027_wp, a5 = 1.061405429_wp, pp = 0.3275911_wp


  !   Integer                     :: global_id_i, jatm, limit, m
  !   Real(Kind=wp)               :: atom_coeffs_i, atom_coeffs_ij, atom_coeffs_j, dr_alpha, &
  !                                  dr_alpha_2, erf_gamma, erfr, inv_mod_r_ij, mod_r_ij, mod_r_ij_2
  !   Real(Kind=wp), Dimension(9) :: stress_temp
  !   Real(Kind=wp), Dimension(3) :: force_temp, force_temp_comp, pos_j

  !   ! initialise potential energy and virial

  !   engcpe_ex = 0.0_wp
  !   vircpe_ex = 0.0_wp

  !   ! initialise stress tensor accumulators

  !   stress_temp = 0.0_wp

  !   ! global identity of iatm

  !   global_id_i = config%ltg(iatm)

  !   ! ignore interaction if the charge is zero

  !   atom_coeffs_i = coeffs(iatm)
  !   If (Abs(atom_coeffs_i) < zero_plus) Return

  !   atom_coeffs_i = atom_coeffs_i * spme_datum%scaling

  !   ! load forces

  !   force_temp(1) = config%parts(iatm)%fxx
  !   force_temp(2) = config%parts(iatm)%fyy
  !   force_temp(3) = config%parts(iatm)%fzz

  !   ! Get neigh%list limit

  !   limit = neigh%list(-1, iatm) - neigh%list(0, iatm)

  !   ! start of primary loop for forces evaluation

  !   Do m = 1, limit

  !     ! atomic index and charge

  !     jatm = neigh%list(neigh%list(0, iatm) + m, iatm)
  !     atom_coeffs_j = coeffs(jatm)

  !     ! interatomic distance

  !     mod_r_ij = dr_j(m)

  !     ! interaction validity and truncation of potential

  !     If (Abs(atom_coeffs_j) > zero_plus .and. mod_r_ij < neigh%cutoff) Then

  !       ! charge product
  !       atom_coeffs_ij = atom_coeffs_j * atom_coeffs_i
  !       ! Squared distance

  !       ! interatomic distance
  !       pos_j = [x_pos(m), y_pos(m), z_pos(m)]

  !       mod_r_ij = dr_j(m)
  !       mod_r_ij_2 = mod_r_ij**2

  !       ! calculate forces

  !       dr_alpha = mod_r_ij * ewld%alpha

  !       ! calculate error function and derivative

  !       If (dr_alpha < 1.0e-2_wp) Then
  !         ! close particles (core-shell units) - small distances limit

  !         dr_alpha_2 = dr_alpha * dr_alpha

  !         erfr = 2.0_wp * atom_coeffs_ij * (ewld%alpha / sqrpi) * &
  !                (1.0_wp + dr_alpha_2 * (-rr3 + dr_alpha_2 * (r10 + dr_alpha_2 * (-r42 + dr_alpha_2 * r216))))

  !         erf_gamma = -4.0_wp * atom_coeffs_ij * (ewld%alpha**3 / sqrpi) * &
  !                     (rr3 + dr_alpha_2 * (-2.0_wp * r10 + dr_alpha_2 * (3.0_wp * r42 - 4.0_wp * dr_alpha_2 * r216)))

  !       Else

  !         ! distant particles - traditional
  !         inv_mod_r_ij = 1.0_wp / mod_r_ij
  !         mod_r_ij_2 = mod_r_ij**2

  !         inv_mod_r_ij = 1.0_wp / mod_r_ij
  !         erfr = atom_coeffs_ij * calc_erf(dr_alpha) * inv_mod_r_ij
  !         erf_gamma = (atom_coeffs_ij * ewld%alpha * calc_erf_deriv(dr_alpha) - erfr) * inv_mod_r_ij**2

  !       End If

  !       ! calculate forces

  !       force_temp_comp = erf_gamma * pos_j

  !       force_temp = force_temp + force_temp_comp

  !       If (jatm <= config%natms) Then

  !         config%parts(jatm)%fxx = config%parts(jatm)%fxx - force_temp_comp(1)
  !         config%parts(jatm)%fyy = config%parts(jatm)%fyy - force_temp_comp(2)
  !         config%parts(jatm)%fzz = config%parts(jatm)%fzz - force_temp_comp(3)

  !       End If

  !       If (jatm <= config%natms .or. global_id_i < config%ltg(jatm)) Then

  !         ! add potential energy and virial

  !         engcpe_ex = engcpe_ex - erfr
  !         vircpe_ex = vircpe_ex - erf_gamma * mod_r_ij_2

  !         ! calculate stress tensor
  !         stress_temp(1:9:3) = stress_temp(1:9:3) + pos_j * force_temp(1)
  !         stress_temp(2:9:3) = stress_temp(2:9:3) + pos_j * force_temp(2)
  !         stress_temp(3:9:3) = stress_temp(3:9:3) + pos_j * force_temp(3)

  !       End If

  !     End If

  !   End Do

  !   ! load back forces

  !   config%parts(iatm)%fxx = force_temp(1)
  !   config%parts(iatm)%fyy = force_temp(2)
  !   config%parts(iatm)%fzz = force_temp(3)

  !   ! complete stress tensor

  !   stress = stress + stress_temp

  ! End Subroutine ewald_excl_forces

  Subroutine ewald_frzn_forces(engcpe_fr, vircpe_fr, stress, ewld, neigh, electro, config, comm)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for calculating corrections to coulombic energy
    !! and forces in a periodic system arising from frozen pairs
    !!
    !! Note: Forces (as well as velocities) on frozen atoms are zeroed at the
    !!       end (and any COM drift removed) but corrections to the stress
    !!       and the virial are important as they feed into the system
    !!       pressure response.  Constant volume ensembles (ensemble < 20)
    !!       need this calculation just once! - controlled by ewld%lf_fce in
    !!       ewald_check<-two_body_forces
    !!
    !! copyright - daresbury laboratory
    !! author    - i.t.todorov december 2015
    !!
    !!-----------------------------------------------------------------------
    Use comms, Only: gsum
    Real(Kind=wp),                 Intent(  Out) :: engcpe_fr, vircpe_fr
    Real(Kind=wp), Dimension(1:9), Intent(InOut) :: stress
    Type(ewald_type),              Intent(InOut) :: ewld
    Type(neighbours_type),         Intent(In   ) :: neigh
    Type(electrostatic_type),      Intent(In   ) :: electro
    Type(configuration_type),      Intent(InOut) :: config
    Type(comms_type),              Intent(InOut) :: comm

    Real(Kind=wp), Parameter :: a1 = 0.254829592_wp, a2 = -0.284496736_wp, a3 = 1.421413741_wp, &
         a4 = -1.453152027_wp, a5 = 1.061405429_wp, pp = 0.3275911_wp

    Character(Len=256)                       :: message
    Integer                                  :: fail, global_id_i, i, ii, j, jj, k, limit, nzfr
    Integer, Allocatable, Dimension(:)       :: l_ind, nz_fr
    Real(Kind=wp)                            :: atom_coeffs_ij, det, erf_gamma, erfr, exp1, &
         mod_r_ij, mod_r_ij_2, rcell(1:9), scl, tt, xrr, &
         xss, yrr, yss, zrr, zss
    Real(Kind=wp), Allocatable, Dimension(:) :: cfr, dr_j, x_pos, xfr, y_pos, yfr, z_pos, zfr
    Real(Kind=wp), Dimension(3)              :: force_temp_comp
    Real(Kind=wp), Dimension(9)              :: stress_temp

    fail = 0
    Allocate (l_ind(1:config%mxatdm), nz_fr(0:comm%mxnode), Stat=fail)
    If (fail > 0) Call error_alloc('l_ind and nz_fr', 'ewald_frzn_forces')

    Call invert(config%cell, rcell, det)

    ! Initialise contributions

    engcpe_fr = 0.0_wp
    vircpe_fr = 0.0_wp

    stress_temp(1) = 0.0_wp
    stress_temp(2) = 0.0_wp
    stress_temp(3) = 0.0_wp
    stress_temp(5) = 0.0_wp
    stress_temp(6) = 0.0_wp
    stress_temp(9) = 0.0_wp

    l_ind = 0; nz_fr = 0
    Do i = 1, config%natms
       If (config%lfrzn(i) > 0 .and. Abs(config%parts(i)%chge) > zero_plus) Then
          nz_fr(comm%idnode + 1) = nz_fr(comm%idnode + 1) + 1
          l_ind(nz_fr(comm%idnode + 1)) = i
       End If
    End Do
    Call gsum(comm, nz_fr)
    nz_fr(0) = Sum(nz_fr(0:comm%idnode)) ! Offset

    scl = r4pie0 / electro%eps
    nzfr = Sum(nz_fr(1:comm%mxnode)) ! Total
    If (nzfr <= 10 * config%mxatms) Then

       Allocate (cfr(1:nzfr), xfr(1:nzfr), yfr(1:nzfr), zfr(1:nzfr), Stat=fail)
       If (fail > 0) Then
          Write (message, '(a,i0)') 'ewald_frzn_forces allocation failure 1'
          Call error(0, message)
       End If

       cfr = 0.0_wp
       xfr = 0.0_wp
       yfr = 0.0_wp
       zfr = 0.0_wp
       Do i = 1, nz_fr(comm%idnode + 1)
          ii = nz_fr(0) + i

          cfr(ii) = config%parts(l_ind(i))%chge
          xfr(ii) = config%parts(l_ind(i))%xxx
          yfr(ii) = config%parts(l_ind(i))%yyy
          zfr(ii) = config%parts(l_ind(i))%zzz
       End Do
       Call gsum(comm, cfr)
       Call gsum(comm, xfr)
       Call gsum(comm, yfr)
       Call gsum(comm, zfr)

       Do i = 1, nz_fr(comm%idnode + 1)
          ii = nz_fr(0) + i

          Do jj = 1, nz_fr(0) ! -, on nodes<comm%idnode
             xrr = xfr(ii) - xfr(jj)
             yrr = yfr(ii) - yfr(jj)
             zrr = zfr(ii) - zfr(jj)

             xss = (rcell(1) * xrr + rcell(4) * yrr + rcell(7) * zrr)
             yss = (rcell(2) * xrr + rcell(5) * yrr + rcell(8) * zrr)
             zss = (rcell(3) * xrr + rcell(6) * yrr + rcell(9) * zrr)

             xss = xss - Anint(xss)
             yss = yss - Anint(yss)
             zss = zss - Anint(zss)

             xrr = (config%cell(1) * xss + config%cell(4) * yss + config%cell(7) * zss)
             yrr = (config%cell(2) * xss + config%cell(5) * yss + config%cell(8) * zss)
             zrr = (config%cell(3) * xss + config%cell(6) * yss + config%cell(9) * zss)

             ! calculate interatomic distance

             mod_r_ij_2 = xrr**2 + yrr**2 + zrr**2

             mod_r_ij = Sqrt(mod_r_ij_2)
             atom_coeffs_ij = cfr(ii) * cfr(jj) * scl

             ! calculate error function and derivative

             exp1 = Exp(-(ewld%alpha * mod_r_ij)**2)
             tt = 1.0_wp / (1.0_wp + pp * ewld%alpha * mod_r_ij)

             erfr = atom_coeffs_ij * &
                  (1.0_wp - tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1) / mod_r_ij

             erf_gamma = -(erfr - 2.0_wp * atom_coeffs_ij * (ewld%alpha / sqrpi) * exp1) / mod_r_ij_2

             force_temp_comp(1) = erf_gamma * xrr
             force_temp_comp(2) = erf_gamma * yrr
             force_temp_comp(3) = erf_gamma * zrr

             ! calculate forces

             config%parts(l_ind(i))%fxx = config%parts(l_ind(i))%fxx - force_temp_comp(1)
             config%parts(l_ind(i))%fyy = config%parts(l_ind(i))%fyy - force_temp_comp(2)
             config%parts(l_ind(i))%fzz = config%parts(l_ind(i))%fzz - force_temp_comp(3)

          End Do

          Do j = i + 1, nz_fr(comm%idnode + 1) ! =, node=comm%idnode (OVERLAP but no SELF)!
             jj = nz_fr(0) + j

             xrr = xfr(ii) - xfr(jj)
             yrr = yfr(ii) - yfr(jj)
             zrr = zfr(ii) - zfr(jj)

             xss = (rcell(1) * xrr + rcell(4) * yrr + rcell(7) * zrr)
             yss = (rcell(2) * xrr + rcell(5) * yrr + rcell(8) * zrr)
             zss = (rcell(3) * xrr + rcell(6) * yrr + rcell(9) * zrr)

             xss = xss - Anint(xss)
             yss = yss - Anint(yss)
             zss = zss - Anint(zss)

             xrr = (config%cell(1) * xss + config%cell(4) * yss + config%cell(7) * zss)
             yrr = (config%cell(2) * xss + config%cell(5) * yss + config%cell(8) * zss)
             zrr = (config%cell(3) * xss + config%cell(6) * yss + config%cell(9) * zss)

             ! calculate interatomic distance

             mod_r_ij_2 = xrr**2 + yrr**2 + zrr**2

             mod_r_ij = Sqrt(mod_r_ij_2)
             atom_coeffs_ij = cfr(ii) * cfr(jj) * scl

             ! calculate error function and derivative

             exp1 = Exp(-(ewld%alpha * mod_r_ij)**2)
             tt = 1.0_wp / (1.0_wp + pp * ewld%alpha * mod_r_ij)

             erfr = atom_coeffs_ij * &
                  (1.0_wp - tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1) / mod_r_ij

             erf_gamma = -(erfr - 2.0_wp * atom_coeffs_ij * (ewld%alpha / sqrpi) * exp1) / mod_r_ij_2

             force_temp_comp(1) = erf_gamma * xrr
             force_temp_comp(2) = erf_gamma * yrr
             force_temp_comp(3) = erf_gamma * zrr

             ! calculate forces

             config%parts(l_ind(i))%fxx = config%parts(l_ind(i))%fxx - force_temp_comp(1)
             config%parts(l_ind(i))%fyy = config%parts(l_ind(i))%fyy - force_temp_comp(2)
             config%parts(l_ind(i))%fzz = config%parts(l_ind(i))%fzz - force_temp_comp(3)

             config%parts(l_ind(j))%fxx = config%parts(l_ind(j))%fxx + force_temp_comp(1)
             config%parts(l_ind(j))%fyy = config%parts(l_ind(j))%fyy + force_temp_comp(2)
             config%parts(l_ind(j))%fzz = config%parts(l_ind(j))%fzz + force_temp_comp(3)

             ! calculate potential energy and virial

             engcpe_fr = engcpe_fr - erfr
             vircpe_fr = vircpe_fr - erf_gamma * mod_r_ij_2

             ! calculate stress tensor

             stress_temp(1) = stress_temp(1) + xrr * force_temp_comp(1)
             stress_temp(2) = stress_temp(2) + xrr * force_temp_comp(2)
             stress_temp(3) = stress_temp(3) + xrr * force_temp_comp(3)
             stress_temp(5) = stress_temp(5) + yrr * force_temp_comp(2)
             stress_temp(6) = stress_temp(6) + yrr * force_temp_comp(3)
             stress_temp(9) = stress_temp(9) + zrr * force_temp_comp(3)
          End Do

          Do jj = nz_fr(0) + nz_fr(comm%idnode + 1) + 1, nzfr ! +, on nodes>comm%idnode
             xrr = xfr(ii) - xfr(jj)
             yrr = yfr(ii) - yfr(jj)
             zrr = zfr(ii) - zfr(jj)

             xss = (rcell(1) * xrr + rcell(4) * yrr + rcell(7) * zrr)
             yss = (rcell(2) * xrr + rcell(5) * yrr + rcell(8) * zrr)
             zss = (rcell(3) * xrr + rcell(6) * yrr + rcell(9) * zrr)

             xss = xss - Anint(xss)
             yss = yss - Anint(yss)
             zss = zss - Anint(zss)

             xrr = (config%cell(1) * xss + config%cell(4) * yss + config%cell(7) * zss)
             yrr = (config%cell(2) * xss + config%cell(5) * yss + config%cell(8) * zss)
             zrr = (config%cell(3) * xss + config%cell(6) * yss + config%cell(9) * zss)

             ! calculate interatomic distance

             mod_r_ij_2 = xrr**2 + yrr**2 + zrr**2

             mod_r_ij = Sqrt(mod_r_ij_2)
             atom_coeffs_ij = cfr(ii) * cfr(jj) * scl

             ! calculate error function and derivative

             exp1 = Exp(-(ewld%alpha * mod_r_ij)**2)
             tt = 1.0_wp / (1.0_wp + pp * ewld%alpha * mod_r_ij)

             erfr = atom_coeffs_ij * &
                  (1.0_wp - tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1) / mod_r_ij

             erf_gamma = -(erfr - 2.0_wp * atom_coeffs_ij * (ewld%alpha / sqrpi) * exp1) / mod_r_ij_2

             force_temp_comp(1) = erf_gamma * xrr
             force_temp_comp(2) = erf_gamma * yrr
             force_temp_comp(3) = erf_gamma * zrr

             ! calculate forces

             config%parts(l_ind(i))%fxx = config%parts(l_ind(i))%fxx - force_temp_comp(1)
             config%parts(l_ind(i))%fyy = config%parts(l_ind(i))%fyy - force_temp_comp(2)
             config%parts(l_ind(i))%fzz = config%parts(l_ind(i))%fzz - force_temp_comp(3)

             ! calculate potential energy and virial

             engcpe_fr = engcpe_fr - erfr
             vircpe_fr = vircpe_fr - erf_gamma * mod_r_ij_2

             ! calculate stress tensor

             stress_temp(1) = stress_temp(1) + xrr * force_temp_comp(1)
             stress_temp(2) = stress_temp(2) + xrr * force_temp_comp(2)
             stress_temp(3) = stress_temp(3) + xrr * force_temp_comp(3)
             stress_temp(5) = stress_temp(5) + yrr * force_temp_comp(2)
             stress_temp(6) = stress_temp(6) + yrr * force_temp_comp(3)
             stress_temp(9) = stress_temp(9) + zrr * force_temp_comp(3)
          End Do
       End Do

       Deallocate (cfr, xfr, yfr, zfr, Stat=fail)
       If (fail > 0) Then
          Write (message, '(a)') 'ewald_frzn_forces deallocation failure 1'
          Call error(0, message)
       End If

    Else

       ! We resort to approximating N*(N-1)/2 interactions
       ! with the short-range one from the two body linked config%cell neigh%list

       Allocate (x_pos(1:neigh%max_list), y_pos(1:neigh%max_list), z_pos(1:neigh%max_list), dr_j(1:neigh%max_list), Stat=fail)
       If (fail > 0) Then
          Write (message, '(a)') 'ewald_frzn_forces allocation failure 2'
          Call error(0, message)
       End If

       Do ii = 1, nz_fr(comm%idnode + 1)
          i = l_ind(nz_fr(comm%idnode + 1))
          global_id_i = config%ltg(ii)

          ! Get neigh%list limit

          limit = neigh%list(-2, i) - neigh%list(-1, i)
          If (limit > 0) Then

             ! calculate interatomic distances

             Do k = 1, limit
                j = neigh%list(neigh%list(-1, i) + k, i)

                x_pos(k) = config%parts(i)%xxx - config%parts(j)%xxx
                y_pos(k) = config%parts(i)%yyy - config%parts(j)%yyy
                z_pos(k) = config%parts(i)%zzz - config%parts(j)%zzz
             End Do

             ! periodic boundary conditions not needed by LC construction
             !
             !           Call images(config%imcon,config%cell,limit,x_pos,y_pos,z_pos)

             ! square of distances

             Do k = 1, limit
                dr_j(k) = Sqrt(x_pos(k)**2 + y_pos(k)**2 + z_pos(k)**2)
             End Do

             Do k = 1, limit
                j = neigh%list(neigh%list(-1, i) + k, i)

                mod_r_ij = dr_j(k)
                If (Abs(config%parts(j)%chge) > zero_plus .and. mod_r_ij < neigh%cutoff) Then
                   atom_coeffs_ij = config%parts(i)%chge * config%parts(j)%chge * scl
                   mod_r_ij_2 = mod_r_ij**2

                   ! calculate error function and derivative

                   exp1 = Exp(-(ewld%alpha * mod_r_ij)**2)
                   tt = 1.0_wp / (1.0_wp + pp * ewld%alpha * mod_r_ij)

                   erfr = atom_coeffs_ij * &
                        (1.0_wp - tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1) / mod_r_ij

                   erf_gamma = -(erfr - 2.0_wp * atom_coeffs_ij * (ewld%alpha / sqrpi) * exp1) / mod_r_ij_2

                   force_temp_comp(1) = erf_gamma * x_pos(k)
                   force_temp_comp(2) = erf_gamma * y_pos(k)
                   force_temp_comp(3) = erf_gamma * z_pos(k)

                   ! calculate forces

                   config%parts(i)%fxx = config%parts(i)%fxx - force_temp_comp(1)
                   config%parts(i)%fyy = config%parts(i)%fyy - force_temp_comp(2)
                   config%parts(i)%fzz = config%parts(i)%fzz - force_temp_comp(3)

                   If (j <= config%natms) Then

                      config%parts(j)%fxx = config%parts(j)%fxx + force_temp_comp(1)
                      config%parts(j)%fyy = config%parts(j)%fyy + force_temp_comp(2)
                      config%parts(j)%fzz = config%parts(j)%fzz + force_temp_comp(3)

                   End If

                   If (j <= config%natms .or. global_id_i < config%ltg(j)) Then

                      ! calculate potential energy and virial

                      engcpe_fr = engcpe_fr - erfr
                      vircpe_fr = vircpe_fr - erf_gamma * mod_r_ij_2

                      ! calculate stress tensor

                      stress_temp(1) = stress_temp(1) + x_pos(k) * force_temp_comp(1)
                      stress_temp(2) = stress_temp(2) + x_pos(k) * force_temp_comp(2)
                      stress_temp(3) = stress_temp(3) + x_pos(k) * force_temp_comp(3)
                      stress_temp(5) = stress_temp(5) + y_pos(k) * force_temp_comp(2)
                      stress_temp(6) = stress_temp(6) + y_pos(k) * force_temp_comp(3)
                      stress_temp(9) = stress_temp(9) + z_pos(k) * force_temp_comp(3)

                   End If
                End If
             End Do

          End If
       End Do

       Deallocate (x_pos, y_pos, z_pos, dr_j, Stat=fail)
       If (fail > 0) Call error_dealloc('position arrays', 'ewald_frzn_forces')

    End If

    ! complete stress tensor

    stress(1) = stress(1) + stress_temp(1)
    stress(2) = stress(2) + stress_temp(2)
    stress(3) = stress(3) + stress_temp(3)
    stress(4) = stress(4) + stress_temp(2)
    stress(5) = stress(5) + stress_temp(5)
    stress(6) = stress(6) + stress_temp(6)
    stress(7) = stress(7) + stress_temp(3)
    stress(8) = stress(8) + stress_temp(6)
    stress(9) = stress(9) + stress_temp(9)

    Deallocate (l_ind, nz_fr, Stat=fail)
    If (fail > 0) Then
       Write (message, '(a)') 'ewald_frzn_forces deallocation failure'
       Call error(0, message)
    End If
  End Subroutine ewald_frzn_forces

!!! Internals

  Subroutine ewald_spme_init(domain, max_atoms, comm, kspace_in, &
    & bspline_in, charge_grid, potential_grid, stress_grid)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to initialise the ewald SPME routines
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!

    Use ewald, Only: ewald_type
    Use bspline, Only: bspline_type, bspline_coeffs_gen
    Use domains, Only: domains_type
    Use parallel_fft, Only: initialize_fft, pfft_indices
    Use comms, Only: comms_type
    Use kspace, Only: setup_kspace
    Type(domains_type),                                Intent(In   ) :: domain
    Integer,                                           Intent(In   ) :: max_atoms
    Type(comms_type),                                  Intent(In   ) :: comm
    Type(kspace_type),                                 Intent(inout) :: kspace_in
    Type(bspline_type),                                Intent(inout) :: bspline_in
    Real(Kind=wp), Allocatable, Dimension(:, :, :),    Intent(  Out) :: charge_grid
    Complex(Kind=wp), Allocatable, Dimension(:, :, :), Intent(  Out) :: potential_grid, stress_grid

    Integer, Dimension(4) :: fail

    fail = 0

    Call setup_kspace(kspace_in, domain, (kspace_in%k_vec_dim), comm)

!!! begin cardinal b-splines set-up

    bspline_in%num_deriv = 2
    Allocate (bspline_in%derivs(3, 0:bspline_in%num_deriv, 1:bspline_in%num_splines, 1:max_atoms), stat=fail(1))
    If (fail(1) > 0) Call error_alloc('bspline_in%derivs', 'ewald_spme_init')

    ! calculate the global b-spline coefficients
    Call bspline_coeffs_gen(kspace_in, bspline_in)

!!! end cardinal b-splines set-up

!!! begin daft set-up

    ! workspace arrays for DaFT

    Allocate (charge_grid   (1:kspace_in%block_fac(1), 1:kspace_in%block_fac(2), 1:kspace_in%block_fac(3)), stat=fail(1))
    Allocate (potential_grid(1:kspace_in%block_fac(1), 1:kspace_in%block_fac(2), 1:kspace_in%block_fac(3)), stat=fail(2))
    Allocate (stress_grid   (1:kspace_in%block_fac(1), 1:kspace_in%block_fac(2), 1:kspace_in%block_fac(3)), stat=fail(3))
    Allocate (pfft_work     (1:kspace_in%block_fac(1), 1:kspace_in%block_fac(2), 1:kspace_in%block_fac(3)), stat=fail(4))
    If (Any(fail > 0)) Call error_alloc('SPME DaFT workspace arrays', 'ewald_spme_init')

!!! end daft set-up

  End Subroutine ewald_spme_init

  Subroutine spme_construct_charge_array(ncalc, ewld, lookup_array, recip_indices, electro, &
    & coeffs, charge_grid)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to construct the charge array for SPME calculations
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins & i.t.todorov & i.j.bush august 2018
    !!
    !!----------------------------------------------------------------------!

    Integer,                           Intent(In   ) :: ncalc
    Type(ewald_type),                  Intent(In   ) :: ewld
    Integer, Dimension(:),             Intent(In   ) :: lookup_array
    Integer, Dimension(:, :),          Intent(In   ) :: recip_indices
    Type(electrostatic_type),          Intent(In   ) :: electro
    Real(Kind=wp), Dimension(:),       Intent(In   ) :: coeffs
    Real(Kind=wp), Dimension(:, :, :), Intent(  Out) :: charge_grid

    Integer               :: atm, i, j, j_hi, j_lo, k, k_hi, k_lo, l, l_hi, l_lo
    Integer, Dimension(3) :: temp
    Real(Kind=wp)         :: atom_coeffs

    charge_grid = 0.0_wp

    ! construct 3d charge array
    ! daft version - use array that holds only the local data

    atom: Do atm = 1, ncalc

      i = lookup_array(atm)
      ! if a particle is charged and in the md cell or in its positive halo
      ! (t(i) >= 0) as the b-splines are negative directionally by propagation

      j_lo = Max(1, recip_indices(1, i) - ewld%kspace%domain_indices(1, 1) - ewld%bspline%num_splines + 3)
      k_lo = Max(1, recip_indices(2, i) - ewld%kspace%domain_indices(2, 1) - ewld%bspline%num_splines + 3)
      l_lo = Max(1, recip_indices(3, i) - ewld%kspace%domain_indices(3, 1) - ewld%bspline%num_splines + 3)
      j_hi = Min(ewld%kspace%domain_indices(1, 2), recip_indices(1, i) + 1) - ewld%kspace%domain_indices(1, 1) + 1
      k_hi = Min(ewld%kspace%domain_indices(2, 2), recip_indices(2, i) + 1) - ewld%kspace%domain_indices(2, 1) + 1
      l_hi = Min(ewld%kspace%domain_indices(3, 2), recip_indices(3, i) + 1) - ewld%kspace%domain_indices(3, 1) + 1

      temp = recip_indices(:, i) - ewld%bspline%num_splines - ewld%kspace%domain_indices(:, 1) + 2

      atom_coeffs = coeffs(i)

      Do l = l_lo, l_hi
        Do k = k_lo, k_hi
          Do j = j_lo, j_hi
            charge_grid(j, k, l) = charge_grid(j, k, l) + atom_coeffs * &
              & ewld%bspline%derivs(1, 0, j - temp(1), i) * &
              & ewld%bspline%derivs(2, 0, k - temp(2), i) * &
              & ewld%bspline%derivs(3, 0, l - temp(3), i)
          End Do
        End Do
      End Do

    End Do atom

  End Subroutine spme_construct_charge_array

  Subroutine spme_construct_potential_grid(ewld, recip_cell, charge_grid, spme_datum, &
    & kernel, potential_grid, stress_contrib)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to construct the potential array for SPME calculations
    !! If given stress_contrib will also calculate the reciprocal contribution
    !! to the stress in the short way
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins & i.t.todorov & i.j.bush august 2018
    !!
    !!----------------------------------------------------------------------!

    Use comms, Only: gsum
    Use constants, Only: pi
    Use parallel_fft, Only: pfft
    Type(ewald_type),                      Intent(In   ) :: ewld
    Real(Kind=wp), Dimension(9),           Intent(In   ) :: recip_cell
    Real(Kind=wp), Dimension(:, :, :),     Intent(In   ) :: charge_grid
    Type(spme_component),                  Intent(In   ) :: spme_datum
    Complex(Kind=wp), External                           :: kernel
    Complex(Kind=wp), Dimension(:, :, :),  Intent(  Out) :: potential_grid
    Real(Kind=wp), Dimension(9), Optional, Intent(  Out) :: stress_contrib

    Complex(Kind=wp)               :: potential_component
    Integer                        :: alpha, beta, j, j_local, jj, k, k_local, kk, l, l_local, ll
    Real(Kind=wp)                  :: bb1, bb2, bb3, cut_off, cut_off_2, f_p_fac, k_vec_2, m, &
                                      pressure_virial, recip_conv_fac
    Real(Kind=wp), Dimension(10)   :: recip_cell_properties
    Real(Kind=wp), Dimension(3, 3) :: recip_pos, stress_temp

!! SPME contribution to the stress
!! Reciprocal lattice vectors
!! Core function to FT
!! bbb(1 to 3) - lengths of cell vectors: a(x,y,z) , b(x,y,z) , c(x,y,z)
!! bbb(4 to 6) - cosines of cell angles: gamma(a,b) , beta(a,c) , alpha(b,c)
!! bbb(7 to 9) - perpendicular cell widths : wx(y,z) , wy(x,z) , wz(x,y)
!! bbb(10)     - cell volume
!! Reciprocal space cut-off
!! Reciprocal space cut-off squared
!! Current Step in K vectors
!! Magnitude of RK vectors
!! Reciprocal convergence factor
!! pi*m/beta
!! Contribution to the potential
!! Virial contribution to the pressure

    recip_conv_fac = 1.0_wp / ewld%alpha

    ! set reciprocal space cutoff

    Call dcell(recip_cell, recip_cell_properties)

    cut_off = 0.5_wp * 1.05_wp * Minval(ewld%kspace%k_vec_dim_real * recip_cell_properties(7:9))
    cut_off_2 = cut_off**2

    ! load charge array into complex array for fft

    potential_grid = Cmplx(charge_grid, Kind=wp)

    ! calculate inverse 3d fft of charge array (in place)

    Call pfft(potential_grid, pfft_work, ewld%kspace%context, 1)

    ! initialise temporary stress tensor

    stress_temp = 0.0_wp

    ! calculate convolution of charge array with gaussian function
    ! daft version - only loop over the local stuff

    Do l_local = 1, ewld%kspace%block_fac(3)
      l = ewld%kspace%index(3, l_local)

      ll = l - 1
      If (2 * ll > ewld%kspace%k_vec_dim(3)) ll = ll - ewld%kspace%k_vec_dim(3)

      recip_pos(:, 3) = Real(ll, wp) * recip_cell(3:9:3)
      bb3 = ewld%bspline%norm2(3, l)

      Do k_local = 1, ewld%kspace%block_fac(2)
        k = ewld%kspace%index(2, k_local)

        kk = k - 1
        If (2 * kk > ewld%kspace%k_vec_dim(2)) kk = kk - ewld%kspace%k_vec_dim(2)

        recip_pos(:, 2) = recip_pos(:, 3) + Real(kk, wp) * recip_cell(2:9:3)
        bb2 = bb3 * ewld%bspline%norm2(2, k)

        Do j_local = 1, ewld%kspace%block_fac(1)
          j = ewld%kspace%index(1, j_local)

          jj = j - 1
          If (2 * jj > ewld%kspace%k_vec_dim(1)) jj = jj - ewld%kspace%k_vec_dim(1)

          recip_pos(:, 1) = recip_pos(:, 2) + Real(jj, wp) * recip_cell(1:9:3)
          bb1 = bb2 * ewld%bspline%norm2(1, j)

          k_vec_2 = Dot_product(recip_pos(:, 1), recip_pos(:, 1))

          If (k_vec_2 <= cut_off_2) Then

            m = Sqrt(k_vec_2)
            f_p_fac = pi * Sqrt(k_vec_2) * recip_conv_fac

            potential_component = kernel(bb1, potential_grid(j_local, k_local, l_local), &
              & f_p_fac, ewld%alpha, spme_datum%pot_order)

            ! By L'Hopital's rule, m=0 does not contribute to stress
            If (Present(stress_contrib) .and. k_vec_2 > 1.0e-6_wp) Then

              pressure_virial = Real(stress_kernel(bb1, potential_grid(j_local, k_local, l_local), &
                & f_p_fac, ewld%alpha, spme_datum%pot_order)  &
                & * Conjg(potential_grid(j_local, k_local, l_local)), wp) * pi * recip_conv_fac

              Do alpha = 1, 3
                Do beta = 1, 3
                  stress_temp(beta, alpha) = stress_temp(beta, alpha) + &
                    & recip_pos(alpha, 1) * recip_pos(beta, 1) * pressure_virial
                End Do
              End Do
            End If
          Else
            potential_component = (0.0_wp, 0.0_wp)
          End If

          potential_grid(j_local, k_local, l_local) = potential_component

        End Do
      End Do
    End Do

    If (Present(stress_contrib)) stress_contrib = Reshape(stress_temp, [9])

    Call pfft(potential_grid, pfft_work, ewld%kspace%context, -1)

  End Subroutine spme_construct_potential_grid

  Subroutine spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, recip_cell, &
    & recip_indices, potential_grid, per_part_step, energies, forces)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to calculate the per-particle energy and force
    !! contributions from the SPME Coulombic potential
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Use comms, Only: gsum
    Use domains, Only: exchange_grid
    Type(ewald_type),                     Intent(In   ) :: ewld
    Type(electrostatic_type),             Intent(In   ) :: electro
    Type(comms_type),                     Intent(inout) :: comm
    Type(domains_type),                   Intent(In   ) :: domain
    Type(configuration_type),             Intent(In   ) :: config
    Real(Kind=wp), Dimension(:),          Intent(In   ) :: coeffs
    Real(Kind=wp), Dimension(9),          Intent(In   ) :: recip_cell
    Integer, Dimension(:, :),             Intent(In   ) :: recip_indices
    Complex(Kind=wp), Allocatable, Dimension(:, :, :), Intent(In   ) :: potential_grid
    Logical,                              Intent(In   ) :: per_part_step
    Real(Kind=wp), Dimension(0:),         Intent(  Out) :: energies
    Real(Kind=wp), Dimension(:, :),       Intent(  Out) :: forces

    Character(Len=256)                                   :: message
    Integer                                              :: fail, i, j, jj, k, kk, l, ll
    Integer, Dimension(3, 2), Save                       :: extended_domain
    Integer, Save                                        :: mxspl2_old = -1
    Real(Kind=wp)                                        :: atom_coeffs, energy_total
    Real(Kind=wp), Allocatable, Dimension(:, :, :), Save :: extended_potential_grid
    Real(Kind=wp), Dimension(1:2)                        :: energy_temp
    Real(Kind=wp), Dimension(3)                          :: curr_force_temp, force_total, &
                                                            recip_kmax
    Real(Kind=wp), Dimension(3, 1:3)                     :: force_temp
    Real(Kind=wp), Dimension(3, 3)                       :: recip_cell_mat
    Real(Kind=wp), Dimension(ewld%bspline%num_splines)   :: bspline_d0_x, bspline_d0_y, &
                                                            bspline_d0_z, bspline_d1_x, &
                                                            bspline_d1_y, bspline_d1_z

!! List of per-particle energy contributions
!! Per-particle forces
!! Grid containing back FT'd potential
!! Coefficients such as charges or potential
!! Reciprocal lattice vectors
!! Reciprocal grid locations of charge centres
!! Whether to perform per-particle measurements
!! Grid with extended halo splines
!! Temporary energy components
!! Total sum of atomic energy
!! Temporary matrix used in force calculation
!! Total sum of all forces
!! Temporary force vec
!! In matrix form
!! Size of extended grid with halo splines

    recip_cell_mat = Reshape(recip_cell, [3, 3])
    recip_kmax = Matmul(recip_cell_mat, ewld%kspace%k_vec_dim_real)

    ! Exchange grid
    If (ewld%bspline%num_spline_padded .ne. mxspl2_old) Then
      mxspl2_old = ewld%bspline%num_spline_padded
      extended_domain(:, 1) = ewld%kspace%domain_indices(:, 1) - ewld%bspline%num_spline_padded
      extended_domain(:, 2) = ewld%kspace%domain_indices(:, 2) + ewld%bspline%num_spline_padded &
        & - ewld%bspline%num_splines

      If (Allocated(extended_potential_grid)) Then
        Deallocate (extended_potential_grid, stat=fail)
        If (fail /= 0) Call error_dealloc('extended_potential_grid', 'spme_calc_force_energy')
      End If

      Allocate (extended_potential_grid( &
        & extended_domain(1, 1):extended_domain(1, 2), &
        & extended_domain(2, 1):extended_domain(2, 2), &
        & extended_domain(3, 1):extended_domain(3, 2)), stat=fail)
      If (fail /= 0) Call error_alloc('extended_potential_grid', 'spme_calc_force_energy')
    End If

    Call exchange_grid(potential_grid, ewld%kspace%domain_indices(:,1), ewld%kspace%domain_indices(:,2), &
      & extended_potential_grid, extended_domain(:,1), extended_domain(:,2), domain, comm)

    If (Any(Minval(recip_indices(:, 1:config%natms), dim=2) + 2 - ewld%bspline%num_splines < extended_domain(:, 1)) .or. &
        Any(Maxval(recip_indices(:, 1:config%natms), dim=2) + 1 > extended_domain(:, 2))) Then
      Write (message, '(A)') 'Atoms beyond box bounds, unstable system'
      Call error(0, message)
    End If

    ! Zero accumulators, energies and forces
    energies(0) = 0.0_wp
    forces = 0.0_wp
    force_total = 0.0_wp

    ! Calculate per-particle contributions
    atom: Do i = 1, config%natms

      energy_total = 0.0_wp
      curr_force_temp = 0.0_wp
      atom_coeffs = coeffs(i)

      bspline_d0_x = ewld%bspline%derivs(1, 0, :, i)
      bspline_d0_y = ewld%bspline%derivs(2, 0, :, i)
      bspline_d0_z = ewld%bspline%derivs(3, 0, :, i)

      bspline_d1_x = ewld%bspline%derivs(1, 1, :, i)
      bspline_d1_y = ewld%bspline%derivs(2, 1, :, i)
      bspline_d1_z = ewld%bspline%derivs(3, 1, :, i)

      Do l = 1, ewld%bspline%num_splines
        ll = recip_indices(3, i) + 1 - ewld%bspline%num_splines + l

        energy_temp(2) = atom_coeffs * bspline_d0_z(l)

        force_temp(1, 3) = atom_coeffs * bspline_d0_z(l)
        force_temp(2, 3) = atom_coeffs * bspline_d0_z(l)
        force_temp(3, 3) = atom_coeffs * bspline_d1_z(l)

        Do k = 1, ewld%bspline%num_splines
          kk = recip_indices(2, i) + 1 - ewld%bspline%num_splines + k

          energy_temp(1) = energy_temp(2) * bspline_d0_y(k)

          force_temp(1, 2) = force_temp(1, 3) * bspline_d0_y(k)
          force_temp(2, 2) = force_temp(2, 3) * bspline_d1_y(k)
          force_temp(3, 2) = force_temp(3, 3) * bspline_d0_y(k)

          Do j = 1, ewld%bspline%num_splines
            jj = recip_indices(1, i) + 1 - ewld%bspline%num_splines + j

            force_temp(1, 1) = force_temp(1, 2) * bspline_d1_x(j) * extended_potential_grid(jj, kk, ll) * recip_kmax(1)
            force_temp(2, 1) = force_temp(2, 2) * bspline_d0_x(j) * extended_potential_grid(jj, kk, ll) * recip_kmax(2)
            force_temp(3, 1) = force_temp(3, 2) * bspline_d0_x(j) * extended_potential_grid(jj, kk, ll) * recip_kmax(3)

            ! Sum force contributions
            force_total = force_total - force_temp(:, 1)
            curr_force_temp = curr_force_temp + force_temp(:, 1)
            ! energy_total now holds omega_j * 2piV
            energy_total = energy_total + energy_temp(1) * bspline_d0_x(j) * extended_potential_grid(jj, kk, ll)

          End Do
        End Do
      End Do

      ! Add to total accumulators
      energies(0) = energies(0) + energy_total
      If (per_part_step) energies(i) = energy_total
      forces(:, i) = forces(:, i) - curr_force_temp

    End Do atom

    ! Correct for CoM term
    Call gsum(comm, force_total)

    force_total = force_total / Real(config%megatm, wp)
    ! Remove CoM
    Do i = 1, config%natms
      forces(:, i) = (forces(:, i) - force_total)
    End Do

  End Subroutine spme_calc_force_energy

  Subroutine spme_calc_stress(ewld, electro, comm, domain, config, coeffs, &
    & recip_cell, recip_indices, stress_grid, stress_out)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to calculate the per-particle energy and force
    !! contributions from the SPME Coulombic potential
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Use comms, Only: gsum
    Use constants, Only: twopi
    Use domains, Only: exchange_grid
    Type(ewald_type),                     Intent(In   ) :: ewld
    Type(electrostatic_type),             Intent(In   ) :: electro
    Type(comms_type),                     Intent(inout) :: comm
    Type(domains_type),                   Intent(In   ) :: domain
    Type(configuration_type),             Intent(In   ) :: config
    Real(Kind=wp), Dimension(:),          Intent(In   ) :: coeffs
    Real(Kind=wp), Dimension(9),          Intent(In   ) :: recip_cell
    Integer, Dimension(:, :),             Intent(In   ) :: recip_indices
    Complex(Kind=wp), Dimension(:, :, :), Intent(In   ) :: stress_grid
    Real(Kind=wp), Dimension(:, 0:),      Intent(  Out) :: stress_out

    Integer                                              :: alpha, beta, fail, i, j, jj, k, kk, l, &
                                                            ll
    Integer, Dimension(3, 2), Save                       :: extended_domain
    Integer, Save                                        :: mxspl2_old = -1
    Real(Kind=wp)                                        :: atom_coeffs, c_fac
    Real(Kind=wp), Allocatable, Dimension(:, :, :), Save :: extended_stress_grid
    Real(Kind=wp), Dimension(3)                          :: recip_kmax
    Real(Kind=wp), Dimension(3, 3)                       :: recip_cell_mat
    Real(Kind=wp), Dimension(3, 3, 0:3)                  :: stress_temp
    Real(Kind=wp), Dimension(ewld%bspline%num_splines)   :: bspline_d0_x, bspline_d0_y, &
                                                            bspline_d0_z, bspline_d1_x, &
                                                            bspline_d1_y, bspline_d1_z, &
                                                            bspline_d2_x, bspline_d2_y, &
                                                            bspline_d2_z

!! Output stress
!! Grid containing back FT'd stress_contrib
!! Coefficients such as charges or potentials
!! Reciprocal lattice vectors
!! Reciprocal grid locations of charge centres
!! Grid with extended halo splines
!! Diag, off-diag
!! In matrix form
!! Size of extended grid with halo splines

    recip_cell_mat = Reshape(recip_cell, [3, 3])
    recip_kmax = Matmul(recip_cell_mat, ewld%kspace%k_vec_dim_real)

    c_fac = 2.0_wp * twopi * ewld%alpha
    c_fac = 1.0_wp / c_fac

    ! Exchange grid
    If (ewld%bspline%num_spline_padded .ne. mxspl2_old) Then
      mxspl2_old = ewld%bspline%num_spline_padded
      extended_domain(:, 1) = ewld%kspace%domain_indices(:, 1) - ewld%bspline%num_spline_padded
      extended_domain(:, 2) = ewld%kspace%domain_indices(:, 2) + ewld%bspline%num_spline_padded &
        & - ewld%bspline%num_splines

      If (Allocated(extended_stress_grid)) Then
        Deallocate (extended_stress_grid, stat=fail)
        If (fail /= 0) Call error_dealloc('extended_stress_grid', 'spme_calc_stress')
      End If
      Allocate (extended_stress_grid( &
        & extended_domain(1, 1):extended_domain(1, 2), &
        & extended_domain(2, 1):extended_domain(2, 2), &
        & extended_domain(3, 1):extended_domain(3, 2)), stat=fail)
      If (fail /= 0) Call error_alloc('extended_stress_grid', 'spme_calc_stress')
    End If

    Call exchange_grid(stress_grid, ewld%kspace%domain_indices(:,1), ewld%kspace%domain_indices(:,2), &
      & extended_stress_grid, extended_domain(:,1), extended_domain(:,2), domain, comm)

    ! Zero accumulator
    stress_out(:, 0) = 0.0_wp

    ! Calculate per-particle contributions
    atom: Do i = 1, config%natms

      stress_temp = 0.0_wp
      atom_coeffs = coeffs(i)

      bspline_d0_x = ewld%bspline%derivs(1, 0, :, i)
      bspline_d0_y = ewld%bspline%derivs(2, 0, :, i)
      bspline_d0_z = ewld%bspline%derivs(3, 0, :, i)

      bspline_d1_x = ewld%bspline%derivs(1, 1, :, i)
      bspline_d1_y = ewld%bspline%derivs(2, 1, :, i)
      bspline_d1_z = ewld%bspline%derivs(3, 1, :, i)

      bspline_d2_x = ewld%bspline%derivs(1, 2, :, i)
      bspline_d2_y = ewld%bspline%derivs(2, 2, :, i)
      bspline_d2_z = ewld%bspline%derivs(3, 2, :, i)

      Do l = 1, ewld%bspline%num_splines
        ll = recip_indices(3, i) + 1 - ewld%bspline%num_splines + l

        stress_temp(1, 1, 3) = atom_coeffs * bspline_d0_z(l)
        stress_temp(2, 1, 3) = atom_coeffs * bspline_d0_z(l)
        stress_temp(3, 1, 3) = atom_coeffs * bspline_d1_z(l)
        stress_temp(2, 2, 3) = atom_coeffs * bspline_d0_z(l)
        stress_temp(3, 2, 3) = atom_coeffs * bspline_d1_z(l)
        stress_temp(3, 3, 3) = atom_coeffs * bspline_d2_z(l)

        Do k = 1, ewld%bspline%num_splines
          kk = recip_indices(2, i) + 1 - ewld%bspline%num_splines + k

          stress_temp(1, 1, 2) = stress_temp(1, 1, 3) * bspline_d0_y(k)
          stress_temp(2, 1, 2) = stress_temp(2, 1, 3) * bspline_d1_y(k)
          stress_temp(3, 1, 2) = stress_temp(3, 1, 3) * bspline_d0_y(k)
          stress_temp(2, 2, 2) = stress_temp(2, 2, 3) * bspline_d2_y(k)
          stress_temp(3, 2, 2) = stress_temp(3, 2, 3) * bspline_d1_y(k)
          stress_temp(3, 3, 2) = stress_temp(3, 3, 3) * bspline_d0_y(k)

          Do j = 1, ewld%bspline%num_splines
            jj = recip_indices(1, i) + 1 - ewld%bspline%num_splines + j

            stress_temp(1, 1, 1) = stress_temp(1, 1, 2) * bspline_d2_x(j) * extended_stress_grid(jj, kk, ll)
            stress_temp(2, 1, 1) = stress_temp(2, 1, 2) * bspline_d1_x(j) * extended_stress_grid(jj, kk, ll)
            stress_temp(3, 1, 1) = stress_temp(3, 1, 2) * bspline_d1_x(j) * extended_stress_grid(jj, kk, ll)
            stress_temp(2, 2, 1) = stress_temp(2, 2, 2) * bspline_d0_x(j) * extended_stress_grid(jj, kk, ll)
            stress_temp(3, 2, 1) = stress_temp(3, 2, 2) * bspline_d0_x(j) * extended_stress_grid(jj, kk, ll)
            stress_temp(3, 3, 1) = stress_temp(3, 3, 2) * bspline_d0_x(j) * extended_stress_grid(jj, kk, ll)

            Do beta = 1, 3
              Do alpha = beta, 3
                stress_temp(alpha, beta, 0) = stress_temp(alpha, beta, 0) - &
                  & (stress_temp(alpha, beta, 1) * recip_kmax(alpha) * recip_kmax(beta) * c_fac)
              End Do
            End Do

          End Do
        End Do
      End Do

      stress_temp(1, 2, 0) = stress_temp(2, 1, 0)
      stress_temp(1, 3, 0) = stress_temp(3, 1, 0)
      stress_temp(2, 3, 0) = stress_temp(3, 2, 0)

      stress_out(:, 0) = stress_out(:, 0) + Reshape(stress_temp(:, :, 0), [9])
      stress_out(:, i) = Reshape(stress_temp(:, :, 0), [9])

    End Do atom
  End Subroutine spme_calc_stress

!!! Kernels

  Function potential_kernel(B_m, pot, pi_m_over_a, conv_factor, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Kernel for calculating energy and forces for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    ! use spme,      only : f_p
    Real(Kind=wp)    :: B_m
    Complex(Kind=wp) :: pot
    Real(Kind=wp)    :: pi_m_over_a, conv_factor
    Integer          :: pot_order
    Complex(Kind=wp) :: potential_kernel

    potential_kernel = B_m * pot * f_p(pi_m_over_a, pot_order)

  End Function potential_kernel

  Function stress_kernel(B_m, pot, pi_m_over_a, conv_factor, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Kernel for calculating energy and forces for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    ! use spme,      only : f_p, f_p_d
    Use constants, Only: pi
    Real(Kind=wp)    :: B_m
    Complex(Kind=wp) :: pot
    Real(Kind=wp)    :: pi_m_over_a, conv_factor
    Integer          :: pot_order
    Complex(Kind=wp) :: stress_kernel

    Real(Kind=wp) :: energy

    If (pi_m_over_a > 1.0e-6) Then
      energy = f_p(pi_m_over_a, pot_order)
      ! (            1/m              )
      stress_kernel = B_m * pot * f_p_d(pi_m_over_a, energy, pot_order) * pi / pi_m_over_a / conv_factor
    Else
      stress_kernel = (0.0_wp, 0.0_wp)
    End If
  End Function stress_kernel

  Function f_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Nth order f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins november 2018
    !! based on  - i.j.bush igf.f90 november 2018
    !!----------------------------------------------------------------------!
    Use constants, Only: sqrpi
    Use numerics, Only: calc_erfc, calc_exp_int, calc_inv_gamma_1_2
    Use spme, Only: f_1, f_2, f_4, f_6, f_12
    Real(kind=wp), Intent(In   ) :: x
    Integer,       Intent(In   ) :: pot_order
    Real(kind=wp)                :: f_p

    Integer       :: curr_pot_order, p_work
    Real(kind=wp) :: base_integ, curr_xp, exp_xsq, x_2, x_fac, xp

    Select Case (pot_order)
    Case (1)
      If (x > 1.0e-6_wp) Then
        f_p = f_1(x)
      Else
        f_p = 0.0_wp
      End If

    Case (2)
      If (x > 1.0e-6_wp) Then
        f_p = f_2(x)
      Else
        f_p = 0.0_wp
      End If
    Case (4)
      f_p = f_4(x)
    Case (6)
      f_p = f_6(x)
    Case (12)
      f_p = f_12(x)
    Case (:0)
      Call error(0, 'Invalid pot order in f_p')
    Case default

      x_2 = x**2
      exp_xsq = Exp(-x_2)
      p_work = 2 - pot_order

      If (Mod(p_work, 2) == 0) Then
        ! even integrals base is I( 0, x )
        base_integ = 0.5_wp * sqrpi * calc_erfc(x)
        curr_pot_order = 0
        xp = 1.0_wp
      Else
        If (p_work > 0) Then
          ! positive odd integrals base is I( 1, x )
          base_integ = 0.5_wp * exp_xsq
          curr_pot_order = 1
          xp = x
        Else
          ! negative odd integrals, base is I( -1, x ), which is 0.5 * E1( x * x )
          ! where e1 is the first order exponential integral
          base_integ = -0.5_wp * calc_exp_int(-x_2)
          curr_pot_order = -1
          xp = 1.0_wp / x
        End If
      End If

      f_p = base_integ

      If (curr_pot_order == p_work) Then
        If (x < 1.0e-6_wp) Then
          f_p = 0.0_wp ! if p < 3 && x is small
          Return
        End If

        Continue
      Else If (curr_pot_order > p_work) Then
        ! recurse down
        x_fac = 1.0_wp / x_2
        curr_xp = xp / x
        Do curr_pot_order = curr_pot_order, p_work + 1, -2
          ! f_p = 2.0_wp * ( f_p - 0.5_wp * x**(curr_pot_order - 1) * exp_xsq ) / real( curr_pot_order - 1, wp )
          f_p = 2.0_wp * (f_p - 0.5_wp * curr_xp * exp_xsq) / Real(curr_pot_order - 1, wp)
          curr_xp = curr_xp * x_fac
        End Do
      Else
        ! recurse up
        x_fac = x_2
        curr_xp = xp * x
        Do curr_pot_order = curr_pot_order, p_work - 1, 2
          ! not tested !!!!! 5/11/18
          ! f_p = 0.5_wp * x ** ( p_now + 1 ) ) * exp_xsq + 0.5_wp * ( curr_pot_order + 1 ) * f_p
          f_p = 0.5_wp * (curr_xp * exp_xsq + Real(curr_pot_order + 1, wp) * f_p)
          curr_xp = curr_xp * x_fac
        End Do
      End If

      f_p = 2.0_wp * x**(pot_order - 3) * calc_inv_gamma_1_2(pot_order) * f_p
    End Select

  End Function f_p

  Function f_p_d(x, energy, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Use constants, Only: inv_gamma_1_2
    Real(Kind=wp) :: x, energy
    Integer       :: pot_order
    Real(Kind=wp) :: f_p_d

    f_p_d = (Real(pot_order - 3, wp) / x * energy) - ((2.0_wp / x) * inv_gamma_1_2(pot_order) * Exp(-(x**2)))

  End Function f_p_d

  Function g_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! General g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Use constants, Only: rsqrpi
    Use numerics, Only: factorial
    Use spme, Only: g_1, g_2, g_6, g_12
    Real(Kind=wp), Intent(In   ) :: x
    Integer,       Intent(In   ) :: pot_order
    Real(Kind=wp)                :: g_p

    Integer       :: i
    Real(Kind=wp) :: den, num, x_2, x_curr

    Select Case (pot_order)
    Case (1)
      g_p = g_1(x)
    Case (2)
      g_p = g_2(x)
    Case (6)
      g_p = g_6(x)
    Case (12)
      g_p = g_12(x)
    Case (:0)
      Call error(0, 'Invalid pot order in g_p')
    Case default

      x_2 = x**2

      If (Mod(pot_order, 2) == 0) Then !Even orders

        g_p = 0.0_wp
        Do i = 0, (pot_order / 2) - 1
          g_p = g_p + Exp(-factorial(i)) * x_2**i
        End Do

        g_p = g_p * Exp(-x_2)

      Else ! Odd orders

        x_curr = x
        g_p = 0.0_wp
        den = 1.0_wp
        num = 1.0_wp

        Do i = 1, pot_order - 1, 2
          num = 2.0_wp * num
          den = den / Real(i, wp)
          g_p = g_p + x_curr * num * den
          x_curr = x_curr * x_2
        End Do
        g_p = g_p * Exp(-x_2) * rsqrpi
        g_p = g_p + calc_erfc(x)

      End If

      g_p = g_p
    End Select

  End Function g_p

  Function g_p_d(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Use constants, Only: inv_gamma_1_2
    Real(Kind=wp), Intent(In   ) :: x
    Integer,       Intent(In   ) :: pot_order
    Real(Kind=wp)                :: g_p_d

    g_p_d = 2.0_wp * inv_gamma_1_2(pot_order) * x**(pot_order - 1) * Exp(-x**2)

  End Function g_p_d

End Module ewald_spole
