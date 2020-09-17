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
  Use bspline,         Only: bspline_splines_gen
  Use comms,           Only: comms_type,&
                             gcheck,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: pi,&
                             r4pie0,&
                             sqrpi,&
                             zero_plus
  Use domains,         Only: domains_type
  Use electrostatic,   Only: electrostatic_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc
  Use ewald,           Only: ewald_type
  Use ewald_general,   Only: ewald_spme_init,&
                             spme_calc_force_energy,&
                             spme_calc_stress,&
                             spme_construct_charge_array,&
                             spme_construct_potential_grid_gen,&
                             stress_kernel
  Use kinds,           Only: wp
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: dcell,&
                             invert
  Use parallel_fft,    Only: pfft
  Use spme,            Only: spme_component
  Use statistics,      Only: calculate_stress,&
                             stats_type
  Use timer,           Only: timer_type

  Implicit None

  Private

  Public ::  ewald_real_forces_coul
  Public ::  ewald_spme_forces_coul
  Public ::  ewald_excl_forces, ewald_frzn_forces

Contains

  Subroutine ewald_real_forces_coul(electro, spme_datum, neigh, config, stats, iatm, x_pos, y_pos, z_pos, mod_dr_ij, &
    & engcpe_rl, vircpe_rl)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using tabulated pots for ewald's method
    !!
    !! note: Real space terms
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith august 1998
    !! amended   - i.t.todorov april 2015
    !! amended   - j. wilkins september 2018
    !! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    !!
    !!-----------------------------------------------------------------------
    Type(electrostatic_type),                   Intent(In   ) :: electro
    Type(spme_component),                       Intent(In   ) :: spme_datum
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(configuration_type),                   Intent(InOut) :: config
    Type(stats_type),                           Intent(InOut) :: stats
    Integer,                                    Intent(In   ) :: iatm
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: x_pos, y_pos, z_pos, mod_dr_ij
    Real(Kind=wp),                              Intent(  Out) :: engcpe_rl, vircpe_rl

    Integer                     :: global_id_i, global_id_j, jatm, m, nearest_sample_index
    Real(Kind=wp)               :: atom_coeffs_i, difference, e_comp, erf_gamma, mod_r_ij, prefac
    Real(Kind=wp), Dimension(9) :: stress_temp, stress_temp_comp
    Real(Kind=wp), Dimension(3) :: force_temp, force_temp_comp, points, pos_j
    Real(Kind=wp), Dimension(2) :: temp

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
        If (nearest_sample_index == 0) points(1) = points(1) * mod_r_ij
        temp(1) = points(1) + (points(2) - points(1)) * difference
        temp(2) = points(2) + (points(3) - points(2)) * (difference - 1.0_wp)
        erf_gamma = prefac * (temp(1) + (temp(2) - temp(1)) * difference * 0.5_wp)
        ! erf_gamma = prefac * electro%erfc_deric%calc(mod_r_ij)

        ! calculate forces ( dU * r/||r|| )

        force_temp_comp = erf_gamma * pos_j
        force_temp = force_temp + force_temp_comp

#ifndef HALF_HALO
        If (jatm <= config%natms .or. global_id_i < global_id_j .or. stats%collect_pp) Then
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - force_temp_comp(1)
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - force_temp_comp(2)
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - force_temp_comp(3)

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

          ! calculate components of G
          nearest_sample_index = Int(mod_r_ij * electro%erfc%recip_spacing)
          difference = mod_r_ij * electro%erfc%recip_spacing - Real(nearest_sample_index, wp)
          points = electro%erfc%table(nearest_sample_index:nearest_sample_index + 2)
          If (nearest_sample_index == 0) points(1) = points(1) * mod_r_ij
          temp(1) = points(1) + (points(2) - points(1)) * difference
          temp(2) = points(2) + (points(3) - points(2)) * (difference - 1.0_wp)
          e_comp = prefac * (temp(1) + (temp(2) - temp(1)) * difference * 0.5_wp)

#ifndef HALF_HALO
        End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
        If (jatm <= config%natms .or. global_id_i < global_id_j) Then
#endif /* HALF_HALO */
          !e_comp = prefac * electro%erfc%calc(mod_r_ij)

          ! calculate interaction energy
          engcpe_rl = engcpe_rl + e_comp
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

#ifndef HALF_HALO
        End If
#endif /* HALF_HALO */

        If (stats%collect_pp) Then
          stress_temp_comp = calculate_stress(pos_j, force_temp_comp)
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + e_comp * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + e_comp * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */
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

  Subroutine ewald_spme_forces_coul(ewld, spme_datum, domain, config, comm, coeffs, stats, &
    & engcpe_rc, vircpe_rc)
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
    Type(ewald_type),            Intent(inout) :: ewld
    Type(spme_component),        Intent(inout) :: spme_datum
    Type(domains_type),          Intent(In   ) :: domain
    Type(configuration_type),    Intent(inout) :: config
    Type(comms_type),            Intent(inout) :: comm
    Real(kind=wp), Dimension(:), Intent(In   ) :: coeffs
    Type(stats_type),            Intent(inout) :: stats
    Real(kind=wp),               Intent(  Out) :: engcpe_rc, vircpe_rc

    Integer                                                 :: dim, i
    Integer, Allocatable, Dimension(:)                      :: to_calc
    Integer, Dimension(4)                                   :: fail
    Logical                                                 :: llspl
    Real(kind=wp)                                           :: det, eng, rvolm, scale
    Real(kind=wp), Allocatable, Dimension(:)                :: Q_abc
    Real(kind=wp), Allocatable, Dimension(:, :)             :: F_abc, S_abc
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

    If (Any(Abs(coeffs) > zero_plus)) Then
      Continue
    Else
      Return
    End If

    If (ewld%newjob_spme) Then
      Call ewald_spme_init(domain, config%mxatms, comm, ewld%kspace, &
        & ewld%bspline)
      ewld%newjob_spme = .false.
    End If

    fail = 0
    Allocate (ewld%kspace%recip_coords(3, config%nlast), stat=fail(1))
    Allocate (ewld%kspace%recip_indices(3, config%nlast), stat=fail(2))
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

    llspl = .true.
    Do i = 1, config%nlast
      Do dim = 1, 3
        ewld%kspace%recip_coords(dim, i) = ewld%kspace%k_vec_dim_real(dim) * ( &
          & rcell(dim) * config%parts(i)%xxx + &
          & rcell(dim + 3) * config%parts(i)%yyy + &
          & rcell(dim + 6) * config%parts(i)%zzz + 0.5_wp)
      End Do

      ! if not dd bound in kmax grid space when .not.llvnl = (ewld%bspline%num_spline_pad == ewld%bspline)

      If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines .and. i <= config%natms) &
        & llspl = llspl .and. ( &
        & All(ewld%kspace%recip_coords(:, i) > ewld%kspace%domain_bounds(:, 1)) .and. &
        & All(ewld%kspace%recip_coords(:, i) < ewld%kspace%domain_bounds(:, 2)))

      ! detect if a particle is charged and in the md cell or in its positive halo
      ! (coords(i) >= -zero_plus) as the b-splines are negative directionally by propagation
      If (All(ewld%kspace%recip_coords(:, i) >= -zero_plus) .and. Abs(coeffs(i)) > zero_plus) Then
        to_calc(0) = to_calc(0) + 1
        to_calc(to_calc(0)) = i
      End If
    End Do

    ewld%kspace%recip_indices = Int(ewld%kspace%recip_coords)

    ! check for breakage of llspl when .not.llvnl = (ewld%bspline%num_spline_pad == ewld%bspline)

    ewld%bspline%num_spline_padded = ewld%bspline%num_spline_pad

    If (ewld%bspline%num_spline_pad == ewld%bspline%num_splines) Then
      Call gcheck(comm, llspl)
      If (.not. llspl) ewld%bspline%num_spline_padded = ewld%bspline%num_splines + 1
    End If

    ! construct b-splines for atoms

    Call bspline_splines_gen(config%nlast, ewld%kspace%recip_coords, ewld%bspline)

    Deallocate (ewld%kspace%recip_coords, stat=fail(1))
    If (fail(1) > 0) Call error_dealloc('recip_coords', 'ewald_spme_forces')

    Call spme_construct_charge_array(to_calc(0), ewld, to_calc(1:), ewld%kspace%recip_indices, coeffs, ewld%kspace%charge_grid)

    If (.not. stats%collect_pp) Then

      ! If we don't need per-particle data, we can use the old method of getting the stress (cheaper)
      Call spme_construct_potential_grid_coul(ewld, rcell, ewld%kspace%charge_grid, ewld%kspace%potential_grid, s_abc(:, 0))
      Call spme_calc_force_energy(ewld, comm, domain, config, coeffs, &
        & rcell, ewld%kspace%recip_indices, ewld%kspace%potential_grid, stats%collect_pp, q_abc, f_abc)

    Else

      Call spme_construct_potential_grid_coul(ewld, rcell, ewld%kspace%charge_grid, ewld%kspace%potential_grid, s_abc(:, 0))
      Call spme_construct_potential_grid_gen(ewld, rcell, ewld%kspace%charge_grid, &
        spme_datum, stress_kernel, ewld%kspace%stress_grid)

      Call spme_calc_force_energy(ewld, comm, domain, config, coeffs, &
        & rcell, ewld%kspace%recip_indices, ewld%kspace%potential_grid, stats%collect_pp, q_abc, f_abc)
      Call spme_calc_stress(ewld, comm, domain, config, coeffs, &
        & rcell, ewld%kspace%recip_indices, ewld%kspace%stress_grid, s_abc)

    End If

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

    Deallocate (ewld%kspace%recip_indices, stat=fail(1))
    ! deallocate (ewld%bspline%derivs, stat=fail(2))
    Deallocate (to_calc, stat=fail(3))
    Deallocate (Q_abc, F_abc, S_abc, stat=fail(4))
    If (Any(fail > 0)) Call error_dealloc('output_arrays', 'ewald_spme_forces')

  End Subroutine ewald_spme_forces_coul

  Subroutine ewald_excl_forces(iatm, xxt, yyt, zzt, rrt, engcpe_ex, vircpe_ex, stress, &
                               neigh, ewld, spme_datum, config)

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
    ! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                    Intent(In   ) :: iatm
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Real(Kind=wp), Dimension(1:9),              Intent(InOut) :: stress
    Real(Kind=wp),                              Intent(  Out) :: vircpe_ex, engcpe_ex
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: rrt, zzt, yyt, xxt
    Type(ewald_type),                           Intent(In   ) :: ewld
    Type(spme_component),                       Intent(In   ) :: spme_datum
    Type(configuration_type),                   Intent(InOut) :: config

    Real(Kind=wp), Parameter :: a1 = 0.254829592_wp, a2 = -0.284496736_wp, a3 = 1.421413741_wp, &
                                a4 = -1.453152027_wp, a5 = 1.061405429_wp, pp = 0.3275911_wp, &
                                r10 = 0.1_wp, r216 = 1.0_wp / 216.0_wp, r42 = 1.0_wp / 42.0_wp, &
                                rr3 = 1.0_wp / 3.0_wp

    Integer       :: idi, jatm, limit, m
    Real(Kind=wp) :: alpr, alpr2, chgea, chgprd, egamma, erfr, exp1, fix, fiy, fiz, fx, fy, fz, &
                     rrr, rsq, strs1, strs2, strs3, strs5, strs6, strs9, tt

    ! initialise potential energy and virial

    engcpe_ex = 0.0_wp
    vircpe_ex = 0.0_wp

    ! initialise stress tensor accumulators

    strs1 = 0.0_wp
    strs2 = 0.0_wp
    strs3 = 0.0_wp
    strs5 = 0.0_wp
    strs6 = 0.0_wp
    strs9 = 0.0_wp

    ! global identity of iatm

    idi = config%ltg(iatm)

    ! ignore interaction if the charge is zero

    chgea = config%parts(iatm)%chge

    If (Abs(chgea) > zero_plus) Then

      chgea = chgea * spme_datum%scaling

      ! load forces

      fix = config%parts(iatm)%fxx
      fiy = config%parts(iatm)%fyy
      fiz = config%parts(iatm)%fzz

      ! Get neigh%list limit

      limit = neigh%list(-1, iatm) - neigh%list(0, iatm)

      ! start of primary loop for forces evaluation

      Do m = 1, limit

        ! atomic index and charge

        jatm = neigh%list(neigh%list(0, iatm) + m, iatm)
        chgprd = config%parts(jatm)%chge

        ! interatomic distance

        rrr = rrt(m)

        ! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < neigh%cutoff) Then

          ! charge product

          chgprd = chgprd * chgea

          ! Squared distance

          rsq = rrr**2

          ! calculate forces

          alpr = rrr * ewld%alpha
          alpr2 = alpr * alpr

          ! calculate error function and derivative

          If (alpr < 1.0e-2_wp) Then

            ! close particles (core-shell units) - small distances limit

            erfr = 2.0_wp * chgprd * (ewld%alpha / sqrpi) * &
                   (1.0_wp + alpr2 * (-rr3 + alpr2 * (r10 + alpr2 * (-r42 + alpr2 * r216))))

            egamma = -4.0_wp * chgprd * (ewld%alpha**3 / sqrpi) * &
                     (rr3 + alpr2 * (-2.0_wp * r10 + alpr2 * (3.0_wp * r42 - 4.0_wp * alpr2 * r216)))

          Else

            ! distant particles - traditional

            exp1 = Exp(-(ewld%alpha * rrr)**2)
            tt = 1.0_wp / (1.0_wp + pp * ewld%alpha * rrr)

            erfr = chgprd * &
                   (1.0_wp - tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1) / rrr

            egamma = -(erfr - 2.0_wp * chgprd * (ewld%alpha / sqrpi) * exp1) / rsq

          End If

          ! calculate forces

          fx = egamma * xxt(m)
          fy = egamma * yyt(m)
          fz = egamma * zzt(m)

          fix = fix + fx
          fiy = fiy + fy
          fiz = fiz + fz

#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */

            config%parts(jatm)%fxx = config%parts(jatm)%fxx - fx
            config%parts(jatm)%fyy = config%parts(jatm)%fyy - fy
            config%parts(jatm)%fzz = config%parts(jatm)%fzz - fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
          If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

            ! add potential energy and virial

            engcpe_ex = engcpe_ex - erfr
            vircpe_ex = vircpe_ex - egamma * rsq

            ! add stress tensor

            strs1 = strs1 + xxt(m) * fx
            strs2 = strs2 + xxt(m) * fy
            strs3 = strs3 + xxt(m) * fz
            strs5 = strs5 + yyt(m) * fy
            strs6 = strs6 + yyt(m) * fz
            strs9 = strs9 + zzt(m) * fz

#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */

        End If

      End Do

      ! load back forces

      config%parts(iatm)%fxx = fix
      config%parts(iatm)%fyy = fiy
      config%parts(iatm)%fzz = fiz

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
    !! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    !!
    !!-----------------------------------------------------------------------
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

#ifndef HALF_HALO
              If (j <= config%natms) Then
#endif /* HALF_HALO */

                config%parts(j)%fxx = config%parts(j)%fxx + force_temp_comp(1)
                config%parts(j)%fyy = config%parts(j)%fyy + force_temp_comp(2)
                config%parts(j)%fzz = config%parts(j)%fzz + force_temp_comp(3)

#ifndef HALF_HALO
              End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
              If (j <= config%natms .or. global_id_i < config%ltg(j)) Then
#endif /* HALF_HALO */

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

#ifndef HALF_HALO
              End If
#endif /* HALF_HALO */
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

  Subroutine spme_construct_potential_grid_coul(ewld, recip_cell, charge_grid, potential_grid, stress_contrib)

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

    Type(ewald_type),                     Intent(InOut) :: ewld
    Real(Kind=wp), Dimension(9),          Intent(In   ) :: recip_cell
    Real(Kind=wp), Dimension(:, :, :),    Intent(In   ) :: charge_grid
    Complex(Kind=wp), Dimension(:, :, :), Intent(  Out) :: potential_grid
    Real(Kind=wp), Dimension(9),          Intent(  Out) :: stress_contrib

    Complex(Kind=wp)               :: potential_component
    Integer                        :: alpha, beta, j, j_local, jj, k, k_local, kk, l, l_local, ll
    Real(Kind=wp)                  :: bb1, bb2, bb3, cut_off, cut_off_2, f_p_fac, k_vec_2, m, &
                                      pressure_virial, recip_conv_fac, test_fac
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

    recip_conv_fac = pi / ewld%alpha
    test_fac = (1.0e-6_wp / recip_conv_fac)**2 ! f_p_fac > 1e-6
    ! set reciprocal space cutoff

    Call dcell(recip_cell, recip_cell_properties)

    cut_off = 0.5_wp * 1.05_wp * Minval(ewld%kspace%k_vec_dim_real * recip_cell_properties(7:9))
    cut_off_2 = cut_off**2

    ! load charge array into complex array for fft

    potential_grid = Cmplx(charge_grid, Kind=wp)

    ! calculate inverse 3d fft of charge array (in place)

    Call pfft(ewld%kspace%potential_grid, ewld%kspace%pfft_work, ewld%kspace%context, 1)

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

          If (k_vec_2 <= cut_off_2 .and. k_vec_2 > test_fac) Then

            m = Sqrt(k_vec_2)
            f_p_fac = Sqrt(k_vec_2) * recip_conv_fac

            potential_component = bb1 * potential_grid(j_local, k_local, l_local) * &
                 & Exp(-(f_p_fac**2)) / (sqrpi * f_p_fac**2)

            ! By L'Hopital's rule, m=0 does not contribute to stress

            pressure_virial = Real(potential_component * &
                                   (-2.0_wp * ((1.0_wp + f_p_fac**2) / k_vec_2)) * &
                                   Conjg(potential_grid(j_local, k_local, l_local)), wp)

            Do alpha = 1, 3
              Do beta = 1, 3
                stress_temp(beta, alpha) = stress_temp(beta, alpha) + &
                     & recip_pos(alpha, 1) * recip_pos(beta, 1) * pressure_virial
              End Do
            End Do
          Else
            potential_component = (0.0_wp, 0.0_wp)
          End If

          potential_grid(j_local, k_local, l_local) = potential_component

        End Do
      End Do
    End Do

    stress_contrib = Reshape(stress_temp, [9])

    Call pfft(potential_grid, ewld%kspace%pfft_work, ewld%kspace%context, -1)

  End Subroutine spme_construct_potential_grid_coul

End Module ewald_spole
