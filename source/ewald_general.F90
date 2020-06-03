Module ewald_general
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

  Public ::  ewald_real_forces_gen
  Public ::  ewald_spme_forces_gen
  Public ::  ewald_spme_init
  Public ::  spme_calc_force_energy
  Public ::  spme_construct_charge_array
  Public ::  spme_construct_potential_grid_gen
  Public ::  spme_calc_stress
  Public ::  stress_kernel

  !> Temp X,Y,Z Scaled Coords (U/mu)
  Real(kind=wp), Dimension(:, :), Allocatable       :: recip_coords
  !> Indices to avoid type conversion
  Integer, Dimension(:, :), Allocatable       :: recip_indices
  !> temporary workspace for parallel fft
  Complex(Kind=wp), Dimension(:, :, :), Allocatable :: pfft_work

Contains

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
    use spme, only : g_p, g_p_d

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

  Subroutine ewald_spme_forces_gen(ewld, spme_datum, electro, domain, config, comm, coeffs, stats, &
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
        & ewld%bspline, charge_grid, potential_grid, stress_grid, pfft_work)
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
      Call spme_construct_potential_grid_gen(ewld, rcell, charge_grid, spme_datum, &
        & potential_kernel, potential_grid, s_abc(:, 0))
      Call stop_timer(tmr, 'Potential')
      Call start_timer(tmr, 'ForceEnergy')
      Call spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, potential_grid, stats%collect_pp, q_abc, f_abc)
      Call stop_timer(tmr, 'ForceEnergy')

    Else

      Call spme_construct_potential_grid_gen(ewld, rcell, charge_grid, spme_datum, &
        & potential_kernel, potential_grid)
      Call spme_construct_potential_grid_gen(ewld, rcell, charge_grid, spme_datum, &
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

  End Subroutine ewald_spme_forces_gen

!!! Internals

  Subroutine ewald_spme_init(domain, max_atoms, comm, kspace_in, &
    & bspline_in, charge_grid, potential_grid, stress_grid, pfft_array)
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
    !> temporary workspace for parallel fft
    Complex(Kind=wp), Dimension(:, :, :), Allocatable, Intent(  Out) :: pfft_array


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
    Allocate (pfft_array     (1:kspace_in%block_fac(1), 1:kspace_in%block_fac(2), 1:kspace_in%block_fac(3)), stat=fail(4))
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
    Integer               :: q
    Integer               :: atm, i, j, j_hi, j_lo, k, k_hi, k_lo, l, l_hi, l_lo
    Integer, Dimension(3) :: temp
    Real(Kind=wp), Dimension(2) :: factor
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
      if (any([j_hi - j_lo, k_hi - k_lo, l_hi - l_lo] < 0)) cycle

      temp = recip_indices(:, i) - ewld%bspline%num_splines - ewld%kspace%domain_indices(:, 1) + 2
      atom_coeffs = coeffs(i)

      Do l = l_lo, l_hi
        factor(1) = atom_coeffs * ewld%bspline%derivs(3, 0, l - temp(3), i)
        Do k = k_lo, k_hi
          factor(2) = factor(2) * ewld%bspline%derivs(2, 0, k - temp(2), i)
          Do j = j_lo, j_hi
            charge_grid(j, k, l) = charge_grid(j, k, l) + factor(2) * ewld%bspline%derivs(1, 0, j - temp(1), i)
          End Do
        End Do
      End Do

    End Do atom

  End Subroutine spme_construct_charge_array

  Subroutine spme_construct_potential_grid_gen(ewld, recip_cell, charge_grid, spme_datum, &
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

  End Subroutine spme_construct_potential_grid_gen

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
    use spme,      only : f_p
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
    use spme,      only : f_p, f_p_d
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

End Module ewald_general
