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
  Use comms,           Only : comms_type
  Use configuration,   Only : configuration_type
  Use domains,         Only : domains_type
  Use electrostatic,   Only : electrostatic_type
  Use errors_warnings, Only : error, error_alloc, error_dealloc
  Use ewald,           Only : ewald_type
  Use kinds,           Only : wp, wi
  Use kspace,          Only : kspace_type
  Use neighbours,      Only : neighbours_type
  Use numerics,        Only : erfcgen, invert, dcell, calc_erf, calc_erf_deriv, calc_erfc
  Use numerics,        Only : invert, dcell
  Use parallel_fft,    Only : initialize_fft, pfft, pfft_indices
  Use constants,       Only : sqrpi,twopi,zero_plus,r4pie0
  Use timer,           Only : start_timer, stop_timer
  Use spme,            Only : spme_component
  Implicit None

  Private
  
  Public ::  ewald_real_forces_coul, ewald_real_forces_gen
  Public ::  ewald_spme_forces
  Public ::  ewald_excl_forces, ewald_frzn_forces

  !> Temp X,Y,Z Scaled Coords (U/mu)
  real( kind = wp ),    dimension( :,: ),     allocatable       :: recip_coords
  !> Indices to avoid type conversion
  integer,              dimension( :,: ),     allocatable       :: recip_indices
  !> temporary workspace for parallel fft
  Complex( Kind = wp ), Dimension( :,:,: ), Allocatable :: pfft_work

Contains

  Subroutine ewald_real_forces_coul(alpha,spme_datum,neigh,config,iatm,x_pos,y_pos,z_pos,mod_dr_ij, &
    & engcpe_rl,vircpe_rl,stress)

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
    use constants, only : rt2
    implicit none

    Type( neighbours_type ),                              Intent( In    ) :: neigh
    Type( configuration_type ),                           Intent( InOut ) :: config
    Type( spme_component ),                               Intent( In    ) :: spme_datum
    Real( Kind = wp ),                                    Intent( In    ) :: alpha
    Integer,                                              Intent( In    ) :: iatm                                         !! Current atom
    Real( Kind = wp ),     Dimension( 1:neigh%max_list ), Intent( In    ) :: x_pos,y_pos,z_pos,mod_dr_ij                  !! Atoms positions (neighbours, not global) and inter-particle separations
    Real( Kind = wp ),                                    Intent(   Out ) :: engcpe_rl,vircpe_rl                          !! Energy and virial for the Real component
    Real( Kind = wp ),     Dimension( 1:9 ),              Intent( InOut ) :: stress                                       !! Stress tensor

    Real( Kind = wp )                                                     :: atom_coeffs_i, atom_coeffs_ij                !! Atom storage of coeffs !!Dimension(size(coeffs,1))
    Real( Kind = wp ),     Dimension( 9 )                                 :: stress_temp                                  !! Tempeorary stress tensor
    Real( Kind = wp ),     Dimension( 3 )                                 :: force_temp, force_temp_comp                  !! Temporary force vectors
    Real( Kind = wp ),     Dimension( 3 )                                 :: pos_j                                        !! Position of ion j
    Real( Kind = wp ),     Dimension( 3 )                                 :: norm_pos_j
    Real( Kind = wp )                                                     :: g_fac, e_comp                                !! g_p & d/dr[g_p]
    Real( Kind = wp )                                                     :: erf_gamma                                    !! Q*g_p
    Real( Kind = wp )                                                     :: prefac                                       !! Coeffs*inv_mod_r_ij**n
    Real( Kind = wp )                                                     :: mod_r_ij, inv_mod_r_ij, alpha_r              !! Inter-particle distances
    Integer                                                               :: global_id_i
    Integer                                                               :: m,jatm

    ! initialise accumulators

    engcpe_rl=0.0_wp
    vircpe_rl=0.0_wp
    stress_temp = 0.0_wp
    force_temp = 0.0_wp
    
    ! global identity of iatm

    global_id_i = config%ltg(iatm)

    atom_coeffs_i = config%parts(iatm)%chge*spme_datum%scaling

    ! ignore interaction if the coeffs or scaling are zero
    if (abs(atom_coeffs_i) < zero_plus) return

    ! start of primary loop for forces evaluation

    Do m=1,neigh%list(0,iatm)

      ! atomic index and charge

      jatm=neigh%list(m,iatm)

      ! interatomic distance
      mod_r_ij=mod_dr_ij(m)
      prefac=config%parts(jatm)%chge

      ! interaction validity and truncation of potential
      if (abs(prefac) > zero_plus .and. mod_r_ij < neigh%cutoff) then

        pos_j = [x_pos(m),y_pos(m),z_pos(m)]
        alpha_r = mod_r_ij * alpha
        inv_mod_r_ij = 1.0_wp / mod_r_ij
      
        ! Complete prefactor
        prefac = atom_coeffs_i*prefac*inv_mod_r_ij

        ! calculate components of G
        e_comp = prefac * calc_erfc(alpha_r)

        ! Because function is g_p(ar)/(r^n)
        ! => -n*(g/r^(n + 1) + a(dg/dr))
        erf_gamma = prefac * ( e_comp * inv_mod_r_ij + rt2*exp(-alpha_r**2)*alpha )
        
        ! calculate forces ( dU * r/||r|| )

        force_temp_comp = erf_gamma*pos_j*inv_mod_r_ij
        force_temp = force_temp + force_temp_comp


        if (jatm <= config%natms .or. global_id_i < config%ltg(jatm)) then
          if (jatm <= config%natms) then
            
            config%parts(jatm)%fxx=config%parts(jatm)%fxx-force_temp_comp(1)
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-force_temp_comp(2)
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-force_temp_comp(3)
            
          end if

          ! calculate interaction energy
          engcpe_rl = engcpe_rl + e_comp

          ! calculate virial ( F.r )

          vircpe_rl = vircpe_rl - erf_gamma*mod_r_ij

          ! calculate stress tensor

          stress_temp(1:9:3) = stress_temp(1:9:3) + pos_j*force_temp_comp(1)
          stress_temp(2:9:3) = stress_temp(2:9:3) + pos_j*force_temp_comp(2)
          stress_temp(3:9:3) = stress_temp(3:9:3) + pos_j*force_temp_comp(3)
        end if

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx=config%parts(iatm)%fxx + force_temp(1)
    config%parts(iatm)%fyy=config%parts(iatm)%fyy + force_temp(2)
    config%parts(iatm)%fzz=config%parts(iatm)%fzz + force_temp(3)

    ! complete stress tensor

    stress = stress + stress_temp

  End Subroutine ewald_real_forces_coul
  
  Subroutine ewald_real_forces_gen(alpha,spme_datum,neigh,config,iatm,coeffs,x_pos,y_pos,z_pos,mod_dr_ij, &
    & engcpe_rl,vircpe_rl,stress)

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
    implicit none

    Type( neighbours_type ),                              Intent( In    ) :: neigh
    Type( configuration_type ),                           Intent( InOut ) :: config
    Type( spme_component ),                               Intent( In    ) :: spme_datum
    Real( Kind = wp ),                                    Intent( In    ) :: alpha
    Integer,                                              Intent( In    ) :: iatm                                         !! Current atom
    Real( Kind = wp ),     Dimension( : ),                Intent( In    ) :: coeffs                                       !! Coulomb charges/ multipole coeffs, etc.
    Real( Kind = wp ),     Dimension( 1:neigh%max_list ), Intent( In    ) :: x_pos,y_pos,z_pos,mod_dr_ij                  !! Atoms positions (neighbours, not global) and inter-particle separations
    Real( Kind = wp ),                                    Intent(   Out ) :: engcpe_rl,vircpe_rl                          !! Energy and virial for the Real component
    Real( Kind = wp ),     Dimension( 1:9 ),              Intent( InOut ) :: stress                                       !! Stress tensor

    Real( Kind = wp )                                                     :: atom_coeffs_i, atom_coeffs_ij                !! Atom storage of coeffs !!Dimension(size(coeffs,1))
    Real( Kind = wp ),     Dimension( 9 )                                 :: stress_temp                                  !! Tempeorary stress tensor
    Real( Kind = wp ),     Dimension( 3 )                                 :: force_temp, force_temp_comp                  !! Temporary force vectors
    Real( Kind = wp ),     Dimension( 3 )                                 :: pos_j                                        !! Position of ion j
    Real( Kind = wp ),     Dimension( 3 )                                 :: norm_pos_j
    Real( Kind = wp )                                                     :: g_fac, e_comp                                !! g_p & d/dr[g_p]
    Real( Kind = wp )                                                     :: erf_gamma                                    !! Q*g_p
    Real( Kind = wp )                                                     :: prefac                                       !! Coeffs*inv_mod_r_ij**n
    Real( Kind = wp )                                                     :: mod_r_ij, inv_mod_r_ij, alpha_r              !! Inter-particle distances
    Integer                                                               :: global_id_i
    Integer                                                               :: m,jatm

    ! initialise accumulators

    engcpe_rl=0.0_wp
    vircpe_rl=0.0_wp
    stress_temp = 0.0_wp
    force_temp = 0.0_wp
    
    ! global identity of iatm

    global_id_i = config%ltg(iatm)

    atom_coeffs_i = coeffs(iatm)*spme_datum%scaling

    ! ignore interaction if the coeffs or scaling are zero
    if (abs(atom_coeffs_i) < zero_plus) return

    ! start of primary loop for forces evaluation

    Do m=1,neigh%list(0,iatm)

      ! atomic index and charge

      jatm=neigh%list(m,iatm)

      ! interatomic distance
      mod_r_ij=mod_dr_ij(m)
      prefac = coeffs(jatm)

      ! interaction validity and truncation of potential
      if (abs(prefac) > zero_plus .and. mod_r_ij < neigh%cutoff) then

        pos_j = [x_pos(m),y_pos(m),z_pos(m)]
        alpha_r = mod_r_ij * alpha
        inv_mod_r_ij = 1.0_wp / mod_r_ij
      
        ! Complete prefactor
        prefac = atom_coeffs_i*prefac*inv_mod_r_ij**spme_datum%pot_order

        ! calculate components of G
        g_fac = g_p(alpha_r,spme_datum%pot_order)

        e_comp = prefac * g_fac

        ! Because function is g_p(ar)/(r^n)
        ! => -n*(g/r^(n + 1) + a(dg/dr))
        erf_gamma = prefac * ( g_p_d(alpha_r,spme_datum%pot_order)*alpha + &
          & spme_datum%pot_order * g_fac * inv_mod_r_ij )

        ! calculate forces ( dU * r/||r|| )

        force_temp_comp = erf_gamma*pos_j*inv_mod_r_ij
        force_temp = force_temp + force_temp_comp


        if (jatm <= config%natms .or. global_id_i < config%ltg(jatm)) then
          if (jatm <= config%natms) then
            
            config%parts(jatm)%fxx=config%parts(jatm)%fxx-force_temp_comp(1)
            config%parts(jatm)%fyy=config%parts(jatm)%fyy-force_temp_comp(2)
            config%parts(jatm)%fzz=config%parts(jatm)%fzz-force_temp_comp(3)
            
          end if

          ! calculate interaction energy
          engcpe_rl = engcpe_rl + e_comp

          ! calculate virial ( F.r )

          vircpe_rl = vircpe_rl - erf_gamma*mod_r_ij

          ! calculate stress tensor

          stress_temp(1:9:3) = stress_temp(1:9:3) + pos_j*force_temp_comp(1)
          stress_temp(2:9:3) = stress_temp(2:9:3) + pos_j*force_temp_comp(2)
          stress_temp(3:9:3) = stress_temp(3:9:3) + pos_j*force_temp_comp(3)
        end if

      end if

    end do

    ! load back forces

    config%parts(iatm)%fxx=config%parts(iatm)%fxx + force_temp(1)
    config%parts(iatm)%fyy=config%parts(iatm)%fyy + force_temp(2)
    config%parts(iatm)%fzz=config%parts(iatm)%fzz + force_temp(3)

    ! complete stress tensor

    stress = stress + stress_temp

  End Subroutine ewald_real_forces_gen

  subroutine ewald_spme_forces(ewld,spme_datum,electro,domain,config,comm,coeffs,nstep,engcpe_rc,vircpe_rc,stress)
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
    use comms,          only : gsum, gcheck, gsync
    use constants,      only : twopi, pi, sqrpi, zero_plus
    use bspline,        only : bspline_splines_gen

    implicit none

!!! Inputs and Outputs
    type( ewald_type ),                                  intent( inout ) :: ewld
    type( spme_component ),                                   intent( inout ) :: spme_datum
    type( electrostatic_type ),                               intent( in    ) :: electro
    type( domains_type ),                                     intent( in    ) :: domain
    type( configuration_type ),                               intent( inout ) :: config
    type( comms_type ),                                       intent( inout ) :: comm
    real( kind = wp ),     dimension(:),                      intent( in    ) :: coeffs              !! Coefficients such as charges or pot params
    integer,                                                  intent( in    ) :: nstep               !! Number of steps taken since calculation start
    real( kind = wp ),                                        intent(   out ) :: engcpe_rc,vircpe_rc !! Energy and virial of Coulomb interaction
    real( kind = wp ),     dimension(1:9),                    intent( inout ) :: stress              !! Output stress tensor

    integer,               dimension( : ), allocatable                        :: to_calc             !! List of points to calculate
    logical                                                                   :: llspl               !! Unknown? Does this want to be saved?
    integer                                                                   :: i, dim              !! Loop counters
    real( kind = wp ),     dimension(9)                                       :: rcell               !! Reciprocal lattice vectors

!!! Data constants and intermediate variables
    real( kind = wp ), dimension(9) :: stress_temp                                         !! Temporary stress tensor
    real( kind = wp )  :: det                                                              !! Determinant of inverse matrix
    real( kind = wp )  :: rvolm                                                            !! Reciprocal volume
    real( kind = wp )  :: scale                                                            !! Coulomb factor / epsq?
    real( kind = wp )  :: eng                                                              !! Energy contribution

    real( kind = wp ),    dimension( :,:,: ), allocatable, save :: charge_grid
    complex( kind = wp ), dimension( :,:,: ), allocatable, save :: potential_grid
    complex( kind = wp ), dimension( :,:,: ), allocatable, save :: stress_grid
    real( kind = wp ),    dimension( : ),     allocatable       :: Q_abc                   !! Per-particle Energies
    real( kind = wp ),    dimension( :,: ),   allocatable       :: F_abc                   !!     ""       Forces
    real( kind = wp ),    dimension( :,:,:),  allocatable       :: S_abc                   !!     ""       Stress
    integer, dimension(4) :: fail                                                          !! Ierr
    logical   :: per_part_step
    logical, save :: newjob = .true.

    call start_timer('Setup')
    per_part_step = mod(nstep,ewld%pp_write_freq) == 0 .and. ewld%pp_write_freq > 0

    if (newjob) then
      call ewald_spme_init(domain, config%mxatms, comm, ewld%kspace, &
        & ewld%bspline, charge_grid, potential_grid, stress_grid)
      newjob = .false.
    end if

    if ( all ( abs(coeffs) < zero_plus )) return

    fail=0
    allocate(recip_coords (3,config%mxatms), stat=fail(1))
    allocate(recip_indices(3,config%mxatms), stat=fail(2))
    allocate(to_calc      (0:config%mxatms), stat=fail(3))

    ! If not per-particle only need global sum, else need everything
    if (.not. per_part_step) then
      allocate(Q_abc(0:0), F_abc(3,config%natms), S_abc(3,3,0:0), stat=fail(4))
    else
      allocate(Q_abc(0:config%natms), F_abc(3,config%natms), S_abc(3,3,0:config%natms), stat=fail(4))
    end if
    if (any(fail > 0)) call error_alloc('output_arrays','ewald_spme_forces')

    ! Initialise accumulator

    to_calc(0) = 0

    ! initialise SPME potential energy and virial

    engcpe_rc = 0.0_wp
    vircpe_rc = 0.0_wp

    ! set working parameters

    rvolm=1.0_wp/config%volm

    ! set scaling constant

    scale=pi*sqrpi*ewld%alpha**(spme_datum%pot_order-3)*(0.5_wp*rvolm)*spme_datum%scaling

    ! calculate reciprocal cell

    call invert(config%cell,rcell,det)
    if (abs(det) < 1.0e-6_wp) call error(120)
    call stop_timer('Setup')

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

    call start_timer('Recip')
    llspl=.true.
    do i=1,config%nlast
      do dim = 1, 3
        recip_coords(dim,i) = ewld%kspace%k_vec_dim_real(dim)*( &
          & rcell(dim  )*config%parts(i)%xxx+ &
          & rcell(dim+3)*config%parts(i)%yyy+ &
          & rcell(dim+6)*config%parts(i)%zzz+0.5_wp)
      end do

      ! if not dd bound in kmax grid space when .not.llvnl = (ewld%bspline%num_spline_pad == ewld%bspline)

      if (ewld%bspline%num_spline_pad == ewld%bspline%num_splines .and. i <= config%natms) &
        & llspl = llspl .and. ( &
        & all(recip_coords(:,i) > ewld%kspace%domain_bounds(:,1)) .and. &
        & all(recip_coords(:,i) < ewld%kspace%domain_bounds(:,2)))

      ! detect if a particle is charged and in the md cell or in its positive halo
      ! (coords(i) >= -zero_plus) as the b-splines are negative directionally by propagation
      if (all (recip_coords (:,i) >= -zero_plus) .and. abs(coeffs(i)) > zero_plus) then
        to_calc(0) = to_calc(0) + 1
        to_calc(to_calc(0)) = i
      end if
    end do

    recip_indices = int (recip_coords)
    call stop_timer('Recip')

    call start_timer('BSpline')
    ! check for breakage of llspl when .not.llvnl = (ewld%bspline%num_spline_pad == ewld%bspline)

    ewld%bspline%num_spline_padded=ewld%bspline%num_spline_pad

    if (ewld%bspline%num_spline_pad == ewld%bspline%num_splines) then
      call gcheck(comm,llspl)
      if (.not.llspl) ewld%bspline%num_spline_padded=ewld%bspline%num_splines+1
    end if

    ! construct b-splines for atoms

    call bspline_splines_gen(config%nlast,recip_coords,ewld%bspline)

    deallocate (recip_coords, stat = fail(1))
    if (fail(1) > 0) call error_dealloc('recip_coords','ewald_spme_forces')
    call stop_timer('BSpline')

    call start_timer('Charge')
    call spme_construct_charge_array(to_calc(0),ewld,to_calc(1:),recip_indices, electro, coeffs, charge_grid)
    call stop_timer('Charge')
    if (.not.per_part_step .or. spme_datum%pot_order /= 1) then

      ! If we don't need per-particle data, we can use the old method of getting the stress (cheaper)
      call start_timer('Potential')
      call spme_construct_potential_grid(ewld, rcell, charge_grid, spme_datum, &
        & potential_kernel, potential_grid, s_abc(:,:,0))
      call stop_timer('Potential')
      call start_timer('ForceEnergy')
      call spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, potential_grid, per_part_step, q_abc, f_abc)
      call stop_timer('ForceEnergy')
    else

      call spme_construct_potential_grid(ewld, rcell, charge_grid, spme_datum, &
        & potential_kernel, potential_grid)
      call spme_construct_potential_grid(ewld, rcell, charge_grid, spme_datum, &
        & stress_kernel, stress_grid)

      call spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, potential_grid, per_part_step, q_abc, f_abc)
      call spme_calc_stress(ewld, electro, comm, domain, config, coeffs, &
        & rcell, recip_indices, stress_grid, s_abc)

    end if

    call start_timer('Output')
    ! Rescale to real space
    q_abc = q_abc * scale / real(comm%mxnode)
    f_abc = f_abc * scale * 2.0_wp
    s_abc = s_abc * scale / real(comm%mxnode)

    if (per_part_step) call write_per_part_contribs(config, comm, q_abc, f_abc, s_abc, nstep, spme_datum%pot_order)

    eng = Q_abc(0)

    ! as only looped over local stuff, we need to gsum the eng

    call gsum(comm, eng)

    engcpe_rc = eng + spme_datum%self_interaction

    ! add up forces

    config%parts(1:config%natms)%fxx = config%parts(1:config%natms)%fxx + f_abc(1,:)
    config%parts(1:config%natms)%fyy = config%parts(1:config%natms)%fyy + f_abc(2,:)
    config%parts(1:config%natms)%fzz = config%parts(1:config%natms)%fzz + f_abc(3,:)

    ! Put accumulated stress tensor into stress temp to save on cache problems and translate to linear regime

    stress_temp = reshape(s_abc(:,:,0), [9])

    ! as only looped over local stuff, we need to gsum stress

    call gsum(comm,stress_temp)

    ! scale strs and distribute per node

    stress_temp(1:9:4) = stress_temp(1:9:4) + eng
    stress = stress + stress_temp
    vircpe_rc = -sum(stress_temp(1:9:4))

    deallocate (recip_indices, stat=fail(1))
    ! deallocate (ewld%bspline%derivs, stat=fail(2))
    deallocate (to_calc, stat=fail(3))
    deallocate (Q_abc, F_abc, S_abc, stat=fail(4))
    if (any(fail > 0)) call error_dealloc('output_arrays','ewald_spme_forces')
    call stop_timer('Output')

  end subroutine ewald_spme_forces

  Subroutine ewald_excl_forces(ewld,spme_datum,neigh,electro,config,coeffs,iatm,x_pos,y_pos,z_pos,dr_j, &
    & engcpe_ex,vircpe_ex,stress)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for calculating coulombic energy and force terms
    !! in a periodic system using ewald's method
    !!
    !! Note: exclusion correction terms
    !!       frozen pairs are ignored by default, they are not dealt with here
    !!
    !! copyright - daresbury laboratory
    !! author    - i.t.todorov february 2015
    !!
    !!-----------------------------------------------------------------------

    Integer,                                          Intent( In    ) :: iatm
    Type( neighbours_type ),                          Intent( In    ) :: neigh
    Type( ewald_type ),                          Intent( In    ) :: ewld
    Type( electrostatic_type ),                       Intent( In    ) :: electro
    Type( spme_component ),                           Intent( In    ) :: spme_datum
    Type( configuration_type ),                       Intent( InOut ) :: config
    Real( Kind = wp ), Dimension( : ),                Intent( In    ) :: coeffs
    Real( Kind = wp ), Dimension( 1:neigh%max_list ), Intent( In    ) :: x_pos,y_pos,z_pos,dr_j
    Real( Kind = wp ),                                Intent(   Out ) :: engcpe_ex,vircpe_ex
    Real( Kind = wp ), Dimension( 1:9 ),              Intent( InOut ) :: stress
    Real( Kind = wp )                                                 :: atom_coeffs_i, atom_coeffs_j, atom_coeffs_ij
    Real( Kind = wp ), Parameter                                      :: a1 =  0.254829592_wp
    Real( Kind = wp ), Parameter                                      :: a2 = -0.284496736_wp
    Real( Kind = wp ), Parameter                                      :: a3 =  1.421413741_wp
    Real( Kind = wp ), Parameter                                      :: a4 = -1.453152027_wp
    Real( Kind = wp ), Parameter                                      :: a5 =  1.061405429_wp
    Real( Kind = wp ), Parameter                                      :: pp =  0.3275911_wp
    Real( Kind = wp ), Parameter                                      :: rr3  = 1.0_wp/3.0_wp
    Real( Kind = wp ), Parameter                                      :: r10  = 0.1_wp
    Real( Kind = wp ), Parameter                                      :: r42  = 1.0_wp/42.0_wp
    Real( Kind = wp ), Parameter                                      :: r216 = 1.0_wp/216.0_wp
    Real( Kind = wp )                                                 :: erfr
    Real( Kind = wp )                                                 :: exp1
    Real( Kind = wp )                                                 :: tt
    Real( Kind = wp ), dimension(9)                                   :: stress_temp
    Real( Kind = wp ), dimension(3)                                   :: pos_j
    Real( Kind = wp ), dimension(3)                                   :: force_temp, force_temp_comp
    Real( kind = wp )                                                 :: inv_mod_r_ij, mod_r_ij_2, mod_r_ij
    Real( Kind = wp )                                                 :: dr_alpha
    Real( Kind = wp )                                                 :: dr_alpha_2
    Real( Kind = wp )                                                 :: erf_gamma

    Integer                                                           :: limit
    Integer                                                           :: global_id_i
    Integer                                                           :: jatm
    Integer                                                           :: m

    ! initialise potential energy and virial

    engcpe_ex=0.0_wp
    vircpe_ex=0.0_wp

    ! initialise stress tensor accumulators

    stress_temp = 0.0_wp

    ! global identity of iatm

    global_id_i=config%ltg(iatm)

    ! ignore interaction if the charge is zero

    atom_coeffs_i = coeffs(iatm)
    if (abs(atom_coeffs_i) < zero_plus) return

    atom_coeffs_i = atom_coeffs_i*spme_datum%scaling

    ! load forces

    force_temp(1)=config%parts(iatm)%fxx
    force_temp(2)=config%parts(iatm)%fyy
    force_temp(3)=config%parts(iatm)%fzz

    ! Get neigh%list limit

    limit=neigh%list(-1,iatm)-neigh%list(0,iatm)

    ! start of primary loop for forces evaluation

    Do m=1,limit

      ! atomic index and charge

      jatm=neigh%list(neigh%list(0,iatm)+m,iatm)
      atom_coeffs_j=coeffs(jatm)

      ! interatomic distance

      mod_r_ij=dr_j(m)

      ! interaction validity and truncation of potential

      If (Abs(atom_coeffs_j) > zero_plus .and. mod_r_ij < neigh%cutoff) Then

        ! charge product
        atom_coeffs_ij=atom_coeffs_j*atom_coeffs_i
        ! Squared distance

        ! interatomic distance
        pos_j = [x_pos(m),y_pos(m),z_pos(m)]

        mod_r_ij=dr_j(m)
        inv_mod_r_ij = 1.0_wp / mod_r_ij
        mod_r_ij_2=mod_r_ij**2

        ! calculate forces

        dr_alpha =mod_r_ij*ewld%alpha

        ! calculate error function and derivative

        If (dr_alpha < 1.0e-2_wp) Then
          ! close particles (core-shell units) - small distances limit

          dr_alpha_2=dr_alpha*dr_alpha

          erfr=2.0_wp*atom_coeffs_ij*(ewld%alpha/sqrpi) * &
            (1.0_wp+dr_alpha_2*(-rr3+dr_alpha_2*(r10+dr_alpha_2*(-r42+dr_alpha_2*r216))))

          erf_gamma=-4.0_wp*atom_coeffs_ij*(ewld%alpha**3/sqrpi) * &
            (rr3+dr_alpha_2*(-2.0_wp*r10+dr_alpha_2*(3.0_wp*r42-4.0_wp*dr_alpha_2*r216)))

        Else

          ! distant particles - traditional

          erfr = atom_coeffs_ij * calc_erf(dr_alpha) * inv_mod_r_ij
          erf_gamma = ( atom_coeffs_ij * ewld%alpha*calc_erf_deriv(dr_alpha) - erfr ) * inv_mod_r_ij**2

        End If

        ! calculate forces

        force_temp_comp = erf_gamma * pos_j

        force_temp = force_temp + force_temp_comp

        If (jatm <= config%natms) Then

          config%parts(jatm)%fxx=config%parts(jatm)%fxx-force_temp_comp(1)
          config%parts(jatm)%fyy=config%parts(jatm)%fyy-force_temp_comp(2)
          config%parts(jatm)%fzz=config%parts(jatm)%fzz-force_temp_comp(3)

        End If

        If (jatm <= config%natms .or. global_id_i < config%ltg(jatm)) Then

          ! add potential energy and virial

          engcpe_ex = engcpe_ex - erfr
          vircpe_ex = vircpe_ex - erf_gamma*mod_r_ij_2

          ! calculate stress tensor
          stress_temp(1:9:3) = stress_temp(1:9:3) + pos_j*force_temp(1)
          stress_temp(2:9:3) = stress_temp(2:9:3) + pos_j*force_temp(2)
          stress_temp(3:9:3) = stress_temp(3:9:3) + pos_j*force_temp(3)

        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx=force_temp(1)
    config%parts(iatm)%fyy=force_temp(2)
    config%parts(iatm)%fzz=force_temp(3)

    ! complete stress tensor

    stress = stress + stress_temp

  End Subroutine ewald_excl_forces

  Subroutine ewald_frzn_forces(engcpe_fr,vircpe_fr,stress,ewld,neigh,electro,config,comm)

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
    Use comms, only : gsum
    implicit none
    Real( Kind = wp  ),                   Intent(   Out ) :: engcpe_fr,vircpe_fr
    Real( Kind = wp  ), Dimension( 1:9 ), Intent( InOut ) :: stress
    Type( ewald_type ),              Intent( InOut ) :: ewld
    Type( neighbours_type ),              Intent( In    ) :: neigh
    Type( electrostatic_type ), Intent( In    ) :: electro
    Type( comms_type ),                   Intent( InOut ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config

    Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
    Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
    Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
    Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
    Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
    Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

    Real( Kind = wp ), dimension(9) :: stress_temp
    Real( Kind = wp ), dimension(3) :: force_temp, force_temp_comp

    Integer           :: fail,i,j,k,ii,jj,global_id_i,nzfr,limit
    Real( Kind = wp ) :: scl,det,rcell(1:9),xrr,yrr,zrr,mod_r_ij,mod_r_ij_2, &
      atom_coeffs_ij,erfr,erf_gamma,exp1,tt

    Real( Kind = wp ) :: xss,yss,zss
    Integer,           Dimension( : ), Allocatable :: l_ind,nz_fr
    Real( Kind = wp ), Dimension( : ), Allocatable :: cfr,xfr,yfr,zfr
    Real( Kind = wp ), Dimension( : ), Allocatable :: x_pos,y_pos,z_pos,dr_j
    Character( Len = 256 ) :: message

    fail=0
    Allocate (l_ind(1:config%mxatdm),nz_fr(0:comm%mxnode), Stat=fail)
    If (fail > 0) call error_alloc('l_ind and nz_fr','ewald_frzn_forces')

    Call invert(config%cell,rcell,det)

    ! Initialise contributions

    engcpe_fr=0.0_wp
    vircpe_fr=0.0_wp

    stress_temp(1) = 0.0_wp
    stress_temp(2) = 0.0_wp
    stress_temp(3) = 0.0_wp
    stress_temp(5) = 0.0_wp
    stress_temp(6) = 0.0_wp
    stress_temp(9) = 0.0_wp

    l_ind=0 ; nz_fr=0
    Do i=1,config%natms
      If (config%lfrzn(i) > 0 .and. Abs(config%parts(i)%chge) > zero_plus) Then
        nz_fr(comm%idnode+1)=nz_fr(comm%idnode+1)+1
        l_ind(nz_fr(comm%idnode+1))=i
      End If
    End Do
    Call gsum(comm, nz_fr)
    nz_fr(0) = Sum(nz_fr(0:comm%idnode)) ! Offset

    scl=r4pie0/electro%eps
    nzfr = Sum(nz_fr(1:comm%mxnode))     ! Total
    If (nzfr <= 10*config%mxatms) Then

      Allocate (cfr(1:nzfr),xfr(1:nzfr),yfr(1:nzfr),zfr(1:nzfr), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a,i0)') 'ewald_frzn_forces allocation failure 1'
        Call error(0,message)
      End If

      cfr=0.0_wp
      xfr=0.0_wp
      yfr=0.0_wp
      zfr=0.0_wp
      Do i=1,nz_fr(comm%idnode+1)
        ii=nz_fr(0)+i

        cfr(ii)=config%parts(l_ind(i))%chge
        xfr(ii)=config%parts(l_ind(i))%xxx
        yfr(ii)=config%parts(l_ind(i))%yyy
        zfr(ii)=config%parts(l_ind(i))%zzz
      End Do
      Call gsum(comm, cfr)
      Call gsum(comm, xfr)
      Call gsum(comm, yfr)
      Call gsum(comm, zfr)

      Do i=1,nz_fr(comm%idnode+1)
        ii=nz_fr(0)+i

        Do jj=1,nz_fr(0) ! -, on nodes<comm%idnode
          xrr=xfr(ii)-xfr(jj)
          yrr=yfr(ii)-yfr(jj)
          zrr=zfr(ii)-zfr(jj)

          xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
          yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
          zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

          xss=xss-Anint(xss)
          yss=yss-Anint(yss)
          zss=zss-Anint(zss)

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          mod_r_ij_2=xrr**2+yrr**2+zrr**2

          mod_r_ij=Sqrt(mod_r_ij_2)
          atom_coeffs_ij=cfr(ii)*cfr(jj)*scl

          ! calculate error function and derivative

          exp1 =Exp(-(ewld%alpha*mod_r_ij)**2)
          tt   =1.0_wp/(1.0_wp+pp*ewld%alpha*mod_r_ij)

          erfr=atom_coeffs_ij * &
            (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/mod_r_ij

          erf_gamma=-(erfr-2.0_wp*atom_coeffs_ij*(ewld%alpha/sqrpi)*exp1)/mod_r_ij_2

          force_temp_comp(1)= erf_gamma*xrr
          force_temp_comp(2) = erf_gamma*yrr
          force_temp_comp(3) = erf_gamma*zrr

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-force_temp_comp(1)
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-force_temp_comp(2)
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-force_temp_comp(3)

          ! redundant calculations copying

          ! If (ewld%lf_cp) Then
          !    ewld%ffx(l_ind(i))=ewld%ffx(l_ind(i))-force_temp_comp(1)
          !    ewld%ffy(l_ind(i))=ewld%ffy(l_ind(i))-force_temp_comp(2)
          !    ewld%ffz(l_ind(i))=ewld%ffz(l_ind(i))-force_temp_comp(3)
          ! End If

          ! infrequent calculations copying

          ! If (ewld%l_cp) Then
          !    ewld%fcx(l_ind(i))=ewld%fcx(l_ind(i))-force_temp_comp(1)
          !    ewld%fcy(l_ind(i))=ewld%fcy(l_ind(i))-force_temp_comp(2)
          !    ewld%fcz(l_ind(i))=ewld%fcz(l_ind(i))-force_temp_comp(3)
          ! End If
        End Do

        Do j=i+1,nz_fr(comm%idnode+1) ! =, node=comm%idnode (OVERLAP but no SELF)!
          jj=nz_fr(0)+j

          xrr=xfr(ii)-xfr(jj)
          yrr=yfr(ii)-yfr(jj)
          zrr=zfr(ii)-zfr(jj)

          xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
          yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
          zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

          xss=xss-Anint(xss)
          yss=yss-Anint(yss)
          zss=zss-Anint(zss)

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          mod_r_ij_2=xrr**2+yrr**2+zrr**2

          mod_r_ij=Sqrt(mod_r_ij_2)
          atom_coeffs_ij=cfr(ii)*cfr(jj)*scl

          ! calculate error function and derivative

          exp1 =Exp(-(ewld%alpha*mod_r_ij)**2)
          tt   =1.0_wp/(1.0_wp+pp*ewld%alpha*mod_r_ij)

          erfr=atom_coeffs_ij * &
            (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/mod_r_ij

          erf_gamma=-(erfr-2.0_wp*atom_coeffs_ij*(ewld%alpha/sqrpi)*exp1)/mod_r_ij_2

          force_temp_comp(1)= erf_gamma*xrr
          force_temp_comp(2) = erf_gamma*yrr
          force_temp_comp(3) = erf_gamma*zrr

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-force_temp_comp(1)
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-force_temp_comp(2)
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-force_temp_comp(3)

          config%parts(l_ind(j))%fxx=config%parts(l_ind(j))%fxx+force_temp_comp(1)
          config%parts(l_ind(j))%fyy=config%parts(l_ind(j))%fyy+force_temp_comp(2)
          config%parts(l_ind(j))%fzz=config%parts(l_ind(j))%fzz+force_temp_comp(3)

          ! redundant calculations copying

          ! If (ewld%lf_cp) Then
          !    ewld%ffx(l_ind(i))=ewld%ffx(l_ind(i))-force_temp_comp(1)
          !    ewld%ffy(l_ind(i))=ewld%ffy(l_ind(i))-force_temp_comp(2)
          !    ewld%ffz(l_ind(i))=ewld%ffz(l_ind(i))-force_temp_comp(3)

          !    ewld%ffx(l_ind(j))=ewld%ffx(l_ind(j))+force_temp_comp(1)
          !    ewld%ffy(l_ind(j))=ewld%ffy(l_ind(j))+force_temp_comp(2)
          !    ewld%ffz(l_ind(j))=ewld%ffz(l_ind(j))+force_temp_comp(3)
          ! End If

          ! infrequent calculations copying

          ! If (ewld%l_cp) Then
          !    ewld%fcx(l_ind(i))=ewld%fcx(l_ind(i))-force_temp_comp(1)
          !    ewld%fcy(l_ind(i))=ewld%fcy(l_ind(i))-force_temp_comp(2)
          !    ewld%fcz(l_ind(i))=ewld%fcz(l_ind(i))-force_temp_comp(3)

          !    ewld%fcx(l_ind(j))=ewld%fcx(l_ind(j))+force_temp_comp(1)
          !    ewld%fcy(l_ind(j))=ewld%fcy(l_ind(j))+force_temp_comp(2)
          !    ewld%fcz(l_ind(j))=ewld%fcz(l_ind(j))+force_temp_comp(3)
          ! End If

          ! calculate potential energy and virial

          engcpe_fr = engcpe_fr - erfr
          vircpe_fr = vircpe_fr - erf_gamma*mod_r_ij_2

          ! calculate stress tensor

          stress_temp(1) = stress_temp(1) + xrr*force_temp_comp(1)
          stress_temp(2) = stress_temp(2) + xrr*force_temp_comp(2)
          stress_temp(3) = stress_temp(3) + xrr*force_temp_comp(3)
          stress_temp(5) = stress_temp(5) + yrr*force_temp_comp(2)
          stress_temp(6) = stress_temp(6) + yrr*force_temp_comp(3)
          stress_temp(9) = stress_temp(9) + zrr*force_temp_comp(3)
        End Do

        Do jj=nz_fr(0)+nz_fr(comm%idnode+1)+1,nzfr ! +, on nodes>comm%idnode
          xrr=xfr(ii)-xfr(jj)
          yrr=yfr(ii)-yfr(jj)
          zrr=zfr(ii)-zfr(jj)

          xss=(rcell(1)*xrr+rcell(4)*yrr+rcell(7)*zrr)
          yss=(rcell(2)*xrr+rcell(5)*yrr+rcell(8)*zrr)
          zss=(rcell(3)*xrr+rcell(6)*yrr+rcell(9)*zrr)

          xss=xss-Anint(xss)
          yss=yss-Anint(yss)
          zss=zss-Anint(zss)

          xrr=(config%cell(1)*xss+config%cell(4)*yss+config%cell(7)*zss)
          yrr=(config%cell(2)*xss+config%cell(5)*yss+config%cell(8)*zss)
          zrr=(config%cell(3)*xss+config%cell(6)*yss+config%cell(9)*zss)

          ! calculate interatomic distance

          mod_r_ij_2=xrr**2+yrr**2+zrr**2

          mod_r_ij=Sqrt(mod_r_ij_2)
          atom_coeffs_ij=cfr(ii)*cfr(jj)*scl

          ! calculate error function and derivative

          exp1 =Exp(-(ewld%alpha*mod_r_ij)**2)
          tt   =1.0_wp/(1.0_wp+pp*ewld%alpha*mod_r_ij)

          erfr=atom_coeffs_ij * &
            (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/mod_r_ij

          erf_gamma=-(erfr-2.0_wp*atom_coeffs_ij*(ewld%alpha/sqrpi)*exp1)/mod_r_ij_2

          force_temp_comp(1) = erf_gamma*xrr
          force_temp_comp(2) = erf_gamma*yrr
          force_temp_comp(3) = erf_gamma*zrr

          ! calculate forces

          config%parts(l_ind(i))%fxx=config%parts(l_ind(i))%fxx-force_temp_comp(1)
          config%parts(l_ind(i))%fyy=config%parts(l_ind(i))%fyy-force_temp_comp(2)
          config%parts(l_ind(i))%fzz=config%parts(l_ind(i))%fzz-force_temp_comp(3)

          ! redundant calculations copying

          ! If (ewld%lf_cp) Then
          !    ewld%ffx(l_ind(i))=ewld%ffx(l_ind(i))-force_temp_comp(1)
          !    ewld%ffy(l_ind(i))=ewld%ffy(l_ind(i))-force_temp_comp(2)
          !    ewld%ffz(l_ind(i))=ewld%ffz(l_ind(i))-force_temp_comp(3)
          ! End If

          ! infrequent calculations copying

          ! If (ewld%l_cp) Then
          !    ewld%fcx(l_ind(i))=ewld%fcx(l_ind(i))-force_temp_comp(1)
          !    ewld%fcy(l_ind(i))=ewld%fcy(l_ind(i))-force_temp_comp(2)
          !    ewld%fcz(l_ind(i))=ewld%fcz(l_ind(i))-force_temp_comp(3)
          ! End If

          ! calculate potential energy and virial

          engcpe_fr = engcpe_fr - erfr
          vircpe_fr = vircpe_fr - erf_gamma*mod_r_ij_2

          ! calculate stress tensor

          stress_temp(1) = stress_temp(1) + xrr*force_temp_comp(1)
          stress_temp(2) = stress_temp(2) + xrr*force_temp_comp(2)
          stress_temp(3) = stress_temp(3) + xrr*force_temp_comp(3)
          stress_temp(5) = stress_temp(5) + yrr*force_temp_comp(2)
          stress_temp(6) = stress_temp(6) + yrr*force_temp_comp(3)
          stress_temp(9) = stress_temp(9) + zrr*force_temp_comp(3)
        End Do
      End Do

      Deallocate (cfr,xfr,yfr,zfr, Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_forces deallocation failure 1'
        Call error(0,message)
      End If

    Else

      ! We resort to approximating N*(N-1)/2 interactions
      ! with the short-range one from the two body linked config%cell neigh%list

      Allocate (x_pos(1:neigh%max_list),y_pos(1:neigh%max_list),z_pos(1:neigh%max_list),dr_j(1:neigh%max_list), Stat=fail)
      If (fail > 0) Then
        Write(message,'(a)') 'ewald_frzn_forces allocation failure 2'
        Call error(0,message)
      End If

      Do ii=1,nz_fr(comm%idnode+1)
        i=l_ind(nz_fr(comm%idnode+1))
        global_id_i=config%ltg(ii)

        ! Get neigh%list limit

        limit=neigh%list(-2,i)-neigh%list(-1,i)
        If (limit > 0) Then

          ! calculate interatomic distances

          Do k=1,limit
            j=neigh%list(neigh%list(-1,i)+k,i)

            x_pos(k)=config%parts(i)%xxx-config%parts(j)%xxx
            y_pos(k)=config%parts(i)%yyy-config%parts(j)%yyy
            z_pos(k)=config%parts(i)%zzz-config%parts(j)%zzz
          End Do

          ! periodic boundary conditions not needed by LC construction
          !
          !           Call images(config%imcon,config%cell,limit,x_pos,y_pos,z_pos)

          ! square of distances

          Do k=1,limit
            dr_j(k)=Sqrt(x_pos(k)**2+y_pos(k)**2+z_pos(k)**2)
          End Do

          Do k=1,limit
            j=neigh%list(neigh%list(-1,i)+k,i)

            mod_r_ij=dr_j(k)
            If (Abs(config%parts(j)%chge) > zero_plus .and. mod_r_ij < neigh%cutoff) Then
              atom_coeffs_ij=config%parts(i)%chge*config%parts(j)%chge*scl
              mod_r_ij_2=mod_r_ij**2

              ! calculate error function and derivative

              exp1 =Exp(-(ewld%alpha*mod_r_ij)**2)
              tt   =1.0_wp/(1.0_wp+pp*ewld%alpha*mod_r_ij)

              erfr=atom_coeffs_ij * &
                (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/mod_r_ij

              erf_gamma=-(erfr-2.0_wp*atom_coeffs_ij*(ewld%alpha/sqrpi)*exp1)/mod_r_ij_2

              force_temp_comp(1) = erf_gamma*x_pos(k)
              force_temp_comp(2) = erf_gamma*y_pos(k)
              force_temp_comp(3) = erf_gamma*z_pos(k)

              ! calculate forces

              config%parts(i)%fxx=config%parts(i)%fxx-force_temp_comp(1)
              config%parts(i)%fyy=config%parts(i)%fyy-force_temp_comp(2)
              config%parts(i)%fzz=config%parts(i)%fzz-force_temp_comp(3)

              ! redundant calculations copying

              ! If (ewld%lf_cp) Then
              !    ewld%ffx(i)=ewld%ffx(i)-force_temp_comp(1)
              !    ewld%ffy(i)=ewld%ffy(i)-force_temp_comp(2)
              !    ewld%ffz(i)=ewld%ffz(i)-force_temp_comp(3)
              ! End If

              ! ! infrequent calculations copying

              ! If (ewld%l_cp) Then
              !    ewld%fcx(i)=ewld%fcx(i)-force_temp_comp(1)
              !    ewld%fcy(i)=ewld%fcy(i)-force_temp_comp(2)
              !    ewld%fcz(i)=ewld%fcz(i)-force_temp_comp(3)
              ! End If

              If (j <= config%natms) Then

                config%parts(j)%fxx=config%parts(j)%fxx+force_temp_comp(1)
                config%parts(j)%fyy=config%parts(j)%fyy+force_temp_comp(2)
                config%parts(j)%fzz=config%parts(j)%fzz+force_temp_comp(3)

                ! redundant calculations copying

                ! If (ewld%lf_cp) Then
                !    ewld%ffx(j)=ewld%ffx(j)+force_temp_comp(1)
                !    ewld%ffy(j)=ewld%ffy(j)+force_temp_comp(2)
                !    ewld%ffz(j)=ewld%ffz(j)+force_temp_comp(3)
                ! End If

                ! ! infrequent calculations copying

                ! If (ewld%l_cp) Then
                !    ewld%fcx(j)=ewld%fcx(j)+force_temp_comp(1)
                !    ewld%fcy(j)=ewld%fcy(j)+force_temp_comp(2)
                !    ewld%fcz(j)=ewld%fcz(j)+force_temp_comp(3)
                ! End If

              End If

              If (j <= config%natms .or. global_id_i < config%ltg(j)) Then

                ! calculate potential energy and virial

                engcpe_fr = engcpe_fr - erfr
                vircpe_fr = vircpe_fr - erf_gamma*mod_r_ij_2

                ! calculate stress tensor

                stress_temp(1) = stress_temp(1) + x_pos(k)*force_temp_comp(1)
                stress_temp(2) = stress_temp(2) + x_pos(k)*force_temp_comp(2)
                stress_temp(3) = stress_temp(3) + x_pos(k)*force_temp_comp(3)
                stress_temp(5) = stress_temp(5) + y_pos(k)*force_temp_comp(2)
                stress_temp(6) = stress_temp(6) + y_pos(k)*force_temp_comp(3)
                stress_temp(9) = stress_temp(9) + z_pos(k)*force_temp_comp(3)

              End If
            End If
          End Do

        End If
      End Do

      Deallocate (x_pos,y_pos,z_pos,dr_j, Stat=fail)
      if (fail > 0) call error_dealloc('position arrays','ewald_frzn_forces')

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

    ! redundant calculations copying

    ! If (ewld%lf_cp) Then
    !    ewld%ef_fr=engcpe_fr
    !    ewld%vf_fr=vircpe_fr

    !    ewld%sf_fr(1) = stress_temp(1)
    !    ewld%sf_fr(2) = stress_temp(2)
    !    ewld%sf_fr(3) = stress_temp(3)
    !    ewld%sf_fr(4) = stress_temp(2)
    !    ewld%sf_fr(5) = stress_temp(5)
    !    ewld%sf_fr(6) = stress_temp(6)
    !    ewld%sf_fr(7) = stress_temp(3)
    !    ewld%sf_fr(8) = stress_temp(6)
    !    ewld%sf_fr(9) = stress_temp(9)
    ! End If

    ! infrequent calculations copying

    ! If (ewld%l_cp) Then
    !    ewld%e_fr=engcpe_fr
    !    ewld%v_fr=vircpe_fr

    !    ewld%s_fr(1) = stress_temp(1)
    !    ewld%s_fr(2) = stress_temp(2)
    !    ewld%s_fr(3) = stress_temp(3)
    !    ewld%s_fr(4) = stress_temp(2)
    !    ewld%s_fr(5) = stress_temp(5)
    !    ewld%s_fr(6) = stress_temp(6)
    !    ewld%s_fr(7) = stress_temp(3)
    !    ewld%s_fr(8) = stress_temp(6)
    !    ewld%s_fr(9) = stress_temp(9)
    ! End If

    Deallocate (l_ind,nz_fr, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'ewald_frzn_forces deallocation failure'
      Call error(0,message)
    End If
  End Subroutine ewald_frzn_forces

!!! Internals

  subroutine ewald_spme_init(domain, max_atoms, comm, kspace, &
    & bspline_in, charge_grid, potential_grid, stress_grid)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to initialise the ewald SPME routines
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!

    use ewald,          only : ewald_type
    use bspline,        only : bspline_type, bspline_coeffs_gen
    use domains,        only : domains_type
    use parallel_fft,   only : initialize_fft, pfft_indices
    use comms,          only : comms_type
    use kspace,         only : setup_kspace
    implicit none


    Type( kspace_type ),                                      Intent ( inout ) :: kspace
    Type( bspline_type ),                                     Intent ( inout ) :: bspline_in
    Type( domains_type ),                                     Intent ( in    ) :: domain
    Type( comms_type ),                                       Intent ( in    ) :: comm

    Integer,                                                  Intent ( in    ) :: max_atoms
    Real( Kind = wp ),       Dimension( :,:,: ), Allocatable, Intent (   out ) :: charge_grid
    complex( Kind = wp ),    Dimension( :,:,: ), Allocatable, Intent (   out ) :: potential_grid
    complex( Kind = wp ),    Dimension( :,:,: ), Allocatable, Intent (   out ) :: stress_grid
    Integer, dimension(4) :: fail
    
    fail = 0

    print*, "hi"
    call setup_kspace(kspace, domain, kspace%k_vec_dim)
    
!!! begin cardinal b-splines set-up

    bspline_in%num_deriv = 2
    allocate(bspline_in%derivs(3,0:bspline_in%num_deriv,1:bspline_in%num_splines,1:max_atoms), stat=fail(1))
    if (fail(1) > 0) call error_alloc('bspline_in%derivs','ewald_spme_init')
    
    ! calculate the global b-spline coefficients
    call bspline_coeffs_gen(kspace, bspline_in)

!!! end cardinal b-splines set-up

!!! begin daft set-up

    ! set up the parallel fft and useful related quantities

    Call initialize_fft( 3, kspace%k_vec_dim, &
      [ domain%nx, domain%ny, domain%nz ], [ domain%idx, domain%idy, domain%idz ],   &
      [ kspace%block_x, kspace%block_y, kspace%block_z ],               &
      comm%comm, kspace%context )

    Call pfft_indices( kspace%k_vec_dim(1), kspace%block_x, domain%idx, domain%nx, kspace%index_x )
    Call pfft_indices( kspace%k_vec_dim(2), kspace%block_y, domain%idy, domain%ny, kspace%index_y )
    Call pfft_indices( kspace%k_vec_dim(3), kspace%block_z, domain%idz, domain%nz, kspace%index_z )

    ! workspace arrays for DaFT

    allocate ( charge_grid   ( 1:kspace%block_x, 1:kspace%block_y, 1:kspace%block_z ), stat = fail(1) )
    allocate ( potential_grid( 1:kspace%block_x, 1:kspace%block_y, 1:kspace%block_z ), stat = fail(2) )
    allocate ( stress_grid   ( 1:kspace%block_x, 1:kspace%block_y, 1:kspace%block_z ), stat = fail(3) )
    allocate ( pfft_work     ( 1:kspace%block_x, 1:kspace%block_y, 1:kspace%block_z ), stat = fail(4) )
    if (any(fail > 0)) call error_alloc('SPME DaFT workspace arrays','ewald_spme_init')

!!! end daft set-up

  end subroutine ewald_spme_init

  subroutine spme_construct_charge_array(ncalc,ewld,lookup_array,recip_indices,electro, &
    & coeffs,charge_grid)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to construct the charge array for SPME calculations
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins & i.t.todorov & i.j.bush august 2018
    !!
    !!----------------------------------------------------------------------!
    implicit none

    Type( ewald_type ),                Intent( in    ) :: ewld
    Type( electrostatic_type ),             Intent( in    ) :: electro
    Real( Kind = wp ), Dimension(:),        Intent( in    ) :: coeffs
    Integer,           Dimension(:,:),      Intent( in    ) :: recip_indices
    Integer,           Dimension(:),        Intent( in    ) :: lookup_array
    Integer,                                Intent( in    ) :: ncalc
    Real( Kind = wp ), Dimension(:,:,:),    Intent(   out ) :: charge_grid

    Real( Kind = wp ) :: atom_coeffs
    Integer, Dimension(3,2) :: spline_bounds
    Integer, Dimension(3) :: temp
    Integer, Dimension(3) :: current_derivs
    Integer :: j_lo,j_hi,k_lo,k_hi,l_lo,l_hi
    Integer :: j,k,l,jj,kk,ll
    Integer :: atm, i

    charge_grid = 0.0_wp

    ! construct 3d charge array
    ! daft version - use array that holds only the local data

    atom:do atm=1,ncalc

      i = lookup_array(atm)
      ! if a particle is charged and in the md cell or in its positive halo
      ! (t(i) >= 0) as the b-splines are negative directionally by propagation

      j_lo = max(1, recip_indices(1,i) - ewld%kspace%domain_indices(1,1) - ewld%bspline%num_splines + 3)
      k_lo = max(1, recip_indices(2,i) - ewld%kspace%domain_indices(2,1) - ewld%bspline%num_splines + 3)
      l_lo = max(1, recip_indices(3,i) - ewld%kspace%domain_indices(3,1) - ewld%bspline%num_splines + 3)
      j_hi = min(ewld%kspace%domain_indices(1,2), recip_indices(1,i) + 1) - ewld%kspace%domain_indices(1,1) + 1
      k_hi = min(ewld%kspace%domain_indices(2,2), recip_indices(2,i) + 1) - ewld%kspace%domain_indices(2,1) + 1
      l_hi = min(ewld%kspace%domain_indices(3,2), recip_indices(3,i) + 1) - ewld%kspace%domain_indices(3,1) + 1

      temp = recip_indices(:,i) - ewld%bspline%num_splines - ewld%kspace%domain_indices(:,1) + 2

      atom_coeffs = coeffs(i)

      do l = l_lo, l_hi
        do k = k_lo, k_hi
          do j = j_lo, j_hi
            charge_grid(j,k,l) = charge_grid(j,k,l) + atom_coeffs * &
              & ewld%bspline%derivs(1,0,j-temp(1),i)* &
              & ewld%bspline%derivs(2,0,k-temp(2),i)* &
              & ewld%bspline%derivs(3,0,l-temp(3),i)
          end do
        end do
      end do

    end do atom

  end subroutine spme_construct_charge_array

  subroutine spme_construct_potential_grid(ewld, recip_cell, charge_grid, spme_datum, &
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

    use comms, only : gsum
    use constants, only : twopi, pi
    use parallel_fft, only : pfft
    implicit none

    type( ewald_type ),                Intent ( in    )           :: ewld
    type( spme_component ),                 Intent ( in    )           :: spme_datum
    Real( Kind = wp ),    Dimension(3,3),   Intent (   out ), optional :: stress_contrib        !! SPME contribution to the stress
    complex( Kind = wp ), Dimension(:,:,:), Intent (   out )           :: potential_grid
    Real( Kind = wp ),    Dimension(:,:,:), Intent ( in    )           :: charge_grid
    Real( Kind = wp ),    Dimension(9),     Intent ( in    )           :: recip_cell            !! Reciprocal lattice vectors
    complex ( Kind = wp ),    external                                 :: kernel                !! Core function to FT
    Real( Kind = wp ),    Dimension(10)                                :: recip_cell_properties !! bbb(1 to 3) - lengths of cell vectors: a(x,y,z) , b(x,y,z) , c(x,y,z)
    !! bbb(4 to 6) - cosines of cell angles: gamma(a,b) , beta(a,c) , alpha(b,c)
    !! bbb(7 to 9) - perpendicular cell widths : wx(y,z) , wy(x,z) , wz(x,y)
    !! bbb(10)     - cell volume
    Real( Kind = wp )                                                  :: cut_off               !! Reciprocal space cut-off
    Real( Kind = wp )                                                  :: cut_off_2             !! Reciprocal space cut-off squared
    Real( Kind = wp ),    Dimension(3,3)                               :: stress_temp
    Real( Kind = wp )                                                  :: bb1,bb2,bb3
    Real( Kind = wp ),    Dimension(3,3)                               :: recip_pos             !! Current Step in K vectors
    Real( Kind = wp )                                                  :: k_vec_2               !! Magnitude of RK vectors
    Real( Kind = wp )                                                  :: recip_conv_fac        !! Reciprocal convergence factor
    Real( Kind = wp )                                                  :: f_p_fac               !! pi*m/beta
    complex( Kind = wp )                                               :: potential_component   !! Contribution to the potential
    Real( Kind = wp )                                                  :: pressure_virial       !! Virial contribution to the pressure
    Integer                                                            :: j,k,l, jj,kk,ll, j_local, k_local, l_local
    Integer                                                            :: alpha,beta

    recip_conv_fac = 1.0_wp/ewld%alpha

    ! set reciprocal space cutoff

    call dcell(recip_cell,recip_cell_properties)

    cut_off=0.5_wp*1.05_wp*minval(ewld%kspace%k_vec_dim_real*recip_cell_properties(7:9))
    cut_off_2=cut_off**2

    ! load charge array into complex array for fft

    potential_grid=cmplx(charge_grid , Kind = wp)

    ! calculate inverse 3d fft of charge array (in place)

    call pfft(potential_grid,pfft_work,ewld%kspace%context,1)

    ! initialise temporary stress tensor

    stress_temp = 0.0_wp

    ! calculate convolution of charge array with gaussian function
    ! daft version - only loop over the local stuff

    do l_local=1,ewld%kspace%block_z
      l=ewld%kspace%index_z(l_local)

      ll = l-1
      if (2*ll > ewld%kspace%k_vec_dim(3)) ll = ll - ewld%kspace%k_vec_dim(3)

      recip_pos(:,3) = Real(ll,wp)*recip_cell(3:9:3)
      bb3 = ewld%bspline%norm2(3,l)

      do k_local=1,ewld%kspace%block_y
        k=ewld%kspace%index_y(k_local)

        kk = k-1
        if (2*kk > ewld%kspace%k_vec_dim(2)) kk = kk - ewld%kspace%k_vec_dim(2)

        recip_pos(:,2) = recip_pos(:,3) + Real(kk,wp)*recip_cell(2:9:3)
        bb2 = bb3*ewld%bspline%norm2(2,k)

        do j_local=1,ewld%kspace%block_x
          j=ewld%kspace%index_x(j_local)

          jj = j-1                                                               
          if (2*jj > ewld%kspace%k_vec_dim(1)) jj = jj - ewld%kspace%k_vec_dim(1)

          recip_pos(:,1) = recip_pos(:,2) + Real(jj,wp)*recip_cell(1:9:3)
          bb1 = bb2 * ewld%bspline%norm2(1,j)

          k_vec_2 = dot_product(recip_pos(:,1),recip_pos(:,1))

          if ( k_vec_2 <= cut_off_2 ) then

            f_p_fac = pi*sqrt(k_vec_2)*recip_conv_fac

            potential_component = kernel(bb1, potential_grid(j_local,k_local,l_local), &
              & f_p_fac, ewld%alpha, spme_datum%pot_order)

            ! By L'Hopital's rule, m=0 does not contribute to stress
            if (present(stress_contrib) .and. k_vec_2 > 1.0e-6_wp) then

              pressure_virial = Real( stress_kernel(bb1, potential_grid(j_local,k_local,l_local), &
                & f_p_fac, ewld%alpha, spme_datum%pot_order)  &
                & * conjg(potential_grid(j_local,k_local,l_local)),wp ) * f_p_fac/k_vec_2
              do alpha = 1,3
                do beta = 1,3
                  stress_temp(beta,alpha) = stress_temp(beta,alpha) + &
                    & recip_pos(alpha,1)*recip_pos(beta,1)*pressure_virial
                end do
              end do
            end if
          else
            potential_component = (0.0_wp, 0.0_wp)
          end if

          potential_grid(j_local,k_local,l_local) = potential_component

        end do
      end do
    end do

    if (present(stress_contrib)) stress_contrib = stress_temp

    call pfft(potential_grid,pfft_work,ewld%kspace%context,-1)

  end subroutine spme_construct_potential_grid

  subroutine spme_calc_force_energy(ewld, electro, comm, domain, config, coeffs, recip_cell, &
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
    use comms,  only : gsum
    use constants,  only : twopi
    use domains, only : exchange_grid

    implicit none

    type ( ewald_type ),                   Intent ( in    ) :: ewld
    type ( electrostatic_type ),                Intent ( in    ) :: electro
    type ( comms_type ),                        Intent ( inout ) :: comm
    type ( domains_type ),                      Intent ( in    ) :: domain
    type ( configuration_type ),                Intent ( in    ) :: config

    Real( Kind = wp ),    Dimension( 0: ),     Intent (   out ) :: energies                !! List of per-particle energy contributions
    Real( Kind = wp ),    Dimension( :,: ),    Intent (   out ) :: forces                  !! Per-particle forces
    complex( Kind = wp ), Dimension( :,:,: ),  Intent ( in    ) :: potential_grid          !! Grid containing back FT'd potential
    Real( Kind = wp ),    Dimension( : ),      Intent ( in    ) :: coeffs                  !! Coefficients such as charges or potential
    Real( Kind = wp ),    Dimension( 9 ),      Intent ( in    ) :: recip_cell              !! Reciprocal lattice vectors
    Integer,              Dimension( :,: ),    Intent ( in    ) :: recip_indices           !! Reciprocal grid locations of charge centres
    logical,                                   Intent ( in    ) :: per_part_step           !! Whether to perform per-particle measurements
    Real( Kind = wp )                                           :: atom_coeffs
    Integer :: l_mp, curr, l_curr

    Real( Kind = wp ),    Dimension( :,:,: ), allocatable, save :: extended_potential_grid !! Grid with extended halo splines
    Real( Kind = wp ),    Dimension( 1:3 )                      :: energy_temp             !! Temporary energy components
    Real( Kind = wp )                                           :: energy_total            !! Total sum of atomic energy
    Real( Kind = wp ),    Dimension( 3,1:3 )                    :: force_temp              !! Temporary matrix used in force calculation
    Real( Kind = wp ),    Dimension( 3 )                        :: force_total             !! Total sum of all forces
    Real( Kind = wp ),    Dimension( 3 )                        :: curr_force_temp         !! Temporary force vec
    Real( Kind = wp ),    Dimension( 3,3 )                      :: recip_cell_mat          !! In matrix form
    Real( Kind = wp ),    Dimension( 3 )                        :: recip_kmax

    Integer, Dimension( 3, 2 ), save :: extended_domain                                   !! Size of extended grid with halo splines
    Integer, Dimension( 3, 2 )       :: spline_bounds

    Integer :: i,j,k,l,jj,kk,ll

    Integer :: fail
    Integer, save :: mxspl2_old = -1

    recip_cell_mat = reshape(recip_cell,[3,3])
    recip_kmax = matmul(recip_cell_mat, ewld%kspace%k_vec_dim_real)

    ! Exchange grid
    if (ewld%bspline%num_spline_padded .ne. mxspl2_old) then
      mxspl2_old = ewld%bspline%num_spline_padded
      extended_domain(:,1) = ewld%kspace%domain_indices(:,1) - ewld%bspline%num_spline_padded
      extended_domain(:,2) = ewld%kspace%domain_indices(:,2) + ewld%bspline%num_spline_padded &
        & - ewld%bspline%num_splines

      if (allocated(extended_potential_grid)) then
        deallocate (extended_potential_grid, stat=fail)
        if (fail > 0) call error_dealloc('extended_potential_grid','spme_calc_force_energy')
      end if

      allocate (extended_potential_grid( &
        & extended_domain(1,1):extended_domain(1,2), &
        & extended_domain(2,1):extended_domain(2,2), &
        & extended_domain(3,1):extended_domain(3,2) ), stat=fail)
      if (fail > 0) call error_alloc('extended_potential_grid','spme_calc_force_energy')
    end if

    call exchange_grid( &
      & ewld%kspace%domain_indices(1,1) , ewld%kspace%domain_indices(1,2) , &
      & ewld%kspace%domain_indices(2,1) , ewld%kspace%domain_indices(2,2) , &
      & ewld%kspace%domain_indices(3,1) , ewld%kspace%domain_indices(3,2) , Real(potential_grid, wp) , &
      & extended_domain(1,1), extended_domain(2,1), extended_domain(3,1), &
      & extended_domain(1,2), extended_domain(2,2), extended_domain(3,2), extended_potential_grid, domain, comm )

    ! Zero accumulators, energies and forces
    energies(0) = 0.0_wp
    forces   = 0.0_wp
    force_total = 0.0_wp

    ! Calculate per-particle contributions

    atom:do i = 1, config%natms

      energy_total = 0.0_wp
      curr_force_temp = 0.0_wp
      atom_coeffs = coeffs(i)

      do l = 1, ewld%bspline%num_splines
        ll = recip_indices(3,i) + 1 - ewld%bspline%num_splines + l

        energy_temp(3)  = atom_coeffs * ewld%bspline%derivs(3,0,l,i)

        force_temp(1,3) = atom_coeffs * ewld%bspline%derivs(3,0,l,i)
        force_temp(2,3) = atom_coeffs * ewld%bspline%derivs(3,0,l,i)
        force_temp(3,3) = atom_coeffs * ewld%bspline%derivs(3,1,l,i)

        do k = 1, ewld%bspline%num_splines
          kk = recip_indices(2,i) + 1  - ewld%bspline%num_splines + k

          energy_temp(2)  = energy_temp(3) * ewld%bspline%derivs(2,0,k,i)

          ! force_temp(:,2) = force_temp(:,3) * ewld%bspline%derivs(2,current_derivs,k,i)
          force_temp(1,2) = force_temp(1,3) * ewld%bspline%derivs(2,0,k,i)
          force_temp(2,2) = force_temp(2,3) * ewld%bspline%derivs(2,1,k,i)
          force_temp(3,2) = force_temp(3,3) * ewld%bspline%derivs(2,0,k,i)

          do j = 1, ewld%bspline%num_splines
            jj = recip_indices(1,i) + 1  - ewld%bspline%num_splines + j

            ! force_temp(:,1) = force_temp(:,2) * ewld%bspline%derivs(1,current_derivs,j,i) * &
            !   & extended_potential_grid(jj,kk,ll) * recip_kmax
            force_temp(1,1) = force_temp(1,2) * ewld%bspline%derivs(1,1,j,i) * & 
              & extended_potential_grid(jj,kk,ll) * recip_kmax(1)
            force_temp(2,1) = force_temp(2,2) * ewld%bspline%derivs(1,0,j,i) * & 
              & extended_potential_grid(jj,kk,ll) * recip_kmax(2)
            force_temp(3,1) = force_temp(3,2) * ewld%bspline%derivs(1,0,j,i) * & 
              & extended_potential_grid(jj,kk,ll) * recip_kmax(3)

            ! Sum force contributions
            force_total = force_total - force_temp(:,1)
            curr_force_temp = curr_force_temp + force_temp(:,1)
            ! energy_total now holds omega_j * 2piV
            energy_total = energy_total + &
              & energy_temp(2) * ewld%bspline%derivs(1,0,j,i) * extended_potential_grid(jj,kk,ll) !energy_temp(1)

          end do
        end do
      end do

      ! Add to total accumulators
      energies(0) = energies(0) + energy_total
      if (per_part_step) energies(i) = energy_total
      forces(:,i) = forces(:,i) - curr_force_temp

    end do atom

    ! Correct for CoM term
    call gsum(comm, force_total)

    force_total = force_total / Real(config%megatm, wp)
    ! Remove CoM
    do i = 1, config%natms
      forces(:,i) = (forces(:,i) - force_total)
    end do

  end subroutine spme_calc_force_energy

  subroutine spme_calc_stress(ewld, electro, comm, domain, config, coeffs, &
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
    use comms,  only : gsum
    use constants,  only : twopi
    use domains, only : exchange_grid

    implicit none

    type( ewald_type ),                   Intent ( in    ) :: ewld
    type( electrostatic_type ),                Intent ( in    ) :: electro
    type( comms_type ),                        Intent ( inout ) :: comm
    type( domains_type ),                      Intent ( in    ) :: domain
    type( configuration_type ),                Intent ( in    ) :: config

    Real( Kind = wp ),    Dimension( :,:,0: ), Intent (   out ) :: stress_out              !! Output stress
    complex( Kind = wp ), Dimension( :,:,: ),  Intent ( in    ) :: stress_grid             !! Grid containing back FT'd stress_contrib
    Real( Kind = wp ),    Dimension( : ),      Intent ( in    ) :: coeffs                  !! Coefficients such as charges or potentials
    Real( Kind = wp ),    Dimension( 9 ),      Intent ( in    ) :: recip_cell              !! Reciprocal lattice vectors
    Integer,              Dimension( :,: ),    Intent ( in    ) :: recip_indices           !! Reciprocal grid locations of charge centres
    Real( Kind = wp ) :: atom_coeffs
    Integer :: l_mp, curr, l_curr

    Real( Kind = wp ),    Dimension( :,:,: ), allocatable, save  :: extended_stress_grid    !! Grid with extended halo splines
    Real( Kind = wp ),    Dimension( 3,3,0:3 )                   :: stress_temp             !! Diag, off-diag
    Real( Kind = wp ),    Dimension( 3,3 )                       :: recip_cell_mat          !! In matrix form
    Real( Kind = wp ),    Dimension( 3 )                         :: recip_kmax

    Integer, Dimension( 3, 2 ), save :: extended_domain                                   !! Size of extended grid with halo splines
    Integer, Dimension( 3, 2 ) :: spline_bounds

    Integer :: i,j,k,l,jj,kk,ll
    Integer :: alpha, beta
    Integer, save :: mxspl2_old = -1
    Integer :: fail

    recip_cell_mat = reshape(recip_cell,[3,3])
    recip_kmax = matmul(recip_cell_mat, ewld%kspace%k_vec_dim_real)

    ! Exchange grid
    if (ewld%bspline%num_spline_padded .ne. mxspl2_old) then
      mxspl2_old = ewld%bspline%num_spline_padded
      extended_domain(:,1) = ewld%kspace%domain_indices(:,1) - ewld%bspline%num_spline_padded
      extended_domain(:,2) = ewld%kspace%domain_indices(:,2) + ewld%bspline%num_spline_padded - ewld%bspline%num_splines

      if (allocated(extended_stress_grid)) then
        deallocate (extended_stress_grid, stat=fail)
        if (fail > 0) call error_dealloc('extended_stress_grid','spme_calc_stress')
      end if
      allocate (extended_stress_grid( &
        & extended_domain(1,1):extended_domain(1,2), &
        & extended_domain(2,1):extended_domain(2,2), &
        & extended_domain(3,1):extended_domain(3,2) ), stat=fail)
      if (fail > 0) call error_alloc('extended_stress_grid','spme_calc_stress')
    end if

    call exchange_grid( &
      & ewld%kspace%domain_indices(1,1) , ewld%kspace%domain_indices(1,2) , &
      & ewld%kspace%domain_indices(2,1) , ewld%kspace%domain_indices(2,2) , &
      & ewld%kspace%domain_indices(3,1) , ewld%kspace%domain_indices(3,2) , Real(stress_grid, wp) , &
      & extended_domain(1,1), extended_domain(2,1), extended_domain(3,1), &
      & extended_domain(1,2), extended_domain(2,2), extended_domain(3,2), extended_stress_grid, domain, comm  )

    ! Zero accumulator
    stress_out(:,:,0) = 0.0_wp

    ! Calculate per-particle contributions

    atom:do i = 1, config%natms

      stress_temp = 0.0_wp
      atom_coeffs = coeffs(i)

      spline_bounds(:,1) = recip_indices(:,i) - ewld%bspline%num_splines + 2
      spline_bounds(:,2) = recip_indices(:,i) + 1

      do l = 1, ewld%bspline%num_splines
        ll = recip_indices(3,i) + 1 - ewld%bspline%num_splines + l

        stress_temp(1,1,3) = atom_coeffs * ewld%bspline%derivs(3,0,l,i)
        stress_temp(2,1,3) = atom_coeffs * ewld%bspline%derivs(3,0,l,i)
        stress_temp(3,1,3) = atom_coeffs * ewld%bspline%derivs(3,1,l,i)
        stress_temp(2,2,3) = atom_coeffs * ewld%bspline%derivs(3,0,l,i)
        stress_temp(3,2,3) = atom_coeffs * ewld%bspline%derivs(3,1,l,i)
        stress_temp(3,3,3) = atom_coeffs * ewld%bspline%derivs(3,2,l,i)

        do k = 1, ewld%bspline%num_splines
          kk = recip_indices(2,i) + 1 - ewld%bspline%num_splines + k

          stress_temp(1,1,2) = stress_temp(1,1,3) * ewld%bspline%derivs(2,0,k,i)
          stress_temp(2,1,2) = stress_temp(2,1,3) * ewld%bspline%derivs(2,1,k,i)
          stress_temp(3,1,2) = stress_temp(3,1,3) * ewld%bspline%derivs(2,0,k,i)
          stress_temp(2,2,2) = stress_temp(2,2,3) * ewld%bspline%derivs(2,2,k,i)
          stress_temp(3,2,2) = stress_temp(3,2,3) * ewld%bspline%derivs(2,1,k,i)
          stress_temp(3,3,2) = stress_temp(3,3,3) * ewld%bspline%derivs(2,0,k,i)

          do j = 1, ewld%bspline%num_splines
            jj = recip_indices(1,i) + 1 - ewld%bspline%num_splines + j

            stress_temp(1,1,1) = stress_temp(1,1,2) * ewld%bspline%derivs(1,2,j,i) * extended_stress_grid(jj,kk,ll)
            stress_temp(2,1,1) = stress_temp(2,1,2) * ewld%bspline%derivs(1,1,j,i) * extended_stress_grid(jj,kk,ll)
            stress_temp(3,1,1) = stress_temp(3,1,2) * ewld%bspline%derivs(1,1,j,i) * extended_stress_grid(jj,kk,ll)
            stress_temp(2,2,1) = stress_temp(2,2,2) * ewld%bspline%derivs(1,0,j,i) * extended_stress_grid(jj,kk,ll)
            stress_temp(3,2,1) = stress_temp(3,2,2) * ewld%bspline%derivs(1,0,j,i) * extended_stress_grid(jj,kk,ll)
            stress_temp(3,3,1) = stress_temp(3,3,2) * ewld%bspline%derivs(1,0,j,i) * extended_stress_grid(jj,kk,ll)

            do beta = 1,3
              do alpha = beta,3
                stress_temp(alpha,beta,0) = stress_temp(alpha,beta,0) - stress_temp(alpha,beta,1) &
                  & * recip_kmax(alpha) * recip_kmax(beta)
              end do
            end do

            stress_temp(1,2,0) = stress_temp(2,1,0)
            stress_temp(1,3,0) = stress_temp(3,1,0)
            stress_temp(2,3,0) = stress_temp(3,2,0)

          end do
        end do
      end do

      stress_out(:,:,0) = stress_out(:,:,0) + stress_temp(:,:,0)
      stress_out(:,:,i) = stress_temp(:,:,0)

    end do atom

  end subroutine spme_calc_stress

!!! Kernels

  function potential_kernel(B_m, pot, pi_m_over_a, conv_factor, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Kernel for calculating energy and forces for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    ! use spme,      only : f_p
    implicit none
    complex ( Kind = wp ) :: potential_kernel
    complex ( Kind = wp ) :: pot
    Real ( Kind = wp ) :: B_m
    Real ( Kind = wp ) :: pi_m_over_a
    Real ( Kind = wp ) :: conv_factor
    Integer :: pot_order

    potential_kernel = B_m * pot * f_p(pi_m_over_a, pot_order)

  end function potential_kernel

  function stress_kernel(B_m, pot, pi_m_over_a, conv_factor, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Kernel for calculating energy and forces for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    ! use spme,      only : f_p, f_p_d
    use constants, only : pi
    implicit none
    Complex ( Kind = wp ) :: stress_kernel

    Real ( Kind = wp ) :: B_m
    Complex ( Kind = wp ) :: pot
    Real ( Kind = wp ) :: pi_m_over_a
    Real ( Kind = wp ) :: conv_factor
    Integer :: pot_order

    Real ( Kind = wp ) :: energy
    Real ( Kind = wp ) :: x

    energy = f_p(pi_m_over_a,pot_order)
    stress_kernel = B_m * pot * f_p_d(pi_m_over_a, energy, pot_order) !* pi / (sqrt(mod_kvec_2)*conv_factor)

  end function stress_kernel

  subroutine write_per_part_contribs(config, comm, energies, forces, stresses, nstep, pot_ref)
    !!----------------------------------------------------------------------!
    !!
    !! Write out per-particle contributions to energy, force, stress, etc
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use mpi,       only : mpi_offset_Kind, mpi_mode_wronly, mpi_info_null, mpi_mode_create, mpi_comm_self
    use io, only : io_type, io_get_parameters, io_set_parameters, io_init, io_open, io_close, &
      & io_finalize, io_write_sorted_file, io_base_comm_not_set, io_allocation_error, &
      & io_unknown_write_option, io_unknown_write_level, io_write_sorted_mpiio, io_delete, io_write_batch
    use io, only : io_histord, io_restart, io_history
    use comms, only : gsync, gsum

    implicit none
    Type( configuration_type ),           Intent ( in    )  :: config    !! Atom details
    Type( comms_type ),                   Intent ( inout )  :: comm      !! Communicator
    Real( Kind = wp ), Dimension(0:),     Intent ( in    )  :: energies  !! Per-particle energies
    Real( Kind = wp ), Dimension(:,1:),   Intent ( in    )  :: forces    !!     ""       forces
    Real( Kind = wp ), Dimension(:,:,0:), Intent ( in    )  :: stresses  !!     ""       stresses
    Integer,                              Intent ( in    )  :: nstep     !! Steps since calculation start
    Integer,                              Intent ( in    )  :: pot_ref


    Type( io_type )  :: my_io                              !! Use our own IO job for now because passing through will be hell

    Real( Kind = wp ), Dimension(:), allocatable :: dummy !! Don't like this, but quick cheat?

    Integer, Parameter                                     :: record_size = 73 !! default record size (apparently)
    Integer( Kind = mpi_offset_Kind )                      :: rec_mpi_io
    Integer                                                :: energy_force_handle !! File handles
    Integer                                                :: io_write !! Write state
    Integer                                                :: batsz
    character(len=record_size)                             :: record
    character, Dimension(record_size,10)                   :: buffer
    character                                              :: lf
    character( len = 40 )                                  :: filename
    Integer                                                :: i, jj
    Integer                                                :: ierr

    call gsync(comm)

    ! Force MPIIO write for now
    io_write = 0
    ! Call io_get_parameters( user_method_write      = io_write )
    Call io_get_parameters( my_io, user_buffer_size_write = batsz, user_line_feed = lf )

    ! Write current time-step to character string
    allocate(dummy(config%natms), stat=ierr)
    if (ierr .ne. 0) call error_alloc('dummy','write_per_part_contribs')
    dummy = 0.0_wp

    write(filename,'("PPCONT",2("_",i0))') pot_ref, nstep

    call io_init( my_io, record_size )

    rec_mpi_io = int(0,mpi_offset_Kind)
    jj=0
    if (comm%idnode == 0) then

      call io_set_parameters( my_io, user_comm = mpi_comm_self )
      call io_delete( my_io, filename, comm ) ! sort existence issues
      call io_open( my_io, io_write, mpi_comm_self, trim(filename), mpi_mode_wronly + mpi_mode_create, energy_force_handle )

      jj=jj+1
      Write(record, Fmt='(a72,a1)') "Energy and force contributions on a per-particle basis",lf
      buffer(:,jj) = [(record(i:i),i=1,record_size)]
      Write(record, Fmt='(a72,a1)') config%cfgname(1:72),lf
      buffer(:,jj) = [(record(i:i),i=1,record_size)]
      jj=jj+1
      Write(record, Fmt='(3i10,42X,a1)') config%imcon,config%megatm,nstep,lf
      buffer(:,jj) = [(record(i:i),i=1,record_size)]

      If (config%imcon > 0) Then
        Do i = 0, 2
          jj=jj+1
          Write(record, Fmt='(3f20.10,a12,a1)') &
            config%cell( 1 + i * 3 ), config%cell( 2 + i * 3 ), config%cell( 3 + i * 3 ), Repeat( ' ', 12 ), lf
          buffer(:,jj) = [(record(i:i),i=1,record_size)]
        End Do
      End If

      call io_write_batch( my_io, energy_force_handle, rec_mpi_io, jj, buffer )

      Call io_close( my_io, energy_force_handle )

    end if

    call gsync(comm)

    call io_set_parameters( my_io, user_comm = comm%comm )
    call io_open( my_io, io_write, comm%comm, trim(filename), mpi_mode_wronly, energy_force_handle ) ! Io sorted mpiio, per-particle contrib

    rec_mpi_io = int(jj,mpi_offset_Kind)
    ! Only write E&F (r/v in write_sorted...) hence 1
    ! Need to skip 0th element (accumulator/total)
    call io_write_sorted_file( my_io, energy_force_handle, 2, io_history, rec_mpi_io, config%natms,      &
      config%ltg, config%atmnam, dummy, dummy, energies(1:config%natms), &
      forces(1,1:config%natms), forces(2,1:config%natms), forces(3,1:config%natms), &
      & stresses(1,1,1:config%natms), stresses(2,2,1:config%natms), stresses(3,3,1:config%natms), &
      & stresses(1,2,1:config%natms), stresses(1,3,1:config%natms), stresses(2,3,1:config%natms), ierr)

    if ( ierr /= 0 ) then
      select case( ierr )
      case( io_base_comm_not_set )
        call error( 1050 )
      case( io_allocation_error )
        call error( 1053 )
      case( io_unknown_write_option )
        call error( 1056 )
      case( io_unknown_write_level )
        call error( 1059 )
      end select
    end if
    call io_close( my_io, energy_force_handle )

    call gsync(comm)

    call io_finalize(my_io)

    deallocate(dummy, stat=ierr)
    if ( ierr > 0 ) call error_dealloc('dummy','write_per_part_contribs')

  end subroutine write_per_part_contribs

  function f_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Nth order f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins november 2018
    !! based on  - i.j.bush igf.f90 november 2018
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi
    use numerics,  only : calc_erfc, calc_exp_int, calc_inv_gamma_1_2
    use spme,      only : f_1, f_2, f_4, f_6, f_12
    implicit none
    real ( kind = wp ), intent ( in     ) :: x
    integer, intent ( in    ) :: pot_order
    integer :: curr_pot_order, p_work
    real ( kind = wp ) :: base_integ
    real ( kind = wp ) :: x_fac, x_2
    real ( kind = wp ) :: exp_xsq
    real ( kind = wp ) :: xp, curr_xp
    real ( kind = wp ) :: f_p

    select case (pot_order)
    case (1)
      if ( x > 1.0e-6_wp ) then
        f_p = f_1(x)
      else
        f_p = 0.0_wp
      end if
      
    case (2)
      if ( x > 1.0e-6_wp ) then
        f_p = f_2(x)
      else
        f_p = 0.0_wp
      end if
     case (4)
       f_p = f_4(x)
     case (6)
       f_p = f_6(x)
     case (12)
       f_p = f_12(x)
     case (:0)
       call error(0,'Invalid pot order in f_p')
    case default
      
      x_2 = x**2
      exp_xsq = exp( -x_2 )
      p_work = 2 - pot_order

      if( mod( p_work, 2 ) == 0 ) then
        ! even integrals base is I( 0, x )
        base_integ = 0.5_wp * sqrpi * calc_erfc( x )
        curr_pot_order = 0
        xp = 1.0_wp
      else
        if( p_work > 0 ) then
          ! positive odd integrals base is I( 1, x )
          base_integ = 0.5_wp * exp_xsq
          curr_pot_order = 1
          xp = x
        else
          ! negative odd integrals, base is I( -1, x ), which is 0.5 * E1( x * x )
          ! where e1 is the first order exponential integral
          base_integ = - 0.5_wp * calc_exp_int( -x_2 )
          curr_pot_order = -1
          xp = 1.0_wp / x
        end if
      end if

      f_p = base_integ

      if ( curr_pot_order == p_work ) then
        if ( x < 1.0e-6_wp ) then
          f_p = 0.0_wp ! if p < 3 && x is small
          return
        end if

        continue
      else if( curr_pot_order > p_work ) then
        ! recurse down
        x_fac = 1.0_wp / x_2
        curr_xp = xp / x
        do curr_pot_order = curr_pot_order, p_work+1, -2
          ! f_p = 2.0_wp * ( f_p - 0.5_wp * x**(curr_pot_order - 1) * exp_xsq ) / real( curr_pot_order - 1, wp )
          f_p = 2.0_wp * ( f_p - 0.5_wp * curr_xp * exp_xsq ) / real( curr_pot_order - 1, wp )
          curr_xp = curr_xp * x_fac
        end do
      else
        ! recurse up
        x_fac = x_2
        curr_xp = xp * x
        do curr_pot_order = curr_pot_order, p_work-1, 2
          ! not tested !!!!! 5/11/18
          ! f_p = 0.5_wp * x ** ( p_now + 1 ) ) * exp_xsq + 0.5_wp * ( curr_pot_order + 1 ) * f_p
          f_p = 0.5_wp * ( curr_xp * exp_xsq + real( curr_pot_order + 1, wp ) * f_p )
          curr_xp = curr_xp * x_fac
        end do
      end if

      f_p = 2.0_wp * x**(pot_order-3) * calc_inv_gamma_1_2(pot_order) * f_p
    end select
    
  end function f_p

  function f_p_d(x, energy, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use constants, only : inv_gamma_1_2
    implicit none
    Real ( Kind = wp ) :: f_p_d
    Real ( kind = wp ) :: energy
    Real ( Kind = wp ) :: x
    Integer            :: pot_order

    f_p_d = ( Real(pot_order - 3,wp)/x * energy) - ((2.0_wp/x)*inv_gamma_1_2(pot_order)*exp(-(x**2)) )

  end function f_p_d

  function g_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! General g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use constants, only : rsqrpi
    use numerics,  only : factorial
    use spme,      only : g_1,g_2,g_6,g_12
    implicit none
    Real ( Kind = wp ) :: g_p
    Real ( Kind = wp ), intent ( in    ) :: x
    Integer,            intent ( in    ) :: pot_order
    Real ( Kind = wp ) :: x_2, x_curr
    Real ( Kind = wp ) :: num, den
    Integer :: i

    select case (pot_order)
    case (1)
      g_p = g_1(x)
    case (2)
      g_p = g_2(x)
    case (6)
      g_p = g_6(x)
    case (12)
      g_p = g_12(x)
    case (:0)
      call error(0,'Invalid pot order in g_p')
    case default
      
      x_2 = x**2

      if (mod(pot_order,2) == 0) then !Even orders

        g_p = 0.0_wp
        do i = 0, (pot_order/2)-1
          g_p = g_p + exp(-factorial(i))*x_2**i
        end do

        g_p = g_p * exp(-x_2)

      else ! Odd orders

        x_curr = x
        g_p = 0.0_wp
        den = 1.0_wp
        num = 1.0_wp

        do i = 1, pot_order-1,2
          num = 2.0_wp * num
          den = den / real(i,wp)
          g_p = g_p + x_curr*num*den
          x_curr = x_curr*x_2
        end do
        g_p = g_p * exp(-x_2) * rsqrpi
        g_p = g_p + calc_erfc(x)

      end if

      g_p = g_p
    end select
  
  end function g_p

  function g_p_d(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use constants, only : inv_gamma_1_2
    implicit none
    Real ( Kind = wp ) :: g_p_d

    Real ( Kind = wp ), intent(in) :: x
    Integer,            intent(in) :: pot_order
    g_p_d = 2.0_wp*inv_gamma_1_2(pot_order) * x**(pot_order-1) * exp(-x**2)

  end function g_p_d

End Module ewald_spole
