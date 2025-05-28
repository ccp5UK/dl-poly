Module statistics
!>
!> dl_poly_4 module declaring global simulation property variables and
!> arrays
!
!> copyright - daresbury laboratory
!> author    - i.t.todorov february 2016
!> refactoring:
!           - a.m.elena march-october 2018
!           - j.madge march-october 2018
!           - a.b.g.chalk march-october 2018
!           - i.scivetti march-october 2018
!> contrib  - i.scivetti September 2019. Required changes to allow EVB simulations
!

  Use comms,           Only: &
                             Spread_tag, comm_self, comms_type, gcheck, girecv, gmax, gmin, gsend, gsum, &
                             gsync, gtime, gwait, mode_create, mode_wronly, offset_kind, &
                             gatherv_scatterv_index_arrays, ggatherv, gscatterv, root_id, &
                             gscatter, gbcast
                             
  Use configuration,   Only: configuration_type
  Use constants,       Only: boltz,&
                             engunit,&
                             eu_ev,&
                             eu_kcpm,&
                             eu_kjpm,&
                             pi,&
                             prsunt,&
                             tenunt,&
                             zero_plus,&
                             voigt_6x6,&
                             voigt_flat_3x3
  Use currents,        Only: current_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc,&
                             info,&
                             warning
  Use filename,        Only: FILE_STATS, FILE_COR, FILE_HEATFLUX, &
                             file_type
  Use flow_control,    Only: RESTART_KEY_OLD, flow_type
  Use hash,            Only: hash_table, &
                             MAX_KEY
  Use io,              Only: &
                             io_allocation_error, io_base_comm_not_set, io_close, io_delete, &
                             io_finalize, io_get_parameters, io_history, io_init, io_open, &
                             io_set_parameters, io_type, io_unknown_write_level, &
                             io_unknown_write_option, io_write_batch, io_write_sorted_file
  Use integrators,     Only: trapezium_rule, simpsons_rule, integrator
  
  Use kinds,           Only: STR_LEN,&
                             li,&
                             wi,&
                             wp
  Use numerics,        Only: dcell,&
                             invert,&
                             pbcshfrc,&
                             pbcshfrl,&
                             pbcshift,&
                             shellsort,&
                             shellsort2,&
                             in_range,&
                             scaling_matrix
  Use rigid_bodies,    Only: rigid_bodies_type
  Use site,            Only: site_type
  Use thermostat,      Only: CONSTRAINT_NONE,&
                             CONSTRAINT_SEMI_ORTHORHOMBIC,&
                             CONSTRAINT_SURFACE_TENSION,&
                             DPD_NULL,&
                             thermostat_type
  Use timer,           Only: start_timer,&
                             stop_timer,&
                             timer_type
  Use z_density,       Only: z_density_collect,&
                             z_density_type
  Use correlators,     Only: correlator, correlator_buffer_type, indices_buffer_type
  Use units,           Only: to_out_units

  Implicit None

  Private

  Integer, Parameter          :: MAX_CORRELATION_NAME_LENGTH = 24
  Character(Len=8), Parameter :: stpval_names(1:27) = (/'eng_tot ', 'temp_tot', 'eng_cfg ', 'eng_src ', &
                                                        'eng_cou ', 'eng_bnd ', 'eng_ang ', 'eng_dih ', & 
                                                        'eng_tet ', 'eng_pv  ', 'temp_rot', 'vir_cfg ', &
                                                        'vir_src ', 'vir_cou ', 'vir_bnd ', 'vir_ang ', &
                                                        'vir_con ', 'vir_tet ', 'volume  ', 'temp_shl', &
                                                        'eng_shl ', 'vir_shl ', 'alpha   ', 'beta    ', &
                                                        'gamma   ', 'vir_pmf ', 'press   '/) 

  Integer, Parameter, Public :: NOT_DISTRIBUTED_OBSERVABLE = 0,   &
                                PER_ATOM_OBSERVABLE = 1, &
                                PER_RIGID_OBSERVABLE = 2
  !> correlation observables, and interface
  Type, Abstract, Public :: observable
    Integer          :: component = 0
    Character(Len=2) :: component_name = ""
    Contains
      Procedure(get_value), Deferred         :: value
      Procedure(get_name),  Deferred         :: name
      Procedure(get_id),    Deferred, NoPass :: id
      Procedure,                      NoPass :: distributed => not_distributed
      Procedure                              :: subtype => no_subtype
  End Type observable

  !> stores correlations <AB\>
  Type, Public :: correlation_data
    !> Common observable A.
    Class(observable), Allocatable :: A
    !> Common observable B.
    Class(observable), Allocatable :: B
    !> All (in per atom case) correlations <AB\>.
    Type(correlator), Allocatable :: correlators(:)
    !> If distributed the local indices (else (/0/)).
    Integer, Allocatable           :: indices(:)
    !> If distributed the global indices (else (/0/)).
    Integer, Allocatable           :: indices_global(:)
    !> Step frequency to update the correlation.
    Integer                        :: freq = 0
    !> Number of valid indices in indices and indices_global.
    Integer                        :: indices_used = 0

    Contains
      Procedure, Public :: is_per_atom
      Procedure, Public :: is_per_rigid
  End Type 

  Type, Public :: statistic_accumulator
    Real(Kind=wp)                            :: mu = 0.0_wp, mu_old = 0.0_wp, &
                                                ss = 0.0_wp, var_tmp = 0.0_wp, &
                                                var = 0.0_wp
    Integer                                  :: initialised = 0, window = 0, &
                                                stack_pos = 0
    Real(Kind=wp), Allocatable, Dimension(:) :: stack
    Contains 
      Procedure update_statistic
  End Type

  Type, Public :: stats_type

    Integer(Kind=wi)                   :: numacc = 0, &
                                         natms0 = 0
    Integer(Kind=wi)                   :: mxnstk
    !> Max stack size for rolling averages
    Integer(Kind=wi)                   :: mxstak = 1
    !> Frequency of STATIS output
    Integer(Kind=wi)                   :: intsta = 100
    !> Whether file open
    Logical                            :: statis_file_open = .false.
    !> Whether file is YAML style
    Logical                            :: file_yaml = .false.
    !> Whether stats has been set up
    Logical                            :: newjob = .true.
    !> Whether any bond, angle, etc. analysis
    Logical                            :: lpana = .false.
    Real(Kind=wp)                      :: consv = 0.0_wp, shlke = 0.0_wp, engke = 0.0_wp, &
                                         engrot = 0.0_wp, engcpe = 0.0_wp, engsrp = 0.0_wp, &
                                         engter = 0.0_wp, engtbp = 0.0_wp, engfbp = 0.0_wp, &
                                         engshl = 0.0_wp, engtet = 0.0_wp, engbnd = 0.0_wp, &
                                         engang = 0.0_wp, engdih = 0.0_wp, enginv = 0.0_wp, &
                                         engfld = 0.0_wp, engcon = 0.0_wp, engpmf = 0.0_wp
    Real(Kind=wp)                      :: stptmp = 0.0_wp, stpprs = 0.0_wp, stpvol = 0.0_wp, &
                                         stpcfg = 0.0_wp, stpeng = 0.0_wp, stpeth = 0.0_wp, &
                                         stpvir = 0.0_wp
    Real(Kind=wp)                      :: virtot = 0.0_wp, vircom = 0.0_wp, vircpe = 0.0_wp, &
                                         virsrp = 0.0_wp, virshl = 0.0_wp, virter = 0.0_wp, &
                                         virtbp = 0.0_wp, virfbp = 0.0_wp, vircon = 0.0_wp, &
                                         virpmf = 0.0_wp, virtet = 0.0_wp, virbnd = 0.0_wp, &
                                         virang = 0.0_wp, virdih = 0.0_wp, virinv = 0.0_wp, &
                                         virfld = 0.0_wp, virdpd = 0.0_wp
    Real(Kind=wp)                      :: strtot(1:9) = 0.0_wp, strkin(1:9) = 0.0_wp, strknf(1:9) = 0.0_wp, &
                                         strknt(1:9) = 0.0_wp, strcom(1:9) = 0.0_wp, strcon(1:9) = 0.0_wp, &
                                         strpmf(1:9) = 0.0_wp, stress(1:9) = 0.0_wp, strdpdr(1:9) = 0.0_wp, &
                                         strdpdd(1:9) = 0.0_wp
    Real(Kind=wp)                      :: clin(1:9) = 0.0_wp
    ! constraints accumulators
    Real(Kind=wp), Public              :: passcnq(1:5) = (/ & ! QUENCHING per call
                                         0.0_wp, & ! cycles counter
                                         0.0_wp, & ! access counter
                                         0.0_wp, & ! average cycles
                                         999999999.0_wp, & ! minimum cycles : ~Huge(1)
                                         0.0_wp/) ! maximum cycles
    Real(Kind=wp), Public              :: passcon(1:5, 1:2, 1:2) = Reshape((/ & ! dim::1-shake, dim:1:-per-call
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp, & ! dim::1-shake, dim:2:-per-tst
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp, & ! dim::2-rattle, dim:1:-per-call
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp, & ! dim::2-rattle, dim:2:-per-tst
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp/), (/5, 2, 2/))
    Real(Kind=wp), Public              :: passpmq(1:5) = (/ & ! QUENCHING per call
                                         0.0_wp, & ! cycles counter
                                         0.0_wp, & ! access counter
                                         0.0_wp, & ! average cycles
                                         999999999.0_wp, & ! minimum cycles : ~Huge(1)
                                         0.0_wp/) ! maximum cycles
    Real(Kind=wp), Public              :: passpmf(1:5, 1:2, 1:2) = Reshape((/ & ! dim::1-shake, dim:1:-per-call
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp, & ! dim::1-shake, dim:2:-per-tst
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp, & ! dim::2-rattle, dim:1:-per-call
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp, & ! dim::2-rattle, dim:2:-per-tst
                                                                   0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp/), (/5, 2, 2/))
    Real(Kind=wp), Public              :: passshl(1:5) = (/ &
                                         0.0_wp, & ! cycles counter
                                         0.0_wp, & ! access counter
                                         0.0_wp, & ! average cycles
                                         999999999.0_wp, & ! minimum cycles : ~Huge(1)
                                         0.0_wp/) ! maximum cycles
    !> Skips, elements are as follows
    !>
    !> - 1 skips counter
    !> - 2 access counter
    !> - 3 average skips
    !> - 4 minimum skips ~Huge(1)
    !> - 5 maximum skips
    Real(Kind=wp), Public              :: neighskip(1:5) = [0.0_wp, 0.0_wp, 0.0_wp, &
                                                    999999999.0_wp, 0.0_wp]
    Real(Kind=wp), Public              :: passmin(1:5) = [ &
                                         0.0_wp, & ! cycles counter
                                         0.0_wp, & ! access counter
                                         0.0_wp, & ! average cycles
                                         999999999.0_wp, & ! minimum cycles : ~Huge(1)
                                         0.0_wp] ! maximum cycles
    Type(current_type)                 :: cur
    Logical                            :: calculate_correlations = .false., per_atom_correlations = .false.
    Integer                            :: cor_deport_buffer = 0, cor_dump_freq = 0, number_of_correlations = 0,&
                                          next_cor = 1, currents_correlations = 0
    Type(statistic_accumulator), Allocatable :: accumulators(:)
    Type(correlation_data), Allocatable :: correlations(:)     
    Real(Kind=wp), Allocatable         :: xin(:), yin(:), zin(:)
    Real(Kind=wp), Allocatable         :: xto(:), yto(:), zto(:), rsd(:)
    Real(Kind=wp), Allocatable         :: stpval(:), stpvl0(:), sumval(:), ssqval(:)
    Real(Kind=wp), Allocatable         :: zumval(:), ravval(:), stkval(:, :)
    Integer, Allocatable               :: found(:), found0(:)
    Integer, Allocatable               :: lsi0(:), lsa0(:), lsa00(:), ltg0(:)
    Real(Kind=wp), Allocatable         :: xin0(:), yin0(:), zin0(:)
    Real(Kind=wp), Allocatable         :: xto0(:), yto0(:), zto0(:)
    Real(Kind=wp), Allocatable         :: stpval0(:), stpvl00(:), sumval0(:), ssqval0(:)
    Real(Kind=wp), Allocatable         :: zumval0(:), ravval0(:), stkval0(:, :)

    Integer                            :: pp_eng_str_frequency = 0
    Integer                            :: born_frequency = 0
    Integer                            :: mom_dens_frequency = 0

    !> store spot heat flux
    Real(Kind=wp)                      :: heat_flux(1:3) = 0.0_wp

    !> spot strain tensor
    Real(Kind=wp)                      :: strain(1:9) = 0.0_wp
    Real(Kind=wp), Allocatable         :: inv_ref_scaling_matrix(:, :)
    !> Accumulate for strain tensor
    Type(statistic_accumulator)        :: strain_accum(1:9)

    !> store spot momentum density
    Real(Kind=wp),    Allocatable :: momentum_density(:, :)
    !> atom types to collect momentum density for
    Integer,          Allocatable :: mom_dens_types(:)
    Character(Len=8), Allocatable :: mom_dens_names(:)

    !> Store for per-particle energy data
    Real(Kind=wp), Allocatable         :: pp_energy(:)

    !> Store for per-particle stress data
    Real(Kind=wp), Allocatable         :: pp_stress(:, :)

    !> Store for per-particle current virial contribution, natms, KPOINTS, 1:3
    Complex(Kind=wp), Allocatable      :: pp_cur_virial(:,:,:)

    !> Store for per-particle k-dependent stress tensor, natms, KPOINTS, 1:6
    Complex(Kind=wp), Allocatable      :: pp_cur_stress(:,:,:)
    
    !> Whether this step is a step to collect per-particle data
    Logical :: collect_pp_eng_str = .false.
    Logical :: collect_born = .false.

    !> Stores correlation_data by correlation type
    !>  e.g. stress_xy-v_y. Per atom cors are stored
    !>  together.
    Type(hash_table), Private :: cor_table
    Logical :: elastic_constants = .false.

    !> Store for per-particle Born term (d^2/dr^2 U - 1/r d/dr U) r_i r_j r_k r_l / r^2 using Voigt notation
    Real(Kind=wp)                      :: born_term(1:21) = 0.0_wp
    !> Accumulate the born term
    Type(statistic_accumulator)        :: born_term_accum(1:21)
    !> Flags for computed born components in vdw
    Logical                            :: born_calculate(1:21) = .false.

    !> Stats for stress tensor average
    Type(statistic_accumulator)        :: stress_accum(1:9)
    !> Whether to report properties (temperature, pressure, stresses, viscosity) in DPD units
    Logical :: dpd_units = .false. 

    Logical :: rigid_body_correlations = .false.

  Contains
    Private

    Procedure, Public :: init              => allocate_statistics_arrays
    Procedure, Public :: init_connect      => allocate_statistics_connect
    Procedure, Public :: init_correlations => init_correlations_table
    Procedure, Public :: init_correlator   => allocate_correlator
    Procedure, Public :: init_born_calculate
    Procedure, Public :: clean_connect     => deallocate_statistics_connect
    Procedure, Public :: calculate_strain
    Procedure, Public :: update_stress
    Procedure, Public :: setup_pp_collection
    Procedure, Public :: pp_result
    Procedure, Public :: correlator_deport
    Procedure, Public :: correlator_recieve
    Procedure, Public :: dump_correlations 
    Procedure, Public :: revive_correlations
    Procedure, Public :: check_collection_frequencies
    Procedure, Public :: calculate_stress_energy_current
    Procedure, Public :: setup_momentum_density
    Procedure, Public :: correlator_reindex

    Procedure         :: allocate_per_particle_arrays
    Procedure         :: deallocate_per_particle_arrays
    Final :: cleanup
  End Type

  Abstract Interface 
    !> Kernal for selecting data for correlation
    Function get_value(t, config, rigid, stats, index) Result(v)
        Import observable, configuration_type, rigid_bodies_type, stats_type, wp
        Class(observable),                     Intent(In   ) :: t
        Type(configuration_type),              Intent(InOut) :: config
        Type(rigid_bodies_type),               Intent(InOut) :: rigid
        Type(stats_type),                      Intent(InOut) :: stats
        Integer,       Optional,               Intent(In   ) :: index
        !> Observable value, complex to handle Real and Complex
        Complex(Kind=wp) :: v
    End Function get_value

    !> utility to get name of observable (i.e. for i/o)
    Function get_name(t, with_component) Result(v)
        Import observable, MAX_CORRELATION_NAME_LENGTH
        Class(observable), Intent(In   )           :: t
        Logical,           Intent(In   ), Optional :: with_component
        Character(Len=MAX_CORRELATION_NAME_LENGTH) :: v
    End Function get_name

    !> utility to get numerical id of observable (i.e. for revive)
    Function get_id() Result(v)
      Integer :: v
    End Function get_id

  End Interface

  Type, Extends(observable), Public :: observable_velocity
  Contains
      Procedure         :: value        => velocity_value
      Procedure         :: name         => velocity_name
      Procedure, NoPass :: id           => velocity_id
      Procedure, NoPass :: distributed  => distributed_per_atom
  End Type

  Type, Extends(observable), Public :: observable_stress
  Contains
      Procedure         :: value        => stress_value
      Procedure         :: name         => stress_name
      Procedure, NoPass :: id           => stress_id
  End Type

  Type, Extends(observable), Public :: observable_heat_flux
  Contains
      Procedure         :: value        => heat_flux_value
      Procedure         :: name         => heat_flux_name
      Procedure, NoPass :: id           => heat_flux_id
  End Type

  Type, Extends(observable), Public :: observable_statis
  Contains
      Procedure         :: value        => statis_value
      Procedure         :: name         => statis_name
      Procedure, NoPass :: id           => statis_id
  End Type

  !> Rigid body observable types
  Integer, Parameter :: RIGID_VELOCITY = 1, &
                        RIGID_POSITION = 2, &
                        RIGID_ORIENTATIONAL_VELOCITY = 3

  Type, Extends(observable), Public :: observable_rigid
    Integer :: rigid_type = 0
  Contains
      Procedure         :: value        => rigid_value
      Procedure         :: name         => rigid_name
      Procedure, NoPass :: id           => rigid_id
      Procedure, NoPass :: distributed  => distributed_per_rigid
      Procedure         :: subtype      => rigid_subtype
  End Type

  Integer, Parameter, Public :: L_MOM_CURRENT = 1, T_MOM_CURRENT = 2, ENG_CURRENT = 3, &
                        K_DENSITY = 4, ENG_DENSITY = 5, K_STRESS = 6

  Type, Extends(observable), Public :: observable_currents
    Integer          :: current_type = 0, kpoint = 0, atom_type = 0
    Character(Len=8) :: atom_type_name = ""
  Contains
      Procedure         :: value        => current_value
      Procedure         :: name         => current_name
      Procedure, NoPass :: id           => current_id
  End Type

  
  Public :: calculate_stress
  Public :: calculate_viscosity
  Public :: calculate_heat_flux
  Public :: calculate_mom_density
  Public :: calculate_thermal_conductivity
  Public :: statistics_collect
  Public :: statistics_connect_frames
  Public :: statistics_connect_set
  Public :: write_per_part_contribs
  Public :: write_header
  Public :: statistics_result
  Public :: correlation_result
  Public :: character_to_observable
  Public :: id_component_to_observable
  Public :: set_currents_observable
  Public :: update_statistic

  Interface write_yaml_vector
    Module Procedure write_real_yaml_vector
    Module Procedure write_char_yaml_vector
  End Interface write_yaml_vector
Contains

  Subroutine write_real_yaml_vector(file_unit, name, values, indent)
    Integer,          Intent(In   ) :: file_unit
    Character(Len=*), Intent(In   ) :: name
    Real(Kind=wp),    Intent(In   ) :: values(:)
    Integer,          Intent(In   ) :: indent

    Write (file_unit, '(a,*(g16.8,","))', advance="no") &
      Repeat(" ",indent)//name//": [", values(1:Size(values)-1)
    Write (file_unit, '(g16.8,"]")') values(Size(values))
  End Subroutine write_real_yaml_vector

  Subroutine write_char_yaml_vector(file_unit, name, values, indent)
    Integer,          Intent(In   ) :: file_unit
    Character(Len=*), Intent(In   ) :: name
    Character(Len=*), Intent(In   ) :: values(:)
    Integer,          Intent(In   ) :: indent

    Write (file_unit, '(a,*(a," , "))', advance="no") &
    Repeat(" ",indent)//name//": [", values(1:Size(values)-1)
    Write (file_unit, '(a,"]")') values(Size(values))
  End Subroutine write_char_yaml_vector

  Subroutine check_collection_frequencies(stats, comm)
    Class(stats_type), Intent(InOut) :: stats
    Class(comms_type), Intent(InOut) :: comm

    Call gbcast(comm, stats%pp_eng_str_frequency, root_id)
    Call gbcast(comm, stats%born_frequency, root_id)
    
  End Subroutine check_collection_frequencies

  Subroutine allocate_statistics_arrays(stats, mxrgd, mxatms, mxatdm, mxatype, variable_cell)
    Class(stats_type), Intent(InOut)   :: stats
    Integer,           Intent(In   )   :: mxrgd, mxatms, mxatdm, mxatype
    Logical,           Intent(In   )   :: variable_cell

    Integer                            :: mxnstk, mxstak, nxatms, i
    Integer,           Dimension(1:6)  :: fail
 
    fail = 0

    If (mxrgd > 0) Then
      nxatms = mxatms
    Else
      nxatms = mxatdm
    End If

    mxnstk = stats%mxnstk
    mxstak = stats%mxstak

    Allocate (stats%xin(1:nxatms), stats%yin(1:nxatms), stats%zin(1:nxatms), Stat=fail(1))
    Allocate (stats%xto(1:mxatdm), stats%yto(1:mxatdm), stats%zto(1:mxatdm), stats%rsd(1:mxatdm), Stat=fail(2))
    Allocate (stats%stpval(0:mxnstk), stats%stpvl0(0:mxnstk), stats%sumval(0:mxnstk), stats%ssqval(0:mxnstk), Stat=fail(3))
    Allocate (stats%zumval(0:mxnstk), stats%ravval(0:mxnstk), stats%stkval(1:mxstak, 0:mxnstk), Stat=fail(4))
    
    Allocate (stats%accumulators(0:mxnstk), Stat=fail(5))
    Do i = 0, mxnstk
      Allocate(stats%accumulators(i)%stack(1:mxstak))
      stats%accumulators(i)%window = mxstak
    End Do

    If (.not. Allocated(stats%mom_dens_names)) Then
      Allocate(stats%mom_dens_names(0))
      Allocate(stats%mom_dens_types(0))
    End If
    Allocate(stats%momentum_density(1:Size(stats%mom_dens_names), 1:3), Stat=fail(6))

    If (Any(fail > 0)) Call error_alloc("allocate_statistics_arrays", "statistics")

    stats%xin = 0.0_wp; stats%yin = 0.0_wp; stats%zin = 0.0_wp
    stats%xto = 0.0_wp; stats%yto = 0.0_wp; stats%zto = 0.0_wp; stats%rsd = 0.0_wp

    stats%stpval = 0.0_wp; stats%stpvl0 = 0.0_wp; stats%sumval = 0.0_wp; stats%ssqval = 0.0_wp
    stats%zumval = 0.0_wp; stats%ravval = 0.0_wp; stats%stkval = 0.0_wp

    If (stats%elastic_constants) Then
      Do i = 1, 21
        Allocate(stats%born_term_accum(i)%stack(1:mxstak))
        stats%born_term_accum(i)%window = mxstak
      End Do
      Do i = 1, 9
        Allocate(stats%stress_accum(i)%stack(1:mxstak))
        stats%stress_accum(i)%window = mxstak
      End Do
    End If

    If (variable_cell) Then
      Do i = 1, 9
        Allocate(stats%strain_accum(i)%stack(1:mxstak))
        stats%strain_accum(i)%window = mxstak
      End Do
    End If

  End Subroutine allocate_statistics_arrays

  Subroutine setup_momentum_density(stats, sites)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine to setup momentum density calculations for 
    !  user selected atom types.
    !
    ! author    - h.l.devereux, Nov 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type), Intent(InOut) :: stats
    Type(site_type),   Intent(In   ) :: sites

    Integer                :: i, loc

    Do i = 1, Size(stats%mom_dens_names)
      loc = Findloc(sites%site_name, Trim(stats%mom_dens_names(i)), 1)
      If (loc == 0) Then
        Call error(0, "Could not find atom type "//Trim(stats%mom_dens_names(i))//" for momentum_density")
      End If
      stats%mom_dens_types(i) = loc
    End Do

  End Subroutine

  Subroutine setup_pp_collection(stats, config, flow)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine to check if pp data is computed this step
    !  if so also allocating pp_arrays and setting up switches for
    !  force routines 
    !
    ! author    - h.l.devereux
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type),         Intent(InOut) :: stats
    Type(configuration_type),  Intent(In   ) :: config
    Type(flow_type),           Intent(In   ) :: flow

    If (stats%intsta > 0 .and. flow%step >= flow%equil_steps) Then
      If (stats%pp_eng_str_frequency > 0) Then
        ! If pp data is required, and to calc stats (per-particle data used) AND not equilibration
        If (Mod(flow%step, stats%pp_eng_str_frequency) == 0) Then
          stats%collect_pp_eng_str = .true.
#ifndef HALF_HALO
          Call stats%allocate_per_particle_arrays(config%natms)
#else /* HALF_HALO */
          Call stats%allocate_per_particle_arrays(config%mxatms)
#endif /* HALF_HALO */
        End If
      End If

      If (stats%born_frequency > 0) Then
        If (Mod(flow%step, stats%born_frequency) == 0) Then
          stats%collect_born = .true.
        End If
      End If
    End If
  End Subroutine setup_pp_collection

  Subroutine pp_result(stats, sites, config, comm, flow, files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine to check if pp data was computed this step
    !  if so resetting switches, and if io is passed writing on 
    !  root process to io_file_heatflux
    !
    ! author    - h.l.devereux
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type),         Intent(InOut)           :: stats
    Type(site_type),           Intent(In   )           :: sites
    Type(configuration_type),  Intent(In   )           :: config
    Type(comms_type),          Intent(InOut)           :: comm
    Type(flow_type),           Intent(In   )           :: flow
    Type(file_type),           Intent(InOut), Optional :: files(:)

    Logical :: write = .false.
    Integer :: i

    If (Present(files)) Then 
      write = .true.
    End If
    If (stats%collect_pp_eng_str) Then
      stats%heat_flux = calculate_heat_flux(stats, config, comm)

      If (flow%heat_flux .and. comm%idnode == 0 .and. Mod(flow%step, stats%intsta) == 0) Then
        If (write) Then
          If (flow%step == stats%intsta) Then
            Open (Newunit=files(FILE_HEATFLUX)%unit_no, File=Trim(files(FILE_HEATFLUX)%filename),Status='replace')
          Else
            Open (Newunit=files(FILE_HEATFLUX)%unit_no, File=Trim(files(FILE_HEATFLUX)%filename), Position='append')
          End If
          If (Size(stats%mom_dens_types) > 0) Then
            Write (files(FILE_HEATFLUX)%unit_no, '(I8.1, 1X, 8(G19.12, 1X))', advance='no') flow%step, stats%stptmp, config%volm, &
              stats%heat_flux
            Do i = 1, Size(stats%mom_dens_types)
              Write (files(FILE_HEATFLUX)%unit_no, '(A, 1X, 3(G19.12, 1X))', advance='no') &
                Trim(sites%unique_atom(i)), stats%momentum_density(i, :)
            End Do
          Else
            Write (files(FILE_HEATFLUX)%unit_no, '(I8.1, 1X, 5(G19.12, 1X))') flow%step, stats%stptmp, config%volm, stats%heat_flux
          End If
          Close (files(FILE_HEATFLUX)%unit_no)
        End If
      End If

      If (flow%write_per_particle) Then
        Call write_per_part_contribs(config, comm, stats%pp_energy, stats%pp_stress, flow%step)
      End If

      Call stats%deallocate_per_particle_arrays()
    End If

    stats%collect_pp_eng_str = .false.
    stats%collect_born = .false.

  End Subroutine pp_result

  Subroutine allocate_per_particle_arrays(stats, natms)
    Class(stats_type), Intent(InOut) :: stats
    Integer,           Intent(In   ) :: natms

    Integer :: fail

    If (.not. Allocated(stats%pp_energy)) Then
      Allocate(stats%pp_energy(natms), stat=fail)
      If (fail > 0) Call error_alloc("stats%pp_energy", "statistics")
    End If

    If (.not. Allocated(stats%pp_stress)) Then
      Allocate(stats%pp_stress(9, natms), stat=fail)
      If (fail > 0) Call error_alloc("stats%pp_stress", "statistics")
    End If

    If ((.not. Allocated(stats%pp_cur_virial)) .and. stats%cur%k_energy_stress_current_on) Then
      Allocate(stats%pp_cur_virial(1:natms, 1:stats%cur%nkpoints, 1:3), stat=fail)
      If (fail > 0) Call error_alloc("stats%pp_cur_virial", "statistics")
      stats%pp_cur_virial = Cmplx(0.0_wp, 0.0_wp, Kind=wp)
    End If

    If ((.not. Allocated(stats%pp_cur_stress)) .and. stats%cur%k_energy_stress_current_on) Then
      Allocate(stats%pp_cur_stress(1:natms, 1:stats%cur%nkpoints, 1:6), stat=fail)
      If (fail > 0) Call error_alloc("stats%pp_cur_stress", "statistics")
      stats%pp_cur_stress = Cmplx(0.0_wp, 0.0_wp, Kind=wp)
    End If

    stats%pp_energy = 0.0_wp
    stats%pp_stress = 0.0_wp

  End Subroutine allocate_per_particle_arrays

  Subroutine deallocate_per_particle_arrays(stats)
    Class(stats_type), Intent(InOut) :: stats

    Integer :: fail

    Deallocate(stats%pp_energy, stat=fail)
    If (fail > 0) Call error_dealloc("stats%pp_energy", "statistics")
    Deallocate(stats%pp_stress, stat=fail)
    If (fail > 0) Call error_dealloc("stats%pp_stress", "statistics")

    If (stats%cur%k_energy_stress_current_on) Then
      Deallocate(stats%pp_cur_virial, stat=fail)
      If (fail > 0) Call error_dealloc("stats%pp_cur_virial", "statistics")
      Deallocate(stats%pp_cur_stress, stat=fail)
      If (fail > 0) Call error_dealloc("stats%pp_cur_stress", "statistics")
    End If

  End Subroutine deallocate_per_particle_arrays

  Subroutine allocate_statistics_connect(stats, mxatdm)
    Class(stats_type), Intent(InOut) :: stats
    Integer,           Intent(InOut) :: mxatdm

    Integer                 :: mxstak
    Integer, Dimension(1:6) :: fail

    mxstak = stats%mxstak
    fail = 0

    Allocate (stats%found(1:mxatdm), stats%found0(1:mxatdm), Stat=fail(1))
    Allocate (stats%lsi0(1:mxatdm), stats%lsa0(1:mxatdm), stats%ltg0(1:mxatdm), Stat=fail(2))
    Allocate (stats%xin0(1:mxatdm), stats%yin0(1:mxatdm), stats%zin0(1:mxatdm), Stat=fail(3))
    Allocate (stats%xto0(1:mxatdm), stats%yto0(1:mxatdm), stats%zto0(1:mxatdm), Stat=fail(4))
    Allocate (stats%stpval0(1:2 * mxatdm), stats%stpvl00(1:2 * mxatdm), stats%sumval0(1:2 * mxatdm), stats%ssqval0(1:2 * mxatdm), &
              Stat=fail(5))
    Allocate (stats%zumval0(1:2 * mxatdm), stats%ravval0(1:2 * mxatdm), stats%stkval0(1:mxstak, 1:2 * mxatdm), Stat=fail(6))

    If (Any(fail > 0)) Call error_alloc("allocate_statistics_connect", "statistics")

  End Subroutine allocate_statistics_connect

  Subroutine deallocate_statistics_connect(stats)
    Class(stats_type), Intent(InOut) :: stats

    Integer, Dimension(1:6) :: fail

    fail = 0

    Deallocate (stats%found, stats%found0, Stat=fail(1))
    Deallocate (stats%lsi0, stats%lsa0, stats%ltg0, Stat=fail(2))
    Deallocate (stats%xin0, stats%yin0, stats%zin0, Stat=fail(3))
    Deallocate (stats%xto0, stats%yto0, stats%zto0, Stat=fail(4))
    Deallocate (stats%stpval0, stats%stpvl00, stats%sumval0, stats%ssqval0, Stat=fail(5))
    Deallocate (stats%zumval0, stats%ravval0, stats%stkval0, Stat=fail(6))

    If (Any(fail > 0)) Call error_dealloc("deallocate_statistics_connect", "statistics")

  End Subroutine deallocate_statistics_connect

  Subroutine cleanup(stats)
    Type(stats_type), Intent(InOut) :: stats
    If (Allocated(stats%xin)) Then
      Deallocate (stats%xin)
    End If
    If (Allocated(stats%yin)) Then
      Deallocate (stats%yin)
    End If
    If (Allocated(stats%zin)) Then
      Deallocate (stats%zin)
    End If

    If (Allocated(stats%xto)) Then
      Deallocate (stats%xto)
    End If
    If (Allocated(stats%yto)) Then
      Deallocate (stats%yto)
    End If
    If (Allocated(stats%zto)) Then
      Deallocate (stats%zto)
    End If
    If (Allocated(stats%rsd)) Then
      Deallocate (stats%rsd)
    End If

    If (Allocated(stats%stpval)) Then
      Deallocate (stats%stpval)
    End If
    If (Allocated(stats%stpvl0)) Then
      Deallocate (stats%stpvl0)
    End If
    If (Allocated(stats%sumval)) Then
      Deallocate (stats%sumval)
    End If
    If (Allocated(stats%ssqval)) Then
      Deallocate (stats%ssqval)
    End If

    If (Allocated(stats%zumval)) Then
      Deallocate (stats%zumval)
    End If
    If (Allocated(stats%ravval)) Then
      Deallocate (stats%ravval)
    End If
    If (Allocated(stats%stkval)) Then
      Deallocate (stats%stkval)
    End If

    If (Allocated(stats%found)) Then
      Deallocate (stats%found)
    End If
    If (Allocated(stats%found0)) Then
      Deallocate (stats%found0)
    End If

    If (Allocated(stats%lsi0)) Then
      Deallocate (stats%lsi0)
    End If
    If (Allocated(stats%lsa0)) Then
      Deallocate (stats%lsa0)
    End If
    If (Allocated(stats%lsa00)) Then
      Deallocate (stats%lsa00)
    End If
    If (Allocated(stats%ltg0)) Then
      Deallocate (stats%ltg0)
    End If

    If (Allocated(stats%xin0)) Then
      Deallocate (stats%xin0)
    End If
    If (Allocated(stats%yin0)) Then
      Deallocate (stats%yin0)
    End If
    If (Allocated(stats%zin0)) Then
      Deallocate (stats%zin0)
    End If

    If (Allocated(stats%xto0)) Then
      Deallocate (stats%xto0)
    End If
    If (Allocated(stats%yto0)) Then
      Deallocate (stats%yto0)
    End If
    If (Allocated(stats%zto0)) Then
      Deallocate (stats%zto0)
    End If

    If (Allocated(stats%stpval0)) Then
      Deallocate (stats%stpval0)
    End If
    If (Allocated(stats%stpvl00)) Then
      Deallocate (stats%stpvl00)
    End If
    If (Allocated(stats%sumval0)) Then
      Deallocate (stats%sumval0)
    End If
    If (Allocated(stats%ssqval0)) Then
      Deallocate (stats%ssqval0)
    End If

    If (Allocated(stats%zumval0)) Then
      Deallocate (stats%zumval0)
    End If
    If (Allocated(stats%ravval0)) Then
      Deallocate (stats%ravval0)
    End If
    If (Allocated(stats%stkval0)) Then
      Deallocate (stats%stkval0)
    End If

    If (Allocated(stats%accumulators)) Then 
      Deallocate(stats%accumulators)
    End If

    If (Allocated(stats%inv_ref_scaling_matrix)) Then
      Deallocate(stats%inv_ref_scaling_matrix)
    End If
    
  End Subroutine cleanup

  Subroutine init_correlations_table(stats, currents_cors)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine to initialise the cor_table.
    !
    ! author    - h.l.devereux 2023
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type), Intent(InOut) :: stats
    Integer,           Intent(In   ) :: currents_cors
    Integer :: s

    s = 2*stats%number_of_correlations+2*stats%currents_correlations*currents_cors
    Call stats%cor_table%init(s+1)
    Allocate(stats%correlations(s))
  End Subroutine init_correlations_table

  Logical Function is_per_atom(cor_data)
    Class(correlation_data), Intent(In   ) :: cor_data
    is_per_atom = cor_data%A%distributed() == PER_ATOM_OBSERVABLE .or. &
      cor_data%B%distributed() == PER_ATOM_OBSERVABLE
  End Function is_per_atom

  Logical Function is_per_rigid(cor_data)
    Class(correlation_data), Intent(In   ) :: cor_data
    is_per_rigid = cor_data%A%distributed() == PER_RIGID_OBSERVABLE .or. &
    cor_data%B%distributed() == PER_RIGID_OBSERVABLE
  End Function is_per_rigid

  Logical Function is_not_distributed(cor_data)
    Class(correlation_data), Intent(In   ) :: cor_data
    is_not_distributed = cor_data%A%distributed() == NOT_DISTRIBUTED_OBSERVABLE .and. &
      cor_data%B%distributed() == NOT_DISTRIBUTED_OBSERVABLE
  End Function is_not_distributed

  Subroutine allocate_correlator(stats, per_atom, per_rigid, config, rigid, comm, &
    blocks, points, window, freq, A, B)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for allocating a particular correlation.
    !
    ! author    - h.l.devereux 2023
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type),         Intent(InOut) :: stats
    Class(configuration_type), Intent(InOut) :: config
    Class(rigid_bodies_type),  Intent(InOut) :: rigid
    Class(comms_type),         Intent(InOut) :: comm
    Logical,                   Intent(In   ) :: per_atom
    Logical,                   Intent(In   ) :: per_rigid
    Integer,                   Intent(In   ) :: blocks, points, window, freq
    Class(observable),         Intent(In   ) :: A, B

    Class(correlation_data), Allocatable :: cor_data
    Integer                              :: i
    Integer,                 Allocatable :: indices(:)
    Character(Len=MAX_KEY)               :: correlation_name

    Allocate(indices(0))

    If (per_atom) Then 
      Do i = 1, config%natms
        indices = [indices, i]
      End Do
    Else If (per_rigid) Then
      Do i = 1,  Size(rigid%index_local,2)
        ! Correlator follows the lead atom.
        If (rigid%index_local(1, i) /= 0 .and. rigid%index_local(1, i) < config%natms) Then
          indices = [indices, i]
        End If
      End Do
      If (size(indices) > 0) Then
        stats%rigid_body_correlations = .true.
      End If
    Else If (comm%idnode /= root_id) Then
      Return
    End If

    If (Size(indices) == 0) Then
      indices = [indices, 0]
    End If

    Write (correlation_name, '(a)') Trim(A%name())//"-"//Trim(B%name())

    Allocate(cor_data)
    cor_data%A = A
    cor_data%B = B
    Allocate(cor_data%indices(1:Size(indices)))
    Allocate(cor_data%indices_global(1:Size(indices)))
    Allocate(cor_data%correlators(1:Size(indices)))

    cor_data%indices_used = 0
    Do i = 1, Size(indices)
      Call cor_data%correlators(i)%init(blocks, points, window)
      If (per_atom) Then
        cor_data%indices(i) = indices(i)
        cor_data%indices_global(i) = config%ltg(indices(i))
        cor_data%indices_used = cor_data%indices_used + 1
      Else If (per_rigid) Then
        cor_data%indices(i) = indices(i)
        cor_data%indices_global(i) = rigid%list(1, indices(i))
        cor_data%indices_used = cor_data%indices_used + 1
      Else
        cor_data%indices(i) = 0
        cor_data%indices_global(i) = 0
      End If
    End Do
    cor_data%freq = freq

    If (stats%cor_table%contains(correlation_name)) Then
      Call error(0, "duplicate correlation")
    End If
    Call stats%cor_table%set(correlation_name, stats%next_cor)
    stats%correlations(stats%next_cor) = cor_data
    stats%next_cor = stats%next_cor + 1
  End Subroutine allocate_correlator

  Subroutine write_yaml_correlation(file_unit, name, &
    points_cor, points, blocks, window, freq, steps, value)
    Integer,                     Intent(InOut) :: file_unit
    Integer,                     Intent(In   ) :: points_cor, &
                                                  points, blocks, window, &
                                                  freq
    Character(Len=*),            Intent(In   ) :: name
    Real(Kind=wp), Dimension(:), Intent(In   ) :: steps, value

    If (points_cor > 1) Then
      Write (file_unit, '(a)') "    "//Trim(name)//":"
      Write (file_unit, '(a)')    "        parameters:"
      Write (file_unit, '(a,i0)') "            points_per_block: ", points
      Write (file_unit, '(a,i0)') "            number_of_blocks: ", blocks
      Write (file_unit, '(a,i0)') "            window_size: ", window
      Write(file_unit,  '(a,i0)') "            update_frequency: ", freq
      Write(file_unit, '(a,*(g16.8,","))',advance="no") "        lags: [", steps(1:points_cor-1)
      Write(file_unit, '(g16.8,a)') steps(points_cor), "]"
      Write(file_unit, '(a,*(g16.8,","))',advance="no") "        value: [", value(1:points_cor-1)
      Write(file_unit, '(g16.8,a)') value(points_cor), "]"
    End If
  End Subroutine write_yaml_correlation

  Subroutine gather_correlation(cor_data, points, blocks, window, count, types, names, &
    dt, correlation_name, file_unit, comm)
    Type(correlation_data),                         Intent(InOut) :: cor_data
    Integer,                                        Intent(In   ) :: points, blocks, window, &
                                                                     count, types(:)
    Integer,                                        Intent(InOut) :: file_unit
    Character(Len=*),                               Intent(In   ) :: names(:)
    Real(Kind=wp),                                  Intent(In   ) :: dt
    Character(Len=MAX_CORRELATION_NAME_LENGTH*2+1), Intent(In   ) :: correlation_name
    Type(comms_type),                               Intent(InOut) :: comm

    Real(Kind=wp), Allocatable :: cor_accumulator(:,:), &
                                  correlation(:),       &
                                  flat_correlation(:),  &
                                  timesteps(:)
    Integer,       Allocatable :: type_counts(:)
    Integer                    :: flat_dim, points_cor, atom, j, &
                                  max_points_cor, min_points_cor

    flat_dim = count*points*blocks

    Allocate(cor_accumulator(1:count, 1:points*blocks))
    Allocate(correlation(1:points*blocks))
    Allocate(flat_correlation(1:flat_dim))
    Allocate(type_counts(1:count))
    Allocate(timesteps(1:points*blocks))

    type_counts = 0
    flat_correlation = 0.0_wp
    correlation = 0.0_wp
    cor_accumulator = 0.0_wp
    timesteps = 0.0_wp
    points_cor = 0
    min_points_cor = 0
    max_points_cor = 0
    ! accumulate average
    Do j = 1, cor_data%indices_used
      If (cor_data%correlators(j)%count_updated == 0) Then
        ! no data was seen in this correlator, distinct from 0
        ! correlation case
        Cycle
      End If
      atom = cor_data%indices(j)

      If (types(atom) == 0) Then
        Cycle
      End If
      correlation = 0.0_wp
      Call cor_data%correlators(j)%get_correlation(correlation, timesteps, points_cor)
      cor_accumulator(types(atom),:) = &
        cor_accumulator(types(atom),:) + correlation
      type_counts(types(atom)) = type_counts(types(atom)) + 1
    End Do
    ! now collect onto root
    points_cor = points_cor - 1
    max_points_cor = points_cor
    Call gmax(comm, max_points_cor)
    min_points_cor = points_cor
    Call gmin(comm, min_points_cor)
    If (max_points_cor /= min_points_cor) Then
      Call error(0, "differing number of correlated points between processors for same correlation");
    End If

    flat_correlation = Reshape(cor_accumulator,(/flat_dim/))
    Call gsum(comm, flat_correlation)
    flat_correlation = flat_correlation
    cor_accumulator = Reshape(flat_correlation,(/count,points*blocks/))
    ! for later averaging
    Call gsum(comm,type_counts)
    timesteps = timesteps * cor_data%freq
    If (comm%idnode == root_id .and. points_cor > 1) Then
      Do j = 1,count
        cor_accumulator(j,:) = cor_accumulator(j,:) / (1+type_counts(j))
        Call write_yaml_correlation(file_unit, Trim(names(j))//"-"//Trim(correlation_name), &
          points_cor, points, blocks, window, cor_data%freq, timesteps, &
          cor_accumulator(j,:))
      End Do
    End If
  End Subroutine gather_correlation

  Subroutine correlation_result(stats, rigid, comm, files, config, sites, dt)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for printing a correlation result
    ! (values, lags and derived data) into io_file_cor
    !
    ! author    - h.l.devereux 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Class(stats_type),        Intent(InOut)       :: stats
    Type(rigid_bodies_type),  Intent(InOut)       :: rigid
    Type(comms_type),         Intent(InOut)       :: comm
    Type(file_type),          Intent(InOut)       :: files(:)
    Type(configuration_type), Intent(In   )       :: config
    Type(site_type),          Intent(In   )       :: sites
    Real(Kind=wp),            Intent(In   )       :: dt

    Integer                                        :: i, j, flat_dim, &
                                                      file_unit, atom, points, window, blocks, &
                                                      points_cor, max_points_cor, min_points_cor, &
                                                      cor_index, freq
    Real(Kind=wp), Allocatable                     :: cor_accumulator(:,:), correlation(:), &
                                                      flat_correlation(:), timesteps(:), visc(:), &
                                                      therm_cond(:), k_visc(:)
    Integer,       Allocatable                     :: type_counts(:)
    Character(Len=MAX_CORRELATION_NAME_LENGTH*2+1) :: correlation_name
    Character(Len=2), Dimension(1:3)               :: components_vector
    Character(Len=2), Dimension(1:9)               :: components_matrix
    Real(Kind=wp)                                  :: conv
    Character(Len=STR_LEN)                         :: units_visc, units_therm
    Character(Len=MAX_KEY), Allocatable            :: cor_keys(:)
    Type(correlation_data)                         :: cor_data
    Logical                                        :: per_atom, per_rigid
    Type(observable_currents)                      :: oc
    
                        
    If (stats%calculate_correlations .eqv. .false.) Then 
      Return 
    End If

    components_vector = (/ 'x', 'y', 'z' /)
    components_matrix = (/'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'/)

    file_unit = files(FILE_COR)%unit_no
    If (comm%idnode == root_id) Then

      Open(Newunit=file_unit,File=Trim(files(FILE_COR)%filename),Status='replace')

      Write (file_unit,'(a)') "%YAML 1.2"
      Write (file_unit,'(a)') "---"
      Write (file_unit, '(a,a,a)') "title: '", Trim(config%cfgname), "'"

      Call calculate_viscosity(stats, dt, visc)
      Call calculate_thermal_conductivity(stats, dt, units_therm, therm_cond)

      If (Allocated(therm_cond) .or. Allocated(visc)) Then
        Write (file_unit, '(a)')       "observables:"
      End if

      If (Allocated(visc)) Then
        Call to_out_units(1.0_wp, "internal_m", conv, units_visc)
        Allocate(k_visc(1:Size(visc)))
        k_visc = visc / ((conv*config%totmas) / stats%accumulators(19)%mu)
        Write (file_unit, '(a)')             "      viscosity:"
        Write (file_unit, '(a,g16.8)')       "            value: ", Sum(visc) / Real(Size(visc), Kind=wp)
        If (Size(visc) > 1) Then
          Call write_real_yaml_vector(file_unit, "components", visc, 12)
        End If
        Write (file_unit, '(a)')             "            units: Katm ps "
        Write (file_unit, '(a)')             "      kinematic-viscosity:"
        Write (file_unit, '(a,g16.8)')       "            value: ", Sum(k_visc) / Real(Size(k_visc), Kind=wp)
        If (Size(visc) > 1) Then
          Call write_real_yaml_vector(file_unit, "components", k_visc, 12)
        End If
        Write (file_unit, '(a)')             "            units: Katm ps / ("//Trim(units_visc)//" / Ang^3)"
        Deallocate(visc)
      End If

      If (stats%elastic_constants) Then
        Call elastic_constants_result(stats, config%natms, file_unit)
      End If

      If (Allocated(therm_cond)) Then
        Write (file_unit, '(a)')       "      thermal-conductivity:"
        Write (file_unit, '(a,g16.8)') "            value: ", Sum(therm_cond) / Real(Size(therm_cond), Kind=wp)
        If (Size(therm_cond) > 1) Then
          Call write_real_yaml_vector(file_unit, "components", therm_cond, 12)
        End If
        Write (file_unit, '(a,a)')     "            units: ", Trim(units_therm)
        Deallocate(therm_cond)
      End If

      Write (file_unit, '(a)') "correlations:"
    End If

    Call stats%cor_table%get_keys(cor_keys)

    Do i = 1, Size(cor_keys)
      If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
      Call stats%cor_table%get(cor_keys(i), cor_index)
      Associate(cor_data => stats%correlations(cor_index))
        freq = stats%correlations(cor_index)%freq
        correlation_name = cor_keys(i)
        per_atom = stats%correlations(cor_index)%A%distributed() == PER_ATOM_OBSERVABLE .or. &
                   stats%correlations(cor_index)%B%distributed() == PER_ATOM_OBSERVABLE
        per_rigid = stats%correlations(cor_index)%A%distributed() == PER_RIGID_OBSERVABLE .or. &
                       stats%correlations(cor_index)%B%distributed() == PER_RIGID_OBSERVABLE

        If (Allocated(correlation)) Then
          Deallocate(correlation)
        End If
        If (Allocated(cor_accumulator)) Then
          Deallocate(cor_accumulator)
        End If
        If (Allocated(flat_correlation)) Then
          Deallocate(flat_correlation)
        End If
        If (Allocated(type_counts)) Then
          Deallocate(type_counts)
        End If
        If (Allocated(timesteps)) Then
          Deallocate(timesteps)
        End If

        points = stats%correlations(cor_index)%correlators(1)%points_per_block
        blocks = stats%correlations(cor_index)%correlators(1)%number_of_blocks
        window = stats%correlations(cor_index)%correlators(1)%window_size
        If (per_atom) Then
          Call gather_correlation(stats%correlations(cor_index), points, blocks, window, sites%mxatyp, &
            config%ltype, sites%unique_atom, dt, correlation_name, file_unit, comm)
        Else If (per_rigid) Then
          Call gather_correlation(stats%correlations(cor_index), points, blocks, window, rigid%unique_types, &
            rigid%type, rigid%type_name, dt, correlation_name, file_unit, comm)
        Else If (comm%idnode == root_id) Then
          flat_dim = points*blocks

          Allocate(correlation(1:points*blocks))
          Allocate(timesteps(1:blocks*points))

          correlation = 0.0_wp
          Call stats%correlations(cor_index)%correlators(1)%get_correlation(correlation, timesteps, points_cor)
          timesteps = timesteps * freq
          points_cor = points_cor - 1
          Call write_yaml_correlation(file_unit, correlation_name, &
            points_cor, points, blocks, window, freq, timesteps, correlation)
        End If
      End Associate
    End Do

    If (comm%idnode == root_id) Then
      Close(file_unit)
    End If

  End Subroutine correlation_result
    
  Subroutine statistics_collect(config, rigid, lsim, leql, nsteql, lmsd, keyres, degfre, degshl, &
                                degrot, nstep, tstep, time, tmst, mxatdm, stats, thermo, zdensity, &
                                sites, files, comm, ff, tmr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating periodic data during the
    ! molecular dynamics simulation and computing the rolling averages
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith & i.t.todorov march 2016
    ! contrib   - a.m.elena february 2017
    ! contrib   - i.t.todorov february 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    ! contrib   - i.t.todorov july 2019 - RSD as the true displacement
    !           - m.a.seaton october 2024 - DPD units
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Logical,                  Intent(In   ) :: lsim, leql
    Integer,                  Intent(In   ) :: nsteql
    Logical,                  Intent(In   ) :: lmsd
    Integer,                  Intent(In   ) :: keyres
    Integer(Kind=li),         Intent(In   ) :: degfre, degshl, degrot
    Integer,                  Intent(In   ) :: nstep
    Real(Kind=wp),            Intent(In   ) :: tstep, time
    Real(Kind=wp),            Intent(InOut) :: tmst
    Integer(Kind=wi),         Intent(In   ) :: mxatdm
    Type(stats_type),         Intent(InOut) :: stats
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(site_type),          Intent(In   ) :: sites
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm
    Integer,                  Intent(In   ) :: ff
    Type(timer_type),         Intent(InOut) :: tmr

    Character(Len=100)         :: fmtt, sunits
    Character(Len=STR_LEN)     :: message
    Integer                    :: fail, i, iadd, j, k, kstak, cor_index, strend
    Logical                    :: ffpass, l_tmp
    Real(Kind=wp)              :: celprp(1:10), h_z, sclnv1, sclnv2, stpcns, stpipv, stprot, &
                                  stpshl, zistk, prsunt0, tenunt0, boltz0, btmp(1:9), rtmp
    Complex(Kind=wp)              observable_a, observable_b
    Real(Kind=wp), Allocatable :: amsd(:), xxt(:), yyt(:), zzt(:)

    Character(Len=MAX_KEY), Allocatable, Dimension(:) :: cor_keys
    Type(correlation_data), Pointer :: cor_data

#ifdef CHRONO 
    Call start_timer(tmr, "Statistics Collect")
#endif

    ffpass = ff == 1

    fail = 0

    ! sort out conversion factors for pressure, stress, tension etc.

    prsunt0 = Merge(1.0_wp, prsunt, stats%dpd_units)
    tenunt0 = Merge(1.0_wp, tenunt, stats%dpd_units)
    boltz0 = Merge(1.0_wp, boltz, stats%dpd_units)

    Allocate (amsd(1:sites%mxatyp), Stat=fail)
    If (fail > 0) Call error_alloc('amsd', 'statistics_collect')

    If (.not. Allocated(stats%inv_ref_scaling_matrix)) Then
      ! User has not set the reference cell, take it as the inital cell
      Allocate(stats%inv_ref_scaling_matrix(1:3, 1:3))
      Call scaling_matrix(config%cell, stats%inv_ref_scaling_matrix)
      Call invert(Reshape(stats%inv_ref_scaling_matrix, [9]), btmp, rtmp)
      stats%inv_ref_scaling_matrix = Transpose(Reshape(btmp, [3, 3]))
    End If

    ! open statistics file and put header
    If (stats%intsta > 0 .and. stats%newjob .and. comm%idnode == 0 .and. ffpass) Then
      stats%newjob = .false.

      ! If the keyres = RESTART_KEY_OLD is the file old (does it exist)?

      l_tmp = .false.
      If (keyres == RESTART_KEY_OLD) Inquire (File=files(FILE_STATS)%filename, Exist=l_tmp)

      If (.not. l_tmp) Then
        Open (Newunit=files(FILE_STATS)%unit_no, File=files(FILE_STATS)%filename, Status='replace')
        stats%statis_file_open = .true.

        If (stats%file_yaml) Then
          Write (files(FILE_STATS)%unit_no, '(a)') "%YAML 1.2"
          Write (files(FILE_STATS)%unit_no, '(a)') "---"
        End If
        Write (files(FILE_STATS)%unit_no, '(a,a,a)') "title: '", config%cfgname, "'"

        If (Abs(engunit - eu_ev) <= zero_plus) Then
          sunits = "electron Volts"
        Else If (Abs(engunit - eu_kcpm) <= zero_plus) Then
          sunits = "kcal/mol"
        Else If (Abs(engunit - eu_kjpm) <= zero_plus) Then
          sunits = "kjoule/mol"
        Else If (Abs(engunit - 1.0_wp) <= zero_plus) Then
          If (stats%dpd_units) Then
            sunits = "DPD units"
          Else
            sunits = "DL_POLY Internal UNITS (10 J/mol)"
          End If
        Else If (Abs(engunit - boltz) <= zero_plus) Then
          sunits = "Kelvin/Boltzmann"
        Else ! once in a blue moon
          sunits = "Unknown"
        End If
        Write (files(FILE_STATS)%unit_no, '(a,a)') "energy unitS: ", Trim(sunits)
        If (stats%file_yaml) Then
          Write (files(FILE_STATS)%unit_no, '(a,a)') "labels: "
          Write (files(FILE_STATS)%unit_no, '(2x,a4,*(a,", "))', advance="no") "- [ ", &
            'step', 'time', 'Total Extended System Energy', 'System Temperature', &
            'Configurational Energy', 'Short Range Potential Energy', 'Electrostatic Energy', &
            'Chemical Bond Energy', 'Valence Angle And 3-Body Potential Energy', &
            'Dihedral Inversion And 4-Body Potential Energy', &
            'Tethering Energy', 'Enthalpy (Total Energy + Pv)', 'Rotational Temperature', 'Total Virial', &
            'Short-Range Virial', 'Electrostatic Virial', 'Bond Virial', 'Valence Angle And 3-Body Virial', &
            'Constraint Bond Virial', 'Tethering Virial', 'Volume', 'Core-Shell Temperature', &
            'Core-Shell Potential Energy', 'Core-Shell Virial', 'Md Cell Angle Α', &
            'Md Cell Angle Β', 'Md Cell Angle Gamma', 'Pmf Constraint Virial', 'Pressure', &
            'External Degree Of Freedom', 'stress xx', 'stress xy', 'stress xz', 'stress yx', &
            'stress yy', 'stress yz', 'stress zx', 'stress zy', 'stress zz'
          If (thermo%key_dpd /= DPD_NULL) &
            Write (files(FILE_STATS)%unit_no, '(*(a,", "))', advance="no") 'config stress xx', 'config stress xy', &
              'config stress xz', 'config stress yx', 'config stress yy', 'config stress yz', 'config stress zx', &
              'config stress zy', 'config stress zz', 'dissipative stress xx', 'dissipative stress xy', 'dissipative stress xz', &
              'dissipative stress yx', 'dissipative stress yy', 'dissipative stress yz', 'dissipative stress zx', &
              'dissipative stress zy', 'dissipative stress zz', 'random stress xx', 'random stress xy', 'random stress xz', &
              'random stress yx', 'random stress yy', 'random stress yz', 'random stress zx', 'random stress zy', &
              'random stress zz', 'kinetic stress xx', 'kinetic stress xy', 'kinetic stress xz', 'kinetic stress yx', &
              'kinetic stress yy', 'kinetic stress yz', 'kinetic stress zx', 'kinetic stress zy', 'kinetic stress zz'
          Do i = 1, sites%ntype_atom - 1
            Write (files(FILE_STATS)%unit_no, '(a)', advance="no") "amsd "//sites%unique_atom(i)//", "
          End Do
          Write (files(FILE_STATS)%unit_no, '(a)', advance="no") "amsd "//sites%unique_atom(sites%ntype_atom)
          If (thermo%variable_cell) Then
            Write (files(FILE_STATS)%unit_no, '(", ",*(a,", "))', advance="no") "cell A1", "cell A2", "cell A3", &
              "cell B1", "cell B2", "cell B3", "cell C1", "cell C2", "cell C3"
            Write (files(FILE_STATS)%unit_no, '(a)', advance="no") "pV"

            If (thermo%iso /= CONSTRAINT_NONE) Then
              Write (files(FILE_STATS)%unit_no, '(a)', advance="no") ",h_z, A_z"
              If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
                Write (files(FILE_STATS)%unit_no, '(a)', advance="no") ",gamma_x, gamma_y"
              End If
            End If
          End If
          Write (files(FILE_STATS)%unit_no, '(a2)') " ]"
          Write (files(FILE_STATS)%unit_no, '(a,a)') "timesteps: "
        End If
      End If
    End If

    ! instantaneous properties of system

    ! system energy
    ! Configurational energy has been defined in subroutine w_calculateorces within drivers.F90
    ! In the case of EVB calculations, the configurational energy is recomputed via diagonalisation
    ! of the EVB matrix (subroutine evb.F90)

    ! Configurational stats%stpcfg energy has been defined in subroutine calculate_forces within drivers.F90

    stats%stpeng = stats%stpcfg + stats%engke + stats%engrot

    ! energy + conserved quantity (for true ensembles)

    stpcns = stats%stpeng + stats%consv

    ! rotational temperature

    stprot = 2.0_wp * (stats%engrot) / (boltz0 * Max(1.0_wp, Real(degrot, wp)))

    ! core-shell units temperature

    stpshl = 2.0_wp * (stats%shlke) / (boltz0 * Max(1.0_wp, Real(degshl, wp)))

    ! system temperature

    stats%stptmp = 2.0_wp * (stats%engke + stats%engrot) / (boltz0 * Real(degfre, wp))

    ! system virial, stats%virtot has been computed in calculate_forces

    stats%stpvir = stats%virtot + stats%vircon + stats%virpmf + stats%vircom + stats%virdpd

    ! system volume

    stats%stpvol = config%volm

    ! system pressure

    stats%stpprs = (2.0_wp * stats%engke - stats%stpvir) / (3.0_wp * stats%stpvol)

    ! system PV

    stpipv = stats%stpprs * stats%stpvol

    ! system enthalpy

    If (thermo%variable_cell) Then ! P_target*V_instantaneous
      stats%stpeth = stats%stpeng + (thermo%press + sum(thermo%stress(1:9:4))/3.0_wp) * stats%stpvol
    Else ! for thermo%variable_cell=.false. V_instantaneous=V_target
      stats%stpeth = stats%stpeng + stpipv ! and there is only P_instantaneous
    End If

    Call dcell(config%cell, celprp)

    stats%strain = stats%calculate_strain(config)

    ! store current values in statistics array

    stats%stpval(0) = stats%consv / engunit
    stats%stpval(1) = stpcns / engunit
    stats%stpval(2) = stats%stptmp
    stats%stpval(3) = stats%stpcfg / engunit
    stats%stpval(4) = (stats%engsrp + stats%engter) / engunit
    stats%stpval(5) = stats%engcpe / engunit
    stats%stpval(6) = stats%engbnd / engunit
    stats%stpval(7) = (stats%engang + stats%engtbp) / engunit
    stats%stpval(8) = (stats%engdih + stats%enginv + stats%engfbp) / engunit
    stats%stpval(9) = stats%engtet / engunit
    stats%stpval(10) = stats%stpeth / engunit
    stats%stpval(11) = stprot
    stats%stpval(12) = stats%stpvir / engunit
    stats%stpval(13) = (stats%virsrp + stats%virter) / engunit
    stats%stpval(14) = stats%vircpe / engunit
    stats%stpval(15) = stats%virbnd / engunit
    stats%stpval(16) = (stats%virtbp + stats%virang) / engunit
    stats%stpval(17) = stats%vircon / engunit
    stats%stpval(18) = stats%virtet / engunit
    stats%stpval(19) = stats%stpvol
    stats%stpval(20) = stpshl
    stats%stpval(21) = stats%engshl / engunit
    stats%stpval(22) = stats%virshl / engunit
    stats%stpval(23) = Acos(celprp(6)) * 180.0_wp / pi
    stats%stpval(24) = Acos(celprp(5)) * 180.0_wp / pi
    stats%stpval(25) = Acos(celprp(4)) * 180.0_wp / pi
    stats%stpval(26) = stats%virpmf / engunit
    stats%stpval(27) = stats%stpprs * prsunt0

    iadd = 27

    ! iadd = iadd + 1 ! for the stpval(0)!!! Thus to account for in printing
    ! pressure tensor (derived for the stress tensor)

    Do i = 1, 9
      stats%stpval(iadd + i) = stats%strtot(i) * prsunt0 / stats%stpvol
    End Do
    iadd = iadd + 9

    If (stats%cur%on .and. Mod(nstep, stats%intsta) == 0 .and. nstep >= nsteql) Then
      If (stats%cur%k_energy_stress_current_on) Then
        Call stats%cur%compute(config, time, comm, sites, stats%pp_energy, &
          stats%pp_cur_virial, stats%pp_cur_stress)
      Else
        Call stats%cur%compute(config, time, comm, sites, stats%pp_energy)
      End If
    End If

    If (stats%mom_dens_frequency > 0 .and. nstep > 0) Then
      If (Mod(nstep, stats%mom_dens_frequency) == 0) Then
        Do i = 1, Size(stats%mom_dens_types)
          stats%momentum_density(i, :) = calculate_mom_density(stats, i, config, comm)
        End Do
      End If
    End If

    ! separated stress tensors for DPD calculations:
    ! conservative/configurational, dissipative, random, kinetic

    If (thermo%key_dpd/=DPD_NULL) Then
      Do i = 1, 9
        stats%stpval(iadd + i) = (stats%strcon(i) + stats%strpmf(i) + stats%stress(i) + stats%strcom(i)) * prsunt0/stats%stpvol
      End Do
      iadd = iadd + 9
      Do i = 1, 9
        stats%stpval(iadd + i) = stats%strdpdd(i) * prsunt0 / stats%stpvol
      End Do
      iadd = iadd + 9
      Do i = 1, 9
        stats%stpval(iadd + i) = stats%strdpdr(i) * prsunt0 / stats%stpvol
      End Do
      iadd = iadd + 9
      Do i = 1, 9
        stats%stpval(iadd + i) = stats%strkin(i) * prsunt0 / stats%stpvol
      End Do
      iadd = iadd + 9
    End If

    ! mean squared displacements per species, dependent on
    ! particle displacements from initial positions (at t=0)

    amsd = 0.0_wp ! initialise

    If (nstep == nsteql + 1) Then ! re-initialise
      Do i = 1, config%natms
        stats%xto(i) = 0.0_wp
        stats%yto(i) = 0.0_wp
        stats%zto(i) = 0.0_wp
      End Do
    End If
    If (nstep > 0) Then
      If (lsim) Then ! real dynamics is happening
        Do i = 1, config%natms

          stats%xto(i) = stats%xto(i) + config%vxx(i) * tstep
          stats%yto(i) = stats%yto(i) + config%vyy(i) * tstep
          stats%zto(i) = stats%zto(i) + config%vzz(i) * tstep

        End Do

        If (stats%calculate_correlations) Then
          If ((.not. leql) .or. nstep >= nsteql) Then
            Call stats%cor_table%get_keys(cor_keys)
            Do j = 1, Size(cor_keys)
              If (.not. stats%cor_table%contains(cor_keys(j))) Cycle
              Call stats%cor_table%get(cor_keys(j), cor_index)
              Associate(cor_data => stats%correlations(cor_index))
                If (nstep > 0 .and. &
                    Mod(nstep, cor_data%freq) == 0) Then
                  If (cor_data%indices_used == 0 &
                      .and. comm%idnode == root_id) Then
                    ! not distributed correlation
                    observable_a = cor_data%A%value(config, rigid, stats)
                    observable_b = cor_data%B%value(config, rigid, stats)
                    Call cor_data%correlators(1)%update(observable_a, observable_b)
                  Else
                    Do i = 1, cor_data%indices_used
                      observable_a = &
                        cor_data%A%value(config, rigid, stats, cor_data%indices(i))
                      observable_b = &
                        cor_data%B%value(config, rigid, stats, cor_data%indices(i))
                      Call cor_data%correlators(i)%update(observable_a, observable_b)
                    End Do
                  End If
                End If
              End Associate
            End Do
          End If
        End If

      Else ! HISTORY is replayed
        Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail)
        If (fail > 0) Call error_alloc("atomic positions", "statistics_collect")
        Do i = 1, config%natms
          xxt(i) = config%parts(i)%xxx
          yyt(i) = config%parts(i)%yyy
          zzt(i) = config%parts(i)%zzz
        End Do
        Call pbcshfrc(config%imcon, config%cell, config%natms, xxt, yyt, zzt)
        Call pbcshfrc(config%imcon, stats%clin, config%natms, stats%xin, stats%yin, stats%zin)
        Do i = 1, config%natms
          stats%xin(i) = xxt(i) - stats%xin(i)
          stats%yin(i) = yyt(i) - stats%yin(i)
          stats%zin(i) = zzt(i) - stats%zin(i)
        End Do
        Deallocate (xxt, yyt, zzt, Stat=fail)
        If (fail > 0) Call error_dealloc("atomic positions", "statistics_collect")

        Call pbcshfrl(config%imcon, config%cell, config%natms, stats%xin, stats%yin, stats%zin)
        Do i = 1, config%natms
          stats%xto(i) = stats%xto(i) + stats%xin(i)
          stats%yto(i) = stats%yto(i) + stats%yin(i)
          stats%zto(i) = stats%zto(i) + stats%zin(i)
        End Do
      End If

      Do i = 1, config%natms
        stats%rsd(i) = Sqrt(stats%xto(i)**2 + stats%yto(i)**2 + stats%zto(i)**2)

        k = config%ltype(i)
        amsd(k) = amsd(k) + stats%rsd(i)**2
      End Do
      Call gsum(comm, amsd(1:sites%ntype_atom))
    End If

    If (lmsd) Then
      Do i = 1, config%natms
        j = 2 * i
        stats%stpval(iadd + j - 1) = stats%rsd(i)**2
        stats%stpval(iadd + j) = config%vxx(i)**2 + config%vyy(i)**2 + config%vzz(i)**2
      End Do
      iadd = iadd + 2 * mxatdm
    End If

! Calculate true displacements from original position in RSD,
! rather than keep the RMSD=Sqrt(MSD)

    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'statistics_collect allocation failure 2'
      Call error(0, message)
    End If
    Do i = 1, config%natms
      xxt(i) = config%parts(i)%xxx - stats%xin(i)
      yyt(i) = config%parts(i)%yyy - stats%yin(i)
      zzt(i) = config%parts(i)%zzz - stats%zin(i)
    End Do
    Call pbcshift(config%imcon, config%cell, config%natms, xxt, yyt, zzt)
    Do i = 1, config%natms
      stats%rsd(i) = xxt(i)**2 + yyt(i)**2 + zzt(i)**2
    End Do
    Do i = 1, config%natms
      stats%rsd(i) = Sqrt(stats%rsd(i))
    End Do
    Deallocate (xxt, yyt, zzt, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'statistics_collect deallocation failure 2'
      Call error(0, message)
    End If

    Do k = 1, sites%ntype_atom
      If (sites%num_type_nf(k) > zero_plus) Then
        stats%stpval(iadd + k) = amsd(k) / sites%num_type_nf(k)
      Else
        stats%stpval(iadd + k) = 0.0_wp
      End If
    End Do
    iadd = iadd + sites%ntype_atom

    If (thermo%variable_cell) Then

      ! cell parameters

      Do i = 1, 9
        stats%stpval(iadd + i) = config%cell(i)
      End Do
      iadd = iadd + 9

      ! instantaneous PV

      stats%stpval(iadd + 1) = stpipv / engunit
      iadd = iadd + 1

      If (thermo%iso /= CONSTRAINT_NONE) Then
        h_z = celprp(9)

        stats%stpval(iadd + 1) = h_z
        stats%stpval(iadd + 2) = stats%stpvol / h_z
        iadd = iadd + 2

        If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
          stats%stpval(iadd + 1) = -h_z * (stats%strtot(1) - (thermo%press + thermo%stress(1))) * tenunt
          stats%stpval(iadd + 2) = -h_z * (stats%strtot(5) - (thermo%press + thermo%stress(5))) * tenunt
          iadd = iadd + 2
        End If
      End If
    End If

    ! write statistics file
    If (stats%intsta > 0) Then 
      If (comm%idnode == 0 .and. Mod(nstep, stats%intsta) == 0 .and. ffpass) Then
        If (.not. stats%statis_file_open) Then
          Open (Newunit=files(FILE_STATS)%unit_no, File=files(FILE_STATS)%filename, Position='append')
          stats%statis_file_open = .true.
        End If

        strend = Merge (36, 72, thermo%key_dpd==DPD_NULL)
        If (lmsd) Then
          If (stats%file_yaml) Then
            Write (fmtt, '(a,i0,a)') '(2x,a4,i0,",",', iadd + 1 - 2 * mxatdm, '(g16.8,","),g16.8,a2)'
            Write (files(FILE_STATS)%unit_no, fmt=Trim(fmtt)) "- [ ", nstep, time, &
              stats%stpval(1:27), stats%stpval(0), stats%stpval(28:strend), &
              stats%stpval(strend + 1 + 2 * mxatdm:iadd), ' ]'
          Else
            Write (files(FILE_STATS)%unit_no, '(i10,1p,e14.6,0p,i10,/, (1p,5e14.6))') &
              nstep, time, iadd + 1 - 2 * mxatdm, stats%stpval(1:27), stats%stpval(0), stats%stpval(28:strend), &
              stats%stpval(strend + 1 + 2 * mxatdm:iadd)
          End If
        Else
          If (stats%file_yaml) Then
            Write (fmtt, '(a,i0,a)') '(2x,a4,i0,",",', iadd + 1, '(g16.8,","),g16.8,a2)'
            Write (files(FILE_STATS)%unit_no, fmt=Trim(fmtt)) "- [ ", nstep, time, &
              stats%stpval(1:27), stats%stpval(0), stats%stpval(28:iadd), ' ]'
          Else
            Write (files(FILE_STATS)%unit_no, '(i10,1p,e14.6,0p,i10,/, (1p,5e14.6))') &
              nstep, time, iadd + 1, stats%stpval(1:27), stats%stpval(0), stats%stpval(28:iadd)
          End If
        End If

      End If

      If (stats%cor_dump_freq > 0 .and. nstep >=stats%intsta) Then
        If (Mod(nstep, stats%cor_dump_freq) == 0) Then
          Call correlation_result(stats, rigid, comm, files, config, sites, thermo%tstep)
        End If
      End If 
    End If
    ! check on number of variables for stack

    If (iadd > stats%mxnstk) Call error(170)

    ! No totals for timestep zero

    If (nstep /= 0) Then

      Do i = 0, stats%mxnstk
        Call stats%accumulators(i)%update_statistic(stats%stpval(i), nstep)   
      End Do

      If (stats%born_frequency > 0) Then
        If (Mod(nstep, stats%born_frequency) == 0) Then
          Call gsum(comm, stats%born_term)
          Do i = 1, 21
            Call stats%born_term_accum(i)%update_statistic(stats%born_term(i), nstep)
          End Do
          Do i = 1, 9
            Call stats%stress_accum(i)%update_statistic(stats%strtot(i)/stats%stpvol, nstep)
          End Do
        End If
      End If

      If (thermo%variable_cell) Then
        Do i = 1, 9
          Call stats%strain_accum(i)%update_statistic(stats%strain(i), nstep)
        End Do
      End If

      ! current stack value

      kstak = Mod(nstep - 1, stats%mxstak) + 1

      ! subtract old stack value from the stack average

      If (nstep > stats%mxstak) Then
        Do i = 0, stats%mxnstk
          stats%zumval(i) = stats%zumval(i) - stats%stkval(kstak, i)
        End Do
      End If

      ! store quantities in stack and update the stack average

      Do i = 0, stats%mxnstk
        stats%stkval(kstak, i) = stats%stpval(i)
        stats%zumval(i) = stats%zumval(i) + stats%stpval(i)
      End Do

      ! calculate rolling averages

      zistk = Real(Min(stats%mxstak, nstep), wp)

      Do i = 0, stats%mxnstk
        stats%ravval(i) = stats%zumval(i) / zistk
      End Do

      ! accumulate totals over steps

      If ((.not. leql) .or. nstep > nsteql) Then
        stats%numacc = stats%numacc + 1
        sclnv2 = 1.0_wp / Real(stats%numacc, wp)
        sclnv1 = Real(stats%numacc - 1, wp) / Real(stats%numacc, wp)

        ! average squared sum and sum (keep in this order!!!)

        If (nstep == nsteql + 1 .or. ((.not. leql) .and. nstep == 1)) stats%stpvl0 = stats%stpval
        stats%stpval = stats%stpval - stats%stpvl0
        Do i = 0, stats%mxnstk
          stats%ssqval(i) = sclnv1 * (stats%ssqval(i) + sclnv2 * (stats%stpval(i) - stats%sumval(i))**2)

          !stats%sumval has to be shifted back tostats%sumval+stpvl0 in statistics_result
          ! when averaging is printed since stpval is only shifted back and forth
          ! which does not affect the fluctuations Sqrtstats%ssqval) only their accuracy

          stats%sumval(i) = sclnv1 * stats%sumval(i) + sclnv2 * stats%stpval(i)
        End Do
        stats%stpval = stats%stpval + stats%stpvl0
      End If
    End If

    ! z-density collection

    If (zdensity%l_collect) Then
      If (((.not. leql) .or. nstep >= nsteql) .and. &
          Mod(nstep, zdensity%frequency) == 0) Call z_density_collect(zdensity, config)
    End If

    ! Catch time of starting statistical averages

    If (((.not. leql) .or. nstep == nsteql) .and. tmst < tstep) tmst = time

    Deallocate (amsd, Stat=fail)
    If (fail > 0) Call error_alloc("amsd", "statistics_collect")

#ifdef CHRONO 
    Call stop_timer(tmr, "Statistics Collect")
#endif

  End Subroutine statistics_collect

  Subroutine statistics_connect_frames(config, megatm, mxatdm, lmsd, dpd, stats, domain, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange exchange of data between neighbouring
    ! domains/nodes in order to reconnect some statistical information
    ! between replayed frames of history
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer(Kind=wi), Intent(In) :: megatm
    Integer(Kind=wi), Intent(In) :: mxatdm
    Logical, Intent(In) :: lmsd, dpd
    Type(stats_type), Intent(InOut) :: stats
    Type(domains_type), Intent(In) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type), Intent(InOut) :: comm

    Integer :: icyc, nres
    Character(Len=STR_LEN) :: message

    stats%found = 0; icyc = 0; nres = 1
    Do While (icyc <= Max(domain%nx, domain%ny, domain%nz) / 2 .and. nres > 0)
      Call match_compress_spread_sort(-1, mxatdm, stats%lsa00) ! -x direction spread
      Call match_compress_spread_sort(1, mxatdm, stats%lsa00) ! +x direction spread

      Call match_compress_spread_sort(-2, mxatdm, stats%lsa00) ! -y direction spread
      Call match_compress_spread_sort(2, mxatdm, stats%lsa00) ! +y direction spread

      Call match_compress_spread_sort(-3, mxatdm, stats%lsa00) ! -z direction spread
      Call match_compress_spread_sort(3, mxatdm, stats%lsa00) ! +z direction spread

      Call match_compress_spread_sort(0, mxatdm, stats%lsa00) ! no spreading

      nres = stats%natms0
      Call gsum(comm, nres)
      If (nres > 0) Then
        nres = Merge(0, Sum(stats%found(1:config%natms)), config%natms > 0)
        Call gsum(comm, nres)
        If (nres /= megatm) icyc = icyc + 1
      End If
    End Do

    If (nres > 0) Then
      Write (message, '(a)') ' particles dynamics properties will be corrupted'
      Call warning(message, .true.)
    End If

  Contains

    Subroutine match_compress_spread_sort(mdir, mxatdm, lsa00)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! dl_poly_4 routine to simplify the repetition of the procedures above
      !
      ! copyright - daresbury laboratory
      ! author    - i.t.todorov february 2016
      ! refactoring:
      !           - a.m.elena march-october 2018
      !           - j.madge march-october 2018
      !           - a.b.g.chalk march-october 2018
      !           - i.scivetti march-october 2018
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer(Kind=wi),              Intent(In   ) :: mdir, mxatdm
    Integer(Kind=wi), Allocatable, Intent(InOut) :: lsa00(:)

    Integer(Kind=wi) :: fail, i, i0, j, j0, kk

! +/-1,+/-2,+/-3,0 is the direction of spread

      ! Search for matches

      stats%found0 = 0
      If (stats%natms0 > 0) Then
        i0 = 1
        Do i = 1, config%natms
          Do While (i0 <= stats%natms0)
            If (config%lsa(i) < stats%lsa0(i0)) Then
              Exit
            Else If (config%lsa(i) == stats%lsa0(i0)) Then
              If (stats%found(config%lsi(i)) > 0) Then ! ghost arrival
                stats%found0(stats%lsi0(i0)) = 1 ! erase at compression
              Else ! new arrival to claim
                stats%found(config%lsi(i)) = 1; stats%found0(stats%lsi0(i0)) = 1

                stats%xin(config%lsi(i)) = stats%xin0(stats%lsi0(i0))
                stats%yin(config%lsi(i)) = stats%yin0(stats%lsi0(i0))
                stats%zin(config%lsi(i)) = stats%zin0(stats%lsi0(i0))

                stats%xto(config%lsi(i)) = stats%xto0(stats%lsi0(i0))
                stats%yto(config%lsi(i)) = stats%yto0(stats%lsi0(i0))
                stats%zto(config%lsi(i)) = stats%zto0(stats%lsi0(i0))

                If (lmsd) Then
                  j = Merge (72, 36, dpd) + 2 * config%lsi(i)
                  j0 = 2 * stats%lsi0(i0)
                  stats%stpvl0(j - 1) = stats%stpvl00(j0 - 1)
                  stats%stpvl0(j) = stats%stpvl00(j0)
                  stats%stpval(j - 1) = stats%stpval0(j0 - 1)
                  stats%stpval(j) = stats%stpval0(j0)
                  stats%zumval(j - 1) = stats%zumval0(j0 - 1)
                  stats%zumval(j) = stats%zumval0(j0)
                  stats%ravval(j - 1) = stats%ravval0(j0 - 1)
                  stats%ravval(j) = stats%ravval0(j0)
                  stats%ssqval(j - 1) = stats%ssqval0(j0 - 1)
                  stats%ssqval(j) = stats%ssqval0(j0)
                  stats%sumval(j - 1) = stats%sumval0(j0 - 1)
                  stats%sumval(j) = stats%sumval0(j0)
                  Do kk = 1, stats%mxstak
                    stats%stkval(kk, j - 1) = stats%stkval0(kk, j0 - 1)
                    stats%stkval(kk, j) = stats%stkval0(kk, j0)
                  End Do
                End If
              End If
            End If
            i0 = i0 + 1 ! move along
          End Do
        End Do
      End If

      ! Invalidate lazies and deallocate at last use

      If (mdir == 0) Then
        i = 1
        Do i0 = 1, stats%natms0
          Do While (lsa00(i) /= 0 .and. i < mxatdm)
            If (stats%lsa0(i0) < lsa00(i)) Then
              Exit
            Else If (stats%lsa0(i0) == lsa00(i)) Then
              stats%found0(stats%lsi0(i0)) = 1 ! erase at compression
            End If
            i = i + 1 ! move along
          End Do
        End Do

        Deallocate (lsa00, Stat=fail)
        If (fail > 0) Then
          Write (message, '(a)') 'match_compress_spread_sort deallocation failure'
          Call error(0, message)
        End If
      End If

      ! Compress remainder

      i0 = 1
      Do While (i0 <= stats%natms0 .and. stats%natms0 > 0)
        If (stats%found0(i0) == 0) Then ! Not claimed
          If (config%ixyz(i0) == 0) Then ! pass along
            If (mdir == -1) Then ! -x to do
              config%ixyz(i0) = 333 ! all since b4 0
            Else If (mdir == 1) Then !  x to do
              config%ixyz(i0) = 331 ! not back to -1
            Else If (mdir == -2) Then ! -y to do
              config%ixyz(i0) = 332 ! not back to  1
            Else If (mdir == 2) Then !  y to do
              config%ixyz(i0) = 313 ! not back to -2
            Else If (mdir == -3) Then ! -z to do
              config%ixyz(i0) = 323 ! not back to  2
            Else If (mdir == 3) Then !  z to do
              config%ixyz(i0) = 133 ! not back to -3
            Else If (mdir == 0) Then ! end of cycle to do
              config%ixyz(i0) = 233 ! not back to  3
            Else ! abort
              Call error(160)
            End If
          End If

          i0 = i0 + 1 ! Increase lower bound marker
        Else ! claimed, to erase entry
          If (stats%found0(stats%natms0) == 0) Then ! try to refill with the last unclaimed entry
            config%ixyz(i0) = config%ixyz(stats%natms0)
            If (config%ixyz(i0) == 0) Then ! pass along
              If (mdir == -1) Then ! -x to do
                config%ixyz(i0) = 333 ! all since b4 0
              Else If (mdir == 1) Then !  x to do
                config%ixyz(i0) = 331 ! not back to -1
              Else If (mdir == -2) Then ! -y to do
                config%ixyz(i0) = 332 ! not back to  1
              Else If (mdir == 2) Then !  y to do
                config%ixyz(i0) = 313 ! not back to -2
              Else If (mdir == -3) Then ! -z to do
                config%ixyz(i0) = 323 ! not back to  2
              Else If (mdir == 3) Then !  z to do
                config%ixyz(i0) = 133 ! not back to -3
              Else If (mdir == 0) Then ! end of cycle to do
                config%ixyz(i0) = 233 ! not back to  3
              Else ! abort
                Call error(160)
              End If
            End If
            stats%ltg0(i0) = stats%ltg0(stats%natms0)

            stats%xin0(i0) = stats%xin0(stats%natms0)
            stats%yin0(i0) = stats%yin0(stats%natms0)
            stats%zin0(i0) = stats%zin0(stats%natms0)

            stats%xto0(i0) = stats%xto0(stats%natms0)
            stats%yto0(i0) = stats%yto0(stats%natms0)
            stats%zto0(i0) = stats%zto0(stats%natms0)

            If (lmsd) Then
              j = 2 * i0
              j0 = 2 * stats%natms0
              stats%stpvl00(j - 1) = stats%stpvl00(j0 - 1)
              stats%stpvl00(j) = stats%stpvl00(j0)
              stats%stpval0(j - 1) = stats%stpval0(j0 - 1)
              stats%stpval0(j) = stats%stpval0(j0)
              stats%zumval0(j - 1) = stats%zumval0(j0 - 1)
              stats%zumval0(j) = stats%zumval0(j0)
              stats%ravval0(j - 1) = stats%ravval0(j0 - 1)
              stats%ravval0(j) = stats%ravval0(j0)
              stats%ssqval0(j - 1) = stats%ssqval0(j0 - 1)
              stats%ssqval0(j) = stats%ssqval0(j0)
              stats%sumval0(j - 1) = stats%sumval0(j0 - 1)
              stats%sumval0(j) = stats%sumval0(j0)
              Do kk = 1, stats%mxstak
                stats%stkval0(kk, j - 1) = stats%stkval0(kk, j0 - 1)
                stats%stkval0(kk, j) = stats%stkval0(kk, j0)
              End Do
            End If

            i0 = i0 + 1 ! increase lower bound marker if entry is refilled
          End If

          !! Erase upper holdings in either case
          !
          !          ixyz(stats%natms0) = 0
          !         stats%ltg0(stats%natms0) = 0
          !
          !          stats%xin0(stats%natms0) = 0
          !          stats%yin0(stats%natms0) = 0
          !          stats%zin0(stats%natms0) = 0
          !
          !         stats%xto0(stats%natms0) = 0
          !         stats%yto0(stats%natms0) = 0
          !         stats%zto0(stats%natms0) = 0
          !
          !          If (lmsd) Then
          !             j0=2*stats%natms0
          !            stats%stpvl00(j0-1)=0.0_wp
          !            stats%stpvl00(j0  )=0.0_wp
          !            stats%stpval0(j0-1)=0.0_wp
          !            stats%stpval0(j0  )=0.0_wp
          !            stats%zumval0(j0-1)=0.0_wp
          !            stats%zumval0(j0  )=0.0_wp
          !            stats%ravval0(j0-1)=0.0_wp
          !            stats%ravval0(j0  )=0.0_wp
          !            stats%ssqval0(j0-1)=0.0_wp
          !            stats%ssqval0(j0  )=0.0_wp
          !            stats%sumval0(j0-1)=0.0_wp
          !            stats%sumval0(j0  )=0.0_wp
          !             Do kk=1,mxstak
          !               stats%stkval0(kk,j0-1)=0.0_wp
          !               stats%stkval0(kk,j0  )=0.0_wp
          !             End Do
          !          End If
          stats%natms0 = stats%natms0 - 1 ! Decrease upper bound marker
        End If
      End Do

      ! Allocate and initialise at first use
      ! Detect unknown lazies sort them lsa like

      If (mdir == -1) Then
        fail = 0
        Allocate (lsa00(1:mxatdm), Stat=fail)
        If (fail > 0) Then
          Write (message, '(a)') 'match_compress_spread_sort allocation failure'
          Call error(0, message)
        End If
        lsa00 = 0

        i = 0
        Do i0 = 1, stats%natms0
          If (config%ixyz(i0) == 333) Then
            i = i + 1
            lsa00(i) = stats%ltg0(i0)
          End If
        End Do
        Call shellsort(i, lsa00)
      End If

      ! Spread atom data in the mdir direction

      If (mdir /= 0) Call statistics_connect_spread(config, mdir, mxatdm, lmsd, dpd, stats, domain, comm)

      ! Sort past frame remainder of global atom indices

      !    lsi0=0 ; lsa0=0
      Do i0 = 1, stats%natms0
        stats%lsi0(i0) = i0
        stats%lsa0(i0) = stats%ltg0(i0)
      End Do
      Call shellsort2(stats%natms0, stats%lsi0, stats%lsa0)

    End Subroutine match_compress_spread_sort

  End Subroutine statistics_connect_frames

  Subroutine statistics_connect_set(config, rcut, mxatdm, lmsd, dpd, stats, domain, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange exchange of data between neighbouring
    ! domains/nodes in order to reconnect some statistical information
    ! between replayed frames of history
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(In) :: rcut
    Integer(Kind=wi), Intent(In) :: mxatdm
    Logical, Intent(In) :: lmsd, dpd
    Type(stats_type), Intent(InOut) :: stats
    Type(domains_type), Intent(In) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type), Intent(InOut) :: comm

    Real(Kind=wp) :: cut

    Integer           :: nlx, nly, nlz, i, i0, kk, strend
    Real(Kind=wp) :: det, celprp(1:10), rcell(1:9), x, y, z, &
                     xdc, ydc, zdc, cwx, cwy, cwz, ecwx, ecwy, ecwz

    If (comm%mxnode > 1) Then

      ! Define cut

      cut = rcut + 1.0e-6_wp

      Call dcell(config%cell, celprp)
      Call invert(config%cell, rcell, det)

      ! calculate link cell dimensions per node

      nlx = Int(celprp(7) / (cut * domain%nx_real))
      nly = Int(celprp(8) / (cut * domain%ny_real))
      nlz = Int(celprp(9) / (cut * domain%nz_real))

      ! Get the total number of link-cells in MD cell per direction

      xdc = Real(nlx * domain%nx, wp)
      ydc = Real(nly * domain%ny, wp)
      zdc = Real(nlz * domain%nz, wp)

      ! link-cell widths in reduced space

      cwx = 1.0_wp / xdc
      cwy = 1.0_wp / ydc
      cwz = 1.0_wp / zdc

      ! Distance from the - edge of this domain

      ecwx = Nearest((-0.5_wp + cwx) + Real(domain%idx, wp) * domain%nx_recip, +1.0_wp) + zero_plus
      ecwy = Nearest((-0.5_wp + cwy) + Real(domain%idy, wp) * domain%ny_recip, +1.0_wp) + zero_plus
      ecwz = Nearest((-0.5_wp + cwz) + Real(domain%idz, wp) * domain%nz_recip, +1.0_wp) + zero_plus

      ! Distance from the + edge of this domain with a possible
      ! extension strip for the one linked cell per domain scenario

      cwx = Nearest((-0.5_wp - cwx) + Real(domain%idx + 1, wp) * domain%nx_recip, -1.0_wp) - &
            zero_plus - Merge(cwx * 1.0e-10_wp, 0.0_wp, nlx == 1)
      cwy = Nearest((-0.5_wp - cwy) + Real(domain%idy + 1, wp) * domain%ny_recip, -1.0_wp) - &
            zero_plus - Merge(cwy * 1.0e-10_wp, 0.0_wp, nly == 1)
      cwz = Nearest((-0.5_wp - cwz) + Real(domain%idz + 1, wp) * domain%nz_recip, -1.0_wp) - &
            zero_plus - Merge(cwz * 1.0e-10_wp, 0.0_wp, nlz == 1)

      config%ixyz(1:mxatdm) = 0 ! Initialise move (former halo) indicator
      Do i = 1, config%natms
        x = rcell(1) * config%parts(i)%xxx + rcell(4) * config%parts(i)%yyy + rcell(7) * config%parts(i)%zzz
        y = rcell(2) * config%parts(i)%xxx + rcell(5) * config%parts(i)%yyy + rcell(8) * config%parts(i)%zzz
        z = rcell(3) * config%parts(i)%xxx + rcell(6) * config%parts(i)%yyy + rcell(9) * config%parts(i)%zzz

        If (x <= ecwx) config%ixyz(i) = config%ixyz(i) + 1
        If (x >= cwx) config%ixyz(i) = config%ixyz(i) + 2

        If (y <= ecwy) config%ixyz(i) = config%ixyz(i) + 10
        If (y >= cwy) config%ixyz(i) = config%ixyz(i) + 20

        If (z <= ecwz) config%ixyz(i) = config%ixyz(i) + 100
        If (z >= cwz) config%ixyz(i) = config%ixyz(i) + 200
      End Do

      config%lsi = 0; config%lsa = 0 ! This is a must, unfortunately
      Do i = 1, config%natms
        config%lsi(i) = i
        config%lsa(i) = config%ltg(i)
      End Do
      Call shellsort2(config%natms, config%lsi, config%lsa)

      stats%natms0 = config%natms
      stats%ltg0(1:stats%natms0) = config%ltg(1:stats%natms0) !;stats%ltg0(stats%natms0+1: ) = 0
      stats%lsa0(1:stats%natms0) = config%lsa(1:stats%natms0) !; lsa0(stats%natms0+1: ) = 0
      stats%lsi0(1:stats%natms0) = config%lsi(1:stats%natms0) !; lsi0(stats%natms0+1: ) = 0

      stats%xin0(1:stats%natms0) = stats%xin(1:stats%natms0) !; stats%xin0(stats%natms0+1: ) = 0 ; stats%xin = 0.0_wp
      stats%yin0(1:stats%natms0) = stats%yin(1:stats%natms0) !; stats%yin0(stats%natms0+1: ) = 0 ; stats%yin = 0.0_wp
      stats%zin0(1:stats%natms0) = stats%zin(1:stats%natms0) !; stats%zin0(stats%natms0+1: ) = 0 ; stats%zin = 0.0_wp

      stats%xto0(1:stats%natms0) = stats%xto(1:stats%natms0) !;stats%xto0(stats%natms0+1: ) = 0
      stats%yto0(1:stats%natms0) = stats%yto(1:stats%natms0) !;stats%yto0(stats%natms0+1: ) = 0
      stats%zto0(1:stats%natms0) = stats%zto(1:stats%natms0) !;stats%zto0(stats%natms0+1: ) = 0

      strend = Merge (72, 36, dpd)
      If (lmsd) Then
        i0 = 2 * stats%natms0
        stats%stpvl00(1:i0) = stats%stpvl0(strend+1:strend + i0) !;stats%stpvl00(i0+1: )=0.0_wp
        stats%stpval0(1:i0) = stats%stpval(strend+1:strend + i0) !;stats%stpval0(i0+1: )=0.0_wp
        stats%zumval0(1:i0) = stats%zumval(strend+1:strend + i0) !;stats%zumval0(i0+1: )=0.0_wp
        stats%ravval0(1:i0) = stats%ravval(strend+1:strend + i0) !;stats%ravval0(i0+1: )=0.0_wp
        stats%ssqval0(1:i0) = stats%ssqval(strend+1:strend + i0) !;stats%ssqval0(i0+1: )=0.0_wp
        stats%sumval0(1:i0) = stats%sumval(strend+1:strend + i0) !;stats%sumval0(i0+1: )=0.0_wp
        Do kk = 1, stats%mxstak
          stats%stkval0(kk, 1:i0) = stats%stkval(kk, 1:i0) !;stats%stkval0(kk,i0+1: )=0.0_wp
        End Do
      End If
    End If

  End Subroutine statistics_connect_set

  Subroutine statistics_connect_spread(config, mdir, mxatdm, lmsd, dpd, stats, domain, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to spread atomic and topological data of particles
    ! leaving this domain
    !
    ! NOTE: When executing on one node we need not get here at all!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Integer(Kind=wi),         Intent(In   ) :: mdir, mxatdm
    Logical,                  Intent(In   ) :: lmsd, dpd
    Type(stats_type),         Intent(InOut) :: stats
    Type(domains_type),       Intent(In   ) :: domain
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)                   :: message
    Integer                                  :: fail, i, iblock, imove, ix, iy, iz, j, jdnode, jj, &
                                                jmove, jxyz, kdnode, keep, kk, kmove, kx, kxyz, &
                                                ky, kz, l, newatm, send
    Logical                                  :: move, safe, stay
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer

    If (comm%mxnode == 1) Return

    fail = 0
    Allocate (buffer(1:config%mxbfss), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'statistics_connect_spread allocation failure'
      Call error(0, message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock = config%mxbfss / 2

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! kxyz - corrected halo reduction factor particles haloing both +&- sides
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0; ky = 0; kz = 0
    If (mdir == -1) Then ! Direction -x
      kx = 1
      jxyz = 1
      kxyz = Merge(3, jxyz, domain%nx <= 2)

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir == 1) Then ! Direction +x
      kx = 1
      jxyz = 2
      kxyz = Merge(3, jxyz, domain%nx <= 2)

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky = 1
      jxyz = 10
      kxyz = Merge(30, jxyz, domain%ny <= 2)

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir == 2) Then ! Direction +y
      ky = 1
      jxyz = 20
      kxyz = Merge(30, jxyz, domain%ny <= 2)

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz = 1
      jxyz = 100
      kxyz = Merge(300, jxyz, domain%nz <= 2)

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir == 3) Then ! Direction +z
      kz = 1
      jxyz = 200
      kxyz = Merge(300, jxyz, domain%nz <= 2)

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(160)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! buffer(1) and buffer(iblock+1) contain the actual number of
    ! particles to get transferred, imove and jmove are the lengths of
    ! the buffers

    imove = 1
    jmove = 1

    ! Initialise how many particles are to be kept and sent

    keep = 0
    send = 0

    ! Initialise array overflow flags

    safe = .true.

    ! LOOP OVER ALL PREVIOUS FRAME'S PARTICLES ON THIS NODE

    Do i = 1, stats%natms0

      ! particle designated directions

      ix = Mod(config%ixyz(i), 10) ! [0,1,2,3]
      iy = Mod(config%ixyz(i) - ix, 100) ! [0,10,20,30]
      iz = Mod(config%ixyz(i) - (ix + iy), 1000) ! [0,100,200,300]

      ! Filter the move index for the selected direction

      j = ix * kx + iy * ky + iz * kz

      ! If the particle is scheduled to be sent in the selected
      ! direction then indicate it in move

      move = .false.
      If (j == jxyz .or. (j > jxyz .and. Mod(j, 3) == 0)) Then
        move = (comm%idnode /= jdnode) ! but don't move it if back to itself

        ! reduce particle move index (ixyz) using the corrected halo reduction
        ! factor when the particle is sent to both +&- sides

        config%ixyz(i) = config%ixyz(i) - Merge(kxyz, jxyz, j /= jxyz)
      End If
      stay = (config%ixyz(i) /= 0) ! decide on keeping it when to be sent elsewhere

      If (stay) Then ! keep it
        keep = keep + 1

        ! retain config indexing and move indexing arrays

        stats%ltg0(keep) = stats%ltg0(i)
        config%ixyz(keep) = config%ixyz(i)

        ! retain initial positions

        stats%xin0(keep) = stats%xin0(i)
        stats%yin0(keep) = stats%yin0(i)
        stats%zin0(keep) = stats%zin0(i)

        ! retain final displacements

        stats%xto0(keep) = stats%xto0(i)
        stats%yto0(keep) = stats%yto0(i)
        stats%zto0(keep) = stats%zto0(i)

        If (lmsd) Then
          jj = 2 * i
          j = 2 * keep
          stats%stpvl00(j - 1) = stats%stpvl00(jj - 1)
          stats%stpvl00(j) = stats%stpvl00(jj)
          stats%stpval0(j - 1) = stats%stpval0(jj - 1)
          stats%stpval0(j) = stats%stpval0(jj)
          stats%zumval0(j - 1) = stats%zumval0(jj - 1)
          stats%zumval0(j) = stats%zumval0(jj)
          stats%ravval0(j - 1) = stats%ravval0(jj - 1)
          stats%ravval0(j) = stats%ravval0(jj)
          stats%ssqval0(j - 1) = stats%ssqval0(jj - 1)
          stats%ssqval0(j) = stats%ssqval0(jj)
          stats%sumval0(j - 1) = stats%sumval0(jj - 1)
          stats%sumval0(j) = stats%sumval0(jj)
          Do kk = 1, stats%mxstak
            stats%stkval0(kk, j - 1) = stats%stkval0(kk, jj - 1)
            stats%stkval0(kk, j) = stats%stkval0(kk, jj)
          End Do
        End If
      End If

      If (move) Then ! copy it
        send = send + 1
        If (imove + 8 <= iblock) Then ! If safe to proceed

          ! pack config indexing and move indexing arrays

          buffer(imove + 1) = Real(stats%ltg0(i), wp)
          buffer(imove + 2) = Real(config%ixyz(i), wp)

          ! pack initial positions

          buffer(imove + 3) = stats%xin0(i)
          buffer(imove + 4) = stats%yin0(i)
          buffer(imove + 5) = stats%zin0(i)

          ! pack final displacements

          buffer(imove + 6) = stats%xto0(i)
          buffer(imove + 7) = stats%yto0(i)
          buffer(imove + 8) = stats%zto0(i)
        Else
          safe = .false.
        End If
        imove = imove + 8

        ! pack MSD arrays

        If (lmsd) Then
          If (imove + 2 * (6 + stats%mxstak) <= iblock) Then
            jj = 2 * i
            buffer(imove + 1) = stats%stpvl00(jj - 1)
            buffer(imove + 2) = stats%stpvl00(jj)
            buffer(imove + 3) = stats%stpval0(jj - 1)
            buffer(imove + 4) = stats%stpval0(jj)
            buffer(imove + 5) = stats%zumval0(jj - 1)
            buffer(imove + 6) = stats%zumval0(jj)
            buffer(imove + 7) = stats%ravval0(jj - 1)
            buffer(imove + 8) = stats%ravval0(jj)
            buffer(imove + 9) = stats%ssqval0(jj - 1)
            buffer(imove + 10) = stats%ssqval0(jj)
            buffer(imove + 11) = stats%sumval0(jj - 1)
            buffer(imove + 12) = stats%sumval0(jj)
            Do kk = 1, stats%mxstak
              l = 12 + 2 * kk
              buffer(imove + l - 1) = stats%stkval0(kk, jj - 1)
              buffer(imove + l) = stats%stkval0(kk, jj)
            End Do
          Else
            safe = .false.
          End If
          imove = imove + 2 * (6 + stats%mxstak)
        End If
      End If

    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm, safe)
    If (.not. safe) Call error(163)

    ! record of number of atoms for transfer

    buffer(1) = Real(send, wp)

    ! exchange information on buffer sizes

    Call girecv(comm, jmove, kdnode, Spread_tag)
    Call gsend(comm, imove, jdnode, Spread_tag)
    Call gwait(comm)

    ! exchange buffers between nodes (this is a MUST)

    Call girecv(comm, buffer(iblock + 1:iblock + jmove), kdnode, Spread_tag)
    Call gsend(comm, buffer(1:imove), jdnode, Spread_tag)
    Call gwait(comm)

    ! check arrays can cope with incoming atom numbers

    kmove = iblock + 1
    jmove = Nint(buffer(kmove))

    ! Test for overloading and collect how many are to really be accepted

    imove = 0
    Do i = 1, jmove
      l = Nint(buffer(kmove + 1))
      If (All(stats%ltg0(1:stats%natms0) /= l)) imove = imove + 1
      kmove = kmove + 8
      If (lmsd) kmove = kmove + 2 * (6 + stats%mxstak)
    End Do

    stats%natms0 = keep + imove

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe = (stats%natms0 <= mxatdm)
    Call gcheck(comm, safe)
    If (.not. safe) Call error(164)

    ! load transferred data

    kmove = iblock + 1 ! restore kmove
    newatm = keep ! restore newatm
    Do i = 1, jmove
      If (imove /= jmove) Then
        l = Nint(buffer(kmove + 1))
        If (Any(stats%ltg0(1:keep) == l)) Then
          kmove = kmove + 8
          If (lmsd) kmove = kmove + 2 * (6 + stats%mxstak)
          Cycle
        End If
      End If

      newatm = newatm + 1

      ! unpack config indexing, site and move indexing arrays

      stats%ltg0(newatm) = Nint(buffer(kmove + 1))
      config%ixyz(newatm) = Nint(buffer(kmove + 2))

      ! unpack initial positions arrays

      stats%xin0(newatm) = buffer(kmove + 3)
      stats%yin0(newatm) = buffer(kmove + 4)
      stats%zin0(newatm) = buffer(kmove + 5)

      ! unpack initial positions arrays

      stats%xto0(newatm) = buffer(kmove + 6)
      stats%yto0(newatm) = buffer(kmove + 7)
      stats%zto0(newatm) = buffer(kmove + 8)

      kmove = kmove + 8

      ! unpack MSD arrays

      If (lmsd) Then
        jj = 2 * newatm
        stats%stpvl00(jj - 1) = buffer(kmove + 1)
        stats%stpvl00(jj) = buffer(kmove + 2)
        stats%stpval0(jj - 1) = buffer(kmove + 3)
        stats%stpval0(jj) = buffer(kmove + 4)
        stats%zumval0(jj - 1) = buffer(kmove + 5)
        stats%zumval0(jj) = buffer(kmove + 6)
        stats%ravval0(jj - 1) = buffer(kmove + 7)
        stats%ravval0(jj) = buffer(kmove + 8)
        stats%ssqval0(jj - 1) = buffer(kmove + 9)
        stats%ssqval0(jj) = buffer(kmove + 10)
        stats%sumval0(jj - 1) = buffer(kmove + 11)
        stats%sumval0(jj) = buffer(kmove + 12)
        Do kk = 1, stats%mxstak
          l = 12 + 2 * kk
          stats%stkval0(kk, jj - 1) = buffer(kmove + l - 1)
          stats%stkval0(kk, jj) = buffer(kmove + l)
        End Do

        kmove = kmove + 2 * (6 + stats%mxstak)
      End If
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'statistics_connect_spread deallocation failure'
      Call error(0, message)
    End If

  End Subroutine statistics_connect_spread

  Subroutine write_header(dpd_units)
    Logical, Intent(In) :: dpd_units
    Character(Len=STR_LEN), Dimension(5) :: messages

    Write (messages(1), '("#",a)') Repeat('-', 130)
    If (dpd_units) Then
      Write (messages(2), '("#",9x,a4,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
        'step', 'eng_tot', 'temp_tot', 'eng_cfg', 'eng_src', 'eng_cou', 'eng_bnd', 'eng_ang', 'eng_dih', 'eng_tet'
      Write (messages(3), '("#",2x,a11,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
        'time[dpd_t]', ' eng_pv', 'temp_rot', 'vir_cfg', 'vir_src', 'vir_cou', 'vir_bnd', 'vir_ang', 'vir_con', 'vir_tet'
      Write (messages(4), '("#",2x,a11,5x,a7,4x,a8,5x,a7,5x,a7,4x,a8,5x,a7,4x,a8,5x,a7,7x,a5)') &
        'cpu     [s]', 'volume', 'temp_shl', 'eng_shl', 'vir_shl', 'alpha[o]', 'beta[o]', 'gamma[o]', 'vir_pmf', 'press'
    Else
      Write (messages(2), '("#",9x,a4,5x,a7,1x,a11,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
        'step', 'eng_tot', 'temp_tot[K]', 'eng_cfg', 'eng_src', 'eng_cou', 'eng_bnd', 'eng_ang', 'eng_dih', 'eng_tet'
      Write (messages(3), '("#",5x,a8,5x,a7,1x,a11,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
        'time[ps]', ' eng_pv', 'temp_rot[K]', 'vir_cfg', 'vir_src', 'vir_cou', 'vir_bnd', 'vir_ang', 'vir_con', 'vir_tet'
      Write (messages(4), '("#",5x,a8,5x,a7,1x,a11,5x,a7,5x,a7,4x,a8,5x,a7,4x,a8,5x,a7,7x,a5)') &
        'cpu  [s]', 'volume', 'temp_shl[K]', 'eng_shl', 'vir_shl', 'alpha[o]', 'beta[o]', 'gamma[o]', 'vir_pmf', 'press'
    End If
    Write (messages(5), '("#",a)') Repeat('-', 130)
    Call info(messages, 5, .true.)

  End Subroutine write_header

  Subroutine statistics_result(config, minim, lmsd, &
                               nstrun, keyshl, megcon, megpmf, &
                               nstep, time, tmst, &
                               mxatdm, neigh_uncond_update, stats, &
                               rigid, thermo, sites, comm, files, tmr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for writing simulation summary
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith & i.t.todorov november 2016
    ! contrib   - m.a.seaton june 2014
    ! contrib   - a.b.g.chalk january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
    Logical,                  Intent(In   ) :: minim, lmsd
    Integer(Kind=wi),         Intent(In   ) :: nstrun, keyshl, megcon, megpmf, nstep
    Real(Kind=wp),            Intent(In   ) :: time, tmst
    Integer(Kind=wi),         Intent(In   ) :: mxatdm
    Logical,                  Intent(InOut) :: neigh_uncond_update
    Type(stats_type),         Intent(InOut) :: stats
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(site_type),          Intent(In   ) :: sites
    Type(comms_type),         Intent(InOut) :: comm
    Type(file_type),          Intent(InOut) :: files(:)
    Type(timer_type),         Intent(InOut) :: tmr

    Character(Len=STR_LEN)               :: message
    Character(Len=STR_LEN), Dimension(5) :: messages
    Integer                              :: i, iadd, mxnstk, strend
    Logical                              :: check
    Real(Kind=wp)                        :: avvol, dc, h_z, srmsd, timelp, tmp, tx, ty, prsunt0, tenunt0

    mxnstk = stats%mxnstk

    Call info('', .true.)

    ! sort out units for pressure, stress, tension etc.

    prsunt0 = Merge(1.0_wp, prsunt, stats%dpd_units)
    tenunt0 = Merge(1.0_wp, tenunt, stats%dpd_units)

    ! VNL skipping statistics

    If (neigh_uncond_update .and. nstep > 0) Then

      Write (message, '(a,f7.2,2(a,i4))') &
        '# VNL skipping run statistics - skips per timestep: average ', stats%neighskip(3), &
        ' minimum ', Nint(Merge(stats%neighskip(4), stats%neighskip(5), stats%neighskip(4) < stats%neighskip(5))), &
        ' maximum ', Nint(stats%neighskip(5))
      Call info(message, .true.)
    End If

    ! minimisation convergence statistics

    If (minim) Then
      Write (message, '(a,f7.2,2(a,i4))') &
        '# minimisation run statistics - cycles per call: average ', stats%passmin(3), &
        ' minimum ', Nint(stats%passmin(4)), ' maximum ', Nint(stats%passmin(5))
      Call info(message, .true.)
    End If

    ! shell relaxation convergence statistics

    If (keyshl == 2) Then
      Write (message, '(a,f7.2,2(a,i4))') &
        '# shell relaxation run statistics - cycles per timestep: average ', stats%passshl(3), &
        ' minimum ', Nint(stats%passshl(4)), ' maximum ', Nint(stats%passshl(5))
      Call info(message, .true.)
    End If

    ! bond constraints iterative cycles statistics

    If (megcon > 0) Then
      Call gmax(comm, stats%passcon(3:5, 1, 1)); Call gmax(comm, stats%passcon(3:5, 2, 1))
      If (stats%passcon(3, 1, 1) > 0.0_wp) Then
        Write (message, '(2(a,f5.2),4(a,i3))') &
          '# constraints shake  run statistics - cycles per call/timestep: average ', &
          stats%passcon(3, 1, 1), ' / ', stats%passcon(3, 2, 1), &
          ' minimum ', Nint(stats%passcon(4, 1, 1)), ' / ', Nint(stats%passcon(4, 2, 1)), &
          ' maximum ', Nint(stats%passcon(5, 1, 1)), ' / ', Nint(stats%passcon(5, 2, 1))
        Call info(message, .true.)
      End If

      Call gmax(comm, stats%passcon(3:5, 1, 2)); Call gmax(comm, stats%passcon(3:5, 2, 2))
      If (stats%passcon(3, 1, 2) > 0.0_wp) Then
        Write (message, '(2(a,f5.2),4(a,i3))') &
          '# constraints rattle  run statistics - cycles per call/timestep: average ', &
          stats%passcon(3, 1, 1), ' / ', stats%passcon(3, 2, 1), &
          ' minimum ', Nint(stats%passcon(4, 1, 2)), ' / ', Nint(stats%passcon(4, 2, 2)), &
          ' maximum ', Nint(stats%passcon(5, 1, 2)), ' / ', Nint(stats%passcon(5, 2, 2))
        Call info(message, .true.)
      End If
    End If

    ! PMF constraints iterative cycles statistics

    If (megpmf > 0) Then
      Call gmax(comm, stats%passpmf(3:5, 1, 1)); Call gmax(comm, stats%passpmf(3:5, 2, 1))
      If (stats%passpmf(3, 1, 1) > 0.0_wp) Then
        Write (message, '(2(a,f5.2),4(a,i3))') &
          '# PMFs shake  run statistics - cycles per call/timestep: average ', &
          stats%passpmf(3, 1, 1), ' / ', stats%passpmf(3, 2, 1), &
          ' minimum ', Nint(stats%passpmf(4, 1, 1)), ' / ', Nint(stats%passpmf(4, 2, 1)), &
          ' maximum ', Nint(stats%passpmf(5, 1, 1)), ' / ', Nint(stats%passpmf(5, 2, 1))
        Call info(message, .true.)
      End If

      Call gmax(comm, stats%passpmf(3:5, 1, 2)); Call gmax(comm, stats%passpmf(3:5, 2, 2))
      If (stats%passpmf(3, 1, 2) > 0.0_wp) Then
        Write (message, '(2(a,f5.2),4(a,i3))') &
          '# PMFs rattle  run statistics - cycles per call/timestep: average ', &
          stats%passpmf(3, 1, 2), ' / ', stats%passpmf(3, 2, 2), &
          ' minimum ', Nint(stats%passpmf(4, 1, 2)), ' / ', Nint(stats%passpmf(4, 2, 2)), &
          ' maximum ', Nint(stats%passpmf(5, 1, 2)), ' / ', Nint(stats%passpmf(5, 2, 2))
        Call info(message, .true.)
      End If
    End If

    ! Get elapsed time

    Call gtime(timelp)

    ! Get simulation time for averages

    If (stats%numacc == 0) Then
      tmp = 0.0_wp
    Else
      tmp = time - tmst
    End If

    ! Report termination

    If ((nstep == 0 .and. nstrun == 0) .or. stats%numacc == 0) Then
      Write (message, '(a)') '# dry run terminated'
    Else
      If (stats%dpd_units) Then
        Write (message, '(2(a,i9,a,f10.3),a)') '# run terminated after ', nstep, &
          ' steps (', time, ' dpd_t), final averages calculated over', stats%numacc, &
          ' steps (', tmp, ' dpd_t)'
      Else
        Write (message, '(2(a,i9,a,f10.3),a)') '# run terminated after ', nstep, &
          ' steps (', time, ' ps), final averages calculated over', stats%numacc, &
          ' steps (', tmp, ' ps)'
      End If
    End If
    Call info(message, .true.)

    ! safe average volume and cell

    avvol = config%volm

    ! If dry/static/minimisation run - NO AVERAGES
    ! Print pressure tensor and jump to possible RDF and Z-Density

    If (nstep == 0 .and. nstrun == 0) Then
      !iadd = 27 + 2 * Merge(mxatdm, 0, lmsd) + sites%ntype_atom
      iadd = 27

      If (comm%idnode == 0) Then
        If (stats%dpd_units) Then
          Write (message, '(a)') 'Pressure tensor  (katms):'
        Else
          Write (message, '(a)') 'Pressure tensor  (dpd_p):'
        End If
        Call info(message, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4)') stats%stpval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%stpval(iadd + 1) + &
                                                         stats%stpval(iadd + 5) + stats%stpval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
      End If

      If (thermo%variable_cell) Then
        Write (messages(1), '(a)') NEW_LINE('A')//'Strain tensor  (angstroms): '
        Do i = 0, 2
          Write (messages(2+i), '(2x,1p,3e12.4)') stats%strain(1+i), stats%strain(4+i), stats%strain(7+i)
        End Do
        Write(messages(5), '(2x,a,1p,e12.4,a)') 'trace/3  ', Sum(stats%strain(1:9:4)) / 3.0_wp, NEW_LINE('A')
        Call info(messages, 5, .true.)
      End If

      Call gtime(timelp)

      Write (message, '("time elapsed since job start: ", f12.3, " sec")') timelp
      Call info(message, .true.)
      Return
    End If

    ! If still running in the pure equilibration regime - NO AVERAGES
    If (stats%numacc /= 0) Then
      ! shift back statistical averages as from statistics_collect

      Do i = 0, stats%mxnstk
        stats%sumval(i) = stats%sumval(i) + stats%stpvl0(i)
      End Do

      ! calculate final fluctuations

      Do i = 0, stats%mxnstk
        stats%ssqval(i) = Sqrt(stats%ssqval(i))
      End Do

      ! average volume

      avvol = stats%sumval(19)

      ! final averages and fluctuations
      Call write_header(stats%dpd_units)

      Write (messages(1), '("#",i13,1p,9e12.4)') stats%numacc, stats%sumval(1:9)
      Write (messages(2), '("#",f13.5,1p,9e12.4)') tmp, stats%sumval(10:18)
      Write (messages(3), '("#",0p,f13.3,1p,9e12.4)') timelp, stats%sumval(19:27)
      Write (messages(4), '("#",a)') ''
      Call info(messages, 4, .true.)

      Write (messages(1), '("#",5x,a8,1p,9e12.4)') ' r.m.s. ', stats%ssqval(1:9)
      Write (messages(2), '("#",5x,a8,1p,9e12.4)') 'fluctu- ', stats%ssqval(10:18)
      Write (messages(3), '("#",5x,a8,1p,9e12.4)') 'ations  ', stats%ssqval(19:27)
      Write (messages(4), '("#",a)') Repeat('-', 130)
      Call info(messages, 4, .true.)

      ! Some extra information - conserved quantity=extended ensemble energy

      Write (messages(1), "(a)") "Extended energy:"
      Write (messages(2), "(a,e12.4)") " average: ", stats%sumval(0)
      Write (messages(3), "(a,e12.4)") " r.m.s. fluctuations:  ", stats%ssqval(0)
      Call info(messages, 3, .true.)

      ! Some extra information - <P*V> term - only matters for NP/sT ensembles

      strend = Merge (37, 73, thermo%key_dpd==DPD_NULL)

      If (thermo%variable_cell) Then
        Write (messages(1), "(a)") NEW_LINE('a')//"<P*V> term:"
        Write (messages(2), "(a,e12.4)") " average: ", stats%sumval(strend + sites%ntype_atom + 2 * Merge(mxatdm, 0, lmsd))
        Write (messages(3), "(a,e12.4)") " r.m.s. fluctuations:  ", &
          stats%ssqval(strend + sites%ntype_atom + 2 * Merge(mxatdm, 0, lmsd))
        Call info(messages, 3, .true.)
      End If

      Write (messages(1), "('#',130('-'))")
      Write (messages(2), '(a)') ''
      Call info(messages, 2, .true.)

      ! Move at the end of the default 27 quantities

      iadd = 27

      ! print out average pressure tensor

      If (comm%idnode == 0) Then
        Write (messages(1), '(a)') 'Pressure tensor:'
        If (stats%dpd_units) Then
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (dpd_p)', 'r.m.s. fluctuations'
        Else
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)', 'r.m.s. fluctuations'
        End If
        Call info(messages, 2, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4,5x,3e12.4)') stats%sumval(i + 1:i + 3), stats%ssqval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%sumval(iadd + 1) + &
                                                         stats%sumval(iadd + 5) + stats%sumval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
        Call info('', .true.)
      End If

      ! print out the average strain tensor
      If (thermo%variable_cell .and. comm%idnode == 0) Then
        Write (messages(1), '(a)') NEW_LINE('A')//'Strain tensor: '
        Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average (angstroms)', 'r.m.s. fluctuations'
        Do i = 0, 2
          Write (messages(3+i), '(2x,1p,3e12.4,5x,3e12.4)') stats%strain_accum(1+i)%mu, stats%strain_accum(4+i)%mu, &
            stats%strain_accum(7+i)%mu, Sqrt(stats%strain_accum(1+i)%var), Sqrt(stats%strain_accum(4+i)%var), &
            Sqrt(stats%strain_accum(7+i)%var)
        End Do
        Call info(messages, 5, .true.)
        Write(message, '(2x,a,1p,e12.4,a)') 'trace/3  ', &
          (stats%strain_accum(1)%mu+stats%strain_accum(5)%mu+stats%strain_accum(9)%mu) / 3.0_wp, NEW_LINE('A')
        Call info(message, .true.)
      End If

      iadd = iadd + 9

      ! print out the separate contributions to pressure tensor if using DPD

      If (thermo%key_dpd/=DPD_NULL .and. comm%idnode == 0) Then
        Write (messages(1), '(a)') 'Pressure tensor (conservative contributions):'
        If (stats%dpd_units) Then
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (dpd_p)', 'r.m.s. fluctuations'
        Else
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)', 'r.m.s. fluctuations'
        End If
        Call info(messages, 2, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4,5x,3e12.4)') stats%sumval(i + 1:i + 3), stats%ssqval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%sumval(iadd + 1) + &
                                                         stats%sumval(iadd + 5) + stats%sumval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
        Call info('', .true.)
        iadd = iadd + 9
        Write (messages(1), '(a)') 'Pressure tensor (dissipative contributions):'
        If (stats%dpd_units) Then
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (dpd_p)', 'r.m.s. fluctuations'
        Else
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)', 'r.m.s. fluctuations'
        End If
        Call info(messages, 2, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4,5x,3e12.4)') stats%sumval(i + 1:i + 3), stats%ssqval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%sumval(iadd + 1) + &
                                                         stats%sumval(iadd + 5) + stats%sumval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
        Call info('', .true.)
        iadd = iadd + 9
        Write (messages(1), '(a)') 'Pressure tensor (random contributions):'
        If (stats%dpd_units) Then
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (dpd_p)', 'r.m.s. fluctuations'
        Else
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)', 'r.m.s. fluctuations'
        End If
        Call info(messages, 2, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4,5x,3e12.4)') stats%sumval(i + 1:i + 3), stats%ssqval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%sumval(iadd + 1) + &
                                                         stats%sumval(iadd + 5) + stats%sumval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
        Call info('', .true.)
        iadd = iadd + 9
        Write (messages(1), '(a)') 'Pressure tensor (kinetic contributions):'
        If (stats%dpd_units) Then
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (dpd_p)', 'r.m.s. fluctuations'
        Else
          Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)', 'r.m.s. fluctuations'
        End If
        Call info(messages, 2, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4,5x,3e12.4)') stats%sumval(i + 1:i + 3), stats%ssqval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%sumval(iadd + 1) + &
                                                         stats%sumval(iadd + 5) + stats%sumval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
        Call info('', .true.)
        iadd = iadd + 9
      Else If (thermo%key_dpd/=DPD_NULL .and. comm%idnode /= 0) Then
        iadd = iadd + 36
      End If
      
      If (lmsd) iadd = iadd + 2 * mxatdm

      ! Write out estimated diffusion coefficients

      Write (messages(1), '(a)') 'Approximate 3D Diffusion Coefficients and square root of MSDs:'
      If (stats%dpd_units) Then
        Write (messages(2), '(6x,a4,3x,a18,4x,a17)') 'atom', 'DC (dpd_l^2/dpd_t)', 'Sqrt[MSD] (dpd_l)'
      Else
        Write (messages(2), '(6x,a4,2x,a19,6x,a15)') 'atom', 'DC (10^-9 m^2 s^-1)', 'Sqrt[MSD] (Ang)'
      End If
      Call info(messages, 2, .true.)

      Do i = 1, sites%ntype_atom
        If (sites%num_type_nf(i) > zero_plus) Then
          dc = Merge(1.0_wp, 10.0_wp, stats%dpd_units) * (stats%ravval(iadd + i) - stats%sumval(iadd + i)) / &
               (3.0_wp * Real(stats%numacc - Min(stats%mxstak, stats%numacc - 1), wp) * thermo%tstep)
          If (dc < 1.0e-10_wp) dc = 0.0_wp

          srmsd = Sqrt(stats%ravval(iadd + i))
          Write (message, '(2x,a8,1p,2(8x,e13.4))') sites%unique_atom(i), dc, srmsd
        Else
          Write (message, '(2x,a8,1p,2(8x,e13.4))') sites%unique_atom(i), 0.0_wp, 0.0_wp
        End If
        Call info(message, .true.)
      End Do
      Call info('', .true.)

      iadd = iadd + sites%ntype_atom

      ! Write out mean cell vectors for npt/nst

      If (thermo%variable_cell) Then

        If (comm%idnode == 0) Then
          If (stats%dpd_units) Then
            Write (message, '(a32,33x,a19)') 'Average cell vectors    (dpd_l) ', 'r.m.s. fluctuations'
          Else
            Write (message, '(a32,33x,a19)') 'Average cell vectors     (Angs) ', 'r.m.s. fluctuations'
          End If
          Call info(message, .true.)

          Do i = iadd, iadd + 6, 3
            Write (message, '(3f20.10,5x,1p,3e12.4)') stats%sumval(i + 1:i + 3), stats%ssqval(i + 1:i + 3)
            Call info(message, .true.)
          End Do
        End If

        iadd = iadd + 9

        ! PV term used above

        iadd = iadd + 1

        If (thermo%iso /= CONSTRAINT_NONE) Then
          h_z = stats%sumval(iadd + 1)

          If (stats%dpd_units) Then
            Write (message, "('Average surface area, fluctuations & mean estimate (dpd_l^2)')")
          Else
            Write (message, "('Average surface area, fluctuations & mean estimate (Angs^2)')")
          End If
          Call info(message, .true.)
          Write (message, '(1p,3e12.4)') stats%sumval(iadd + 2), stats%ssqval(iadd + 2), avvol / h_z
          Call info(message, .true.)

          iadd = iadd + 2

          If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
            tx = -h_z * (stats%sumval(29) / prsunt0 - (thermo%press + thermo%stress(1))) * tenunt0
            ty = -h_z * (stats%sumval(30) / prsunt0 - (thermo%press + thermo%stress(5))) * tenunt0
            If (stats%dpd_units) Then
              Write (message, "('Average surface tension, fluctuations & mean estimate in x (dpd_f/dpd_l)')")
            Else
              Write (message, "('Average surface tension, fluctuations & mean estimate in x (dyn/cm)')")
            End If
            Call info(message, .true.)
            Write (message, '(1p,3e12.4)') stats%sumval(iadd + 1), stats%ssqval(iadd + 1), tx
            Call info(message, .true.)
            If (stats%dpd_units) Then
              Write (message, "('Average surface tension, fluctuations & mean estimate in y (dpd_f/dpd_l)')")
            Else
              Write (message, "('Average surface tension, fluctuations & mean estimate in y (dyn/cm)')")
            End If
            Call info(message, .true.)
            Write (message, '(1p,3e12.4)') stats%sumval(iadd + 2), stats%ssqval(iadd + 2), ty
            Call info(message, .true.)

            iadd = iadd + 2
          End If
        End If

      End If

      ! Write out remaining registers

      check = .false.
      Do i = iadd + 1, stats%mxnstk
        If (Abs(stats%sumval(i)) > zero_plus .or. Abs(stats%ssqval(i)) > zero_plus) check = .true.
      End Do

      If (check) Then
        Write (messages(1), "('Remaining non-zero statistics registers:')")
        Write (messages(2), "(4x,'Register',7x,'Average value',8x,'r.m.s. fluc.')")
        Call info(messages, 2, .true.)
      End If

      If (comm%idnode == 0) Then
        Do i = iadd + 1, mxnstk
          If (Abs(stats%sumval(i)) > zero_plus .or. Abs(stats%ssqval(i)) > zero_plus) Then
            Write (message, '(2x,i10,2f20.10)') i, stats%sumval(i), stats%ssqval(i)
            Call info(message, .true.)
          End If
        End Do
      End If

    End If

    If (stats%calculate_correlations) Then
      Call correlation_result(stats,rigid,comm,files,config,sites, thermo%tstep)
    End If
    ! print final time check

    Call gtime(timelp)

    Write (message, '("time elapsed since job start: ", f12.3, " sec")') timelp
    Call info(message, .true.)

  End Subroutine statistics_result

  Pure Function calculate_stress(r, f)

    Real(Kind=wp), Dimension(3), Intent(In   ) :: r, f
    Real(Kind=wp), Dimension(9)                :: calculate_stress

    calculate_stress(1:9:3) = r * f(1)
    calculate_stress(2:9:3) = r * f(2)
    calculate_stress(3:9:3) = r * f(3)

  End Function calculate_stress

  Function calculate_strain(stats, config) Result(strain)
    !!----------------------------------------------------------------------!
    !!
    !! Calculate the strain tensor from the scaling matrix,
    !!
    !! \eps = 1/2 * (H_0^(-1, T) H^T H H_0^-1 - I)
    !! 
    !! H_0 (ref) is expected to be the ensemble average of
    !! the scaling matrix H, or the reference box.
    !!
    !! CF G. Clavier et al., 2017, Molecular Simulation
    !!
    !! author    - h.l.devereux jun 2024
    !!
    !!----------------------------------------------------------------------!
    Class(stats_type),                              Intent(InOut) :: stats
    Class(configuration_type),                      Intent(In   ) :: config
    Real(Kind=wp), Dimension(1:3, 1:3) :: eps, h
    Real(Kind=wp)                      :: det, strain(1:9)

    Call scaling_matrix(config%cell, h)
    
    eps = MatMul(Transpose(stats%inv_ref_scaling_matrix), MatMul(Transpose(h), MatMul(h, stats%inv_ref_scaling_matrix)))
    eps(1, 1) = eps(1, 1) - 1.0_wp
    eps(2, 2) = eps(2, 2) - 1.0_wp
    eps(3, 3) = eps(3, 3) - 1.0_wp
    eps = eps * 0.5_wp
    strain = Reshape(eps, (/9/))
  End Function

  Subroutine update_stress(t, s)
    Class(stats_type)            :: t
    Real(kind=wp), Intent(In   ) :: s(9)

    t%stress(1) = t%stress(1) + s(1)
    t%stress(2) = t%stress(2) + s(2)
    t%stress(3) = t%stress(3) + s(3)
    t%stress(4) = t%stress(4) + s(2)
    t%stress(5) = t%stress(5) + s(4)
    t%stress(6) = t%stress(6) + s(5)
    t%stress(7) = t%stress(7) + s(3)
    t%stress(8) = t%stress(8) + s(5)
    t%stress(9) = t%stress(9) + s(6)

  End Subroutine update_stress

  Function calculate_mom_density(stats, atype, config, comm) Result(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for calculating momentum density for a given atom
    ! type.
    !
    !
    ! author    - h.l.devereux Nov 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Use comms, Only: gsum
    Type(stats_type),         Intent(In   ) :: stats
    Integer,                  Intent(In   ) :: atype
    Type(configuration_type), Intent(In   ) :: config
    Type(comms_type),         Intent(InOut) :: comm
    Real(Kind=wp), Dimension(3)             :: j
    
    Integer :: iatm

    j = 0.0_wp
    Do iatm = 1, config%natms
      If (config%ltype(iatm) == atype) Then
        j = j + config%weight(iatm) * [config%vxx(iatm), config%vyy(iatm), config%vzz(iatm)]
      End If
    End Do
    Call gsum(comm, j)
    j = j / stats%stpvol
  End Function

  Function calculate_heat_flux(stats, config, comm) Result(heat_flux)
    Use comms, Only: gsum
    Type(stats_type),         Intent(In   ) :: stats
    Type(configuration_type), Intent(In   ) :: config
    Type(comms_type),         Intent(InOut) :: comm
    Real(Kind=wp), Dimension(3)             :: heat_flux

    Integer                     :: iatm
    Real(Kind=wp), Dimension(3) :: e_v, S_v, velocity

    !! Per-particle energy * velocity
    !! Per-particle stress * velocity
    e_v = 0.0_wp
    S_v = 0.0_wp
    Do iatm = 1, config%natms
      velocity = [config%vxx(iatm), config%vyy(iatm), config%vzz(iatm)]
      !      Σ    (        P              +                                  K                           ) *     V
      e_v = e_v + (stats%pp_energy(iatm) + 0.5_wp * config%weight(iatm) * Dot_product(velocity, velocity)) * velocity
      S_v = S_v + Matmul(Reshape(stats%pp_stress(:, iatm), [3, 3]), velocity)
    End Do
    Call gsum(comm, e_v)
    Call gsum(comm, S_v)

    heat_flux = (e_v + S_v) / (engunit * config%volm)

  End Function calculate_heat_flux

  Subroutine calculate_stress_energy_current(stats, config, iatm, jatm, rij, r_rsq, gamma)
    Class(stats_type),        Intent(InOut) :: stats
    Type(configuration_type), Intent(In   ) :: config
    Integer,                  Intent(In   ) :: iatm, jatm
    Real(Kind=wp),            Intent(In   ) :: r_rsq, gamma, rij(1:3)

    Complex(Kind=wp) :: ikdotrij, vir_pre, pkij, cur_str(1:6), cur_vir(1:3)
    Real(Kind=wp)    :: vi(3), vj(3)
    Integer          :: b, kpoint

    If (.not. stats%cur%k_energy_stress_current_on) Return

    vi = (/config%vxx(iatm), config%vyy(iatm), config%vzz(iatm)/)
    vj = (/config%vxx(jatm), config%vyy(jatm), config%vzz(jatm)/)
    Do kpoint = 1, config%k%n
      cur_vir = Cmplx(0.0_wp, 0.0_wp, Kind=wp)
      cur_str = Cmplx(0.0_wp, 0.0_wp, Kind=wp)
      ikdotrij = Cmplx(0.0_wp, 1.0_wp, Kind=wp) * Dot_product(config%k%r(:, kpoint), rij)
      If (ikdotrij /= Cmplx(0.0_wp, 0.0_wp, Kind=wp)) Then
        pkij = (1-Exp(-1.0_wp*ikdotrij))/(ikdotrij)
        vir_pre = r_rsq*gamma*pkij
        Do b = 1, 3
          cur_vir = cur_vir + vir_pre*(vi(b)+vj(b))*rij(b)
        End Do

        cur_str(1) = wi*(vi(1)**2)-0.5*(rij(1)**2*r_rsq)*pkij
        cur_str(2) = wi*(vi(1)*vi(2))-0.5*(rij(1)*rij(2)*r_rsq)*pkij
        cur_str(3) = wi*(vi(1)*vi(3))-0.5*(rij(1)*rij(3)*r_rsq)*pkij
        cur_str(4) = wi*(vi(2)**2)-0.5*(rij(2)**2*r_rsq)*pkij
        cur_str(5) = wi*(vi(2)*vi(3))-0.5*(rij(1)*rij(3)*r_rsq)*pkij
        cur_str(6) = wi*(vi(3)**2)-0.5*(rij(3)**2*r_rsq)*pkij
      End If
      stats%pp_cur_virial(iatm, kpoint, :) = stats%pp_cur_virial(iatm, kpoint, :) +&
        cur_vir
      stats%pp_cur_stress(iatm, kpoint, :) = stats%pp_cur_stress(iatm, kpoint, :) +&
        cur_str
    End Do
  End Subroutine calculate_stress_energy_current

  Subroutine calculate_viscosity(stats, dt, viscosity)
    Type(stats_type),  Intent(InOut)              :: stats
    Real(Kind=wp),     Intent(In   )              :: dt
    Real(Kind=wp),     Intent(  Out), Allocatable :: viscosity(:)

    Real(Kind=wp),     Allocatable  :: correlation(:)
    Class(integrator), Allocatable  :: inter
    Integer                         :: i, freq

    Character(Len=2), Dimension(1:6), Parameter :: components = (/"xy", "xz", "yx", "yz", "zx", "zy"/)
    Real(Kind=wp)                  :: boltz0

    boltz0 = Merge(1.0_wp, boltz, stats%dpd_units)

    Allocate(simpsons_rule::inter)
    Allocate(viscosity(0))

    Do i = 1, Size(components)
      Call get_correlation_value(stats, "stress_"//components(i)//"-stress_"//components(i), correlation)
      If (Allocated(correlation)) Then
        freq = get_correlation_frequency(stats, "stress_"//components(i)//"-stress_"//components(i))
        viscosity = [viscosity, inter%integrate_uniform(correlation(:), dt*freq)]
        Deallocate(correlation)
      End If
    End Do

    If (Size(viscosity) > 0) Then
      viscosity = prsunt * ( stats%accumulators(19)%mu / (boltz0*stats%accumulators(2)%mu) ) * viscosity
    Else
      Deallocate(viscosity)
    End If

  End Subroutine calculate_viscosity

  Subroutine init_born_calculate(stats, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for checking which born terms to calculate in
    ! vdw. If a user correlation conmensurate with it is present, that
    ! born tensor entry will be calculated as well.
    !
    ! author    - h.l.devereux July 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type),              Intent(InOut) :: stats
    Type(comms_type),               Intent(InOut) :: comm

    Character(Len=1) :: symbols(1:3) = (/"x", "y", "z"/), a, b, c, d
    Integer          :: v
    If (comm%idnode == root_id) Then
      Do v = 1, 21
        a = symbols(voigt_6x6(v, 1))
        b = symbols(voigt_6x6(v, 2))
        c = symbols(voigt_6x6(v, 3))
        d = symbols(voigt_6x6(v, 4))
        If (stats%cor_table%contains("stress_"//a//b//"-stress_"//c//d)) Then
          stats%born_calculate(v) = .true.
        End If
      End Do
    End If
    Call gbcast(comm, stats%born_calculate, root_id)
  End Subroutine init_born_calculate

  Subroutine elastic_constants_result(stats, megatm, file_unit)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for writing out elastic constants. Up to 21
    !   (independent) components are possible. All possible components
    !   calculable from user stress correlations are outputted in Voigt
    !   order (see constants.F90:voigt_6x6): C1111, C1122, C1133, ..., 
    !   C1212.
    !
    !   Stress fluctuation method: e.g. G. Clavier, et al., 
    !   Molecular Simulation, 2017, 
    !   https://doi.org/10.1080/08927022.2017.1313418
    !
    ! author    - h.l.devereux July 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stats_type),              Intent(InOut) :: stats
    Integer,                       Intent(In   ) :: megatm
    Integer,                       Intent(In   ) :: file_unit

    Character(Len=1)               :: symbols(1:3) = (/"x", "y", "z"/), a, b, c, d
    Integer                        :: v, i, j, k, l
    Real(Kind=wp)                  :: stress_prefactor, kinetic_term, cijkl, del
    Real(Kind=wp),    Allocatable  :: correlation(:), elasticity(:)
    Character(Len=6), Allocatable  :: component_names(:)           

    stress_prefactor = prsunt * stats%accumulators(19)%mu/(boltz*stats%accumulators(2)%mu) ! V / (kbT)
    kinetic_term = prsunt * 2.0_wp*boltz*stats%accumulators(2)%mu * megatm / stats%accumulators(19)%mu ! NKbT/V
    Allocate(elasticity(0))
    ! Debug crashes if component_names is allocated as 0
    ! The Character(Len=6) is left undefined.
    Allocate(component_names(1))
    Do v = 1, 21
      i = voigt_6x6(v, 1)
      a = symbols(i)
      j = voigt_6x6(v, 2)
      b = symbols(j)
      k = voigt_6x6(v, 3)
      c = symbols(k)
      l = voigt_6x6(v, 4)
      d = symbols(l)
      Call get_correlation_value(stats, "stress_"//a//b//"-stress_"//c//d, correlation)
      If (Allocated(correlation)) Then
        del = 0.0_wp
        If (i == k .and. j == l) del = del + 1.0_wp
        If (i == l .and. j == k) del = del + 1.0_wp
        cijkl = prsunt*stats%born_term_accum(v)%mu/ stats%accumulators(19)%mu &
              - stress_prefactor*(correlation(1) - &
              stats%stress_accum((i-1)*3+j)%mu*stats%stress_accum((k-1)*3+l)%mu)+ &
              kinetic_term * del
        elasticity = [elasticity, cijkl]
        component_names = [component_names, Trim("C_"//a//b//c//d)]
        Deallocate(correlation)
      End If
    End Do
    If (Size(elasticity) > 0) Then
      Write (file_unit, '(a)') "      elasticity_tensor:"
      Call write_char_yaml_vector(file_unit, "components", component_names(2:), 12)
      Call write_real_yaml_vector(file_unit, "values", elasticity, 12)
      Write (file_unit, '(a)') "            units: Katm"
    End If

  End Subroutine

  Subroutine calculate_thermal_conductivity(stats, dt, units, therm_cond)
    Type(stats_type),       Intent(InOut)              :: stats
    Real(Kind=wp),          Intent(In   )              :: dt
    Character(Len=STR_LEN), Intent(  Out)              :: units
    Real(Kind=wp),          Intent(  Out), Allocatable :: therm_cond(:)
    
    Real(Kind=wp)                   :: conv, boltz0
    Real(Kind=wp),     Allocatable  :: correlation(:)
    Class(integrator), Allocatable  :: inter
    Integer                         :: i, freq

    Character(Len=1), Dimension(1:3), Parameter :: components = (/"x", "y", "z"/)
    
    boltz0 = Merge(1.0_wp, boltz, stats%dpd_units)

    Allocate(simpsons_rule::inter)
    Allocate(therm_cond(0))

    Do i = 1, Size(components)
      Call get_correlation_value(stats, "heat_flux_"//components(i)//"-heat_flux_"//components(i), correlation)
      If (Allocated(correlation)) Then
        freq = get_correlation_frequency(stats, "heat_flux_"//components(i)//"-heat_flux_"//components(i))
        therm_cond = [therm_cond, inter%integrate_uniform(correlation(:), dt*freq)]
        Deallocate(correlation)
      End If
    End Do

    If (Size(therm_cond) > 0) Then
      Call to_out_units(1.0_wp, "internal_e", conv, units)
      ! already divided through by volume
      therm_cond = stats%accumulators(19)%mu / ( (stats%accumulators(2)%mu**2) * (boltz0/engunit)) * therm_cond
      If (stats%dpd_units) Then
        units = Trim(units)//" / (dpd_t dpd_l dpd_temp)"
      Else
        units = Trim(units)//" / (ps Ang K)"
      End If
    Else
      Deallocate(therm_cond)
    End If
  End Subroutine calculate_thermal_conductivity

  Subroutine write_per_part_contribs(config, comm, energies, stresses, nstep) !, forces
    !!----------------------------------------------------------------------!
    !!
    !! Write out per-particle contributions to energy, force, stress, etc
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Type(configuration_type),        Intent(In   ) :: config
    Type(comms_type),                Intent(InOut) :: comm
    Real(Kind=wp), Dimension(1:),    Intent(In   ) :: energies
    Real(Kind=wp), Dimension(:, 1:), Intent(In   ) :: stresses
    Integer,                         Intent(In   ) :: nstep

    Integer, Parameter :: record_size = 73

    Character                                :: lf
    Character(len=40)                        :: filename
    Character(len=record_size)               :: record
    Character, Dimension(record_size, 10)    :: buffer
    Integer                                  :: batsz, energy_force_handle, i, ierr, io_write, jj
    Integer(Kind=offset_Kind)                :: rec_mpi_io
    Real(Kind=wp)                            :: velocity(3)
    Real(Kind=wp), Allocatable, Dimension(:) :: dummy
    Type(io_type)                            :: my_io

!! Atom details
!! Communicator
!! Per-particle energies
! Real( Kind = wp ), Dimension(:,1:),   Intent ( In    )  :: forces    !!     ""       forces
!!     ""       stresses
!! Steps since calculation start
!! Use our own IO job for now because passing through will be hell
!! Don't like this, but quick cheat?
!! default record size (apparently)
!! File handles
!! Write state

    Call gsync(comm)

    ! Force MPIIO write for now
    io_write = 0
    ! Call io_get_parameters( user_method_write      = io_write )
    Call io_get_parameters(my_io, user_buffer_size_write=batsz, user_line_feed=lf)

    ! Write current time-step to character string
    Allocate (dummy(config%natms), stat=ierr)
    If (ierr .ne. 0) Call error_alloc('dummy', 'write_per_part_contribs')
    dummy = engunit

    Write (filename, '("PPCONT",("_",i0))') nstep

    Call io_init(my_io, record_size)

    rec_mpi_io = Int(0, offset_Kind)
    jj = 0
    If (comm%idnode == 0) Then

      Call io_set_parameters(my_io, user_comm=comm_self)
      Call io_delete(my_io, filename, comm) ! sort existence issues
      Call io_open(my_io, io_write, comm_self, Trim(filename), mode_wronly + mode_create, energy_force_handle)

      jj = jj + 1
      Write (record, Fmt='(a72,a1)') "Energy and force contributions on a per-particle basis", lf
      buffer(:, jj) = [(record(i:i), i=1, record_size)]
      Write (record, Fmt='(a72,a1)') config%cfgname(1:72), lf
      buffer(:, jj) = [(record(i:i), i=1, record_size)]
      jj = jj + 1
      Write (record, Fmt='(3i10,42X,a1)') config%imcon, config%megatm, nstep, lf
      buffer(:, jj) = [(record(i:i), i=1, record_size)]

      If (config%imcon > 0) Then
        Do i = 0, 2
          jj = jj + 1
          Write (record, Fmt='(3f20.10,a12,a1)') &
            config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3), Repeat(' ', 12), lf
          buffer(:, jj) = [(record(i:i), i=1, record_size)]
        End Do
      End If

      Call io_write_batch(my_io, energy_force_handle, rec_mpi_io, jj, buffer)

      Call io_close(my_io, energy_force_handle)

    End If

    Call gsync(comm)

    Do i = 1, config%natms
      velocity = [config%vxx(i), config%vyy(i), config%vzz(i)]
      dummy(i) = 0.5_wp * config%weight(i) * Dot_product(velocity, velocity) / engunit
    End Do

    Call io_set_parameters(my_io, user_comm=comm%comm)
    Call io_open(my_io, io_write, comm%comm, Trim(filename), mode_wronly, energy_force_handle) ! Io sorted mpiio, per-particle contrib

    rec_mpi_io = Int(jj, offset_Kind)
    ! Only write E&F (r/v in write_sorted...) hence 1
    ! Need to skip 0th element (accumulator/total)
    Call io_write_sorted_file(my_io, energy_force_handle, 2, io_history, rec_mpi_io, config%natms, &
      config%ltg, config%atmnam, dummy, config%weight(1:config%natms), energies(1:config%natms) / engunit, &
      & stresses(1, 1:config%natms) * prsunt, stresses(2, 1:config%natms) * prsunt, stresses(3, 1:config%natms) * prsunt, &
      & stresses(4, 1:config%natms) * prsunt, stresses(5, 1:config%natms) * prsunt, stresses(6, 1:config%natms) * prsunt, &
      & stresses(7, 1:config%natms) * prsunt, stresses(8, 1:config%natms) * prsunt, stresses(9, 1:config%natms) * prsunt, ierr)
    ! forces(1,1:config%natms), forces(2,1:config%natms), forces(3,1:config%natms), &

    Select Case (ierr)
    Case (0)
      Continue
    Case (io_base_comm_not_set)
      Call error(1050)
    Case (io_allocation_error)
      Call error(1053)
    Case (io_unknown_write_option)
      Call error(1056)
    Case (io_unknown_write_level)
      Call error(1059)
    End Select
    Call io_close(my_io, energy_force_handle)

    Call gsync(comm)

    Call io_finalize(my_io)

    Deallocate (dummy, stat=ierr)
    If (ierr > 0) Call error_dealloc('dummy', 'write_per_part_contribs')

  End Subroutine write_per_part_contribs

  Subroutine get_correlation_value(stats, name, correlation, atom)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for getting a correlation from a
    ! given name (see get_correlation_index). correlation is
    ! unallocated if the correlation was not found.
    !
    ! author    - h.l.devereux 2023
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stats_type),             Intent(InOut) :: stats
    Character(Len=*),             Intent(In   ) :: name
    Real(Kind=wp),   Allocatable, Intent(InOut) :: correlation(:)
    Integer,         Optional,    Intent(In   ) :: atom

    Integer                    :: cor_index, atom_index
    Integer                    :: points, blocks
    Real(Kind=wp), Allocatable :: timesteps(:)

    If (stats%cor_table%contains(name)) Then
      Call stats%cor_table%get(name, cor_index)
      atom_index = 1
      If (Present(atom)) Then
        atom_index = atom
      End If
      points = stats%correlations(cor_index)%correlators(atom_index)%points_per_block
      blocks = stats%correlations(cor_index)%correlators(atom_index)%number_of_blocks
      If (points-1 > 1) Then
        Allocate(correlation(1:points*blocks))
        Allocate(timesteps(1:points*blocks))
        Call stats%correlations(cor_index)%correlators(atom_index)%get_correlation(correlation, timesteps, points)
      End If
    End If
  End Subroutine

  Integer Function get_correlation_frequency(stats, name)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 function to get a correlations update frequency.
    !
    ! author    - h.l.devereux August 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stats_type),             Intent(InOut) :: stats
    Character(Len=*),             Intent(In   ) :: name

    Integer :: cor_index    

    If (stats%cor_table%contains(name)) Then
      Call stats%cor_table%get(name, cor_index)
      get_correlation_frequency = stats%correlations(cor_index)%freq
    Else
      Call error(0, "Cannot get frequency for non-existent correlation "//name)
    End If
  End Function get_correlation_frequency

  Subroutine correlator_reindex(stats, config, rigid)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 for updating correlator local indices for per-atom
    ! correlations after deport/receipt.
    !
    ! author    - h.l.devereux 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type),                       Intent(InOut)  :: stats
    Type(configuration_type),                Intent(InOut)  :: config
    Type(rigid_bodies_type),                 Intent(InOut)  :: rigid

    Integer                         :: i, j, new_index, cor_index
    Type(correlation_data), Pointer :: cor_data
    Character(Len=MAX_KEY), Allocatable, Dimension(:) :: cor_keys

    If (stats%per_atom_correlations .and. stats%calculate_correlations) Then
      Call stats%cor_table%get_keys(cor_keys)
      Do j = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(j))) Cycle
        Call stats%cor_table%get(cor_keys(j), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          If (is_per_atom(cor_data)) Then
            Do i = 1, cor_data%indices_used
              new_index = Findloc(cor_data%indices_global, config%ltg(i), 1)
              cor_data%indices(new_index) = i
            End Do
          Else If (is_per_rigid(cor_data)) Then
            Do i = 1, Size(rigid%list, 2)
              new_index = Findloc(cor_data%indices_global, rigid%list(1, i), 1)
              ! Rigid bodies may be shared. But only need one cor.
              If (new_index /= 0) Then
                cor_data%indices(new_index) = i
              End If
            End Do
          End If
        End Associate
      End Do
    End If
  End Subroutine correlator_reindex

  Subroutine correlator_recieve(stats, buffer, buffer_index, dep_type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 for recieving correlations defined on a per-atom basis
    ! see also correlator_deport.
    !
    ! author    - h.l.devereux 2023
    !
    ! data packed as 
    !   deportations count (if 0, skips) - always written if per-atom correlations calculated
    !   global atom index                - dependent on deportations count
    !   A observable code                         |
    !   A component                               |
    !   B observable code                         |
    !   B component                               |
    !   update_freq                               |
    !   flat correlator data (see correlator.F90) |
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Class(stats_type),                       Intent(InOut)  :: stats
      Real(Kind=wp),           Dimension(:),   Intent(InOut)  :: buffer
      Integer,                                 Intent(InOut)  :: buffer_index
      Integer,                                 Intent(In   )  :: dep_type

      Class(observable), Allocatable  :: A, B
      Integer                         :: iA, iB, gid, &
                                         window, blocks, points, deportations, d, &
                                         new_index, s, c_a, c_b, freq, cor_index, n, &
                                         type_a, type_b
      Character(Len=MAX_KEY)          :: name
      Type(correlation_data), Pointer :: cor_data

      If (stats%per_atom_correlations .and. stats%calculate_correlations) Then

        s = buffer_index
        buffer_index = buffer_index + 1
        deportations = INT(buffer(buffer_index))

        If (deportations < 1) Then
          buffer_index = s+1
          Return
        End If

        Do d = 1, deportations
          buffer_index = buffer_index + 1
          gid = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          iA = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          c_a = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          iB = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          c_b = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          freq = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          type_a = INT(buffer(buffer_index))
          buffer_index = buffer_index + 1
          type_b = INT(buffer(buffer_index))

          Call id_component_to_observable(iA, c_a, A, type_a)
          Call id_component_to_observable(iB, c_b, B, type_b)
          name = Trim(A%name())//"-"//Trim(B%name())
          Call stats%cor_table%get(name, cor_index)
          Associate(cor_data => stats%correlations(cor_index))
            new_index = cor_data%indices_used+1
            If (new_index > Size(cor_data%indices)) Then
              ! A spare entry does not exist, extend by 1.
              cor_data%correlators = [cor_data%correlators, cor_data%correlators(Size(cor_data%correlators))]
              cor_data%indices = [cor_data%indices, cor_data%indices(Size(cor_data%indices))]
              cor_data%indices_global = [cor_data%indices_global, cor_data%indices_global(Size(cor_data%indices_global))]
            End If

            cor_data%indices_global(new_index) = gid
            blocks = cor_data%correlators(1)%number_of_blocks
            points = cor_data%correlators(1)%points_per_block
            window = cor_data%correlators(1)%window_size
            Call cor_data%correlators(new_index)%recieve_buffer(buffer,buffer_index)
            cor_data%freq = freq
            cor_data%indices_used = cor_data%indices_used + 1
          End Associate
        End Do

        buffer_index = s + deportations*stats%cor_deport_buffer + 1

      End If

  End Subroutine correlator_recieve

  Subroutine correlator_deport(stats, buffer, gid, buffer_index, dep_type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 for deporting correlations defined on a per-atom basis
    ! see also correlator_receive.
    !
    ! author    - h.l.devereux 2023
    !
    ! data packed as 
    !   deportations count (if 0, skips) - always written if per-atom correlations calculated
    !   global atom index                - dependent on deportations count
    !   A observable code                         |
    !   A component                               |
    !   B observable code                         |
    !   B component                               |
    !   update frequency                          |
    !   flat correlator data (see correlator.F90) |
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stats_type),                    Intent(InOut) :: stats
    Real(Kind=wp), Dimension(:),          Intent(InOut) :: buffer
    Integer,                              Intent(In   ) :: gid, dep_type
    Integer,                              Intent(InOut) :: buffer_index

    Integer                             :: i, A, B, c_a, c_b, s, n, &
                                           cor_index, location, &
                                           type_a, type_b
    Integer,                Allocatable :: locations(:), cor_indices(:)
    Character(Len=MAX_KEY), Allocatable :: cor_keys(:)
    Type(correlation_data), Pointer     :: cor_data

    If (stats%per_atom_correlations .and. stats%calculate_correlations) Then

      s = buffer_index
      Allocate(locations(0))
      Allocate(cor_indices(0))

      Call stats%cor_table%get_keys(cor_keys)

      Do i = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
        Call stats%cor_table%get(cor_keys(i), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          If (cor_data%A%distributed() /= dep_type) Cycle
          location = Findloc(cor_data%indices_global, gid, 1)
          ! RB entries may be "deported" multiple times due to sharing.
          If (location /= 0 .and. location <= cor_data%indices_used) Then
            locations = [locations, location]
            cor_indices = [cor_indices, cor_index]
          End If
        End Associate
      End Do

      buffer_index = buffer_index + 1 
      buffer(buffer_index) = Size(cor_indices)
      If (Size(cor_indices) == 0) Return

      Do i = 1, Size(cor_indices)
        Associate(cor_data => stats%correlations(cor_indices(i)))
          A = cor_data%A%id()
          c_a = cor_data%A%component
          B = cor_data%B%id()
          c_b = cor_data%B%component
          buffer_index = buffer_index + 1
          buffer(buffer_index) = gid
          buffer_index = buffer_index + 1
          buffer(buffer_index) = A
          buffer_index = buffer_index + 1
          buffer(buffer_index) = c_a
          buffer_index = buffer_index + 1
          buffer(buffer_index) = B
          buffer_index = buffer_index + 1
          buffer(buffer_index) = c_b
          buffer_index = buffer_index + 1
          buffer(buffer_index) = cor_data%freq
          buffer_index = buffer_index + 1
          buffer(buffer_index) = cor_data%A%subtype()
          buffer_index = buffer_index + 1
          buffer(buffer_index) = cor_data%B%subtype()

          location = locations(i)
          Call cor_data%correlators(location)%deport_buffer(buffer,buffer_index,.true.)
          ! remove old correlator
          n = cor_data%indices_used-1
          If (location /= n+1) Then
            ! replace hole with off-end element
            cor_data%correlators(location) = cor_data%correlators(n+1)
            cor_data%indices(location) = cor_data%indices(n+1)
            cor_data%indices_global(location) = cor_data%indices_global(n+1)
          End If
          cor_data%indices_used = cor_data%indices_used - 1
        End Associate
      End Do
      buffer_index = s + Size(cor_indices)*stats%cor_deport_buffer + 1
    End If

  End Subroutine correlator_deport

  !!!!!!!! correlators revive !!!!!!!!!

  Subroutine dump_correlations(stats, comm, config, unit)
    Class(stats_type),              Intent(InOut) :: stats
    Type(comms_type),               Intent(InOut) :: comm
    Type(configuration_type),       Intent(InOut) :: config
    Integer,                        Intent(In   ) :: unit

    Real(Kind=wp),               Allocatable :: data_buffer(:), local_buffer(:)
    Integer,                     Allocatable :: local_ids(:), ids_buffer(:)
    Type(correlator_buffer_type)             :: packed_correlators
    Type(indices_buffer_type)                :: packed_ids
    Integer                                  :: i, atom, idx, buffer_size, correlations, &
                                                buffer_index, local_correlations, local_buffer_size, &
                                                n_local_cor, cor_index, &
                                                attributes = 9, header = 4
    Character(Len=MAX_KEY),      Allocatable :: cor_keys(:)

  
    n_local_cor = stats%number_of_correlations
    Call gsum(comm, n_local_cor)
    If (n_local_cor == 0) Return

    ! determine total buffer sizes needed for root

    buffer_size = 0
    correlations = 0
    local_correlations = 0
    local_buffer_size = 0
    If (stats%calculate_correlations) Then

      Call stats%cor_table%get_keys(cor_keys)

      Do i = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
        Call stats%cor_table%get(cor_keys(i), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          correlations = correlations + cor_data%indices_used
          buffer_size = buffer_size + cor_data%indices_used*(cor_data%correlators(1)%buffer_size-header)
        End Associate
      End Do
      ! id_a, id_b, component_a, component_b, atom, buffer_size
      Allocate(local_ids(1:attributes*correlations))
      Allocate(local_buffer(1:buffer_size))
      ! collect local buffer first
      buffer_index = 0
      idx = 1
      Do i = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
        Call stats%cor_table%get(cor_keys(i), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          Do atom = 1, cor_data%indices_used
            local_ids((idx-1)*attributes+1) = cor_data%A%id()
            local_ids((idx-1)*attributes+2) = cor_data%A%component
            local_ids((idx-1)*attributes+3) = cor_data%B%id()
            local_ids((idx-1)*attributes+4) = cor_data%B%component
            local_ids((idx-1)*attributes+5) = cor_data%freq
            local_ids((idx-1)*attributes+6) = cor_data%A%subtype()
            local_ids((idx-1)*attributes+7) = cor_data%B%subtype()
            If (cor_data%indices(1) /= 0) Then
              local_ids((idx-1)*attributes+8) = cor_data%indices_global(atom)
            Else
              local_ids((idx-1)*attributes+8) = 0
            End If
            local_ids((idx-1)*attributes+9) = cor_data%correlators(1)%buffer_size-header
            Call cor_data%correlators(atom)%deport_buffer(local_buffer,buffer_index)
            idx = idx + 1
          End Do
        End Associate
      End Do
      local_correlations = correlations
      local_buffer_size = buffer_size
    End If

    If (local_correlations == 0) Then
      ! ggatherv will call Size, so must at least
      !   have a dummy allocation
      If (.not. Allocated(local_ids)) Allocate(local_ids(0))
      If (.not. Allocated(local_buffer)) Allocate(local_buffer(0))
    End If

    ! now globally sum sizes
    Call gsum(comm,buffer_size)
    Call gsum(comm,correlations)

    Call packed_ids%initialise(comm, attributes*correlations)
    Call packed_correlators%initialise(comm, buffer_size)

    If (comm%idnode == root_id) Then

      Allocate(data_buffer(1:buffer_size))
      Allocate(ids_buffer(1:attributes*correlations))
      data_buffer = 0.0_wp
      ids_buffer = 0

    End If

    Call gatherv_scatterv_index_arrays(comm, &
    attributes*local_correlations, &
    packed_ids%mpi%counts, &
    packed_ids%mpi%displ & 
    )

    Call ggatherv(comm, local_ids, &
      packed_ids%mpi%counts, &
      packed_ids%mpi%displ, &
      packed_ids%buffer)

    Call gatherv_scatterv_index_arrays(comm, &
      local_buffer_size, &
      packed_correlators%mpi%counts, &
      packed_correlators%mpi%displ & 
    )

    Call ggatherv(comm, local_buffer, &
      packed_correlators%mpi%counts, &
      packed_correlators%mpi%displ, &
      packed_correlators%buffer)

    ! dump on root

    If (comm%idnode == root_id) Then
      buffer_index = 1
      Do i = 1, correlations
        Write (unit) packed_ids%buffer((i-1)*attributes+1), &
          packed_ids%buffer((i-1)*attributes+2), &
          packed_ids%buffer((i-1)*attributes+3), &
          packed_ids%buffer((i-1)*attributes+4), &
          packed_ids%buffer((i-1)*attributes+5), &
          packed_ids%buffer((i-1)*attributes+6), &
          packed_ids%buffer((i-1)*attributes+7), &
          packed_ids%buffer((i-1)*attributes+8), &
          packed_ids%buffer((i-1)*attributes+9), &
          packed_correlators%buffer(buffer_index:(buffer_index+packed_ids%buffer((i-1)*attributes+attributes)-1))
        buffer_index = buffer_index + packed_ids%buffer((i-1)*attributes+attributes)

      End Do

    End If

    Call packed_ids%finalise()
    Call packed_correlators%finalise()

  End Subroutine dump_correlations

  Subroutine revive_correlations(stats, comm, config, unit, keyio, no_advance, format)
    Class(stats_type),              Intent(InOut) :: stats
    Type(comms_type),               Intent(InOut) :: comm
    Type(configuration_type),       Intent(InOut) :: config
    Integer,                        Intent(In   ) :: unit
    Logical,                        Intent(In   ) :: no_advance
    Character(Len=40),              Intent(In   ) :: format
    Integer,                        Intent(InOut) :: keyio

    Real(Kind=wp),               Allocatable :: data_buffer(:), local_buffer(:)
    Integer,                     Allocatable :: local_ids(:), ids_buffer(:), &
                                                offests_buffer(:), sizes_buffer(:)
    Type(correlator_buffer_type)             :: packed_correlators
    Type(indices_buffer_type)                :: packed_ids
    Integer                                  :: i, idx, buffer_size, correlations, &
                                                buffer_index, j, packed_index, &
                                                A, B, atom, offset, local_correlations, &
                                                local_buffer_size, n_local_cor, &
                                                c_a, c_b, cor_index, type_a, type_b, &
                                                attributes = 9, header = 4
    Character(Len=MAX_KEY),      Allocatable :: cor_keys(:)
    Type(correlation_data),      Pointer     :: cor_data
  
    n_local_cor = stats%number_of_correlations
    Call gsum(comm, n_local_cor)
    If (n_local_cor == 0) Return
    local_correlations = 0
    local_buffer_size = 0
    buffer_size = 0
    correlations = 0
    If (stats%calculate_correlations) Then     
      Call stats%cor_table%get_keys(cor_keys)

      Do i = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
        Call stats%cor_table%get(cor_keys(i), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          correlations = correlations + Size(cor_data%indices)
          buffer_size = buffer_size + Size(cor_data%indices)*(cor_data%correlators(1)%buffer_size-header)
        End Associate
      End Do                                            
      Allocate(local_ids(1:attributes*correlations))
      local_correlations = correlations
      ! determine total buffer sizes needed for root

      ! collect local id buffers
      idx = 1
      Do i = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
        Call stats%cor_table%get(cor_keys(i), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          Do atom = 1, Size(cor_data%indices)
            local_ids((idx-1)*attributes+1) = cor_data%A%id()
            local_ids((idx-1)*attributes+2) = cor_data%A%component
            local_ids((idx-1)*attributes+3) = cor_data%B%id()
            local_ids((idx-1)*attributes+4) = cor_data%B%component
            local_ids((idx-1)*attributes+5) = cor_data%freq
            local_ids((idx-1)*attributes+6) = cor_data%A%subtype()
            local_ids((idx-1)*attributes+7) = cor_data%B%subtype()
            If (cor_data%indices(1) /= 0) Then
              local_ids((idx-1)*attributes+8) = cor_data%indices_global(atom)
            Else
              local_ids((idx-1)*attributes+8) = 0
            End If
            local_ids((idx-1)*attributes+9) = cor_data%correlators(1)%buffer_size-header
            idx = idx + 1
          End Do
        End Associate
      End Do
      Allocate(local_buffer(1:buffer_size))
      local_buffer_size = buffer_size
    End If

    If (local_correlations == 0) Then
      ! ggatherv will call Size, so must at least
      !   have a dummy allocation
      Allocate(local_ids(0))
      Allocate(local_buffer(0))
    End If

    Call gsum(comm,buffer_size)
    Call gsum(comm,correlations)

    Call packed_ids%initialise(comm, attributes*correlations)
    Call packed_correlators%initialise(comm, buffer_size)

    If (comm%idnode == root_id) Then

      Allocate(data_buffer(1:buffer_size))
      Allocate(sizes_buffer(1:correlations))
      Allocate(ids_buffer(1:attributes*correlations))
      Allocate(offests_buffer(1:correlations))
    
    End If

    Call gatherv_scatterv_index_arrays(comm, &
    attributes*local_correlations, &
    packed_ids%mpi%counts, &
    packed_ids%mpi%displ & 
    )

    Call ggatherv(comm, local_ids, &
      packed_ids%mpi%counts, &
      packed_ids%mpi%displ, &
      packed_ids%buffer)

    If (comm%idnode == root_id) Then
      ! root reads the data as is
      buffer_index = 1
      Do i = 1, correlations
          offests_buffer(i) = buffer_index
          If (no_advance) Then
            Read (Unit=unit, IOStat=keyio, Fmt=format, Advance = 'No') &
              ids_buffer((i-1)*attributes+1), ids_buffer((i-1)*attributes+2), ids_buffer((i-1)*attributes+3), &
              ids_buffer((i-1)*attributes+4), ids_buffer((i-1)*attributes+5), ids_buffer((i-1)*attributes+6), &
              ids_buffer((i-1)*attributes+7), ids_buffer((i-1)*attributes+8), ids_buffer((i-1)*attributes+9), &
              data_buffer(buffer_index:(buffer_index+ids_buffer((i-1)*attributes+attributes)-1))
          Else
            Read (Unit=unit, IOStat=keyio) &
              ids_buffer((i-1)*attributes+1), ids_buffer((i-1)*attributes+2), ids_buffer((i-1)*attributes+3), &
              ids_buffer((i-1)*attributes+4), ids_buffer((i-1)*attributes+5), ids_buffer((i-1)*attributes+6), &
              ids_buffer((i-1)*attributes+7), ids_buffer((i-1)*attributes+8), ids_buffer((i-1)*attributes+9), &
              data_buffer(buffer_index:(buffer_index+ids_buffer((i-1)*attributes+attributes)-1))
          End If
          sizes_buffer(i) = ids_buffer((i-1)*attributes+attributes)
          buffer_index = buffer_index + sizes_buffer(i)
      End Do

    End If

    If (comm%idnode == root_id) Then

      ! now setup arrays for scattering
      packed_correlators%buffer = data_buffer

      buffer_index = 1

      Do i = 1, correlations

        A = packed_ids%buffer((i-1)*attributes+1)
        c_a = packed_ids%buffer((i-1)*attributes+2)
        B = packed_ids%buffer((i-1)*attributes+3)
        c_b = packed_ids%buffer((i-1)*attributes+4)
        type_a = packed_ids%buffer((i-1)*attributes+6)
        type_b = packed_ids%buffer((i-1)*attributes+7)
        atom = packed_ids%buffer((i-1)*attributes+8)
        packed_index = -1

        ! find where this data should be placed
        Do j = 1, correlations
          If (ids_buffer((j-1)*attributes+1) == A   .and. &
              ids_buffer((j-1)*attributes+2) == c_a .and. &
              ids_buffer((j-1)*attributes+3) == B   .and. &
              ids_buffer((j-1)*attributes+4) == c_b .and. &
              ids_buffer((j-1)*attributes+6) == type_a .and. &
              ids_buffer((j-1)*attributes+7) == type_b .and. &
              ids_buffer((j-1)*attributes+8) == atom) Then

            ! found
            packed_index = j
            exit

          End If
        End Do

        If (packed_index == -1) Then
          Call error(0, "correlator not found in revive")
        End If

        offset = offests_buffer(packed_index)
        buffer_size = sizes_buffer(packed_index)

        packed_correlators%buffer(buffer_index:(buffer_index+buffer_size-1)) = &
          data_buffer(offset:(offset+buffer_size-1))

        buffer_index = buffer_index + buffer_size
      
      End Do

    End If

    ! data packed, deport

    Call gatherv_scatterv_index_arrays(comm, &
    local_buffer_size, &
    packed_correlators%mpi%counts, &
    packed_correlators%mpi%displ & 
    )

    Call gscatterv(comm, packed_correlators%buffer, &
      packed_correlators%mpi%counts, &
      packed_correlators%mpi%displ, &
      local_buffer, root_id)

    buffer_index = 0
    If (stats%calculate_correlations) Then  
      Do i = 1, Size(cor_keys)
        If (.not. stats%cor_table%contains(cor_keys(i))) Cycle
        Call stats%cor_table%get(cor_keys(i), cor_index)
        Associate(cor_data => stats%correlations(cor_index))
          Do j = 1, Size(cor_data%correlators)
            Call cor_data%correlators(j)%recieve_buffer(local_buffer,buffer_index)
          End Do
          Call stats%cor_table%set(cor_keys(i), cor_index)
        End Associate
      End Do
    End If

    Call packed_ids%finalise()
    Call packed_correlators%finalise()

  End Subroutine revive_correlations

  !!!!!!!! observables !!!!!!!! 

    !> True if the observable name refers to a values in rigid_bodies_type
  Pure Logical Function is_rigid_observable(observable_name)
    Character(Len=MAX_CORRELATION_NAME_LENGTH), Intent(In   ) :: observable_name

    Select Case (Trim(observable_name))
    Case ("rigid_body_position")
      is_rigid_observable = .true.
    Case ("rigid_body_velocity")
      is_rigid_observable = .true.
    Case ("rigid_body_omega")
      is_rigid_observable = .true.
    Case ("rbp")
      is_rigid_observable = .true.
    Case ("rbv")
      is_rigid_observable = .true.
    Case ("rbo")
      is_rigid_observable = .true.
    Case Default
      is_rigid_observable = .false.
    End Select
  End Function is_rigid_observable

  Subroutine set_rigid_observable_name(observable_name, o)
    Character(Len=MAX_CORRELATION_NAME_LENGTH), Intent(In   ) :: observable_name
    Class(observable_rigid),                 Intent(InOut) :: o

    If (Trim(observable_name) == "rigid_body_position" .or. &
        Trim(observable_name) == "rbp") Then
      o%rigid_type = RIGID_POSITION
    Else If (Trim(observable_name) == "rigid_body_velocity" .or. &
             Trim(observable_name) == "rbv") Then
      o%rigid_type = RIGID_VELOCITY
    Else If (Trim(observable_name) == "rigid_body_omega" .or. &
             Trim(observable_name) == "rbo") Then
      o%rigid_type = RIGID_ORIENTATIONAL_VELOCITY
    End If

  End Subroutine set_rigid_observable_name


  Integer Function component_symbol_to_index(s)
    Character(Len=*), Intent(In   ) :: s

    Character(Len=2), Dimension(1:3) :: components_vector
    Character(Len=2), Dimension(1:9) :: components_matrix

    components_vector = (/ 'x', 'y', 'z' /)
    components_matrix = (/'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'/)

    component_symbol_to_index = FindLoc(components_vector, s, 1)
    If (component_symbol_to_index == 0) Then
      component_symbol_to_index = FindLoc(components_matrix, s, 1)
    End If

    If (component_symbol_to_index == 0) Then
      Call error(0, "no index for component symbol "//s)
    End If
  End Function component_symbol_to_index

  !> True if the observable name refers to a value in current_type
  Logical Function is_currents_observable(observable_name)
    Character(Len=MAX_CORRELATION_NAME_LENGTH), Intent(In   ) :: observable_name

    Integer                                                     :: i
    Character(Len=MAX_CORRELATION_NAME_LENGTH), Dimension(1:12) :: current_names = &
      (/"longitudinal_current", "transverse_current  ", "energy_current      ", &
        "kdensity            ", "energy_density      ", "kstress             ", &
        "lc                  ", "tc                  ", "ec                  ", &
        "kd                  ", "edc                 ", "ks                  "/)

    Do i = 1, Size(current_names)
      If (Trim(current_names(i)) == Trim(observable_name)) Then
        is_currents_observable = .true.
        return
      End If
    End Do
    is_currents_observable = .false.
  End Function is_currents_observable

  Subroutine set_currents_observable(observable_name, o, k, atom_type, component_name, component, site_name)
    Character(Len=MAX_CORRELATION_NAME_LENGTH), Intent(In   )           :: observable_name
    Class(observable_currents),                 Intent(InOut)           :: o
    Integer,                                    Intent(In   ), Optional :: k, atom_type, component
    Character(Len=2),                           Intent(In   ), Optional :: component_name
    Character(Len=8),                           Intent(In   ), Optional :: site_name

    If (Trim(observable_name) == "longitudinal_current" .or. &
        Trim(observable_name) == "lc") Then
      o%current_type = L_MOM_CURRENT
    Else If (Trim(observable_name) == "transverse_current" .or. &
             Trim(observable_name) == "tc") Then
      o%current_type = T_MOM_CURRENT
    Else If (Trim(observable_name) == "energy_current" .or. &
             Trim(observable_name) == "ec") Then
      o%current_type = ENG_CURRENT
    Else If (Trim(observable_name) == "kdensity" .or. &
             Trim(observable_name) == "kd") Then
      o%current_type = K_DENSITY
    Else If (Trim(observable_name) == "energy_density" .or. &
             Trim(observable_name) == "edc") Then
      o %current_type = ENG_DENSITY
    Else If (Trim(observable_name) == "kstress" .or. &
             Trim(observable_name) == "ks") Then
      o%current_type = K_STRESS
    End If

    If (Present(k)) o%kpoint = k
    If (Present(atom_type)) o%atom_type = atom_type
    If (Present(component_name)) o%component_name = component_name
    If (Present(component)) o%component = component
    If (Present(site_name)) o%atom_type_name = site_name

  End Subroutine set_currents_observable

  Subroutine character_to_observable(c, o)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 for parsing correlation observable options.
    !   format is expected to be NAME_COMPONENT for observables
    !   with components or simple NAME for those without.
    !
    ! author    - h.l.devereux 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*),                Intent(In   ) :: c
    Class(observable), Allocatable,  Intent(  Out) :: o
    
    Logical                                    :: success
    Integer                                    :: i
    Integer                                    :: component
    Character(Len=2)                           :: component_sym
    Character(Len=STR_LEN)                     :: msg, component_name
    Character(Len=MAX_CORRELATION_NAME_LENGTH) :: observable_name
    Type(observable_currents)                  :: oc
    Type(observable_rigid)                     :: or

    ! Check if in stpval
    Do i = 1, Size(stpval_names)
      If (Index(stpval_names(i), Trim(c)) > 0) Then
        Allocate(observable_statis::o)
        o%component = i
        o%component_name = stpval_names(i)
        Return
      End If
    End Do

    i = Scan(c, '_', .true.)

    If (in_range(i, (/1, Len(c)-1/))) Then
      component_name = c(i+1:Len(c))
      observable_name = c(1:i-1)
      component_sym = component_name(1:Min(Len(component_name),2))
      component = component_symbol_to_index(component_sym)
    Else If (c == "kd" .or. c == "kdensity") Then
      Call set_currents_observable(c, oc)
      Allocate(observable_currents::o)
      oc%current_type = K_DENSITY
      o = oc
      Return
    Else
      Write (msg, ('(a)')) "correlation without component, please specify a component with _, got: "//Trim(c)
      Call error(0, msg)
    End If

    success = .false.
    If (observable_name == velocity_name(observable_velocity(), .false.) .or. observable_name == "v") Then 
      If (Len(Trim(component_sym)) /= 1) Then
        Write (msg, ('(a)')) "velocity requires components x, y, or z. Got: "//Trim(c)
        Call error(0, msg)
      End If
      Allocate(observable_velocity::o)
      o%component = component
      o%component_name = component_sym
      success = .true.
    Else If (observable_name == stress_name(observable_stress(), .false.) .or. observable_name == "s") Then
      If (Len(Trim(component_sym)) /= 2) Then
        Write (msg, ('(a)')) "stress requires components xx, xy, xz, yx, yy, yz, zx, zy, or, zz. Got: "//Trim(c)
        Call error(0, msg)
      End If
      Allocate(observable_stress::o)
      o%component = component
      o%component_name = component_sym
      success = .true.
    Else If (observable_name == heat_flux_name(observable_heat_flux(), .false.) .or. observable_name == "hf") Then
      If (Len(Trim(component_sym)) /= 1) Then
        Write (msg, ('(a)')) "heat_flux requires components x, y, or z. Got: "//Trim(c)
        Call error(0, msg)
      End If
      Allocate(observable_heat_flux::o)
      o%component = component
      o%component_name = component_sym
      success = .true.
    Else If (is_currents_observable(observable_name)) Then
      Call set_currents_observable(observable_name, oc)
      Allocate(observable_currents::o)
      o = oc
      o%component = component
      o%component_name = component_sym
      success = .true.
    Else If (is_rigid_observable(observable_name)) Then
      Call set_rigid_observable_name(observable_name, or)
      Allocate(observable_rigid::o)
      o = or
      o%component = component
      o%component_name = component_sym
      success = .true.
    End If


    If (success .eqv. .false.) Then
      Call error(0,"correlation observable could not be allocated from character: "//c)
    End If

  End Subroutine character_to_observable

  Subroutine id_component_to_observable(id, component, o, subtype)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 for parsing an id and component into an observable.
    !  Used for correlation reciept between processors.
    !
    ! author    - h.l.devereux 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer,                        Intent(In   )  :: id, component, subtype
    Class(observable), Allocatable, Intent(  Out)  :: o

    Logical                          :: success
    Character(Len=100)               :: msg
    Character(Len=2), Dimension(1:3) :: components_vector
    Character(Len=2), Dimension(1:9) :: components_matrix
    Type(observable_rigid)           :: or

    components_vector = (/ 'x', 'y', 'z' /)
    components_matrix = (/'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'/)
    
    success = .false.
    If (id == velocity_id()) Then
      Allocate(observable_velocity::o)
      o%component_name = components_vector(component)
      success = .true.
    Else If (id == stress_id()) Then
      Allocate(observable_stress::o)
      o%component_name = components_matrix(component)
      success = .true.
    Else If (id == heat_flux_id()) Then
      Allocate(observable_heat_flux::o)
      o%component_name = components_vector(component)
      success = .true.
    Else If (id == current_id()) Then
      Allocate(observable_currents::o)
      o%component_name = components_vector(component)
      success = .true.
    Else If (id == rigid_id()) Then
      Allocate(observable_rigid::o)
      or%component_name = components_vector(component)
      or%rigid_type = subtype
      o = or
      success = .true.
    End If

    If (success .eqv. .false.) Then
      Write(msg,'(a,i0)') "correlation observable could not be allocated from internal code: ", id
      Call error(0,msg)
    End If

    o%component = component
  End Subroutine id_component_to_observable

  !> Observable has no subtype.
  Integer(Kind=wi) Function no_subtype(t)
    Class(observable), Intent(In   ) :: t
    no_subtype = 0_wi
  End Function no_subtype

  !> Observable not distributed.
  Integer(Kind=wi) Function not_distributed()
    not_distributed = NOT_DISTRIBUTED_OBSERVABLE
  End Function not_distributed

  !> Observable is distributed per atom.
  Integer(Kind=wi) Function distributed_per_atom()
    distributed_per_atom = PER_ATOM_OBSERVABLE
  End Function distributed_per_atom

  !> Observable is distributed per rigid body.
  Integer(Kind=wi) Function distributed_per_rigid()
    distributed_per_rigid = PER_RIGID_OBSERVABLE
  End Function distributed_per_rigid
  
  !!!!!!!!!! observable_velocity !!!!!!!!!!

  Function velocity_value(t, config, rigid, stats, index) Result(v)
    Class(observable_velocity),            Intent(In   ) :: t
    Type(configuration_type),              Intent(InOut) :: config
    Type(rigid_bodies_type),               Intent(InOut) :: rigid
    Type(stats_type),                      Intent(InOut) :: stats
    Integer,                    Optional,  Intent(In   ) :: index

    Complex(Kind=wp)       :: v
    Real(Kind=wp)          :: velocity(1:3)

    If (.not. Present(index)) Then
      Call error(0,message="no atom index specified for velocity correlator")
    End If

    velocity = (/config%vxx(index), config%vyy(index), config%vzz(index)/)
    v = Cmplx(velocity(t%component), Kind=wp)
    
  End Function velocity_value

  Function velocity_name(t, with_component) Result(v)
      Class(observable_velocity), Intent(In   )           :: t
      Logical,                    Intent(In   ), Optional :: with_component
      Character(Len=MAX_CORRELATION_NAME_LENGTH)  :: v
      v = 'velocity_'//t%component_name
      If (Present(with_component)) Then
        If (.not. with_component) Then
          v = 'velocity'
        End If
      End If
  End Function velocity_name

  Function velocity_id() Result(v)
    Integer :: v
    v = 0
  End Function velocity_id

  !!!!!!!!!! observable stress !!!!!!!!!!

  Function stress_value(t, config, rigid, stats, index) Result(v)
    Class(observable_stress),             Intent(In   ) :: t
    Type(configuration_type),             Intent(InOut) :: config
    Type(rigid_bodies_type),              Intent(InOut) :: rigid
    Type(stats_type),                     Intent(InOut) :: stats
    Integer,                    Optional, Intent(In   ) :: index
    Complex(Kind=wp)       :: v

    v = Cmplx(stats%strtot(t%component) / stats%stpvol, Kind=wp)

  End Function stress_value

  Function stress_name(t, with_component) Result(v)
      Class(observable_stress), Intent(In   )           :: t
      Logical,                  Intent(In   ), Optional :: with_component
      Character(Len=MAX_CORRELATION_NAME_LENGTH)  :: v
      v = 'stress_'//t%component_name
      If (Present(with_component)) Then
        If (.not. with_component) Then
          v = 'stress'
        End If
      End If
  End Function stress_name

  Function stress_id() Result(v)
    Integer :: v
    v = 1
  End Function stress_id

  !!!!!!!!!! observable heat_flux !!!!!!!!!!

  Function heat_flux_value(t, config, rigid, stats, index) Result(v)
    Class(observable_heat_flux),          Intent(In   ) :: t
    Type(configuration_type),             Intent(InOut) :: config
    Type(rigid_bodies_type),              Intent(InOut) :: rigid
    Type(stats_type),                     Intent(InOut) :: stats
    Integer,                    Optional, Intent(In   ) :: index
    
    Complex(Kind=wp)       :: v

    v = Cmplx(stats%heat_flux(t%component), Kind=wp)

  End Function heat_flux_value

  Function heat_flux_name(t, with_component) Result(v)
      Class(observable_heat_flux), Intent(In   )           :: t
      Logical,                     Intent(In   ), Optional :: with_component
      Character(Len=MAX_CORRELATION_NAME_LENGTH)     :: v
      v = 'heat_flux_'//t%component_name
      If (Present(with_component)) Then
        If (.not. with_component) Then
          v = 'heat_flux'
        End If
      End If
  End Function heat_flux_name

  Function heat_flux_id() Result(v)
    Integer :: v
    v = 2
  End Function heat_flux_id

  !!!!!!!!!! observable statis !!!!!!!!!!

  Function statis_value(t, config, rigid, stats, index) Result(v)
    Class(observable_statis),             Intent(In   ) :: t
    Type(configuration_type),             Intent(InOut) :: config
    Type(rigid_bodies_type),              Intent(InOut) :: rigid
    Type(stats_type),                     Intent(InOut) :: stats
    Integer,                    Optional, Intent(In   ) :: index
    
    Complex(Kind=wp)          :: v
    Character(Len=STR_LEN)    :: msg

    If (t%component < 1 .or. t%component > Size(stpval_names)) Then
      Write(msg, '(a, i0)') "Correlating non-existant/unsupported statis component ", t%component
      Call error(0, msg)
    End If
    v = Cmplx(stats%stpval(t%component), 0.0_wp, kind=wp)

  End Function statis_value

  Function statis_name(t, with_component) Result(v)
      Class(observable_statis), Intent(In   )           :: t
      Logical,                  Intent(In   ), Optional :: with_component

      Character(Len=MAX_CORRELATION_NAME_LENGTH) :: v
      v = Trim(stpval_names(t%component))
  End Function statis_name

  Function statis_id() Result(v)
    Integer :: v
    v = 3
  End Function statis_id

  !!!!!!!!!! observable currents !!!!!!!!!!

  Function current_value(t, config, rigid, stats, index) Result(v)
    Class(observable_currents),                         Intent(In   ) :: t
    Type(configuration_type),                           Intent(InOut) :: config
    Type(rigid_bodies_type),                            Intent(InOut) :: rigid
    Type(stats_type),                                   Intent(InOut) :: stats
    Integer,                    Optional,               Intent(In   ) :: index
    
    Complex(Kind=wp)       :: v

    If (.not. in_range(t%kpoint, (/1, stats%cur%nkpoints/))) Then
      Call error(0, "invalid kpoint in observable currents")
    End If

    If (.not. in_range(t%atom_type, (/1, Size(stats%cur%density_jlk, 3)/))) Then
      Call error(0, "invalid atom type in observable currents")
    End If

    Select Case (t%current_type)
      Case (L_MOM_CURRENT)
       v = stats%cur%longitudinal_jlk(t%kpoint, t%component, t%atom_type)
      Case (T_MOM_CURRENT)
        v = stats%cur%transverse_jlk(t%kpoint, t%component, t%atom_type)
      Case (ENG_CURRENT)
        v = stats%cur%energy_jlk(t%kpoint, t%component, t%atom_type)
      Case (K_DENSITY)
        v = stats%cur%density_jlk(t%kpoint, 1, t%atom_type)
      Case (ENG_DENSITY)
        v = stats%cur%energy_density_jlk(t%kpoint, t%component, t%atom_type)
      Case (K_STRESS)
        v = stats%cur%stress_jlk(t%kpoint, voigt_flat_3x3(t%component), t%atom_type)
      Case Default
        Call error(0, "invalid type for observable currents")
    End Select

  End Function current_value

  Function current_name(t, with_component) Result(v)
      Class(observable_currents), Intent(In   )           :: t
      Logical,                    Intent(In   ), Optional :: with_component

      Character(Len=MAX_CORRELATION_NAME_LENGTH) :: v
      Character(Len=4)                           :: kpoint
      Logical                                    :: comp

      If (Present(with_component)) Then
        comp = with_component
      Else
        comp = .true.
      End If

      Write (kpoint, '(i0)') t%kpoint
      Select Case (t%current_type)
        Case (L_MOM_CURRENT)
          If (comp) Then
            v = Trim(t%atom_type_name)//"-lc_"//Trim(kpoint)//"_"//Trim(t%component_name)
          Else
            v = "lc"
          End If
        Case (T_MOM_CURRENT)
          If (comp) Then
            v =  Trim(t%atom_type_name)//"-tc_"//Trim(kpoint)//"_"//Trim(t%component_name)
          Else
            v = "tc"
          End If
        Case (ENG_CURRENT)
          If (comp) Then
            v =  Trim(t%atom_type_name)//"-ec_"//Trim(kpoint)//"_"//Trim(t%component_name)
          Else
            v = "ec"
          End If
        Case (K_DENSITY)
          If (comp) Then
            v =  Trim(t%atom_type_name)//"-kd_"//Trim(kpoint)
          Else
            v = "kd"
          End If
        Case (ENG_DENSITY)
          If (comp) Then
            v =  Trim(t%atom_type_name)//"-edc_"//Trim(kpoint)//"_"//Trim(t%component_name)
          Else
            v = "edc"
          End If
        Case (K_STRESS)
          If (comp) Then
            v =  Trim(t%atom_type_name)//"-ks_"//Trim(kpoint)//"_"//Trim(t%component_name)
          Else
            v = "ks"
          End If
      End Select
  End Function current_name

  Function current_id() Result(v)
    Integer :: v
    v = 4
  End Function current_id

  !!!!!!!!!! observable_rigid !!!!!!!!!!

  Function rigid_value(t, config, rigid, stats, index) Result(v)
    Class(observable_rigid),               Intent(In   ) :: t
    Type(configuration_type),              Intent(InOut) :: config
    Type(rigid_bodies_type),               Intent(InOut) :: rigid
    Type(stats_type),                      Intent(InOut) :: stats
    Integer,                    Optional,  Intent(In   ) :: index

    Complex(Kind=wp)       :: v
    Character(Len=STR_LEN) :: msg
    Real(Kind=wp)          :: values(1:3)

    If (.not. Present(index)) Then
      Call error(0,message="no rigid body index specified for rigid body velocity correlator")
    End If

    If (t%rigid_type == RIGID_POSITION) Then
      values = (/rigid%xxx(index), rigid%yyy(index), rigid%zzz(index)/)
    Else If (t%rigid_type == RIGID_VELOCITY) Then
      values = (/rigid%vxx(index), rigid%vyy(index), rigid%vzz(index)/)
    Else If (t%rigid_type == RIGID_ORIENTATIONAL_VELOCITY) Then
      values = (/rigid%oxx(index), rigid%oyy(index), rigid%ozz(index)/)
    End If
    v = Cmplx(values(t%component), Kind=wp)

  End Function rigid_value

  Function rigid_name(t, with_component) Result(v)
      Class(observable_rigid), Intent(In   )   :: t
      Logical,                  Intent(In   ), Optional :: with_component
      Character(Len=MAX_CORRELATION_NAME_LENGTH)  :: v
      If (t%rigid_type == RIGID_POSITION) Then
        v = 'rb_position_'//t%component_name
      Else If (t%rigid_type == RIGID_VELOCITY) Then
        v = 'rb_velocity_'//t%component_name
      Else If (t%rigid_type == RIGID_ORIENTATIONAL_VELOCITY) Then
        v = 'rb_omega_'//t%component_name
      End If

  End Function rigid_name

  Function rigid_id() Result(v)
    Integer :: v
    v = 5
  End Function rigid_id

  Integer(Kind=wi) Function rigid_subtype(t)
    Class(observable_rigid), Intent(In   ) :: t
    rigid_subtype = t%rigid_type
  End Function

  Subroutine update_statistic(stat, v, step)
    Class(statistic_accumulator), Intent(InOut) :: stat
    Real(Kind=wp)                               :: v
    Integer                                     :: step

    stat%mu_old = stat%mu
    stat%stack_pos = Mod((step - 1),stat%window) + 1

    If (stat%initialised == stat%window) Then 
      stat%mu = (stat%window * stat%mu_old - stat%stack(stat%stack_pos) + v) / stat%window
      stat%ss = stat%ss - stat%stack(stat%stack_pos)**2 + v*v
      ! also multiply by n / (n-1) for unbiased, mu already divided by n
      stat%var = stat%ss / (stat%window-1) - (stat%mu**2)  * stat%window / (stat%window - 1)
    Else 
      stat%mu = (stat%initialised * stat%mu + v) / (stat%initialised + 1)
      stat%var_tmp = stat%var_tmp + (v - stat%mu_old) * (v -stat%mu)
      ! account for bias, / n-1
      stat%var = stat%var_tmp / Max(1,stat%initialised)
      stat%ss = stat%ss + v*v
      stat%initialised = stat%initialised + 1
    End If

    stat%stack(stat%stack_pos) = v

  End Subroutine update_statistic

End Module statistics