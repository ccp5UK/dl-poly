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
                             Spread_tag, comm_self, comms_type, gcheck, girecv, gmax, gsend, gsum, &
                             gsync, gtime, gwait, mode_create, mode_wronly, offset_kind, &
                             gatherv_scatterv_index_arrays, ggatherv, gscatterv, root_id, &
                             gscatter
                             
  Use configuration,   Only: configuration_type
  Use constants,       Only: boltz,&
                             engunit,&
                             eu_ev,&
                             eu_kcpm,&
                             eu_kjpm,&
                             pi,&
                             prsunt,&
                             tenunt,&
                             zero_plus
  Use currents,        Only: current_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc,&
                             info,&
                             warning
  Use filename,        Only: FILE_STATS, FILE_COR, &
                             file_type
  Use flow_control,    Only: RESTART_KEY_OLD
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
                             shellsort2
  Use site,            Only: site_type
  Use thermostat,      Only: CONSTRAINT_NONE,&
                             CONSTRAINT_SEMI_ORTHORHOMBIC,&
                             CONSTRAINT_SURFACE_TENSION,&
                             thermostat_type
  Use timer,           Only: timer_type, stop_timer, start_timer
  Use z_density,       Only: z_density_collect,&
                             z_density_type
  Use correlators,     Only: correlator, correlator_buffer_type, indices_buffer_type
  Use units,           Only: to_out_units 

  Implicit None

  Private

  Integer, Parameter :: MAX_CORRELATION_NAME_LENGTH = 16

  ! correlation observables, and interface
  Type, Abstract, Public :: observable
    Contains
      Procedure(get_value),     Deferred :: value
      Procedure(get_dimension), Deferred :: dimension
      Procedure(get_name),      Deferred :: name
      Procedure(get_id)  ,      Deferred :: id
      Procedure(is_per_atom),   Deferred :: per_atom
  End Type observable

  Type, Public :: observable_holder
    Class(observable), Allocatable :: observable
  End Type

  Type, Public :: correlation
    Class(observable), Allocatable :: A
    Class(observable), Allocatable :: B
    ! 0 indicates not tracking an atom but a global property
    ! >= 1 indicates a local atom's index, or global resp.
    Integer                        :: atom, atom_global
  End Type 

  Type, Public :: correlator_holder
    Type(correlator)               :: correlator
    Type(correlation)              :: correlation
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
                                         strpmf(1:9) = 0.0_wp, stress(1:9) = 0.0_wp, strdpd(1:9) = 0.0_wp
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
    Integer                            :: number_of_correlations = 0, max_buffer_per_atom = 0
    Type(correlator_holder), Allocatable :: correlations(:)
    Type(correlation),       Allocatable :: unique_correlations(:)
    Integer,                 Allocatable :: unique_correlation_params(:)
    Integer                              :: cor_dump_freq = 0 
    ! Integer, Allocatable               :: number_of_blocks(:), points_per_block(:), window_size(:)
    ! Integer, Allocatable               :: min_distance(:), dim_left(:), dim_right(:)
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

    !> store spot heat flux
    Real(Kind=wp)                      :: heat_flux(1:3) = [0.0_wp, 0.0_wp, 0.0_wp]

    !> Store for per-particle energy data
    Real(Kind=wp), Allocatable         :: pp_energy(:)

    !> Store for per-particle stress data
    Real(Kind=wp), Allocatable         :: pp_stress(:, :)

    !> Whether per-particle information is needed
    Logical :: require_pp = .false.

    !> Whether this step is a step to collect per-particle data
    Logical :: collect_pp = .false.

    !> Whether heat_flux is correlated
    Logical :: correlating_heat_flux = .false.

  Contains
    Private

    Procedure, Public :: init              => allocate_statistics_arrays
    Procedure, Public :: init_connect      => allocate_statistics_connect
    Procedure, Public :: init_correlations => allocate_correlations_arrays
    Procedure, Public :: init_correlator   => allocate_correlator
    Procedure, Public :: clean_connect     => deallocate_statistics_connect
    Procedure, Public :: update_stress
    Procedure, Public, Pass :: allocate_per_particle_arrays
    Procedure, Public, Pass :: deallocate_per_particle_arrays
    Procedure, Public :: correlator_deport
    Procedure, Public :: correlator_recieve
    Procedure, Public :: dump_correlations 
    Procedure, Public :: revive_correlations 
    Procedure, Public :: reindex_correlators
    Final :: cleanup
  End Type

  Abstract Interface 
    ! Kernal for selecting data for correlation
    Subroutine get_value(t, config, stats, v, atom)
        Import observable, configuration_type, stats_type, wp
        Class(observable),          Intent(In   ) :: t
        Type(configuration_type),   Intent(InOut) :: config
        Type(stats_type),           Intent(InOut) :: stats
        Real(Kind=wp), Allocatable, Intent(InOut) :: v(:)
        Integer,       Optional,    Intent(In   ) :: atom
    End Subroutine get_value

    ! utility to get size of observable
    Function get_dimension(t) Result(v)
        Import observable
        Class(observable), Intent(In   ) :: t
        Integer                          :: v
    End Function get_dimension

    ! utility to get name of observable (i.e. for i/o)
    Function get_name(t) Result(v)
        Import observable, MAX_CORRELATION_NAME_LENGTH
        Class(observable), Intent(In   ) :: t
        Character(Len=MAX_CORRELATION_NAME_LENGTH)                 :: v
    End Function get_name

    ! utility to get numerical id of observable (i.e. for revive)
    Function get_id(t) Result(v)
      Import observable
      Class(observable), Intent(In   ) :: t
      Integer                          :: v
    End Function get_id

    Function is_per_atom(t) Result(v)
      Import observable
      Class(observable), Intent(In   ) :: t
      Logical                          :: v
    End Function is_per_atom
  End Interface

  Interface operator (==)
    Module Procedure is_equal
  End Interface

  Type, Extends(observable), Public :: observable_velocity
  Contains
      Procedure :: value     => velocity_value
      Procedure :: dimension => velocity_dimension
      Procedure :: name      => velocity_name
      Procedure :: id        => velocity_id
      Procedure :: per_atom  => velocity_per_atom
  End Type

  Type, Extends(observable), Public :: observable_stress
  Contains
      Procedure :: value     => stress_value
      Procedure :: dimension => stress_dimension
      Procedure :: name      => stress_name
      Procedure :: id        => stress_id
      Procedure :: per_atom  => stress_per_atom
  End Type
  Type, Extends(observable), Public :: observable_heat_flux
  Contains
      Procedure :: value     => heat_flux_value
      Procedure :: dimension => heat_flux_dimension
      Procedure :: name      => heat_flux_name
      Procedure :: id        => heat_flux_id
      Procedure :: per_atom  => heat_flux_per_atom
  End Type
  
  Public :: calculate_stress
  Public :: calculate_viscosity
  Public :: calculate_heat_flux
  Public :: calculate_thermal_conductivity
  Public :: statistics_collect
  Public :: statistics_connect_frames
  Public :: statistics_connect_set
  Public :: write_per_part_contribs
  Public :: write_header
  Public :: statistics_result
  Public :: correlation_result
  Public :: character_to_observable
  Public :: code_to_observable
Contains

  Subroutine allocate_statistics_arrays(stats, mxrgd, mxatms, mxatdm)
    Class(stats_type), Intent(InOut)   :: stats
    Integer,           Intent(In   )   :: mxrgd, mxatms, mxatdm

    Integer                            :: mxnstk, mxstak, nxatms
    Integer,           Dimension(1:4)  :: fail
 
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
    If (Any(fail > 0)) Call error_alloc("allocate_statistics_arrays", "statistics")

    stats%xin = 0.0_wp; stats%yin = 0.0_wp; stats%zin = 0.0_wp
    stats%xto = 0.0_wp; stats%yto = 0.0_wp; stats%zto = 0.0_wp; stats%rsd = 0.0_wp

    stats%stpval = 0.0_wp; stats%stpvl0 = 0.0_wp; stats%sumval = 0.0_wp; stats%ssqval = 0.0_wp
    stats%zumval = 0.0_wp; stats%ravval = 0.0_wp; stats%stkval = 0.0_wp

  End Subroutine allocate_statistics_arrays

  Subroutine allocate_per_particle_arrays(stats, natms)
    Class(stats_type), Intent(InOut) :: stats
    Integer,           Intent(In   ) :: natms

    Integer :: fail

    If (.not. Allocated(stats%pp_energy)) Then
      Allocate (stats%pp_energy(natms), stat=fail)
      If (fail > 0) Call error_alloc("stats%pp_energy", "statistics")
    End If

    If (.not. Allocated(stats%pp_stress)) Then
      Allocate (stats%pp_stress(9, natms), stat=fail)
      If (fail > 0) Call error_alloc("stats%pp_stress", "statistics")
    End If

    stats%pp_energy = 0.0_wp
    stats%pp_stress = 0.0_wp
    stats%collect_pp = .true.

  End Subroutine allocate_per_particle_arrays

  Subroutine deallocate_per_particle_arrays(stats)
    Class(stats_type), Intent(InOut) :: stats

    Integer :: fail

    Deallocate (stats%pp_energy, stat=fail)
    If (fail > 0) Call error_dealloc("stats%pp_energy", "statistics")
    Deallocate (stats%pp_stress, stat=fail)
    If (fail > 0) Call error_dealloc("stats%pp_stress", "statistics")

    stats%collect_pp = .false.

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

    !Call deallocate_correlations_arrays()
  End Subroutine cleanup

  Subroutine allocate_correlations_arrays(stats)
    Class(stats_type), Intent(InOut) :: stats
    Integer,           Dimension(1)  :: fail

    Allocate(stats%correlations(1:stats%number_of_correlations), Stat = fail(1))

    If (Any(fail > 0)) Call error_alloc("allocate_correlations_arrays", "statistics")

  End Subroutine allocate_correlations_arrays

  Subroutine allocate_correlator(stats, atom, global, blocks, points, window, &
    A, B, correlator_index)
    Class(stats_type),        Intent(InOut) :: stats
    Integer,                  Intent(In   ) :: atom, global, blocks, points, window, &
                                               correlator_index
    Class(observable),        Intent(In   ) :: A, B
    Integer                                 :: dim_left, dim_right
    
    stats%correlations(correlator_index)%correlation%A = A
    stats%correlations(correlator_index)%correlation%B = B
    stats%correlations(correlator_index)%correlation%atom = atom
    stats%correlations(correlator_index)%correlation%atom_global = global


    dim_left = stats%correlations(correlator_index)%correlation%A%dimension()
    dim_right = stats%correlations(correlator_index)%correlation%B%dimension()

    Call stats%correlations(correlator_index)%correlator%init(blocks, points, window, dim_left, dim_right)
  
  End Subroutine allocate_correlator

  Subroutine correlation_result(stats, comm, files, config, sites, nstep, time)
    
    Class(stats_type),        Intent(InOut)       :: stats
    Type(comms_type),         Intent(InOut)       :: comm
    Type(file_type),          Intent(InOut)       :: files(:)
    Type(configuration_type), Intent(In   )       :: config
    Type(site_type),          Intent(In   )       :: sites
    Integer,                  Intent(In   )       :: nstep
    Real(Kind=wp),            Intent(In   )       :: time
    Integer                                       :: i, tau, j, k, flat_dim, l, r, &
                                                     file_unit, atom
    Real(Kind=wp), Allocatable                    :: cor_accumulator(:,:,:,:), correlation(:,:,:)
    Real(Kind=wp), Allocatable                    :: flat_correlation(:)
    Real(Kind=wp), Allocatable                    :: timesteps(:)    
    Integer,       Allocatable                    :: type_counts(:)
    Character(Len=MAX_CORRELATION_NAME_LENGTH*2+1)  :: correlation_name
    Character(Len=MAX_CORRELATION_NAME_LENGTH)    :: char_left, char_right, &
                                                     component_left, component_right
    Character(Len=2), Dimension(1:3)              :: components_vector
    Character(Len=2), Dimension(1:9)              :: components_matrix
    Real(Kind=wp)                                 :: t, dt, visc, therm_cond, conv
    Integer                                       :: points, window, blocks,&
                                                     dim_left, dim_right
    Character(Len=STR_LEN)                        :: units
                                              
                        
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
      Write (file_unit, '(a)') "correlations:"
    
    End If

    Do i = 1, Size(stats%unique_correlations)

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

      char_left = stats%unique_correlations(i)%A%name()
      char_right = stats%unique_correlations(i)%B%name()
      correlation_name = Trim(char_left)//'-'//Trim(char_right)

      If (stats%unique_correlations(i)%atom > 0) Then

        Do j = 1, stats%number_of_correlations

          If (stats%correlations(j)%correlation%A == stats%unique_correlations(i)%A .and. &
              stats%correlations(j)%correlation%B == stats%unique_correlations(i)%B ) Then

            ! found a match, assume all same (they should be)
            
            points = stats%correlations(j)%correlator%points_per_block
            blocks = stats%correlations(j)%correlator%number_of_blocks
            dim_left = stats%correlations(j)%correlator%left_dim
            dim_right = stats%correlations(j)%correlator%right_dim
            window = stats%correlations(j)%correlator%window_size

            Allocate(cor_accumulator(1:sites%mxatyp, 1:points*blocks, &
              1:dim_left, 1:dim_right ))

            Allocate(correlation(1:points*blocks,1:dim_left,1:dim_right))

            flat_dim = sites%mxatyp*points*blocks*dim_left*dim_right

            Allocate(flat_correlation(1:flat_dim))
            Allocate(type_counts(1:sites%mxatyp))

            type_counts = 0
            flat_correlation = 0.0
            correlation = 0.0
            cor_accumulator = 0.0
            
            ! data allocated, now exit to accumulate 

            exit
          End If

        End Do

        ! accumulate average

        Do j = 1, stats%number_of_correlations

          If (stats%correlations(j)%correlator%count_updated == 0) Then
            ! no data was seen in this correlator, distinct from 0 
            ! correlation case
            Cycle
          End If
          
          If (stats%correlations(j)%correlation%A == stats%unique_correlations(i)%A .and. &
              stats%correlations(j)%correlation%B == stats%unique_correlations(i)%B ) Then

            atom = stats%correlations(j)%correlation%atom

            correlation = 0.0
            Call stats%correlations(j)%correlator%get_correlation(correlation)

            cor_accumulator(config%ltype(atom),:,:,:) = &
              cor_accumulator(config%ltype(atom),:,:,:) + correlation

            type_counts(config%ltype(atom)) = type_counts(config%ltype(atom)) + 1

          End If

        End Do

        ! now collect onto root

        flat_correlation = Reshape(cor_accumulator,(/flat_dim/))
        Call gsum(comm, flat_correlation)
        flat_correlation = flat_correlation 
        cor_accumulator = Reshape(flat_correlation,(/sites%mxatyp,points*blocks,&
         dim_left,dim_right/))

        ! for later averaging
        Call gsum(comm,type_counts)
        
        If (comm%idnode == root_id) Then

          Do j = 1,sites%mxatyp

            Write (file_unit, '(*(a))') "    - name: [", correlation_name, ", ", &
                                                        Trim(sites%unique_atom(j)), "]"

            Write (file_unit, '(a)') "      parameters:"
            Write (file_unit, '(a,i0)') "            points_per_block: ", &
                    points
            Write (file_unit, '(a,i0)') "            number_of_blocks: ", &
                    blocks
            Write (file_unit, '(a,i0)') "            window_size: ", &
                    window

            t = 0.0
            dt = time / nstep

            If (Allocated(timesteps)) Then
              Deallocate(timesteps)
            End If

            Allocate(timesteps(1:points*blocks))

            Do tau = 1,points*blocks
              timesteps(tau) = t
              t = t + dt
            End Do  
            
            
            Write(file_unit, '(a,*(g16.8,","))',advance="no") "      lags: [", timesteps(1:Size(timesteps)-1)
            
            Write(file_unit, '(g16.8,a)') timesteps(Size(timesteps)), "]"

            cor_accumulator(j,:,:,:) = cor_accumulator(j,:,:,:) / (1+type_counts(j))

            Write (file_unit, '(a)') "      components: "

            Do l = 1,dim_left
              
              If (dim_left == 3) Then
                component_left = Trim(components_vector(l))
              Else If (dim_left == 9) Then 
                component_left = Trim(components_matrix(l))
              End If

              Do r = 1,dim_right

                If (dim_right == 3) Then
                  component_right = Trim(components_vector(r))
                Else If (dim_right == 9) Then 
                  component_right = Trim(components_matrix(r))
                End If

                Write (file_unit, '(a,a,a)',advance="no") "           ", &
                  Trim(char_left)//'_'//Trim(component_left)//"-"//Trim(char_right)//'_'//Trim(component_right), ": "

                  Write(file_unit, '(a,*(g16.8,","))',advance="no") "[", cor_accumulator(j,1:Size(cor_accumulator,2)-1,l,r)
                  
                  Write(file_unit, '(g16.8,a)') cor_accumulator(j,Size(cor_accumulator,2),l,r), "]"

              End Do
            End Do  
          End Do
        End If

      Else

        If (comm%idnode == root_id) Then

          k = 0
          Do j = 1, stats%number_of_correlations

            If (stats%correlations(j)%correlation%A == stats%unique_correlations(i)%A .and. &
                stats%correlations(j)%correlation%B == stats%unique_correlations(i)%B ) Then

              ! found a match (should be unique)
              k = j
              
              points = stats%correlations(j)%correlator%points_per_block
              blocks = stats%correlations(j)%correlator%number_of_blocks
              dim_left = stats%correlations(j)%correlator%left_dim
              dim_right = stats%correlations(j)%correlator%right_dim
              window = stats%correlations(j)%correlator%window_size

              Allocate(correlation(1:points*blocks,1:dim_left,1:dim_right))

              flat_dim = points*blocks*dim_left*dim_right

              Allocate(flat_correlation(1:flat_dim))

              flat_correlation = 0.0
              correlation = 0.0
            
              exit
            End If
            
          End Do

          If (k == 0) Then
            Call error(0,"correlation not found for correlation result")
          End If

          Call stats%correlations(k)%correlator%get_correlation(correlation)

          Write (file_unit, '(*(a))') "    - name: [", correlation_name, ", global]"

          Write (file_unit, '(a)') "      parameters:"
          Write (file_unit, '(a,i0)') "            points_per_block: ", &
                  points
          Write (file_unit, '(a,i0)') "            number_of_blocks: ", &
                  blocks
          Write (file_unit, '(a,i0)') "            window_size: ", &
                  window

          t = 0.0
          dt = time / nstep

          If (Allocated(timesteps)) Then
            Deallocate(timesteps)
          End If

          Allocate(timesteps(1:points*blocks))

          Do tau = 1,points*blocks
            timesteps(tau) = t
            t = t + dt
          End Do  

          If (char_left == stress_name(observable_stress()) .and. &
              char_right == stress_name(observable_stress())) Then

            visc = calculate_viscosity(stats, correlation, dt)

            Write (file_unit, '(a)') "      derived:"
            Write (file_unit, '(a)') "            viscosity:"
            Write (file_unit, '(a,g16.8)') "                  value: ", visc
            Write (file_unit, '(a)') "                  units: Katm ps "

            correlation = correlation * prsunt * prsunt

          Else If(char_left == heat_flux_name(observable_heat_flux()) .and. &
                  char_right == heat_flux_name(observable_heat_flux())) Then

            therm_cond = calculate_thermal_conductivity(stats, correlation, dt, units)

            Write (file_unit, '(a)') "      derived:"
            Write (file_unit, '(a)') "            thermal-conductivity:"
            Write (file_unit, '(a,g16.8)') "                  value: ", therm_cond
            Write (file_unit, '(a,a)') "                  units: ", units
            
            Call to_out_units(1.0_wp, "internal_e", conv)

            correlation = correlation / conv

          End If
          
          Write(file_unit, '(a,*(g16.8,","))',advance="no") "      lags: [", timesteps(1:Size(timesteps)-1)
          
          Write(file_unit, '(g16.8,a)') timesteps(Size(timesteps)), "]"

          Write (file_unit, '(a)') "      components: "

          Do l = 1,dim_left
              
            If (dim_left == 3) Then
              component_left = Trim(components_vector(l))
            Else If (dim_left == 9) Then 
              component_left = Trim(components_matrix(l))
            End If

            Do r = 1,dim_right

              If (dim_right == 3) Then
                component_right = Trim(components_vector(r))
              Else If (dim_right == 9) Then 
                component_right = Trim(components_matrix(r))
              End If              

              Write (file_unit, '(a,a,a)',advance="no") "           ", &
              Trim(char_left)//'_'//Trim(component_left)//"-"//Trim(char_right)//'_'//Trim(component_right), ": "

                Write(file_unit, '(a,*(g16.8,","))',advance="no") "[", correlation(1:Size(correlation,1)-1,l,r)
                
                Write(file_unit, '(g16.8,a)') correlation(Size(correlation,1),l,r), "]"

            End Do
          End Do  
          
        End If
        
      End If

    End Do

    If (comm%idnode == root_id) Then
      Close(file_unit)
    End If

  End Subroutine correlation_result
    

  Subroutine statistics_collect(config, lsim, leql, nsteql, lmsd, keyres, degfre, degshl, &
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
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(InOut) :: config
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
    Integer                    :: fail, i, iadd, j, k, kstak
    Logical                    :: ffpass, l_tmp
    Real(Kind=wp)              :: celprp(1:10), h_z, sclnv1, sclnv2, stpcns, stpipv, stprot, &
                                  stpshl, zistk
    Real(Kind=wp), Allocatable :: amsd(:), xxt(:), yyt(:), zzt(:)
    Real(Kind=wp), Allocatable :: observable_a(:), observable_b(:)

#ifdef CHRONO
    Call start_timer(tmr, "Stats")
#endif

    ffpass = ff == 1

    fail = 0

    Allocate (amsd(1:sites%mxatyp), Stat=fail)
    If (fail > 0) Call error_alloc('amsd', 'statistics_collect')

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
          sunits = "DL_POLY Internal UNITS (10 J/mol)"
        Else If (Abs(engunit - boltz) <= zero_plus) Then
          sunits = "Kelvin/Boltzmann"
        Else ! once in a blue moon
          sunits = "DPD (Unknown)"
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

    stprot = 2.0_wp * (stats%engrot) / (boltz * Max(1.0_wp, Real(degrot, wp)))

    ! core-shell units temperature

    stpshl = 2.0_wp * (stats%shlke) / (boltz * Max(1.0_wp, Real(degshl, wp)))

    ! system temperature

    stats%stptmp = 2.0_wp * (stats%engke + stats%engrot) / (boltz * Real(degfre, wp))

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
    stats%stpval(27) = stats%stpprs * prsunt

    iadd = 27

    ! iadd = iadd + 1 ! for the stpval(0)!!! Thus to account for in printing
    ! pressure tensor (derived for the stress tensor)

    Do i = 1, 9
      stats%stpval(iadd + i) = stats%strtot(i) * prsunt / stats%stpvol
    End Do
    iadd = iadd + 9

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

        If ((.not. leql) .or. nstep >= nsteql) Then
      
          Do j = 1, stats%number_of_correlations

            i = stats%correlations(j)%correlation%atom

            If (Allocated(observable_a)) Then 
              Deallocate(observable_a)
            End If

            If (Allocated(observable_b)) Then 
              Deallocate(observable_b)
            End If

            If (i > 0) Then
                
              Call stats%correlations(j)%correlation%A%value(config, stats, observable_a, i)
              Call stats%correlations(j)%correlation%B%value(config, stats, observable_b, i)

              Call stats%correlations(j)%correlator%update(observable_a,observable_b)

            Else If (comm%idnode == root_id) Then

              Call stats%correlations(j)%correlation%A%value(config, stats, observable_a)
              Call stats%correlations(j)%correlation%B%value(config, stats, observable_b)

              Call stats%correlations(j)%correlator%update(observable_a,observable_b)
            
            End If

          End Do

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

        If (lmsd) Then
          If (stats%file_yaml) Then
            Write (fmtt, '(a,i0,a)') '(2x,a4,i0,",",', iadd + 1 - 2 * mxatdm, '(g16.8,","),g16.8,a2)'
            Write (files(FILE_STATS)%unit_no, fmt=Trim(fmtt)) "- [ ", nstep, time, &
              stats%stpval(1:27), stats%stpval(0), stats%stpval(28:36), &
              stats%stpval(37 + 2 * mxatdm:iadd), ' ]'
          Else
            Write (files(FILE_STATS)%unit_no, '(i10,1p,e14.6,0p,i10,/, (1p,5e14.6))') &
              nstep, time, iadd + 1 - 2 * mxatdm, stats%stpval(1:27), stats%stpval(0), stats%stpval(28:36), &
              stats%stpval(37 + 2 * mxatdm:iadd)
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

      If (stats%cor_dump_freq > 0) Then
        If (Mod(nstep, stats%cor_dump_freq) == 0) Then
          Call correlation_result(stats, comm, files, config, sites, nstep, time)
        End If
      End If 
    End If
    ! check on number of variables for stack

    If (iadd > stats%mxnstk) Call error(170)

    ! No totals for timestep zero

    If (nstep /= 0) Then

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
    Call stop_timer(tmr, "Stats")
#endif

  End Subroutine statistics_collect

  Subroutine statistics_connect_frames(config, megatm, mxatdm, lmsd, stats, domain, comm)

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
    Logical, Intent(In) :: lmsd
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
                  j = 36 + 2 * config%lsi(i)
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

      If (mdir /= 0) Call statistics_connect_spread(config, mdir, mxatdm, lmsd, stats, domain, comm)

      ! Sort past frame remainder of global atom indices

      !    lsi0=0 ; lsa0=0
      Do i0 = 1, stats%natms0
        stats%lsi0(i0) = i0
        stats%lsa0(i0) = stats%ltg0(i0)
      End Do
      Call shellsort2(stats%natms0, stats%lsi0, stats%lsa0)

    End Subroutine match_compress_spread_sort

  End Subroutine statistics_connect_frames

  Subroutine statistics_connect_set(config, rcut, mxatdm, lmsd, stats, domain, comm)

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
    Logical, Intent(In) :: lmsd
    Type(stats_type), Intent(InOut) :: stats
    Type(domains_type), Intent(In) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type), Intent(InOut) :: comm

    Real(Kind=wp) :: cut

    Integer           :: nlx, nly, nlz, i, i0, kk
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

      If (lmsd) Then
        i0 = 2 * stats%natms0
        stats%stpvl00(1:i0) = stats%stpvl0(37:36 + i0) !;stats%stpvl00(i0+1: )=0.0_wp
        stats%stpval0(1:i0) = stats%stpval(37:36 + i0) !;stats%stpval0(i0+1: )=0.0_wp
        stats%zumval0(1:i0) = stats%zumval(37:36 + i0) !;stats%zumval0(i0+1: )=0.0_wp
        stats%ravval0(1:i0) = stats%ravval(37:36 + i0) !;stats%ravval0(i0+1: )=0.0_wp
        stats%ssqval0(1:i0) = stats%ssqval(37:36 + i0) !;stats%ssqval0(i0+1: )=0.0_wp
        stats%sumval0(1:i0) = stats%sumval(37:36 + i0) !;stats%sumval0(i0+1: )=0.0_wp
        Do kk = 1, stats%mxstak
          stats%stkval0(kk, 1:i0) = stats%stkval(kk, 1:i0) !;stats%stkval0(kk,i0+1: )=0.0_wp
        End Do
      End If
    End If

  End Subroutine statistics_connect_set

  Subroutine statistics_connect_spread(config, mdir, mxatdm, lmsd, stats, domain, comm)

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
    Logical,                  Intent(In   ) :: lmsd
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

  Subroutine write_header()
    Character(Len=STR_LEN), Dimension(5) :: messages

    Write (messages(1), '(a)') Repeat('-', 130)
    Write (messages(2), '(9x,a4,5x,a7,1x,a11,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
      'step', 'eng_tot', 'temp_tot[K]', 'eng_cfg', 'eng_src', 'eng_cou', 'eng_bnd', 'eng_ang', 'eng_dih', 'eng_tet'
    Write (messages(3), '(5x,a8,5x,a7,1x,a11,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
      'time[ps]', ' eng_pv', 'temp_rot[K]', 'vir_cfg', 'vir_src', 'vir_cou', 'vir_bnd', 'vir_ang', 'vir_con', 'vir_tet'
    Write (messages(4), '(5x,a8,5x,a7,1x,a11,5x,a7,5x,a7,4x,a8,5x,a7,4x,a8,5x,a7,7x,a5)') &
      'cpu  [s]', 'volume', 'temp_shl[K]', 'eng_shl', 'vir_shl', 'alpha[o]', 'beta[o]', 'gamma[o]', 'vir_pmf', 'press'
    Write (messages(5), '(a)') Repeat('-', 130)
    Call info(messages, 5, .true.)

  End Subroutine write_header

  Subroutine statistics_result(config, minim, lmsd, &
                               nstrun, keyshl, megcon, megpmf, &
                               nstep, time, tmst, &
                               mxatdm, neigh_uncond_update, stats, thermo, sites, comm, files, tmr)

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
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(site_type),          Intent(In   ) :: sites
    Type(comms_type),         Intent(InOut) :: comm
    Type(file_type),          Intent(InOut) :: files(:)
    Type(timer_type),         Intent(InOut) :: tmr

    Character(Len=STR_LEN)               :: message
    Character(Len=STR_LEN), Dimension(5) :: messages
    Integer                              :: i, iadd, mxnstk
    Logical                              :: check
    Real(Kind=wp)                        :: avvol, dc, h_z, srmsd, timelp, tmp, tx, ty

    mxnstk = stats%mxnstk

    Call info('', .true.)

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
      Write (message, '(2(a,i9,a,f10.3),a)') '# run terminated after ', nstep, &
        ' steps (', time, ' ps), final averages calculated over', stats%numacc, &
        ' steps (', tmp, ' ps)'
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
        Write (message, '(a)') 'Pressure tensor  (katms):'
        Call info(message, .true.)

        Do i = iadd, iadd + 6, 3
          Write (message, '(2x,1p,3e12.4)') stats%stpval(i + 1:i + 3)
          Call info(message, .true.)
        End Do

        Write (message, '(2x,a,1p,e12.4)') 'trace/3  ', (stats%stpval(iadd + 1) + &
                                                         stats%stpval(iadd + 5) + stats%stpval(iadd + 9)) / 3.0_wp
        Call info(message, .true.)
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
      Call write_header()

      Write (messages(1), '(i13,1p,9e12.4)') stats%numacc, stats%sumval(1:9)
      Write (messages(2), '(f13.5,1p,9e12.4)') tmp, stats%sumval(10:18)
      Write (messages(3), '(0p,f13.3,1p,9e12.4)') timelp, stats%sumval(19:27)
      Write (messages(4), '(a)') ''
      Call info(messages, 4, .true.)

      Write (messages(1), '(6x,a8,1p,9e12.4)') ' r.m.s. ', stats%ssqval(1:9)
      Write (messages(2), '(6x,a8,1p,9e12.4)') 'fluctu- ', stats%ssqval(10:18)
      Write (messages(3), '(6x,a8,1p,9e12.4)') 'ations  ', stats%ssqval(19:27)
      Write (messages(4), '(a)') Repeat('-', 130)
      Call info(messages, 4, .true.)

      ! Some extra information - conserved quantity=extended ensemble energy

      Write (message, "(a,1p,e12.4,5x,a,1p,e12.4)") &
        "Extended energy:       ", stats%sumval(0), &
        " r.m.s. fluctuations:  ", stats%ssqval(0)
      Call info(message, .true.)

      ! Some extra information - <P*V> term - only matters for NP/sT ensembles

      If (thermo%variable_cell) Then
        Write (message, "(a,1p,e12.4,5x,a,1p,e12.4)") &
          "<P*V> term:            ", stats%sumval(37 + sites%ntype_atom + 2 * Merge(mxatdm, 0, lmsd)), &
          " r.m.s. fluctuations:  ", stats%ssqval(37 + sites%ntype_atom + 2 * Merge(mxatdm, 0, lmsd))
        Call info(message, .true.)
      End If

      Write (messages(1), "(130('-'))")
      Write (messages(2), '(a)') ''
      Call info(messages, 2, .true.)

      ! Move at the end of the default 27 quantities

      iadd = 27

      ! print out average pressure tensor

      If (comm%idnode == 0) Then
        Write (messages(1), '(a)') 'Pressure tensor:'
        Write (messages(2), '(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)', 'r.m.s. fluctuations'
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

      iadd = iadd + 9

      If (lmsd) iadd = iadd + 2 * mxatdm

      ! Write out estimated diffusion coefficients

      Write (messages(1), '(a)') 'Approximate 3D Diffusion Coefficients and square root of MSDs:'
      Write (messages(2), '(6x,a4,2x,a19,6x,a15)') 'atom', 'DC (10^-9 m^2 s^-1)', 'Sqrt[MSD] (Ang)'
      Call info(messages, 2, .true.)

      Do i = 1, sites%ntype_atom
        If (sites%num_type_nf(i) > zero_plus) Then
          dc = 10.0_wp * (stats%ravval(iadd + i) - stats%sumval(iadd + i)) / &
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
          Write (message, '(a32,33x,a19)') 'Average cell vectors     (Angs) ', 'r.m.s. fluctuations'
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

          Write (message, "('Average surface area, fluctuations & mean estimate (Angs^2)')")
          Call info(message, .true.)
          Write (message, '(1p,3e12.4)') stats%sumval(iadd + 2), stats%ssqval(iadd + 2), avvol / h_z
          Call info(message, .true.)

          iadd = iadd + 2

          If (Any(thermo%iso == [CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC])) Then
            tx = -h_z * (stats%sumval(29) / prsunt - (thermo%press + thermo%stress(1))) * tenunt
            ty = -h_z * (stats%sumval(30) / prsunt - (thermo%press + thermo%stress(5))) * tenunt
            Write (message, "('Average surface tension, fluctuations & mean estimate in x (dyn/cm)')")
            Call info(message, .true.)
            Write (message, '(1p,3e12.4)') stats%sumval(iadd + 1), stats%ssqval(iadd + 1), tx
            Call info(message, .true.)
            Write (message, "('Average surface tension, fluctuations & mean estimate in y (dyn/cm)')")
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
      Call correlation_result(stats,comm,files,config,sites, nstep, time)
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

  Function calculate_viscosity(stats, correlation, dt) Result(visc)
    Type(stats_type),  Intent(In   ) :: stats
    Real(Kind=wp),     Intent(In   ) :: correlation(:,:,:)
    Real(Kind=wp),     Intent(In   ) :: dt
    Real(Kind=wp)                    :: visc
    Class(integrator), Allocatable   :: inter 

    ! combine xy, yz, zx

    Allocate(simpsons_rule::inter)

    visc = inter%integrate_uniform( (correlation(:,2,2) + correlation(:,6,6) + correlation(:,7,7))/3.0_wp, dt)

    ! in internal pressure-time units

    visc = ( stats%stpvol / (boltz*stats%stptmp) ) * visc
    
    ! now in Katm - internal time 

    visc = visc * prsunt 

  End Function calculate_viscosity

  Function calculate_thermal_conductivity(stats, correlation, dt, units) Result(therm_cond)
    Type(stats_type),       Intent(In   ) :: stats
    Real(Kind=wp),          Intent(In   ) :: correlation(:,:,:)
    Real(Kind=wp),          Intent(In   ) :: dt
    Character(Len=STR_LEN), Intent(  Out) :: units
    Real(Kind=wp)                         :: therm_cond, conv
    Class(integrator),      Allocatable   :: inter 
                                        

    ! z component
    
    Allocate(simpsons_rule::inter)

    therm_cond = inter%integrate_uniform(correlation(:,3,3),dt)

    ! careful, DL_POLY already divides by volume in calculate_heat_flux

    Call to_out_units(1.0_wp, "internal_e", conv, units)

    therm_cond = stats%stpvol / (stats%stptmp * stats%stptmp * (boltz/engunit)) * therm_cond

    units = Trim(units)//" / (ps Ang K)"
    
  End Function calculate_thermal_conductivity


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

  Subroutine correlator_recieve(this, config, newatm, buffer, buffer_index)
      Class(stats_type),                       Intent(InOut)  :: this
      Type(configuration_type),                Intent(InOut)  :: config
      Integer,                                 Intent(In   )  :: newatm
      Real(Kind=wp),           Dimension(:),   Intent(InOut)  :: buffer
      Integer,                                 Intent(InOut)  :: buffer_index
      Type(correlator_holder), Allocatable                    :: tmp_cors(:)
      Class(observable),       Allocatable                    :: A, B
      Integer                                                 :: i, iA, iB, global_index, local_index, &
                                                                 window, blocks, points, dim_left, dim_right, &
                                                                 jA, jB, new_index, s

      If (this%per_atom_correlations .and. this%calculate_correlations) Then

        ! data packed as 
          ! global atom index
          ! A observable code
          ! B observable code
          ! flat correlator data (see correlator.F90) 

        s = buffer_index
        buffer_index = buffer_index + 1
        global_index = INT(buffer(buffer_index))

        buffer_index = buffer_index + 1
        iA = INT(buffer(buffer_index))

        buffer_index = buffer_index + 1
        iB = INT(buffer(buffer_index))

        Call code_to_observable(iA,A)
        Call code_to_observable(iB,B)

        Allocate(tmp_cors(1:(Size(this%correlations)+1)))

        local_index = newatm

        Do i = 1, Size(this%correlations)

          tmp_cors(i) = this%correlations(i)

        End Do

        dim_left = A%dimension()
        dim_right = B%dimension()
        
        ! obtain correlation parameters
        Do i = 1, Size(this%unique_correlations)

          jA = this%unique_correlations(i)%A%id()
          jB = this%unique_correlations(i)%B%id()

          If (iA == jA .and. iB == jB) Then
            blocks = this%unique_correlation_params((i-1)*3+1)
            points = this%unique_correlation_params((i-1)*3+2)
            window = this%unique_correlation_params((i-1)*3+3)
            exit
          End If

        End Do

        ! initialise new correlator
        new_index = Size(tmp_cors)
        Call code_to_observable(iA,tmp_cors(new_index)%correlation%A)
        Call code_to_observable(iB,tmp_cors(new_index)%correlation%B)
        tmp_cors(new_index)%correlation%atom = local_index
        tmp_cors(new_index)%correlation%atom_global = global_index
        Call tmp_cors(new_index)%correlator%init(blocks,points,window,dim_left,dim_right)

        Call move_alloc(tmp_cors,this%correlations)

        ! recieve from buffer
        Call this%correlations(new_index)%correlator%recieve_buffer(buffer,buffer_index)

        buffer_index = s + this%max_buffer_per_atom

        this%number_of_correlations = Size(this%correlations)
 
        Call reindex_correlators(this, config)

      End If

  End Subroutine correlator_recieve

  Subroutine correlator_deport(this, config, buffer, atom_index, buffer_index)
    Class(stats_type),                    Intent(InOut) :: this
    Type(configuration_type),             Intent(InOut) :: config
    Real(Kind=wp), Dimension(:),          Intent(InOut) :: buffer
    Integer,                              Intent(In)    :: atom_index
    Integer,                              Intent(InOut) :: buffer_index
    Type(correlator_holder), Allocatable                :: tmp_cors(:)
    Integer                                             :: i, A, B, j, s, deportations

    If (this%per_atom_correlations .and. this%calculate_correlations) Then


      ! data packed as 
        ! global atom index
        ! A observable code
        ! B observable code
        ! flat correlator data (see correlator.F90) 
    
      j = 0
      deportations = 0
      ! deport to buffer
      
      s = buffer_index

      Do i = 1, Size(this%correlations)
        If (this%correlations(i)%correlation%atom_global == config%ltg(atom_index)) Then
          A = this%correlations(i)%correlation%A%id()
          B = this%correlations(i)%correlation%B%id()
          buffer_index = buffer_index + 1
          buffer(buffer_index) = config%ltg(atom_index)
          buffer_index = buffer_index + 1
          buffer(buffer_index) = A
          buffer_index = buffer_index + 1
          buffer(buffer_index) = B 
          Call this%correlations(i)%correlator%deport_buffer(buffer,buffer_index,.true.)
          deportations = deportations + 1
        Else 
          j = j + 1
        End If
      End Do
      buffer_index = s + deportations*this%max_buffer_per_atom
      ! remove old correlator

      Allocate(tmp_cors(1:j))

      j = 1
      Do i = 1, Size(this%correlations)
        If (this%correlations(i)%correlation%atom_global /= config%ltg(atom_index)) Then
          tmp_cors(j) = this%correlations(i)
          j = j + 1
        End If
      End Do

      Call move_alloc(tmp_cors,this%correlations)

      this%number_of_correlations = Size(this%correlations)

      Call reindex_correlators(this, config)

    End If

  End Subroutine correlator_deport

  Subroutine reindex_correlators(this, config)
    Class(stats_type),                    Intent(InOut)     :: this
    Type(configuration_type),             Intent(In   )     :: config
    Integer                                                 :: i, j, k, &
                                                               cor_global, cor_local, iA, iB
    Logical                                                 :: found
    Logical,                 Allocatable                    :: cors(:)
    Type(correlator_holder), Allocatable                    :: tmp_cors(:)


    k = 0

    Allocate(cors(1:this%number_of_correlations))
    cors = .false.

    Do i = 1,this%number_of_correlations

      cor_global = this%correlations(i)%correlation%atom_global
      cor_local = this%correlations(i)%correlation%atom
      
      found = .false.
      If (cor_local > 0) Then
        If(cor_global /= config%ltg(cor_local)) Then
          Do j = 1, config%natms
            If (config%ltg(j) == cor_global) Then
              this%correlations(i)%correlation%atom = j
              cors(i) = .true.
              k = k + 1
              Exit
            End If
          End Do
        Else 
          cors(i) = .true.
            k = k + 1
        End If
      Else
        cors(i) = .true.
        k = k + 1
      End if

    End Do

    ! check all correlations are correct, ignore any 
    !  that are not tracking a local atom

    Allocate(tmp_cors(1:k))

    k = 1
    Do i = 1, this%number_of_correlations
      If (cors(i)) Then
        tmp_cors(k) = this%correlations(i)
        iA = tmp_cors(k)%correlation%A%id()
        iB = tmp_cors(k)%correlation%B%id()
        k = k + 1
      End If
    End Do

    Call move_alloc(tmp_cors,this%correlations)

    this%number_of_correlations = Size(this%correlations)


  End Subroutine reindex_correlators 

  !!!!!!!! correlators revive !!!!!!!!!

  Subroutine dump_correlations(this, comm, config, unit)
    Class(stats_type),              Intent(InOut) :: this
    Type(comms_type),               Intent(InOut) :: comm
    Type(configuration_type),       Intent(InOut) :: config
    Integer,                        Intent(In   ) :: unit

    Real(Kind=wp),    Allocatable                 :: data_buffer(:), local_buffer(:)
    Integer,          Allocatable                 :: local_ids(:), ids_buffer(:)
    Type(correlator_buffer_type)                  :: packed_correlators
    Type(indices_buffer_type)                     :: packed_ids
    Integer                                       :: i, buffer_size, correlations, &
                                                     buffer_index, A, B, &
                                                     local_correlations, local_buffer_size
    Real(Kind=wp)                                 :: write_sum
  
    ! determine total buffer sizes needed for root

    buffer_size = 0
    correlations = 0
    local_correlations = 0
    local_buffer_size = 0
    If (this%calculate_correlations) Then

      Do i = 1, Size(this%correlations)
        correlations = correlations + 1
        buffer_size = buffer_size + this%correlations(i)%correlator%buffer_size-3
      End Do

      Allocate(local_ids(1:4*correlations))
      Allocate(local_buffer(1:buffer_size))

      ! collect local buffer first
      buffer_index = 0
      Do i = 1, Size(this%correlations)
        A = this%correlations(i)%correlation%A%id()
        B = this%correlations(i)%correlation%B%id()
        local_ids((i-1)*4+1) = A
        local_ids((i-1)*4+2) = B
        If (this%correlations(i)%correlation%atom > 0) Then
          local_ids((i-1)*4+3) = this%correlations(i)%correlation%atom_global
        Else
          local_ids((i-1)*4+3) = 0
        End If
        local_ids((i-1)*4+4) = this%correlations(i)%correlator%buffer_size-3
        Call this%correlations(i)%correlator%deport_buffer(local_buffer,buffer_index)
      End Do

      local_correlations = correlations
      local_buffer_size = buffer_size

    End If

    If (local_correlations == 0) Then
      ! ggatherv will call Size, so must at least
      !   have a dummy allocation
      Allocate(local_ids(0))
      Allocate(local_buffer(0))
    End If

    ! now globally sum sizes
    Call gsum(comm,buffer_size)
    Call gsum(comm,correlations)

    Call packed_ids%initialise(comm, 4*correlations)
    Call packed_correlators%initialise(comm, buffer_size)

    If (comm%idnode == root_id) Then

      Allocate(data_buffer(1:buffer_size))
      Allocate(ids_buffer(1:4*correlations))

    End If

    Call gatherv_scatterv_index_arrays(comm, &
    4*local_correlations, &
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
      write_sum = 0.0
      Do i = 1, correlations
        Write (unit) packed_ids%buffer((i-1)*4+1), &
          packed_ids%buffer((i-1)*4+2), &
          packed_ids%buffer((i-1)*4+3), &
          packed_ids%buffer((i-1)*4+4), &
          packed_correlators%buffer(buffer_index:(buffer_index+packed_ids%buffer((i-1)*4+4)-1))

        write_sum = write_sum + Sum(packed_correlators%buffer(buffer_index:(buffer_index+packed_ids%buffer((i-1)*4+4)-1)))
        buffer_index = buffer_index + packed_ids%buffer((i-1)*4+4)

      End Do

    End If

    Call packed_ids%finalise()
    Call packed_correlators%finalise()

  End Subroutine dump_correlations

  Subroutine revive_correlations(this, comm, config, unit, keyio, no_advance, format)
    Class(stats_type),              Intent(InOut) :: this
    Type(comms_type),               Intent(InOut) :: comm
    Type(configuration_type),       Intent(InOut) :: config
    Integer,                        Intent(In   ) :: unit
    Logical,                        Intent(In   ) :: no_advance
    Character(Len=40),              Intent(In   ) :: format
    Integer,                        Intent(InOut) :: keyio

    Real(Kind=wp),    Allocatable                 :: data_buffer(:), local_buffer(:)
    Integer,          Allocatable                 :: local_ids(:), ids_buffer(:), &
                                                     offests_buffer(:), sizes_buffer(:)
    Type(correlator_buffer_type)                  :: packed_correlators
    Type(indices_buffer_type)                     :: packed_ids
    Integer                                       :: i, buffer_size, correlations, &
                                                     buffer_index, j, packed_index, &
                                                     A, B, atom, offset, local_correlations, &
                                                     local_buffer_size

    local_correlations = 0
    local_buffer_size = 0
    buffer_size = 0
    correlations = 0
    If (this%calculate_correlations) Then                                                 
      Allocate(local_ids(1:4*Size(this%correlations)))

      local_correlations = Size(this%correlations)

      ! determine total buffer sizes needed for root

      ! collect local id buffers
      Do i = 1, Size(this%correlations)
        correlations = correlations + 1
        buffer_size = buffer_size + this%correlations(i)%correlator%buffer_size-3
        A = this%correlations(i)%correlation%A%id()
        B = this%correlations(i)%correlation%B%id()
        local_ids((i-1)*4+1) = A
        local_ids((i-1)*4+2) = B
        If (this%correlations(i)%correlation%atom > 0) Then
          local_ids((i-1)*4+3) = config%ltg(this%correlations(i)%correlation%atom)
        Else
          local_ids((i-1)*4+3) = 0
        End If
        local_ids((i-1)*4+4) = this%correlations(i)%correlator%buffer_size-3

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

    Call packed_ids%initialise(comm, 4*correlations)
    Call packed_correlators%initialise(comm, buffer_size)

    If (comm%idnode == root_id) Then

      Allocate(data_buffer(1:buffer_size))
      Allocate(sizes_buffer(1:correlations))
      Allocate(ids_buffer(1:4*correlations))
      Allocate(offests_buffer(1:correlations))
    
    End If

    Call gatherv_scatterv_index_arrays(comm, &
    4*local_correlations, &
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
              ids_buffer((i-1)*4+1), ids_buffer((i-1)*4+2), ids_buffer((i-1)*4+3), &
              ids_buffer((i-1)*4+4), &
              data_buffer(buffer_index:(buffer_index+ids_buffer((i-1)*4+4)-1))
          Else
            Read (Unit=unit, IOStat=keyio) &
              ids_buffer((i-1)*4+1), ids_buffer((i-1)*4+2), ids_buffer((i-1)*4+3), &
              ids_buffer((i-1)*4+4), &
              data_buffer(buffer_index:(buffer_index+ids_buffer((i-1)*4+4)-1))
          End If
          sizes_buffer(i) = ids_buffer((i-1)*4+4)
          buffer_index = buffer_index + sizes_buffer(i)

      End Do

    End If

    If (comm%idnode == root_id) Then

      ! now setup arrays for scattering
      packed_correlators%buffer = data_buffer

      buffer_index = 1

      Do i = 1, correlations

        A = packed_ids%buffer((i-1)*4+1)
        B = packed_ids%buffer((i-1)*4+2)
        atom = packed_ids%buffer((i-1)*4+3)

        packed_index = -1

        ! find where this data should be placed
        Do j = 1, correlations
          If (ids_buffer((j-1)*4+1) == A .and. &
              ids_buffer((j-1)*4+2) == B .and. &
              ids_buffer((j-1)*4+3) == atom) Then

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
    If (this%calculate_correlations) Then  
      Do i = 1, Size(this%correlations)
        Call this%correlations(i)%correlator%recieve_buffer(local_buffer,buffer_index)
      End Do
    End If

    Call packed_ids%finalise()
    Call packed_correlators%finalise()

  End Subroutine revive_correlations

  !!!!!!!! observables !!!!!!!! 

  Function is_equal(left,right) result (b)
    Class(observable), Intent(In   ) :: left
    Class(observable), Intent(In   ) :: right
    Logical                          :: b
    Integer                          :: lid, rid

    lid = left%id()
    rid = right%id()

    b = lid == rid      
  End Function is_equal

  Subroutine character_to_observable(c, o)
    Character(Len=*),                Intent(In   ) :: c
    Class(observable), Allocatable,  Intent(  Out) :: o
    Logical                                        :: success

    success = .false.
    If (c == velocity_name(observable_velocity()) .or. c == "v") Then 
      Allocate(observable_velocity::o)
      success = .true.
    Else If (c == stress_name(observable_stress()) .or. c == "s") Then
      Allocate(observable_stress::o)
      success = .true.
    Else If (c == heat_flux_name(observable_heat_flux()) .or. c == "hf") Then
      Allocate(observable_heat_flux::o)
      success = .true.
    End If

    If (success .eqv. .false.) Then
      Call error(0,"correlation observable could not be allocated from character: "//c)
    End If

  End Subroutine character_to_observable

  Subroutine code_to_observable(c, o)
    Integer,                        Intent(In   )  :: c
    Class(observable), Allocatable, Intent(  Out)  :: o
    Logical                                        :: success
    Character(Len=100)                             :: msg
    
    success = .false.
    If (c == velocity_id(observable_velocity())) Then 
      Allocate(observable_velocity::o)
      success = .true.
    Else If (c == stress_id(observable_stress())) Then
      Allocate(observable_stress::o)
      success = .true.
    Else If (c == heat_flux_id(observable_heat_flux())) Then
      Allocate(observable_heat_flux::o)
      success = .true.
    End If

    If (success .eqv. .false.) Then
      Write(msg,'(a,i0)') "correlation observable could not be allocated from internal code: ", c
      Call error(0,msg)
    End If

  End Subroutine code_to_observable
  
  !!!!!!!!!! observable_velocity !!!!!!!!!!

  Subroutine velocity_value(t, config, stats, v, atom)
    Class(observable_velocity),   Intent(In   ) :: t
    Type(configuration_type),     Intent(InOut) :: config
    Type(stats_type),             Intent(InOut) :: stats
    Real(Kind=wp), Allocatable,   Intent(InOut) :: v(:)
    Integer, Optional,            Intent(In   ) :: atom
    
    Integer                                     :: d

    d = velocity_dimension(t)

    Allocate(v(1:d))

    If (Present(atom)) Then
      v(1) = config%vxx(atom)
      v(2) = config%vyy(atom)
      v(3) = config%vzz(atom)
    Else
      Call error(0,message="no atom index specified")
    End If
    
  End Subroutine velocity_value

  Function velocity_dimension(t) Result(v)
      Class(observable_velocity), Intent(In   ) :: t
      Integer                                   :: v
      v = 3
  End Function velocity_dimension

  Function velocity_name(t) Result(v)
      Class(observable_velocity), Intent(In   )   :: t
      Character(Len=MAX_CORRELATION_NAME_LENGTH)  :: v
      v = 'velocity'
  End Function velocity_name

  Function velocity_id(t) Result(v)
    Class(observable_velocity),   Intent(In   ) :: t
    Integer                                     :: v
    v = 0
  End Function velocity_id

  Function velocity_per_atom(t) Result(v)
    Class(observable_velocity), Intent(In   ) :: t
    Logical                                   :: v
    v = .true.
  End Function velocity_per_atom

  !!!!!!!!!! observable stress !!!!!!!!!!

  Subroutine stress_value(t, config, stats, v, atom)
    Class(observable_stress),     Intent(In   ) :: t
    Type(configuration_type),     Intent(InOut) :: config
    Type(stats_type),             Intent(InOut) :: stats
    Real(Kind=wp), Allocatable,   Intent(InOut) :: v(:)
    Integer, Optional,            Intent(In   ) :: atom
    
    Integer                                     :: d, i

    d = stress_dimension(t)

    Allocate(v(1:d))

    Do i = 1, 9
      v(i) = stats%strtot(i) / stats%stpvol
    End Do

  End Subroutine stress_value

  Function stress_dimension(t) Result(v)
      Class(observable_stress), Intent(In   ) :: t
      Integer                                 :: v
      v = 9
  End Function stress_dimension

  Function stress_name(t) Result(v)
      Class(observable_stress), Intent(In   )     :: t
      Character(Len=MAX_CORRELATION_NAME_LENGTH)  :: v
      v = 'stress'
  End Function stress_name

  Function stress_id(t) Result(v)
    Class(observable_stress), Intent(In   ) :: t
    Integer                                 :: v
    v = 1
  End Function stress_id

  Function stress_per_atom(t) Result(v)
    Class(observable_stress), Intent(In   ) :: t
    Logical                                 :: v
    v = .false.
  End Function stress_per_atom

    !!!!!!!!!! observable heat_flux !!!!!!!!!!

  Subroutine heat_flux_value(t, config, stats, v, atom)
    Class(observable_heat_flux),     Intent(In   ) :: t
    Type(configuration_type),        Intent(InOut) :: config
    Type(stats_type),                Intent(InOut) :: stats
    Real(Kind=wp), Allocatable,      Intent(InOut) :: v(:)
    Integer, Optional,               Intent(In   ) :: atom
    
    Integer                                        :: d, i

    d = heat_flux_dimension(t)

    Allocate(v(1:d))

    v = stats%heat_flux

  End Subroutine heat_flux_value

  Function heat_flux_dimension(t) Result(v)
      Class(observable_heat_flux),   Intent(In   ) :: t
      Integer                                      :: v
      v = 3
  End Function heat_flux_dimension

  Function heat_flux_name(t) Result(v)
      Class(observable_heat_flux), Intent(In   )     :: t
      Character(Len=MAX_CORRELATION_NAME_LENGTH)     :: v
      v = 'heat_flux'
  End Function heat_flux_name

  Function heat_flux_id(t) Result(v)
    Class(observable_heat_flux), Intent(In   ) :: t
    Integer                                    :: v
    v = 2
  End Function heat_flux_id

  Function heat_flux_per_atom(t) Result(v)
    Class(observable_heat_flux), Intent(In   ) :: t
    Logical                                    :: v
    v = .false.
  End Function heat_flux_per_atom

End Module statistics