Module new_control
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle new style, logical consistent, control file in dlpoly
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins april 2020
  !! contrib - a.m.elena april 2022
  !!-----------------------------------------------------------------------

  Use angles,                   Only: angles_type
  Use angular_distribution,     Only: adf_type
  Use bonds,                    Only: bonds_type
  Use bspline,                  Only: MAX_SPLINES,&
                                      MIN_SPLINES
  Use comms,                    Only: comms_type
  Use configuration,            Only: IMCON_NOPBC,&
                                      IMCON_SLAB,&
                                      configuration_type
  Use constants,                Only: pi,&
                                      tenunt,&
                                      zero_plus
  Use constraints,              Only: constraints_type
  Use control_parameter_module, Only: DATA_BOOL,&
                                      DATA_FLOAT,&
                                      DATA_INT,&
                                      DATA_OPTION,&
                                      DATA_STRING,&
                                      DATA_VECTOR3,&
                                      DATA_VECTOR6,&
                                      control_parameter,&
                                      parameters_hash_table
  Use coord,                    Only: coord_type
  Use core_shell,               Only: core_shell_type
  Use defects,                  Only: defects_type
  Use development,              Only: development_type
  Use dihedrals,                Only: dihedrals_type
  Use electrostatic,            Only: ELECTROSTATIC_COULOMB,&
                                      ELECTROSTATIC_COULOMB_FORCE_SHIFT,&
                                      ELECTROSTATIC_COULOMB_REACTION_FIELD,&
                                      ELECTROSTATIC_DDDP,&
                                      ELECTROSTATIC_EWALD,&
                                      ELECTROSTATIC_NULL,&
                                      electrostatic_type
  Use errors_warnings,          Only: check_print_level,&
                                      error,&
                                      error_units,&
                                      info,&
                                      set_print_level,&
                                      warning
  Use ewald,                    Only: ewald_type
  Use filename,                 Only: &
                                      FILE_CONFIG, FILE_FIELD, FILE_HISTORF, FILE_HISTORY, &
                                      FILE_MSD, FILE_OUTPUT, FILE_RDF, FILE_REVCON, FILE_REVIVE, &
                                      FILE_REVOLD, FILE_STATS, FILE_TABANG, FILE_TABBND, &
                                      FILE_TABDIH, FILE_TABEAM, FILE_TABINV, FILE_TABVDW, file_type
  Use flow_control,             Only: DFTB,&
                                      RESTART_KEY_CLEAN,&
                                      RESTART_KEY_NOSCALE,&
                                      RESTART_KEY_OLD,&
                                      RESTART_KEY_SCALE,&
                                      flow_type
  Use four_body,                Only: four_body_type
  Use greenkubo,                Only: greenkubo_type
  Use hash,                     Only: MAX_KEY
  Use impacts,                  Only: impact_type
  Use inversions,               Only: inversions_type
  Use io,                       Only: &
                                      IO_READ_DIRECT, IO_READ_MASTER, IO_READ_MPIIO, &
                                      IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
                                      IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT, &
                                      IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, &
                                      io_get_parameters, io_set_parameters, io_type
  Use kim,                      Only: kim_type
  Use kinds,                    Only: STR_LEN,&
                                      wp
  Use metal,                    Only: metal_type
  Use minimise,                 Only: MIN_DISTANCE,&
                                      MIN_ENERGY,&
                                      MIN_FORCE,&
                                      minimise_type
  Use mpole,                    Only: POLARISATION_CHARMM,&
                                      POLARISATION_DEFAULT,&
                                      mpole_type
  Use msd,                      Only: msd_type
  Use neighbours,               Only: neighbours_type
  Use numerics,                 Only: dcell,&
                                      seed_type
  Use parse,                    Only: get_line,&
                                      get_word,&
                                      lower_case,&
                                      word_2_real
  Use plumed,                   Only: plumed_type
  Use pmf,                      Only: pmf_type
  Use rdfs,                     Only: rdf_type
  Use rsds,                     Only: rsd_type
  Use statistics,               Only: stats_type
  Use tersoff,                  Only: tersoff_type
  Use thermostat,               Only: &
                                      CONSTRAINT_NONE, CONSTRAINT_SEMI_ORTHORHOMBIC, &
                                      CONSTRAINT_SURFACE_AREA, CONSTRAINT_SURFACE_TENSION, &
                                      DPD_FIRST_ORDER, DPD_NULL, DPD_SECOND_ORDER, &
                                      ENS_NPT_BERENDSEN, ENS_NPT_BERENDSEN_ANISO, &
                                      ENS_NPT_LANGEVIN, ENS_NPT_LANGEVIN_ANISO, ENS_NPT_MTK, &
                                      ENS_NPT_MTK_ANISO, ENS_NPT_NOSE_HOOVER, &
                                      ENS_NPT_NOSE_HOOVER_ANISO, ENS_NVE, ENS_NVT_ANDERSON, &
                                      ENS_NVT_BERENDSEN, ENS_NVT_EVANS, ENS_NVT_GENTLE, &
                                      ENS_NVT_LANGEVIN, ENS_NVT_LANGEVIN_INHOMO, &
                                      ENS_NVT_NOSE_HOOVER, thermostat_type
  Use three_body,               Only: threebody_type
  Use timer,                    Only: timer_type
  Use trajectory,               Only: TRAJ_KEY_COMPRESSED,&
                                      TRAJ_KEY_COORD,&
                                      TRAJ_KEY_COORD_VEL,&
                                      TRAJ_KEY_COORD_VEL_FORCE,&
                                      trajectory_type
  Use ttm,                      Only: ttm_type
  Use units,                    Only: atomic_units,&
                                      convert_units,&
                                      current_units => out_units,&
                                      hartree_units,&
                                      internal_units,&
                                      set_out_units,&
                                      set_timestep,&
                                      si_units,&
                                      units_scheme
  Use vdw,                      Only: MIX_FENDER_HALSEY,&
                                      MIX_FUNCTIONAL,&
                                      MIX_HALGREN,&
                                      MIX_HOGERVORST,&
                                      MIX_LORENTZ_BERTHELOT,&
                                      MIX_NULL,&
                                      MIX_TANG_TOENNIES,&
                                      MIX_WALDMAN_HAGLER,&
                                      vdw_type
  Use z_density,                Only: z_density_type

  Implicit None

  Private

  !> Max number of params, increase to reduce hash collisions
  Integer, Parameter :: PARAMS_TABLE_SIZE = 500

  Real(kind=wp), Parameter :: REAL_TOL = 1.0e-6_wp
!  Character(Len=7), Parameter :: FMT_REAL="1p g12.5e2"

  Public :: initialise_control
  Public :: read_new_control
  Public :: read_io
  Public :: read_ttm
  Public :: read_ensemble
  Public :: read_forcefield
  Public :: read_devel
  Public :: read_structure_analysis
  Public :: read_units
  Public :: read_bond_analysis
  Public :: read_run_parameters
  Public :: read_system_parameters
  Public :: write_parameters

  ! For unit testing
  Public :: parse_control_file
Contains

  Subroutine read_new_control(control_file, params, comm, can_parse)
    !!-----------------------------------------------------------------------
    !!
    !! Set read in a new_style control file
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type(file_type),             Intent(InOut) :: control_file
    Type(parameters_hash_table), Intent(  Out) :: params
    Type(comms_type),            Intent(InOut) :: comm
    Logical,                     Intent(  Out) :: can_parse

    Integer :: ierr

    Call initialise_control(params)
    Open (newunit=control_file%unit_no, file=control_file%filename, status='old', action='read', iostat=ierr)
    If (ierr /= 0) Call error(0, 'CONTROL file not found')

    can_parse = try_parse(control_file%unit_no, comm)

    ! Possibly old style
    If (.not. can_parse) Then
      Call control_file%close ()
      Return
    End If

    Call parse_control_file(control_file%unit_no, params, comm)
    Call control_file%close ()

  End Subroutine read_new_control

  Subroutine read_devel(params, devel, tmr, seed)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(development_type),      Intent(InOut) :: devel
    Type(timer_type),            Intent(InOut) :: tmr
    Type(seed_type),             Intent(InOut) :: seed

    Integer                     :: print_level
    Real(kind=wp), Dimension(3) :: vtmp

    Call params%retrieve('print_level', print_level)
    Call set_print_level(print_level)
    Call params%retrieve('unsafe_comms', devel%l_fast)

    Call params%retrieve('output_energy', devel%l_eng)
    Call params%retrieve('io_write_ascii_revive', devel%l_rout)
    Call params%retrieve('io_read_ascii_revold', devel%l_rin)
    Call params%retrieve('initial_minimum_separation', devel%r_dis)

    Call params%retrieve('timer_depth', tmr%max_depth)
    Call params%retrieve('timer_per_mpi', tmr%proc_detail)

    Call params%retrieve('time_job', tmr%job)
    If (tmr%job < 0.0_wp) tmr%job = Huge(1.0_wp)

    Call params%retrieve('time_close', tmr%clear_screen)
    If (tmr%clear_screen < 0.0_wp) tmr%clear_screen = 0.01_wp * tmr%job

    If (params%is_set('random_seed')) Then
      Call params%retrieve('random_seed', vtmp(1:3))
      Call seed%init(Nint(vtmp(1:3)))
    End If

    Call params%retrieve('dftb_test', devel%test_dftb_library)

  End Subroutine read_devel

  Subroutine read_io(params, io_data, files, comm)
    !!-----------------------------------------------------------------------
    !!
    !! Read in the io parameters
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type(parameters_hash_table),   Intent(In   ) :: params
    Type(io_type),                 Intent(InOut) :: io_data
    Type(file_type), Dimension(:), Intent(InOut) :: files
    Type(comms_type),              Intent(InOut) :: comm

    Integer, Parameter :: MAX_BATCH_SIZE = 10000000, MAX_BUFFER_SIZE = 100000

    Character(Len=STR_LEN) :: curr_option, message
    Integer                :: io_read, io_write, itmp
    Logical                :: ltmp
    Real(kind=wp)          :: rtmp

    Call params%retrieve('io_read_method', curr_option)

    Select Case (curr_option)
    Case ('mpiio')
      io_read = IO_READ_MPIIO

    Case ('direct')
      io_read = IO_READ_DIRECT

    Case ('master')
      io_read = IO_READ_MASTER

    Case Default
      Call bad_option('io_read_method', curr_option)

    End Select

    Call io_set_parameters(io_data, user_method_read=io_read)

    Select Case (io_read)
    Case (IO_READ_MPIIO, IO_READ_DIRECT)
      ! Need to calculate number of readers
      Call params%retrieve('io_read_readers', itmp)

      If (itmp == 0) Then
        rtmp = Min(Real(comm%mxnode, wp), 2.0_wp * Real(comm%mxnode, wp)**0.5_wp)
        itmp = 2**Int(Nearest(Log(rtmp) / Log(2.0_wp), +1.0_wp))
        Do While (Mod(comm%mxnode, itmp) /= 0)
          itmp = itmp - 1
        End Do
      Else If (itmp < comm%mxnode) Then
        Do While (Mod(comm%mxnode, itmp) /= 0)
          itmp = itmp - 1
        End Do
      Else If (itmp > comm%mxnode) Then
        rtmp = Min(Real(comm%mxnode, wp), 2.0_wp * Real(comm%mxnode, wp)**0.5_wp)
        itmp = 2**Int(Nearest(Log(rtmp) / Log(2.0_wp), +1.0_wp))
        Do While (Mod(comm%mxnode, itmp) /= 0)
          itmp = itmp - 1
        End Do
      Else
        Call error(0, 'Cannot have negative number of I/O readers')
      End If

      ! the number of readers is now ready to set

      Call io_set_parameters(io_data, user_n_io_procs_read=itmp)

      ! Sort read batch size
      ! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
      ! Note zero or negative values indicate use the default

      Call params%retrieve('io_read_batch_size', itmp)

      Select Case (itmp)
      Case (0)
        Call io_get_parameters(io_data, user_batch_size_read=itmp)

      Case (1:)
        itmp = Min(itmp, MAX_BATCH_SIZE)
        Call io_set_parameters(io_data, user_batch_size_read=itmp)

      Case Default
        Call error(0, 'Cannot have negative I/O read batch size')

      End Select

    Case (IO_READ_MASTER)
      Write (message, '(a,i10)') 'I/O readers (enforced) ', 1
      Call info(message, .true.)

    Case Default
      rtmp = Min(Real(comm%mxnode, wp), 2.0_wp * Real(comm%mxnode, wp)**0.5_wp)
      itmp = 2**Int(Nearest(Log(rtmp) / Log(2.0_wp), +1.0_wp))
      ! the number of readers is now ready to set
      Call io_set_parameters(io_data, user_n_io_procs_read=itmp)

    End Select

    ! Get read buffer size

    Call params%retrieve('io_read_buffer_size', itmp)

    Select Case (itmp)
    Case (0)
      Call io_get_parameters(io_data, user_buffer_size_read=itmp)

    Case (1:99)
      Call io_set_parameters(io_data, user_buffer_size_read=100)

    Case (100:MAX_BUFFER_SIZE)
      Call io_set_parameters(io_data, user_buffer_size_read=itmp)

    Case (MAX_BUFFER_SIZE + 1:)
      Call io_set_parameters(io_data, user_buffer_size_read=MAX_BUFFER_SIZE)

    Case Default
      Call error(0, 'Negative read buffer size not valid, min = 1')

    End Select

    ! Get parallel read error checking

    If (io_read /= IO_READ_MASTER) Then
      Call params%retrieve('io_read_error_check', ltmp)
      Call io_set_parameters(io_data, user_error_check=ltmp)

    End If

    ! Get write settings

    Call params%retrieve('io_write_method', curr_option)
    Call params%retrieve('io_write_sorted', ltmp)

    Select Case (curr_option)
    Case ('mpiio')
      If (ltmp) Then
        io_write = IO_WRITE_SORTED_MPIIO
      Else
        io_write = IO_WRITE_UNSORTED_MPIIO
      End If

    Case ('direct')
      If (ltmp) Then
        io_write = IO_WRITE_SORTED_DIRECT
      Else
        io_write = IO_WRITE_UNSORTED_DIRECT
      End If

    Case ('master')

      If (ltmp) Then
        io_write = IO_WRITE_SORTED_MASTER
      Else
        io_write = IO_WRITE_UNSORTED_MASTER
      End If

    Case Default
      Call bad_option('io_write_method', curr_option)

    End Select

    ! the write method and type are now ready to set

    Call io_set_parameters(io_data, user_method_write=io_write)

    Select Case (io_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT)
      Call params%retrieve('io_write_writers', itmp)

      If (itmp == 0) Then
        rtmp = Min(Real(comm%mxnode, wp), 8.0_wp * Real(comm%mxnode, wp)**0.5_wp)
        itmp = 2**Int(Nearest(Log(rtmp) / Log(2.0_wp), +1.0_wp))
        Do While (Mod(comm%mxnode, itmp) /= 0)
          itmp = itmp - 1
        End Do

      Else If (itmp < comm%mxnode) Then
        Do While (Mod(comm%mxnode, itmp) /= 0)
          itmp = itmp - 1
        End Do

      Else If (itmp > comm%mxnode) Then
        rtmp = Min(Real(comm%mxnode, wp), 8.0_wp * Real(comm%mxnode, wp)**0.5_wp)
        itmp = 2**Int(Nearest(Log(rtmp) / Log(2.0_wp), +1.0_wp))
        Do While (Mod(comm%mxnode, itmp) /= 0)
          itmp = itmp - 1
        End Do

      Else
        Call error(0, 'Cannot have negative number of I/O writers')

      End If

      ! the number of writers is now ready to set

      Call io_set_parameters(io_data, user_n_io_procs_write=itmp)

      Call params%retrieve('io_write_batch_size', itmp)

      Select Case (itmp)
      Case (0)
        Call io_get_parameters(io_data, user_batch_size_write=itmp)

      Case (1:)
        itmp = Min(itmp, MAX_BATCH_SIZE)
        Call io_set_parameters(io_data, user_batch_size_write=itmp)

      Case Default
        Call error(0, 'Cannot have negative I/O write batch size')

      End Select

    Case (IO_WRITE_UNSORTED_MASTER, IO_WRITE_SORTED_MASTER)
      Continue

    Case Default
      rtmp = Min(Real(comm%mxnode, wp), 8.0_wp * Real(comm%mxnode, wp)**0.5_wp)
      itmp = 2**Int(Nearest(Log(rtmp) / Log(2.0_wp), +1.0_wp))
      ! the number of writers is now ready to set
      Call io_set_parameters(io_data, user_n_io_procs_write=itmp)

    End Select

    Call params%retrieve('io_write_buffer_size', itmp)

    Select Case (itmp)
    Case (0)
      Call io_get_parameters(io_data, user_buffer_size_write=itmp)

    Case (1:99)
      Call io_set_parameters(io_data, user_buffer_size_write=100)

    Case (100:MAX_BUFFER_SIZE)
      Call io_set_parameters(io_data, user_buffer_size_write=itmp)

    Case (MAX_BUFFER_SIZE + 1:)
      Call io_set_parameters(io_data, user_buffer_size_write=MAX_BUFFER_SIZE)

    Case Default
      Call error(0, 'Negative write buffer size not valid, min = 1')
    End Select

    ! switch error checking flag for writing

    If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
      Call params%retrieve('io_write_error_check', ltmp)
      Call io_set_parameters(io_data, user_error_check=ltmp)
    End If

    Call params%retrieve('io_file_output', curr_option)
    files(FILE_OUTPUT)%filename = curr_option
    Call params%retrieve('io_file_config', curr_option)
    files(FILE_CONFIG)%filename = curr_option
    Call params%retrieve('io_file_field', curr_option)
    files(FILE_FIELD)%filename = curr_option
    Call params%retrieve('io_file_statis', curr_option)
    files(FILE_STATS)%filename = curr_option
    Call params%retrieve('io_file_history', curr_option)
    files(FILE_HISTORY)%filename = curr_option
    Call params%retrieve('io_file_historf', curr_option)
    files(FILE_HISTORF)%filename = curr_option
    Call params%retrieve('io_file_revive', curr_option)
    files(FILE_REVIVE)%filename = curr_option
    Call params%retrieve('io_file_revcon', curr_option)
    files(FILE_REVCON)%filename = curr_option
    Call params%retrieve('io_file_revold', curr_option)
    files(FILE_REVOLD)%filename = curr_option
    Call params%retrieve('io_file_rdf', curr_option)
    files(FILE_RDF)%filename = curr_option
    Call params%retrieve('io_file_msd', curr_option)
    files(FILE_MSD)%filename = curr_option
    Call params%retrieve('io_file_tabbnd', curr_option)
    files(FILE_TABBND)%filename = curr_option
    Call params%retrieve('io_file_tabang', curr_option)
    files(FILE_TABANG)%filename = curr_option
    Call params%retrieve('io_file_tabdih', curr_option)
    files(FILE_TABDIH)%filename = curr_option
    Call params%retrieve('io_file_tabinv', curr_option)
    files(FILE_TABINV)%filename = curr_option
    Call params%retrieve('io_file_tabvdw', curr_option)
    files(FILE_TABVDW)%filename = curr_option
    Call params%retrieve('io_file_tabeam', curr_option)
    files(FILE_TABEAM)%filename = curr_option

  End Subroutine read_io

  Subroutine read_ensemble(params, thermo, nvdw, ttm_active)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(thermostat_type),       Intent(InOut) :: thermo
    Integer,                     Intent(In   ) :: nvdw
    Logical,                     Intent(In   ) :: ttm_active

    Character(Len=STR_LEN) :: option
    Logical                :: ltmp

    Call params%retrieve('ensemble', option, required=.true.)
    If (ttm_active .and. option /= 'ttm') Call error(0, 'TTM requested, ensemble not ttm')
    Select Case (option)
    Case ('nve', 'pmf')
      thermo%ensemble = ENS_NVE

    Case ('nvt')

      Call params%retrieve('ensemble_method', option, required=.true.)

      Select Case (option)
      Case ('evans')
        thermo%ensemble = ENS_NVT_EVANS

      Case ('langevin')

        thermo%ensemble = ENS_NVT_LANGEVIN

        Call params%retrieve('ensemble_thermostat_friction', thermo%chi, .true.)

      Case ('andersen')

        thermo%ensemble = ENS_NVT_ANDERSON

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_thermostat_softness', thermo%soft)

      Case ('berendsen')

        thermo%ensemble = ENS_NVT_BERENDSEN

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)

      Case ('hoover', 'nose', 'nose-hoover')

        thermo%ensemble = ENS_NVT_NOSE_HOOVER

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)

      Case ('gentle', 'gst')

        thermo%ensemble = ENS_NVT_GENTLE

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_thermostat_friction', thermo%gama)

      Case ('ttm')

        thermo%ensemble = ENS_NVT_LANGEVIN_INHOMO
        Call params%retrieve('ttm_e-phonon_friction', thermo%chi_ep)
        Call params%retrieve('ttm_e-stopping_friction', thermo%chi_es, .true.)
        Call params%retrieve('ttm_e-stopping_velocity', thermo%vel_es2)
        thermo%vel_es2 = thermo%vel_es2 * thermo%vel_es2 ! square of cutoff velocity for inhomogeneous Langevin thermostat and ttm

      Case ('dpd')
        thermo%ensemble = ENS_NVE

        Call info('Ensemble: NVT dpd (Dissipative Particle Dynamics)', .true.)

        Call params%retrieve('ensemble_dpd_order', option)
        Select Case (option)
        Case ('first', '1')
          thermo%key_dpd = DPD_FIRST_ORDER
        Case ('second', '2')
          thermo%key_dpd = DPD_SECOND_ORDER
        Case default
          Call bad_option('ensemble_dpd_order', option)
        End Select

        ! ALLOCATE DPD ARRAYS
        Call thermo%init_dpd(nvdw + 1) ! Account for extra one in bounds?!

        Call params%retrieve('ensemble_dpd_drag', thermo%gamdpd(0))

      Case default
        Call bad_option('NVT ensemble_method', option)
      End Select

    Case ('npt')

      thermo%variable_cell = .true.

      Call params%retrieve('ensemble_method', option)
      Select Case (option)
      Case ('langevin')

        thermo%ensemble = ENS_NPT_LANGEVIN
        thermo%l_langevin = .true.

        Call params%retrieve('ensemble_thermostat_friction', thermo%chi, .true.)
        Call params%retrieve('ensemble_barostat_friction', thermo%tai, .true.)

      Case ('berendsen')

        thermo%ensemble = ENS_NPT_BERENDSEN

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case ('hoover', 'nose', 'nose-hoover')

        thermo%ensemble = ENS_NPT_NOSE_HOOVER

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case ('mtk')

        thermo%ensemble = ENS_NPT_MTK

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case default
        Call bad_option('NPT ensemble_method', option)
      End Select

    Case ('nst')

      thermo%variable_cell = .true.
      thermo%anisotropic_pressure = .true.

      Call params%retrieve('ensemble_method', option)
      Select Case (option)
      Case ('langevin')

        thermo%ensemble = ENS_NPT_LANGEVIN_ANISO
        thermo%l_langevin = .true.

        Call params%retrieve('ensemble_thermostat_friction', thermo%chi, .true.)
        Call params%retrieve('ensemble_barostat_friction', thermo%tai, .true.)

      Case ('berendsen')

        thermo%ensemble = ENS_NPT_BERENDSEN_ANISO

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case ('hoover', 'nose', 'nose-hoover')

        thermo%ensemble = ENS_NPT_NOSE_HOOVER_ANISO

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case ('mtk')

        thermo%ensemble = ENS_NPT_MTK_ANISO

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case default
        Call bad_option('NST ensemble method', option)
      End Select

      ! Semi isotropic ensembles

      Call params%retrieve('ensemble_semi_isotropic', option)
      Call params%retrieve('ensemble_semi_orthorhombic', ltmp)

      Select Case (option)
      Case ('off')
        Continue
      Case ('area')
        thermo%iso = CONSTRAINT_SURFACE_AREA
      Case ('surface', 'tension')

        thermo%iso = CONSTRAINT_SURFACE_TENSION

        Call params%retrieve('ensemble_tension', thermo%tension, required=.true.)

        thermo%tension = thermo%tension / tenunt

        If (ltmp) Then
          thermo%iso = CONSTRAINT_SEMI_ORTHORHOMBIC
        End If

      Case ('ortho', 'orthorhombic')

        thermo%iso = CONSTRAINT_SURFACE_TENSION
        If (ltmp) Then
          thermo%iso = CONSTRAINT_SEMI_ORTHORHOMBIC
        End If
      Case default
        Call bad_option('ensemble_semi_isotropic', option)
      End Select

    Case default
      Call bad_option('ensemble', option)
    End Select

  End Subroutine read_ensemble

  Subroutine read_bond_analysis(params, flow, bond, angle, dihedral, inversion, max_grid_analysis)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(flow_type),             Intent(InOut) :: flow
    Type(bonds_type),            Intent(InOut) :: bond
    Type(angles_type),           Intent(InOut) :: angle
    Type(dihedrals_type),        Intent(InOut) :: dihedral
    Type(inversions_type),       Intent(InOut) :: inversion
    Integer,                     Intent(  Out) :: max_grid_analysis

    Real(kind=wp), Parameter :: minimum_bond_anal_length = 2.5_wp

    Integer :: itmp, itmp2
    Logical :: ltmp

    Call params%retrieve('analyse_bonds', flow%analyse_bond)
    Call params%retrieve('analyse_angles', flow%analyse_ang)
    Call params%retrieve('analyse_dihedrals', flow%analyse_dih)
    Call params%retrieve('analyse_inversions', flow%analyse_inv)
    Call params%retrieve('analyse_all', ltmp)
    If (ltmp) Then
      flow%analyse_bond = .true.
      flow%analyse_ang = .true.
      flow%analyse_dih = .true.
      flow%analyse_inv = .true.
    End If

    ! Set global
    Call params%retrieve('analyse_frequency', itmp)
    Call params%retrieve('analyse_num_bins', itmp2)

    If (flow%analyse_bond) Then
      Call params%retrieve('analyse_max_dist', bond%rcut)
      If (bond%rcut < minimum_bond_anal_length) Then
        bond%rcut = minimum_bond_anal_length
        Call warning('Bond less than minimum length, setting to minimum (2.5 Ang)', .true.)
      End If
      bond%bin_pdf = itmp2
      If (params%is_set('analyse_num_bins_bonds')) &
        Call params%retrieve('analyse_num_bins_bonds', bond%bin_pdf)

      Call params%retrieve('analyse_frequency_bonds', flow%freq_bond)
    Else
      bond%bin_pdf = -1
      flow%freq_bond = -1
    End If

    If (flow%analyse_ang) Then
      angle%bin_adf = itmp2
      If (params%is_set('analyse_num_bins_angles')) &
        Call params%retrieve('analyse_num_bins_angles', angle%bin_adf)
      Call params%retrieve('analyse_frequency_angles', flow%freq_angle)
    Else
      angle%bin_adf = -1
      flow%freq_angle = -1
    End If

    If (flow%analyse_dih) Then
      dihedral%bin_adf = itmp2
      If (params%is_set('analyse_num_bins_dihedrals')) &
        Call params%retrieve('analyse_num_bins_dihedrals', dihedral%bin_adf)
      Call params%retrieve('analyse_frequency_dihedrals', flow%freq_dihedral)
    Else
      dihedral%bin_adf = -1
      flow%freq_dihedral = -1
    End If

    If (flow%analyse_inv) Then
      inversion%bin_adf = itmp2
      If (params%is_set('analyse_num_bins_inversions')) &
        Call params%retrieve('analyse_num_bins_inversions', inversion%bin_adf)

      Call params%retrieve('analyse_frequency_inversions', flow%freq_inversion)
    Else
      inversion%bin_adf = -1
      flow%freq_inversion = -1
    End If

    If (Any([flow%analyse_bond, flow%analyse_ang, flow%analyse_dih, flow%analyse_inv])) Then
      max_grid_analysis = Max(bond%bin_pdf, angle%bin_adf, dihedral%bin_adf, inversion%bin_adf)
    End If

    flow%freq_bond = Max(1, itmp, flow%freq_bond)
    flow%freq_angle = Max(1, itmp, flow%freq_angle)
    flow%freq_dihedral = Max(1, itmp, flow%freq_dihedral)
    flow%freq_inversion = Max(1, itmp, flow%freq_inversion)

  End Subroutine read_bond_analysis

  Subroutine read_structure_analysis(params, stats, msd_data, rdf, vaf, zden, adf, coords, traj, defect, displacement)

    Type(parameters_hash_table), Intent(In   ) :: params
    Type(stats_type),            Intent(InOut) :: stats
    Type(msd_type),              Intent(InOut) :: msd_data
    Type(rdf_type),              Intent(InOut) :: rdf
    Type(greenkubo_type),        Intent(InOut) :: vaf
    Type(z_density_type),        Intent(InOut) :: zden
    Type(adf_type),              Intent(InOut) :: adf
    Type(coord_type),            Intent(InOut) :: coords
    Type(trajectory_type),       Intent(InOut) :: traj
    Type(defects_type),          Intent(InOut) :: defect(:)
    Type(rsd_type),              Intent(InOut) :: displacement

    Character(Len=STR_LEN) :: option
    Logical                :: ltmp

    ! MSD
    Call params%retrieve('msd_calculate', msd_data%l_msd)

    If (msd_data%l_msd) Then
      Call params%retrieve('msd_start', msd_data%start)
      Call params%retrieve('msd_frequency', msd_data%freq)
    Else If (params%is_any_set([Character(13) :: 'msd_frequency', 'msd_start'])) Then
      Call warning('msd_start or msd_frequency found without msd_calculate')
    End If

    ! VAF
    Call params%retrieve('vaf_calculate', vaf%l_collect)

    If (vaf%l_collect) Then
      Call params%retrieve('vaf_frequency', vaf%freq)
      If (vaf%freq <= 0) vaf%freq = 50
      Call params%retrieve('vaf_binsize', vaf%binsize)
      If (vaf%binsize <= 0) vaf%binsize = Merge(2 * vaf%freq, 100, vaf%freq >= 100)
      Call params%retrieve('vaf_print', vaf%l_print)

      vaf%samp = Ceiling(Real(vaf%binsize, wp) / Real(vaf%freq, wp))

    Else If (params%is_any_set([Character(13) :: 'vaf_frequency', 'vaf_binsize', 'vaf_print'])) Then
      Call warning('vaf_print, vaf_frequency or vaf_binsize found without vaf_calculate')
    End If

    ! RDF
    Call params%retrieve('rdf_calculate', rdf%l_collect)
    ! Necessary for some reason
    Call params%retrieve('rdf_binsize', rdf%rbin)

    If (rdf%l_collect) Then
      Call params%retrieve('rdf_error_analysis', option)

      Select Case (option)
      Case ('jackknife')
        rdf%l_errors_jack = .true.
      Case ('block')
        rdf%l_errors_block = .true.
      Case ('off')
        Continue
      Case default
        Call bad_option('rdf_error_analysis', option)
      End Select

      If (rdf%l_errors_jack .or. rdf%l_errors_block) &
           & Call params%retrieve('rdf_error_analysis_blocks', rdf%num_blocks)

      Call params%retrieve('rdf_frequency', rdf%freq)
      Call params%retrieve('rdf_print', rdf%l_print)

    Else If (params%is_any_set([Character(13) :: 'rdf_frequency', 'rdf_binsize', 'rdf_print'])) Then
      Call warning('rdf_print, rdf_frequency or rdf_binsize found without rdf_calculate')
    End If

    Call params%retrieve('print_probability_distribution', stats%lpana)
    If (stats%lpana .and. .not. rdf%l_print) Then
      Call warning('RDF printing triggered due to a PDA printing request', .true.)
      rdf%l_print = stats%lpana
    End If

    ! ZDen

    Call params%retrieve('zden_calculate', zden%l_collect)

    If (zden%l_collect) Then

      Call params%retrieve('zden_frequency', zden%frequency)
      Call params%retrieve('zden_binsize', zden%bin_width)
      Call params%retrieve('zden_print', zden%l_print)

    Else If (params%is_any_set([Character(14) :: 'zden_frequency', 'zden_binsize', 'zden_print'])) Then
      Call warning('zden_print, zden_frequency or zden_binsize found without zden_calculate')
    End If

    ! ADF

    Call params%retrieve('adf_calculate', adf%adfon)

    If (adf%adfon) Then

      Call params%retrieve('adf_frequency', adf%interval)
      Call params%retrieve('adf_precision', adf%prec)

    Else If (params%is_any_set([Character(13) :: 'adf_frequency', 'adf_precision'])) Then
      Call warning('adf_frequency or adf_precision found without adf_calculate')
    End If

    ! Coord

    Call params%retrieve('coord_calculate', coords%coordon)

    If (coords%coordon) Then
      Call params%retrieve('coord_start', coords%coordstart)
      Call params%retrieve('coord_interval', coords%coordinterval)
      Call params%retrieve('coord_ops', option)
      Select Case (option)
      Case ("icoord")
        coords%coordops = 0
      Case ("ccoord")
        coords%coordops = 1
      Case ("full")
        coords%coordops = 2
      Case Default
        Call bad_option("coord_ops", option)
      End Select

    Else If (params%is_any_set([Character(14) :: 'coord_start', 'coord_interval', 'coord_ops'])) Then
      Call warning('coord_start, coord_interval or coord_ops found without coord_calculate')
    End If

    ! Trajectory

    Call params%retrieve('traj_calculate', traj%ltraj)

    If (traj%ltraj) Then
      Call params%retrieve('traj_start', traj%start)
      Call params%retrieve('traj_interval', traj%freq)
      Call params%retrieve('traj_key', option)

      Select Case (option)
      Case ('pos')
        traj%key = TRAJ_KEY_COORD
      Case ('pos-vel')
        traj%key = TRAJ_KEY_COORD_VEL
      Case ('pos-vel-force')
        traj%key = TRAJ_KEY_COORD_VEL_FORCE
      Case ('compressed')
        traj%key = TRAJ_KEY_COMPRESSED
      Case default
        Call bad_option('traj_key', option)
      End Select

      ! Need to dealias
      Call traj%init((traj%key - 1), (traj%freq), (traj%start))

    Else If (params%is_any_set([Character(13) :: 'traj_start', 'traj_interval', 'traj_key'])) Then
      Call warning('traj_start, traj_interval or traj_key found without traj_calculate')
    End If

    Call params%retrieve('defects_calculate', defect(1)%ldef)

    If (defect(1)%ldef) Then

      Call params%retrieve('defects_start', defect(1)%nsdef)
      Call params%retrieve('defects_interval', defect(1)%isdef)
      Call params%retrieve('defects_distance', defect(1)%rdef)
      ! if (defects%rdef < )

      defect(1)%newjob = .true.
      ! Name REFERENCE and DEFECTS files
      defect(1)%reffile = 'REFERENCE'
      defect(1)%deffile = 'DEFECTS'

      Call params%retrieve('defects_backup', ltmp)
      If (ltmp) Then
        defect(2)%ldef = .true.
        defect(2)%nsdef = defect(1)%nsdef
        defect(2)%isdef = defect(1)%isdef
        defect(2)%rdef = defect(1)%rdef
        defect(2)%newjob = .true.
        defect(2)%reffile = 'REFERENCE1'
        defect(2)%deffile = 'DEFECTS1'
      End If

    Else If (params%is_any_set([Character(16) :: 'defects_start', 'defects_interval', 'defects_distance'])) Then
      Call warning('defects_start, defects_interval or defects_distance found without defects_calculate')
    End If

    Call params%retrieve('displacements_calculate', displacement%lrsd)

    If (displacement%lrsd) Then

      Call params%retrieve('displacements_start', displacement%nsrsd)
      Call params%retrieve('displacements_interval', displacement%isrsd)
      Call params%retrieve('displacements_distance', displacement%rrsd)
      If (displacement%rrsd < 0.15_wp) Then
        displacement%rrsd = 0.15_wp
        Call warning('Displacement_distance too small, reset to 0.15 ang')
      End If
    Else If (params%is_any_set([Character(22) :: 'displacements_start', 'displacements_interval', 'displacements_distance'])) Then
      Call warning('displacements_start, displacements_interval or displacements_distance found without displacements_calculate')
    End If

  End Subroutine read_structure_analysis

  Subroutine read_ttm(params, ttm)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(ttm_type),              Intent(  Out) :: ttm

    Character(Len=STR_LEN) :: option
    Logical                :: ltmp

    Call params%retrieve('ttm_calculate', ttm%l_ttm)
    If (.not. ttm%l_ttm) Return

    Call params%retrieve('ttm_num_ion_cells', ttm%ntsys(3))
    Call params%retrieve('ttm_num_elec_cells', ttm%eltsys)
    Call params%retrieve('ttm_metal', ttm%ismetal)

    Call params%retrieve('ttm_heat_cap_model', option)
    Select Case (option)
    Case ('constant')
      ttm%cetype = 0
      Call params%retrieve('ttm_heat_cap', ttm%ce0)
    Case ('tanh')
      ttm%cetype = 1
      Call params%retrieve('ttm_heat_cap', ttm%sh_A)
      Call params%retrieve('ttm_temp_term', ttm%sh_B)

      If (ttm%sh_A <= zero_plus .or. ttm%sh_B <= zero_plus) &
        Call error(0, 'Electronic specific heat not fully specified')
    Case ('linear')
      ttm%cetype = 2
      Call params%retrieve('ttm_heat_cap', ttm%Cemax)
      Call params%retrieve('ttm_fermi_temp', ttm%Tfermi)

      If (ttm%Tfermi <= zero_plus .or. ttm%Cemax <= zero_plus) &
        Call error(0, 'Electronic specific heat not fully specified')
    Case ('tabulated')
      ttm%cetype = 3
    Case default
      Call bad_option('ttm_heat_cap_model', option)
    End Select

    Call params%retrieve('ttm_dens_model', option)
    Select Case (option)
    Case ('constant')
      ttm%ttmdyndens = .false.
      Call params%retrieve('ttm_dens', ttm%cellrho, .true.)
      If (ttm%cellrho <= zero_plus) Call error(0, 'Bad ttm_dens (<= 0)')
      ttm%rcellrho = 1.0_wp / ttm%cellrho

      ! Rescale
      ttm%sh_A = ttm%sh_A * ttm%cellrho
      ttm%Cemax = ttm%Cemax * ttm%cellrho
      ttm%epc_to_chi = convert_units(1.0e-12_wp * ttm%rcellrho / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')

    Case ('dynamic')
      ttm%ttmdyndens = .true.
      ttm%CeType = ttm%CeType + 4
      ttm%epc_to_chi = convert_units(1.0e-12_wp / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')

    Case default
      Call bad_option('ttm_dens_model', option)
    End Select

    If (ttm%ismetal) Then
      ttm%detype = 0

      Call params%retrieve('ttm_elec_cond_model', option, required=.true.)
      Select Case (option)
      Case ('infinite')
        ttm%ketype = 0
      Case ('constant')
        ttm%ketype = 1
        Call params%retrieve('ttm_elec_cond', ttm%ka0)
        If (ttm%ka0 <= zero_plus) &
          Call error(0, 'Electronic thermal conductivity not fully specified')
      Case ('drude')
        ttm%ketype = 2
        Call params%retrieve('ttm_elec_cond', ttm%ka0)

        If (ttm%ka0 <= zero_plus) &
          Call error(0, 'Electronic thermal conductivity not fully specified')
      Case ('tabulated')
        ttm%ketype = 3
      Case default
        Call bad_option('ttm_elec_cond_model', option)
      End Select
    Else
      ttm%ketype = 0

      Call params%retrieve('ttm_diff_model', option, required=.true.)
      Select Case (option)
      Case ('constant')
        ttm%detype = 1
        Call params%retrieve('ttm_diff', ttm%diff0)
        If (ttm%diff0 <= zero_plus) &
          Call error(0, 'Thermal diffusivity of non-metal not specified')
      Case ('recip', 'reciprocal')
        ttm%detype = 2
        Call params%retrieve('ttm_diff', ttm%diff0)
        Call params%retrieve('ttm_fermi_temp', ttm%Tfermi)
        If (ttm%diff0 <= zero_plus .or. ttm%Tfermi <= zero_plus) &
          Call error(0, 'Thermal diffusivity of non-metal not specified')
      Case ('tabulated')
        ttm%detype = 3
      Case default
        Call bad_option('ttm_diff_model', option)
      End Select

    End If

    Call params%retrieve('ttm_variable_ep', option, required=.true.)
    Select Case (option)
    Case ('homo')
      ttm%gvar = 1
    Case ('hetero')
      ttm%gvar = 2
    Case default
      Call bad_option('ttm_variable_ep', option)
    End Select

    Call params%retrieve('ttm_com_correction', option)
    Select Case (option)
    Case ('full')
      ttm%ttmthvel = .true.
      ttm%ttmthvelz = .false.
    Case ('zdir')
      ttm%ttmthvel = .true.
      ttm%ttmthvelz = .true.
    Case ('off')
      ttm%ttmthvel = .false.
      ttm%ttmthvelz = .false.
    Case default
      Call bad_option('ttm_com_correction', option)
    End Select

    Call params%retrieve('ttm_redistribute', ttm%redistribute)

    Call params%retrieve('ttm_min_atoms', ttm%amin)
    ttm%amin = Max(ttm%amin, 1) ! minimum number of atoms for ttm ionic temperature cell

    Call params%retrieve('ttm_stopping_power', ttm%dedx)

    Call params%retrieve('ttm_spatial_dist', option)
    Select Case (option)
    Case ('gaussian')
      ttm%sdepoType = 1
      Call params%retrieve('ttm_spatial_sigma', ttm%sig)
      Call params%retrieve('ttm_spatial_cutoff', ttm%sigmax)

    Case ('flat')
      ttm%sdepoType = 2

    Case ('laser')
      Call params%retrieve('ttm_laser_type', option)
      Call params%retrieve('ttm_fluence', ttm%fluence)
      Call params%retrieve('ttm_penetration_depth', ttm%pdepth)

      Select Case (option)
      Case ('flat')
        ttm%sdepoType = 2

      Case ('exponential')
        ttm%sdepoType = 3

      Case default
        Call bad_option('ttm_laser_type', option)
      End Select

    Case default
      Call bad_option('ttm_spatial_dist', option)
    End Select

    Call params%retrieve('ttm_temporal_dist', option)
    Select Case (option)
    Case ('gaussian')
      ttm%tdepotype = 1
      Call params%retrieve('ttm_temporal_duration', ttm%tdepo)
      Call params%retrieve('ttm_temporal_cutoff', ttm%tcdepo)

    Case ('exponential')
      ttm%tdepotype = 2
      Call params%retrieve('ttm_temporal_duration', ttm%tdepo)
      Call params%retrieve('ttm_temporal_cutoff', ttm%tcdepo)

    Case ('delta')
      ttm%tdepotype = 3

    Case ('square')
      ttm%tdepotype = 4
      Call params%retrieve('ttm_temporal_duration', ttm%tdepo)
      If (ttm%tdepo <= zero_plus) Then
        ttm%tdepoType = 3
      End If

    Case default
      Call bad_option('ttm_temporal_dist', option)
    End Select

    Call params%retrieve('ttm_boundary_condition', option)
    Call params%retrieve('ttm_boundary_xy', ltmp)
    Select Case (option)
    Case ('periodic')
      ttm%bctypee = 1

    Case ('dirichlet')
      If (ltmp) Then
        ttm%bcTypeE = 4
      Else
        ttm%bcTypeE = 2
      End If

    Case ('neumann')
      ttm%bcTypeE = 3
    Case ('robin')

      Call params%retrieve('ttm_boundary_heat_flux', ttm%fluxout)

      If (ltmp) Then
        ttm%bcTypeE = 6
      Else
        ttm%bcTypeE = 5
      End If

    Case default
      Call bad_option('ttm_boundary_condition', option)
    End Select

    Call params%retrieve('ttm_time_offset', ttm%ttmoffset)
    Call params%retrieve('ttm_oneway', ttm%oneway)
    Call params%retrieve('ttm_stats_frequency', ttm%ttmstats)
    Call params%retrieve('ttm_traj_frequency', ttm%ttmtraj)

  End Subroutine read_ttm

  Subroutine read_units(params)
    Type(parameters_hash_table), Intent(In   ) :: params

    Character(Len=STR_LEN) :: option
    Real(kind=wp)          :: test
    Type(units_scheme)     :: out_units

    out_units = current_units

    If (params%is_set("io_units_scheme")) Then
      Call params%retrieve("io_units_scheme", option)

      Select Case (option)
      Case ('internal')

        out_units = internal_units

      Case ('si')

        out_units = si_units

      Case ('atomic')

        out_units = atomic_units

      Case ('hartree')

        out_units = hartree_units

      End Select
    End If

    If (params%is_set("io_units_length")) &
      Call params%retrieve("io_units_length", out_units%length)
    If (params%is_set("io_units_time")) &
      Call params%retrieve("io_units_time", out_units%time)
    If (params%is_set("io_units_mass")) &
      Call params%retrieve("io_units_mass", out_units%mass)
    If (params%is_set("io_units_charge")) &
      Call params%retrieve("io_units_charge", out_units%charge)
    If (params%is_set("io_units_energy")) &
      Call params%retrieve("io_units_energy", out_units%energy)
    If (params%is_set("io_units_pressure")) &
      Call params%retrieve("io_units_pressure", out_units%pressure)
    If (params%is_set("io_units_force")) &
      Call params%retrieve("io_units_force", out_units%force)
    If (params%is_set("io_units_velocity")) &
      Call params%retrieve("io_units_velocity", out_units%velocity)
    If (params%is_set("io_units_power")) &
      Call params%retrieve("io_units_power", out_units%power)
    If (params%is_set("io_units_surface_tension")) &
      Call params%retrieve("io_units_surface_tension", out_units%surf_ten)
    If (params%is_set("io_units_emf")) &
      Call params%retrieve("io_units_emf", out_units%emf)

    ! Initialise timestep unit
    Call params%retrieve('timestep', test, required=.true.)
    Call set_timestep(test)

    Call set_out_units(out_units)

  End Subroutine read_units

  Subroutine read_forcefield(params, neigh, config, xhi, yhi, zhi, flow, vdws, electro, ewld, mpoles, cshell, met, &
                             kim_Data, bond, threebody, fourbody, tersoffs)
    !!-----------------------------------------------------------------------
    !!
    !! Read forcefield
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !! contrib - a.m.elena march 2021 - we always need a cell if NO_PBC or SLAB
    !!-----------------------------------------------------------------------
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(neighbours_type),       Intent(InOut) :: neigh
    Type(configuration_type),    Intent(InOut) :: config
    Real(Kind=wp),               Intent(In   ) :: xhi, yhi, zhi
    Type(flow_type),             Intent(InOut) :: flow
    Type(vdw_type),              Intent(InOut) :: vdws
    Type(electrostatic_type),    Intent(InOut) :: electro
    Type(ewald_type),            Intent(InOut) :: ewld
    Type(mpole_type),            Intent(InOut) :: mpoles
    Type(core_shell_type),       Intent(InOut) :: cshell
    Type(metal_type),            Intent(InOut) :: met
    Type(kim_type),              Intent(In   ) :: kim_data
    Type(bonds_type),            Intent(In   ) :: bond
    Type(threebody_type),        Intent(InOut) :: threebody
    Type(four_body_type),        Intent(InOut) :: fourbody
    Type(tersoff_type),          Intent(In   ) :: tersoffs

    Real(kind=wp), Parameter :: minimum_rcut = 1.0_wp

    Character(Len=STR_LEN)       :: message, option
    Integer                      :: SPLINE_LIMITS
    Logical                      :: metals_on, vdws_on
    Real(Kind=wp)                :: cut, rtmp, tol, tol1
    Real(Kind=wp), Dimension(10) :: cell_properties

    ! ---------------- SUBCELLING ----------------------------------------------

    Call params%retrieve('subcell_threshold', neigh%pdplnc)
    If (neigh%pdplnc < 1.0_wp) Call error(0, 'subcell_threshold less than minimum (1.0)')

    ! ---------------- VDW SETUP -----------------------------------------------

    Call params%retrieve('vdw_method', option)
    Select Case (option)
    Case ('tabulated')
      Continue
    Case ('direct')
      vdws%l_direct = .true.
    Case ('off')
      vdws%no_vdw = .true.
    Case ('ewald')
      ewld%vdw = .true.
    Case default
      Call bad_option('vdw_method', option)
    End Select
    If (vdws%max_vdw <= 0) vdws%no_vdw = .true.

    vdws_on = .not. vdws%no_vdw
    metals_on = met%max_metal > 0

    Call params%retrieve('vdw_mix_method', option)
    Select Case (option)
    Case ('off')
      vdws%mixing = MIX_NULL
    Case ('lorentz-bethelot')
      vdws%mixing = MIX_LORENTZ_BERTHELOT
    Case ('fender-halsey')
      vdws%mixing = MIX_FENDER_HALSEY
    Case ('hogervorst')
      vdws%mixing = MIX_HOGERVORST
    Case ('halgren')
      vdws%mixing = MIX_HALGREN
    Case ('waldman-hagler')
      vdws%mixing = MIX_WALDMAN_HAGLER
    Case ('tang-tonnies')
      vdws%mixing = MIX_TANG_TOENNIES
    Case ('functional')
      vdws%mixing = MIX_FUNCTIONAL
    Case default
      Call bad_option('vdw_mix_method', option)
    End Select

    Call params%retrieve('vdw_force_shift', vdws%l_force_shift)

    ! ---------------- METALS --------------------------------------------------

    Call params%retrieve('metal_direct', met%l_direct)

    Call params%retrieve('metal_sqrtrho', met%l_emb)

    ! ---------------- CUTOFF --------------------------------------------------

    Call params%retrieve('cutoff', neigh%cutoff, required=.true.)
    If (neigh%cutoff < minimum_rcut) Then
      neigh%cutoff = minimum_rcut
      Call warning('neighbour cutoff less than minimum_cutoff (1.0 Ang), setting to minimum_cutoff', .true.)
    End If

    Call params%retrieve('padding', neigh%padding)
    If (neigh%padding < 0.0_wp) Then
      neigh%padding = 0.0_wp
      Call warning('Bad padding value, reset to 0.0', .true.)
    End If
    flow%reset_padding = params%is_set('padding')

    If (vdws_on) Then
      Call params%retrieve('vdw_cutoff', rtmp)
      If (vdws%cutoff > REAL_TOL .and. vdws%cutoff < rtmp .and. Abs(vdws%cutoff - rtmp) > REAL_TOL) Then
        Write (message, '(a,1p,e12.4)') 'VdW cutoff set by bounds to (Angs): ', vdws%cutoff
        Call warning(message, .true.)

      Else
        vdws%cutoff = rtmp
      End If

      If (vdws%cutoff < minimum_rcut .and. vdws%cutoff > 0.0_wp) Then
        vdws%cutoff = neigh%cutoff
        Call warning('VdW cutoff less than minimum cutoff, setting to global cutoff', .true.)
      End If
      !  there is no  specification of cutoff for vdw in the control

      If (vdws%cutoff < 0.0_wp) Then
        vdws%cutoff = neigh%cutoff
        Call info('# VdW cutoff set to global cutoff', .true.)
      End If

    Else
      vdws%cutoff = 0.0_wp
    End If

    If (metals_on .and. met%rcut < REAL_TOL) Then
      met%rcut = Max(met%rcut, neigh%cutoff, vdws%cutoff)
      Call warning('metal_cutoff not set, setting to max of global cutoff and vdw cutoff', .true.)
    End If

    neigh%cutoff = Max(neigh%cutoff, vdws%cutoff, met%rcut, kim_data%cutoff, bond%rcut, &
                       2.0_wp * tersoffs%cutoff + REAL_TOL)

    ! Sort met%rcut=neigh%cutoff if metal interactions are in play, even if
    ! they are defined by EAM since met%rcut can be /= neigh%cutoff in such
    ! instances, this can break the NLAST check in metal_ld_set_halo

    If (metals_on) met%rcut = neigh%cutoff

    ! Sort vdws%cutoff=neigh%cutoff if VDW interactions are in play

    If (vdws_on .and. vdws%cutoff > neigh%cutoff) Then
      vdws%cutoff = neigh%cutoff
      Call warning('VdW cutoff greater than global cutoff, setting to global cutoff', .true.)
    End If

    If (threebody%mxtbp > 0 .and. threebody%cutoff < REAL_TOL) Then
      Call warning('three body cutoff not set, setting to half of global cutoff', .true.)
      threebody%cutoff = 0.5_wp * neigh%cutoff
    End If

    If (fourbody%max_four_body > 0 .and. fourbody%cutoff < REAL_TOL) Then
      Call warning('four body cutoff not set, setting to half of global cutoff', .true.)
      fourbody%cutoff = 0.5_wp * neigh%cutoff
    End If

    !  sort out cell for NON-PBC and slab...
    cut = neigh%cutoff + 1e-6_wp

    If (config%imcon == IMCON_NOPBC) Then
      config%cell = 0.0_wp
      config%cell(1) = Max(2.0_wp * xhi + cut, 3.0_wp * cut, config%cell(1))
      config%cell(5) = Max(2.0_wp * yhi + cut, 3.0_wp * cut, config%cell(5))
      config%cell(9) = Max(2.0_wp * zhi + cut, 3.0_wp * cut, config%cell(9))
    Else If (config%imcon == IMCON_SLAB) Then
      config%cell(9) = Max(2.0_wp * zhi + cut, 3.0_wp * cut, config%cell(9))
    End If

    ! ---------------- POLARISATION --------------------------------------------

    Call params%retrieve('polarisation_model', option)
    Select Case (option)
    Case ('charmm')
      If (cshell%mxshl > 0 .and. mpoles%max_mpoles > 0) mpoles%key = POLARISATION_CHARMM
      If (mpoles%max_mpoles == 0 .or. cshell%mxshl == 0) &
        Call error(0, 'CHARMM polarisation selected with no mpoles or shells')

      electro%lecx = .true.
      Call params%retrieve('polarisation_thole', mpoles%thole)

    Case ('default')
      mpoles%key = POLARISATION_DEFAULT
    Case default
      Call bad_option('polarisation_model', option)
    End Select

    If (cshell%mxshl > 0) Then
      Call params%retrieve('rlx_tol', cshell%rlx_tol(1))
      Call params%retrieve('rlx_cgm_step', cshell%rlx_tol(2))
      If (cshell%rlx_tol(1) < 1.0_wp) Call error(0, 'Relaxed shell CGM tolerance < 1.0')
    End If

    ! ---------------- ELECTROSTATICS ------------------------------------------

    Call params%retrieve('coul_method', option)

    Select Case (option)
    Case ('off')
      electro%no_elec = .true.
      ! reinitialise multipolar electrostatics indicators
      mpoles%max_mpoles = 0
      mpoles%max_order = 0
      mpoles%key = POLARISATION_DEFAULT
      electro%key = ELECTROSTATIC_NULL
    Case ('ewald', 'spme')
      electro%key = ELECTROSTATIC_EWALD
      ewld%active = .true.

    Case ('dddp')
      electro%key = ELECTROSTATIC_DDDP

    Case ('pairwise')

      electro%key = ELECTROSTATIC_COULOMB

    Case ('force_shifted')

      electro%key = ELECTROSTATIC_COULOMB_FORCE_SHIFT

    Case ('reaction_field')

      electro%key = ELECTROSTATIC_COULOMB_REACTION_FIELD

    Case default

      Call bad_option('coul_method', option)

    End Select

    ! If it's been forcibly set by polarisation_model
    If (.not. electro%lecx) Call params%retrieve('coul_extended_exclusion', electro%lecx)

    Call params%retrieve('coul_dielectric_constant', electro%eps)

    If (params%is_set([Character(14) :: 'coul_damping', 'coul_precision'])) Then
      Call error(0, 'Both damping and precision set')

    Else If (params%is_set('coul_damping')) Then
      Call params%retrieve('coul_damping', electro%damping)

    Else If (params%is_set('coul_precision')) Then
      Call params%retrieve('coul_precision', rtmp)
      rtmp = Max(Min(rtmp, 0.5_wp), 1.0e-20_wp)
      tol = Sqrt(Abs(Log(rtmp * neigh%cutoff)))
      electro%damping = Sqrt(Abs(Log(rtmp * neigh%cutoff * tol))) / neigh%cutoff

    End If

    If (electro%damping > zero_plus) Then
      Call info('Fennell damping applied', .true.)
      If (neigh%cutoff < 12.0_wp) Call warning(7, neigh%cutoff, 12.0_wp, 0.0_wp)
    End If
    If (ewld%active) Then

      Call params%retrieve('ewald_nsplines', ewld%bspline%num_splines)
      ! Check splines in min
      SPLINE_LIMITS = MIN_SPLINES + mpoles%max_order
      If (ewld%bspline%num_splines < SPLINE_LIMITS) Then
        Write (message, '(a,i0,a,i0,a)') &
          "Number of bsplines (", ewld%bspline%num_splines, ") less than minimum permitted (", &
          SPLINE_LIMITS, ") -- Resetting"
        Call Warning(message, .true.)
        ewld%bspline%num_splines = Max(ewld%bspline%num_splines, SPLINE_LIMITS)
      End If

      ! Check splines even
      If (Mod(ewld%bspline%num_splines, 2) /= 0) &
        ewld%bspline%num_splines = ewld%bspline%num_splines + 1

      !Check splines in max
      If (ewld%bspline%num_splines > MAX_SPLINES) Then
        Write (message, '(a,i0,a,i0,a)') &
          "Number of bsplines (", ewld%bspline%num_splines, ") greater than maxiumum permitted (", &
          MAX_SPLINES, ")"
        Call error(0, message)
      End If

      Call dcell(config%cell, cell_properties)

      If (params%is_set([Character(15) :: 'ewald_precision', 'ewald_alpha'])) Then

        Call error(0, 'Cannot specify both precision and manual ewald parameters')

      Else If (params%is_set('ewald_alpha')) Then

        Call params%retrieve('ewald_alpha', ewld%alpha)
        If (params%is_set([Character(18) :: 'ewald_kvec', 'ewald_kvec_spacing'])) Then

          Call error(0, 'Cannot specify both explicit k-vec grid and k-vec spacing')
        Else If (params%is_set('ewald_kvec')) Then

          Call params%retrieve('ewald_kvec', ewld%kspace%k_vec_dim_cont)
        Else

          Call params%retrieve('ewald_kvec_spacing', rtmp)
          ewld%kspace%k_vec_dim_cont = Nint(rtmp / cell_properties(7:9))
        End If

        ! Sanity check for ill defined ewald sum parameters 1/8*2*2*2 == 1
        tol = ewld%alpha * Real(Product(ewld%kspace%k_vec_dim_cont), wp)
        If (Int(tol) < 1) Call error(9)

      Else

        Call params%retrieve('ewald_precision', ewld%precision)

        tol = Sqrt(Abs(Log(ewld%precision * neigh%cutoff)))
        ewld%alpha = Sqrt(Abs(Log(ewld%precision * neigh%cutoff * tol))) / neigh%cutoff
        tol1 = Sqrt(-Log(ewld%precision * neigh%cutoff * (2.0_wp * tol * ewld%alpha)**2))

        ewld%kspace%k_vec_dim_cont = 2 * Nint(0.25_wp + cell_properties(7:9) * ewld%alpha * tol1 / pi)

      End If
    End If

  End Subroutine read_forcefield

  Subroutine read_run_parameters(params, flow, thermo, stats, read_config_indices)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(flow_type),             Intent(InOut) :: flow
    Type(thermostat_type),       Intent(InOut) :: thermo
    Type(stats_type),            Intent(InOut) :: stats
    Logical,                     Intent(  Out) :: read_config_indices

    Logical :: ltmp

    Call params%retrieve('timestep', thermo%tstep, required=.true.)
    If (thermo%tstep < zero_plus) Call error(0, 'Timestep too small')
    Call params%retrieve('timestep_variable', thermo%lvar)

    If (thermo%key_dpd /= DPD_NULL .and. thermo%lvar) Then
      thermo%lvar = .false.
      Call warning('Variable timestep unavalable in DPD themostats', .true.)
    Else If (thermo%lvar) Then
      Call params%retrieve('timestep_variable_min_dist', thermo%mndis)
      Call params%retrieve('timestep_variable_max_dist', thermo%mxdis)
      Call params%retrieve('timestep_variable_max_delta', thermo%mxstp)

      If (thermo%mxstp > zero_plus) Then
        thermo%tstep = Min(thermo%tstep, thermo%mxstp)
      Else
        thermo%mxstp = Huge(1.0_wp)
      End If
      If (thermo%mxdis < 2.5_wp * thermo%mndis .or. thermo%mndis <= 0.0_wp) Then
        Call warning(140, thermo%mndis, thermo%mxdis, 0.0_wp)
        Call error(518)
      End If
    End If

    Call params%retrieve('stack_size', stats%mxstak)

    Call params%retrieve('ignore_config_indices', read_config_indices)

    Call params%retrieve('strict_checks', flow%strict)

    Call params%retrieve('time_run', flow%run_steps, .true.)

    Call params%retrieve('time_equilibration', flow%equil_steps)

    Call params%retrieve('record_equilibration', ltmp)
    flow%equilibration = .not. ltmp

    Call params%retrieve('data_dump_frequency', flow%freq_restart)

    Call params%retrieve('print_frequency', flow%freq_output)

    Call params%retrieve('stats_frequency', stats%intsta)

    Call params%retrieve('print_topology_info', flow%print_topology)

    Call params%retrieve('currents_calculate', stats%cur%on)

  End Subroutine read_run_parameters

  Subroutine read_system_parameters(params, flow, config, thermo, impa, minim, plume, cons, pmf, ttm_active)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(flow_type),             Intent(InOut) :: flow
    Type(configuration_type),    Intent(InOut) :: config
    Type(thermostat_type),       Intent(InOut) :: thermo
    Type(impact_type),           Intent(InOut) :: impa
    Type(minimise_type),         Intent(InOut) :: minim
    Type(plumed_type),           Intent(InOut) :: plume
    Type(constraints_type),      Intent(InOut) :: cons
    Type(pmf_type),              Intent(In   ) :: pmf
    Logical,                     Intent(In   ) :: ttm_active

    Character(Len=STR_LEN)      :: option
    Logical                     :: ltmp, stat
    Real(Kind=wp)               :: rtmp
    Real(Kind=wp), Dimension(6) :: vtmp
    Type(control_parameter)     :: param

    Call params%retrieve('title', option)
    config%sysname = option(1:72)

    ! ---------------- PHYSICAL PROPERTIES -------------------------------------

    Call params%retrieve('temperature', thermo%temp)

    If (params%num_set([Character(22) :: 'pressure_tensor', 'pressure_hydrostatic', 'pressure_perpendicular']) > 1) Then
      Call error(0, 'Multiple pressure specifications')
      ! else if (params%num_set([Character(22) :: 'pressure_tensor', 'pressure_hydrostatic', 'pressure_perpendicular']) == 0) then
      !   call error(0, 'No pressure specification')
    End If

    thermo%stress = 0.0_wp
    If (params%is_set('pressure_tensor')) Then
      Call params%retrieve('pressure_tensor', vtmp)
      thermo%stress(1) = vtmp(1)
      thermo%stress(5) = vtmp(2)
      thermo%stress(9) = vtmp(3)
      thermo%stress(2) = vtmp(4)
      thermo%stress(4) = vtmp(4)
      thermo%stress(3) = vtmp(5)
      thermo%stress(7) = vtmp(5)
      thermo%stress(6) = vtmp(6)
      thermo%stress(8) = vtmp(6)

    Else If (params%is_set('pressure_perpendicular')) Then
      Call params%retrieve('pressure_perpendicular', vtmp(1:3))
      thermo%stress(1:9:4) = vtmp(1:3)
    Else
      Call params%retrieve('pressure_hydrostatic', rtmp)
      thermo%stress(1:9:4) = rtmp
    End If

    If (.not. thermo%anisotropic_pressure .or. thermo%iso /= CONSTRAINT_NONE) Then
      thermo%press = Sum(thermo%stress(1:9:4)) / 3.0_wp
      thermo%stress = 0.0_wp
    End If

    Call params%retrieve('fixed_com', config%l_vom)
    If (.not. config%l_vom .and. .not. ttm_active) Then
      Call info('"no fixed_com" option auto-switched on - COM momentum removal will be abandoned', .true.)
      Call warning('this may lead to a build up of the COM momentum ' &
                   //'and a manifestation of the "flying ice-cube" effect', .true.)
    End If

    ! ---------------- INITIALISATION ------------------------------------------

    Call params%retrieve('restart', option)
    Select Case (option)
    Case ('rescale')
      flow%restart_key = RESTART_KEY_SCALE
    Case ('noscale')
      flow%restart_key = RESTART_KEY_NOSCALE
    Case ('continue')
      flow%restart_key = RESTART_KEY_OLD
    Case ('clean')
      flow%restart_key = RESTART_KEY_CLEAN
    Case default
      Call bad_option('restart', option)
    End Select

    If (config%levcfg == 0 .and. flow%restart_key /= RESTART_KEY_CLEAN) Then
      Call warning('CONFIG contains positions only, forcing clean restart')
      flow%restart_key = RESTART_KEY_CLEAN
    End If

    ! ---------------- EQULIBRATION --------------------------------------------

    Call params%retrieve('reset_temperature_interval', thermo%freq_zero)
    thermo%l_zero = thermo%freq_zero > 0
    If (thermo%freq_zero == 0) thermo%freq_zero = flow%equil_steps + 1

    Call params%retrieve('regauss_frequency', thermo%freq_tgaus)
    thermo%l_tgaus = thermo%freq_tgaus > 0
    If (thermo%freq_tgaus == 0) thermo%freq_tgaus = flow%equil_steps + 1

    Call params%retrieve('rescale_frequency', thermo%freq_tscale)
    thermo%l_tscale = thermo%freq_tscale > 0
    If (thermo%freq_tscale == 0) thermo%freq_tscale = flow%equil_steps + 1

    If (params%is_set('equilibration_force_cap')) Then
      flow%force_cap = .true.
      Call params%retrieve('equilibration_force_cap', config%fmax)
    End If

    Call params%retrieve('minimisation_criterion', option)
    minim%minimise = option /= 'off'

    Call params%get('minimisation_tolerance', param)
    Read (param%val, *) minim%tolerance

    Call params%retrieve('minimisation_frequency', minim%freq)
    Call params%retrieve('minimisation_step_length', minim%step_length)

    Select Case (option)
    Case ('off')
      Continue
    Case ('force')
      minim%key = MIN_FORCE

      If (.not. params%is_set('minimisation_tolerance')) Then
        param%units = "internal_f"
      End If

      param%internal_units = "internal_f"
      minim%tolerance = convert_units(minim%tolerance, param%units, param%internal_units, stat)
      If (.not. stat) Call error_units(param%units, param%internal_units, 'Minimisation tolerance wrong type')

      If (minim%tolerance < 1.0_wp .or. minim%tolerance > 1000.0_wp) &
        minim%tolerance = 50.0_wp

    Case ('energy')
      minim%key = MIN_ENERGY

      If (.not. params%is_set('minimisation_tolerance')) Then
        param%units = "internal_e"
      End If

      param%internal_units = "internal_e"
      minim%tolerance = convert_units(minim%tolerance, param%units, param%internal_units, stat)
      If (.not. stat) Call error_units(param%units, param%internal_units, 'Minimisation tolerance wrong type')

      If (minim%tolerance < zero_plus .or. minim%tolerance > 0.01_wp) &
        minim%tolerance = 0.005_wp

    Case ('distance')
      minim%key = MIN_DISTANCE

      If (.not. params%is_set('minimisation_tolerance')) Then
        param%units = "internal_l"
      End If

      param%internal_units = "internal_l"
      minim%tolerance = convert_units(minim%tolerance, param%units, param%internal_units, stat)
      If (.not. stat) Call error_units(param%units, param%internal_units, 'Minimisation tolerance wrong type')

      If (minim%tolerance < REAL_TOL .or. minim%tolerance > 0.1_wp) &
        minim%tolerance = 0.005_wp

    Case default
      Call bad_option('minimisation_criterion', option)
    End Select

    If (minim%freq == 0) minim%freq = flow%equil_steps + 1

    ! ---------------- PSEUDO-THERMOSTAT ---------------------------------------

    Call params%retrieve('pseudo_thermostat_method', option)
    If (option /= 'off') Then

      thermo%l_stochastic_boundaries = .true.
      Call info('pseudo thermostat attached to MD cell boundary', .true.)

      Select Case (option)
      Case ('langevin-direct')
        thermo%key_pseudo = 0
      Case ('langevin')
        thermo%key_pseudo = 1
      Case ('gaussian')
        thermo%key_pseudo = 2
      Case ('direct')
        thermo%key_pseudo = 3
      Case default
        Call bad_option('pseudo_thermostat_method', option)
      End Select

      Call params%retrieve("pseudo_thermostat_width", thermo%width_pseudo)
      If (thermo%width_pseudo < 2.0_wp) Then
        thermo%width_pseudo = 2.0_wp
        Call info('thermostat thickness insufficient - reset to 2 Angs', .true.)
      End If

      Call params%retrieve("pseudo_thermostat_temperature", thermo%temp_pseudo)
      thermo%temp_pseudo = Max(0.0_wp, thermo%temp_pseudo)

    End If

    ! --------------- CONSTRAINTS ----------------------------------------------

    If (cons%mxcons > 0 .or. pmf%mxpmf > 0) Then
      Call params%retrieve('shake_max_iter', cons%max_iter_shake)
      Call params%retrieve('shake_tolerance', cons%tolerance)
    End If

    ! --------------- PLUMED ---------------------------------------------------

    Call params%retrieve('plumed', plume%l_plumed)

    If (plume%l_plumed) Then
      Call params%retrieve('plumed_input', option)
      plume%input = option(1:125)
      Call params%retrieve('plumed_log', option)
      plume%logfile = option(1:125)
      Call params%retrieve('plumed_precision', plume%prec)
      Call params%retrieve('plumed_restart', ltmp)

      If (ltmp) Then
        plume%restart = 1
      Else
        plume%restart = 0
      End If
    End If

    ! --------------- IMPACT ---------------------------------------------------

    impa%active = params%is_any_set([Character(17) :: 'impact_part_index', 'impact_energy', 'impact_time'])
    If (impa%active) Then

      Call params%retrieve('impact_part_index', impa%imd)
      Call params%retrieve('impact_time', impa%tmd)
      Call params%retrieve('impact_energy', impa%emd)
      Call params%retrieve('impact_direction', vtmp(1:3))
      impa%vmx = vtmp(1)
      impa%vmy = vtmp(2)
      impa%vmz = vtmp(3)
    End If

    ! --------------- PER-PARTICLE ---------------------------------------------

    Call params%retrieve('heat_flux', flow%heat_flux)
    Call params%retrieve('write_per_particle', flow%write_per_particle)

    ! --------------- EXPANSION ------------------------------------------------

    Call params%retrieve('nfold', vtmp(1:3))
    config%l_exp = Any(Nint(vtmp(1:3)) > 1)
    If (config%l_exp) Then
      config%nx = Nint(vtmp(1))
      config%ny = Nint(vtmp(2))
      config%nz = Nint(vtmp(3))
    End If

    Call params%retrieve('replay', flow%simulation)
    flow%simulation = .not. flow%simulation
    Call params%retrieve('replay_calculate_forces', flow%replay_recalculate_forces)

  End Subroutine read_system_parameters

  Subroutine write_parameters(io_data, files, neigh, config, link_cell, flow, stats, thermo, ttm, mpoles, vdws, &
                              electro, cshell, ewld, met, impa, minim, plume, cons, pmf, bond, angle, dihedral, inversion, &
                              msd_data, rdf, vaf, zdensity, adf, coords, defect, traj, displacement)
    Type(io_type),            Intent(In   ) :: io_data
    Type(file_type),          Intent(In   ) :: files(:)
    Type(neighbours_type),    Intent(In   ) :: neigh
    Type(configuration_type), Intent(In   ) :: config
    Integer, Dimension(3),    Intent(In   ) :: link_cell
    Type(flow_type),          Intent(In   ) :: flow
    Type(stats_type),         Intent(In   ) :: stats
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(ttm_type),           Intent(In   ) :: ttm
    Type(mpole_type),         Intent(In   ) :: mpoles
    Type(vdw_type),           Intent(In   ) :: vdws
    Type(electrostatic_type), Intent(In   ) :: electro
    Type(core_shell_type),    Intent(In   ) :: cshell
    Type(ewald_type),         Intent(In   ) :: ewld
    Type(metal_type),         Intent(In   ) :: met
    Type(impact_type),        Intent(In   ) :: impa
    Type(minimise_type),      Intent(In   ) :: minim
    Type(plumed_type),        Intent(In   ) :: plume
    Type(constraints_type),   Intent(In   ) :: cons
    Type(pmf_type),           Intent(In   ) :: pmf
    Type(bonds_type),         Intent(In   ) :: bond
    Type(angles_type),        Intent(In   ) :: angle
    Type(dihedrals_type),     Intent(In   ) :: dihedral
    Type(inversions_type),    Intent(In   ) :: inversion
    Type(msd_type),           Intent(In   ) :: msd_data
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(greenkubo_type),     Intent(In   ) :: vaf
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(adf_type),           Intent(In   ) :: adf
    Type(coord_type),         Intent(In   ) :: coords
    Type(defects_type),       Intent(In   ) :: defect(:)
    Type(trajectory_type),    Intent(In   ) :: traj
    Type(rsd_type),           Intent(In   ) :: displacement

    Character(Len=80) :: banner(6)

    Write (banner(1), '(a)') ''
    Write (banner(2), '(a)') '#'//Repeat('*', 79)
    Write (banner(3), '(a4,a72,a4)') '#** ', 'title:'//Repeat(' ', 66), ' ***'
    Write (banner(4), '(a4,a72,a4)') '#** ', config%sysname, ' ***'
    Write (banner(5), '(a)') '#'//Repeat('*', 79)
    Write (banner(6), '(a)') ''
    Call info(banner, 6, .true.)

    If (check_print_level(1)) Call write_io(io_data, files)
    If (check_print_level(1)) Call write_system_parameters(flow, config, stats, thermo, impa, minim, plume, cons, pmf)
    If (check_print_level(1)) Call write_ensemble(thermo)
    If (check_print_level(1)) Call write_forcefield(link_cell, neigh, vdws, electro, ewld, mpoles, cshell, met)
    If (ttm%l_ttm .and. check_print_level(1)) Call write_ttm(thermo, ttm)
    If (check_print_level(1)) Call write_bond_analysis(stats, flow, bond, angle, dihedral, inversion)
    If (check_print_level(1)) &
      Call write_structure_analysis(stats, msd_data, rdf, vaf, zdensity, adf, coords, traj, defect, displacement)
    Call info('', .true.)

  End Subroutine write_parameters

  Subroutine write_io(io_data, files)
    Type(io_type),   Intent(In   ) :: io_data
    Type(file_type), Intent(In   ) :: files(:)

    Character(len=256) :: message

    !    Integer :: io_read, io_write, procs_write, procs_read, batch_write, batch_read

    Call info('', .true.)
    Call info('I/O Parameters: ', .true.)

    Select Case (io_data%method_read)
    Case (IO_READ_MPIIO)
      Call info('  I/O read method: parallel by using MPI-I/O', .true.)

    Case (IO_READ_DIRECT)
      Call info('  I/O read method: parallel by using direct access', .true.)

    Case (IO_READ_MASTER)
      Call info('  I/O read method: serial by using a single master process', .true.)

    End Select

    Write (message, '(a,i0.1)') '  I/O readers: ', io_data%n_io_procs_read
    Call info(message, .true.)

    Select Case (io_data%method_write)
    Case (IO_READ_MPIIO, IO_READ_DIRECT)
      Write (message, '(a,i0.1)') '  I/O read batch size (assumed): ', io_data%batch_size_read
      Call info(message, .true., level=3)
    End Select

    Write (message, '(a,i0.1)') '  I/O read buffer size: ', io_data%buffer_size_read
    Call info(message, .true., level=3)

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_MPIIO)
      Call info('  I/O write method: parallel by using MPI-I/O', .true.)

    Case (IO_WRITE_SORTED_DIRECT, IO_WRITE_UNSORTED_DIRECT)
      Call info('  I/O write method: parallel by using direct access', .true.)
      Call warning('  in parallel this I/O write method has portability issues', .true.)

    Case (IO_WRITE_SORTED_MASTER, IO_WRITE_UNSORTED_MASTER)
      Call info('  I/O write method: serial by using a single master process', .true.)

    End Select

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MASTER, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT)
      Call info('  I/O write type: data sorting ON', .true., level=3)

    Case (IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT)
      Call info('  I/O write type: data sorting OFF', .true.)

    End Select

    Write (message, '(a,i0.1)') '  I/O writers: ', io_data%n_io_procs_write
    Call info(message, .true.)

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_SORTED_DIRECT, &
          IO_WRITE_UNSORTED_DIRECT)
      Write (message, '(a,i10)') '  I/O write batch size (assumed): ', io_data%batch_size_write
      Call info(message, .true., level=3)
    End Select

    Write (message, '(a,i10)') '  I/O write buffer size: ', io_data%buffer_size_write
    Call info(message, .true., level=3)

    If (io_data%global_error_check) Then
      Call info('  I/O parallel error checking ON', .true.)
    Else
      Call info('  I/O parallel error checking OFF', .true., level=3)
    End If

    Call info('', .true.)
    Call info('File outputs: ', .true.)
    If (files(FILE_OUTPUT)%filename /= "OUTPUT") Call info('  OUTPUT file: '//files(FILE_OUTPUT)%filename, .true.)
    If (files(FILE_CONFIG)%filename /= "CONFIG") Call info('  CONFIG file: '//files(FILE_CONFIG)%filename, .true.)
    If (files(FILE_FIELD)%filename /= "FIELD") Call info('  FIELD file: '//files(FILE_FIELD)%filename, .true.)
    If (files(FILE_STATS)%filename /= "STATIS") Call info('  STATIS file: '//files(FILE_STATS)%filename, .true.)
    If (files(FILE_HISTORY)%filename /= "HISTORY") Call info('  HISTORY file: '//files(FILE_HISTORY)%filename, .true.)
    If (files(FILE_HISTORF)%filename /= "HISTORF") Call info('  HISTORF file: '//files(FILE_HISTORF)%filename, .true.)
    If (files(FILE_REVIVE)%filename /= "REVIVE") Call info('  REVIVE file: '//files(FILE_REVIVE)%filename, .true.)
    If (files(FILE_REVCON)%filename /= "REVCON") Call info('  REVCON file: '//files(FILE_REVCON)%filename, .true.)
    If (files(FILE_REVOLD)%filename /= "REVOLD") Call info('  REVOLD file: '//files(FILE_REVOLD)%filename, .true.)
    If (files(FILE_RDF)%filename /= 'RDFDAT') Call info('  RDF file: '//files(FILE_RDF)%filename, .true.)
    If (files(FILE_MSD)%filename /= 'MSDTMP') Call info('  MSD file: '//files(FILE_MSD)%filename, .true.)
    If (files(FILE_TABBND)%filename /= 'TABBND') Call info('  TABBND file: '//files(FILE_TABBND)%filename, .true.)
    If (files(FILE_TABANG)%filename /= 'TABANG') Call info('  TABANG file: '//files(FILE_TABANG)%filename, .true.)
    If (files(FILE_TABDIH)%filename /= 'TABDIH') Call info('  TABDIH file: '//files(FILE_TABDIH)%filename, .true.)
    If (files(FILE_TABINV)%filename /= 'TABINV') Call info('  TABINV file: '//files(FILE_TABINV)%filename, .true.)
    If (files(FILE_TABVDW)%filename /= 'TABVDW') Call info('  TABVDW file: '//files(FILE_TABVDW)%filename, .true.)
    If (files(FILE_TABEAM)%filename /= 'TABEAM') Call info('  TABEAM file: '//files(FILE_TABEAM)%filename, .true.)

  End Subroutine write_io

  Subroutine write_bond_analysis(stats, flow, bond, angle, dihedral, inversion)
    Type(stats_type),      Intent(In   ) :: stats
    Type(flow_type),       Intent(In   ) :: flow
    Type(bonds_type),      Intent(In   ) :: bond
    Type(angles_type),     Intent(In   ) :: angle
    Type(dihedrals_type),  Intent(In   ) :: dihedral
    Type(inversions_type), Intent(In   ) :: inversion

    Character(Len=STR_LEN) :: messages(4)

    Call info('', .true.)
    If (.not. Any([flow%analyse_bond, flow%analyse_ang, flow%analyse_dih, flow%analyse_inv])) Then
      Call info('# No intramolecular distribution collection requested', .true., level=3)
    Else
      Call info('Intramolecular distribution collection requested for:', .true.)

      If (flow%analyse_bond) Then
        Write (messages(1), '(a)') '  Bonds:'
        Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_bond
        Write (messages(3), '(a, i0.1)') '  -- Num samples  (points): ', bond%bin_pdf
        Write (messages(4), '(a, 1p g12.5e2)') '  -- Cutoff         (Angs): ', bond%rcut
        Call info(messages, 4, .true.)
      End If

      If (flow%analyse_ang) Then
        Write (messages(1), '(a)') '  Angles:'
        Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_angle
        Write (messages(3), '(a, i0.1)') '  -- Num samples  (points): ', angle%bin_adf
        Call info(messages, 3, .true.)
      End If

      If (flow%analyse_dih) Then
        Write (messages(1), '(a)') '  Dihedrals:'
        Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_dihedral
        Write (messages(3), '(a, i0.1)') '  -- Num samples  (points): ', dihedral%bin_adf
        Call info(messages, 3, .true.)
      End If

      If (flow%analyse_inv) Then
        Write (messages(1), '(a)') '  Inversions:'
        Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_inversion
        Write (messages(3), '(a, i0.1)') '  -- Num samples  (points): ', inversion%bin_adf
        Call info(messages, 3, .true.)
      End If

    End If

    If (stats%lpana) Then
      Call info('# Probability distribution analysis printing requested', .true.)
    Else
      Call info('# No probability distribution analysis printing requested', .true., level=3)
    End If

  End Subroutine write_bond_analysis

  Subroutine write_structure_analysis(stats, msd_data, rdf, vaf, zdensity, adf, coords, traj, defect, displacement)
    Type(stats_type),      Intent(In   ) :: stats
    Type(msd_type),        Intent(In   ) :: msd_data
    Type(rdf_type),        Intent(InOut) :: rdf
    Type(greenkubo_type),  Intent(In   ) :: vaf
    Type(z_density_type),  Intent(InOut) :: zdensity
    Type(adf_type),        Intent(In   ) :: adf
    Type(coord_type),      Intent(In   ) :: coords
    Type(trajectory_type), Intent(In   ) :: traj
    Type(defects_type),    Intent(In   ) :: defect(:)
    Type(rsd_type),        Intent(In   ) :: displacement

    Character(Len=STR_LEN) :: messages(4)

    Call info('', .true.)

    If (rdf%l_collect) Then
      Write (messages(1), '(a)') 'RDF collection requested:'
      Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', rdf%freq
      Write (messages(3), '(a, 1p g12.5e2)') '  -- Bin size       (Angs): ', rdf%rbin
      Call info(messages, 3, .true.)

      If (rdf%l_print) Then
        Call info('  -- RDF printing requested', .true., level=3)
      Else If (stats%lpana) Then
        Call info('  -- RDF printing triggered due to a PDA printing request', .true.)
      Else
        Call info('  -- No RDF printing requested', .true.)
      End If

      If (rdf%max_rdf == 0) Then
        Call info('  -- No RDF pairs specified in FIELD', .true., level=3)
      Else
        Call info('  -- RDF pairs specified in FIELD', .true.)
      End If
    Else
      Call info('No RDF analysis requested', .true., level=3)
    End If

    If (zdensity%l_collect) Then
      Write (messages(1), '(a)') 'Z-density profiles requested:'
      Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', zdensity%frequency
      Write (messages(3), '(a, 1p g12.5e2)') '  -- Bin size       (Angs): ', zdensity%bin_width
      Call info(messages, 3, .true.)

      If (zdensity%l_print) Then
        Call info('  -- Z-density printing requested', .true., level=3)
      Else
        Call info('  -- No Z-density printing requested', .true.)
      End If

    Else
      Call info('No Z-density analysis requested', .true., level=3)
    End If

    If (vaf%samp > 0) Then
      Write (messages(1), '(a)') 'VAF profiles requested:'
      Write (messages(2), '(a, i0.1)') '  -- Collect every (steps): ', vaf%freq
      Write (messages(3), '(a, 1p g12.5e2)') '  -- Bin size       (Angs): ', vaf%binsize
      Call info(messages, 3, .true.)

      If (vaf%l_print) Then
        Call info('  -- VAF printing requested', .true., level=3)
      Else
        Call info('  -- No VAF printing requested', .true.)
      End If

      If (vaf%l_average) Then
        Call info('  -- Time-averaged VAF profile', .true.)
      Else
        Call info('  -- Instantaneous VAF profiles', .true.)
      End If

    Else
      Call info('No VAF analysis requested', .true., level=3)
    End If

    If (msd_data%l_msd) Then
      Write (messages(1), '(a)') 'MSDTMP profiles requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', msd_data%start
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', msd_data%freq
      Call info(messages, 3, .true.)
    Else
      Call info('No MSD analysis requested', .true., level=3)
    End If

    If (traj%ltraj) Then
      Write (messages(1), '(a)') 'Trajectory recording requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', traj%start
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', traj%freq
      Select Case (traj%key)
      Case (TRAJ_KEY_COORD)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COORD'
      Case (TRAJ_KEY_COORD_VEL)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COORD, VEL'
      Case (TRAJ_KEY_COORD_VEL_FORCE)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COORD, VEL, FORCE'
      Case (TRAJ_KEY_COMPRESSED)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COMPRESSED'
      End Select

      Call info(messages, 4, .true.)
    Else
      Call info('No trajectory recording requested', .true., level=3)
    End If

    If (defect(1)%ldef) Then
      Write (messages(1), '(a)') 'Defects analysis requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', defect(1)%nsdef
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', defect(1)%isdef
      Write (messages(4), '(a,1p g12.5e2)') '  -- Distance condition (Angs): ', defect(1)%rdef
      Call info(messages, 4, .true.)
      If (defect(2)%ldef) Call info('DEFECTS1 file option: ON', .true.)
    Else
      Call info('No defects analysis requested', .true., level=3)
    End If

    If (displacement%lrsd) Then
      Write (messages(1), '(a)') 'Displacements analysis requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', displacement%nsrsd
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', displacement%isrsd
      Write (messages(4), '(a,1p g12.5e2)') '  -- Distance condition (Angs): ', displacement%rrsd
      Call info(messages, 4, .true.)
    Else
      Call info('No displacements analysis requested', .true., level=3)
    End If

    If (coords%coordon) Then
      Write (messages(1), '(a)') 'Coordination analysis requested:'
      Write (messages(2), '(a)') '  -- display to be implemented'
      Call info(messages, 2, .true.)
    Else
      Call info('No coordination analysis requested', .true., level=3)
    End If

    If (adf%adfon) Then
      Write (messages(1), '(a)') 'Angular analysis requested:'
      Write (messages(2), '(a)') '  -- display to be implemented'
      Call info(messages, 2, .true.)
    Else
      Call info('No angular analysis requested', .true., level=3)
    End If

  End Subroutine write_structure_analysis

  Subroutine write_system_parameters(flow, config, stats, thermo, impa, minim, plume, cons, pmf)
    Type(flow_type),          Intent(In   ) :: flow
    Type(configuration_type), Intent(In   ) :: config
    Type(stats_type),         Intent(In   ) :: stats
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(impact_type),        Intent(In   ) :: impa
    Type(minimise_type),      Intent(In   ) :: minim
    Type(plumed_type),        Intent(In   ) :: plume
    Type(constraints_type),   Intent(In   ) :: cons
    Type(pmf_type),           Intent(In   ) :: pmf

    Character(Len=STR_LEN) :: message, messages(5)
    Integer            :: itmp

    Call info('', .true.)
    Call info("System parameters: ", .true.)

    If (.not. thermo%lvar) Then
      Write (message, '(a,1p g12.5e2)') '  Fixed simulation timestep (ps): ', thermo%tstep
      Call info(message, .true.)
    Else
      Write (messages(1), '(a,1p g12.5e2)') '  Variable simulation timestep (ps): ', thermo%tstep
      Write (messages(2), '(a)') '  Controls for variable timestep:'
      Write (messages(3), '(a,1p g12.5e2)') '  -- Minimum distance Dmin (Angs): ', thermo%mndis
      Write (messages(4), '(a,1p g12.5e2)') '  -- Maximum distance Dmax (Angs): ', thermo%mxdis
      Call info(messages, 4, .true.)
      If (thermo%mxstp > zero_plus .and. thermo%mxstp < 1.0e10_wp) Then
        Write (message, '(a,1p g12.5e2)') '  -- Timestep ceiling max step (ps): ', thermo%mxstp
        Call info(message, .true.)
      End If
    End If

    Write (message, '(a,i10)') '  Run Duration (steps): ', flow%run_steps
    If (flow%run_steps > 0) Call info(message, .true.)

    If (flow%equil_Steps > 0) Then
      Write (message, '(a,i10)') '  Equilibration period (steps): ', flow%equil_steps
      Call info(message, .true.)

      If (.not. flow%equilibration) Call info('  -- Equilibration included in overall averages', .true.)

      If (thermo%freq_tgaus > 0) Then
        Call info('  -- Temperature regaussing: ON', .true.)
        Write (message, '(a,i10)') '  -- Temperature regaussing interval (steps): ', thermo%freq_tgaus
        Call info(message, .true.)
      End If

      If (thermo%freq_tscale > 0) Then
        Call info('  -- Temperature scaling: ON', .true.)
        Write (message, '(a,i10)') '  -- Temperature scaling interval (steps): ', thermo%freq_tscale
        Call info(message, .true.)
      End If

      If (flow%force_cap) Then
        Write (messages(1), '(a)') '  -- Force capping: ON'
        Write (messages(2), '(a,1p g12.5e2)') '  -- Force capping limit (kT/Angs): ', config%fmax
        Call info(messages, 2, .true.)
      End If

      If (minim%minimise) Then
        Write (messages(1), '(a)') '  -- Minimisation: ON'
        Select Case (minim%key)
        Case (MIN_FORCE)
          Write (messages(2), '(a,a8)') '   -- Minimisation criterion:         ', "Force"
        Case (MIN_ENERGY)
          Write (messages(2), '(a,a8)') '   -- Minimisation criterion:         ', "Energy"
        Case (MIN_DISTANCE)
          Write (messages(2), '(a,a8)') '   -- Minimisation criterion:         ', "Distance"
        End Select

        Write (messages(3), '(a,i10)') '   -- Minimisation frequency (steps): ', minim%freq
        Write (messages(4), '(a,1p g12.5e2)') '   -- Minimisation tolerance:         ', minim%tolerance
        Write (messages(5), '(a,1p g12.5e2)') '   -- Minimisation CGM step:          ', minim%step_length

        Call info(messages, 5, .true.)
      End If
    End If

    If (stats%mxstak > 0) Then
      Write (message, '(a,i10)') '  Rolling averages length (steps): ', stats%mxstak
      Call info(message, .true.)
    End If

    If (Any([flow%freq_restart > 0, flow%freq_output > 0, stats%intsta > 0,.not. flow%print_topology])) Then
      Call info('  File output info:', .true.)
    End If

    If (flow%freq_restart > 0) Then
      Write (message, '(a,i10)') '  -- Restart dumping interval (steps): ', flow%freq_restart
      Call info(message, .true.)
    End If

    If (flow%freq_output > 0) Then
      Write (message, '(a,i10)') '  -- Data printing interval (steps): ', flow%freq_output
      Call info(message, .true.)
    End If

    If (stats%intsta > 0) Then
      Write (message, '(a,i10)') '  -- Statistics file interval (steps): ', stats%intsta
      Call info(message, .true.)
    End If

    If (.not. flow%print_topology) Then
      Call info('  -- Not printing extended FIELD topology in OUTPUT', .true., level=3)
    Else
      Call info('  -- Printing extended FIELD topology in OUTPUT: ON', .true.)
    End If

    If (stats%cur%on) Call info("  Computing currents: ON", .true.)
    If (flow%heat_flux) Call info("  Computing heat flux: ON", .true.)
    If (flow%write_per_particle) Call info("  Per-particle info: ON", .true.)

    Select Case (flow%restart_key)
    Case (RESTART_KEY_SCALE)
      Call info('  Scaled restart requested (starting a new simulation)', .true.)
    Case (RESTART_KEY_NOSCALE)
      Call info('  Unscaled restart requested (starting a new simulation)', .true.)
    Case (RESTART_KEY_OLD)
      Call info('  Restart requested (continuing an old simulation)', .true.)
      Call warning('  Timestep from REVOLD overides specification in CONTROL', .true.)
    Case (RESTART_KEY_CLEAN)
      If (config%levcfg /= 0) Then
        Call info('  Clean start requested, discarding CONFIG velocities', .true.)
      End If
    End Select

    Write (message, '(a,1p g12.5e2)') '  Simulation temperature (K):  ', thermo%temp
    Call info(message, .true.)

    If (thermo%l_stochastic_boundaries) Then
      Call info('  Pseudo-thermostat attached to MD cell boundary: ON', .true.)

      Select Case (thermo%key_pseudo)
      Case (0)
        Call info('  -- Thermostat control: Langevin + direct temperature scaling', .true.)
      Case (1)
        Call info('  -- Thermostat control: Langevin temperature scaling', .true.)
      Case (2)
        Call info('  -- Thermostat control: gaussian temperature scaling', .true.)
      Case (3)
        Call info('  -- Thermostat control: direct temperature scaling', .true.)
      End Select

      Write (message, '(a,1p g12.5e2)') '  -- Thermostat thickness (Angs): ', thermo%width_pseudo
      Call info(message, .true.)
      Write (message, '(a,1p g12.5e2)') '  -- Thermostat temperature (K): ', thermo%temp_pseudo
      Call info(message, .true.)
    End If

    If (thermo%press > 0.0_wp) Then
      Write (message, '(3A,1p g12.5e2)') '  Simulation pressure (', 'katms', '):', convert_units(thermo%press, 'internal_p', 'katm')
      Call info(message, .true.)
    Else If (Any(thermo%stress > 0.0_wp)) Then
      Write (messages(1), '(3A)') '  Simulation pressure (', 'katms', '):'
      Write (messages(2), '(2X, 3(1p g12.5e2))') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp=1, 3)
      Write (messages(3), '(2X, 3(1p g12.5e2))') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp=4, 6)
      Write (messages(4), '(2X, 3(1p g12.5e2))') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp=7, 9)
      Call info(messages, 4, .true.)
    End If

    If (cons%mxcons > 0 .or. pmf%mxpmf > 0) Then
      Write (messages(1), '(a,i0.1)') '  Iterations for SHAKE/RATTLE: ', cons%max_iter_shake
      Write (messages(2), '(a,1p g12.5e2)') '  Tolerance for SHAKE/RATTLE (Angs): ', cons%tolerance
      Call info(messages, 2, .true.)
    End If

    If (config%l_exp) Then
      Write (message, '(a,9x,3i5)') '  System expansion: ', config%nx, config%ny, config%nz
      Call info(message, .true.)
    End If

    If (config%dvar > 1.0_wp) Then
      Write (message, '(a,9x,1p g12.5e2)') '  Permitted density variance (%): ', &
        convert_units(config%dvar - 1.0_wp, '', '%')
      Call info(message, .true.)
    End If

    If (impa%active) Then
      Call info('  Impact calculation: ON', .true.)
      Write (messages(1), '(a)') ''
      Write (messages(2), '(a,i10)') '  -- Particle (index): ', impa%imd
      Write (messages(3), '(a,i10)') '  -- Timestep (steps): ', impa%tmd
      Write (messages(4), '(a,1p g12.5e2)') '  -- Energy   (keV): ', impa%emd
      Write (messages(5), '(a,3(1p g12.5e2))') '  -- v-r(x,y,z): ', impa%vmx, impa%vmy, impa%vmz

      Call info(messages, 5, .true.)

    End If

    If (plume%l_plumed) Then
      Call info('  Plumed calculation: ON', .true.)
      Write (messages(1), '(a)') ''
      Write (messages(2), '(a,a)') '  -- PLUMED input: ', plume%input
      Write (messages(3), '(a,a)') '  -- PLUMED log: ', plume%logfile
      Write (messages(4), '(a,i0)') '  -- PLUMED precision: ', plume%prec
      Write (messages(5), '(a,i0)') '  -- PLUMED restart: ', plume%restart
      Call info(messages, 5, .true.)
    End If

  End Subroutine write_system_parameters

  Subroutine write_forcefield(link_cell, neigh, vdws, electro, ewld, mpoles, cshell, met)
    Integer, Dimension(3),    Intent(In   ) :: link_cell
    Type(neighbours_type),    Intent(In   ) :: neigh
    Type(vdw_type),           Intent(In   ) :: vdws
    Type(electrostatic_type), Intent(In   ) :: electro
    Type(ewald_type),         Intent(In   ) :: ewld
    Type(mpole_type),         Intent(In   ) :: mpoles
    Type(core_shell_type),    Intent(In   ) :: cshell
    Type(metal_type),         Intent(In   ) :: met

    Character(Len=STR_LEN) :: message, messages(4)

    Call info('', .true.)
    Call info('Forcefield Parameters: ', .true.)

    ! ---------------- CUTOFF --------------------------------------------------

    Write (message, '(a,3i6)') "  Final link-cell decomposition (x,y,z): ", link_cell
    Call info(message, .true., level=1)

    Write (message, '(a,1p g12.5e2)') '  Real space cutoff (Angs): ', neigh%cutoff
    Call info(message, .true.)
    Write (message, '(a,1p g12.5e2)') '  Cutoff padding (Angs): ', neigh%padding
    Call info(message, .true.)

    ! ---------------- SUBCELLING ----------------------------------------------

    Write (message, '(a,1p g12.5e2)') '  Subcelling threshold density: ', neigh%pdplnc
    Call info(message, .true.)

    ! ---------------- VDW SETUP -----------------------------------------------

    If (vdws%l_direct) Then
      Call info('  VdWs: Direct calculation', .true.)
    Else If (vdws%no_vdw) Then
      Call info('  VdWs: Disabled', .true.)
      ! else if (ewld%vdw) then
      !   Call info('VdWs - Ewald', .true.)
    Else
      Call info('  VdWs: Tabulated', .true.)
    End If

    If (.not. vdws%no_vdw) Then
      Write (message, '(a,1p g12.5e2)') '  -- VdW cutoff (Angs): ', vdws%cutoff
      Call info(message, .true.)

      If (vdws%mixing /= MIX_NULL) Then
        Call info('  -- Vdw cross terms mixing opted (for undefined mixed potentials)', .true.)
        Call info('    mixing is limited to potentials of the same type only', .true.)
        Call info('    mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)', .true.)

        Select Case (vdws%mixing)
        Case (MIX_LORENTZ_BERTHELOT)
          Call info('  -- Mixing scheme: Lorentz-Berthelot - e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2', .true.)

        Case (MIX_FENDER_HALSEY)
          Call info('  -- Mixing scheme: Fender-Halsey - e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2', .true.)

        Case (MIX_HOGERVORST)
          Call info('  -- Mixing scheme: Hogervorst (good hope) - ' &
                    //'e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)', .true.)

        Case (MIX_HALGREN)
          Call info('  -- Mixing scheme: Halgren HHG - ' &
                    //'e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)', .true.)

        Case (MIX_WALDMAN_HAGLER)
          Call info('  -- Mixing scheme: Waldman-Hagler - ' &
                    //'e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ;s_ij=[(s_i^6+s_j^6)/2]^(1/6)', .true.)

        Case (MIX_TANG_TOENNIES)
          Call info('  -- Mixing scheme: Tang-Toennies - ' &
                    //' e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13', .true.)
          Call info(Repeat(' ', 43)//'s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)', .true.)

        Case (MIX_FUNCTIONAL)
          Call info('  -- Mixing scheme: Functional - ' &
                    //'e_ij=3 * (e_i*e_j)^(1/2) * (s_i*s_j)^3 / ' &
                    //'SUM_L=0^2{[(s_i^3+s_j^3)^2 / (4*(s_i*s_j)^L)]^(6/(6-2L))}', .true.)
          Call info(Repeat(' ', 40)//'s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}', .true.)
        End Select
      End If

      If (vdws%l_force_shift) Call info('  -- VdW force-shifting: ON', .true.)
    End If

    ! ---------------- METALS --------------------------------------------------

    If (met%l_direct) Call info('  Metal direct: ON', .true.)
    If (met%l_emb) Call info('  Metal sqrtrho: ON', .true.)

    ! ---------------- POLARISATION --------------------------------------------

    If (mpoles%max_mpoles > 0) Then
      Select Case (mpoles%key)
      Case (POLARISATION_CHARMM)
        Call info('  Polarisation method: CHARMM')
      Case (POLARISATION_DEFAULT)
        Call info('  Polarisation method: Default')
      End Select
    End If

    If (cshell%mxshl > 0) Then
      Write (message, '(a,1p g12.5e2)') '  -- Relaxed shell model CGM tolerance (Force): ', cshell%rlx_tol(1)
      Call info(message, .true.)
      If (cshell%rlx_tol(2) > 0.0_wp) Then
        Write (message, '(a,1p g12.5e2)') '  -- Relaxed shell model CGM step (Angs): ', cshell%rlx_tol(2)
        Call info(message, .true.)
      End If
    End If

    ! ---------------- ELECTROSTATICS ------------------------------------------

    Select Case (electro%key)
    Case (ELECTROSTATIC_NULL)
      Call info('  Electrostatics: Disabled', .true.)
    Case (ELECTROSTATIC_EWALD)

      Call info('  Electrostatics: Smooth Particle Mesh Ewald', .true.)

      If (ewld%precision > 0.0_wp) Then
        Write (message, '(a,1p g12.5e2)') '  -- Ewald sum precision: ', ewld%precision
        Call info(message, .true.)
      End If

      Write (messages(1), '(a,1p g12.5e2)') '  -- Ewald convergence parameter (Angs^-1): ', ewld%alpha
      Write (messages(2), '(a,3i5)') '  -- Ewald kmax1 kmax2 kmax3   (x2): ', ewld%kspace%k_vec_dim_cont
      If (Any(ewld%kspace%k_vec_dim /= ewld%kspace%k_vec_dim_cont)) Then
        Write (messages(3), '(a,3i5)') '  -- DaFT adjusted kmax values (x2): ', ewld%kspace%k_vec_dim
        Write (messages(4), '(a,i0.1)') '  -- B-spline interpolation order: ', ewld%bspline%num_splines
        Call info(messages, 4, .true.)
      Else
        Write (messages(3), '(a,i0.1)') '  -- B-spline interpolation order: ', ewld%bspline%num_splines
        Call info(messages, 3, .true.)
      End If

    Case (ELECTROSTATIC_DDDP)
      Call info('  Electrostatics: Distance Dependent Dielectric', .true.)

    Case (ELECTROSTATIC_COULOMB)

      Call info('  Electrostatics: Coulombic Potential', .true.)

    Case (ELECTROSTATIC_COULOMB_FORCE_SHIFT)
      Call info('  Electrostatics: Force-Shifted Coulombic Potential', .true.)

    Case (ELECTROSTATIC_COULOMB_REACTION_FIELD)
      Call info('  Electrostatics: Reaction Field', .true.)
    End Select

    If (electro%key /= ELECTROSTATIC_NULL) Then
      If (Abs(electro%eps - 1.0_wp) > zero_plus) Then
        Write (message, '(a,1p g12.5e2)') '  -- Relative dielectric constant: ', electro%eps
        Call info(message, .true.)
      End If

      !! Fix with electro merge
      If (electro%key /= ELECTROSTATIC_EWALD .and. electro%damping > zero_plus) Then
        Call info('  -- Fennell damping applied', .true.)
        Write (message, '(a,1p g12.5e2)') '  -- Damping parameter (A^-1): ', electro%damping
        Call info(message, .true.)
      End If

      If (electro%lecx) Then
        Call info('  -- Extended Coulombic eXclusion: ON', .true.)
      Else
        Call info('  -- Extended Coulombic eXclusion: OFF', .true.)
      End If

    End If

  End Subroutine write_forcefield

  Subroutine write_ensemble(thermo)
    Type(thermostat_type), Intent(In   ) :: thermo

    Character(Len=STR_LEN)               :: message
    Character(Len=STR_LEN), Dimension(4) :: messages

    Call info('', .true.)
    Call info('Thermostat details:', .true.)

    ensembles:Select Case(thermo%ensemble)
  Case (ENS_NVE)
    Select Case (thermo%key_dpd)
    Case (DPD_NULL)
      Call info('  Ensemble: NVE (Microcanonical)', .true.)
      Exit ensembles
    Case (DPD_FIRST_ORDER)

      Call info('  Ensemble: NVT dpd (Dissipative Particle Dynamics)', .true.)
      Call info("  Ensemble type: Shardlow's first order splitting (S1)", .true.)

    Case (DPD_SECOND_ORDER)

      Call info('  Ensemble: NVT dpd (Dissipative Particle Dynamics)', .true.)
      Call info("  Ensemble type: Shardlow's first order splitting (S2)", .true.)

    End Select

    If (thermo%gamdpd(0) > zero_plus) Then
      Write (message, '(a,1p g12.5e2)') '  -- Drag coefficient (Dalton/ps): ', thermo%gamdpd(0)
      Call info(message, .true.)
    End If

  Case (ENS_NVT_EVANS)

    Call info('  Ensemble: NVT Evans (Isokinetic)', .true.)
    Call info('  Gaussian temperature constraints in use', .true.)

  Case (ENS_NVT_LANGEVIN)

    Call info('  Ensemble: NVT Langevin (Stochastic Dynamics)', .true.)
    Write (message, '(a,1p g12.5e2)') '  -- Thermostat friction (ps^-1): ', thermo%chi
    Call info(message, .true.)

  Case (ENS_NVT_ANDERSON)

    Write (messages(1), '(a)') '  Ensemble: NVT Andersen'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Softness (dimensionless): ', thermo%soft

    Call info(messages, 3, .true.)

  Case (ENS_NVT_BERENDSEN)

    Call info('  Ensemble: NVT Berendsen', .true.)
    Write (message, '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Call info(message, .true.)
    Call warning('If you plan to use the Berendsen thermostat, please read https://doi.org/10.1021/acs.jctc.8b00446', .true.)

  Case (ENS_NVT_NOSE_HOOVER)

    Call info('  Ensemble: NVT Nose-Hoover', .true.)
    Write (message, '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Call info(message, .true.)

  Case (ENS_NVT_GENTLE)

    Write (messages(1), '(a)') '  Ensemble: NVT gentle stochastic thermostat'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Friction on thermostat  (ps^-1): ', thermo%gama
    Call info(messages, 3, .true.)

  Case (ENS_NVT_LANGEVIN_INHOMO)

    Write (messages(1), '(a)') '  Ensemble: NVT inhomogeneous Langevin (Stochastic Dynamics)'
    Write (messages(2), '(a,1p g12.5e2)') '  -- e-phonon friction (ps^-1): ', thermo%chi_ep
    Write (messages(3), '(a,1p g12.5e2)') '  -- e-stopping friction (ps^-1): ', thermo%chi_es
    Write (messages(4), '(a,1p g12.5e2)') '  -- e-stopping velocity (ang ps^-1): ', thermo%vel_es2
    Call info(messages, 4, .true.)

  Case (ENS_NPT_LANGEVIN)

    Write (messages(1), '(a)') '  Ensemble: NPT isotropic Langevin (Stochastic Dynamics)'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat friction (ps^-1): ', thermo%chi
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat friction (ps^-1): ', thermo%tai
    Call info(messages, 3, .true.)

  Case (ENS_NPT_BERENDSEN)

    Write (messages(1), '(a)') '  Ensemble: NPT isotropic Berendsen'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ', thermo%tau_p
    Call info(messages, 3, .true.)

  Case (ENS_NPT_NOSE_HOOVER)

    Write (messages(1), '(a)') '  Ensemble: NPT isotropic Nose-Hoover (Melchionna)'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ', thermo%tau_p
    Call info(messages, 3, .true.)

  Case (ENS_NPT_MTK)

    Write (messages(1), '(a)') '  Ensemble: NPT isotropic Martyna-Tuckerman-Klein'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ', thermo%tau_p
    Call info(messages, 3, .true.)

  Case (ENS_NPT_LANGEVIN_ANISO)

    Write (messages(1), '(a)') '  Ensemble: NPT anisotropic Langevin (Stochastic Dynamics)'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat friction (ps^-1): ', thermo%chi
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat friction (ps^-1): ', thermo%tai
    Call info(messages, 3, .true.)

  Case (ENS_NPT_BERENDSEN_ANISO)

    Write (messages(1), '(a)') '  Ensemble: NPT anisotropic Berendsen'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ', thermo%tau_p
    Call info(messages, 3, .true.)

  Case (ENS_NPT_NOSE_HOOVER_ANISO)

    Write (messages(1), '(a)') '  Ensemble: NPT anisotropic Nose-Hoover (Melchionna)'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ', thermo%tau_p
    Call info(messages, 3, .true.)

  Case (ENS_NPT_MTK_ANISO)

    Write (messages(1), '(a)') '  Ensemble: NPT anisotropic Martyna-Tuckerman-Klein'
    Write (messages(2), '(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ', thermo%tau_t
    Write (messages(3), '(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ', thermo%tau_p
    Call info(messages, 3, .true.)

    End Select ensembles

    ! Semi isotropic ensembles

    Select Case (thermo%iso)
    Case (CONSTRAINT_NONE)
      Continue
    Case (CONSTRAINT_SURFACE_AREA)
      Call info('  Semi-isotropic barostat: constant normal pressure (Pn) &', .true.)
      Call info('         (N-Pn-A-T)      : constant surface area (A)', .true.)
    Case (CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC)

      If (thermo%tension > 0.0_wp) Then
        Write (messages(1), '(a)') '  Semi-isotropic barostat: constant normal pressure (Pn) &'
        Write (messages(2), '(a)') '       (N-Pn-gamma-T)    : constant surface tension (gamma)'
        Write (messages(3), '(a,e11.4)') &
          '  -- Simulation surface tension ('//'dyn/cm'//'): ', convert_units(thermo%tension, 'internal_f/internal_l', 'dyn/cm')
        Call info(messages, 3, .true.)
      Else
        Call info('  Semi-isotropic barostat: orthorhombic MD cell constraints', .true.)
      End If

      If (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) Then
        Call info('  Semi-isotropic barostat: semi-orthorhombic MD cell constraints', .true.)
      End If

    End Select

    If (Any(thermo%iso == [CONSTRAINT_SURFACE_AREA, CONSTRAINT_SURFACE_TENSION])) Then
      Call info('  -- Semi-isotropic ensembles are only correct for infinite ', .true.)
      Call info('       interfaces placed perpendicularly to the z axis', .true.)
    End If

  End Subroutine write_ensemble

  Subroutine write_ttm(thermo, ttm)
    Type(thermostat_type), Intent(In   ) :: thermo
    Type(ttm_type),        Intent(In   ) :: ttm

    Character(Len=STR_LEN) :: message, messages(3)

    Call info('', .true.)
    Call info('TTM Parameters: ', .true.)

    Write (messages(1), '(a,3(1x,i8.1))') '  Ionic temperature grid size (x,y,z): ', ttm%ntsys(1:3)
    Write (messages(2), '(a,3(1x,1p g12.5e2))') '  Temperature grid size (x,y,z): ', ttm%delx, ttm%dely, ttm%delz
    Write (messages(3), '(a,1p g12.5e2)') '  Average number of atoms/cell: ', ttm%sysrho * ttm%volume
    Call info(messages, 3, .true.)
    Write (message, '(a,3(1x,i8.1))') '  Electronic temperature grid size (x,y,z): ', &
      ttm%eltsys(1:3)
    Call info(message, .true.)

    If (ttm%ismetal) Then
      Call info('  Electronic subsystem: Metal (thermal conductivity required)', .true.)
    Else
      Call info('  Electronic subsystem: Non-Metal (thermal diffusivity required)', .true.)
    End If

    If (ttm%ttmdyndens) Then
      Call info('  Dynamic atomic density: ON', .true.)
    Else
      Write (message, '(a,1p g12.5e2)') '  Atomic density (A^-3): ', ttm%cellrho
      Call info(message, .true.)
    End If

    Select Case (ttm%cetype)
    Case (0) !Constant
      Call info('  Electronic specific heat capacity: Constant', .true.)
      Write (message, '(a,1p g12.5e2)') '  -- Electronic S.H.C. (kB/atom): ', ttm%Ce0 / ttm%cellrho
      Call info(message, .true.)

    Case (4) !Constant
      Call info('  Electronic specific heat capacity: Constant', .true.)
      Write (message, '(a,1p g12.5e2)') '  -- Electronic S.H.C. (kB/atom): ', ttm%Ce0
      Call info(message, .true.)

    Case (1) !tanh
      Write (messages(1), '(a)') '  Electronic specific heat capacity: Hyperbolic tangent'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Constant term A (kB/atom): ', ttm%sh_A / ttm%cellrho
      Write (messages(3), '(a,1p g12.5e2)') '  -- Temperature term B (K^-1): ', ttm%sh_B * 1.0e4_wp
      Call info(messages, 3, .true.)

    Case (5) !tanh
      Write (messages(1), '(a)') '  Electronic specific heat capacity: Hyperbolic tangent'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Constant term A (kB/atom): ', ttm%sh_A
      Write (messages(3), '(a,1p g12.5e2)') '  -- Temperature term B (K^-1): ', ttm%sh_B * 1.0e4_wp
      Call info(messages, 3, .true.)

    Case (2, 6) !linear
      Write (messages(1), '(a)') '  Electronic specific heat capacity: Linear'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Max. electronic S.H.C. (kB/atom): ', ttm%Cemax
      Write (messages(3), '(a,1p g12.5e2)') '  -- Fermi temperature (K):', ttm%Tfermi
      Call info(messages, 3, .true.)

    Case (3) !tabulated
      Call info('  Electronic specific heat capacity: Tabulated', .true.)

    End Select

    Select Case (ttm%ketype)
    Case (0) ! Infinite
      Call info('  Electronic thermal conductivity: Infinity', .true.)

    Case (1) ! Constant
      Write (messages(1), '(a)') '  Electronic thermal conductivity: Constant'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Electronic T.C. (W m^-1 K^-1): ', &
        convert_units(ttm%Ka0, 'k_b/ps/A', 'W/m/K')
      Call info(messages, 2, .true.)

    Case (2) ! Drude
      Write (messages(1), '(a)') '  Electronic thermal conductivity: Drude model'
      Write (messages(2), '(a,1p g12.5e2)') '  -- T.C. at system temp. (W m^-1 K^-1): ', &
        convert_units(ttm%Ka0, 'k_b/ps/A', 'W/m/K')
      Call info(messages, 2, .true.)

    Case (3) ! Tabulated

      Call info('  Electronic thermal conductivity: Tabulated', .true.)

    End Select

    Select Case (ttm%detype)
    Case (0) ! Off
      Continue
    Case (1) ! Constant
      Write (messages(1), '(a)') '  Electronic thermal diffusivity: Constant'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Electronic T.D. (m^2 s^-1): ', convert_units(ttm%diff0, 'ang^2/ps', 'm^2/s')
      Call info(messages, 2, .true.)

    Case (2) ! Recip

      Write (messages(1), '(a)') '  Electronic thermal diffusivity: Reciprocal'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Datum electronic T.D. (m^2 s^-1): ', &
        convert_units(ttm%diff0, 'ang^2/ps', 'm^2/s') / thermo%temp
      Write (messages(3), '(a,1p g12.5e2)') '  -- Fermi temperature (K): ', ttm%Tfermi
      Call info(messages, 3, .true.)

    Case (3) ! Tabulated

      Call info('  Electronic thermal diffusivity: Tabulated', .true.)

    End Select

    Write (message, '(a,i8)') '  Mininum number of atoms for ionic cells: ', ttm%amin
    Call info(message, .true.)

    If (ttm%redistribute) Then
      Call info('  Energy redistribution: ON', .true.)
    End If

    Write (message, '(a,1p g12.5e2)') '  Elec. stopping power (eV/nm): ', convert_units(ttm%dedx, 'e.V/ang', 'e.V/nm')
    Call info(message, .true.)

    If (ttm%fluence < zero_plus) Then
      Select Case (ttm%sdepoType)
      Case (1)
        Write (messages(1), '(a)') '   Spatial energy deposition: Gaussian'
        Write (messages(2), '(a,1p g12.5e2)') '  -- Sigma of distribution (nm): ', convert_units(ttm%sig, 'ang', 'nm')
        Write (messages(3), '(a,1p g12.5e2)') '  -- Distribution cutoff (nm): ', convert_units(ttm%sigmax * ttm%sig, 'ang', 'nm')
        Call info(messages, 3, .true.)

      Case (2)
        Call info('  Spatial energy deposition: Homogeneous', .true.)
      End Select

    Else
      Select Case (ttm%sdepoType)
      Case (2)
        Write (messages(1), '(a)') '  Spatial energy deposition: Homogeneous Laser'
        Write (messages(2), '(a,1p g12.5e2)') &
          '  -- Absorbed fluence (mJ cm^-2): ', convert_units(ttm%fluence, 'e.V/ang^2', 'mJ/cm^2')
        Write (messages(3), '(a,1p g12.5e2)') '  -- Penetration depth (nm): ', convert_units(ttm%pdepth, 'ang', 'nm')
        Call info(messages, 3, .true.)

      Case (3)

        Write (messages(1), '(a)') '  Spatial energy deposition: Z-exponential decaying Laser '
        Write (messages(2), '(a,1p g12.5e2)') &
          '  -- Absorbed fluence at surface (mJ cm^-2): ', convert_units(ttm%fluence, 'e.V/ang^2', 'mJ/cm^2')
        Write (messages(3), '(a,1p g12.5e2)') '  -- Penetration depth (nm): ', convert_units(ttm%pdepth, 'ang', 'nm')
        Call info(messages, 3, .true.)
      End Select

    End If

    Select Case (ttm%tdepotype)
    Case (1)
      Write (messages(1), '(a)') '  Temporal energy deposition: Gaussian'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Sigma of distribution (ps): ', ttm%tdepo
      Write (messages(3), '(a,1p g12.5e2)') '  -- Distribution cutoff (ps): ', 2.0_wp * ttm%tcdepo * ttm%tdepo
      Call info(messages, 3, .true.)

    Case (2)
      Write (messages(1), '(a)') '  Temporal energy deposition: Decaying Exponential'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Tau of distribution (ps): ', ttm%tdepo
      Write (messages(3), '(a,1p g12.5e2)') '  -- Distribution cutoff (ps): ', ttm%tcdepo * ttm%tdepo
      Call info(messages, 3, .true.)

    Case (3)
      Call info('  Temporal energy deposition: Dirac delta', .true.)

    Case (4)
      Write (messages(1), '(a)') '  Temporal energy deposition: Square pulse'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Pulse duration (ps): ', ttm%tdepo
      Call info(messages, 2, .true.)

    End Select

    Select Case (ttm%gvar)
    Case (1)
      Write (messages(1), '(a)') '  Variable electron-phonon coupling: Homogeneous'
    Case (2)
      Write (messages(1), '(a)') '  Variable electron-phonon coupling: Heterogeneous'
    End Select
    Write (messages(2), '(a)') '    (overrides value given for ensemble, required tabulated stopping terms in g.dat file)'
    Call info(messages, 2, .true.)

    Select Case (ttm%bcTypeE)
    Case (1)
      Call info('  Electronic temperature boundary conditions: Periodic', .true.)

    Case (2)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Dirichlet'
      Write (messages(2), '(a)') '  -- Boundaries set to system temperature'
      Call info(messages, 2, .true.)

    Case (3)
      Call info('  Electronic temperature boundary conditions: Neumann', .true.)

    Case (4)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Dirichlet (XY), Neumann (Z)'
      Write (messages(2), '(a)') '  -- XY boundaries set to system temperature'
      Call info(messages, 2, .true.)

    Case (5)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Robin'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Temperature leakage (%): ', convert_units(ttm%fluxout, '', '%')
      Call info(messages, 2, .true.)

    Case (6)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Robin (XY), Neumann (Z):'
      Write (messages(2), '(a,1p g12.5e2)') '  -- XY Temperature leakage (%): ', convert_units(ttm%fluxout, '', '%')
      Call info(messages, 3, .true.)

    End Select

    Write (message, '(a,1p g12.5e2)') '  Electron-ion coupling offset (ps): ', ttm%ttmoffset
    Call info(message, .true.)

    If (ttm%oneway) Call info('  One-way electron-phonon coupling: ON', .true.)

    If (ttm%ttmstats > 0) Then
      Write (messages(1), '(a)') '  TTM statistics file: ON'
      Write (messages(2), '(a,i0.1)') '  -- TTM statistics file interval (steps): ', ttm%ttmstats
      Call info(messages, 2, .true.)
    End If

    If (ttm%ttmtraj > 0) Then
      Write (messages(1), '(a)') '  TTM trajectory (temperature profile) file: ON'
      Write (messages(2), '(a,i0.1)') '  -- TTM trajectory file interval (steps): ', ttm%ttmtraj
      Call info(messages, 2, .true.)
    End If

  End Subroutine write_ttm

  Subroutine initialise_control(table)
    !!-----------------------------------------------------------------------
    !!
    !! Initialise all control parameters to their default
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type(parameters_hash_table), Intent(  Out) :: table

    Call table%init(PARAMS_TABLE_SIZE)

    run_properties:block
      Call table%set('title', control_parameter( &
                     key='title', &
                     name='Title', &
                     val='', &
                     description="Run title", &
                     data_type=DATA_STRING))

      Call table%set('simulation_method', control_parameter( &
                     key="simulation_method", &
                     name="Simulation Method ", &
                     val="md", &
                     description="Set simulation method, options: MD, EVB, FFS", &
                     data_type=DATA_OPTION))

      Call table%set("random_seed", control_parameter( &
                     key="random_seed", &
                     name="Random seed", &
                     val="1 2 3", &
                     description="Set random seed", &
                     data_type=DATA_VECTOR3))

      Call table%set("density_variance", control_parameter( &
                     key="density_variance", &
                     name="Expected density variance", &
                     val="0.0", &
                     units="%", &
                     internal_units="", &
                     description="Set expected density variance for determining maximum array sizes", &
                     data_type=DATA_FLOAT))

      Call table%set("data_dump_frequency", control_parameter( &
                     key="data_dump_frequency", &
                     name="Data dumping", &
                     val="1000", &
                     units="steps", &
                     internal_units="steps", &
                     description="Set data dumping frequency for restarts", &
                     data_type=DATA_FLOAT))

      Call table%set("subcell_threshold", control_parameter( &
                     key="subcell_threshold", &
                     name="Subcelling threshold density", &
                     val="50.0", &
                     description="Set subcelling threshold density for setting minimum particles per link-cell", &
                     data_type=DATA_FLOAT))

      run_times:block
        Call table%set("time_run", control_parameter( &
                       key="time_run", &
                       name="Calculation run length", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Set duration of simulation (inc. equilibration)", &
                       data_type=DATA_FLOAT))

        Call table%set("time_equilibration", control_parameter( &
                       key="time_equilibration", &
                       name="Equilibration run length", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Set equilibration duration", &
                       data_type=DATA_FLOAT))

        Call table%set("time_job", control_parameter( &
                       key="time_job", &
                       name="Calculation job length", &
                       val="-1.0", &
                       units="hr", &
                       internal_units="s", &
                       description="Set total job time before attempted safe closure", &
                       data_type=DATA_FLOAT))

        Call table%set("time_close", control_parameter( &
                       key="time_close", &
                       name="Calculation close length", &
                       val="-1.0", &
                       units="min", &
                       internal_units="s", &
                       description="Estimated closure time for finite-time jobs", &
                       data_type=DATA_FLOAT))
      End block run_times

      statistics:block
        Call table%set("stats_frequency", control_parameter( &
                       key="stats_frequency", &
                       name="Stats Print Frequency", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Set frequency of stats sampling to statis file", &
                       data_type=DATA_FLOAT))

        Call table%set("stack_size", control_parameter( &
                       key="stack_size", &
                       name="Rolling average stack size", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Set rolling average stack to n timesteps", &
                       data_type=DATA_FLOAT))

        Call table%set("record_equilibration", control_parameter( &
                       key="record_equilibration", &
                       name="Record equilibration", &
                       val="off", &
                       description="Include equilibration in output", &
                       data_type=DATA_BOOL))

        Call table%set("print_probability_distribution", control_parameter( &
                       key="print_probability_distribution", &
                       name="Print probability distribution", &
                       val="off", &
                       description="Calculate and print probability distribution (enforces RDF print)", &
                       data_type=DATA_BOOL))

        analysis_options:block
          Call table%set("analyse_all", control_parameter( &
                         key="analyse_all", &
                         name="Analyse all", &
                         val="off", &
                         description="Enable analysis for all bonds, angles, dihedrals and inversions", &
                         data_type=DATA_BOOL))

          Call table%set("analyse_angles", control_parameter( &
                         key="analyse_angles", &
                         name="Analyse angles", &
                         val="off", &
                         description="Enable analysis for all angles", &
                         data_type=DATA_BOOL))

          Call table%set("analyse_bonds", control_parameter( &
                         key="analyse_bonds", &
                         name="Analyse bonds", &
                         val="off", &
                         description="Enable analysis for all bonds", &
                         data_type=DATA_BOOL))

          Call table%set("analyse_dihedrals", control_parameter( &
                         key="analyse_dihedrals", &
                         name="Analyse dihedrals", &
                         val="off", &
                         description="Enable analysis for all dihedrals", &
                         data_type=DATA_BOOL))

          Call table%set("analyse_inversions", control_parameter( &
                         key="analyse_inversions", &
                         name="Analyse inversions", &
                         val="off", &
                         description="Enable analysis for all inversions", &
                         data_type=DATA_BOOL))

          Call table%set("analyse_frequency", control_parameter( &
                         key="analyse_frequency", &
                         name="Analysis frequency", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set global frequency of data analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("analyse_frequency_bonds", control_parameter( &
                         key="analyse_frequency_bonds", &
                         name="Analysis frequency for bonds", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of bonds data analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("analyse_frequency_angles", control_parameter( &
                         key="analyse_frequency_angles", &
                         name="Analysis frequency for angles", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of angles data analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("analyse_frequency_dihedrals", control_parameter( &
                         key="analyse_frequency_dihedrals", &
                         name="Analysis frequency for dihedrals", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of dihedrals data analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("analyse_frequency_inversions", control_parameter( &
                         key="analyse_frequency_inversions", &
                         name="Analysis frequency for inversions", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of inversions data analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("analyse_max_dist", control_parameter( &
                         key="analyse_max_dist", &
                         name="Max analyse distance", &
                         val="2.0", &
                         units="ang", &
                         internal_units="internal_l", &
                         description="Set cutoff for bonds analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("analyse_num_bins", control_parameter( &
                         key="analyse_num_bins", &
                         name="Analysis number of bins", &
                         val="-1", &
                         description="Set global number of bins to be used in bonding analysis", &
                         data_type=DATA_INT))

          Call table%set("analyse_num_bins_bonds", control_parameter( &
                         key="analyse_num_bins_bonds", &
                         name="Analysis number of bins for bonds", &
                         val="0", &
                         description="Set number of bins to be used in bond analysis", &
                         data_type=DATA_INT))

          Call table%set("analyse_num_bins_angles", control_parameter( &
                         key="analyse_num_bins_angles", &
                         name="Analysis number of bins for angles", &
                         val="0", &
                         description="Set number of bins to be used in angle analysis", &
                         data_type=DATA_INT))

          Call table%set("analyse_num_bins_dihedrals", control_parameter( &
                         key="analyse_num_bins_dihedrals", &
                         name="Analysis number of bins for dihedrals", &
                         val="0", &
                         description="Set number of bins to be used in dihedral analysis", &
                         data_type=DATA_INT))

          Call table%set("analyse_num_bins_inversions", control_parameter( &
                         key="analyse_num_bins_inversions", &
                         name="Analysis number of bins for inversions", &
                         val="0", &
                         description="Set number of bins to be used in inversion analysis", &
                         data_type=DATA_INT))
        End block analysis_options

        msd:block
          Call table%set("msd_calculate", control_parameter( &
                         key="msd_calculate", &
                         name="MSD calculating", &
                         val="off", &
                         description="Enable calculation of MSD", &
                         data_type=DATA_BOOL))

          Call table%set("msd_print", control_parameter( &
                         key="msd_print", &
                         name="MSD printing", &
                         val="off", &
                         description="Enable printing of MSD", &
                         data_type=DATA_BOOL))

          Call table%set("msd_start", control_parameter( &
                         key="msd_start", &
                         name="MSD Start calculating", &
                         val="0", &
                         units="steps", &
                         internal_units="steps", &
                         description="Start timestep for dumping MSD configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("msd_frequency", control_parameter( &
                         key="msd_frequency", &
                         name="MSD calculation interval", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Interval between dumping MSD configurations", &
                         data_type=DATA_FLOAT))
        End block msd

        traj:block
          Call table%set("traj_calculate", control_parameter( &
                         key="traj_calculate", &
                         name="Trajectory calculating", &
                         val="off", &
                         description="Enable calculation of trajectory", &
                         data_type=DATA_BOOL))

          Call table%set("traj_key", control_parameter( &
                         key="traj_key", &
                         name="Trajectory key ", &
                         val="pos", &
                         description="Set trajectory output, options: pos, pos-vel, pos-vel-force, compressed", &
                         data_type=DATA_OPTION))

          Call table%set("traj_start", control_parameter( &
                         key="traj_start", &
                         name="trajectory Start calculating", &
                         val="0", &
                         units="steps", &
                         internal_units="steps", &
                         description="Start timestep for dumping trajectory configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("traj_interval", control_parameter( &
                         key="traj_interval", &
                         name="trajectory calculation interval", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Interval between dumping trajectory configurations", &
                         data_type=DATA_FLOAT))
        End block traj

        defects:block
          Call table%set("defects_calculate", control_parameter( &
                         key="defects_calculate", &
                         name="Defects calculating", &
                         val="off", &
                         description="Enable calculation of defects", &
                         data_type=DATA_BOOL))

          Call table%set("defects_start", control_parameter( &
                         key="defects_start", &
                         name="Defects Start calculating", &
                         val="0", &
                         units="steps", &
                         internal_units="steps", &
                         description="Start timestep for dumping defects configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("defects_interval", control_parameter( &
                         key="defects_interval", &
                         name="Defects calculation interval", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Interval between dumping defects configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("defects_distance", control_parameter( &
                         key="defects_distance", &
                         name="Defects distance condition ", &
                         val="0.75", &
                         units="ang", &
                         internal_units="internal_l", &
                         description="Set cutoff for deviation to be considered by defects as interstitial", &
                         data_type=DATA_FLOAT))

          Call table%set("defects_backup", control_parameter( &
                         key="defects_backup", &
                         name="Defects backup ", &
                         val="off", &
                         description="Enable defects backup", &
                         data_type=DATA_BOOL))

        End block defects

        displacements:block
          Call table%set("displacements_calculate", control_parameter( &
                         key="displacements_calculate", &
                         name="Displacements calculating", &
                         val="off", &
                         description="Enable calculation of displacements", &
                         data_type=DATA_BOOL))

          Call table%set("displacements_start", control_parameter( &
                         key="displacements_start", &
                         name="Displacements Start calculating", &
                         val="0", &
                         units="steps", &
                         internal_units="steps", &
                         description="Start timestep for dumping displacements configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("displacements_interval", control_parameter( &
                         key="displacements_interval", &
                         name="Displacements calculation interval", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Interval between dumping displacements configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("displacements_distance", control_parameter( &
                         key="displacements_distance", &
                         name="Displacements distance condition ", &
                         val="0.75", &
                         units="ang", &
                         internal_units="internal_l", &
                         description="Set cutoff for qualifying as displacement", &
                         data_type=DATA_FLOAT))

        End block displacements

        coord:block
          Call table%set("coord_calculate", control_parameter( &
                         key="coord_calculate", &
                         name="Coordination calculating", &
                         val="off", &
                         description="Enable calculation of coordination", &
                         data_type=DATA_BOOL))

          Call table%set("coord_ops", control_parameter( &
                         key="coord_ops", &
                         name="Coord ops", &
                         val="icoord", &
                         description="Set Coordops, options: icoord: only dumps the coordination of each atom at i; "// &
                         "CCOORD: only dumps coordination of each atom at i; "// &
                         "FULL: dumps the coordination of each atom every j steps", &
                         data_type=DATA_OPTION))

          Call table%set("coord_start", control_parameter( &
                         key="coord_start", &
                         name="Coordination start calculating", &
                         val="0", &
                         units="steps", &
                         internal_units="steps", &
                         description="Start timestep for dumping coordination configurations", &
                         data_type=DATA_FLOAT))

          Call table%set("coord_interval", control_parameter( &
                         key="coord_interval", &
                         name="Coordination calculation interval", &
                         val="100", &
                         units="steps", &
                         internal_units="steps", &
                         description="Interval between dumping coordination configurations", &
                         data_type=DATA_FLOAT))
        End block coord

        adf:block
          Call table%set("adf_calculate", control_parameter( &
                         key="adf_calculate", &
                         name="ADF calculating", &
                         val="off", &
                         description="Enable calculation of ADF", &
                         data_type=DATA_BOOL))

          Call table%set("adf_frequency", control_parameter( &
                         key="adf_frequency", &
                         name="ADF Sampling Frequency", &
                         val="100", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of ADF sampling", &
                         data_type=DATA_FLOAT))

          Call table%set("adf_precision", control_parameter( &
                         key="adf_precision", &
                         name="ADF Precision", &
                         val="0.0", &
                         description="Set precision of angular distribution bins in ADF analysis", &
                         data_type=DATA_FLOAT))
        End block adf

        rdf:block
          Call table%set("rdf_calculate", control_parameter( &
                         key="rdf_calculate", &
                         name="RDF calculating", &
                         val="off", &
                         description="Enable calculation of RDF", &
                         data_type=DATA_BOOL))

          Call table%set("rdf_print", control_parameter( &
                         key="rdf_print", &
                         name="RDF printing", &
                         val="on", &
                         description="Enable printing of RDF", &
                         data_type=DATA_BOOL))

          Call table%set("rdf_frequency", control_parameter( &
                         key="rdf_frequency", &
                         name="RDF Sampling Frequency", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of RDF sampling", &
                         data_type=DATA_FLOAT))

          Call table%set("rdf_binsize", control_parameter( &
                         key="rdf_binsize", &
                         name="RDF number of bins", &
                         val="0.05", &
                         units="ang", &
                         internal_units="ang", &
                         description="Set number of bins to be used in RDF analysis", &
                         data_type=DATA_FLOAT))

          Call table%set("rdf_error_analysis", control_parameter( &
                         key="rdf_error_analysis", &
                         name="RDF Error Analysis", &
                         val="off", &
                         description="Enable RDF error analysis, options: Off, Jackknife, Block", &
                         data_type=DATA_OPTION))

          Call table%set("rdf_error_analysis_blocks", control_parameter( &
                         key="rdf_error_analysis_blocks", &
                         name="Num RDF Error Analysis blocks", &
                         val="1", &
                         description="Set number of RDF error analysis blocks", &
                         data_type=DATA_INT))
        End block rdf

        zden:block
          Call table%set("zden_calculate", control_parameter( &
                         key="zden_calculate", &
                         name="ZDen calculating", &
                         val="off", &
                         description="Enable calculation of ZDen", &
                         data_type=DATA_BOOL))

          Call table%set("zden_print", control_parameter( &
                         key="zden_print", &
                         name="ZDen printing", &
                         val="on", &
                         description="Enable printing of ZDen", &
                         data_type=DATA_BOOL))

          Call table%set("zden_frequency", control_parameter( &
                         key="zden_frequency", &
                         name="ZDen Sampling Frequency", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of ZDen sampling", &
                         data_type=DATA_FLOAT))

          Call table%set("zden_binsize", control_parameter( &
                         key="zden_binsize", &
                         name="ZDen number of bins", &
                         val="0.05", &
                         units="ang", &
                         internal_units="ang", &
                         description="Set number of bins to be used in ZDen analysis", &
                         data_type=DATA_FLOAT))
        End block zden

        vaf:block
          Call table%set("vaf_calculate", control_parameter( &
                         key="vaf_calculate", &
                         name="VAF calculating", &
                         val="off", &
                         description="Enable calculation of VAF", &
                         data_type=DATA_BOOL))

          Call table%set("vaf_print", control_parameter( &
                         key="vaf_print", &
                         name="VAF printing", &
                         val="on", &
                         description="Enable printing of VAF", &
                         data_type=DATA_BOOL))

          Call table%set("vaf_frequency", control_parameter( &
                         key="vaf_frequency", &
                         name="VAF Sampling Frequency", &
                         val="1", &
                         units="steps", &
                         internal_units="steps", &
                         description="Set frequency of VAF sampling", &
                         data_type=DATA_FLOAT))

          Call table%set("vaf_binsize", control_parameter( &
                         key="vaf_binsize", &
                         name="VAF number of bins", &
                         val="0", &
                         description="Set number of bins to be used in VAF analysis", &
                         data_type=DATA_INT))

          Call table%set("vaf_averaging", control_parameter( &
                         key="vaf_averaging", &
                         name="VAF Time Averaging", &
                         val="on", &
                         description="Ignore time-averaging of VAF, "// &
                         "report all calculated VAF to VAFDAT files and final profile to OUTPUT", &
                         data_type=DATA_BOOL))
        End block vaf

        Call table%set("currents_calculate", control_parameter( &
                       key="currents_calculate", &
                       name="Calculate currents", &
                       val="off", &
                       description="Enable calculation of currents", &
                       data_type=DATA_BOOL))

      End block statistics

      Call table%set("heat_flux", control_parameter( &
                     key="heat_flux", &
                     name="Calculate heat flux", &
                     val="off", &
                     description="Enable calculation of heat flux", &
                     data_type=DATA_BOOL))

      Call table%set("write_per_particle", control_parameter( &
                     key="write_per_particle", &
                     name="Dump per-particle interaction information to file", &
                     val="off", &
                     description="Enable dumping of per-particle information", &
                     data_type=DATA_BOOL))
      per_particle:block

      End block per_particle

    End block run_properties

    io:block
      Call table%set("print_frequency", control_parameter( &
                     key="print_frequency", &
                     name="Results Print Frequency", &
                     val="0", &
                     units="steps", &
                     internal_units="steps", &
                     description="Set frequency of printing results to output", &
                     data_type=DATA_FLOAT))

      units:block
        Call table%set("io_units_scheme", control_parameter( &
                       key="io_units_scheme", &
                       name="I/O units scheme", &
                       val="internal", &
                       description="Set I/O units scheme, options: internal, si, atomic (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_length", control_parameter( &
                       key="io_units_length", &
                       name="I/O units length", &
                       val="internal_l", &
                       description="Set I/O units for length (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_time", control_parameter( &
                       key="io_units_time", &
                       name="I/O units time", &
                       val="internal_t", &
                       description="Set I/O units for time (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_mass", control_parameter( &
                       key="io_units_mass", &
                       name="I/O units mass", &
                       val="internal_m", &
                       description="Set I/O units for mass (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_charge", control_parameter( &
                       key="io_units_charge", &
                       name="I/O units charge", &
                       val="internal_q", &
                       description="Set I/O units for charge (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_energy", control_parameter( &
                       key="io_units_energy", &
                       name="I/O units energy", &
                       val="internal_e", &
                       description="Set I/O units for energy (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_pressure", control_parameter( &
                       key="io_units_pressure", &
                       name="I/O units pressure", &
                       val="internal_p", &
                       description="Set I/O units for pressure (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_force", control_parameter( &
                       key="io_units_force", &
                       name="I/O units force", &
                       val="internal_f", &
                       description="Set I/O units for force (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_velocity", control_parameter( &
                       key="io_units_velocity", &
                       name="I/O units velocity", &
                       val="internal_v", &
                       description="Set I/O units for velocity (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_power", control_parameter( &
                       key="io_units_power", &
                       name="I/O units power", &
                       val="internal_e/internal_t", &
                       description="Set I/O units for power (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_surface_tension", control_parameter( &
                       key="io_units_surface_tension", &
                       name="I/O units surface tension", &
                       val="internal_f/internal_l", &
                       description="Set I/O units for surface tension (*unused*)", &
                       data_type=DATA_OPTION))

        Call table%set("io_units_emf", control_parameter( &
                       key="io_units_emf", &
                       name="I/O units emf", &
                       val="internal_e/internal_q", &
                       description="Set I/O units for electromotive force (*unused*)", &
                       data_type=DATA_OPTION))

      End block units

      io_read:block
        Call table%set("io_read_method", control_parameter( &
                       key="io_read_method", &
                       name="I/O read method", &
                       val="mpiio", &
                       description="Set I/O read method, possible read methods: "// &
                       "mpiio, direct, master", &
                       data_type=DATA_OPTION))

        Call table%set("io_read_readers", control_parameter( &
                       key="io_read_readers", &
                       name="Num parallel I/O readers", &
                       val="0", &
                       units="(Automatic)", &
                       description="Set number of parallel I/O readers", &
                       data_type=DATA_INT))

        Call table%set("io_read_batch_size", control_parameter( &
                       key="io_read_batch_size", &
                       name="I/O reader batch size", &
                       val="0", &
                       units="(Automatic)", &
                       description="Set I/O reader batch size", &
                       data_type=DATA_INT))

        Call table%set("io_read_buffer_size", control_parameter( &
                       key="io_read_buffer_size", &
                       name="I/O reader buffer size", &
                       val="0", &
                       units="(Automatic)", &
                       description="Set I/O reader buffer size", &
                       data_type=DATA_INT))

        Call table%set("io_read_error_check", control_parameter( &
                       key="io_read_error_check", &
                       name="I/O error check on read", &
                       val="off", &
                       description="Enable extended error checking on read", &
                       data_type=DATA_BOOL))

        Call table%set('io_read_ascii_revold', control_parameter( &
                       key='io_read_ascii_revold', &
                       name='Read plain-text REVOLD', &
                       val='off', &
                       data_type=DATA_BOOL, &
                       description="Read human-readable (ASCII) REVOLD file"))

      End block io_read

      io_write:block
        Call table%set("io_write_method", control_parameter( &
                       key="io_write_method", &
                       name="I/O write method", &
                       val="mpiio", &
                       description="Set I/O write method, possible write methods: "// &
                       "mpiio, direct, master", &
                       data_type=DATA_OPTION))

        Call table%set("io_write_writers", control_parameter( &
                       key="io_write_writers", &
                       name="Num parallel I/O writers", &
                       val="0", &
                       units="(Automatic)", &
                       description="Set number of parallel I/O writers", &
                       data_type=DATA_INT))

        Call table%set("io_write_batch_size", control_parameter( &
                       key="io_write_batch_size", &
                       name="I/O writer batch size", &
                       val="0", &
                       units="(Automatic)", &
                       description="Set I/O writer batch size", &
                       data_type=DATA_INT))

        Call table%set("io_write_buffer_size", control_parameter( &
                       key="io_write_buffer_size", &
                       name="I/O writer buffer size", &
                       val="0", &
                       units="(Automatic)", &
                       description="Set I/O writer buffer size", &
                       data_type=DATA_INT))

        Call table%set("io_write_sorted", control_parameter( &
                       key="io_write_sorted", &
                       name="I/O sorted writing", &
                       val="on", &
                       description="Enable sorted output for atomic data", &
                       data_type=DATA_BOOL))

        Call table%set("io_write_error_check", control_parameter( &
                       key="io_write_error_check", &
                       name="I/O error check on write", &
                       val="off", &
                       description="Enable extended error checking on write", &
                       data_type=DATA_BOOL))

        Call table%set('io_write_ascii_revive', control_parameter( &
                       key='io_write_ascii_revive', &
                       name='Write plain-text REVIVE', &
                       val='off', &
                       data_type=DATA_BOOL, &
                       description="Write REVIVE as a human-readable (ASCII) file"))

      End block io_write

      io_file:block
        Call table%set("io_file_output", control_parameter( &
                       key="io_file_output", &
                       name="Output filepath", &
                       val="OUTPUT", &
                       description="Set output filepath, special options: SCREEN, NONE", &
                       data_type=DATA_STRING))

        Call table%set("io_file_config", control_parameter( &
                       key="io_file_config", &
                       name="Config filepath", &
                       val="CONFIG", &
                       description="Set input configuration filepath", &
                       data_type=DATA_STRING))

        Call table%set("io_file_field", control_parameter( &
                       key="io_file_field", &
                       name="Field filepath", &
                       val="FIELD", &
                       description="Set input field filepath", &
                       data_type=DATA_STRING))

        Call table%set("io_file_statis", control_parameter( &
                       key="io_file_statis", &
                       name="Statistics filepath", &
                       val="STATIS", &
                       description="Set output statistics filepath, special options: NONE", &
                       data_type=DATA_STRING))

        Call table%set("io_file_history", control_parameter( &
                       key="io_file_history", &
                       name="History filepath", &
                       val="HISTORY", &
                       description="Set output history filepath, special options: NONE", &
                       data_type=DATA_STRING))

        Call table%set("io_file_historf", control_parameter( &
                       key="io_file_historf", &
                       name="Historf filepath", &
                       val="HISTORF", &
                       description="Set output historf filepath, special options: NONE", &
                       data_type=DATA_STRING))

        Call table%set("io_file_revive", control_parameter( &
                       key="io_file_revive", &
                       name="Revive filepath", &
                       val="REVIVE", &
                       description="Set output revive filepath, special options: NONE", &
                       data_type=DATA_STRING))

        Call table%set("io_file_revold", control_parameter( &
                       key="io_file_revold", &
                       name="Revold filepath", &
                       val="REVOLD", &
                       description="Set output revold filepath, special options: NONE", &
                       data_type=DATA_STRING))

        Call table%set("io_file_revcon", control_parameter( &
                       key="io_file_revcon", &
                       name="Revcon filepath", &
                       val="REVCON", &
                       description="Set output revcon filepath, special options: NONE", &
                       data_type=DATA_STRING))

        Call table%set('io_file_rdf', control_parameter( &
                       key='io_file_rdf', &
                       name='rdf filepath', &
                       val='RDFDAT', &
                       description='Set output RDF filepath, special options: NONE', &
                       data_type=DATA_STRING))

        Call table%set('io_file_msd', control_parameter( &
                       key='io_file_msd', &
                       name='msd filepath', &
                       val='MSDTMP', &
                       description='Set output MSD filepath, special options: NONE', &
                       data_type=DATA_STRING))

        Call table%set('io_file_tabbnd', control_parameter( &
                       key='io_file_tabbnd', &
                       name='tabbnd filepath', &
                       val='TABBND', &
                       description='Set input TABBND filepath', &
                       data_type=DATA_STRING))

        Call table%set('io_file_tabang', control_parameter( &
                       key='io_file_tabang', &
                       name='tabang filepath', &
                       val='TABANG', &
                       description='Set input TABANG filepath', &
                       data_type=DATA_STRING))

        Call table%set('io_file_tabdih', control_parameter( &
                       key='io_file_tabdih', &
                       name='tabdih filepath', &
                       val='TABDIH', &
                       description='Set input TABDIH filepath', &
                       data_type=DATA_STRING))

        Call table%set('io_file_tabinv', control_parameter( &
                       key='io_file_tabinv', &
                       name='tabinv filepath', &
                       val='TABINV', &
                       description='Set input TABINV filepath', &
                       data_type=DATA_STRING))

        Call table%set('io_file_tabvdw', control_parameter( &
                       key='io_file_tabvdw', &
                       name='tabvdw filepath', &
                       val='TABVDW', &
                       description='Set input TABVDW filepath', &
                       data_type=DATA_STRING))

        Call table%set('io_file_tabeam', control_parameter( &
                       key='io_file_tabeam', &
                       name='tabeam filepath', &
                       val='TABEAM', &
                       description='Set input TABEAM filepath', &
                       data_type=DATA_STRING))

      End block io_file

      Call table%set('output_energy', control_parameter( &
                     key='output_energy', &
                     name='Output Final Energy', &
                     val='off', &
                     data_type=DATA_BOOL, &
                     description="Output final energy e_tot in output file"))

      Call table%set("ignore_config_indices", control_parameter( &
                     key="ignore_config_indices", &
                     name="Ignore Config Indices", &
                     val="off", &
                     description="Ignore indices as defined in CONFIG and use read order instead", &
                     data_type=DATA_BOOL))

      Call table%set("print_topology_info", control_parameter( &
                     key="print_topology_info", &
                     name="Print topology information ", &
                     val="off", &
                     description="Print topology information in output file", &
                     data_type=DATA_BOOL))

      Call table%set("print_level", control_parameter( &
                     key="print_level", &
                     name="Print level", &
                     val="1", &
                     description="Disable unnecessary printing, levels: 0 - silent, 1 - quiet, 2 - standard, 3 - full", &
                     data_type=DATA_INT))

      Call table%set("timer_depth", control_parameter( &
                     key="timer_depth", &
                     name="Timer print level", &
                     val="4", &
                     description="Do not display timers beyond this many levels in final timer output", &
                     data_type=DATA_INT))

      Call table%set("timer_per_mpi", control_parameter( &
                     key="timer_per_mpi", &
                     name="Per Process timing", &
                     val="off", &
                     description="Write timings for each MPI process individually", &
                     data_type=DATA_BOOL))

    End block io

    simulation_properties:block

      timestep:block
        Call table%set("timestep", control_parameter( &
                       key="timestep", &
                       name="Timestep", &
                       val="0.0", &
                       units="internal_t", &
                       internal_units="internal_t", &
                       description="Set calculation timestep or initial timestep for variable timestep calculations", &
                       data_type=DATA_FLOAT))

        Call table%set("timestep_variable", control_parameter( &
                       key="timestep_variable", &
                       name="Variable timestep", &
                       val="off", &
                       description="Enable variable timestep", &
                       data_type=DATA_BOOL))

        Call table%set("timestep_variable_min_dist", control_parameter( &
                       key="timestep_variable_min_dist", &
                       name="Variable timestep minimum distance", &
                       val="0.03", &
                       units="ang", &
                       internal_units="internal_l", &
                       description="Set minimum permissible distance for variable timestep", &
                       data_type=DATA_FLOAT))

        Call table%set("timestep_variable_max_dist", control_parameter( &
                       key="timestep_variable_max_dist", &
                       name="Variable timestep maximum distance", &
                       val="0.1", &
                       units="ang", &
                       internal_units="internal_l", &
                       description="Set maximum permissible distance for variable timestep", &
                       data_type=DATA_FLOAT))

        Call table%set("timestep_variable_max_delta", control_parameter( &
                       key="timestep_variable_max_delta", &
                       name="Variable timestep max delta", &
                       val="0.0", &
                       units="internal_t", &
                       internal_units="internal_t", &
                       description="Set maximum timestep delta for variable timestep", &
                       data_type=DATA_FLOAT))
      End block timestep

      ensemble:block
        Call table%set("ensemble", control_parameter( &
                       key="ensemble", &
                       name="Ensemble constraints", &
                       val="NVE", &
                       description="Set ensemble constraints, options: NVE, PMF, NVT, NPT, NST", &
                       data_type=DATA_OPTION))

        Call table%set("ensemble_method", control_parameter( &
                       key="ensemble_method", &
                       name="Ensemble method", &
                       val="", &
                       description="Set ensemble method, options "// &
                       "NVT: Evans, Langevin, Andersen, Berendsen, Hoover, gentle, ttm, dpds1, dpds2. "// &
                       "NP|ST: Langevin, Berendsen, Hoover, MTK.", &
                       data_type=DATA_OPTION))

        Call table%set("ensemble_thermostat_coupling", control_parameter( &
                       key="ensemble_thermostat_coupling", &
                       name="Thermostat coupling", &
                       val="0.0", &
                       units="ps", &
                       internal_units="ps", &
                       description="Set thermostat relaxation/decorrelation" &
                       //" times (use ensemble_thermostat_friction for langevin)", &
                       data_type=DATA_FLOAT))

        Call table%set("ensemble_dpd_order", control_parameter( &
                       key="ensemble_dpd_order", &
                       name="Ensemble DPD order", &
                       val="off", &
                       description="Set dpd method, options: off, first, second", &
                       data_type=DATA_OPTION))

        Call table%set("ensemble_dpd_drag", control_parameter( &
                       key="ensemble_dpd_drag", &
                       name="Ensemble DPD drag coefficient", &
                       val="0.0", &
                       units="Da/ps", &
                       internal_units="Da/ps", &
                       description="Set DPD drag coefficient", &
                       data_type=DATA_FLOAT))

        Call table%set("ensemble_thermostat_friction", control_parameter( &
                       key="ensemble_thermostat_friction", &
                       name="Thermostat Friction", &
                       val="0.0", &
                       units="ps^-1", &
                       internal_units="ps^-1", &
                       description="Set thermostat friction for langevin and gentle stochastic thermostats", &
                       data_type=DATA_FLOAT))

        Call table%set("ensemble_thermostat_softness", control_parameter( &
                       key="ensemble_thermostat_softness", &
                       name="Ensemble thermostat softness", &
                       val="0.0", &
                       units="", &
                       internal_units="", &
                       description="Set thermostat softness for Andersen thermostat", &
                       data_type=DATA_FLOAT))

        Call table%set("ensemble_barostat_coupling", control_parameter( &
                       key="ensemble_barostat_coupling", &
                       name="Barostat coupling", &
                       val="0.0", &
                       units="ps", &
                       internal_units="ps", &
                       description="Set barostat relaxation/decorrelation times (use ensemble_barostat_friction for langevin)", &
                       data_type=DATA_FLOAT))

        Call table%set("ensemble_barostat_friction", control_parameter( &
                       key="ensemble_barostat_friction", &
                       name="Barostat Friction", &
                       val="0.0", &
                       units="ps^-1", &
                       internal_units="ps^-1", &
                       description="Set barostat friction", &
                       data_type=DATA_FLOAT))

        Call table%set("ensemble_semi_isotropic", control_parameter( &
                       key="ensemble_semi_isotropic", &
                       name="Ensemble semi-isotropic constraint", &
                       val="off", &
                       description="Enable semi-isotropic barostat constraints, options: area, tension, orthorhombic", &
                       data_type=DATA_OPTION))

        Call table%set("ensemble_semi_orthorhombic", control_parameter( &
                       key="ensemble_semi_orthorhombic", &
                       name="Ensemble semi-orthorhombic constraint", &
                       val="off", &
                       description="Enable semi-orthorhombic barostat constraints", &
                       data_type=DATA_BOOL))

        Call table%set("ensemble_tension", control_parameter( &
                       key="ensemble_tension", &
                       name="Constrained surface tension", &
                       val="0.0", &
                       units="N/m", &
                       internal_units="dyn/cm", &
                       description="Set tension in NPngT calctulation", &
                       data_type=DATA_FLOAT))
      End block ensemble

      system_properties:block
        Call table%set("pressure_tensor", control_parameter( &
                       key="pressure_tensor", &
                       name="Pressure tensor", &
                       val="0.0 0.0 0.0 0.0 0.0 0.0", &
                       units="katm", &
                       internal_units="internal_p", &
                       description="Set the target pressure tensor for NsT calculations", &
                       data_type=DATA_VECTOR6))

        Call table%set("pressure_hydrostatic", control_parameter( &
                       key="pressure_hydrostatic", &
                       name="Hydrostatic pressure", &
                       val="0.0", &
                       units="katm", &
                       internal_units="internal_p", &
                       description="Set the target hydrostatic pressure (1/3Tr[P]) for NPT calculations", &
                       data_type=DATA_FLOAT))

        Call table%set("pressure_perpendicular", control_parameter( &
                       key="pressure_perpendicular", &
                       name="Perpendicular pressure", &
                       val="0.0 0.0 0.0", &
                       units="katm", &
                       internal_units="internal_p", &
                       description="Set the target pressure as x, y, z perpendicular to cell faces for NPT calculations", &
                       data_type=DATA_VECTOR3))

        Call table%set("temperature", control_parameter( &
                       key="temperature", &
                       name="Initial/Target temperature", &
                       val="0.0", &
                       units="K", &
                       internal_units="K", &
                       description="Set the initial temperature or target temperature (for thermostats)", &
                       data_type=DATA_FLOAT))
      End block system_properties

      pseudo_thermostat:block
        Call table%set("pseudo_thermostat_method", control_parameter( &
                       key="pseudo_thermostat_method", &
                       name="Pseudo thermostat method", &
                       val="off", &
                       description="Set pseudo thermostat method," &
                       //" possible options: Off, Langevin-Direct, Langevin, Gauss, Direct", &
                       data_type=DATA_OPTION))

        Call table%set("pseudo_thermostat_width", control_parameter( &
                       key="pseudo_thermostat_width", &
                       name="Pseudo thermostat width", &
                       val="2.0", &
                       units="ang", &
                       internal_units="internal_l", &
                       description="Set the width of thermostatted boundaries for pseudo thermostats", &
                       data_type=DATA_FLOAT))

        Call table%set("pseudo_thermostat_temperature", control_parameter( &
                       key="pseudo_thermostat_temperature", &
                       name="Pseudo thermostat temperature", &
                       val="0.0", &
                       units="K", &
                       internal_units="K", &
                       description="Set the temperature of the pseudo thermostat", &
                       data_type=DATA_FLOAT))
      End block pseudo_thermostat

      impact:block
        Call table%set("impact_part_index", control_parameter( &
                       key="impact_part_index", &
                       name="Impact particle index", &
                       val="0", &
                       description="Set particle index for impact simulations", &
                       data_type=DATA_INT))

        Call table%set("impact_time", control_parameter( &
                       key="impact_time", &
                       name="Impact time", &
                       val="0.0", &
                       units="internal_t", &
                       internal_units="steps", &
                       description="Set time for impact in impact simulations", &
                       data_type=DATA_FLOAT))

        Call table%set("impact_energy", control_parameter( &
                       key="impact_energy", &
                       name="Impact energy", &
                       val="0.0", &
                       units="ke.V", &
                       internal_units="ke.V", &
                       description="Set impact energy for impact simulations", &
                       data_type=DATA_FLOAT))

        Call table%set("impact_direction", control_parameter( &
                       key="impact_direction", &
                       name="Impact direction", &
                       val="1.0 1.0 1.0", &
                       units="", &
                       internal_units="", &
                       description="Direction vector for impact simulations", &
                       data_type=DATA_VECTOR3))
      End block impact

      ttm:block

        Call table%set("ttm_calculate", control_parameter( &
                       key="ttm_calculate", &
                       name="TTM Active", &
                       val="off", &
                       description="Enable calculation of two-temperature model", &
                       data_type=DATA_BOOL))

        Call table%set("ttm_num_ion_cells", control_parameter( &
                       key="ttm_num_ion_cells", &
                       name="Number of TTM ionic cells", &
                       val="10", &
                       description="Set number of coarse-grained ion temperature cells (CIT)", &
                       data_type=DATA_INT))

        Call table%set("ttm_num_elec_cells", control_parameter( &
                       key="ttm_num_elec_cells", &
                       name="Number of TTM electronic cells", &
                       val="50 50 50", &
                       description="Set number of coarse-grained electronic temperature cells (CET)", &
                       data_type=DATA_VECTOR3))

        Call table%set("ttm_metal", control_parameter( &
                       key="ttm_metal", &
                       name="TTM Metallic", &
                       val="off", &
                       description="Specifies parameters for metallic system are required for two-temperature model"// &
                       ", i.e. thermal conductivity", &
                       data_type=DATA_BOOL))

        Call table%set("ttm_heat_cap_model", control_parameter( &
                       key="ttm_heat_cap_model", &
                       name="TTM Heat Capacity model", &
                       val="", &
                       description="Sets model for specific heat capacity in TTM, options: const, linear, tabulated, tanh", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_heat_cap", control_parameter( &
                       key="ttm_heat_cap", &
                       name="TTM Heat Capacity", &
                       val="0.0", &
                       units="internal_e/internal_m/K", &
                       internal_units="internal_e/internal_m/K", &
                       description="Sets constant, scale or maximum heat capcity in TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_temp_term", control_parameter( &
                       key="ttm_temp_term", &
                       name="TTM Tanh temperature term", &
                       val="0.0", &
                       units="K^-1", &
                       internal_units="K^-1", &
                       description="Set Fermi temperature in TTM, for tanh", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_fermi_temp", control_parameter( &
                       key="ttm_fermi_temp", &
                       name="TTM Fermi temperature", &
                       val="0.0", &
                       units="K", &
                       internal_units="K", &
                       description="Set Fermi temperature in TTM, for linear", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_elec_cond_model", control_parameter( &
                       key="ttm_elec_cond_model", &
                       name="TTM Electonic conductivity model", &
                       val="", &
                       description="Set electronic conductivity model in TTM, options: Infinite, constant, drude, tabulated", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_elec_cond", control_parameter( &
                       key="ttm_elec_cond", &
                       name="TTM Electronic conductivity", &
                       val="0.0", &
                       units="W/m/K", &
                       internal_units="k_b/ps/Ang", &
                       description="Set electronic conductivity in TTM ", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_diff_model", control_parameter( &
                       key="ttm_diff_model", &
                       name="TTM Diffusion model", &
                       val="", &
                       description="Set diffusion model in TTM, options: constant, recip, tabulated", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_diff", control_parameter( &
                       key="ttm_diff", &
                       name="TTM Thermal diffusivity", &
                       val="0.0", &
                       units="m^2/s", &
                       internal_units="ang^2/ps", &
                       description="Set TTM thermal diffusivity", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_dens_model", control_parameter( &
                       key="ttm_dens_model", &
                       name="TTM Density model", &
                       val="", &
                       description="Set density model in TTM, options are: constant, dynamic", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_dens", control_parameter( &
                       key="ttm_dens", &
                       name="TTM Density", &
                       val="0.0", &
                       units="ang^-3", &
                       internal_units="internal_l^-3", &
                       description="Set constant density in TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_min_atoms", control_parameter( &
                       key="ttm_min_atoms", &
                       name="TTM minimum atoms", &
                       val="0", &
                       description="Minimum number of atoms needed per ionic temperature cell", &
                       data_type=DATA_INT))

        Call table%set("ttm_stopping_power", control_parameter( &
                       key="ttm_stopping_power", &
                       name="TTM Electron stopping power", &
                       val="0.0", &
                       units="e.V/nm", &
                       internal_units="e.V/ang", &
                       description="Electronic stopping power of projectile entering electronic system ", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_spatial_dist", control_parameter( &
                       key="ttm_spatial_dist", &
                       name="TTM Spatial deposition distribution", &
                       val="", &
                       description="Set the spatial distribution of TTM, options: flat, gaussian, flat-laser, exp-laser", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_spatial_sigma", control_parameter( &
                       key="ttm_spatial_sigma", &
                       name="TTM Spatial sigma", &
                       val="1.0", &
                       units="nm", &
                       internal_units="internal_l", &
                       description="Set the sigma for spatial distributions of TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_spatial_cutoff", control_parameter( &
                       key="ttm_spatial_cutoff", &
                       name="TTM Spatial cutoff", &
                       val="5.0", &
                       units="nm", &
                       internal_units="nm", &
                       description="Set the cutoff for spatial distributions of TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_fluence", control_parameter( &
                       key="ttm_fluence", &
                       name="TTM laser absorbed energy", &
                       val="0.0", &
                       units="mJ/cm^2", &
                       internal_units="e.V/ang^2", &
                       description="Initial energy deposition into electronic system by laser for TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_penetration_depth", control_parameter( &
                       key="ttm_penetration_depth", &
                       name="TTM laser penetration depth", &
                       val="0.0", &
                       units="nm", &
                       internal_units="ang", &
                       description="Set laser penetration depth for TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_laser_type", control_parameter( &
                       key="ttm_laser_type", &
                       name="TTM Deposition type", &
                       val="flat", &
                       description="Set laser deposition type. options: flat, exponential", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_temporal_dist", control_parameter( &
                       key="ttm_temporal_dist", &
                       name="TTM Temporal distribution", &
                       val="", &
                       description="Set temporal distribution for TTM, options: gaussian, exponential, delta, square", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_temporal_duration", control_parameter( &
                       key="ttm_temporal_duration", &
                       name="TTM Duration", &
                       val="0.001", &
                       units="ps", &
                       internal_units="ps", &
                       description="Set duration of energy deposition for TTM (gaussian, exponential, square)", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_temporal_cutoff", control_parameter( &
                       key="ttm_temporal_cutoff", &
                       name="TTM Time Cutoff", &
                       val="5.0", &
                       units="ps", &
                       internal_units="ps", &
                       description="Set temporal cutoff for TTM (gaussian, exponential)", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_variable_ep", control_parameter( &
                       key="ttm_variable_ep", &
                       name="TTM Variable electron-phonon coupling", &
                       val="", &
                       description="Set electron-phonon coupling for TTM, options: homo, hetero", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_boundary_condition", control_parameter( &
                       key="ttm_boundary_condition", &
                       name="TTM Boundary conditions", &
                       val="", &
                       description="Set boundary conditions for TTM, options: periodic, dirichlet, neumann, robin", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_boundary_xy", control_parameter( &
                       key="ttm_boundary_xy", &
                       name="TTM Neumann Boundary in Z", &
                       val="off", &
                       description="Fix Neumann (zero-flux) boundary in Z", &
                       data_type=DATA_BOOL))

        Call table%set("ttm_boundary_heat_flux", control_parameter( &
                       key="ttm_boundary_heat_flux", &
                       name="TTM Heat flux", &
                       val="96", &
                       units="%", &
                       internal_units="", &
                       description="Set boundary heat flux in Robin boundaries for TTM", &
                       data_type=DATA_BOOL))

        Call table%set("ttm_time_offset", control_parameter( &
                       key="ttm_time_offset", &
                       name="TTM time offset", &
                       val="0.0", &
                       units="ps", &
                       internal_units="ps", &
                       description="Set electron-ion coupling offset for TTM", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_oneway", control_parameter( &
                       key="ttm_oneway", &
                       name="TTM oneway", &
                       val="off", &
                       description="Enable one-way electron-phonon coupling "// &
                       "when electronic temperature is greater than ionic temperature", &
                       data_type=DATA_BOOL))

        Call table%set("ttm_stats_frequency", control_parameter( &
                       key="ttm_stats_frequency", &
                       name="Frequency of TTM statis output", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Frequency of write to TTM PEAK_E and PEAK_I", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_traj_frequency", control_parameter( &
                       key="ttm_traj_frequency", &
                       name="Frequency of TTM trajectory output", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Frequency of write to TTM LATS_E and LATS_I", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_com_correction", control_parameter( &
                       key="ttm_com_correction", &
                       name="TTM CoM correction", &
                       val="full", &
                       description="Apply inhomogeneous Langevin thermostat to" &
                       //" selected directions in TTM, options: full, zdir, off", &
                       data_type=DATA_OPTION))

        Call table%set("ttm_redistribute", control_parameter( &
                       key="ttm_redistribute", &
                       name="", &
                       val="", &
                       description="Redistribute electronic energy in newly-deactivated" &
                       //" temperature cells to nearest active neighbours", &
                       data_type=DATA_BOOL))

        Call table%set("ttm_e-phonon_friction", control_parameter( &
                       key="ttm_e-phonon_friction", &
                       name="TTM Electron-phonon friction", &
                       val="0.0", &
                       units="ps^-1", &
                       internal_units="ps^-1", &
                       description="Set TTM electron-phonon friction", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_e-stopping_friction", control_parameter( &
                       key="ttm_e-stopping_friction", &
                       name="TTM Electron-stopping friction", &
                       val="0.0", &
                       units="ps^-1", &
                       internal_units="ps^-1", &
                       description="Set TTM electron-stopping friction", &
                       data_type=DATA_FLOAT))

        Call table%set("ttm_e-stopping_velocity", control_parameter( &
                       key="ttm_e-stopping_velocity", &
                       name="TTM Electron-stopping velocity", &
                       val="0.0", &
                       units="ang/ps", &
                       internal_units="internal_v", &
                       description="Set TTM electron-stopping velocity", &
                       data_type=DATA_FLOAT))

      End block ttm

      integrator_tolerances:block

        Call table%set("rlx_cgm_step", control_parameter( &
                       key="rlx_cgm_step", &
                       name="Relaxed shell CGM step", &
                       val="-1.0", &
                       units="ang", &
                       internal_units="internal_l", &
                       description="Set CGM stepping for relaxed shell model", &
                       data_type=DATA_FLOAT))

        Call table%set("rlx_tol", control_parameter( &
                       key="rlx_tol", &
                       name="Relaxed shell force tolerance", &
                       val="1.0", &
                       units="internal_f", &
                       internal_units="internal_f", &
                       description="Set force tolerance for relaxed shell model", &
                       data_type=DATA_FLOAT))

        Call table%set("shake_max_iter", control_parameter( &
                       key="shake_max_iter", &
                       name="Max SHAKE/RATTLE iterations", &
                       val="250", &
                       description="Set maximum number of SHAKE/RATTLE iterations", &
                       data_type=DATA_INT))

        Call table%set("shake_tolerance", control_parameter( &
                       key="shake_tolerance", &
                       name="SHAKE/RATTLE tolerance", &
                       val="1e-6", &
                       units="ang", &
                       internal_units="internal_l", &
                       description="Set accepted SHAKE/RATTLE tolerance", &
                       data_type=DATA_FLOAT))

        ! call table%set("quaternion_max_iter", control_parameter( &
        !      key = "quaternion_max_iter", &
        !      name = "Max Leapfrog Verlet FIQA iterations", &
        !      val = "100", &
        !      description = "Set max iterations for FIQA in Leapfrog Verlet", &
        !      data_type = DATA_INT))

        ! call table%set("quaternion_tol", control_parameter( &
        !      key = "quaternion_tol", &
        !      name = "FIQA tolerance", &
        !      val = "1e-8", &
        !      units = "", &
        !      internal_units = "", &
        !      description = "Set accepted FIQA tolerance in Leapfrog Verlet", &
        !      data_type = DATA_FLOAT))
      End block integrator_tolerances

      Call table%set("dftb", control_parameter( &
                     key="dftb", &
                     name="Enable DFTB", &
                     val="off", &
                     description="Enable DFTB", &
                     data_type=DATA_BOOL))

      Call table%set("fixed_com", control_parameter( &
                     key="fixed_com", &
                     name="Fix centre of mass", &
                     val="on", &
                     description="Remove net centre of mass momentum", &
                     data_type=DATA_BOOL))

    End block simulation_properties

    equilibration_properties:block
      Call table%set("reset_temperature_interval", control_parameter( &
                     key="reset_temperature_interval", &
                     name="Reset temperature time", &
                     val="-1", &
                     units="steps", &
                     internal_units="steps", &
                     description="Interval between temperature zeroing during equilibration for minimisation", &
                     data_type=DATA_FLOAT))

      Call table%set("regauss_frequency", control_parameter( &
                     key="regauss_frequency", &
                     name="Regauss frequency", &
                     val="-1", &
                     units="steps", &
                     internal_units="steps", &
                     description="Set the frequency of temperature regaussing", &
                     data_type=DATA_FLOAT))

      Call table%set("rescale_frequency", control_parameter( &
                     key="rescale_frequency", &
                     name="Rescale frequency", &
                     val="-1", &
                     units="steps", &
                     internal_units="steps", &
                     description="Set the frequency of temperature rescaling", &
                     data_type=DATA_FLOAT))

      Call table%set("equilibration_force_cap", control_parameter( &
                     key="equilibration_force_cap", &
                     name="Equilibration force cap", &
                     val="1000.0", &
                     units="k_b.temp/ang", &
                     internal_units="k_b.temp/ang", &
                     description="Set force cap clamping maximum force during equilibration", &
                     data_type=DATA_FLOAT))

      minimisation:block
        Call table%set("minimisation_criterion", control_parameter( &
                       key="minimisation_criterion", &
                       name="Minimisation criterion", &
                       val="off", &
                       description="Set minimisation criterion, options: off, force, energy, distance", &
                       data_type=DATA_OPTION))

        Call table%set("minimisation_tolerance", control_parameter( &
                       key="minimisation_tolerance", &
                       name="Minimisation tolerance", &
                       val="0.0", &
                       description="Set minimisation tolerance, units: determined by criterion", &
                       data_type=DATA_FLOAT))

        Call table%set("minimisation_step_length", control_parameter( &
                       key="minimisation_step_length", &
                       name="Minimisation step length", &
                       val="-1.0", &
                       units="ang", &
                       internal_units="internal_l", &
                       description="Set minimisation tolerance", &
                       data_type=DATA_FLOAT))

        Call table%set("minimisation_frequency", control_parameter( &
                       key="minimisation_frequency", &
                       name="Minimisation frequency", &
                       val="0", &
                       units="steps", &
                       internal_units="steps", &
                       description="Set minimisation frequency", &
                       data_type=DATA_FLOAT))

      End block minimisation

    End block equilibration_properties

    initialisation_parameters:block
      Call table%set('initial_minimum_separation', control_parameter( &
                     key='initial_minimum_separation', &
                     name='Initial minimum separation', &
                     val='-1.0', &
                     units='internal_l', &
                     internal_units='internal_l', &
                     data_type=DATA_FLOAT, &
                     description='Turn on the check on minimum separation distance between VNL pairs at re/start'))

      Call table%set("restart", control_parameter( &
                     key="restart", &
                     name="Restart", &
                     val="clean", &
                     description="Set restart settings, possible options: Clean, Continue, Rescale, Noscale", &
                     data_type=DATA_OPTION))

      Call table%set("nfold", control_parameter( &
                     key="nfold", &
                     name="N Fold", &
                     val="1 1 1", &
                     description="Expand cell before running", &
                     data_type=DATA_VECTOR3))

    End block initialisation_parameters

    forcefields:block

      global:block
        Call table%set("cutoff", control_parameter( &
                       key="cutoff", &
                       name="Real-space global cutoff", &
                       val="1.0", &
                       units="internal_l", &
                       internal_units="internal_l", &
                       description="Set the global cutoff for real-speace potentials", &
                       data_type=DATA_FLOAT))

        Call table%set("padding", control_parameter( &
                       key="padding", &
                       name="Set VNL padding", &
                       val="0.0", &
                       units="internal_l", &
                       internal_units="internal_l", &
                       description="Set padding for sizing of Verlet neighbour lists", &
                       data_type=DATA_FLOAT))
      End block global

      coulomb:block
        Call table%set("coul_damping", control_parameter( &
                       key="coul_damping", &
                       name="Electrostatics Fennell damping", &
                       val="0.0", &
                       units="1/Ang", &
                       internal_units="1/internal_l", &
                       description="Calculate electrostatics using Fennell damping (Ewald-like) with given alpha", &
                       data_type=DATA_FLOAT))

        Call table%set("coul_dielectric_constant", control_parameter( &
                       key="coul_dielectric_constant", &
                       name="Dielectric Constant", &
                       val="1.0", &
                       description="Set dielectric constant relative to vacuum", &
                       data_type=DATA_FLOAT))

        Call table%set("coul_extended_exclusion", control_parameter( &
                       key="coul_extended_exclusion", &
                       name="Extended Exclusion", &
                       val="off", &
                       description="Enable extended coulombic exclusion affecting intra-molecular interactions", &
                       data_type=DATA_BOOL))

        Call table%set("coul_method", control_parameter( &
                       key="coul_method", &
                       name="Electrostatics method", &
                       val="off", &
                       description="Set method for electrostatics method, "// &
                       "options: off, spme, dddp, pairwise, reaction_field, force_shifted", &
                       data_type=DATA_OPTION))

        Call table%set("coul_precision", control_parameter( &
                       key="coul_precision", &
                       name="Electrostatics Fennell precision", &
                       val="0.0", &
                       description="Calculate electrostatics using Fennell damping (Ewald-like) with given precision", &
                       data_type=DATA_FLOAT))

        ewald:block
          Call table%set("ewald_precision", control_parameter( &
                         key="ewald_precision", &
                         name="Ewald precision", &
                         val="1.0e-6", &
                         description="Set Ewald parameters to calculate within given precision for Ewald calculations", &
                         data_type=DATA_FLOAT))

          Call table%set("ewald_alpha", control_parameter( &
                         key="ewald_alpha", &
                         name="Ewald alpha", &
                         val="0.0", &
                         units="ang^-1", &
                         internal_units="internal_l^-1", &
                         description="Set real-recip changeover location for Ewald calculations", &
                         data_type=DATA_FLOAT))

          Call table%set("ewald_kvec", control_parameter( &
                         key="ewald_kvec", &
                         name="Ewald k-vector samples", &
                         val="0 0 0", &
                         units="", &
                         internal_units="", &
                         description="Set number of k-space samples for Ewald calculations", &
                         data_type=DATA_VECTOR3))

          Call table%set("ewald_kvec_spacing", control_parameter( &
                         key="ewald_kvec_spacing", &
                         name="Ewald k-vector spacing", &
                         val="0.0", &
                         units="ang^-1", &
                         internal_units="internal_l^-1", &
                         description="Calculate k-vector samples for an even sampling of given spacing in Ewald calculations", &
                         data_type=DATA_FLOAT))

          Call table%set("ewald_nsplines", control_parameter( &
                         key="ewald_nsplines", &
                         name="Number of B-Splines", &
                         val="8", &
                         description="Set number of B-Splines for Ewald SPME calculations, min=3", &
                         data_type=DATA_INT))
        End block ewald

      End block coulomb

      polarisation:block
        Call table%set("polarisation_model", control_parameter( &
                       key="polarisation_model", &
                       name="Polarisation model", &
                       val="default", &
                       description="Enable polarisation, options: default, CHARMM", &
                       data_type=DATA_OPTION))

        Call table%set("polarisation_thole", control_parameter( &
                       key="polarisation_thole", &
                       name="Polarisation", &
                       val="1.3", &
                       description="Set global atomic damping factor", &
                       data_type=DATA_FLOAT))
      End block polarisation

      metal:block
        Call table%set('metal_direct', control_parameter( &
                       key='metal_direct', &
                       name='Direct metallic force calculation', &
                       val='off', &
                       description="Enable direct (non-tabulated) calculation of metallic forces", &
                       data_type=DATA_BOOL))

        Call table%set("metal_sqrtrho", control_parameter( &
                       key="metal_sqrtrho", &
                       name="SqrtRho metallic interpolation", &
                       val="off", &
                       description="Enable metal sqrtrho interpolation option for EAM embeding function in TABEAM", &
                       data_type=DATA_BOOL))
      End block metal

      vdw:block
        Call table%set("vdw_method", control_parameter( &
                       key="vdw_method", &
                       name="VdW Method", &
                       val="tabulated", &
                       description="Set method for Van der Waal's calculations, options: off, direct, tabulated, ewald", &
                       data_type=DATA_OPTION))

        Call table%set("vdw_cutoff", control_parameter( &
                       key="vdw_cutoff", &
                       name="VdW cutoff", &
                       val="-1.0", &
                       units="internal_l", &
                       internal_units="internal_l", &
                       description="Set cut-off for Van der Waal's potentials", &
                       data_type=DATA_FLOAT))

        Call table%set('vdw_mix_method', control_parameter( &
                       key='vdw_mix_method', &
                       name='VdW mixing method', &
                       val='off', &
                       description="Enable VdW mixing, possible mixing schemes: Off, "// &
                       "Lorentz-Berthelot, Fender-Hasley, Hogervorst, Waldman-Hagler, Tang-Toennies, Functional", &
                       data_type=DATA_OPTION))

        Call table%set('vdw_force_shift', control_parameter( &
                       key='vdw_force_shift', &
                       name='VdW force shifting', &
                       val='off', &
                       description="Enable force shift corrections to Van der Waals' forces", &
                       data_type=DATA_BOOL))

      End block vdw

    End block forcefields

    plumed:block
      Call table%set("plumed", control_parameter( &
                     key="plumed", &
                     name="Enable Plumed", &
                     val="off", &
                     description="Enabled plumed dynamics", &
                     data_type=DATA_BOOL))

      Call table%set("plumed_input", control_parameter( &
                     key="plumed_input", &
                     name="Plumed input", &
                     val="", &
                     description="Set plumed input file", &
                     data_type=DATA_STRING))

      Call table%set("plumed_log", control_parameter( &
                     key="plumed_log", &
                     name="Plumed log", &
                     val="", &
                     description="Set plumed log file", &
                     data_type=DATA_STRING))

      Call table%set("plumed_precision", control_parameter( &
                     key="plumed_precision", &
                     name="Plumed precision", &
                     val="8", &
                     description="Set plumed numerical precision (4=single, 8=double)", &
                     data_type=DATA_INT))

      Call table%set("plumed_restart", control_parameter( &
                     key="plumed_restart", &
                     name="Is Plumed restart", &
                     val="on", &
                     description="Restart plumed dynamics", &
                     data_type=DATA_BOOL))

    End block plumed

    miscellaneous:block

      Call table%set("strict_checks", control_parameter( &
                     key="strict_checks", &
                     name="Enable strict", &
                     val="on", &
                     description="Enforce strict checks such as: "// &
                     "good system cutoff, particle index contiguity, disable non-error warnings, minimisation information", &
                     data_type=DATA_BOOL))

      Call table%set("unsafe_comms", control_parameter( &
                     key="unsafe_comms", &
                     name="Disable parallel safety", &
                     val="off", &
                     description="Do not ensure checks of logicals are enforced in parallel", &
                     data_type=DATA_BOOL))

      Call table%set("dftb_test", control_parameter( &
                     key="dftb_test", &
                     name="Run dftb tests", &
                     val="off", &
                     description="Do not perform a DLPOLY run, instead run dftb tests", &
                     data_type=DATA_BOOL))

    End block miscellaneous

    replay:block

      Call table%set("replay", control_parameter( &
                     key="replay", &
                     name="Replay history", &
                     val="off", &
                     description="Don't perform simulation, instead replay history", &
                     data_type=DATA_BOOL))

      Call table%set("replay_calculate_forces", control_parameter( &
                     key="replay_calculate_forces", &
                     name="Replay history recalculate forces", &
                     val="on", &
                     description="On history replay recalculate forces", &
                     data_type=DATA_BOOL))

    End block replay

  End Subroutine initialise_control

  Function try_parse(ifile, comm) Result(can_parse)
    !!-----------------------------------------------------------------------
    !!
    !! Attempt to detect if a control file is new or old style
    !! Only based on first keyword (usually title)
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins may 2020
    !!-----------------------------------------------------------------------
    Integer,          Intent(In   ) :: ifile
    Type(comms_type), Intent(InOut) :: comm
    Logical                         :: can_parse

    Character(Len=STR_LEN) :: input, key, units, val
    Logical                :: line_read

    Do
      key = ''; val = ''; units = ''
      Call get_line(line_read, ifile, input, comm)
      If (.not. line_read) Exit
      ! Trim comments
      If (Scan(input, "#!") > 0) input = input(:Scan(input, "#!") - 1)
      Call get_word(input, key)
      ! Skip blanks
      If (key == "") Cycle
      Call lower_case(key)
      can_parse = key == "title"
      ! can_parse = params%in(key)
      Exit
    End Do
    Rewind (ifile)

  End Function try_parse

  Subroutine parse_control_file(ifile, params, comm)
    !!-----------------------------------------------------------------------
    !!
    !! Read a control file filling the params table
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------

    Integer,                     Intent(In   ) :: ifile
    Type(parameters_hash_table), Intent(InOut) :: params
    Type(comms_type),            Intent(InOut) :: comm

    Character(Len=STR_LEN)  :: input, key, units, val
    Logical                 :: line_read
    Type(control_parameter) :: param

    Do
      key = ''; val = ''; units = ''
      Call get_line(line_read, ifile, input, comm)
      If (.not. line_read) Exit
      ! Trim comments
      If (Scan(input, "#!") > 0) input = input(:Scan(input, "#!") - 1)
      Call get_word(input, key)
      ! Skip blanks
      If (key == "") Cycle
      Call lower_case(key)
      If (.not. params%in(key)) Call error(0, 'Unrecognised key '//Trim(key))
      Call params%get(key, param)
      If (param%set) Call error(0, 'Param '//Trim(key)//' already set')
      Call read_control_param(input, param, ifile, comm)
      param%set = .true.
      Call params%set(key, param)
    End Do

    Call params%fix()
  End Subroutine parse_control_file

  Subroutine read_control_param(input, param, ifile, comm)
    Character(Len=*),        Intent(InOut) :: input
    Type(control_parameter), Intent(InOut) :: param
    Integer,                 Intent(In   ) :: ifile
    Type(comms_type),        Intent(InOut) :: comm

    Character(Len=MAX_KEY) :: unit
    Character(Len=STR_LEN) :: tmp, val
    Integer                :: i, test_int
    Logical                :: line_read
    Real(kind=wp)          :: test_real

    Select Case (param%data_type)
    Case (DATA_INT)
      Call get_word(input, val)
      test_int = Nint(word_2_real(val))
      param%val = val
    Case (DATA_FLOAT)
      Call get_word(input, val)
      test_real = word_2_real(val)
      param%val = val
      Call get_word(input, unit)
      param%units = unit
    Case (DATA_VECTOR3, DATA_VECTOR6)
      ! Handle Multiline
      i = Index(input, '&')
      tmp = ''
      Do While (i > 0)
        tmp = Trim(tmp)//input(:i - 1)
        If (input(i + 1:) /= '') Call error(0, 'Unexpected junk '//Trim(input)//' after key '//Trim(param%key))
        Call get_line(line_read, ifile, input, comm)
        i = Index(input, '&')
      End Do
      tmp = Adjustl(Trim(tmp)//input)

      If (tmp(1:1) /= '[') Call error(0, '')
      i = Index(tmp, ']', back=.true.)

      ! Cut off braces
      input = tmp(i + 1:)
      tmp = tmp(2:i - 1)

      param%val = ""
      Do While (tmp /= '')
        Call get_word(tmp, val)
        ! Check all valid reals
        test_real = word_2_real(val)
        param%val = Trim(param%val)//' '//val
      End Do

      Call get_word(input, unit)
      param%units = unit

    Case (DATA_OPTION, DATA_BOOL)
      Call get_word(input, val)
      If (val == "") Call error(0, 'Missing option for '//Trim(param%key))
      Call lower_case(val)
      param%val = val
      input = ""
    Case (DATA_STRING)
      param%val = Adjustl(input)
      input = ""
    Case Default
      Call error(0, 'Unknown data type while parsing '//Trim(param%key))
    End Select

    If (input /= '') Call error(0, "Unexpected junk "//Trim(input)//" after key "//Trim(param%key))

  End Subroutine read_control_param

  Subroutine bad_option(key, option)
    Character(Len=*), Intent(In   ) :: key, option

    Call info('Unrecognised option '//Trim(option)//' for key '//Trim(key), .true.)
    Call error(3)

  End Subroutine bad_option

End Module new_control
