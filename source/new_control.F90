Module new_control

  Use angles,               Only: angles_type
  Use angular_distribution, Only: adf_type
  Use bspline,              Only: MIN_SPLINES, MAX_SPLINES
  Use bonds,                Only: bonds_type
  Use comms,                Only: comms_type,&
       gcheck
  Use configuration,        Only: configuration_type,&
       IMCON_NOPBC,&
       IMCON_SLAB
  Use constants,            Only: pi,&
       prsunt,&
       tenunt,&
       zero_plus
  Use constraints,          Only: constraints_type
  Use coord,                Only: coord_type
  Use core_shell,           Only: core_shell_type
  Use defects,              Only: defects_type
  Use development,          Only: development_type
  Use dihedrals,            Only: dihedrals_type
  Use electrostatic,        Only: ELECTROSTATIC_COULOMB,&
       ELECTROSTATIC_COULOMB_FORCE_SHIFT,&
       ELECTROSTATIC_COULOMB_REACTION_FIELD,&
       ELECTROSTATIC_DDDP,&
       ELECTROSTATIC_EWALD,&
       ELECTROSTATIC_NULL,&
       ELECTROSTATIC_POISSON,&
       electrostatic_type
  Use errors_warnings,      Only: error,&
       info,&
       warning,&
       set_print_level,&
       check_print_level
  Use ewald,                Only: ewald_type
  Use filename,             Only: FILE_CONFIG,&
       FILE_CONTROL,&
       FILE_FIELD,&
       FILE_HISTORF,&
       FILE_HISTORY,&
       FILE_OUTPUT,&
       FILE_REVCON,&
       FILE_REVIVE,&
       FILE_REVOLD,&
       FILE_STATS,&
       FILE_RDF,&
       FILE_MSD,&
       FILE_TABBND,&
       FILE_TABANG,&
       FILE_TABDIH,&
       FILE_TABINV,&
       FILE_TABVDW,&
       FILE_TABEAM,&
       file_type
  Use flow_control,         Only: RESTART_KEY_CLEAN,&
       RESTART_KEY_NOSCALE,&
       RESTART_KEY_OLD,&
       RESTART_KEY_SCALE,&
       flow_type, &
       DFTB, &
       MD_STD
  Use greenkubo,            Only: greenkubo_type
  Use impacts,              Only: impact_type
  Use inversions,           Only: inversions_type
  Use io,                   Only: &
       IO_READ_DIRECT, IO_READ_MASTER, IO_READ_MPIIO, IO_READ_NETCDF, &
       IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
       IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_NETCDF, &
       IO_WRITE_UNSORTED_DIRECT, IO_WRITE_UNSORTED_MASTER, &
       IO_WRITE_UNSORTED_MPIIO, io_get_parameters, io_nc_compiled, &
       io_nc_set_real_precision, io_set_parameters, io_type
  Use kim,                  Only: kim_type
  Use kinds,                Only: dp,&
       sp,&
       wi,&
       wp
  Use metal,                Only: metal_type
  Use minimise,             Only: MIN_DISTANCE,&
       MIN_ENERGY,&
       MIN_FORCE,&
       MIN_NULL,&
       minimise_type
  Use mpole,                Only: POLARISATION_CHARMM,&
       POLARISATION_DEFAULT,&
       mpole_type
  Use msd,                  Only: msd_type
  Use neighbours,           Only: neighbours_type
  Use netcdf_wrap,          Only: netcdf_param
  Use numerics,             Only: dcell,&
       invert,&
       seed_type
  Use parse,                Only: get_line,&
       get_word,&
       lower_case,&
       strip_blanks,&
       word_2_real
  Use plumed,               Only: plumed_type
  Use pmf,                  Only: pmf_type
  Use poisson,              Only: poisson_type
  Use rdfs,                 Only: rdf_type
  Use rigid_bodies,         Only: rigid_bodies_type
  Use rsds,                 Only: rsd_type
  Use statistics,           Only: stats_type
  Use tersoff,              Only: tersoff_type
  Use thermostat,           Only: &
       CONSTRAINT_NONE, CONSTRAINT_SEMI_ORTHORHOMBIC, &
       CONSTRAINT_SURFACE_AREA, CONSTRAINT_SURFACE_TENSION, &
       DPD_FIRST_ORDER, DPD_NULL, DPD_SECOND_ORDER, ENS_NPT_BERENDSEN, &
       ENS_NPT_BERENDSEN_ANISO, ENS_NPT_LANGEVIN, &
       ENS_NPT_LANGEVIN_ANISO, ENS_NPT_MTK, ENS_NPT_MTK_ANISO, &
       ENS_NPT_NOSE_HOOVER, ENS_NPT_NOSE_HOOVER_ANISO, ENS_NVE, &
       ENS_NVT_ANDERSON, ENS_NVT_BERENDSEN, ENS_NVT_EVANS, &
       ENS_NVT_GENTLE, ENS_NVT_LANGEVIN, ENS_NVT_LANGEVIN_INHOMO, &
       ENS_NVT_NOSE_HOOVER, thermostat_type
  Use timer,                Only: timer_type
  Use trajectory,           Only: trajectory_type,&
                                  TRAJ_KEY_COORD,&
                                  TRAJ_KEY_COORD_VEL,&
                                  TRAJ_KEY_COORD_VEL_FORCE,&
                                  TRAJ_KEY_COMPRESSED
  Use ttm,                  Only: ttm_type
  Use vdw,                  Only: MIX_FENDER_HALSEY,&
       MIX_FUNCTIONAL,&
       MIX_HALGREN,&
       MIX_HOGERVORST,&
       MIX_LORENTZ_BERTHELOT,&
       MIX_NULL,&
       MIX_TANG_TOENNIES,&
       MIX_WALDMAN_HAGLER,&
       vdw_type
  Use z_density,            Only: z_density_type
  Use comms, only: comms_type
  Use hash, only: hash_table, MAX_KEY, STR_LEN
  Use parse, only: get_word, get_line, lower_case
  Use filename, Only : file_type,FILE_CONTROL,FILE_OUTPUT,FILE_CONFIG,FILE_FIELD, &
       FILE_STATS,FILE_HISTORY,FILE_HISTORF,FILE_REVIVE,FILE_REVCON, &
       FILE_REVOLD
  Use errors_warnings, only: error, info, warning, error_alloc, error_dealloc
  Use io,     Only : io_set_parameters,io_type,        &
       io_get_parameters,        &
       io_nc_set_real_precision, &
       io_nc_compiled,           &
       IO_READ_MPIIO,            &
       IO_READ_DIRECT,           &
       IO_READ_MASTER,           &
       IO_READ_NETCDF,           &
       IO_WRITE_UNSORTED_MPIIO,  &
       IO_WRITE_UNSORTED_DIRECT, &
       IO_WRITE_UNSORTED_MASTER, &
       IO_WRITE_SORTED_MPIIO,    &
       IO_WRITE_SORTED_DIRECT,   &
       IO_WRITE_SORTED_NETCDF,   &
       IO_WRITE_SORTED_MASTER
  Use three_body,      Only: threebody_type
  Use four_body,       Only: four_body_type
  Use units, only : convert_units, set_timestep, units_scheme, internal_units, set_out_units
  Use control_parameter_module, only : control_parameter, parameters_hash_table, &
       DATA_INT, DATA_FLOAT, DATA_STRING, DATA_BOOL, DATA_OPTION, DATA_VECTOR3, DATA_VECTOR6

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
contains

  Subroutine read_new_control(control_file, params, comm, can_parse)
    !!-----------------------------------------------------------------------
    !!
    !! Set read in a new_style control file
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), intent(   Out ) :: params
    Type( file_type ), Intent( InOut ) :: control_file
    Type( comms_type ), Intent( InOut ) :: comm
    Logical, Intent(   Out) :: can_parse


    Call initialise_control(params)
    open(newunit=control_file%unit_no, file=control_file%filename, status='old', action='read')
    can_parse = try_parse(control_file%unit_no, params, comm)

    ! Possibly old style
    if (.not. can_parse) then
      call control_file%close()
      return
    end if

    Call parse_control_file(control_file%unit_no, params, comm)
    call control_file%close()

  end Subroutine read_new_control

  Subroutine read_devel(params, devel, tmr, seed)
    Type( parameters_hash_table ), Intent( In    ) :: params
    Type( development_type ),      Intent( InOut ) :: devel
    Type( timer_type),             Intent( InOut ) :: tmr
    Type( seed_type),              Intent( InOut ) :: seed
    Real(kind=wp), Dimension(3) :: vtmp
    Character(len=STR_LEN) :: option
    Character(len=STR_LEN) :: word
    Integer :: print_level

    Call params%retrieve('print_level', print_level)
    Call set_print_level(print_level)
    Call params%retrieve('unsafe_comms', devel%l_fast)

    call params%retrieve('output_energy', devel%l_eng)
    call params%retrieve('io_write_ascii_revive', devel%l_rout)
    call params%retrieve('io_read_ascii_revold', devel%l_rin)
    call params%retrieve('initial_minimum_separation', devel%r_dis)

    call params%retrieve('timer_depth', tmr%max_depth)
    call params%retrieve('timer_per_mpi', tmr%proc_detail)

    call params%retrieve('time_job', tmr%job)
    if (tmr%job < 0.0_wp) tmr%job = Huge(1.0_wp)

    call params%retrieve('time_close', tmr%clear_screen)
    if (tmr%clear_screen < 0.0_wp) tmr%clear_screen = 0.01_wp * tmr%job

    if (params%is_set('random_seed')) then
      call params%retrieve('random_seed', vtmp(1:3))
      call seed%init(nint(vtmp(1:3)))
    end if

    call params%retrieve('dftb_test', devel%test_dftb_library)

  end Subroutine read_devel

  Subroutine read_io(params, io_data, netcdf, files, comm)
    !!-----------------------------------------------------------------------
    !!
    !! Read in the io parameters
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), Intent( In    ) :: params
    Type( io_type ), Intent( InOut ) :: io_data
    Type( netcdf_param ), Intent( InOut ) :: netcdf
    Type( file_type ), Dimension(:), Intent( InOut ) :: files
    Type( comms_type ), Intent( InOut ) :: comm
    Character(Len=STR_LEN) :: curr_option
    Character(Len=STR_LEN) :: message

    Integer, Parameter :: MAX_BATCH_SIZE = 10000000, MAX_BUFFER_SIZE = 100000

    Integer :: io_read, io_write
    Real(kind=wp) :: rtmp
    Integer :: itmp
    Logical :: ltmp

    call params%retrieve('io_read_method', curr_option)

    Select Case(curr_option)
    Case ( 'mpiio' )
      io_read = IO_READ_MPIIO

    Case ( 'direct' )
      io_read = IO_READ_DIRECT

    Case ( 'netcdf' )
      io_read = IO_READ_NETCDF

    Case ( 'master' )
      io_read = IO_READ_MASTER

    Case Default
      Call bad_option('io_read_method', curr_option)

    End Select

    Call io_set_parameters(io_data, user_method_read = io_read)

    Select Case (io_read)
    Case (IO_READ_MPIIO, IO_READ_DIRECT, IO_READ_NETCDF)
      ! Need to calculate number of readers
      call params%retrieve('io_read_readers', itmp)

      If (itmp == 0) then
        rtmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
        itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
        Do While ( Mod( comm%mxnode, itmp ) /= 0 )
          itmp = itmp - 1
        End Do
      else if (itmp < comm%mxnode) then
        Do While ( Mod( comm%mxnode, itmp ) /= 0 )
          itmp = itmp - 1
        End Do
      else if (itmp > comm%mxnode) then
        rtmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
        itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
        Do While ( Mod( comm%mxnode, itmp ) /= 0 )
          itmp = itmp - 1
        End Do
      else
        Call error(0, 'Cannot have negative number of I/O readers')
      end If


      ! the number of readers is now ready to set

      Call io_set_parameters(io_data, user_n_io_procs_read = itmp)

      ! Sort read batch size
      ! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
      ! Note zero or negative values indicate use the default

      call params%retrieve('io_read_batch_size', itmp)

      Select Case (itmp)
      Case(0)
        Call io_get_parameters(io_data, user_batch_size_read = itmp )

      Case(1:)
        itmp = Min( itmp, MAX_BATCH_SIZE )
        Call io_set_parameters(io_data, user_batch_size_read = itmp )

      Case Default
        Call error(0, 'Cannot have negative I/O read batch size')

      End Select

    Case(IO_READ_MASTER)
      Write(message,'(a,i10)') 'I/O readers (enforced) ', 1
      Call info(message,.true.)

    Case Default
      rtmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
      itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
      ! the number of readers is now ready to set
      Call io_set_parameters(io_data, user_n_io_procs_read = itmp )

    end Select

    ! Get read buffer size

    call params%retrieve('io_read_buffer_size', itmp)

    Select Case(itmp)
    Case(0)
      Call io_get_parameters(io_data, user_buffer_size_read = itmp )

    Case(1:99)
      Call io_set_parameters(io_data, user_buffer_size_read = 100 )

    Case(100:MAX_BUFFER_SIZE)
      Call io_set_parameters(io_data, user_buffer_size_read = itmp )

    Case(MAX_BUFFER_SIZE+1:)
      Call io_set_parameters(io_data, user_buffer_size_read = MAX_BUFFER_SIZE )

    Case Default
      Call error(0, 'Negative read buffer size not valid, min = 1')

    End Select

    ! Get parallel read error checking

    if (io_read /= IO_READ_MASTER) then
      call params%retrieve('io_read_error_check', ltmp)
      Call io_set_parameters(io_data, user_error_check = ltmp )

    End If

    ! Get write settings

    call params%retrieve('io_write_method', curr_option)
    call params%retrieve('io_write_sorted', ltmp)

    Select Case (curr_option)
    Case ( 'mpiio' )
      if (ltmp) then
        io_write = IO_WRITE_SORTED_MPIIO
      else
        io_write = IO_WRITE_UNSORTED_MPIIO
      end if

    Case ( 'direct' )
      if (ltmp) then
        io_write = IO_WRITE_SORTED_DIRECT
      else
        io_write = IO_WRITE_UNSORTED_DIRECT
      end if

    Case ( 'netcdf' )
      io_write = IO_WRITE_SORTED_NETCDF

      call params%retrieve('io_read_netcdf_format', curr_option)

      Select Case (curr_option)
      Case('amber', '32bit', '32-bit')
        ! Use 32-bit quantities in output for real numbers
        Call io_nc_set_real_precision( sp, netcdf, itmp )
      Case('64-bit', '64bit')
        ! Use 64-bit quantities in output for real numbers
        Call io_nc_set_real_precision( dp, netcdf, itmp )
      Case Default
        Call error(3)
      End Select

    Case ( 'master' )

      if (ltmp) then
        io_write = IO_WRITE_SORTED_MASTER
      else
        io_write = IO_WRITE_UNSORTED_MASTER
      end if

    Case Default
      call bad_option('io_write_method', curr_option)

    End Select


    ! the write method and type are now ready to set

    Call io_set_parameters(io_data, user_method_write = io_write )

    Select Case (io_write)
    Case ( IO_WRITE_SORTED_NETCDF, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT )
      call params%retrieve('io_write_writers', itmp)

      if (itmp == 0) then
        rtmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
        itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
        Do While ( Mod( comm%mxnode, itmp ) /= 0 )
          itmp = itmp - 1
        End Do

      else if (itmp < comm%mxnode) then
        Do While ( Mod( comm%mxnode, itmp ) /= 0 )
          itmp = itmp - 1
        End Do

      else if (itmp > comm%mxnode) then
        rtmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
        itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
        Do While ( Mod( comm%mxnode, itmp ) /= 0 )
          itmp = itmp - 1
        End Do

      else
        Call error(0, 'Cannot have negative number of I/O writers')

      End if

      ! the number of writers is now ready to set

      Call io_set_parameters(io_data, user_n_io_procs_write = itmp )

      call params%retrieve('io_write_batch_size', itmp)

      Select Case (itmp)
      Case(0)
        Call io_get_parameters(io_data, user_batch_size_write = itmp )

      Case(1:)
        itmp = Min( itmp, MAX_BATCH_SIZE )
        Call io_set_parameters(io_data, user_batch_size_write = itmp )

      Case Default
        Call error(0, 'Cannot have negative I/O write batch size')

      end Select

    Case (IO_WRITE_UNSORTED_MASTER, IO_WRITE_SORTED_MASTER)
      Continue

    Case Default
      rtmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
      itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
      ! the number of writers is now ready to set
      Call io_set_parameters(io_data, user_n_io_procs_write = itmp )

    End Select

    call params%retrieve('io_write_buffer_size', itmp)

    Select Case(itmp)
    Case(0)
      Call io_get_parameters(io_data, user_buffer_size_write = itmp )

    Case(1:99)
      Call io_set_parameters(io_data, user_buffer_size_write = 100 )

    Case(100:MAX_BUFFER_SIZE)
      Call io_set_parameters(io_data, user_buffer_size_write = itmp )

    Case(MAX_BUFFER_SIZE+1:)
      Call io_set_parameters(io_data, user_buffer_size_write = MAX_BUFFER_SIZE )

    Case Default
      Call error(0, 'Negative write buffer size not valid, min = 1')
    end Select

    ! switch error checking flag for writing

    If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
      call params%retrieve('io_write_error_check', ltmp)
      Call io_set_parameters(io_data, user_error_check = ltmp )
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
    Type( parameters_hash_table ), intent( In    ) :: params
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Integer, Intent( In    ) :: nvdw
    Logical, Intent( In    ) :: ttm_active
    Character(Len=STR_LEN) :: option
    Logical :: ltmp

    call params%retrieve('ensemble', option, required = .true.)
    if (ttm_active .and. option /= 'ttm') call error(0, 'TTM requested, ensemble not ttm')
    select case(option)
    case ('nve', 'pmf')
      thermo%ensemble = ENS_NVE

    case ('nvt')

      call params%retrieve('ensemble_method', option, required = .true.)

      select case(option)
      case ('evans')
        thermo%ensemble = ENS_NVT_EVANS

      case ('langevin')

        thermo%ensemble = ENS_NVT_LANGEVIN

        Call params%retrieve('ensemble_thermostat_friction', thermo%chi, .true.)

      case ('andersen')

        thermo%ensemble = ENS_NVT_ANDERSON

        Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        Call params%retrieve('ensemble_thermostat_softness', thermo%soft)

      case ('berendsen')

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

        Call info('Ensemble: NVT dpd (Dissipative Particle Dynamics)',.true.)

        call params%retrieve('ensemble_dpd_order', option)
        select case (option)
        case ('first', '1')
          thermo%key_dpd = DPD_FIRST_ORDER
        case ('second', '2')
          thermo%key_dpd = DPD_SECOND_ORDER
        case default
          call bad_option('ensemble_dpd_order', option)
        end select

        ! ALLOCATE DPD ARRAYS
        call thermo%init_dpd(nvdw + 1) ! Account for extra one in bounds?!

        call params%retrieve('ensemble_dpd_drag', thermo%gamdpd(0))

      case default
        call bad_option('NVT ensemble_method', option)
      end select

    case ('npt')

      thermo%variable_cell = .true.

      call params%retrieve('ensemble_method', option)
      select case (option)
      case ('langevin')

        thermo%ensemble = ENS_NPT_LANGEVIN
        thermo%l_langevin = .true.

        call params%retrieve('ensemble_thermostat_friction', thermo%chi, .true.)
        call params%retrieve('ensemble_barostat_friction', thermo%tai, .true.)

      case ('berendsen')

        thermo%ensemble = ENS_NPT_BERENDSEN

        call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      case ('hoover', 'nose', 'nose-hoover')


        thermo%ensemble = ENS_NPT_NOSE_HOOVER

        call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      case ('mtk')

        thermo%ensemble = ENS_NPT_MTK

        call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      case default
        call bad_option('NPT ensemble_method', option)
      end select

    case ('nst')

      thermo%variable_cell = .true.
      thermo%anisotropic_pressure = .true.

      call params%retrieve('ensemble_method', option)
      select case (option)
      case ('langevin')

        thermo%ensemble = ENS_NPT_LANGEVIN_ANISO
        thermo%l_langevin = .true.

        call params%retrieve('ensemble_thermostat_friction', thermo%chi, .true.)
        call params%retrieve('ensemble_barostat_friction', thermo%tai, .true.)

      case ('berendsen')

        thermo%ensemble = ENS_NPT_BERENDSEN_ANISO

        call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      case ('hoover', 'nose', 'nose-hoover')

        thermo%ensemble = ENS_NPT_NOSE_HOOVER_ANISO

        call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      Case ('mtk')

        thermo%ensemble = ENS_NPT_MTK_ANISO

        call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
        call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

      case default
        call bad_option('NST ensemble method', option)
      end select

      ! Semi isotropic ensembles

      call params%retrieve('ensemble_semi_isotropic', option)
      call params%retrieve('ensemble_semi_orthorhombic', ltmp)

      select case (option)
      case ('off')
        continue
      case ('area')
        thermo%iso = CONSTRAINT_SURFACE_AREA
      case ('surface', 'tension')

        thermo%iso = CONSTRAINT_SURFACE_TENSION

        call params%retrieve('ensemble_tension', thermo%tension, required=.true.)

        thermo%tension=thermo%tension/tenunt

        if (ltmp) then
          thermo%iso = CONSTRAINT_SEMI_ORTHORHOMBIC
        end if

      case ('ortho', 'orthorhombic')

        thermo%iso = CONSTRAINT_SURFACE_TENSION
        if (ltmp) then
          thermo%iso = CONSTRAINT_SEMI_ORTHORHOMBIC
        end if
      case default
        call bad_option('ensemble_semi_isotropic', option)
      end select

    case default
      call bad_option('ensemble', option)
    end select

  end Subroutine read_ensemble

  Subroutine read_bond_analysis(params, flow, bond, angle, dihedral, inversion, max_grid_analysis)
    Type(parameters_hash_table ), Intent(In   ) :: params
    Type(flow_type),              Intent(InOut) :: flow
    Type(bonds_type),             Intent(InOut) :: bond
    Type(angles_type),            Intent(InOut) :: angle
    Type(dihedrals_type),         Intent(InOut) :: dihedral
    Type(inversions_type),        Intent(InOut) :: inversion
    Integer,                      Intent(  Out) :: max_grid_analysis
    Logical :: ltmp
    Integer :: itmp, itmp2

    Real( kind = wp ), Parameter :: minimum_bond_anal_length = 2.5_wp

    call params%retrieve('analyse_bonds', flow%analyse_bond)
    call params%retrieve('analyse_angles', flow%analyse_ang)
    call params%retrieve('analyse_dihedrals', flow%analyse_dih)
    call params%retrieve('analyse_inversions', flow%analyse_inv)
    call params%retrieve('analyse_all', ltmp)
    if (ltmp) then
      flow%analyse_bond = .true.
      flow%analyse_ang = .true.
      flow%analyse_dih = .true.
      flow%analyse_inv = .true.
    end if

    ! Set global
    call params%retrieve('analyse_frequency', itmp)
    call params%retrieve('analyse_num_bins', itmp2)

    if (flow%analyse_bond) then
      call params%retrieve('analyse_max_dist', bond%rcut)
      if (bond%rcut < minimum_bond_anal_length) then
        bond%rcut = minimum_bond_anal_length
        call warning('Bond less than minimum length, setting to minimum (2.5 Ang)', .true.)
      end if
      bond%bin_pdf = itmp2
      if (params%is_set('analyse_num_bins_bonds')) &
           call params%retrieve('analyse_num_bins_bonds', bond%bin_pdf)

      call params%retrieve('analyse_frequency_bonds', flow%freq_bond)
    else
      bond%bin_pdf = -1
      flow%freq_bond = -1
    end if

    if (flow%analyse_ang) then
      angle%bin_adf = itmp2
      if (params%is_set('analyse_num_bins_angles')) &
           call params%retrieve('analyse_num_bins_angles', angle%bin_adf)
      call params%retrieve('analyse_frequency_angles', flow%freq_angle)
    else
      angle%bin_adf = -1
      flow%freq_angle = -1
    end if

    if (flow%analyse_dih) then
      dihedral%bin_adf = itmp2
      if (params%is_set('analyse_num_bins_dihedrals')) &
           call params%retrieve('analyse_num_bins_dihedrals', dihedral%bin_adf)
      call params%retrieve('analyse_frequency_dihedrals', flow%freq_dihedral)
    else
      dihedral%bin_adf = -1
      flow%freq_dihedral = -1
    end if

    if (flow%analyse_inv) then
      inversion%bin_adf = itmp2
      if (params%is_set('analyse_num_bins_inversions')) &
           call params%retrieve('analyse_num_bins_inversions', inversion%bin_adf)

      call params%retrieve('analyse_frequency_inversions', flow%freq_inversion)
    else
      inversion%bin_adf = -1
      flow%freq_inversion = -1
    end if

    if (any([flow%analyse_bond, flow%analyse_ang, flow%analyse_dih, flow%analyse_inv])) then
      max_grid_analysis = Max(bond%bin_pdf, angle%bin_adf, dihedral%bin_adf, inversion%bin_adf)
    end if

    flow%freq_bond = max(1, itmp, flow%freq_bond)
    flow%freq_angle = max(1, itmp, flow%freq_angle)
    flow%freq_dihedral = max(1, itmp, flow%freq_dihedral)
    flow%freq_inversion = max(1, itmp, flow%freq_inversion)

  End Subroutine read_bond_analysis

  Subroutine read_structure_analysis(params, stats, msd_data, rdf, vaf, zden, adf, coords, traj, defect, displacement)

    Type(parameters_hash_table), Intent(In   ) :: params
    Type(stats_type),            Intent(InOut) :: stats
    Type(msd_type),              Intent(InOut) :: msd_data
    Type(greenkubo_type),        Intent(InOut) :: vaf
    Type(rdf_type),              Intent(InOut) :: rdf
    Type(z_density_type),        Intent(InOut) :: zden
    Type(adf_type),              Intent(InOut) :: adf
    Type(coord_type),            intent(InOut) :: coords
    Type(trajectory_type),       Intent(InOut) :: traj
    Type(defects_type),          Intent(InOut) :: defect(:)
    Type(rsd_type),              Intent(InOut) :: displacement

    Character(Len=STR_LEN) :: option
    Logical :: ltmp

    ! MSD
    call params%retrieve('msd_calculate', msd_data%l_msd)

    if (msd_data%l_msd) then
      call params%retrieve('msd_start', msd_data%start)
      call params%retrieve('msd_frequency', msd_data%freq)
    else if (params%is_any_set([Character(13) :: 'msd_frequency', 'msd_start'])) then
      Call warning('msd_start or msd_frequency found without msd_calculate')
    end if

    ! VAF
    call params%retrieve('vaf_calculate', vaf%l_collect)

    if (vaf%l_collect) then
      call params%retrieve('vaf_frequency', vaf%freq)
      If (vaf%freq <= 0) vaf%freq=50
      call params%retrieve('vaf_binsize', vaf%binsize)
      If (vaf%binsize <= 0) vaf%binsize=Merge(2*vaf%freq,100,vaf%freq >= 100)
      call params%retrieve('vaf_print', vaf%l_print)

      vaf%samp = Ceiling(Real(vaf%binsize,wp)/Real(vaf%freq,wp))

    else if (params%is_any_set([Character(13) :: 'vaf_frequency', 'vaf_binsize', 'vaf_print'])) then
      Call warning('vaf_print, vaf_frequency or vaf_binsize found without vaf_calculate')
    end if

    ! RDF
    call params%retrieve('rdf_calculate', rdf%l_collect)
    ! Necessary for some reason
    call params%retrieve('rdf_binsize', rdf%rbin)

    if (rdf%l_collect) then
      call params%retrieve('rdf_error_analysis', option)

      select case (option)
      case ('jackknife')
        rdf%l_errors_jack = .TRUE.
      case ('block')
        rdf%l_errors_block = .TRUE.
      case ('off')
        continue
      case default
        call bad_option('rdf_error_analysis', option)
      end select

      if (rdf%l_errors_jack .or. rdf%l_errors_block) &
           & call params%retrieve('rdf_error_analysis_blocks', rdf%num_blocks)

      call params%retrieve('rdf_frequency', rdf%freq)
      call params%retrieve('rdf_print', rdf%l_print)

    else if (params%is_any_set([Character(13) :: 'rdf_frequency', 'rdf_binsize', 'rdf_print'])) then
      Call warning('rdf_print, rdf_frequency or rdf_binsize found without rdf_calculate')
    end if

    Call params%retrieve('print_probability_distribution', stats%lpana)
    if (stats%lpana .and. .not. rdf%l_print) then
      Call warning('RDF printing triggered due to a PDA printing request', .true.)
      rdf%l_print = stats%lpana
    end if

    ! ZDen

    call params%retrieve('zden_calculate', zden%l_collect)

    if (zden%l_collect) then

      call params%retrieve('zden_frequency', zden%frequency)
      call params%retrieve('zden_binsize', zden%bin_width)
      call params%retrieve('zden_print', zden%l_print)

    else if (params%is_any_set([Character(14) :: 'zden_frequency', 'zden_binsize', 'zden_print'])) then
      Call warning('zden_print, zden_frequency or zden_binsize found without zden_calculate')
    end if

    ! ADF

    call params%retrieve('adf_calculate', adf%adfon)

    if (adf%adfon) then

      call params%retrieve('adf_frequency', adf%interval)
      call params%retrieve('adf_precision', adf%prec)

    else if (params%is_any_set([Character(13) :: 'adf_frequency', 'adf_precision'])) then
      Call warning('adf_frequency or adf_precision found without adf_calculate')
    end if

    ! Coord

    call params%retrieve('coord_calculate', coords%coordon)

    if (coords%coordon) then
      call params%retrieve('coord_start', coords%coordstart)
      call params%retrieve('coord_interval', coords%coordinterval)
      call params%retrieve('coord_ops', coords%coordops)

    else if (params%is_any_set([Character(14) :: 'coord_start', 'coord_interval', 'coord_ops'])) then
      Call warning('coord_start, coord_interval or coord_ops found without coord_calculate')
    end if

    ! Trajectory

    call params%retrieve('traj_calculate', traj%ltraj)

    if (traj%ltraj) then
      call params%retrieve('traj_start', traj%start)
      call params%retrieve('traj_interval', traj%freq)
      call params%retrieve('traj_key', option)

      select case (option)
      case ('pos')
        traj%key = TRAJ_KEY_COORD
      case ('pos-vel')
        traj%key = TRAJ_KEY_COORD_VEL
      case ('pos-vel-force')
        traj%key = TRAJ_KEY_COORD_VEL_FORCE
      case ('compressed')
        traj%key = TRAJ_KEY_COMPRESSED
      case default
        call bad_option('traj_key', option)
      end select

      ! Need to dealias
      call traj%init((traj%key), (traj%freq), (traj%start))

    else if (params%is_any_set([Character(13) :: 'traj_start', 'traj_interval', 'traj_key'])) then
      Call warning('traj_start, traj_interval or traj_key found without traj_calculate')
    end if

    call params%retrieve('defects_calculate', defect(1)%ldef)

    if (defect(1)%ldef) then

      call params%retrieve('defects_start', defect(1)%nsdef)
      call params%retrieve('defects_interval', defect(1)%isdef)
      call params%retrieve('defects_distance', defect(1)%rdef)
      ! if (defects%rdef < )

      defect(1)%newjob = .true.
      ! Name REFERENCE and DEFECTS files
      defect(1)%reffile = 'REFERENCE'
      defect(1)%deffile = 'DEFECTS'

      call params%retrieve('defects_backup', ltmp)
      if (ltmp) then
        defect(2)%ldef = .true.
        defect(2)%nsdef = defect(1)%nsdef
        defect(2)%isdef = defect(1)%isdef
        defect(2)%rdef =  defect(1)%rdef
        defect(2)%newjob = .true.
        defect(2)%reffile = 'REFERENCE1'
        defect(2)%deffile = 'DEFECTS1'
      end if

    else if (params%is_any_set([Character(16) :: 'defects_start', 'defects_interval', 'defects_distance'])) then
      Call warning('defects_start, defects_interval or defects_distance found without defects_calculate')
    end if

    call params%retrieve('displacements_calculate', displacement%lrsd)

    if (displacement%lrsd) then

      call params%retrieve('displacements_start', displacement%nsrsd)
      call params%retrieve('displacements_interval', displacement%isrsd)
      call params%retrieve('displacements_distance', displacement%rrsd)
      if (displacement%rrsd < 0.15_wp) then
        displacement%rrsd = 0.15_wp
        call warning('Displacement_distance too small, reset to 0.15 ang')
      end if
    else if (params%is_any_set([Character(22) :: 'displacements_start', 'displacements_interval', 'displacements_distance'])) then
      Call warning('displacements_start, displacements_interval or displacements_distance found without displacements_calculate')
    end if

  end Subroutine read_structure_analysis

  Subroutine read_ttm(params, ttm)
    Type( parameters_hash_table ), Intent( In    ) :: params
    Type(ttm_type),           Intent(   Out ) :: ttm
    Character(Len=STR_LEN) :: option
    Logical :: ltmp


    call params%retrieve('ttm_calculate', ttm%l_ttm)
    if (.not. ttm%l_ttm) return

    call params%retrieve('ttm_num_ion_cells', ttm%ntsys(3))
    call params%retrieve('ttm_num_elec_cells', ttm%eltsys)
    call params%retrieve('ttm_metal', ttm%ismetal)
    call params%retrieve('ttm_dens_model', option)

    call params%retrieve('ttm_heat_cap_model', option)
    select case (option)
    case ('constant')
      ttm%cetype = 0
      call params%retrieve('ttm_heat_cap', ttm%ce0)
    case ('tanh')
      ttm%cetype = 1
      call params%retrieve('ttm_heat_cap', ttm%sh_A)
      call params%retrieve('ttm_temp_term', ttm%sh_B)

      If (ttm%sh_A <= zero_plus .or. ttm%sh_B <= zero_plus) &
           Call error(0, 'Electronic specific heat not fully specified')
    case ('linear')
      ttm%cetype = 2
      call params%retrieve('ttm_heat_cap', ttm%Cemax)
      call params%retrieve('ttm_fermi_temp', ttm%Tfermi)

      If (ttm%Tfermi <= zero_plus .or. ttm%Cemax <= zero_plus) &
           Call error(0, 'Electronic specific heat not fully specified')
    case ('tabulated')
      ttm%cetype = 3
    case default
      call bad_option('ttm_heat_cap_model', option)
    end select

    select case (option)
    case ('constant')
      ttm%ttmdyndens = .false.
      call params%retrieve('ttm_dens', ttm%cellrho, .true.)
      If (ttm%cellrho <= zero_plus) call error(0, 'Bad ttm_dens (<= 0)')
      ttm%rcellrho = 1.0_wp / ttm%cellrho

      ! Rescale
      ttm%sh_A = ttm%sh_A * ttm%cellrho
      ttm%Cemax = ttm%Cemax * ttm%cellrho
      ttm%epc_to_chi = convert_units(1.0e-12_wp * ttm%rcellrho / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')

    case ('dynamic')
      ttm%ttmdyndens = .true.
      ttm%CeType = ttm%CeType + 4
      ttm%epc_to_chi = convert_units(1.0e-12_wp / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')

    case default
      call bad_option('ttm_dens_model', option)
    end select

    if (ttm%ismetal) then
      ttm%detype = 0

      call params%retrieve('ttm_elec_cond_model', option, required=.true.)
      select case (option)
      case ('infinite')
        ttm%ketype = 0
      case ('constant')
        ttm%ketype = 1
        call params%retrieve('ttm_elec_cond', ttm%ka0)
        If (ttm%ka0 <= zero_plus) &
             Call error(0, 'Electronic thermal conductivity not fully specified')
      case ('drude')
        ttm%ketype = 2
        call params%retrieve('ttm_elec_cond', ttm%ka0)

        If (ttm%ka0 <= zero_plus) &
             Call error(0, 'Electronic thermal conductivity not fully specified')
      case ('tabulated')
        ttm%ketype = 3
      case default
        call bad_option('ttm_elec_cond_model', option)
      end select
    else
      ttm%ketype = 0

      call params%retrieve('ttm_diff_model', option, required=.true.)
      select case (option)
      case ('constant')
        ttm%detype = 1
        call params%retrieve('ttm_diff', ttm%diff0)
        if (ttm%diff0 <= zero_plus) &
             Call error(0, 'Thermal diffusivity of non-metal not specified')
      case ('recip', 'reciprocal')
        ttm%detype = 2
        call params%retrieve('ttm_diff', ttm%diff0)
        call params%retrieve('ttm_fermi_temp', ttm%Tfermi)
        if (ttm%diff0 <= zero_plus .or. ttm%Tfermi <= zero_plus) &
             Call error(0, 'Thermal diffusivity of non-metal not specified')
      case ('tabulated')
        ttm%detype = 3
      case default
        call bad_option('ttm_diff_model', option)
      end select

    end if

    call params%retrieve('ttm_variable_ep', option, required=.true.)
    select case (option)
    case ('homo')
      ttm%gvar = 1
    case ('hetero')
      ttm%gvar = 2
    case default
      call bad_option('ttm_variable_ep', option)
    end select

    call params%retrieve('ttm_com_correction', option)
    select case (option)
    case ('full')
      ttm%ttmthvel = .true.
      ttm%ttmthvelz = .false.
    case ('zdir')
      ttm%ttmthvel = .true.
      ttm%ttmthvelz = .true.
    case ('off')
      ttm%ttmthvel = .false.
      ttm%ttmthvelz = .false.
    case default
      call bad_option('ttm_com_correction', option)
    end select

    call params%retrieve('ttm_redistribute', ttm%redistribute)


    call params%retrieve('ttm_min_atoms', ttm%amin)
    ttm%amin = Max(ttm%amin, 1) ! minimum number of atoms for ttm ionic temperature cell

    call params%retrieve('ttm_stopping_power', ttm%dedx)

    call params%retrieve('ttm_spatial_dist', option)
    select case (option)
    case ('gaussian')
      ttm%sdepoType = 1
      call params%retrieve('ttm_spatial_sigma', ttm%sig)
      call params%retrieve('ttm_spatial_cutoff', ttm%sigmax)

    case ('flat')
      ttm%sdepoType = 2

    case ('laser')
      call params%retrieve('ttm_laser_type', option)
      call params%retrieve('ttm_fluence', ttm%fluence)
      call params%retrieve('ttm_penetration_depth', ttm%pdepth)

      select case (option)
      case ('flat')
        ttm%sdepoType = 2

      case ('exponential')
        ttm%sdepoType = 3

      case default
        call bad_option('ttm_laser_type', option)
      end select

    case default
      call bad_option('ttm_spatial_dist', option)
    end select

    call params%retrieve('ttm_temporal_dist', option)
    select case (option)
    case ('gaussian')
      ttm%tdepotype = 1
      call params%retrieve('ttm_temporal_duration', ttm%tdepo)
      call params%retrieve('ttm_temporal_cutoff', ttm%tcdepo)

    case ('exponential')
      ttm%tdepotype = 2
      call params%retrieve('ttm_temporal_duration', ttm%tdepo)
      call params%retrieve('ttm_temporal_cutoff', ttm%tcdepo)

    case ('delta')
      ttm%tdepotype = 3

    case ('square')
      ttm%tdepotype = 4
      call params%retrieve('ttm_temporal_duration', ttm%tdepo)
      If (ttm%tdepo <= zero_plus) Then
        ttm%tdepoType = 3
      End If

    case default
      call bad_option('ttm_temporal_dist', option)
    end select

    call params%retrieve('ttm_boundary_condition', option)
    call params%retrieve('ttm_boundary_xy', ltmp)
    select case (option)
    case ('periodic')
      ttm%bctypee = 1

    case ('dirichlet')
      if (ltmp) then
        ttm%bcTypeE = 4
      else
        ttm%bcTypeE = 2
      end if

    case ('neumann')
      ttm%bcTypeE = 3
    case ('robin')

      call params%retrieve('ttm_boundary_heat_flux', ttm%fluxout)

      if (ltmp) then
        ttm%bcTypeE = 6
      else
        ttm%bcTypeE = 5
      end if

    case default
      call bad_option('ttm_boundary_condition', option)
    end select

    call params%retrieve('ttm_time_offset', ttm%ttmoffset)
    call params%retrieve('ttm_oneway', ttm%oneway)
    call params%retrieve('ttm_stats_frequency', ttm%ttmstats)
    call params%retrieve('ttm_traj_frequency', ttm%ttmtraj)

  end Subroutine read_ttm

  Subroutine read_units(params)
    Type( parameters_hash_table ), intent( In    ) :: params
    Type( units_scheme ) :: out_units
    Character(Len=STR_LEN) :: option
    Real(kind = wp) :: test

    call params%retrieve("io_units_scheme", option)

    select case (option)
    case ('internal')
      out_units = internal_units

    case ('si')
      out_units%length   = 'm'
      out_units%time     = 's'
      out_units%mass     = 'kg'
      out_units%charge   = 'C'
      out_units%energy   = 'J'
      out_units%pressure = 'Pa'
      out_units%force    = 'N'
      out_units%velocity = 'm/s'
      out_units%power    = 'W'
      out_units%surf_ten = 'N/m'
      out_units%emf      = 'V'

    case ('atomic')
      out_units%length   = 'ang'
      out_units%time     = 'ps'
      out_units%mass     = 'amu'
      out_units%charge   = 'q_e'
      out_units%energy   = 'e.V'
      out_units%pressure = 'GPa'
      out_units%force    = 'e.V/ang'
      out_units%velocity = 'ang/ps'
      out_units%power    = 'e.V/ps'
      out_units%surf_ten = 'e.V/ang^2'
      out_units%emf      = 'e.V/q_e'

    case ('hartree')
      out_units%length   = 'bohr'
      out_units%time     = 'aut'
      out_units%mass     = 'm_e'
      out_units%charge   = 'q_e'
      out_units%energy   = 'Ha'
      out_units%pressure = 'Ha/bohr^3'
      out_units%force    = 'Ha/bohr'
      out_units%velocity = 'auv'
      out_units%power    = 'Ha/aut'
      out_units%surf_ten = 'Ha/bohr^2'
      out_units%emf      = 'Ha/q_e'

    end select

    if (params%is_set("io_units_length")) &
         call params%retrieve("io_units_length", out_units%length)
    if (params%is_set("io_units_time")) &
         call params%retrieve("io_units_time", out_units%time)
    if (params%is_set("io_units_mass")) &
         call params%retrieve("io_units_mass", out_units%mass)
    if (params%is_set("io_units_charge")) &
         call params%retrieve("io_units_charge", out_units%charge)
    if (params%is_set("io_units_energy")) &
         call params%retrieve("io_units_energy", out_units%energy)
    if (params%is_set("io_units_pressure")) &
         call params%retrieve("io_units_pressure", out_units%pressure)
    if (params%is_set("io_units_force")) &
         call params%retrieve("io_units_force", out_units%force)
    if (params%is_set("io_units_velocity")) &
         call params%retrieve("io_units_velocity", out_units%velocity)
    if (params%is_set("io_units_power")) &
         call params%retrieve("io_units_power", out_units%power)
    if (params%is_set("io_units_surface_tension")) &
         call params%retrieve("io_units_surface_tension", out_units%surf_ten)
    if (params%is_set("io_units_emf")) &
         call params%retrieve("io_units_emf", out_units%emf)

    ! Check units validity
    test = convert_units(1.0_wp, out_units%length, internal_units%length)
    test = convert_units(1.0_wp, out_units%time, internal_units%time)
    test = convert_units(1.0_wp, out_units%mass, internal_units%mass)
    test = convert_units(1.0_wp, out_units%charge, internal_units%charge)
    test = convert_units(1.0_wp, out_units%energy, internal_units%energy)
    test = convert_units(1.0_wp, out_units%pressure, internal_units%pressure)
    test = convert_units(1.0_wp, out_units%force, internal_units%force)
    test = convert_units(1.0_wp, out_units%velocity, internal_units%velocity)
    test = convert_units(1.0_wp, out_units%power, internal_units%power)
    test = convert_units(1.0_wp, out_units%surf_ten, internal_units%surf_ten)
    test = convert_units(1.0_wp, out_units%emf, internal_units%emf)

    ! Initialise timestep unit
    call params%retrieve('timestep', test, .true.)
    call set_timestep(test)

    call set_out_units(out_units)

  end Subroutine read_units

  Subroutine read_forcefield(params, neigh, config, xhi, yhi, zhi, flow, vdws, electro, ewld, mpoles, cshell, met, &
       kim_Data, bond, threebody, fourbody, tersoffs)
    !!-----------------------------------------------------------------------
    !!
    !! Read forcefield
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(neighbours_type),       Intent(InOut) :: neigh
    Type(configuration_type),    Intent(In   ) :: config
    Type(flow_type),             Intent(InOut) :: flow
    Type(mpole_type),            Intent(InOut) :: mpoles
    Type(vdw_type),              Intent(InOut) :: vdws
    Type(electrostatic_type),    Intent(InOut) :: electro
    Type(core_shell_type),       Intent(InOut) :: cshell
    Type(ewald_type),            Intent(InOut) :: ewld
    Type(metal_type),            Intent(InOut) :: met
    Type(kim_type),              Intent(In   ) :: kim_data
    Type(bonds_type),            Intent(In   ) :: bond
    Type(threebody_type),        Intent(InOut) :: threebody
    Type(four_body_type),        Intent(InOut) :: fourbody
    Type(tersoff_type),          Intent(In   ) :: tersoffs
    Real(Kind=wp),               Intent(In   ) :: xhi, yhi, zhi

    Real(kind=wp), Parameter :: minimum_rcut = 1.0_wp

    Real(Kind=wp), Dimension(9) :: cell
    Real(Kind=wp), Dimension(10):: cell_properties
    Real(Kind=wp) :: cut, tol, tol1, rtmp
    Integer :: SPLINE_LIMITS

    Logical :: metals_on, vdws_on
    Character(Len=STR_LEN) :: message
    Character(Len=STR_LEN) :: option

    ! ---------------- SUBCELLING ----------------------------------------------

    call params%retrieve('subcell_threshold', neigh%pdplnc)
    if (neigh%pdplnc < 1.0_wp) call error(0, 'subcell_threshold less than minimum (1.0)')

    ! ---------------- VDW SETUP -----------------------------------------------

    call params%retrieve('vdw_method', option)
    select case(option)
    case ('tabulated')
      continue
    case ('direct')
      vdws%l_direct = .true.
    case ('off')
      vdws%no_vdw = .true.
    case ('ewald')
      ewld%vdw = .true.
    case default
      call bad_option('vdw_method', option)
    end select
    if (vdws%max_vdw <= 0) vdws%no_vdw = .true.

    vdws_on = .not. vdws%no_vdw
    metals_on = met%max_metal > 0

    call params%retrieve('vdw_mix_method', option)
    select case (option)
    case ('off')
      vdws%mixing = MIX_NULL
    case ('lorentz-bethelot')
      vdws%mixing = MIX_LORENTZ_BERTHELOT
    case ('fender-halsey')
      vdws%mixing = MIX_FENDER_HALSEY
    case ('hogervorst')
      vdws%mixing = MIX_HOGERVORST
    case ('halgren')
      vdws%mixing = MIX_HALGREN
    case ('waldman-hagler')
      vdws%mixing = MIX_WALDMAN_HAGLER
    case ('tang-tonnies')
      vdws%mixing = MIX_TANG_TOENNIES
    case ('functional')
      vdws%mixing = MIX_FUNCTIONAL
    case default
      call bad_option('vdw_mix_method', option)
    end select

    call params%retrieve('vdw_force_shift', vdws%l_force_shift)

    ! ---------------- METALS --------------------------------------------------

    call params%retrieve('metal_direct', met%l_direct)

    call params%retrieve('metal_sqrtrho', met%l_emb)

    ! ---------------- CUTOFF --------------------------------------------------

    call params%retrieve('cutoff', neigh%cutoff)
    if (neigh%cutoff < minimum_rcut) then
      neigh%cutoff = minimum_rcut
      call warning('neighbour cutoff less than minimum_cutoff (1.0 Ang), setting to minimum_cutoff', .true.)
    end if

    call params%retrieve('padding', neigh%padding)
    if (neigh%padding < 0.0_wp) then
      neigh%padding = 0.0_wp
      call warning('Bad padding value, reset to 0.0', .true.)
    end if
    flow%reset_padding = params%is_set('padding')

    if (vdws_on) then
      call params%retrieve('vdw_cutoff', rtmp)

      if (vdws%cutoff > REAL_TOL .and. vdws%cutoff < rtmp .and. Abs(vdws%cutoff - rtmp) > REAL_TOL) then
        Write (message, '(a,1p,e12.4)') 'VdW cutoff set by bounds to (Angs): ', vdws%cutoff
        Call warning(message, .true.)

      else
        vdws%cutoff = rtmp
      End If

      if (vdws%cutoff < minimum_rcut) then
        vdws%cutoff = neigh%cutoff
        call warning('VdW cutoff less than minimum cutoff, setting to global cutoff', .true.)
      end if

    else
      vdws%cutoff = 0.0_wp
    end if

    if (metals_on .and. met%rcut < REAL_TOL) then
      met%rcut=Max(met%rcut, neigh%cutoff,vdws%cutoff)
      call warning('metal_cutoff not set, setting to max of global cutoff and vdw cutoff', .true.)
    end if

    neigh%cutoff = Max(neigh%cutoff, vdws%cutoff, met%rcut, kim_data%cutoff, bond%rcut, &
         2.0_wp * tersoffs%cutoff + REAL_TOL)

    ! Sort met%rcut=neigh%cutoff if metal interactions are in play, even if
    ! they are defined by EAM since met%rcut can be /= neigh%cutoff in such
    ! instances, this can break the NLAST check in metal_ld_set_halo

    if (metals_on) met%rcut = neigh%cutoff

    ! Sort vdws%cutoff=neigh%cutoff if VDW interactions are in play

    If (vdws_on .and. vdws%cutoff > neigh%cutoff) Then
      vdws%cutoff = neigh%cutoff
      Call warning('VdW cutoff greater than global cutoff, setting to global cutoff', .true.)
    End If

    if (threebody%mxtbp > 0 .and. threebody%cutoff < REAL_TOL) then
      call warning('three body cutoff not set, setting to half of global cutoff', .true.)
      threebody%cutoff = 0.5_wp * neigh%cutoff
    end if

    If (fourbody%max_four_body > 0 .and. fourbody%cutoff < REAL_TOL) then
      call warning('four body cutoff not set, setting to half of global cutoff', .true.)
      fourbody%cutoff = 0.5_wp * neigh%cutoff
    end If

    ! ---------------- POLARISATION --------------------------------------------

    call params%retrieve('polarisation_model', option)
    select case (option)
    case ('charmm')
      if (cshell%mxshl > 0 .and. mpoles%max_mpoles > 0) mpoles%key = POLARISATION_CHARMM
      if (mpoles%max_mpoles == 0 .or. cshell%mxshl == 0) &
           call error(0, 'CHARMM polarisation selected with no mpoles or shells')

      electro%lecx = .true.
      call params%retrieve('polarisation_thole', mpoles%thole)

    case ('default')
      mpoles%key = POLARISATION_DEFAULT
    case default
      call bad_option('polarisation_model', option)
    end select

    if (cshell%mxshl > 0) then
      call params%retrieve('rlx_tol', cshell%rlx_tol(1))
      call params%retrieve('rlx_cgm_step', cshell%rlx_tol(2))
      if (cshell%rlx_tol(1) < 1.0_wp) call error(0, 'Relaxed shell CGM tolerance < 1.0')
    end if

    ! ---------------- ELECTROSTATICS ------------------------------------------

    call params%retrieve('coul_method', option)

    select case (option)
    case ('off')
      electro%no_elec = .true.
      ! reinitialise multipolar electrostatics indicators
      mpoles%max_mpoles = 0
      mpoles%max_order = 0
      mpoles%key = POLARISATION_DEFAULT
      electro%key = ELECTROSTATIC_NULL
    case ('ewald')
      electro%key = ELECTROSTATIC_EWALD
      ewld%active = .true.

    case ('dddp')
      electro%key = ELECTROSTATIC_DDDP

    case ('pairwise')

      electro%key = ELECTROSTATIC_COULOMB

    case ('force_shifted')

      electro%key = ELECTROSTATIC_COULOMB_FORCE_SHIFT

    case ('reaction_field')

      electro%key = ELECTROSTATIC_COULOMB_REACTION_FIELD

    case default

      call bad_option('coul_method', option)

    end select

    ! If it's been forcibly set by polarisation_model
    if (.not. electro%lecx) call params%retrieve('coul_extended_exclusion', electro%lecx)

    call params%retrieve('coul_dielectric_constant', electro%eps)

    if (params%is_set([Character(14) :: 'coul_damping', 'coul_precision'])) then
      call error(0, 'Both damping and precision set')

    else if (params%is_set('coul_damping')) then
      call params%retrieve('coul_damping', electro%damping)

    else if (params%is_set('coul_precision')) then
      call params%retrieve('coul_precision', rtmp)
      rtmp = Max(Min(rtmp, 0.5_wp), 1.0e-20_wp)
      tol = Sqrt(Abs(Log(rtmp * neigh%cutoff)))
      electro%damping = Sqrt(Abs(Log(rtmp * neigh%cutoff * tol))) / neigh%cutoff

    end if

    If (electro%damping > zero_plus) Then
      Call info('Fennell damping applied', .true.)
      If (neigh%cutoff < 12.0_wp) Call warning(7, neigh%cutoff, 12.0_wp, 0.0_wp)
    End If

    if (ewld%active) then

      cell = config%cell
      cut = neigh%cutoff + 1e-6_wp

      if (config%imcon == IMCON_NOPBC) then
        cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
        cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
        cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

        cell(2) = 0.0_wp
        cell(3) = 0.0_wp
        cell(4) = 0.0_wp
        cell(6) = 0.0_wp
        cell(7) = 0.0_wp
        cell(8) = 0.0_wp
      else if (config%imcon == IMCON_SLAB) then
        cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))
      End If

      call params%retrieve('ewald_nsplines', ewld%bspline%num_splines)
      ! Check splines in min
      SPLINE_LIMITS = MIN_SPLINES + mpoles%max_order
      If (ewld%bspline%num_splines < SPLINE_LIMITS) Then
        Write(message, '(a,i0,a,i0,a)') &
             "Number of bsplines (", ewld%bspline%num_splines,") less than minimum permitted (", &
             SPLINE_LIMITS, ") -- Resetting"
        Call Warning(message, .true.)
        ewld%bspline%num_splines = Max(ewld%bspline%num_splines, SPLINE_LIMITS)
      End If

      ! Check splines even
      if (Mod(ewld%bspline%num_splines, 2) /= 0) &
           ewld%bspline%num_splines = ewld%bspline%num_splines + 1

      !Check splines in max
      If (ewld%bspline%num_splines > MAX_SPLINES) Then
        Write(message, '(a,i0,a,i0,a)') &
             "Number of bsplines (", ewld%bspline%num_splines,") greater than maxiumum permitted (", &
             MAX_SPLINES, ")"
        Call error(0, message)
      End If

      Call dcell(cell,cell_properties)

      if (params%is_set([Character(15) :: 'ewald_precision', 'ewald_alpha'])) then

        call error(0, 'Cannot specify both precision and manual ewald parameters')

      else if (params%is_set('ewald_alpha')) then

        call params%retrieve('ewald_alpha', ewld%alpha)
        if (params%is_set([character(18) :: 'ewald_kvec', 'ewald_kvec_spacing'])) then

          call error(0, 'Cannot specify both explicit k-vec grid and k-vec spacing')
        else if (params%is_set('ewald_kvec')) then

          call params%retrieve('ewald_kvec', ewld%kspace%k_vec_dim_cont)
        else

          call params%retrieve('ewald_kvec_spacing', rtmp)
          ewld%kspace%k_vec_dim_cont = Nint(rtmp / cell_properties(7:9))
        end if

        ! Sanity check for ill defined ewald sum parameters 1/8*2*2*2 == 1
        tol=ewld%alpha*real(product(ewld%kspace%k_vec_dim_cont), wp)
        If (Int(tol) < 1) Call error(9)

      else

        call params%retrieve('ewald_precision', ewld%precision)

        tol = Sqrt(Abs(Log(ewld%precision*neigh%cutoff)))
        ewld%alpha = Sqrt(Abs(Log(ewld%precision*neigh%cutoff*tol)))/neigh%cutoff
        tol1 = Sqrt(-Log(ewld%precision*neigh%cutoff*(2.0_wp*tol*ewld%alpha)**2))

        ewld%kspace%k_vec_dim_cont = 2*Nint(0.25_wp + cell_properties(7:9)*ewld%alpha*tol1/pi)

      end if
    end if

  End Subroutine read_forcefield

  Subroutine read_run_parameters(params, flow, thermo, stats, read_config_indices)
    Type(parameters_hash_table), Intent(In   ) :: params
    Type(flow_type),             Intent(InOut) :: flow
    Type(thermostat_type),       Intent(InOut) :: thermo
    Type(stats_type),            Intent(InOut) :: stats
    Logical,                     Intent(  Out) :: read_config_indices

    Logical :: ltmp

    CAll params%retrieve('timestep', thermo%tstep)
    if (thermo%tstep < zero_plus) call error(0, 'Timestep too small')
    call params%retrieve('timestep_variable', thermo%lvar)

    If (thermo%key_dpd /= DPD_NULL .and. thermo%lvar) then
      thermo%lvar = .false.
      Call warning('Variable timestep unavalable in DPD themostats', .true.)
    Else If (thermo%lvar) Then
      call params%retrieve('timestep_variable_min_dist', thermo%mndis)
      call params%retrieve('timestep_variable_max_dist', thermo%mxdis)
      call params%retrieve('timestep_variable_max_delta', thermo%mxstp)

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

    call params%retrieve('stack_size', stats%mxstak)

    call params%retrieve('ignore_config_indices', read_config_indices)

    call params%retrieve('strict_checks', flow%strict)

    call params%retrieve('time_run', flow%run_steps, .true.)

    call params%retrieve('time_equilibration', flow%equil_steps)

    call params%retrieve('record_equilibration', ltmp)
    flow%equilibration = .not. ltmp

    call params%retrieve('data_dump_frequency', flow%freq_restart)

    call params%retrieve('print_frequency', flow%freq_output)

    call params%retrieve('stats_frequency', stats%intsta)

    call params%retrieve('print_topology_info', flow%print_topology)

    call params%retrieve('currents_calculate', stats%cur%on)

  end Subroutine read_run_parameters

  Subroutine read_system_parameters(params, flow, config, thermo, impa, minim, plume, cons, pmf, ttm_active)
    Type( parameters_hash_table ), Intent(In   ) :: params
    Type( flow_type ), Intent(InOut) :: flow
    Type( configuration_type ), Intent(InOut) :: config
    Type( thermostat_type ), Intent(InOut) :: thermo
    Type( impact_type ), Intent(InOut) :: impa
    Type( minimise_type ), Intent(InOut) :: minim
    Type( plumed_type ), Intent(InOut) :: plume
    Type( constraints_type ), Intent(InOut) :: cons
    Type( pmf_type ), Intent(In    ) :: pmf
    Logical, Intent(In   ) :: ttm_active

    Character(Len=STR_LEN) :: option
    Type( control_parameter ) :: param
    Real(Kind=wp), Dimension(6) :: vtmp
    Real(Kind=wp) :: rtmp
    Logical :: ltmp

    Call params%retrieve('title', option)
    config%sysname = option(1:72)

    ! ---------------- PHYSICAL PROPERTIES -------------------------------------

    call params%retrieve('temperature', thermo%temp)

    if (params%num_set([Character(22) :: 'pressure_tensor', 'pressure_hydrostatic', 'pressure_perpendicular']) > 1) then
      call error(0, 'Multiple pressure specifications')
    ! else if (params%num_set([Character(22) :: 'pressure_tensor', 'pressure_hydrostatic', 'pressure_perpendicular']) == 0) then
    !   call error(0, 'No pressure specification')
    end if

    thermo%stress = 0.0_wp
    if (params%is_set('pressure_tensor')) then
      call params%retrieve('pressure_tensor', vtmp)
      thermo%stress(1) = vtmp(1)
      thermo%stress(5) = vtmp(2)
      thermo%stress(9) = vtmp(3)
      thermo%stress(2) = vtmp(4)
      thermo%stress(4) = vtmp(4)
      thermo%stress(3) = vtmp(5)
      thermo%stress(7) = vtmp(5)
      thermo%stress(6) = vtmp(6)
      thermo%stress(8) = vtmp(6)

    else if (params%is_set('pressure_perpendicular')) then
      call params%retrieve('pressure_perpendicular', vtmp(1:3))
      thermo%stress(1:9:4) = vtmp(1:3)
    else
      call params%retrieve('pressure_hydrostatic', rtmp)
      thermo%stress(1:9:4) = rtmp
    end if

    if (.not. thermo%anisotropic_pressure .or. thermo%iso /= CONSTRAINT_NONE) then
      thermo%press = sum(thermo%stress(1:9:4)) / 3.0_wp
      thermo%stress = 0.0_wp
    end if

    call params%retrieve('fixed_com', config%l_vom)
    if (.not. config%l_vom .and. .not. ttm_active) then
      Call info('"no fixed_com" option auto-switched on - COM momentum removal will be abandoned', .true.)
      Call warning('this may lead to a build up of the COM momentum ' &
           //'and a manifestation of the "flying ice-cube" effect', .true.)
    end if

    ! ---------------- INITIALISATION ------------------------------------------

    call params%retrieve('restart', option)
    select case (option)
    case ('rescale')
      flow%restart_key = RESTART_KEY_SCALE
    case ('noscale')
      flow%restart_key = RESTART_KEY_NOSCALE
    case ('continue')
      flow%restart_key = RESTART_KEY_OLD
    case ('clean')
      flow%restart_key = RESTART_KEY_CLEAN
    case default
      call bad_option('restart', option)
    end select

    if (config%levcfg == 0 .and. flow%restart_key /= RESTART_KEY_CLEAN) then
      Call warning('CONFIG contains positions only, forcing clean restart')
      flow%restart_key = RESTART_KEY_CLEAN
    end if

    ! ---------------- EQULIBRATION --------------------------------------------

    call params%retrieve('reset_temperature_interval', thermo%freq_zero)
    thermo%l_zero = thermo%freq_zero > 0
    If (thermo%freq_zero == 0) thermo%freq_zero = flow%equil_steps + 1

    call params%retrieve('regauss_frequency', thermo%freq_tgaus)
    thermo%l_tgaus = thermo%freq_tgaus > 0
    If (thermo%freq_tgaus == 0) thermo%freq_tgaus = flow%equil_steps + 1

    call params%retrieve('rescale_frequency', thermo%freq_tscale)
    thermo%l_tscale = thermo%freq_tscale > 0
    If (thermo%freq_tscale == 0) thermo%freq_tscale = flow%equil_steps + 1

    if (params%is_set('equilibration_force_cap')) then
      flow%force_cap = .true.
      call params%retrieve('equilibration_force_cap', config%fmax)
    end if

    call params%retrieve('minimisation_criterion', option)
    minim%minimise = option /= 'off'

    call params%get('minimisation_tolerance', param)
    read(param%val, *) minim%tolerance

    call params%retrieve('minimisation_frequency', minim%freq)
    call params%retrieve('minimisation_step_length', minim%step_length)

    select case (option)
    case ('off')
      continue
    case ('force')
      minim%key = MIN_FORCE
      param%internal_units = "internal_f"
      minim%tolerance = convert_units(minim%tolerance, param%units, param%internal_units)
      If (minim%tolerance < 1.0_wp .or. minim%tolerance > 1000.0_wp) &
           minim%tolerance = 50.0_wp

    case ('energy')
      minim%key = MIN_ENERGY
      param%internal_units = "internal_e"
      minim%tolerance = convert_units(minim%tolerance, param%units, param%internal_units)
      If (minim%tolerance < zero_plus .or. minim%tolerance > 0.01_wp) &
           minim%tolerance = 0.005_wp

    case ('distance')
      minim%key = MIN_DISTANCE
      param%internal_units = "internal_l"
      minim%tolerance = convert_units(minim%tolerance, param%units, param%internal_units)
      If (minim%tolerance < REAL_TOL .or. minim%tolerance > 0.1_wp) &
           minim%tolerance = 0.005_wp

    case default
      call bad_option('minimisation_criterion', option)
    end select

    If (minim%freq == 0) minim%freq = flow%equil_steps + 1

    ! ---------------- PSEUDO-THERMOSTAT ---------------------------------------

    call params%retrieve('pseudo_thermostat_method', option)
    if (option /= 'off') then

      thermo%l_stochastic_boundaries = .true.
      Call info('pseudo thermostat attached to MD cell boundary', .true.)

      select case (option)
      case ('langevin-direct')
        thermo%key_pseudo = 0
      case ('langevin')
        thermo%key_pseudo = 1
      case ('gaussian')
        thermo%key_pseudo = 2
      case ('direct')
        thermo%key_pseudo = 3
      case default
        call bad_option('pseudo_thermostat_method', option)
      end select

      Call params%retrieve("pseudo_thermostat_width", thermo%width_pseudo)
      if (thermo%width_pseudo < 2.0_wp) then
        thermo%width_pseudo = 2.0_wp
        Call info('thermostat thickness insufficient - reset to 2 Angs', .true.)
      end if

      call params%retrieve("pseudo_thermostat_temperature", thermo%temp_pseudo)
      thermo%temp_pseudo = max(0.0_wp, thermo%temp_pseudo)

    end if

    ! --------------- CONSTRAINTS ----------------------------------------------

    If (cons%mxcons > 0 .or. pmf%mxpmf > 0) Then
      call params%retrieve('shake_max_iter', cons%max_iter_shake)
      call params%retrieve('shake_tolerance', cons%tolerance)
    End If

    ! --------------- PLUMED ---------------------------------------------------

    call params%retrieve('plumed', plume%l_plumed)

    if (plume%l_plumed) then
      call params%retrieve('plumed_input', option)
      plume%input = option(1:125)
      call params%retrieve('plumed_log', option)
      plume%logfile = option(1:125)
      call params%retrieve('plumed_precision', plume%prec)
      call params%retrieve('plumed_restart', ltmp)

      if (ltmp) then
        plume%restart = 1
      else
        plume%restart = 0
      end if
    end if

    ! --------------- IMPACT ---------------------------------------------------

    impa%active = params%is_any_set([Character(17) :: 'impact_part_index', 'impact_energy', 'impact_time'])
    if (impa%active) then

      call params%retrieve('impact_part_index', impa%imd)
      call params%retrieve('impact_time', impa%tmd)
      call params%retrieve('impact_energy', impa%emd)
      call params%retrieve('impact_direction', vtmp(1:3))
      impa%vmx = vtmp(1)
      impa%vmy = vtmp(2)
      impa%vmz = vtmp(3)
    end if

    ! --------------- EXPANSION ------------------------------------------------

    call params%retrieve('nfold', vtmp(1:3))
    config%l_exp = any(nint(vtmp(1:3)) > 1)
    if (config%l_exp) then
      config%nx = nint(vtmp(1))
      config%ny = nint(vtmp(2))
      config%nz = nint(vtmp(3))
    end if

  end Subroutine read_system_parameters

  Subroutine write_parameters(io_data, netcdf, files, neigh, config, link_cell, flow, stats, thermo, ttm, mpoles, vdws, &
       electro, cshell, ewld, met, impa, minim, plume, cons, pmf, bond, angle, dihedral, inversion, &
       msd_data, rdf, vaf, zdensity, adf, coords, defect, traj, displacement)
    Type(io_type),              Intent(In)    :: io_data
    Type(netcdf_param),         Intent(In)    :: netcdf
    Type(file_type),            Intent(In)    :: files(:)
    Type(neighbours_type),      Intent(In)    :: neigh
    Type(configuration_type),   Intent(In)    :: config
    Type(flow_type),            Intent(In)    :: flow
    Type(stats_type),           Intent(In)    :: stats
    Type(thermostat_type),      Intent(In)    :: thermo
    Type(ttm_type),             Intent(In)    :: ttm
    Type(mpole_type),           Intent(In)    :: mpoles
    Type(vdw_type),             Intent(In)    :: vdws
    Type(electrostatic_type),   Intent(In)    :: electro
    Type(core_shell_type),      Intent(In)    :: cshell
    Type(ewald_type),           Intent(In)    :: ewld
    Type(metal_type),           Intent(In)    :: met
    Type(plumed_type),          Intent(In)    :: plume
    Type(impact_type),          Intent(In)    :: impa
    Type(minimise_type),        Intent(In)    :: minim
    Type(constraints_type),     Intent(In)    :: cons
    Type(pmf_type),             Intent(In)    :: pmf
    Type(bonds_type),           Intent(In)    :: bond
    Type(angles_type),          Intent(In)    :: angle
    Type(dihedrals_type),       Intent(In)    :: dihedral
    Type(inversions_type),      Intent(In)    :: inversion
    Type(msd_type),             Intent(In)    :: msd_data
    Type(rdf_type),             Intent(InOut) :: rdf
    Type(greenkubo_type),       Intent(In)    :: vaf
    Type(z_density_type),       Intent(InOut) :: zdensity
    Type(adf_type),             Intent(In)    :: adf
    Type(coord_type),           Intent(In)    :: coords
    Type(trajectory_type),      Intent(In)    :: traj
    Type(defects_type),         Intent(In)    :: defect(:)
    Type(rsd_type),             Intent(In)    :: displacement
    Integer, Dimension(3),      Intent(In)    :: link_cell
    Character(Len=80)                         :: banner(6)

    Write (banner(1), '(a)') ''
    Write (banner(2), '(a)') Repeat('*', 80)
    Write (banner(3), '(a4,a72,a4)') '*** ', 'title:'//Repeat(' ', 66), ' ***'
    Write (banner(4), '(a4,a72,a4)') '*** ', config%sysname, ' ***'
    Write (banner(5), '(a)') Repeat('*', 80)
    Write (banner(6), '(a)') ''
    Call info(banner, 6, .true.)

    if (check_print_level(1)) Call write_io(io_data, netcdf, files)
    if (check_print_level(1)) Call write_system_parameters(flow, config, stats, thermo, impa, minim, plume, cons, pmf)
    if (check_print_level(1)) Call write_ensemble(thermo)
    If (check_print_level(1)) Call write_forcefield(link_cell, neigh, vdws, electro, ewld, mpoles, cshell, met)
    if (ttm%l_ttm .and. check_print_level(1)) Call write_ttm(thermo, ttm)
    if (check_print_level(1)) Call write_bond_analysis(stats, config, flow, bond, angle, dihedral, inversion)
    if (check_print_level(1)) &
         Call write_structure_analysis(stats, msd_data, rdf, vaf, zdensity, adf, coords, traj, defect, displacement)
    Call info('', .true.)

  end Subroutine write_parameters

  Subroutine write_io(io_data, netcdf, files)
    Type(io_type), Intent(In) :: io_data
    Type(netcdf_param), Intent(In) :: netcdf
    Type(file_type), Intent(In) :: files(:)

    Character(len=256) :: message
    !    Integer :: io_read, io_write, procs_write, procs_read, batch_write, batch_read

    Call info('', .true.)
    Call info('I/O Parameters: ', .true.)

    select case (io_data%method_read)
    Case (IO_READ_MPIIO)
      Call info('  I/O read method: parallel by using MPI-I/O',.true.)

    Case (IO_READ_DIRECT)
      Call info('  I/O read method: parallel by using direct access',.true.)

    Case (IO_READ_NETCDF)
      Call info('  I/O read method: parallel by using netCDF',.true.)

    Case (IO_READ_MASTER)
      Call info('  I/O read method: serial by using a single master process',.true.)

    End Select

    Write(message,'(a,i0.1)') '  I/O readers set to ', io_data%n_io_procs_read
    Call info(message,.true.)

    Select Case (io_data%method_write)
    Case (IO_READ_MPIIO, IO_READ_DIRECT, IO_READ_NETCDF)
      Write(message,'(a,i0.1)') '  I/O read batch size (assumed) ', io_data%batch_size_read
      Call info(message,.true., level=3)
    end Select

    Write(message,'(a,i0.1)') '  I/O read buffer size set to ', io_data%buffer_size_read
    Call info(message,.true., level=3)

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_MPIIO)
      Call info('  I/O write method: parallel by using MPI-I/O',.true.)

    Case (IO_WRITE_SORTED_DIRECT, IO_WRITE_UNSORTED_DIRECT)
      Call info('  I/O write method: parallel by using direct access',.true.)
      Call warning('  in parallel this I/O write method has portability issues',.true.)

    Case (IO_WRITE_SORTED_NETCDF)
#ifdef NETCDF
      Select Case(netcdf%ncp)
      Case( sp )
        Call info('  I/O write method: parallel by using netCDF in the amber-like/32-bit format',.true.)
      Case ( dp )
        Call info('  I/O write method: parallel by using netCDF in 64-bit format',.true.)
      end Select
#else
      call error(0, 'NetCDF not compiled')
#endif

    Case (IO_WRITE_SORTED_MASTER, IO_WRITE_UNSORTED_MASTER)
      Call info('  I/O write method: serial by using a single master process',.true.)

    End Select

    Select Case (io_data%method_write)
    Case(IO_WRITE_SORTED_MASTER, IO_WRITE_SORTED_NETCDF, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT)
      Call info('  I/O write type: data sorting ON',.true., level=3)

    Case(IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT)
      Call info('  I/O write type: data sorting OFF',.true.)

    End Select

    Write(message,'(a,i0.1)') '  I/O writers set to ', io_data%n_io_procs_write
    Call info(message,.true.)

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_SORTED_DIRECT, &
         IO_WRITE_UNSORTED_DIRECT, IO_WRITE_SORTED_NETCDF)
      Write(message,'(a,i10)') '  I/O write batch size (assumed) ', io_data%batch_size_write
      Call info(message,.true., level=3)
    end Select

    Write(message,'(a,i10)') '  I/O write buffer size set to ', io_data%buffer_size_write
    Call info(message,.true., level=3)


    If (io_data%global_error_check) Then
      Call info('  I/O parallel error checking ON',.true.)
    Else
      Call info('  I/O parallel error checking OFF',.true., level=3)
    End If

    Call info('', .true.)
    Call info('File outputs: ',.true.)
    if (files(FILE_OUTPUT)%filename /= "OUTPUT") Call info('  OUTPUT file is '//files(FILE_OUTPUT)%filename,.true.)
    if (files(FILE_CONFIG)%filename /= "CONFIG") Call info('  CONFIG file is '//files(FILE_CONFIG)%filename,.true.)
    if (files(FILE_FIELD)%filename /= "FIELD") Call info('  FIELD file is '//files(FILE_FIELD)%filename,.true.)
    if (files(FILE_STATS)%filename /= "STATIS") Call info('  STATIS file is '//files(FILE_STATS)%filename,.true.)
    if (files(FILE_HISTORY)%filename /= "HISTORY") Call info('  HISTORY file is '//files(FILE_HISTORY)%filename,.true.)
    if (files(FILE_HISTORF)%filename /= "HISTORF") Call info('  HISTORF file is '//files(FILE_HISTORF)%filename,.true.)
    if (files(FILE_REVIVE)%filename /= "REVIVE") Call info('  REVIVE file is '//files(FILE_REVIVE)%filename,.true.)
    if (files(FILE_REVCON)%filename /= "REVCON") Call info('  REVCON file is '//files(FILE_REVCON)%filename,.true.)
    if (files(FILE_REVOLD)%filename /= "REVOLD") Call info('  REVOLD file is '//files(FILE_REVOLD)%filename,.true.)
    if (files(FILE_RDF)%filename /= 'RDFDAT') Call info('  RDF file is '//files(FILE_RDF)%filename,.true.)
    if (files(FILE_MSD)%filename /= 'MSDTMP') Call info('  MSD file is '//files(FILE_MSD)%filename,.true.)
    if (files(FILE_TABBND)%filename /= 'TABBND') Call info('  TABBND file is '//files(FILE_TABBND)%filename,.true.)
    if (files(FILE_TABANG)%filename /= 'TABANG') Call info('  TABANG file is '//files(FILE_TABANG)%filename,.true.)
    if (files(FILE_TABDIH)%filename /= 'TABDIH') Call info('  TABDIH file is '//files(FILE_TABDIH)%filename,.true.)
    if (files(FILE_TABINV)%filename /= 'TABINV') Call info('  TABINV file is '//files(FILE_TABINV)%filename,.true.)
    if (files(FILE_TABVDW)%filename /= 'TABVDW') Call info('  TABVDW file is '//files(FILE_TABVDW)%filename,.true.)
    if (files(FILE_TABEAM)%filename /= 'TABEAM') Call info('  TABEAM file is '//files(FILE_TABEAM)%filename,.true.)

  end Subroutine write_io

  Subroutine write_bond_analysis(stats, config, flow, bond, angle, dihedral, inversion)
    Type(stats_type),         Intent(In) :: stats
    Type(configuration_type), Intent(In) :: config
    Type(flow_type),          Intent(In) :: flow
    Type(bonds_type),         Intent(In) :: bond
    Type(angles_type),        Intent(In) :: angle
    Type(dihedrals_type),     Intent(In) :: dihedral
    Type(inversions_type),    Intent(In) :: inversion

    Character(Len=256) :: messages(4)

    Call info('', .true.)
    If (.not. any([flow%analyse_bond, flow%analyse_ang, flow%analyse_dih, flow%analyse_inv])) Then
      Call info('No intramolecular distribution collection requested', .true., level=3)
    Else
      Call info('Intramolecular distribution collection requested for:', .true.)

      if (flow%analyse_bond) then
        Write(messages(1), '(a)') '  Bonds:'
        Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_bond
        Write(messages(3), '(a, i0.1)') '  -- Num samples  (points): ', bond%bin_pdf
        Write(messages(4), '(a, 1p g12.5e2)') '  -- Cutoff         (Angs): ', bond%rcut
        Call info(messages, 4, .true.)
      end if

      if (flow%analyse_ang) then
        Write(messages(1), '(a)') '  Angles:'
        Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_angle
        Write(messages(3), '(a, i0.1)') '  -- Num samples  (points): ', angle%bin_adf
        Call info(messages, 3, .true.)
      end if

      if (flow%analyse_dih) then
        Write(messages(1), '(a)') '  Dihedrals:'
        Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_dihedral
        Write(messages(3), '(a, i0.1)') '  -- Num samples  (points): ', dihedral%bin_adf
        Call info(messages, 3, .true.)
      end if

      if (flow%analyse_inv) then
        Write(messages(1), '(a)') '  Inversions:'
        Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', flow%freq_inversion
        Write(messages(3), '(a, i0.1)') '  -- Num samples  (points): ', inversion%bin_adf
        Call info(messages, 3, .true.)
      end if

    End If

    If (stats%lpana) Then
      Call info('Probability distribution analysis printing requested', .true.)
    Else
      Call info('No probability distribution analysis printing requested', .true., level=3)
    End If

  end Subroutine write_bond_analysis

  Subroutine write_structure_analysis(stats, msd_data, rdf, vaf, zdensity, adf, coords, traj, defect, displacement)
    Type(stats_type),         Intent(In)    :: stats
    Type(msd_type),           Intent(In)    :: msd_data
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(greenkubo_type),     Intent(In)    :: vaf
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(adf_type),           Intent(In)    :: adf
    Type(coord_type),         Intent(In)    :: coords
    Type(trajectory_type),    Intent(In)    :: traj
    Type(defects_type),       Intent(In)    :: defect(:)
    Type(rsd_type),           Intent(In)    :: displacement

    Character(Len=256) :: messages(4)

    Call info('', .true.)

    If (rdf%l_collect) Then
      Write(messages(1), '(a)') 'RDF collection requested:'
      Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', rdf%freq
      Write(messages(3), '(a, 1p g12.5e2)') '  -- Bin size       (Angs): ', rdf%rbin
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
      Write(messages(1), '(a)') 'Z-density profiles requested:'
      Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', zdensity%frequency
      Write(messages(3), '(a, 1p g12.5e2)') '  -- Bin size       (Angs): ', zdensity%bin_width
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
      Write(messages(1), '(a)') 'VAF profiles requested:'
      Write(messages(2), '(a, i0.1)') '  -- Collect every (steps): ', vaf%freq
      Write(messages(3), '(a, 1p g12.5e2)') '  -- Bin size       (Angs): ', vaf%binsize
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
      Write(messages(1), '(a)') 'MSDTMP profiles requested:'
      Write(messages(2), '(a,i0.1)') '  -- File start: ', msd_data%start
      Write(messages(3), '(a,i0.1)') '  -- File interval: ', msd_data%freq
      Call info(messages, 3, .true.)
    Else
      Call info('No MSD analysis requested', .true., level=3)
    End If

    if (traj%ltraj) then
      Write (messages(1), '(a)') 'Trajectory recording requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', traj%start
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', traj%freq
      Select Case (traj%key)
      Case(TRAJ_KEY_COORD)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COORD'
      Case(TRAJ_KEY_COORD_VEL)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COORD, VEL'
      Case(TRAJ_KEY_COORD_VEL_FORCE)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COORD, VEL, FORCE'
      Case(TRAJ_KEY_COMPRESSED)
        Write (messages(4), '(a)') '  -- Trajectory file detail: COMPRESSED'
      end Select

      Call info(messages, 4, .true.)
    Else
      Call info('No trajectory recording requested', .true., level=3)
    end if

    if (defect(1)%ldef) then
      Write (messages(1), '(a)') 'Defects analysis requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', defect(1)%nsdef
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', defect(1)%isdef
      Write (messages(4), '(a,1p g12.5e2)') '  -- Distance condition (Angs): ', defect(1)%rdef
      Call info(messages, 4, .true.)
      if (defect(2)%ldef) Call info('DEFECTS1 file option: ON', .true.)
    Else
      Call info('No defects analysis requested', .true., level=3)
    end if

    if (displacement%lrsd) then
      Write (messages(1), '(a)') 'Displacements analysis requested:'
      Write (messages(2), '(a,i0.1)') '  -- File start: ', displacement%nsrsd
      Write (messages(3), '(a,i0.1)') '  -- File interval: ', displacement%isrsd
      Write (messages(4), '(a,1p g12.5e2)') '  -- Distance condition (Angs): ', displacement%rrsd
      Call info(messages, 4, .true.)
    Else
      Call info('No displacements analysis requested', .true., level=3)
    end if

  End Subroutine write_structure_analysis

  Subroutine write_system_parameters(flow, config, stats, thermo, impa, minim, plume, cons, pmf)
    Type(flow_type),          Intent(In) :: flow
    Type(stats_type),         Intent(In) :: stats
    Type(configuration_type), Intent(In) :: config
    Type(thermostat_type),    Intent(In) :: thermo
    Type(impact_type),        Intent(In) :: impa
    Type(minimise_type),      Intent(In) :: minim
    Type(plumed_type),        Intent(In) :: plume
    Type(constraints_type),   Intent(In) :: cons
    Type(pmf_type),           Intent(In) :: pmf

    Character(Len=256) :: message, messages(5)
    Integer :: itmp

    Call info('', .true.)
    Call info("System parameters: ", .true.)

    If (.not. thermo%lvar) then
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
      end If
    end If

    Write (message, '(a,i10)') '  Run Duration (steps): ', flow%run_steps
    if (flow%run_steps > 0) Call info(message, .true.)

    if (flow%equil_Steps > 0) then
      Write (message, '(a,i10)') '  Equilibration period (steps): ', flow%equil_steps
      Call info(message, .true.)

      if (.not. flow%equilibration) Call info('  -- Equilibration included in overall averages', .true.)

      if (thermo%freq_tgaus > 0) then
        Call info('  -- Temperature regaussing: ON', .true.)
        Write (message, '(a,i10)') '  -- Temperature regaussing interval (steps): ', thermo%freq_tgaus
        Call info(message, .true.)
      end if

      if (thermo%freq_tscale > 0) then
        Call info('  -- Temperature scaling: ON', .true.)
        Write (message, '(a,i10)') '  -- Temperature scaling interval (steps): ', thermo%freq_tscale
        Call info(message, .true.)
      end if

      if (flow%force_cap) then
        Write (messages(1), '(a)') '  -- Force capping: ON'
        Write (messages(2), '(a,1p g12.5e2)') '  -- Force capping limit (kT/Angs): ', config%fmax
        Call info(messages, 2, .true.)
      end if

      if (minim%minimise) then
        Write (messages(1), '(a)') '  -- Minimisation: ON'
        Select Case(minim%key)
        Case(MIN_FORCE)
          Write (messages(2), '(a,a8)')     '   -- Minimisation criterion:         ', "Force"
        Case(MIN_ENERGY)
          Write (messages(2), '(a,a8)')     '   -- Minimisation criterion:         ', "Energy"
        Case(MIN_DISTANCE)
          Write (messages(2), '(a,a8)')     '   -- Minimisation criterion:         ', "Distance"
        end Select

        Write (messages(3), '(a,i10)')      '   -- Minimisation frequency (steps): ', minim%freq
        Write (messages(4), '(a,1p g12.5e2)') '   -- Minimisation tolerance:         ', minim%tolerance
        Write (messages(5), '(a,1p g12.5e2)') '   -- Minimisation CGM step:          ', minim%step_length

        Call info(messages, 5, .true.)
      end if
    end if

    if (stats%mxstak > 0) then
      Write(message, '(a,i10)') '  Rolling averages length (steps): ', stats%mxstak
      Call info(message, .true.)
    end if

    if (Any([flow%freq_restart > 0, flow%freq_output > 0, stats%intsta > 0, .not. flow%print_topology])) then
      Call info('  File output info:', .true.)
    end if

    if (flow%freq_restart > 0) Then
      Write (message, '(a,i10)') '  -- Restart dumping interval (steps): ', flow%freq_restart
      Call info(message, .true.)
    end if

    if (flow%freq_output > 0) then
      Write (message, '(a,i10)') '  -- Data printing interval (steps): ', flow%freq_output
      Call info(message, .true.)
    end if

    if (stats%intsta > 0) then
      Write (message, '(a,i10)') '  -- Statistics file interval (steps): ', stats%intsta
      Call info(message, .true.)
    end if

    if (.not. flow%print_topology) then
      Call info('  -- Not printing extended FIELD topology in OUTPUT', .true., level=3)
    else
      Call info('  -- Printing extended FIELD topology in OUTPUT: ON', .true.)
    end if

    if (stats%cur%on) Call info("  Computing currents: ON", .true.)

    select case (flow%restart_key)
    case (RESTART_KEY_SCALE)
      Call info('  Scaled restart requested (starting a new simulation)', .true.)
    case (RESTART_KEY_NOSCALE)
      Call info('  Unscaled restart requested (starting a new simulation)', .true.)
    case (RESTART_KEY_OLD)
      Call info('  Restart requested (continuing an old simulation)', .true.)
      Call warning('  Timestep from REVOLD overides specification in CONTROL', .true.)
    case (RESTART_KEY_CLEAN)
      if (config%levcfg /= 0) then
        Call info('  Clean start requested, discarding CONFIG velocities', .true.)
      end if
    end select

    Write (message, '(a,1p g12.5e2)') '  Simulation temperature (K):  ', thermo%temp
    call info(message, .true.)

    if (thermo%l_stochastic_boundaries) then
      Call info('  Pseudo-thermostat attached to MD cell boundary: ON', .true.)

      select case (thermo%key_pseudo)
      case (0)
        Call info('  -- Thermostat control: Langevin + direct temperature scaling', .true.)
      case (1)
        Call info('  -- Thermostat control: Langevin temperature scaling', .true.)
      case (2)
        Call info('  -- Thermostat control: gaussian temperature scaling', .true.)
      case (3)
        Call info('  -- Thermostat control: direct temperature scaling', .true.)
      end select

      Write (message, '(a,1p g12.5e2)') '  -- Thermostat thickness (Angs): ', thermo%width_pseudo
      Call info(message, .true.)
      Write (message, '(a,1p g12.5e2)') '  -- Thermostat temperature (K): ', thermo%temp_pseudo
      Call info(message, .true.)
    end if

    If (thermo%press > 0.0_wp) then
      write(message, '(3A,1p g12.5e2)') '  Simulation pressure (','katms','):', convert_units(thermo%press, 'internal_p', 'katm')
      Call info(message, .true.)
    Else If (any(thermo%stress > 0.0_wp)) Then
      write (messages(1), '(3A)') '  Simulation pressure (','katms','):'
      Write (messages(2), '(2X, 3(1p g12.5e2))') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp = 1,3)
      Write (messages(3), '(2X, 3(1p g12.5e2))') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp = 4,6)
      Write (messages(4), '(2X, 3(1p g12.5e2))') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp = 7,9)
      Call info(messages, 4, .true.)
    end If

    If (cons%mxcons > 0 .or. pmf%mxpmf > 0) Then
      Write (messages(1), '(a,i0.1)') '  Iterations for SHAKE/RATTLE: ', cons%max_iter_shake
      Write (messages(2), '(a,1p g12.5e2)') '  Tolerance for SHAKE/RATTLE (Angs): ', cons%tolerance
      Call info(messages, 2, .true.)
    end If

    if (config%l_exp) then
      Write (message, '(a,9x,3i5)') '  System expansion: ', config%nx, config%ny, config%nz
      Call info(message, .true.)
    end if

    if (config%dvar > 1.0_wp) then
      Write (message, '(a,9x,1p g12.5e2)') '  Permitted density variance (%): ', &
           convert_units(config%dvar - 1.0_wp, '', '%')
      Call info(message, .true.)
    end if

    if (impa%active) then
      Call info('  Impact calculation: ON', .true.)
      Write (messages(1), '(a)') ''
      Write (messages(2), '(a,i10)') '  -- Particle (index): ', impa%imd
      Write (messages(3), '(a,i10)') '  -- Timestep (steps): ', impa%tmd
      Write (messages(4), '(a,1p g12.5e2)') '  -- Energy   (keV): ', impa%emd
      Write (messages(5), '(a,3(1p g12.5e2))') '  -- v-r(x,y,z): ', impa%vmx, impa%vmy, impa%vmz

      Call info(messages, 5, .true.)

    end if

  end Subroutine write_system_parameters

  Subroutine write_forcefield(link_cell, neigh, vdws, electro, ewld, mpoles, cshell, met)
    Type(neighbours_type),       Intent(In) :: neigh
    Type(mpole_type),            Intent(In) :: mpoles
    Type(vdw_type),              Intent(In) :: vdws
    Type(electrostatic_type),    Intent(In) :: electro
    Type(core_shell_type),       Intent(In) :: cshell
    Type(ewald_type),            Intent(In) :: ewld
    Type(metal_type),            Intent(In) :: met

    Integer, Dimension(3),       Intent(In) :: link_cell
    Character(Len=256) :: message, messages(4)

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
    call info(message, .true.)

    ! ---------------- VDW SETUP -----------------------------------------------

    if (vdws%l_direct) then
      Call info('  VdWs: Direct calculation', .true.)
    else if (vdws%no_vdw) then
      Call info('  VdWs: Disabled', .true.)
      ! else if (ewld%vdw) then
      !   Call info('VdWs - Ewald', .true.)
    else
      Call info('  VdWs: Tabulated', .true.)
    end if

    if (.not. vdws%no_vdw) then
      Write (message, '(a,1p g12.5e2)') '  -- VdW cutoff (Angs): ', vdws%cutoff
      call info(message, .true.)

      if (vdws%mixing /= MIX_NULL) then
        Call info('  -- Vdw cross terms mixing opted (for undefined mixed potentials)', .true.)
        Call info('    mixing is limited to potentials of the same type only', .true.)
        Call info('    mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)', .true.)

        Select Case (vdws%mixing)
        case (MIX_LORENTZ_BERTHELOT)
          Call info('  -- Mixing scheme: Lorentz-Berthelot - e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2', .true.)

        case (MIX_FENDER_HALSEY)
          Call info('  -- Mixing scheme: Fender-Halsey - e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2', .true.)

        case (MIX_HOGERVORST)
          Call info('  -- Mixing scheme: Hogervorst (good hope) - ' &
               //'e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)', .true.)

        case (MIX_HALGREN)
          Call info('  -- Mixing scheme: Halgren HHG - ' &
               //'e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)', .true.)

        case (MIX_WALDMAN_HAGLER)
          Call info('  -- Mixing scheme: Waldman-Hagler - ' &
               //'e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ;s_ij=[(s_i^6+s_j^6)/2]^(1/6)', .true.)

        case (MIX_TANG_TOENNIES)
          Call info('  -- Mixing scheme: Tang-Toennies - ' &
               //' e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13', .true.)
          Call info(Repeat(' ', 43)//'s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)', .true.)

        case (MIX_FUNCTIONAL)
          Call info('  -- Mixing scheme: Functional - ' &
               //'e_ij=3 * (e_i*e_j)^(1/2) * (s_i*s_j)^3 / ' &
               //'SUM_L=0^2{[(s_i^3+s_j^3)^2 / (4*(s_i*s_j)^L)]^(6/(6-2L))}', .true.)
          Call info(Repeat(' ', 40)//'s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}', .true.)
        end Select
      end if

      if (vdws%l_force_shift) Call info('  -- VdW force-shifting: ON', .true.)
    end if

    ! ---------------- METALS --------------------------------------------------

    if (met%l_direct) Call info('  Metal direct: ON', .true.)
    if (met%l_emb) Call info('  Metal sqrtrho: ON', .true.)

    ! ---------------- POLARISATION --------------------------------------------

    if (mpoles%max_mpoles > 0) then
      Select Case(mpoles%key)
      Case(POLARISATION_CHARMM)
        Call info('  Polarisation method: CHARMM')
      Case(POLARISATION_DEFAULT)
        Call info('  Polarisation method: Default')
      end Select
    end if

    if (cshell%mxshl > 0) then
      Write (message, '(a,1p g12.5e2)') '  -- Relaxed shell model CGM tolerance (Force): ', cshell%rlx_tol(1)
      Call info(message, .true.)
      if (cshell%rlx_tol(2) > 0.0_wp) then
        Write (message, '(a,1p g12.5e2)') '  -- Relaxed shell model CGM step (Angs): ', cshell%rlx_tol(2)
        Call info(message, .true.)
      end if
    end if

    ! ---------------- ELECTROSTATICS ------------------------------------------

    Select Case (electro%key)
    case (ELECTROSTATIC_NULL)
      Call info('  Electrostatics: Disabled', .true.)
    case (ELECTROSTATIC_EWALD)

      Call info('  Electrostatics: Smooth Particle Mesh Ewald', .true.)

      if (ewld%precision > 0.0_wp) then
        Write(message, '(a,1p g12.5e2)') '  -- Ewald sum precision: ', ewld%precision
        Call info(message, .true.)
      end if

      Write (messages(1), '(a,1p g12.5e2)') '  -- Ewald convergence parameter (Angs^-1): ', ewld%alpha
      Write (messages(2), '(a,3i5)') '  -- Ewald kmax1 kmax2 kmax3   (x2): ', ewld%kspace%k_vec_dim_cont
      if (any(ewld%kspace%k_vec_dim /= ewld%kspace%k_vec_dim_cont)) then
        Write (messages(3), '(a,3i5)') '  -- DaFT adjusted kmax values (x2): ', ewld%kspace%k_vec_dim
        Write (messages(4), '(a,i0.1)') '  -- B-spline interpolation order: ', ewld%bspline%num_splines
        Call info(messages, 4, .true.)
      Else
        Write (messages(3), '(a,i0.1)') '  -- B-spline interpolation order: ', ewld%bspline%num_splines
        Call info(messages, 3, .true.)
      end If

    case (ELECTROSTATIC_DDDP)
      Call info('  Electrostatics: Distance Dependent Dielectric', .true.)

    case (ELECTROSTATIC_COULOMB)

      Call info('  Electrostatics: Coulombic Potential', .true.)

    case (ELECTROSTATIC_COULOMB_FORCE_SHIFT)
      Call info('  Electrostatics: Force-Shifted Coulombic Potential', .true.)

    case (ELECTROSTATIC_COULOMB_REACTION_FIELD)
      Call info('  Electrostatics: Reaction Field', .true.)
    end select

    if (electro%key /= ELECTROSTATIC_NULL) then
      if (abs(electro%eps - 1.0_wp) > zero_plus) then
        Write (message, '(a,1p g12.5e2)') '  -- Relative dielectric constant: ', electro%eps
        Call info(message, .true.)
      end if

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

    end if

  end Subroutine write_forcefield

  Subroutine write_ensemble(thermo)
    Type( thermostat_type ), Intent(In) :: thermo
    Character(Len=STR_LEN) :: message
    Character(Len=STR_LEN), dimension(4) :: messages

    Call info('', .true.)
    Call info('Thermostat details:', .true.)

    ensembles:select case(thermo%ensemble)
    case (ENS_NVE)
      Select Case (thermo%key_dpd)
      case (DPD_NULL)
        Call info('  Ensemble: NVE (Microcanonical)',.true.)
        exit ensembles
      Case(DPD_FIRST_ORDER)

        Call info('  Ensemble: NVT dpd (Dissipative Particle Dynamics)',.true.)
        Call info("  Ensemble type: Shardlow's first order splitting (S1)",.true.)

      Case (DPD_SECOND_ORDER)

        Call info('  Ensemble: NVT dpd (Dissipative Particle Dynamics)',.true.)
        Call info("  Ensemble type: Shardlow's first order splitting (S2)",.true.)

      end select

      If (thermo%gamdpd(0) > zero_plus) Then
        Write(message,'(a,1p g12.5e2)') '  -- Drag coefficient (Dalton/ps): ', thermo%gamdpd(0)
        Call info(message,.true.)
      End If

    case (ENS_NVT_EVANS)

      Call info('  Ensemble: NVT Evans (Isokinetic)',.true.)
      Call info('  Gaussian temperature constraints in use',.true.)

    case (ENS_NVT_LANGEVIN)

      Call info('  Ensemble: NVT Langevin (Stochastic Dynamics)',.true.)
      Write(message,'(a,1p g12.5e2)') '  -- Thermostat friction (ps^-1): ', thermo%chi
      Call info(message,.true.)

    case (ENS_NVT_ANDERSON)

      Write(messages(1),'(a)') '  Ensemble: NVT Andersen'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Softness (dimensionless): ',thermo%soft

      Call info(messages,3,.true.)

    case (ENS_NVT_BERENDSEN)

      Call info('  Ensemble: NVT Berendsen',.true.)
      Write(message,'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Call info(message,.true.)
      Call warning('If you plan to use the Berendsen thermostat, please read https://doi.org/10.1021/acs.jctc.8b00446', .true.)

    case (ENS_NVT_NOSE_HOOVER)

      Call info('  Ensemble: NVT Nose-Hoover',.true.)
      Write(message,'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Call info(message,.true.)

    Case (ENS_NVT_GENTLE)

      Write(messages(1),'(a)') '  Ensemble: NVT gentle stochastic thermostat'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Friction on thermostat  (ps^-1): ',thermo%gama
      Call info(messages,3,.true.)

    Case (ENS_NVT_LANGEVIN_INHOMO)

      Write (messages(1), '(a)') '  Ensemble: NVT inhomogeneous Langevin (Stochastic Dynamics)'
      Write (messages(2), '(a,1p g12.5e2)') '  -- e-phonon friction (ps^-1): ', thermo%chi_ep
      Write (messages(3), '(a,1p g12.5e2)') '  -- e-stopping friction (ps^-1): ', thermo%chi_es
      Write (messages(4), '(a,1p g12.5e2)') '  -- e-stopping velocity (ang ps^-1): ', thermo%vel_es2
      Call info(messages, 4, .true.)


    Case (ENS_NPT_LANGEVIN)

      Write(messages(1),'(a)') '  Ensemble: NPT isotropic Langevin (Stochastic Dynamics)'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat friction (ps^-1): ',thermo%chi
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat friction (ps^-1): ',thermo%tai
      Call info(messages,3,.true.)

    Case (ENS_NPT_BERENDSEN)


      Write(messages(1),'(a)') '  Ensemble: NPT isotropic Berendsen'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ',thermo%tau_p
      Call info(messages,3,.true.)

    Case (ENS_NPT_NOSE_HOOVER)


      Write(messages(1),'(a)') '  Ensemble: NPT isotropic Nose-Hoover (Melchionna)'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ',thermo%tau_p
      Call info(messages,3,.true.)

    Case (ENS_NPT_MTK)


      Write(messages(1),'(a)') '  Ensemble: NPT isotropic Martyna-Tuckerman-Klein'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ',thermo%tau_p
      Call info(messages,3,.true.)


    Case (ENS_NPT_LANGEVIN_ANISO)

      Write(messages(1),'(a)') '  Ensemble: NPT anisotropic Langevin (Stochastic Dynamics)'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat friction (ps^-1): ',thermo%chi
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat friction (ps^-1): ',thermo%tai
      Call info(messages,3,.true.)

    Case (ENS_NPT_BERENDSEN_ANISO)

      Write(messages(1),'(a)') '  Ensemble: NPT anisotropic Berendsen'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ',thermo%tau_p
      Call info(messages,3,.true.)

    Case (ENS_NPT_NOSE_HOOVER_ANISO)

      Write(messages(1),'(a)') '  Ensemble: NPT anisotropic Nose-Hoover (Melchionna)'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ',thermo%tau_p
      Call info(messages,3,.true.)

    Case (ENS_NPT_MTK_ANISO)

      Write(messages(1),'(a)') '  Ensemble: NPT anisotropic Martyna-Tuckerman-Klein'
      Write(messages(2),'(a,1p g12.5e2)') '  -- Thermostat relaxation time (ps): ',thermo%tau_t
      Write(messages(3),'(a,1p g12.5e2)') '  -- Barostat relaxation time (ps): ',thermo%tau_p
      Call info(messages,3,.true.)

    end select ensembles

    ! Semi isotropic ensembles

    select case (thermo%iso)
    case (CONSTRAINT_NONE)
      continue
    case (CONSTRAINT_SURFACE_AREA)
      Call info('  Semi-isotropic barostat: constant normal pressure (Pn) &',.true.)
      Call info('         (N-Pn-A-T)      : constant surface area (A)',.true.)
    case (CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC)

      if (thermo%tension > 0.0_wp) then
        Write(messages(1),'(a)') '  Semi-isotropic barostat: constant normal pressure (Pn) &'
        Write(messages(2),'(a)') '       (N-Pn-gamma-T)    : constant surface tension (gamma)'
        Write(messages(3),'(a,e11.4)') &
             '  -- Simulation surface tension ('//'dyn/cm'//'): ', convert_units(thermo%tension, 'internal_f/internal_l', 'dyn/cm')
        Call info(messages,3,.true.)
      else
        Call info('  Semi-isotropic barostat: orthorhombic MD cell constraints',.true.)
      end if

      if (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) then
        Call info('  Semi-isotropic barostat: semi-orthorhombic MD cell constraints',.true.)
      end if

    end select

    If (Any(thermo%iso == [CONSTRAINT_SURFACE_AREA,CONSTRAINT_SURFACE_TENSION])) Then
      Call info('  -- Semi-isotropic ensembles are only correct for infinite ', .true.)
      Call info('       interfaces placed perpendicularly to the z axis', .true.)
    End If

  end Subroutine write_ensemble

  Subroutine write_ttm(thermo, ttm)
    Type(thermostat_type),      Intent(In) :: thermo
    Type( ttm_type ), Intent(In) :: ttm
    Character(Len=256) :: message, messages(3)

    Call info('', .true.)
    Call info('TTM Parameters: ', .true.)

    Write (messages(1), '(a,3(1x,i8.1))') '  Ionic temperature grid size (x,y,z): ', ttm%ntsys(1:3)
    Write (messages(2), '(a,3(1x,1p g12.5e2))') '  Temperature grid size (x,y,z): ', ttm%delx, ttm%dely, ttm%delz
    Write (messages(3), '(a,1p g12.5e2)') '  Average number of atoms/cell: ', ttm%sysrho * ttm%volume
    Call info(messages, 3, .true.)
    Write (message, '(a,3(1x,i8.1))') '  Electronic temperature grid size (x,y,z): ', &
         ttm%eltsys(1:3)
    Call info(message, .true.)

    if (ttm%ismetal) then
      Call info('  Electronic subsystem: Metal (thermal conductivity required)', .true.)
    else
      Call info('  Electronic subsystem: Non-Metal (thermal diffusivity required)', .true.)
    end if

    if (ttm%ttmdyndens) then
      Call info('  Dynamic atomic density: ON', .true.)
    else
      Write (message, '(a,1p g12.5e2)') '  Atomic density (A^-3): ', ttm%cellrho
      Call info(message, .true.)
    end if

    select case (ttm%cetype)
    case (0) !Constant
      Call info('  Electronic specific heat capacity: Constant', .true.)
      Write (message, '(a,1p g12.5e2)') '  -- Electronic S.H.C. (kB/atom): ', ttm%Ce0 / ttm%cellrho
      Call info(message, .true.)

    case (4) !Constant
      Call info('  Electronic specific heat capacity: Constant', .true.)
      Write (message, '(a,1p g12.5e2)') '  -- Electronic S.H.C. (kB/atom): ', ttm%Ce0
      Call info(message, .true.)

    case (1) !tanh
      Write (messages(1), '(a)') '  Electronic specific heat capacity: Hyperbolic tangent'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Constant term A (kB/atom): ', ttm%sh_A / ttm%cellrho
      Write (messages(3), '(a,1p g12.5e2)') '  -- Temperature term B (K^-1): ', ttm%sh_B * 1.0e4_wp
      Call info(messages, 3, .true.)

    case (5) !tanh
      Write (messages(1), '(a)') '  Electronic specific heat capacity: Hyperbolic tangent'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Constant term A (kB/atom): ', ttm%sh_A
      Write (messages(3), '(a,1p g12.5e2)') '  -- Temperature term B (K^-1): ', ttm%sh_B * 1.0e4_wp
      Call info(messages, 3, .true.)

    case (2, 6) !linear
      Write (messages(1), '(a)') '  Electronic specific heat capacity: Linear'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Max. electronic S.H.C. (kB/atom): ', ttm%Cemax
      Write (messages(3), '(a,1p g12.5e2)') '  -- Fermi temperature (K):', ttm%Tfermi
      Call info(messages, 3, .true.)

    case (3) !tabulated
      Call info('  Electronic specific heat capacity: Tabulated', .true.)

    end select

    select case (ttm%ketype)
    case (0) ! Infinite
      Call info('  Electronic thermal conductivity: Infinity', .true.)

    case (1) ! Constant
      Write (messages(1), '(a)') '  Electronic thermal conductivity: Constant'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Electronic T.C. (W m^-1 K^-1): ', &
           convert_units(ttm%Ka0, 'k_b/ps/A', 'W/m/K')
      Call info(messages, 2, .true.)

    case (2) ! Drude
      Write (messages(1), '(a)') '  Electronic thermal conductivity: Drude model'
      Write (messages(2), '(a,1p g12.5e2)') '  -- T.C. at system temp. (W m^-1 K^-1): ', &
           convert_units(ttm%Ka0, 'k_b/ps/A', 'W/m/K')
      Call info(messages, 2, .true.)

    case (3) ! Tabulated

      Call info('  Electronic thermal conductivity: Tabulated', .true.)

    end select

    select case (ttm%detype)
    case (0) ! Off
      continue
    case (1) ! Constant
      Write (messages(1), '(a)') '  Electronic thermal diffusivity: Constant'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Electronic T.D. (m^2 s^-1): ', convert_units(ttm%diff0, 'ang^2/ps', 'm^2/s')
      Call info(messages, 2, .true.)

    case (2) ! Recip

      Write (messages(1), '(a)') '  Electronic thermal diffusivity: Reciprocal'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Datum electronic T.D. (m^2 s^-1): ', &
           convert_units(ttm%diff0, 'ang^2/ps', 'm^2/s') / thermo%temp
      Write (messages(3), '(a,1p g12.5e2)') '  -- Fermi temperature (K): ', ttm%Tfermi
      Call info(messages, 3, .true.)

    case (3) ! Tabulated

      Call info('  Electronic thermal diffusivity: Tabulated', .true.)

    end select


    Write (message, '(a,i8)') '  Mininum number of atoms for ionic cells: ', ttm%amin
    Call info(message, .true.)

    if (ttm%redistribute) then
      Call info('  Energy redistribution: ON', .true.)
    End If

    Write (message, '(a,1p g12.5e2)') '  Elec. stopping power (eV/nm): ', convert_units(ttm%dedx, 'e.V/ang', 'e.V/nm')
    Call info(message, .true.)

    if (ttm%fluence < zero_plus) then
      select case (ttm%sdepoType)
      case (1)
        Write (messages(1), '(a)') '   Spatial energy deposition: Gaussian'
        Write (messages(2), '(a,1p g12.5e2)') '  -- Sigma of distribution (nm): ', convert_units(ttm%sig, 'ang', 'nm')
        Write (messages(3), '(a,1p g12.5e2)') '  -- Distribution cutoff (nm): ', convert_units(ttm%sigmax * ttm%sig, 'ang', 'nm')
        Call info(messages, 3, .true.)

      case (2)
        Call info('  Spatial energy deposition: Homogeneous', .true.)
      end select

    else
      select case (ttm%sdepoType)
      case (2)
        Write (messages(1), '(a)') '  Spatial energy deposition: Homogeneous Laser'
        Write (messages(2), '(a,1p g12.5e2)') &
             '  -- Absorbed fluence (mJ cm^-2): ', convert_units(ttm%fluence, 'e.V/ang^2', 'mJ/cm^2')
        Write (messages(3), '(a,1p g12.5e2)') '  -- Penetration depth (nm): ', convert_units(ttm%pdepth, 'ang', 'nm')
        Call info(messages, 3, .true.)

      case (3)

        Write (messages(1), '(a)') '  Spatial energy deposition: Z-exponential decaying Laser '
        Write (messages(2), '(a,1p g12.5e2)') &
             '  -- Absorbed fluence at surface (mJ cm^-2): ', convert_units(ttm%fluence, 'e.V/ang^2', 'mJ/cm^2')
        Write (messages(3), '(a,1p g12.5e2)') '  -- Penetration depth (nm): ', convert_units(ttm%pdepth, 'ang', 'nm')
        Call info(messages, 3, .true.)
      end select

    end if

    select case (ttm%tdepotype)
    case (1)
      Write (messages(1), '(a)') '  Temporal energy deposition: Gaussian'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Sigma of distribution (ps): ', ttm%tdepo
      Write (messages(3), '(a,1p g12.5e2)') '  -- Distribution cutoff (ps): ', 2.0_wp * ttm%tcdepo * ttm%tdepo
      Call info(messages, 3, .true.)

    case (2)
      Write (messages(1), '(a)') '  Temporal energy deposition: Decaying Exponential'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Tau of distribution (ps): ', ttm%tdepo
      Write (messages(3), '(a,1p g12.5e2)') '  -- Distribution cutoff (ps): ', ttm%tcdepo * ttm%tdepo
      Call info(messages, 3, .true.)

    case (3)
      Call info('  Temporal energy deposition: Dirac delta', .true.)

    case (4)
      Write (messages(1), '(a)') '  Temporal energy deposition: Square pulse'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Pulse duration (ps): ', ttm%tdepo
      Call info(messages, 2, .true.)

    end select

    Select Case (ttm%gvar)
    Case (1)
      Write (messages(1), '(a)') '  Variable electron-phonon coupling: Homogeneous'
    Case (2)
      Write (messages(1), '(a)') '  Variable electron-phonon coupling: Heterogeneous'
    End Select
    Write (messages(2), '(a)') '    (overrides value given for ensemble, required tabulated stopping terms in g.dat file)'
    Call info(messages, 2, .true.)

    select case(ttm%bcTypeE)
    case(1)
      Call info('  Electronic temperature boundary conditions: Periodic', .true.)

    Case (2)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Dirichlet'
      Write (messages(2), '(a)') '  -- Boundaries set to system temperature'
      Call info(messages, 2, .true.)

    case (3)
      Call info('  Electronic temperature boundary conditions: Neumann', .true.)

    case (4)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Dirichlet (XY), Neumann (Z)'
      Write (messages(2), '(a)') '  -- XY boundaries set to system temperature'
      Call info(messages, 2, .true.)

    Case (5)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Robin'
      Write (messages(2), '(a,1p g12.5e2)') '  -- Temperature leakage (%): ', convert_units(ttm%fluxout, '', '%')
      Call info(messages, 2, .true.)

    case (6)
      Write (messages(1), '(a)') '  Electronic temperature boundary conditions: Robin (XY), Neumann (Z):'
      Write (messages(2), '(a,1p g12.5e2)') '  -- XY Temperature leakage (%): ', convert_units(ttm%fluxout, '', '%')
      Call info(messages, 3, .true.)

    end select

    Write (message, '(a,1p g12.5e2)') '  Electron-ion coupling offset (ps): ', ttm%ttmoffset
    Call info(message, .true.)

    if (ttm%oneway) Call info('  One-way electron-phonon coupling: ON', .true.)

    if (ttm%ttmstats > 0) then
      Write (messages(1), '(a)') '  TTM statistics file: ON'
      Write (messages(2), '(a,i0.1)') '  -- TTM statistics file interval (steps): ', ttm%ttmstats
      Call info(messages, 2, .true.)
    end if

    if (ttm%ttmtraj > 0) then
      Write (messages(1), '(a)') '  TTM trajectory (temperature profile) file: ON'
      Write (messages(2), '(a,i0.1)') '  -- TTM trajectory file interval (steps): ', ttm%ttmtraj
      Call info(messages, 2, .true.)
    end if

  End Subroutine write_ttm

  Subroutine initialise_control(table)
    !!-----------------------------------------------------------------------
    !!
    !! Initialise all control parameters to their default
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), intent(   Out ) :: table

    call table%init(PARAMS_TABLE_SIZE)

    run_properties: block
      call table%set('title', control_parameter( &
           key = 'title', &
           name = 'Title', &
           val = '', &
           description = "Run title", &
           data_type = DATA_STRING))

      call table%set('simulation_method', control_parameter( &
           key = "simulation_method", &
           name = "Simulation Method ", &
           val = "md", &
           description = "Set trajectory output, options: MD, EVB, FFS", &
           data_type = DATA_OPTION))

      call table%set("random_seed", control_parameter( &
           key = "random_seed", &
           name = "Random seed", &
           val = "1 2 3", &
           description = "Set random seed", &
           data_type = DATA_VECTOR3))

      call table%set("density_variance", control_parameter( &
           key = "density_variance", &
           name = "Expected density variance", &
           val = "0.0", &
           units = "%", &
           internal_units = "", &
           description = "Set expected density variance for determining maximum array sizes", &
           data_type = DATA_FLOAT))

      call table%set("data_dump_frequency", control_parameter( &
           key = "data_dump_frequency", &
           name = "Data dumping", &
           val = "1000", &
           units = "steps", &
           internal_units = "steps", &
           description = "Set data dumping frequency", &
           data_type = DATA_FLOAT))

      call table%set("subcell_threshold", control_parameter( &
           key = "subcell_threshold", &
           name = "Subcelling threshold density", &
           val = "50.0", &
           description = "Set subcelling threshold density", &
           data_type = DATA_FLOAT))

      run_times: block
        call table%set("time_run", control_parameter( &
             key = "time_run", &
             name = "Calculation run length", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Set calculation run length", &
             data_type = DATA_FLOAT))

        call table%set("time_equilibration", control_parameter( &
             key = "time_equilibration", &
             name = "Equilibration run length", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Set equilibration run length", &
             data_type = DATA_FLOAT))

        call table%set("time_job", control_parameter( &
             key = "time_job", &
             name = "Calculation job length", &
             val = "-1.0", &
             units = "hr", &
             internal_units = "s", &
             description = "Set total job time before attempted safe closure", &
             data_type = DATA_FLOAT))

        call table%set("time_close", control_parameter( &
             key = "time_close", &
             name = "Calculation close length", &
             val = "-1.0", &
             units = "min", &
             internal_units = "s", &
             description = "Estimated closure time for finite-time jobs", &
             data_type = DATA_FLOAT))
      end block run_times

      statistics: block
        call table%set("stats_frequency", control_parameter( &
             key = "stats_frequency", &
             name = "Stats Print Frequency", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Set frequency of stats sampling", &
             data_type = DATA_FLOAT))

        call table%set("stack_size", control_parameter( &
             key = "stack_size", &
             name = "Rolling average stack size", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Set rolling average stack to n timesteps", &
             data_type = DATA_FLOAT))

        call table%set("record_equilibration", control_parameter( &
             key = "record_equilibration", &
             name = "Record equilibration", &
             val = "off", &
             description = "Include equilibration in outputs", &
             data_type = DATA_BOOL))

        call table%set("print_per_particle_contrib", control_parameter( &
             key = "print_per_particle_contrib", &
             name = "Per-particle contributions", &
             val = "off", &
             description = "Calculate and print per-particle contributions to energy, force and stress to file every stats step", &
             data_type = DATA_BOOL))

        call table%set("print_probability_distribution", control_parameter( &
             key = "print_probability_distribution", &
             name = "Print probability distribution", &
             val = "off", &
             description = "Calculate and print probability distribution (enforces RDF print)", &
             data_type = DATA_BOOL))

        analysis_options: block
          call table%set("analyse_all", control_parameter( &
               key = "analyse_all", &
               name = "Analyse all", &
               val = "off", &
               description = "Enable analysis for all bonds, angles, dihedrals and inversions", &
               data_type = DATA_BOOL))

          call table%set("analyse_angles", control_parameter( &
               key = "analyse_angles", &
               name = "Analyse angles", &
               val = "off", &
               description = "Enable analysis for all angles", &
               data_type = DATA_BOOL))

          call table%set("analyse_bonds", control_parameter( &
               key = "analyse_bonds", &
               name = "Analyse bonds", &
               val = "off", &
               description = "Enable analysis for all bonds", &
               data_type = DATA_BOOL))

          call table%set("analyse_dihedrals", control_parameter( &
               key = "analyse_dihedrals", &
               name = "Analyse dihedrals", &
               val = "off", &
               description = "Enable analysis for all dihedrals", &
               data_type = DATA_BOOL))

          call table%set("analyse_inversions", control_parameter( &
               key = "analyse_inversions", &
               name = "Analyse inversions", &
               val = "off", &
               description = "Enable analysis inversions", &
               data_type = DATA_BOOL))

          call table%set("analyse_frequency", control_parameter( &
               key = "analyse_frequency", &
               name = "Analysis frequency", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of analysis data", &
               data_type = DATA_FLOAT))

          call table%set("analyse_frequency_bonds", control_parameter( &
               key = "analyse_frequency_bonds", &
               name = "Analysis frequency for bonds", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of bonds data analysis", &
               data_type = DATA_FLOAT))

          call table%set("analyse_frequency_angles", control_parameter( &
               key = "analyse_frequency_angles", &
               name = "Analysis frequency for angles", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of angles data analysis", &
               data_type = DATA_FLOAT))

          call table%set("analyse_frequency_dihedrals", control_parameter( &
               key = "analyse_frequency_dihedrals", &
               name = "Analysis frequency for dihedrals", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of dihedrals data analysis", &
               data_type = DATA_FLOAT))


          call table%set("analyse_frequency_inversions", control_parameter( &
               key = "analyse_frequency_inversions", &
               name = "Analysis frequency for inversions", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of inversions data analysis", &
               data_type = DATA_FLOAT))

          call table%set("analyse_max_dist", control_parameter( &
               key = "analyse_max_dist", &
               name = "Max analyse distance", &
               val = "2.0", &
               units = "ang", &
               internal_units = "internal_l", &
               description = "Set cutoff for bonds analysis", &
               data_type = DATA_FLOAT))

          call table%set("analyse_num_bins", control_parameter( &
               key = "analyse_num_bins", &
               name = "Analysis number of bins", &
               val = "-1", &
               description = "Set number of bins to be used in bonding analysis", &
               data_type = DATA_INT))

          call table%set("analyse_num_bins_bonds", control_parameter( &
               key = "analyse_num_bins_bonds", &
               name = "Analysis number of bins for bonds", &
               val = "0", &
               description = "Set number of bins to be used in bond analysis", &
               data_type = DATA_INT))

          call table%set("analyse_num_bins_angles", control_parameter( &
               key = "analyse_num_bins_angles", &
               name = "Analysis number of bins for angles", &
               val = "0", &
               description = "Set number of bins to be used in angle analysis", &
               data_type = DATA_INT))

          call table%set("analyse_num_bins_dihedrals", control_parameter( &
               key = "analyse_num_bins_dihedrals", &
               name = "Analysis number of bins for dihedrals", &
               val = "0", &
               description = "Set number of bins to be used in dihedral analysis", &
               data_type = DATA_INT))

          call table%set("analyse_num_bins_inversions", control_parameter( &
               key = "analyse_num_bins_inversions", &
               name = "Analysis number of bins for inversions", &
               val = "0", &
               description = "Set number of bins to be used in inversion analysis", &
               data_type = DATA_INT))
        end block analysis_options

        msd: block
          call table%set("msd_calculate", control_parameter( &
               key = "msd_calculate", &
               name = "MSD calculating", &
               val = "off", &
               description = "Enable calculation of MSD", &
               data_type = DATA_BOOL))

          call table%set("msd_print", control_parameter( &
               key = "msd_print", &
               name = "MSD printing", &
               val = "off", &
               description = "Enable printing of MSD", &
               data_type = DATA_BOOL))

          call table%set("msd_start", control_parameter( &
               key = "msd_start", &
               name = "MSD Start calculating", &
               val = "0", &
               units = "steps", &
               internal_units = "steps", &
               description = "Start timestep for dumping MSD configurations", &
               data_type = DATA_FLOAT))

          call table%set("msd_frequency", control_parameter( &
               key = "msd_frequency", &
               name = "MSD calculation interval", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Interval between dumping MSD configurations", &
               data_type = DATA_FLOAT))
        end block msd

        traj: block
          call table%set("traj_calculate", control_parameter( &
               key = "traj_calculate", &
               name = "Trajectory calculating", &
               val = "off", &
               description = "Enable calculation of trajectory", &
               data_type = DATA_BOOL))

          call table%set("traj_key", control_parameter( &
               key = "traj_key", &
               name = "Trajectory key ", &
               val = "pos", &
               description = "Set trajectory output, options: pos, pos-vel, pos-vel-force, compressed", &
               data_type = DATA_OPTION))

          call table%set("traj_start", control_parameter( &
               key = "traj_start", &
               name = "trajectory Start calculating", &
               val = "0", &
               units = "steps", &
               internal_units = "steps", &
               description = "Start timestep for dumping trajectory configurations", &
               data_type = DATA_FLOAT))

          call table%set("traj_interval", control_parameter( &
               key = "traj_interval", &
               name = "trajectory calculation interval", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Interval between dumping trajectory configurations", &
               data_type = DATA_FLOAT))
        end block traj

        defects: block
          call table%set("defects_calculate", control_parameter( &
               key = "defects_calculate", &
               name = "Defects calculating", &
               val = "off", &
               description = "Enable calculation of defects", &
               data_type = DATA_BOOL))

          call table%set("defects_start", control_parameter( &
               key = "defects_start", &
               name = "Defects Start calculating", &
               val = "0", &
               units = "steps", &
               internal_units = "steps", &
               description = "Start timestep for dumping defects configurations", &
               data_type = DATA_FLOAT))

          call table%set("defects_interval", control_parameter( &
               key = "defects_interval", &
               name = "Defects calculation interval", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Interval between dumping defects configurations", &
               data_type = DATA_FLOAT))

          call table%set("defects_distance", control_parameter( &
               key = "defects_distance", &
               name = "Defects distance condition ", &
               val = "0.75", &
               units = "ang", &
               internal_units = "internal_l", &
               description = "Set defects condition", &
               data_type = DATA_FLOAT))

          call table%set("defects_backup", control_parameter( &
               key = "defects_backup", &
               name = "Defects backup ", &
               val = "off", &
               description = "Enable defects backup", &
               data_type = DATA_BOOL))

        end block defects

        displacements: block
          call table%set("displacements_calculate", control_parameter( &
               key = "displacements_calculate", &
               name = "Displacements calculating", &
               val = "off", &
               description = "Enable calculation of displacements", &
               data_type = DATA_BOOL))

          call table%set("displacements_start", control_parameter( &
               key = "displacements_start", &
               name = "Displacements Start calculating", &
               val = "0", &
               units = "steps", &
               internal_units = "steps", &
               description = "Start timestep for dumping displacements configurations", &
               data_type = DATA_FLOAT))

          call table%set("displacements_interval", control_parameter( &
               key = "displacements_interval", &
               name = "Displacements calculation interval", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Interval between dumping displacements configurations", &
               data_type = DATA_FLOAT))

          call table%set("displacements_distance", control_parameter( &
               key = "displacements_distance", &
               name = "Displacements distance condition ", &
               val = "0.75", &
               units = "ang", &
               internal_units = "internal_l", &
               description = "Set displacements condition", &
               data_type = DATA_FLOAT))

        end block displacements

        coord: block
          call table%set("coord_calculate", control_parameter( &
               key = "coord_calculate", &
               name = "Coordination calculating", &
               val = "off", &
               description = "Enable calculation of Coordination", &
               data_type = DATA_BOOL))

          call table%set("coord_ops", control_parameter( &
               key = "coord_ops", &
               name = "Coord ops", &
               val = "0", &
               description = "Set Coordops", &
               data_type = DATA_INT))

          call table%set("coord_start", control_parameter( &
               key = "coord_start", &
               name = "Coordination Start calculating", &
               val = "0", &
               units = "steps", &
               internal_units = "steps", &
               description = "Start timestep for dumping Coordination configurations", &
               data_type = DATA_FLOAT))

          call table%set("coord_interval", control_parameter( &
               key = "coord_interval", &
               name = "Coordination calculation interval", &
               val = "100", &
               units = "steps", &
               internal_units = "steps", &
               description = "Interval between dumping Coordination configurations", &
               data_type = DATA_FLOAT))
        end block coord

        adf: block
          call table%set("adf_calculate", control_parameter( &
               key = "adf_calculate", &
               name = "ADF calculating", &
               val = "off", &
               description = "Enable calculation of Adf", &
               data_type = DATA_BOOL))

          call table%set("adf_frequency", control_parameter( &
               key = "adf_frequency", &
               name = "ADF Sampling Frequency", &
               val = "100", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of Adf sampling", &
               data_type = DATA_FLOAT))

          call table%set("adf_precision", control_parameter( &
               key = "adf_precision", &
               name = "ADF Precision", &
               val = "0.0", &
               description = "Set precision in Adf analysis", &
               data_type = DATA_FLOAT))
        end block adf

        rdf: block
          call table%set("rdf_calculate", control_parameter( &
               key = "rdf_calculate", &
               name = "RDF calculating", &
               val = "off", &
               description = "Enable calculation of RDF", &
               data_type = DATA_BOOL))

          call table%set("rdf_print", control_parameter( &
               key = "rdf_print", &
               name = "RDF printing", &
               val = "on", &
               description = "Enable printing of RDF", &
               data_type = DATA_BOOL))

          call table%set("rdf_frequency", control_parameter( &
               key = "rdf_frequency", &
               name = "RDF Sampling Frequency", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of RDF sampling", &
               data_type = DATA_FLOAT))

          call table%set("rdf_binsize", control_parameter( &
               key = "rdf_binsize", &
               name = "RDF number of bins", &
               val = "0.05", &
               units = "ang", &
               internal_units = "ang", &
               description = "Set number of bins to be used in RDF analysis", &
               data_type = DATA_FLOAT))

          call table%set("rdf_error_analysis", control_parameter( &
               key = "rdf_error_analysis", &
               name = "RDF Error Analysis", &
               val = "off", &
               description = "Enable RDF error analysis, options: Off, Jackknife, Block", &
               data_type = DATA_OPTION))

          call table%set("rdf_error_analysis_blocks", control_parameter( &
               key = "rdf_error_analysis_blocks", &
               name = "Num RDF Error Analysis blocks", &
               val = "1", &
               description = "Set number of RDF error analysis blocks", &
               data_type = DATA_INT))
        end block rdf

        zden: block
          call table%set("zden_calculate", control_parameter( &
               key = "zden_calculate", &
               name = "Zden calculating", &
               val = "off", &
               description = "Enable calculation of Zden", &
               data_type = DATA_BOOL))

          call table%set("zden_print", control_parameter( &
               key = "zden_print", &
               name = "Zden printing", &
               val = "on", &
               description = "Enable printing of Zden", &
               data_type = DATA_BOOL))

          call table%set("zden_frequency", control_parameter( &
               key = "zden_frequency", &
               name = "ZDen Sampling Frequency", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of ZDen sampling", &
               data_type = DATA_FLOAT))

          call table%set("zden_binsize", control_parameter( &
               key = "zden_binsize", &
               name = "ZDen number of bins", &
               val = "0.05", &
               units = "ang", &
               internal_units = "ang", &
               description = "Set number of bins to be used in ZDen analysis", &
               data_type = DATA_FLOAT))
        end block zden

        vaf: block
          call table%set("vaf_calculate", control_parameter( &
               key = "vaf_calculate", &
               name = "VAF calculating", &
               val = "off", &
               description = "Enable calculation of VAF", &
               data_type = DATA_BOOL))

          call table%set("vaf_print", control_parameter( &
               key = "vaf_print", &
               name = "VAF printing", &
               val = "on", &
               description = "Enable printing of VAF", &
               data_type = DATA_BOOL))

          call table%set("vaf_frequency", control_parameter( &
               key = "vaf_frequency", &
               name = "VAF Sampling Frequency", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of VAF sampling", &
               data_type = DATA_FLOAT))

          call table%set("vaf_binsize", control_parameter( &
               key = "vaf_binsize", &
               name = "VAF number of bins", &
               val = "0", &
               description = "Set number of bins to be used in VAF analysis", &
               data_type = DATA_INT))

          call table%set("vaf_averaging", control_parameter( &
               key = "vaf_averaging", &
               name = "VAF Time Averaging", &
               val = "on", &
               description = "Ignore time-averaging of VAF, "// &
               "report all calculated VAF to VAFDAT files and final profile to OUTPUT", &
               data_type = DATA_BOOL))
        end block vaf

        call table%set("currents_calculate", control_parameter( &
             key = "currents_calculate", &
             name = "Calculate currents", &
             val = "off", &
             description = "Enable calculation of currents", &
             data_type = DATA_BOOL))

      end block statistics

    end block run_properties

    io: block
      call table%set("print_frequency", control_parameter( &
           key = "print_frequency", &
           name = "Results Print Frequency", &
           val = "0", &
           units = "steps", &
           internal_units = "steps", &
           description = "Set frequency of results sampling", &
           data_type = DATA_FLOAT))

      units: block
        call table%set("io_units_scheme", control_parameter( &
             key = "io_units_scheme", &
             name = "I/O units scheme", &
             val = "internal", &
             description = "Set I/O units scheme, options: internal, si, atomic", &
             data_type = DATA_OPTION))

        call table%set("io_units_length", control_parameter( &
             key = "io_units_length", &
             name = "I/O units length", &
             val = "internal_l", &
             description = "Set I/O units for length", &
             data_type = DATA_OPTION))

        call table%set("io_units_time", control_parameter( &
             key = "io_units_time", &
             name = "I/O units time", &
             val = "internal_t", &
             description = "Set I/O units for time", &
             data_type = DATA_OPTION))

        call table%set("io_units_mass", control_parameter( &
             key = "io_units_mass", &
             name = "I/O units mass", &
             val = "internal_m", &
             description = "Set I/O units for mass", &
             data_type = DATA_OPTION))

        call table%set("io_units_charge", control_parameter( &
             key = "io_units_charge", &
             name = "I/O units charge", &
             val = "internal_q", &
             description = "Set I/O units for charge", &
             data_type = DATA_OPTION))

        call table%set("io_units_energy", control_parameter( &
             key = "io_units_energy", &
             name = "I/O units energy", &
             val = "internal_e", &
             description = "Set I/O units for energy", &
             data_type = DATA_OPTION))

        call table%set("io_units_pressure", control_parameter( &
             key = "io_units_pressure", &
             name = "I/O units pressure", &
             val = "internal_p", &
             description = "Set I/O units for pressure", &
             data_type = DATA_OPTION))

        call table%set("io_units_force", control_parameter( &
             key = "io_units_force", &
             name = "I/O units force", &
             val = "internal_f", &
             description = "Set I/O units for force", &
             data_type = DATA_OPTION))

        call table%set("io_units_velocity", control_parameter( &
             key = "io_units_velocity", &
             name = "I/O units velocity", &
             val = "internal_v", &
             description = "Set I/O units for velocity", &
             data_type = DATA_OPTION))

        call table%set("io_units_power", control_parameter( &
             key = "io_units_power", &
             name = "I/O units power", &
             val = "internal_e/internal_t", &
             description = "Set I/O units for power", &
             data_type = DATA_OPTION))

        call table%set("io_units_surface_tension", control_parameter( &
             key = "io_units_surface_tension", &
             name = "I/O units surface tension", &
             val = "internal_f/internal_l", &
             description = "Set I/O units for surface tension", &
             data_type = DATA_OPTION))

        call table%set("io_units_emf", control_parameter( &
             key = "io_units_emf", &
             name = "I/O units emf", &
             val = "internal_e/internal_q", &
             description = "Set I/O units for electromotive force", &
             data_type = DATA_OPTION))

      end block units

      io_read: block
        call table%set("io_read_method", control_parameter( &
             key = "io_read_method", &
             name = "I/O read method", &
             val = "mpiio", &
             description = "Set I/O read method, possible read methods: "//&
             "mpiio, direct, netcdf, master", &
             data_type = DATA_OPTION))

        call table%set("io_read_readers", control_parameter( &
             key = "io_read_readers", &
             name = "Num parallel I/O readers", &
             val = "0", &
             units = "(Automatic)", &
             description = "Set number of parallel I/O readers", &
             data_type = DATA_INT))

        call table%set("io_read_batch_size", control_parameter( &
             key = "io_read_batch_size", &
             name = "I/O reader batch size", &
             val = "0", &
             units = "(Automatic)", &
             description = "Set I/O reader batch size", &
             data_type = DATA_INT))

        call table%set("io_read_buffer_size", control_parameter( &
             key = "io_read_buffer_size", &
             name = "I/O reader buffer size", &
             val = "0", &
             units = "(Automatic)", &
             description = "Set I/O reader buffer size", &
             data_type = DATA_INT))

        call table%set("io_read_error_check", control_parameter( &
             key = "io_read_error_check", &
             name = "I/O error check on read", &
             val = "off", &
             description = "Enable extended error checking on read", &
             data_type = DATA_BOOL))

        call table%set('io_read_ascii_revold', control_parameter( &
             key = 'io_read_ascii_revold', &
             name = 'Read plain-text REVOLD', &
             val = 'off', &
             data_type = DATA_BOOL, &
             description = "Read human-readable (ASCII) REVOLD file"))

      end block io_read

      io_write: block
        call table%set("io_write_method", control_parameter( &
             key = "io_write_method", &
             name = "I/O write method", &
             val = "mpiio", &
             description = "Set I/O write method, possible write methods: "//&
             "mpiio, direct, netcdf, master", &
             data_type = DATA_OPTION))

        call table%set("io_write_writers", control_parameter( &
             key = "io_write_writers", &
             name = "Num parallel I/O writers", &
             val = "0", &
             units = "(Automatic)", &
             description = "Set number of parallel I/O writers", &
             data_type = DATA_INT))

        call table%set("io_write_batch_size", control_parameter( &
             key = "io_write_batch_size", &
             name = "I/O writer batch size", &
             val = "0", &
             units = "(Automatic)", &
             description = "Set I/O writer batch size", &
             data_type = DATA_INT))

        call table%set("io_write_buffer_size", control_parameter( &
             key = "io_write_buffer_size", &
             name = "I/O writer buffer size", &
             val = "0", &
             units = "(Automatic)", &
             description = "Set I/O writer buffer size", &
             data_type = DATA_INT))

        call table%set("io_write_sorted", control_parameter( &
             key = "io_write_sorted", &
             name = "I/O sorted writing", &
             val = "on", &
             description = "Enable sorted output for atomic data", &
             data_type = DATA_BOOL))

        call table%set("io_write_error_check", control_parameter( &
             key = "io_write_error_check", &
             name = "I/O error check on write", &
             val = "off", &
             description = "Enable extended error checking on write", &
             data_type = DATA_BOOL))

        call table%set("io_write_netcdf_format", control_parameter( &
             key = "io_write_netcdf_format", &
             name = "Netcdf Format", &
             val = "64bit", &
             description = "Set netcdf write format, options: amber, 32bit, 32-bit, 64-bit, 64bit", &
             data_type = DATA_OPTION))

        call table%set('io_write_ascii_revive', control_parameter( &
             key = 'io_write_ascii_revive', &
             name = 'Write plain-text REVIVE', &
             val = 'off', &
             data_type = DATA_BOOL, &
             description = "Write REVIVE as a human-readable (ASCII) file"))

      end block io_write

      io_file: block
        call table%set("io_file_output", control_parameter( &
             key = "io_file_output", &
             name = "Output filepath", &
             val = "OUTPUT", &
             description = "Set output filepath, special options: SCREEN, NONE", &
             data_type = DATA_STRING))

        call table%set("io_file_config", control_parameter( &
             key = "io_file_config", &
             name = "Config filepath", &
             val = "CONFIG", &
             description = "Set input configuration filepath", &
             data_type = DATA_STRING))

        call table%set("io_file_field", control_parameter( &
             key = "io_file_field", &
             name = "Field filepath", &
             val = "FIELD", &
             description = "Set input field filepath", &
             data_type = DATA_STRING))

        call table%set("io_file_statis", control_parameter( &
             key = "io_file_statis", &
             name = "Statistics filepath", &
             val = "STATIS", &
             description = "Set output statistics filepath, special options: NONE", &
             data_type = DATA_STRING))

        call table%set("io_file_history", control_parameter( &
             key = "io_file_history", &
             name = "History filepath", &
             val = "HISTORY", &
             description = "Set output history filepath, special options: NONE", &
             data_type = DATA_STRING))

        call table%set("io_file_historf", control_parameter( &
             key = "io_file_historf", &
             name = "Historf filepath", &
             val = "HISTORF", &
             description = "Set output historf filepath, special options: NONE", &
             data_type = DATA_STRING))

        call table%set("io_file_revive", control_parameter( &
             key = "io_file_revive", &
             name = "Revive filepath", &
             val = "REVIVE", &
             description = "Set output revive filepath, special options: NONE", &
             data_type = DATA_STRING))

        call table%set("io_file_revold", control_parameter( &
             key = "io_file_revold", &
             name = "Revold filepath", &
             val = "REVOLD", &
             description = "Set output revold filepath, special options: NONE", &
             data_type = DATA_STRING))

        call table%set("io_file_revcon", control_parameter( &
             key = "io_file_revcon", &
             name = "Revcon filepath", &
             val = "REVCON", &
             description = "Set output revcon filepath, special options: NONE", &
             data_type = DATA_STRING))

        call table%set('io_file_rdf', control_parameter( &
             key = 'io_file_rdf', &
             name = 'rdf filepath', &
             val = 'RDFDAT', &
             description = 'Set output RDF filepath, special options: NONE', &
             data_type = DATA_STRING))

        call table%set('io_file_msd', control_parameter( &
             key = 'io_file_msd', &
             name = 'msd filepath', &
             val = 'MSDTMP', &
             description = 'Set output MSD filepath, special options: NONE', &
             data_type = DATA_STRING))

        call table%set('io_file_tabbnd', control_parameter( &
             key = 'io_file_tabbnd', &
             name = 'tabbnd filepath', &
             val = 'TABBND', &
             description = 'Set input TABBND filepath', &
             data_type = DATA_STRING))

        call table%set('io_file_tabang', control_parameter( &
             key = 'io_file_tabang', &
             name = 'tabang filepath', &
             val = 'TABANG', &
             description = 'Set input TABANG filepath', &
             data_type = DATA_STRING))

        call table%set('io_file_tabdih', control_parameter( &
             key = 'io_file_tabdih', &
             name = 'tabdih filepath', &
             val = 'TABDIH', &
             description = 'Set input TABDIH filepath', &
             data_type = DATA_STRING))

        call table%set('io_file_tabinv', control_parameter( &
             key = 'io_file_tabinv', &
             name = 'tabinv filepath', &
             val = 'TABINV', &
             description = 'Set input TABINV filepath', &
             data_type = DATA_STRING))

        call table%set('io_file_tabvdw', control_parameter( &
             key = 'io_file_tabvdw', &
             name = 'tabvdw filepath', &
             val = 'TABVDW', &
             description = 'Set input TABVDW filepath', &
             data_type = DATA_STRING))

        call table%set('io_file_tabeam', control_parameter( &
             key = 'io_file_tabeam', &
             name = 'tabeam filepath', &
             val = 'TABEAM', &
             description = 'Set input TABEAM filepath', &
             data_type = DATA_STRING))

      end block io_file

      call table%set('output_energy', control_parameter( &
           key = 'output_energy', &
           name= 'Output Final Energy', &
           val = 'off', &
           data_type = DATA_BOOL, &
           description = "Output final energy e_tot in output file"))

      call table%set("ignore_config_indices", control_parameter( &
           key = "ignore_config_indices", &
           name = "Ignore Config Indices", &
           val = "off", &
           description = "Ignore indices as defined in CONFIG and use read order instead", &
           data_type = DATA_BOOL))

      call table%set("print_topology_info", control_parameter( &
           key = "print_topology_info", &
           name = "Print topology information ", &
           val = "off", &
           description = "Print topology information in output file", &
           data_type = DATA_BOOL))

      call table%set("print_level", control_parameter( &
           key = "print_level", &
           name = "Print level", &
           val = "1", &
           description = "Disable unnecessary printing, levels: 0 - silent, 1 - quiet, 2 - standard, 3 - full", &
           data_type = DATA_INT))

      call table%set("timer_depth", control_parameter( &
           key = "timer_depth", &
           name = "Timer print level", &
           val = "4", &
           description = "Do not display timers beyond so this depth", &
           data_type = DATA_INT))

      call table%set("timer_per_mpi", control_parameter( &
           key = "timer_per_mpi", &
           name = "Per Process timing", &
           val = "off", &
           description = "Time each MPI process individually", &
           data_type = DATA_BOOL))

    end block io

    simulation_properties: block

      timestep: block
        call table%set("timestep", control_parameter( &
             key = "timestep", &
             name = "Timestep", &
             val = "0.0", &
             units = "internal_t", &
             internal_units = "internal_t", &
             description = "Set calculation timestep or initial timestep for variable timestep calculations", &
             data_type = DATA_FLOAT))

        call table%set("timestep_variable", control_parameter( &
             key = "timestep_variable", &
             name = "Variable timestep", &
             val = "off", &
             description = "Enable variable timestep", &
             data_type = DATA_BOOL))

        call table%set("timestep_variable_min_dist", control_parameter( &
             key = "timestep_variable_min_dist", &
             name = "Variable timestep minimum distance", &
             val = "0.03", &
             units = "ang", &
             internal_units = "internal_l", &
             description = "Set minimum permissible distance for variable timestep", &
             data_type = DATA_FLOAT))

        call table%set("timestep_variable_max_dist", control_parameter( &
             key = "timestep_variable_max_dist", &
             name = "Variable timestep maximum distance", &
             val = "0.1", &
             units = "ang", &
             internal_units = "internal_l", &
             description = "Set maximum permissible distance for variable timestep", &
             data_type = DATA_FLOAT))

        call table%set("timestep_variable_max_delta", control_parameter( &
             key = "timestep_variable_max_delta", &
             name = "Variable timestep max delta", &
             val = "0.0", &
             units = "internal_t", &
             internal_units = "internal_t", &
             description = "Set maximum timestep delta for variable timestep", &
             data_type = DATA_FLOAT))
      end block timestep

      ensemble: block
        call table%set("ensemble", control_parameter( &
             key = "ensemble", &
             name = "Ensemble constraints", &
             val = "NVE", &
             description = "Set ensemble constraints, options: NVE, PMF, NVT, NPT, NST",&
             data_type = DATA_OPTION))

        call table%set("ensemble_method", control_parameter( &
             key = "ensemble_method", &
             name = "Ensemble method", &
             val = "", &
             description = "Set ensemble method, options "// &
             "NVT: Evans, Langevin, Andersen, Berendsen, Hoover, gentle, ttm, dpds1, dpds2. "// &
             "NP|ST: Langevin, Berendsen, Hoover, MTK.", &
             data_type = DATA_OPTION))

        call table%set("ensemble_thermostat_coupling", control_parameter( &
             key = "ensemble_thermostat_coupling", &
             name = "Thermostat coupling", &
             val = "0.0", &
             units = "ps", &
             internal_units = "ps", &
             description = "Set thermostat relaxation/decorrelation times (use ensemble_thermostat_friction for langevin)", &
             data_type = DATA_FLOAT))

        call table%set("ensemble_dpd_order", control_parameter( &
             key = "ensemble_dpd_order", &
             name = "Ensemble DPD order", &
             val = "off", &
             description = "Set dpd method, options: off, first, second", &
             data_type = DATA_OPTION))

        call table%set("ensemble_dpd_drag", control_parameter( &
             key = "ensemble_dpd_drag", &
             name = "Ensemble DPD drag coefficient", &
             val = "0.0", &
             units = "Da/ps", &
             internal_units = "Da/ps", &
             description = "Set DPD drag coefficient", &
             data_type = DATA_FLOAT))

        call table%set("ensemble_thermostat_friction", control_parameter( &
             key = "ensemble_thermostat_friction", &
             name = "Thermostat Friction", &
             val = "0.0", &
             units = "ps^-1", &
             internal_units = "ps^-1", &
             description = "Set thermostat friction for langevin and gentle stochastic thermostats", &
             data_type = DATA_FLOAT))

        call table%set("ensemble_thermostat_softness", control_parameter( &
             key = "ensemble_thermostat_softness", &
             name = "Ensemble thermostat softness", &
             val = "0.0", &
             units = "", &
             internal_units = "", &
             description = "Set thermostat softness for Andersen thermostat", &
             data_type = DATA_FLOAT))

        call table%set("ensemble_barostat_coupling", control_parameter( &
             key = "ensemble_barostat_coupling", &
             name = "Barostat coupling", &
             val = "0.0", &
             units = "ps", &
             internal_units = "ps", &
             description = "Set barostat relaxation/decorrelation times (use ensemble_barostat_friction for langevin)", &
             data_type = DATA_FLOAT))

        call table%set("ensemble_barostat_friction", control_parameter( &
             key = "ensemble_barostat_friction", &
             name = "Barostat Friction", &
             val = "0.0", &
             units = "ps^-1", &
             internal_units = "ps^-1", &
             description = "Set barostat friction", &
             data_type = DATA_FLOAT))

        call table%set("ensemble_semi_isotropic", control_parameter( &
             key = "ensemble_semi_isotropic", &
             name = "Ensemble semi-isotropic constraint", &
             val = "off", &
             description = "Enable semi-isotropic barostat constraints, options: area, tension, orthorhombic", &
             data_type = DATA_OPTION))

        call table%set("ensemble_semi_orthorhombic", control_parameter( &
             key = "ensemble_semi_orthorhombic", &
             name = "Ensemble semi-orthorhombic constraint", &
             val = "off", &
             description = "Enable semi-orthorhombic barostat constraints", &
             data_type = DATA_BOOL))

        call table%set("ensemble_tension", control_parameter( &
             key = "ensemble_tension", &
             name = "Constrained surface tension", &
             val = "0.0", &
             units = "N/m", &
             internal_units = "dyn/cm", &
             description = "Set tension in NPngT calctulation", &
             data_type = DATA_FLOAT))
      end block ensemble

      system_properties: block
        call table%set("pressure_tensor", control_parameter( &
             key = "pressure_tensor", &
             name = "Pressure tensor", &
             val = "0.0 0.0 0.0 0.0 0.0 0.0", &
             units = "katm", &
             internal_units = "internal_p", &
             description = "Set the target pressure tensor for NsT calculations", &
             data_type = DATA_VECTOR6))

        call table%set("pressure_hydrostatic", control_parameter( &
             key = "pressure_hydrostatic", &
             name = "Hydrostatic pressure", &
             val = "0.0", &
             units = "katm", &
             internal_units = "internal_p", &
             description = "Set the target hydrostatic pressure (1/3Tr[P]) for NPT calculations", &
             data_type = DATA_FLOAT))

        call table%set("pressure_perpendicular", control_parameter( &
             key = "pressure_perpendicular", &
             name = "Perpendicular pressure", &
             val = "0.0 0.0 0.0", &
             units = "katm", &
             internal_units = "internal_p", &
             description = "Set the target pressure as x, y, z perpendicular to cell faces for NPT calculations", &
             data_type = DATA_VECTOR3))

        call table%set("temperature", control_parameter( &
             key = "temperature", &
             name = "Initial/Target temperature", &
             val = "0.0", &
             units = "K", &
             internal_units = "K", &
             description = "Set the initial temperature or target temperature (for thermostats)", &
             data_type = DATA_FLOAT))
      end block system_properties

      pseudo_thermostat: block
        call table%set("pseudo_thermostat_method", control_parameter( &
             key = "pseudo_thermostat_method", &
             name = "Pseudo thermostat method", &
             val = "off", &
             description = "Set pseudo thermostat method, possible options: Off, Langevin-Direct, Langevin, Gauss, Direct", &
             data_type = DATA_OPTION))

        call table%set("pseudo_thermostat_width", control_parameter( &
             key = "pseudo_thermostat_width", &
             name = "Pseudo thermostat width", &
             val = "2.0", &
             units = "ang", &
             internal_units = "internal_l", &
             description = "Set the width of thermostatted boundaries for pseudo thermostats", &
             data_type = DATA_FLOAT))

        call table%set("pseudo_thermostat_temperature", control_parameter( &
             key = "pseudo_thermostat_temperature", &
             name = "Pseudo thermostat temperature", &
             val = "0.0", &
             units = "K", &
             internal_units = "K", &
             description = "Set the temperature of the pseudo thermostat", &
             data_type = DATA_FLOAT))
      end block pseudo_thermostat

      impact: block
        call table%set("impact_part_index", control_parameter( &
             key = "impact_part_index", &
             name = "Impact particle index", &
             val = "0", &
             description = "Set particle index for impact simulations", &
             data_type = DATA_INT))

        call table%set("impact_time", control_parameter( &
             key = "impact_time", &
             name = "Impact time", &
             val = "0.0", &
             units = "internal_t", &
             internal_units = "steps", &
             description = "Set time for impact in impact simulations", &
             data_type = DATA_FLOAT))

        call table%set("impact_energy", control_parameter( &
             key = "impact_energy", &
             name = "Impact energy", &
             val = "0.0", &
             units = "ke.V", &
             internal_units = "ke.V", &
             description = "Set impact energy for impact simulations", &
             data_type = DATA_FLOAT))

        call table%set("impact_direction", control_parameter( &
             key = "impact_direction", &
             name = "Impact direction", &
             val = "1.0 1.0 1.0", &
             units = "", &
             internal_units = "", &
             description = "Direction vector for impact simulations", &
             data_type = DATA_VECTOR3))
      end block impact

      ttm: block

        call table%set("ttm_calculate", control_parameter( &
             key = "ttm_calculate", &
             name = "TTM Active", &
             val = "off", &
             description = "Enable calculation of two-temperature model", &
             data_type = DATA_BOOL))

        call table%set("ttm_num_ion_cells", control_parameter( &
             key = "ttm_num_ion_cells", &
             name = "Number of TTM ionic cells", &
             val = "10", &
             description = "Set number of coarse-grained ion temperature cells (CIT)", &
             data_type = DATA_INT))

        call table%set("ttm_num_elec_cells", control_parameter( &
             key = "ttm_num_elec_cells", &
             name = "Number of TTM electronic cells", &
             val = "50 50 50", &
             description = "Set number of coarse-grained electronic temperature cells (CET)", &
             data_type = DATA_VECTOR3))

        call table%set("ttm_metal", control_parameter( &
             key = "ttm_metal", &
             name = "TTM Metallic", &
             val = "off", &
             description = "Specifies parameters for metallic system are required for two-temperature model"// &
             ", i.e. thermal conductivity", &
             data_type = DATA_BOOL))

        call table%set("ttm_heat_cap_model", control_parameter( &
             key = "ttm_heat_cap_model", &
             name = "TTM Heat Capacity model", &
             val = "", &
             description = "Sets model for specific heat capacity in TTM, options: const, linear, tabulated, tanh", &
             data_type = DATA_OPTION))

        call table%set("ttm_heat_cap", control_parameter( &
             key = "ttm_heat_cap", &
             name = "TTM Heat Capacity", &
             val = "0.0", &
             units = "internal_e/internal_m/K", &
             internal_units = "internal_e/internal_m/K", &
             description = "Sets constant, scale or maximum heat capcity in TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_temp_term", control_parameter( &
             key = "ttm_temp_term", &
             name = "TTM Tanh temperature term", &
             val = "0.0", &
             units = "K^-1", &
             internal_units = "K^-1", &
             description = "Set Fermi temperature in TTM, for tanh", &
             data_type = DATA_FLOAT))

        call table%set("ttm_fermi_temp", control_parameter( &
             key = "ttm_fermi_temp", &
             name = "TTM Fermi temperature", &
             val = "0.0", &
             units = "K", &
             internal_units = "K", &
             description = "Set Fermi temperature in TTM, for linear", &
             data_type = DATA_FLOAT))

        call table%set("ttm_elec_cond_model", control_parameter( &
             key = "ttm_elec_cond_model", &
             name = "TTM Electonic conductivity model", &
             val = "", &
             description = "Set electronic conductivity model in TTM, options: Infinite, constant, drude, tabulated", &
             data_type = DATA_OPTION))

        call table%set("ttm_elec_cond", control_parameter( &
             key = "ttm_elec_cond", &
             name = "TTM Electronic conductivity", &
             val = "0.0", &
             units = "W/m/K", &
             internal_units = "k_b/ps/Ang", &
             description = "Set electronic conductivity in TTM ", &
             data_type = DATA_FLOAT))

        call table%set("ttm_diff_model", control_parameter( &
             key = "ttm_diff_model", &
             name = "TTM Diffusion model", &
             val = "", &
             description = "Set diffusion model in TTM, options: constant, recip, tabulated", &
             data_type = DATA_OPTION))

        call table%set("ttm_diff", control_parameter( &
             key = "ttm_diff", &
             name = "TTM Thermal diffusivity", &
             val = "0.0", &
             units = "m^2/s", &
             internal_units = "ang^2/ps", &
             description = "Set TTM thermal diffusivity", &
             data_type = DATA_FLOAT))

        call table%set("ttm_dens_model", control_parameter( &
             key = "ttm_dens_model", &
             name = "TTM Density model", &
             val = "", &
             description = "Set density model in TTM, options are: constant, dynamic", &
             data_type = DATA_OPTION))

        call table%set("ttm_dens", control_parameter( &
             key = "ttm_dens", &
             name = "TTM Density", &
             val = "0.0", &
             units = "ang^-3", &
             internal_units = "internal_l^-3", &
             description = "Set constant density in TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_min_atoms", control_parameter( &
             key = "ttm_min_atoms", &
             name = "TTM minimum atoms", &
             val = "0", &
             description = "Minimum number of atoms needed per ionic temperature cell", &
             data_type = DATA_INT))

        call table%set("ttm_stopping_power", control_parameter( &
             key = "ttm_stopping_power", &
             name = "TTM Electron stopping power", &
             val = "0.0", &
             units = "e.V/nm", &
             internal_units = "e.V/ang", &
             description = "Electronic stopping power of projectile entering electronic system ", &
             data_type = DATA_FLOAT))

        call table%set("ttm_spatial_dist", control_parameter( &
             key = "ttm_spatial_dist", &
             name = "TTM Spatial deposition distribution", &
             val = "", &
             description = "Set the spatial distribution of TTM, options: flat, gaussian, flat-laser, exp-laser", &
             data_type = DATA_OPTION))

        call table%set("ttm_spatial_sigma", control_parameter( &
             key = "ttm_spatial_sigma", &
             name = "TTM Spatial sigma", &
             val = "1.0", &
             units = "nm", &
             internal_units = "internal_l", &
             description = "Set the sigma for spatial distributions of TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_spatial_cutoff", control_parameter( &
             key = "ttm_spatial_cutoff", &
             name = "TTM Spatial cutoff", &
             val = "5.0", &
             units = "nm", &
             internal_units = "nm", &
             description = "Set the cutoff for spatial distributions of TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_fluence", control_parameter( &
             key = "ttm_fluence", &
             name = "TTM laser absorbed energy", &
             val = "0.0", &
             units = "mJ/cm^2", &
             internal_units = "e.V/ang^2", &
             description = "Initial energy deposition into electronic system by laser for TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_penetration_depth", control_parameter( &
             key = "ttm_penetration_depth", &
             name = "TTM laser penetration depth", &
             val = "0.0", &
             units = "nm", &
             internal_units = "ang", &
             description = "Set laser penetration depth for TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_laser_type", control_parameter( &
             key = "ttm_laser_type", &
             name = "TTM Deposition type", &
             val = "flat", &
             description = "Set laser deposition type. options: flat, exponential", &
             data_type = DATA_OPTION))

        call table%set("ttm_temporal_dist", control_parameter( &
             key = "ttm_temporal_dist", &
             name = "TTM Temporal distribution", &
             val = "", &
             description = "Set temporal distribution for TTM, options: gaussian, exponential, delta, square", &
             data_type = DATA_OPTION))

        call table%set("ttm_temporal_duration", control_parameter( &
             key = "ttm_temporal_duration", &
             name = "TTM Duration", &
             val = "0.001", &
             units = "ps", &
             internal_units = "ps", &
             description = "Set duration of energy deposition for TTM (gaussian, exponential, square)", &
             data_type = DATA_FLOAT))

        call table%set("ttm_temporal_cutoff", control_parameter( &
             key = "ttm_temporal_cutoff", &
             name = "TTM Time Cutoff", &
             val = "5.0", &
             units = "ps", &
             internal_units = "ps", &
             description = "Set temporal cutoff for TTM (gaussian, exponential)", &
             data_type = DATA_FLOAT))

        call table%set("ttm_variable_ep", control_parameter( &
             key = "ttm_variable_ep", &
             name = "TTM Variable electron-phonon coupling", &
             val = "", &
             description = "Set electron-phonon coupling for TTM, options: homo, hetero", &
             data_type = DATA_OPTION))

        call table%set("ttm_boundary_condition", control_parameter( &
             key = "ttm_boundary_condition", &
             name = "TTM Boundary conditions", &
             val = "", &
             description = "Set boundary conditions for TTM, options: periodic, dirichlet, neumann, robin", &
             data_type = DATA_OPTION))

        call table%set("ttm_boundary_xy", control_parameter( &
             key = "ttm_boundary_xy", &
             name = "TTM Neumann Boundary in Z", &
             val = "off", &
             description = "Fix Neumann (zero-flux) boundary in Z", &
             data_type = DATA_BOOL))

        call table%set("ttm_boundary_heat_flux", control_parameter( &
             key = "ttm_boundary_heat_flux", &
             name = "TTM Heat flux", &
             val = "96", &
             units = "%", &
             internal_units = "", &
             description = "Set boundary heat flux in Robin boundaries for TTM", &
             data_type = DATA_BOOL))

        call table%set("ttm_time_offset", control_parameter( &
             key = "ttm_time_offset", &
             name = "TTM time offset", &
             val = "0.0", &
             units = "ps", &
             internal_units = "ps", &
             description = "Set electron-ion coupling offset for TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_oneway", control_parameter( &
             key = "ttm_oneway", &
             name = "TTM oneway", &
             val = "off", &
             description = "Enable one-way electron-phonon coupling "// &
             "when electronic temperature is greater than ionic temperature", &
             data_type = DATA_BOOL))

        call table%set("ttm_stats_frequency", control_parameter( &
             key = "ttm_stats_frequency", &
             name = "Frequency of TTM statis output", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Frequency of write to TTM PEAK_E and PEAK_I", &
             data_type = DATA_FLOAT))

        call table%set("ttm_traj_frequency", control_parameter( &
             key = "ttm_traj_frequency", &
             name = "Frequency of TTM trajectory output", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Frequency of write to TTM LATS_E and LATS_I", &
             data_type = DATA_FLOAT))

        call table%set("ttm_com_correction", control_parameter( &
             key = "ttm_com_correction", &
             name = "TTM CoM correction", &
             val = "full", &
             description = "Apply inhomogeneous Langevin thermostat to selected directions in TTM, options: full, zdir, off", &
             data_type = DATA_OPTION))

        call table%set("ttm_redistribute", control_parameter( &
             key = "ttm_redistribute", &
             name = "", &
             val = "", &
             description = "Redistribute electronic energy in newly-deactivated temperature cells to nearest active neighbours", &
             data_type = DATA_BOOL))

        call table%set("ttm_e-phonon_friction", control_parameter( &
             key = "ttm_e-phonon_friction", &
             name = "TTM Electron-phonon friction", &
             val = "0.0", &
             units = "ps^-1", &
             internal_units = "ps^-1", &
             description = "Set TTM electron-phonon friction", &
             data_type = DATA_FLOAT))

        call table%set("ttm_e-stopping_friction", control_parameter( &
             key = "ttm_e-stopping_friction", &
             name = "TTM Electron-stopping friction", &
             val = "0.0", &
             units = "ps^-1", &
             internal_units = "ps^-1", &
             description = "Set TTM electron-stopping friction", &
             data_type = DATA_FLOAT))

        call table%set("ttm_e-stopping_velocity", control_parameter( &
             key = "ttm_e-stopping_velocity", &
             name = "TTM Electron-stopping velocity", &
             val = "0.0", &
             units = "ang/ps", &
             internal_units = "internal_v", &
             description = "Set TTM electron-stopping velocity", &
             data_type = DATA_FLOAT))

      end block ttm

      integrator_tolerances: block

        call table%set("rlx_cgm_step", control_parameter( &
             key = "rlx_cgm_step", &
             name = "Relaxed shell CGM step", &
             val = "-1.0", &
             units = "ang", &
             internal_units = "internal_l", &
             description = "Set CGM stepping for relaxed shell model", &
             data_type = DATA_FLOAT))

        call table%set("rlx_tol", control_parameter( &
             key = "rlx_tol", &
             name = "Relaxed shell force tolerance", &
             val = "1.0", &
             units = "internal_f", &
             internal_units = "internal_f", &
             description = "Set force tolerance for relaxed shell model", &
             data_type = DATA_FLOAT))

        call table%set("shake_max_iter", control_parameter( &
             key = "shake_max_iter", &
             name = "Max SHAKE/RATTLE iterations", &
             val = "250", &
             description = "Set maximum number of SHAKE/RATTLE iterations", &
             data_type = DATA_INT))

        call table%set("shake_tolerance", control_parameter( &
             key = "shake_tolerance", &
             name = "SHAKE/RATTLE tolerance", &
             val = "1e-6", &
             units = "ang", &
             internal_units = "internal_l", &
             description = "Set accepted SHAKE/RATTLE tolerance", &
             data_type = DATA_FLOAT))

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
      end block integrator_tolerances

      call table%set("dftb", control_parameter( &
           key = "dftb", &
           name = "Enable DFTB", &
           val = "off", &
           description = "Enable DFTB", &
           data_type = DATA_BOOL))

      call table%set("fixed_com", control_parameter( &
           key = "fixed_com", &
           name = "Fix centre of mass", &
           val = "on", &
           description = "Remove net centre of mass momentum", &
           data_type = DATA_BOOL))

    end block simulation_properties

    equilibration_properties: block
      call table%set("reset_temperature_interval", control_parameter( &
           key = "reset_temperature_interval", &
           name = "Reset temperature time", &
           val = "-1", &
           units = "steps", &
           internal_units = "steps", &
           description = "Interval between temperature zeroing during equilibration for minimisation", &
           data_type = DATA_FLOAT))

      call table%set("regauss_frequency", control_parameter( &
           key = "regauss_frequency", &
           name = "Regauss frequency", &
           val = "-1", &
           units = "steps", &
           internal_units = "steps", &
           description = "Set the frequency of temperature regaussing", &
           data_type = DATA_FLOAT))

      call table%set("rescale_frequency", control_parameter( &
           key = "rescale_frequency", &
           name = "Rescale frequency", &
           val = "-1", &
           units = "steps", &
           internal_units = "steps", &
           description = "Set the frequency of temperature rescaling", &
           data_type = DATA_FLOAT))

      call table%set("equilibration_force_cap", control_parameter( &
           key = "equilibration_force_cap", &
           name = "Equilibration force cap", &
           val = "1000.0", &
           units = "k_b.temp/ang", &
           internal_units = "k_b.temp/ang", &
           description = "Set force cap clamping maximum force during equilibration", &
           data_type = DATA_FLOAT))

      minimisation: block
        call table%set("minimisation_criterion", control_parameter( &
             key = "minimisation_criterion", &
             name = "Minimisation criterion", &
             val = "off", &
             description = "Set minimisation criterion, options: off, force, energy, distance", &
             data_type = DATA_OPTION))

        call table%set("minimisation_tolerance", control_parameter( &
             key = "minimisation_tolerance", &
             name = "Minimisation tolerance", &
             val = "0.0", &
             description = "Set minimisation tolerance, units: determined by criterion", &
             data_type = DATA_FLOAT))

        call table%set("minimisation_step_length", control_parameter( &
             key = "minimisation_step_length", &
             name = "Minimisation step length", &
             val = "-1.0", &
             units = "ang", &
             internal_units = "internal_l", &
             description = "Set minimisation tolerance, units: unknown", &
             data_type = DATA_FLOAT))

        call table%set("minimisation_frequency", control_parameter( &
             key = "minimisation_frequency", &
             name = "Minimisation frequency", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Set minimisation frequency", &
             data_type = DATA_FLOAT))

      end block minimisation

    end block equilibration_properties

    initialisation_parameters: block
      call table%set('initial_minimum_separation', control_parameter( &
           key = 'initial_minimum_separation', &
           name = 'Initial minimum separation', &
           val = '-1.0', &
           units = 'internal_l', &
           internal_units = 'internal_l', &
           data_type = DATA_FLOAT, &
           description = 'Turn on the check on minimum separation distance between VNL pairs at re/start'))

      call table%set("restart", control_parameter( &
           key = "restart", &
           name = "Restart", &
           val = "clean", &
           description = "Set restart settings, possible options: Clean, Continue, Rescale, Noscale", &
           data_type = DATA_OPTION))

      call table%set("nfold", control_parameter( &
           key = "nfold", &
           name = "N Fold", &
           val = "1 1 1", &
           description = "Expand cell before running", &
           data_type = DATA_VECTOR3))


    end block initialisation_parameters

    forcefields: block

      global: block
        call table%set("cutoff", control_parameter( &
             key = "cutoff", &
             name = "Real-space global cutoff", &
             val = "1.0", &
             units = "internal_l", &
             internal_units = "internal_l", &
             description = "Set the global cutoff for real-speace potentials", &
             data_type = DATA_FLOAT))

        call table%set("padding", control_parameter( &
             key = "padding", &
             name = "Set VNL padding", &
             val = "0.0", &
             units = "internal_l", &
             internal_units = "internal_l", &
             description = "Set padding for sizing of Verlet neighbour lists", &
             data_type = DATA_FLOAT))
      end block global

      coulomb: block
        call table%set("coul_damping", control_parameter( &
             key = "coul_damping", &
             name = "Electrostatics Fennell damping", &
             val = "0.0", &
             units = "1/Ang", &
             internal_units = "1/internal_l", &
             description = "Calculate electrostatics using Fennell damping (Ewald-like) with given alpha", &
             data_type = DATA_FLOAT))

        call table%set("coul_dielectric_constant", control_parameter( &
             key = "coul_dielectric_constant", &
             name = "Dielectric Constant", &
             val = "1.0", &
             description = "Set dielectric constant relative to vacuum", &
             data_type = DATA_FLOAT))

        call table%set("coul_extended_exclusion", control_parameter( &
             key = "coul_extended_exclusion", &
             name = "Extended Exclusion", &
             val = "off", &
             description = "Enable extended coulombic exclusion affecting intra-molecular interactions", &
             data_type = DATA_BOOL))

        call table%set("coul_method", control_parameter( &
             key = "coul_method", &
             name = "Electrostatics method", &
             val = "off", &
             description = "Set method for electrostatics method, "// &
             "options: off, spme, dddp, pairwise, reaction_field, force_shifted", &
             data_type = DATA_OPTION))

        call table%set("coul_precision", control_parameter( &
             key = "coul_precision", &
             name = "Electrostatics Fennell precision", &
             val = "0.0", &
             description = "Calculate electrostatics using Fennell damping (Ewald-like) with given precision", &
             data_type = DATA_FLOAT))

        ewald: block
          call table%set("ewald_precision", control_parameter( &
               key = "ewald_precision", &
               name = "Ewald precision", &
               val = "1.0e-6", &
               description = "Set Ewald parameters to calculate within given precision for Ewald calculations", &
               data_type = DATA_FLOAT))

          call table%set("ewald_alpha", control_parameter( &
               key = "ewald_alpha", &
               name = "Ewald alpha", &
               val = "0.0", &
               units = "ang^-1", &
               internal_units = "internal_l^-1", &
               description = "Set real-recip changeover location for Ewald calculations", &
               data_type = DATA_FLOAT))

          call table%set("ewald_kvec", control_parameter( &
               key = "ewald_kvec", &
               name = "Ewald k-vector samples", &
               val = "0 0 0", &
               units = "", &
               internal_units = "", &
               description = "Set number of k-space samples for Ewald calculations", &
               data_type = DATA_VECTOR3))

          call table%set("ewald_kvec_spacing", control_parameter( &
               key = "ewald_kvec_spacing", &
               name = "Ewald k-vector spacing", &
               val = "0.0", &
               units = "ang^-1", &
               internal_units = "internal_l^-1", &
               description = "Calculate k-vector samples for an even sampling of given spacing in Ewald calculations", &
               data_type = DATA_FLOAT))

          call table%set("ewald_nsplines", control_parameter( &
               key = "ewald_nsplines", &
               name = "Number of B-Splines", &
               val = "8", &
               description = "Set number of B-Splines for Ewald SPME calculations, min=3", &
               data_type = DATA_INT))
        end block ewald


      end block coulomb

      polarisation: block
        call table%set("polarisation_model", control_parameter( &
             key = "polarisation_model", &
             name = "Polarisation model", &
             val = "default", &
             description = "Enable polarisation, options: default, CHARMM", &
             data_type = DATA_OPTION))

        call table%set("polarisation_thole", control_parameter( &
             key = "polarisation_thole", &
             name = "Polarisation", &
             val = "1.3", &
             description = "Set global atomic damping factor", &
             data_type = DATA_FLOAT))
      end block polarisation

      metal: block
        call table%set('metal_direct', control_parameter( &
             key = 'metal_direct', &
             name = 'Direct metallic force calculation', &
             val = 'off', &
             description = "Enable direct (non-tabulated) calculation of metallic forces", &
             data_type = DATA_BOOL))

        call table%set("metal_sqrtrho", control_parameter( &
             key = "metal_sqrtrho", &
             name = "SqrtRho metallic interpolation", &
             val = "off", &
             description = "Enable metal sqrtrho interpolation option for EAM embeding function in TABEAM", &
             data_type = DATA_BOOL))
      end block metal

      vdw: block
        call table%set("vdw_method", control_parameter( &
             key = "vdw_method", &
             name = "VdW Method", &
             val = "tabulated", &
             description = "Set method for Van der Waal's calculations, options: off, direct, tabulated, ewald", &
             data_type = DATA_OPTION))

        call table%set("vdw_cutoff", control_parameter( &
             key = "vdw_cutoff", &
             name = "VdW cutoff", &
             val = "0.0", &
             units = "internal_l", &
             internal_units = "internal_l", &
             description = "Set cut-off for Van der Waal's potentials", &
             data_type = DATA_FLOAT))

        call table%set('vdw_mix_method', control_parameter( &
             key = 'vdw_mix_method', &
             name = 'VdW mixing method', &
             val = 'off', &
             description = "Enable VdW mixing, possible mixing schemes: Off, "// &
             "Lorentz-Berthelot, Fender-Hasley, Hogervorst, Waldman-Hagler, Tang-Toennies, Functional", &
             data_type = DATA_OPTION))

        call table%set('vdw_force_shift', control_parameter( &
             key = 'vdw_force_shift', &
             name = 'VdW force shifting', &
             val = 'off', &
             description = "Enable force shift corrections to Van der Waals' forces", &
             data_type = DATA_BOOL))

      end block vdw

    end block forcefields

    plumed: block
      call table%set("plumed", control_parameter( &
           key = "plumed", &
           name = "Enable Plumed", &
           val = "off", &
           description = "Enabled plumed dynamics", &
           data_type = DATA_BOOL))

      call table%set("plumed_input", control_parameter( &
           key = "plumed_input", &
           name = "Plumed input", &
           val = "", &
           description = "Set plumed input file", &
           data_type = DATA_STRING))

      call table%set("plumed_log", control_parameter( &
           key = "plumed_log", &
           name = "Plumed log", &
           val = "", &
           description = "Set plumed log file", &
           data_type = DATA_STRING))

      call table%set("plumed_precision", control_parameter( &
           key = "plumed_precision", &
           name = "Plumed precision", &
           val = "1.0", &
           description = "Set plumed precision", &
           data_type = DATA_FLOAT))

      call table%set("plumed_restart", control_parameter( &
           key = "plumed_restart", &
           name = "Is Plumed restart", &
           val = "on", &
           description = "Restart plumed dynamics", &
           data_type = DATA_BOOL))

    end block plumed

    miscellaneous: block

      call table%set("strict_checks", control_parameter( &
           key = "strict_checks", &
           name = "Enable strict", &
           val = "on", &
           description = "Enforce strict checks such as: "// &
           "good system cutoff, particle index contiguity, disable non-error warnings, minimisation information", &
           data_type = DATA_BOOL))

      call table%set("unsafe_comms", control_parameter( &
           key = "unsafe_comms", &
           name = "Disable parallel safety", &
           val = "off", &
           description = "Do not ensure checks of logicals are enforced in parallel", &
           data_type = DATA_BOOL))


      call table%set("dftb_test", control_parameter( &
           key = "dftb_test", &
           name = "Run dftb tests", &
           val = "off", &
           description = "Do not perform a DLPOLY run, instead run dftb tests", &
           data_type = DATA_BOOL))

    end block miscellaneous

  end Subroutine initialise_control

  Function try_parse(ifile, params, comm) result(can_parse)
    !!-----------------------------------------------------------------------
    !!
    !! Attempt to detect if a control file is new or old style
    !! Only based on first keyword (usually title)
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins may 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), intent( In    ) :: params
    Type( comms_type ), Intent ( InOut ) :: comm
    Integer, Intent( In    ) :: ifile
    Logical :: line_read
    Logical :: can_parse
    Character(Len=STR_LEN) :: input, key, val, units

    do
      key = ''; val = ''; units = ''
      Call get_line(line_read,ifile,input,comm)
      if (.not. line_read) exit
      ! Trim comments
      if (scan(input, "#!") > 0) input = input(:scan(input, "#!")-1)
      call get_word(input, key)
      ! Skip blanks
      if (key == "") cycle
      Call lower_case(key)
      can_parse = params%in(key)
      exit
    end do
    rewind(ifile)

  end Function try_parse

  Subroutine parse_control_file(ifile, params, comm)
    !!-----------------------------------------------------------------------
    !!
    !! Read a control file filling the params table
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------

    Type( parameters_hash_table ), intent( InOut ) :: params
    Type( control_parameter ) :: param
    Type( comms_type ), Intent ( InOut ) :: comm
    Integer, Intent( In    ) :: ifile

    Character(Len=STR_LEN) :: input, key, val, units
    Logical :: line_read

    do
      key = ''; val = ''; units = ''
      Call get_line(line_read,ifile,input,comm)
      if (.not. line_read) exit
      ! Trim comments
      if (scan(input, "#!") > 0) input = input(:scan(input, "#!")-1)
      call get_word(input, key)
      ! Skip blanks
      if (key == "") cycle
      Call lower_case(key)
      if (.not. params%in(key)) call error(0, 'Unrecognised key '//trim(key))
      Call params%get(key, param)
      if (param%set) call error(0, 'Param '//trim(key)//' already set')
      Call read_control_param(input, param, ifile, comm)
      param%set = .true.
      call params%set(key, param)
    end do

    call params%fix()
  end Subroutine parse_control_file

  Subroutine read_control_param(input, param, ifile, comm)
    Type( control_parameter ), Intent ( InOut ) :: param
    Integer, Intent( In    ) :: ifile
    Type( comms_type ), Intent ( InOut ) :: comm
    Character(Len=*), Intent( InOut ) :: input
    Character(Len=STR_LEN) :: tmp, val
    Character(Len=MAX_KEY) :: unit
    Logical :: line_read
    Integer :: test_int, i
    Real(kind=wp) :: test_real

    select case (param%data_type)
    case (DATA_INT)
      call get_word(input, val)
      test_int = nint(word_2_real(val))
      param%val = val
    case (DATA_FLOAT)
      call get_word(input, val)
      test_real = word_2_real(val)
      param%val = val
      call get_word(input, unit)
      param%units = unit
    case (DATA_VECTOR3, DATA_VECTOR6)
      ! Handle Multiline
      i = index(input, '&')
      tmp = ''
      do while (i > 0)
        tmp = trim(tmp) // input(:i-1)
        if (input(i+1:) /= '') call error(0, 'Unexpected junk '//trim(input)//' after key '//trim(param%key))
        call get_line(line_read, ifile, input, comm)
        i = index(input, '&')
      end do
      tmp = adjustl(trim(tmp) // input)

      if (tmp(1:1) /= '[') call error(0, '')
      i = index(tmp, ']', back=.true.)

      ! Cut off braces
      input = tmp(i+1:)
      tmp = tmp(2:i-1)

      param%val = ""
      do while (tmp /= '')
        call get_word(tmp, val)
        ! Check all valid reals
        test_real = word_2_real(val)
        param%val = trim(param%val)//' '//val
      end do

      call get_word(input, unit)
      param%units = unit

    case (DATA_OPTION, DATA_BOOL)
      call get_word(input, val)
      if (val == "") call error(0, 'Missing option for '//trim(param%key))
      call lower_case(val)
      param%val = val
      input = ""
    case (DATA_STRING)
      param%val = adjustl(input)
      input = ""
    case Default
      call error(0, 'Unknown data type while parsing '//trim(param%key))
    end select

    if (input /= '') call error(0, "Unexpected junk "//trim(input)//" after key "//trim(param%key))

  End Subroutine read_control_param

  Subroutine bad_option(key, option)
    Character(Len=*), Intent(In) :: key, option

    Call info('Unrecognised option '//trim(option)//' for key '//trim(key),.true.)
    Call error(3)

  end Subroutine bad_option

End Module new_control
