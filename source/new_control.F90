Module new_control

  Use angles,               Only: angles_type
  Use angular_distribution, Only: adf_type
  Use bonds,                Only: bonds_type
  Use comms,                Only: comms_type,&
       gcheck
  Use configuration,        Only: configuration_type
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
       warning
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
       file_type
  Use flow_control,         Only: RESTART_KEY_CLEAN,&
       RESTART_KEY_NOSCALE,&
       RESTART_KEY_OLD,&
       RESTART_KEY_SCALE,&
       flow_type, &
       DFTB, &
       MD
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
  Use langevin,             Only: langevin_allocate_arrays
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
  Use trajectory,           Only: trajectory_type
  Use ttm,                  Only: ttm_type
  Use vdw,                  Only: MIX_FENDER_HASLEY,&
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
  Use units, only : convert_units, set_timestep
  Use control_parameter_module, only : control_parameter, parameters_hash_table, &
       DATA_INT, DATA_FLOAT, DATA_STRING, DATA_BOOL, DATA_OPTION, DATA_VECTOR3, DATA_VECTOR6
  Implicit None


  !> Max number of params, increase to reduce hash collisions
  Integer, Parameter :: PARAMS_TABLE_SIZE = 500

  Public :: read_new_control
  Public :: initialise_control

  Public :: parse_file

  Public :: bad_option, read_ensemble

contains

  Subroutine read_new_control(control_file, params, comm)
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

    Call initialise_control(params)
    Call parse_file(control_file%unit_no, params, comm)

  end Subroutine read_new_control

  Subroutine setup_file_io(params, io, netcdf, files, comm)
    !!-----------------------------------------------------------------------
    !!
    !! Read in the io parameters
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), Intent( In    ) :: params
    Type( io_type ), Intent( InOut ) :: io
    Type( netcdf_param ), Intent( InOut ) :: netcdf
    Type( file_type ), Dimension(:), Intent( InOut ) :: files
    Type( comms_type ), Intent( InOut ) :: comm
    Type( control_parameter ) :: curr_param
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
       Call info('I/O read method: parallel by using MPI-I/O',.true.)

    Case ( 'direct' )
       io_read = IO_READ_DIRECT
       Call info('I/O read method: parallel by using direct access',.true.)

    Case ( 'netcdf' )
       io_read = IO_READ_NETCDF
       Call info('I/O read method: parallel by using netCDF',.true.)

    Case ( 'master' )
       io_read = IO_READ_MASTER
       Call info('I/O read method: serial by using a single master process',.true.)

    Case Default
       Call bad_option('io_read_method', curr_option)

    End Select

    Call io_set_parameters(io, user_method_read = io_read)

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
          Write(message,'(a,i10)') 'I/O readers (assumed) ',itmp
          Call info(message,.true.)
       else if (itmp < comm%mxnode) then
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O readers set to ',itmp
          Call info(message,.true.)
       else if (itmp > comm%mxnode) then
          rtmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O readers (enforced) ',itmp
          Call info(message,.true.)
       else
          Call error(0, 'Cannot have negative number of I/O readers')
       end If


       ! the number of readers is now ready to set

       Call io_set_parameters(io, user_n_io_procs_read = itmp)

       ! Sort read batch size
       ! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
       ! Note zero or negative values indicate use the default

       call params%retrieve('io_read_batch_size', itmp)

       Select Case (itmp)
       Case(0)
          Call io_get_parameters(io, user_batch_size_read = itmp )
          Write(message,'(a,i10)') 'I/O read batch size (assumed) ', itmp
          Call info(message,.true.)

       Case(1:)
          itmp = Min( itmp, MAX_BATCH_SIZE )
          Call io_set_parameters(io, user_batch_size_read = itmp )
          Write(message,'(a,i10)') 'I/O read batch size set to ', itmp
          Call info(message,.true.)

       Case Default
          Call error(0, 'Cannot have negative I/O read batch size')

       End Select

    Case(IO_READ_MASTER)
       Write(message,'(a,i10)') 'I/O readers (enforced) ', 1
       Call info(message,.true.)

    Case Default
       rtmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
       itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
       Write(message,'(a,i10)') 'I/O readers (enforced) ', itmp
       Call info(message,.true.)
       ! the number of readers is now ready to set
       Call io_set_parameters(io, user_n_io_procs_read = itmp )

    end Select

    ! Get read buffer size

    call params%retrieve('io_read_buffer_size', itmp)

    Select Case(itmp)
    Case(0)
       Call io_get_parameters(io, user_buffer_size_read = itmp )
       Write(message,'(a,i10)') 'I/O read buffer size (assumed) ',itmp
       Call info(message,.true.)

    Case(1:99)
       Call io_set_parameters(io, user_buffer_size_read = 100 )
       Write(message,'(a,i10)') 'I/O read buffer size set to ',100
       Call info(message,.true.)

    Case(100:MAX_BUFFER_SIZE)
       Call io_set_parameters(io, user_buffer_size_read = itmp )
       Write(message,'(a,i10)') 'I/O read buffer size set to ',itmp
       Call info(message,.true.)

    Case(MAX_BUFFER_SIZE+1:)
       Call io_set_parameters(io, user_buffer_size_read = MAX_BUFFER_SIZE )
       Write(message,'(a,i10)') 'I/O read buffer size set to ', MAX_BUFFER_SIZE
       Call info(message,.true.)

    Case Default
       Call error(0, 'Negative read buffer size not valid, min = 1')

    End Select

    ! Get parallel read error checking

    if (io_read /= IO_READ_MASTER) then
       call params%retrieve('io_read_error_check', ltmp)

       If (.not. ltmp) Then
          Call info('I/O parallel read error checking off',.true.)
       Else
          Call info('I/O parallel read error checking on',.true.)
       End If
       Call io_set_parameters(io, user_error_check = ltmp )
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
       Call info('I/O write method: parallel by using MPI-I/O',.true.)

    Case ( 'direct' )
       if (ltmp) then
          io_write = IO_WRITE_SORTED_DIRECT
       else
          io_write = IO_WRITE_UNSORTED_DIRECT
       end if
       Call info('I/O write method: parallel by using direct access',.true.)
       Call warning('in parallel this I/O write method has portability issues',.true.)

    Case ( 'netcdf' )
       io_write = IO_WRITE_SORTED_NETCDF

       call params%retrieve('io_read_netcdf_format', curr_option)

       Select Case (curr_option)
       Case('amber', '32bit', '32-bit')
          ! Use 32-bit quantities in output for real numbers
          Call info('I/O write method: parallel by using netCDF in the amber-like/32-bit format',.true.)
          Call io_nc_set_real_precision( sp, netcdf, itmp )
       Case('64-bit', '64bit')
          ! Use 64-bit quantities in output for real numbers
          Call info('I/O write method: parallel by using netCDF in 64-bit format',.true.)
          Call io_nc_set_real_precision( dp, netcdf, itmp )
       Case Default
          Call info('io_write_netcdf_format '//curr_option//' unrecognised option.',.true.)
          Call error(3)
       End Select

    Case ( 'master' )

       if (ltmp) then
          io_write = IO_WRITE_SORTED_MASTER
       else
          io_write = IO_WRITE_UNSORTED_MASTER
       end if

       Call info('I/O write method: serial by using a single master process',.true.)
    Case Default
       call bad_option('io_write_method', curr_option)

    End Select


    Select Case (io_write)
    Case(IO_WRITE_SORTED_MASTER, IO_WRITE_SORTED_NETCDF, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT)
       Call info('I/O write type: data sorting on',.true.)

    Case(IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT)
       Call info('I/O write type: data sorting off',.true.)

    Case Default
       call bad_option('io_write_method', curr_option)

    End Select

    ! the write method and type are now ready to set

    Call io_set_parameters(io, user_method_write = io_write )

    Select Case (io_write)
    Case ( IO_WRITE_SORTED_NETCDF, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT )
       call params%retrieve('io_write_writers', itmp)

       if (itmp == 0) then
          rtmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O writers (assumed) ',itmp
          Call info(message,.true.)

       else if (itmp < comm%mxnode) then
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O writers set to ',itmp
          Call info(message,.true.)

       else if (itmp > comm%mxnode) then
          rtmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O writers (enforced) ',itmp
          Call info(message,.true.)

       else
          Call error(0, 'Cannot have negative number of I/O writers')

       End if

       ! the number of writers is now ready to set

       Call io_set_parameters(io, user_n_io_procs_write = itmp )

       call params%retrieve('io_write_batch_size', itmp)

       Select Case (itmp)
       Case(0)
          Call io_get_parameters(io, user_batch_size_write = itmp )
          Write(message,'(a,i10)') 'I/O write batch size (assumed) ', itmp
          Call info(message,.true.)

       Case(1:)
          itmp = Min( itmp, MAX_BATCH_SIZE )
          Call io_set_parameters(io, user_batch_size_write = itmp )
          Write(message,'(a,i10)') 'I/O write batch size set to ', itmp
          Call info(message,.true.)

       Case Default
          Call error(0, 'Cannot have negative I/O write batch size')

       end Select

    Case (IO_WRITE_UNSORTED_MASTER, IO_WRITE_SORTED_MASTER)
       Write(message,'(a,i10)') 'I/O writers (enforced) ',1
       Call info(message,.true.)

    Case Default
       rtmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
       itmp = 2**Int(Nearest( Log(rtmp)/Log(2.0_wp) , +1.0_wp ))
       Write(message,'(a,i10)') 'I/O writers (enforced) ',itmp
       Call info(message,.true.)
       ! the number of writers is now ready to set
       Call io_set_parameters(io, user_n_io_procs_write = itmp )

    End Select

    call params%retrieve('io_write_buffer_size', itmp)

    Select Case(itmp)
    Case(0)
       Call io_get_parameters(io, user_buffer_size_write = itmp )
       Write(message,'(a,i10)') 'I/O write buffer size (assumed) ',itmp
       Call info(message,.true.)

    Case(1:99)
       Call io_set_parameters(io, user_buffer_size_write = 100 )
       Write(message,'(a,i10)') 'I/O write buffer size set to ',100
       Call info(message,.true.)

    Case(100:MAX_BUFFER_SIZE)
       Call io_set_parameters(io, user_buffer_size_write = itmp )
       Write(message,'(a,i10)') 'I/O write buffer size set to ',itmp
       Call info(message,.true.)

    Case(MAX_BUFFER_SIZE+1:)
       Call io_set_parameters(io, user_buffer_size_write = MAX_BUFFER_SIZE )
       Write(message,'(a,i10)') 'I/O write buffer size set to ', MAX_BUFFER_SIZE
       Call info(message,.true.)

    Case Default
       Call error(0, 'Negative write buffer size not valid, min = 1')
    end Select

    ! switch error checking flag for writing

    If (io_write /= IO_WRITE_UNSORTED_MASTER .and. io_write /= IO_WRITE_SORTED_MASTER) Then
       call params%retrieve('io_write_error_check', ltmp)
       If (.not. ltmp) Then
          Call info('I/O parallel read error checking off',.true.)
       Else
          Call info('I/O parallel read error checking on',.true.)
       End If
       Call io_set_parameters(io, user_error_check = ltmp )
    End If

    call params%retrieve('io_file_output', curr_option)
    if (curr_option /= '') Call info('OUTPUT file is '//files(FILE_OUTPUT)%filename,.true.)
    call params%retrieve('io_file_config', curr_option)
    if (curr_option /= '') Call info('CONFIG file is '//files(FILE_CONFIG)%filename,.true.)
    call params%retrieve('io_file_field', curr_option)
    if (curr_option /= '') Call info('FIELD file is '//files(FILE_FIELD)%filename,.true.)
    call params%retrieve('io_file_statis', curr_option)
    if (curr_option /= '') Call info('STATIS file is '//files(FILE_STATS)%filename,.true.)
    call params%retrieve('io_file_history', curr_option)
    if (curr_option /= '') Call info('HISTORY file is '//files(FILE_HISTORY)%filename,.true.)
    call params%retrieve('io_file_historf', curr_option)
    if (curr_option /= '') Call info('HISTORF file is '//files(FILE_HISTORF)%filename,.true.)
    call params%retrieve('io_file_revive', curr_option)
    if (curr_option /= '') Call info('REVIVE file is '//files(FILE_REVIVE)%filename,.true.)
    call params%retrieve('io_file_revcon', curr_option)
    if (curr_option /= '') Call info('REVCON file is '//files(FILE_REVCON)%filename,.true.)
    call params%retrieve('io_file_revold', curr_option)
    if (curr_option /= '') Call info('REVOLD file is '//files(FILE_REVOLD)%filename,.true.)

  End Subroutine setup_file_io

  Subroutine read_ensemble(params, thermo)
    Type( parameters_hash_table ), intent( In    ) :: params
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Character(Len=STR_LEN) :: message
    Character(Len=STR_LEN) :: option
    Character(Len=STR_LEN), dimension(4) :: messages
    Logical :: ltmp

    call params%retrieve('ensemble', option, required = .true.)

    select case(option)
    case ('nve', 'pmf')
       thermo%ensemble = ENS_NVE
       Call info('Ensemble : NVE (Microcanonical)',.true.)

    case ('nvt')

       call params%retrieve('ensemble_method', option, required = .true.)

       select case(option)
       case ('evans')
          thermo%ensemble = ENS_NVT_EVANS

          Call info('Ensemble : NVT Evans (Isokinetic)',.true.)
          Call info('Gaussian temperature constraints in use',.true.)

       case ('langevin')

          thermo%ensemble = ENS_NVT_LANGEVIN

          Call params%retrieve('ensemble_thermostat_friction', thermo%chi)

          Call info('Ensemble : NVT Langevin (Stochastic Dynamics)',.true.)
          Write(message,'(a,1p,e12.4)') 'thermostat friction (ps^-1)', thermo%chi
          Call info(message,.true.)

       case ('andersen')

          thermo%ensemble = ENS_NVT_ANDERSON

          Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          Call params%retrieve('ensemble_thermostat_softness', thermo%soft)

          Write(messages(1),'(a)') 'Ensemble : NVT Andersen'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'softness (dimensionless)',thermo%soft
          if (thermo%soft > 1) call error(0, 'Andersen softness greater than 1 '//message)

          Call info(messages,3,.true.)

       case ('berendsen')

          thermo%ensemble = ENS_NVT_BERENDSEN

          Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)

          Call info('Ensemble : NVT Berendsen',.true.)
          Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Call info(message,.true.)

       Case ('hoover', 'nose', 'nose-hoover')

          thermo%ensemble = ENS_NVT_NOSE_HOOVER

          Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)

          Call info('Ensemble : NVT Nose-Hoover',.true.)
          Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Call info(message,.true.)

       Case ('gentle', 'gst')

          thermo%ensemble = ENS_NVT_GENTLE

          Call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          Call params%retrieve('ensemble_thermostat_friction', thermo%gama)

          Write(messages(1),'(a)') 'Ensemble : NVT gentle stochastic thermostat'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'friction on thermostat  (ps^-1) ',thermo%gama
          Call info(messages,3,.true.)

       Case ('ttm')

          thermo%ensemble = ENS_NVT_LANGEVIN_INHOMO

          Call params%retrieve('ttm_e-phonon_friction', thermo%chi_ep)
          Call params%retrieve('ttm_e-stopping_friction', thermo%chi_es)
          Call params%retrieve('ttm_e-stopping_velocity', thermo%vel_es2)

          Write (messages(1), '(a)') 'Ensemble : NVT inhomogeneous Langevin (Stochastic Dynamics)'
          Write (messages(2), '(a,1p,e12.4)') 'e-phonon friction (ps^-1) ', thermo%chi_ep
          Write (messages(3), '(a,1p,e12.4)') 'e-stopping friction (ps^-1) ', thermo%chi_es
          Write (messages(4), '(a,1p,e12.4)') 'e-stopping velocity (A ps^-1) ', thermo%vel_es2
          Call info(messages, 4, .true.)

       Case ('dpd')
          thermo%ensemble = ENS_NVE

          Call info('Ensemble : NVT dpd (Dissipative Particle Dynamics)',.true.)

          call params%retrieve('ensemble_dpd_order', option)
          select case (option)
          case ('first', '1')
             thermo%key_dpd = DPD_FIRST_ORDER
             Call info("Ensemble type : Shardlow's first order splitting (S1)",.true.)
          case ('second', '2')
             thermo%key_dpd = DPD_SECOND_ORDER
             Call info("Ensemble type : Shardlow's first order splitting (S2)",.true.)
          case default
             call bad_option('ensemble_dpd_order', option)

          end select


          call params%retrieve('ensemble_dpd_drag', thermo%gamdpd(0))

          If (thermo%gamdpd(0) > zero_plus) Then
             Write(message,'(a,1p,e12.4)') 'drag coefficient (Dalton/ps) ', thermo%gamdpd(0)
             Call info(message,.true.)
          End If

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

          call params%retrieve('ensemble_thermostat_friction', thermo%chi)
          call params%retrieve('ensemble_barostat_friction', thermo%tai)

          Write(messages(1),'(a)') 'Ensemble : NPT isotropic Langevin (Stochastic Dynamics)'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat friction (ps^-1)',thermo%chi
          Write(messages(3),'(a,1p,e12.4)') 'barostat friction (ps^-1)',thermo%tai
          Call info(messages,3,.true.)

       case ('berendsen')

          thermo%ensemble = ENS_NPT_BERENDSEN

          call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

          Write(messages(1),'(a)') 'Ensemble : NPT isotropic Berendsen'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
          Call info(messages,3,.true.)

       case ('hoover', 'nose', 'nose-hoover')


          thermo%ensemble = ENS_NPT_NOSE_HOOVER

          call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

          Write(messages(1),'(a)') 'Ensemble : NPT isotropic Nose-Hoover (Melchionna)'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
          Call info(messages,3,.true.)

       case ('mtk')

          thermo%ensemble = ENS_NPT_MTK

          call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

          Write(messages(1),'(a)') 'Ensemble : NPT isotropic Martyna-Tuckerman-Klein'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
          Call info(messages,3,.true.)

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

          call params%retrieve('ensemble_thermostat_friction', thermo%chi)
          call params%retrieve('ensemble_barostat_friction', thermo%tai)

          Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Langevin (Stochastic Dynamics)'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat friction (ps^-1)',thermo%chi
          Write(messages(3),'(a,1p,e12.4)') 'barostat friction (ps^-1)',thermo%tai
          Call info(messages,3,.true.)

       case ('berendsen')

          thermo%ensemble = ENS_NPT_BERENDSEN_ANISO

          call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

          Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Berendsen'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
          Call info(messages,3,.true.)

       case ('hoover', 'nose', 'nose-hoover')

          thermo%ensemble = ENS_NPT_NOSE_HOOVER_ANISO

          call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

          Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Nose-Hoover (Melchionna)'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
          Call info(messages,3,.true.)

       Case ('mtk')

          thermo%ensemble = ENS_NPT_MTK_ANISO

          call params%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          call params%retrieve('ensemble_barostat_coupling', thermo%tau_p)

          Write(messages(1),'(a)') 'Ensemble : NPT anisotropic Martyna-Tuckerman-Klein'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'barostat relaxation time (ps) ',thermo%tau_p
          Call info(messages,3,.true.)

       case default
          call bad_option('NST ensemble method', option)
       end select

       ! Semi isotropic ensembles

       call params%retrieve('ensemble_semi_isotropic', option)
       select case (option)
       case ('off')
          continue
       case ('area')
          thermo%iso = CONSTRAINT_SURFACE_AREA
          Call info('semi-isotropic barostat : constant normal pressure (Pn) &',.true.)
          Call info('       (N-Pn-A-T)       : constant surface area (A)',.true.)
       case ('tens')

          thermo%iso = CONSTRAINT_SURFACE_TENSION

          call params%retrieve('ensemble_tension', thermo%tension, required=.true.)

          Write(messages(1),'(a)') 'semi-isotropic barostat : constant normal pressure (Pn) &'
          Write(messages(2),'(a)') '     (N-Pn-gamma-T)     : constant surface tension (gamma)'
          Write(messages(3),'(a,1p,e11.4)') 'sumulation surface tension (dyn/cm)', thermo%tension
          Call info(messages,3,.true.)
          thermo%tension=thermo%tension/tenunt

          call params%retrieve('ensemble_semi_orthorhombic', ltmp)
          if (ltmp) then
             thermo%iso = CONSTRAINT_SEMI_ORTHORHOMBIC
             Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
          end if

       case ('ortho', 'orthorhombic')

          thermo%iso = CONSTRAINT_SURFACE_TENSION
          Call info('semi-isotropic barostat : orthorhombic MD cell constraints',.true.)
          call params%retrieve('ensemble_semi_orthorhombie', ltmp)
          if (ltmp) then
             thermo%iso = CONSTRAINT_SEMI_ORTHORHOMBIC
             Call info('semi-isotropic barostat : semi-orthorhombic MD cell constraints',.true.)
          end if
       case default
          call bad_option('ensemble_semi_isotropic', option)
       end select
       If (Any(thermo%iso == [CONSTRAINT_SURFACE_AREA,CONSTRAINT_SURFACE_TENSION])) Then
          Call warning('semi-isotropic ensembles are only correct for infinite' &
               //'interfaces placed perpendicularly to the z axis',.true.)
       End If

    case default
       call bad_option('ensemble', option)
    end select

  end Subroutine read_ensemble

  Subroutine read_structure_analysis(params, msd, rdf, vaf, zden)

    Type( parameters_hash_table ), intent( In    ) :: params
    Type( msd_type ), Intent( InOut ) :: msd
    Type( greenkubo_type ), Intent( InOut ) :: vaf
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( z_density_type ), Intent( InOut ) :: zden
    Real(kind=wp) :: rtmp
    Character(Len=STR_LEN) :: option

    ! VAF
    call params%retrieve('vaf_calculate', vaf%l_collect)

    if (vaf%l_collect) then
       call params%retrieve('vaf_frequency', rtmp)
       call params%retrieve('vaf_binsize', vaf%binsize)
       call params%retrieve('vaf_print', vaf%l_print)

       vaf%freq = nint(rtmp)
       If (vaf%freq <= 0) vaf%freq=50
       If (vaf%binsize <= 0) vaf%binsize=Merge(2*vaf%freq,100,vaf%freq >= 100)

       vaf%samp = Ceiling(Real(vaf%binsize,wp)/Real(vaf%freq,wp))

    else if (params%is_any_set([Character(13) :: 'vaf_frequency', 'vaf_binsize', 'vaf_print'])) then
       Call warning('vaf_print, vaf_frequency or vaf_binsize found without vaf_calculate')
    end if

    ! RDF
    call params%retrieve('rdf_calculate', rdf%l_collect)

    if (rdf%l_collect) then
       call params%retrieve('rdf_error_analysis', option)

       select case (option)
       case ('jack')
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

       call params%retrieve('rdf_frequency', rtmp)
       call params%retrieve('rdf_binsize', rdf%rbin)
       call params%retrieve('rdf_print', rdf%l_print)

       rdf%freq = nint(rtmp)

    else if (params%is_any_set([Character(13) :: 'rdf_frequency', 'rdf_binsize', 'rdf_print'])) then
       Call warning('rdf_print, rdf_frequency or rdf_binsize found without rdf_calculate')
    end if

    ! ZDen

    call params%retrieve('zden_calculate', zden%l_collect)

    if (zden%l_collect) then

       call params%retrieve('zden_frequency', rtmp)
       call params%retrieve('zden_binsize', zden%bin_width)
       call params%retrieve('zden_print', zden%l_print)

       zden%frequency = nint(rtmp)

    else if (params%is_any_set([Character(13) :: 'zden_frequency', 'zden_binsize', 'zden_print'])) then
       Call warning('zden_print, zden_frequency or zden_binsize found without zden_calculate')
    end if


  end Subroutine read_structure_analysis

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

      call table%set("random_seed", control_parameter( &
           key = "random_seed", &
           name = "Random seed", &
           val = "1 2 3", &
           description = "Set random seed", &
           data_type = DATA_VECTOR3))

      call table%set("density_variance", control_parameter( &
           key = "density_variance", &
           name = "Expected density variance", &
           val = "1.0", &
           units = "%", &
           internal_units = "%", &
           description = "Set expected density variance for determining maximum array sizes", &
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
             val = "-1", &
             units = "hr", &
             internal_units = "s", &
             description = "Set total job time before attempted safe closure", &
             data_type = DATA_FLOAT))

        call table%set("time_close", control_parameter( &
             key = "time_close", &
             name = "Calculation close length", &
             val = "1", &
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
             description = "Set rolling average stack to n timesteps", &
             data_type = DATA_INT))

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

          call table%set("analyse_frequency", control_parameter( &
               key = "analyse_frequency", &
               name = "Analysis frequency", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Set frequency of analysis data", &
               data_type = DATA_FLOAT))

          call table%set("analyse_inversions", control_parameter( &
               key = "analyse_inversions", &
               name = "Analyse inversions", &
               val = "off", &
               description = "Enable analysis inversions", &
               data_type = DATA_BOOL))

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
               key = "msd", &
               name = "MSD calculating", &
               val = "off", &
               description = "Enable calculation of MSD", &
               data_type = DATA_BOOL))

          call table%set("msd_print", control_parameter( &
               key = "msd", &
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

          call table%set("msd_interval", control_parameter( &
               key = "msd_interval", &
               name = "MSD calculation interval", &
               val = "1", &
               units = "steps", &
               internal_units = "steps", &
               description = "Interval between dumping MSD configurations", &
               data_type = DATA_FLOAT))
        end block msd

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
               val = "off", &
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
               val = "0", &
               description = "Set number of bins to be used in RDF analysis", &
               data_type = DATA_INT))

          call table%set("rdf_error_analysis", control_parameter( &
               key = "rdf_error_analysis", &
               name = "RDF Error Analysis", &
               val = "off", &
               description = "Enable RDF error analysis, options: Off, Jack, Block", &
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
               val = "off", &
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
               val = "0", &
               description = "Set number of bins to be used in ZDen analysis", &
               data_type = DATA_INT))
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
               val = "off", &
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
             val = "ON", &
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
             description = "Set netcdf write format, options: 'amber', '32bit', '32-bit', '64-bit', '64bit'", &
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
             val = "SCREEN", &
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
           val = "ON", &
           description = "Print topology information in output file", &
           data_type = DATA_BOOL))

      call table%set("print_level", control_parameter( &
           key = "print_level", &
           name = "Print level", &
           val = "1", &
           description = "Disable unnecessary printing, levels: 0 - silent, 1 - quiet, 2 - standard, 3 - full", &
           data_type = DATA_INT))

      call table%set("timer_level", control_parameter( &
           key = "timer_level", &
           name = "Timer level", &
           val = "4", &
           description = "Do not display timers beyond so much depth", &
           data_type = DATA_INT))

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
             key = "variable_timestep_min_dist", &
             name = "Variable timestep minimum distance", &
             val = "0.03", &
             units = "internal_l", &
             internal_units = "internal_l", &
             description = "Set minimum permissible distance for variable timestep", &
             data_type = DATA_FLOAT))

        call table%set("timestep_variable_max_dist", control_parameter( &
             key = "variable_timestep_max_dist", &
             name = "Variable timestep maximum distance", &
             val = "0.1", &
             units = "internal_l", &
             internal_units = "internal_l", &
             description = "Set maximum permissible distance for variable timestep", &
             data_type = DATA_FLOAT))

        call table%set("timestep_variable_max_delta", control_parameter( &
             key = "variable_timestep_max_delta", &
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
             description = "Set ensemble constraints, options: NVE, PMF, NVT, NPT, NST, NPnAT, NPng³",&
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
             val = "", &
             description = "Set dpd method, options: first, second", &
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
             data_type = DATA_BOOL))

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
             units = "internal_p", &
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
             internal_units = "W/m/K", &
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
             internal_units = "m^2/s", &
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
             internal_units = "e.V/nm", &
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
             internal_units = "nm", &
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
             internal_units = "mJ/cm^2", &
             description = "Initial energy deposition into electronic system by laser for TTM", &
             data_type = DATA_FLOAT))

        call table%set("ttm_penetration_depth", control_parameter( &
             key = "ttm_penetration_depth", &
             name = "TTM laser penetration depth", &
             val = "0.0", &
             units = "nm", &
             internal_units = "nm", &
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
             key = "ttm_statis_frequency", &
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
             units = "A/ps", &
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
             units = "", &
             internal_units = "", &
             description = "Set accepted SHAKE/RATTLE tolerance", &
             data_type = DATA_FLOAT))

        call table%set("quaternion_max_iter", control_parameter( &
             key = "quaternion_max_iter", &
             name = "Max Leapfrog Verlet FIQA iterations", &
             val = "100", &
             description = "Set max iterations for FIQA in Leapfrog Verlet", &
             data_type = DATA_INT))

        call table%set("quaternion_tol", control_parameter( &
             key = "quaternion_tol", &
             name = "FIQA tolerance", &
             val = "1e-8", &
             units = "", &
             internal_units = "", &
             description = "Set accepted FIQA tolerance in Leapfrog Verlet", &
             data_type = DATA_FLOAT))
      end block integrator_tolerances

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
             key = "minimisation_step_length", &
             name = "Minimisation step length", &
             val = "0", &
             units = "steps", &
             internal_units = "steps", &
             description = "Set minimisation frequency", &
             data_type = DATA_FLOAT))

      end block minimisation

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
           val = "0", &
           units = "steps", &
           internal_units = "steps", &
           description = "Interval between temperature zeroing during equilibration for minimisation", &
           data_type = DATA_FLOAT))

      call table%set("regauss_frequency", control_parameter( &
           key = "regauss_frequency", &
           name = "Regauss frequency", &
           val = "0", &
           units = "steps", &
           internal_units = "steps", &
           description = "Set the frequency of temperature regaussing", &
           data_type = DATA_FLOAT))

      call table%set("rescale_frequency", control_parameter( &
           key = "rescale_frequency", &
           name = "Rescale frequency", &
           val = "0", &
           units = "steps", &
           internal_units = "steps", &
           description = "Set the frequency of temperature rescaling", &
           data_type = DATA_FLOAT))

      call table%set("equilibration_force_cap", control_parameter( &
           key = "equilibration_force_cap", &
           name = "Equilibration force cap", &
           val = "0.0", &
           units = "N", &
           internal_units = "internal_f", &
           description = "Set force cap clamping maximum force during equilibration", &
           data_type = DATA_FLOAT))
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
    end block initialisation_parameters

    forcefields: block

      global: block
        call table%set("cutoff", control_parameter( &
             key = "cutoff", &
             name = "Real-space global cutoff", &
             val = "0.0", &
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
               units = "1/Ang", &
               internal_units = "1/internal_l", &
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
               units = "1/ang", &
               internal_units = "1/internal_l", &
               description = "Calculate k-vector samples for an even sampling of given spacing in Ewald calculations", &
               data_type = DATA_FLOAT))

          call table%set("ewald_nsplines", control_parameter( &
               key = "ewald_nsplines", &
               name = "Number of B-Splines", &
               val = "12", &
               description = "Set number of B-Splines for Ewald SPME calculations", &
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
             key = "polarisation", &
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
             "Lorentz-Berthelot, Fender-Hasley/Halsey, Hogervorst, Waldman-Hagler, Tang-Toennies, Functional", &
             data_type = DATA_OPTION))

        call table%set('vdw_force_shift', control_parameter( &
             key = 'vdw_force_shift', &
             name = 'VdW force shifting', &
             val = 'off', &
             description = "Enable force shift corrections to Van der Waals' forces", &
             data_type = DATA_BOOL))

      end block vdw

    end block forcefields

    call table%set("unsafe", control_parameter( &
         key = "unsafe", &
         name = "Disable strict", &
         val = "off", &
         description = "Ignore strict checks such as; "// &
         "good system cutoff, particle âindex contiguity, disable non-error warnings, minimisation information", &
         data_type = DATA_BOOL))

    !     call table%set("disable_linked_cell", control_parameter( &
    !          key = "disable_linked_cell", &
    !          name = "", &
    !          val = "", &
    !          units = "", &
    !          internal_units = "", &
    !          description = "", &
    !          data_type = DATA_))

  end Subroutine initialise_control

  Subroutine parse_file(ifile, params, comm)
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
       call get_word(input, key)
       ! Skip comments and blanks
       if (key(1:1) == "#" .or. key(1:1) == "!" .or. key(1:1) == " ") cycle
       Call lower_case(key)
       if (.not. params%in(key)) call error(0, 'Unrecognised key '//trim(key))
       Call params%get(key, param)
       if (param%set) call error(0, 'Param '//trim(key)//' already set')
       Call read_control_param(input, param, ifile, comm)
       param%set = .true.
       call params%set(key, param)
    end do

    call params%fix()
  end Subroutine parse_file

  Subroutine read_control_param(input, param, ifile, comm)
    Type( control_parameter ), Intent ( InOut ) :: param
    Integer, Intent( In    ) :: ifile
    Type( comms_type ), Intent ( InOut ) :: comm
    Character(Len=*), Intent( InOut ) :: input
    Character(Len=STR_LEN) :: tmp, val, unit, junk
    Logical :: line_read
    Integer :: test_int, i
    Real(kind=wp) :: test_real
    Integer :: ierr

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
       tmp = tmp(:i)
       input = tmp(i+1:)

       test_int = 0
       do i = 1, len_trim(tmp)
          if (tmp(i:i) == '[') then
             test_int = test_int + 1
             tmp(i:i) = ' '
          else if (tmp(i:i) == ']') then
             test_int = test_int - 1
             tmp(i:i) = ' '
          end if
       end do
       if (test_int /= 0) call error(0, 'Unmatched brackets in input '//trim(tmp))

       do while (tmp /= '')
          call get_word(tmp, val)
          ! Check all valid reals
          test_real = word_2_real(val)
          param%val = trim(param%val)//' '//val
       end do

       call get_word(input, unit)
       param%units = unit

       ! call get_word(input, val)
       ! if (val(1:1) /= '[') call error(0, 'Expected vector in key '//trim(param%key)//' received '//trim(input))
       ! test_int = 1
       ! ! Skip [ on read
       ! if (trim(val(2:)) /= '') then
       !    test_real = word_2_real(val(2:))
       !    param%val = trim(val(2:))
       ! end if

       ! do
       !    call get_word(input, val)
       !    ! Handle multiline
       !    do while(trim(input) == '')
       !       call get_line(line_read, ifile, input, comm)
       !       if (.not. line_read) call error(0, 'End of file while parsing '//trim(param%key))
       !    end do

       !    i = index(val, ']')
       !    if (i > 1) then ! Val]
       !       test_real = word_2_real(val(1:i-1))
       !       param%val = trim(param%val) // ' ' // val(1:i-1)
       !       test_int = test_int - 1
       !    else if (i == 1) then ! Independent ]
       !       test_int = test_int - 1
       !    else
       !       test_real = word_2_real(val)
       !       param%val = trim(param%val) // ' ' // val
       !    end if
       !    if (test_int == 0) exit
       ! end do

       ! ! Handle multiline
       ! do while(trim(input) == '')
       !    call get_line(line_read, ifile, input, comm)
       !    if (.not. line_read) call error(0, 'End of file while parsing '//trim(param%key))
       ! end do
       ! call get_word(input, unit)
       ! param%units = unit

    case (DATA_OPTION, DATA_BOOL)
       call get_word(input, val)
       call lower_case(val)
       param%val = val
       input = ""
    case (DATA_STRING)
       param%val = adjustl(input)
       input = ""
    case Default
       call error(0, 'Unknown data type while parsing '//trim(param%key))
    end select

    if (trim(input) /= '') call error(0, "Unexpected junk "//trim(input)//" after key "//trim(param%key))

  End Subroutine read_control_param

  Subroutine bad_option(case, option)
    Character(Len=*), Intent(In) :: case, option

    Call info(trim(case)//' '//trim(option)//' unrecognised option.',.true.)
    Call error(3)

  end Subroutine bad_option


End Module new_control
