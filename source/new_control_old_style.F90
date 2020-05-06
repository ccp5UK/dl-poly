module new_control_old_style
  ! Temporary module to allow testing and direct comparison with old style control parsing

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
  Use hash, only: MAX_KEY, STR_LEN
  Use parse, only: get_word, get_line, lower_case
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
  Use control_parameter_module, only : parameters_hash_table, control_parameter, &
       DATA_INT, DATA_FLOAT, DATA_STRING, DATA_BOOL, DATA_OPTION, DATA_VECTOR3, DATA_VECTOR6
  Use new_control, only : bad_option, read_ensemble, read_structure_analysis, &
       read_bond_analysis, write_ensemble
  Implicit None

  Private
  Public :: setup_file_io
  Public :: scan_new_control_pre_old
  Public :: scan_new_control_output_old
  Public :: scan_new_control_old
  Public :: read_new_control_old

contains

  Subroutine scan_new_control_pre_old(params, image_convention, density_variance)
    !!-----------------------------------------------------------------------
    !!
    !! Scan control for use in set_bounds
    !! This function is stupid.
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), intent( In    ) :: params
    Integer, Intent( InOut ) :: image_convention
    Real(kind=wp), Intent(   Out) :: density_variance
    Logical :: ltmp

    call params%retrieve('slab', ltmp)
    call params%retrieve('density_variance', density_variance)

    if (ltmp) image_convention = 6

  end Subroutine scan_new_control_pre_old

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

  Subroutine scan_new_control_output_old(params, files)
    Type( parameters_hash_table ), intent( In    ) :: params
    Type( file_type ), Intent( InOut ), Dimension(:) :: files

    Call params%retrieve('io_file_output', files(FILE_OUTPUT)%filename)
    Call params%retrieve('io_file_config', files(FILE_CONFIG)%filename)
    Call params%retrieve('io_file_field',  files(FILE_FIELD)%filename)
    Call params%retrieve('io_file_statis',  files(FILE_STATS)%filename)
    Call params%retrieve('io_file_history',files(FILE_HISTORY)%filename)
    Call params%retrieve('io_file_historf',files(FILE_HISTORF)%filename)
    Call params%retrieve('io_file_revive', files(FILE_REVIVE)%filename)
    Call params%retrieve('io_file_revcon', files(FILE_REVCON)%filename)
    Call params%retrieve('io_file_revold', files(FILE_REVOLD)%filename)

  end Subroutine scan_new_control_output_old

  Subroutine read_new_control_old(params, lfce, impa, ttm, dfcts, rigid, &
       rsdc, cshell, cons, pmf, stats, thermo, green, devel, plume, msd_data, met, &
       pois, bond, angle, dihedral, inversion, zdensity, neigh, vdws, &
       rdf, minim, mpoles, electro, ewld, seed, traj, files, tmr, config, flow, crd, adf)
    !!-----------------------------------------------------------------------
    !!
    !! Read remaining contents of control file -- DEPRECATED
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------

    Type( parameters_hash_table ), Intent( In    ) :: params
    Logical,                  Intent(  Out) :: lfce
    Type(impact_type),        Intent(  Out) :: impa
    Type(ttm_type),           Intent(InOut) :: ttm
    Type(defects_type),       Intent(InOut) :: dfcts(:)
    Type(rigid_bodies_type),  Intent(InOut) :: rigid
    Type(rsd_type),           Intent(InOut) :: rsdc
    Type(core_shell_type),    Intent(InOut) :: cshell
    Type(constraints_type),   Intent(InOut) :: cons
    Type(pmf_type),           Intent(InOut) :: pmf
    Type(stats_type),         Intent(InOut) :: stats
    Type(thermostat_type),    Intent(InOut) :: thermo
    Type(greenkubo_type),     Intent(InOut) :: green
    Type(development_type),   Intent(InOut) :: devel
    Type(plumed_type),        Intent(InOut) :: plume
    Type(msd_type),           Intent(InOut) :: msd_data
    Type(metal_type),         Intent(InOut) :: met
    Type(poisson_type),       Intent(InOut) :: pois
    Type(bonds_type),         Intent(InOut) :: bond
    Type(angles_type),        Intent(In   ) :: angle
    Type(dihedrals_type),     Intent(In   ) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(minimise_type),      Intent(InOut) :: minim
    Type(mpole_type),         Intent(InOut) :: mpoles
    Type(electrostatic_type), Intent(InOut) :: electro
    Type(ewald_type),         Intent(InOut) :: ewld
    Type(seed_type),          Intent(InOut) :: seed
    Type(trajectory_type),    Intent(InOut) :: traj
    Type(file_type),          Intent(InOut) :: files(:)
    Type(timer_type),         Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(flow_type),          Intent(InOut) :: flow
    Type(coord_type),         Intent(InOut) :: crd
    Type(adf_type),           Intent(InOut) :: adf

    Logical :: limp, lforc
    Character(Len=256) :: message, messages(9)
    Character(Len=80) :: banner(9)

    Character(Len=STR_LEN) :: option
    Type(control_parameter) :: param
    Real(kind=wp), dimension(6) :: vtmp
    Real(kind=wp) :: rtmp, tol
    Integer :: itmp
    Logical :: ltmp

    call params%retrieve('title', option)
    config%sysname = option(1:72)
    Write (banner(1), '(a)') ''
    Write (banner(2), '(a)') Repeat('*', 80)
    Write (banner(3), '(a4,a72,a4)') '*** ', 'title:'//Repeat(' ', 66), ' ***'
    Write (banner(4), '(a4,a72,a4)') '*** ', config%sysname, ' ***'
    Write (banner(5), '(a)') Repeat('*', 80)
    Write (banner(6), '(a)') ''
    Call info(banner, 6, .true.)

    Call info('simulation control parameters', .true.)


    call params%retrieve('output_energy', devel%l_eng)
    call params%retrieve('io_write_ascii_revive', devel%l_rout)
    call params%retrieve('io_read_ascii_revold', devel%l_rin)

    devel%l_dis = params%is_set('initial_minimum_separation')
    call params%retrieve('initial_minimum_separation', devel%r_dis)

    call params%retrieve('unit_test', devel%run_unit_tests)
    if (devel%run_unit_tests) then
       devel%unit_test%configuration = .true.
    end if

    call params%retrieve('vdw_method', option)
    select case(option)
    case ('tabulated')
       continue
    case ('direct')
       vdws%l_direct = .true.
    case ('off')
      Call info('vdw potential terms switched off', .true.)
    case ('ewald')
       ! Add when merged
       !  ewld%vdw = .true.
    case default
       call bad_option('vdw_method', option)
    end select

    Write (message, '(a,1p,e12.4)') 'real space cutoff (Angs) ', neigh%cutoff
    Call info(message, .true.)
    Write (message, '(a,1p,e12.4)') 'cutoff padding (Angs) ', neigh%padding
    Call info(message, .true.)

    Write (message, '(a,1p,e12.4)') 'vdw cutoff (Angs) ', vdws%cutoff
    call info(message, .true.)

    call params%retrieve('nfold', vtmp(1:3))
    config%l_exp = any(nint(vtmp(1:3)) > 1)
    if (config%l_exp) then
       config%nx = nint(vtmp(1))
       config%ny = nint(vtmp(2))
       config%nz = nint(vtmp(3))
       Write (message, '(a,9x,3i5)') 'system expansion opted', config%nx, config%ny, config%nz
       Call info(message, .true.)
     end if

    call params%retrieve('vdw_mix_method', option)
    if (option /= 'off') then

       Call info('vdw cross terms mixing opted (for undefined mixed potentials)', .true.)
       Call info('mixing is limited to potentials of the same type only', .true.)
       Call info('mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)', .true.)

       select case (option)
       case ('off')
          continue

       case ('lorentz-bethelot')
          vdws%mixing = MIX_LORENTZ_BERTHELOT
          Call info('type of mixing selected - Lorentz-Berthelot :: e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2', .true.)

       case ('fender-hasley')
          vdws%mixing = MIX_FENDER_HASLEY
          Call info('type of mixing selected - Fender-Halsey :: e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2', .true.)

       case ('hogervorst')
          vdws%mixing = MIX_HOGERVORST
          Call info('type of mixing selected - Hogervorst (good hope) :: ' &
               //'e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)', .true.)

       case ('halgren')
          vdws%mixing = MIX_HALGREN
          Call info('type of mixing selected - Halgren HHG :: ' &
               //'e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)', .true.)

       case ('waldman-hagler')
          vdws%mixing = MIX_WALDMAN_HAGLER
          Call info('type of mixing selected - WaldmanHagler :: ' &
               //'e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ;s_ij=[(s_i^6+s_j^6)/2]^(1/6)', .true.)

       case ('tang-tonnies')
          vdws%mixing = MIX_TANG_TOENNIES
          Call info('type of mixing selected - Tang-Toennies :: ' &
               //' e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13', .true.)
          Call info(Repeat(' ', 43)//'s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)', .true.)

       case ('functional')
          vdws%mixing = MIX_FUNCTIONAL
          Call info('type of mixing selected - Functional :: ' &
               //'e_ij=3 * (e_i*e_j)^(1/2) * (s_i*s_j)^3 / ' &
               //'SUM_L=0^2{[(s_i^3+s_j^3)^2 / (4*(s_i*s_j)^L)]^(6/(6-2L))}', .true.)
          Call info(Repeat(' ', 40)//'s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}', .true.)

       case default
          call bad_option('vdw_mix_method', option)

       end select
    end if

    call params%retrieve('vdw_force_shift', vdws%l_force_shift)
    if (vdws%l_force_shift) Call info('vdw force-shifting option on', .true.)


    call params%retrieve('metal_direct', met%l_direct)
    if (met%l_direct) Call info('metal direct option on', .true.)

    call params%retrieve('metal_sqrtrho', met%l_emb)
    if (met%l_emb) Call info('metal sqrtrho option on', .true.)

    limp = params%is_any_set([Character(17) :: 'impact_part_index', 'impact_energy', 'impact_time'])
    if (limp) then
       call params%retrieve('impact_part_index', impa%imd)
       call params%retrieve('impact_time', impa%tmd)
       call params%retrieve('impact_energy', impa%emd)
       call params%retrieve('impact_direction', vtmp(1:3))
       impa%vmx = vtmp(1)
       impa%vmy = vtmp(2)
       impa%vmz = vtmp(3)

       Write (messages(1), '(a)') ''
       Write (messages(2), '(a,i10)') 'particle (index)', impa%imd
       Write (messages(3), '(a,i10)') 'timestep (steps)', impa%tmd
       Write (messages(4), '(a,1p,e12.4)') 'energy   (keV)  ', impa%emd
       Write (messages(5), '(a,1p,3e12.4)') 'v-r(x,y,z)      ', impa%vmx, impa%vmy, impa%vmz

       Call info(messages, 5, .true.)
    end if

    if (params%is_set('random_seed')) then
       call params%retrieve('random_seed', vtmp(1:3))
       call seed%init(nint(vtmp(1:3)))
       Write (message, '(a,3i5)') 'radomisation seeds supplied: ', seed%seed(1:3)
       Call info(message, .true.)
    end if

    call read_ensemble(params, thermo, ttm%l_ttm)
    call write_ensemble(thermo)
    If (thermo%l_langevin) Call langevin_allocate_arrays(thermo, config%mxatms)

    call params%retrieve('temperature', thermo%temp)
    Write (message, '(a,1p,e12.4)') 'simulation temperature (K)  ', thermo%temp
    call info(message, .true.)

    call params%retrieve('reset_temperature_interval', thermo%freq_zero)
    thermo%l_zero = thermo%freq_zero > 0
    If (thermo%freq_zero == 0) thermo%freq_zero = flow%equil_steps + 1

    if (params%num_set([Character(22) :: 'pressure_tensor', 'pressure_hydrostatic', 'pressure_perpendicular']) > 1) then
       call error(0, 'Multiple pressure specifications')
    ! else if (params%num_set([Character(22) :: 'pressure_tensor', 'pressure_hydrostatic', 'pressure_perpendicular']) == 0) then
    !    call error(0, 'No pressure specification')
    end if

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

    Call info('simulation pressure tensor (katms):', .true.)
    Write (messages(1), '(3f20.10)') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp = 1,3)
    Write (messages(2), '(3f20.10)') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp = 4,6)
    Write (messages(3), '(3f20.10)') (convert_units(thermo%stress(itmp), 'internal_p', 'katm'), itmp = 7,9)
    Call info(messages, 3, .true.)

    if (.not. thermo%anisotropic_pressure) then
       thermo%press = sum(thermo%stress(1:9:4)) / 3.0_wp
       thermo%stress = 0.0_wp
    end if

    call params%retrieve('currents_calculate', stats%cur%on)
    if (stats%cur%on) Call info("computing currents is on!", .true.)

    call params%retrieve('restart', option)
    select case (option)
    case ('rescale')
       flow%restart_key = RESTART_KEY_SCALE
       Call info('scaled restart requested (starting a new simulation)', .true.)
    case ('noscale')
       flow%restart_key = RESTART_KEY_NOSCALE
       Call info('unscaled restart requested (starting a new simulation)', .true.)
    case ('continue')
       flow%restart_key = RESTART_KEY_OLD
       Call info('restart requested (continuing an old simulation)', .true.)
       Call warning('timestep from REVOLD overides specification in CONTROL', .true.)
    case ('clean')
       flow%restart_key = RESTART_KEY_CLEAN
    case default
       call bad_option('restart', option)
    end select

    call params%retrieve('timestep', thermo%tstep, .true.)
    call params%retrieve('timestep_variable', thermo%lvar)

    If (thermo%key_dpd == DPD_NULL .and. thermo%lvar) then
        thermo%lvar = .false.
        Call warning('variable timestep unavalable in DPD themostats', .true.)
        Write (message, '(a,1p,e12.4)') 'fixed simulation timestep (ps) ', thermo%tstep
        Call info(message, .true.)
    Else If (thermo%lvar) Then
       call params%retrieve('timestep_variable_min_dist', thermo%mndis)
       call params%retrieve('timestep_variable_max_dist', thermo%mxdis)
       call params%retrieve('timestep_variable_max_delta', thermo%mxstp)

       Write (messages(1), '(a,1p,e12.4)') 'variable simulation timestep (ps) ', thermo%tstep
       Write (messages(2), '(a)') 'controls for variable timestep:'
       Write (messages(3), '(2x,a,1p,e12.4)') 'minimum distance Dmin (Angs) ', thermo%mndis
       Write (messages(4), '(2x,a,1p,e12.4)') 'maximum distance Dmax (Angs) ', thermo%mxdis
       Call info(messages, 4, .true.)

       If (thermo%mxstp > zero_plus) Then
          Write (message, '(a,1p,e12.4)') 'timestep ceiling mxstp (ps) ', thermo%mxstp
          Call info(message, .true.)
          thermo%tstep = Min(thermo%tstep, thermo%mxstp)
       Else
          thermo%mxstp = Huge(1.0_wp)
       End If
       If (thermo%mxdis < 2.5_wp * thermo%mndis .or. thermo%mndis <= 0.0_wp) Then
          Call warning(140, thermo%mndis, thermo%mxdis, 0.0_wp)
          Call error(518)
       End If
    Else
      Write (message, '(a,1p,e12.4)') 'fixed simulation timestep (ps) ', thermo%tstep
      Call info(message, .true.)
    End If

    call params%retrieve('time_run', flow%run_steps, .true.)
    Write (message, '(a,i10)') 'selected number of timesteps ', flow%run_steps
    if (flow%run_steps > 0) Call info(message, .true.)

    call params%retrieve('time_equilibration', flow%equil_steps)
    Write (message, '(a,i10)') 'equilibration period (steps) ', flow%equil_steps
    if (flow%equil_steps > 0) Call info(message, .true.)

    call params%retrieve('record_equilibration', ltmp)
    flow%equilibration = .not. ltmp
    if (ltmp) Call info('equilibration included in overall averages', .true.)

    call params%retrieve('pseudo_thermostat_method', option)
    if (option /= 'off') then

       thermo%l_stochastic_boundaries = .true.
       Call info('pseudo thermostat attached to MD cell boundary', .true.)

       select case (option)
       case ('langevin-direct')
          thermo%key_pseudo = 0
          Call info('thermostat control: Langevin + direct temperature scaling', .true.)
       case ('langevin')
          thermo%key_pseudo = 1
          Call info('thermostat control: Langevin temperature scaling', .true.)
       case ('gaussian')
          thermo%key_pseudo = 2
          Call info('thermostat control: gaussian temperature scaling', .true.)
       case ('direct')
          thermo%key_pseudo = 3
          Call info('thermostat control: direct temperature scaling', .true.)
       case default
          call bad_option('pseudo_thermostat_method', option)
       end select

       Call params%retrieve("pseudo_thermostat_width", thermo%width_pseudo)
       if (thermo%width_pseudo < 2.0_wp) then
          thermo%width_pseudo = 2.0_wp
          Call info('thermostat thickness insufficient - reset to 2 Angs', .true.)
       end if
       Write (message, '(a,1p,e12.4)') 'thermostat thickness (Angs) ', thermo%width_pseudo
       Call info(message, .true.)

       call params%retrieve("pseudo_thermostat_temperature", thermo%temp_pseudo)
       thermo%temp_pseudo = max(0.0_wp, thermo%temp_pseudo)
       Write (message, '(a,1p,e12.4)') 'thermostat temperature (K) ', thermo%temp_pseudo
       Call info(message, .true.)

    end if

    call params%retrieve('minimisation_criterion', option)
    minim%minimise = option /= 'off'

    call params%retrieve('minimisation_tolerance', minim%tolerance)
    call params%get('minimisation_tolerance', param)

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
       If (minim%tolerance < 1.0e-6_wp .or. minim%tolerance > 0.1_wp) &
            minim%tolerance = 0.005_wp

    case default
       call bad_option('minimisation_criterion', option)
    end select

    if (minim%minimise) then
       Write (messages(1), '(a)') 'minimisation option on (during equilibration)'
       Write (messages(2), '(a,a8)') 'minimisation criterion        ', trim(option)
       Write (messages(3), '(a,i10)') 'minimisation frequency (steps)', minim%freq
       Write (messages(4), '(a,1p,e12.4)') 'minimisation tolerance        ', minim%tolerance
       Write (messages(5), '(a,1p,e12.4)') 'minimisation CGM step         ', minim%step_length

       Call info(messages, 5, .true.)
    end if
    If (minim%freq == 0) minim%freq = flow%equil_steps + 1

    call params%retrieve('regauss_frequency', thermo%freq_tgaus)
    if (thermo%freq_tgaus > 0) then
       thermo%l_tgaus = .true.
       Write (message, '(a,i10)') 'temperature regaussing interval ', thermo%freq_tgaus
       Call info(message, .true.)
    end if
    If (thermo%freq_tgaus == 0) thermo%freq_tgaus = flow%equil_steps + 1

    call params%retrieve('rescale_frequency', thermo%freq_tscale)

    if (thermo%freq_tscale > 0) then
       thermo%l_tscale = .true.
       Call info('temperature scaling on (during equilibration)', .true.)
       Write (message, '(a,i10)') 'temperature scaling interval (steps)', thermo%freq_tscale
       Call info(message, .true.)
    end if
    If (thermo%freq_tscale == 0) thermo%freq_tscale = flow%equil_steps + 1

    call params%retrieve('polarisation_model', option)
    select case (option)
    case ('charmm')
       if (mpoles%max_mpoles == 0 .or. cshell%mxshl == 0) &
            call error(0, 'CHARMM polarisation selected with no mpoles or shells')

       electro%lecx = .true.
       call params%retrieve('polarisation_thole', mpoles%thole)

       Write (message, '(a,f5.2)') &
            'CHARMM polarisation scheme selected with optional atomic thole dumping of ', mpoles%thole
       Call info(message, .true.)
    case ('default')
       ! Handled in scan control
    case default
       call bad_option('polarisation_model', option)
    end select

    call params%retrieve('density_variance', rtmp)
    Write (message, '(a,1p,e12.4)') 'density variation allowance (%) ', rtmp
    Call info(message, .true.)

    call params%retrieve('coul_method', option)
    lforc = option /= 'off'

    select case (option)
    case ('off')
       electro%key = ELECTROSTATIC_NULL
       Call info('Electrostatics switched off!!!', .true.)
    case ('ewald')
       electro%key = ELECTROSTATIC_EWALD

       Call info('Electrostatics : Smooth Particle Mesh Ewald', .true.)

       if (params%is_set('ewald_precision')) then
          call params%retrieve('ewald_precision', rtmp)
          Write (message, '(a,1p,e12.4)') 'Ewald sum precision ', rtmp
          Call info(message, .true.)
       end if

       Write (messages(1), '(a,1p,e12.4)') 'Ewald convergence parameter (A^-1) ', electro%alpha
       Write (messages(2), '(a,3i5)') 'Ewald kmax1 kmax2 kmax3   (x2) ', ewld%fft_dim_a1, ewld%fft_dim_b1, ewld%fft_dim_c1
       If (ewld%fft_dim_a /= ewld%fft_dim_a1 .or. ewld%fft_dim_b /= ewld%fft_dim_b1 .or. ewld%fft_dim_c /= ewld%fft_dim_c1) Then
          Write (messages(3), '(a,3i5)') 'DaFT adjusted kmax values (x2) ', ewld%fft_dim_a, ewld%fft_dim_b, ewld%fft_dim_c
          Write (messages(4), '(a,1p,i5)') 'B-spline interpolation order ', ewld%bspline
          Call info(messages, 4, .true.)
       Else
          Write (messages(3), '(a,1p,i5)') 'B-spline interpolation order ', ewld%bspline
          Call info(messages, 3, .true.)
       end If

    case ('dddp')
       electro%key = ELECTROSTATIC_DDDP
       Call info('Electrostatics : Distance Dependent Dielectric', .true.)

    case ('pairwise')

       electro%key = ELECTROSTATIC_COULOMB
       Call info('Electrostatics : Coulombic Potential', .true.)

    case ('force_shifted')

       electro%key = ELECTROSTATIC_COULOMB_FORCE_SHIFT
       Call info('Electrostatics : Force-Shifted Coulombic Potential', .true.)

    case ('reaction_field')

       electro%key = ELECTROSTATIC_COULOMB_REACTION_FIELD
       Call info('Electrostatics : Reaction Field', .true.)

    case default

       call bad_option('coul_method', option)

    end select

    call params%retrieve('coul_dielectric_constant', electro%eps)

    if (params%is_set([Character(14) :: 'coul_damping', 'coul_precision'])) then
       call error(0, 'Both damping and precision set')

    else if (params%is_set('coul_damping')) then
       call params%retrieve('coul_damping', electro%alpha)
       Write (message, '(a,1p,e12.4)') 'damping parameter (A^-1) ', electro%alpha
       Call info(message, .true.)

    else if (params%is_set('coul_precision')) then
       call params%retrieve('coul_precision', rtmp)
       Write (message, '(a,1p,e12.4)') 'precision parameter ', rtmp
       Call info(message, .true.)
       rtmp = Max(Min(rtmp, 0.5_wp), 1.0e-20_wp)
       tol = Sqrt(Abs(Log(rtmp * neigh%cutoff)))
       electro%alpha = Sqrt(Abs(Log(rtmp * neigh%cutoff * tol))) / neigh%cutoff
       Write (message, '(a,1p,e12.4)') 'damping parameter (A^-1) derived ', electro%alpha
       Call info(message, .true.)

    end if

    If (electro%key /= ELECTROSTATIC_EWALD .and. electro%alpha > zero_plus) Then
       Call info('Fennell damping applied', .true.)
       If (neigh%cutoff < 12.0_wp) Call warning(7, neigh%cutoff, 12.0_wp, 0.0_wp)
    End If

    ! If it's been forcibly set by polarisation_model
    if (.not. electro%lecx) call params%retrieve('coul_extended_exclusion', electro%lecx)
    If (electro%key /= ELECTROSTATIC_NULL) Then
      If (electro%lecx) Then
        Call info('Extended Coulombic eXclusion : YES', .true.)
      Else
        Call info('Extended Coulombic eXclusion : NO', .true.)
      End If
    End If

    if (params%is_set('equilibration_force_cap')) then
       flow%force_cap = .true.
       call params%retrieve('equilibration_force_cap', config%fmax)
       Write (messages(1), '(a)') 'force capping on (during equilibration)'
       Write (messages(2), '(a,1p,e12.4)') 'force capping limit (kT/Angs)', config%fmax
       Call info(messages, 2, .true.)
    end if

    if (.not. flow%strict) then
       Call info('no strict option on', .true.)
       Write (messages(1), '(a)') '*** It skips printing inessential information in OUTPUT such as many       ***'
       Write (messages(2), '(a)') '*** warnings, FIELD digested information and full iteration cycles         ***'
       Write (messages(3), '(a)') '*** information from CGM based routines!  However, it also assumes some,   ***'
       Write (messages(4), '(a)') '*** deemed safe, defaults for some specified as well as unspecified by the ***'
       Write (messages(5), '(a)') '*** user options, that may or may not be needed for the simulation to run! ***'
       Write (messages(6), '(a)') '*** The defaults are deemed to deliver safer passage as well as optimal    ***'
       Write (messages(7), '(a)') '*** performance without sacrificing on accuracy!  While it may, by chance, ***'
       Write (messages(8), '(a)') '*** help to pass previously failing runs it may as well lead to a run      ***'
       Write (messages(9), '(a)') '*** failure without warnings!  Beware, avoid usage if uncertain!           ***'
       Call info(messages, 9, .true.)
    end if

    call params%retrieve('print_topology_info', flow%print_topology)
    if (.not. flow%print_topology) &
         Call info('no topology option on (avoids printing extended FIELD topology in OUTPUT)', .true.)

    call params%retrieve('vaf_averaging', green%l_average)

    call params%retrieve('fixed_com', config%l_vom)
    if (.not. config%l_vom .and. .not. ttm%l_ttm) then
       Call info('"no fixed_com" option auto-switched on - COM momentum removal will be abandoned', .true.)
       Call warning('this may lead to a build up of the COM momentum ' &
            //'and a manifestation of the "flying ice-cube" effect', .true.)
    end if

    if (cshell%mxshl > 0) then
       call params%retrieve('rlx_tol', cshell%rlx_tol(1))
       call params%retrieve('rlx_cgm_step', cshell%rlx_tol(2))
       if (cshell%rlx_tol(1) < 1.0_wp) call error(0, 'Relaxed shell CGM tolerance < 1.0')


       Write (message, '(a,1p,e12.4)') 'relaxed shell model CGM tolerance ', cshell%rlx_tol(1)
       Call info(message, .true.)
       if (cshell%rlx_tol(2) > 0.0_wp) then
          Write (message, '(a,1p,e12.4)') 'relaxed shell model CGM step ', cshell%rlx_tol(2)
          Call info(message, .true.)
       end if
    end if

    call params%retrieve('shake_max_iter', cons%max_iter_shake)
    call params%retrieve('shake_tolerance', cons%tolerance)
    If (cons%mxcons > 0 .or. pmf%mxpmf > 0) Then
      Write (messages(1), '(a,i10)') 'iterations for shake/rattle ', cons%max_iter_shake
      Write (messages(2), '(a,1p,e12.4)') 'tolerance for shake/rattle (Angs) ', cons%tolerance
      Call info(messages, 2, .true.)
    End If

    ttm_block: if (ttm%l_ttm) then

       Write (messages(1), '(a,3(1x,i8))') 'ionic temperature grid size (x,y,z):', ttm%ntsys(1:3)
       Write (messages(2), '(a,3(1x,f8.4))') 'temperature grid size (x,y,z):', ttm%delx, ttm%dely, ttm%delz
       Write (messages(3), '(a,f10.4)') 'average number of atoms/cell: ', ttm%sysrho * ttm%volume
       Call info(messages, 3, .true.)
       Write (message, '(a,3(1x,i8))') 'electronic temperature grid size (x,y,z):', &
            ttm%eltsys(1:3)
       Call info(message, .true.)

       if (ttm%ismetal) then
          Call info('electronic subsystem represents metal: thermal conductivity required', .true.)
       else
          Call info('electronic subsystem represents non-metal: thermal diffusivity required', .true.)
       end if

       call params%retrieve('ttm_dens_model', option)
       select case (option)
       case ('constant')
          ttm%ttmdyndens = .false.
          call params%retrieve('ttm_dens', ttm%cellrho, .true.)

          Write (message, '(a,f10.4)') 'user-specified atomic density (A^-3) ', ttm%cellrho
          Call info(message, .true.)

          If (ttm%cellrho <= zero_plus) Then
             ttm%cellrho = ttm%sysrho
             call error(0, 'Bad ttm_dens (<= 0)')
          Else If (ttm%cellrho > zero_plus) Then
             ttm%rcellrho = 1.0_wp / ttm%cellrho
          Else
             ttm%rcellrho = 0.0_wp
          End If

          ttm%epc_to_chi = convert_units(1.0e-12_wp * ttm%rcellrho / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')

       case ('dynamic')
          ttm%ttmdyndens = .true.
          ttm%CeType = ttm%CeType + 4
          ttm%epc_to_chi = convert_units(1.0e-12_wp / 3.0_wp, 'J.m^-3.K^-1', 'k_B/internal_l^3')
          Call info('dynamic calculations of average atomic density in active ionic cells', .true.)

       case default
          call bad_option('ttm_dens_model', option)
       end select

       select case (ttm%cetype)
       case (0, 4) !Constant
          call params%retrieve('ttm_heat_cap', ttm%ce0)
          Call info('electronic specific heat capacity set to constant value', .true.)
          Write (message, '(a,1p,e12.4)') 'electronic s.h.c. (kB/atom) ', ttm%Ce0
          Call info(message, .true.)

       case (1, 5) !tanh
          call params%retrieve('ttm_heat_cap', ttm%sh_A)
          call params%retrieve('ttm_temp_term', ttm%sh_B)

          Write (messages(1), '(a)') 'electronic specific heat capacity set to hyperbolic tangent function'
          Write (messages(2), '(a,1p,e12.4)') 'constant term A (kB/atom) ', ttm%sh_A
          Write (messages(3), '(a,1p,e12.4)') 'temperature term B (K^-1) ', ttm%sh_B
          Call info(messages, 3, .true.)

          If (ttm%sh_A <= zero_plus .or. ttm%sh_B <= zero_plus) &
               Call error(0, 'Electronic specific heat not fully specified')
          ttm%sh_B = ttm%sh_B * 1.0e-4_wp

       case (2, 6) !linear
          call params%retrieve('ttm_heat_cap', ttm%Cemax)
          call params%retrieve('ttm_fermi_temp', ttm%Tfermi)

          Write (messages(1), '(a)') 'electronic specific heat capacity set to linear function up to Fermi temperature'
          Write (messages(2), '(a,1p,e12.4)') 'max. electronic s.h.c. (kB/atom) ', ttm%Cemax
          Write (messages(3), '(a,1p,e12.4)') 'Fermi temperature (K)', ttm%Tfermi
          Call info(messages, 3, .true.)

          If (ttm%Tfermi <= zero_plus .or. ttm%Cemax <= zero_plus) &
               Call error(0, 'Electronic specific heat not fully specified')
       case (3) !tabulated
          Call info('electronic volumetric heat capacity given as tabulated function of temperature', .true.)

       end select

       If (ttm%cetype < 4) then
          ttm%sh_A = ttm%sh_A * ttm%cellrho
          ttm%Cemax = ttm%Cemax * ttm%cellrho
       end If

       select case (ttm%ketype)
       case (0) ! Infinite
          Call info('electronic thermal conductivity set to infinity', .true.)

       case (1) ! Constant
          call params%retrieve('ttm_elec_cond', ttm%ka0)
          Write (messages(1), '(a)') 'electronic thermal conductivity set to constant value'
          Write (messages(2), '(a,1p,e12.4)') 'electronic t.c. (W m^-1 K^-1) ', &
               convert_units(ttm%Ka0, 'k_b/ps/A', 'W/m/K')
          Call info(messages, 2, .true.)

          If (ttm%ka0 <= zero_plus) &
               Call error(0, 'Electronic thermal conductivity not fully specified')

       case (2) ! Drude
          call params%retrieve('ttm_elec_cond', ttm%ka0)
          Write (messages(1), '(a)') 'electronic thermal conductivity set to drude model'
          Write (messages(2), '(a,1p,e12.4)') 't.c. at system thermo%temp. (W m^-1 K^-1) ', &
               convert_units(ttm%Ka0, 'k_b/ps/A', 'W/m/K')
          Call info(messages, 2, .true.)

          If (ttm%ka0 <= zero_plus) &
               Call error(0, 'Electronic thermal conductivity not fully specified')

       case (3) ! Tabulated
          Write (messages(1), '(a)') 'electronic thermal conductivity given as tabulated function of temperature:'
          Write (messages(2), '(a)') 'uses ionic or system temperature to calculate cell conductivity value'
          Write (messages(3), '(a)') 'for thermal diffusion equation'
          Call info(messages, 3, .true.)

       end select

       select case (ttm%detype)
       case (0) ! Off
          continue
       case (1) ! Constant
          call params%retrieve('ttm_diff', ttm%diff0)
          Write (messages(1), '(a)') 'electronic thermal diffusivity set to constant value'
          Write (messages(2), '(a,1p,e12.4)') 'electronic t.d. (m^2 s^-1) ', convert_units(ttm%diff0, 'ang^2/ps', 'm^2/s')
          Call info(messages, 2, .true.)

          if (ttm%diff0 <= zero_plus) &
               Call error(0, 'Thermoal diffusivity of non-metal not specified')

       case (2) ! Recip

          call params%retrieve('ttm_diff', ttm%diff0)
          call params%retrieve('ttm_fermi_temp', ttm%Tfermi)

          Write (messages(1), '(a)') 'electronic thermal diffusivity set to reciprocal function up to Fermi temperature'
          Write (messages(2), '(a,1p,e12.4)') 'datum electronic t.d. (m^2 s^-1) ', convert_units(ttm%diff0, 'ang^2/ps', 'm^2/s')
          Write (messages(3), '(a,1p,e12.4)') 'Fermi temperature (K) ', ttm%Tfermi
          Call info(messages, 3, .true.)

          ! Scaled by system temperature
          ttm%diff0 = thermo%temp*ttm%diff0

       case (3) ! Tabulated

          Call info('electronic thermal diffusivity given as tabulated function of temperature', .true.)

       end select


       call params%retrieve('ttm_min_atoms', ttm%amin)
       Write (message, '(a,1p,i8)') 'min. atom no. for ionic cells ', ttm%amin
       Call info(message, .true.)

       call params%retrieve('ttm_redistribute', ttm%redistribute)
       if (ttm%redistribute) then
          Write (messages(1), '(a)') 'redistributing energy from deactivated electronic cells into active neighbours'
          Write (messages(2), '(a)') '(requires at least one electronic temperature cell beyond ionic cells)'
          Call info(messages, 2, .true.)
       End If

       call params%retrieve('ttm_stopping_power', ttm%dedx)
       Write (message, '(a,1p,e12.4)') 'elec. stopping power (eV/nm) ', convert_units(ttm%dedx, 'e.V/ang', 'e.V/nm')
       Call info(message, .true.)

       call params%retrieve('ttm_spatial_dist', option)
       select case (option)
       case ('gaussian')
          ttm%sdepoType = 1
          call params%retrieve('ttm_spatial_sigma', ttm%sig)
          call params%retrieve('ttm_spatial_cutoff', ttm%sigmax)
          Write (messages(1), '(a)') 'initial gaussian spatial energy deposition in electronic system'
          Write (messages(2), '(a,1p,e12.4)') 'thermo%ttm%sigma of distribution (nm) ', convert_units(ttm%sig, 'ang', 'nm')
          Write (messages(3), '(a,1p,e12.4)') 'distribution cutoff (nm) ', convert_units(ttm%sigmax * ttm%sig, 'ang', 'nm')
          Call info(messages, 3, .true.)

       case ('flat')
          ttm%sdepoType = 2
          Call info('initial homogeneous (flat) spatial energy deposition in electronic system', .true.)

       case ('laser')
          call params%retrieve('ttm_laser_type', option)
          call params%retrieve('ttm_fluence', ttm%fluence)
          call params%retrieve('ttm_penetration_depth', ttm%pdepth)

          select case (option)
          case ('flat')
             ttm%sdepoType = 2
             Write (messages(1), '(a)') 'initial homogeneous (flat) spatial ' &
                  //'energy deposition in electronic system due to laser'
             Write (messages(2), '(a,1p,e12.4)') &
                  'absorbed ttm%fluence (mJ cm^-2) ', convert_units(ttm%fluence, 'e.V/ang^2', 'mJ/cm^2')
             Write (messages(3), '(a,1p,e12.4)') 'penetration depth (nm) ', convert_units(ttm%pdepth, 'ang', 'nm')

          case ('exponential')
             ttm%sdepoType = 3

             Write (messages(1), '(a)') 'initial xy-homogeneous, z-exponential ' &
                  //'decaying spatial energy deposition in electronic system due to laser'
             Write (messages(2), '(a,1p,e12.4)') &
                  'absorbed ttm%fluence at surface (mJ cm^-2) ', convert_units(ttm%fluence, 'e.V/ang^2', 'mJ/cm^2')
             Write (messages(3), '(a,1p,e12.4)') 'penetration depth (nm) ', convert_units(ttm%pdepth, 'ang', 'nm')
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

          Write (messages(1), '(a)') 'gaussian temporal energy deposition in electronic system'
          Write (messages(2), '(a,1p,e12.4)') 'thermo%ttm%sigma of distribution (ps) ', ttm%tdepo
          Write (messages(3), '(a,1p,e12.4)') 'distribution cutoff (ps) ', 2.0_wp * ttm%tcdepo * ttm%tdepo
          Call info(messages, 3, .true.)

       case ('exponential')
          ttm%tdepotype = 2
          call params%retrieve('ttm_temporal_duration', ttm%tdepo)
          call params%retrieve('ttm_temporal_cutoff', ttm%tcdepo)

          Write (messages(1), '(a)') 'decaying exponential temporal energy deposition in electronic system'
          Write (messages(2), '(a,1p,e12.4)') 'tau of distribution (ps) ', ttm%tdepo
          Write (messages(3), '(a,1p,e12.4)') 'distribution cutoff (ps) ', ttm%tcdepo * ttm%tdepo
          Call info(messages, 3, .true.)
       case ('delta')
          ttm%tdepotype = 3
          Call info('dirac delta temporal energy deposition in electronic system', .true.)

       case ('square')
          ttm%tdepotype = 4
          call params%retrieve('ttm_temporal_duration', ttm%tdepo)
          If (ttm%tdepo <= zero_plus) Then
             ttm%tdepoType = 3
             Call info('square pulse temporal energy deposition in electronic' &
                  //'system of zero duration: being treated as dirac delta')
          Else
             Write (messages(1), '(a)') 'square pulse temporal energy deposition in electronic system'
             Write (messages(2), '(a,1p,e12.4)') 'pulse duration (ps) ', ttm%tdepo
             Call info(messages, 2, .true.)
          End If

       case default
          call bad_option('ttm_temporal_dist', option)
       end select

       Select Case (ttm%gvar)
       Case (1)
          Write (messages(1), '(a)') 'variable electron-phonon coupling values to be applied homogeneously'
       Case (2)
          Write (messages(1), '(a)') 'variable electron-phonon coupling values to be applied heterogeneously'
       End Select
       Write (messages(2), '(a)') '(overrides value given for ensemble, required tabulated stopping'
       Write (messages(3), '(a)') 'terms in g.dat file)'
       Call info(messages, 3, .true.)

       call params%retrieve('ttm_boundary_condition', option)
       call params%retrieve('ttm_boundary_xy', ltmp)
       select case (option)
       case ('periodic')
          ttm%bctypee = 1
          Call info('electronic temperature boundary conditions set as periodic', .true.)

       case ('dirichlet')
          if (ltmp) then
             ttm%bcTypeE = 4
             Write (messages(1), '(a)') 'electronic temperature boundary conditions set as dirichlet (xy), neumann (z):'
             Write (messages(2), '(a)') 'system temperature at x and y boundaries'
             Write (messages(3), '(a)') 'zero energy flux at z boundaries'
             Call info(messages, 3, .true.)
          else
             ttm%bcTypeE = 2
             Write (messages(1), '(a)') 'electronic temperature boundary conditions set as dirichlet:'
             Write (messages(2), '(a)') 'setting boundaries to system temperature'
             Call info(messages, 2, .true.)
          end if

       case ('neumann')
          ttm%bcTypeE = 3
          Write (messages(1), '(a)') 'electronic temperature boundary conditions set as neumann:'
          Write (messages(2), '(a)') 'zero energy flux at boundaries'
          Call info(messages, 2, .true.)
       case ('robin')

          call params%retrieve('ttm_boundary_heat_flux', ttm%fluxout)

          if (ltmp) then
             ttm%bcTypeE = 6
             Write (messages(1), '(a)') 'electronic temperature boundary conditions set as robin (xy), neumann (z):'
             Write (messages(2), '(a,1p,e11.4)') 'temperature leakage at x and y boundaries of ', ttm%fluxout
             Write (messages(3), '(a)') 'zero energy flux at z boundaries'
             Call info(messages, 3, .true.)
          else
             ttm%bcTypeE = 5
             Write (messages(1), '(a)') 'electronic temperature boundary conditions set as robin:'
             Write (messages(2), '(a,1p,e11.4)') 'temperature leakage at boundaries of ', ttm%fluxout
             Call info(messages, 2, .true.)
          end if

       case default
          call bad_option('ttm_boundary_condition', option)
       end select

       call params%retrieve('ttm_time_offset', ttm%ttmoffset)
       Write (message, '(a,1p,e12.4)') 'electron-ion coupling offset (ps) ', ttm%ttmoffset
       Call info(message, .true.)

       call params%retrieve('ttm_oneway', ttm%oneway)
       if (ttm%oneway) Call info('one-way electron-phonon coupling option switched on', .true.)

       call params%retrieve('ttm_stats_frequency', ttm%ttmstats)
       if (ttm%ttmstats > 0) then
          Write (messages(1), '(a)') 'ttm statistics file option on'
          Write (messages(2), '(a,i10)') 'ttm statistics file interval ', ttm%ttmstats
          Call info(messages, 2, .true.)
       end if

       call params%retrieve('ttm_traj_frequency', ttm%ttmtraj)
       if (ttm%ttmtraj > 0) then
          Write (messages(1), '(a)') 'ttm trajectory (temperature profile) file option on'
          Write (messages(2), '(a,i10)') 'ttm trajectory file interval', ttm%ttmtraj
          Call info(messages, 2, .true.)
       end if

    end if ttm_block

    ! Others set in scan
    call params%retrieve('analyse_frequency', itmp)
    call params%retrieve('analyse_frequency_bonds', flow%freq_bond)
    flow%freq_bond = max(1, itmp, flow%freq_bond)

    call params%retrieve('analyse_frequency_angles', flow%freq_angle)
    flow%freq_angle = max(1, itmp, flow%freq_angle)
    call params%retrieve('analyse_frequency_dihedrals', flow%freq_dihedral)
    flow%freq_dihedral = max(1, itmp, flow%freq_dihedral)
    call params%retrieve('analyse_frequency_inversions', flow%freq_inversion)
    flow%freq_inversion = max(1, itmp, flow%freq_inversion)

    call read_structure_analysis(params, msd_data, rdf, green, zdensity, adf, crd, traj, dfcts, rsdc)

    If (rdf%l_collect) Then
        Write (messages(1), '(a)') 'rdf collection requested:'
        Write (messages(2), '(2x,a,i10)') 'rdf collection interval ', rdf%freq
        Write (messages(3), '(2x,a,1p,e12.4)') 'rdf binsize (Angstroms) ', rdf%rbin
        Call info(messages, 3, .true.)

      If (rdf%l_print) Then
        Call info('rdf printing requested', .true.)
      Else
        If (stats%lpana) Then
          Call info('rdf printing triggered due to a PDA printing request', .true.)
          rdf%l_print = stats%lpana
        Else
          Call info('no rdf printing requested', .true.)
        End If
      End If

      If (rdf%max_rdf == 0) Then
        Call info('no rdf pairs specified in FIELD', .true.)
      Else
        Call info('rdf pairs specified in FIELD', .true.)
      End If

      If ((.not. rdf%l_collect) .or. rdf%max_rdf == 0) Then
        Call info('rdf routines not to be activated', .true.)
        rdf%l_collect = .false.
        rdf%l_print = .false.
      End If
    End If

    If (zdensity%l_collect) Then
       Write (messages(1), '(a)') 'z-density profiles requested:'
       Write (messages(2), '(2x,a,i10)') 'z-density collection interval ', zdensity%frequency
       Write (messages(3), '(2x,a,1p,e12.4)') 'z-density binsize (Angstroms) ', rdf%rbin
       Call info(messages, 3, .true.)

       If (zdensity%l_print) Then
        Call info('z-density printing requested', .true.)
      Else
        Call info('no z-density printing requested', .true.)
      End If

      If (.not. zdensity%l_collect) Then
        Call info('z-density routines not to be activated', .true.)
        zdensity%l_print = .false.
      End If
    End If

    If (green%samp > 0) Then
       Write (messages(1), '(a)') 'vaf profiles requested:'
       Write (messages(2), '(2x,a,i10)') 'vaf collection interval ', green%freq
       Write (messages(3), '(2x,a,i10)') 'vaf binsize  ', green%binsize
       Call info(messages, 3, .true.)

      If (green%l_print) Then
        Call info('vaf printing requested', .true.)
      Else
        Call info('no vaf printing requested', .true.)
      End If

      If (green%l_average) Then
        Call info('time-averaged vaf profile', .true.)
      Else
        Call info('instantaneous vaf profiles', .true.)
      End If
    End If

    if (msd_data%l_msd) then
       Write (messages(1), '(a)') 'MSDTMP file option on'
       Write (messages(2), '(2x,a,i10)') 'MSDTMP file start ', msd_data%start
       Write (messages(3), '(2x,a,i10)') 'MSDTMP file interval ', msd_data%freq
       Call info(messages, 3, .true.)
    end if

    if (traj%ltraj) then
       Write (messages(1), '(a)') 'trajectory file option on'
       Write (messages(2), '(2x,a,i10)') 'trajectory file start ', traj%start
       Write (messages(3), '(2x,a,i10)') 'trajectory file interval ', traj%freq
       Write (messages(4), '(2x,a,i10)') 'trajectory file info key ', traj%key
       Call info(messages, 4, .true.)
    end if

    if (dfcts(1)%ldef) then
       Write (messages(1), '(a)') 'defects file option on'
       Write (messages(2), '(2x,a,i10)') 'defects file start ', dfcts(1)%nsdef
       Write (messages(3), '(2x,a,i10)') 'defects file interval ', dfcts(1)%isdef
       Write (messages(4), '(2x,a,1p,e12.4)') 'defects distance condition (Angs) ', dfcts(1)%rdef
       Call info(messages, 4, .true.)
       if (dfcts(2)%ldef) Call info('defects1 file option on', .true.)
    end if

    if (rsdc%lrsd) then
       Write (messages(1), '(a)') 'displacements file option on'
       Write (messages(2), '(2x,a,i10)') 'DISPDAT file start ', rsdc%nsrsd
       Write (messages(3), '(2x,a,i10)') 'DISPDAT file interval ', rsdc%isrsd
       Write (messages(4), '(2x,a,1p,e12.4)') 'DISPDAT distance condition (Angs) ', rsdc%rrsd
       Call info(messages, 4, .true.)
     end if

    if (params%is_set('stack_size')) then
        Write (message, '(a,i10)') 'data stacking interval (steps) ', stats%mxstak
        Call info(message, .true.)
     end if

     call params%retrieve('print_frequency', flow%freq_output)
     if (flow%freq_output > 0) then
        Write (message, '(a,i10)') 'data printing interval (steps) ', flow%freq_output
        Call info(message, .true.)
     end if

     call params%retrieve('stats_frequency', stats%intsta)
     if (stats%intsta > 0) then
        Write (message, '(a,i10)') 'statistics file interval ', stats%intsta
        Call info(message, .true.)
     end if

     call params%retrieve('time_depth', tmr%max_depth)
     call params%retrieve('time_per_mpi', tmr%proc_detail)

     call params%retrieve('data_dump_frequency', flow%freq_restart)
     Write (messages(1), '(a,i10)') 'data dumping interval (steps) ', flow%freq_restart

     call params%retrieve('subcell_threshold', neigh%pdplnc)
     neigh%pdplnc = Max(neigh%pdplnc, 1.0_wp) ! disallow any less than 1
     Write (messages(2), '(a,1p,e12.4)') 'subcelling threshold density ', neigh%pdplnc


     call params%retrieve('time_job', tmr%job)
     if (tmr%job < 0.0_wp) then
        tmr%job = Huge(1.0_wp)
        Write (messages(3), '(a)') 'Unlimited job time'
     else
        Write (messages(3), '(a,1p,e12.4)') 'allocated job run time (s) ', tmr%job
     end if

     call params%retrieve('time_close', tmr%clear_screen)
     if (tmr%clear_screen < 0.0_wp) then
        tmr%clear_screen = 0.01_wp * tmr%job
        if (tmr%job < Huge(1.0_wp) - 1.0_wp) then
           Write (messages(4), '(a,1p,e12.4)') 'allocated job close time (s) ', tmr%clear_screen
        else
           Write (messages(4), '(a,1p,e12.4)') 'Unlimited job close time'
        end if
     else
        Write (messages(4), '(a,1p,e12.4)') 'allocated job close time (s) ', tmr%clear_screen
     end if
     Call info(messages, 4, .true.)

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


    If (flow%restart_key == RESTART_KEY_CLEAN) Then
      Call info('clean start requested', .true.)
    Else If (config%levcfg == 0) Then
      Call warning(200, 0.0_wp, 0.0_wp, 0.0_wp)
      flow%restart_key = RESTART_KEY_CLEAN
    End If

    If (config%l_vom .and. ttm%l_ttm) Then
      Call info('no vom option off - COM momentum removal will be used', .true.)
      Call warning('this may lead to incorrect dynamic behaviour for' &
                   //'two-temperature model: COM momentum removal recommended', .true.)
    End If

    If (.not. flow%simulation) Then
      If (lfce) Then
        Write (messages(1), '(a)') '*** HISTORF will be replayed with full force recalculation ***'
        Write (messages(2), '(a)') '*** There is no actual dynamics/integration!!!             ***'
        Call info(messages, 2, .true.)
      Else
        Write (messages(1), '(a)') '*** HISTORY will be replayed (no actual simulation) ***'
        Write (messages(2), '(a)') '*** with structural properties will be recalculated ***'
        Call info(messages, 2, .true.)
        ! abort if there's no structural property to recalculate
        If (.not. (rdf%l_collect .or. zdensity%l_collect .or. dfcts(1)%ldef .or. &
                   msd_data%l_msd .or. rsdc%lrsd .or. (config%mxgana > 0))) Then
          Call error(580)
        End If
      End If

      If (flow%restart_key /= RESTART_KEY_CLEAN) Then
        flow%restart_key = RESTART_KEY_CLEAN ! Force clean restart
        Call info('clean start enforced', .true.)
      End If
    End If

  end Subroutine read_new_control_old

  Subroutine scan_new_control_old(params, rcter, max_rigid, imcon,imc_n,cell,xhi,yhi,zhi,mxgana, &
       l_n_r,lzdn,l_ind,nstfce,ttm,cshell,stats,thermo,green,msd_data,met, &
       pois,bond,angle,dihedral,inversion,zdensity,neigh,vdws,tersoffs,rdf,mpoles, &
       electro,ewld,kim_data,flow)
    !!-----------------------------------------------------------------------
    !!
    !! Scan the contents of control file -- DEPRECATED
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------

    Type( parameters_hash_table ), intent( In    ) :: params
    Type( ttm_type ), Intent( InOut ) :: ttm
    Logical,           Intent(   Out ) :: l_n_r,lzdn,l_ind
    Integer,           Intent( In    ) :: max_rigid,imcon
    Integer,           Intent( InOut ) :: imc_n
    Integer,           Intent(   Out ) :: mxgana,nstfce
    Real( Kind = wp ), Intent( In    ) :: xhi,yhi,zhi,rcter
    Real( Kind = wp ), Intent( InOut ) :: cell(1:9)
    Type( core_shell_type ), Intent (   In  )   :: cshell
    Type( stats_type ), Intent( InOut ) :: stats
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( greenkubo_type ), Intent( InOut ) :: green
    Type( msd_type ), Intent( InOut ) :: msd_data
    Type( metal_type ), Intent( InOut ) :: met
    Type( poisson_type ), Intent( InOut ) :: pois
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( vdw_type ), Intent( InOut ) :: vdws
    Type( tersoff_type ), Intent( In    )  :: tersoffs
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( electrostatic_type ), Intent( InOut ) :: electro
    Type( ewald_type ), Intent( InOut ) :: ewld
    Type( kim_type), Intent( InOut ) :: kim_data
    Type( flow_type ), Intent( InOut ) :: flow

    Real( kind = wp ), Parameter :: minimum_rcut = 1.0_wp, &
         minimum_bin_size = 0.05_wp, &
         minimum_bond_anal_length = 2.5_wp

    Real(kind=wp), dimension(3) :: vtmp
    Real(kind=wp), dimension(10) :: celprp
    Real(kind=wp) :: cut, eps0, fac, tol, tol1
    Character(Len=STR_LEN) :: option
    Real(Kind=wp) :: rtmp
    Logical :: ltmp
    Logical :: la_bnd,la_ang,la_dih,la_inv,lelec

    nstfce = -1

    call params%retrieve('timestep', thermo%tstep, required=.true.)
    call set_timestep(thermo%tstep)

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

    call params%retrieve('vdw_method', option)
    vdws%no_vdw = option == 'off' .or. vdws%max_vdw <= 0

    call params%retrieve('vdw_cutoff', vdws%cutoff)
    if (vdws%cutoff < minimum_rcut) then
       vdws%cutoff = neigh%cutoff
       call warning('vdw_cutoff less than global cutoff, setting to global cutoff', .true.)
    end if
    if (met%max_metal > 0 .and. met%rcut < 1.0e-6_wp) then
       met%rcut=Max(neigh%cutoff,vdws%cutoff)
       call warning('metal_cutoff not set, setting to max of global cutoff and vdw cutoff', .true.)
    end if

    call params%retrieve('ensemble_dpd_order', option)
    select case (option)
    case ('off')
       thermo%key_dpd = DPD_NULL
    case ('first', '1')
       thermo%key_dpd = DPD_FIRST_ORDER
    case ('second', '2')
       thermo%key_dpd = DPD_SECOND_ORDER
    case default
       call bad_option('ensemble_dpd_order', option)
    end select

    call params%retrieve('stack_size', stats%mxstak)

    call params%retrieve('msd_calculate', msd_data%l_msd)

    ! Velocity autocorrelation function

    call params%retrieve('vaf_calculate', green%l_collect)

    if (green%l_collect) then
       call params%retrieve('vaf_frequency', rtmp)
       call params%retrieve('vaf_binsize', green%binsize)
       call params%retrieve('vaf_print', green%l_print)

       green%freq = nint(rtmp)
       If (green%freq <= 0) green%freq=50
       If (green%binsize <= 0) green%binsize=Merge(2*green%freq,100,green%freq >= 100)

       green%samp = Ceiling(Real(green%binsize,wp)/Real(green%freq,wp))

    else if (params%is_any_set([Character(13) :: 'vaf_frequency', 'vaf_binsize', 'vaf_print'])) then
       Call warning('vaf_print, vaf_frequency or vaf_binsize found without vaf_calculate')
    end if

    call params%retrieve('zden_calculate', lzdn)

    call params%retrieve('rdf_calculate', rdf%l_collect)
    l_n_r = .not. rdf%l_collect

    call params%retrieve('rdf_binsize', rdf%rbin)
    call params%retrieve('zden_binsize', zdensity%bin_width)
    if (rdf%rbin < zero_plus) then
       rdf%rbin = minimum_bin_size
       Call warning('RDF too small, resetting to default (0.05)', .true.)
    End If


    if (params%is_set('polarisation_model')) then
       call params%retrieve('polarisation_model', option)
       select case (option)
       case ('charmm')
          if (cshell%mxshl > 0 .and. mpoles%max_mpoles > 0) mpoles%key = POLARISATION_CHARMM
       case ('default')
          mpoles%key = POLARISATION_DEFAULT
       case default
          call bad_option('polarisation_model', option)
       end select
    end if

    call params%retrieve('coul_method', option)
    lelec = option /= 'off'
    electro%no_elec = option == 'off'
    ! reinitialise multipolar electrostatics indicators

    If (electro%no_elec) Then
       mpoles%max_mpoles = 0
       mpoles%max_order = 0
       mpoles%key = POLARISATION_DEFAULT
    End If

    ! Set ewald params
    if (option == "ewald") then

       !! Old method
       call params%retrieve('ewald_nsplines', ewld%bspline)
       Call dcell(cell,celprp)

       cut = neigh%cutoff + 1e-6_wp

       if (imcon == 0) then
          cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
          cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
          cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

          cell(2) = 0.0_wp
          cell(3) = 0.0_wp
          cell(4) = 0.0_wp
          cell(6) = 0.0_wp
          cell(7) = 0.0_wp
          cell(8) = 0.0_wp
       else if (imcon == 6) then
          cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))
       End If

       if (params%is_set([Character(15) :: 'ewald_precision', 'ewald_alpha'])) then

          call error(0, 'Cannot specify both precision and manual ewald parameters')

       else if (params%is_set('ewald_alpha')) then

          call params%retrieve('ewald_alpha', electro%alpha)
          if (params%is_set([Character(18) :: 'ewald_kvec', 'ewald_kvec_spacing'])) then

             call error(0, 'Cannot specify both explicit k-vec grid and k-vec spacing')
          else if (params%is_set('ewald_kvec')) then

             call params%retrieve('ewald_kvec', vtmp(1:3))
             ewld%fft_dim_a1 = Nint(vtmp(1))
             ewld%fft_dim_b1 = Nint(vtmp(2))
             ewld%fft_dim_c1 = Nint(vtmp(3))
          else

             call params%retrieve('ewald_kvec_spacing', rtmp)
             ewld%fft_dim_a1 = Nint(rtmp / celprp(7))
             ewld%fft_dim_b1 = Nint(rtmp / celprp(8))
             ewld%fft_dim_c1 = Nint(rtmp / celprp(9))
          end if

          ! Sanity check for ill defined ewald sum parameters 1/8*2*2*2 == 1
          tol=electro%alpha*real(ewld%fft_dim_a1 * ewld%fft_dim_b1 * ewld%fft_dim_c1, wp)

          If (Int(tol) < 1) Call error(9)

       else

          call params%retrieve('ewald_precision', eps0)

          tol = Sqrt(Abs(Log(eps0*neigh%cutoff)))
          electro%alpha = Sqrt(Abs(Log(eps0*neigh%cutoff*tol)))/neigh%cutoff
          tol1 = Sqrt(-Log(eps0*neigh%cutoff*(2.0_wp*tol*electro%alpha)**2))

          fac = 1.0_wp
          If (imcon == 4 .or. imcon == 5 .or. imcon == 7) fac = 2.0_wp**(1.0_wp/3.0_wp)

          ewld%fft_dim_a1 = 2 * Nint(0.25_wp + fac * celprp(7) * electro%alpha * tol1 / pi)
          ewld%fft_dim_b1 = 2 * Nint(0.25_wp + fac * celprp(8) * electro%alpha * tol1 / pi)
          ewld%fft_dim_c1 = 2 * Nint(0.25_wp + fac * celprp(9) * electro%alpha * tol1 / pi)

       end if

       !! Remember to reset after merge

       ! ewld%active = .true.
       ! cut = neigh%cutoff + 1e-6_wp

       ! if (imcon == 0) then
       !    cell(1) = Max(2.0_wp*xhi+cut,3.0_wp*cut,cell(1))
       !    cell(5) = Max(2.0_wp*yhi+cut,3.0_wp*cut,cell(5))
       !    cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))

       !    cell(2) = 0.0_wp
       !    cell(3) = 0.0_wp
       !    cell(4) = 0.0_wp
       !    cell(6) = 0.0_wp
       !    cell(7) = 0.0_wp
       !    cell(8) = 0.0_wp
       ! else if (imcon == 6) then
       !    cell(9) = Max(2.0_wp*zhi+cut,3.0_wp*cut,cell(9))
       ! End If

       ! call params%retrieve('ewald_nsplines', ewld%bspline%num_splines)
       ! Call dcell(cell,celprp)

       ! if (params%is_set(['ewald_precision', 'ewald_alpha'])) then

       !    call error(0, 'Cannot specify both precision and manual ewald parameters')

       ! else if (params%is_set('ewald_alpha')) then

       !    call params%retrieve('ewald_alpha', ewld%alpha)
       !    if (params%is_set(['ewald_kvec', 'ewald_kvec_spacing'])) then

       !       call error(0, 'Cannot specify both explicit k-vec grid and k-vec spacing')
       !    else if (params%is_set('ewald_kvec')) then

       !       call params%retrieve('ewald_kvec', ewld%kspace%k_vec_dim_cont)
       !    else

       !       call params%retrieve('ewald_kvec_spacing', rtmp)
       !       ewld%kspace%k_vec_dim_cont = Nint(rtmp / celprp(7:9))
       !    end if

       !    ! Sanity check for ill defined ewald sum parameters 1/8*2*2*2 == 1
       !    tol=ewld%alpha*real(product(ewld%kspace%k_vec_dim_cont), wp)
       !    If (Int(tol) < 1) Call error(9)

       ! else

       !    call params%retreive('ewald_precision', eps0)

       !    tol = Sqrt(Abs(Log(eps0*neigh%cutoff)))
       !    ewld%alpha = Sqrt(Abs(Log(eps0*neigh%cutoff*tol)))/neigh%cutoff
       !    tol1 = Sqrt(-Log(eps0*neigh%cutoff*(2.0_wp*tol*ewld%alpha)**2))

       !    fac = 1.0_wp
       !    If (imcon == 4 .or. imcon == 5 .or. imcon == 7) fac = 2.0_wp**(1.0_wp/3.0_wp)

       !    ewld%kspace%k_vec_dim_cont = 2*Nint(0.25_wp + fac*celprp(7:9)*ewld%alpha*tol1/pi)

       ! end if

    end if

    call params%retrieve('ignore_config_indices', l_ind)
    if (l_ind) Call info('no index (reading in CONFIG) option on',.true.)

    call params%retrieve('strict_checks', flow%strict)

    call read_bond_analysis(params, flow, bond, angle, dihedral, inversion)

    ! call params%retrieve('analyse_bonds', la_bnd)
    ! call params%retrieve('analyse_angles', la_ang)
    ! call params%retrieve('analyse_dihedrals', la_dih)
    ! call params%retrieve('analyse_inversions', la_inv)
    ! call params%retrieve('analyse_all', ltmp)
    ! if (ltmp) then
    !    la_bnd = .true.
    !    la_ang = .true.
    !    la_dih = .true.
    !    la_inv = .true.
    ! end if

    ! if (la_bnd) then
    !    call params%retrieve('analyse_max_dist', bond%rcut)
    !    if (bond%rcut < minimum_bond_anal_length) then
    !       bond%rcut = minimum_bond_anal_length
    !       write(message, '(A, f5.2, A, f5.2)') 'Bond rcut (',bond%rcut,') less than minimum cutoff ('&
    !            , minimum_bond_anal_length, ') setting to global cutoff'
    !       call warning(message, .true.)
    !    end if
    ! end if

    ! ! Set global
    ! call params%retrieve('analyse_num_bins', bond%bin_pdf)
    ! angle%bin_adf = bond%bin_pdf
    ! dihedral%bin_adf = bond%bin_pdf
    ! inversion%bin_adf = bond%bin_pdf

    ! if (params%is_set('analyse_num_bins_bonds')) call params%retrieve('analyse_num_bins_bonds', bond%bin_pdf)
    ! if (params%is_set('analyse_num_bins_angles')) call params%retrieve('analyse_num_bins_angles', angle%bin_adf)
    ! if (params%is_set('analyse_num_bins_dihedrals')) call params%retrieve('analyse_num_bins_dihedrals', dihedral%bin_adf)
    ! if (params%is_set('analyse_num_bins_inversions')) call params%retrieve('analyse_num_bins_inversions', inversion%bin_adf)


    if (thermo%ensemble == ENS_NVT_LANGEVIN_INHOMO) then
       ttm%l_ttm = .true.

       call params%retrieve('ttm_num_ion_cells', ttm%ntsys(3))
       call params%retrieve('ttm_num_elec_cells', ttm%eltsys)
       call params%retrieve('ttm_metal', ttm%ismetal)

       call params%retrieve('ttm_heat_cap_model', option)
       select case (option)
       case ('constant')
          ttm%cetype = 0
       case ('tanh')
          ttm%cetype = 1
       case ('linear')
          ttm%cetype = 2
       case ('tabulated')
          ttm%cetype = 3
       case default
          call bad_option('ttm_heat_cap_model', option)
       end select

       if (ttm%ismetal) then
          ttm%detype = 0

          call params%retrieve('ttm_elec_cond_model', option, required=.true.)
          select case (option)
          case ('infinite')
             ttm%ketype = 0
          case ('constant')
             ttm%ketype = 1
          case ('drude')
             ttm%ketype = 2
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
          case ('recip', 'reciprocal')
             ttm%detype = 2
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

    end if

    call params%retrieve('dftb', ltmp)
    if (ltmp) flow%simulation_method = DFTB

    neigh%cutoff=Max(neigh%cutoff,vdws%cutoff,met%rcut,kim_data%cutoff,bond%rcut, 2.0_wp*rcter+1.0e-6_wp)
    ! If KIM model requires
    If (kim_data%padding_neighbours_required) Then
       If (neigh%cutoff < kim_data%influence_distance) Then
          neigh%cutoff = kim_data%influence_distance
       End If
    End If


    ! Reset vdws%cutoff, met%rcut and neigh%cutoff when only tersoff potentials are opted for
    If (tersoffs%max_ter > 0 .and. electro%no_elec .and. vdws%no_vdw .and. met%max_metal <= 0 .and. .not. rdf%l_collect) Then
       vdws%cutoff=0.0_wp
       met%rcut=0.0_wp
       If (.not.flow%strict) Then
          If (max_rigid == 0) Then ! compensate for Max(Size(RBs))>vdws%cutoff
             neigh%cutoff=2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp
          Else
             neigh%cutoff=Max(neigh%cutoff,2.0_wp*Max(bond%rcut,rcter)+1.0e-6_wp)
          End If
       End If
    End If

    ! no SPME electrostatics is specified but neigh%cutoff is still needed for
    ! domain decompositioning and link-celling
    ! It is needed for the rest of the types of electrostatic
    !! Remember to set when ewald merged
    !    If (.not. ewld%active .and. .not. lrcut .and. lelec) Call error(382)

    ! Set cutoff^2 now that cutoff won't change
    !!neigh%cutoff_2 = neigh%cutoff**2

  end Subroutine scan_new_control_old

end module new_control_old_style
