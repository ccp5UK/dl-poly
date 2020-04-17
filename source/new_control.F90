Module new_control

  Use comms, only: comms_type
  Use hash, only: parameters_hash_table, control_parameter, STR_LEN, &
       & DATA_INT, DATA_FLOAT, DATA_STRING, DATA_BOOL, DATA_OPTION, DATA_VECTOR3, DATA_VECTOR6
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

  Subroutine read_new_control(file, params, comm)
    Type( parameters_hash_table ), intent(   Out ) :: params
    Type( file_type ), Intent( InOut ) :: control_file
    Type( comms_type ), Intent( InOut ) :: comm

    Call initialise_control(params)
    Call parse_file(control_file%unit_no, params, comm)

  end Subroutine read_new_control

  Subroutine setup_file_io(params, io, netcdf, files, comm)
    Type( params ), Intent( In    ) :: params
    Type( io_type ), Intent( InOut ) :: io
    Type( netcdf_param ), Intent( InOut ) :: netcdf
    Type( file_type ), Dimension(:), Intent( InOut ) :: files
    Type( comms_type ), Intent( InOut ) :: comm
    Character(Len=MAX_LEN) :: curr_val, curr_unit
    Type( control_parameter ) :: curr_param
    Integer :: io_read, io_write
    Integer :: itmp
    Logical :: ltmp

    curr_val = params%get_val('io_read_method', default = 'mpiio')
    Call lower_case(curr_val)

    Select Case(Trim(curr_val%val))
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
       Call info('io_read_method '//curr_val//' unrecognised option.',.true.)
       Call error(3)

    End Select

    Call io_set_parameters(io, user_method_read = io_read)

    Select Case (io_read)
    Case (IO_READ_MPIIO, IO_READ_DIRECT, IO_READ_NETCDF)
       ! Need to calculate number of readers
       curr_val = params%get_val('io_readers', '0')
       itmp = parse_int(curr_val)

       Select Case (itmp)
       Case(0)
          tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O readers (assumed) ',itmp
          Call info(message,.true.)
       Case(1:comm%mxnode)
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O readers set to ',itmp
          Call info(message,.true.)
       Case(comm%mxnode+1:)
          tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O readers (enforced) ',itmp
          Call info(message,.true.)
       Case Default
          Call error(0, 'Cannot have negative number of I/O readers')
       End Select

       ! the number of readers is now ready to set

       Call io_set_parameters(io, user_n_io_procs_read = itmp)

       ! Sort read batch size
       ! 1 <= batch <= MAX_BATCH_SIZE, default 2000000
       ! Note zero or negative values indicate use the default

       curr_val = params%get_val('io_read_batch_size', '0')
       itmp = parse_int(curr_val)

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
       tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
       itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
       Write(message,'(a,i10)') 'I/O readers (enforced) ', itmp
       Call info(message,.true.)
       ! the number of readers is now ready to set
       Call io_set_parameters(io, user_n_io_procs_read = itmp )

    end Select

    ! Get read buffer size

    curr_val = params%get_val('io_read_buffer_size', '0')
    itmp = parse_int(curr_val)

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
       curr_val = params%get_val('io_read_error_check', 'N')
       ltmp = parse_bool(curr_val)
       If (.not. ltmp) Then
          Call info('I/O parallel read error checking off',.true.)
       Else
          Call info('I/O parallel read error checking on',.true.)
       End If
       Call io_set_parameters(io, user_error_check = l_tmp )
    End If

    ! Get write settings

    curr_val = params%get_val('io_write_method', default = 'mpiio')
    Call lower_case(curr_val)
    ltmp = parse_bool(params%get_val('io_write_sorted', default = 'ON'))

    Select Case (Trim(curr_val))
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

       curr_val = params%get_val('io_write_netcdf_format', '64bit')
       Call lower_case(curr_val)

       Select Case (Trim(curr_val))
       Case('amber', '32bit', '32-bit')
          ! Use 32-bit quantities in output for real numbers
          Call info('I/O write method: parallel by using netCDF in the amber-like/32-bit format',.true.)
          Call io_nc_set_real_precision( sp, netcdf, err_r )
       Case('64-bit', '64bit')
          ! Use 64-bit quantities in output for real numbers
          Call info('I/O write method: parallel by using netCDF in 64-bit format',.true.)
          Call io_nc_set_real_precision( dp, netcdf, err_r )
          record1=' '
          record1=word(1:Len_Trim(word)+1) // record ! back up
          record=record1
       Case Default
          Call info('io_write_netcdf_format '//curr_val//' unrecognised option.',.true.)
          Call error(3)
       End Select

    Case ( 'master' ) Then

       if (ltmp) then
          io_write = IO_WRITE_SORTED_MASTER
       else
          io_write = IO_WRITE_UNSORTED_MASTER
       end if

       Call info('I/O write method: serial by using a single master process',.true.)
    Case Default
       Call info('io_write_method '//curr_val//' unrecognised option.',.true.)
       Call error(3)

    End Select


    Select Case (io_write)
    Case(IO_WRITE_SORTED_MASTER, IO_WRITE_SORTED_NETCDF, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT)
       Call info('I/O write type: data sorting on',.true.)

    Case(IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT)
       Call info('I/O write type: data sorting off',.true.)

    Case Default
       Call info('io_write_method '//curr_val//' unrecognised option.',.true.)
       Call error(3)

    End Select

    ! the write method and type are now ready to set

    Call io_set_parameters(io, user_method_write = io_write )

    Select Case (io_write)
    Case ( IO_WRITE_SORTED_NETCDF, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT )
       curr_val = params%get_val('io_writers', default = '0')
       itmp = parse_int(curr_val)

       Select Case( itmp )
       Case(0)
          tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O writers (assumed) ',itmp
          Call info(message,.true.)

       Case(1:comm%mxnode)
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O writers set to ',itmp
          Call info(message,.true.)

       Case(comm%mxnode+1:)
          tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
          itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
          Do While ( Mod( comm%mxnode, itmp ) /= 0 )
             itmp = itmp - 1
          End Do
          Write(message,'(a,i10)') 'I/O writers (enforced) ',itmp
          Call info(message,.true.)

       Case Default
          Call error(0, 'Cannot have negative number of I/O writers')

       End Select

       ! the number of writers is now ready to set

       Call io_set_parameters(io, user_n_io_procs_write = itmp )

       curr_val = params%get_val('io_write_batch_size', default = '0')
       itmp = parse_int(curr_val)

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
       tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
       itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
       Write(message,'(a,i10)') 'I/O writers (enforced) ',itmp
       Call info(message,.true.)
       ! the number of writers is now ready to set
       Call io_set_parameters(io, user_n_io_procs_write = itmp )

    End Select

    curr_val = params%get_val('io_write_buffer_size', default = '0')
    itmp = parse_int(curr_val)

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
       curr_val = params%get_val('io_write_error_check', 'N')
       ltmp = parse_bool(curr_val)
       If (.not. ltmp) Then
          Call info('I/O parallel read error checking off',.true.)
       Else
          Call info('I/O parallel read error checking on',.true.)
       End If
       Call io_set_parameters(io, user_error_check = l_tmp )
    End If

    curr_val = params%get_val('io_output', '')
    if (curr_val /= '') Call info('OUTPUT file is '//files(FILE_OUTPUT)%filename,.true.)
    curr_val = params%get_val('io_config', '')
    if (curr_val /= '') Call info('CONFIG file is '//files(FILE_CONFIG)%filename,.true.)
    curr_val = params%get_val('io_field', '')
    if (curr_val /= '') Call info('FIELD file is '//files(FILE_FIELD)%filename,.true.)
    curr_val = params%get_val('io_statis', '')
    if (curr_val /= '') Call info('STATIS file is '//files(FILE_STATS)%filename,.true.)
    curr_val = params%get_val('io_history', '')
    if (curr_val /= '') Call info('HISTORY file is '//files(FILE_HISTORY)%filename,.true.)
    curr_val = params%get_val('io_historf', '')
    if (curr_val /= '') Call info('HISTORF file is '//files(FILE_HISTORF)%filename,.true.)
    curr_val = params%get_val('io_revive', '')
    if (curr_val /= '') Call info('REVIVE file is '//files(FILE_REVIVE)%filename,.true.)
    curr_val = params%get_val('io_revcon', '')
    if (curr_val /= '') Call info('REVCON file is '//files(FILE_REVCON)%filename,.true.)
    curr_val = params%get_val('io_revold', '')
    if (curr_val /= '') Call info('REVOLD file is '//files(FILE_REVOLD)%filename,.true.)

  End Subroutine setup_file_io

  Subroutine initialise_control(table)
    Type( parameters_hash_table ), intent(   Out ) :: table

    call table%init(256)

    call table%set('title', control_parameter( &
         key = 'title', &
         name = 'Title', &
         val = '', &
         description = "Run title", &
         data_type = DATA_STRING))

    call table%set('output_energy', control_parameter( &
         key = 'output_energy', &
         name= 'Output Final Energy', &
         val = 'OFF', &
         data_type = DATA_BOOL, &
         description = "Output final energy e_tot in output file"))

    call table%set('write_ascii_revive', control_parameter( &
         key = 'write_ascii_revive', &
         name = 'Write plain-text REVIVE', &
         val = 'OFF', &
         data_type = DATA_BOOL, &
         description = "Write REVIVE as a human-readable (ASCII) file"))

    call table%set('read_ascii_revold', control_parameter( &
         key = 'read_ascii_revold', &
         name = 'Read plain-text REVOLD', &
         val = 'OFF', &
         data_type = DATA_BOOL, &
         description = "Read human-readable (ASCII) REVOLD file"))

    call table%set('initial_minimum_separation', control_parameter( &
         key = 'initial_minimum_separation', &
         name = 'Initial minimum separation', &
         val = '-1.0', &
         units = 'internal_l', &
         internal_units = 'internal_l', &
         data_type = DATA_FLOAT, &
         description = 'Turn on the check on minimum separation distance between VNL pairs at re/start'))

    call table%set('metal_direct', control_parameter( &
         key = 'metal_direct', &
         name = 'Direct metallic force calculation', &
         val = 'OFF', &
         description = "Enable direct (non-tabulated) calculation of metallic forces", &
         data_type = DATA_BOOL))

    call table%set("metal_sqrtrho", control_parameter( &
         key = "metal_sqrtrho", &
         name = "SqrtRho metallic interpolation", &
         val = "OFF", &
         description = "Enable metal sqrtrho interpolation option for EAM embeding function in TABEAM", &
         data_type = DATA_BOOL))

    call table%set("slab", control_parameter( &
         key = "slab", &
         name = "Slab", &
         val = "OFF", &
         description = "Enable slab mechanics", &
         data_type = DATA_BOOL))

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

    call table%set("random_seed", control_parameter( &
         key = "random_seed", &
         name = "Random seed", &
         val = "1 2 3", &
         description = "Set random seed", &
         data_type = DATA_VECTOR3))

    call table%set("temperature", control_parameter( &
         key = "temperature", &
         name = "Initial/Target temperature", &
         val = "0.0", &
         units = "K", &
         internal_units = "K", &
         description = "Set the initial temperature or target temperature (for thermostats)", &
         data_type = DATA_FLOAT))

    call table%set("reset_temperature_interval", control_parameter( &
         key = "reset_temperature_interval", &
         name = "Reset temperature time", &
         val = "0", &
         units = "steps", &
         internal_units = "steps", &
         description = "Interval between temperature zeroing during equilibration for minimisation", &
         data_type = DATA_FLOAT))

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
         val = "0.0", &
         units = "katm", &
         internal_units = "internal_p", &
         description = "Set the target pressure as x, y, z perpendicular to cell faces for NPT calculations", &
         data_type = DATA_VECTOR3))

    call table%set("restart", control_parameter( &
         key = "restart", &
         name = "Restart", &
         val = "OFF", &
         description = "Set restart settings, possible options: Off, Continue, Rescale, Noscale", &
         data_type = DATA_OPTION))

    call table%set("timestep", control_parameter( &
         key = "timestep", &
         name = "Timestep", &
         val = "0.0", &
         units = "internal_t", &
         internal_units = "internal_t", &
         description = "Set calculation timestep or initial timestep for variable timestep calculations", &
         data_type = DATA_FLOAT))

    call table%set("variable_timestep", control_parameter( &
         key = "variable_timestep", &
         name = "Variable timestep", &
         val = "OFF", &
         description = "Enable variable timestep", &
         data_type = DATA_BOOL))

    call table%set("variable_timestep_min_dist", control_parameter( &
         key = "variable_timestep_min_dist", &
         name = "Variable timestep minimum distance", &
         val = "0.03", &
         units = "internal_l", &
         internal_units = "internal_l", &
         description = "Set minimum permissible distance for variable timestep", &
         data_type = DATA_FLOAT))

    call table%set("variable_timestep_max_dist", control_parameter( &
         key = "variable_timestep_max_dist", &
         name = "Variable timestep maximum distance", &
         val = "0.1", &
         units = "internal_l", &
         internal_units = "internal_l", &
         description = "Set maximum permissible distance for variable timestep", &
         data_type = DATA_FLOAT))

    call table%set("variable_timestep_max_delta", control_parameter( &
         key = "variable_timestep_max_delta", &
         name = "Variable timestep max delta", &
         val = "0.0", &
         units = "internal_t", &
         internal_units = "internal_t", &
         description = "Set maximum timestep delta for variable timestep", &
         data_type = DATA_FLOAT))

    call table%set("run_time", control_parameter( &
         key = "run_time", &
         name = "Calculation run length", &
         val = "0", &
         units = "steps", &
         internal_units = "steps", &
         description = "Set calculation run length", &
         data_type = DATA_FLOAT))

    call table%set("equilibration_time", control_parameter( &
         key = "equilibration_time", &
         name = "Equilibration run length", &
         val = "0", &
         units = "steps", &
         internal_units = "steps", &
         description = "Set equilibration run length", &
         data_type = DATA_FLOAT))

    call table%set("record_equilibration", control_parameter( &
         key = "record_equilibration", &
         name = "Record equilibration", &
         val = "OFF", &
         description = "Include equilibration in outputs", &
         data_type = DATA_BOOL))

    call table%set("pseudo_thermostat_method", control_parameter( &
         key = "pseudo_thermostat_method", &
         name = "Pseudo thermostat method", &
         val = "Langevin-Direct", &
         description = "Set pseudo thermostat method, possible options: Langevin-Direct, Langevin, Gauss, Direct", &
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
         val = "1.0", &
         units = "K", &
         internal_units = "K", &
         description = "Set the temperature of the pseudo thermostat", , &
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

    call table%set("polarisation_model", control_parameter( &
         key = "polarisation_model", &
         name = "Polarisation model", &
         val = "CHARMM", &
         description = "Enable polarisation, options: CHARMM", &
         data_type = DATA_OPTION))

    call table%set("polarisation_thole", control_parameter( &
         key = "polarisation", &
         name = "Polarisation", &
         val = "1.3", &
         description = "Set global atomic damping factor", &
         data_type = DATA_FLOAT))

    call table%set("ensemble", control_parameter( &
         key = "ensemble", &
         name = "Ensemble constraints", &
         val = "NVE", &
         description = "Set ensemble constraints, options: NVE, NVT, NPT, NST, NPnAT, NPng³", &
         data_type = DATA_OPTION))

    call table%set("ensemble_method", control_parameter( &
         key = "ensemble_method", &
         name = "Ensemble method", &
         val = "", &
         description = "Set ensemble method, options "// &
         "NVT: Evans, Langevin, Andersen, Berendsen, Hoover, gentle, ttm, dpd. "// &
         "NP|ST: Langevin, Berendsen, Hoover, MTK.", &
         data_type = DATA_OPTION))

    call table%set("ensemble_thermostat_coupling", control_parameter( &
         key = "ensemble_thermostat_coupling", &
         name = "Thermostat coupling", &
         val = "0.0", &
         units = "ps", &
         internal_units = "ps", &
         description = "Set thermostat relaxation/decorrelation times (inverse units for langevin)", &
         data_type = DATA_FLOAT))

    call table%set("ensemble_barostat_coupling", control_parameter( &
         key = "ensemble_barostat_coupling", &
         name = "Barostat coupling", &
         val = "0.0", &
         units = "ps", &
         internal_units = "ps", &
         description = "Set barostat relaxation/decorrelation times (inverse units for langevin)", &
         data_type = DATA_FLOAT))

    call table%set("ensemble_semi_isotropic", control_parameter( &
         key = "ensemble_semi_isotropic", &
         name = "Ensemble semi-isotropic constraint", &
         val = "OFF", &
         description = "Enable semi-isotropic barostat constraints", &
         data_type = DATA_BOOL))

    call table%set("ensemble_orthorhombic", control_parameter( &
         key = "ensemble_orthorhombic", &
         name = "Ensemble orthorhombic constraint", &
         val = "OFF", &
         description = "Enable orthorhombic barostat constraints", &
         data_type = DATA_BOOL))

    call table%set("ensemble_tension", control_parameter( &
         key = "ensemble_tension", &
         name = "Constrained surface tension", &
         val = "0.0", &
         units = "N/m", &
         internal_units = "internal_f/internal_l", &
         description = "Set tension in NPngT calctulation", &
         data_type = DATA_FLOAT))

    call table%set("density_variance", control_parameter( &
         key = "density_variance", &
         name = "Expected density variance", &
         val = "1.0", &
         units = "%", &
         internal_units = "%", &
         description = "Set expected density variance for determining maximum array sizes", &
         data_type = DATA_FLOAT))

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
         val = '', &
         description = "Enable VdW mixing, possible mixing schemes: "// &
         "Lorentz-Berthelot, Fender-Hasley/Halsey, Hogervorst, Waldman-Hagler, Tang-Toennies, Functional", &
         data_type = DATA_OPTION))

    call table%set('vdw_force_shift', control_parameter( &
         key = 'vdw_force_shift', &
         name = 'VdW force shifting', &
         val = 'OFF', &
         description = "Enable force shift corrections to Van der Waals' forces", &
         data_type = DATA_BOOL))

    call table%set("print_per_particle_contrib", control_parameter( &
         key = "print_per_particle_contrib", &
         name = "Per-particle contributions", &
         val = "OFF", &
         description = "Calculate and print per-particle contributions to energy, force and stress to file every stats step", &
         data_type = DATA_BOOL))

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
         description = "Set number of B-Splines for Ewald SPME calculations"))

    call table%set("coul_method", control_parameter( &
         key = "coul_method", &
         name = "Electrostatics method", &
         val = "off", &
         description = "Set method for electrostatics method, options: off, ewald, dddp, pairwise, reaction_field, force_shifted", &
         data_type = DATA_OPTION))

    call table%set("coul_damping", control_parameter( &
         key = "coul_damping", &
         name = "Electrostatics Fennell damping", &
         val = "0.0", &
         units = "1/Ang", &
         internal_units = "1/internal_l", &
         description = "Calculate electrostatics using Fennell damping (Ewald-like) with given alpha", &
         data_type = DATA_FLOAT))

    call table%set("coul_precision", control_parameter( &
         key = "coul_precision", &
         name = "Electrostatics Fennell precision", &
         val = "0.0", &
         description = "Calculate electrostatics using Fennell damping (Ewald-like) with given precision", &
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
         val = "OFF", &
         description = "Enable extended coulombic exclusion affecting intra-molecular interactions", &
         data_type = DATA_BOOL))

    call table%set("equilibration_force_cap", control_parameter( &
         key = "equilibration_force_cap", &
         name = "Equilibration force cap", &
         val = "0.0", &
         units = "N", &
         internal_units = "internal_f", &
         description = "Set force cap clamping maximum force during equilibration", &
         data_type = DATA_FLOAT))

    call table%set("ignore_config_indices", control_parameter( &
         key = "ignore_config_indices", &
         name = "Ignore Config Indices", &
         val = "OFF", &
         description = "Ignore indices as defined in CONFIG and use read order instead", &
         data_type = DATA_BOOL))

    call table%set("unsafe", control_parameter( &
         key = "unsafe", &
         name = "Disable strict", &
         val = "OFF", &
         description = "Ignore strict checks such as; "// &
         "good system cutoff, particle âindex contiguity, disable non-error warnings, minimisation information", &
         data_type = DATA_BOOL))

    call table%set("print_topology_info", control_parameter( &
         key = "print_topology_info", &
         name = "Print topology information ", &
         val = "ON", &
         description = "Print topology information in output file", &
         data_type = DATA_BOOL))

    call table%set("vaf_averaging", control_parameter( &
         key = "vaf_averaging", &
         name = "VAF Time Averaging", &
         val = "ON", &
         description = "ignore time-averaging of VAF, report all calculated VAF to VAFDAT files and final profile to OUTPUT", &
         data_type = DATA_BOOL))

    call table%set("fixed_com", control_parameter( &
         key = "fixed_com", &
         name = "Fix centre of mass", &
         val = "ON", &
         description = "Remove net centre of mass momentum", &
         data_type = DATA_BOOL))

    !     call table%set("disable_linked_cell", control_parameter( &
    !          key = "disable_linked_cell", &
    !          name = "", &
    !          val = "", &
    !          units = "", &
    !          internal_units = "", &
    !          description = "", &
    !          data_type = DATA_))

    call table%set("rlx_cgm_step", control_parameter( &
         key = "rlx_cgm_step", &
         name = "Relaxed shell CGM step", &
         val = "0", &
         units = "", &
         internal_units = "", &
         description = "Set CGM stepping for relaxed shell model", &
         data_type = DATA_INT))

    call table%set("rlx_tol", control_parameter( &
         key = "rlx_tol", &
         name = "Relaxed shell force tolerance", &
         val = "1", &
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

    call table%set("shake_tol", control_parameter( &
         key = "shake_tol", &
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

    call table%set("ttm_num_ion_cells", control_parameter( &
         key = "ttm_num_ion_cells", &
         name = "Number of TTM ionic cells", &
         val = "10", &
         description = "Set number of coarse-grained ion temperature cells (CIT)", &
         data_type = DATA_INT))

    call table%set("ttm_num_elec_cells", control_parameter( &
         key = "ttm_num_elec_cells", &
         name = "Number of TTM electronic cells", &
         val = "20 20 20", &
         description = "Set number of coarse-grained electronic temperature cells (CET)", &
         data_type = DATA_VECTOR3))

    call table%set("ttm_metal", control_parameter( &
         key = "ttm_metal", &
         name = "TTM Metallic", &
         val = "OFF", &
         description = "Specifies parameters for metallic system are required for two-temperature model, i.e. thermal conductivity", &
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

    call table%set("ttm_fermi_temp", control_parameter( &
         key = "ttm_fermi_temp", &
         name = "TTM Fermi temperature", &
         val = "0.0", &
         units = "1/K", &
         internal_units = "1/K", &
         description = "Set Fermi temperature in TTM", &
         data_type = DATA_FLOAT))

    ! call table%set("ttm_temp_gradient", control_parameter( &
    !      key = "ttm_temp_gradient", &
    !      name = "", &
    !      val = "", &
    !      units = "", &
    !      internal_units = "", &
    !      description = "", &
    !      data_type = DATA_))

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
         units = "J/s/m/K", &
         internal_units = "J/s/m/K", &
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
         units = "", &
         internal_units = "", &
         description = "", &
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
         units = "1/ang^3", &
         internal_units = "1/internal_l^3", &
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
         internal_units = "internal_e/internal_l", &
         description = "Electronic stopping power of projectile entering electronic system ", &
         data_type = DATA_FLOAT))

    call table%set("ttm_spatial_dist", control_parameter( &
         key = "ttm_spatial_dist", &
         name = "TTM Spatial deposition distribution", &
         val = "", &
         description = "Set the spatial distribution of TTM, options: flat, gaussian, laser", &
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

    call table%set("ttm_absorbed_e", control_parameter( &
         key = "ttm_absorbed_e", &
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

    call table%set("ttm_deposition_type", control_parameter( &
         key = "ttm_deposition_type", &
         name = "TTM Deposition type", &
         val = "zdep", &
         description = "Set laser deposition type. options: zdep", &
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
         val = "OFF", &
         description = "Fix Neumann (zero-flux) boundary in Z", &
         data_type = DATA_BOOL))

    call table%set("ttm_boundary_heat_flux", control_parameter( &
         key = "ttm_boundary_heat_flux", &
         name = "TTM Heat flux", &
         val = "96", &
         units = "%", &
         internal_units = "", &
         description = "Set boundary heat flux in Robin boundaries for TTM", &
         data_type = DATA_BOOL)

    call table%set("ttm_time_offset", control_parameter( &
         key = "ttm_time_offset", &
         name = "TTM time offset", &
         val = "0.0", &
         units = "ps", &
         internal_units = "ps", &
         description = "Set initial electronic-ionic coupling via thermostat for TTM", &
         data_type = DATA_FLOAT))

    call table%set("ttm_oneway", control_parameter( &
         key = "ttm_oneway", &
         name = "TTM oneway", &
         val = "OFF", &
         description = "Enable one-way electron-phonon coupling when electronic temperature is greater than ionic temperature", &
         data_type = DATA_BOOL))

    call table%set("ttm_statis_frequency", control_parameter( &
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
         val = "Full", &
         description = "Apply inhomogeneous Langevin thermostat to selected directions in TTM, options: full, zdir, off", &
         data_type = DATA_OPTION))

    call table%set("ttm_redistribute", control_parameter( &
         key = "ttm_redistribute", &
         name = "", &
         val = "", &
         description = "Redistribute electronic energy in newly-deactivated temperature cells to nearest active neighbours", &
         data_type = DATA_BOOL))

    call table%set("msd_start", control_parameter( &
         key = "mst_start", &
         name = "MSD Start calculating", &
         val = "0", &
         units = "steps", &
         internal_units = "steps", &
         description = "Start timestep for dumping MSD configurations", &
         data_type = DATA_FLOAT))

    call table%set("msd_interval", control_parameter( &
         key = "mst_interval", &
         name = "MSD calculation interval", &
         val = "1", &
         units = "steps", &
         internal_units = "steps", &
         description = "Interval between dumping MSD configurations", &
         data_type = DATA_FLOAT))

    call table%set("analyse_all", control_parameter( &
         key = "analyse_all", &
         name = "Analyse all", &
         val = "OFF", &
         description = "Enable analysis for all bonds, angles, dihedrals and inversions", &
         data_type = DATA_BOOL))

    call table%set("analyse_bonds", control_parameter( &
         key = "analyse_bonds", &
         name = "Analyse bonds", &
         val = "OFF", &
         description = "Enable analysis for all bonds", &
         data_type = DATA_BOOL))

    call table%set("analyse_angles", control_parameter( &
         key = "analyse_angles", &
         name = "Analyse angles", &
         val = "OFF", &
         description = "Enable analysis for all angles", &
         data_type = DATA_BOOL))

    call table%set("analyse_dihedrals", control_parameter( &
         key = "analyse_dihedrals", &
         name = "Analyse dihedrals", &
         val = "OFF", &
         description = "Enable analysis for all dihedrals", &
         data_type = DATA_BOOL))

    call table%set("analyse_inversions", control_parameter( &
         key = "analyse_inversions", &
         name = "Analyse inversions", &
         val = "OFF", &
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

    call table%set("analyse_num_bins", control_parameter( &
         key = "analyse_num_bins", &
         name = "Analysis number of bins", &
         val = "0", &
         description = "Set number of bins to be used in bonding analysis", &
         data_type = DATA_INT))

    call table%set("analyse_max_dist", control_parameter( &
         key = "analyse_max_dist", &
         name = "Max analyse distance", &
         val = "2.0", &
         units = "ang", &
         internal_units = "internal_l", &
         description = "Set cutoff for bonds analysis", &
         data_type = DATA_FLOAT))

    call table%set("rdf_frequency", control_parameter( &
         key = "rdf_frequency", &
         name = "RDF Sampling Frequency", &
         val = "0", &
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

    call table%set("zden_frequency", control_parameter( &
         key = "zden_frequency", &
         name = "ZDen Sampling Frequency", &
         val = "0", &
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

    call table%set("vaf_frequency", control_parameter( &
         key = "vaf_frequency", &
         name = "VAF Sampling Frequency", &
         val = "0", &
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

    call table%set("print_frequency", control_parameter( &
         key = "print_frequency", &
         name = "Results Print Frequency", &
         val = "0", &
         units = "steps", &
         internal_units = "steps", &
         description = "Set frequency of results sampling", &
         data_type = DATA_FLOAT))

    call table%set("stack_size", control_parameter( &
         key = "stack_size", &
         name = "Rolling average stack size", &
         val = "0", &
         description = "Set rolling average stack to n timesteps", &
         data_type = DATA_INT))

    call table%set("stats_frequency", control_parameter( &
         key = "stats_frequency", &
         name = "Stats Print Frequency", &
         val = "0", &
         units = "steps", &
         internal_units = "steps", &
         description = "Set frequency of stats sampling", &
         data_type = DATA_FLOAT))

  end Subroutine initialise_control

  Subroutine parse_file(ifile, params, comm)

    Type( parameters_hash_table ), intent(   Out ) :: params
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
       Call lower_case(key)
       Call params%get(key, param)
       Call read_control_param(input, param, file, comm)
       call table%set(key, param)
    end do

  end Subroutine parse_file

  Subroutine read_control_param(input, param, ifile, comm)
    Type( control_parameter ), Intent ( InOut ) :: param
    Integer, Intent( In    ) :: ifile
    Type( comms_type ), Intent ( InOut ) :: comm
    Character(Len=*), Intent( InOut ) :: input
    Character(Len=STR_LEN) :: val, unit, junk
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
    case (DATA_VECTOR3. DATA_VECTOR6)
       call get_word(input, val)
       if (val(1:1) /= '[') call error(0, 'Expected vector in key '//trim(param%key)//' received '//trim(input))
       test_int = 1
       ! Skip [ on read
       if (trim(val(2:)) /= '') then
          test_real = word_2_real(val(2:))
          param%val = trim(val(2:))
       end if

       do
          call get_word(input, val)
          ! Handle multiline
          do while(trim(input) == '')
             call get_line(line_read, ifile, input, comm)
             if (.not. line_read) call error(0, 'End of file while parsing '//trim(param%key))
          end do

          i = index(val, ']')
          if (i > 1) then ! Val]
             test_real = word_2_real(val(1:i-1))
             param%val = trim(param%val) // " " // val(1:i-1)
             exit
          else if (i == 1) then ! Independent ]
             exit
          else
             test_real = word_2_real(val)
             param%val = trim(param%val) // " " // val
          end if

       end do

       ! Handle multiline
       do while(trim(input) == '')
          call get_line(line_read, ifile, input, comm)
          if (.not. line_read) call error(0, 'End of file while parsing '//trim(param%key))
       end do
       call get_word(input, unit)
       param%units = unit

    case (DATA_OPTION, DATA_BOOL)
       call lower_case(val)
       param%val = val
    case (DATA_STRING)
       param%val = val
    case Default
       call error(0, 'Unknown data type while parsing '//trim(param%key))
    end select

    if (trim(input) /= '') call error(0, "Unexpected junk "//trim(input)//" after key "//trim(param%key))

  End Subroutine read_control_param

  Function parse_vector3(input)
    Character(Len=*), Value :: input
    Character(Len=STR_LEN) :: parse
    Real(kind=wp), dimension(3) :: parse_vector3

    Integer :: i

    do i = 1, 3
       call get_word(input, parse)
       parse_vector3(i) = word_2_real(parse)
    end do

  End Function parse_vector3

  Function parse_int(input)
    Character(Len=*), Intent( In    ) :: input
    Integer :: parse_int

    read(input, '(i0)') parse_int

  End Function parse_int

  Function parse_bool(input)
    Character(Len=*), Intent( In    ) :: input
    Logical :: parse_bool

    Select Case(trim(input))
    Case('on', 'y')
       parse_bool = .true.
    Case Default
       parse_bool = .false.
    end Select

  End Function parse_bool

End Module new_control
