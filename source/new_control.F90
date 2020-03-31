Module new_control

  Use comms, only: comms_type
  Use hash, only: parameters_hash_table, control_parameter, MAX_LEN
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

  Subroutine read_new_control(file, params)
    Type( parameters_hash_table ), intent(   Out ) :: params
    Type( file_type ), Intent( InOut ) :: control_file

    Call parse_file(control_file%unit_no, params)

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

       End If

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

  Function parse_int(input)
    Character(Len=*), Intent( In    ) :: input
    Integer :: parse_int

    read(input, '(i0)') parse_int

  End Function parse_int

  Function parse_bool(input)
    Character(Len=*), Intent( In    ) :: input
    Logical :: parse_bool

    Call lower_case(input)

    Select Case(trim(input))
    Case('on', 'y')
       parse_bool = .true.
    Case Default
       parse_bool = .false.
    end Select

  End Function parse_bool

  Subroutine parse_file(file, params)

    Type( parameters_hash_table ), intent(   Out ) :: params
    Type( control_parameter ) :: param
    Integer, Intent( In    ) :: file
    Character(Len=MAX_LEN) :: input, key, val, units
    Logical :: line_read

    call table%init(256)

    do
       key = ''; val = ''; units = ''
       Call get_line(line_read,file,input,comm)
       if (.not. line_read) exit
       call get_word(input, key)
       Call lower_case(key)
       ! Special case for title -- Read as one val
       if (trim(key) == 'title') then
          val = adjustl(input)
       else
          call get_word(input, val)
          if (trim(input) /= '') call get_word(input, units)
       end if
       param = control_parameter(key, val, units)
       call table%set(key, param)
    end do

    call table%keyvals

  end Subroutine parse_file



End Module new_control
