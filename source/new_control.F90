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
  Use units, only : convert_units, set_timestep
  Implicit None

  Integer, Parameter, Public :: DATA_INT=1, DATA_FLOAT=2, DATA_STRING=3, DATA_BOOL=4, DATA_OPTION=5, DATA_VECTOR3=6, DATA_VECTOR6=7

  Type, Private, Extends(hash_table) :: parameters_hash_table
   contains
     Generic, Public  :: get => get_int, get_double, get_complex, get_param
     Procedure, Private :: get_param
     Generic, Public  :: retrieve => retrieve_option_or_string, retrieve_float, retrieve_vector3, retrieve_vector6, retrieve_int
     Procedure, Pass(table), Private :: retrieve_option_or_string, retrieve_float, retrieve_vector3, retrieve_vector6, retrieve_int
  end type parameters_hash_table


  Type, Public :: control_parameter
     !! Type containing breakdown of control parameter
     !> Internal key name
     Character(Len=30) :: key = ""
     !> Long name to be printed on high print_level
     Character(Len=STR_LEN) :: name = ""
     !> Current value -- Initialise to default
     Character(Len=STR_LEN) :: val = ""
     !> User specified units
     Character(Len=STR_LEN) :: units = ""
     !> Units to be converted to internally
     Character(Len=STR_LEN) :: internal_units = ""
     !> Information to be printed with help
     Character(Len=STR_LEN) :: description = ""
     !> Control parameter data type (int, float, vector3, vector6, string, bool, option)
     Character(Len=10) :: data_type = 0
     !> Is value set
     Logical :: set = .false.
  End Type control_parameter

contains

  Subroutine read_new_control(file, params, comm)
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
    table%can_overwrite = .false.

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
    Character(Len=256) :: curr_option
    Integer :: io_read, io_write
    Integer :: itmp
    Logical :: ltmp

    call params%retrieve('io_read_method', curr_option)

    Select Case(Trim(curr_option))
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
       Call info('io_read_method '//curr_option//' unrecognised option.',.true.)
       Call error(3)

    End Select

    Call io_set_parameters(io, user_method_read = io_read)

    Select Case (io_read)
    Case (IO_READ_MPIIO, IO_READ_DIRECT, IO_READ_NETCDF)
       ! Need to calculate number of readers
       call params%retrieve('io_read_readers', itmp)

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
       tmp = Min( Real(comm%mxnode,wp), 2.0_wp*Real(comm%mxnode,wp)**0.5_wp )
       itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
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
       call param%retrieve('io_read_error_check', ltmp)

       If (.not. ltmp) Then
          Call info('I/O parallel read error checking off',.true.)
       Else
          Call info('I/O parallel read error checking on',.true.)
       End If
       Call io_set_parameters(io, user_error_check = l_tmp )
    End If

    ! Get write settings

    call params%retrieve('io_write_method', curr_option)
    call params%retrieve('io_write_sorted', ltmp)

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

       call current_params%retrieve('io_read_netcdf_format', curr_option)

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
       call params%retrieve('io_write_writers', itmp)

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
       tmp = Min( Real(comm%mxnode,wp), 8.0_wp*Real(comm%mxnode,wp)**0.5_wp )
       itmp = 2**Int(Nearest( Log(tmp)/Log(2.0_wp) , +1.0_wp ))
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
       Call io_set_parameters(io, user_error_check = l_tmp )
    End If

    call params%retrieve('io_file_output', curr_val)
    if (trim(curr_val) /= '') Call info('OUTPUT file is '//files(FILE_OUTPUT)%filename,.true.)
    call params%retrieve('io_file_config', curr_val)
    if (trim(curr_val) /= '') Call info('CONFIG file is '//files(FILE_CONFIG)%filename,.true.)
    call params%retrieve('io_file_field', curr_val)
    if (trim(curr_val) /= '') Call info('FIELD file is '//files(FILE_FIELD)%filename,.true.)
    call params%retrieve('io_file_statis', curr_val)
    if (trim(curr_val) /= '') Call info('STATIS file is '//files(FILE_STATS)%filename,.true.)
    call params%retrieve('io_file_history', curr_val)
    if (trim(curr_val) /= '') Call info('HISTORY file is '//files(FILE_HISTORY)%filename,.true.)
    call params%retrieve('io_file_historf', curr_val)
    if (trim(curr_val) /= '') Call info('HISTORF file is '//files(FILE_HISTORF)%filename,.true.)
    call params%retrieve('io_file_revive', curr_val)
    if (trim(curr_val) /= '') Call info('REVIVE file is '//files(FILE_REVIVE)%filename,.true.)
    call params%retrieve('io_file_revcon', curr_val)
    if (trim(curr_val) /= '') Call info('REVCON file is '//files(FILE_REVCON)%filename,.true.)
    call params%retrieve('io_file_revold', curr_val)
    if (trim(curr_val) /= '') Call info('REVOLD file is '//files(FILE_REVOLD)%filename,.true.)

  End Subroutine setup_file_io

  Subroutine scan_new_control_pre(params, image_convention, density_variance)
    !!-----------------------------------------------------------------------
    !!
    !! Scan control for use in set_bounds
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), intent( In    ) :: params
    Integer, Intent( InOut ) :: image_convention
    Real(kind=wp), Intent(   Out) :: permitted_density_variation
    Logical :: ltmp

    call params%retrieve('slab', ltmp)
    call params%retrieve('density_variance', density_variance)

    if (ltmp) image_convention = 6

  end Subroutine scan_new_control_pre

  Subroutine scan_new_control_output(params, files)
    Type( parameters_hash_table ), intent( In    ) :: table
    Type( file_type ), Intent( InOut ) :: files(:)

    Call params%retrieve('io_file_output', files(FILE_OUTPUT)%filename)
    Call params%retrieve('io_file_config', files(FILE_CONFIG)%filename)
    Call params%retrieve('io_file_field',  files(FILE_FIELD)%filename)
    Call params%retrieve('io_file_stats',  files(FILE_STATS)%filename)
    Call params%retrieve('io_file_history',files(FILE_HISTORY)%filename)
    Call params%retrieve('io_file_historf',files(FILE_HISTORF)%filename)
    Call params%retrieve('io_file_revive', files(FILE_REVIVE)%filename)
    Call params%retrieve('io_file_revcon', files(FILE_REVCON)%filename)
    Call params%retrieve('io_file_revold', files(FILE_REVOLD)%filename)

  end Subroutine scan_new_control_output

  Subroutine scan_new_control(rcter,max_rigid,imcon,imc_n,cell,xhi,yhi,zhi,mxgana, &
    l_n_r,lzdn,l_ind,nstfce,ttm,cshell,stats,thermo,green,devel,msd_data,met, &
    pois,bond,angle,dihedral,inversion,zdensity,neigh,vdws,tersoffs,rdf,mpoles, &
    electro,ewld,kim_data,files,flow,comm)

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
    Type( development_type ), Intent( InOut ) :: devel
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
    Type( file_type ), Intent( InOut ) :: files(:)
    Type( flow_type ), Intent( InOut ) :: flow
    Type( comms_type ), Intent( InOut ) :: comm

    Character(Len=50) :: option
    Logical :: ltmp
    Logical :: la_ana,la_bnd,la_ang,la_dih,la_inv, &
         lrcut,lrvdw,lrmet,lelec,lvdw,lmet,l_n_m,lter,l_exp

    call params%retrieve('cutoff', neigh%cutoff)
    lrcut = neigh%cutoff > zero_plus

    call param%retrieve('padding', neigh%padding)
    if (neigh%padding < 0.0_wp) neigh%padding = 0.0_wp
    flow%reset_padding = .true.

    call param%retrieve('vdw_cutoff', vdws%cutoff)
    lrvdw = vdws%cutoff > zero_plus

    call param%retrieve('rdf_binsize', rdf%rbin)
    rdf%rbin = abs(rdf%rbin)

    call param%retrieve('zden_binsize', zdensity%bin_width)
    zdensity%bin_width = abs(zdensity%bin_width)

    call param%retrieve('ensemble', option, required = .true.)

    select case(trim(option))
    case ('nve', 'pmf')
       thermo%ensemble = NVE
       Call info('Ensemble : NVE (Microcanonical)',.true.)

    case ('nvt')

       call param%retrieve('ensemble_method', option, required = .true.)

       select case(trim(option))
       case ('evans')
          thermo%ensemble = ENS_NVT_EVANS

          Call info('Ensemble : NVT Evans (Isokinetic)',.true.)
          Call info('Gaussian temperature constraints in use',.true.)

       case ('langevin')

          thermo%ensemble = ENS_NVT_LANGEVIN

          Call param%retrieve('ensemble_thermostat_friction', thermo%chi)

          Call info('Ensemble : NVT Langevin (Stochastic Dynamics)',.true.)
          Write(message,'(a,1p,e12.4)') 'thermostat friction (ps^-1)', thermo%chi
          Call info(message,.true.)

       case ('andersen')

          thermo%ensemble = ENS_NVT_ANDERSON

          Call param%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          Call param%retrieve('ensemble_thermostat_softness', theremo%soft)

          Write(messages(1),'(a)') 'Ensemble : NVT Andersen'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'softness (dimensionless)',thermo%soft
          if (thermo%soft > 1) call error(0, 'Andersen softness greater than 1 '//message)

          Call info(messages,3,.true.)

       case ('berendsen')

          thermo%ensemble = ENS_NVT_BERENDSEN

          Call param%retrieve('ensemble_thermostat_coupling', thermo%tau_t)

          Call info('Ensemble : NVT Berendsen',.true.)
          Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Call info(message,.true.)

       Case ('hoover', 'nose', 'nose-hoover')

          thermo%ensemble = ENS_NVT_NOSE_HOOVER

          Call param%retrieve('ensemble_thermostat_coupling', thermo%tau_t)

          Call info('Ensemble : NVT Nose-Hoover',.true.)
          Write(message,'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Call info(message,.true.)

       Case ('gentle', 'gst')

          thermo%ensemble = ENS_NVT_GENTLE

          Call param%retrieve('ensemble_thermostat_coupling', thermo%tau_t)
          Call param%retrieve('ensemble_thermostat_friction', thermo%gama)

          Write(messages(1),'(a)') 'Ensemble : NVT gentle stochastic thermostat'
          Write(messages(2),'(a,1p,e12.4)') 'thermostat relaxation time (ps) ',thermo%tau_t
          Write(messages(3),'(a,1p,e12.4)') 'friction on thermostat  (ps^-1) ',thermo%gama
          Call info(messages,3,.true.)

       Case ('ttm')


          ! I hate TTM parameters. I refuse (maybe later)

       Case ('dpd')
          thermo%ensemble = ENS_NVE

          Call info('Ensemble : NVT dpd (Dissipative Particle Dynamics)',.true.)

          call params%retrieve('ensemble_dpd_order', option)
          select case (trim(option))
          case ('first', '1')
             thermo%key_dpd = DPD_FIRST_ORDER
             Call info("Ensemble type : Shardlow's first order splitting (S1)",.true.)
          case ('second', '2')
             thermo%key_dpd = DPD_SECOND_ORDER
             Call info("Ensemble type : Shardlow's first order splitting (S2)",.true.)
          case default
             call error(0, 'Unrecognised option for ensemble_dpd_order '//trim(option))
          end select


          call params%retrieve('ensemble_dpd_drag', thermo%gamdpd(0))

          If (thermo%gamdpd(0) > zero_plus) Then
             Write(message,'(a,1p,e12.4)') 'drag coefficient (Dalton/ps) ', thermo%gamdpd(0)
             Call info(message,.true.)
          End If

       case default
          call error(0, 'Unrecognised option for nvt ensemble_method '//trim(option))
       end select

    case ('npt')

       thermo%variable_cell = .true.

       call params%retrieve('ensemble_method', option)
       select case (trim(option))
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
          call error(0, 'Unrecognised option for npt ensemble_method '//trim(option))
       end select

    case ('nst')

       thermo%variable_cell = .true.
       thermo%anisotropic_pressure = .true.

       call params%retrieve('ensemble_method', option)
       select case (trim(option))
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
          call error(0, 'Unrecognised option for nst ensemble_method '//trim(option))
       end select

       ! Semi isotropic ensembles

       call params%retrieve('ensemble_semi_isotropic', option)
       select case (trim(option))
       case ('OFF')
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

          call params%retrieve('ensemble_semi_orthorhombie', ltmp)
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
          call error(0, 'Unrecognised option for ensemble '//trim(option))
       end select
       If (Any(thermo%iso == [CONSTRAINT_SURFACE_AREA,CONSTRAINT_SURFACE_TENSION])) Then
          Call warning('semi-isotropic ensembles are only correct for infinite' &
               //'interfaces placed perpendicularly to the z axis',.true.)
       End If

    case default
       call error(0, 'Unrecognised option for ensemble '//trim(option))
    end select


  end Subroutine scan_new_control



  Subroutine initialise_control(table)
    !!-----------------------------------------------------------------------
    !!
    !! Initialise all control parameters to their default
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------
    Type( parameters_hash_table ), intent(   Out ) :: table

    call table%init(256)
    table%can_overwrite = .true.

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

    call table%set("io_read_error_check", control_parameter( &
         key = "io_read_error_check", &
         name = "I/O error check on read", &
         vall = "OFF", &
         description = "Enable extended error checking on read", &
         data_type = DATA_BOOL))

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
         vall = "OFF", &
         description = "Enable extended error checking on write", &
         data_type = DATA_BOOL))

    call table%set("io_write_netcdf_format", control_parameter( &
         key = "io_write_netcdf_format", &
         name = "Netcdf Format", &
         val = "64bit", &
         description = "Set netcdf write format, options: 'amber', '32bit', '32-bit', '64-bit', '64bit'", &
         data_type = DATA_OPTION))

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
         description = "Set ensemble constraints, options: NVE, PMF, NVT, NPT, NST, NPnAT, NPng³", &
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
         val = "OFF", &
         description = "Enable semi-isotropic barostat constraints, options: area, tension, orthorhombic", &
         data_type = DATA_BOOL))

    call table%set("ensemble_semi_orthorhombic", control_parameter( &
         key = "ensemble_semi_orthorhombic", &
         name = "Ensemble semi-orthorhombic constraint", &
         val = "OFF", &
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
    !!-----------------------------------------------------------------------
    !!
    !! Read a control file filling the params table
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2020
    !!-----------------------------------------------------------------------

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
       if (param%set) call error(0, 'Param '//key//' already set')
       Call read_control_param(input, param, file, comm)
       param%set = .true.
       call table%set(key, param)
    end do

  end Subroutine parse_file

  Subroutine read_control_param(input, param, ifile, comm)
    Type( control_parameter ), Intent ( InOut ) :: param
    Integer, Intent( In    ) :: ifile
    Type( comms_type ), Intent ( InOut ) :: comm
    Character(Len=*), Intent( InOut ) :: input
    Character(Len=STR_LEN) :: tmp, val, unit, junk
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
          if (trim(input(i+1:) /= '')) call error(0, 'Unexpected junk '//trim(input)//' after key '//trim(param%key))
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

       do while (trim(tmp) /= '')
          call get_word(tmp, val)
          ! Check all valid reals
          test_real = word_2_real(val)
          test%val = trim(test%val)//' '//val
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
       call lower_case(val)
       param%val = val
    case (DATA_STRING)
       param%val = val
    case Default
       call error(0, 'Unknown data type while parsing '//trim(param%key))
    end select

    if (trim(input) /= '') call error(0, "Unexpected junk "//trim(input)//" after key "//trim(param%key))

  End Subroutine read_control_param

  Subroutine retrieve_option_or_string(table, key, output, required)
    Type( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Character(Len=STR_LEN), Intent( Out    ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    output = param%val

  end Subroutine retrieve_option_or_string

  Subroutine retrieve_float(table, key, output, required)
    Type( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Character(Len=STR_LEN), Intent( Out    ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    val = param%val
    call get_word(val, parse)
    output = word_2_real(parse)
    output = convert_units(output, param%units, param%internal_units)

  End Subroutine retrieve_float

  Subroutine retrieve_vector3(table, key, output, required)
    Type( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Real(kind=wp), dimension(3), Intent( Out    ) :: output

    Integer :: i

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    val = param%val

    do i = 1, 3
       call get_word(val, parse)
       if (trim(parse) == "") exit
       output(i) = word_2_real(parse)
       output(i) = convert_units(output(i), param%units, param%internal_units)
    end do
    if (i /= 3) call error(0, "Bad length vector 3")

  End Subroutine retrieve_vector3

  Subroutine retrieve_vector6(table, key, output, required)
    Type( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse
    Type( control_parameter ) :: param
    Real(kind=wp), dimension(9) :: tmp
    Logical, Intent( In    ), Optional :: required
    Real(kind=wp), dimension(6), Intent( Out    ) :: output

    Integer :: i

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    val = param%val

    do i = 1, 9
       call get_word(val, parse)
       if (trim(parse) == "") exit
       tmp(i) = word_2_real(parse)
       tmp(i) = convert_units(tmp(i), param%units, param%internal_units)
    end do

    select case(i)
    case(7)
       output = tmp(1:6)
    case(9)
       output = [tmp(1), tmp(5), tmp(9), tmp(2), tmp(3), tmp(4)]
    case default
       call error(0, "Bad length vector 6")
    end select


  End Subroutine retrieve_vector6

  Subroutine retrieve_int(table, key, output, required)
    Type( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Integer, Intent(   Out ) :: output

    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if
    read(param%val, '(i0)') output

  End Subroutine retrieve_int

  Subroutine retrieve_bool(table, key, output, required)
    Type( parameters_hash_table ) :: table
    Character(Len=*), Intent( In    ) :: key
    Character(Len=STR_LEN) :: parse
    Type( control_parameter ) :: param
    Logical, Intent( In    ), Optional :: required
    Logical, Intent(   Out ) :: output



    call table%get(key, param)
    if (present(required)) then
       if (required .and. .not. param%set) call error(0, 'Necessary parameter '//trim(key)//' not set')
    end if

    Select Case(trim(param%val))
    Case('on', 'y')
       output = .true.
    Case Default
       output = .false.
    end Select

  End Subroutine retrieve_bool

  Subroutine get_param( table, key, val, default )
    Class( parameters_hash_table ), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type(control_parameter), Intent(   Out ) :: val
    Type(control_parameter), Intent( In    ), Optional :: default
    Class( * ), Allocatable :: stuff

    stuff = table%get_cont(key, default)

    Select Type( stuff )
    Type is ( control_param )
       val = stuff
    Class Default
       Call error('Trying to get control_param from a not control_param')
    End Select

  End Subroutine get_param

End Module new_control
