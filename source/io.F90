Module io

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module implementing the sorted parallel I/O strategy.  See
  ! http://www.hpcx.ac.uk/research/hpc/technical_reports/HPCxTR0806.pdf
  ! for more details which describes the ideas behind the method.
  !
  ! Note this module does NOT deal with any headers or extra "meta-data"
  ! that the dump may require - this deals only with writing data about
  ! the atoms themselves.
  !
  ! As a general point no return codes from MPI routines are checked because
  ! the default behaviour is not to return but to abort on error ...
  ! And nobody bothers to change that ( famous last words ).
  !
  ! copyright - daresbury laboratory
  ! author    - i.j.bush april 2011
  ! amended   - m.a.seaton november 2011
  ! amended   - a.m.elena september 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp, li
  Use comms, Only : comms_type,  wp_mpi
  Use errors_warnings, Only : error
  Use timer, Only : start_timer_path, stop_timer_path
  Use particle, Only: corePart
#ifdef SERIAL
  Use mpi_api
#else
  Use mpi
#endif
  Use netcdf_wrap, Only : netcdf_param,netcdf_desc,netcdf_set_real_precision, &
    netcdf_get_real_precision,netcdf_set_def,netcdf_compiled, &
    netcdf_get_file_real_precision,netcdf_get_att, &
    netcdf_get_var,netcdf_put_var,netcdf_close,netcdf_open, &
    netcdf_get_def,netcdf_get_dim,netcdf_create

  Implicit None

  Integer, Parameter, Public     :: recsz = 73 ! default record size
  ! MPI-I/O representation
  Character( Len = 6 ), Parameter :: datarep = 'native'
  ! Public interface

  ! Values to tell the module quite which file we are writing
  Integer, Parameter, Public :: IO_RESTART  = 1
  Integer, Parameter, Public :: IO_HISTORY  = 2
  Integer, Parameter, Public :: IO_HISTORD  = 21
  Integer, Parameter, Public :: IO_MSDTMP   = 3

  ! Values to indicate what I/O method is to be used
  Integer, Parameter, Public :: IO_NO_METHOD             = -1

  Integer, Parameter, Public :: IO_WRITE_UNSORTED_MPIIO  = 0
  Integer, Parameter, Public :: IO_WRITE_UNSORTED_DIRECT = 1
  Integer, Parameter, Public :: IO_WRITE_UNSORTED_MASTER = 2
  Integer, Parameter, Public :: IO_WRITE_SORTED_MPIIO    = 10
  Integer, Parameter, Public :: IO_WRITE_SORTED_DIRECT   = 11
  Integer, Parameter, Public :: IO_WRITE_SORTED_MASTER   = 12
  Integer, Parameter, Public :: IO_WRITE_SORTED_NETCDF   = 13

  Integer, Parameter, Public :: IO_READ_MPIIO   = 0
  Integer, Parameter, Public :: IO_READ_DIRECT  = 1
  Integer, Parameter, Public :: IO_READ_MASTER  = 2
  Integer, Parameter, Public :: IO_READ_NETCDF  = 3

  !Values to set which particle subset we're using
  Integer, Parameter, Public :: IO_SUBSET_POSITIONS = 1
  Integer, Parameter, Public :: IO_SUBSET_FORCES    = 2
  Integer, Parameter, Public :: IO_SUBSET_CHARGES   = 4

  ! Error codes
  Integer, Parameter, Public :: IO_BASE_COMM_NOT_SET     = 1
  Integer, Parameter, Public :: IO_ALLOCATION_ERROR      = 2
  Integer, Parameter, Public :: IO_DEALLOCATION_ERROR    = 3
  Integer, Parameter, Public :: IO_UNKNOWN_WRITE_OPTION  = 4
  Integer, Parameter, Public :: IO_UNKNOWN_WRITE_LEVEL   = 5

  Public :: io_init
  Public :: io_finalize
  Public :: io_open
  Public :: io_close
  Public :: io_delete
  Public :: io_set_parameters
  Public :: io_get_parameters
  Public :: io_write_sorted_file
  Public :: io_write_record
  Public :: io_write_batch
  Public :: io_read_batch
  Public :: io_nc_compiled
  Public :: io_nc_create
  Public :: io_nc_get_dim
  Public :: io_nc_set_def
  Public :: io_nc_put_var
  Public :: io_get_var
  Public :: io_nc_get_var
  Public :: io_nc_get_att
  Public :: io_nc_set_real_precision
  Public :: io_nc_get_real_precision
  Public :: io_nc_get_file_real_precision

  ! Private entities
  Private

  ! Default values
  Integer  , Parameter :: default_method_write      = IO_WRITE_SORTED_MPIIO
  Integer  , Parameter :: default_method_read       = IO_READ_MPIIO
  Integer  , Parameter :: default_batch_size_write  = 2000000
  Integer  , Parameter :: default_batch_size_read   = 2000000
  Integer  , Parameter :: default_n_io_procs_write  = 8
  Integer  , Parameter :: default_n_io_procs_read   = 2
  Integer  , Parameter :: default_buffer_size_write = 20000
  Integer  , Parameter :: default_buffer_size_read  = 20000
  Character, Parameter :: default_lf                = new_line('A')

  ! Range of FORTRAN unit numbers we shall consider when trying to find a unit number
  Integer, Parameter :: LOW_HANDLE  = 10
  Integer, Parameter :: HIGH_HANDLE = 999

  Integer, Parameter :: UNUSED       = 0
  Integer, Parameter :: FILE_MPI     = 1
  Integer, Parameter :: FILE_FORTRAN = 2
  Integer, Parameter :: FILE_NETCDF  = 3

  Integer, Parameter :: READ_ONLY  = 1
  Integer, Parameter :: WRITE_ONLY = 2

  Type file_data
    Character( Len = 32 )  :: name        = 'UNDEFINED'
    Integer                :: method      = UNUSED
    Integer                :: file_handle = -1
    Integer                :: action      = -1
    Type( netcdf_desc )    :: desc
  End Type file_data
  Type, Public :: io_type
    ! The record size
    Integer:: rec_size = recsz

    ! The record length for FORTRAN write
    Character( Len = 6 ):: forma    = '73a'
    ! Parameters for method.  Default to sensible values where possible
    Integer  :: base_comm          = MPI_COMM_NULL                   ! Top level MPI Communicator
    Integer  :: method_write       = default_method_write            ! How to perform the I/O
    Integer  :: method_read        = default_method_read             ! How to perform the I/O
    Integer  :: n_io_procs_write   = default_n_io_procs_write        ! Number of I/O processors
    Integer  :: n_io_procs_read    = default_n_io_procs_read         ! Number of I/O processors
    Integer  :: batch_size_write   = default_batch_size_write        ! Number of atoms in each batch
    Integer  :: batch_size_read    = default_batch_size_read         ! Number of atoms in each batch
    Integer  :: buffer_size_write  = default_buffer_size_write       ! Number of lines in each buffer
    Integer  :: buffer_size_read   = default_buffer_size_read        ! Number of lines in each buffer
    Logical  :: global_error_check = .false.                         ! Whether to check across all processors for errors.
    Character:: lf    =  default_lf                     ! The line feed character to use
    ! These indices depend on the write level, so can't parameterise them
    Integer :: Q_IND
    Integer :: W_IND
    Integer :: D_IND
    ! The handle for the MPI derived type used in the I/O
    Integer :: rec_type
    Type( file_data ), Dimension( LOW_HANDLE:HIGH_HANDLE ):: known_files
  End Type io_type
  ! Parameters to make indexing of combined arrays more legible
  Integer, Parameter :: RX_IND = 1
  Integer, Parameter :: RY_IND = 2
  Integer, Parameter :: RZ_IND = 3
  Integer, Parameter :: VX_IND = 4
  Integer, Parameter :: VY_IND = 5
  Integer, Parameter :: VZ_IND = 6
  Integer, Parameter :: FX_IND = 7
  Integer, Parameter :: FY_IND = 8
  Integer, Parameter :: FZ_IND = 9

  ! "Extra" data over and above the basics that needs to be written for each file.
  Integer, Parameter :: N_RESTART_DATA = 0
  Integer, Parameter :: N_HISTORY_DATA = 3
  Integer, Parameter :: N_HISTORD_DATA = 1
  Integer, Parameter :: N_MSDTMP_DATA  = 2






  Interface io_nc_put_var
    Module Procedure io_nc_put_var_rwp_0d
    Module Procedure io_nc_put_var_rwp_1d
    Module Procedure io_nc_put_var_rwp_2d
    Module Procedure io_nc_put_var_int_0d
    Module Procedure io_nc_put_var_int_1d
    Module Procedure io_nc_put_var_chr_2d
  End Interface io_nc_put_var

  Interface io_nc_get_var
    Module Procedure io_nc_get_var_rwp_0d
    Module Procedure io_nc_get_var_rwp_1d
    Module Procedure io_nc_get_var_rwp_2d
    Module Procedure io_nc_get_var_int_0d
    Module Procedure io_nc_get_var_int_1d
    Module Procedure io_nc_get_var_chr_1d
    Module Procedure io_nc_get_var_chr_2d
  End Interface io_nc_get_var

  Interface io_nc_get_att
    Module Procedure io_nc_get_att_int
    Module Procedure io_nc_get_att_chr
  End Interface io_nc_get_att

  Interface io_write_sorted_file
    Module Procedure io_write_sorted_file_parts
    Module Procedure io_write_sorted_file_parts_subset
    Module Procedure io_write_sorted_file_arrays
  End Interface io_write_sorted_file

Contains

  Subroutine io_init( io,this_rec_size )

    ! Initialise the sorted I/O method

    Type(io_type), Intent( InOut ) :: io
    Integer, Intent( In    ) :: this_rec_size
    Integer :: ierr

    io%rec_size = this_rec_size

    ! Create the derived type for I/O if MPI-I/O is being used anywhere
    If ( ( io%method_write == IO_WRITE_UNSORTED_MPIIO .or. io%method_write == IO_WRITE_SORTED_MPIIO ) .or. &
      io%method_read  == IO_READ_MPIIO ) Then
      Call MPI_TYPE_CONTIGUOUS( io%rec_size, MPI_CHARACTER, io%rec_type, ierr )
      Call MPI_TYPE_COMMIT( io%rec_type, ierr )
    End If

    ! If using direct access set up the format
    If ( ( io%method_write == IO_WRITE_UNSORTED_DIRECT .or. io%method_write == IO_WRITE_SORTED_DIRECT ) .or. &
      io%method_read  == IO_READ_DIRECT ) Then
      io%forma =' '
      Write( io%forma, "('(',i0,'a)')" ) io%rec_size
    End If

  End Subroutine io_init

  Subroutine io_finalize(io)
    Type(io_type), Intent( InOut ) :: io

    ! Finalize the sorted I/O method for MPI-I/O

    Integer :: ierr
    ! Free the MPI derived type if MPI-I/O is used anywhere
    If ( ( io%method_write == IO_WRITE_UNSORTED_MPIIO .or. io%method_write == IO_WRITE_SORTED_MPIIO ) .or. &
      io%method_read  == IO_READ_MPIIO ) Then
      Call MPI_TYPE_FREE( io%rec_type, ierr )
    End If

    ! Free the communicator - note after a comm_free the communicator
    ! is, by the standard, set to MPI_COMM_NULL
    Call MPI_COMM_FREE( io%base_comm, ierr )

  End Subroutine io_finalize

  Subroutine io_nc_create( param, comm, file_name, title, n )

    ! Create a netCDF file.  Note it will overwrite any existing file.

    Type( netcdf_param ), Intent( In    ) :: param
    Integer             , Intent( In    ) :: comm
    Character( Len = * ), Intent( In    ) :: file_name
    Character( Len = * ), Intent( In    ) :: title
    Integer             , Intent( In    ) :: n

    Type( netcdf_desc ) :: desc

    Call netcdf_create( Trim( file_name ), desc, comm, MPI_INFO_NULL )
    Call netcdf_set_def(title, n, param, desc )
    Call netcdf_close( desc )

  End Subroutine io_nc_create

  Subroutine io_open(io, method, comm, file_name, flags, file_handle )

    ! Open a file for access
    ! Arguments are:
    ! METHOD     : The I/O method, needed to get the correct intent:
    !              RD or WR as only one is allowed at a time/call here!
    ! COMM       : The communicator within which the I/O will occur
    ! FILE_NAME  : the file name
    ! FLAGS      : The flags with which the file shall be opened
    ! FILE_HANDLE: The returned file handle

    Type(io_type), Intent( InOut ) :: io
    Integer             , Intent( In    ) :: method
    Integer             , Intent( In    ) :: comm
    Character( Len = * ), Intent( In    ) :: file_name
    Integer             , Intent( In    ) :: flags
    Integer             , Intent(   Out ) :: file_handle

    Type( netcdf_desc ) :: desc

    Integer :: fh,ierr

    Logical :: do_read, do_write

    do_read  = method == IO_READ_MPIIO  .or. &
      method == IO_READ_DIRECT .or. &
      method == IO_READ_NETCDF
    do_write = method == IO_WRITE_UNSORTED_MPIIO  .or. method == IO_WRITE_SORTED_MPIIO  .or. &
      method == IO_WRITE_UNSORTED_DIRECT .or. method == IO_WRITE_SORTED_DIRECT .or. &
      method == IO_WRITE_SORTED_NETCDF

    Call get_file_handle(io, file_handle )
    If ( file_handle == LOW_HANDLE - 1 ) Then
      file_handle = -1
      Return
    End If

    If      ( ( method == IO_WRITE_UNSORTED_MPIIO .or. method == IO_WRITE_SORTED_MPIIO ) .or. &
      ( method == IO_READ_MPIIO ) ) Then

      ! using MPI-I/O
      Call MPI_FILE_OPEN( comm, file_name, flags, MPI_INFO_NULL, fh, ierr )

      ! And set up how it is to be accessed
      Call MPI_FILE_SET_VIEW( fh, Int(0,MPI_OFFSET_KIND), io%rec_type, io%rec_type, &
        datarep, MPI_INFO_NULL, ierr )

      If ( ierr == 0 ) Then
        io%known_files( file_handle )%method      = FILE_MPI
        io%known_files( file_handle )%file_handle = fh
        io%known_files( file_handle )%name        = file_name
        If ( do_read ) Then
          io%known_files( file_handle )%action = read_only
        Else If ( do_write ) Then
          io%known_files( file_handle )%action = write_only
        End If
      Else
        file_handle = -1
      End If

    Else If ( ( method == IO_WRITE_UNSORTED_DIRECT .or. method == IO_WRITE_SORTED_DIRECT ) .or. &
      ( method == IO_READ_DIRECT ) ) Then

      ! using parallel direct access
      If ( file_handle /= LOW_HANDLE - 1 ) Then

        If ( do_read ) Then
          Open(Unit=file_handle, File=file_name, Form='formatted', Access='direct', Action = 'read' , Recl=io%rec_size)
        Else If ( do_write ) Then
          Open(Unit=file_handle, File=file_name, Form='formatted', Access='direct', Action = 'write', Recl=io%rec_size)
        End If

        io%known_files( file_handle )%method      = FILE_FORTRAN
        io%known_files( file_handle )%file_handle = file_handle
        io%known_files( file_handle )%name        = file_name
        If ( do_read ) Then
          io%known_files( file_handle )%action = read_only
        Else If ( do_write ) Then
          io%known_files( file_handle )%action = write_only
        End If

      Else

        file_handle = -1

      End If

    Else If ( ( method == IO_WRITE_SORTED_NETCDF ) .or. &
      ( method == IO_READ_NETCDF ) ) Then

      If ( file_handle /= LOW_HANDLE - 1 ) Then

        Call netcdf_open( Trim( file_name ), desc, comm, MPI_INFO_NULL )
        Call netcdf_get_def( desc )

        io%known_files( file_handle )%method      = FILE_NETCDF
        io%known_files( file_handle )%file_handle = file_handle
        io%known_files( file_handle )%desc        = desc
        io%known_files( file_handle )%name        = file_name
        If ( do_read ) Then
          io%known_files( file_handle )%action = read_only
        Else If ( do_write ) Then
          io%known_files( file_handle )%action = write_only
        End If

      Else

        file_handle = -1

      End If

    Else

      file_handle = -1

    End If

  End Subroutine io_open

  Subroutine io_close(io, io_file_handle )

    ! Close a file
    ! Arguments are:
    ! FILE_HANDLE: The file handle

    Type(io_type), Intent( InOut ) :: io
    Integer, Intent( In    ) :: io_file_handle

    Integer :: file_handle
    Integer :: ierr

    If ( io%known_files( io_file_handle )%method /= UNUSED ) Then
      file_handle = io%known_files( io_file_handle )%file_handle

      ! using MPI-I/O
      If ( io%known_files( io_file_handle )%method == FILE_MPI ) Then
        Call MPI_FILE_CLOSE( file_handle, ierr )

        ! using parallel direct access
      Else If ( io%known_files( io_file_handle )%method == FILE_FORTRAN ) Then
        ! Leave in sync
        Call MPI_BARRIER( io%base_comm, ierr )
        Close(Unit=file_handle)

        ! using netCDF
      Else If ( io%known_files( io_file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_close( io%known_files( io_file_handle )%desc )
      End If

      io%known_files( io_file_handle )%name        = 'UNDEFINED'
      io%known_files( io_file_handle )%method      = UNUSED
      io%known_files( io_file_handle )%file_handle = -1
      io%known_files( io_file_handle )%action      = -1
    End If

  End Subroutine io_close

  Subroutine io_nc_get_dim(io, what, io_file_handle, val )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Integer             , Intent(   Out ) :: val

    Call netcdf_get_dim( what, io%known_files( io_file_handle )%desc, val )

  End Subroutine io_nc_get_dim

  Subroutine io_delete(io, file_name, comm )

    ! Delete a file - note that it is NOT an error to try to
    ! delete a non-existent file, it is silently ignored
    ! THIS ROUTINE BLOCKS ACROSS ALL THE COMM.
    ! Arguments are:
    ! FILE_NAME  : the file name

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: file_name
    Type(comms_type), Intent( InOut ) :: comm

    Integer :: file_handle, ierr

    Logical :: exist

    ! Only 1 need do it and avoid contention.
    If ( comm%idnode == 0 ) Then

      Inquire( File = file_name, Exist = exist )

      ! If it does delete it
      If ( exist ) Then

        Call get_file_handle(io, file_handle )

        Open(Unit=file_handle, File=file_name)
        Close(Unit=file_handle, Status='delete')

      End If

    End If

    ! Make sure that proc 0 has finished deleting the file before
    ! everybody else carriers on
    Call MPI_BARRIER( io%base_comm, ierr )

  End Subroutine io_delete

  Subroutine io_set_parameters(io, user_comm, user_method_write, user_method_read, user_n_io_procs_write, user_n_io_procs_read, &
      user_batch_size_write, user_batch_size_read, user_buffer_size_write, user_buffer_size_read,                           &
      user_line_feed, user_error_check )

    ! Set the options for the module
    ! Arguments are:
    ! USER_COMM             : The communicator that the data for this run spans
    ! USER_METHOD_WRITE     : Control how the writing is performed
    ! USER_METHOD_READ      : Control how the reading is performed
    ! USER_N_IO_PROCS_WRITE : The number of writing processors.  A non-positive number
    !                         means use the default number of processors.
    ! USER_N_IO_PROCS_READ  : The number of reading processors.  A non-positive number
    !                         means use the default number of processors.
    ! USER_BATCH_SIZE_WRITE : The maximum number of atoms that the writing processors
    !                         deal with at any one time.  A non-positive number means use the
    !                         default
    ! USER_BATCH_SIZE_READ  : The maximum number of atoms that the reading processors
    !                         deal with at any one time.  A non-positive number means use the
    !                         default
    ! USER_BUFFER_SIZE_WRITE: The maximum number of records that will be written to disk in
    !                         any 1 write transaction
    ! USER_BUFFER_SIZE_READ : The maximum number of records that will be read from disk in
    !                         any 1 read transaction
    ! USER_LINE_FEED        : The Character to use as a line feed character
    ! USER_ERROR_CHECK      : If true perform error checking across all the processors ( the default ).
    !                         If false perform error checking only locally.  This latter option
    !                         is less safe but has a much lower communication cost.

    Type(io_type), Intent( InOut ) :: io
    Integer  , Intent( In    ), Optional :: user_comm
    Integer  , Intent( In    ), Optional :: user_method_write
    Integer  , Intent( In    ), Optional :: user_method_read
    Integer  , Intent( In    ), Optional :: user_n_io_procs_write
    Integer  , Intent( In    ), Optional :: user_n_io_procs_read
    Integer  , Intent( In    ), Optional :: user_batch_size_write
    Integer  , Intent( In    ), Optional :: user_batch_size_read
    Integer  , Intent( In    ), Optional :: user_buffer_size_write
    Integer  , Intent( In    ), Optional :: user_buffer_size_read
    Character, Intent( In    ), Optional :: user_line_feed
    Logical  , Intent( In    ), Optional :: user_error_check

    Integer :: ierr
    ! Set the communicator
    If ( Present( user_comm ) ) Then
      ! Use a copy of the communicator to avoid potential message clashes with what
      ! is happening outside this module - e.g. if asynchronous messages are in
      ! flight while the I/O is occurring
      Call MPI_COMM_DUP( user_comm, io%base_comm, ierr )
    End If

    ! The writing method
    If ( Present( user_method_write ) ) Then
      io%method_write = user_method_write
    End If

    ! The reading method
    If ( Present( user_method_read ) ) Then
      io%method_read = user_method_read
    End If

    ! Set the number of writing processors
    If ( Present( user_n_io_procs_write ) ) Then
      If ( user_n_io_procs_write > 0 ) Then
        io%n_io_procs_write = user_n_io_procs_write
      Else
        io%n_io_procs_write = default_n_io_procs_write
      End If
    End If

    ! Set the number of reading processors
    If ( Present( user_n_io_procs_read ) ) Then
      If ( user_n_io_procs_read > 0 ) Then
        io%n_io_procs_read = user_n_io_procs_read
      Else
        io%n_io_procs_read = default_n_io_procs_read
      End If
    End If

    ! Set the batch size
    If ( Present( user_batch_size_write ) ) Then
      If ( user_batch_size_write > 0 ) Then
        io%batch_size_write = user_batch_size_write
      Else
        io%batch_size_write = default_batch_size_write
      End If
    End If

    If ( Present( user_batch_size_read ) ) Then
      If ( user_batch_size_read > 0 ) Then
        io%batch_size_read = user_batch_size_read
      Else
        io%batch_size_read = default_batch_size_read
      End If
    End If

    ! Set the buffer size
    If ( Present( user_buffer_size_write ) ) Then
      If ( user_buffer_size_write > 0 ) Then
        io%buffer_size_write = user_buffer_size_write
      Else
        io%buffer_size_write = default_buffer_size_write
      End If
    End If

    If ( Present( user_buffer_size_read ) ) Then
      If ( user_buffer_size_read > 0 ) Then
        io%buffer_size_read = user_buffer_size_read
      Else
        io%buffer_size_read = default_buffer_size_read
      End If
    End If

    ! Set the line feed
    If ( Present( user_line_feed ) ) Then
      io%lf = user_line_feed
    End If

    ! Set the error checking method
    If ( Present( user_error_check ) ) Then
      io%global_error_check = user_error_check
    End If

  End Subroutine io_set_parameters

  Subroutine io_get_parameters(io, user_comm, user_method_write, user_method_read, user_n_io_procs_write, user_n_io_procs_read, &
      user_batch_size_write, user_batch_size_read, user_buffer_size_write, user_buffer_size_read,                           &
      user_line_feed, user_error_check )

    ! Return the options currently in use by the module
    ! Arguments are:
    ! USER_COMM             : The communicator that the data for this run spans
    ! USER_METHOD_WRITE     : The method used for writing
    ! USER_METHOD_READ      : The method used for reading
    ! USER_N_IO_PROCS_WRITE : The number of writing processors.  A non-positive number
    !                         means use the default number of processors.
    ! USER_N_IO_PROCS_READ  : The number of reading processors.  A non-positive number
    !                         means use the default number of processors.
    ! USER_BATCH_SIZE_WRITE : The maximum number of atoms that the writing processors
    !                         deal with at any one time.  A non-positive number means use the
    !                         default
    ! USER_BATCH_SIZE_READ  : The maximum number of atoms that the reading processors
    !                         deal with at any one time.  A non-positive number means use the
    !                         default
    ! USER_LINE_FEED        : The character being used for a line feed
    ! USER_ERROR_CHECK      : If true perform error checking across all the processors ( the default ).
    !                         If false perform error checking only locally.  This latter option
    !                         is less safe but has a much lower communication cost.

    Type(io_type), Intent( InOut ) :: io
    Integer  , Intent(   Out ), Optional :: user_comm
    Integer  , Intent(   Out ), Optional :: user_method_write
    Integer  , Intent(   Out ), Optional :: user_method_read
    Integer  , Intent(   Out ), Optional :: user_n_io_procs_write
    Integer  , Intent(   Out ), Optional :: user_n_io_procs_read
    Integer  , Intent(   Out ), Optional :: user_batch_size_write
    Integer  , Intent(   Out ), Optional :: user_batch_size_read
    Integer  , Intent(   Out ), Optional :: user_buffer_size_write
    Integer  , Intent(   Out ), Optional :: user_buffer_size_read
    Character, Intent(   Out ), Optional :: user_line_feed
    Logical  , Intent(   Out ), Optional :: user_error_check

    ! Get the communicator
    If ( Present( user_comm ) ) Then
      user_comm = io%base_comm
    End If

    ! Get the writing method
    If ( Present( user_method_write ) ) Then
      user_method_write = io%method_write
    End If

    ! Get the reading method
    If ( Present( user_method_read ) ) Then
      user_method_read = io%method_read
    End If

    ! Get the number of I/O processors
    If ( Present( user_n_io_procs_write ) ) Then
      user_n_io_procs_write = io%n_io_procs_write
    End If

    ! Get the number of I/O processors
    If ( Present( user_n_io_procs_read ) ) Then
      user_n_io_procs_read = io%n_io_procs_read
    End If

    ! Get the batch size
    If ( Present( user_batch_size_write ) ) Then
      user_batch_size_write = io%batch_size_write
    End If

    If ( Present( user_batch_size_read ) ) Then
      user_batch_size_read = io%batch_size_read
    End If

    ! Get the buffer size
    If ( Present( user_buffer_size_write ) ) Then
      user_buffer_size_write = io%buffer_size_write
    End If

    If ( Present( user_buffer_size_read ) ) Then
      user_buffer_size_read = io%buffer_size_read
    End If

    ! Get the line feed
    If ( Present( user_line_feed ) ) Then
      user_line_feed = io%lf
    End If

    ! Get the error checking method
    If ( Present( user_error_check ) ) Then
      user_error_check = io%global_error_check
    End If

  End Subroutine io_get_parameters

  Subroutine io_write_sorted_file_parts(io, file_handle, write_level, write_options, first_record,   &
      n_atoms, global_indices, atom_name, weight, rsd, &
      parts, vx, vy, vz,                      &
      error )

    ! Write a single atomic configuration using MPI-I/O to the file corresponding to FILE_HANDLE.
    ! Arguments are:
    ! FILE_HANDLE   : The file handle
    ! WRITE_LEVEL   : How detailed the output will be, i.e. if velocities and forces are written
    ! WRITE_OPTIONS : Identifies quite which version of the closely related files is being written
    ! FIRST_RECORD  : Where the first output will be in the file
    ! N_ATOMS       : The number of atoms stored by this processor
    ! GLOBAL_INDICES: The global index number of the atom
    ! ATOM_NAME     : The atoms' names
    ! WEIGHT        : The atoms' weights
    ! CHARGE        : The atoms' charges
    ! RSD           : The atoms' RMS displacement (??)
    ! RX, RY, RZ    : The components of the atoms' positions
    ! VX, VY, VZ    : The components of the atoms' velocities
    ! FX, FY, FZ    : The components of the forces on the atoms
    ! ERROR         : Return code.  Non-zero indicates an error

    Type(io_type), Intent( InOut ) :: io
    Integer                             , Intent( In    ) :: file_handle
    Integer                             , Intent( In    ) :: write_level
    Integer                             , Intent( In    ) :: write_options
    Integer( Kind = MPI_OFFSET_KIND )   , Intent( In    ) :: first_record
    Integer                             , Intent( In    ) :: n_atoms
    Integer             , Dimension( : ), Intent( In    ) :: global_indices
    Character( Len = * ), Dimension( : ), Intent( In    ) :: atom_name
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: weight
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: rsd
    Type( corePart )    , Dimension( : ), Intent( In    ) :: parts
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vx
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vy
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vz
    Integer                             , Intent(   Out ) :: error

    Type( netcdf_desc ) :: desc

    Real( Kind = wp ), Dimension( :, : ), Allocatable :: local_data
    Real( Kind = wp ), Dimension( :, : ), Allocatable :: gathered_data
    Real( Kind = wp ), Dimension( :, : ), Allocatable :: sorted_data

    Integer, Dimension( : ), Allocatable :: global_index_rank
    Integer, Dimension( : ), Allocatable :: local_global_indices
    Integer, Dimension( : ), Allocatable :: gathered_global_indices
    Integer, Dimension( : ), Allocatable :: sorted_indices
    Integer, Dimension( : ), Allocatable :: n_reorg

    Integer :: io_comm, io_gather_comm
    Integer :: actual_io_procs
    Integer :: size_local
    Integer :: bottom_batch, top_batch
    Integer :: local_bottom, local_top
    Integer :: this_batch_size
    Integer :: n_gathered
    Integer :: tot_atoms
    Integer :: me_in_io, n_me
    Integer :: itmp
    Integer :: iter

    Logical :: do_io

    Character( Len = 1 ), Dimension( :, : ), Allocatable :: local_name
    Character( Len = 1 ), Dimension( :, : ), Allocatable :: gathered_name
    Character( Len = 1 ), Dimension( :, : ), Allocatable :: sorted_name
    Integer :: ierr
    ! Ever the optimist
    error = 0

    ! Check we have a communicator
    If ( .not. ok(io, io%base_comm /= MPI_COMM_NULL, io%base_comm ) ) Then
      error = IO_BASE_COMM_NOT_SET
      Return
    End If

    ! Create the I/O communicators, IO_COMM which contains only those procs
    ! which will perform I/O, and IO_GATHER_COMM, which contains at rank 0 the
    ! proc which will perform the I/O for all other procs in the give instantiation
    ! of IO_GATHER_COMM.  As a processor can not ( easily ) determine if it is
    ! a member of a given communicator or not ( if it's not in a communicator
    ! it has no data about it ) also return DO_IO to indicate whether this
    ! processor will perform I/O, i.e.  if it is a member of IO_COMM.
    Call split_io_comm( io%base_comm, io%n_io_procs_write, io_comm, io_gather_comm, do_io )
    If ( do_io ) Then
      Call MPI_COMM_RANK( io_comm, me_in_io, ierr )
      Call MPI_COMM_SIZE( io_comm, actual_io_procs, ierr )
    End If

    ! Allocate required data structures for stuff local to this processor,
    ! including those required for ranking the data
    ! Put as much of the data as possible into one array to avoid latencies
    ! in MPI later
    ! Quite the size of this depends on quite what we are writing and with
    ! how much detail.
    ! In an attempt to catch bugs set things that are not used by the certain
    ! options to silly values

    io%W_IND = -10000
    io%Q_IND = -10000
    io%D_IND = -10000

    If ( write_options /= IO_MSDTMP ) Then
      ! Determine the amount of data required for each atom due to the detail
      ! of writing that is requested
      Select Case( write_level )
      Case( 0 )
        ! Coordinates only
        size_local = 1 * 3
      Case( 1 )
        ! Coordinates and velocities
        size_local = 2 * 3
      Case( 2 )
        ! Coordinates, velocities and forces
        size_local = 3 * 3
      Case Default
        size_local = -10000
        error = IO_UNKNOWN_WRITE_LEVEL
      End Select
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_UNKNOWN_WRITE_LEVEL
        Return
      End If
    Else
      size_local = 0
    End If

    ! and what extras are required due to the exact file type.
    io%W_IND = 0
    io%Q_IND = 0
    io%D_IND = 0
    Select Case( write_options )
    Case( IO_RESTART )
      size_local = N_RESTART_DATA + size_local
    Case( IO_HISTORY )
      size_local = N_HISTORY_DATA + size_local
      Select Case( write_level )
      Case( 0 )
        io%W_IND = 4
        io%Q_IND = 5
        io%D_IND = 6
      Case( 1 )
        io%W_IND = 7
        io%Q_IND = 8
        io%D_IND = 9
      Case( 2 )
        io%W_IND = 10
        io%Q_IND = 11
        io%D_IND = 12
      End Select
    Case( IO_HISTORD )
      size_local = N_HISTORD_DATA + size_local
      io%D_IND = 4
    Case( IO_MSDTMP )
      size_local = N_MSDTMP_DATA + size_local
      io%W_IND = 1
      io%Q_IND = 2
    Case Default
      error = IO_UNKNOWN_WRITE_OPTION
    End Select
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_UNKNOWN_WRITE_OPTION
      Return
    End If

    Allocate ( local_data( 1:size_local, 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    ! MPI only knows about Character( Len = 1 ) so will put the names in such an array
    Allocate ( local_name( 1:Len( atom_name ), 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    Allocate ( local_global_indices( 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If

    ! Allocate stuff required only on the I/O processors
    ! Note it's not a standard to pass an unallocated array unless the
    ! dummy argument is allocatable.  Here this is not the case hence allocate
    ! to size zero on procs that will not need the array.
    If ( do_io ) Then
      Allocate ( n_reorg( 0:actual_io_procs - 1 ), Stat = error )
    Else
      Allocate ( n_reorg( 0:-1 ), Stat = error )
    End If
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If

    ! Rank the atoms on this proc
    Allocate ( global_index_rank( 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    Call rank( global_indices( 1:n_atoms ), global_index_rank, error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      Return
    End If

    ! Sort the local data using the ranks and copy into arrays to be used for gathering
    Call sort_local( write_level, write_options, global_index_rank, &
      global_indices, parts, vx, vy, vz, weight, rsd, atom_name, &
      local_global_indices, local_data, local_name )

    Deallocate ( global_index_rank , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If

    ! Find max number of atoms on any proc
    tot_atoms = get_tot_atoms( n_atoms, io%base_comm )

    ! Output the data in batches

    bottom_batch   = 1
    top_batch      = Min( io%batch_size_write, tot_atoms )

    ! For netCDF close the file and reopen it so that the communicator
    ! associated with the file only contains those processors which will actually
    ! do the I/O.  This seems to avoid some problems on the Cray XT series.
    If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
      desc = io%known_files( file_handle )%desc
      Call netcdf_close( desc )
    End If

    If ( do_io ) Then
      If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_open( Trim( io%known_files( file_handle )%name ), desc, io_comm, MPI_INFO_NULL )
        Call netcdf_get_def( desc )
      End If
    End If

    Do

      ! Work out how large a set of GLOBAL indices we can consider at once.
      ! This is controlled by the batch size, which in turn controls the
      ! maximum amount of LOCAL memory we can allocate.  Thus iterate to
      ! get a loose bound on the range which maximises the batch but stays within
      ! the limit, and then use that.  We don't want to try to hard on this to get
      ! a tight limit as otherwise the extra comms could outweigh any performance
      ! gain we get from larger batches.
      ! Note the earlier method had to make the most pessimistic assumption -
      ! that all atoms for this batch belonged to one I/O processor.  This
      ! is particularly poor when there are many writers, partially because
      ! the number of atoms written by each processor will be, on average,
      ! small, partially because the messages used to construct the sorted
      ! arrays will be short.
      this_batch_size = Min( io%batch_size_write, tot_atoms )

      iter = 0
      Do
        iter = iter + 1

        top_batch = Min( bottom_batch + this_batch_size - 1, tot_atoms )

        ! Find where the local data relevant for atoms with global index BOTTOM_BATCH:TOP_BATCH
        ! reside in the local arrays.  As sorted can use binary search
        If ( bottom_batch /= 1 .or. top_batch /= tot_atoms ) Then
          Call batch_limits_set( local_global_indices, bottom_batch, top_batch, &
            local_bottom, local_top )
        Else
          local_bottom  = 1
          local_top     = n_atoms
        End If

        !IJB for testing - if uncommented return behaviour to same as older code.
        !!$       n_gathered = this_batch_size
        !!$       Exit

        ! For this batch size what is the biggest number of atoms ant processor
        ! needs to hold ? First find out within an io_gather_group how many
        ! atoms the I/O processor for that group will hold
        Call MPI_ALLREDUCE( Max( local_top - local_bottom + 1, 0 ), n_gathered, 1, &
          MPI_INTEGER, MPI_SUM, io_gather_comm, ierr )
        ! If we are covering all the remaining atoms can now exit
        If ( top_batch == tot_atoms ) Then
          Exit
        End If
        ! Then max over the I/O processors and broadcast it back over the io_gather_group.
        If ( do_io ) Then
          Call MPI_ALLREDUCE( n_gathered, itmp, 1, &
            MPI_INTEGER, MPI_MAX, io_comm, ierr )
          n_gathered = itmp
        End If
        Call MPI_BCAST( n_gathered, 1, MPI_INTEGER, 0, io_gather_comm, ierr )
        ! Check if the required memory to hold the atoms is too much
        If ( n_gathered > io%batch_size_write ) Then
          ! If it is the last iteration was the last acceptable, so go back
          ! to that and exit.
          this_batch_size = this_batch_size / 2
          top_batch = Min( bottom_batch + this_batch_size - 1, tot_atoms )
          Call batch_limits_set( local_global_indices, bottom_batch, top_batch, &
            local_bottom, local_top )
          Call MPI_ALLREDUCE( local_top - local_bottom + 1, n_gathered, 1, &
            MPI_INTEGER, MPI_SUM, io_gather_comm, ierr )
          Exit
        Else
          ! The batch is till smaller than the limit.  Try again.
          this_batch_size = this_batch_size * 2
        End If
      End Do
      this_batch_size = n_gathered
      error = ierr

      ! Gather the data onto the I/O processors.  Note allocation to full
      ! size only occurs on the I/O processor, as that is where it is needed,
      ! but to keep to the standard alloc to zero on the other procs
      If ( do_io ) Then
        Allocate ( gathered_global_indices( 1:this_batch_size ), Stat = error )
      Else
        Allocate ( gathered_global_indices( 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      If ( do_io ) Then
        Allocate ( gathered_data( 1:Size( local_data, Dim = 1 ), 1:this_batch_size ), &
          Stat = error )
      Else
        Allocate ( gathered_data( 1:Size( local_data, Dim = 1 ), 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      If ( do_io ) Then
        Allocate ( gathered_name( 1:Size( local_name, Dim = 1 ), 1:this_batch_size ), &
          Stat = error )
      Else
        Allocate ( gathered_name( 1:Size( local_name, Dim = 1 ), 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      Call gather_data( io_gather_comm, local_bottom, local_top,               &
        local_global_indices   , local_data   , local_name,    &
        gathered_global_indices, gathered_data, gathered_name, &
        n_gathered, error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        Return
      End If

      itmp=0 ; error=0
      ! Data now on the I/O procs - so they only do work now.
      IO_PROCS_ONLY: If ( do_io ) Then

        ! Sort the data currently on the I/O processors into order of increasing ( local values of
        ! the ) global indices
        Call io_sort( n_gathered, gathered_global_indices, gathered_data, gathered_name, error )

        ! Avoid deadlock problems on error

        If ( ok(io,  error == 0, io_comm ) ) Then

          ! Work out how much data is on each IO processor
          ! after the reorganization in the global sort
          Call how_much_after_reorg( io_comm, n_gathered, n_reorg )

          ! Now allocate the data for use after the global sort
          n_me = n_reorg( me_in_io )
          Allocate ( sorted_indices( 1:n_me ), Stat = error )
          If ( ok(io,  error == 0, io_comm ) ) Then
            Allocate ( sorted_data( 1:size_local, 1:n_me ), Stat = error )
            If ( ok(io,  error == 0, io_comm ) ) Then
              Allocate ( sorted_name( 1:Len( atom_name ), 1:n_me ), Stat = error )
              If ( ok(io,  error == 0, io_comm ) ) Then
                ! And sort the data across the I/O processors
                Call global_sort( io_comm, &
                  n_gathered, gathered_global_indices, gathered_data, gathered_name, &
                  n_reorg   , sorted_indices         , sorted_data  , sorted_name,   &
                  error )

                If ( ok(io,  error == 0, io_comm ) ) Then
                  ! Finally write the damn thing !!
                  Call config_out(io, io_comm, write_level, write_options, file_handle, first_record, &
                    n_me, sorted_indices, sorted_data, sorted_name, error )
                End If

                Deallocate ( sorted_name , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
                Deallocate ( sorted_data , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
                Deallocate ( sorted_indices , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
              Else
                error = IO_ALLOCATION_ERROR
                Exit
              End If
            Else
              error = IO_ALLOCATION_ERROR
              Exit
            End If
          Else
            error = IO_ALLOCATION_ERROR
            Exit
          End If
        End If

      Else

        error = 0

      End If IO_PROCS_ONLY

      Deallocate ( gathered_name , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If
      Deallocate ( gathered_data , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If
      Deallocate ( gathered_global_indices , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      bottom_batch = top_batch + 1

      If ( bottom_batch > tot_atoms ) Then
        Exit
      End If

    End Do

    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      If ( io%global_error_check ) Then
        Call MPI_ALLREDUCE( error, itmp, 1, MPI_INTEGER, MPI_MAX,  io%base_comm, ierr )
        error = itmp
      End If
      Return
    End If

    Deallocate ( local_global_indices , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If
    Deallocate ( local_name , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If
    Deallocate ( local_data , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If

    ! For netCDF reopen the file in its original state
    If ( do_io ) Then
      If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_close( desc )
      End If
    End If

    If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
      Call netcdf_open( Trim( io%known_files( file_handle )%name ), io%known_files( file_handle )%desc, &
        io%base_comm, MPI_INFO_NULL )
      Call netcdf_get_def( io%known_files( file_handle )%desc )
    End If

    ! Free comms
    Call free_io_comm( do_io, io_comm, io_gather_comm )

    ! Leave in sync
    Call MPI_BARRIER( io%base_comm, ierr )

  Contains

    Subroutine sort_local( write_level, write_options, global_index_rank,                    &
        global_indices, parts, vx, vy, vz, weight, rsd, atom_name, &
        local_global_indices, local_data, local_name )

      ! Copy the disjoint arrays into contiguous data structures to avoid latencies when using MPI
      ! At the same time sort them using the ranking array.
      ! Arguments are:
      ! WRITE_LEVEL         : How detailed the output will be
      ! WRITE_OPTIONS       : Identifies exactly the file we are writing
      ! GLOBAL_INDEX_RANK   : The ranking of the global indices of the atoms
      ! GLOBAL_INDICES      : The global indices
      ! ATOM_NAME           : The atoms' names
      ! WEIGHT              : The atoms' weights
      ! CHARGE              : The atoms' charges
      ! RMS                 : The atoms' displacement from its position at t=0
      ! RX, RY, RZ          : The components of the atoms' positions
      ! VX, VY, VZ          : The components of the atoms' velocities
      ! FX, FY, FZ          : The components of the forces on the atoms
      ! LOCAL_GLOBAL_INDICES: The output sorted global indices
      ! LOCAL_DATA          : The output single sorted array for real quantities
      ! LOCAL_NAME          : The output atom names as a 2d array of   Character( Len = 1 )

      Integer                                , Intent( In    ) :: write_level
      Integer                                , Intent( In    ) :: write_options
      Integer             , Dimension( :    ), Intent( In    ) :: global_index_rank
      Integer             , Dimension( :    ), Intent( In    ) :: global_indices
      Character( Len = * ), Dimension( :    ), Intent( In    ) :: atom_name
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: weight
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: rsd
      Type( corePart )    , Dimension( :    ), Intent( In    ) :: parts
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vx
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vy
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vz
      Integer             , Dimension( :    ), Intent(   Out ) :: local_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: local_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: local_name

      Integer :: n_atoms
      Integer :: i, j

      n_atoms = Size( local_data, Dim = 2 )

      ! Sort the global indices
      Do i = 1, n_atoms
        local_global_indices( global_index_rank( i ) ) = global_indices( i )
      End Do

      ! Pack and sort the name array
      Do i = 1, n_atoms
        Do j = 1, Len( atom_name )
          local_name( j, global_index_rank( i ) ) = atom_name( i )( j:j )
        End Do
      End Do

      ! Pack and sort the real arrays
      Do i = 1, n_atoms
        If ( write_options /= IO_MSDTMP ) Then
          ! First be basics - coords, velocities and forces
          local_data( RX_IND, global_index_rank( i ) ) = parts( i )%xxx
          local_data( RY_IND, global_index_rank( i ) ) = parts( i )%yyy
          local_data( RZ_IND, global_index_rank( i ) ) = parts( i )%zzz
          If ( write_level > 0 ) Then
            local_data( VX_IND, global_index_rank( i ) ) = vx( i )
            local_data( VY_IND, global_index_rank( i ) ) = vy( i )
            local_data( VZ_IND, global_index_rank( i ) ) = vz( i )
          End If
          If ( write_level > 1 ) Then
            local_data( FX_IND, global_index_rank( i ) ) = parts( i )%fxx
            local_data( FY_IND, global_index_rank( i ) ) = parts( i )%fyy
            local_data( FZ_IND, global_index_rank( i ) ) = parts( i )%fzz
          End If
        End If

        ! Now whatever extra may be required
        Select Case( write_options )
        Case( IO_HISTORY )
          local_data( io%Q_IND, global_index_rank( i ) ) = parts(i)%chge
          local_data( io%W_IND, global_index_rank( i ) ) = weight( i )
          local_data( io%D_IND, global_index_rank( i ) ) = rsd( i )
        Case( IO_HISTORD )
          local_data( io%D_IND, global_index_rank( i ) ) = rsd( i )
        Case( IO_MSDTMP )
          local_data( io%W_IND, global_index_rank( i ) ) = weight( i )
          local_data( io%Q_IND, global_index_rank( i ) ) = parts(i)%chge
        End Select
      End Do

    End Subroutine sort_local

    Subroutine rank( a, rank_array, error )

      ! Rank a 1d integer array using heapsort
      ! Arguments:
      ! A         : Array to rank
      ! RANK_ARRAY: Returned ranking array

      Integer, Dimension( : ), Intent( In    ) :: a
      Integer, Dimension( : ), Intent(   Out ) :: rank_array
      Integer                                  :: error

      Integer, Dimension( : ), Allocatable :: work

      Integer :: temp_scal

      Integer :: n
      Integer :: index
      Integer :: p, q
      Integer :: i, j

      n = Size( a )

      ! zero sized arrays correction

      If ( n < 1 ) Return

      Allocate ( work( 1:n ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      ! First use heapsort to generate an indexing array.

      Do i = 1, n
        work( i ) = i
      End Do
      p = n / 2 + 1
      q = n
      Do While ( q /= 1 .or. p /= 1 )
        If ( p > 1 ) Then
          p          = p - 1
          index      = work( p )
          temp_scal  = a( index )
        Else
          index      = work( q )
          temp_scal  = a( index )
          work( q )  = work( 1 )
          q = q - 1
        End If
        If ( q == 1 .and. p == 1 ) Then
          work( 1 ) = index
        Else
          i = p
          j = p + p
          Do While ( j <= q )
            If ( j < q ) Then
              If ( a( work( j ) ) < a( work( j + 1 ) ) ) Then
                j = j + 1
              End If
            End If
            If ( temp_scal < a( work( j ) ) ) Then
              work( i ) = work( j )
              i = j
              j = j + j
            Else
              j = q + 1
            End If
          End Do
          work( i ) = index
        End If
      End Do

      ! Now turn the Indexing array into a ranking array

      Do i = 1, n
        rank_array( work( i ) ) = i
      End Do

      Deallocate ( work , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

    End Subroutine rank

    Subroutine gather_data( io_gather_comm, local_bottom, local_top,               &
        local_global_indices   , local_data   , local_name,    &
        gathered_global_indices, gathered_data, gathered_name, &
        n_gathered, error )

      ! Gather the local data onto the rank 0 processor in the gather comm
      ! Arguments are:
      ! IO_GATHER_COMM         : The communicator to gather across
      ! LOCAL_BOTTOM           : First index of the local arrays to use
      ! LOCAL_TOP              : Last index of the local arrays to use
      ! LOCAL_GLOBAL_INDICES   : The global indices of the atoms on this proc
      ! LOCAL_DATA             : The data about the atoms on this proc
      ! LOCAL_NAME             : The names of the atoms on this proc
      ! GATHERED_GLOBAL_INDICES: The gathered global indices
      ! GATHERED_DATA          : The gathered atomic data
      ! GATHERED_NAME          : The gathered atomic names
      ! N_GATHERED             : Amount of data gathered onto the I/O proc
      ! ERROR                  : Zero on successful return

      Integer                                , Intent( In    ) :: io_gather_comm
      Integer                                , Intent( In    ) :: local_bottom
      Integer                                , Intent( In    ) :: local_top
      Integer             , Dimension(    : ), Intent( In    ) :: local_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: local_data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: local_name
      Integer             , Dimension(    : ), Intent(   Out ) :: gathered_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: gathered_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: gathered_name
      Integer                                , Intent(   Out ) :: n_gathered
      Integer                                , Intent(   Out ) :: error

      Integer, Dimension( : ), Allocatable :: n_to_gather
      Integer, Dimension( : ), Allocatable :: displs_for_gather

      Integer :: n_atoms, n_data
      Integer :: n_procs_gather, me_in_gather
      Integer :: this_size
      Integer :: i

      error = 0

      ! Number of atoms in this batch
      n_atoms = local_top - local_bottom + 1
      If ( local_top < local_bottom ) Then
        n_atoms = 0
      End If

      ! Amount of data associated with each atom
      n_data = Size( local_data, Dim = 1 )

      Call MPI_COMM_SIZE( io_gather_comm, n_procs_gather, ierr )
      Call MPI_COMM_RANK( io_gather_comm, me_in_gather  , ierr )

      Allocate ( n_to_gather( 0:n_procs_gather - 1 ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      Allocate ( displs_for_gather( 0:n_procs_gather - 1 ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      ! How many atoms are on each processor in this batch ?
      Call MPI_ALLGATHER( n_atoms, 1, MPI_INTEGER, n_to_gather, 1, MPI_INTEGER, io_gather_comm, ierr )

      ! And hence the displacements for the gather of the actual data
      displs_for_gather( 0 ) = 0
      Do i = 1, Ubound( displs_for_gather, Dim = 1 )
        displs_for_gather( i ) = displs_for_gather( i - 1 ) + n_to_gather( i - 1 )
      End Do

      ! And now the gather have to be careful if gathering zero
      ! atoms from this processor
      If ( n_atoms /= 0 ) Then

        this_size = n_atoms * Size( local_name, Dim = 1 )

        Call MPI_GATHERV( local_name( :, local_bottom:local_top ), this_size, MPI_CHARACTER, &
          gathered_name, n_to_gather * Size( local_name, Dim = 1 ),          &
          displs_for_gather * Size( local_name, Dim = 1 ), MPI_CHARACTER,    &
          0, io_gather_comm, ierr )

        this_size = n_atoms

        Call MPI_GATHERV( local_global_indices( local_bottom:local_top ), this_size, MPI_INTEGER,&
          gathered_global_indices, n_to_gather, displs_for_gather, MPI_INTEGER,  &
          0, io_gather_comm, ierr )

        this_size = n_atoms * Size( local_data, Dim = 1 )

        Call MPI_GATHERV( local_data( :, local_bottom:local_top ), this_size, wp_mpi, &
          gathered_data, n_to_gather * n_data,                        &
          displs_for_gather * n_data, wp_mpi,                         &
          0, io_gather_comm, ierr )

      Else

        this_size = n_atoms * Size( local_name, Dim = 1 )

        Call MPI_GATHERV( local_name, this_size, MPI_CHARACTER,                           &
          gathered_name, n_to_gather * Size( local_name, Dim = 1 ),       &
          displs_for_gather * Size( local_name, Dim = 1 ), MPI_CHARACTER, &
          0, io_gather_comm, ierr )

        this_size = n_atoms

        Call MPI_GATHERV( local_global_indices, this_size, MPI_INTEGER,                         &
          gathered_global_indices, n_to_gather, displs_for_gather, MPI_INTEGER, &
          0, io_gather_comm, ierr )

        this_size = n_atoms * Size( local_data, Dim = 1 )

        Call MPI_GATHERV( local_data, this_size, wp_mpi,       &
          gathered_data, n_to_gather * n_data, &
          displs_for_gather * n_data, wp_mpi,  &
          0, io_gather_comm, ierr )
      End If

      n_gathered = Sum( n_to_gather )

      Deallocate ( displs_for_gather , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( n_to_gather       , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

    End Subroutine gather_data

    Subroutine io_sort( n, global_indices, data, name, error )

      ! Sort the data according to the indices in GLOBAL_INDICES
      ! Arguments are:
      ! N             : Amount of data on this proc
      ! GLOBAL_INDICES: The global indices
      ! DATA          : the atomic data
      ! NAME          : The atomic names
      ! ERROR         : zero on successful return.

      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( InOut ) :: global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( InOut ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( InOut ) :: name
      Integer                                , Intent(   Out ) :: error

      Integer :: i,nd,nn

      Integer, Dimension( : ), Allocatable :: ranks

      Integer,              Dimension( : ),    Allocatable :: gi
      Real( Kind = wp ),    Dimension( :, : ), Allocatable :: dt
      Character( Len = 1 ), Dimension( :, : ), Allocatable :: nm

      error = 0

      If ( n /= 0 ) Then

        Allocate ( ranks( 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If

        Call rank( global_indices( 1:n ), ranks, error )
        If ( error /= 0 ) Then
          Return
        End If

        !         global_indices( [ ranks ] ) = global_indices( 1:n )
        Allocate ( gi( 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          gi( ranks( i ) ) = global_indices( i )
        End Do
        global_indices = gi
        Deallocate ( gi , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        !         data( :, [ ranks ] ) = data( :, 1:n )
        nd = Size( data, Dim = 1 )
        Allocate ( dt( 1:nd, 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          dt( :, ranks( i ) ) = data( :, i )
        End Do
        data = dt
        Deallocate ( dt , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        !         name( :, [ ranks ] ) = name( :, 1:n )
        nn = Size( name, Dim = 1 )
        Allocate ( nm( 1:nn, 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          nm( :, ranks( i ) ) = name( :, i )
        End Do
        name = nm
        Deallocate ( nm , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        Deallocate ( ranks , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

      End If

    End Subroutine io_sort

    Subroutine how_much_after_reorg( io_comm, n, n_reorg )

      ! Work out how many atoms will be held by each I/O processor after the
      ! data reorganization which is part of the global sort.
      ! Arguments are:
      ! IO_COMM: The communicator which the I/O procs span
      ! N      : The number of atoms in this batch for this I/O processor
      ! N_REORG: The number of atoms that are on each processor after the reorganisation

      Integer                 , Intent( In    ) :: io_comm
      Integer                 , Intent( In    ) :: n
      Integer, Dimension( 0: ), Intent(   Out ) :: n_reorg

      Integer :: n_in_io_comm
      Integer :: n_total
      Integer :: n_av, n_left

      Call MPI_COMM_SIZE( io_comm, n_in_io_comm, ierr )

      Call MPI_ALLREDUCE( n, n_total, 1, MPI_INTEGER, MPI_SUM, io_comm, ierr )

      n_av    = n_total / n_in_io_comm
      n_left  = n_total - n_av * n_in_io_comm

      n_reorg( 0:n_in_io_comm - 1 ) = n_av
      n_reorg( 0:n_left - 1       ) = n_reorg( 0:n_left - 1 ) + 1

    End Subroutine how_much_after_reorg

    Subroutine global_sort( io_comm,                                           &
        n      , indices       , data       , name       , &
        n_reorg, sorted_indices, sorted_data, sorted_name, &
        error )

      ! Sort the data according to the indexing in GATHERED_GLOBAL_INDICES which
      ! is spread across the processors spanned by IO_COMM.  The local arrays
      ! are assumed already sorted.
      ! The arguments are:
      ! IO_COMM       : The communicator
      ! N             : The number of atoms in the batch currently on this processor
      ! INDICES       : The atomic indices
      ! DATA          : The atomic data
      ! NAME          : The atoms' names
      ! N_REORG       : The number of atoms on each of the processors after the reorganization
      ! SORTEio%D_INDICES: The final, globally sorted indices
      ! SORTED_DATA   : The final, globally sorted atomic data
      ! SORTED_NAME   : The final, globally sorted atom names
      ! ERROR         : zero on successful return.

      Integer                                , Intent( In    ) :: io_comm
      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( In    ) :: indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: name
      Integer             , Dimension(   0: ), Intent( In    ) :: n_reorg
      Integer             , Dimension(    : ), Intent(   Out ) :: sorted_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: sorted_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: sorted_name
      Integer                                , Intent(   Out ) :: error

      Integer, Dimension( : ), Allocatable :: first_atom
      Integer, Dimension( : ), Allocatable :: to_send
      Integer, Dimension( : ), Allocatable :: displs_send
      Integer, Dimension( : ), Allocatable :: to_recv
      Integer, Dimension( : ), Allocatable :: displs_recv

      Integer :: n_in_io_comm, me_in_io_comm
      Integer :: my_first, first_in_batch
      Integer :: at1
      Integer :: leng, ndat
      Integer :: i, j

      error = 0

      Call MPI_COMM_SIZE( io_comm,  n_in_io_comm, ierr )
      Call MPI_COMM_RANK( io_comm, me_in_io_comm, ierr )

      ! Find the first atom in the batch
      If ( n /= 0 ) Then
        my_first = indices( 1 )
      Else
        my_first = Huge( my_first )
      End If
      Call MPI_ALLREDUCE( my_first, first_in_batch, 1, MPI_INTEGER, MPI_MIN, io_comm, ierr )

      ! Find the first atom on each proc - remember the atoms are being outputted
      ! in one big contiguous block of indices
      Allocate ( first_atom( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      first_atom( 0 ) = first_in_batch
      Do i = 1, n_in_io_comm - 1
        first_atom( i ) = first_atom( i - 1 ) + n_reorg( i - 1 )
      End Do

      ! Now work out how much data this proc will send to each of the other
      ! processors
      Allocate ( to_send( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      to_send = 0
      Do i = 1, n
        at1 = indices( i )
        Do j = 1, n_in_io_comm - 1
          If ( at1 < first_atom( j ) ) Then
            Exit
          End If
        End Do
        to_send( j - 1 ) = to_send( j - 1 ) + 1
      End Do

      ! Displacements for sending in alltoallv later on
      Allocate ( displs_send( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      displs_send( 0 ) = 0
      Do i = 1, n_in_io_comm - 1
        displs_send( i ) = displs_send(  i - 1 ) + to_send( i - 1 )
      End Do

      ! The amount to be received in the alltoallv later on is simply the amount sent
      Allocate ( to_recv( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      Call MPI_ALLTOALL( to_send, 1, MPI_INTEGER, to_recv, 1, MPI_INTEGER, io_comm, ierr )

      ! And now the displacements for receiving
      Allocate ( displs_recv( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      displs_recv( 0 ) = 0
      Do i = 1, n_in_io_comm - 1
        displs_recv( i ) = displs_recv(  i - 1 ) + to_recv( i - 1 )
      End Do

      ! And reorganize the data !! ( at long last ... )
      leng = Size( name, Dim = 1 )
      ndat = Size( data, Dim = 1 )
      Call MPI_ALLTOALLV(        name, leng * to_send, leng * displs_send, MPI_CHARACTER, &
        sorted_name, leng * to_recv, leng * displs_recv, MPI_CHARACTER, &
        io_comm, ierr )
      Call MPI_ALLTOALLV(        data, ndat * to_send, ndat * displs_send, wp_mpi, &
        sorted_data, ndat * to_recv, ndat * displs_recv, wp_mpi, &
        io_comm, ierr )
      Call MPI_ALLTOALLV(        indices, to_send, displs_send, MPI_INTEGER, &
        sorted_indices, to_recv, displs_recv, MPI_INTEGER, &
        io_comm, ierr )

      Deallocate ( displs_recv , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( to_recv     , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( displs_send , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( to_send     , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( first_atom  , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      If ( n_reorg( me_in_io_comm ) /= 0 ) Then

        ! One final sort ....
        Call io_sort( n_reorg( me_in_io_comm ), sorted_indices, sorted_data, sorted_name, error )

      End If

      ! Phew !

    End Subroutine global_sort

    Subroutine config_out(io, io_comm, write_level, write_options, file_handle, first_record, &
        n, indices, data, name, error )

      ! Write the config to disk !
      ! Arguments are:
      ! IO_COMM     : The communicator within which we are working.
      ! WRITE_LEVEL : Controls the amount of data to be written
      ! FILE_HANDLE : The MPI-I/O file handle
      ! FIRST_RECORD: the first record to be written to - i.e. where indices( 1 ) should go
      ! N           : The number of atoms in the batch currently on this processor
      ! INDICES     : The final, globally sorted indices
      ! DATA        : The final, globally sorted atomic data
      ! NAME        : The final, globally sorted atom names
      ! ERROR       : zero on successful return.

      Type(io_type), Intent( InOut ) :: io
      Integer                                , Intent( In    ) :: io_comm
      Integer                                , Intent( In    ) :: write_level
      Integer                                , Intent( In    ) :: write_options
      Integer                                , Intent( In    ) :: file_handle
      Integer( Kind = MPI_OFFSET_KIND )      , Intent( In    ) :: first_record
      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( In    ) :: indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: name
      Integer                                , Intent(   Out ) :: error

      Integer( Kind = MPI_OFFSET_KIND ) :: next_rec

      Integer :: n_buff
      Integer :: in_buffer
      Integer :: ln, frame
      Integer :: i, j

      Character( Len = 1 ), Dimension( :, : ), Allocatable :: buffer

      Character( Len = io%rec_size ) :: line

      error = 0

      If ( io%method_write /= IO_WRITE_SORTED_NETCDF ) Then

        ! The io_buffer
        n_buff = Min( n, io%buffer_size_write )
        Allocate ( buffer( 1:io%rec_size, 1:n_buff ), Stat = error )
        If ( .not. ok(io, error == 0, io_comm ) ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If

        ! Loop over the atoms filling up the buffer.  Once full dump to disk.
        in_buffer = 0
        If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then
          next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
            Int( indices( 1 ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND )
        Else
          next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
            Int( indices( 1 ) - 1, MPI_OFFSET_KIND )
        End If

        Do i = 1, n

          ! Atom name - change structure according to whether
          ! the weight, charge and rsd are stored or not
          ! Also note slight pain due to array of character not being the same as character( Len( Size( array ) )
          ! being extremely careful here as MPI only knows about character( Len = 1 ), so the derived
          ! type can ONLY be an array of char ....
          Select Case( write_options )
          Case( IO_RESTART )
            ! Restart file
            Write( line, '( 8a1, i10, a54, a1 )' ) name( :, i ), indices( i ), Repeat( ' ', 54 ), io%lf
          Case( IO_HISTORY )
            ! History File
            Write( line, '( 8a1, i10, 3(1x,f13.6), a12, a1 )' ) name( :, i ), indices( i ), &
              data( io%W_IND, i ), data( io%Q_IND, i ), data( io%D_IND, i ), Repeat( ' ', 12 ), io%lf
          Case( IO_HISTORD )
            ! Short History File
            Write( line, '( 6a1, 4f7.1, a1 )' ) &
              name( 1:6, i ), data( RX_IND, i ), data( RY_IND, i ), data( RZ_IND, i ), data( io%D_IND, i ), io%lf
          Case( IO_MSDTMP )
            ! MSDTMP file
            Write( line, '( 8a1, i10, 1p, 2e13.4, a8, a1 )' )  name( :, i ), indices( i ), &
              data( io%W_IND, i ), data( io%Q_IND, i ), Repeat( ' ', 8 ), io%lf
          End Select
          in_buffer = in_buffer + 1
          Do j = 1, io%rec_size
            buffer( j, in_buffer ) = line( j:j )
          End Do

          If ( in_buffer == n_buff ) Then
            Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
            in_buffer = 0
            If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                Int(1,MPI_OFFSET_KIND)
            Else
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + Int( indices( i ) - 1, MPI_OFFSET_KIND ) + &
                Int(1,MPI_OFFSET_KIND)
            End If
          End If

          If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then

            ! Atomic coordinates
            Write( line, '( 3g20.10, a12, a1 )' ) &
              data( RX_IND, i ), data( RY_IND, i ), data( RZ_IND, i ), Repeat( ' ', 12 ), io%lf
            in_buffer = in_buffer + 1
            Do j = 1, io%rec_size
              buffer( j, in_buffer ) = line( j:j )
            End Do

            If ( in_buffer == n_buff ) Then
              Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
              in_buffer = 0
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                Int(2,MPI_OFFSET_KIND)
            End If

            ! Velocities, if required
            If ( write_level > 0 ) Then
              Write( line, '( 3g20.10, a12, a1 )' ) &
                data( VX_IND, i ), data( VY_IND, i ), data( VZ_IND, i ), Repeat( ' ', 12 ), io%lf
              in_buffer = in_buffer + 1
              Do j = 1, io%rec_size
                buffer( j, in_buffer ) = line( j:j )
              End Do

              If ( in_buffer == n_buff ) Then
                Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
                in_buffer = 0
                next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                  Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                  Int(3,MPI_OFFSET_KIND)
              End If
            End If

            ! Forces, if required
            If ( write_level > 1 ) Then
              Write( line, '( 3g20.10, a12, a1 )' ) &
                data( FX_IND, i ), data( FY_IND, i ), data( FZ_IND, i ), Repeat( ' ', 12 ), io%lf
              in_buffer = in_buffer + 1
              Do j = 1, io%rec_size
                buffer( j, in_buffer ) = line( j:j )
              End Do

              If ( in_buffer == n_buff ) Then
                Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
                in_buffer = 0
                next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                  Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                  Int(4,MPI_OFFSET_KIND)
              End If
            End If

          End If

        End Do

        If ( in_buffer /= 0 ) Then
          Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
          in_buffer = 0
        End If

        Deallocate ( buffer , Stat = error )
        If ( .not. ok(io, error == 0, io_comm ) ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

      Else

        ! netCDF write
        ln = Size( name, Dim = 1 )
        frame = Int(first_record,Kind(frame))

        Select Case( write_options )
        Case( IO_RESTART )
          ! Restart file
          Call netcdf_put_var( 'atomnames', desc, name   , [ indices( 1 ), frame ], [ ln, n, 1 ] )
          Call netcdf_put_var( 'indices'  , desc, indices, [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_HISTORY )
          ! History File
          If (frame == 1) Then
            Call netcdf_put_var( 'atomnames', desc, name              , [ indices( 1 ), frame ], [ ln, n, 1 ] )
            Call netcdf_put_var( 'indices'  , desc, indices           , [ indices( 1 ), frame ], [     n, 1 ] )
            Call netcdf_put_var( 'masses'   , desc, data( io%W_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
            Call netcdf_put_var( 'charges'  , desc, data( io%Q_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
          End If
          Call netcdf_put_var( 'rsd'      , desc, data( io%D_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_HISTORD )
          ! Short History File
          If (frame == 1) &
            Call netcdf_put_var( 'atomnames', desc, name              , [ indices( 1 ), frame ], [ ln, n, 1 ] )
          Call netcdf_put_var( 'rsd'      , desc, data( io%D_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_MSDTMP )

          ! Impossible for the time being

        End Select

        If ( write_options /= IO_MSDTMP ) Then
          Select Case( write_level )
          Case( 0 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          Case( 1 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'velocities' , desc, data( VX_IND:VZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          Case( 2 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'velocities' , desc, data( VX_IND:VZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'forces'     , desc, data( FX_IND:FZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          End Select
        End If

      End If

    End Subroutine config_out

  End Subroutine io_write_sorted_file_parts

  Subroutine io_write_sorted_file_parts_subset(io, file_handle, write_level, write_options, first_record,   &
      n_atoms, global_indices, atom_name, weight, rsd, &
      parts, vx, vy, vz, subset,               &
      error )



    ! Write a single atomic configuration using MPI-I/O to the file corresponding to FILE_HANDLE.
    ! Arguments are:
    ! FILE_HANDLE   : The file handle
    ! WRITE_LEVEL   : How detailed the output will be, i.e. if velocities and forces are written
    ! WRITE_OPTIONS : Identifies quite which version of the closely related files is being written
    ! FIRST_RECORD  : Where the first output will be in the file
    ! N_ATOMS       : The number of atoms stored by this processor
    ! GLOBAL_INDICES: The global index number of the atom
    ! ATOM_NAME     : The atoms' names
    ! WEIGHT        : The atoms' weights
    ! CHARGE        : The atoms' charges
    ! RSD           : The atoms' RMS displacement (??)
    ! RX, RY, RZ    : The components of the atoms' positions
    ! VX, VY, VZ    : The components of the atoms' velocities
    ! FX, FY, FZ    : The components of the forces on the atoms
    ! ERROR         : Return code.  Non-zero indicates an error

    Type(io_type), Intent( InOut ) :: io
    Integer                             , Intent( In    ) :: file_handle
    Integer                             , Intent( In    ) :: write_level
    Integer                             , Intent( In    ) :: write_options
    Integer( Kind = MPI_OFFSET_KIND )   , Intent( In    ) :: first_record
    Integer                             , Intent( In    ) :: n_atoms
    Integer             , Dimension( : ), Intent( In    ) :: global_indices
    Character( Len = * ), Dimension( : ), Intent( In    ) :: atom_name
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: weight
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: rsd
    Type( corePart )    , Dimension( : ), Intent( In    ) :: parts
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vx
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vy
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vz
    Integer                             , Intent( In    ) :: subset
    Integer                             , Intent(   Out ) :: error

    Type( netcdf_desc ) :: desc

    Real( Kind = wp ), Dimension( :, : ), Allocatable :: local_data
    Real( Kind = wp ), Dimension( :, : ), Allocatable :: gathered_data
    Real( Kind = wp ), Dimension( :, : ), Allocatable :: sorted_data

    Integer, Dimension( : ), Allocatable :: global_index_rank
    Integer, Dimension( : ), Allocatable :: local_global_indices
    Integer, Dimension( : ), Allocatable :: gathered_global_indices
    Integer, Dimension( : ), Allocatable :: sorted_indices
    Integer, Dimension( : ), Allocatable :: n_reorg

    Integer :: io_comm, io_gather_comm
    Integer :: actual_io_procs
    Integer :: size_local
    Integer :: bottom_batch, top_batch
    Integer :: local_bottom, local_top
    Integer :: this_batch_size
    Integer :: n_gathered
    Integer :: tot_atoms
    Integer :: me_in_io, n_me
    Integer :: itmp
    Integer :: iter

    Logical :: do_io

    Character( Len = 1 ), Dimension( :, : ), Allocatable :: local_name
    Character( Len = 1 ), Dimension( :, : ), Allocatable :: gathered_name
    Character( Len = 1 ), Dimension( :, : ), Allocatable :: sorted_name
    Integer :: ierr
    ! Ever the optimist
    error = 0

    ! Check we have a communicator
    If ( .not. ok(io, io%base_comm /= MPI_COMM_NULL, io%base_comm ) ) Then
      error = IO_BASE_COMM_NOT_SET
      Return
    End If

    ! Create the I/O communicators, IO_COMM which contains only those procs
    ! which will perform I/O, and IO_GATHER_COMM, which contains at rank 0 the
    ! proc which will perform the I/O for all other procs in the give instantiation
    ! of IO_GATHER_COMM.  As a processor can not ( easily ) determine if it is
    ! a member of a given communicator or not ( if it's not in a communicator
    ! it has no data about it ) also return DO_IO to indicate whether this
    ! processor will perform I/O, i.e.  if it is a member of IO_COMM.
    Call split_io_comm( io%base_comm, io%n_io_procs_write, io_comm, io_gather_comm, do_io )
    If ( do_io ) Then
      Call MPI_COMM_RANK( io_comm, me_in_io, ierr )
      Call MPI_COMM_SIZE( io_comm, actual_io_procs, ierr )
    End If

    ! Allocate required data structures for stuff local to this processor,
    ! including those required for ranking the data
    ! Put as much of the data as possible into one array to avoid latencies
    ! in MPI later
    ! Quite the size of this depends on quite what we are writing and with
    ! how much detail.
    ! In an attempt to catch bugs set things that are not used by the certain
    ! options to silly values

    io%W_IND = -10000
    io%Q_IND = -10000
    io%D_IND = -10000

    If ( write_options /= IO_MSDTMP ) Then
      ! Determine the amount of data required for each atom due to the detail
      ! of writing that is requested
      Select Case( write_level )
      Case( 0 )
        ! Coordinates only
        size_local = 1 * 3
      Case( 1 )
        ! Coordinates and velocities
        size_local = 2 * 3
      Case( 2 )
        ! Coordinates, velocities and forces
        size_local = 3 * 3
      Case Default
        size_local = -10000
        error = IO_UNKNOWN_WRITE_LEVEL
      End Select
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_UNKNOWN_WRITE_LEVEL
        Return
      End If
    Else
      size_local = 0
    End If

    ! and what extras are required due to the exact file type.
    io%W_IND = 0
    io%Q_IND = 0
    io%D_IND = 0
    Select Case( write_options )
    Case( IO_RESTART )
      size_local = N_RESTART_DATA + size_local
    Case( IO_HISTORY )
      size_local = N_HISTORY_DATA + size_local
      Select Case( write_level )
      Case( 0 )
        io%W_IND = 4
        io%Q_IND = 5
        io%D_IND = 6
      Case( 1 )
        io%W_IND = 7
        io%Q_IND = 8
        io%D_IND = 9
      Case( 2 )
        io%W_IND = 10
        io%Q_IND = 11
        io%D_IND = 12
      End Select
    Case( IO_HISTORD )
      size_local = N_HISTORD_DATA + size_local
      io%D_IND = 4
    Case( IO_MSDTMP )
      size_local = N_MSDTMP_DATA + size_local
      io%W_IND = 1
      io%Q_IND = 2
    Case Default
      error = IO_UNKNOWN_WRITE_OPTION
    End Select
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_UNKNOWN_WRITE_OPTION
      Return
    End If

    Allocate ( local_data( 1:size_local, 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    ! MPI only knows about Character( Len = 1 ) so will put the names in such an array
    Allocate ( local_name( 1:Len( atom_name ), 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    Allocate ( local_global_indices( 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If

    ! Allocate stuff required only on the I/O processors
    ! Note it's not a standard to pass an unallocated array unless the
    ! dummy argument is allocatable.  Here this is not the case hence allocate
    ! to size zero on procs that will not need the array.
    If ( do_io ) Then
      Allocate ( n_reorg( 0:actual_io_procs - 1 ), Stat = error )
    Else
      Allocate ( n_reorg( 0:-1 ), Stat = error )
    End If
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If

    ! Rank the atoms on this proc
    Allocate ( global_index_rank( 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    Call rank( global_indices( 1:n_atoms ), global_index_rank, error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      Return
    End If

    ! Sort the local data using the ranks and copy into arrays to be used for gathering
    Call sort_local( write_level, write_options, global_index_rank, &
      global_indices, parts, subset, vx, vy, vz, weight, rsd, atom_name, &
      local_global_indices, local_data, local_name )

    Deallocate ( global_index_rank , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If

    ! Find max number of atoms on any proc
    tot_atoms = get_tot_atoms( n_atoms, io%base_comm )

    ! Output the data in batches

    bottom_batch   = 1
    top_batch      = Min( io%batch_size_write, tot_atoms )

    ! For netCDF close the file and reopen it so that the communicator
    ! associated with the file only contains those processors which will actually
    ! do the I/O.  This seems to avoid some problems on the Cray XT series.
    If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
      desc = io%known_files( file_handle )%desc
      Call netcdf_close( desc )
    End If

    If ( do_io ) Then
      If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_open( Trim( io%known_files( file_handle )%name ), desc, io_comm, MPI_INFO_NULL )
        Call netcdf_get_def( desc )
      End If
    End If

    Do

      ! Work out how large a set of GLOBAL indices we can consider at once.
      ! This is controlled by the batch size, which in turn controls the
      ! maximum amount of LOCAL memory we can allocate.  Thus iterate to
      ! get a loose bound on the range which maximises the batch but stays within
      ! the limit, and then use that.  We don't want to try to hard on this to get
      ! a tight limit as otherwise the extra comms could outweigh any performance
      ! gain we get from larger batches.
      ! Note the earlier method had to make the most pessimistic assumption -
      ! that all atoms for this batch belonged to one I/O processor.  This
      ! is particularly poor when there are many writers, partially because
      ! the number of atoms written by each processor will be, on average,
      ! small, partially because the messages used to construct the sorted
      ! arrays will be short.
      this_batch_size = Min( io%batch_size_write, tot_atoms )

      iter = 0
      Do
        iter = iter + 1

        top_batch = Min( bottom_batch + this_batch_size - 1, tot_atoms )

        ! Find where the local data relevant for atoms with global index BOTTOM_BATCH:TOP_BATCH
        ! reside in the local arrays.  As sorted can use binary search
        If ( bottom_batch /= 1 .or. top_batch /= tot_atoms ) Then
          Call batch_limits_set( local_global_indices, bottom_batch, top_batch, &
            local_bottom, local_top )
        Else
          local_bottom  = 1
          local_top     = n_atoms
        End If

        !IJB for testing - if uncommented return behaviour to same as older code.
        !!$       n_gathered = this_batch_size
        !!$       Exit

        ! For this batch size what is the biggest number of atoms ant processor
        ! needs to hold ? First find out within an io_gather_group how many
        ! atoms the I/O processor for that group will hold
        Call MPI_ALLREDUCE( Max( local_top - local_bottom + 1, 0 ), n_gathered, 1, &
          MPI_INTEGER, MPI_SUM, io_gather_comm, ierr )
        ! If we are covering all the remaining atoms can now exit
        If ( top_batch == tot_atoms ) Then
          Exit
        End If
        ! Then max over the I/O processors and broadcast it back over the io_gather_group.
        If ( do_io ) Then
          Call MPI_ALLREDUCE( n_gathered, itmp, 1, &
            MPI_INTEGER, MPI_MAX, io_comm, ierr )
          n_gathered = itmp
        End If
        Call MPI_BCAST( n_gathered, 1, MPI_INTEGER, 0, io_gather_comm, ierr )
        ! Check if the required memory to hold the atoms is too much
        If ( n_gathered > io%batch_size_write ) Then
          ! If it is the last iteration was the last acceptable, so go back
          ! to that and exit.
          this_batch_size = this_batch_size / 2
          top_batch = Min( bottom_batch + this_batch_size - 1, tot_atoms )
          Call batch_limits_set( local_global_indices, bottom_batch, top_batch, &
            local_bottom, local_top )
          Call MPI_ALLREDUCE( local_top - local_bottom + 1, n_gathered, 1, &
            MPI_INTEGER, MPI_SUM, io_gather_comm, ierr )
          Exit
        Else
          ! The batch is till smaller than the limit.  Try again.
          this_batch_size = this_batch_size * 2
        End If
      End Do
      this_batch_size = n_gathered
      error = ierr

      ! Gather the data onto the I/O processors.  Note allocation to full
      ! size only occurs on the I/O processor, as that is where it is needed,
      ! but to keep to the standard alloc to zero on the other procs
      If ( do_io ) Then
        Allocate ( gathered_global_indices( 1:this_batch_size ), Stat = error )
      Else
        Allocate ( gathered_global_indices( 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      If ( do_io ) Then
        Allocate ( gathered_data( 1:Size( local_data, Dim = 1 ), 1:this_batch_size ), &
          Stat = error )
      Else
        Allocate ( gathered_data( 1:Size( local_data, Dim = 1 ), 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      If ( do_io ) Then
        Allocate ( gathered_name( 1:Size( local_name, Dim = 1 ), 1:this_batch_size ), &
          Stat = error )
      Else
        Allocate ( gathered_name( 1:Size( local_name, Dim = 1 ), 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      Call gather_data( io_gather_comm, local_bottom, local_top,               &
        local_global_indices   , local_data   , local_name,    &
        gathered_global_indices, gathered_data, gathered_name, &
        n_gathered, error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        Return
      End If

      itmp=0 ; error=0
      ! Data now on the I/O procs - so they only do work now.
      IO_PROCS_ONLY: If ( do_io ) Then

        ! Sort the data currently on the I/O processors into order of increasing ( local values of
        ! the ) global indices
        Call io_sort( n_gathered, gathered_global_indices, gathered_data, gathered_name, error )

        ! Avoid deadlock problems on error

        If ( ok(io,  error == 0, io_comm ) ) Then

          ! Work out how much data is on each IO processor
          ! after the reorganization in the global sort
          Call how_much_after_reorg( io_comm, n_gathered, n_reorg )

          ! Now allocate the data for use after the global sort
          n_me = n_reorg( me_in_io )
          Allocate ( sorted_indices( 1:n_me ), Stat = error )
          If ( ok(io,  error == 0, io_comm ) ) Then
            Allocate ( sorted_data( 1:size_local, 1:n_me ), Stat = error )
            If ( ok(io,  error == 0, io_comm ) ) Then
              Allocate ( sorted_name( 1:Len( atom_name ), 1:n_me ), Stat = error )
              If ( ok(io,  error == 0, io_comm ) ) Then
                ! And sort the data across the I/O processors
                Call global_sort( io_comm, &
                  n_gathered, gathered_global_indices, gathered_data, gathered_name, &
                  n_reorg   , sorted_indices         , sorted_data  , sorted_name,   &
                  error )

                If ( ok(io,  error == 0, io_comm ) ) Then
                  ! Finally write the damn thing !!
                  Call config_out(io, io_comm, write_level, write_options, file_handle, first_record, &
                    n_me, sorted_indices, sorted_data, sorted_name, error )
                End If

                Deallocate ( sorted_name , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
                Deallocate ( sorted_data , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
                Deallocate ( sorted_indices , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
              Else
                error = IO_ALLOCATION_ERROR
                Exit
              End If
            Else
              error = IO_ALLOCATION_ERROR
              Exit
            End If
          Else
            error = IO_ALLOCATION_ERROR
            Exit
          End If
        End If

      Else

        error = 0

      End If IO_PROCS_ONLY

      Deallocate ( gathered_name , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If
      Deallocate ( gathered_data , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If
      Deallocate ( gathered_global_indices , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      bottom_batch = top_batch + 1

      If ( bottom_batch > tot_atoms ) Then
        Exit
      End If

    End Do

    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      If ( io%global_error_check ) Then
        Call MPI_ALLREDUCE( error, itmp, 1, MPI_INTEGER, MPI_MAX,  io%base_comm, ierr )
        error = itmp
      End If
      Return
    End If

    Deallocate ( local_global_indices , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If
    Deallocate ( local_name , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If
    Deallocate ( local_data , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If

    ! For netCDF reopen the file in its original state
    If ( do_io ) Then
      If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_close( desc )
      End If
    End If

    If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
      Call netcdf_open( Trim( io%known_files( file_handle )%name ), io%known_files( file_handle )%desc, &
        io%base_comm, MPI_INFO_NULL )
      Call netcdf_get_def( io%known_files( file_handle )%desc )
    End If

    ! Free comms
    Call free_io_comm( do_io, io_comm, io_gather_comm )

    ! Leave in sync
    Call MPI_BARRIER( io%base_comm, ierr )

  Contains

    Subroutine sort_local( write_level, write_options, global_index_rank,                    &
        global_indices, parts, subset, vx, vy, vz, weight, rsd, atom_name, &
        local_global_indices, local_data, local_name )

      ! Copy the disjoint arrays into contiguous data structures to avoid latencies when using MPI
      ! At the same time sort them using the ranking array.
      ! Arguments are:
      ! WRITE_LEVEL         : How detailed the output will be
      ! WRITE_OPTIONS       : Identifies exactly the file we are writing
      ! GLOBAL_INDEX_RANK   : The ranking of the global indices of the atoms
      ! GLOBAL_INDICES      : The global indices
      ! ATOM_NAME           : The atoms' names
      ! WEIGHT              : The atoms' weights
      ! CHARGE              : The atoms' charges
      ! RMS                 : The atoms' displacement from its position at t=0
      ! RX, RY, RZ          : The components of the atoms' positions
      ! VX, VY, VZ          : The components of the atoms' velocities
      ! FX, FY, FZ          : The components of the forces on the atoms
      ! LOCAL_GLOBAL_INDICES: The output sorted global indices
      ! LOCAL_DATA          : The output single sorted array for real quantities
      ! LOCAL_NAME          : The output atom names as a 2d array of   Character( Len = 1 )

      Integer                                , Intent( In    ) :: write_level
      Integer                                , Intent( In    ) :: write_options
      Integer             , Dimension( :    ), Intent( In    ) :: global_index_rank
      Integer             , Dimension( :    ), Intent( In    ) :: global_indices
      Character( Len = * ), Dimension( :    ), Intent( In    ) :: atom_name
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: weight
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: rsd
      Type( corePart )    , Dimension( :    ), Intent( In    ) :: parts
      Integer             ,                    Intent( In    ) :: subset
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vx
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vy
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vz
      Integer             , Dimension( :    ), Intent(   Out ) :: local_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: local_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: local_name

      Integer :: n_atoms
      Integer :: i, j

      n_atoms = Size( local_data, Dim = 2 )

      ! Sort the global indices
      Do i = 1, n_atoms
        local_global_indices( global_index_rank( i ) ) = global_indices( i )
      End Do

      ! Pack and sort the name array
      Do i = 1, n_atoms
        Do j = 1, Len( atom_name )
          local_name( j, global_index_rank( i ) ) = atom_name( i )( j:j )
        End Do
      End Do

      ! Pack and sort the real arrays
      Do i = 1, n_atoms
        If ( write_options /= IO_MSDTMP ) Then
          ! First be basics - coords, velocities and forces
          If( IAND(subset, IO_SUBSET_POSITIONS) > 0) Then
            local_data( RX_IND, global_index_rank( i ) ) = parts( i )%xxx
            local_data( RY_IND, global_index_rank( i ) ) = parts( i )%yyy
            local_data( RZ_IND, global_index_rank( i ) ) = parts( i )%zzz
          End If
          If ( write_level > 0 ) Then
            local_data( VX_IND, global_index_rank( i ) ) = vx( i )
            local_data( VY_IND, global_index_rank( i ) ) = vy( i )
            local_data( VZ_IND, global_index_rank( i ) ) = vz( i )
          End If
          If ( write_level > 1 .and. IAND(subset, IO_SUBSET_FORCES) > 0 ) Then
            local_data( FX_IND, global_index_rank( i ) ) = parts( i )%fxx
            local_data( FY_IND, global_index_rank( i ) ) = parts( i )%fyy
            local_data( FZ_IND, global_index_rank( i ) ) = parts( i )%fzz
          End If
        End If

        ! Now whatever extra may be required
        Select Case( write_options )
        Case( IO_HISTORY )
          If(IAND(subset, IO_SUBSET_CHARGES) > 0 ) Then
            local_data( io%Q_IND, global_index_rank( i ) ) = parts( i )%chge
          End If
          local_data( io%W_IND, global_index_rank( i ) ) = weight( i )
          local_data( io%D_IND, global_index_rank( i ) ) = rsd( i )
        Case( IO_HISTORD )
          local_data( io%D_IND, global_index_rank( i ) ) = rsd( i )
        Case( IO_MSDTMP )
          local_data( io%W_IND, global_index_rank( i ) ) = weight( i )
          If(IAND(subset, IO_SUBSET_CHARGES) > 0 ) Then
            local_data( io%Q_IND, global_index_rank( i ) ) = parts( i )%chge
          End If
        End Select
      End Do

    End Subroutine sort_local

    Subroutine rank( a, rank_array, error )

      ! Rank a 1d integer array using heapsort
      ! Arguments:
      ! A         : Array to rank
      ! RANK_ARRAY: Returned ranking array

      Integer, Dimension( : ), Intent( In    ) :: a
      Integer, Dimension( : ), Intent(   Out ) :: rank_array
      Integer                                  :: error

      Integer, Dimension( : ), Allocatable :: work

      Integer :: temp_scal

      Integer :: n
      Integer :: index
      Integer :: p, q
      Integer :: i, j

      n = Size( a )

      ! zero sized arrays correction

      If ( n < 1 ) Return

      Allocate ( work( 1:n ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      ! First use heapsort to generate an indexing array.

      Do i = 1, n
        work( i ) = i
      End Do
      p = n / 2 + 1
      q = n
      Do While ( q /= 1 .or. p /= 1 )
        If ( p > 1 ) Then
          p          = p - 1
          index      = work( p )
          temp_scal  = a( index )
        Else
          index      = work( q )
          temp_scal  = a( index )
          work( q )  = work( 1 )
          q = q - 1
        End If
        If ( q == 1 .and. p == 1 ) Then
          work( 1 ) = index
        Else
          i = p
          j = p + p
          Do While ( j <= q )
            If ( j < q ) Then
              If ( a( work( j ) ) < a( work( j + 1 ) ) ) Then
                j = j + 1
              End If
            End If
            If ( temp_scal < a( work( j ) ) ) Then
              work( i ) = work( j )
              i = j
              j = j + j
            Else
              j = q + 1
            End If
          End Do
          work( i ) = index
        End If
      End Do

      ! Now turn the Indexing array into a ranking array

      Do i = 1, n
        rank_array( work( i ) ) = i
      End Do

      Deallocate ( work , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

    End Subroutine rank

    Subroutine gather_data( io_gather_comm, local_bottom, local_top,               &
        local_global_indices   , local_data   , local_name,    &
        gathered_global_indices, gathered_data, gathered_name, &
        n_gathered, error )

      ! Gather the local data onto the rank 0 processor in the gather comm
      ! Arguments are:
      ! IO_GATHER_COMM         : The communicator to gather across
      ! LOCAL_BOTTOM           : First index of the local arrays to use
      ! LOCAL_TOP              : Last index of the local arrays to use
      ! LOCAL_GLOBAL_INDICES   : The global indices of the atoms on this proc
      ! LOCAL_DATA             : The data about the atoms on this proc
      ! LOCAL_NAME             : The names of the atoms on this proc
      ! GATHERED_GLOBAL_INDICES: The gathered global indices
      ! GATHERED_DATA          : The gathered atomic data
      ! GATHERED_NAME          : The gathered atomic names
      ! N_GATHERED             : Amount of data gathered onto the I/O proc
      ! ERROR                  : Zero on successful return

      Integer                                , Intent( In    ) :: io_gather_comm
      Integer                                , Intent( In    ) :: local_bottom
      Integer                                , Intent( In    ) :: local_top
      Integer             , Dimension(    : ), Intent( In    ) :: local_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: local_data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: local_name
      Integer             , Dimension(    : ), Intent(   Out ) :: gathered_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: gathered_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: gathered_name
      Integer                                , Intent(   Out ) :: n_gathered
      Integer                                , Intent(   Out ) :: error

      Integer, Dimension( : ), Allocatable :: n_to_gather
      Integer, Dimension( : ), Allocatable :: displs_for_gather

      Integer :: n_atoms, n_data
      Integer :: n_procs_gather, me_in_gather
      Integer :: this_size
      Integer :: i

      error = 0

      ! Number of atoms in this batch
      n_atoms = local_top - local_bottom + 1
      If ( local_top < local_bottom ) Then
        n_atoms = 0
      End If

      ! Amount of data associated with each atom
      n_data = Size( local_data, Dim = 1 )

      Call MPI_COMM_SIZE( io_gather_comm, n_procs_gather, ierr )
      Call MPI_COMM_RANK( io_gather_comm, me_in_gather  , ierr )

      Allocate ( n_to_gather( 0:n_procs_gather - 1 ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      Allocate ( displs_for_gather( 0:n_procs_gather - 1 ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      ! How many atoms are on each processor in this batch ?
      Call MPI_ALLGATHER( n_atoms, 1, MPI_INTEGER, n_to_gather, 1, MPI_INTEGER, io_gather_comm, ierr )

      ! And hence the displacements for the gather of the actual data
      displs_for_gather( 0 ) = 0
      Do i = 1, Ubound( displs_for_gather, Dim = 1 )
        displs_for_gather( i ) = displs_for_gather( i - 1 ) + n_to_gather( i - 1 )
      End Do

      ! And now the gather have to be careful if gathering zero
      ! atoms from this processor
      If ( n_atoms /= 0 ) Then

        this_size = n_atoms * Size( local_name, Dim = 1 )

        Call MPI_GATHERV( local_name( :, local_bottom:local_top ), this_size, MPI_CHARACTER, &
          gathered_name, n_to_gather * Size( local_name, Dim = 1 ),          &
          displs_for_gather * Size( local_name, Dim = 1 ), MPI_CHARACTER,    &
          0, io_gather_comm, ierr )

        this_size = n_atoms

        Call MPI_GATHERV( local_global_indices( local_bottom:local_top ), this_size, MPI_INTEGER,&
          gathered_global_indices, n_to_gather, displs_for_gather, MPI_INTEGER,  &
          0, io_gather_comm, ierr )

        this_size = n_atoms * Size( local_data, Dim = 1 )

        Call MPI_GATHERV( local_data( :, local_bottom:local_top ), this_size, wp_mpi, &
          gathered_data, n_to_gather * n_data,                        &
          displs_for_gather * n_data, wp_mpi,                         &
          0, io_gather_comm, ierr )

      Else

        this_size = n_atoms * Size( local_name, Dim = 1 )

        Call MPI_GATHERV( local_name, this_size, MPI_CHARACTER,                           &
          gathered_name, n_to_gather * Size( local_name, Dim = 1 ),       &
          displs_for_gather * Size( local_name, Dim = 1 ), MPI_CHARACTER, &
          0, io_gather_comm, ierr )

        this_size = n_atoms

        Call MPI_GATHERV( local_global_indices, this_size, MPI_INTEGER,                         &
          gathered_global_indices, n_to_gather, displs_for_gather, MPI_INTEGER, &
          0, io_gather_comm, ierr )

        this_size = n_atoms * Size( local_data, Dim = 1 )

        Call MPI_GATHERV( local_data, this_size, wp_mpi,       &
          gathered_data, n_to_gather * n_data, &
          displs_for_gather * n_data, wp_mpi,  &
          0, io_gather_comm, ierr )
      End If

      n_gathered = Sum( n_to_gather )

      Deallocate ( displs_for_gather , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( n_to_gather       , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

    End Subroutine gather_data

    Subroutine io_sort( n, global_indices, data, name, error )

      ! Sort the data according to the indices in GLOBAL_INDICES
      ! Arguments are:
      ! N             : Amount of data on this proc
      ! GLOBAL_INDICES: The global indices
      ! DATA          : the atomic data
      ! NAME          : The atomic names
      ! ERROR         : zero on successful return.

      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( InOut ) :: global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( InOut ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( InOut ) :: name
      Integer                                , Intent(   Out ) :: error

      Integer :: i,nd,nn

      Integer, Dimension( : ), Allocatable :: ranks

      Integer,              Dimension( : ),    Allocatable :: gi
      Real( Kind = wp ),    Dimension( :, : ), Allocatable :: dt
      Character( Len = 1 ), Dimension( :, : ), Allocatable :: nm

      error = 0

      If ( n /= 0 ) Then

        Allocate ( ranks( 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If

        Call rank( global_indices( 1:n ), ranks, error )
        If ( error /= 0 ) Then
          Return
        End If

        !         global_indices( [ ranks ] ) = global_indices( 1:n )
        Allocate ( gi( 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          gi( ranks( i ) ) = global_indices( i )
        End Do
        global_indices = gi
        Deallocate ( gi , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        !         data( :, [ ranks ] ) = data( :, 1:n )
        nd = Size( data, Dim = 1 )
        Allocate ( dt( 1:nd, 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          dt( :, ranks( i ) ) = data( :, i )
        End Do
        data = dt
        Deallocate ( dt , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        !         name( :, [ ranks ] ) = name( :, 1:n )
        nn = Size( name, Dim = 1 )
        Allocate ( nm( 1:nn, 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          nm( :, ranks( i ) ) = name( :, i )
        End Do
        name = nm
        Deallocate ( nm , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        Deallocate ( ranks , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

      End If

    End Subroutine io_sort

    Subroutine how_much_after_reorg( io_comm, n, n_reorg )

      ! Work out how many atoms will be held by each I/O processor after the
      ! data reorganization which is part of the global sort.
      ! Arguments are:
      ! IO_COMM: The communicator which the I/O procs span
      ! N      : The number of atoms in this batch for this I/O processor
      ! N_REORG: The number of atoms that are on each processor after the reorganisation

      Integer                 , Intent( In    ) :: io_comm
      Integer                 , Intent( In    ) :: n
      Integer, Dimension( 0: ), Intent(   Out ) :: n_reorg

      Integer :: n_in_io_comm
      Integer :: n_total
      Integer :: n_av, n_left

      Call MPI_COMM_SIZE( io_comm, n_in_io_comm, ierr )

      Call MPI_ALLREDUCE( n, n_total, 1, MPI_INTEGER, MPI_SUM, io_comm, ierr )

      n_av    = n_total / n_in_io_comm
      n_left  = n_total - n_av * n_in_io_comm

      n_reorg( 0:n_in_io_comm - 1 ) = n_av
      n_reorg( 0:n_left - 1       ) = n_reorg( 0:n_left - 1 ) + 1

    End Subroutine how_much_after_reorg

    Subroutine global_sort( io_comm,                                           &
        n      , indices       , data       , name       , &
        n_reorg, sorted_indices, sorted_data, sorted_name, &
        error )

      ! Sort the data according to the indexing in GATHERED_GLOBAL_INDICES which
      ! is spread across the processors spanned by IO_COMM.  The local arrays
      ! are assumed already sorted.
      ! The arguments are:
      ! IO_COMM       : The communicator
      ! N             : The number of atoms in the batch currently on this processor
      ! INDICES       : The atomic indices
      ! DATA          : The atomic data
      ! NAME          : The atoms' names
      ! N_REORG       : The number of atoms on each of the processors after the reorganization
      ! SORTEio%D_INDICES: The final, globally sorted indices
      ! SORTED_DATA   : The final, globally sorted atomic data
      ! SORTED_NAME   : The final, globally sorted atom names
      ! ERROR         : zero on successful return.

      Integer                                , Intent( In    ) :: io_comm
      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( In    ) :: indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: name
      Integer             , Dimension(   0: ), Intent( In    ) :: n_reorg
      Integer             , Dimension(    : ), Intent(   Out ) :: sorted_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: sorted_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: sorted_name
      Integer                                , Intent(   Out ) :: error

      Integer, Dimension( : ), Allocatable :: first_atom
      Integer, Dimension( : ), Allocatable :: to_send
      Integer, Dimension( : ), Allocatable :: displs_send
      Integer, Dimension( : ), Allocatable :: to_recv
      Integer, Dimension( : ), Allocatable :: displs_recv

      Integer :: n_in_io_comm, me_in_io_comm
      Integer :: my_first, first_in_batch
      Integer :: at1
      Integer :: leng, ndat
      Integer :: i, j

      error = 0

      Call MPI_COMM_SIZE( io_comm,  n_in_io_comm, ierr )
      Call MPI_COMM_RANK( io_comm, me_in_io_comm, ierr )

      ! Find the first atom in the batch
      If ( n /= 0 ) Then
        my_first = indices( 1 )
      Else
        my_first = Huge( my_first )
      End If
      Call MPI_ALLREDUCE( my_first, first_in_batch, 1, MPI_INTEGER, MPI_MIN, io_comm, ierr )

      ! Find the first atom on each proc - remember the atoms are being outputted
      ! in one big contiguous block of indices
      Allocate ( first_atom( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      first_atom( 0 ) = first_in_batch
      Do i = 1, n_in_io_comm - 1
        first_atom( i ) = first_atom( i - 1 ) + n_reorg( i - 1 )
      End Do

      ! Now work out how much data this proc will send to each of the other
      ! processors
      Allocate ( to_send( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      to_send = 0
      Do i = 1, n
        at1 = indices( i )
        Do j = 1, n_in_io_comm - 1
          If ( at1 < first_atom( j ) ) Then
            Exit
          End If
        End Do
        to_send( j - 1 ) = to_send( j - 1 ) + 1
      End Do

      ! Displacements for sending in alltoallv later on
      Allocate ( displs_send( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      displs_send( 0 ) = 0
      Do i = 1, n_in_io_comm - 1
        displs_send( i ) = displs_send(  i - 1 ) + to_send( i - 1 )
      End Do

      ! The amount to be received in the alltoallv later on is simply the amount sent
      Allocate ( to_recv( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      Call MPI_ALLTOALL( to_send, 1, MPI_INTEGER, to_recv, 1, MPI_INTEGER, io_comm, ierr )

      ! And now the displacements for receiving
      Allocate ( displs_recv( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      displs_recv( 0 ) = 0
      Do i = 1, n_in_io_comm - 1
        displs_recv( i ) = displs_recv(  i - 1 ) + to_recv( i - 1 )
      End Do

      ! And reorganize the data !! ( at long last ... )
      leng = Size( name, Dim = 1 )
      ndat = Size( data, Dim = 1 )
      Call MPI_ALLTOALLV(        name, leng * to_send, leng * displs_send, MPI_CHARACTER, &
        sorted_name, leng * to_recv, leng * displs_recv, MPI_CHARACTER, &
        io_comm, ierr )
      Call MPI_ALLTOALLV(        data, ndat * to_send, ndat * displs_send, wp_mpi, &
        sorted_data, ndat * to_recv, ndat * displs_recv, wp_mpi, &
        io_comm, ierr )
      Call MPI_ALLTOALLV(        indices, to_send, displs_send, MPI_INTEGER, &
        sorted_indices, to_recv, displs_recv, MPI_INTEGER, &
        io_comm, ierr )

      Deallocate ( displs_recv , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( to_recv     , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( displs_send , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( to_send     , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( first_atom  , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      If ( n_reorg( me_in_io_comm ) /= 0 ) Then

        ! One final sort ....
        Call io_sort( n_reorg( me_in_io_comm ), sorted_indices, sorted_data, sorted_name, error )

      End If

      ! Phew !

    End Subroutine global_sort

    Subroutine config_out(io, io_comm, write_level, write_options, file_handle, first_record, &
        n, indices, data, name, error )

      ! Write the config to disk !
      ! Arguments are:
      ! IO_COMM     : The communicator within which we are working.
      ! WRITE_LEVEL : Controls the amount of data to be written
      ! FILE_HANDLE : The MPI-I/O file handle
      ! FIRST_RECORD: the first record to be written to - i.e. where indices( 1 ) should go
      ! N           : The number of atoms in the batch currently on this processor
      ! INDICES     : The final, globally sorted indices
      ! DATA        : The final, globally sorted atomic data
      ! NAME        : The final, globally sorted atom names
      ! ERROR       : zero on successful return.

      Type(io_type), Intent( InOut ) :: io
      Integer                                , Intent( In    ) :: io_comm
      Integer                                , Intent( In    ) :: write_level
      Integer                                , Intent( In    ) :: write_options
      Integer                                , Intent( In    ) :: file_handle
      Integer( Kind = MPI_OFFSET_KIND )      , Intent( In    ) :: first_record
      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( In    ) :: indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: name
      Integer                                , Intent(   Out ) :: error

      Integer( Kind = MPI_OFFSET_KIND ) :: next_rec

      Integer :: n_buff
      Integer :: in_buffer
      Integer :: ln, frame
      Integer :: i, j

      Character( Len = 1 ), Dimension( :, : ), Allocatable :: buffer

      Character( Len = io%rec_size ) :: line

      error = 0

      If ( io%method_write /= IO_WRITE_SORTED_NETCDF ) Then

        ! The io_buffer
        n_buff = Min( n, io%buffer_size_write )
        Allocate ( buffer( 1:io%rec_size, 1:n_buff ), Stat = error )
        If ( .not. ok(io, error == 0, io_comm ) ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If

        ! Loop over the atoms filling up the buffer.  Once full dump to disk.
        in_buffer = 0
        If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then
          next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
            Int( indices( 1 ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND )
        Else
          next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
            Int( indices( 1 ) - 1, MPI_OFFSET_KIND )
        End If

        Do i = 1, n

          ! Atom name - change structure according to whether
          ! the weight, charge and rsd are stored or not
          ! Also note slight pain due to array of character not being the same as character( Len( Size( array ) )
          ! being extremely careful here as MPI only knows about character( Len = 1 ), so the derived
          ! type can ONLY be an array of char ....
          Select Case( write_options )
          Case( IO_RESTART )
            ! Restart file
            Write( line, '( 8a1, i10, a54, a1 )' ) name( :, i ), indices( i ), Repeat( ' ', 54 ), io%lf
          Case( IO_HISTORY )
            ! History File
            Write( line, '( 8a1, i10, 3(1x,f13.6), a12, a1 )' ) name( :, i ), indices( i ), &
              data( io%W_IND, i ), data( io%Q_IND, i ), data( io%D_IND, i ), Repeat( ' ', 12 ), io%lf
          Case( IO_HISTORD )
            ! Short History File
            Write( line, '( 6a1, 4f7.1, a1 )' ) &
              name( 1:6, i ), data( RX_IND, i ), data( RY_IND, i ), data( RZ_IND, i ), data( io%D_IND, i ), io%lf
          Case( IO_MSDTMP )
            ! MSDTMP file
            Write( line, '( 8a1, i10, 1p, 2e13.4, a8, a1 )' )  name( :, i ), indices( i ), &
              data( io%W_IND, i ), data( io%Q_IND, i ), Repeat( ' ', 8 ), io%lf
          End Select
          in_buffer = in_buffer + 1
          Do j = 1, io%rec_size
            buffer( j, in_buffer ) = line( j:j )
          End Do

          If ( in_buffer == n_buff ) Then
            Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
            in_buffer = 0
            If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                Int(1,MPI_OFFSET_KIND)
            Else
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + Int( indices( i ) - 1, MPI_OFFSET_KIND ) + &
                Int(1,MPI_OFFSET_KIND)
            End If
          End If

          If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then

            ! Atomic coordinates
            Write( line, '( 3g20.10, a12, a1 )' ) &
              data( RX_IND, i ), data( RY_IND, i ), data( RZ_IND, i ), Repeat( ' ', 12 ), io%lf
            in_buffer = in_buffer + 1
            Do j = 1, io%rec_size
              buffer( j, in_buffer ) = line( j:j )
            End Do

            If ( in_buffer == n_buff ) Then
              Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
              in_buffer = 0
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                Int(2,MPI_OFFSET_KIND)
            End If

            ! Velocities, if required
            If ( write_level > 0 ) Then
              Write( line, '( 3g20.10, a12, a1 )' ) &
                data( VX_IND, i ), data( VY_IND, i ), data( VZ_IND, i ), Repeat( ' ', 12 ), io%lf
              in_buffer = in_buffer + 1
              Do j = 1, io%rec_size
                buffer( j, in_buffer ) = line( j:j )
              End Do

              If ( in_buffer == n_buff ) Then
                Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
                in_buffer = 0
                next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                  Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                  Int(3,MPI_OFFSET_KIND)
              End If
            End If

            ! Forces, if required
            If ( write_level > 1 ) Then
              Write( line, '( 3g20.10, a12, a1 )' ) &
                data( FX_IND, i ), data( FY_IND, i ), data( FZ_IND, i ), Repeat( ' ', 12 ), io%lf
              in_buffer = in_buffer + 1
              Do j = 1, io%rec_size
                buffer( j, in_buffer ) = line( j:j )
              End Do

              If ( in_buffer == n_buff ) Then
                Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
                in_buffer = 0
                next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                  Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                  Int(4,MPI_OFFSET_KIND)
              End If
            End If

          End If

        End Do

        If ( in_buffer /= 0 ) Then
          Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
          in_buffer = 0
        End If

        Deallocate ( buffer , Stat = error )
        If ( .not. ok(io, error == 0, io_comm ) ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

      Else

        ! netCDF write
        ln = Size( name, Dim = 1 )
        frame = Int(first_record,Kind(frame))

        Select Case( write_options )
        Case( IO_RESTART )
          ! Restart file
          Call netcdf_put_var( 'atomnames', desc, name   , [ indices( 1 ), frame ], [ ln, n, 1 ] )
          Call netcdf_put_var( 'indices'  , desc, indices, [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_HISTORY )
          ! History File
          If (frame == 1) Then
            Call netcdf_put_var( 'atomnames', desc, name              , [ indices( 1 ), frame ], [ ln, n, 1 ] )
            Call netcdf_put_var( 'indices'  , desc, indices           , [ indices( 1 ), frame ], [     n, 1 ] )
            Call netcdf_put_var( 'masses'   , desc, data( io%W_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
            Call netcdf_put_var( 'charges'  , desc, data( io%Q_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
          End If
          Call netcdf_put_var( 'rsd'      , desc, data( io%D_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_HISTORD )
          ! Short History File
          If (frame == 1) &
            Call netcdf_put_var( 'atomnames', desc, name              , [ indices( 1 ), frame ], [ ln, n, 1 ] )
          Call netcdf_put_var( 'rsd'      , desc, data( io%D_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_MSDTMP )

          ! Impossible for the time being

        End Select

        If ( write_options /= IO_MSDTMP ) Then
          Select Case( write_level )
          Case( 0 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          Case( 1 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'velocities' , desc, data( VX_IND:VZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          Case( 2 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'velocities' , desc, data( VX_IND:VZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'forces'     , desc, data( FX_IND:FZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          End Select
        End If

      End If

    End Subroutine config_out

  End Subroutine io_write_sorted_file_parts_subset


  Subroutine io_write_sorted_file_arrays(io, file_handle, write_level, write_options, first_record,   &
      n_atoms, global_indices, atom_name, weight, charge, rsd, &
      rx, ry, rz, vx, vy, vz, fx, fy, fz,            &
      error )

    ! Write a single atomic configuration using MPI-I/O to the file corresponding to FILE_HANDLE.
    ! Arguments are:
    ! FILE_HANDLE   : The file handle
    ! WRITE_LEVEL   : How detailed the output will be, i.e. if velocities and forces are written
    ! WRITE_OPTIONS : Identifies quite which version of the closely related files is being written
    ! FIRST_RECORD  : Where the first output will be in the file
    ! N_ATOMS       : The number of atoms stored by this processor
    ! GLOBAL_INDICES: The global index number of the atom
    ! ATOM_NAME     : The atoms' names
    ! WEIGHT        : The atoms' weights
    ! CHARGE        : The atoms' charges
    ! RSD           : The atoms' RMS displacement (??)
    ! RX, RY, RZ    : The components of the atoms' positions
    ! VX, VY, VZ    : The components of the atoms' velocities
    ! FX, FY, FZ    : The components of the forces on the atoms
    ! ERROR         : Return code.  Non-zero indicates an error

    Type(io_type), Intent( InOut ) :: io
    Integer                             , Intent( In    ) :: file_handle
    Integer                             , Intent( In    ) :: write_level
    Integer                             , Intent( In    ) :: write_options
    Integer( Kind = MPI_OFFSET_KIND )   , Intent( In    ) :: first_record
    Integer                             , Intent( In    ) :: n_atoms
    Integer             , Dimension( : ), Intent( In    ) :: global_indices
    Character( Len = * ), Dimension( : ), Intent( In    ) :: atom_name
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: weight
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: charge
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: rsd
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: rx
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: ry
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: rz
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vx
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vy
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: vz
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: fx
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: fy
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: fz
    Integer                             , Intent(   Out ) :: error

    Type( netcdf_desc ) :: desc

    Real( Kind = wp ), Dimension( :, : ), Allocatable :: local_data
    Real( Kind = wp ), Dimension( :, : ), Allocatable :: gathered_data
    Real( Kind = wp ), Dimension( :, : ), Allocatable :: sorted_data

    Integer, Dimension( : ), Allocatable :: global_index_rank
    Integer, Dimension( : ), Allocatable :: local_global_indices
    Integer, Dimension( : ), Allocatable :: gathered_global_indices
    Integer, Dimension( : ), Allocatable :: sorted_indices
    Integer, Dimension( : ), Allocatable :: n_reorg

    Integer :: io_comm, io_gather_comm
    Integer :: actual_io_procs
    Integer :: size_local
    Integer :: bottom_batch, top_batch
    Integer :: local_bottom, local_top
    Integer :: this_batch_size
    Integer :: n_gathered
    Integer :: tot_atoms
    Integer :: me_in_io, n_me
    Integer :: itmp
    Integer :: iter

    Logical :: do_io

    Character( Len = 1 ), Dimension( :, : ), Allocatable :: local_name
    Character( Len = 1 ), Dimension( :, : ), Allocatable :: gathered_name
    Character( Len = 1 ), Dimension( :, : ), Allocatable :: sorted_name
    Integer :: ierr
    ! Ever the optimist
    error = 0

    ! Check we have a communicator
    If ( .not. ok(io, io%base_comm /= MPI_COMM_NULL, io%base_comm ) ) Then
      error = IO_BASE_COMM_NOT_SET
      Return
    End If

    ! Create the I/O communicators, IO_COMM which contains only those procs
    ! which will perform I/O, and IO_GATHER_COMM, which contains at rank 0 the
    ! proc which will perform the I/O for all other procs in the give instantiation
    ! of IO_GATHER_COMM.  As a processor can not ( easily ) determine if it is
    ! a member of a given communicator or not ( if it's not in a communicator
    ! it has no data about it ) also return DO_IO to indicate whether this
    ! processor will perform I/O, i.e.  if it is a member of IO_COMM.
    Call split_io_comm( io%base_comm, io%n_io_procs_write, io_comm, io_gather_comm, do_io )
    If ( do_io ) Then
      Call MPI_COMM_RANK( io_comm, me_in_io, ierr )
      Call MPI_COMM_SIZE( io_comm, actual_io_procs, ierr )
    End If

    ! Allocate required data structures for stuff local to this processor,
    ! including those required for ranking the data
    ! Put as much of the data as possible into one array to avoid latencies
    ! in MPI later
    ! Quite the size of this depends on quite what we are writing and with
    ! how much detail.
    ! In an attempt to catch bugs set things that are not used by the certain
    ! options to silly values

    io%W_IND = -10000
    io%Q_IND = -10000
    io%D_IND = -10000

    If ( write_options /= IO_MSDTMP ) Then
      ! Determine the amount of data required for each atom due to the detail
      ! of writing that is requested
      Select Case( write_level )
      Case( 0 )
        ! Coordinates only
        size_local = 1 * 3
      Case( 1 )
        ! Coordinates and velocities
        size_local = 2 * 3
      Case( 2 )
        ! Coordinates, velocities and forces
        size_local = 3 * 3
      Case Default
        size_local = -10000
        error = IO_UNKNOWN_WRITE_LEVEL
      End Select
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_UNKNOWN_WRITE_LEVEL
        Return
      End If
    Else
      size_local = 0
    End If

    ! and what extras are required due to the exact file type.
    io%W_IND = 0
    io%Q_IND = 0
    io%D_IND = 0
    Select Case( write_options )
    Case( IO_RESTART )
      size_local = N_RESTART_DATA + size_local
    Case( IO_HISTORY )
      size_local = N_HISTORY_DATA + size_local
      Select Case( write_level )
      Case( 0 )
        io%W_IND = 4
        io%Q_IND = 5
        io%D_IND = 6
      Case( 1 )
        io%W_IND = 7
        io%Q_IND = 8
        io%D_IND = 9
      Case( 2 )
        io%W_IND = 10
        io%Q_IND = 11
        io%D_IND = 12
      End Select
    Case( IO_HISTORD )
      size_local = N_HISTORD_DATA + size_local
      io%D_IND = 4
    Case( IO_MSDTMP )
      size_local = N_MSDTMP_DATA + size_local
      io%W_IND = 1
      io%Q_IND = 2
    Case Default
      error = IO_UNKNOWN_WRITE_OPTION
    End Select
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_UNKNOWN_WRITE_OPTION
      Return
    End If

    Allocate ( local_data( 1:size_local, 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    ! MPI only knows about Character( Len = 1 ) so will put the names in such an array
    Allocate ( local_name( 1:Len( atom_name ), 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    Allocate ( local_global_indices( 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If

    ! Allocate stuff required only on the I/O processors
    ! Note it's not a standard to pass an unallocated array unless the
    ! dummy argument is allocatable.  Here this is not the case hence allocate
    ! to size zero on procs that will not need the array.
    If ( do_io ) Then
      Allocate ( n_reorg( 0:actual_io_procs - 1 ), Stat = error )
    Else
      Allocate ( n_reorg( 0:-1 ), Stat = error )
    End If
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If

    ! Rank the atoms on this proc
    Allocate ( global_index_rank( 1:n_atoms ), Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_ALLOCATION_ERROR
      Return
    End If
    Call rank( global_indices( 1:n_atoms ), global_index_rank, error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      Return
    End If

    ! Sort the local data using the ranks and copy into arrays to be used for gathering
    Call sort_local( write_level, write_options, global_index_rank, &
      global_indices, rx, ry, rz, vx, vy, vz, fx, fy, fz, weight, charge, rsd, atom_name, &
      local_global_indices, local_data, local_name )

    Deallocate ( global_index_rank , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If

    ! Find max number of atoms on any proc
    tot_atoms = get_tot_atoms( n_atoms, io%base_comm )

    ! Output the data in batches

    bottom_batch   = 1
    top_batch      = Min( io%batch_size_write, tot_atoms )

    ! For netCDF close the file and reopen it so that the communicator
    ! associated with the file only contains those processors which will actually
    ! do the I/O.  This seems to avoid some problems on the Cray XT series.
    If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
      desc = io%known_files( file_handle )%desc
      Call netcdf_close( desc )
    End If

    If ( do_io ) Then
      If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_open( Trim( io%known_files( file_handle )%name ), desc, io_comm, MPI_INFO_NULL )
        Call netcdf_get_def( desc )
      End If
    End If

    Do

      ! Work out how large a set of GLOBAL indices we can consider at once.
      ! This is controlled by the batch size, which in turn controls the
      ! maximum amount of LOCAL memory we can allocate.  Thus iterate to
      ! get a loose bound on the range which maximises the batch but stays within
      ! the limit, and then use that.  We don't want to try to hard on this to get
      ! a tight limit as otherwise the extra comms could outweigh any performance
      ! gain we get from larger batches.
      ! Note the earlier method had to make the most pessimistic assumption -
      ! that all atoms for this batch belonged to one I/O processor.  This
      ! is particularly poor when there are many writers, partially because
      ! the number of atoms written by each processor will be, on average,
      ! small, partially because the messages used to construct the sorted
      ! arrays will be short.
      this_batch_size = Min( io%batch_size_write, tot_atoms )

      iter = 0
      Do
        iter = iter + 1

        top_batch = Min( bottom_batch + this_batch_size - 1, tot_atoms )

        ! Find where the local data relevant for atoms with global index BOTTOM_BATCH:TOP_BATCH
        ! reside in the local arrays.  As sorted can use binary search
        If ( bottom_batch /= 1 .or. top_batch /= tot_atoms ) Then
          Call batch_limits_set( local_global_indices, bottom_batch, top_batch, &
            local_bottom, local_top )
        Else
          local_bottom  = 1
          local_top     = n_atoms
        End If

        !IJB for testing - if uncommented return behaviour to same as older code.
        !!$       n_gathered = this_batch_size
        !!$       Exit

        ! For this batch size what is the biggest number of atoms ant processor
        ! needs to hold ? First find out within an io_gather_group how many
        ! atoms the I/O processor for that group will hold
        Call MPI_ALLREDUCE( Max( local_top - local_bottom + 1, 0 ), n_gathered, 1, &
          MPI_INTEGER, MPI_SUM, io_gather_comm, ierr )
        ! If we are covering all the remaining atoms can now exit
        If ( top_batch == tot_atoms ) Then
          Exit
        End If
        ! Then max over the I/O processors and broadcast it back over the io_gather_group.
        If ( do_io ) Then
          Call MPI_ALLREDUCE( n_gathered, itmp, 1, &
            MPI_INTEGER, MPI_MAX, io_comm, ierr )
          n_gathered = itmp
        End If
        Call MPI_BCAST( n_gathered, 1, MPI_INTEGER, 0, io_gather_comm, ierr )
        ! Check if the required memory to hold the atoms is too much
        If ( n_gathered > io%batch_size_write ) Then
          ! If it is the last iteration was the last acceptable, so go back
          ! to that and exit.
          this_batch_size = this_batch_size / 2
          top_batch = Min( bottom_batch + this_batch_size - 1, tot_atoms )
          Call batch_limits_set( local_global_indices, bottom_batch, top_batch, &
            local_bottom, local_top )
          Call MPI_ALLREDUCE( local_top - local_bottom + 1, n_gathered, 1, &
            MPI_INTEGER, MPI_SUM, io_gather_comm, ierr )
          Exit
        Else
          ! The batch is till smaller than the limit.  Try again.
          this_batch_size = this_batch_size * 2
        End If
      End Do
      this_batch_size = n_gathered
      error = ierr

      ! Gather the data onto the I/O processors.  Note allocation to full
      ! size only occurs on the I/O processor, as that is where it is needed,
      ! but to keep to the standard alloc to zero on the other procs
      If ( do_io ) Then
        Allocate ( gathered_global_indices( 1:this_batch_size ), Stat = error )
      Else
        Allocate ( gathered_global_indices( 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      If ( do_io ) Then
        Allocate ( gathered_data( 1:Size( local_data, Dim = 1 ), 1:this_batch_size ), &
          Stat = error )
      Else
        Allocate ( gathered_data( 1:Size( local_data, Dim = 1 ), 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      If ( do_io ) Then
        Allocate ( gathered_name( 1:Size( local_name, Dim = 1 ), 1:this_batch_size ), &
          Stat = error )
      Else
        Allocate ( gathered_name( 1:Size( local_name, Dim = 1 ), 0:-1 ), Stat = error )
      End If
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      Call gather_data( io_gather_comm, local_bottom, local_top,               &
        local_global_indices   , local_data   , local_name,    &
        gathered_global_indices, gathered_data, gathered_name, &
        n_gathered, error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        Return
      End If

      itmp=0 ; error=0
      ! Data now on the I/O procs - so they only do work now.
      IO_PROCS_ONLY: If ( do_io ) Then

        ! Sort the data currently on the I/O processors into order of increasing ( local values of
        ! the ) global indices
        Call io_sort( n_gathered, gathered_global_indices, gathered_data, gathered_name, error )

        ! Avoid deadlock problems on error

        If ( ok(io,  error == 0, io_comm ) ) Then

          ! Work out how much data is on each IO processor
          ! after the reorganization in the global sort
          Call how_much_after_reorg( io_comm, n_gathered, n_reorg )

          ! Now allocate the data for use after the global sort
          n_me = n_reorg( me_in_io )
          Allocate ( sorted_indices( 1:n_me ), Stat = error )
          If ( ok(io,  error == 0, io_comm ) ) Then
            Allocate ( sorted_data( 1:size_local, 1:n_me ), Stat = error )
            If ( ok(io,  error == 0, io_comm ) ) Then
              Allocate ( sorted_name( 1:Len( atom_name ), 1:n_me ), Stat = error )
              If ( ok(io,  error == 0, io_comm ) ) Then
                ! And sort the data across the I/O processors
                Call global_sort( io_comm, &
                  n_gathered, gathered_global_indices, gathered_data, gathered_name, &
                  n_reorg   , sorted_indices         , sorted_data  , sorted_name,   &
                  error )

                If ( ok(io,  error == 0, io_comm ) ) Then
                  ! Finally write the damn thing !!
                  Call config_out(io, io_comm, write_level, write_options, file_handle, first_record, &
                    n_me, sorted_indices, sorted_data, sorted_name, error )
                End If

                Deallocate ( sorted_name , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
                Deallocate ( sorted_data , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
                Deallocate ( sorted_indices , Stat = error )
                If ( .not. ok(io, error == 0, io_comm ) ) Then
                  error = IO_DEALLOCATION_ERROR
                  Exit
                End If
              Else
                error = IO_ALLOCATION_ERROR
                Exit
              End If
            Else
              error = IO_ALLOCATION_ERROR
              Exit
            End If
          Else
            error = IO_ALLOCATION_ERROR
            Exit
          End If
        End If

      Else

        error = 0

      End If IO_PROCS_ONLY

      Deallocate ( gathered_name , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If
      Deallocate ( gathered_data , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If
      Deallocate ( gathered_global_indices , Stat = error )
      If ( .not. ok(io, error == 0, io%base_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      bottom_batch = top_batch + 1

      If ( bottom_batch > tot_atoms ) Then
        Exit
      End If

    End Do

    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      If ( io%global_error_check ) Then
        Call MPI_ALLREDUCE( error, itmp, 1, MPI_INTEGER, MPI_MAX,  io%base_comm, ierr )
        error = itmp
      End If
      Return
    End If

    Deallocate ( local_global_indices , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If
    Deallocate ( local_name , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If
    Deallocate ( local_data , Stat = error )
    If ( .not. ok(io, error == 0, io%base_comm ) ) Then
      error = IO_DEALLOCATION_ERROR
      Return
    End If

    ! For netCDF reopen the file in its original state
    If ( do_io ) Then
      If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
        Call netcdf_close( desc )
      End If
    End If

    If ( io%known_files( file_handle )%method == FILE_NETCDF ) Then
      Call netcdf_open( Trim( io%known_files( file_handle )%name ), io%known_files( file_handle )%desc, &
        io%base_comm, MPI_INFO_NULL )
      Call netcdf_get_def( io%known_files( file_handle )%desc )
    End If

    ! Free comms
    Call free_io_comm( do_io, io_comm, io_gather_comm )

    ! Leave in sync
    Call MPI_BARRIER( io%base_comm, ierr )

  Contains

    Subroutine sort_local( write_level, write_options, global_index_rank,                    &
        global_indices, rx, ry, rz, vx, vy, vz, fx, fy, fz, weight, charge, rsd, atom_name, &
        local_global_indices, local_data, local_name )

      ! Copy the disjoint arrays into contiguous data structures to avoid latencies when using MPI
      ! At the same time sort them using the ranking array.
      ! Arguments are:
      ! WRITE_LEVEL         : How detailed the output will be
      ! WRITE_OPTIONS       : Identifies exactly the file we are writing
      ! GLOBAL_INDEX_RANK   : The ranking of the global indices of the atoms
      ! GLOBAL_INDICES      : The global indices
      ! ATOM_NAME           : The atoms' names
      ! WEIGHT              : The atoms' weights
      ! CHARGE              : The atoms' charges
      ! RMS                 : The atoms' displacement from its position at t=0
      ! RX, RY, RZ          : The components of the atoms' positions
      ! VX, VY, VZ          : The components of the atoms' velocities
      ! FX, FY, FZ          : The components of the forces on the atoms
      ! LOCAL_GLOBAL_INDICES: The output sorted global indices
      ! LOCAL_DATA          : The output single sorted array for real quantities
      ! LOCAL_NAME          : The output atom names as a 2d array of   Character( Len = 1 )

      Integer                                , Intent( In    ) :: write_level
      Integer                                , Intent( In    ) :: write_options
      Integer             , Dimension( :    ), Intent( In    ) :: global_index_rank
      Integer             , Dimension( :    ), Intent( In    ) :: global_indices
      Character( Len = * ), Dimension( :    ), Intent( In    ) :: atom_name
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: weight
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: charge
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: rsd
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: rx
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: ry
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: rz
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vx
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vy
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: vz
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: fx
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: fy
      Real( Kind = wp )   , Dimension( :    ), Intent( In    ) :: fz
      Integer             , Dimension( :    ), Intent(   Out ) :: local_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: local_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: local_name

      Integer :: n_atoms
      Integer :: i, j

      n_atoms = Size( local_data, Dim = 2 )

      ! Sort the global indices
      Do i = 1, n_atoms
        local_global_indices( global_index_rank( i ) ) = global_indices( i )
      End Do

      ! Pack and sort the name array
      Do i = 1, n_atoms
        Do j = 1, Len( atom_name )
          local_name( j, global_index_rank( i ) ) = atom_name( i )( j:j )
        End Do
      End Do

      ! Pack and sort the real arrays
      Do i = 1, n_atoms
        If ( write_options /= IO_MSDTMP ) Then
          ! First be basics - coords, velocities and forces
          local_data( RX_IND, global_index_rank( i ) ) = rx( i )
          local_data( RY_IND, global_index_rank( i ) ) = ry( i )
          local_data( RZ_IND, global_index_rank( i ) ) = rz( i )
          If ( write_level > 0 ) Then
            local_data( VX_IND, global_index_rank( i ) ) = vx( i )
            local_data( VY_IND, global_index_rank( i ) ) = vy( i )
            local_data( VZ_IND, global_index_rank( i ) ) = vz( i )
          End If
          If ( write_level > 1 ) Then
            local_data( FX_IND, global_index_rank( i ) ) = fx( i )
            local_data( FY_IND, global_index_rank( i ) ) = fy( i )
            local_data( FZ_IND, global_index_rank( i ) ) = fz( i ) 
          End If
        End If

        ! Now whatever extra may be required
        Select Case( write_options )
        Case( IO_HISTORY )
          local_data( io%Q_IND, global_index_rank( i ) ) = charge( i )
          local_data( io%W_IND, global_index_rank( i ) ) = weight( i )
          local_data( io%D_IND, global_index_rank( i ) ) = rsd( i )
        Case( IO_HISTORD )
          local_data( io%D_IND, global_index_rank( i ) ) = rsd( i )
        Case( IO_MSDTMP )
          local_data( io%W_IND, global_index_rank( i ) ) = weight( i )
          local_data( io%Q_IND, global_index_rank( i ) ) = charge( i )
        End Select
      End Do

    End Subroutine sort_local

    Subroutine rank( a, rank_array, error )

      ! Rank a 1d integer array using heapsort
      ! Arguments:
      ! A         : Array to rank
      ! RANK_ARRAY: Returned ranking array

      Integer, Dimension( : ), Intent( In    ) :: a
      Integer, Dimension( : ), Intent(   Out ) :: rank_array
      Integer                                  :: error

      Integer, Dimension( : ), Allocatable :: work

      Integer :: temp_scal

      Integer :: n
      Integer :: index
      Integer :: p, q
      Integer :: i, j

      n = Size( a )

      ! zero sized arrays correction

      If ( n < 1 ) Return

      Allocate ( work( 1:n ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      ! First use heapsort to generate an indexing array.

      Do i = 1, n
        work( i ) = i
      End Do
      p = n / 2 + 1
      q = n
      Do While ( q /= 1 .or. p /= 1 )
        If ( p > 1 ) Then
          p          = p - 1
          index      = work( p )
          temp_scal  = a( index )
        Else
          index      = work( q )
          temp_scal  = a( index )
          work( q )  = work( 1 )
          q = q - 1
        End If
        If ( q == 1 .and. p == 1 ) Then
          work( 1 ) = index
        Else
          i = p
          j = p + p
          Do While ( j <= q )
            If ( j < q ) Then
              If ( a( work( j ) ) < a( work( j + 1 ) ) ) Then
                j = j + 1
              End If
            End If
            If ( temp_scal < a( work( j ) ) ) Then
              work( i ) = work( j )
              i = j
              j = j + j
            Else
              j = q + 1
            End If
          End Do
          work( i ) = index
        End If
      End Do

      ! Now turn the Indexing array into a ranking array

      Do i = 1, n
        rank_array( work( i ) ) = i
      End Do

      Deallocate ( work , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

    End Subroutine rank

    Subroutine gather_data( io_gather_comm, local_bottom, local_top,               &
        local_global_indices   , local_data   , local_name,    &
        gathered_global_indices, gathered_data, gathered_name, &
        n_gathered, error )

      ! Gather the local data onto the rank 0 processor in the gather comm
      ! Arguments are:
      ! IO_GATHER_COMM         : The communicator to gather across
      ! LOCAL_BOTTOM           : First index of the local arrays to use
      ! LOCAL_TOP              : Last index of the local arrays to use
      ! LOCAL_GLOBAL_INDICES   : The global indices of the atoms on this proc
      ! LOCAL_DATA             : The data about the atoms on this proc
      ! LOCAL_NAME             : The names of the atoms on this proc
      ! GATHERED_GLOBAL_INDICES: The gathered global indices
      ! GATHERED_DATA          : The gathered atomic data
      ! GATHERED_NAME          : The gathered atomic names
      ! N_GATHERED             : Amount of data gathered onto the I/O proc
      ! ERROR                  : Zero on successful return

      Integer                                , Intent( In    ) :: io_gather_comm
      Integer                                , Intent( In    ) :: local_bottom
      Integer                                , Intent( In    ) :: local_top
      Integer             , Dimension(    : ), Intent( In    ) :: local_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: local_data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: local_name
      Integer             , Dimension(    : ), Intent(   Out ) :: gathered_global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: gathered_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: gathered_name
      Integer                                , Intent(   Out ) :: n_gathered
      Integer                                , Intent(   Out ) :: error

      Integer, Dimension( : ), Allocatable :: n_to_gather
      Integer, Dimension( : ), Allocatable :: displs_for_gather

      Integer :: n_atoms, n_data
      Integer :: n_procs_gather, me_in_gather
      Integer :: this_size
      Integer :: i

      error = 0

      ! Number of atoms in this batch
      n_atoms = local_top - local_bottom + 1
      If ( local_top < local_bottom ) Then
        n_atoms = 0
      End If

      ! Amount of data associated with each atom
      n_data = Size( local_data, Dim = 1 )

      Call MPI_COMM_SIZE( io_gather_comm, n_procs_gather, ierr )
      Call MPI_COMM_RANK( io_gather_comm, me_in_gather  , ierr )

      Allocate ( n_to_gather( 0:n_procs_gather - 1 ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      Allocate ( displs_for_gather( 0:n_procs_gather - 1 ), Stat = error )
      If ( error /= 0 ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If

      ! How many atoms are on each processor in this batch ?
      Call MPI_ALLGATHER( n_atoms, 1, MPI_INTEGER, n_to_gather, 1, MPI_INTEGER, io_gather_comm, ierr )

      ! And hence the displacements for the gather of the actual data
      displs_for_gather( 0 ) = 0
      Do i = 1, Ubound( displs_for_gather, Dim = 1 )
        displs_for_gather( i ) = displs_for_gather( i - 1 ) + n_to_gather( i - 1 )
      End Do

      ! And now the gather have to be careful if gathering zero
      ! atoms from this processor
      If ( n_atoms /= 0 ) Then

        this_size = n_atoms * Size( local_name, Dim = 1 )

        Call MPI_GATHERV( local_name( :, local_bottom:local_top ), this_size, MPI_CHARACTER, &
          gathered_name, n_to_gather * Size( local_name, Dim = 1 ),          &
          displs_for_gather * Size( local_name, Dim = 1 ), MPI_CHARACTER,    &
          0, io_gather_comm, ierr )

        this_size = n_atoms

        Call MPI_GATHERV( local_global_indices( local_bottom:local_top ), this_size, MPI_INTEGER,&
          gathered_global_indices, n_to_gather, displs_for_gather, MPI_INTEGER,  &
          0, io_gather_comm, ierr )

        this_size = n_atoms * Size( local_data, Dim = 1 )

        Call MPI_GATHERV( local_data( :, local_bottom:local_top ), this_size, wp_mpi, &
          gathered_data, n_to_gather * n_data,                        &
          displs_for_gather * n_data, wp_mpi,                         &
          0, io_gather_comm, ierr )

      Else

        this_size = n_atoms * Size( local_name, Dim = 1 )

        Call MPI_GATHERV( local_name, this_size, MPI_CHARACTER,                           &
          gathered_name, n_to_gather * Size( local_name, Dim = 1 ),       &
          displs_for_gather * Size( local_name, Dim = 1 ), MPI_CHARACTER, &
          0, io_gather_comm, ierr )

        this_size = n_atoms

        Call MPI_GATHERV( local_global_indices, this_size, MPI_INTEGER,                         &
          gathered_global_indices, n_to_gather, displs_for_gather, MPI_INTEGER, &
          0, io_gather_comm, ierr )

        this_size = n_atoms * Size( local_data, Dim = 1 )

        Call MPI_GATHERV( local_data, this_size, wp_mpi,       &
          gathered_data, n_to_gather * n_data, &
          displs_for_gather * n_data, wp_mpi,  &
          0, io_gather_comm, ierr )
      End If

      n_gathered = Sum( n_to_gather )

      Deallocate ( displs_for_gather , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( n_to_gather       , Stat = error )
      If ( error /= 0 ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

    End Subroutine gather_data

    Subroutine io_sort( n, global_indices, data, name, error )

      ! Sort the data according to the indices in GLOBAL_INDICES
      ! Arguments are:
      ! N             : Amount of data on this proc
      ! GLOBAL_INDICES: The global indices
      ! DATA          : the atomic data
      ! NAME          : The atomic names
      ! ERROR         : zero on successful return.

      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( InOut ) :: global_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( InOut ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( InOut ) :: name
      Integer                                , Intent(   Out ) :: error

      Integer :: i,nd,nn

      Integer, Dimension( : ), Allocatable :: ranks

      Integer,              Dimension( : ),    Allocatable :: gi
      Real( Kind = wp ),    Dimension( :, : ), Allocatable :: dt
      Character( Len = 1 ), Dimension( :, : ), Allocatable :: nm

      error = 0

      If ( n /= 0 ) Then

        Allocate ( ranks( 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If

        Call rank( global_indices( 1:n ), ranks, error )
        If ( error /= 0 ) Then
          Return
        End If

        !         global_indices( [ ranks ] ) = global_indices( 1:n )
        Allocate ( gi( 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          gi( ranks( i ) ) = global_indices( i )
        End Do
        global_indices = gi
        Deallocate ( gi , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        !         data( :, [ ranks ] ) = data( :, 1:n )
        nd = Size( data, Dim = 1 )
        Allocate ( dt( 1:nd, 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          dt( :, ranks( i ) ) = data( :, i )
        End Do
        data = dt
        Deallocate ( dt , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        !         name( :, [ ranks ] ) = name( :, 1:n )
        nn = Size( name, Dim = 1 )
        Allocate ( nm( 1:nn, 1:n ), Stat = error )
        If ( error /= 0 ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If
        Do i = 1,n
          nm( :, ranks( i ) ) = name( :, i )
        End Do
        name = nm
        Deallocate ( nm , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

        Deallocate ( ranks , Stat = error )
        If ( error /= 0 ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

      End If

    End Subroutine io_sort

    Subroutine how_much_after_reorg( io_comm, n, n_reorg )

      ! Work out how many atoms will be held by each I/O processor after the
      ! data reorganization which is part of the global sort.
      ! Arguments are:
      ! IO_COMM: The communicator which the I/O procs span
      ! N      : The number of atoms in this batch for this I/O processor
      ! N_REORG: The number of atoms that are on each processor after the reorganisation

      Integer                 , Intent( In    ) :: io_comm
      Integer                 , Intent( In    ) :: n
      Integer, Dimension( 0: ), Intent(   Out ) :: n_reorg

      Integer :: n_in_io_comm
      Integer :: n_total
      Integer :: n_av, n_left

      Call MPI_COMM_SIZE( io_comm, n_in_io_comm, ierr )

      Call MPI_ALLREDUCE( n, n_total, 1, MPI_INTEGER, MPI_SUM, io_comm, ierr )

      n_av    = n_total / n_in_io_comm
      n_left  = n_total - n_av * n_in_io_comm

      n_reorg( 0:n_in_io_comm - 1 ) = n_av
      n_reorg( 0:n_left - 1       ) = n_reorg( 0:n_left - 1 ) + 1

    End Subroutine how_much_after_reorg

    Subroutine global_sort(io_comm,                                           &
        n      , indices       , data       , name       , &
        n_reorg, sorted_indices, sorted_data, sorted_name, &
        error )

      ! Sort the data according to the indexing in GATHERED_GLOBAL_INDICES which
      ! is spread across the processors spanned by IO_COMM.  The local arrays
      ! are assumed already sorted.
      ! The arguments are:
      ! IO_COMM       : The communicator
      ! N             : The number of atoms in the batch currently on this processor
      ! INDICES       : The atomic indices
      ! DATA          : The atomic data
      ! NAME          : The atoms' names
      ! N_REORG       : The number of atoms on each of the processors after the reorganization
      ! SORTEio%D_INDICES: The final, globally sorted indices
      ! SORTED_DATA   : The final, globally sorted atomic data
      ! SORTED_NAME   : The final, globally sorted atom names
      ! ERROR         : zero on successful return.

      Integer                                , Intent( In    ) :: io_comm
      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( In    ) :: indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: name
      Integer             , Dimension(   0: ), Intent( In    ) :: n_reorg
      Integer             , Dimension(    : ), Intent(   Out ) :: sorted_indices
      Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: sorted_data
      Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: sorted_name
      Integer                                , Intent(   Out ) :: error

      Integer, Dimension( : ), Allocatable :: first_atom
      Integer, Dimension( : ), Allocatable :: to_send
      Integer, Dimension( : ), Allocatable :: displs_send
      Integer, Dimension( : ), Allocatable :: to_recv
      Integer, Dimension( : ), Allocatable :: displs_recv

      Integer :: n_in_io_comm, me_in_io_comm
      Integer :: my_first, first_in_batch
      Integer :: at1
      Integer :: leng, ndat
      Integer :: i, j

      error = 0

      Call MPI_COMM_SIZE( io_comm,  n_in_io_comm, ierr )
      Call MPI_COMM_RANK( io_comm, me_in_io_comm, ierr )

      ! Find the first atom in the batch
      If ( n /= 0 ) Then
        my_first = indices( 1 )
      Else
        my_first = Huge( my_first )
      End If
      Call MPI_ALLREDUCE( my_first, first_in_batch, 1, MPI_INTEGER, MPI_MIN, io_comm, ierr )

      ! Find the first atom on each proc - remember the atoms are being outputted
      ! in one big contiguous block of indices
      Allocate ( first_atom( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      first_atom( 0 ) = first_in_batch
      Do i = 1, n_in_io_comm - 1
        first_atom( i ) = first_atom( i - 1 ) + n_reorg( i - 1 )
      End Do

      ! Now work out how much data this proc will send to each of the other
      ! processors
      Allocate ( to_send( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      to_send = 0
      Do i = 1, n
        at1 = indices( i )
        Do j = 1, n_in_io_comm - 1
          If ( at1 < first_atom( j ) ) Then
            Exit
          End If
        End Do
        to_send( j - 1 ) = to_send( j - 1 ) + 1
      End Do

      ! Displacements for sending in alltoallv later on
      Allocate ( displs_send( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      displs_send( 0 ) = 0
      Do i = 1, n_in_io_comm - 1
        displs_send( i ) = displs_send(  i - 1 ) + to_send( i - 1 )
      End Do

      ! The amount to be received in the alltoallv later on is simply the amount sent
      Allocate ( to_recv( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      Call MPI_ALLTOALL( to_send, 1, MPI_INTEGER, to_recv, 1, MPI_INTEGER, io_comm, ierr )

      ! And now the displacements for receiving
      Allocate ( displs_recv( 0:n_in_io_comm - 1 ), Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_ALLOCATION_ERROR
        Return
      End If
      displs_recv( 0 ) = 0
      Do i = 1, n_in_io_comm - 1
        displs_recv( i ) = displs_recv(  i - 1 ) + to_recv( i - 1 )
      End Do

      ! And reorganize the data !! ( at long last ... )
      leng = Size( name, Dim = 1 )
      ndat = Size( data, Dim = 1 )
      Call MPI_ALLTOALLV(        name, leng * to_send, leng * displs_send, MPI_CHARACTER, &
        sorted_name, leng * to_recv, leng * displs_recv, MPI_CHARACTER, &
        io_comm, ierr )
      Call MPI_ALLTOALLV(        data, ndat * to_send, ndat * displs_send, wp_mpi, &
        sorted_data, ndat * to_recv, ndat * displs_recv, wp_mpi, &
        io_comm, ierr )
      Call MPI_ALLTOALLV(        indices, to_send, displs_send, MPI_INTEGER, &
        sorted_indices, to_recv, displs_recv, MPI_INTEGER, &
        io_comm, ierr )

      Deallocate ( displs_recv , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( to_recv     , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( displs_send , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( to_send     , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      Deallocate ( first_atom  , Stat = error )
      If ( .not. ok(io, error == 0, io_comm ) ) Then
        error = IO_DEALLOCATION_ERROR
        Return
      End If

      If ( n_reorg( me_in_io_comm ) /= 0 ) Then

        ! One final sort ....
        Call io_sort( n_reorg( me_in_io_comm ), sorted_indices, sorted_data, sorted_name, error )

      End If

      ! Phew !

    End Subroutine global_sort

    Subroutine config_out(io, io_comm, write_level, write_options, file_handle, first_record, &
        n, indices, data, name, error )

      ! Write the config to disk !
      ! Arguments are:
      ! IO_COMM     : The communicator within which we are working.
      ! WRITE_LEVEL : Controls the amount of data to be written
      ! FILE_HANDLE : The MPI-I/O file handle
      ! FIRST_RECORD: the first record to be written to - i.e. where indices( 1 ) should go
      ! N           : The number of atoms in the batch currently on this processor
      ! INDICES     : The final, globally sorted indices
      ! DATA        : The final, globally sorted atomic data
      ! NAME        : The final, globally sorted atom names
      ! ERROR       : zero on successful return.

      Type(io_type), Intent( InOut ) :: io
      Integer                                , Intent( In    ) :: io_comm
      Integer                                , Intent( In    ) :: write_level
      Integer                                , Intent( In    ) :: write_options
      Integer                                , Intent( In    ) :: file_handle
      Integer( Kind = MPI_OFFSET_KIND )      , Intent( In    ) :: first_record
      Integer                                , Intent( In    ) :: n
      Integer             , Dimension(    : ), Intent( In    ) :: indices
      Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: data
      Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: name
      Integer                                , Intent(   Out ) :: error

      Integer( Kind = MPI_OFFSET_KIND ) :: next_rec

      Integer :: n_buff
      Integer :: in_buffer
      Integer :: ln, frame
      Integer :: i, j

      Character( Len = 1 ), Dimension( :, : ), Allocatable :: buffer

      Character( Len = io%rec_size ) :: line

      error = 0

      If ( io%method_write /= IO_WRITE_SORTED_NETCDF ) Then

        ! The io_buffer
        n_buff = Min( n, io%buffer_size_write )
        Allocate ( buffer( 1:io%rec_size, 1:n_buff ), Stat = error )
        If ( .not. ok(io, error == 0, io_comm ) ) Then
          error = IO_ALLOCATION_ERROR
          Return
        End If

        ! Loop over the atoms filling up the buffer.  Once full dump to disk.
        in_buffer = 0
        If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then
          next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
            Int( indices( 1 ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND )
        Else
          next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
            Int( indices( 1 ) - 1, MPI_OFFSET_KIND )
        End If

        Do i = 1, n

          ! Atom name - change structure according to whether
          ! the weight, charge and rsd are stored or not
          ! Also note slight pain due to array of character not being the same as character( Len( Size( array ) )
          ! being extremely careful here as MPI only knows about character( Len = 1 ), so the derived
          ! type can ONLY be an array of char ....
          Select Case( write_options )
          Case( IO_RESTART )
            ! Restart file
            Write( line, '( 8a1, i10, a54, a1 )' ) name( :, i ), indices( i ), Repeat( ' ', 54 ), io%lf
          Case( IO_HISTORY )
            ! History File
            Write( line, '( 8a1, i10, 3(1x,f13.6), a12, a1 )' ) name( :, i ), indices( i ), &
              data( io%W_IND, i ), data( io%Q_IND, i ), data( io%D_IND, i ), Repeat( ' ', 12 ), io%lf
          Case( IO_HISTORD )
            ! Short History File
            Write( line, '( 6a1, 4f7.1, a1 )' ) &
              name( 1:6, i ), data( RX_IND, i ), data( RY_IND, i ), data( RZ_IND, i ), data( io%D_IND, i ), io%lf
          Case( IO_MSDTMP )
            ! MSDTMP file
            Write( line, '( 8a1, i10, 1p, 2e13.4, a8, a1 )' )  name( :, i ), indices( i ), &
              data( io%W_IND, i ), data( io%Q_IND, i ), Repeat( ' ', 8 ), io%lf
          End Select
          in_buffer = in_buffer + 1
          Do j = 1, io%rec_size
            buffer( j, in_buffer ) = line( j:j )
          End Do

          If ( in_buffer == n_buff ) Then
            Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
            in_buffer = 0
            If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                Int(1,MPI_OFFSET_KIND)
            Else
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + Int( indices( i ) - 1, MPI_OFFSET_KIND ) + &
                Int(1,MPI_OFFSET_KIND)
            End If
          End If

          If ( write_options /= IO_MSDTMP .and. write_options /= IO_HISTORD ) Then

            ! Atomic coordinates
            Write( line, '( 3g20.10, a12, a1 )' ) &
              data( RX_IND, i ), data( RY_IND, i ), data( RZ_IND, i ), Repeat( ' ', 12 ), io%lf
            in_buffer = in_buffer + 1
            Do j = 1, io%rec_size
              buffer( j, in_buffer ) = line( j:j )
            End Do

            If ( in_buffer == n_buff ) Then
              Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
              in_buffer = 0
              next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                Int(2,MPI_OFFSET_KIND)
            End If

            ! Velocities, if required
            If ( write_level > 0 ) Then
              Write( line, '( 3g20.10, a12, a1 )' ) &
                data( VX_IND, i ), data( VY_IND, i ), data( VZ_IND, i ), Repeat( ' ', 12 ), io%lf
              in_buffer = in_buffer + 1
              Do j = 1, io%rec_size
                buffer( j, in_buffer ) = line( j:j )
              End Do

              If ( in_buffer == n_buff ) Then
                Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
                in_buffer = 0
                next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                  Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                  Int(3,MPI_OFFSET_KIND)
              End If
            End If

            ! Forces, if required
            If ( write_level > 1 ) Then
              Write( line, '( 3g20.10, a12, a1 )' ) &
                data( FX_IND, i ), data( FY_IND, i ), data( FZ_IND, i ), Repeat( ' ', 12 ), io%lf
              in_buffer = in_buffer + 1
              Do j = 1, io%rec_size
                buffer( j, in_buffer ) = line( j:j )
              End Do

              If ( in_buffer == n_buff ) Then
                Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
                in_buffer = 0
                next_rec = Int( first_record, MPI_OFFSET_KIND ) + &
                  Int( indices( i ) - 1, MPI_OFFSET_KIND ) * Int( write_level + 2, MPI_OFFSET_KIND ) + &
                  Int(4,MPI_OFFSET_KIND)
              End If
            End If

          End If

        End Do

        If ( in_buffer /= 0 ) Then
          Call io_write_batch(io, file_handle, next_rec, in_buffer, buffer )
          in_buffer = 0
        End If

        Deallocate ( buffer , Stat = error )
        If ( .not. ok(io, error == 0, io_comm ) ) Then
          error = IO_DEALLOCATION_ERROR
          Return
        End If

      Else

        ! netCDF write
        ln = Size( name, Dim = 1 )
        frame = Int(first_record,Kind(frame))

        Select Case( write_options )
        Case( IO_RESTART )
          ! Restart file
          Call netcdf_put_var( 'atomnames', desc, name   , [ indices( 1 ), frame ], [ ln, n, 1 ] )
          Call netcdf_put_var( 'indices'  , desc, indices, [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_HISTORY )
          ! History File
          If (frame == 1) Then
            Call netcdf_put_var( 'atomnames', desc, name              , [ indices( 1 ), frame ], [ ln, n, 1 ] )
            Call netcdf_put_var( 'indices'  , desc, indices           , [ indices( 1 ), frame ], [     n, 1 ] )
            Call netcdf_put_var( 'masses'   , desc, data( io%W_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
            Call netcdf_put_var( 'charges'  , desc, data( io%Q_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
          End If
          Call netcdf_put_var( 'rsd'      , desc, data( io%D_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_HISTORD )
          ! Short History File
          If (frame == 1) &
            Call netcdf_put_var( 'atomnames', desc, name              , [ indices( 1 ), frame ], [ ln, n, 1 ] )
          Call netcdf_put_var( 'rsd'      , desc, data( io%D_IND, 1:n ), [ indices( 1 ), frame ], [     n, 1 ] )
        Case( IO_MSDTMP )

          ! Impossible for the time being

        End Select

        If ( write_options /= IO_MSDTMP ) Then
          Select Case( write_level )
          Case( 0 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          Case( 1 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'velocities' , desc, data( VX_IND:VZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          Case( 2 )
            Call netcdf_put_var( 'coordinates', desc, data( RX_IND:RZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'velocities' , desc, data( VX_IND:VZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
            Call netcdf_put_var( 'forces'     , desc, data( FX_IND:FZ_IND, 1:n ), [  1, indices( 1 ), frame ], [ 3, n, 1 ] )
          End Select
        End If

      End If

    End Subroutine config_out

  End Subroutine io_write_sorted_file_arrays

  Subroutine io_write_record( io,io_file_handle, rec_num, record )

    ! Write a record in a shared file
    ! Arguments are:
    ! IO_FILE_HANDLE: The file_handle
    ! REC_NUM       : The record number
    ! RECORD        : the data to be written

    Type(io_type), Intent( InOut ) :: io
    Integer                          , Intent( In    ) :: io_file_handle
    Integer( Kind = MPI_OFFSET_KIND ), Intent( In    ) :: rec_num
    Character( Len = * )             , Intent( In    ) :: record

    Integer               :: status(1:MPI_STATUS_SIZE)
    Integer :: file_handle, ierr
    Integer :: i,j

    Character( Len = 1 ), Dimension( 1:io%rec_size ) :: line

    ! Get the actual file handle or unit in use, rather than the one given to the outside world
    file_handle = io%known_files( io_file_handle )%file_handle

    ! MPI strictly only knows about Character( Len = 1 ) ...
    j = Min( io%rec_size, Len( record ) )
    Do i = 1, j
      line( i ) = record( i:i )
    End Do

    ! Complete line
    If ( j < io%rec_size ) Then
      Do i = j+1, io%rec_size-1
        line( i ) = ' '
      End Do
      line( io%rec_size ) = io%lf
    End If

    If      ( io%method_write == IO_WRITE_UNSORTED_MPIIO  .or. io%method_write == IO_WRITE_SORTED_MPIIO  ) Then

      ! using MPI-I/O
      Call MPI_FILE_WRITE_AT( file_handle, rec_num, line, 1, io%rec_type, status, ierr )

    Else If ( io%method_write == IO_WRITE_UNSORTED_DIRECT .or. io%method_write == IO_WRITE_SORTED_DIRECT ) Then

      ! using parallel direct access
      Write(Unit=file_handle, Fmt=io%forma, Rec=Int(rec_num,li)+Int(1,li)) line

    End If

  End Subroutine io_write_record

  Subroutine io_write_batch(io, io_file_handle, next_rec, in_buffer, buffer )

    ! Write a batch in a shared file
    ! Arguments are:
    ! IO_FILE_HANDLE: The file_handle
    ! NEXT_REC      : The starting record number
    ! IN_BUFFER     : the number of records in the batch
    ! BUFFER        : the data batch to be written

    Type(io_type), Intent( InOut ) :: io
    Integer,                                 Intent( In    ) :: io_file_handle,in_buffer
    Integer( Kind = MPI_OFFSET_KIND ),       Intent( In    ) :: next_rec
    Character( Len = 1 ), Dimension( :, : ), Intent( In    ) :: buffer

    Integer               :: status(1:MPI_STATUS_SIZE) 
    Integer :: file_handle, ierr
    Integer :: i

    ! Get the actual file handle or unit in use, rather than the one given to the outside world
    file_handle = io%known_files( io_file_handle )%file_handle

    If      ( io%method_write == IO_WRITE_UNSORTED_MPIIO  .or. io%method_write == IO_WRITE_SORTED_MPIIO  ) Then

      ! using MPI-I/O
      Call MPI_FILE_WRITE_AT( file_handle, next_rec, buffer, in_buffer, io%rec_type, status, ierr )

    Else If ( io%method_write == IO_WRITE_UNSORTED_DIRECT .or. io%method_write == IO_WRITE_SORTED_DIRECT ) Then

      ! using parallel direct access
      Write(Unit=file_handle, Fmt=io%forma, Rec=Int(next_rec,li)+Int(1,li)) ( buffer( : , i ) , i = 1 , in_buffer )

    End If

  End Subroutine io_write_batch

  Subroutine io_read_batch(io, io_file_handle, next_rec, in_buffer, buffer, error )

    ! Read a batch in a shared file
    ! Arguments are:
    ! IO_FILE_HANDLE: The file_handle
    ! NEXT_REC      : The starting record number
    ! IN_BUFFER     : the number of records in the batch
    ! BUFFER        : the data batch to be read

    Type(io_type), Intent( InOut ) :: io
    Integer,                                 Intent( In    ) :: io_file_handle,in_buffer
    Integer( Kind = MPI_OFFSET_KIND ),       Intent( In    ) :: next_rec
    Character( Len = 1 ), Dimension( :, : ), Intent(   Out ) :: buffer
    Integer                                , Intent(   Out ) :: error

    Integer :: file_handle
    Integer :: count
    Integer :: i, ierr

    Integer               :: status(1:MPI_STATUS_SIZE) 

    ! Get the actual file handle or unit in use, rather than the one given to the outside world
    file_handle = io%known_files( io_file_handle )%file_handle

    If      ( io%method_read == IO_READ_MPIIO ) Then

      ! using MPI-I/O
      Call MPI_FILE_READ_AT( file_handle, next_rec, buffer, in_buffer, io%rec_type, status, ierr )
      Call MPI_GET_COUNT( status, io%rec_type, count, ierr )
      error = Merge( 0, -1, count == in_buffer )

    Else If ( io%method_read == IO_READ_DIRECT ) Then

      ! using parallel direct access
      Read(Unit=file_handle, Fmt=io%forma, Rec=Int(next_rec,li)+Int(1,li), Iostat=error) &
        ( buffer( : , i ) , i = 1 , in_buffer )

    End If

  End Subroutine io_read_batch

  Subroutine io_nc_set_def(io, param, io_file_handle, title, n )

    Type(io_type), Intent( InOut ) :: io
    Type( netcdf_param ), Intent( InOut ) :: param
    Integer              , Intent( In    ) :: io_file_handle
    Character( Len = * ) , Intent( In    ) :: title
    Integer              , Intent( In    ) :: n

    Call netcdf_set_def( title, n, param, io%known_files( io_file_handle )%desc )

  End Subroutine io_nc_set_def

  Subroutine io_nc_put_var_rwp_0d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Real( Kind = wp )   , Intent( In    ) :: val
    Integer             , Intent( In    ) :: start
    Integer             , Intent( In    ) :: count

    Call netcdf_put_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_put_var_rwp_0d

  Subroutine io_nc_put_var_rwp_1d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                , Intent( In    ) :: what
    Integer                             , Intent( In    ) :: io_file_handle
    Real( Kind = wp )   , Dimension( : ), Intent( In    ) :: val
    Integer             , Dimension( : ), Intent( In    ) :: start
    Integer             , Dimension( : ), Intent( In    ) :: count

    Call netcdf_put_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_put_var_rwp_1d

  Subroutine io_nc_put_var_rwp_2d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                   , Intent( In    ) :: what
    Integer                                , Intent( In    ) :: io_file_handle
    Real( Kind = wp )   , Dimension( :, : ), Intent( In    ) :: val
    Integer             , Dimension( :    ), Intent( In    ) :: start
    Integer             , Dimension( :    ), Intent( In    ) :: count

    Call netcdf_put_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_put_var_rwp_2d

  Subroutine io_nc_put_var_int_0d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Integer             , Intent( In    ) :: val
    Integer             , Intent( In    ) :: start
    Integer             , Intent( In    ) :: count

    Call netcdf_put_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_put_var_int_0d

  Subroutine io_nc_put_var_int_1d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                , Intent( In    ) :: what
    Integer                             , Intent( In    ) :: io_file_handle
    Integer             , Dimension( : ), Intent( In    ) :: val
    Integer             , Dimension( : ), Intent( In    ) :: start
    Integer             , Dimension( : ), Intent( In    ) :: count

    Call netcdf_put_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_put_var_int_1d

  Subroutine io_nc_put_var_chr_2d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                   , Intent( In    ) :: what
    Integer                                , Intent( In    ) :: io_file_handle
    Character           , Dimension( :, : ), Intent( In    ) :: val
    Integer             , Dimension( :    ), Intent( In    ) :: start
    Integer             , Dimension( :    ), Intent( In    ) :: count

    Call netcdf_put_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_put_var_chr_2d

  Subroutine io_get_var(io, what, fh, start, count, x, y, z )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )             , Intent( In    ) :: what
    Integer                          , Intent( In    ) :: fh
    Integer   ,        Dimension( : ), Intent( In    ) :: start
    Integer   ,        Dimension( : ), Intent( In    ) :: count
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: x
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: y
    Real( Kind = wp ), Dimension( : ), Intent(   Out ) :: z

    Real( Kind = wp ), Dimension( :, : ), Allocatable :: buff

    Integer :: to_read
    Integer :: fail
    Integer :: i

    Character( Len = 256 ) :: message

    to_read = count( 2 )

    Allocate (buff( 1:3, 1:to_read ), Stat=fail)
    If (fail /= 0) Then
      Write( message, '(a)') 'allocation failure in io_get_var'
      Call error(0, message)
    End If

    Call io_nc_get_var(io,  what, fh, buff, start, count )

    Do i = 1, to_read
      x( i ) = buff( 1, i )
      y( i ) = buff( 2, i )
      z( i ) = buff( 3, i )
    End Do

    Deallocate (buff, Stat=fail)
    If (fail /= 0) Then
      Write( message, '(a)') 'deallocation failure in io_get_var'
      Call error(0, message)
    End If
  End Subroutine io_get_var

  Subroutine io_nc_get_var_rwp_0d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Real( Kind  = wp )  , Intent(   Out ) :: val
    Integer             , Intent( In    ) :: start
    Integer             , Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_rwp_0d

  Subroutine io_nc_get_var_rwp_1d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                , Intent( In    ) :: what
    Integer                             , Intent( In    ) :: io_file_handle
    Real( Kind = wp )   , Dimension( : ), Intent(   Out ) :: val
    Integer             , Dimension( : ), Intent( In    ) :: start
    Integer             , Dimension( : ), Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_rwp_1d

  Subroutine io_nc_get_var_rwp_2d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                   , Intent( In    ) :: what
    Integer                                , Intent( In    ) :: io_file_handle
    Real( Kind = wp )   , Dimension( :, : ), Intent(   Out ) :: val
    Integer             , Dimension( :    ), Intent( In    ) :: start
    Integer             , Dimension( :    ), Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_rwp_2d

  Subroutine io_nc_get_var_int_0d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Integer             , Intent(   Out ) :: val
    Integer             , Intent( In    ) :: start
    Integer             , Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_int_0d

  Subroutine io_nc_get_var_int_1d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                , Intent( In    ) :: what
    Integer                             , Intent( In    ) :: io_file_handle
    Integer             , Dimension( : ), Intent(   Out ) :: val
    Integer             , Dimension( : ), Intent( In    ) :: start
    Integer             , Dimension( : ), Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_int_1d

  Subroutine io_nc_get_var_chr_1d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                , Intent( In    ) :: what
    Integer                             , Intent( In    ) :: io_file_handle
    Character( Len = * ), Dimension( : ), Intent(   Out ) :: val
    Integer             , Dimension( : ), Intent( In    ) :: start
    Integer             , Dimension( : ), Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_chr_1d

  Subroutine io_nc_get_var_chr_2d(io, what, io_file_handle, val, start, count )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * )                   , Intent( In    ) :: what
    Integer                                , Intent( In    ) :: io_file_handle
    Character( Len = * ), Dimension( :, : ), Intent(   Out ) :: val
    Integer             , Dimension( :    ), Intent( In    ) :: start
    Integer             , Dimension( :    ), Intent( In    ) :: count

    Call netcdf_get_var( what, io%known_files( io_file_handle )%desc, val, start, count )

  End Subroutine io_nc_get_var_chr_2d

  Subroutine io_nc_get_att_int(io, what, io_file_handle, val )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Integer             , Intent(   Out ) :: val

    Call netcdf_get_att( what, io%known_files( io_file_handle )%desc, val )

  End Subroutine io_nc_get_att_int

  Subroutine io_nc_get_att_chr( io,what, io_file_handle, val )

    Type(io_type), Intent( InOut ) :: io
    Character( Len = * ), Intent( In    ) :: what
    Integer             , Intent( In    ) :: io_file_handle
    Character           , Intent(   Out ) :: val

    Call netcdf_get_att( what, io%known_files( io_file_handle )%desc, val )

  End Subroutine io_nc_get_att_chr

  Subroutine io_nc_set_real_precision( k, param, error )

    Integer, Intent( In    ) :: k
    Type( netcdf_param ), Intent( InOut ) :: param
    Integer, Intent(   Out ) :: error

    Call netcdf_set_real_precision( k, param, error )

  End Subroutine io_nc_set_real_precision

  Subroutine io_nc_get_real_precision( param, p, r, error )

    Type( netcdf_param ), Intent( In    ) :: param
    Integer, Intent(   Out ) :: p
    Integer, Intent(   Out ) :: r
    Integer, Intent(   Out ) :: error

    Call netcdf_get_real_precision( param, p, r, error )

  End Subroutine io_nc_get_real_precision

  Subroutine io_nc_get_file_real_precision(io, io_file_handle, p, r, error )

    Type(io_type), Intent( InOut ) :: io
    Integer, Intent( In    ) :: io_file_handle
    Integer, Intent(   Out ) :: p
    Integer, Intent(   Out ) :: r
    Integer, Intent(   Out ) :: error

    Call netcdf_get_file_real_precision( io%known_files( io_file_handle )%desc, p, r, error )

  End Subroutine io_nc_get_file_real_precision

  Subroutine split_io_comm( base_comm, n_io_procs_write, io_comm, io_gather_comm, do_io )

    ! Derive from the base communicator two further ones, IO_COMM and IO_GATHER_COMM.
    ! IO_COMM has N_IO_PROCS_WRITE and contains the processors that will perform I/O.
    ! Each member in IO_COMM is also processor zero in an instantiation of
    ! IO_GATHER_COMM.  This communicator is used to gather ( or, for reads, scatter )
    ! the data onto the I/O processor.  DO_IO indicates if this processor is a member
    ! of IO_COMM.
    ! Arguments are:
    ! BASE_COMM        : The base communicator
    ! N_IO_PROCS_WRITE : The number of processors that will perform I/O
    ! IO_COMM          : The communicator that will contain the I/O processors
    ! IO_GATHER_COMM   : The communicator for gathering onto/scattering from the I/O processors
    ! DO_IO            : true if this processor is will perform I/O

    Integer, Intent( In    ) :: base_comm
    Integer, Intent( In    ) :: n_io_procs_write
    Integer, Intent(   Out ) :: io_comm
    Integer, Intent(   Out ) :: io_gather_comm
    Logical, Intent(   Out ) :: do_io

    Integer :: actual_io_procs
    Integer :: n_base, rank_base
    Integer :: everyth_for_io
    Integer :: colour, key, ierr

    Call MPI_COMM_SIZE( base_comm,    n_base, ierr )
    Call MPI_COMM_RANK( base_comm, rank_base, ierr )

    ! After much though decided not to flag trying to use more I/O procs
    ! than there are in the communicator as an error - it avoids users
    ! potentially having to change their input for different processor counts,
    ! and not annoying users is the primary objective ( so that don't annoy me ).
    ! If more I/O procs asked for than available simply use all that are
    ! available
    actual_io_procs = Min( n_base, n_io_procs_write )

    ! Keep original order of ranks
    key = rank_base

    everyth_for_io = n_base / actual_io_procs

    ! Now create the comms.  Again have to be careful if the number of procs is
    ! not a multiple of the number of io procs.
    If ( Mod( n_base, actual_io_procs ) == 0 ) Then
      ! Exact multiple case

      ! Generate io_comm
      colour = Mod( rank_base, everyth_for_io )

      do_io = ( colour == 0 )

      If ( .not. do_io ) Then
        colour = MPI_UNDEFINED
      End If

      Call MPI_COMM_SPLIT( base_comm, colour, key, io_comm, ierr )

      ! And now the gather comm
      colour = rank_base / everyth_for_io
      Call MPI_COMM_SPLIT( base_comm, colour, key, io_gather_comm, ierr )

    Else
      ! Not an exact multiple.  Give out the extra processors as evenly as possible
      ! to each gather comm

      ! First the I/O comm
      If ( Mod( rank_base, everyth_for_io + 1 ) == 0 ) Then
        colour = 0
      Else
        colour = MPI_UNDEFINED
      End If
      do_io = colour == 0
      Call MPI_COMM_SPLIT( base_comm, colour, key, io_comm, ierr )

      ! And now the gather comm
      colour = rank_base / ( everyth_for_io + 1 )
      Call MPI_COMM_SPLIT( base_comm, colour, key, io_gather_comm, ierr )
    End If

  End Subroutine split_io_comm

  Subroutine free_io_comm( do_io, io_comm, io_gather_comm )

    ! Freeing IO_COMM and IO_GATHER_COMM.

    Logical                  :: do_io
    ! Note we can not set an intent on this - on procs NOT doing I/O it is
    ! intent Out, on procs doing I/O it is InOut
    Integer                  :: io_comm
    Integer, Intent( InOut ) :: io_gather_comm
    Integer :: ierr
    If( do_io ) Then
      Call MPI_COMM_FREE( io_comm, ierr )
    End If
    Call MPI_COMM_FREE( io_gather_comm, ierr )

  End Subroutine free_io_comm

  Subroutine batch_limits_set( local_global_indices, bottom_batch, top_batch, &
      local_bottom, local_top )

    ! Search for where to find the range of atoms BOTTOM_BATCH:TOP_BATCH in the local
    ! data.  Note we probably will not actually have the very atom BOTTOM_BATCH so
    ! have to be a little careful.  Similarly for TOP_BATCH.
    ! The arrays have assumed sorted, so can use a binary search.

    Integer, Dimension( : ), Intent( In    ) :: local_global_indices
    Integer,                 Intent( In    ) :: bottom_batch
    Integer,                 Intent( In    ) :: top_batch
    Integer,                 Intent(   Out ) :: local_bottom
    Integer,                 Intent(   Out ) :: local_top

    ! Set the bottom of the batch
    ! The search returns the first element above number searched for, or
    ! 1 more than the size of the array if the number is bigger than
    ! all elements in the array
    Call binary_search( bottom_batch, local_global_indices, local_bottom )

    ! Set the top of the batch
    Call binary_search( top_batch, local_global_indices, local_top )
    ! Want local_top to be the last atom held, not one above if we do not actually
    ! hold it.  Careful not to go out of bounds
    If ( local_top > Size( local_global_indices ) ) Then
      local_top = Size( local_global_indices )
    Else
      If ( local_global_indices( local_top ) /= top_batch ) Then
        local_top = local_top - 1
      End If
    End If

  End Subroutine batch_limits_set

  Subroutine binary_search( n, a, index )

    ! Binary search for N in the sorted array A.  The result is in INDEX.
    ! if the value is not found INDEX holds the value of the next highest entry
    ! in A, or one bigger than the size of the array if N is bigger than everything
    ! in A.

    Integer                , Intent( In    ) :: n
    Integer, Dimension( : ), Intent( In    ) :: a
    Integer                , Intent(   Out ) :: index

    Integer :: na
    Integer :: top, bottom, middle

    na = Size( a )

    ! zero sized arrays correction

    If ( na < 1 ) Then
      index = 1
      Return
    End If

    If ( a( na ) < n ) Then
      ! Catch over top of array
      index = na + 1
    Else If ( a( 1 ) >= n ) Then
      ! Catch below or at start of array
      index = 1
    Else
      ! Binary search
      bottom = 1
      top    = na
      Do
        If ( top - bottom == 1 ) Then
          Exit
        End If
        middle = ( top + bottom ) / 2
        If ( a( middle ) >= n ) Then
          top = middle
          If ( a( middle ) == n ) Then
            Exit
          End If
        Else
          bottom = middle
        End If
      End Do
      If ( a( bottom ) == n ) Then
        index = bottom
      Else
        index = top
      End If
    End If

  End Subroutine binary_search

  Function ok(io, flag, comm )
    Type(io_type), Intent( InOut ) :: io

    ! Check that FLAG is TRUE across all of the communicator COMM
    ! Arguments are:
    ! FLAG: The flag
    ! COMM: The communicator

    Logical                  :: ok

    Logical, Intent( In    ) :: flag
    Integer, Intent( In    ) :: comm

    Logical :: g_flag
    Integer :: ierr

    If ( io%global_error_check ) Then

      Call MPI_ALLREDUCE( flag, g_flag, 1, MPI_LOGICAL, MPI_LAND, comm, ierr )
      ok = g_flag

    Else

      ok = flag

    End If

  End Function ok

  Function get_tot_atoms( atoms, comm )

    ! Get the maximum number of atoms held by processors in communicator COMM
    ! Arguments are:
    ! ATOMS: The local number of atoms
    ! COMM : The communicator

    Integer                  :: get_tot_atoms

    Integer, Intent( In    ) :: atoms
    Integer, Intent( In    ) :: comm

    Integer :: g_atoms, ierr

    Call MPI_ALLREDUCE( atoms, g_atoms, 1, MPI_INTEGER, MPI_SUM, comm, ierr )

    get_tot_atoms = g_atoms

  End Function get_tot_atoms

  Subroutine get_file_handle( io,file_handle )
    Type(io_type), Intent( InOut ) :: io

    Integer, Intent(   Out ) :: file_handle

    Integer :: i
    logical :: open

    Do i = LOW_HANDLE, HIGH_HANDLE
      If ( io%known_files( i )%method == UNUSED ) Then
        Inquire( unit = i, opened = open )
        If ( .not. open ) Then
          file_handle = i
          Return
        End If
      End If
    End Do

    ! Place a known elephant
    file_handle = LOW_HANDLE - 1

  End Subroutine get_file_handle

  Subroutine io_nc_compiled()

    Call netcdf_compiled()

  End Subroutine io_nc_compiled

End Module io

