Module kim
!> Module providing support for the KIM v2.0 library
!>
!> Copyright - Daresbury Laboratory
!>
!> Author  - J. Madge July 2018
!> Contrib - Y. Afshar Dec 2019
!> Contrib - A.M. Elena Dec 2019
!>
!> Based on KIM v2.0 simulator example in Fortran and the KIM v1 implementation
!> by R. S. Elliot
  Use, Intrinsic :: iso_c_binding,  Only: c_double,&
                                          c_int,&
                                          c_char,& 
                                          c_ptr,&
                                          c_null_ptr,& 
                                          c_funloc,& 
                                          c_loc,& 
                                          c_f_pointer,& 
                                          c_signed_char
  Use comms,                        Only: Export_tag,&
                                          comms_type,&
                                          girecv,&
                                          gsend,&
                                          gwait
  Use errors_warnings,              Only: error,&
                                          info,&
                                          warning
#ifdef KIM
  Use kim_log_module,               Only: kim_pop_default_verbosity,&
                                          kim_push_default_verbosity
  Use kim_log_verbosity_module,     Only: KIM_LOG_VERBOSITY_DEBUG,&
                                          KIM_LOG_VERBOSITY_ERROR
  Use kim_simulator_headers_module, Only: &
                                          KIM_COLLECTION_ITEM_TYPE_PORTABLE_MODEL, &
                                          KIM_COLLECTION_ITEM_TYPE_SIMULATOR_MODEL, &
                                          KIM_DATA_TYPE_INTEGER, &
                                          kim_cache_list_of_item_metadata_files, &
                                          kim_charge_unit_e, kim_collection_item_type_type, &
                                          kim_collections_create, kim_collections_destroy, &
                                          kim_collections_handle_type, kim_compute, &
                                          kim_compute_argument_name_coordinates, &
                                          kim_compute_argument_name_number_of_particles, &
                                          kim_compute_argument_name_partial_energy, &
                                          kim_compute_argument_name_partial_forces, &
                                          kim_compute_argument_name_partial_virial, &
                                          kim_compute_argument_name_particle_contributing, &
                                          kim_compute_argument_name_particle_species_codes, &
                                          kim_compute_arguments_create, &
                                          kim_compute_arguments_destroy, &
                                          kim_compute_arguments_handle_type, &
                                          kim_compute_callback_name_get_neighbor_list, &
                                          kim_data_type_type, kim_energy_unit_amu_a2_per_ps2, &
                                          kim_from_string, kim_get_influence_distance, &
                                          kim_get_item_metadata_file_values, kim_get_item_type, &
                                          kim_get_neighbor_list_values, &
                                          kim_get_number_of_neighbor_lists, &
                                          kim_get_number_of_parameters, &
                                          kim_get_parameter_metadata, &
                                          kim_get_species_support_and_code, &
                                          kim_language_name_fortran, kim_length_unit_a, &
                                          kim_model_create, kim_model_destroy, &
                                          kim_model_handle_type, kim_numbering_one_based, &
                                          kim_set_argument_pointer, kim_set_callback_pointer, &
                                          kim_species_name_type, kim_temperature_unit_k, &
                                          kim_time_unit_ps, kim_to_string, Operator(.eq.)
#endif
  Use kinds,                        Only: wi,&
                                          wp
  Use numerics,                     Only: local_index
  Use particle,                     Only: corepart

  Implicit None

  Private

  !> Flag indicating whether the program has been compiled with the KIM preprocessor flag
#ifdef KIM
  Logical, Parameter, Public :: COMPILED_WITH_KIM = .true.
#else
  Logical, Parameter, Public :: COMPILED_WITH_KIM = .false.
#endif

  !> Name of this source file
  Character(Len=*), Parameter :: FILE_NAME = "kim.F90"

  !> Number of scalars per particle to be sent during communication.
  !> In this case force_x, force_y, > force_z and global_id
  Integer(Kind=wi), Parameter :: span = 4

  !> Neighbour list type containing KIM neighbour list data
  Type :: kim_neighbour_list_type
    Private

    !> Potential cut off
    Real(Kind=c_double)              :: cutoff
    !> Total number of particles
    Integer(Kind=c_int)              :: n_part
    !> Number of neighbours for each particle
    Integer(Kind=c_int), Allocatable :: n_neigh(:)
    !> Neighbour list
    !>
    !> `list(n,part)` is the number of the `n`th neighbour of particle `part`
    Integer(Kind=c_int), Allocatable :: neigh_list(:, :)
  Contains
    Private

    Procedure :: init => kim_neighbour_list_type_init
    Final     :: kim_neighbour_list_type_cleanup
  End Type kim_neighbour_list_type

  !> Neighbour list type which the KIM api will pass to the routine
  !> get_neigh. This type must be C interoperable and so cannot have allocatable
  !> components or type bound procedures.
  Type, Bind(c) :: kim_neighbour_list_pointer_type
    Private

    !> Maximum number of atoms per node; the dimension of n_neigh and second
    !> dimsenion of neigh_list
    Integer(c_int) :: max_atoms
    !> Maximum number of neighbours per particle; the first dimension of neigh_list
    Integer(c_int) :: max_neigh
    !> Pointer to kim_neighbour_list_type cutoff
    Type(c_ptr)    :: cutoff
    !> Pointer to kim_neighbour_list_type n_part
    Type(c_ptr)    :: n_part
    !> Pointer to kim_neighbour_list_type n_neigh
    Type(c_ptr)    :: n_neigh
    !> Pointer to kim_neighbour_list_type neigh_list
    Type(c_ptr)    :: neigh_list
  End Type kim_neighbour_list_pointer_type

  !> Type containing data necessary for the KIM message passing in the routine
  !> kim_share_halo_forces
  Type :: kim_comms_type
    Private

    !> Buffer for message passing
    Real(Kind=wp), Allocatable :: buffer(:)
    !> Half buffer size, used to divide the portion for sending and the portion
    !> for receiving
    Integer(Kind=wi)           :: recv_start
    !> Number of atoms to receive in each direction
    Integer(Kind=wi)           :: n_recv(6)
    !> ids of first and last atoms to send respectively in each direction
    Integer(Kind=wi)           :: first(6), last(6)

  Contains
    Private

    Procedure         :: init => kim_comms_type_init
    Procedure, Public :: set => kim_comms_type_set
    Final             :: kim_comms_type_cleanup
  End Type kim_comms_type

  !> Type containing KIM data
  Type, Public :: kim_type
    Private

    !> Flag indicating whether KIM is in use
    Logical, Public                                    :: active = .false.
    !> Flag indicating if any neighour lists do not have the padding hint and so
    !> require neighbours of non-contributing (padding) particles
    Logical, Public                                    :: padding_neighbours_required = .false.
    !> Name of the KIM model requested
    Character(Len=:, Kind=c_char), Allocatable, Public :: model_name
    !> Largest model cutoff, used to determining cell list and domain size
    Real(Kind=wp), Public                              :: cutoff = 0.0_wp
    !> Model influence distance
    Real(Kind=c_double), Public                        :: influence_distance
    !> Neighbour list type
    Type(kim_neighbour_list_type), Allocatable         :: neigh(:)
    !> Neighbour lsit pointer type
    Type(kim_neighbour_list_pointer_type), Allocatable :: neigh_pointer(:)
    !> Comms type
    Type(kim_comms_type), Public                       :: kcomms

#ifdef KIM
    !> Model handle
    Type(kim_model_handle_type)             :: model_handle
    !> Compute arguments handle
    Type(kim_compute_arguments_handle_type) :: compute_arguments_handle
    !> Model type
    Type(kim_collection_item_type_type)     :: model_type
    !> Number of particles
    Integer(Kind=c_int)                     :: n_particles
    !> Number of model neighbour lists
    Integer(Kind=c_int)                     :: n_lists
    !> Neighbour list padding hints
    Integer(Kind=c_int), Allocatable        :: hints_padding(:)
    ! Species
    !> List of species codes
    Integer(Kind=c_int), Allocatable        :: species_code(:)
    !> List of unique species codes
    Integer(Kind=c_int), Allocatable        :: unique_species_code(:)
    !> Contributing list
    !>
    !> - 1 if particle contributes
    !> - 0 otherwise
    Integer(Kind=c_int), Allocatable        :: contributing(:)
    !> Coordinates
    Real(Kind=c_double), Allocatable        :: coords(:, :)
    !> KIM energy
    Real(Kind=c_double)                     :: energy
    !> KIM forces
    Real(Kind=c_double), Allocatable        :: forces(:, :)
    !> KIM virial
    Real(Kind=c_double)                     :: virial(6)
#endif
  Contains
    Private

    Procedure :: init => kim_type_init
    Final     :: kim_type_cleanup
  End Type kim_type

  Public :: kim_setup
  Public :: kim_cutoff
  Public :: kim_interactions
  Public :: kim_energy_and_forces
  Public :: get_neigh
  Public :: kim_citations

Contains

  !> Interface to DL_POLY error routine
  Subroutine kim_error(message, line)
    !> Error message for printing in the standard output file
    Character(Kind=c_char, Len=*), Intent(In   ) :: message
    Integer(Kind=wi),              Intent(In   ) :: line

    Character(Len=3)              :: line_str
    Character(Len=:), Allocatable :: error_message

!> Line number in the source file where error happened

    Write (line_str, '(i3)') line
    error_message = 'KIM error: '//Trim(message)//' line: '// &
                    Trim(Adjustl(line_str))//' file: '//FILE_NAME
    Call error(0, error_message, .true.)
  End Subroutine kim_error

  !> Interface to DL_POLY warning routine
  Subroutine kim_warning(message, line)
    !> Warning message for printing in the standard output file
    Character(Kind=c_char, Len=*), Intent(In   ) :: message
    Integer(Kind=wi),              Intent(In   ) :: line

    Character(Len=3)              :: line_str
    Character(Len=:), Allocatable :: warning_message

!> Line number in the source file

    Write (line_str, '(i3)') line
    warning_message = 'KIM warning: '//Trim(message)//' line: '// &
                      Trim(Adjustl(line_str))//' file: '//FILE_NAME
    Call warning(warning_message, .true.)
  End Subroutine kim_warning

  !> Initialise KIM types
  Subroutine kim_setup(kim_data, mxatms, mxatdm, megatm, max_list, mxbfxp, mxnode)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data
    !> Extent of particle arrays including halo particles (coordinates, forces, etc.)
    Integer(Kind=wi),       Intent(In   ) :: mxatms
    !> Extent of particle arrays not including halo particles (neighbours, contraints, etc.)
    Integer(Kind=wi),       Intent(In   ) :: mxatdm
    !> Total number of particles on all domains/nodes
    Integer(Kind=wi),       Intent(In   ) :: megatm
    !> Extent of neighbour list arrays
    Integer(Kind=wi),       Intent(In   ) :: max_list
    !> Extent of force arrays
    Integer(Kind=wi),       Intent(In   ) :: mxbfxp
    !> Number of nodes
    Integer(Kind=wi),       Intent(In   ) :: mxnode

    Real(Kind=c_double), Allocatable :: cutoffs(:)
    Integer(Kind=c_int)              :: requested_units_accepted
    Integer(Kind=wi)                 :: list_index, max_atoms
    Integer(Kind=c_int)              :: kerror
    Integer(Kind=c_int)              :: n_parameters
    Integer(Kind=c_int)              :: extent
    Integer(Kind=c_int)              :: parameter_index
    Integer(Kind=wi)                 :: max_len, i
    Integer(Kind=wi)                 :: fail(4)
    Character(Kind=c_char, Len=256)  :: parameter_name
    Character(Kind=c_char, Len=1024) :: parameter_description
    Character(Kind=c_char, Len=256)  :: parameter_string
    Character(Len=68)                :: fmt, banner(7)
    Character(Len=256)               :: message

#ifdef KIM
    Type(kim_data_type_type) :: kim_data_type
#endif

    ! Early return if a KIM model is not being used
    If (.not. kim_data%active) Then
      Return
    End If

#ifdef KIM
    If (kim_data%model_type .eq. &
        KIM_COLLECTION_ITEM_TYPE_PORTABLE_MODEL) Then

      ! Allocate KIM data type arrays
      Call kim_data%init(mxatms)

      ! Initialise comms type
      Call kim_data%kcomms%init(mxbfxp * span, mxnode)

      ! Print KIM model parameters
      Call kim_get_number_of_parameters(kim_data%model_handle, n_parameters)

      Write (banner(1), '(a)') ""
      Write (banner(2), '(a)') '//'//Repeat("=", 30)//' KIM '// &
        Repeat("=", 31)
      Write (banner(3), '(a)') '||'
      Write (banner(5), '(a)') '||'
      Write (banner(6), '(a)') '||'
      Write (banner(7), '(a)') '\\'//Repeat("=", 66)
      If (n_parameters .gt. 0_c_int) Then
        Write (banner(4), '(a, i0, a)') '|| This model has ', n_parameters, &
          ' mutable parameters.'
        Call info(banner, 6, .true.)

        max_len = 18
        Do parameter_index = 1_c_int, n_parameters
          Call kim_get_parameter_metadata( &
            kim_data%model_handle, &
            parameter_index, &
            kim_data_type, &
            extent, &
            parameter_name, &
            parameter_description, &
            kerror)
          If (kerror /= 0_c_int) Then
            Call kim_error('kim_setup, kim_get_parameter_metadata '// &
                           'returned an error.', __LINE__)
          End If
          max_len = Max(max_len, Len_trim(parameter_name))
        End Do

        Write (message, '(a, a, a)') '|| No.     |  Parameter name  ', &
          Repeat(' ', max_len - 18), '|  Data type  |  Extent'
        Call info(message, .true.)

        Write (message, '(a)') '||'//Repeat('=', 66)
        Call info(message, .true.)

        Do parameter_index = 1_c_int, n_parameters
          Call kim_get_parameter_metadata( &
            kim_data%model_handle, &
            parameter_index, &
            kim_data_type, &
            extent, &
            parameter_name, &
            parameter_description, &
            kerror)

          Call kim_to_string(kim_data_type, parameter_string)

          Write (fmt, '(i0)') parameter_index

          If (kim_data_type .eq. KIM_DATA_TYPE_INTEGER) Then
            Write (message, '(a, i0, 7a, i0)') '|| ', parameter_index, &
              Repeat(' ', 8 - Len_trim(fmt)), '| ', Trim(parameter_name), &
              Repeat(' ', Max(17, Len_trim(parameter_name)) - &
                     Len_trim(parameter_name)), &
              '|  "', Trim(parameter_string), '"  | ', extent
          Else
            Write (message, '(a, i0, 7a, i0)') '|| ', parameter_index, &
              Repeat(' ', 8 - Len_trim(fmt)), '| ', Trim(parameter_name), &
              Repeat(' ', Max(17, Len_trim(parameter_name)) - &
                     Len_trim(parameter_name)), &
              '|  "', Trim(parameter_string), '"   | ', extent
          End If
          Call info(message, .true.)
        End Do
        Call info(banner(7), .true.)
      Else
        Write (banner(4), '(a)') '|| This model has No mutable parameters.'

        Call info(banner, 5, .true.)
        Call info(banner(7), .true.)
      End If

      ! Create compute_arguments object
      Call kim_compute_arguments_create( &
        kim_data%model_handle, &
        kim_data%compute_arguments_handle, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_compute_arguments_create, '// &
                       'failed to allocate a new ComputeArguments object.', __LINE__)
      End If

      ! Set number of particles
      kim_data%n_particles = Int(megatm, c_int)

      fail = 0

      ! Allocate neighbour list
      Allocate (kim_data%neigh(kim_data%n_lists), Stat=fail(1))
      ! Allocate neighbour list pointer
      Allocate (kim_data%neigh_pointer(kim_data%n_lists), Stat=fail(2))
      ! Allocate hint arrays
      Allocate (kim_data%hints_padding(kim_data%n_lists), Stat=fail(3))
      ! Get neighbour list cutoffs
      Allocate (cutoffs(kim_data%n_lists), Stat=fail(4))
      If (Any(fail /= 0)) Then
        Call kim_error('kim_setup, allocation failure.', __LINE__)
      End If

      Call kim_get_neighbor_list_values( &
        kim_data%model_handle, &
        cutoffs, &
        kim_data%hints_padding, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_get_neighbor_list_values, '// &
                       'wrong number of neighbour lists.', __LINE__)
      End If

      Do list_index = 1, Int(kim_data%n_lists, wi)
        kim_data%neigh(list_index)%cutoff = cutoffs(list_index)
      End Do

      fail(1) = 0
      Deallocate (cutoffs, Stat=fail(1))
      If (fail(1) /= 0) Then
        Call kim_error('kim_setup, cutoffs deallocation failure.', __LINE__)
      End If

      ! Initialise neighbour list and neighbour list pointer types
      Do list_index = 1, Int(kim_data%n_lists, wi)
        If (kim_data%hints_padding(list_index) == 1_c_int) Then
          max_atoms = mxatdm
        Else
          max_atoms = mxatms
        End If

        ! Initialise neighbour list type
        Call kim_data%neigh(list_index)%init(max_atoms, max_list)

        ! Initialise neighbour list pointer type
        Call kim_neighbour_list_pointer_type_init( &
          kim_data%neigh_pointer(list_index), &
          kim_data%neigh(list_index), &
          max_atoms, &
          max_list)
      End Do

      ! Allocate KIM pointers
      ! Number of particles
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_number_of_particles, &
        kim_data%n_particles, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_argument_pointer, '// &
                       'kim_compute_arguments_set_arugment_pointer, '// &
                       'failed to set arugment pointer for number of particles.', __LINE__)
      End If

      ! Species codes
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_particle_species_codes, &
        kim_data%species_code, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_argument_pointer, '// &
                       'kim_compute_arguments_set_arugment_pointer, '// &
                       'failed to set arugment pointer for species codes.', __LINE__)
      End If

      ! Contributing status
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_particle_contributing, &
        kim_data%contributing, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_argument_pointer, '// &
                       'kim_compute_arguments_set_arugment_pointer '// &
                       'failed to set arugment pointer for contributing status.', __LINE__)
      End If

      ! Coordinates
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_coordinates, &
        kim_data%coords, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_argument_pointer, '// &
                       'kim_compute_arguments_set_arugment_pointer '// &
                       'failed to set arugment pointer for coordinates.', __LINE__)
      End If

      ! Energy
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_partial_energy, &
        kim_data%energy, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_argument_pointer, '// &
                       'kim_compute_arguments_set_arugment_pointer '// &
                       'failed to set arugment pointer for energy.', __LINE__)
      End If

      ! Forces
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_partial_forces, &
        kim_data%forces, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_argument_pointer, '// &
                       'kim_compute_arguments_set_arugment_pointer '// &
                       'failed to set arugment pointer for forces.', __LINE__)
      End If

      ! Virials
      Call kim_set_argument_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_argument_name_partial_virial, &
        kim_data%virial, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_warning('kim_setup, kim_set_argument_pointer, '// &
                         'The selected KIM model does not compute '// &
                         'virials, stress and pressure will be incorrect.', __LINE__)
        ! Set KIM virials to 0 so they will not contribute to the total
        kim_data%virial = 0.0_wp
      End If

      ! Set KIM pointer to neighbour list routine and type
      Call kim_set_callback_pointer( &
        kim_data%compute_arguments_handle, &
        kim_compute_callback_name_get_neighbor_list, &
        kim_language_name_fortran, &
        c_funloc(get_neigh), &
        c_loc(kim_data%neigh_pointer), &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_setup, kim_set_callback_pointer '// &
                       'returned an error.', __LINE__)
      End If
    Else If (kim_data%model_type .eq. &
             KIM_COLLECTION_ITEM_TYPE_SIMULATOR_MODEL) Then
      Call kim_error('kim_setup, currently DL_POLY does not support '// &
                     'KIM Simulator model.', __LINE__)
    Else
      Call kim_error('kim_setup, unknown model type.', __LINE__)
    End If
#endif
  End Subroutine kim_setup

  !> Function to return the cutoff of the kim model 'model_name'.
  !>
  !> This is necessary as the cutoff is required for ensuring domains and cells
  !> are of an appropriate size. Perhaps with a reworking of the way the input
  !> files are read could remedy this clumsy implementation.
  Subroutine kim_cutoff(kim_data)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data

    Real(Kind=c_double), Allocatable :: cutoffs(:)
    Integer(Kind=c_int), Allocatable :: hints_padding(:)
    Integer(Kind=c_int)              :: requested_units_accepted
    Integer(Kind=c_int)              :: kerror
    Integer(Kind=wi)                 :: fail(2)

#ifdef KIM
    Type(kim_collections_handle_type) :: kim_coll

    ! This is the debugging flag to be set manually.
    ! @todo
    ! It can be set with the CMAKE flag as CMAKE_BUILD_TYPE=Debug
    Logical :: debug = .false.
#endif

    If (.not. COMPILED_WITH_KIM) Then
      Call error(0, 'KIM directive found in FIELD, but the program is '// &
                 'not built with openKIM support', .true.)
    End If

#ifdef KIM
    ! Set the the KIM API verbosity mode.
    Call kim_pop_default_verbosity()
    If (debug) Then
      ! The standard debug verbosity. It should not be used for runtime
      ! since there would be too much of information to print.
      Call kim_push_default_verbosity(KIM_LOG_VERBOSITY_DEBUG)
    Else
      ! The standard error verbosity. It only reports when the execution of
      ! some task could not be completed.
      Call kim_push_default_verbosity(KIM_LOG_VERBOSITY_ERROR)
    End If

    ! Create a new KIM API collection object
    Call kim_collections_create(kim_coll, kerror)
    If (kerror /= 0_c_int) Then
      Call kim_error('kim_cutoff, kim_collections_create, unable to '// &
                     'access KIM Collections to find the requested Model.', __LINE__)
    End If

    ! Get the collection item with the requested name. It should be
    ! a KIM portableModel type otherwise it fails
    Call kim_get_item_type( &
      kim_coll, &
      Trim(kim_data%model_name), &
      kim_data%model_type, &
      kerror)
    If (kerror /= 0_c_int) Then
      Call kim_error('kim_cutoff, kim_get_item_type, KIM Model name = `'// &
                     Trim(kim_data%model_name)//'` not found.', __LINE__)
    End If

    ! Destroy a previously created collection.
    Call kim_collections_destroy(kim_coll)

    ! Check if the item type is a KIM portable model
    If (kim_data%model_type .eq. &
        KIM_COLLECTION_ITEM_TYPE_PORTABLE_MODEL) Then
      ! Create KIM object
      ! kim_numbering_one_based - numbering from 1 (fortran style arrays)
      !
      ! DL_POLY internal units:
      ! - length Angstrom
      ! - energy 10J/mol (equivalent KIM energy unit selected here)
      ! - charge electrons
      ! - temperature K
      ! - time ps
      Call kim_model_create( &
        kim_numbering_one_based, &
        kim_length_unit_a, &
        kim_energy_unit_amu_a2_per_ps2, &
        kim_charge_unit_e, &
        kim_temperature_unit_k, &
        kim_time_unit_ps, &
        Trim(kim_data%model_name), &
        requested_units_accepted, &
        kim_data%model_handle, &
        kerror)
      If (requested_units_accepted == 0_c_int) Then
        Call kim_error('kim_cutoff, kim_model_create, the selected '// &
                       'KIM model does not support the DL_POLY internal units.', __LINE__)
      End If
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_cutoff, kim_model_create, Unable to load '// &
                       'KIM Model.', __LINE__)
      End If

      ! Get number of neighbour lists
      Call kim_get_number_of_neighbor_lists(kim_data%model_handle, &
                                            kim_data%n_lists)

      ! Get neighbour list cutoffs and hints

      fail = 0

      ! Allocate cutoffs and hint arrays
      Allocate (hints_padding(kim_data%n_lists), Stat=fail(1))
      Allocate (cutoffs(kim_data%n_lists), Stat=fail(2))
      If (Any(fail /= 0)) Then
        Call kim_error('kim_cutoff, cutoffs, and/or hints_padding '// &
                       'allocation failure.', __LINE__)
      End If

      ! Get neighbour list cutoffs and record maximum
      Call kim_get_neighbor_list_values( &
        kim_data%model_handle, &
        cutoffs, &
        hints_padding, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_cutoff, kim_get_neighbor_list_values, '// &
                       'Wrong number of neighbour lists.', __LINE__)
      End If

      kim_data%cutoff = Real(Maxval(cutoffs), wp)

      ! Check if the padding hint is defined
      If (Any(hints_padding /= 1_c_int)) Then
        Call kim_warning('The selected KIM model requires neighbours '// &
                         'of non-contributing particles. This may significantly '// &
                         'affect the performance.', __LINE__)

        kim_data%padding_neighbours_required = .true.
      End If

      ! Get influence distance
      Call kim_get_influence_distance(kim_data%model_handle, &
                                      kim_data%influence_distance)

      fail = 0
      ! Deallocate temporary variables
      Deallocate (cutoffs, Stat=fail(1))
      Deallocate (hints_padding, Stat=fail(2))
      If (Any(fail /= 0)) Then
        Call kim_error('kim_cutoff, cutoffs, and/or hints_padding '// &
                       'deallocation failure.', __LINE__)
      End If
    Else If (kim_data%model_type .eq. &
             KIM_COLLECTION_ITEM_TYPE_SIMULATOR_MODEL) Then
      Call kim_error('kim_cutoff, currently DL_POLY does not support '// &
                     'KIM Simulator model.', __LINE__)
    Else
      Call kim_error('kim_cutoff, unknown model type.', __LINE__)
    End If
#endif
  End Subroutine kim_cutoff

  !> Function that defines a list of N species, and defines a mapping between
  !> atom types in DL_POLY to the available species in the OpenKIM IM
  !>
  !> This function must be preceded by the DL_POLY defination of the number
  !> of atom types and a 'kim_cutoff' call through using the 'kim_init'
  !> keyword in the FIELD file
  Subroutine kim_interactions(kim_data, ntype_atom, unique_atom, &
                              kim_ntype_atom, kim_species_name)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data
    !> Number of atom types
    Integer(Kind=wi), Intent(In)          :: ntype_atom
    !> Unique atom labels
    Character(Len=8), Intent(In)          :: unique_atom(:)
    !> Number of atom types in the KIM interaction
    Integer(Kind=wi), Intent(In)          :: kim_ntype_atom
    !> Unique atom labels in the KIM interaction
    Character(Len=8), Intent(In)          :: kim_species_name(:)

    Integer(Kind=c_int) :: species_is_supported
    Integer(Kind=c_int) :: species_code
    Integer(Kind=c_int) :: kerror
    Integer(Kind=wi)    :: iatom, jatom
    Integer(Kind=wi)    :: fail

#ifdef KIM
    Type(kim_species_name_type) :: species_name
#endif

    fail = 0

#ifdef KIM
    ! Allocate KIM unique_species_code array
    Allocate (kim_data%unique_species_code(ntype_atom), Stat=fail)
    If (fail /= 0) Then
      Call kim_error('kim_interactions, unique_species_code.'// &
                     'allocation failure.', __LINE__)
    End If

    ! Initialise the species code
    kim_data%unique_species_code = -2_c_int

    ! Loop through the kim interactions atom types
    Do iatom = 1, kim_ntype_atom
      Call kim_from_string(Trim(kim_species_name(iatom)), species_name)

      ! Check model supports the requested species
      Call kim_get_species_support_and_code( &
        kim_data%model_handle, &
        species_name, &
        species_is_supported, &
        species_code, &
        kerror)
      If ((kerror /= 0_c_int) .or. (species_is_supported /= 1_c_int)) Then
        Call kim_error('kim_interactions, '// &
                       'The species `'//Trim(kim_species_name(iatom))// &
                       '` is not supported by KIM API.', __LINE__)
      End If

      Do jatom = 1, ntype_atom
        If (kim_species_name(iatom) == unique_atom(jatom)) Then
          ! Assign the unique code
          kim_data%unique_species_code(jatom) = species_code
          ! Break the inside loop
          Exit
        End If
      End Do
    End Do

    ! Extra check
    If (Count(kim_data%unique_species_code /= -2_c_int) /= &
        Int(kim_ntype_atom, c_int)) Then
      Call kim_error('kim_interactions, There is an undefined '// &
                     'species type passed to `kim_interactions`.', __LINE__)
    End If

#endif
  End Subroutine kim_interactions

  !> Provide KIM citations to reference publication
  Subroutine kim_citations(kim_data, comm)
    !> KIM data type
    Type(kim_type), Target, Intent(In) :: kim_data
    !> Comms data
    Type(comms_type),       Intent(In) :: comm

#ifdef KIM
    Type(kim_collections_handle_type) :: kim_coll
#endif

    Integer                           :: unit_no = -2
    Integer(Kind=c_int)               :: kerror
    Integer(Kind=c_int)               :: extent
    Integer(Kind=c_int)               :: index
    Integer(Kind=c_signed_char)       :: cite_file_raw_data(10000)
    Character(Kind=c_char, Len=2048)  :: cite_file_name
    Character(Kind=c_char, Len=10000) :: cite_file_string

#ifdef KIM
    ! Create a KIM API collection object
    Call kim_collections_create(kim_coll, kerror)
    If (kerror /= 0_c_int) Then
      Call kim_error('kim_citations, kim_collections_create, '// &
                     'unable to access KIM Collections to find the Model.', __LINE__)
    End If

    ! Check if the item type is a KIM portable model
    If (kim_data%model_type .eq. &
        KIM_COLLECTION_ITEM_TYPE_PORTABLE_MODEL) Then
      ! Cache a list of an item's metadata files.
      Call kim_cache_list_of_item_metadata_files( &
        kim_coll, &
        KIM_COLLECTION_ITEM_TYPE_PORTABLE_MODEL, &
        Trim(kim_data%model_name), &
        extent, &
        kerror)
      If (kerror /= 0_c_int) Then
        Call kim_error('kim_citations, '// &
                       'kim_cache_list_of_item_metadata_files, failed to cache '// &
                       'a list of metadata files.', __LINE__)
      End If
    Else If (kim_data%model_type .eq. &
             KIM_COLLECTION_ITEM_TYPE_SIMULATOR_MODEL) Then
      Call kim_error('kim_citations, currently DL_POLY does not '// &
                     'support KIM Simulator model.', __LINE__)
    Else
      Call kim_error('kim_citations, unknown model type.', __LINE__)
    End If

    If (comm%idnode == 0) Then
      Open (Newunit=unit_no, File='kim.cite', Status='replace')

      Do index = 1_c_int, extent
        ! Get the name and content of the metadata files.
        Call kim_get_item_metadata_file_values( &
          kim_coll, &
          index, &
          cite_file_name, &
          cite_file_raw_data, &
          cite_file_string, &
          kerror)
        If (kerror /= 0_c_int) Then
          Call kim_error('kim_citations, '// &
                         'kim_get_item_metadata_file_values, failed to get '// &
                         'the name and content of metadata files.', __LINE__)
        End If

        If (cite_file_name(1:7) == 'kimcite') Then
          Write (unit_no, Fmt='(A)') Trim(cite_file_string)
          Write (unit_no, Fmt='(A)')
        End If
      End Do

      Close (Unit=unit_no)
    End If

    ! Destroy a previously created collection.
    Call kim_collections_destroy(kim_coll)
#endif
  End Subroutine kim_citations

  !> Compute KIM energy and forces
  Subroutine kim_energy_and_forces(kim_data, natms, nlast, parts, neigh_list, &
                                   map, lsite, ltype, lsi, lsa, ltg, site_name, energy_kim, virial_kim, stress, &
                                   comm)
    !> KIM data type
    Type(kim_type), Target,       Intent(InOut) :: kim_data
    !> Number of particles in this domain (excluding the halo)
    Integer(Kind=wi),             Intent(In   ) :: natms
    !> Number of particles in this domain and it's halo
    Integer(Kind=wi),             Intent(In   ) :: nlast
    !> Particles
    Type(corepart), Dimension(:), Intent(InOut) :: parts
    !> DL_POLY neighbour list
    !>
    !> - list(0,n) is the number of neighbours of particle n
    !> - list(m,n) is the 'm'th neighbour of particle n
    Integer(Kind=wi),             Intent(In   ) :: neigh_list(-3:, 1:)
    !> Map of neighbouring domain ids
    Integer(Kind=wi),             Intent(In   ) :: map(1:26)
    !> Site local index to global index
    Integer(Kind=wi),             Intent(In   ) :: lsite(:)
    !> Type of each local atom
    Integer(Kind=wi),             Intent(In   ) :: ltype(:)
    !> Some sort of local to global arrays
    Integer(Kind=wi),             Intent(In   ) :: lsi(:)
    !> Some sort of local to global arrays
    Integer(Kind=wi),             Intent(In   ) :: lsa(:)
    !> Local to global id
    Integer(Kind=wi),             Intent(In   ) :: ltg(:)
    !> Names of each atom
    Character(Len=8),             Intent(In   ) :: site_name(:)
    !> KIM model energy
    Real(Kind=wp),                Intent(Out  ) :: energy_kim
    !> KIM model virial
    Real(Kind=wp),                Intent(Out  ) :: virial_kim
    !> Virial
    Real(Kind=wp),                Intent(InOut) :: stress(1:9)
    !> Comms data
    Type(comms_type),             Intent(InOut) :: comm

    Integer(Kind=wi)    :: atom, list_index
    Integer(Kind=c_int) :: species_code
    Integer(Kind=c_int) :: kerror

#ifdef KIM
    ! Set number of particles
    kim_data%n_particles = Int(nlast, c_int)

    ! Copy coordinates to KIM coordinates array
    kim_data%coords(:, nlast + 1:) = 0.0_c_double
    Do atom = 1, nlast
      kim_data%coords(1, atom) = Real(parts(atom)%xxx, c_double)
      kim_data%coords(2, atom) = Real(parts(atom)%yyy, c_double)
      kim_data%coords(3, atom) = Real(parts(atom)%zzz, c_double)
    End Do

    ! Set contributing status, atoms in the halo do not contribute
    kim_data%contributing = 0_c_int
    kim_data%contributing(1:natms) = 1_c_int

    ! Enter species information of species codes & contributing
    Do atom = 1, nlast
      species_code = kim_data%unique_species_code(ltype(atom))
      kim_data%species_code(atom) = species_code
      If (species_code < 0_c_int) Then
        ! If the current model does not support a specific species, means the
        ! species does not contribute. Here we mark the species and later
        ! correct the contributing array to only 0 and 1
        kim_data%contributing(atom) = -2_c_int
      End If
    End Do

    ! Construct neighbour lists
    Do list_index = 1, Int(kim_data%n_lists, wi)
      Call kim_neighbour_list(kim_data, list_index, nlast, natms, neigh_list)
    End Do

    ! Call KIM API to compute energy and forces
    Call kim_compute( &
      kim_data%model_handle, &
      kim_data%compute_arguments_handle, &
      kerror)
    If (kerror /= 0_c_int) Then
      Call kim_error('kim_energy_and_forces, kim_compute, '// &
                     'returned an error.', __LINE__)
    End If

    ! Retrieve KIM energy and forces from pointers (allocated in kim_setup)
    energy_kim = Real(kim_data%energy, wp)

    Do atom = 1, natms
      parts(atom)%fxx = parts(atom)%fxx + Real(kim_data%forces(1, atom), wp)
      parts(atom)%fyy = parts(atom)%fyy + Real(kim_data%forces(2, atom), wp)
      parts(atom)%fzz = parts(atom)%fzz + Real(kim_data%forces(3, atom), wp)
    End Do

    ! Distribute the partial forces on boundary atoms calculated on this processor
    ! to the appropriate neighbour
    Call kim_share_halo_forces(kim_data, parts, natms, nlast, lsi, lsa, &
                               ltg, map, comm)

    ! Virials (and pressure?)
    ! In OpenKIM virial has 6 components and is stored as
    ! a 6-element vector in the following order: xx, yy, zz, yz, xz, xy.
    stress(1) = stress(1) - Real(kim_data%virial(1), wp)
    stress(2) = stress(2) - Real(kim_data%virial(6), wp)
    stress(3) = stress(3) - Real(kim_data%virial(5), wp)
    stress(4) = stress(4) - Real(kim_data%virial(6), wp)
    stress(5) = stress(5) - Real(kim_data%virial(2), wp)
    stress(6) = stress(6) - Real(kim_data%virial(4), wp)
    stress(7) = stress(7) - Real(kim_data%virial(5), wp)
    stress(8) = stress(8) - Real(kim_data%virial(4), wp)
    stress(9) = stress(9) - Real(kim_data%virial(3), wp)

    virial_kim = Real(kim_data%virial(1) + kim_data%virial(2) + &
                      kim_data%virial(3), wp)
#else
    energy_kim = 0.0_wp
    virial_kim = 0.0_wp
#endif
  End Subroutine kim_energy_and_forces

  !> Prepare the KIM neighbour list used to compute energy and forces
  Subroutine kim_neighbour_list(kim_data, list_index, nlast, natms, list)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data
    !> Label of list to generate
    Integer(Kind=wi), Intent(In) :: list_index
    !> Number of particles in this domain and it's halo
    Integer(Kind=wi), Intent(In) :: nlast
    !> Number of particles in this domain (excluding the halo)
    Integer(Kind=wi), Intent(In) :: natms
    !> DL_POLY neighbour list
    !>
    !> - list(0,n) is the number of neighbours of particle n
    !> - list(m,n) is the 'm'th neighbour of particle n
    Integer(Kind=wi), Intent(In) :: list(-3:, 1:)

    Integer(Kind=wi) :: neighbour
    Logical :: padding_hint, any_non_supporting_species

    Integer(Kind=c_int) :: ipart, jpart
    Integer(Kind=c_int) :: c_natms, c_nlast
    Real(Kind=c_double) :: dr(1:3), r

#ifdef KIM
    !
    ! With the current KIM neighbour list implementation here in DL_POLY
    ! The user can combine the KIM model with any other internal DL_POLY
    ! interaction. While this is very interesting, It puts the
    ! responsibility on the user to make sure including all the additional
    ! interactions plus KIM interactions.
    !
    ! @todo
    ! Add extra check to the DL_POLY to make sure all the interactions are
    ! considered

    ! Check if any species in the simulation is not supported by the current
    ! model
    any_non_supporting_species = Any(kim_data%contributing < 0_c_int)

    padding_hint = kim_data%hints_padding(list_index) == 1_c_int

    c_natms = Int(natms, c_int)

    Associate (n_neigh=>kim_data%neigh(list_index)%n_neigh, &
               kim_list=>kim_data%neigh(list_index)%neigh_list, &
               cutoff=>kim_data%neigh(list_index)%cutoff)

      ! Set number of particles
      kim_data%neigh%n_part = c_natms

      ! Initialise neighbour list and number of neighbours
      kim_list = 0_c_int
      n_neigh = 0_c_int

      ! If all species in the simulation are supported by the current model
      If (.not. any_non_supporting_species) Then

        ! Build KIM neighbour list
        Do ipart = 1_c_int, c_natms
          Do neighbour = 1, list(0, ipart)
            ! Determine the 'neighbour'th neighbour of particles ipart
            jpart = list(neighbour, ipart)

            ! Add this neighbour to the KIM neighbour list of ipart
            n_neigh(ipart) = n_neigh(ipart) + 1_c_int

            kim_list(n_neigh(ipart), ipart) = jpart

            ! Add symmetric entry as KIM requires a full list. If the padding hint
            ! is true, then this is only nessecary if jpart is not a halo atom.
            If ((.not. padding_hint) .or. jpart <= c_natms) Then
              n_neigh(jpart) = n_neigh(jpart) + 1_c_int
              kim_list(n_neigh(jpart), jpart) = ipart
            End If ! padding_hint
          End Do ! neighbour
        End Do ! ipart

        ! Determine padding neighbours of padding particles if required. The
        ! non-padding neighbours of padding particles have been added to the
        ! neighbour list above.
        If (.not. padding_hint) Then
          c_nlast = Int(nlast, c_int)

          Do ipart = c_natms + 1_c_int, c_nlast - 1_c_int
            Do jpart = ipart + 1_c_int, c_nlast
              ! Add to list if the pair ipart-jpart are within the cutoff
              dr(1:3) = kim_data%coords(1:3, jpart) - &
                        kim_data%coords(1:3, ipart)
              r = Sqrt(Dot_product(dr, dr))

              If (r < cutoff) Then
                n_neigh(ipart) = n_neigh(ipart) + 1_c_int
                kim_list(n_neigh(ipart), ipart) = jpart

                ! Add symmetric entry
                n_neigh(jpart) = n_neigh(jpart) + 1_c_int
                kim_list(n_neigh(jpart), jpart) = ipart
              End If ! r < cutoff
            End Do ! jpart
          End Do ! ipart
        End If ! padding_hint

        ! If there is any non-supporting species by model in the simulation
      Else

        ! Build KIM neighbour list
        Do ipart = 1_c_int, c_natms
          ! If the atom is contributing we need its neighbours
          If (kim_data%contributing(ipart) > 0_c_int) Then
            Do neighbour = 1, list(0, ipart)
              ! Determine the 'neighbour'th neighbour of particles ipart
              jpart = list(neighbour, ipart)

              ! If the species is supported by model we add it to the neigbour list
              If (kim_data%contributing(jpart) > -1_c_int) Then
                ! Add this neighbour to the KIM neighbour list of ipart
                n_neigh(ipart) = n_neigh(ipart) + 1_c_int

                kim_list(n_neigh(ipart), ipart) = jpart

                ! Add symmetric entry as KIM requires a full list. If the padding hint
                ! is true, then this is only nessecary if jpart is not a halo atom.
                If ((.not. padding_hint) .or. jpart <= c_natms) Then
                  n_neigh(jpart) = n_neigh(jpart) + 1_c_int
                  kim_list(n_neigh(jpart), jpart) = ipart
                End If ! padding_hint
              End If ! contributing(jpart)
            End Do ! neighbour
          End If ! contributing(ipart)
        End Do ! ipart

        ! Determine padding neighbours of padding particles if required. The
        ! non-padding neighbours of padding particles have been added to the
        ! neighbour list above.
        If (.not. padding_hint) Then
          c_nlast = Int(nlast, c_int)

          Do ipart = c_natms + 1_c_int, c_nlast - 1_c_int
            ! If the species is supported by model we want the neigbour list
            If (kim_data%contributing(ipart) > -1_c_int) Then
              Do jpart = ipart + 1_c_int, c_nlast
                ! If the species is supported by model we check to add it to the neigbour list
                If (kim_data%contributing(jpart) > -1_c_int) Then
                  ! Add to list if the pair ipart-jpart are within the cutoff
                  dr(1:3) = kim_data%coords(1:3, jpart) - &
                            kim_data%coords(1:3, ipart)
                  r = Sqrt(Dot_product(dr, dr))

                  If (r < cutoff) Then
                    n_neigh(ipart) = n_neigh(ipart) + 1_c_int
                    kim_list(n_neigh(ipart), ipart) = jpart

                    ! Add symmetric entry
                    n_neigh(jpart) = n_neigh(jpart) + 1_c_int
                    kim_list(n_neigh(jpart), jpart) = ipart
                  End If ! r < cutoff
                End If ! contributing(jpart)
              End Do ! jpart
            End If ! contributing(ipart)
          End Do ! ipart
        End If ! padding_hint

        ! Correct the array to contributing and non-contributing atoms
        kim_data%contributing = Merge(0_c_int, kim_data%contributing, &
                                      (kim_data%contributing < 0_c_int))

      End If ! any_non_supporting_species

    End Associate
#endif
  End Subroutine kim_neighbour_list

  !> KIM does not assume that the entirity of an atoms force will be calculated
  !> on the domain to which it belongs. It is therefore necessary to send the
  !> partial forces calculated on padding atoms to the appropriate domain.
  Subroutine kim_share_halo_forces(kim_data, parts, natms, nlast, lsi, lsa, ltg, &
                                   map, comm)
    !> KIM data type
    Type(kim_type), Target,       Intent(InOut) :: kim_data
    !> Particles
    Type(corepart), Dimension(:), Intent(InOut) :: parts
    !> Number of particles in this domain (excluding the halo)
    Integer(Kind=wi),             Intent(In   ) :: natms
    !> Number of particles in this domain and it's halo
    Integer(Kind=wi),             Intent(In   ) :: nlast
    !> Some sort of local to global arrays
    Integer(Kind=wi),             Intent(In   ) :: lsi(:), lsa(:)
    !> Local to global id
    Integer(Kind=wi),             Intent(In   ) :: ltg(:)
    !> Map of neighbouring domain ids
    Integer(Kind=wi),             Intent(In   ) :: map(1:26)
    !> Communicator type
    Type(comms_type),             Intent(InOut) :: comm

    !> Direction to send and receive in
    Integer(Kind=wi) :: direction
    !> id of the node to send and receive from respectively
    Integer(Kind=wi) :: destination, source
    !> Local atom number iterator
    Integer(Kind=wi) :: atom
    !> Buffer location
    Integer(Kind=wi) :: buf
    !> Size of buffer to send and receive respectively
    Integer(Kind=wi) :: send_size, recv_size
    !> The position in the buffer where data to receive begins
    Integer(Kind=wi) :: recv_start
    Integer(Kind=wi) :: i

#ifdef KIM
    ! Iterate over each direction to send and receive
    Do direction = 1, 6

      ! Determine the id of the appropriate neighbour
      If (direction == 1) Then
        ! -x
        destination = map(1)
        source = map(2)
      Else If (direction == 2) Then
        ! +x
        destination = map(2)
        source = map(1)
      Else If (direction == 3) Then
        ! -y
        destination = map(3)
        source = map(4)
      Else If (direction == 4) Then
        ! +y
        destination = map(4)
        source = map(3)
      Else If (direction == 5) Then
        ! -z
        destination = map(5)
        source = map(6)
      Else If (direction == 6) Then
        ! +z
        destination = map(6)
        source = map(5)
      End If

      Associate (buffer=>kim_data%kcomms%buffer)
        ! Pack data to send
        buf = 0
        ! Loop over atoms to send
        Do atom = kim_data%kcomms%first(direction), &
          kim_data%kcomms%last(direction)
          ! Pack forces
          buffer(buf + 1) = Real(kim_data%forces(1, atom), wp)
          buffer(buf + 2) = Real(kim_data%forces(2, atom), wp)
          buffer(buf + 3) = Real(kim_data%forces(3, atom), wp)

          ! Pack global atom id
          buffer(buf + 4) = Real(ltg(atom), wp)

          buf = buf + span
        End Do

        ! Determine buffer size to send
        send_size = buf

        ! Exchange buffers
        If (comm%mxnode > 1) Then
          ! Determine buffer size to receive
          recv_size = kim_data%kcomms%n_recv(direction) * span
          recv_start = kim_data%kcomms%recv_start
          If (recv_size > 0) Then
            Call girecv(comm, buffer(recv_start + 1:recv_start + &
                                     recv_size), source, Export_tag)
            Call gsend(comm, buffer(1:send_size), destination, Export_tag)
            Call gwait(comm)
          End If
        Else
          recv_size = send_size
          recv_start = kim_data%kcomms%recv_start
        End If

        ! Unpack received data
        ! Determine where received data begins in the buffer
        buf = Merge(recv_start, 0, comm%mxnode > 1)
        ! Loop number over atoms to receive
        Do i = 1, recv_size / span
          ! Find local index of received atom
          atom = local_index(Nint(buffer(buf + 4)), nlast, lsi, lsa)

          ! Add partial force to local force arrays
          parts(atom)%fxx = parts(atom)%fxx + buffer(buf + 1)
          parts(atom)%fyy = parts(atom)%fyy + buffer(buf + 2)
          parts(atom)%fzz = parts(atom)%fzz + buffer(buf + 3)

          buf = buf + span
        End Do
      End Associate
    End Do
#endif
  End Subroutine kim_share_halo_forces

  !> Neighbour list access function (called by KIM API)
  Subroutine get_neigh(data_object, n_lists, cutoffs, neighbour_list_index, &
                       request, n_neigh, first, kerror) bind(c)
    !> KIM neighbour type address
    Type(c_ptr), Value,    Intent(In)  :: data_object
    !> Number of potential cutoffs in model
    Integer(c_int), Value, Intent(In)  :: n_lists
    !> Array of cutoffs
    Real(c_double),        Intent(In)  :: cutoffs(n_lists)
    !> Which neighbour list to inspect
    Integer(c_int), Value, Intent(In)  :: neighbour_list_index
    !> Number of the particle for which list information is requested
    Integer(c_int), Value, Intent(In)  :: request
    !> Number of neighbours of the requested particle
    Integer(c_int),        Intent(Out) :: n_neigh
    !> Location of the first neighbour of the returned neighbour list
    Type(c_ptr),           Intent(Out) :: first
    !> Error status
    Integer(c_int),        Intent(Out) :: kerror

    Type(kim_neighbour_list_pointer_type), Pointer :: neigh(:)
    Integer(c_int), Pointer                        :: kim_neigh(:)
    Integer(c_int), Pointer                        :: kim_list(:, :)
    Real(c_double), Pointer                        :: kim_cutoff
    Integer(c_int), Pointer                        :: kim_n_part
    Character(Len=80)                              :: message

#ifdef KIM
    ! Assign local KIM neighbour type pointer to address of C data object
    Call c_f_pointer(data_object, neigh, [n_lists])
    Call c_f_pointer(neigh(neighbour_list_index)%n_neigh, kim_neigh, &
                     [neigh(neighbour_list_index)%max_atoms])
    Call c_f_pointer(neigh(neighbour_list_index)%neigh_list, kim_list, &
                     [neigh(neighbour_list_index)%max_neigh, &
                      neigh(neighbour_list_index)%max_atoms])
    Call c_f_pointer(neigh(neighbour_list_index)%cutoff, kim_cutoff)
    Call c_f_pointer(neigh(neighbour_list_index)%n_part, kim_n_part)

    ! Ensure neighbour list cutoff is larger than KIM model cutoff
    If (cutoffs(neighbour_list_index) > kim_cutoff) Then
      Call kim_warning('get_neigh, Neighbour list cutoff too small for '// &
                       'model cutoff.', __LINE__)
      kerror = 1_c_int
      Return
    End If

    ! Ensure requested particle exists
    If ((request > kim_n_part) .or. (request < 1)) Then
      Write (message, '(a, g0)') 'get_neigh, Invalid particle number '// &
        'requested: ', request
      Call kim_warning(message, __LINE__)
      kerror = 1_c_int
      Return
    End If

    ! Set number of neighbours for requested particle
    n_neigh = kim_neigh(request)
    ! Set the location for beginning of the returned neighbour list
    first = c_loc(kim_list(1, request))

    kerror = 0_c_int
#else
    n_neigh = 0_c_int
    first = c_null_ptr
    kerror = 0_c_int
#endif
  End Subroutine get_neigh

  !> Initialise the KIM neighbour list type
  Subroutine kim_neighbour_list_type_init(T, max_atoms, max_list)
    Class(kim_neighbour_list_type)  :: T
    Integer(Kind=wi), Intent(In   ) :: max_atoms, max_list

    Integer(Kind=wi) :: fail

!> Extent of particle arrays. When the padding hint is true this is the
!> maximum number of atoms on a node excluding the halo. When the hint is
!> false this is the maximum number of atoms on a node including the halo.
!> Extent of neighbour list arrays

    fail = 0
    ! Allocate KIM neighbour list type arrays
    Allocate (T%n_neigh(max_atoms), Stat=fail)
    If (fail /= 0) Then
      Call kim_error('kim_neighbour_list_type_init, n_neigh '// &
                     'allocation failure.', __LINE__)
    End If

    Allocate (T%neigh_list(max_list, max_atoms), Stat=fail)
    If (fail /= 0) Then
      Call kim_error('kim_neighbour_list_type_init, neigh_list '// &
                     'allocation failure.', __LINE__)
    End If
  End Subroutine kim_neighbour_list_type_init

  !> Deallocate memory used by the KIM neighbour list type
  Subroutine kim_neighbour_list_type_cleanup(T)
    Type(kim_neighbour_list_type) :: T

    Integer(Kind=wi) :: fail

    fail = 0

    If (Allocated(T%n_neigh)) Then
      Deallocate (T%n_neigh, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_neighbour_list_type_cleanup, n_neigh '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    If (Allocated(T%neigh_list)) Then
      Deallocate (T%neigh_list, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_neighbour_list_type_cleanup, neigh_list '// &
                       'deallocation failure.', __LINE__)
      End If
    End If
  End Subroutine kim_neighbour_list_type_cleanup

  !> Initialise the KIM neighbour list pointer type
  Subroutine kim_neighbour_list_pointer_type_init(neigh_pointer, neigh, &
                                                  max_atoms, max_list)
    Type(kim_neighbour_list_pointer_type), Intent(InOut) :: neigh_pointer
    Type(kim_neighbour_list_type), Target, Intent(In   ) :: neigh
    Integer(Kind=wi),                      Intent(In   ) :: max_atoms, max_list

    neigh_pointer%cutoff = c_loc(neigh%cutoff)
    neigh_pointer%n_part = c_loc(neigh%n_part)
    neigh_pointer%neigh_list = c_loc(neigh%neigh_list)
    neigh_pointer%n_neigh = c_loc(neigh%n_neigh)
    neigh_pointer%max_atoms = Int(max_atoms, c_int)
    neigh_pointer%max_neigh = Int(max_list, c_int)
  End Subroutine kim_neighbour_list_pointer_type_init

  !> Initialise the KIM comms type
  Subroutine kim_comms_type_init(T, buffer_size, mxnode)
    Class(kim_comms_type)           :: T
    Integer(Kind=wi), Intent(In   ) :: buffer_size, mxnode

    Integer(Kind=wi) :: fail

!> Size of the buffer to allocate
!> Number of nodes

    fail = 0

    ! Allocate buffer
    Allocate (T%buffer(buffer_size), Stat=fail)
    If (fail /= 0) Then
      Call kim_error('kim_comms_type_init, buffer '// &
                     'allocation failure.', __LINE__)
    End If

    ! Set begining of receive portion
    T%recv_start = buffer_size / Merge(2, 1, mxnode > 1)
  End Subroutine kim_comms_type_init

  !> Set the ids of atoms to be sent in kim_share_halo_forces
  Subroutine kim_comms_type_set(T, direction, n_recv, first, last)
    Class(kim_comms_type)           :: T
    Integer(Kind=wi), Intent(In   ) :: direction, n_recv
    Integer(Kind=wi)                :: first, last

!> The direction in which data is being sent and received
!> Number of atoms to receive in each direction
!> ids of first and last atoms to send respectively in each direction

    T%n_recv(direction) = n_recv
    T%first(direction) = first
    T%last(direction) = last
  End Subroutine kim_comms_type_set

  !> Deallocate memory from KIM comms type
  Subroutine kim_comms_type_cleanup(T)
    Type(kim_comms_type) :: T

    Integer(Kind=wi) :: fail

    fail = 0

    If (Allocated(T%buffer)) Then
      Deallocate (T%buffer, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_comms_type_cleanup, buffer '// &
                       'deallocation failure.', __LINE__)
      End If
    End If
  End Subroutine kim_comms_type_cleanup

  !> Allocate memory for the KIM type
  Subroutine kim_type_init(T, mxatms)
    Class(kim_type) :: T
    !> Extent of particle arrays including halo particles (coordinates, forces, etc.)
    Integer(Kind=wi), Intent(In) :: mxatms

    Integer(Kind=wi) :: fail(5)
#ifdef KIM
    fail = 0

    Allocate (T%species_code(mxatms), Stat=fail(2))
    Allocate (T%contributing(mxatms), Stat=fail(3))
    Allocate (T%coords(3, mxatms), Stat=fail(4))
    Allocate (T%forces(3, mxatms), Stat=fail(5))

    If (Any(fail /= 0)) Then
      Call kim_error('kim_type_init, allocation failure.', __LINE__)
    End If
#endif
  End Subroutine kim_type_init

  !> Deallocate memory from KIM type
  Subroutine kim_type_cleanup(T)
    Type(kim_type) :: T

    Integer(Kind=c_int) :: kerror
    Integer(Kind=wi)    :: fail

    fail = 0

#ifdef KIM
    If (Allocated(T%species_code)) Then
      Deallocate (T%species_code, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_type_cleanup, species_code '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    If (Allocated(T%unique_species_code)) Then
      Deallocate (T%unique_species_code, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_type_cleanup, unique_species_code '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    If (Allocated(T%forces)) Then
      Deallocate (T%forces, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_type_cleanup, forces '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    If (Allocated(T%coords)) Then
      Deallocate (T%coords, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_type_cleanup, coords '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    If (Allocated(T%neigh)) Then
      Deallocate (T%neigh, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_type_cleanup, neigh '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    If (Allocated(T%neigh_pointer)) Then
      Deallocate (T%neigh_pointer, Stat=fail)
      If (fail /= 0) Then
        Call kim_error('kim_type_cleanup, neigh_pointer '// &
                       'deallocation failure.', __LINE__)
      End If
    End If

    Call kim_compute_arguments_destroy( &
      T%model_handle, &
      T%compute_arguments_handle, &
      kerror)
    If (kerror /= 0_c_int) Then
      Call kim_error('kim_type_cleanup, kim_compute_arguments_destroy '// &
                     'returned an error.', __LINE__)
    End If

    Call kim_model_destroy(T%model_handle)
#endif
  End Subroutine kim_type_cleanup
End Module kim
