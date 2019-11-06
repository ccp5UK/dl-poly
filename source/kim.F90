!> Module providing support for the KIM v2.0 library
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J. Madge July 2018
!>
!> Based on KIM v2.0 simulator example in Fortran and the KIM v1 implementation
!> by R. S. Elliot
Module kim
  Use, Intrinsic :: iso_c_binding, Only : c_double,c_int,c_char,c_ptr, &
    c_null_ptr,c_funloc,c_loc,c_f_pointer
  Use kinds, Only : wi,wp
  Use numerics, Only: local_index
  Use particle, Only : corepart
  Use comms, Only : comms_type,gsend,girecv,gwait,Export_tag
  Use errors_warnings, Only : error,warning
#ifdef KIM
  Use kim_simulator_headers_module, Only : &
    kim_model_create, &
    kim_model_destroy, &
    kim_compute, &
    kim_compute_arguments_create, &
    kim_compute_arguments_destroy, &
    kim_compute_arguments_handle_type, &
    kim_set_argument_pointer, &
    kim_set_callback_pointer, &
    kim_compute_argument_name_coordinates, &
    kim_compute_argument_name_number_of_particles, &
    kim_compute_argument_name_partial_energy, &
    kim_compute_argument_name_partial_forces, &
    kim_compute_argument_name_partial_virial, &
    kim_compute_argument_name_particle_contributing, &
    kim_compute_argument_name_particle_species_codes, &
    kim_compute_callback_name_get_neighbor_list, &
    kim_language_name_fortran, &
    kim_get_species_support_and_code, &
    kim_get_influence_distance, &
    kim_get_number_of_neighbor_lists, &
    kim_get_neighbor_list_values, &
    kim_model_handle_type, &
    kim_numbering_one_based, &
    kim_from_string, &
    kim_species_name_type, &
    kim_length_unit_a, kim_energy_unit_amu_a2_per_ps2, kim_charge_unit_e, &
    kim_temperature_unit_k, kim_time_unit_ps
#endif
  Implicit None

  Private

  !> Flag indicating whether the program has been compiled with the KIM preprocessor flag
#ifdef KIM
  Logical, Parameter, Public :: COMPILED_WITH_KIM = .true.
#else
  Logical, Parameter, Public :: COMPILED_WITH_KIM = .false.
#endif

  !> Name of this source file
  Character(Len=*), Parameter :: FILE_NAME="kim.F90"

  !> Number of scalars per particle to be sent during communication.
  !> In this case force_x, force_y, > force_z and global_id
  Integer(Kind=wi), Parameter :: span = 4

  !> Neighbour list type containing KIM neighbour list data
  Type :: kim_neighbour_list_type
    Private

    !> Potential cut off
    Real(Kind=c_double) :: cutoff
    !> Total number of particles
    Integer(Kind=c_int) :: n_part
    !> Number of neighbours for each particle
    Integer(Kind=c_int), Allocatable :: n_neigh(:)
    !> Neighbour list
    !>
    !> `list(n,part)` is the number of the `n`th neighbour of particle `part`
    Integer(Kind=c_int), Allocatable :: neigh_list(:,:)
  Contains
    Private

    Procedure :: init => kim_neighbour_list_type_init
    Final :: kim_neighbour_list_type_cleanup
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
    Type(c_ptr) :: cutoff
    !> Pointer to kim_neighbour_list_type n_part
    Type(c_ptr) :: n_part
    !> Pointer to kim_neighbour_list_type n_neigh
    Type(c_ptr) :: n_neigh
    !> Pointer to kim_neighbour_list_type neigh_list
    Type(c_ptr) :: neigh_list
  End Type kim_neighbour_list_pointer_type

  !> Type containing data necessary for the KIM message passing in the routine
  !> kim_share_halo_forces
  Type :: kim_comms_type
    Private

    !> Buffer for message passing
    Real(Kind=wp), Allocatable :: buffer(:)
    !> Half buffer size, used to divide the portion for sending and the portion
    !> for receiving
    Integer(Kind=wi) :: recv_start
    !> Number of atoms to receive in each direction
    Integer(Kind=wi) :: n_recv(6)
    !> ids of first and last atoms to send respectively in each direction
    Integer(Kind=wi) :: first(6), last(6)

  Contains
    Private

    Procedure :: init => kim_comms_type_init
    Procedure, Public :: set => kim_comms_type_set
    Final :: kim_comms_type_cleanup
  End Type kim_comms_type

  !> Type containing KIM data
  Type, Public :: kim_type
    Private

    !> Flag indicating whether KIM is in use
    Logical, Public :: active = .false.
    !> Flag indicating if any neighour lists do not have the padding hint and so
    !> require neighbours of non-contributing (padding) particles
    Logical, Public :: padding_neighbours_required = .false.

    !> Name of the KIM model requested
    Character(Len=:, Kind=c_char), Allocatable, Public :: model_name

    !> Largest model cutoff, used to determining cell list and domain size
    Real(Kind=wp), Public :: cutoff = 0.0_wp
    !> Model influence distance
    Real(Kind=c_double), Public :: influence_distance

    !> Neighbour list type
    Type(kim_neighbour_list_type), Allocatable :: neigh(:)
    !> Neighbour lsit pointer type
    Type(kim_neighbour_list_pointer_type), Allocatable :: neigh_pointer(:)
    !> Comms type
    Type(kim_comms_type), Public :: kcomms

#ifdef KIM
    !> Model handle
    Type(kim_model_handle_type) :: model_handle
    !> Compute arguments handle
    Type(kim_compute_arguments_handle_type) :: compute_arguments_handle

    !> Number of particles
    Integer(Kind=c_int) :: n_particles

    !> Number of model neighbour lists
    Integer(Kind=c_int) :: n_lists

    !> Neighbour list padding hints
    Integer(Kind=c_int), Allocatable :: hints_padding(:)

    ! Species
    !> List of species names
    Type(kim_species_name_type), Allocatable :: species_name(:)
    !> List of species codes
    Integer(Kind=c_int), Allocatable :: species_code(:)

    !> Contributing list
    !>
    !> - 1 if particle contributes
    !> - 0 otherwise
    Integer(Kind=c_int), Allocatable :: contributing(:)

    !> Coordinates
    Real(Kind=c_double), Allocatable :: coords(:,:)

    !> KIM energy
    Real(Kind=c_double) :: energy
    !> KIM forces
    Real(Kind=c_double), Allocatable :: forces(:,:)
    !> KIM virial
    Real(Kind=c_double) :: virial(6)
#endif
  Contains
    Private

    Procedure :: init => kim_type_init
    Final :: kim_type_cleanup
  End Type kim_type

  Public :: kim_setup, kim_cutoff
  Public :: kim_energy_and_forces
  Public :: get_neigh

Contains

  !> Interface to DL_POLY error routine
  Subroutine kim_error(message,line)
    Character(Kind=c_char, Len=*), Intent(In) :: message
    Integer(Kind=wi), Intent(In) :: line

    Character(Len=3) :: line_str
    Character(Len=:), Allocatable :: error_message

    Write(line_str,'(i3)') line
    error_message = 'KIM error: ' // Trim(message) // ' line: ' // &
      Trim(Adjustl(line_str)) // ' file: ' // FILE_NAME
    Call error(0,error_message,.true.)
  End subroutine kim_error

  !> Interface to DL_POLY warning routine
  Subroutine kim_warning(message,line)
    Character(Kind=c_char, Len=*), Intent(In) :: message
    Integer(Kind=wi), Intent(In) :: line

    Character(Len=3) :: line_str
    Character(Len=:), Allocatable :: warning_message

    Write(line_str,'(i3)') line
    warning_message = 'KIM warning: ' // Trim(message) // ' line: ' // &
      Trim(Adjustl(line_str)) // ' file: ' // FILE_NAME
    Call warning(warning_message,.true.)
  End subroutine kim_warning

  !> Initialise KIM types
  Subroutine kim_setup(kim_data,mxatms,mxatdm,megatm,max_list,mxbfxp,mxnode)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data
    !> Extent of particle arrays including halo particles (coordinates, forces, etc.)
    Integer(Kind=wi), Intent(In) :: mxatms
    !> Extent of particle arrays not including halo particles (neighbours, contraints, etc.)
    Integer(Kind=wi), Intent(In) :: mxatdm
    !> Total number of particles on all domains/nodes
    Integer(Kind=wi), Intent(In) :: megatm
    !> Extent of neighbour list arrays
    Integer(Kind=wi), Intent(In) :: max_list
    !> Extent of force arrays
    Integer(Kind=wi), Intent(In) :: mxbfxp
    !> Number of nodes
    Integer(Kind=wi), Intent(In) :: mxnode

    Integer(Kind=c_int) :: requested_units_accepted
    Real(Kind=c_double), Allocatable :: cutoffs(:)
    Integer(Kind=wi) :: list_index, max_atoms
    Integer(Kind=c_int) :: kerror

    ! Early return if a KIM model is not being used
    If (.not. kim_data%active) Then
      Return
    End If

#ifdef KIM
    ! Allocate KIM data type arrays
    Call kim_data%init(mxatms)

    ! Initialise comms type
    Call kim_data%kcomms%init(mxbfxp*span,mxnode)

    ! Create KIM object
    ! kim_numbering_one_based - numbering from 1 (fortran style arrays)
    !
    ! DL_POLY internal units:
    ! - length Angstrom
    ! - energy 10J/mol (equivalent KIM energy unit selected here)
    ! - charge electrons
    ! - temperature K
    ! - time ps
    Call kim_model_create(kim_numbering_one_based, &
      kim_length_unit_a, &
      kim_energy_unit_amu_a2_per_ps2, &
      kim_charge_unit_e, &
      kim_temperature_unit_k, &
      kim_time_unit_ps, &
      Trim(kim_data%model_name), &
      requested_units_accepted, &
      kim_data%model_handle, kerror)
    If (requested_units_accepted /= 1) Then
      Call kim_error('kim_model_create, the selected KIM model does not support DL_POLY internal units', __LINE__)
    End If
    If (kerror /= 0) Then
      Call kim_error('kim_model_create', __LINE__)
    End If

    ! Create compute_arguments object
    Call kim_compute_arguments_create(kim_data%model_handle, &
      kim_data%compute_arguments_handle,kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_create', __LINE__)
    End If

    ! Set number of particles
    kim_data%n_particles = megatm

    ! Get influence distance
    Call kim_get_influence_distance(kim_data%model_handle, &
      kim_data%influence_distance)

    ! Get number of neighbour lists
    Call kim_get_number_of_neighbor_lists(kim_data%model_handle, &
      kim_data%n_lists)

    ! Allocate neighbour list, neighbour list pointer
    Allocate(kim_data%neigh(kim_data%n_lists))
    Allocate(kim_data%neigh_pointer(kim_data%n_lists))
    ! Allocate hint arrays
    Allocate(kim_data%hints_padding(kim_data%n_lists))

    ! Get neighbour list cutoffs and hints
    Allocate(cutoffs(kim_data%n_lists))
    Call kim_get_neighbor_list_values(kim_data%model_handle, &
      cutoffs, kim_data%hints_padding, kerror)
    Do list_index = 1, kim_data%n_lists
      kim_data%neigh(list_index)%cutoff = cutoffs(list_index)
    End Do
    Deallocate(cutoffs)
    ! Initialise neighbour list and neighbour list pointer types
    Do list_index = 1, kim_data%n_lists
      If (kim_data%hints_padding(list_index) == 1) Then
        max_atoms = mxatdm
      Else
        max_atoms = mxatms
      End If
      ! Initialise neighbour list type
      Call kim_data%neigh(list_index)%init(max_atoms,max_list)

      ! Initialise neighbour list pointer type
      Call kim_neighbour_list_pointer_type_init(kim_data%neigh_pointer(list_index), &
        kim_data%neigh(list_index),max_atoms,max_list)
    End Do

    If (kerror /= 0) Then
      Call kim_error('kim_get_neighbor_list_values', __LINE__)
    End If

    ! Allocate KIM pointers
    ! Number of particles
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_number_of_particles, &
      kim_data%n_particles, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_set_arugment_pointer for number of particles', __LINE__)
    End If

    ! Species codes
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_particle_species_codes, &
      kim_data%species_code, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_set_arugment_pointer for species codes', __LINE__)
    End If

    ! Contributing status
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_particle_contributing, &
      kim_data%contributing, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_set_arugment_pointer for contributing status', __LINE__)
    End If

    ! Coordinates
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_coordinates, &
      kim_data%coords, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_set_arugment_pointer for coordinates', __LINE__)
    End If

    ! Energy
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_partial_energy, &
      kim_data%energy, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_set_arugment_pointer for energy', __LINE__)
    End If

    ! Forces
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_partial_forces, &
      kim_data%forces, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_compute_arguments_set_arugment_pointer for forces', __LINE__)
    End If

    ! Virials
    Call kim_set_argument_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_argument_name_partial_virial, &
      kim_data%virial, kerror)
    If (kerror /= 0) Then
      Call kim_warning('The selected KIM model does not compute virials, stress and pressure will be incorrect', __LINE__)
      ! Set KIM virials to 0 so they will not contribute to the total
      kim_data%virial = 0.0_wp
    End If

    ! Set KIM pointer to neighbour list routine and type
    Call kim_set_callback_pointer( &
      kim_data%compute_arguments_handle, &
      kim_compute_callback_name_get_neighbor_list, kim_language_name_fortran, &
      C_funloc(get_neigh), C_loc(kim_data%neigh_pointer), kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_set_callback_pointer', __LINE__)
    End If
#endif
  End Subroutine kim_setup

  !> Function to return the cutoff of the kim model 'model_name'.
  !>
  !> This is necessary as the cutoff is required for ensuring domains and cells
  !> are of an appropriate size. Perhaps with a reworking of the way the input
  !> files are read could remedy this clumsy implementation.
  Subroutine kim_cutoff(kim_data)
    Type(kim_type), Target, Intent(InOut) :: kim_data

#ifdef KIM
    Type(kim_model_handle_type) :: model_handle
#endif
    Integer(Kind=c_int) :: requested_units_accepted
    Integer(Kind=c_int) :: kerror
    Integer(Kind=c_int) :: n_lists
    Real(Kind=c_double), Allocatable :: cutoffs(:)
    Integer(Kind=c_int), Allocatable :: hints_padding(:)

    If (COMPILED_WITH_KIM .eqv. .false.) Then
      Call error(0,'KIM directive found in FIELD, but the program is not built with openKIM support',.true.)
    End If

#ifdef KIM
    Call kim_model_create(kim_numbering_one_based, &
      kim_length_unit_a, &
      kim_energy_unit_amu_a2_per_ps2, &
      kim_charge_unit_e, &
      kim_temperature_unit_k, &
      kim_time_unit_ps, &
      Trim(kim_data%model_name), &
      requested_units_accepted, &
      model_handle, kerror)
    If (requested_units_accepted == 0) Then
      Call kim_error('kim_model_create, the selected KIM model does not support DL_POLY internal units', __LINE__)
    End If
    If (kerror /= 0) Then
      Call kim_error('kim_cutoff', __LINE__)
    End If

    ! Get number of neighbour lists
    Call kim_get_number_of_neighbor_lists(model_handle, n_lists)

    ! Get neighbour list cutoffs and record maximum
    Allocate(cutoffs(n_lists), hints_padding(n_lists))
    Call kim_get_neighbor_list_values(model_handle, &
      cutoffs, hints_padding, kerror)
    kim_data%cutoff = Real(maxval(cutoffs), wp)

    ! Check if the padding hint is defined
    If (Any(hints_padding /= 1)) Then
      Call kim_warning('The selected KIM model requires neighbours of ' // &
        'non-contributing particles. This may significantly affect performance', &
        __LINE__)
      kim_data%padding_neighbours_required = .true.
    End If

    ! Deallocate temporary variables
    Deallocate(cutoffs, hints_padding)
    Call kim_model_destroy(model_handle)
#endif
  End Subroutine kim_cutoff

  !> Compute KIM energy and forces
  Subroutine kim_energy_and_forces(kim_data,natms,nlast,parts,neigh_list,map, &
      lsite,lsi,lsa,ltg,site_name,energy_kim,virial_kim,stress,comm)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data
    !> Number of particles in this domain (excluding the halo)
    Integer(Kind=wi), Intent(In) :: natms
    !> Number of particles in this domain and it's halo
    Integer(Kind=wi), Intent(In) :: nlast
    !> Particles
    Type(corepart), Dimension(:), Intent(InOut) :: parts
    !> DL_POLY neighbour list
    !>
    !> - list(0,n) is the number of neighbours of particle n
    !> - list(m,n) is the 'm'th neighbour of particle n
    Integer(Kind=wi), Intent(In) :: neigh_list(-3:,1:)
    !> Map of neighbouring domain ids
    Integer(Kind=wi), Intent(In) :: map(1:26)
    !> Site local index to global index
    Integer(Kind=wi), Intent(In) :: lsite(:)
    !> Some sort of local to global arrays
    Integer(Kind=wi), Intent(In) :: lsi(:),lsa(:)
    !> Local to global id
    Integer(Kind=wi), Intent(In) :: ltg(:)
    !> Names of each atom
    Character(Len=8), Intent(In) :: site_name(:)
    !> KIM model energy
    Real(Kind=wp), Intent(Out) :: energy_kim
    !> KIM model virial
    Real(Kind=wp), Intent(Out) :: virial_kim
    !> Virial
    Real(Kind=wp), Dimension(9), Intent(InOut) :: stress
    !> Comms data
    Type(comms_type), Intent(InOut) :: comm

    Integer(Kind=wi) :: atom, list_index
    Integer(Kind=c_int) :: species_is_supported
    Integer(Kind=c_int) :: kerror

#ifdef KIM
    ! Set number of particles
    kim_data%n_particles = Int(nlast, c_int)

    ! Copy coordinates to KIM coordinates array
    kim_data%coords = 0
    Do atom=1,nlast
      kim_data%coords(1,atom) = Real(parts(atom)%xxx, c_double)
      kim_data%coords(2,atom) = Real(parts(atom)%yyy, c_double)
      kim_data%coords(3,atom) = Real(parts(atom)%zzz, c_double)
    End Do

    ! Enter species information
    Do atom = 1, nlast
      Call kim_from_string(Trim(site_name(lsite(atom))), &
        kim_data%species_name(atom))

      ! Check model supports the requested species
      Call kim_get_species_support_and_code(kim_data%model_handle, &
        kim_data%species_name(atom),species_is_supported, &
        kim_data%species_code(atom),kerror)

      If ( (kerror/=0) .or. (species_is_supported/=1) ) Then
        Call kim_error('Model does not support species', __LINE__)
      End If
    End Do

    ! Set contributing status, atoms in the halo do not contribute
    kim_data%contributing = 0
    kim_data%contributing(1:natms) = 1

    ! Construct neighbour lists
    Do list_index = 1, kim_data%n_lists
      Call kim_neighbour_list(kim_data,list_index,nlast,natms,neigh_list)
    End Do

    ! Call KIM API to compute energy and forces
    Call kim_compute(kim_data%model_handle, &
      kim_data%compute_arguments_handle, kerror)
    If (kerror /= 0) Then
      Call kim_error('kim_model_compute returned an error', __LINE__)
    End If

    ! Retrieve KIM energy and forces from pointers (allocated in kim_setup)
    energy_kim = Real(kim_data%energy, wp)
    Do atom = 1, natms
      parts(atom)%fxx = parts(atom)%fxx + Real(kim_data%forces(1,atom), wp)
      parts(atom)%fyy = parts(atom)%fyy + Real(kim_data%forces(2,atom), wp)
      parts(atom)%fzz = parts(atom)%fzz + Real(kim_data%forces(3,atom), wp)
    End Do
    ! Distribute the partial forces on boundary atoms calculated on this processor
    ! to the appropriate neighbour
    Call kim_share_halo_forces(kim_data,parts,natms,nlast,lsi,lsa,ltg,map,comm)

    ! Virials (and pressure?)
    ! Why the strange numbers?
    stress(1) = stress(1) - Real(kim_data%virial(1), wp)
    stress(2) = stress(2) - Real(kim_data%virial(6), wp)
    stress(3) = stress(3) - Real(kim_data%virial(5), wp)
    stress(4) = stress(4) - Real(kim_data%virial(6), wp)
    stress(5) = stress(5) - Real(kim_data%virial(2), wp)
    stress(6) = stress(6) - Real(kim_data%virial(4), wp)
    stress(7) = stress(7) - Real(kim_data%virial(5), wp)
    stress(8) = stress(8) - Real(kim_data%virial(4), wp)
    stress(9) = stress(9) - Real(kim_data%virial(3), wp)

    virial_kim = Real(kim_data%virial(1)+kim_data%virial(2)+kim_data%virial(3), wp)
#else
    energy_kim = 0.0_wp
    virial_kim = 0.0_wp
#endif
  End Subroutine kim_energy_and_forces

  !> Prepare the KIM neighbour list used to compute energy and forces
  Subroutine kim_neighbour_list(kim_data,list_index,nlast,natms,list)
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
    Integer(Kind=wi), Intent(In) :: list(-3:,1:)

    Integer(Kind=wi) :: ipart, jpart
    Integer(Kind=wi) :: neighbour
    Logical :: padding_hint
    Real(Kind=c_double) :: dr(1:3), r

#ifdef KIM
    padding_hint = kim_data%hints_padding(list_index) == 1

    Associate(n_neigh => kim_data%neigh(list_index)%n_neigh, &
        kim_list => kim_data%neigh(list_index)%neigh_list, &
        cutoff => kim_data%neigh(list_index)%cutoff)

      ! Set number of particles
      kim_data%neigh%n_part = natms

      ! Initialise neighbour list and number of neighbours
      kim_list = 0
      n_neigh = 0

      ! Build KIM neighbour list
      Do ipart = 1, natms

        Do neighbour = 1, list(0,ipart)
          ! Determine the 'neighbour'th neighbour of particles ipart
          jpart = list(neighbour,ipart)

          ! Add this neighbour to the KIM neighbour list of ipart
          n_neigh(ipart) = n_neigh(ipart) + 1
          kim_list(n_neigh(ipart),ipart) = jpart

          ! Add symmetric entry as KIM requires a full list. If the padding hint
          ! is true, then this is only nessecary if jpart is not a halo atom.
          If (padding_hint .eqv. .false. .or. jpart <= natms) Then
            n_neigh(jpart) = n_neigh(jpart) + 1
            kim_list(n_neigh(jpart),jpart) = ipart
          End If
        End Do
      End Do

      ! Determine padding neighbours of padding particles if required. The
      ! non-padding neighbours of padding particles have been added to the
      ! neighbour list above.
      If (padding_hint .eqv. .false.) Then
        Do ipart = natms+1, nlast-1
          Do jpart = ipart+1, nlast

            ! Add to list if the pair ipart-jpart are within the cutoff
            dr(1:3) = kim_data%coords(1:3,jpart) - kim_data%coords(1:3,ipart)
            r = Sqrt(Dot_product(dr,dr))

            If (r < cutoff) Then
              n_neigh(ipart) = n_neigh(ipart) + 1
              kim_list(n_neigh(ipart),ipart) = jpart

              ! Add symmetric entry
              n_neigh(jpart) = n_neigh(jpart) + 1
              kim_list(n_neigh(jpart),jpart) = ipart
            End If

          End Do
        End Do
      End If
    End Associate
#endif
  End Subroutine kim_neighbour_list

  !> KIM does not assume that the entirity of an atoms force will be calculated
  !> on the domain to which it belongs. It is therefore necessary to send the
  !> partial forces calculated on padding atoms to the appropriate domain.
  Subroutine kim_share_halo_forces(kim_data,parts,natms,nlast,lsi,lsa,ltg,map, &
      comm)
    !> KIM data type
    Type(kim_type), Target, Intent(InOut) :: kim_data
    !> Particles
    Type(corepart), Dimension(:), Intent(InOut) :: parts
    !> Number of particles in this domain (excluding the halo)
    Integer(Kind=wi), Intent(In) :: natms
    !> Number of particles in this domain and it's halo
    Integer(Kind=wi), Intent(In) :: nlast
    !> Some sort of local to global arrays
    Integer(Kind=wi), Intent(In) :: lsi(:),lsa(:)
    !> Local to global id
    Integer(Kind=wi), Intent(In) :: ltg(:)
    !> Map of neighbouring domain ids
    Integer(Kind=wi), Intent(In) :: map(1:26)
    !> Communicator type
    Type(comms_type), Intent(InOut) :: comm

    !> Direction to send and receive in
    Integer(Kind=wi) :: direction
    !> id of the node to send and receive from respectively
    Integer(Kind=wi) :: destination,source
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
    Do direction=1,6

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

      Associate(buffer => kim_data%kcomms%buffer)
        ! Pack data to send
        buf=0
        ! Loop over atoms to send
        Do atom=kim_data%kcomms%first(direction),kim_data%kcomms%last(direction)
          ! Pack forces
          buffer(buf+1)=Real(kim_data%forces(1,atom), wp)
          buffer(buf+2)=Real(kim_data%forces(2,atom), wp)
          buffer(buf+3)=Real(kim_data%forces(3,atom), wp)

          ! Pack global atom id
          buffer(buf+4)=Real(ltg(atom), wp)

          buf=buf+span
        End Do

        ! Exchange buffers
        If (comm%mxnode > 1) Then
          ! Determine sizes of buffers to send and receive
          send_size = buf
          recv_size = kim_data%kcomms%n_recv(direction)*span
          recv_start = kim_data%kcomms%recv_start
          If (recv_size > 0) Then
            Call girecv(comm,buffer(recv_start+1:recv_start+recv_size),source,Export_tag)
            Call gsend(comm,buffer(1:send_size),destination,Export_tag)
            Call gwait(comm)
          End  If
        Else
          recv_size = send_size
          recv_start = kim_data%kcomms%recv_start
        End If

        ! Unpack data received
        ! Determine where received data begins in the buffer
        buf = Merge(recv_start,0,comm%mxnode > 1)
        ! Loop number over atoms to receive
        Do i = 1, recv_size/span
          ! Find local index of received atom
          atom = local_index(Nint(buffer(buf+4)),nlast,lsi,lsa)

          ! Add partial force to local force arrays
          parts(atom)%fxx = parts(atom)%fxx + buffer(buf+1)
          parts(atom)%fyy = parts(atom)%fyy + buffer(buf+2)
          parts(atom)%fzz = parts(atom)%fzz + buffer(buf+3)

          buf = buf + span
        End Do
      End Associate
    End Do
#endif
  End Subroutine kim_share_halo_forces

  !> Neighbour list access function (called by KIM API)
  Subroutine get_neigh(data_object,n_lists,cutoffs,neighbour_list_index, &
      request,n_neigh,first,kerror) bind(c)
    !> KIM neighbour type address
    Type(c_ptr), Value, Intent(In) :: data_object
    !> Number of potential cutoffs in model
    Integer(c_int), Value, Intent(In) :: n_lists
    !> Array of cutoffs
    Real(c_double), Intent(In) :: cutoffs(n_lists)
    !> Which neighbour list to inspect
    Integer(c_int), Value, Intent(In) :: neighbour_list_index
    !> Number of the particle for which list information is requested
    Integer(c_int), Value, Intent(In) :: request
    !> Number of neighbours of the requested particle
    Integer(c_int), Intent(out) :: n_neigh
    !> Location of the first neighbour of the returned neighbour list
    Type(c_ptr), Intent(out) :: first
    !> Error status
    Integer(c_int), Intent(out) :: kerror

    Type(kim_neighbour_list_pointer_type), Pointer :: neigh(:)
    Integer(c_int), Pointer :: kim_neigh(:)
    Integer(c_int), Pointer :: kim_list(:,:)
    Real(c_double), Pointer :: kim_cutoff
    Integer(c_int), Pointer :: kim_n_part
    Character(Len=80) :: message

#ifdef KIM
    ! Assign local KIM neighbour type pointer to address of C data object
    Call c_f_pointer(data_object, neigh, [n_lists])
    Call c_f_pointer(neigh(neighbour_list_index)%n_neigh, kim_neigh, &
      [neigh(neighbour_list_index)%max_atoms])
    Call c_f_pointer(neigh(neighbour_list_index)%neigh_list, kim_list, &
      [neigh(neighbour_list_index)%max_neigh,neigh(neighbour_list_index)%max_atoms])
    Call c_f_pointer(neigh(neighbour_list_index)%cutoff, kim_cutoff)
    Call c_f_pointer(neigh(neighbour_list_index)%n_part, kim_n_part)

    ! Ensure neighbour list cutoff is larger than KIM model cutoff
    If (cutoffs(neighbour_list_index) > kim_cutoff) Then
      Call kim_warning('Neighbour list cutoff too small for model cutoff', __LINE__)
      kerror = 1
      Return
    End If

    ! Ensure requested particle exists
    If ( (request > kim_n_part) .or. (request < 1) ) Then
      Write(message,'(a,g0)') 'Invalid particle number requested in get_neigh: ', request
      Call kim_warning(message, __LINE__)
      kerror = 1
      Return
    End If

    ! Set number of neighbours for requested particle
    n_neigh = kim_neigh(request)
    ! Set the location for beginning of the returned neighbour list
    first = c_loc(kim_list(1,request))

    kerror = 0
#else
    n_neigh = 0_c_int
    first = c_null_ptr
    kerror = 0_c_int
#endif
  End Subroutine get_neigh

  !> Initialise the KIM neighbour list type
  Subroutine kim_neighbour_list_type_init(T,max_atoms,max_list)
  Class(kim_neighbour_list_type) :: T
    !> Extent of particle arrays. When the padding hint is true this is the
    !> maximum number of atoms on a node excluding the halo. When the hint is
    !> false this is the maximum number of atoms on a node including the halo.
    Integer(Kind=wi), Intent(In) :: max_atoms
    !> Extent of neighbour list arrays
    Integer(Kind=wi), Intent(In) :: max_list

    ! Allocate KIM neighbour list type arrays
    Allocate(T%n_neigh(max_atoms))
    Allocate(T%neigh_list(max_list,max_atoms))
  End Subroutine kim_neighbour_list_type_init

  !> Deallocate memory used by the KIM neighbour list type
  Subroutine kim_neighbour_list_type_cleanup(T)
    Type(kim_neighbour_list_type) :: T

    If (Allocated(T%n_neigh)) Then
      Deallocate(T%n_neigh)
    End If
    If (Allocated(T%neigh_list)) Then
      Deallocate(T%neigh_list)
    End If
  End Subroutine kim_neighbour_list_type_cleanup

  !> Initialise the KIM neighbour list pointer type
  Subroutine kim_neighbour_list_pointer_type_init(neigh_pointer,neigh,max_atoms, &
      max_list)
    Type(kim_neighbour_list_pointer_type), Intent(InOut) :: neigh_pointer
    Type(kim_neighbour_list_type), Target, Intent(In) :: neigh
    Integer(Kind=wi), Intent(In) :: max_atoms
    Integer(Kind=wi), Intent(In) :: max_list

    neigh_pointer%cutoff = C_loc(neigh%cutoff)
    neigh_pointer%n_part = C_loc(neigh%n_part)
    neigh_pointer%neigh_list = C_loc(neigh%neigh_list)
    neigh_pointer%n_neigh = C_loc(neigh%n_neigh)
    neigh_pointer%max_atoms = max_atoms
    neigh_pointer%max_neigh = max_list
  End Subroutine kim_neighbour_list_pointer_type_init

  !> Initialise the KIM comms type
  Subroutine kim_comms_type_init(T,buffer_size,mxnode)
  Class(kim_comms_type) :: T
    !> Size of the buffer to allocate
    Integer(Kind=wi), Intent(In) :: buffer_size
    !> Number of nodes
    Integer(Kind=wi), Intent(In) :: mxnode

    ! Allocate buffer
    Allocate(T%buffer(buffer_size))

    ! Set begining of receive portion
    T%recv_start = buffer_size/Merge(2,1,mxnode > 1)
  End Subroutine kim_comms_type_init

  !> Set the ids of atoms to be sent in kim_share_halo_forces
  Subroutine kim_comms_type_set(T,direction,n_recv,first,last)
  Class(kim_comms_type) :: T
    !> The direction in which data is being sent and received
    Integer(Kind=wi), Intent(In) :: direction
    !> Number of atoms to receive in each direction
    Integer(Kind=wi), Intent(In) :: n_recv
    !> ids of first and last atoms to send respectively in each direction
    Integer(Kind=wi) :: first, last

    T%n_recv(direction) = n_recv
    T%first(direction) = first
    T%last(direction) = last
  End Subroutine kim_comms_type_set

  !> Deallocate memory from KIM comms type
  Subroutine kim_comms_type_cleanup(T)
    Type(kim_comms_type) :: T

    If (Allocated(T%buffer)) Then
      Deallocate(T%buffer)
    End If
  End Subroutine kim_comms_type_cleanup

  !> Allocate memory for the KIM type
  Subroutine kim_type_init(T,mxatms)
  Class(kim_type) :: T
    !> Extent of particle arrays including halo particles (coordinates, forces, etc.)
    Integer(Kind=wi), Intent(In) :: mxatms

#ifdef KIM
    Allocate(T%species_name(mxatms))
    Allocate(T%species_code(mxatms))

    Allocate(T%contributing(mxatms))

    Allocate(T%coords(3,mxatms))
    Allocate(T%forces(3,mxatms))
#endif
  End Subroutine kim_type_init

  !> Deallocate memory from KIM type
  Subroutine kim_type_cleanup(T)
    Type(kim_type) :: T

    Integer(Kind=c_int) :: kerror

#ifdef KIM
    If (Allocated(T%species_name)) Then
      Deallocate(T%species_name)
    End If
    If (Allocated(T%species_code)) Then
      Deallocate(T%species_code)
    End If

    If (Allocated(T%forces)) Then
      Deallocate(T%forces)
    End If
    If (Allocated(T%coords)) Then
      Deallocate(T%coords)
    End If

    If (Allocated(T%neigh)) Then
      Deallocate(T%neigh)
    End If
    If (Allocated(T%neigh_pointer)) Then
      Deallocate(T%neigh_pointer)
    End If

    Call kim_compute_arguments_destroy(T%model_handle, &
      T%compute_arguments_handle, kerror)

    Call kim_model_destroy(T%model_handle)
#endif
  End Subroutine kim_type_cleanup
End Module kim
