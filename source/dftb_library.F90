Module dftb_library

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for calling DFTB+ 19.2 external library 
  !
  ! copyright - daresbury laboratory
  ! author    - A. Buccheri Nov 2018 - Jan 2020
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef WITH_DFTBP 
  Use, Intrinsic :: iso_fortran_env, Only: output_unit
  !> dl_poly modules
  Use comms,           Only : wp, root_id
  Use configuration,   Only : configuration_type, reallocate, unpack_gathered_coordinates, &
       ordering_indices
  Use errors_warnings, Only : error, warning
  Use asserts,         Only : assert
  Use parse,           Only : lower_case, upper_case, number_of_lines 

  !> DFTB+ API modules 
  Use dftb_mmapi,      Only : TDftbPlus_init, TDftbPlus, TDftbPlusInput, TDftbPlus_destruct!, string

  Implicit none

  Private
  
  !> Type to hold data gathered from geometry DL_POLY processes
  !> for use in DFTB+ 19.2 library
  !
  Type, Public :: dftb_geometry_type
     private 
     !> 'Periodic system' logical
     Logical               :: periodic
     !> 'Coordinates in fractional coordinates' logical
     Logical               :: fractCoord
     !> Atomic coordinates for all atoms in system (Bohr)
     Real(wp), Allocatable :: coords(:,:)
     !> Lattice vectors, stored columnwise 
     Real(wp)              :: lattice(3,3)
     !> Index of species type for each atom in the system  
     Integer,  Allocatable :: species(:)
     !> Map of species index to species label:
     !    speciesName(species(iatom)) = 'element_name'
     !  A water molecule for example:
     !    speciesNames = (/'H','O'/)
     !    species = (/1,2,2/)
     Character(50), Allocatable :: speciesNames(:)

   contains
     procedure :: initialise => initialise_dftb_geometry_type
     procedure :: finalise => finalise_dftb_geometry_type
     procedure :: set_geometry => set_dftb_geometry_type
   
  End type dftb_geometry_type

  !> Unit conversions
  !> E0 to J (== 10 J/mol) ref pg 7 DL POLY Manual 
  Real(wp), Public, Parameter :: EMolecular__J = 1.6605402e-23_wp
  !> Bohr to Angstrom
  Real(wp), Public, Parameter :: Bohr__AA = 0.529177249_wp
  !> Angstrom to Bohr
  Real(wp), Public, Parameter :: AA__Bohr = 1.0_wp / Bohr__AA
  !> Hartree to Joule
  Real(wp), Public, Parameter :: Hartree__J = 4.3597441775e-18_wp

  !> Modular parameters
  Logical :: IO
  !> Default dftb input file name 
  Character(Len=11), Parameter :: dftb_fname = 'dftb_in.hsd'
  !> Suppress IO 
  Logical,   Parameter :: suppressIO = .True.
  
  Public :: run_dftbplus, dftb_geometry_type, convert_unit
  
Contains
  
  !> @brief Perform a DFTB calculation with DFTB+ 19.2
  !!
  !! @param[in]  MDstatus           MD object containing initial, current and final steps
  !! @param[in]  DFTBglobalMpiComm  MPI communicator information in object of DFTB+ type 
  !! @param[inout] geo              Object containing all geometry data
  !! @param[inout] forces           Forces (in units of Ha/B)
  !
  Subroutine run_dftbplus(comm, flow, geo, forces)

    !Arguments 
    Type(comms_type), Intent( In ) :: comm
    Type(flow_type),  Intent( In ) :: flow
    Type(dftb_geometry_type), Intent( In) :: geo
    Real(wp), Allocatable, Intent(InOut)  :: forces(:,:)

#ifdef WITH_DFTBP
    !Local data
    Type(TDftbPlus), Save :: dftb
    Integer,  Allocatable, Save :: prior_global_ordering(:)
    Real(wp), Allocatable, Save :: atomic_charges(:)
    Integer,  Allocatable :: map_prior_to_current(:), map_current_to_prior(:)
    Logical  :: atomic_ordering_changed = .false.
    Real(wp) :: merminEnergy
    
    Call assert(Allocated(forces), &
         'Forces must be allocated upon calling run_dftbplus')
    
    If(flow%step == initial_MD_step) Then
       IO = comm%idnode == root_id 
       Call TDftbPlus_init(dftb, output_unit, comm%comm)
       Call initialise_dftbplus_tree(comm, geo, hsd_tree)
       Call dftb%setupCalculator(hsd_tree)
       Allocate(atomic_charges(config%megatm))
       !TODO(Alex) Implement this correctly
       !Should be a set of unique indices for all atoms in the system
       !Probably need to gather or something 
       Allocate(prior_global_ordering(config%megatm))
       prior_global_ordering = config%ltg
    EndIf
      
    !Note, calculation fails if lattice is not passed (DFTB+ issue) 
    Call dftb%setGeometry(geo%coords, geo%lattice)

    atomic_ordering_changed = All(prior_global_ordering = config%ltg)
    If(atomic_ordering_changed) Then
       Call dftb%setSpeciesAndDependents(geo%speciesNames, geo%species)
       Call ordering_indices(prior_global_ordering, config%ltg, &
            map_current_to_prior, map_prior_to_current)
       atomic_charges = atomic_charges(map_prior_to_current)
       Call dftb%setAtomicCharges(atomic_charges)
    Endif

    Call dftb%getEnergy(merminEnergy)  
    Call dftb%getGradients(forces)
    forces = -1._wp * forces
    Call dftb%getGrossCharges(atomic_charges)
 
    If(flow%step == flow%run_steps) Then
       Call TDftbPlus_destruct(dftb)
       Deallocate(global_ordering)
       Deallocate(atomic_charges)
    Endif
    
#endif
  End Subroutine run_dftbplus

  
  !> @brief Initialise DFTB tree
  !!
  !! Functions to modify nodes rather than overwrite not exposed in DFTB+ API
  !! As such, it's the user's job to correctly set everything not in the geometry block     
  !!
  !! @param[in]   comm        MPI communicator object
  !! @param[in]   geo         Geometry for DFTB+
  !! @param[out]  hsd_tree    DFTB+ Input tree
  !! @param[in]   input_fname Optional DFTB+ input file name 
  !
  Subroutine initialise_dftbplus_tree(comm, geo, hsd_tree, input_fname)
    !Arguments 
    Type(comms_type),         Intent( InOut ) :: comm
    Type(dftb_geometry_type), Intent( In)     :: geo 
    Type(TDftbPlusInput),     Intent(   Out ) :: hsd_tree
    Character(Len=*),         Intent( In    ),   Optional :: input_fname

    !Local data
    !> DFTB+ input file name 
    Character(Len=*) :: fname
    !> Pointers to the parts of the input tree that will be set                                
    Type(fnode), Pointer :: pRoot, pGeo, pOptions, pParserOpts, pDftb, pKpoints
    !> "Does geometry already exist in DTFB+ input?" (== "replace geometry in HSD tree?") 
    Logical  :: replace_geometry
    Character(Len=100) :: error_message
      
    If(Present(input_fname)) Then
       fname = input_fname
    Else
       fname = dftb_fname
    EndIf

    replace_geometry = check_dftb_input_file_for_geom(fname)
    Call check_dftb_input_file_for_latticeopt(fname)
    
    !Read dftb input file into hsd_tree
    Call dftb%getInputFromFile(fname, hsd_tree)

    !Set structure data, retain rest
    Call hsd_tree%getRootNode(pRoot)
    Call setChild(pRoot, "Geometry", pGeo, replace=replace_geometry)
    Call setChildValue(pGeo, "Periodic", geo%periodic)
    Call setChildValue(pGeo, "LatticeVectors", geo%lattice)
    Call setChildValue(pGeo, "TypeNames", geo%speciesNames)

    !See DFTB+ subroutine in lib_type/typegeometryhsd.F90
    Call setChildValue(pGeo, "TypesAndCoordinates", &
         Reshape(geo%species, (/ 1, Size(geo%species) /)), geo%coords)

    !Always compute forces
    Call setChild(pRoot, "Analysis", pAnalysis, replace=.True.)
    Call setChildValue(pAnalysis, "CalculateForces", .True.)
    
    If(suppressIO) Then
       Call setChildValue(pAnalysis, "WriteBandOut", .False.)
       Call setChildValue(pAnalysis, "MullikenAnalysis", .False.)
       
       Call setChild(pRoot, "ParserOptions", pParserOpts, replace=.True.)
       Call setChildValue(pParserOpts, "ParserVersion", 8)
       Call setChildValue(pParserOpts, "WriteDetailedOut", .False.)
       Call setChildValue(pParserOpts, "WriteXMLInput", .False.)
       Call setChildValue(pParserOpts, "WriteHSDInput", .False.)
       Write(*,*) 'Suppressed file IO from DFTB+'
    Endif
      
  End Subroutine initialise_dftbplus_tree
    
  
  !> @brief Check if DFTB+ input file contains "Geometry" node
  !!
  !! Checks that "Geometry" node is not present in dftb_in.hsd
  !! as this is provided by DL_POLY
  !
  Function check_dftb_input_file_for_geom(input_fname) Result(geometry_present)

    Character(Len=*), Intent(In), Optional :: input_fname
    Logical           :: geometry_present
  
    Character(Len=*)  :: fname
    Character(Len=60) :: line
    Integer           :: N,i,noccurrences,ios
    Logical           :: exist
    Character(Len=100):: message 

    If(Present(input_fname)) Then
       fname = input_fname
    Else
       fname = dftb_fname
    Endif

    If(IO) Write(*,*) 'Checking '//Trim(fname)
    noccurrences=0

    Inquire(file=Trim(Adjustl(fname)), exist=exist)

    If(exist) Then
       N=number_of_lines(fname)
       Open(unit=100, file=Trim(Adjustl(fname)), form='formatted', status='old', iostat=ios)
       Do i=1,N
          Read(100, '(A)') line
          !Catch uncommented instances of geometry in fname
          If(Index(lower_case(line), 'geometry') >0 .and. &
             Index(lower_case(line), 'genformat')>0 .and. &
             line(1:1)/='#' ) Then
             noccurrences = noccurrences + 1
          Endif
       Enddo
       Close(100)
       
    Elseif(.Not. exist) Then
       message = 'Error openning: '//Trim(Adjustl(fname))
       Call error(0, message=Trim(Adjustl(message)), master_only=.true.)
    Endif
    
    If(noccurrences > 0) Then
       message = 'Geometry is present in '//Trim(Adjustl(fname))//NEW_LINE('A')//&
            ' This will be overwritten with data from CONFIG'
       Call warning(Trim(Adjustl(message)), master_only=.true.) 
       geometry_present = .True.
    Else
       geometry_present = .False.
    Endif
  
  End Function check_dftb_input_file_for_geom

 
  !> @brief Check if dftb_in.hsd contains lattice relaxation option 
  !!
  !! Checks that 'LatticeOpt = Yes' is not present in dftb_in.hsd
  !! If found, an error is issued and the interface stops
  !!
  !
  Subroutine check_dftb_input_file_for_latticeopt(input_fname)

    !Declarations
    Character(Len=*), Optional, Intent(In) :: input_fname
    Character(Len=20) :: fname 
    Character(Len=60) :: line
    Character(Len=200):: message
    Integer           :: N,i,noccurrences,ios
    Logical           :: exist

    !Main Routine
    If(Present(input_fname)) Then
       fname = input_fname
    Else
       fname = dftb_fname
    Endif
    
    noccurrences=0

    Inquire(file=Trim(Adjustl(fname)), exist=exist)

    If(exist) Then
       N=number_of_lines(fname)
       Open(unit=100, file=Trim(Adjustl(fname)), form='formatted', status='old', iostat=ios)
       Do i = 1,N
          Read(100, '(A)') line
          If(Index(lower_case(line), 'latticeopt') >0 .And. &
             Index(lower_case(line), 'yes')>0 .And. &
             line(1:1)/='#' ) Then
             noccurrences = noccurrences + 1
          Endif
       Enddo
       Close(100)
       
    Elseif(.Not. exist) Then
       message = 'Error openning: '//Trim(Adjustl(fname))
       Call error(0, message=Trim(Adjustl(message)), master_only=.true.)
    Endif

    If(noccurrences > 0) Then
       message = 'LatticeOpt set to yes in '//Trim(Adjustl(fname))//'.'//NEW_LINE('A')//&
            'Lattice optimisation should not be used in conjunction with dl_poly'
       Call error(0, message=Trim(Adjustl(message)), master_only=.true.)
    Endif
  
  End Subroutine check_dftb_input_file_for_latticeopt 

  
  !> @brief Assign DFTB+ boundary condition given DLPOLY setting 
  !!
  !! @param[in]      config           Object containing geometry and forces
  !! @param[inout]   geo              Geometry object
  !
  Subroutine boundary_option(geo, config)
    
    !Subroutine arguments
    Type(configuration_type), Intent(In)    :: config
    Type(dftb_geometry_type), Intent(InOut) :: geo 

    !Local data
    Character(Len=1)  :: boundary
    Character(Len=50) :: error_message

    !(0) = finite boundary
    boundary = 'c'
    !(1,2,3) = cubic, orthorhombic, parallelpiped
    !(6) = xy periodic, z finite  
    If(config%imcon > 0) boundary = 's'

    Select Case (boundary)
    !Periodic (s)upercell   
    Case("S","s") 
       geo%periodic  = .true.
       geo%fracCoord = .false.
    !Periodic supercell. Atoms defined with (f)ractional coordinates 
    Case("F","f") 
       geo%periodic  = .true.
       geo%fracCoord = .true.
    !Finite (c)luster 
    Case("C", "c")
       geo%periodic  = .false.
       geo%fracCoord = .false.
    Case default
       error_message = 'Invalid DFTB+ boundary condition option used: '//boundary
       Call error(0,message=error_message,master_only=.true.)
    End Select
    
  End Subroutine boundary_option

  !> @brief Find the name of each species in the system 
  !!
  !! For example, water would give unique_species_names = (/'H','O'/)
  !! 
  !! @param[in]     atmnam                 Name of each atom in system
  !! @param[inout]  unique_species_names   Name of each species in system 
  !
  Function find_unqiue_species(atmnam), Result(unique_species_names)

    !Arguments
    Character(Len=len_atmnam), Allocatable, Intent(In)    :: atmnam(:)
    Character(Len=3),          Allocatable  :: unique_species_names(:)

    !Local data
    Integer           :: ia,ja,cnt,natms,status 
    Logical           :: already_found

    Allocate(unique_species_names(118))
    natms=Size(atmnam)
    unique_species_names(1)=atmnam(1)
    cnt=1
    already_found=.false.
    
    Do ja=1,natms 
       Do ia=1,cnt
          If( Trim(Adjustl( unique_species_names(ia) )) == &     
              Trim(Adjustl( atmnam(ja) )) )Then
             already_found=.true.
          End If
       End Do
       If(already_found.eqv..false.)Then
          cnt=cnt+1
          unique_species_names(cnt)=atmnam(ja)
       End If
       already_found=.false.
    End Do

    Call reallocate(cnt-Size(unique_species_names), unique_species_names, status)

  End Function find_unqiue_species


  !> Initialise dftb geometry 
  Subroutine initialise_dftb_geometry_type(this, config, all_atom_names)
    Class(dftb_geometry_type),     Intent(InOut) :: this
    Type(configuration_type),      Intent(In)    :: config
    Character(Len=8), Allocatable, Intent(In)    :: all_atom_names(:)  

    Allocate(this%species(config%megatm))
    Allocate(this%coords(3,config%megatm))
    this%species=0
    this%coords=0._wp
    Call boundary_option(this,config)
    !Set labels for unique atomic species in the simulation
    ! - should never change even if all species aren't present at all times
    !TODO(Alex) Don't know if I can assign an allocatable array in this manner
    this%speciesNames = find_unqiue_species(all_atom_names) 
  
  End Subroutine initialise_dftb_geometry_type

  
  !> Finalise dftb geometry 
  Subroutine finalise_dftb_geometry_type(this)
    Class(dftb_geometry_type),     Intent(InOut) :: this
    
    If(Allocated(this%species)) Deallocate(this%species)
    If(Allocated(this%coords))  Deallocate(this%coords)
    If(Allocated(this%speciesNames)) Deallocate(this%speciesNames)
    this%lattice = 0._wp
    this%periodic = .false.
    this%fracCoord = .false.
    
  End Subroutine finalise_dftb_geometry_type

  
  !> @brief Fill DFTB geometry object 
  !!
  !! @param[inout]  this             DFTB+ geometry data
  !! @param[in]     comm             MPI communicator and data
  !! @param[in]     config           Configuration data: config.cell = lattice vectors 
  !! @param[in]     all_coords       1D real array containing atomic coordinates for all atoms
  !! @param[in]     MPI_buffer       MPI distribution arrays used by gather/scatter for all_coords
  !! @param[in]     all_atom_names   Name of each atom in system
  !
  Subroutine set_dftb_geometry_type(this, comm, config, gathered, all_atom_names)
   
    !Subroutine arguments
    Class(dftb_geometry_type),     Intent(InOut) :: this
    Type(comms_type),              Intent(In)    :: comm
    Type(configuration_type),      Intent(In)    :: config
    Type(coordinate_buffer_type),  Intent(In)    :: gathered
    Character(Len=8), Allocatable, Intent(In)    :: all_atom_names(:)

    !Local data
    Integer            :: ia,ja,cnt,natms_local,ip,ix,iy,iz
    Character(Len=100) :: message, rank_label 

    Call assert(Allocated(this%coords), message='coords in dftb_geometry not allocated')
    Call assert(Allocated(this%species),message='species in dftb_geometry not allocated')
    Call assert(Size(this%coords,1) == 3, message='size(coords,1) /= 3')
    Call assert(Size(this%coords,2) == config%megatm, message='size(coords,2) /= megatm')
    message = 'Size(all_atom_names) /= Size(this%speciesNames) in set_dftb_geometry'
    Call assert(Size(all_atom_names) == Size(this%speciesNames), message)

    !Pass packed array of atomic coordinates to geometry object
    cnt=0
    Do ip=1,comm%mxnode
       natms_local=int(gathered%mpi%counts(ip)/3)
       Write(rank_label, '(I6)') ip-1 
       message = 'Local number of atoms on rank '//Trim(Adjustl(rank_label))//&
            ' , config%natms /= int(gathered%mpi%counts(rank+1)/3)' 
       Call assert(natms_local == config%natms, message)
       Do ia=1,natms_local
          cnt=cnt+1
          !mpi_process_offset + xyz offset + local_index
          ix = gathered%mpi%displ(ip) +                  ia
          iy = gathered%mpi%displ(ip) + natms_local   +  ia
          iz = gathered%mpi%displ(ip) + natms_local*2 +  ia
          this%coords(1:3,cnt)= &
               (/gathered%coords(ix),gathered%coords(iy),gathered%coords(iz)/) * AA__Bohr
       End Do
    End Do

    !TODO(Alex) Could use subroutine in configuration module
    !Call unpack_gathered_coordinates(comm, config, gathered, this%coords)
    !this%coords(:,:) = this%coords(:,:) * AA__Bohr

    !Assign unique integer to each species
    Do ia=1,config%megatm
       Do ja=1,Size(this%speciesNames)
          If(Trim(Adjustl(all_atom_names(ia))) ==  &
             Trim(Adjustl(this%speciesNames(ja))) )Then
             this%species(ia)=ja
             Exit
          End If
       End Do
    End Do
   
    If(this%periodic)Then
       !Stored columnwise
       this%lattice(1:3,1) = config%cell(1:3) * AA__Bohr 
       this%lattice(1:3,2) = config%cell(4:6) * AA__Bohr 
       this%lattice(1:3,3) = config%cell(7:9) * AA__Bohr        
    End If
    
  End Subroutine set_dftb_geometry_type

  !> @brief Convert between internal units of DLPOLY and DFTB+, or vice versa
  !!
  !! DFTB+ internal units are 'Hartree atomic units'
  !! DL POLY internal units given on page 7 of the manual (4.08)
  !! Valid quantities for conversion: length, force 
  !! 
  !! @param[in]   from     Character: Program name of unit to convert from
  !! @param[in]   to       Character: Program name of unit to convert to
  !! @param[in]   unit     Character: Physical quantity to convert 
  !! @param[out]  factor   Conversion factor
  !
  Function convert_unit(from,to,unit) Result(factor)
    Character(Len=*),  Intent(In)  :: from,to
    Character(Len=*),  Intent(In)  :: unit
    Real(wp)                       :: factor
    Character(Len=Len(from))       :: from_local
    Character(Len=Len(to))         :: to_local
    Character(Len=Len(unit))       :: unit_local
    Character(Len=50)              :: warning_message

    factor = 0._wp

    !Enforce case consistency
    from_local= Trim(Adjustl(uppercase(from)))
    to_local  = Trim(Adjustl(uppercase(to)))
    unit_local= Trim(Adjustl(lowercase(unit)))

    !No unit conversion required
    If( from_local == to_local ) Then
       factor = 1._wp
       Return
    End If
       
    !DLPOLY to DFTB+
    If(from_local=='DLPOLY' .and. to_local=='DFTB+') Then
       If(unit_local=='length' ) Then
          factor = AA__Bohr
          Return
       Elseif( unit_local=='force' ) Then
          factor = (EMolecular__J/Hartree__J)*Bohr__AA
          Return
       Else
          warning_message = 'Not able to convert between internal units for ('//&
               unit_local//')'//new_line('A')//' Conversion factor will be zeroed'
           Call warning(warning_message,master_only=.true.)
          Return
       End If
    End If

    !DFTB+ to DLPOLY
    If(from_local=='DFTB+' .and. to_local=='DLPOLY') Then  
       If(unit_local=='length' )Then
          factor=Bohr__AA
          Return
       Elseif( unit_local=='force' ) Then
          factor = (Hartree__J/EMolecular__J)*AA__Bohr
          Return
       Else
          warning_message = 'Not able to convert between internal units for ('//&
               unit_local//')'//New_line('A')//' Conversion factor will be zeroed'
          Call warning(warning_message, master_only=.true.)
          Return
       End If
    End If

    !Erroneous inputs for to/from
    warning_message = 'Options not recognised: '//from_local//' and/or '//to_local//&
         New_line('A')//' Conversion factor will be zeroed'
    Call warning(warning_message,master_only=.true.)
    
  End Function convert_unit

#endif
End Module dftb_library
