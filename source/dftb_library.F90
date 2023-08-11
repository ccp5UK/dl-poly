Module dftb_library

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for calling DFTB+ 19.2 external library 
  !
  ! copyright - daresbury laboratory
  ! author    - A. Buccheri Nov 2018 - Mar 2020
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DFTBP 
  Use, Intrinsic :: iso_fortran_env, Only: output_unit
  !> dl_poly modules
  Use comms,           Only : comms_type, wp, root_id
  Use configuration,   Only : configuration_type, reallocate, unpack_gathered_coordinates, &
                              ordering_indices, len_atmnam, coordinate_buffer_type, gather_forces
  Use errors_warnings, Only : error, warning
  Use asserts,         Only : assert
  Use parse,           Only : number_of_lines 
  Use flow_control,    Only : flow_type
  Use kinds,           Only : STR_FILENAME
  
  !> DFTB+ API modules 
  Use dftbp_mmapi,      Only : TDftbPlus_init, TDftbPlus, TDftbPlusInput, TDftbPlus_destruct!, string
  Use dftbp_hsdapi,     Only : fnode, getChild, getChildren, setChild, getChildValue, &
                               setChildValue, dumpHsd
  Implicit none

  Private
  
  !> Type to hold data gathered from geometry DL_POLY processes
  !> for use in DFTB+ 19.2 library
  !
  Type dftb_geometry_type
     !> 'Periodic system' logical
     Logical               :: periodic
     !> 'Coordinates in fractional coordinates' logical
     Logical               :: fracCoord
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
  Character(Len=STR_FILENAME), Parameter :: dftb_fname = 'dftb_in.hsd'
  !> Suppress IO - Change for debugging
  Logical,   Parameter :: suppressIO = .true.

  Interface output_dftb_forces
     Module Procedure output_dftb_forces_from_array, output_dftb_forces_from_config
  End Interface output_dftb_forces
  
  Public :: run_dftbplus, dftb_geometry_type, convert_unit, print_DFTB_geometry_data,&
            output_dftb_forces, print_DFTB_geo_in_xyz
  
Contains
  
  !> @brief Perform a DFTB calculation with DFTB+ 19.2
  !!
  !! @param[inout] comm            Object containing MPI communicator       
  !! @param[in]    flow            Object containing MD step/time data     
  !! @param[inout] geo             Object containing all geometry data
  !! @param[inout] forces          Forces (in units of Ha/B)
  !! @param[inout] atomic_charges  Atomic charges
  !! @param[in] run_app_test       Logical indicating whether or not dftb+ app
  !!                               test is running 
  !
  Subroutine run_dftbplus(comm, flow, geo, forces, atomic_charges, run_app_test)

    Type(comms_type),         Intent( InOut ) :: comm
    Type(flow_type),          Intent( In    ) :: flow
    Type(dftb_geometry_type), Intent( In    ) :: geo
    Real(wp), Allocatable,    Intent( InOut ) :: forces(:,:), atomic_charges(:)
    Logical,  Optional,       Intent( In    ) :: run_app_test
    
    Type(TDftbPlus), Save :: dftb
    Type(TDftbPlusInput)  :: hsd_tree
    Real(wp) :: merminEnergy
    Integer, Allocatable, Save :: prior_species(:)
    Logical :: always_set_species_and_dependents
    
    Call assert(Allocated(forces), &
         message = 'Forces must be allocated upon calling run_dftbplus')
    Call assert(Allocated(atomic_charges), &
         message = 'Atomic charges must be allocated upon calling run_dftbplus')

    If (Present(run_app_test)) Then
       always_set_species_and_dependents = run_app_test
    Else
       always_set_species_and_dependents = .false.
    Endif
    
    If(comm%idnode == root_id) Then
       Write(*,*) 'On MD step: ', flow%step
    Endif
    
    If(flow%step == flow%initial_MD_step) Then
       IO = comm%idnode == root_id
       Call TDftbPlus_init(dftb, output_unit, comm%comm)
       Call initialise_dftbplus_tree(geo, dftb, hsd_tree)
       !Dump hsd tree to fort.001
       !Call dumpHsd(hsd_tree%hsdTree, 001)
       Call dftb%setupCalculator(hsd_tree)
       Allocate(prior_species(Size(geo%species)))
       prior_species = geo%species
    EndIf

    Call dftb%setGeometry(geo%coords, geo%lattice)

    If(Any(prior_species /= geo%species) .or. &
       always_set_species_and_dependents) Then
       !write(*,*) 'Resetting Species and Dependents got called on MD step:',flow%step
       prior_species = geo%species
       Call dftb%setSpeciesAndDependents(geo%speciesNames, geo%species)
    Endif
    
    Call dftb%getEnergy(merminEnergy)  
    Call dftb%getGradients(forces)
    forces = -1._wp * forces
    Call dftb%getGrossCharges(atomic_charges)
    
    If(flow%step == flow%run_steps) Then
       Call TDftbPlus_destruct(dftb)
       Deallocate(prior_species)
    Endif
    
  End Subroutine run_dftbplus

  
  !> @brief Initialise DFTB tree
  !!
  !! Functions to modify nodes rather than overwrite not exposed in DFTB+ API
  !! As such, it's the user's job to correctly set everything not in the geometry block     
  !!
  !! @param[in]   geo         Geometry for DFTB+
  !! @param[inout]dftb        DFTB+ caculation object 
  !! @param[out]  hsd_tree    DFTB+ Input tree
  !! @param[in]   input_fname Optional DFTB+ input file name 
  !
  Subroutine initialise_dftbplus_tree(geo, dftb, hsd_tree, input_fname)

    Type(dftb_geometry_type), Intent( In    ) :: geo
    Type(TDftbPlus),          Intent( InOut ) :: dftb
    Type(TDftbPlusInput),     Intent(   Out ) :: hsd_tree
    Character(Len=*),         Intent( In    ),   Optional :: input_fname

    !> DFTB+ input file name 
    Character(Len=STR_FILENAME) :: fname
    !> Pointers to the parts of the input tree that will be set                                
    Type(fnode), Pointer        :: pRoot, pGeo, pParserOpts, pAnalysis
    !> "Does geometry already exist in DTFB+ input?" (== "replace geometry in HSD tree?") 
    Logical  :: replace_geometry
      
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
       Call setChildValue(pParserOpts, "ParserVersion",          8)
       Call setChildValue(pParserOpts, "WriteDetailedOut", .False.)
       Call setChildValue(pParserOpts, "WriteXMLInput",    .False.)
       Call setChildValue(pParserOpts, "WriteHSDInput",    .False.)
       Call setChildValue(pParserOpts, 'WriteResultsTag',  .False.)
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
  
    Character(Len=STR_FILENAME) :: fname
    Character(Len=60)           :: line
    Integer                     :: N,i,noccurrences,ios
    Logical                     :: exist
    Character(Len=100)          :: message 

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
  !
  Subroutine check_dftb_input_file_for_latticeopt(input_fname)

    !Declarations
    Character(Len=*), Optional, Intent(In) :: input_fname
    Character(Len=STR_FILENAME) :: fname 
    Character(Len=60)           :: line
    Character(Len=200)          :: message
    Integer                     :: N,i,noccurrences,ios
    Logical                     :: exist

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
  Subroutine find_unqiue_species(atmnam, unique_species_names) 

    !Arguments
    Character(Len=len_atmnam),      Intent(In)   :: atmnam(:)
    Character(Len=50), Allocatable, Intent(Out)  :: unique_species_names(:)

    !Local data
    Integer           :: ia,ja,cnt,status 
    Logical           :: already_found

    Allocate(unique_species_names(118))
    unique_species_names(1)=atmnam(1)
    cnt=1
    already_found=.false.
    
    Do ja=1,Size(atmnam)
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

  End Subroutine find_unqiue_species


  !> Initialise dftb geometry
  !!
  !! @param[inout]  this             DFTB+ geometry data
  !! @param[in]     config           Configuration data
  !! @param[in]     all_atom_names   List of all atom names in system
  !!                                 Size(megatm)
  !
  Subroutine initialise_dftb_geometry_type(this, config, all_atom_names)
    Class(dftb_geometry_type),     Intent(InOut) :: this
    Type(configuration_type),      Intent(In)    :: config
    Character(Len=len_atmnam),     Intent(In)    :: all_atom_names(:)  

    Allocate(this%species(config%megatm))
    Allocate(this%coords(3,config%megatm))
    this%species=0
    this%coords=0._wp
    Call boundary_option(this,config)
    !Set labels for unique atomic species in the simulation
    ! - should never change even if all species aren't present at all times
    Call find_unqiue_species(all_atom_names, this%speciesNames)
  
  End Subroutine initialise_dftb_geometry_type

  
  !> Finalise dftb geometry 
  !! @param[inout]  this             DFTB+ geometry data 
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
   
    Class(dftb_geometry_type),     Intent(InOut) :: this
    Type(comms_type),              Intent(In)    :: comm
    Type(configuration_type),      Intent(In)    :: config
    Type(coordinate_buffer_type),  Intent(In)    :: gathered
    Character(Len=8), Allocatable, Intent(In)    :: all_atom_names(:)

    Integer            :: ia,ja
    Character(Len=100) :: message, rank_label 

    Call assert(Allocated(this%coords), message='coords in dftb_geometry not allocated')
    Call assert(Allocated(this%species),message='species in dftb_geometry not allocated')
    Call assert(Size(this%coords,1) == 3, message='size(coords,1) /= 3')
    Call assert(Size(this%coords,2) == config%megatm, message='size(coords,2) /= megatm')
    Call assert(Size(all_atom_names) == Size(this%species), &
         message='Size(all_atom_names) /= Size(this%species) in set_dftb_geometry')

    Write(rank_label, '(I6)') comm%idnode 
    message = 'Local number of atoms on rank '//Trim(Adjustl(rank_label))//&
         ' , config%natms /= int(gathered%mpi%counts(rank+1)/3)' 
    Call assert(config%natms == int(gathered%mpi%counts(comm%idnode+1)/3), &
                message=message)

    !Pass packed array of atomic coordinates to geometry object
    Call unpack_gathered_coordinates(comm, config, gathered, this%coords)
    this%coords(:,:) = this%coords(:,:) * AA__Bohr

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
    Character(Len=50)              :: message

    factor = 0._wp

    !Enforce case consistency
    from_local= Trim(Adjustl(upper_case(from)))
    to_local  = Trim(Adjustl(upper_case(to)))
    unit_local= Trim(Adjustl(lower_case(unit)))

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
          message = 'Not able to convert between internal units for ('//&
               unit_local//')'
           Call error(0, message=trim(adjustl(message)), master_only=.true.)
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
          message = 'Not able to convert between internal units for ('//&
               unit_local//')'
          Call error(0, message=trim(adjustl(message)), master_only=.true.)
          Return
       End If
    End If

    !Erroneous inputs for to/from
    message = 'Options not recognised: '//from_local//' and/or '//to_local
    Call error(0, message=trim(adjustl(message)), master_only=.true.)
    
  End Function convert_unit

  
  !> @brief Output forces from all processes in DFTB+ units/format
  !!
  !! @param[inout] comm      MPI communicator 
  !! @param[in]    flow      MD step/time data 
  !! @param[in]    config    Configuration data  
  !
  Subroutine output_dftb_forces_from_config(comm, flow, config)
    Type(comms_type),         Intent( InOut ) :: comm
    Type(flow_type),          Intent( In    ) :: flow 
    Type(configuration_type), Intent( In    ) :: config
    Real(wp), Allocatable :: forces(:,:)
    Real(wp)              :: unit_factor
    
    Allocate(forces(3,config%megatm))
    Call gather_forces(comm, config, forces)
    unit_factor = convert_unit('DLPOLY','DFTB+','force')
    forces = forces * unit_factor
    Call output_dftb_forces_from_array(comm, flow, forces)
    
  End Subroutine output_dftb_forces_from_config
    

  !> @brief Output forces from all processes in DFTB+ units/format
  !!
  !! @param[inout] comm      MPI communicator  
  !! @param[in]    flow      MD step/time data
  !! @param[in]    forces    Forces on all atoms in system
  !!                         Expected output from DFTB+
  !
  Subroutine output_dftb_forces_from_array(comm, flow, forces)
    Type(comms_type),   Intent( InOut ) :: comm
    Type(flow_type),    Intent( In    ) :: flow 
    Real(wp),           Intent( In    ) :: forces(:,:)
    Integer             :: iatom
    Character(Len = 99) :: file_name
    Character(Len=2)    :: md_lab

    Call assert(Size(forces,1) == 3, message=&
         "forces shape /= (3,megatm)") 
    
    If(comm%idnode == root_id) Then
       Write(md_lab, '(I2)') flow%step
       file_name = 'dl_forces_'//trim(adjustl(md_lab))//'.dat'
       Open(001, file=Trim(Adjustl(file_name)))
       Write(001,'(A,I2)') "Total Forces, running with np=",comm%mxnode
       Do iatom = 1, Size(forces, 2)
          !Format settings taken from DFTB+, line 2968 in mainio.F90                               
          Write(001, "(I5, 3F20.12)") iatom, forces(:, iatom)
       Enddo
       Close(001)
    Endif

  End Subroutine output_dftb_forces_from_array


  !> @brief Print DFTB+ geometry data in xyz format
  !
  !! @param[in]    geo       Object containing all geometry data
  !! @param[in]    md_step   MD step             
  !
  subroutine print_DFTB_geo_in_xyz(geo, md_step)
    Type(dftb_geometry_type),   Intent(In) :: geo
    integer,                    Intent(In) :: md_step

    integer            :: ia,itype
    character(len=6)   :: md_step_lab
    character(len=:),  Allocatable :: fname 
    real(wp)           :: a(3), b(3), c(3)

    Write(md_step_lab,'(I6)') md_step
    fname='structure_'//Trim(Adjustl(md_step_lab))
    fname = fname//'.xyz'
    Open(unit=100,file=Trim(Adjustl(fname)))

    !Line 1
    Write(100, '(1X,I5)') Size(geo%coords,2)  

    !Line 2
    If(    (geo%periodic .eqv. .true. ) .and. (geo%fracCoord .eqv. .false.)) Then
       a(:) = geo%lattice(:,1) * Bohr__AA 
       b(:) = geo%lattice(:,2) * Bohr__AA
       c(:) = geo%lattice(:,3) * Bohr__AA
       Write(100,'(A, 8(F8.5,1X),F8.5,A)') 'Lattice="', a(:), b(:), c(:), '"'
       
    Elseif((geo%Periodic .eqv. .false.) .and. (geo%fracCoord .eqv. .false.)) Then
       Write(100,'(A)') "  "
       
    Elseif( (geo%periodic .eqv. .true. ) .and. (geo%fracCoord .eqv. .true.) ) Then
       Write(*,*) 'Not set up for fractional'
       stop
    Endif

    !Line 3 to n_atoms+3
    do ia=1,Size(geo%coords,2)
       itype = geo%species(ia)
       Write(100, '(A,1X,3(f10.6,1X))') trim(adjustl(geo%speciesNames(itype))), geo%coords(1:3,ia)*Bohr__AA
    enddo

    Close(100)
  end subroutine print_DFTB_geo_in_xyz

  
  !> @brief Print DFTB+ geometry data to file, in the .gen format
  !!
  !! Outputs in DFTB+'s gen format, with the naming convention:
  !! structure_mdstep_label.gen
  !!
  !! @param[in]           geo         Object containing all geometry data
  !! @param[in]           md_step     Current MD step
  !! @param[in,optional]  label       Additional label providing more info in filename       
  !
  Subroutine print_DFTB_geometry_data(geo, md_step, label) 
    Type(dftb_geometry_type),   Intent(In) :: geo
    integer,                    Intent(In) :: md_step
    character(len=*), Optional, Intent(In) :: label 
    integer            :: ia
    character(len=1)   :: boundary
    character(len=6)   :: md_step_lab
    character(len=:), Allocatable :: fname 
    character(len=100) :: header

    If(     (geo%periodic .eqv. .true. ) .and. (geo%fracCoord .eqv. .false.)) Then
       boundary='S'
    Elseif( (geo%periodic .eqv. .true. ) .and. (geo%fracCoord .eqv. .true.) ) Then
       boundary='F'
    Elseif( (geo%Periodic .eqv. .false.) .and. (geo%fracCoord .eqv. .false.)) Then
       boundary='C'
    Endif
    
    Write(md_step_lab,'(I6)') md_step
    fname='structure_'//Trim(Adjustl(md_step_lab))

    If(Present(label))then
       fname = fname//'_'//Trim(Adjustl(label))
    Endif
    
    fname = fname//'.gen'
    Open(unit=100,file=Trim(Adjustl(fname)))
      
    !Lines 1 and 2
    header = '# .gen file created by DLPOLY, at MD step '//md_step_lab//&
             '. Atomic positions in ang'
    Write(100,'(A)')  header
    Write(100,'(1X,I3,1X,A)') Size(geo%coords,2), boundary

    !Line 3
    Do ia=1,Size(geo%speciesNames)
       Write(100,'(A,1X)',advance="NO") Trim(Adjustl(geo%speciesNames(ia)))
    Enddo
    Write(100,*)

    !Rest of file
    Do ia=1,Size(geo%coords,2)
       Write(100,'(I3,I3,1X,3(f10.6,1X))') ia, geo%species(ia), geo%coords(1:3,ia)*Bohr__AA
    Enddo

    If(geo%periodic)Then
       !Origin
       Write(100,'(3(f10.6,1X))') 0._wp, 0._wp, 0._wp
       !Stored columnwise internally, row-wise in .gen
       Do ia=1,3
          Write(100,'(3(f10.6,1X))') geo%lattice(1:3,ia) * Bohr__AA 
       Enddo    
    End If
    
    Close(100)

  End Subroutine print_DFTB_geometry_data

  
  FUNCTION lower_case(s1)  RESULT (s2)
    CHARACTER(*)       :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER            :: i
    
    DO i = 1,LEN(s1)
       ch = s1(i:i)
       IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
       s2(i:i) = ch
    END DO
  END FUNCTION lower_Case

  
  FUNCTION upper_case(s1)  RESULT (s2)
    CHARACTER(*)       :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER            :: i
    
    DO i = 1,LEN(s1)
       ch = s1(i:i)
       IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
       s2(i:i) = ch
    END DO
  END FUNCTION upper_case
  
#endif
End Module dftb_library
