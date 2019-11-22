Module test_configuration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl-poly 4 module: Unit tests for subroutines implemented
  ! in configuration.F90, in the absence of a unit-testing framework
  !
  ! copyright - daresbury laboratory
  ! author    - A. Buccheri january 2020
  !
  ! Some notes/Issues. Jan 2020 
  !  One still requires dummy CONFIG, FIELD and CONTROL files to
  !  run this module.
  !  The CONTROL file only needs to contain 'unit_test' and 'rvdw cutoff'
  !  CONFIG and FIELD must be consistent, and contain > 1000 atoms for
  !  testing with 4 cores (else DD will kill the prgram in setup)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use iso_fortran_env, Only : error_unit
  Use kinds, Only : wp
  Use comms, Only : comms_type, root_id, mpi_distribution_type, &
        gatherv_scatterv_index_arrays, op_sum, gallreduce
  Use configuration, Only : configuration_type, gather_coordinates, len_atmnam, &
       coordinate_buffer_type, unpack_gathered_coordinates, gather_atomic_names,&
       distribute_forces
  Use asserts, Only : assert
  Use mpi, Only : MPI_COMM_WORLD, mpi_comm_rank, mpi_barrier
  Implicit None

  Private
  Public :: run_configuration_tests

  !Data for tests 
  Integer, Parameter :: n_atoms = 35
  
  Real(wp), DIMENSION(3, n_atoms) :: position= reshape( &
   (/2.442746018,        0.4025324114E-01,     1.836669280,    &
     1.455907528,         2.350702857,        0.7792575923E-01,&
    -1.071892774,         1.963028995,         3.628125339,    &
   -0.9628280811,         2.989241235,         2.322539321,    &
    0.8130007294,         3.785895632,        0.6008188321,    &
    0.3753060852,         2.013460291,         4.292511636,    &
     2.869103650,         2.554322249,         4.744985661,    &
     3.466028318,        0.5070745127,         2.993257274,    &
     1.676507250,         1.339620571,         1.288168533,    &
     2.115025043,       -0.4892942532E-02,     7.234723249,    &
     1.432829998,         2.375768905,         5.461013556,    &
    -1.147549122,         1.983400427,         8.840267355,    &
   -0.9455381456,         3.095401656,         7.698082405,    &
    0.9844029117,         3.779679892,         6.079565408,    &
    0.2462016987,         1.866273849,         9.640345970,    &
     2.804215918,         2.631374597,         10.16790585,    &
     3.011460904,        0.6656878862,         8.433434549,    &
     1.309060241,         1.170890820,         6.532670282,    &
     2.097689296,        0.7660428222E-01,     12.74562716,    &  
     1.296317611,         2.400231488,         10.77277417,    &
    -1.195496140,         2.159294427,         14.45198593,    &
    -1.105227753,         3.218526343,         13.18332872,    &
    0.6949791977,         3.763828973,         11.32474898,    &
    0.2842851838,         1.895859430,         15.08695423,    &
     2.932551923,         2.635157703,         15.64953749,    &
     3.031796457,        0.5969663469,         14.05064549,    &
     1.351620029,         1.313643155,         11.99743161,    &
     2.086607600,        0.2310361326E-01,     18.17518538,    &
     1.385412938,         2.389960312,         16.21011575,    &
    -1.106759713,         2.049205353,         19.77022701,    &
   -0.9321799496,         3.011469826,         18.45831071,    &
    0.8793036619,         3.889146335,         16.77030058,    &
    0.2779189123,         1.733074931,         20.54423275,    &
     2.749777264,         2.604636095,         21.14705411,    &
     2.994398395,        0.5774762656,         19.42164302/),  (/3,n_atoms/) )

  Character(Len=len_atmnam),Parameter :: si = "Si"
  Character(Len=len_atmnam),Parameter :: o = "O"
  Character(Len=len_atmnam), Dimension(n_atoms) :: species = [ &
       si,si,si,o,o,o,o,o,o,&
       si,si,si,o,o,o,o,o,o,&
       si,si,si,o,o,o,o,o,o,&
       si,si,si,o,o,o,o,o]
  
Contains

  ! Runs all unit tests for configuration module 
  Subroutine run_configuration_tests()
    Type(comms_type)           :: comm
    Type(configuration_type)   :: config
    Logical :: passed(3)
    Integer :: i, n_processes, np, ierr

    !Test coordinates are gathered from all processes and broadcast to all
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, n_processes, ierr)
    Call initialise_comm(n_processes, comm)
    Call initialise_config(comm, position, species, config)

    passed(1) = test_gather_coordinates_1(comm, config)
    passed(2) = test_gather_atomic_names(comm, config)
    passed(3) = test_distribute_forces(comm, config)
    
    If(comm%idnode == root_id) Then
       Write(*,*) "Summary of unit tests for configuration.F90"
       Do i = 1, Size(passed)
          Write(*,*) "Test ",i,": ",passed(i)
       Enddo
    Endif
    
    Call finalise_config(config)
    Call finalise_comm(comm)

  End Subroutine run_configuration_tests

  
  !Initialise a comm object, with required data only
  Subroutine initialise_comm(n_processes, comm)
    Integer,          Intent(In)  :: n_processes
    Type(comms_type), Intent(InOut) :: comm
    Integer                       :: ierr
    comm%comm = MPI_COMM_WORLD
    comm%mxnode = n_processes 
    Call MPI_COMM_RANK(comm%comm, comm%idnode, ierr)
  End Subroutine initialise_comm

  
  Subroutine finalise_comm(comm)
    Type(comms_type), Intent(InOut) :: comm
    comm%comm = 0
    comm%mxnode = 0
    comm%idnode = 0
  End Subroutine finalise_comm

   
  !Initialise a config object for testing for required member data, only
  Subroutine initialise_config(comm, position, species, config)
    Type(comms_type),          Intent(In)  :: comm
    Real(wp),                  Intent(In ) :: position(3, n_atoms)
    Character(Len=len_atmnam), Intent(In)  :: species(n_atoms)
    Type(configuration_type),  Intent(InOut) :: config
    Integer                                :: istart,istop, ia,i
    
    Call assert(size(position, 1) == Size(species), &
         "Size(position,1) /= Size(species)")

    config%megatm = n_atoms
    !Alternative to domain decomposition
    !Call distribute_loop(comm, config%megatm, istart, istop)
    !Doesn't see to be working, so hard-coded my own distributions
    !for up to 4 cores
    Call distribute_loop(comm, istart, istop)
      
    ! Initialise required parts of config
    config%natms = istop - istart +1    
    Allocate(config%parts(config%natms))
    Allocate(config%atmnam(config%natms))

    i = 0
    Do ia = istart, istop
       i=i+1
       config%parts(i)%xxx = position(1,ia)
       config%parts(i)%yyy = position(2,ia)
       config%parts(i)%zzz = position(3,ia)
       config%atmnam(i) = species(ia) 
    Enddo
    
  End Subroutine initialise_config

  Subroutine finalise_config(config)
    Type(configuration_type), Intent(InOut) :: config
    Deallocate(config%parts)
    Deallocate(config%atmnam)
  End Subroutine finalise_config


  ! Distribute coordindates, gather into packed 1D array.
  ! Unpack array and compare to initial coordinates 
  Function test_gather_coordinates_1(comm, config) Result(passed)
    Type(comms_type),            Intent(InOut) :: comm
    Type(configuration_type),    Intent( In )  :: config
    Logical                                    :: passed
    
    Type(coordinate_buffer_type) :: gathered
    Logical, Parameter :: to_master_only = .false.
    Real(wp), Allocatable :: coords(:,:)

    Integer :: ia,i
    
    !Gather coordinates in to packed array
    Call gathered%initialise(comm, config%megatm)         
    Call gather_coordinates(comm, config, to_master_only, gathered)
    Call assert(Size(gathered%coords) == config%megatm*3, &
         "Size(gathered&coords) /= config%megatm*3")
   
    !Unpack into same shape as position
    Allocate(coords(3, config%megatm))
    Call unpack_gathered_coordinates(comm, config, gathered, coords)
 
    passed = all(position == coords)
    If(.not. passed) Then
       Write(error_unit,'(/,1x,a)') "test_gather_coordinates_1: position /= coords"
    Endif
    
    Call gathered%finalise()
  End Function test_gather_coordinates_1


  ! Initialisation CONFIG distributed species in this module
  ! Gather atomic names and compare to species (i.e. source data)
  Function test_gather_atomic_names(comm, config) Result(passed)
    Type(comms_type),            Intent(InOut) :: comm
    Type(configuration_type),    Intent( In )  :: config
    Logical                                    :: passed
    
    Character(Len=len_atmnam),   Allocatable   :: atmnam(:) 
    Logical, Parameter :: to_master_only = .false.

    Allocate(atmnam(config%megatm))
    !Local subset of atmnam in config
    Call gather_atomic_names(comm, config, to_master_only, atmnam)
    
    passed = all(atmnam == species)
    If(.not. passed) Then
       Write(error_unit,'(/,1x,a)') "test_gather_atomic_names: atmnam /= species"
    Endif

    Deallocate(atmnam)
  End Function test_gather_atomic_names


  ! Works on coordinate data as well as forces (no unit conversion) 
  ! Distribute positions, then compare to coordinates distributed
  ! upon initialisation of config
  Function test_distribute_forces(comm, config) Result(passed)
    Type(comms_type),            Intent(InOut) :: comm
    Type(configuration_type),    Intent( In )  :: config
    Logical                                    :: passed

    Type(mpi_distribution_type) :: mpi
    Type(configuration_type)    :: config2
    Real(wp), Allocatable :: coords(:,:), diff(:)
    Integer :: i, istart, istop
    Integer :: passed_local(comm%mxnode), passed_per_process(comm%mxnode)

    !Local coords in config%parts(:)%xxx, etc, of size config%natms
    !Construct mpi index arrays without gathering coordinates into buffer
    Call mpi%init(comm)
    Call gatherv_scatterv_index_arrays(comm, config%natms*3, mpi%counts, mpi%displ)

    
    !Initialise new config object and fill with position
    config2%megatm = n_atoms
    Call distribute_loop(comm, istart, istop)
    config2%natms = istop - istart +1    
    Allocate(config2%parts(config2%natms))

    !Distribute position(3,megatm) into config2%parts(:)%fxx, etc
    Allocate(coords(3, config2%megatm), source=position)
    Call distribute_forces(comm, mpi, coords, config2)

    
    !Check data 
    Allocate(diff(config%natms))
    diff = (config%parts(:)%xxx - config2%parts(:)%fxx) + &
           (config%parts(:)%yyy - config2%parts(:)%fyy) + &
           (config%parts(:)%zzz - config2%parts(:)%fzz)
    
    passed_local = 0 
    If(Any(diff >= 1.e-10)) passed_local(comm%idnode+1) = 1
    !Easier than summing logicals 
    Call gallreduce(comm, passed_local, passed_per_process, op_sum)
    
    If(All(passed_per_process == 0)) Then
       passed = .true.
    Else
       passed = .false.
       If(comm%idnode == root_id) Then
          write(*,*) 'config%parts did not agree on these processes'
          Do i = 1,comm%mxnode
             If(passed_per_process(i)>0)Then
                Write(*,*) "process ",i-1
             Endif
          Enddo
       Endif
    Endif
       
  End Function test_distribute_forces
    
  
  !> Stupid means of loop distribution for testing purposes
  subroutine distribute_loop(comm, istart, istop)
    Type(comms_type), Intent(In) :: comm
    !integer, intent(in)  :: n
    integer, intent(out) :: istart, istop

    if(comm%mxnode == 1)Then
       istart= 1
       istop= 35
       return
    endif
    
    if(comm%mxnode == 2)Then
       if(comm%idnode==0)Then
          istart=1
          istop= 18
       elseif(comm%idnode==1)Then
          istart= 19
          istop = 35
       endif
       return
    endif
    
    if(comm%mxnode == 3)Then
       if(comm%idnode==0)Then
          istart = 1
          istop= 12
       elseif(comm%idnode==1)Then
          istart= 13
          istop = 24
       elseif(comm%idnode==2)Then
          istart= 25
          istop = 35
       endif
    endif
    
    if(comm%mxnode == 4)Then
       if(comm%idnode==0)Then
          istart = 1
          istop= 8
       elseif(comm%idnode==1)Then
          istart = 9 
          istop = 16
       elseif(comm%idnode==2)Then
          istart= 17
          istop = 24
       elseif(comm%idnode==3)Then
          istart= 25
          istop = 35
       endif
    endif
    
  end subroutine distribute_loop

  
End Module test_configuration
