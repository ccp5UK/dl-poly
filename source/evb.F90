Module evb
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module that declares arrays and variables and defines 
  ! subroutines for the computation of the Empirical Valence Bond method 
  ! to simulate reactive dynamics via the coupling of non-reactive 
  ! Force-Fields (FFs)
  !
  ! copyright - daresbury laboratory
  ! author    - i.scivetti march 2019
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use angles,                   Only: angles_type,&
                                      ANGLE_TAB
  Use bonds,                    Only: bonds_type,&
                                      BOND_TAB
  Use comms,                    Only: comms_type,&
                                      gcheck,&
                                      grecv,&
                                      gsend,&
                                      gsum,&
                                      WriteConf_tag 
  Use configuration,            Only: configuration_type
  Use constants,                Only: engunit,&
                                      eu_ev,&
                                      eu_kcpm,&
                                      eu_kjpm,&
                                      boltz
  Use constraints,              Only: constraints_type
  Use core_shell,               Only: core_shell_type
  Use dihedrals,                Only: dihedrals_type,&
                                      DIHEDRAL_TAB
  Use errors_warnings,          Only: error,&
                                      info
  Use external_field,           Only: external_field_type,& 
                                      FIELD_ELECTRIC,& 
                                      FIELD_MAGNETIC 
  Use filename,                 Only: file_type,&  
                                      FILE_POPEVB,&  
                                      FILE_SETEVB  
  Use flow_control,             Only: flow_type,& 
                                      RESTART_KEY_OLD
  Use four_body,                Only: four_body_type
  Use inversions,               Only: inversions_type,&
                                      INVERSION_TAB
  Use kinds,                    Only: wp, &
                                      wi
  Use metal,                    Only: metal_type
  Use numerics,                 Only: invert
  Use particle,                 Only: corePart


  Use parse,                    Only: get_line,&
                                      get_word,&
                                      lower_case,&
                                      clean_string, &
                                      word_2_real,&
                                      strip_blanks
  Use pmf,                      Only: pmf_type

  Use rigid_bodies,             Only: rigid_bodies_type
  Use site,                     Only: site_type
  Use statistics,               Only: stats_type
  Use tersoff,                  Only: tersoff_forces,&
                                      tersoff_type
  Use tethers,                  Only: tethers_type
  Use thermostat,               Only: thermostat_type

  Use three_body,               Only: threebody_type


  Implicit None
  Private

  Type, Public ::  evb_type
!> Type of EVB simulation: standard (0) or multi-component (1). It is set to standard by default
  Integer(Kind = wi)              :: typsim = 0
!> Number of atoms of the EVB reactive unit 
  Integer(Kind = wi), Allocatable :: num_at(:)
!> Number of molecules that make the EVB site (this might depend on the FF)
  Integer(Kind = wi), Allocatable :: FFunits(:)
!> Flag to activate the printing of EVB population    
  Logical                         :: population = .False.
!> Flag for opening EVB population file     
  Logical                         :: population_file_open = .False.
!> Flag for newjob
  Logical                         :: newjob = .True.
!> EVB ionic force
  Real(Kind = wp), Allocatable    :: force(:,:)
!> Energy shift for force fields, convenient to model asymmetries in the potential energy surface (PES)
  Real(Kind = wp), Allocatable    :: eshift(:)
!> EVB Energy matrix
  Real(Kind = wp), Allocatable    :: ene_matrix(:,:)
!> EVB Force matrix
  Real(Kind = wp), Allocatable    :: force_matrix(:,:)
!> Matrix fro computation of the EVB stress
  Real(Kind = wp), Allocatable    :: stress_matrix(:,:)
!> Matrix with the derivative of the EVB energy with respect to each component of the lattice vectors
  Real(Kind = wp), Allocatable    :: dE_dh (:,:)
!> EVB stress tensor
  Real(Kind = wp), Allocatable    :: stress(:,:)
!> EVB eigenvalues
  Real(Kind = wp), Allocatable    :: eigval(:)
!> Energy for each force field, including any potential shift 
  Real(Kind = wp), Allocatable    :: eneFF(:)
!> EVB eigenvectors 
  Real(Kind = wp), Allocatable    :: psi(:,:)
!> Coupling parameters between fields
  Integer                         :: maxfitparam=4        !  maximum number of fitting parameters
  Character(Len = 5), Allocatable :: typcoupl(:,:)        !  Matrix with the specification of the type of coupling   
  Real(Kind = wp),    Allocatable :: couplparam(:,:,:)    !  Matrix of coupling parameters
  Real(Kind = wp),    Allocatable :: grad_coupl(:,:)      !  Matrix for the gradient of coupling parameters
!> Working arrays for diagonalization
  Real(Kind = wp), Allocatable    :: work(:)
  Integer(Kind = wi), Allocatable :: ifail(:), iwork(:)
!> Number of intramolecular interactions corresponding only to the EVB site
  Integer(Kind = wi), Allocatable :: num_bond(:), num_angle(:), num_dihedral(:)
  Integer(Kind = wi), Allocatable :: num_inversion(:)
!> Flag in case of zero coupling
  Logical                         :: no_coupling = .False. 

  Contains
    Private
    Procedure, Public :: init => allocate_evb_arrays
    Final             :: cleanup
  End Type evb_type

  Public :: read_evb_settings  
  Public :: evb_check_configs
  Public :: evb_check_constraints
  Public :: evb_check_intramolecular 
  Public :: evb_check_intermolecular
  Public :: evb_check_external
  Public :: evb_check_intrinsic
  Public :: evb_pes
  Public :: evb_population
  Public :: evb_setzero
  Public :: evb_merge_stochastic
Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! EVB subroutines have been divided in five different categories
!
! 1) Allocation/Deallocation of EVB related variables
!    - allocate_evb_arrays
!    - cleanup
!
! 2) Reading EVB settings from file
!    - read_evb_settings
!
! 3) Checking consistency/correctness for the initial configurations,
!    intrinsic properties and force-field parameters
!    - evb_check_intrinsic
!    - evb_intrinsic_error
!    - evb_check_configs
!    - evb_check_constraints
!    - evb_constraint_error
!    - evb_check_intermolecular
!    - evb_intermolecular_error
!    - evb_check_intramolecular
!    - evb_intramolecular_number_error
!    - evb_check_external
!
! 4) Computation of EVB energy, ionic forces and stress tensor (and virial)
!    - evb_pes
!    - evb_energy
!    - evb_energy_couplings
!    - evb_force
!    - evb_stress
!    - evb_population
!    - evb_setzero
!
! 5) Stochastic dynamics
!    - evb_merge_stochastic
!
! The purpose of each subroutine is described below
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine allocate_evb_arrays(evb,num_ff)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to allocate EVB related variables 
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Class(evb_type), Intent(InOut)     :: evb
    Integer(Kind = wi), Intent (In   ) :: num_ff
    
    Integer :: fail(1:22)
  
    Allocate(evb%eshift(1:num_ff),             Stat=fail(1))
    Allocate(evb%ene_matrix(num_ff,num_ff),    Stat=fail(2))
    Allocate(evb%psi(num_ff,num_ff),           Stat=fail(3))
    Allocate(evb%eigval(num_ff),               Stat=fail(4))
    Allocate(evb%ifail(num_ff),                Stat=fail(5))
    Allocate(evb%iwork(5*num_ff),              Stat=fail(6))
    Allocate(evb%work(8*num_ff),               Stat=fail(7))
    Allocate(evb%eneFF(num_ff),                Stat=fail(8)) 
    Allocate(evb%force(3,num_ff),              Stat=fail(9)) 
    Allocate(evb%force_matrix(num_ff,num_ff),  Stat=fail(10))
    Allocate(evb%dE_dh(3,3),                   Stat=fail(11))
    Allocate(evb%stress(3,3),                  Stat=fail(12))
    Allocate(evb%stress_matrix(num_ff,num_ff), Stat=fail(13))
    Allocate(evb%couplparam(num_ff,num_ff,evb%maxfitparam), Stat=fail(14))
    Allocate(evb%typcoupl(num_ff,num_ff),                   Stat=fail(15))
    Allocate(evb%grad_coupl(num_ff,num_ff),    Stat=fail(16))
    Allocate(evb%FFunits(num_ff),                Stat=fail(17))
    Allocate(evb%num_at(num_ff),               Stat=fail(18))
    Allocate(evb%num_bond(num_ff),             Stat=fail(19))
    Allocate(evb%num_angle(num_ff),            Stat=fail(20))
    Allocate(evb%num_dihedral(num_ff),         Stat=fail(21))
    Allocate(evb%num_inversion(num_ff),        Stat=fail(22))
  
    If (Any(fail /= 0 ))Then
      Call error(0,' error - allocation failure in evb -> allocate_evb_arrays')
    End If
  
  End Subroutine allocate_evb_arrays

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to De-allocate EVB related variables 
    ! allocated in allocate_evb_array
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type) :: T

    If (Allocated(T%eshift)) Then
      Deallocate(T%eshift)
    End If
    If (Allocated(T%ene_matrix)) Then
      Deallocate(T%ene_matrix)
    End If
    If (Allocated(T%psi)) Then
      Deallocate(T%psi)
    End If
    If (Allocated(T%eigval)) Then
      Deallocate(T%eigval)
    End If
    If (Allocated(T%ifail)) Then
      Deallocate(T%ifail)
    End If
    If (Allocated(T%iwork)) Then
      Deallocate(T%iwork)
    End If
    If (Allocated(T%work)) Then
      Deallocate(T%work)
    End If
    If (Allocated(T%eneFF)) Then
      Deallocate(T%eneFF)
    End If
    If (Allocated(T%force_matrix)) Then
      Deallocate(T%force_matrix)
    End If
    If (Allocated(T%force)) Then
      Deallocate(T%force)
    End If
    If (Allocated(T%dE_dh)) Then
      Deallocate(T%dE_dh)
    End If
    If (Allocated(T%stress)) Then
      Deallocate(T%stress)
    End If
    If (Allocated(T%stress_matrix)) Then
      Deallocate(T%stress_matrix)
    End If
    If (Allocated(T%couplparam)) Then
      Deallocate(T%couplparam)
    End If
    If (Allocated(T%typcoupl)) Then
      Deallocate(T%typcoupl)
    End If
    If (Allocated(T%grad_coupl)) Then
      Deallocate(T%grad_coupl)
    End If
    If (Allocated(T%FFunits)) Then
      Deallocate(T%FFunits)
    End If
    If (Allocated(T%num_at)) Then
      Deallocate(T%num_at)
    End If
    If (Allocated(T%num_bond)) Then
      Deallocate(T%num_bond)
    End If
    If (Allocated(T%num_angle)) Then
      Deallocate(T%num_angle)
    End If
    If (Allocated(T%num_dihedral)) Then
      Deallocate(T%num_dihedral)
    End If
    If (Allocated(T%num_inversion)) Then
      Deallocate(T%num_inversion)
    End If

  End Subroutine cleanup        


  Subroutine read_evb_settings(evb, flow, sites, files, comm) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to read setting and parameter for EVB simulations 
    ! from the SETEVB file. After reading, few checks are performed to verify 
    ! the correctness of the input values. Coupling parameters are printed to
    ! OUTPUT
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),   Intent(InOut) :: evb
    Type(flow_type),  Intent(InOut) :: flow
    Type(site_type),  Intent(In   ) :: sites(:)
    Type(file_type),  Intent(InOut) :: files(:)
    Type(comms_type), Intent(InOut) :: comm
  
  
    Logical               :: carry,safe
    Character(Len = 200)  :: record
    Character(Len = 40 )  :: word
  
    Integer               :: ncoupl, icoupl
    Integer               :: m, i, j, k, kparam, ieshift
    Logical               :: FFmolflag, couplerror, eshifterror
  
    Logical, Allocatable  :: couplflag(:,:) 
    Logical, Allocatable  :: eshiftflag(:) 
  
    Character(Len = 256)  :: message, messages(5)
    Character(Len = 256)  :: evbunit
  
    Call info(' ',.True.)
    Call info(' Reading EVB settings from SETEVB file...',.True.) 
    Call info(' ',.True.)
  
    ! Check if units are the same for all FIELD files
    ! To this purpose, field ffunit was added to site_type 
    Do m=1, flow%NUM_FF-1
      If(Abs(sites(m)%ffunit-sites(m+1)%ffunit) >= epsilon(sites(m)%ffunit)) Then
        Write(message, '(1x,2(a,i2),a)') 'error- Units used in FF', m, ' **DIFFER** from those specified in FF ', m+1, &
                                         '. Units MUST be the same for all FIELD files' 
        Call error(0,message)
      End If
    End Do
  
    ! Allocate matrix for checking all possible couplings have been set 
    Allocate(couplflag(flow%NUM_FF,flow%NUM_FF))
    Allocate(eshiftflag(flow%NUM_FF)) 
 
    ! Initialization of flags for control over reading input EVB settings 
    couplflag  = .False.
    eshiftflag = .False.
  
    FFmolflag   = .False.
    couplerror  = .False.
    eshifterror = .False.
  
    !Initialise coupling terms and energy shifts
    evb%couplparam= 0.0_wp
    evb%eshift=0.0_wp
  
    ! Set the total number of coupling elements for later check
    ncoupl=(flow%NUM_FF-1)*flow%NUM_FF/2 
  
    ! initialise counters
    icoupl=0
    ieshift=0
  
    ! Set safe flag
    safe = .True.
  
    ! Open the SETEVB file with EVB settings
  
    If (comm%idnode == 0)Then
      Inquire(File=files(FILE_SETEVB)%filename, Exist=safe)
    End If  
    
    Call gcheck(comm,safe,"enforce")
    
    If (.not.safe) Then
      Write(message,'(1x,a)') 'error - File SETEVB (settings for EVB coupling terms) not found'
      Call error(0,message)
    Else
      If (comm%idnode == 0) Then
        Open(Newunit=files(FILE_SETEVB)%unit_no, File=files(FILE_SETEVB)%filename,Status='old')
      End If
    End If
    Call get_line(safe,files(FILE_SETEVB)%unit_no,record,comm)  
  
    If (safe) Then
      carry = .True.
      Do While (carry)
        Call get_line(safe,files(FILE_SETEVB)%unit_no,record,comm)
        If (.not.safe) Exit
        Call lower_case(record)
        Call get_word(record,word)

        If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
  
  ! Read setting to print (or not) EVB population during the course of the MD sumilation
        Else If (word(1:6) == 'evbpop') Then
          evb%population = .True.
  
  ! Read the number of force fields that describe the EVB unit  
        Else If (word(1:10) == 'evbffunits') Then
          FFmolflag=.True.
          Do k=1, flow%NUM_FF         
            Call get_word(record,word)
            If (word(1:2) == '  ') Then
              Write(message,'(1x,2a,i2,a)') 'error - In file SETEVB, incomplete line after the "evbFFunits" key. ', &
                                            'Value of evbFFunits for FF', k, ' is missing. See manual for correct syntax '
              Call error(0,message)
            Else
              evb%FFunits(k)= Nint(word_2_real(word,0.0_wp))
              If( evb%FFunits(k) <= 0 .Or. evb%FFunits(k) > sites(k)%ntype_mol)Then
                Write(message,'(1x,2(a,i2),2a,i2)') 'error - In file SETEVB, wrong input for index', k,&
                                                    ' of the evbFFunits list (corresponding to FF', k, &
                                                    '). This value should be > 1 and <= MOLECULES, ', &
                                                    'where MOLECULES for this FF is equal to ', sites(k)%ntype_mol
                Call error(0,message)                      
              End If        
            End If
          End Do

  ! Read all coupling elements
  
        Else If(word(1:8) == 'evbcoupl') Then
          Call get_word(record,word)
          i = Nint(word_2_real(word,0.0_wp))
  
          Call get_word(record,word)
          j = Nint(word_2_real(word,0.0_wp))
          If (i< 1 .Or. i > flow%NUM_FF .Or. j < 1 .Or. j > flow%NUM_FF) Then
            Write(message,'(1x,a)') 'error - In file SETEVB, either wrong or incomplete list for the two FFs to be & 
                                     & coupled after option "evbcoupl". ACTION: check syntaxis for this option.' 
            Call error(0,message)  
          End If        
  
          If(i==j)Then
            Write(message,'(1x,a,i2,a)') 'error - In file SETEVB, wrong labelling following a "evbcoupl" option: FF ', i,&
                                         ' CANNOT couple to itself'
            Call error(0,message)
          End If

          If(couplflag(i,j))Then
            Write(message,'(1x,2(a,i2),a)') 'error- In SETEVB file, the coupling between field ',i, ' and ', j, & 
                                             ' cannot be defined more than once. ACTION: remove duplication'
            Call error(0,message)
          Else  
            couplflag(i,j)=.True.
            couplflag(j,i)=.True.
          End If
  
          ! Read type of functional form for the coupling term
          Call get_word(record,evb%typcoupl(i,j))
          If (word(1:3) == '   ' .Or. word(1:1) == '!' .Or. word(1:1) == '#') Then
            couplerror=.True.
          Else 
            evb%typcoupl(j,i)= evb%typcoupl(i,j) 
          End If        
  
          ! Read parameters for coupling. The number of parameters depend on the type, as follows 
          If (evb%typcoupl(i,j) == 'const') Then
            kparam=1      
            Do k=1, kparam        
              Call get_word(record,word)
              If (word(1:3) == '   ' .Or. word(1:1) == '!' .Or. word(1:1) == '#') Then
                couplerror=.True.
              Else
                evb%couplparam(i,j,k)= word_2_real(word)
                evb%couplparam(j,i,k)=evb%couplparam(i,j,k)
              End If
            End Do
          Else If (evb%typcoupl(i,j) == 'gauss') Then
            kparam=evb%maxfitparam
            Do k=1, kparam        
              Call get_word(record,word)
              If (word(1:3) == '   ' .Or. word(1:1) == '!' .Or. word(1:1) == '#') Then
                couplerror=.True.
              Else
                evb%couplparam(i,j,k)= word_2_real(word)
                evb%couplparam(j,i,k)= evb%couplparam(i,j,k)
              End If
            End Do
            If(Abs(evb%couplparam(i,j,3)) < epsilon(evb%couplparam(i,j,3)) )Then
              Write(messages(1),'(1x,a)')          'error - Wrong input parameter A3 for the gauss coupling between' 
              Write(messages(2),'(1x,2(a,i2),2a)') 'FFs', i, ' and ', j , '. This value MUST be different from zero.'
              Write(messages(3),'(1x,a)')          'See manual for details in the functional form of gauss coupling.'
              Call info(messages, 3, .True.)
              Call error(0)
            End If        
          Else 
            couplerror=.True.
          End If  

          If(couplerror)Then
            Write(messages(1),'(1x,a)') 'error - Wrong input data for coupling parameters in file SETEVB. The correct format is: ' 
            Write(messages(2),'(1x,a)') '        evbcoupl       i       j       type      list with (real) coupling-parameters'
            Write(messages(3),'(1x,a)') '        where i and j are the fields to be coupled' 
            Write(messages(4),'(1x,a)') '        If type=const, 1 coupling-parameter is needed'
            Write(messages(5),'(1x,a)') '        If type=gauss, 4 coupling-parameters are needed'
            Call info(messages,5,.True.)
            Call error(0)
          End If        
  
          icoupl=icoupl+1
  
  ! Read the energy shift for each force field
  
        Else If (word(1:8) == 'evbshift') Then
          Call get_word(record,word)
          i = Nint(word_2_real(word,0.0_wp))
          If (i < 1 .Or. i > flow%NUM_FF) Then
            eshifterror = .True.        
          End If
  
          ! Read evb%eshift(i) 
          Call get_word(record,word)
          If(word(1:3) == '   ' .Or. word(1:1) == '!' .Or. word(1:1) == '#') Then
            eshifterror = .True. 
          Else
             If(.Not. eshifterror)Then
               evb%eshift(i)= word_2_real(word)
             End If
          End If 

          If(eshifterror)Then
            Write(messages(1),'(1x,a)') 'error - Wrong input data for energy shift in file SETEVB. The correct format is: ' 
            Write(messages(2),'(1x,a)') '        evbshift       i      E0'
            Write(messages(3),'(1x,a)') '        '
            Write(messages(4),'(1x,a)') '        where i is the field to be shifted by E_i0' 
            Write(messages(5),'(1x,a)') '        This energy shift E_i0 should be in units of the FIELD file that model the state i' 
            Call info(messages,5,.True.)
            Call error(0)
          Else
            If(eshiftflag(i))Then
              Write(message,'(1x,a,i2,a)') 'error- In SETEVB files, energy shift for field ',i,' cannot be defined more than once'
              Call info(message,.True.)
              Call error(0,message)
            Else
              eshiftflag(i)=.True.
            End If        
          End If        

          ieshift=ieshift+1
  
        Else If (word(1:6) == 'finish') Then
          carry=.false.
  
        Else 
          Call strip_blanks(record)
          Write(message,'(1x,3a)') 'error - unknown directive ', word(1:Len_Trim(word)+1), 'found in SETEVB file'
          call error(0,message)
  
        End If
  
      End Do
    End If
  
    If(icoupl  /= ncoupl)Then
      Write(message,'(1x,a)') 'error - In file SETEVB there are missing specification for "evbcoupl". '&
                             &'All combinations, even if they are zero, shall be specified.'        
      Call error(0,message)
    End If        
    If(ieshift /= flow%NUM_FF)Then
      Write(message,'(1x,a)') 'error - In file SETEVB there are missing specification for "evbshift". All fields '&
                             &'must be specified'        
      Call error(0,message)
    End If        
    If(.Not. FFmolflag)Then
      Write(message,'(1x,a)') 'error - Flag evbFFunits not found in the SETEVB file.'    
      Call error(0,message) 
    End If        
  
    Deallocate(eshiftflag, couplflag)
 
    ! Close SETEVB file 
    If (comm%idnode == 0)Then
      Close(files(FILE_SETEVB)%unit_no)   
    End If  
  
    ! For standard EVB (evb%tysim=0), only one reactive EVB site is allowed. 
    ! Such a reactive site can be well composed of one or many fragments, depending on the chemical state
    If(evb%typsim == 0)Then
      Do m=1, flow%NUM_FF  
        Do i=1, evb%FFunits(m)
         If(sites(m)%num_mols(i) /= 1 )Then
           write(messages(1),'(1x,a)') 'error - For standard EVB simulations, only a single EVB reactive site is allowed'
           write(messages(2),'(1x,a)') 'Thus, NUMMOL MUST be equal to 1 for the first evbFFunits types of molecular'
           write(messages(3),'(1x,a)') 'units/fragments that define the EVB site (evbFFunits is set in the SETEVB file)'
           Call info(messages,3,.True.)
           Call error(0)
         End If        
        End Do
      End Do
    End If 
  
    ! Initialise the number of atoms that belong to the EVB site (calculated below)
    evb%num_at    = 0
  
    ! Calculate the total number of atoms being part of the EVB reactive site(s)
    Do m=1, flow%NUM_FF  
      Do i=1, evb%FFunits(m)
        evb%num_at(m)=evb%num_at(m)+sites(m)%num_mols(i)*sites(m)%num_site(i)
      End Do
    End Do
  
    ! Check if the number of EVB atoms is the same for all FFs
    Do m=1, flow%NUM_FF-1  
      If(evb%num_at(m) /= evb%num_at(m+1))Then
        write(messages(1),'(1x,2(a,i2))') 'error - The total number of EVB atoms differ between FF', m, 'and FF', m+1
        write(messages(2),'(1x,a)')       'This means that the number of EVB atoms is different for different FIELD files'
        write(messages(3),'(1x,a)')       'ACTION: check the number of EVB atoms is the same for all FIELD files'
        Call info(messages,3,.True.)
        Call error(0)
      End If        
    End Do
  
  ! Recover the label consistent witht the engunit factor   
  
    If (abs(engunit-eu_ev) <= epsilon(eu_ev)) Then
      evbunit = 'eV'
    Else If (abs(engunit-eu_kcpm) <= epsilon(eu_kcpm)) Then
      evbunit = 'kcal/mol'
    Else If (abs(engunit-eu_kjpm) <= epsilon(eu_kjpm)) Then
      evbunit = 'kJ/mol'
    Else If (abs(engunit-1.0_wp) <= epsilon(1.0_wp)) Then
      evbunit = '10 J/mol'
    Else If (abs(engunit-boltz) <= epsilon(boltz)) Then
      evbunit = 'Kelvin/Boltzmann'
    End If
     
  ! Print info about coupling  
    Write(messages(1),'(1x,a)') '------------'
    Write(messages(2),'(1x,a)') 'EVB settings'
    Write(messages(3),'(1x,a)') '------------'
    Call info(messages,3,.True.) 
    Write(messages(1),'(1x,a)') '1) Coupling terms between Force Fields.'
    Write(messages(2),'(1x,a)') 'The adopted functional form for the coupling between fields i and j (C_{ij})' 
    Write(messages(3),'(1x,a)') 'can be of two types: constant (const) or Gaussian (gauss)'
    Write(messages(4),'(1x,a)') 'See manual for details on the functional forms'
    Write(messages(5),'(1x,a)') 'Details for coupling terms are summarised in the following table:'
    Call info(messages,5,.True.)
    Write(messages(1),'(1x,a)')               '______________________________________________________________________'
    Write(messages(2),'(1x,a,5x,a,5x,a)')           'Force Field pair' , 'Type' , 'Parameters ('//trim(evbunit)//')'
    Write(messages(3),'(5x,a,5x,a,10x,a,4(15x,a))')             'i','j', '    '  , 'A1', 'A2','A3','A4'
    Write(messages(4),'(1x,a)')               '______________________________________________________________________'
    Call info(messages,4,.True.)
    Do i=1,flow%NUM_FF
     Do j=i+1,flow%NUM_FF
       If(evb%typcoupl(i,j)=='const')Then
         Write(message,'(2(4x,i2),10x,a,10x,1(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%couplparam(i,j,k),k=1,1)
         Call info(message,.True.)
       ElseIf(evb%typcoupl(i,j)=='gauss')Then
         Write(message,'(2(4x,i2),10x,a,10x,4(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%couplparam(i,j,k),k=1,4)
         Call info(message,.True.)
       End If  
     End Do
    EnD Do
    Write(message,'(1x,a)') '______________________________________________________________________'
    Call info(message,.True.)
    Call info(' ',.True.)
    Write(messages(1),'(1x,a)') '2) Energy shift for each Force Fields.'
    Write(messages(2),'(1x,a)') 'In addition to the interactions described in each FIELD file, one'
    Write(messages(3),'(1x,a)') 'might want to shift the potential by a certain amount'
    Write(messages(4),'(1x,a)') 'Energy shifts are reported in the following table:'
    Call info(messages,4,.True.)
    Write(messages(1),'(1x,a)') '_______________________________________________________'
    Write(messages(2),'(1x,a,10x,a)') 'Force Field', 'Energy shift ('//trim(evbunit)//')'
    Write(messages(3),'(1x,a)') '_______________________________________________________'
    Call info(messages,3,.True.)
    Do i=1,flow%NUM_FF
      Write(message,'(4x,i2,22x,E12.5)') i , evb%eshift(i)
      Call info(message,.True.)
    EnD Do
    Write(message,'(1x,a)')     '_______________________________________________________'
    Call info(message,.True.)
    Call info(' ',.True.)
    
    ! Convert coupling EVB parameters and energy shifts to internal units
    evb%eshift     = engunit*evb%eshift
    evb%couplparam = engunit*evb%couplparam  
  
    ! If There is no coupling defined, e.i. all type of coupling are constants and equal to zero
    ! set the flag evb%no_coupling = .True. and avoid the execution of evb_setzero 
    icoupl=0
  
    Do i = 1, flow%NUM_FF
       Do j = i+1, flow%NUM_FF
         If( (abs(evb%couplparam(i,j,1)) <= epsilon(evb%couplparam(i,j,1))) .And. &
                  evb%typcoupl(i,j) == 'const'  )Then
            icoupl = icoupl + 1
         End If       
       End Do
    End Do
  
    If( icoupl == ncoupl )Then
      evb%no_coupling=.True.
    End If        

    Call info(' Reading finished!',.True.) 
    Call info(' ',.True.)

  End Subroutine read_evb_settings
          

  Subroutine evb_check_intrinsic(evb, sites, config, flow, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the intrinsic properties for
    ! each atom. Since atomic labels and charges might change for the EVB part 
    ! between different FFs (different chamistry), check is set as follows:
    ! -labels and charges are checked to be the same only for non-EVB part
    ! -masses are checked to be the same for all atoms in the system, as there
    !  is mass conservation for the reactive process
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          
    Type(evb_type),           Intent(In   ) :: evb
    Type(site_type),          Intent(In   ) :: sites(:)
    Type(configuration_type), Intent(In   ) :: config(:)
    Type(flow_type),          Intent(In   ) :: flow
    Type(comms_type),         Intent(InOut) :: comm

    Integer                :: m, i 
    Integer                :: mol1, mol2   
  
    Logical                :: loopi
    Logical                :: lmass, lchg, ltag


    Integer                :: st1num(2), st2num(2)
    Integer                :: list(2)

    Character(Len =  8)    :: st1lab(2), st2lab(2)
    Integer                :: mol(2)
  
    Character(Len = 256)   :: labelint 
     
    Call info(' ',.True.)
    Call info(' EVB checking for consistency in the intrinsic properties of atoms & 
              & between the FIELD fies',.True.) 
    Call info(' ',.True.)
  
    ! Initialise flags
    lmass=.False. 
    lchg=.False. 
    ltag=.False.

    Do m = 1, flow%NUM_FF-1
      loopi=.True.
      i=1
      Do While ( i <= config(m)%natms .And. loopi)
        If(Abs(config(m)%weight(i)-config(m+1)%weight(i)) > &
           epsilon(config(m)%weight(i)))Then
          labelint= 'Mass'   
          lmass=.True.
        End If
        If(config(m)%ltg(i) > evb%num_at(m))Then
          If(config(m)%ltg(i) /= config(m+1)%ltg(i))Then
            ltag = .True.      
            labelint= 'Name'   
          Else If(Abs(config(m)%parts(i)%chge - config(m+1)%parts(i)%chge) > &
                  epsilon(config(m)%parts(i)%chge))Then
            labelint= 'Charge'   
            lchg=.True.   
          End If       
        End If
        If(lmass .Or. lchg .Or. ltag)Then
          list(1)  =config(m)%ltg(i)
          st1lab(1)=config(m)%atmnam(i)
          Call obtain_sites_from_list(1,m,sites,list,st1num,mol1)
          mol(1)=mol1   
          list(1)  =config(m+1)%ltg(i)
          st2lab(2)=config(m+1)%atmnam(i)
          Call obtain_sites_from_list(1,m+1,sites,list,st2num,mol2)
          mol(2)=mol2   
          loopi=.False.
        End If
        i=i+1
      End Do
      Call evb_intrinsic_error(labelint, m, comm, lmass, ltag, lchg, st1num, st2num, st1lab, st2lab, mol)
    End Do

  End Subroutine evb_check_intrinsic


  Subroutine evb_intrinsic_error(string, field, comm, lmass, ltag, lchg, st1num, st2num, st1lab, st2lab, mol)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intrinsic) to print an 
    ! error message and abort if there was an inconsitency found.
    ! 
    ! Communication is needed to detect error in all processors.
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

    Character(Len = 256),  Intent(InOut) :: string
    Integer,               Intent(In   ) :: field 
    Type(comms_type),      Intent(InOut) :: comm
    Logical,               Intent(InOut) :: lmass 
    Logical,               Intent(InOut) :: ltag
    Logical,               Intent(InOut) :: lchg
    Integer,               Intent(InOut) :: st1num(:) 
    Integer,               Intent(InOut) :: st2num(:) 
    Character(Len = 8 ),   Intent(InOut) :: st1lab(:) 
    Character(Len = 8 ),   Intent(InOut) :: st2lab(:) 
    Integer,               Intent(InOut) :: mol(:) 
  
    Character(Len = 256)   :: messages(5)
    Integer                :: jdnode
    Logical                :: siterror
  
    siterror = .False.
     
    Write(messages(3),'(1x,3a)')     'ACTION: check settings in FIELD files. If, in principle, the above element numbers are &
                                     &inconsistent with' 
    Write(messages(4),'(1x,3a)')     'the amount of atomic sites listed for this molecular unit in the FIELD file, please consider &
                                     &the specification'
    Write(messages(5),'(1x,3a)')     'for the repetiton number of each site (if present, repetitions are set following the &
                                     &values for charges)'
 
    If(comm%idnode == 0)Then
      Do jdnode=0,comm%mxnode-1
        If(jdnode>0)Then
          Call grecv(comm, string, jdnode, WriteConf_tag)     
          Call grecv(comm, ltag,   jdnode, WriteConf_tag)     
          Call grecv(comm, lchg,   jdnode, WriteConf_tag)     
          Call grecv(comm, lmass,  jdnode, WriteConf_tag)     
          Call grecv(comm, st1num, jdnode, WriteConf_tag)
          Call grecv(comm, st2num, jdnode, WriteConf_tag)
          Call grecv(comm, st1lab, jdnode, WriteConf_tag)
          Call grecv(comm, st2lab, jdnode, WriteConf_tag)
          Call grecv(comm, mol,    jdnode, WriteConf_tag)
        End If
        If(lmass .Or. ltag .Or. lchg) Then
          Write(messages(1),'(1x,5a,i6,2(a,i2))') 'error- ', trim(string), ' for site ', trim(st1lab(1)), ' corresponding to the ',&
                                                   st1num(1), ' element of molecular unit ', mol(1), ' in the FF ', field 
          siterror=.True.
          If(lmass .Or. lchg)Then                 
            Write(messages(2),'(1x,a,i6,2(a,i2))')  '**DIFFERS** from the value assigned to the equivalent site: element ', &
                                                    st2num(1), ' in molecular unit', mol(2), ' of FF ', field+1     
          Else If(ltag)Then                 
            Write(messages(2),'(1x,4a,i6,2(a,i2))') '**DIFFERS** from the name ', trim(st2lab(1)), 'assigned to the equivalent ', & 
                                                    'site: element ', st2num(1), ' in molecular unit ', mol(2), ' of FF ', field+1     
          End If                              
        End If

        If(siterror)Then
          Call info(messages,5,.True.)
          Call error(0)
        End if  
      End Do
    Else  
      Call gsend(comm, string, 0, WriteConf_tag)     
      Call gsend(comm, ltag,   0, WriteConf_tag)     
      Call gsend(comm, lchg,   0, WriteConf_tag)     
      Call gsend(comm, lmass,  0, WriteConf_tag)     
      Call gsend(comm, st1num, 0, WriteConf_tag)     
      Call gsend(comm, st2num, 0, WriteConf_tag)     
      Call gsend(comm, st1lab, 0, WriteConf_tag)
      Call gsend(comm, st2lab, 0, WriteConf_tag)
      Call gsend(comm, mol,    0, WriteConf_tag)     
    End If 
  
  End Subroutine evb_intrinsic_error


  Subroutine evb_check_configs(config,flow,comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check all CONFIG files (one CONFIG file per FF)
    ! have the same 
    !
    ! - number of atoms
    ! - coordinates
    ! - simulation cell
    ! - symmetry (image convention) 
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(In   ) :: config(:)
    Type(flow_type),          Intent(In   ) :: flow
    Type(comms_type),         Intent(InOut) :: comm
 
  
    Integer                :: i, j, m
    Character(Len = 256)   :: message
    Character(Len = 256)   :: messages(3)
    Character(Len = 1)     :: coord, coord0
    Character(Len = 6)     :: string1, string2
    Logical                :: carry
  
    Integer                :: jdnode, stnum
  
    Real(Kind = wp)        :: cell(flow%NUM_FF,3,3)
  
    Call info(' ',.True.)
    Call info(' EVB check for consistency of atomic coordinates and supercell dimensions between different CONFIG files....',.True.) 
    Call info(' ',.True.)
  
    ! Comparison between different CONFIG files: check the value imcon (image convention/symmetry)
    Do m = 1, flow%NUM_FF-1 
      If(config(m)%imcon /= config(m+1)%imcon)then
        If(m == 1) Then       
          write(string1,'(a1)') ' ' 
        Else
          write(string1,'(i2)') m      
        End If        
        write(string2,'(i2)') m+1
        ! In case an inconsistency has been found, complain and abort       
        Write(message,'(1x,5a)')  'error - Value for imcon *DIFFERS* between ', &
                                  'CONFIG'//adjustl(trim(string1)), ' and ','CONFIG'//adjustl(trim(string2)), ' files'     
        Call error(0,message)
      End If 
    End Do       
  
    ! Comparison between different CONFIG files: check number of atoms
    Do m = 1, flow%NUM_FF-1
      If(config(m)%megatm /= config(m+1)%megatm)Then
        Write(messages(1),'(1x,2(a,i2))') 'error- Different number of atoms for FFs', m, ' and', m+1
        Write(messages(2),'(1x,a)')       'ACTION: check settings for all FIELD/CONFIG files and make sure &
                                          &the number of atoms is the same'
        Call info(messages, 2, .True.)
        Call error(0)
      End If        
    End Do  
  
    ! Comparison between different CONFIG files: check simulation cell 
    ! Arrange config(1)%cell in a more convenient index notation, stored in cell
    Do m = 1, flow%NUM_FF
      Do i=1,3
        Do j=1,3
          cell(m,i,j)=config(m)%cell(3*(i-1)+j)
        End Do
      End Do
    End Do
  
    Do i=1,3
      Do j=1,3
        Do m = 1, flow%NUM_FF-1 
          If(abs(cell(m,i,j)-cell(m+1,i,j))>= epsilon(cell(m,i,j)))then
           If(m == 1) Then       
             write(string1,'(a1)') ' ' 
           Else
             write(string1,'(i2)') m      
           End If        
           write(string2,'(i2)') m+1
           ! In case an inconsistency has been found, complain and abort       
           Write(messages(1),'(1x,a,2i2,a)') 'error- Component (', i, j,') of the simulation cell differs between'
           Write(messages(2),'(1x,4a)')      'CONFIG'//adjustl(trim(string1)), ' and ','CONFIG'//adjustl(trim(string2)), ' files'     
           Call info(messages,2,.True.)
           Call error(0) 
          End If 
        End Do       
      End Do
    End Do
  
    ! Comparison of ionic position between all CONFIG files
    ! set initial flags for comparion
    coord0='w'
    coord=coord0
    carry=.True.
    
    Do i = 1 , config(1)%natms
       Do m = 1, flow%NUM_FF-1
         If(carry)Then
           ! check x-coordinate of atom i      
           If(abs(config(m)%parts(i)%xxx - config(m+1)%parts(i)%xxx) >= epsilon(config(m)%parts(i)%xxx))then
             coord='x'
             stnum=config(m)%ltg(i)
           End If
           ! check y-coordinate of atom i      
           If(abs(config(m)%parts(i)%yyy - config(m+1)%parts(i)%yyy) >= epsilon(config(m)%parts(i)%yyy))then
             coord='y'
             stnum=config(m)%ltg(i)
           End If
           ! check z-coordinate of atom i      
           If(abs(config(m)%parts(i)%zzz - config(m+1)%parts(i)%zzz) >= epsilon(config(m)%parts(i)%zzz))then
             coord='z'
             stnum=config(m)%ltg(i)
           End If
           If(coord /= coord0)Then
             carry=.False.      
             If(m == 1) Then
               write(string1,'(a1)') ' '
             Else
               write(string1,'(i2)') m      
             End If
             write(string2,'(i2)') m+1 
           End If        
         End If
       End Do
     End Do
  
    ! Since atoms are divided in different domains and processors we need to pass the flags (computed above) to the master node
     If(comm%idnode == 0)Then
       Do jdnode=0,comm%mxnode-1
         If(jdnode>0)Then
           Call grecv(comm,stnum,jdnode,WriteConf_tag)
           Call grecv(comm,carry,jdnode,WriteConf_tag)
           Call grecv(comm,coord,jdnode,WriteConf_tag)
           Call grecv(comm,string1,jdnode,WriteConf_tag)
           Call grecv(comm,string2,jdnode,WriteConf_tag)
         End If
         If(.Not. carry)Then
           ! In case an inconsistency has been found, complain and abort       
           Write(messages(1),'(1x,a,i3,3a, i6)') 'error- Process ', jdnode, ' found that the ', coord, '-coordinate for atom', stnum
           Write(messages(2),'(1x,6a)')          '**DIFFERS** between ', 'CONFIG'//adjustl(trim(string1)), ' and ', &   
                                                 'CONFIG'//adjustl(trim(string2))     
           Write(messages(3),'(1x,a)')           'ACTION: Fix the error. Atomic coordinates should not vary between CONFIG files'
           Call info(messages,3,.True.)
           Call error(0) 
         End If
       End Do
     Else  
         Call gsend(comm,stnum, 0,WriteConf_tag)     
         Call gsend(comm,carry, 0,WriteConf_tag)
         Call gsend(comm,coord, 0,WriteConf_tag)     
         Call gsend(comm,string1, 0,WriteConf_tag)     
         Call gsend(comm,string2, 0,WriteConf_tag)     
     End If  
  
  
    Call info(' Check Passed!',.True.) 
    Call info(' ',.True.)

  End Subroutine evb_check_configs


  Subroutine evb_check_constraints(config,cons,cshell,tether,sites,flow,rigid,comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the definition of constraints
    ! between atoms (constraints are defined for each FIELD file). 
    ! Check is carried out over frozen atoms, rigid bond constraint, rigid bodies, 
    ! core shells and tether sites. 
    ! Whatever constraints one defines, they MUST NOT differ between FIELD files.
    ! In other words, a group of atoms cannot be a rigid body for one FIELD 
    ! and a flexible molecule for abother FIELD. The same logic applies to other 
    ! constraints. Even though the EVB could be computed independently, different 
    ! contraint settings will lead to a dynamics that depends on the field, which 
    ! is conceptually wrong.
    ! At first sight, one might consider there is a lot of code repetition below. 
    ! However, this is a consequence of the different formats and types related to 
    ! the settings for each constraint. There might be ways to improve/simplify
    ! the coding, but I (i.scivetti) could not find a simpler way. 
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(In   ) :: config(:)
    Type(constraints_type),   Intent(In   ) :: cons(:)
    Type(core_shell_type),    Intent(In   ) :: cshell(:)
    Type(tethers_type),       Intent(InOut) :: tether(:)
    Type(site_type),          Intent(In   ) :: sites(:)
    Type(flow_type),          Intent(In   ) :: flow
    Type(rigid_bodies_type),  Intent(In   ) :: rigid(:)        
    Type(comms_type),         Intent(InOut) :: comm
  
    Integer                :: m, i, j, l, k
    Integer                :: mol1, mol2   
  
    Character(Len = 256)   :: labelint 
  
    Logical                :: loopi, loopj, loopk
    Logical                :: sitemiss, siteparam
    Integer                :: listdim


    Integer, Allocatable   :: st1num(:), st2num(:)
    Integer, Allocatable   :: list(:)

    Character(Len =  8)    :: st1lab(2), st2lab(2)
    Integer                :: mol(2)
  
    Real(Kind = wp)        :: param(2) 
     
    Call info(' ',.True.)
    Call info(' EVB checking for consistency in the constraints (rigid bodies,& 
              & bond constraints, frozen atoms, core-shells and tether sites)',.True.) 
    Call info(' ',.True.)
  
    
    ! Initialise flags
    sitemiss=.False.
    siteparam=.False.


    ! Check consistency in the number of rigid bodies
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    Do m = 1, flow%NUM_FF-1
       If(rigid(m)%total /= rigid(m+1)%total)Then
         labelint='Total number of rigid body units'
         Call evb_constraint_error(labelint, m)
       End If        
     End Do  
 
    listdim=0
    If(rigid(m)%total /= 0 )Then
 
    ! If there are rigid bodies, set listdim equal to rigid(m)%max_list
      listdim = rigid(1)%max_list

    ! Check all fields have exaclty the same number of maximum number of atoms for a rigid body unit.
    ! If not, print an error
      Do m = 1, flow%NUM_FF-1
        If(rigid(m)%max_list /= rigid(m+1)%max_list)Then
          labelint='maximum number of atoms for a rigid body'       
          Call evb_constraint_error(labelint, m)
        End If        
      End Do
    Else
      listdim=2       
    End If     

    Allocate(st1num(listdim), st2num(listdim),list(listdim))

    !Initilization of vectors
    mol=0
    st1num=0
    st2num=0
    list=0

    If(rigid(m)%total /= 0 )Then
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While ( i <= rigid(m)%max_rigid .And. loopi)
          loopk=.True.
          If(rigid(m)%list(-1,i) /= rigid(m+1)%list(-1,i))Then
            Do l = 1, rigid(m)%list(-1,i)
              list(l) = rigid(m)%list(l,i)
            End Do
            Call obtain_sites_from_list(rigid(m)%list(-1,i),m,sites,list,st1num,mol1)
            mol(1)=mol1 
            Do l = 1, rigid(m+1)%list(-1,i)
              list(l) = rigid(m+1)%list(l,i)
            End Do
            Call obtain_sites_from_list(rigid(m+1)%list(-1,i),m+1,sites,list,st2num,mol2)
            mol(2)=mol2 
            listdim=rigid(m)%list(-1,i)  
            siteparam=.True.
            loopi=.False.
          Else        
            k=1
            Do While ( k <= rigid(m)%list(-1,i) .And. loopk)
              If(rigid(m)%list(k,i) /= rigid(m+1)%list(k,i))Then
                Do l = 1, rigid(m)%list(-1,i)
                  list(l) = rigid(m)%list(l,i)
                End Do
                Call obtain_sites_from_list(rigid(m)%list(-1,i),m,sites,list,st1num,mol1)
                mol(1)=mol1 
                listdim=rigid(m)%list(-1,i)  
                sitemiss=.True.
                loopi=.False.
                loopk=.False.
              End If
              k=k+1
            End Do
          End If
        i=i+1
        End Do
        labelint='Rigid-body'
        Call evb_constraint_error(labelint, m, comm, sitemiss, siteparam, st1num, st2num, st1lab, st2lab, mol, listdim)
      End Do
    End If

  
    ! Check consistency in the frozen nature of ions 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Do m = 1, flow%NUM_FF-1
      loopi=.True.
      i=1
      Do While ( i <= config(m)%natms .And. loopi)
        If(config(m)%lfrzn(i) /= config(m+1)%lfrzn(i))Then
          list(1)  =config(m)%ltg(i)
          st1lab(1)=config(m)%atmnam(i)
          Call obtain_sites_from_list(1,m,sites,list,st1num,mol1)
          mol(1)=mol1   
          list(1)  =config(m+1)%ltg(i)
          st2lab(1)=config(m+1)%atmnam(i)
          Call obtain_sites_from_list(1,m+1,sites,list,st2num,mol2)
          mol(2)=mol2   
          siteparam=.True.
          loopi=.False.
        End If
        i=i+1
      End Do
      labelint='Frozen'
      Call evb_constraint_error(labelint, m, comm, sitemiss, siteparam, st1num, st2num, st1lab, st2lab, mol,1)
    End Do
  
    ! Check consistency in the definition of rigid bond constraints 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check that the number of bonds constraints is the same
    Do m = 1, flow%NUM_FF-1
      If( cons(m)%megcon /= cons(m+1)%megcon )Then
        labelint='Total number of bond-constraints'
        Call evb_constraint_error(labelint, m)
      End If  
    End Do
  
    ! Check that bond constraint specification does not change beetween files
    ! Perform the checking only if the number of bond contraints is different from zero 
    If(cons(1)%megcon /= 0)Then
  
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While (i <= cons(m)%ntcons .And. loopi)
          loopj=.True.
          j=1
            Do While (j <= cons(m+1)%ntcons .And. loopj)
                If((cons(m)%listcon(1,i) == cons(m+1)%listcon(1,j)) .Or. (cons(m)%listcon(1,i) == cons(m+1)%listcon(2,j)) )Then
                  If((cons(m)%listcon(2,i) == cons(m+1)%listcon(1,j)) .Or. (cons(m)%listcon(2,i) == cons(m+1)%listcon(2,j)) )Then
                    loopj=.False.
                    param(1)=cons(m)%prmcon(cons(m)%listcon(0,i))
                    param(2)=cons(m+1)%prmcon(cons(m+1)%listcon(0,j))
                    If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                      loopi=.False.
                      siteparam =.True.
                      Do k=1,2
                       list(k)=cons(m)%listcon(k,i)
                      End Do
                      Call obtain_sites_from_list(2,m,sites,list,st1num,mol1)
                      mol(1)=mol1
                      Do k=1,2
                       list(k)=cons(m+1)%listcon(k,j)
                      End Do
                      Call obtain_sites_from_list(2,m+1,sites,list,st2num,mol2)
                      mol(2)=mol2
                    End If  
                  End If  
                End If
              j=j+1
            End do
    
            If(loopj)Then
              loopi=.False.
              sitemiss=.True.
              Do k=1,2
               list(k)=cons(m)%listcon(k,i)
              End Do
              Call obtain_sites_from_list(2,m,sites,list,st1num,mol1)
              mol(1)=mol1
            End If        
    
            i=i+1
        End Do
  
        labelint='Bond-constraints'
        Call evb_constraint_error(labelint, m, comm, sitemiss, siteparam, st1num, st2num, st1lab, st2lab, mol, 2)
    
      End Do
  
    End If
  
  
    ! Check consistency in the definition of core-shell units 
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    ! Check that the number of core-shell units is the same for all FIELD files
    Do m = 1, flow%NUM_FF-1
      If( cshell(m)%megshl /= cshell(m+1)%megshl )Then
        labelint='Total number of core-shell units'
        Call evb_constraint_error(labelint, m)
      End If  
    End Do
  
    ! Check that core-shell specification does not change beetween files
    ! Perform the checking only if the number of core-shell units is different from zero 
  
    If(cshell(1)%megshl /= 0)Then
  
    ! Initialise flags
    sitemiss=.False.
    siteparam=.False.
  
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While (i <= cshell(m)%ntshl .And. loopi)
          loopj=.True.
          j=1
          Do While (j <= cshell(m+1)%ntshl .And. loopj)
            If((cshell(m)%listshl(1,i) == cshell(m+1)%listshl(1,j)) .Or. &
               (cshell(m)%listshl(1,i) == cshell(m+1)%listshl(2,j)) )Then
              If((cshell(m)%listshl(2,i) == cshell(m+1)%listshl(1,j)) .Or. &
                 (cshell(m)%listshl(2,i) == cshell(m+1)%listshl(2,j)) )Then
                loopj=.False.
                Do l=1,2
                  param(1)=cshell(m)%prmshl(l,cshell(m)%listshl(0,i))
                  param(2)=cshell(m+1)%prmshl(l,cshell(m+1)%listshl(0,j))
                  If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                    loopi=.False.
                    siteparam =.True. 
                    Do k=1,2
                      list(k)=cshell(m)%listshl(k,i)
                    End Do
                    Call obtain_sites_from_list(2,m,sites,list,st1num,mol1)
                    mol(1)=mol1
                    Do k=1,2
                      list(k)=cshell(m+1)%listshl(k,j)
                    End Do
                    Call obtain_sites_from_list(2,m+1,sites,list,st2num,mol2)
                    mol(2)=mol2
                  End If
                End Do
               End If  
            End If
            j=j+1
          End do
    
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            Do k=1,2
              list(k)=cshell(m)%listshl(k,i)
            End Do
            Call obtain_sites_from_list(2,m,sites,list,st1num,mol1)
            mol(1)=mol1
           End If        
        i=i+1
      End Do
  
      labelint='Core-shell'
      Call evb_constraint_error(labelint, m, comm, sitemiss, siteparam, st1num, st2num, st1lab, st2lab, mol, 2)
    
      End Do
  
    End If
  
    ! Check consistency in the definition of tether potentials 
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the total number of tethered sites
    Do m = 1, flow%NUM_FF
      tether(m)%megteth=tether(m)%ntteth
     Call gsum(comm, tether(m)%megteth)
    End Do
  
    Do m = 1, flow%NUM_FF-1
      If( tether(m)%megteth /= tether(m+1)%megteth )Then
        labelint='Total number of tether units'
        Call evb_constraint_error(labelint, m)
      End If  
    End Do
  
    ! Second, check that tether specification does not change beetween files
    ! Perform the checking only if the number of core-shell units is different from zero 
  
    If(tether(1)%megteth /= 0)Then
    ! Initialise flags
    sitemiss=.False.
    siteparam=.False.
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While (i <= tether(m)%ntteth .And. loopi)
          loopj=.True.
          j=1
          Do While (j <= tether(m+1)%ntteth .And. loopj)
            If((tether(m)%listtet(1,i) == tether(m+1)%listtet(1,j)))Then
              loopj=.False.
              If(tether(m)%keytet(tether(m)%listtet(0,i)) /= tether(m+1)%keytet(tether(m+1)%listtet(0,j)))Then
                loopi=.False.
                siteparam=.True.    
                list(1)=tether(m)%listtet(1,i)
                Call obtain_sites_from_list(1,m,sites,list,st1num,mol1)
                mol(1)=mol1
                list(1)=tether(m+1)%listtet(1,j)
                Call obtain_sites_from_list(1,m+1,sites,list,st2num,mol2)
                mol(2)=mol2  
              Else 
                Do l=1,tether(m)%mxpteth  
                  param(1)= tether(m)%prmtet(l,tether(m)%listtet(0,i))
                  param(2)= tether(m+1)%prmtet(l,tether(m+1)%listtet(0,j))
                  If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                    loopi=.False.
                    siteparam =.True. 
                    list(1)=tether(m)%listtet(1,i)
                    Call obtain_sites_from_list(1,m,sites,list,st1num,mol1)
                    mol(1)=mol1
                    list(1)=tether(m+1)%listtet(1,j)
                    Call obtain_sites_from_list(1,m+1,sites,list,st2num,mol2)
                    mol(2)=mol2
                  End If
                End Do
              End If  
            End If
            j=j+1
          End do
    
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            list(1)=tether(m)%listtet(1,i)
            Call obtain_sites_from_list(1,m,sites,list,st1num,mol1)
            mol(1)=mol1
          End If        
          i=i+1
        End Do
        labelint='Tether'
        Call evb_constraint_error(labelint, m, comm, sitemiss, siteparam, st1num, st2num, st1lab, st2lab, mol, 2)
      End Do
    End If
  

    Deallocate(list, st2num, st1num)

    ! Print and exit 
    Call info(' Check Passed!',.True.) 
    Call info(' ',.True.)
  
  End Subroutine evb_check_constraints


  Subroutine obtain_sites_from_list(dimlist,field,sites,list,st,mol) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_constraints) to obtain the atomic
    ! site(s) and corresponding molecule given the index(es) assigned to the atom(s) 
    ! from the input list
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    Integer,         Intent(In   ) :: dimlist  
    Integer,         Intent(In   ) :: field   
    Type(site_type), Intent(In   ) :: sites(:)
    Integer,         Intent(In   ) :: list(:)
    Integer,         Intent(  Out) :: st(:) 
    Integer,         Intent(  Out) :: mol   
  
    Integer :: i, j, k
    Integer :: atnum, nsite
    Integer :: num, den       
  
    Logical :: flag
 
    Do i=1,dimlist  ! The input list has dimension dimlist (single atom, pair of atoms, etc)
      flag=.True.
      atnum=list(i)
  
      Do k=sites(field)%ntype_mol-1,1,-1
       If(flag)Then
         nsite=0
         Do j=k,1,-1
           nsite=nsite+sites(field)%num_mols(j)*sites(field)%num_site(j)
         End Do
         If(nsite < atnum)Then
           num=atnum-nsite
           den=sites(field)%num_site(k+1)
           st(i) = mod(num,den)
           If(st(i)==0)Then
             st(i)=sites(field)%num_site(k+1)
           End If        
           mol=k+1
           flag=.False.
         End If  
        End If
      End Do
  
      If(flag)Then  ! If the above search failed it means 
        num=atnum
        den=sites(field)%num_site(1)      
        st(i)=mod(num,den)
        If(st(i)==0)Then
          st(i)=sites(field)%num_site(1)
        End If
        mol=1     
      End If        
  
    End Do

  End subroutine  obtain_sites_from_list


  Subroutine evb_constraint_error(string, field, comm, sitemiss, siteparam, st1num, st2num, st1lab, st2lab, mol, listdim)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_constraints) to print an 
    ! error message and abort if either
    ! 
    ! - the number for a type of constraint is different between FIELD files 
    !   (no optional variable present when subroutine is called)
    ! - the contraint unit has been defined in one FIELD but not found in other (sitemiss=.True.)
    ! - the spefication for the constraint unit changes between FIELD files (siteparam=.True.)
    ! 
    ! Communication is needed to detect error in all processors.
    ! The format of the output message depends on the type of interactions 
    ! 
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

    Character(Len = 256),  Intent(In   )           :: string
    Integer,               Intent(In   )           :: field 
    Type(comms_type),      Intent(InOut), Optional :: comm
    Logical,               Intent(InOut), Optional :: sitemiss 
    Logical,               Intent(InOut), Optional :: siteparam
    Integer,               Intent(InOut), Optional :: st1num(:) 
    Integer,               Intent(InOut), Optional :: st2num(:) 
    Character(Len = 8 ),   Intent(InOut), Optional :: st1lab(:) 
    Character(Len = 8 ),   Intent(InOut), Optional :: st2lab(:) 
    Integer,               Intent(InOut), Optional :: mol(:) 
    Integer,               Intent(In   ), Optional :: listdim
  
    Character(Len = 256)   :: messages(4)
    Integer                :: jdnode, option
    Logical                :: siterror
    Integer                :: i
  
    siterror = .False.
     
    Write(messages(1),'(1x,a)')           'error-'
    Write(messages(4),'(1x,3a)')          'ACTION: chek settings for ' , trim(string), ' in FIELD files'
  
    If(present(comm))Then
  
      If(string(1:6)=='Tether')Then
        option=0
      Else If(string(1:6)=='Frozen')Then  
        option=1
      Else If(string(1:15) == 'Bond-constraint' .Or. &
              string(1:10) == 'Core-shell')Then
        option=2        
      Else If(string(1:10) == 'Rigid-body')Then  
        option=3
      End If        
  
      If(comm%idnode == 0)Then
        Do jdnode=0,comm%mxnode-1
          If(jdnode>0)Then
            Call grecv(comm, sitemiss  , jdnode, WriteConf_tag)
            Call grecv(comm, siteparam , jdnode, WriteConf_tag)
            Call grecv(comm, st1num    , jdnode, WriteConf_tag)
            Call grecv(comm, st2num    , jdnode, WriteConf_tag)
            Call grecv(comm, st1lab    , jdnode, WriteConf_tag)
            Call grecv(comm, st2lab    , jdnode, WriteConf_tag)
            Call grecv(comm, mol       , jdnode, WriteConf_tag)
          End If
      
          
          If(sitemiss) Then
            If(option == 0)Then
              Write(messages(2),'(1x,2a,2(i4,a),2(i2,a))')   trim(string),' unit for atomic site', st1num(1), &
                                                             ' of molecule type ', mol(1), ' (set in  FF ', field, ')'
            Else If(option == 2)Then
              Write(messages(2),'(1x,a,2(a,i4),2(a,i2),a)')  trim(string),' unit between atomic sites', st1num(1), ' and ',&
                                                             st1num(2), ' of molecule type ', mol(1), ' (set in  FF ', field, ')'
            Else If(option == 3)Then
              Write(messages(2),*)                           trim(string),' unit composed of atomic sites', &
                                                             (st1num(i), i = 1, listdim) , ' of molecule type ', &
                                                             mol(1), ' (set in  FF ', field, ')'
            End If                                         

            Write(messages(3),'(1x,a,i2,a)')                 'could not find its equivalent in FF ', field+1, '.' 
            siterror=.True.
      
          Else If(siteparam) Then
            If(option == 0)Then
              Write(messages(2),'(1x,2a,i4,a,2(i2,a))')      trim(string),' unit specification for atomic site ', st1num(1), & 
                                                             ' of molecule type ', mol(1), ' (set in FF ', field,') **DIFFERS**'
              Write(messages(3),'(1x,a,i4,a,2(i2,a))')       'from the specification for atomic site', & 
                                                              st2num(1), ' of molecule type ', mol(2), ' (set in FF ', field+1,')'
  
            ElseIf(option == 1)Then
              Write(messages(2),'(1x,4a,2(i2,a))')           trim(string),' unit specification for atomic site ', trim(st1lab(1)), & 
                                                             ' of molecule type ', mol(1),        & 
                                                             ' (set in FF ', field,') **DIFFERS** from'
              Write(messages(3),'(1x,3a,2(i2,a))')           'the specification for atomic site ', trim(st2lab(1)),&
                                                             ' of molecule type ', mol(2), ' (set in FF ', field+1,')'
            Else If(option == 2)Then
              Write(messages(2),'(1x,2a,2(i6,a),2(i2,a))')   trim(string),' unit specification between atomic sites ', st1num(1), & 
                                                             ' and ' , st1num(2) ,' of molecule type ', mol(1),        & 
                                                             ' (set in FF ', field,') **DIFFERS** from'
              Write(messages(3),'(1x,3a,2(i6,a),2(i2,a))')   'the ', trim(string), ' specification between atomic sites  ', &
                                                             st2num(1), ' and ', st2num(2),           &
                                                             ' of molecule type ', mol(2), ' (set in FF ', field+1,')'
            Else If(option == 3)Then
              Write(messages(2),*)                           trim(string),' unit composed of atomic sites', &
                                                             (st1num(i), i = 1, listdim) 
              Write(messages(3),'(2(a,i2),a,2(a,i2))')       ' of molecule type', mol(1), ' (set in  FF ', field, &
                                                             ') **DIFFERS** from the what it should be the corresponding unit ',  &
                                                             'defined in the molecule type ', mol(2), ' of FF ', field+1
  
            End If                                        
            siterror=.True.
          End If
          If(siterror)Then
            Call info(messages,4,.True.)
            Call error(0)
          End if  
        End Do
      Else  
        Call gsend(comm, sitemiss , 0, WriteConf_tag)     
        Call gsend(comm, siteparam, 0, WriteConf_tag)     
        Call gsend(comm, st1num   , 0, WriteConf_tag)     
        Call gsend(comm, st2num   , 0, WriteConf_tag)     
        Call gsend(comm, st1lab   , 0, WriteConf_tag)
        Call gsend(comm, st2lab   , 0, WriteConf_tag)
        Call gsend(comm, mol      , 0, WriteConf_tag)     
      End If 
  
    Else
  
      Write(messages(2),'(1x,2a,i2)')    trim(string), ' set in FIELD file', field
      Write(messages(3),'(1x,a,i2)')     '**DIFFERS** from the number of sites in FIELD file ', field+1
      Call info(messages,4,.True.)
      Call error(0)
  
    End If  


  End Subroutine evb_constraint_error


  Subroutine evb_check_intermolecular(evb, flow, sites, tersoff, met, threebody, fourbody)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the definition of intermolecular
    ! interactions between FIELDS files. Check is carried out for tersoff and metallic
    ! potentials, as well as for three and four-body interactions. Since EVB is 
    ! designed to describe reactions (bond breaking and formation) this subroutine also
    ! checks that NONE of the EVB atoms interact via in ANY of this type of intermolecular 
    ! interactions. Thus, the present check applies to all those atoms that are not part of
    ! the EVB reactive site.
    !
    ! For example, EVB could be used to model the bond breaking of a molecule 
    ! supported by a metallic substrate, but metallic interactions will ONLY manifest when 
    ! considering the interactions between atoms of the substrate. Molecule and surface might
    ! interact via vdW forces.
    ! 
    ! For the case of four body interactions, I (i.scivetti) have also decided to prevent 
    ! any possible calculations with EVB, mainly beacuse we have never tested 
    ! this type of potential, to date (Feb 2020).
    !
    ! At first sight, and similarly to subroutine evb_check_constraints, 
    ! one might consider there is a lot of code repetition. 
    ! However, this is a consequence of the different formats and types related to 
    ! the settings for each interaction. 
    !
    ! copyright - daresbury laboratory
    ! copyright - daresbury laboratory
    ! author    - i.scivetti January 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),         Intent(InOut) :: evb
    Type(flow_type),        Intent(In   ) :: flow
    Type(site_type),        Intent(In   ) :: sites(:)
    Type(tersoff_type),     Intent(InOut) :: tersoff(:)
    Type(metal_type),       Intent(InOut) :: met(:)
    Type(threebody_type),   Intent(InOut) :: threebody(:)
    Type(four_body_type),   Intent(InOut) :: fourbody(:)
  
  
    Integer               :: m, i, j, k, l, il, ic
    Integer               :: nprter
  
    Character(Len = 256)  :: labelint 
  
    Logical               :: loopi, loopj
    Logical               :: sitemiss, siteparam
  
    Character(Len = 8)    :: stlab(3)
  
    Real(Kind = wp)       :: param(2) 
  
    Character(Len = 256)  :: messages(5)
  
    ! Initialise flags
    sitemiss=.False.
    siteparam=.False.
  
    Call info(' ',.True.)
    Call info(' EVB checking of intermolecular interations for the non-EVB part....', .True.) 
    Call info(' ',.True.)
  
    ! Check consistency in the definition of METAL potentials
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check that the number of metallic pair interactions are the same for all FIELD files
    Do m = 1, flow%NUM_FF-1
      If( met(m)%n_potentials /= met(m+1)%n_potentials )Then
        labelint='Total number of METAL pair interactions'
        Call evb_intermolecular_error(labelint, m)
      End If  
    End Do
   
    ! Perform the checking only if the number of metallic pair interactions is different from zero 
    If(met(1)%n_potentials /= 0)Then
      ! Check there is no EVB atom that interacts via metallic forces
      Do m = 1, flow%NUM_FF         
        Do i=1,evb%num_at(m)
          Do j= 1, met(m)%n_potentials 
            Do k = 1, 2
              If( trim(sites(m)%site_name(i)) == trim(met(m)%labunit(k,j)) )Then
                write(messages(1),'(1x,3a,i2,a)') 'error- Site ', trim(sites(m)%site_name(i)),  &
                                                  ' belongs to the EVB site in FF ', m,            &
                                                  ' and MUST NOT be part of metallic interactions'    
                Call error(0, messages(1))
              End If        
            End Do
          End Do
        End Do
      End do
     
      ! Check that metallic specification does not change beetween files
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While (i <= met(m)%n_potentials .And. loopi)
          loopj=.True.
          j=1
          Do While (j <= met(m)%n_potentials .And. loopj)
            If((met(m)%labunit(1,i) == met(m+1)%labunit(1,j)) .Or. &
               (met(m)%labunit(1,i) == met(m+1)%labunit(2,j)) )Then
              Do k=1,2
                If(met(m)%labunit(1,i) == met(m+1)%labunit(k,j)) il = k
              End Do
              Do k=1,2
                If( (met(m)%labunit(2,i) == met(m+1)%labunit(k,j)) .And. k /= il )Then
                  loopj=.False.
                  If(met(m)%labunit(3,i) /= met(m+1)%labunit(3,j) ) Then 
                    loopi=.False.
                    siteparam =.True.
                    stlab(1) =  met(m)%labunit(1,i)
                    stlab(2) =  met(m)%labunit(2,i)
                  Else
                    Do l=1,met(m)%max_param
                      param(1)=met(m)%prm(l,i)
                      param(2)=met(m+1)%prm(l,j)
                      If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                        loopi=.False.
                        siteparam =.True.
                        stlab(1) =  met(m)%labunit(1,i)
                        stlab(2) =  met(m)%labunit(2,i) 
                      End If
                    End Do
                  End If  
                End If 
              End Do
            End If
            j=j+1
          End do
      
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            stlab(1)=   met(m)%labunit(1,i)
            stlab(2)=   met(m)%labunit(2,i)
          End If        
    
          i=i+1
    
        End Do
        labelint='Metal'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab) 
      End Do
    End If
  
    ! Check consistency in the definition of TERSOFF potentials
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check that the number of metallic pair interactions are the same for all FIELD files
    Do m = 1, flow%NUM_FF-1
      If( tersoff(m)%n_potential /= tersoff(m+1)%n_potential)Then
        labelint='Total number of Tersoff units'
        Call evb_intermolecular_error(labelint, m)
      End If  
    End Do
  
    ! Perform the checking only if the number of tersoff interactions is different from zero 
    If(tersoff(1)%n_potential /= 0)Then
    ! Check there is no EVB atoms that interacts via Tersoff forces
      Do m = 1, flow%NUM_FF         
        Do i=1,evb%num_at(m)
          Do j= 1, tersoff(m)%n_potential 
              If( trim(sites(m)%site_name(i)) == trim(tersoff(m)%labunit(1,j)) )Then
                write(messages(1),'(1x,3a,i2,a)') 'PROBLEMS: Site ', trim(sites(m)%site_name(i)),  &
                                                  ' belongs to the EVB site in FF ', m,            &
                                                  ' and MUST NOT be part of Tersoff interactions'    
                Call error(0, messages(1))
              End If        
          End Do
        End Do
      End do
      
      ! Check that Tersoff specification for each atomic unit does not change beetween files
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While (i <= tersoff(m)%n_potential .And. loopi)
          loopj=.True.
          j=1
          Do While (j <= tersoff(m)%n_potential .And. loopj)
            If((tersoff(m)%labunit(1,i) == tersoff(m+1)%labunit(1,j)))Then
              If( tersoff(m)%labunit(2,i) == tersoff(m+1)%labunit(2,j) )Then
                loopj=.False.      
                Do l=1,tersoff(m)%max_param
                  param(1)=tersoff(m)%param(l,i)
                  param(2)=tersoff(m+1)%param(l,j)
                  If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                    loopi=.False.
                    siteparam =.True.
                    stlab(1) =  tersoff(m)%labunit(1,i)
                  End If
                End Do
              End If
            End If
            j=j+1
          End do
    
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            stlab(1)= tersoff(m)%labunit(1,i)
          End If        
  
          i=i+1
  
        End Do
  
        labelint='Tersoff-unit'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab) 
    
      End Do
  
      ! Further checking in case the type of Tersoff potential is "tersoff"
      ! Here we chack consistency in the pairwise interction terms of the potential. 
      If(tersoff(1)%key_pot == 1)Then
        nprter = (tersoff(1)%max_ter * (tersoff(1)%max_ter + 1)) / 2
        Do m = 1, flow%NUM_FF-1
          loopi=.True.
          i=1
          Do While (i <= nprter .And. loopi)
            loopj=.True.
            j=1
            Do While (j <= nprter .And. loopj)
              If((tersoff(m)%labpair(1,i) == tersoff(m+1)%labpair(1,j)) .Or. &
                 (tersoff(m)%labpair(1,i) == tersoff(m+1)%labpair(2,j)) )Then
                Do k=1,2
                  If(tersoff(m)%labpair(1,i) == tersoff(m+1)%labpair(k,j)) il = k
                End Do
                Do k=1,2
                  If( (tersoff(m)%labpair(2,i) == tersoff(m+1)%labpair(k,j)) .And. k /= il )Then
                    loopj=.False.
                    Do l=1,2
                      param(1)=tersoff(m)%param2(i,l)
                      param(2)=tersoff(m+1)%param2(j,l)
                      If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                        loopi=.False.
                        siteparam =.True.
                        stlab(1) = tersoff(m)%labpair(1,i)
                        stlab(2) = tersoff(m)%labpair(2,i) 
                      End If
                    End Do
                  End If  
                End Do
              End If
             j=j+1
            End do
      
            If(loopj)Then
              loopi=.False.
              sitemiss=.True.
              stlab(1)= tersoff(m)%labpair(1,i)
              stlab(2)= tersoff(m)%labpair(2,i)
            End If        
    
            i=i+1
    
          End Do
    
        labelint='Tersoff-pair'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab) 
      
        End Do
    
      End If
    
    End If
  
  
    ! Check consistency in the definition of TBP potentials
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check that the number of TBP pair interactions are the same for all FIELD files
    Do m = 1, flow%NUM_FF-1
      If( threebody(m)%ntptbp /= threebody(m+1)%ntptbp )Then
        labelint='Total number of TBP pair interactions'
        Call evb_constraint_error(labelint, m)
      End If  
    End Do
  
    ! Perform the checking only if the number of TBP interactions set is different from zero 
    If(threebody(1)%ntptbp /= 0)Then
    
      ! Check there is no EVB atom that interacts via TBP
      Do m = 1, flow%NUM_FF         
        Do i = 1, evb%num_at(m)
          Do j= 1, threebody(m)%ntptbp
            Do k = 1, 3
              If( trim(sites(m)%site_name(i)) == trim(threebody(m)%labunit(k,j)) )Then
                write(messages(1),'(1x,3a,i2,a)') 'PROBLEMS: Site ', trim(sites(m)%site_name(i)),  &
                                                  ' belongs to the EVB site in FF ', m,            &
                                                  ' and MUST NOT be part of TBP interactions'    
                Call error(0, messages(1))
              End If        
            End Do
          End Do
        End Do
      End do
    
      ! Check that TBP specification does not change beetween files
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        Do While (i <= threebody(m)%ntptbp .And. loopi)
          loopj=.True.
          j=1
          Do While (j <= threebody(m)%ntptbp .And. loopj)
            ic=0
            Do k=1,3
             If(threebody(m)%labunit(k,i) == threebody(m+1)%labunit(k,j))Then
              ic=ic+1
             End If 
            End Do
            If(ic == 3) loopj=.False.
     
            If(.not. loopj)Then
              If(threebody(m)%labunit(4,i) /= threebody(m+1)%labunit(4,j))Then
                loopi=.False.
                siteparam =.True.
                stlab(1) =  threebody(m)%labunit(1,i)
                stlab(2) =  threebody(m)%labunit(2,i) 
              Else
                Do l=1,threebody(m)%mxptbp-1
                  param(1)=threebody(m)%prmtbp(l,i)
                  param(2)=threebody(m+1)%prmtbp(l,j)
                  If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                    loopi=.False.
                    siteparam =.True.
                    Do k=1,3
                      stlab(k) =  threebody(m)%labunit(k,i)
                    End Do
                  End If
                End Do
              End If  
            End If
            j=j+1
          End do
      
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            Do k=1,3
             stlab(k)=   threebody(m)%labunit(k,i)
            End Do
          End If        
    
          i=i+1
    
        End Do
    
        labelint='TBP'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab) 
      
      End Do
    
    End If
  
    ! Check for four-body potentials. If there is any interaction found, kill the job. 
    ! This has to be reviewed once we have a test for fbp 
    Do m=1, flow%NUM_FF-1
      If(fourbody(m)%max_four_body /= 0)Then
        Write(messages(1),'(1x,a)') 'error - EVB  simulations with four body potentials are currently not implemented'      
        Call info(messages,1,.True.)
        Call error(0) 
      End If
    End Do

    ! Print and exit 
    Call info(' Check Passed!',.True.) 
    Call info(' ',.True.)
 
  End Subroutine evb_check_intermolecular


  Subroutine evb_intermolecular_error(string, field, sitemiss, siteparam, stlab)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intermolecular) to print an 
    ! error message and abort if either
    ! 
    ! - the number for a type of intermolecular interactions is different between FIELD files 
    !   (no optional variable present)
    ! - the interaction unit has been defined in one FIELD but not found in other (sitemiss=.True.)
    ! - the spefication for the intermolecular interation unit changes between FIELD files (siteparam=.True.)
    !  
    ! The format of the output message depends on the type of interactions 
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    Character(Len = 256),  Intent(In   )           :: string
    Integer,               Intent(In   )           :: field 
    Logical,               Intent(InOut), Optional :: sitemiss 
    Logical,               Intent(InOut), Optional :: siteparam
    Character(Len = 8),    Intent(InOut), Optional :: stlab(:) 
  
    Character(Len = 256) :: messages(4)
    Integer                :: option
    Logical                :: siterror
  
    siterror = .False. 
    Write(messages(1),'(1x,a)')        'error-'
    Write(messages(4),'(1x,3a)')       'ACTION: chek settings for ' , trim(string), &
                                       ' interactions, as they MUST be the same for all FIELD files'
  
    If(present(stlab))Then
  
      If(string(1:3)=='TBP')Then
        option=0
      Else If(string(1:5)=='Metal')Then  
        option=1
      Else If(string(1:13) == 'Tersoff-unit')Then
        option=2        
      Else If(string(1:13) == 'Tersoff-pair')Then
        option=3        
      End If        
  
          
      If(sitemiss) Then
        If(option == 0)Then
          Write(messages(2),'(1x,8a,i2,a)')           trim(string),' unit ', &
                                                      trim(stlab(1)), '-', trim(stlab(2)), '-', trim(stlab(3)),&
                                                      ' (set in  FF ', field, ')'
        Else If(option == 1)Then
          Write(messages(2),'(1x,6a,i2,a)')           trim(string),' pair unit ', trim(stlab(1)), '-', &
                                                      trim(stlab(2)), ' (set in  FF ', field, ')'
        Else If(option == 2)Then
          Write(messages(2),'(1x,4a,i2,a)')           trim(string), ' ', trim(stlab(1)),' (set in  FF ', field, ')'
  
        Else If(option == 3)Then
          Write(messages(2),'(1x,6a,i2,a)')           trim(string),' ', trim(stlab(1)),'-', trim(stlab(1)) ,&
                                                      ' (set in  FF ', field, ')'
        End If                                     
  
        Write(messages(3),'(a,i2,a)')                 ' could not find its equivalent in FF ', field+1, '.' 
        siterror=.True.
     
      Else If(siteparam) Then
          If(option == 0)Then
            Write(messages(2),'(1x,8a,i2,a)')         trim(string),' specification/parameters for unit ', &
                                                      trim(stlab(1)), '-', trim(stlab(2)), '-', trim(stlab(3)), &
                                                      ' (set in FF ', field,')' 
            Write(messages(3),'(1a,i2)')              ' **DIFFERS** from the specification set in FF ', field+1
  
          Else If(option == 1)Then
            Write(messages(2),'(1x,6a,i2,a)')         trim(string),' specification/parameters for pair ',&
                                                      trim(stlab(1)), '-', trim(stlab(2)), ' (set in FF ', field,')' 
            Write(messages(3),'(2a,i2)')              ' **DIFFERS** from the specification for the same pair ', & 
                                                      'set in FF ', field+1
          Else If(option == 2)Then
            Write(messages(2),'(1x,4a,i2,a)')         trim(string),' specification/parameters ',&
                                                      trim(stlab(1)), ' (set in FF ', field,')' 
            Write(messages(3),'(a,i2)')               ' **DIFFERS** from the specification set in FF ', field+1
        
          Else If(option == 3)Then
            Write(messages(2),'(1x,6a,i2,a)')         trim(string),' specification/parameters for pair ',&
                                                      trim(stlab(1)), '-', trim(stlab(2)), ' (set in FF ', field,')' 
            Write(messages(3),'(2a,i2)')              ' **DIFFERS** from the specification for the same pair ', & 
                                                      'set in FF ', field+1
          End If                                        
          siterror=.True.
      End If
  
      If(siterror)Then
          Call info(messages,4,.True.)
          Call error(0)
      End if  
  
    Else
  
      Write(messages(2),'(1x,2a,i2)')    trim(string), ' set in FIELD file', field
      Write(messages(3),'(a,i2)')        ' **DIFFERS** from the number set in FIELD file ' , field+1
      Call info(messages,4,.True.)
      Call error(0)
  
    End If  

  End Subroutine evb_intermolecular_error



  Subroutine evb_check_intramolecular(evb, flow, sites, bond, angle, dihedral, inversion)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the definition of intramolecular
    ! interactions between FIELDS files ONLY for the non-reactive part of the system.
    ! In fact, functional forms nor parameters for interactions SHOULD change
    ! between atoms of the non-reacitve part of the system.
    !
    ! Check is carried out for bond, angle, dihedral and inversion potentials.
    ! As part of the checking, the number of intramolecular interactions of each type
    ! is also computed for the EVB site.
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti February 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),             Intent(InOut) :: evb
    Type(flow_type),            Intent(In   ) :: flow
    Type(site_type),            Intent(In   ) :: sites(:)
    Type(bonds_type),           Intent(In   ) :: bond(:)
    Type(angles_type),          Intent(In   ) :: angle(:)
    Type(dihedrals_type),       Intent(In   ) :: dihedral(:)
    Type(inversions_type),      Intent(In   ) :: inversion(:)
  
    Integer              :: m, i
  
    Character(Len = 256) :: labelint 
  
    ! Temporary array to count total number of interactions)
    Integer, Allocatable :: numtot(:)
  
    Allocate(numtot(flow%NUM_FF))
  
    Call info(' ',.True.)
    Call info(' EVB checking of intramolecular interations for the non-reactive part of the system', .True.)
    Call info(' ',.True.)
  
    ! Initialise the number for each type of interaction at EVB site 
    evb%num_bond       = 0 
    evb%num_angle      = 0
    evb%num_dihedral   = 0
    evb%num_inversion  = 0
  
    ! Check consistency in the definition of bonds
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    labelint= 'Bond'
    numtot=0

    ! For each FIELD, calculate the number of bonds for the EVB site and the total number of bonds for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+bond(m)%num(i)
        If(i <= evb%FFunits(m))Then   
          evb%num_bond(m) = evb%num_bond(m)+bond(m)%num(i)
        End If        
      End Do
    End Do 
  
    ! Compare if the total number of bonds for the non-reactive part is the same for all FFs (if not, print error)
    Do m=1, flow%NUM_FF-1
      If( (numtot(m)-evb%num_bond(m)) /= (numtot(m+1)-evb%num_bond(m+1)) )Then
        Call evb_intramolecular_number_error(labelint, m)
      End If        
    End Do
  
    ! For the non-reactive part, check that specification bonds does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 2, bond(m)%max_param, evb%num_bond(m), evb%num_bond(m+1), &
                               numtot(m), numtot(m+1), bond(m)%key, bond(m+1)%key, bond(m)%lst, bond(m+1)%lst, &
                               bond(m)%param, bond(m+1)%param, bond(m)%num, sites(m)%ntype_mol, BOND_TAB)
    End Do                      
  
  
    ! Check consistency in the definition of angles 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    labelint= 'Angle'
    numtot=0
  
    ! For each FIELD, calculate the number of angles for the EVB site and the total number of bonds for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+angle(m)%num(i)
        If(i <= evb%FFunits(m))Then    
          evb%num_angle(m) = evb%num_angle(m)+angle(m)%num(i)
        End If
      End Do
    End Do 
  
    ! Compare if the total number of angles for the non-reactive part is the same for all FFs (if not, print error)
    Do m=1, flow%NUM_FF-1  
      If( (numtot(m)-evb%num_angle(m)) /= (numtot(m+1)-evb%num_angle(m+1)) )Then
        Call evb_intramolecular_number_error(labelint, m)
      End If        
    End Do
  
    ! For the non-reactive part, check that specification angles does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 3, angle(m)%max_param, evb%num_angle(m), evb%num_angle(m+1), &
                                  numtot(m), numtot(m+1), angle(m)%key, angle(m+1)%key, angle(m)%lst, angle(m+1)%lst, &
                                  angle(m)%param, angle(m+1)%param, angle(m)%num, sites(m)%ntype_mol, ANGLE_TAB)
    End Do                      
  
  
    ! Check consistency in the definition of dihedrals 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    labelint= 'Dihedral'
    numtot=0
  
    ! For each FIELD, calculate the number of dihedrals for the EVB site and the total number of bonds for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+dihedral(m)%num(i)
        If(i <= evb%FFunits(m))Then    
          evb%num_dihedral(m) = evb%num_dihedral(m)+dihedral(m)%num(i)
        End If
      End Do
    End Do 
  
    ! Compare if the total number of dihedrals for the non-reactive part is the same for all FFs (if not, print error)
    Do m=1, flow%NUM_FF-1  
      If( (numtot(m)-evb%num_dihedral(m)) /= (numtot(m+1)-evb%num_dihedral(m+1)) )Then
        Call evb_intramolecular_number_error(labelint, m)
      End If        
    End Do
  
    ! For the non-reactive part, check that specification dihedrals does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 4, dihedral(m)%max_param, evb%num_dihedral(m), evb%num_dihedral(m+1), &
                                  numtot(m), numtot(m+1), dihedral(m)%key, dihedral(m+1)%key, dihedral(m)%lst, dihedral(m+1)%lst, &
                                  dihedral(m)%param, dihedral(m+1)%param, dihedral(m)%num, sites(m)%ntype_mol, DIHEDRAL_TAB)
    End Do                      
  
  
    ! Check consistency in the definition of inversions 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    labelint= 'Inversion'
    numtot=0
  
    ! For each FIELD, calculate the number of inversions for the EVB site and the total number of bonds for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+inversion(m)%num(i)
        If(i <= evb%FFunits(m))Then    
          evb%num_inversion(m) = evb%num_inversion(m)+inversion(m)%num(i)
        End If
      End Do
    End Do 
  
    ! Compare if the total number of inversions for the non-reactive part is the same for all FFs (if not, print error)
    Do m=1, flow%NUM_FF-1  
      If( (numtot(m)-evb%num_inversion(m)) /= (numtot(m+1)-evb%num_inversion(m+1)) )Then
        Call evb_intramolecular_number_error(labelint, m)
      End If        
    End Do
  
    ! For the non-reactive part, check that specification inversions does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 4, inversion(m)%max_param, evb%num_inversion(m), evb%num_inversion(m+1), &
                               numtot(m), numtot(m+1), inversion(m)%key, inversion(m+1)%key, inversion(m)%lst, inversion(m+1)%lst, &
                               inversion(m)%param, inversion(m+1)%param, inversion(m)%num, sites(m)%ntype_mol, INVERSION_TAB)
    End Do                      
  
    ! Deallocation of temparary array numtot
    Deallocate(numtot)
  
  End Subroutine evb_check_intramolecular


  Subroutine compare_intramolecular(labelint, field, listdim, maxparam, evb_num1, evb_num2, &
                                 numtot1, numtot2, key1, key2, list1, list2,             &
                                 param1, param2, numxmol, ntype_mol, keyexcl) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intramolecular) to compare intramolecular
    ! interactions between FIELDS and print an error message (and abort) if either
    ! 
    ! - the interaction unit has been defined in one FIELD but not found in other (llist=.True.)
    ! - the spefication for the intermolecular interation unit changes between FIELD files (lparam or lkey =.True.)
    !  
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len = 256), Intent(In   ) :: labelint 
    Integer(Kind = wi),   Intent(In   ) :: field, listdim, maxparam
    Integer(Kind = wi),   Intent(In   ) :: evb_num1, evb_num2
    Integer(Kind = wi),   Intent(In   ) :: numtot1 , numtot2
    Integer(Kind = wi),   Intent(In   ) :: key1(:) , key2(:)   
    Integer(Kind = wi),   Intent(In   ) :: list1(:,:), list2(:,:)   
    Integer(Kind = wi),   Intent(In   ) :: numxmol(:)
    Integer(Kind = wi),   Intent(In   ) :: ntype_mol
    Integer(Kind = wi),   Intent(In   ) :: keyexcl 
    Real(Kind = wp),      Intent(In   ) :: param1(:,:), param2(:,:) 

    Integer                :: i, j, k, mol
    Logical                :: llist, lparam, lkey
    Character(Len = 256)   :: messages(6)

    ! Initialise flags
    llist=.False.
    lparam=.False.
    lkey=.False.

    i=evb_num1+1
    j=evb_num2+1

    Do While (i <= numtot1 .And. j <= numtot2)
      Do k=1,listdim
        If( list1(k,i) /= list2(k,j))Then
          llist=.True.       
        End If           
      End Do
      If(.Not. llist)Then
        If( key1(i) /= key2(j) ) Then
          lkey=.True.
        End If        
        If(key1(i) /= keyexcl)Then
          Do k=1,maxparam
            If( abs(param1(k,i)-param2(k,j)) >= epsilon(param1(k,i)) ) Then
              lparam=.True. 
            End If
          End Do
        End If
      End If

      If(lkey .Or. llist .Or. lparam)Then

        Call obtain_molecule_from_intra(i, numxmol, ntype_mol, mol)

        write(messages(1), '(1x,a)' )        'error-'
        write(messages(2), * )                trim(labelint),' interaction between atoms', (list1(k,i), k = 1, listdim) 
        write(messages(3),'(1x,a,2(i2,a))')  'is defined in FF ', field, & 
                                             ' for molecular unit', mol, ' (see the corresponding FIELD file).' 

        If(lparam)Then
          write(messages(4),'(1x,a)')        'Parameters for this unit *DO NOT* agree with the values set'
        End If                      
        If(lkey)Then
          write(messages(4),'(1x,a)')        'Key-type of interaction for this unit *DIFFERS* from the type specified'
        End If                      
        If(llist)Then
          write(messages(4),'(1x,a)')        'This interaction unit CANNOT be found '
        End If                      

        write(messages(5),'(1x,a,i2,3a)')    'in FF', field+1, '. ACTION: find unit and verify it is equally defined for ', &
                                             'all FIELD files, as settings for ', trim(labelint)
        write(messages(6),'(1x,a)')          'interactions for the non-reactive part of the system must be the same for all FFs'
        Call info(messages,6,.True.)
        call error(0)
      End If

      i=i+1
      j=j+1
    End Do

  End Subroutine compare_intramolecular


  Subroutine obtain_molecule_from_intra(num_unit, numxmol, ntype_mol, mol)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to compare_intramolecular) to obtain the molecular unit
    ! given the index of an intramolecular interaction
    !  
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Integer(Kind = wi),   Intent(In   ) :: num_unit
    Integer(Kind = wi),   Intent(In   ) :: numxmol(:)
    Integer(Kind = wi),   Intent(In   ) :: ntype_mol
    Integer(Kind = wi),   Intent(InOut) :: mol
 
    Integer :: j, k
    Integer :: intnum 
 
    Logical :: flag
 
    flag=.True.
 
    Do k=ntype_mol-1,1,-1
      If(flag)Then
        intnum=0
        Do j=k,1,-1
          intnum=intnum+numxmol(j)
        End Do
        If(intnum < num_unit)Then
          mol=k+1
          flag=.False.
        End If  
      End If
    End Do
 

  End Subroutine obtain_molecule_from_intra        


  Subroutine evb_intramolecular_number_error(labelint, field)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intramolecular) to print
    ! and arror message if the number of intramolecular interactions (of a
    ! given type) for the non-reactive part of the system is not the same 
    ! for all FIELD files
    !  
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    Character(Len = 256), Intent(In   ) :: labelint 
    Integer(Kind = wi),   Intent(In   ) :: field
    
    
    Character(Len = 256) :: messages(5)
  
    write(messages(1),'(1x,3a)')        'error- The total number of intramolecular ', trim(labelint),&
                                        ' interactions between atoms that' 
    write(messages(2),'(1x,a,i2,a,i2)') 'are NOT part of the EVB reactive fragment **DIFFERS** between FF ',&
                                         field, ' and FF ', field+1 
    write(messages(3),'(1x,3a)')        'ACTION: check settings for ', trim(labelint), ' interactions and make sure' 
    write(messages(4),'(1x,a)')         'they are equally defined in all FIELD files.'
    write(messages(5),'(1x,a)')         'Settings for the non-reactive part of the system MUST be the same for all FIELD files'
    Call info(messages,5,.True.)
    Call error(0)

  End Subroutine evb_intramolecular_number_error


  Subroutine evb_check_external(evb, flow, ext_field) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that aborts the execution of a standard EVB runs 
    ! if either electric of magnetic fields are present. 
    !  
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),             Intent(InOut) :: evb
    Type(flow_type),            Intent(In   ) :: flow
    Type(external_field_type),  Intent(In   ) :: ext_field(:)
  
    Integer :: m
    Character(Len = 256) :: message 
  
    ! For standard EVB (evb%tysim=0), only one reactive EVB site is allowed. 
    ! Such a reactive site can be well composed of one or many fragments, depending on the chemical state 
    If(evb%typsim == 0)Then
  
      ! Kill the run for standard EVB with external electric or magnetic fields
      Do m=1, flow%NUM_FF
        If (ext_field(m)%key == FIELD_ELECTRIC .Or. ext_field(m)%key == FIELD_MAGNETIC ) Then
          Write(message, '(1x,a)') 'error - Standard EVB formalism is not valid for external Electric or Magnetic fields'
          Call error(0,message)
        End If
      End Do
  
    End If 

  End Subroutine evb_check_external




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines for the computation of the EVB problem 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Subroutine evb_pes(evb,flow,config,stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that calls for the computation of 
    ! - EVB energy
    ! - Gradients of the EVB coupling terms
    ! - EVB forces 
    ! - EVB stress
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti Febrauary 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    Type(evb_type), Intent(InOut) :: evb
    Type(flow_type), Intent(In   ) :: flow
    Type(configuration_type), Intent(InOut) :: config(:)
    Type(stats_type), Intent(InOut) :: stat(:)
  
  
    !Compute EVB Energy 
    Call evb_energy(evb,flow,stat)
  
    !Compute gradients for the coupling terms of the EVB matrix
    Call grad_evb_couplings(flow,evb)
  
    !Compute EVB forces
    Call evb_force(evb,flow,config)
  
    !Compute EVB stress
    Call evb_stress(evb,flow,config,stat)


  End Subroutine evb_pes


  Subroutine evb_energy(evb,flow,stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes the EVB energy. Steps:
    !
    ! - buildig the EVB matrix. Conformational energies for each FF are set in the
    !   diagonal terms; off diagonal terms are computed via suroutine 
    !   evb_energy_couplings
    ! - Diagonalization of the EVB matrix via the BLAS subroutine dsyevx
    ! - Check if diagonalization was successful
    ! - EVB eigenvectors are assigned to evb%psi
    ! - Lowest EVB eigenvalue evb%eigval is the conformational energy 
    !
    ! This procedure is repeated at each time step
    !
    ! #ifdef condition is to allow compilation in serial
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti August 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    Type(evb_type),   Intent(InOut) :: evb
    Type(flow_type),  Intent(In   ) :: flow
    Type(stats_type), Intent(InOut) :: stat(:)
  
    Integer(Kind = wi) :: m,k           ! Indices for matrix elements
    Integer(Kind = wi) :: mevb, evbinfo ! Working intergers

    Character(Len = 256)   :: messages(5)

#ifdef EVB
    !Initialise matrix elements
    evb%ene_matrix = 0.0_wp
  
    ! Build EVB matrix 
    !!!!!!!!!!!!!!!!!!
    ! Diagonal elements
    Do m=1,flow%NUM_FF 
      evb%eneFF(m) = stat(m)%stpcfg+evb%eshift(m)
      evb%ene_matrix(m,m)= evb%eneFF(m)
    End Do
  
    !Off-diagonal terms
    Do m=1,flow%NUM_FF
      Do k=m+1,flow%NUM_FF
        ! Compute coupling, element (k,m) of the EVB matrix   
        Call evb_energy_couplings(m,k,evb)
        ! Make EVB matrix symmetric, this is a MUST to ensure Hermiticity
        evb%ene_matrix(k,m) =evb%ene_matrix(m,k)
      End Do
    End Do
    
    !Diagonalization
    !!!!!!!!!!!!!!!!
    call dsyevx( 'V', 'I', 'U', flow%NUM_FF, evb%ene_matrix, flow%NUM_FF, -1., 0., 1, flow%NUM_FF, &
                  1.d-30, mevb, evb%eigval, evb%psi, flow%NUM_FF, evb%work, 8*flow%NUM_FF,     &
                  evb%iwork, evb%ifail, evbinfo)
 
    !Check that diagonalization was successful
    If(evbinfo /= 0)Then
      If(evbinfo < 0)Then
        write(messages(1),'(1x,a,i2,2a)')  'error- EVB diagonalization failed as the', evbinfo, '-th argument ',&
                                           'has and illegal value.'   
      Else If(evbinfo > 0)Then
        write(messages(1),'(1x,a,i2,2a)')  'error- EVB diagonalization failed as eigenvector ', evbinfo, ' failed ',&
                                           'to converge'   
      End If

      write(messages(3),'(1x,a)')          'ACTION 1: check for incorrect input parameters for coupling terms in SETEVB'
      write(messages(4),'(1x,a)')          'ACTION 2: check for incorrect input parameters in FIELD and CONFIG files. It'
      write(messages(5),'(1x,a)')          'might be helful to run standard MD simulations for each FFs separately'

      Call info(messages,5,.True.)
      Call error(0)
    End If        

    !Assign eigenvalues to conformational energy
    Do m=1,flow%NUM_FF
      stat(m)%stpcfg=evb%eigval(1)
    End Do
#endif    

  End subroutine evb_energy


  Subroutine evb_energy_couplings(m,k,evb)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes the coupling (off-diagonal) terms 
    ! of the EVB matrix 
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti August 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer       , Intent(In   ) :: m,k
    Type(evb_type), Intent(InOut) :: evb
  
    Real(Kind = wp)           :: ediff, elimit
    Real(Kind = wp)           :: A(evb%maxfitparam)
    Character(Len = 256)      :: messages(4)
    
    ! Copy evb coupling parameters to array A, for the sake of clarity in coding the formulae
    A(:)  = evb%couplparam(m,k,:)
    ! Define energy difference, which is the reaction corrdinate between fields m and k
    ediff = evb%eneFF(m)-evb%eneFF(k)
   
    If(evb%typcoupl(m,k)=='const')Then
      evb%ene_matrix(m,k) = A(1)
    Else If (evb%typcoupl(m,k)=='gauss')Then
      elimit=Abs((ediff-A(2))/A(3))
      If(elimit < 700.0_wp)Then     ! Here we use the sensible limit of 700 for the argument to compute exp()
        evb%ene_matrix(m,k) =  A(1)*exp(-((ediff-A(2))/A(3))**2)+A(4)
      Else
        write(messages(1),'(1x,a)')          'error- Argument for the exponential term of the gauss coupling'
        write(messages(2),'(1x,2(a,i2),a)')  'between FFs', m, ' and ', k, ' is beyond the range of computation'  
        write(messages(3),'(1x,a)')          '(-700.0 < argument < 700.0). Either one energy term has diverged or'
        write(messages(4),'(1x,a)')          'there is a problem with the input parameters set in SETEVB.'
        Call info(messages, 4, .True.)
        Call error(0)
      End If  
    End If        
  
  
    End Subroutine evb_energy_couplings
  
    Subroutine grad_evb_couplings(flow,evb)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes the derivative of the coupling terms 
    ! with respect to the energy difference ediff= E_{FF,m}-E_{FF,k}
    ! The resulting matrix evb%grad_coupl is then used to compute EFV forces and
    ! stress tensor
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti August 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      Type(flow_type),  Intent(In   ) :: flow
      Type(evb_type  ), Intent(InOut) :: evb

      Integer                         :: m,k
      Real(Kind = wp)                 :: ediff, elimit
      Real(Kind = wp)                 :: A(evb%maxfitparam)
      Character(Len = 256)            :: messages(4)
    
      evb%grad_coupl=0.0d0
    
      Do m=1,flow%NUM_FF
        Do k=m+1,flow%NUM_FF
        ! Copy evb coupling parameters to array A, for the sake of clarity in coding the formulae
          A(:)  = evb%couplparam(m,k,:)
        ! Define energy difference, which is the reaction corrdinate between fields m and k
          ediff = evb%eneFF(m)-evb%eneFF(k)
          If(evb%typcoupl(m,k)=='const')Then
            evb%grad_coupl(m,k) = 0.0d0 
          Else If (evb%typcoupl(m,k)=='gauss')Then
            elimit=Abs((ediff-A(2))/A(3))
            If(elimit < 700.0_wp)Then     ! Here we use the sensible limit of 700 for the argument to compute exp()
              evb%grad_coupl(m,k) = -2.0_wp * (A(1)/A(3)**2)*(ediff-A(2))*exp(-((ediff-A(2))/A(3))**2)
            Else
              write(messages(1),'(1x,a)')          'error- Argument for the exponential term of the gauss coupling'
              write(messages(2),'(1x,2(a,i2),a)')  'between FFs', m, ' and ', k, ' is beyond the range of computation'  
              write(messages(3),'(1x,a)')          '(-700.0 < argument < 700.0). Either one energy term has diverged or'
              write(messages(4),'(1x,a)')          'there is a problem with the input parameters set in SETEVB.'
              Call info(messages, 4, .True.)
              Call error(0)
            End If
          End If       
    
        End Do
      End Do 


  End Subroutine grad_evb_couplings

  
  Subroutine evb_force(evb,flow,config)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes for the EVB force over each ion of the
    ! system. 
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti September 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    Type(evb_type),           Intent(InOut) :: evb
    Type(flow_type),          Intent(In   ) :: flow
    Type(configuration_type), Intent(InOut) :: config(:)
  
  
    Integer(Kind = wi) :: m,k         ! Indices for matrix elements
    Integer(Kind = wi) :: i           ! Indices for atoms
    Integer(Kind = wi) :: alpha       ! Indices for cartesian coordinates x, y and z
  
    Integer(Kind = wi) :: natms
    Real(Kind = wp)    :: fdiff
  
    !Maximum number of atoms. 
    !In evb_check_config, config%natms has been checked to be the same for all CONFIG files
    natms=config(1)%natms 
  
    Do i=1, natms ! Copy forces from config%parts(i) to evb%force. x->1, y->2, y->3. We do this for each force field.

      Do m=1,flow%NUM_FF
        evb%force(1,m)=config(m)%parts(i)%fxx
        evb%force(2,m)=config(m)%parts(i)%fyy
        evb%force(3,m)=config(m)%parts(i)%fzz
      End Do

      Do alpha=1,3 !For particle i, loop over the three coordinates
        evb%force_matrix=0.0_wp
        ! Diagonal elements
        Do m=1,flow%NUM_FF
          evb%force_matrix(m,m)=evb%force(alpha,m)
        End Do
        ! Off-diagonal elements     
        Do m=1,flow%NUM_FF
          Do k=m+1,flow%NUM_FF
            fdiff=evb%force(alpha,m)-evb%force(alpha,k)       
            evb%force_matrix(m,k)=fdiff*evb%grad_coupl(m,k) 
            !Make matrix symmetric
            evb%force_matrix(k,m)=evb%force_matrix(m,k)
          End Do
        End Do
  
        !Set the first column of evb%force to zero   
        evb%force(alpha,1)=0.0_wp
  
        ! Matrix multiplication (evb%psi)^T*evb%force_matrix*(evb%psi) to compute evb%force(j,1)
        ! Store EVB force in column 1 of evb%force
        Do m=1,flow%NUM_FF
          Do k=1,flow%NUM_FF
            evb%force(alpha,1)=evb%force(alpha,1)+evb%psi(m,1)*evb%force_matrix(m,k)*evb%psi(k,1)
          End Do
        End Do
      End Do !Finish loop over coordinates
  
    ! Once we have compute the the three force components, we copy evb-force to the configurational force of each FF
      Do m=1,flow%NUM_FF
        config(m)%parts(i)%fxx=evb%force(1,1)
        config(m)%parts(i)%fyy=evb%force(2,1)
        config(m)%parts(i)%fzz=evb%force(3,1)
      End Do
  
    End Do

   End subroutine evb_force


  Subroutine evb_stress(evb,flow,config,stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes for the EVB stress tensor and 
    ! EVB configurational virial 
    ! 
    ! In principle, it is possible to have an EVB decomposition of the virial into separate 
    ! contributions for varios type of interactions (e.g. angles, bonds, dihedrals, 
    ! coulombic, etc). However, i.scivetti has decided not to compute these contributions
    ! as their implementation would entail a significant amount of changes to the code, and
    ! it would be compltely irrelevant to describing the dynamics. Instead, the total 
    ! configurational virial is computed using the trace of the EVB stress tensor
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti September 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    Type(evb_type),           Intent(InOut) :: evb
    Type(flow_type),          Intent(In   ) :: flow
    Type(configuration_type), Intent(InOut) :: config(:)
    Type(stats_type),         Intent(InOut) :: stat(:)
  
    ! Working arrays and variables
    Integer(Kind = wi) :: m,k                ! Indices for matrix elements
    Integer(Kind = wi) :: alpha, beta, gama  ! Indices for coordinates x, y, z

    Real(Kind = wp)    :: invcell(1:9)
    Real(Kind = wp)    :: h(3,3)                 ! Matrix describing the simulation cell 
    Real(Kind = wp)    :: det                    ! Volume of h  
    Real(Kind = wp)    :: invh(3,3)              ! Inverse of h
    Real(Kind = wp)    :: strff(flow%NUM_FF,3,3) ! Temporary working array to compute stress
     
  
    !Arrange stat(m)%stress in a more convenient index notation and store in strff. We do that for every field. 
    Do m=1,flow%NUM_FF 
      Do alpha=1,3
        Do beta=1,3
          strff(m,alpha,beta)=stat(m)%stress(3*(alpha-1)+beta)
        End Do
      End Do
    End Do
  
    !Arrange config(1)%cell in a more convenient index notation, stored in matrix h 
    Do alpha =1, 3
      Do beta =1, 3
        h(alpha,beta)=config(1)%cell(3*(alpha-1)+beta) 
      End Do
    End Do
       
    !Compute the inverse of the lattice vectors   
     Call invert(config(1)%cell,invcell,det)
  
    !Arrange invcell in a more convenient index notation, stored in invh
    Do alpha=1,3
      Do beta=1,3
        invh(alpha,beta)=invcell(3*(alpha-1)+beta) 
      End Do
    End Do
  
    Do alpha=1,3
      Do beta=1,3

      ! Build the alpha,beta component of the evb%dE_dh matrix

        ! Diagonal elements
        Do m=1,flow%NUM_FF
          evb%stress_matrix(m,m)= 0.0_wp
          Do gama=1,3  
            evb%stress_matrix(m,m)= evb%stress_matrix(m,m)+strff(m,alpha,gama)*invh(beta,gama)
          End Do
        End Do
        ! Off-diagonal elements     
        Do m=1,flow%NUM_FF
          Do k=m+1,flow%NUM_FF
            evb%stress_matrix(m,k)=0.0_wp
            Do gama=1,3
              evb%stress_matrix(m,k)=evb%stress_matrix(m,k)+(strff(m,alpha,gama)-strff(k,alpha,gama))*invh(beta,gama)
            End Do
             ! Compute coupling stress for element (k,m) of the matrix
             evb%stress_matrix(m,k) = evb%stress_matrix(m,k)*evb%grad_coupl(m,k)
             !  Make matrix symmetric
            evb%stress_matrix(k,m)=evb%stress_matrix(m,k)
          End Do 
        End Do
        
        ! Initialise dE_dh to zero
        evb%dE_dh(alpha,beta)=0.0_wp
    
        ! Matrix multiplication (evb%psi)^T*(evb%stress_matrix)*(evb%psi) to compute stat(1)%stress(j)
        Do m=1,flow%NUM_FF
          Do k=1,flow%NUM_FF
            evb%dE_dh(alpha,beta)=evb%dE_dh(alpha,beta)+evb%psi(m,1)*evb%stress_matrix(m,k)*evb%psi(k,1)
          End Do
        End Do

      End Do
    End Do 
  
    !Obtain EVB stress tensor
    Do alpha=1,3
      Do beta=alpha,3
       ! Set stress to zero
        evb%stress(alpha,beta)=0.0_wp
        Do gama=1,3
          evb%stress(alpha,beta)=evb%stress(alpha,beta)+evb%dE_dh(alpha,gama)*h(beta,gama)
        End Do
        ! Make stress tensor symmetric
        If(alpha /= beta)Then
          evb%stress(beta,alpha)=evb%stress(alpha,beta)
        End If  
      End Do
    End Do 
  
    !Copy evb%stress to stat(m)%stress for each FIELD
    Do alpha=1,3
      Do beta=1,3
        Do m=1,flow%NUM_FF  
          stat(m)%stress(3*(alpha-1)+beta)=evb%stress(alpha,beta)
        End Do  
      End Do
    End Do
  
  
    Do m=1,flow%NUM_FF
       stat(m)%virtot=-(stat(m)%stress(1)+stat(m)%stress(5)+stat(m)%stress(9))
    End Do

  End subroutine evb_stress


  Subroutine evb_population(evb,flow,files,comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to compute the population of EVB states over time
    ! Results are written in file POPEVB ater equilibration, only if the flag
    ! evbpop is present in the SETEVB file
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti January 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),   Intent(InOut) :: evb
    Type(flow_type),  Intent(In   ) :: flow
    Type(file_type),  Intent(InOut) :: files(:)
    Type(comms_type), Intent(InOut) :: comm

    Integer(Kind = wi) :: ff
    Logical            :: l_tmp
  
    ! open EVB file file and write header
    If (evb%newjob .And. comm%idnode == 0) Then
  
      evb%newjob = .false.
      l_tmp=.false.
      
      ! If the flow%restart_key = RESTART_KEY_OLD is the file old (does it exist)?
      If (flow%restart_key == RESTART_KEY_OLD)Then
        Inquire(File=files(FILE_POPEVB)%filename, Exist=l_tmp)
      End If
  
      If (.not.l_tmp) Then
        Open(Newunit=files(FILE_POPEVB)%unit_no,File=files(FILE_POPEVB)%filename, Status='Replace')
        write(files(FILE_POPEVB)%unit_no,'(a,17x,a)') '# Time (ps)','Weights of each FFs in the total EVB state (FF 1, 2, etc)'
      Else
        Open(Newunit=files(FILE_POPEVB)%unit_no,File=files(FILE_POPEVB)%filename, Position='append')      
      End If  
  
    End If
  
    ! Only print after equilibration
    If(flow%step>flow%equil_steps+1)Then
       If (comm%idnode == 0)Then
          write(files(FILE_POPEVB)%unit_no,*) flow%time, (evb%psi(ff,1)**2, ff=1,flow%NUM_FF) 
       End If
    End If

  End subroutine evb_population


  Subroutine evb_setzero(flow,stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that sets to zero most all the decomposed terms 
    ! of the conformational virial and energy.  
    ! 
    ! copyright - daresbury laboratory
    ! author    - i.scivetti January 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(flow_type),  Intent(In   ) :: flow
    Type(stats_type), Intent(InOut) :: stat(:)
  
    Integer(Kind = wi) :: ff
  
    Do ff=1,flow%NUM_FF
    !Energy components
      stat(ff)%engcpe = 0.0_wp
      stat(ff)%engsrp = 0.0_wp 
      stat(ff)%engter = 0.0_wp 
      stat(ff)%engtbp = 0.0_wp 
      stat(ff)%engfbp = 0.0_wp 
      stat(ff)%engshl = 0.0_wp 
      stat(ff)%engtet = 0.0_wp 
      stat(ff)%engbnd = 0.0_wp 
      stat(ff)%engang = 0.0_wp 
      stat(ff)%engdih = 0.0_wp 
      stat(ff)%enginv = 0.0_wp
    !Virial components
      stat(ff)%vircpe = 0.0_wp 
      stat(ff)%virsrp = 0.0_wp
      stat(ff)%virter = 0.0_wp
      stat(ff)%virtbp = 0.0_wp
      stat(ff)%virfbp = 0.0_wp
      stat(ff)%virshl = 0.0_wp
      stat(ff)%virtet = 0.0_wp 
      stat(ff)%virbnd = 0.0_wp
      stat(ff)%virang = 0.0_wp
      stat(ff)%virdih = 0.0_wp 
      stat(ff)%virinv = 0.0_wp
    End Do 

  End subroutine evb_setzero

  Subroutine evb_merge_stochastic(flow,config,stat,rigid,thermo,cshell,cons,pmf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that copies various fields of variable type of the force-field 1
    ! to the rest of the force fields. This copy is needed only when using stochastic features: 
    !
    ! -the regauss option for equilibration
    ! -the pseudo thermostat for higly non-equilibrium dynamics.
    !
    !  
    ! copyright - daresbury laboratory
    ! author    - i.scivetti January 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(flow_type),          Intent(In   ) :: flow
    Type(configuration_type), Intent(InOut) :: config(:)
    Type(stats_type),         Intent(InOut) :: stat(:)
    Type(rigid_bodies_type),  Intent(InOut) :: rigid(:)
    Type(thermostat_type),    Intent(InOut) :: thermo(:)
    Type(core_shell_type),    Intent(InOut) :: cshell(:)
    Type(constraints_type),   Intent(InOut) :: cons(:)
    Type(pmf_type),           Intent(InOut) :: pmf(:)
  
    Integer(Kind = wi) :: m
  
    Do m=2,flow%NUM_FF
  
      config(m)%vxx=config(1)%vxx
      config(m)%vyy=config(1)%vyy
      config(m)%vzz=config(1)%vzz
      
      config(m)%parts(:)%fxx=config(1)%parts(:)%fxx
      config(m)%parts(:)%fyy=config(1)%parts(:)%fyy
      config(m)%parts(:)%fzz=config(1)%parts(:)%fzz
      
      config(m)%parts(:)%xxx=config(1)%parts(:)%xxx
      config(m)%parts(:)%yyy=config(1)%parts(:)%yyy
      config(m)%parts(:)%zzz=config(1)%parts(:)%zzz
      
      rigid(m)=rigid(1)
      thermo(m)=thermo(1)
      cshell(m)=cshell(1)
      cons(m)=cons(1)
      pmf(m) = pmf(1)
      stat(m)=stat(1)
    
    End Do 

  End subroutine evb_merge_stochastic

End Module evb        
