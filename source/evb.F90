Module evb
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module that declares arrays and variables and defines
  ! subroutines for the computation of the Empirical Valence Bond (EVB) method
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
  Use vdw,                      Only: vdw_type

  Implicit None
  Private

Type, Public ::  evb_type
!> Type of EVB simulation: standard EVB (0). It is set to standard by default. Plan: Free-Energy-Perturbation (=1) and Multi-component EVB (=2).
  Integer(Kind = wi)              :: typsim = 0
!> Number of sites that correspond to the EVB reactive unit (FF dependent)
  Integer(Kind = wi), Allocatable :: num_site(:)
!> Number of atoms of the EVB reactive unit
  Integer(Kind = wi), Allocatable :: num_at(:)
!> Number of types-of-molecules of the FIELD file that describe the EVB site (FF dependent)
  Integer(Kind = wi), Allocatable :: typemols(:)

!> Coupling related variables for computing the coupling between fields
  Integer                         :: maxparam=4         !  maximum number of parameters used for the functional forms of coupling terms
  Character(Len = 5), Allocatable :: typcoupl(:,:)       !  Matrix with the type of functional form (const or gauss) for coupling terms
  Real(Kind = wp),    Allocatable :: coupl_param(:,:,:)  !  Matrix with input parameters for the functional forms for the
                                                         !  computations of the coupling terms

!> Energy shift for force fields, convenient to model asymmetries in the potential energy surface (PES)
  Real(Kind = wp), Allocatable    :: eshift(:)

!> Energy for each force field, including any potential shift (working array to build ene_matrix)
  Real(Kind = wp), Allocatable    :: eneFF(:)

!> EVB Energy matrix
  Real(Kind = wp), Allocatable    :: ene_matrix(:,:)

!> Working arrays for diagonalization
  Real(Kind = wp), Allocatable    :: work(:)
  Integer(Kind = wi), Allocatable :: ifail(:), iwork(:)
!> EVB eigenvalues
  Real(Kind = wp), Allocatable    :: eigval(:)
!> EVB eigenvectors
  Real(Kind = wp), Allocatable    :: psi(:,:)

! Limit for the maximum absolute value of the argument to compute exp()
  Real(Kind = wp)                 :: elimit

!> Matrix for coupling terms
  Real(Kind = wp),    Allocatable :: coupl(:,:)

!> Matrix for the gradient of coupling terms
  Real(Kind = wp),    Allocatable :: grad_coupl(:,:)

!> EVB Force matrix
  Real(Kind = wp), Allocatable    :: force_matrix(:,:)
!> EVB ionic force (working array)
  Real(Kind = wp), Allocatable    :: force(:,:)

!> Matrix for computation of the EVB stress
  Real(Kind = wp), Allocatable    :: stress_matrix(:,:)
!> Matrix with the derivative of the EVB energy with respect to each component of the lattice vectors
  Real(Kind = wp), Allocatable    :: dE_dh (:,:)
!> EVB stress tensor (working array)
  Real(Kind = wp), Allocatable    :: stress(:,:)

!> Number of intramolecular interactions corresponding only to the EVB site
  Integer(Kind = wi), Allocatable :: num_bond(:), num_angle(:), num_dihedral(:)
  Integer(Kind = wi), Allocatable :: num_inversion(:)
!> Flag in case of zero coupling
  Logical                         :: no_coupling = .False.
!> Flag to activate the printing of EVB population
  Logical                         :: population = .False.
!> Flag for opening EVB population file POPEVB
  Logical                         :: population_file_open = .False.
!> Flag for newjob
  Logical                         :: newjob = .True.

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
  Public :: evb_check_vdw
  Public :: evb_merge_stochastic
  Public :: evb_pes
  Public :: evb_population
  Public :: evb_prevent
  Public :: evb_setzero
  Public :: print_evb_banner
Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! EVB subroutines can be divided in five different categories
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
!    - evb_check_vdw
!    - evb_check_configs
!    - evb_check_constraints
!    - evb_check_intermolecular
!    - evb_check_intramolecular
!    - evb_check_external
!    - evb_prevent
!
!    with all the auxiliary subroutines (see bellow)
!
! 4) Computation of EVB energy, ionic forces and stress tensor (and virial)
!    - evb_pes
!    - evb_couplings
!    - evb_energy
!    - evb_force
!    - evb_stress
!    - evb_population
!    - evb_setzero
!
!    with all the auxiliary subroutines (see bellow)
!
! 5) Stochastic dynamics, when either the option regauss of a pseudo thermostat is used
!    - evb_merge_stochastic
!
! The purpose of each subroutine is described below
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Subroutine allocate_evb_arrays(evb, num_ff)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to allocate EVB arrays of evb_type variables
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Class(evb_type), Intent(InOut)     :: evb
    Integer(Kind = wi), Intent (In   ) :: num_ff

    Integer :: fail(1:24)

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
    Allocate(evb%coupl_param(num_ff,num_ff,evb%maxparam), Stat=fail(14))
    Allocate(evb%typcoupl(num_ff,num_ff),      Stat=fail(15))
    Allocate(evb%coupl(num_ff,num_ff),         Stat=fail(16))
    Allocate(evb%grad_coupl(num_ff,num_ff),    Stat=fail(17))
    Allocate(evb%typemols(num_ff),             Stat=fail(18))
    Allocate(evb%num_at(num_ff),               Stat=fail(19))
    Allocate(evb%num_site(num_ff),             Stat=fail(20))
    Allocate(evb%num_bond(num_ff),             Stat=fail(21))
    Allocate(evb%num_angle(num_ff),            Stat=fail(22))
    Allocate(evb%num_dihedral(num_ff),         Stat=fail(23))
    Allocate(evb%num_inversion(num_ff),        Stat=fail(24))

    If (Any(fail /= 0 ))Then
      Call error(0,' error - allocation failure in evb -> allocate_evb_arrays')
    End If

  End Subroutine allocate_evb_arrays

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to De-allocate those EVB related arrays
    ! allocated in subroutine allocate_evb_array
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
    If (Allocated(T%coupl_param)) Then
      Deallocate(T%coupl_param)
    End If
    If (Allocated(T%typcoupl)) Then
      Deallocate(T%typcoupl)
    End If
    If (Allocated(T%coupl)) Then
      Deallocate(T%coupl)
    End If
    If (Allocated(T%grad_coupl)) Then
      Deallocate(T%grad_coupl)
    End If
    If (Allocated(T%typemols)) Then
      Deallocate(T%typemols)
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
    ! from the SETEVB file. Few checks are also performed to verify
    ! the correctness of the input values. If any inconsistency is found,
    ! execution is aborted with an error message. Coupling parameters
    ! and energy shifts are printed to OUTPUT
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
    Call info(' Reading EVB settings from the SETEVB file and checking correctness...',.True.)
    Call info(' ',.True.)

    ! Check if units are the same for all FIELD files
    ! To this purpose, field ffunit was added to site_type (site.F90). Values for sites(m)%ffunit were set
    ! equal to engunit in subroutine read_field when field was read.
    Do m=1, flow%NUM_FF-1
      If(Abs(sites(m)%ffunit-sites(m+1)%ffunit) >= epsilon(sites(m)%ffunit)) Then
        Write(message, '(1x,2(a,i2),a)') 'error - Units used in FF', m, ' **DIFFER** from those specified in FF ', m+1, &
                                         '. Units MUST be the same for all FIELD files'
        Call error(0,message)
      End If
    End Do

    ! Allocate matrix for checking all possible couplings have been set
    Allocate(couplflag(flow%NUM_FF,flow%NUM_FF))
    Allocate(eshiftflag(flow%NUM_FF))

    ! Initialise the number of atoms that belong to the EVB site (calculated below)
    evb%num_at    = 0

    ! Initialise the number of sites that belong to the EVB site
    evb%num_site    = 0

    ! Initialization of flags for control over the reading of EVB input settings
    couplflag  = .False.
    eshiftflag = .False.

    FFmolflag   = .False.
    couplerror  = .False.
    eshifterror = .False.

    !Initialise arrays of coupling parameters and energy shifts
    evb%coupl_param= 0.0_wp
    evb%eshift=0.0_wp

    ! Set the total number of coupling terms to be considered
    ncoupl=(flow%NUM_FF-1)*flow%NUM_FF/2

    ! initialise counters
    icoupl=0
    ieshift=0

    ! Set safe flag for reading
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

        ! For each chemical state, read the number of the type-of-molecules that describe the EVB site.
        ! For example, let us assume a EVB site in vacuum with two possible chemical states. In chemical state 1 the EVB site is a
        ! single molecule (described by one type-of-molecule in the FIELD file). In chemical state 2, the system
        ! breaks in two fragments described by 2 types-of-molecules (FIELD2 file).
        ! Thence, the specification for evbtypemols in the SETEVB file SHOULD be:
        !
        ! evbtypemols    1     2
        !
        ! The number of values specified after the option "evbtypemols" MUST be equal to flow%NUM_FF, otherwise the run stops.
        ! Values for evbtypemols should not be larger than the corresponding value for option "molecules" specified in each FIELD file.
        !
        ! By construction, any other type-of-molecule describing the non-EVB part MUST be specified in the FIELD file ONLY after
        ! the specification of all "evbtypemols" type-of-molecules that describe the EVB site. To exemplify this point, lets
        ! consider the EVB site already described and assume a set of N_h2o surrounding non-reactive water molecules, described by
        ! a water-type-of-molecule (with nummols n_h2o). The specification for "evbtypemols" should be the same as above.
        !
        Else If (word(1:11) == 'evbtypemols') Then
          FFmolflag=.True.
          Do k=1, flow%NUM_FF
            Call get_word(record,word)
            If (word(1:2) == '  ') Then
              Write(message,'(1x,2a,i2,a)') 'error - In file SETEVB, incomplete line after the "evbtypemols" key. ', &
                                            'Value of evbtypemols for FF', k, ' is missing. See manual for correct syntax '
              Call error(0,message)
            Else
              evb%typemols(k)= Nint(word_2_real(word,0.0_wp))

              ! Check the value makes sense, otherwise print error and abort
              If( evb%typemols(k) <= 0 .Or. evb%typemols(k) > sites(k)%ntype_mol)Then
                Write(message,'(1x,2(a,i2),2a,i2)') 'error - In file SETEVB, wrong input for index', k,&
                                                    ' of the evbtypemols list (corresponding to FF', k, &
                                                    '). This value should be > 1 and <= MOLECULES, ', &
                                                    'where MOLECULES for this FF is equal to ', sites(k)%ntype_mol
                Call error(0,message)
              End If
            End If
          End Do

        ! Read all coupling elements
        ! All the possible ncoupl pairs of couplings MUST be specified, otherwise the run stops (see counter icoupl)
        ! The syntax for the spectication of each coupling term is the following, depending on the functional type
        !
        !             "Force Field pair of indexes"    "Type of interaction"   "Parameters (in the units of the FF)"
        ! evbcoupl          m            k                     const                   A1
        ! evbcoupl          m            k                     gauss                   A1   A2   A3   A4
        !
        ! If specification has misspellings, or there is a wrong syntax, the simulation stops (hopefully) with
        ! the error printed
        !
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
            Write(message,'(1x,2(a,i2),a)') 'error - In SETEVB file, the coupling between field ',i, ' and ', j, &
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
                evb%coupl_param(i,j,k)= word_2_real(word)
                evb%coupl_param(j,i,k)=evb%coupl_param(i,j,k)
              End If
            End Do
          Else If (evb%typcoupl(i,j) == 'gauss') Then
            kparam=evb%maxparam
            Do k=1, kparam
              Call get_word(record,word)
              If (word(1:3) == '   ' .Or. word(1:1) == '!' .Or. word(1:1) == '#') Then
                couplerror=.True.
              Else
                evb%coupl_param(i,j,k)= word_2_real(word)
                evb%coupl_param(j,i,k)= evb%coupl_param(i,j,k)
              End If
            End Do
            If(Abs(evb%coupl_param(i,j,3)) < epsilon(evb%coupl_param(i,j,3)) )Then
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
        ! Force field energies can be shifted by a given amount of energy, specified via evbshift.
        ! All energy shifts should be specified, even though shifts are zero. The syntax is:
        !
        ! evbshift      Field       DeltaE
        !
        ! If there are mispellings, or missing input the run stops.

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
              Write(message,'(1x,a,i2,a)') 'error - In SETEVB file, energy shift for field ',i,' cannot be defined more than once'
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
      Write(message,'(1x,a)') 'error - In file SETEVB there are missing specification for "evbshift". All fields &
                             &must be specified'
      Call error(0,message)
    End If
    If(.Not. FFmolflag)Then
      Write(message,'(1x,a)') 'error - Flag evbtypemols not found in the SETEVB file.'
      Call error(0,message)
    End If

    Deallocate(eshiftflag, couplflag)

    ! Close SETEVB file
    If (comm%idnode == 0)Then
      Close(files(FILE_SETEVB)%unit_no)
    End If

    ! For standard EVB (evb%tysim=0), only one reactive EVB site is allowed.
    ! Such a reactive site can be well composed of one or many fragments, depending on the chemical state.
    ! However, the nummols for the first evbtypemols type-of-molecules MUST be equal to 1.
    If(evb%typsim == 0)Then
      Do m=1, flow%NUM_FF
        Do i=1, evb%typemols(m)
         If(sites(m)%num_mols(i) /= 1 )Then
           Write(messages(1),'(1x,a)') 'error - For standard EVB simulations, only a single EVB reactive site is allowed'
           Write(messages(2),'(1x,a)') 'Thus, NUMMOL MUST be equal to 1 for the first evbtypemols types-of-molecules'
           Write(messages(3),'(1x,a)') 'that define the EVB site (evbtypemols is set in the SETEVB file)'
           Call info(messages,3,.True.)
           Call error(0)
         End If
        End Do
      End Do
    End If


    ! Calculate the total number of atoms being part of the EVB reactive site
    Do m=1, flow%NUM_FF
      Do i=1, evb%typemols(m)
        evb%num_at(m)=evb%num_at(m)+sites(m)%num_mols(i)*sites(m)%num_site(i)
      End Do
    End Do


    ! Calculate the total sites being part of the EVB reactive site(s)
    Do m=1, flow%NUM_FF
      Do i=1, evb%typemols(m)
        evb%num_site(m)=evb%num_site(m)+sites(m)%num_site(i)
      End Do
    End Do

    ! Check if the number of EVB atoms is the same for all FFs. Print an error message and about otherwise.
    Do m=1, flow%NUM_FF-1
      If(evb%num_at(m) /= evb%num_at(m+1))Then
        Write(messages(1),'(1x,2(a,i2))') 'error - The total number of EVB atoms differ between FF', m, ' and FF ', m+1
        Write(messages(2),'(1x,a)')       'This means that the number of EVB atoms is different for different FIELD files'
        Write(messages(3),'(1x,a)')       'ACTION: check the number of EVB atoms, as well as the settings for evbtypemols '
        Write(messages(4),'(1x,a)')       'in SETEVB and the structure of type-of-molecules in the FIELD files'
        Call info(messages,4,.True.)
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
    Write(messages(2),'(1x,a)') 'The adopted functional form for the coupling between fields m and k (C_{mk})'
    Write(messages(3),'(1x,a)') 'can be of two types: constant (const) or Gaussian (gauss)'
    Write(messages(4),'(1x,a)') '(See manual for details of the Gaussian functional form)'
    Write(messages(5),'(1x,a)') 'Details for coupling terms are summarised in the following table:'
    Call info(messages,5,.True.)
    Write(messages(1),'(1x,a)')               '______________________________________________________________________&
                                             &_____________________________'
    Write(messages(2),'(1x,a,5x,a,5x,a)')           'Force Field pair' , 'Type' , 'Parameters ('//trim(evbunit)//')'
    Write(messages(3),'(5x,a,5x,a,10x,a,4(15x,a))')             'm','k', '    '  , 'A1', 'A2','A3','A4'
    Write(messages(4),'(1x,a)')               '______________________________________________________________________&
                                             &_____________________________'
    Call info(messages,4,.True.)
    Do i=1,flow%NUM_FF
     Do j=i+1,flow%NUM_FF
       If(evb%typcoupl(i,j)=='const')Then
         Write(message,'(2(4x,i2),10x,a,10x,1(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%coupl_param(i,j,k),k=1,1)
         Call info(message,.True.)
       ElseIf(evb%typcoupl(i,j)=='gauss')Then
         Write(message,'(2(4x,i2),10x,a,10x,4(E12.5,5x))') i , j, evb%typcoupl(i,j), (evb%coupl_param(i,j,k),k=1,4)
         Call info(message,.True.)
       End If
     End Do
    EnD Do
    Write(message,'(1x,a)')               '______________________________________________________________________&
                                             &_____________________________'
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
    evb%coupl_param = engunit*evb%coupl_param

    ! If There is no coupling defined, e.i. all type of coupling are constants and equal to zero
    ! set the flag evb%no_coupling = .True. and this avoids the execution of evb_setzero in
    ! calculate_forces (drivers.F90)
    icoupl=0

    Do i = 1, flow%NUM_FF
       Do j = i+1, flow%NUM_FF
         If( (abs(evb%coupl_param(i,j,1)) <= epsilon(evb%coupl_param(i,j,1))) .And. &
                  evb%typcoupl(i,j) == 'const'  )Then
            icoupl = icoupl + 1
         End If
       End Do
    End Do

    If( icoupl == ncoupl )Then
      evb%no_coupling=.True.
    End If

    ! Set the limit for the maximum value of an argument to be computed by the exponential function.
    evb%elimit=700.0_wp

  End Subroutine read_evb_settings


  Subroutine evb_prevent(lplumed, lkim, lttm, thermo_dpd)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that aborts the execution of a standard EVB runs
    ! if the following are activated
    ! - plumed
    ! - kim
    ! - dpd
    ! - ttm
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti August 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,     Intent(In   ) :: lplumed
    Logical,     Intent(In   ) :: lkim
    Logical,     Intent(In   ) :: lttm
    Integer,     Intent(In   ) :: thermo_dpd

    Logical              :: flag
    Character(Len = 256) :: appex, base
    Character(Len = 256) :: message

    flag=.False.
    base = 'error - Standard EVB formalism is not valid for '


    If(lplumed)Then
      flag=.True.
      appex='Plumed'
    Else If(lkim)Then
      flag=.True.
      appex='KIM'
    Else If(lttm)Then
      flag=.True.
      appex='TTM'
   Else If(thermo_dpd/=0)Then
      flag=.True.
      appex='DPD'
   End If

    If(flag)Then
      Call info(' ', .True.)
      Write(message, '(1x,a,1x,2a)') Trim(base), Trim(appex), ' simulations. Please reconsider the settings'
      Call error(0,message)
    End If

  End Subroutine evb_prevent


  Subroutine evb_check_intrinsic(evb, sites, config, flow, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the intrinsic properties for
    ! each atom specified in the FIELD files. Since atomic labels and charges
    ! might change for the EVB part between different FFs (different chemistry),
    ! intrisic checking is set as follows:
    !
    ! - labels and charges are checked to be the same only for the non-EVB part
    ! - masses are checked to be the same for all atoms in the system, as there
    !   is mass conservation despite of the reactive process
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
    Call info(' intrinsic properties of atoms...',.True.)

    ! Initialise flags
    lmass=.False.
    lchg=.False.
    ltag=.False.

    Do m = 1, flow%NUM_FF-1
      loopi=.True.
      i=1
      Do While ( i <= config(m)%natms .And. loopi)
        ! Check masses for all atoms
        If(Abs(config(m)%weight(i)-config(m+1)%weight(i)) > &
           epsilon(config(m)%weight(i)))Then
          labelint= 'Mass'
          lmass=.True.
        End If
        ! For all atoms of the non-EVB part
        If(config(m)%ltg(i) > evb%num_at(m))Then
          ! Check labels
          If(config(m)%ltg(i) /= config(m+1)%ltg(i))Then
            ltag = .True.
            labelint= 'Name'
          ! Check charges
          Else If(Abs(config(m)%parts(i)%chge - config(m+1)%parts(i)%chge) > &
                  epsilon(config(m)%parts(i)%chge))Then
            labelint= 'Charge'
            lchg=.True.
          End If
        End If
        ! If an inconsistency was found, mark the atomic site and type-of-molecule
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
      ! Due to domain decomposition, inconsistencies do not need to manifest for all processors.
      ! Subroutine evb_intrinsic_error checks if an error was found for any processors.
      Call evb_intrinsic_error(labelint, m, comm, lmass, ltag, lchg, st1num, st2num, st1lab, st2lab, mol)
    End Do

  End Subroutine evb_check_intrinsic


  Subroutine evb_intrinsic_error(string, field, comm, lmass, ltag, lchg, st1num, st2num, st1lab, st2lab, mol)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intrinsic) to print an
    ! error message and abort if there was an inconsitency found.
    !
    ! Communication is needed to detect for errors in all processors.
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

    ! These are common error messages
    Write(messages(3),'(1x,3a)')     'ACTION: check settings in FIELD files. If, in principle, the above element numbers are &
                                     &inconsistent with'
    Write(messages(4),'(1x,3a)')     'the amount of atomic sites listed for this type-of-molecule in the FIELD file, &
                                     &please consider the specification'
    Write(messages(5),'(1x,3a)')     'for the repetiton number of each site (if present, repetitions are set following the &
                                     &values for charges)'

    If(comm%idnode == 0)Then
      Do jdnode=0,comm%mxnode-1
        If(jdnode>0)Then
          ! In node 0, receive the following variables from the rest of the nodes
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
                                                   st1num(1), ' element of type-of-molecule ', mol(1), ' in the FF ', field
          siterror=.True.
          If(lmass .Or. lchg)Then
            ! Error from inconsistency in mass or charge
            Write(messages(2),'(1x,a,i6,2(a,i2))')  '**DIFFERS** from the value assigned to the equivalent site: element ', &
                                                    st2num(1), ' in type-of-molecule', mol(2), ' of FF ', field+1
          Else If(ltag)Then
            ! Error from inconsistency in the labelling
            Write(messages(2),'(1x,4a,i6,2(a,i2))') '**DIFFERS** from the name ', trim(st2lab(1)), 'assigned to the equivalent ', &
                                                    'site: element ', st2num(1), ' in type-of-molecule ', mol(2), ' of FF ', field+1
          End If
        End If

        If(siterror)Then
          Call info(messages,5,.True.)
          Call error(0)
        End if
      End Do
    Else
      ! If the node is not node 0, send the following variables to node 0
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


  Subroutine evb_check_configs(config, flow, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check all CONFIG files (one CONFIG file per FF)
    ! have the same:
    !
    ! - number of atoms
    ! - coordinates
    ! - simulation cell dimensions
    ! - symmetry (image convention)
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(configuration_type), Intent(In   ) :: config(:)
    Type(flow_type),          Intent(In   ) :: flow
    Type(comms_type),         Intent(InOut) :: comm


    Integer                :: alpha, beta
    Integer                :: m, i
    Character(Len = 256)   :: message
    Character(Len = 256)   :: messages(3)
    Character(Len = 1)     :: coord, coord0
    Character(Len = 6)     :: string1, string2
    Logical                :: carry

    Integer                :: jdnode, stnum

    Real(Kind = wp)        :: cell(flow%NUM_FF,3,3)

    Call info(' ',.True.)
    Call info(' atomic coordinates, symmetry and dimensions for the simulation cell...',.True.)

    ! Comparison between different CONFIG files: check the value imcon (image convention/symmetry)
    Do m = 1, flow%NUM_FF-1
      If(config(m)%imcon /= config(m+1)%imcon)then
        If(m == 1) Then
          Write(string1,'(a1)') ' '
        Else
          Write(string1,'(i2)') m
        End If
        Write(string2,'(i2)') m+1
        ! In case an inconsistency has been found, complain and abort
        Write(message,'(1x,5a)')  'error - Value for imcon *DIFFERS* between ', &
                                  'CONFIG'//adjustl(trim(string1)), ' and ','CONFIG'//adjustl(trim(string2)), ' files'
        Call error(0,message)
      End If
    End Do

    ! Comparison between different CONFIG files: check number of atoms
    Do m = 1, flow%NUM_FF-1
      If(config(m)%megatm /= config(m+1)%megatm)Then
        Write(messages(1),'(1x,2(a,i2))') 'error - Different number of atoms for FFs', m, ' and', m+1
        Write(messages(2),'(1x,a)')       'ACTION: check settings for all FIELD/CONFIG files and make sure &
                                          &the number of atoms is the same'
        Call info(messages, 2, .True.)
        Call error(0)
      End If
    End Do

    ! Comparison between different CONFIG files: check simulation cell dimensions
    ! To this purpose, arrange config(m)%cell in a more convenient index notation, stored in cell
    Do m = 1, flow%NUM_FF
      Do alpha=1,3
        Do beta=1,3
          cell(m,alpha,beta)=config(m)%cell(3*(alpha-1)+beta)
        End Do
      End Do
    End Do

    ! Compare cell dimensions
    Do alpha = 1, 3
      Do beta = 1, 3
        Do m = 1, flow%NUM_FF-1
          If(abs(cell(m,alpha,beta)-cell(m+1,alpha,beta))>= epsilon(cell(m,alpha,beta)))then
           If(m == 1) Then
             Write(string1,'(a1)') ' '
           Else
             Write(string1,'(i2)') m
           End If
           Write(string2,'(i2)') m+1
           ! In case an inconsistency has been found, complain and abort
           Write(messages(1),'(1x,a,2i2,a)') 'error - Component (', alpha, beta,') of the simulation cell differs between'
           Write(messages(2),'(1x,4a)')      'CONFIG'//adjustl(trim(string1)), ' and ','CONFIG'//adjustl(trim(string2)), ' files'
           Call info(messages,2,.True.)
           Call error(0)
          End If
        End Do
      End Do
    End Do

    ! Comparison of ionic positions between all CONFIG files
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
               Write(string1,'(a1)') ' '
             Else
               Write(string1,'(i2)') m
             End If
             Write(string2,'(i2)') m+1
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
           Write(messages(1),'(1x,3a,i7,6a)')    'error -  ', coord, '-coordinate for atom', stnum, &
                                                 ' **DIFFERS** between ', 'CONFIG'//adjustl(trim(string1)),&
                                                 ' and ', 'CONFIG'//adjustl(trim(string2))
           Write(messages(2),'(1x,a)')           'ACTION: Fix the error. Atomic coordinates should not vary between CONFIG files'
           Call info(messages,2,.True.)
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

  End Subroutine evb_check_configs


  Subroutine evb_check_constraints(evb, config, cons, cshell, tether, sites, flow, rigid, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the definition of constraints
    ! between atoms (constraints are defined in each FIELD file).
    ! Check is carried out over frozen atoms, rigid bond constraint, rigid bodies,
    ! core shells and tether sites.
    !
    ! Whatever constraints one defines, they MUST NOT differ between FIELD files.
    ! In other words, a group of atoms cannot be a rigid body for one FIELD
    ! and a flexible molecule for another FIELD. The same logic applies to other
    ! constraints. Even though the EVB could be computed independently, different
    ! contraint settings will lead to a dynamics that depends on the field, which
    ! is conceptually wrong.
    !
    ! To check consistency between all FFs, the process starts by comparing FF1 and FF2,
    ! continues with FF2 with FF3, FF3 with FF4,....., and finishes FF(NFF-1) with FF(NFF).
    !
    ! At first sight, one might consider there is a lot of code repetition below.
    ! However, this is a consequence of having different formats for the specification
    ! of different constraint. There might be ways to improve/simplify the coding.
    ! This is left for further work.
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),           Intent(In   ) :: evb
    Type(configuration_type), Intent(In   ) :: config(:)
    Type(constraints_type),   Intent(In   ) :: cons(:)
    Type(core_shell_type),    Intent(In   ) :: cshell(:)
    Type(tethers_type),       Intent(InOut) :: tether(:)
    Type(site_type),          Intent(In   ) :: sites(:)
    Type(flow_type),          Intent(In   ) :: flow
    Type(rigid_bodies_type),  Intent(In   ) :: rigid(:)
    Type(comms_type),         Intent(InOut) :: comm

    Integer                :: m, i, j, l, k
    Integer                :: field
    Integer                :: mol1, mol2

    Character(Len = 256)   :: labelint

    Logical                :: loopi, loopj, loopk
    Logical                :: siteevb, sitemiss, siteparam
    Integer                :: listdim


    Integer, Allocatable   :: st1num(:), st2num(:)
    Integer, Allocatable   :: list(:)

    Character(Len =  8)    :: st1lab(2), st2lab(2)
    Integer                :: mol(2)

    Real(Kind = wp)        :: param(2)

    Call info(' ',.True.)
    Call info(' constraints (if any)...',.True.)

    ! Initialise flags
    siteevb   = .False.   ! To check if the constraint unit is part of the EVB site (rigid bodies are not allowed)
    sitemiss  = .False.   ! To check if a given constraint unit is defined in all FFs
    siteparam = .False.   ! To check if the parameters parameters with the constrained
                          ! unit is consistently defined for all FFs


    ! Check consistency in the total number of rigid bodies
    ! (there must be the same number of rigid bodies defined for all chemical states)

    ! Loop over FFs
    Do m = 1, flow%NUM_FF-1
      If(rigid(m)%total /= rigid(m+1)%total)Then
        labelint='Total number of rigid body units'
        Call evb_constraint_error(labelint, m)
      End If
    End Do

    ! Define the mximum dimension for the list of atoms (listdim) that are part of the constraints.
    ! All the constraints but rigid bodies contain two are less atoms per constraint. If
    ! rigid bodies are defined, then set listdim to the maximum number of atoms found for a rigid body.

    If(rigid(1)%total /= 0 )Then
      listdim = rigid(1)%max_list
      ! Check all fields have exaclty the same number of maximum number of atoms for a rigid body unit.
      ! If not, print an error message and abort.
      Do m = 1, flow%NUM_FF-1
        If(rigid(m)%max_list /= rigid(m+1)%max_list)Then
          labelint='maximum number of atoms for a rigid body'
          Call evb_constraint_error(labelint, m)
        End If
      End Do
    Else
      ! If no rigid body is present, set listdim to 2, which is the dimension for the list
      ! of atoms corresponding to any constraint except rigid bodies
      listdim=2
    End If

    ! Allocate the following arrays according to the computed value for listdim
    Allocate(st1num(listdim), st2num(listdim),list(listdim))

    !Initilization of vectors
    mol=0
    st1num=0
    st2num=0
    list=0

    !! Rigid bodies
    !!
    !! Check consistency in the definition rigid bodies, if defined.
    !!
    !! DL_POLY will complain and abort the execution if a rigid body unit:
    !!
    !! - is not found in all FIELD files
    !! - does not contain the same number of atoms in all FIELD files
    !! - is part of the EVB site.
    !!
    !! Specification for rigid bodies should follow the same sequence for all FIELD files.
    !! If the specification ordering is altered, the run stops.

    If(rigid(1)%total /= 0 )Then
      ! Loop over FFs
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        ! Loop over rigid body units
        Do While ( i <= rigid(m)%max_rigid .And. loopi)
          ! For this constraint, identify constituent sites and the type-of-molecule in FF(m)
          Do l = 1, rigid(m)%list(-1,i)
            list(l) = rigid(m)%list(l,i)
          End Do
          Call obtain_sites_from_list(rigid(m)%list(-1,i),m,sites,list,st1num,mol1)
          mol(1)=mol1
          ! For this constraint, identify constituent sites and the type-of-molecule in FF(m+1)
          Do l = 1, rigid(m+1)%list(-1,i)
            list(l) = rigid(m+1)%list(l,i)
          End Do
          Call obtain_sites_from_list(rigid(m+1)%list(-1,i),m+1,sites,list,st2num,mol2)
          mol(2)=mol2

          ! Check if the rigid body is part of the EVB site for both FFs
          If(mol1 <= evb%typemols(m) .Or. mol2 <= evb%typemols(m+1))Then
            listdim=rigid(m)%list(-1,i)
            siteevb=.True.
            loopi=.False.
            field=m
            If(mol2 <= evb%typemols(m+1))Then
              st1num = st2num
              mol1  = mol2
              field = m+1
            End If
          Else
            ! If rigid body is not part of the EVB site, continue with checking.
            ! Does the rigid body contain the same amount of atoms in both FFs?
            If(rigid(m)%list(-1,i) /= rigid(m+1)%list(-1,i))Then
              listdim=rigid(m)%list(-1,i)
              siteparam=.True.
              loopi=.False.
              field = m
            Else
              loopk=.True.
              k=1
              ! Is the rigid body composed of the same atoms for both FFs?
              Do While ( k <= rigid(m)%list(-1,i) .And. loopk)
                If(rigid(m)%list(k,i) /= rigid(m+1)%list(k,i))Then
                  listdim=rigid(m)%list(-1,i)
                  sitemiss=.True.
                  loopi=.False.
                  loopk=.False.
                  field = m
                End If
                k=k+1
              End Do
            End If
          End If
          i=i+1
        End Do
        ! Call for checking if there was an error for any processor
        labelint='Rigid-body'
        Call evb_constraint_error(labelint, field, comm, st1num, st2num, st1lab, st2lab, mol, listdim, sitemiss, siteparam, siteevb)
      End Do
    End If

    !! Frozen ions
    !!
    !! Check if ions are frozen or not.
    !! A frozen ion should be frozen in all FIELD files.
    !! In principle one might have a frozen ion in the EVB part, as long as
    !! NOT all EVB atoms are frozen, obviously! This might be convenient for
    !! equilibrating large molecules with a rather small EVB region

    ! Loop over FFs
    Do m = 1, flow%NUM_FF-1
      loopi=.True.
      i=1
      ! Loop over the free/frozen nature of units
      Do While ( i <= config(m)%natms .And. loopi)
      ! If an ion is frozen in one FFs but free in the other, identify the corresponding site and type-of molecule. Exit loop
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
      ! Call for checking if there was an error for any processor
      labelint='Frozen'
      Call evb_constraint_error(labelint, m, comm, st1num, st2num, st1lab, st2lab, mol,1, sitemiss, siteparam, siteevb)
    End Do


    !! Bond constraints
    !!
    !! Check consistency in the definition of these units, only if bond constraints are defined.
    !!
    !! DL_POLY will complain and abort the execution if a bond constraint unit:
    !!
    !! - is not found in all FIELD files
    !! - is not composed of the same atoms
    !! - does not have the same constraint distance in all FIELD files
    !!
    !! In contrast to rigid bodies, the sequence of definition can be altered between FIELD files.
    !! Bond constraints can be part of the EVB site. This might be convenient to descibe a large
    !! molecule with a relatively small EVB region.


    ! Check that the number of bonds constraints is the same for all FIELD files
    ! Loop over FFs
    Do m = 1, flow%NUM_FF-1
      If( cons(m)%megcon /= cons(m+1)%megcon )Then
        labelint='Total number of bond-constraints'
        Call evb_constraint_error(labelint, m)
      End If
    End Do

    ! Check that bond constraint specification does not change beetween files
    ! Perform the checking only if the number of bond contraints is different from zero
    If(cons(1)%megcon /= 0)Then

      ! Loop over FFs
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        ! Loop over bond-constraints units of FF(m)
        Do While (i <= cons(m)%ntcons .And. loopi)
          loopj=.True.
          j=1
          ! Loop over bond-constraints units of FF(m+1)
          Do While (j <= cons(m+1)%ntcons .And. loopj)
            If((cons(m)%listcon(1,i) == cons(m+1)%listcon(1,j)) .Or. (cons(m)%listcon(1,i) == cons(m+1)%listcon(2,j)) )Then
              If((cons(m)%listcon(2,i) == cons(m+1)%listcon(1,j)) .Or. (cons(m)%listcon(2,i) == cons(m+1)%listcon(2,j)) )Then

                ! If the rigid bond unit is found in FF-(m+1), check if bond constraint distances are consistent
                ! If distances are different, identify sites and type-of-molecules and exit loop.
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

          ! If rigid bond unit is not found in FF-(m+1), identify sites and type-of-molecules.
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
        ! Call for checking if there was an error for any processor
        labelint='Bond-constraints'
        Call evb_constraint_error(labelint, m, comm, st1num, st2num, st1lab, st2lab, mol, listdim, sitemiss, siteparam, siteevb)
      End Do
    End If

    !! Core-shells
    !!
    !! Check consistency in the definition of these units, only if core-shells are defined.
    !!
    !! DL_POLY will complain and abort the execution if a core-shell unit:
    !!
    !! - is not found in all FIELD files
    !! - is described with different parameters in different FFs
    !!
    !! The sequence of definition of core-shells can be altered between FIELD files.
    !! Core-shells can be part of the EVB region.

    ! Check that the number of core-shell units is the same for all FIELD files
    ! Loop over FFs
    Do m = 1, flow%NUM_FF-1
      If( cshell(m)%megshl /= cshell(m+1)%megshl )Then
        labelint='Total number of core-shell units'
        Call evb_constraint_error(labelint, m)
      End If
    End Do

    ! Perform the checking only if the number of core-shell units is different from zero
    If(cshell(1)%megshl /= 0)Then
      ! Loop over FFs
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        ! Loop over core-shell units of FF(m)
        Do While (i <= cshell(m)%ntshl .And. loopi)
          loopj=.True.
          j=1

          ! Loop over core-shell units of FF(m+1)
          Do While (j <= cshell(m+1)%ntshl .And. loopj)
            If((cshell(m)%listshl(1,i) == cshell(m+1)%listshl(1,j)) .Or. &
               (cshell(m)%listshl(1,i) == cshell(m+1)%listshl(2,j)) )Then
              If((cshell(m)%listshl(2,i) == cshell(m+1)%listshl(1,j)) .Or. &
                 (cshell(m)%listshl(2,i) == cshell(m+1)%listshl(2,j)) )Then

                ! If the core-shell defined in FF(m) was found in FF(m+1), check they are both
                ! defined with the same parameters and leave loop j.
                loopj=.False.
                Do l=1,2
                  param(1)=cshell(m)%prmshl(l,cshell(m)%listshl(0,i))
                  param(2)=cshell(m+1)%prmshl(l,cshell(m+1)%listshl(0,j))

                  ! If assigned core-shell parameters change between FFs, identify
                  ! constituent sites and type-of-molecule. Exit loop i.
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

          ! If core-shell unit of FF(m) was not found in FF(m+1), identify atomic sites and type-of-molecule
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

      ! Call for checking if there was an error for any processor
      labelint='Core-shell'
      Call evb_constraint_error(labelint, m, comm, st1num, st2num, st1lab, st2lab, mol, listdim, sitemiss, siteparam, siteevb)

      End Do

    End If

    !! Tether units
    !!
    !! Check consistency in the definition of these units, only if they are defined.
    !!
    !! DL_POLY will complain and abort the execution if a tether unit:
    !!
    !! - is not found in all FIELD files
    !! - is described with different parameters in different FFs
    !!
    !! The sequence of definition for tether units can be altered between FIELD files.
    !! Tehter can be part of the EVB region.


    ! Compute the total number of tethered sites
    ! Loop over FFs
    Do m = 1, flow%NUM_FF
      tether(m)%megteth=tether(m)%ntteth
     Call gsum(comm, tether(m)%megteth)
    End Do

    ! Check that the number of tether units is the same for all FIELD files
    ! Loop over FFs
    Do m = 1, flow%NUM_FF-1
      If( tether(m)%megteth /= tether(m+1)%megteth )Then
        labelint='Total number of tether units'
        Call evb_constraint_error(labelint, m)
      End If
    End Do

    ! Perform the checking only if the number of core-shell units is different from zero
    If(tether(1)%megteth /= 0)Then
    ! Initialise flags
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        ! Loop over tether units of FF(m)
        Do While (i <= tether(m)%ntteth .And. loopi)
          loopj=.True.
          j=1
          ! Loop over tether units of FF(m+1)
          Do While (j <= tether(m+1)%ntteth .And. loopj)
            If((tether(m)%listtet(1,i) == tether(m+1)%listtet(1,j)))Then
              loopj=.False.

              ! If the tether defined in FF(m) was found in FF(m+1), check they are both
              ! defined with the same parameters and leave loop j.
              ! If assigned core-shell parameters and/or interaction key change between FFs, identify
              ! constituent sites and type-of-molecule. Exit loop i.

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

          ! If tether unit of FF(m) was not found in FF(m+1), identify atomic sites and type-of-molecule
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            list(1)=tether(m)%listtet(1,i)
            Call obtain_sites_from_list(1,m,sites,list,st1num,mol1)
            mol(1)=mol1
          End If
          i=i+1
        End Do

        ! Call for checking if there was an error for any processor
        labelint='Tether'
        Call evb_constraint_error(labelint, m, comm, st1num, st2num, st1lab, st2lab, mol, listdim, sitemiss, siteparam)
      End Do
    End If

    ! Dealocate working arrays
    Deallocate(list, st2num, st1num)

  End Subroutine evb_check_constraints


  Subroutine obtain_sites_from_list(dimlist, field, sites, list, st, mol)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_constraints and evb_check_intrinsic)
    ! to obtain the atomic site(s) and corresponding type-of-molecule given the index(es)
    ! assigned to the atom(s) via the input list.
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

    ! The input list has dimension dimlist (single atom, pair of atoms, rigid bodies)
    Do i=1,dimlist
      flag=.True.
      atnum=list(i)

      ! Search backwards
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

      ! If the above search failed it means that the atoms belongs to a site of the firts type-of-molecule
      If(flag)Then
        num=atnum
        den=sites(field)%num_site(1)
        st(i)=mod(num,den)
        If(st(i)==0)Then
          st(i)=sites(field)%num_site(1)
        End If
        mol=1
      End If

    End Do

  End Subroutine  obtain_sites_from_list


  Subroutine evb_constraint_error(string, field, comm, st1num, st2num, st1lab, st2lab, mol, listdim, sitemiss, siteparam, siteevb)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_constraints) to print contraint related errors.
    ! Constraints are rigid bodies, tethers, bond-constraints, core shells and frozen ions. This
    ! subroutine to print an error message and abort if either
    !
    ! - the number for a given type of constraint is different between FIELD files
    !   (no optional variable present when subroutine is called)
    ! - the constraint unit is defined in one FIELD but not found in other (via input sitemiss=.True.)
    ! - the specification for the constraint unit changes between FIELD files (via input siteparam=.True.)
    ! - a rigid unit has been found in the EVB region (via siteevb=.True.)
    !
    ! Error messages depend on the constraint, which is identified via the input string
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
    Integer,               Intent(InOut), Optional :: st1num(:)
    Integer,               Intent(InOut), Optional :: st2num(:)
    Character(Len = 8 ),   Intent(InOut), Optional :: st1lab(:)
    Character(Len = 8 ),   Intent(InOut), Optional :: st2lab(:)
    Integer,               Intent(InOut), Optional :: mol(:)
    Integer,               Intent(In   ), Optional :: listdim
    Logical,               Intent(InOut), Optional :: sitemiss
    Logical,               Intent(InOut), Optional :: siteparam
    Logical,               Intent(InOut), Optional :: siteevb

    Character(Len = 256)   :: messages(4)
    Character(Len = 256)   :: message1, message2
    Integer                :: jdnode, option
    Logical                :: siterror
    Integer                :: i

    siterror = .False.

    Write(messages(1),'(1x,a)')           'error -'
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
            Call grecv(comm, siteevb   , jdnode, WriteConf_tag)
            Call grecv(comm, sitemiss  , jdnode, WriteConf_tag)
            Call grecv(comm, siteparam , jdnode, WriteConf_tag)
            Call grecv(comm, st1num    , jdnode, WriteConf_tag)
            Call grecv(comm, st2num    , jdnode, WriteConf_tag)
            Call grecv(comm, st1lab    , jdnode, WriteConf_tag)
            Call grecv(comm, st2lab    , jdnode, WriteConf_tag)
            Call grecv(comm, mol       , jdnode, WriteConf_tag)
          End If

          If(siteevb)Then
            If(option == 3)Then
              Write(message1,*)                              trim(string),' unit composed of atomic sites', &
                                                             (st1num(i), i = 1, listdim)
              Write(message2,'(3x,2(a,i3),a)')              'is found in type-of-molecule ', mol(1), ' (FF ', field, ').'

              messages(2) = trim(message1)//trim(message2)
            End If

            Write(messages(3),'(1x,3a)')                    'No ', trim(string), ' MUST BE defined for the EVB site (see values for&
                                                            & evbtypemols in the SETEVB file)'
            siterror=.True.

          Else If(sitemiss) Then
            If(option == 0)Then
              Write(messages(2),'(1x,2a,2(i4,a),2(i2,a))')   trim(string),' unit for atomic site', st1num(1), &
                                                             ' of type-of-molecule ', mol(1), ' (set in  FF ', field, ')'
            Else If(option == 2)Then
              Write(messages(2),'(1x,a,2(a,i4),2(a,i2),a)')  trim(string),' unit between atomic sites', st1num(1), ' and ',&
                                                             st1num(2), ' of type-of-molecule ', mol(1), ' (set in  FF ', field, ')'
            Else If(option == 3)Then
              Write(message1,*)                              trim(string),' unit composed of atomic sites', &
                                                             (st1num(i), i = 1, listdim)
              Write(message2,'(3x,2(a,i3),a)')              'of type-of-molecule ', mol(1), ' (set in  FF ', field, ')'
              messages(2) = trim(message1)//trim(message2)
              Write(messages(3),'(1x,a,i2,a)')              'could either not find its equivalent in FF ', field+1, ' or the &
                                                            &specification ordering for units is not consistent between these FFs'

            End If

            If(option /= 3)Then
              Write(messages(3),'(1x,a,i2,a)')              'could not find its equivalent in FF ', field+1, '.'
            End If

            siterror=.True.

          Else If(siteparam) Then
            If(option == 0)Then
              Write(messages(2),'(1x,2a,i4,a,2(i2,a))')     trim(string),' unit specification for atomic site ', st1num(1), &
                                                            ' of type-of-molecule ', mol(1), ' (set in FF ', field,') **DIFFERS**'
              Write(messages(3),'(1x,a,i4,a,2(i2,a))')      'from the specification for the corresponding atomic site', &
                                                            st2num(1), ' of type-of-molecule ', mol(2), ' (set in FF ', field+1,')'

            ElseIf(option == 1)Then
              Write(messages(2),'(1x,4a,2(i2,a))')          trim(string),' unit specification for atomic site ', trim(st1lab(1)), &
                                                            ' of type-of-molecule ', mol(1),        &
                                                            ' (set in FF ', field,') **DIFFERS** from'
              Write(messages(3),'(1x,3a,2(i2,a))')          'the specification for the corresponding atomic site ',&
                                                            trim(st2lab(1)),' of type-of-molecule ', mol(2),&
                                                            ' (set in FF ', field+1,')'
            Else If(option == 2)Then
              Write(messages(2),'(1x,2a,2(i6,a),2(i2,a))')  trim(string),' unit specification between atomic sites ', st1num(1), &
                                                            ' and ' , st1num(2) ,' of type-of-molecule ', mol(1),        &
                                                            ' (set in FF ', field,') **DIFFERS** from'
              Write(messages(3),'(1x,3a,2(i6,a),2(i2,a))')  'the ', trim(string), ' specification between atomic sites  ', &
                                                            st2num(1), ' and ', st2num(2),           &
                                                            ' of type-of-molecule ', mol(2), ' (set in FF ', field+1,')'
            Else If(option == 3)Then
              Write(messages(2),*)                          trim(string),' unit composed of atomic sites', &
                                                            (st1num(i), i = 1, listdim)
              Write(messages(3),'(2(a,i2),a,2(a,i2))')      ' of type-of-molecule', mol(1), ' (set in  FF ', field, &
                                                            ') **DIFFERS** from the what it should be the corresponding unit ',  &
                                                            'defined in type-of-molecule ', mol(2), ' of FF ', field+1

            End If
            siterror=.True.
          End If
          If(siterror)Then
            Call info(messages,4,.True.)
            Call error(0)
          End if
        End Do
      Else
        Call gsend(comm, siteevb  , 0, WriteConf_tag)
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


  Subroutine evb_check_vdw(evb, flow, sites, vdws)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the definition of vdW
    ! interactions between FIELD files. Check is carried out for pairs of vdW
    ! interactions that ONLY include atoms belonging to the non-reactive part of the system.
    ! In fact, the vdW interactions between atoms of the non-reactive site should remain the
    ! independently of the chemical state.
    !
    ! To check consistency between all FFs, the process starts by comparing FF1 and FF2,
    ! continues with FF2 with FF3, FF3 with FF4,....., and finishes FF(NFF) with FF1.
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti January 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(evb_type),         Intent(InOut) :: evb
    Type(flow_type),        Intent(In   ) :: flow
    Type(site_type),        Intent(In   ) :: sites(:)
    Type(vdw_type),         Intent(InOut) :: vdws(:)

    Integer               :: m, m2, k, i, j, l
    Integer               :: icvdw0, icvdw
    Integer               :: n_vdw_max

    Logical               :: lnovdw, vdwerror
    Logical               :: loopi, loopj

    Logical, Allocatable  :: vdw_noEVB(:,:)

    Character(Len = 256)  :: labelint

    Character(Len = 8)    :: stlab(2)

    Real(Kind = wp)       :: param(2)

    Character(Len = 256)  :: messages(3)

    ! Initialize values
    n_vdw_max = 0

    ! Initialise counters
    icvdw0=0

    ! Initialise flags
    lnovdw=.False.
    vdwerror=.False.

    Call info(' ',.True.)
    Call info(' vdW interactions for the non-reactive part of the system...', .True.)

    ! Check to find if there are vdW interactions specified in each FIELD file.
    ! In case there are no vdw Interactions for any of the FIELD files, counter icvdw
    ! is zero and checking is not done
    Do m=1,flow%NUM_FF
      If(vdws(m)%n_vdw == 0)Then
       icvdw0=icvdw0+1
      Else
        ! Identify the largest amount of vdw interactions between all FIELD files
        ! Due to the formation of new bonds, the amount of vdW interctions will be
        ! generally different for the involved FFs.
        If(vdws(m)%n_vdw > n_vdw_max)Then
          n_vdw_max=vdws(m)%n_vdw
        End If
       End If
    End Do

    If(icvdw0 == flow%NUM_FF)Then
     lnovdw=.True.
    Else
     ! Allocate logical array to indicate which vdW interactions only include pair
     ! of atoms outside the reactive region (see below)
     Allocate(vdw_noEVB(flow%NUM_FF,n_vdw_max))
     vdw_noEVB=.False.
    End If

    ! Perform the checking only if there are vdW interactions
    If(.Not. lnovdw)Then

      ! Identify those vdW interaction for which involved atoms are NOT part of the EVB region
      Do m = 1, flow%NUM_FF
        Do k = 1, vdws(m)%n_vdw
          icvdw=0
          Do j =1, 2
            i= evb%num_site(m)+1
            loopi=.True.
            Do While (i <= sites(m)%max_site .And. loopi)
              If(vdws(m)%labpair(j,k) == sites(m)%site_name(i))Then
                loopi=.False.
                icvdw=icvdw+1
              End If
            i=i+1
            End Do
          End Do
          ! If the two sites of vdws(m)%labpair(:,k) are part of the non-EVB region, set the element of
          ! vdw_noEVB to .true.
          If(icvdw==2)Then
            vdw_noEVB(m,k) = .True.
          End If
        End Do
      End Do


      ! Check consistency of the vdW interactions between FFs only for interacting ions that are not part of the EVB site
      ! Loop over the FFs
      Do m = 1, flow%NUM_FF
        ! Here FF(N_FF) and FF1 are also compared
        If(m==flow%NUM_FF)Then
          m2=1
        Else
          m2=m+1
        End If

        loopi=.True.
        i=1
        ! Loop over vdW interactions for FF(m)
        Do While (i <= vdws(m)%n_vdw .And. loopi)
          If(vdw_noEVB(m,i))Then
            loopj=.True.
            j=1
            ! Loop over vdW interactions for FF(m2)
            Do While (j <= vdws(m2)%n_vdw  .And. loopj)
              If(vdw_noEVB(m2,j))Then
                If((vdws(m)%labpair(1,i) == vdws(m2)%labpair(1,j)) .Or. (vdws(m)%labpair(1,i) == vdws(m2)%labpair(2,j)) )Then
                  If((vdws(m)%labpair(2,i) == vdws(m2)%labpair(1,j)) .Or. (vdws(m)%labpair(2,i) == vdws(m2)%labpair(2,j)))Then
                    loopj=.False.
                    ! If the vdW pair defined in FF(m) was found in FF(m2), check interaction key and parameters are consistent
                    Do l=1,2
                       stlab(l) =  vdws(m)%labpair(l,i)
                    End Do
                    If(vdws(m)%ltp(i) /= vdws(m2)%ltp(j))Then
                      vdwerror=.True.
                      loopi = .False.
                      labelint= 'key'
                    Else
                      Do k=1, vdws(m)%max_param
                         param(1)= vdws(m)%param(k,i)
                         param(2)=vdws(m2)%param(k,j)
                         If(abs(param(1)-param(2)) .gt. epsilon(param(1)))Then
                           vdwerror=.True.
                           loopi=.False.
                           labelint='parameters'
                         End If
                      End Do
                    End If
                  End If
                End If
              End If
              j=j+1
            End do

            ! If the vdW pair defined in FF(m) was found in FF(m2)....
            If(loopj)Then
              vdwerror=.True.
              loopi=.False.
              labelint= 'missed'
              Do l=1,2
                stlab(l) = vdws(m)%labpair(l,i)
              End Do
            End If

          End If
          i=i+1
        End Do

        ! Print message and abort if errors were found. Since all processors read all vdW interactions
        ! from all FFs, communication is NOT needed here.
        If(vdwerror)Then
          Write(messages(1), '(1x,5a,i2)') 'error - vdW interaction between sites ', trim(stlab(1)), ' and ', trim(stlab(2)),&
                                           ' in FF ', m
          If(labelint(1:6)=='missed')Then
            Write(messages(2), '(1x,a,i2)')  'has not been found in FF', m2
          Else
            Write(messages(2), '(1x,4a,i2)') 'has a different setting for the ', trim(labelint), ' with respect to the same ', &
                                             'interaction in FF ', m2
          End If
          Write(messages(3),'(1x,a)')     'ACTION: chek consistency of settings in vdw interactions for sites that are not &
                                          &part of the reactive region'
          Call info(messages, 3, .True.)
          Call error(0)
        End If

      End Do

      If (Allocated(vdw_noEVB)) Then
        Deallocate(vdw_noEVB)
      End If

    End If

  End Subroutine evb_check_vdw


  Subroutine evb_check_intermolecular(evb, flow, sites, tersoff, met, threebody, fourbody)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to check consistency in the definition of intermolecular
    ! interactions between FIELD files. Check is carried out for:
    !
    ! - Tersoff potentials
    ! - metallic potentials
    ! - three-body potentials (TBPs)
    ! - four-body potentials (FBPs)
    !
    ! Since EVB is designed to describe reactions (bond breaking and formation)
    ! the above interactions are only meaningful for the non-reactive part of the system. Thus,
    ! this subroutine first check that NONE of the EVB atoms interact via in ANY of these types of
    ! intermolecular interactions.
    !
    ! For example, EVB could be used to model the bond breaking/formation of a molecule
    ! supported by a metallic substrate. Metallic interactions will ONLY manifest when
    ! considering the interactions between atoms of the substrate. Molecule and surface might
    ! interact via vdW forces.
    !
    ! To check consistency between all FFs, the process starts by comparing FF1 and FF2,
    ! continues with FF2 with FF3, FF3 with FF4,....., and finishes FF(NFF-1) with FF(NFF).
    !
    ! For the case of four-body interactions, I (i.scivetti) have decided to prevent the EVB runs,
    ! mainly because four-body interactions were never tested to date (Feb 2020).
    !
    ! At first sight, and similarly to subroutine evb_check_constraints,
    ! one might consider there is a lot of code repetition.
    ! However, this is a consequence of having different formats related to
    ! the settings for the different intermolecular interactions.
    !
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


    Call info(' ',.True.)
    Call info(' inter-molecular interactions for the non-reactive part of the system...', .True.)

    ! Initialise flags
    sitemiss=.False.
    siteparam=.False.

    !! Metallic interactions
    !!
    !! Check consistency in the definition of metallic, only if they are defined in the FIELD files.
    !!
    !! DL_POLY will complain and abort the execution if a given metallic pair interactions:
    !!
    !! - is not found in all FIELD files
    !! - is described with different parameters in different FFs
    !!
    !! The sequence in the definition of tersoff interactions can be altered between FIELD files.
    !! Metallic interactions MUST NOT include atoms of the EVB region.

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
      ! Loop over FFs
      Do m = 1, flow%NUM_FF
        ! Loop over EVB atoms
        Do i=1,evb%num_at(m)
          ! Loop over metallic units
          Do j= 1, met(m)%n_potentials
            ! Loop over each member of the unit
            Do k = 1, 2
              If( trim(sites(m)%site_name(i)) == trim(met(m)%labunit(k,j)) )Then
                Write(messages(1),'(1x,3a,i2,a)') 'error - Site ', trim(sites(m)%site_name(i)),  &
                                                  ' belongs to the EVB site in FF ', m,            &
                                                  ' and MUST NOT interact via metallic forces'
                Call error(0, messages(1))
              End If
            End Do
          End Do
        End Do
      End do

      ! Check that metallic specification does not change beetween files
      ! Loop over FF
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        ! Loop over metallic interaction of FF(m)
        Do While (i <= met(m)%n_potentials .And. loopi)
          loopj=.True.
          j=1
          ! Loop over metallic interaction of FF(m+1)
          Do While (j <= met(m)%n_potentials .And. loopj)
            If((met(m)%labunit(1,i) == met(m+1)%labunit(1,j)) .Or. &
               (met(m)%labunit(1,i) == met(m+1)%labunit(2,j)) )Then
              Do k=1,2
                If(met(m)%labunit(1,i) == met(m+1)%labunit(k,j)) il = k
              End Do
              ! Check over the interacting pair sites
              Do k=1,2
                ! Check the label of each pair
                If( (met(m)%labunit(2,i) == met(m+1)%labunit(k,j)) .And. k /= il )Then
                  loopj=.False.
                  If(met(m)%labunit(3,i) /= met(m+1)%labunit(3,j) ) Then
                    loopi=.False.
                    siteparam =.True.
                    stlab(1) =  met(m)%labunit(1,i)
                    stlab(2) =  met(m)%labunit(2,i)
                  Else
                    ! Check interaction parameters
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

          ! If the metallic unit of FF(m) was not found in FF(m+1), identify the unit
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            stlab(1)=   met(m)%labunit(1,i)
            stlab(2)=   met(m)%labunit(2,i)
          End If

          i=i+1

        End Do
        ! Print erros and abort if there was an inconsistency found in the specification of metallic interactions
        labelint='Metal'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab)
      End Do
    End If

    !! Tersoff interactions
    !!
    !! Check consistency in the definition of Tersoff potentials, only if they are defined in the FIELD files.
    !!
    !! DL_POLY will complain and abort the execution if a given Tersoff atomic unit:
    !!
    !! - is not found in all FIELD files
    !! - is described with different parameters in different FFs
    !!
    !! For key "tersoff", checking is also conducted to compare the correction of the pairwise interaction
    !! between the Tersoff units.
    !!
    !! The sequence in the definition of Tersoff interactions can be altered between FIELD files.
    !! Atoms of the EVB region MUST NOT interact via Tersoff forces.

    ! Check that the number of Tersoff interactions are the same for all FIELD files
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
                Write(messages(1),'(1x,3a,i2,a)') 'PROBLEMS: Site ', trim(sites(m)%site_name(i)),  &
                                                  ' belongs to the EVB site in FF ', m,            &
                                                  ' and MUST NOT interact via Tersoff forces'
               Call error(0, messages(1))
              End If
          End Do
        End Do
      End do

      ! Check that Tersoff specification for each atomic unit does not change beetween files
      ! Loop over FFs
      Do m = 1, flow%NUM_FF-1
        loopi=.True.
        i=1
        ! Loop over Tersoff units for FF(m)
        Do While (i <= tersoff(m)%n_potential .And. loopi)
          loopj=.True.
          j=1
          ! Loop over Tersoff units for FF(m+1)
          Do While (j <= tersoff(m+1)%n_potential .And. loopj)
            If((tersoff(m)%labunit(1,i) == tersoff(m+1)%labunit(1,j)))Then
              If( tersoff(m)%labunit(2,i) == tersoff(m+1)%labunit(2,j) )Then
                loopj=.False.
                ! If Tersoff atomic unit and key were coincide, check input parameters
                Do l=1,tersoff(m)%max_param
                  param(1)=tersoff(m)%param(l,i)
                  param(2)=tersoff(m+1)%param(l,j)
                  ! If any of the values are inconsistent, exit loop i
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

          ! If Tersoff atomic unit of FF(m) was not found in FF(m+1), set sitemiss=.True. and exit loop.
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            stlab(1)= tersoff(m)%labunit(1,i)
          End If

          i=i+1

        End Do

        ! Print errors and abort if there was an inconsistency found in the specification of Tersoff units
        labelint='Tersoff-unit'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab)

      End Do

      ! Further checking in case the key of Tersoff potential is "tersoff"
      ! Check consistency in the pairwise interaction terms of the potential
      If(tersoff(1)%key_pot == 1)Then
        ! Check all the combinations are counted via definition of nprter
        nprter = (tersoff(1)%max_ter * (tersoff(1)%max_ter + 1)) / 2

        ! Loop over FFs
        Do m = 1, flow%NUM_FF-1
          loopi=.True.
          i=1
          ! Loop over pairwise interactions defined in FF(m)
          Do While (i <= nprter .And. loopi)
            loopj=.True.
            j=1
            ! Loop over pairwise interactions defined in FF(m+1)
            Do While (j <= nprter .And. loopj)
              If((tersoff(m)%labpair(1,i) == tersoff(m+1)%labpair(1,j)) .Or. &
                 (tersoff(m)%labpair(1,i) == tersoff(m+1)%labpair(2,j)) )Then
                Do k=1,2
                  If(tersoff(m)%labpair(1,i) == tersoff(m+1)%labpair(k,j)) il = k
                End Do
                ! If pair was found check the interaction parameters
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
        ! Print erros and abort if there was an inconsistency found in the specification of
        ! Tersoff pairwise interactions only for key "tersoff"
        labelint='Tersoff-pair'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab)
        End Do
      End If
    End If

    !! Three-body potentials (TBP)
    !!
    !! Check consistency in the definition of TBP, only if they are defined in the FIELD files.
    !!
    !! DL_POLY will complain and abort the execution if a given TBP unit:
    !!
    !! - is not found in all FIELD files
    !! - is described with different parameters in different FFs
    !!
    !! The sequence in the definition of Tersoff interactions can be altered between FIELD files.
    !! Atoms of the EVB region MUST NOT interact via TBP.

    ! Check that the number of defined units for TBP interactions are the same for all FIELD files
    Do m = 1, flow%NUM_FF-1
      If( threebody(m)%ntptbp /= threebody(m+1)%ntptbp )Then
        labelint='Total number of TBP pair interactions'
        Call evb_intermolecular_error(labelint, m)
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
                Write(messages(1),'(1x,3a,i2,a)') 'PROBLEMS: Site ', trim(sites(m)%site_name(i)),  &
                                                  ' belongs to the EVB site in FF ', m,            &
                                                  ' and MUST NOT interact via TBP forces'
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
        ! Loop over TBP for FF(m)
        Do While (i <= threebody(m)%ntptbp .And. loopi)
          loopj=.True.
          j=1
          ! Loop over TBP for FF(m+1)
          Do While (j <= threebody(m+1)%ntptbp .And. loopj)
            ic=0
            ! Loop over the three atoms to check if TBP units are the same
            Do k=1,3
             If(threebody(m)%labunit(k,i) == threebody(m+1)%labunit(k,j))Then
              ic=ic+1
             End If
            End Do
            If(ic == 3) loopj=.False.

            If(.not. loopj)Then
                    ! Check if TBP key is the same. If not, identify unit and leave loop
              If(threebody(m)%labunit(4,i) /= threebody(m+1)%labunit(4,j))Then
                loopi=.False.
                siteparam =.True.
                stlab(1) =  threebody(m)%labunit(1,i)
                stlab(2) =  threebody(m)%labunit(2,i)
              Else
                      ! Check if TBP parameters are the same. If not identify unit and leave loop
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

          ! If TBP unit has not been found in FF(m+1)
          If(loopj)Then
            loopi=.False.
            sitemiss=.True.
            Do k=1,3
             stlab(k) = threebody(m)%labunit(k,i)
            End Do
          End If

          i=i+1
        End Do

        ! Print errors and abort if there was an inconsistency found in the specification of Tersoff units.
        labelint='TBP'
        Call evb_intermolecular_error(labelint, m, sitemiss, siteparam, stlab)

      End Do

    End If

    !! Four-body potentials (TBP)
    !!
    !! As we have no test for this feature, it was decided to prevent the execution of DL_POLY is FBP are present
    !!
    Do m=1, flow%NUM_FF-1
      If(fourbody(m)%max_four_body /= 0)Then
        Write(messages(1),'(1x,a)') 'error - EVB  simulations with four-body potentials are currently not implemented'
        Call info(messages,1,.True.)
        Call error(0)
      End If
    End Do

  End Subroutine evb_check_intermolecular


  Subroutine evb_intermolecular_error(string, field, sitemiss, siteparam, stlab)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intermolecular) to print an
    ! error message and abort if either
    !
    ! - the number for a type of intermolecular interactions is different between FIELD files
    ! - the interaction unit is defined in one FIELD but not found in other (sitemiss=.True.)
    ! - the specification for the intermolecular interation unit changes between FIELD files (siteparam=.True.)
    !
    ! The format of the output error message depends on the intermolecular interactions
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
    Write(messages(1),'(1x,a)')        'error -'
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
            Write(messages(2),'(1x,4a,i2,a)')         trim(string),' specification/parameters for ',&
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
    ! interactions between FIELD files ONLY between atoms of the non-reactive part of the system.
    ! In fact, functional forms and parameters for intramolecular interactions SHOULD NOT change
    ! between atoms of the non-reacitve part of the system.
    !
    ! Check is carried out for bond, angle, dihedral and inversion potentials.
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

    ! Temporary array to count total number of interactions
    Integer, Allocatable :: numtot(:)

    ! Allocate array for the total number of intra-molecular interactions
    Allocate(numtot(flow%NUM_FF))

    Call info(' ',.True.)
    Call info(' intra-molecular interactions for the non-reactive part of the system...', .True.)

    ! Initialise the number for each type of interaction at EVB site
    evb%num_bond       = 0
    evb%num_angle      = 0
    evb%num_dihedral   = 0
    evb%num_inversion  = 0

    !! Bonds

    labelint= 'Bond'
    numtot=0

    ! For each FIELD, calculate the number of bonds for the EVB site and the total number of bonds for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+bond(m)%num(i)
        If(i <= evb%typemols(m))Then
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

    ! For the non-reactive part, check that specification for bonds does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 2, bond(m)%max_param, evb%num_bond(m), evb%num_bond(m+1), &
                               numtot(m), numtot(m+1), bond(m)%key, bond(m+1)%key, bond(m)%lst, bond(m+1)%lst, &
                               bond(m)%param, bond(m+1)%param, bond(m)%num, sites(m)%ntype_mol, BOND_TAB)
    End Do


    !! Angles
    labelint= 'Angle'
    numtot=0

    ! For each FIELD, calculate the number of angles for the EVB site and the total number of angles for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+angle(m)%num(i)
        If(i <= evb%typemols(m))Then
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

    ! For the non-reactive part, check that specification for angles does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 3, angle(m)%max_param, evb%num_angle(m), evb%num_angle(m+1), &
                                  numtot(m), numtot(m+1), angle(m)%key, angle(m+1)%key, angle(m)%lst, angle(m+1)%lst, &
                                  angle(m)%param, angle(m+1)%param, angle(m)%num, sites(m)%ntype_mol, ANGLE_TAB)
    End Do


    !! Dihedrals
    labelint= 'Dihedral'
    numtot=0

    ! For each FIELD, calculate the number of dihedrals for the EVB site and the total number of dihedrals for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+dihedral(m)%num(i)
        If(i <= evb%typemols(m))Then
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

    ! For the non-reactive part, check that specification for dihedrals does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 4, dihedral(m)%max_param, evb%num_dihedral(m), evb%num_dihedral(m+1), &
                                  numtot(m), numtot(m+1), dihedral(m)%key, dihedral(m+1)%key, dihedral(m)%lst, dihedral(m+1)%lst, &
                                  dihedral(m)%param, dihedral(m+1)%param, dihedral(m)%num, sites(m)%ntype_mol, DIHEDRAL_TAB)
    End Do


    !! Inversion
    labelint= 'Inversion'
    numtot=0

    ! For each FIELD, calculate the number of inversions for the EVB site and the total number of inversions for the whole system
    Do m= 1, flow%NUM_FF
      Do i= 1, sites(m)%ntype_mol
        numtot(m) = numtot(m)+inversion(m)%num(i)
        If(i <= evb%typemols(m))Then
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

    ! For the non-reactive part, check that specification for inversions does not change beetween FIELD files
    Do m = 1, flow%NUM_FF-1
      Call compare_intramolecular(labelint, m, 4, inversion(m)%max_param, evb%num_inversion(m), evb%num_inversion(m+1), &
                               numtot(m), numtot(m+1), inversion(m)%key, inversion(m+1)%key, inversion(m)%lst, inversion(m+1)%lst, &
                               inversion(m)%param, inversion(m+1)%param, inversion(m)%num, sites(m)%ntype_mol, INVERSION_TAB)
    End Do

    ! Deallocation of temporary array numtot
    Deallocate(numtot)

  End Subroutine evb_check_intramolecular


  Subroutine evb_intramolecular_number_error(labelint, field)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intramolecular) to print
    ! and error message and abort the execution if the number of intramolecular
    ! interactions (of a given type) for the non-reactive part of the system
    ! is not preserved for all FIELD files
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti December 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Character(Len = 256), Intent(In   ) :: labelint
    Integer(Kind = wi),   Intent(In   ) :: field


    Character(Len = 256) :: messages(5)

    Write(messages(1),'(1x,3a)')        'error - The total number of intramolecular ', trim(labelint),&
                                        ' interactions between atoms that'
    Write(messages(2),'(1x,a,i2,a,i2)') 'are NOT part of the EVB reactive fragment **DIFFERS** between FF ',&
                                         field, ' and FF ', field+1
    Write(messages(3),'(1x,3a)')        'ACTION: check settings for ', trim(labelint), ' interactions and make sure'
    Write(messages(4),'(1x,a)')         'they are equally defined in all FIELD files.'
    Write(messages(5),'(1x,a)')         'Settings for the non-reactive part of the system MUST be the same for all FIELD files'
    Call info(messages,5,.True.)
    Call error(0)

  End Subroutine evb_intramolecular_number_error


  Subroutine compare_intramolecular(labelint, field, listdim, maxparam, evb_num1, evb_num2, &
                                 numtot1, numtot2, key1, key2, list1, list2,             &
                                 param1, param2, numxmol, ntype_mol, keyexcl)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine (auxiliary to evb_check_intramolecular) to compare intramolecular
    ! interactions between FIELD files and print an error message (and abort) if either:
    !
    ! - the interaction unit is defined in one FIELD but not found in other
    ! - the specification for the intermolecular interation unit changes between FIELD files
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

        Write(messages(1), '(1x,a)' )        'error -'
        Write(messages(2), * )                trim(labelint),' interaction between atoms', (list1(k,i), k = 1, listdim)
        Write(messages(3),'(1x,a,2(i2,a))')  'is defined in FF ', field, &
                                             ' for type-of-molecule', mol, ' (see the corresponding FIELD file).'

        If(lparam)Then
          Write(messages(4),'(1x,a)')        'Parameters for this unit *DO NOT* agree with the values set'
        End If
        If(lkey)Then
          Write(messages(4),'(1x,a)')        'Key-type of interaction for this unit *DIFFERS* from the type specified'
        End If
        If(llist)Then
          Write(messages(4),'(1x,a)')        'This interaction unit CANNOT be found '
        End If

        Write(messages(5),'(1x,a,i2,3a)')    'in FF', field+1, '. ACTION: find unit and verify it is equally defined for ', &
                                             'all FIELD files, as settings for ', trim(labelint)
        Write(messages(6),'(1x,a)')          'interactions for the non-reactive part of the system must be the same for all FFs'
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
    ! dl_poly_4 subroutine (auxiliary to compare_intramolecular) to obtain the
    ! type-of-molecule given the index of an intramolecular interaction
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


  Subroutine evb_pes(evb, flow, config, stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that calls for the computation of:
    !
    ! - EVB coupling terms and their derivatives
    ! - EVB energy
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

    Integer(Kind = wi) :: m

    ! Shift configurational energies
    Do m=1,flow%NUM_FF
      evb%eneFF(m) = stat(m)%stpcfg+evb%eshift(m)
    End Do

    !Compute coupling terms and their gradients
    Call evb_couplings(flow, evb)

    !Compute EVB Energy
    Call evb_energy(evb, flow, stat)

    !Compute EVB forces
    Call evb_force(evb, flow, config)

    !Compute EVB stress
    Call evb_stress(evb, flow, config, stat)

  End Subroutine evb_pes


  Subroutine evb_couplings(flow, evb)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes:
    !
    ! - matrix of coupling terms (evb%coupl), used to compute the EVB matrix
    ! - matrix with the derivatives of coupling terms (evb%grad_coupl), used
    !   to EVB forces and stress tensor
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti August 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Type(flow_type),  Intent(In   ) :: flow
      Type(evb_type  ), Intent(InOut) :: evb

      Integer                         :: m,k
      Real(Kind = wp)                 :: ediff, earg
      Real(Kind = wp)                 :: A(evb%maxparam)

      ! Initialise coupling matrices to zero
       evb%coupl=0.0d0
       evb%grad_coupl=0.0d0

      Do m=1,flow%NUM_FF
        Do k=m+1,flow%NUM_FF
        ! Copy evb coupling parameters to array A, for the sake of clarity in coding the formulae
          A(:)  = evb%coupl_param(m,k,:)
        ! Define energy difference, which is the reaction corrdinate between fields m and k
          ediff = evb%eneFF(m)-evb%eneFF(k)
          If(evb%typcoupl(m,k)=='const')Then
            evb%coupl(m,k) = A(1)
            evb%grad_coupl(m,k) = 0.0_wp
          Else If (evb%typcoupl(m,k)=='gauss')Then
            earg=Abs((ediff-A(2))/A(3))**2
            If(earg < evb%elimit)Then     ! Here we use the sensible limit of 700 for the argument to compute exp()
              evb%coupl(m,k) = A(1)*exp(-((ediff-A(2))/A(3))**2)+A(4)
              evb%grad_coupl(m,k) = -2.0_wp * (A(1)/A(3)**2)*(ediff-A(2))*exp(-((ediff-A(2))/A(3))**2)
            Else
              Call evb_couplings_error(m,k,evb)
            End If
          End If
          ! Make both matrices symmetric
          evb%coupl(k,m)      = evb%coupl(m,k)
          evb%grad_coupl(k,m) = evb%grad_coupl(m,k)

        End Do
      End Do

  End Subroutine evb_couplings


  Subroutine evb_couplings_error(m,k,evb)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to print an error message and abort if any of the
    ! coupling terms exhibit numerical instabilities
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti August 2019
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer(Kind = wi), Intent(In   ) :: m, k
    Type(evb_type  ),   Intent(InOut) :: evb

    Character(Len = 256)   :: messages(4)

    Write(messages(1),'(1x,a)')          'error - Argument for the exponential term of the "gauss" coupling'
    Write(messages(2),'(1x,2(a,i2),a)')  'between FFs', m, ' and ', k, ' is beyond the range of computation'
    Write(messages(3),'(1x,2(f6.1,a))')   -evb%elimit, ' < argument < ', evb%elimit, '. Either the energy of one FF '
    Write(messages(4),'(1x,a)')          'has diverged or the coupling parameters set in SETEVB are not correct.'
    Call info(messages, 4, .True.)
    Call error(0)

  End Subroutine evb_couplings_error


  Subroutine evb_energy(evb, flow, stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes the EVB energy. Steps:
    !
    ! - buildig the EVB matrix. Shifted conformational energies for each FF
    !   are set in the diagonal elements; off diagonal terms are set equal
    !   to evb%coupl, previously computed via suroutine evb_couplings.
    ! - Diagonalization of the EVB matrix via the BLAS subroutine dsyevx
    ! - Check if diagonalization was successful
    ! - EVB eigenvectors are assigned to evb%psi
    ! - Lowest eigenvalue evb%eigval is the EVB energy
    !
    ! This procedure is repeated at each time step
    !
    ! #ifdef condition is to allow serial compilation with blas subroutine dsyevx
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


#ifdef EVB
    !Initialise matrix elements
    evb%ene_matrix = 0.0_wp

    !! Build EVB matrix
    ! Diagonal elements
    Do m=1,flow%NUM_FF
      evb%ene_matrix(m,m)= evb%eneFF(m)
    End Do

    !Off-diagonal terms
    Do m=1,flow%NUM_FF
      Do k=m+1,flow%NUM_FF
        evb%ene_matrix(m,k) = evb%coupl(m,k)
        evb%ene_matrix(k,m) = evb%coupl(k,m)
      End Do
    End Do

    !!Diagonalization
    call dsyevx( 'V', 'I', 'U', flow%NUM_FF, evb%ene_matrix, flow%NUM_FF, -1., 0., 1, flow%NUM_FF, &
                  1.d-30, mevb, evb%eigval, evb%psi, flow%NUM_FF, evb%work, 8*flow%NUM_FF,     &
                  evb%iwork, evb%ifail, evbinfo)

    !Check that diagonalization was successful
    If(evbinfo /= 0)Then
      Call evb_diag_error(evbinfo)
    End If

    !Assign eigenvalues to conformational energy
    Do m=1,flow%NUM_FF
      stat(m)%stpcfg=evb%eigval(1)
    End Do
#endif

  End Subroutine evb_energy


  Subroutine evb_diag_error(evbinfo)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to print an error and abort if the diagonalization
    ! was unsuccessful
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti Febraury 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer(Kind = wi), Intent(In   ) :: evbinfo

    Character(Len = 256)   :: messages(4)


    If(evbinfo < 0)Then
      Write(messages(1),'(1x,a,i3,a)')   'error - EVB diagonalization failed as argument ', -evbinfo, &
                                         ' has an illegal value.'
    Else If(evbinfo > 0)Then
      Write(messages(1),'(1x,a,i2,2a)')  'error - EVB diagonalization failed as eigenvector ', evbinfo, ' failed ',&
                                         'to converge'
    End If

    Write(messages(2),'(1x,a)')          'ACTION 1: check for incorrect input parameters for coupling terms in SETEVB'
    Write(messages(3),'(1x,a)')          'ACTION 2: check for incorrect input parameters in FIELD and CONFIG files. It'
    Write(messages(4),'(1x,a)')          'might be helful to run standard MD simulations for each FFs separately'

    Call info(messages,4,.True.)
    Call error(0)

  End Subroutine evb_diag_error


  Subroutine evb_force(evb, flow, config)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes for the EVB force over each ion of the
    ! system. Steps:
    !
    ! - build force matrix for each ion and coordinate, using matrix evb%grad_coupl
    !   computed previously by subroutine evb_couplings.
    ! - compute force evb%force using the EVB eigenvector.
    ! - copy evb%force to the force components of config(m)%parts(i)
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti March 2020
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
      End Do

    ! Once we have compute the the three force components, we copy evb-force to the configurational force of each FF
      Do m=1,flow%NUM_FF
        config(m)%parts(i)%fxx=evb%force(1,1)
        config(m)%parts(i)%fyy=evb%force(2,1)
        config(m)%parts(i)%fzz=evb%force(3,1)
      End Do

    End Do

   End Subroutine evb_force


  Subroutine evb_stress(evb, flow, config, stat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that computes for the EVB stress tensor and
    ! EVB configurational virial. Steps:
    !
    ! - build stress matrix evb%dE_dh using matrix evb%grad_coupl (computed previously
    !   by subroutine evb_couplings).
    ! - compute stress tensor evb%stress using the EVB eigenvector.
    ! - copy evb%stress to the stat(m)%stress
    ! - compute the total EVB virial
    !
    ! In principle, it is possible to have an EVB decomposition of the virial into separate
    ! contributions for various type of interactions (e.g. angles, bonds, dihedrals,
    ! coulombic, etc). However, i.scivetti has decided not to compute these contributions
    ! as their implementation would entail a significant amount of changes to the code, and
    ! it would be completely irrelevant to describing the dynamics. Instead, the total
    ! configurational virial is computed using the trace of the EVB stress tensor
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti March 2020
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

  End Subroutine evb_stress


  Subroutine evb_population(evb, flow, files, comm)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to compute the population of EVB states at each timestep
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
        Write(files(FILE_POPEVB)%unit_no,'(a,17x,a)') '# Time (ps)','Weights of each FFs in the total EVB state (FF 1, 2, etc)'
      Else
        Open(Newunit=files(FILE_POPEVB)%unit_no,File=files(FILE_POPEVB)%filename, Position='append')
      End If

    End If

    ! Only print after equilibration
    If(flow%step>flow%equil_steps+1)Then
       If (comm%idnode == 0)Then
         ! The specified writing format below allows to write a variable number of FFs
         Write(files(FILE_POPEVB)%unit_no,'(*(e16.6))') flow%time, (evb%psi(ff,1)**2, ff=1,flow%NUM_FF)
       End If
    End If

  End Subroutine evb_population


  Subroutine evb_setzero(flow, stat)
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

  End Subroutine evb_setzero

  Subroutine evb_merge_stochastic(flow, config, stat, rigid, thermo, cshell, cons, pmf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine that copies various variable type from the force-field 1
    ! to the rest of the force fields. This copy is needed only when using the following features:
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

  End Subroutine evb_merge_stochastic

  Subroutine print_evb_banner(ffnum)
    Integer,          Intent(In   ) :: ffnum

    Character(Len=*), Parameter :: fmt1 = '(a)', fmt2  = '(a,i2)'

    Character(Len=256) :: banner(7)

    Write (banner(1), fmt1)  '*******************************************************************'
    Write (banner(2), fmt1)  '*******************************************************************'
    Write (banner(3), fmt1)  ' Simulation with the Empirical Valence Bond (EVB) method        '
    Write (banner(4), fmt1)  '  '
    Write (banner(5), fmt2)  ' Number of Force-Fields (FF) to be coupled ', ffnum
    Write (banner(6), fmt1)  '*******************************************************************'
    Write (banner(7), fmt1)  '*******************************************************************'
    Call info(banner, 7, .true., level=-1)
  End Subroutine print_evb_banner

End module evb

