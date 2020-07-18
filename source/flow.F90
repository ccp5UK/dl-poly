Module flow_control

!> Module to define flow settings
!>
!> Copyright - Daresbury Laboratory
!
!> Author - a.m.elena March 2018
!> refactoring:
!>           - a.m.elena march-october 2018
!>           - j.madge march-october 2018
!>           - a.b.g.chalk march-october 2018
!>           - i.scivetti march-october 2018
!> contrib   - i.scivetti Aug 2019 - Add read_simtype

  Use kinds,           Only : wi,&
                              wp
  Use parse,           Only : get_line,&
                              get_word,&
                              lower_case,&
                              word_2_real

  Use errors_warnings, Only : error,&
                              info

  Use comms,           Only : comms_type,gcheck

  Implicit None

  Private

  ! Simulation restart keys
  !> Clean start (New simulation?)
  Integer(Kind=wi), Parameter, Public :: RESTART_KEY_CLEAN = 0
  !> Continue an old simulation
  Integer(Kind=wi), Parameter, Public :: RESTART_KEY_OLD = 1
  !> Scaled restart
  Integer(Kind=wi), Parameter, Public :: RESTART_KEY_SCALE = 2
  !> Unscaled restart
  Integer(Kind=wi), Parameter, Public :: RESTART_KEY_NOSCALE = 3
  !> Simulation type keys: STD
  Integer(Kind=wi), Parameter, Public :: MD = 1
  !> Simulation type keys: EVB
  Integer(Kind=wi), Parameter, Public :: EmpVB = 2
  !> Simulation type keys: FFS
  Integer(Kind=wi), Parameter, Public :: FFS = 3
  !> Simulation type keys: DFTB+
  Integer(Kind=wi), Parameter, Public :: DFTB = 4

  !> Type containing program flow data
  Type, Public :: flow_type
    Private

    !> Check if is first time we call build_book_intra
    Logical, Public          :: newjob_build_book = .true.
    Logical, Public          :: oldjob_shared_units = .false.

    ! STDOUT printing control
    !> Number of print events before starting a new 'page'
    Integer(Kind=wi)         :: npage = 8
    !> Current number of print events
    Integer(Kind=wi), Public :: lines = 0
    !> Check if first call of md_vv or calculate_forces
    Logical, Public          :: newjob = .true.
    !> Strict mode
    Logical, Public          :: strict
    !> Topology printing switch
    Logical, Public          :: print_topology
    !> Force capping switch
    Logical, Public          :: force_cap
    !> Equilibration state flag
    Logical, Public          :: equilibration
    !> Full simulation (_i.e._ not replay) flag
    Logical, Public          :: simulation
    !> Book keeping flag
    Logical, Public          :: book
    !> Excluded interactions flag
    Logical, Public :: exclusions

    !> Write per-particle information
    Logical, Public :: write_per_particle

    !> Calculate heat-flux
    Logical, Public :: heat_flux

    !> Restart key
    Integer(Kind=wi), Public :: restart_key
    !> Current simulation step
    Integer(kind=wi), Public :: step
    !> Current simulation time (step * timestep)
    Real(Kind=wp), Public    :: time
    !> Starting time, non-zero if job is restarted
    Real(Kind=wp), Public    :: start_time
    !> Number of production steps
    Integer(Kind=wi), Public :: run_steps
    !> Number of equilibration steps
    Integer(Kind=wi), Public :: equil_steps
    !> Data printing interval (in steps)
    Integer(Kind=wi), Public :: freq_output
    !> Bond distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_bond
    !> Angle distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_angle
    !> Dihedral distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_dihedral
    !> Inversion distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_inversion
    !> Restart files creation period (in steps
    Integer(Kind=wi), Public :: freq_restart
    !> Reset padding flag
    Logical, Public :: reset_padding
    !> Type of Simulation we perform
    Integer, Public          :: simulation_method = MD
    !> MD step that DL_POLY starts at
    Integer, Public          :: initial_md_step = 0
    !> Define number of force-fields to be coupled
    Integer(Kind = wi), Public :: NUM_FF
    !> Flags to kill EVB if reading "evb" settings failed
    Logical,            Public :: evbfail = .False.

    Logical, Public          :: l_vdw = .False.

  Contains
    Procedure, Public :: new_page => flow_type_new_page
    Procedure, Public :: line_printed => flow_type_line_printed
  End Type flow_type

  Public :: read_simtype

Contains

  Pure Function flow_type_new_page(T) Result(new_page)
    Class(flow_type), Intent(In   ) :: T
    Logical                         :: new_page

    new_page = Mod(T%lines, T%npage) == 0
  End Function flow_type_new_page

  Subroutine flow_type_line_printed(T)
    Class(flow_type), Intent(InOut) :: T

    T%lines = T%lines + 1
  End Subroutine flow_type_line_printed


  Subroutine read_simtype(flow,comm) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to read option evb from CONTROL file. If present,
    ! the number of fields to be coupled is assigned to flow%NUM_FF. If there 
    ! is an error in the specification for this option evb or the CONTROL file is 
    ! not found, variables are temporarily assigned and DL_POLY will print an 
    ! error message and abort later once OUTPUT file has been opened
    !
    ! copyright - daresbury laboratory
    ! author    - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(flow_type) , Intent(  Out) :: flow
    Type(comms_type), Intent(InOut) :: comm
  
    Logical                :: carry,safe, stdtype
    Character( Len = 200 ) :: record
    Character( Len = 40  ) :: word
  
    Integer       :: unit_no
  
    ! Set safe flag
      safe    =.true.
      stdtype =.true.
  
    ! Open CONTROL file
    If (comm%idnode == 0) Inquire(File='CONTROL', Exist=safe)
    Call gcheck(comm,safe,"enforce")
    If (.not.safe) Then
      ! If CONTROL file is not found, set the following variables and return.
      ! DL_POLY will later abort by printing an error message
      ! At this stage we cannot stop the execution because OUTPUT has not been opened yet
      flow%simulation_method=MD
      flow%NUM_FF = 1 
      Return
    Else
      ! If CONTROL file is found, proceed
      If (comm%idnode == 0) Then
        Open(Newunit=unit_no, File='CONTROL',Status='old')
      End If
    End If
  
    Call get_line(safe,unit_no,record,comm)  
  
    If (safe) Then
      carry = .true.
      Do While (carry)
  
        Call get_line(safe,unit_no,record,comm)
        If (.not.safe) Exit
        Call lower_case(record)
        Call get_word(record,word)
  
        ! read EVB option. If for any reason there was a typo with the option "evb", DL_POLY will complain later
        If (word(1:3) == 'evb') Then
          flow%simulation_method=EmpVB
          stdtype =.false.
          Call get_word(record,word)
          flow%NUM_FF = Nint(word_2_real(word,0.0_wp))
          If(flow%NUM_FF <= 1)Then
          ! If there is no value assigned or a wrong input value has been set, activate evbfail
            flow%NUM_FF  = 1
            flow%evbfail = .True.         
          End If        
        Else If (word(1:6) == 'finish') Then
          carry=.false.
        End If
      End Do
  
      If(stdtype)Then
      ! Set standard option
        flow%simulation_method=MD
        flow%NUM_FF = 1 
      End If         
  
    End If
  
    ! Close CONTROL file
    If (comm%idnode == 0) Close(unit_no)    

  End Subroutine read_simtype
          

End Module flow_control
