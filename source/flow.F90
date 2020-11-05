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
  !> Simulation type keys: MD_STD
  Integer(Kind=wi), Parameter, Public :: MD_STD = 1
  !> Simulation type keys: EVB
  Integer(Kind=wi), Parameter, Public :: EmpVB = 2
  !> Simulation type keys: FFS
  Integer(Kind=wi), Parameter, Public :: FFS = 3
  !> Simulation type keys: DFTB+
  Integer(Kind=wi), Parameter, Public :: DFTB = 4

  Integer(Kind=wi), Parameter, Public :: MAX_FF = 3

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
    Logical, Public          :: strict = .true.
    !> Topology printing switch
    Logical, Public          :: print_topology = .true.
    !> Force capping switch
    Logical, Public          :: force_cap = .false.
    !> Equilibration state flag
    Logical, Public          :: equilibration = .true.
    !> Full simulation (_i.e._ not replay) flag
    Logical, Public          :: simulation = .true.
    !> Replay calculate forces
    Logical, Public          :: replay_recalculate_forces = .false.
    !> Book keeping flag
    Logical, Public          :: book = .true.
    !> Excluded interactions flag
    Logical, Public :: exclusions

    !> Write per-particle information
    Logical, Public :: write_per_particle

    !> Calculate heat-flux
    Logical, Public :: heat_flux

    !> Restart key
    Integer(Kind=wi), Public :: restart_key = RESTART_KEY_CLEAN
    !> Current simulation step
    Integer(kind=wi), Public :: step =0
    !> Current simulation time (step * timestep)
    Real(Kind=wp), Public    :: time = 0.0_wp
    !> Starting time, non-zero if job is restarted
    Real(Kind=wp), Public    :: start_time = 0.0_wp
    !> Number of production steps
    Integer(Kind=wi), Public :: run_steps = 0
    !> Number of equilibration steps
    Integer(Kind=wi), Public :: equil_steps = 0
    !> Data printing interval (in steps)
    Integer(Kind=wi), Public :: freq_output = 100
    !Analyse
    !> Analyse bonds
    Logical, Public :: analyse_bond = .false.
    !> Analyse angles
    Logical, Public :: analyse_ang = .false.
    !> Analyse dihedrals
    Logical, Public :: analyse_dih = .false.
    !> Analyse inversions
    Logical, Public :: analyse_inv = .false.
    !> Bond distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_bond = 0
    !> Angle distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_angle = 0
    !> Dihedral distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_dihedral = 0
    !> Inversion distribution calculation period (in steps)
    Integer(Kind=wi), Public :: freq_inversion = 0
    !> Restart files creation period (in steps)
    Integer(Kind=wi), Public :: freq_restart = 1000
    !> Reset padding flag
    Logical, Public :: reset_padding
    !> Type of Simulation we are performing
    Integer, Public          :: simulation_method = MD_STD
    !> MD step that DL_POLY starts at
    Integer, Public          :: initial_md_step = 0
    !> Define number of force-fields to be coupled
    Integer(Kind = wi), Public :: num_ff = 1
    Logical, Public          :: l_vdw = .False.

  Contains
    Procedure, Public :: new_page => flow_type_new_page
    Procedure, Public :: line_printed => flow_type_line_printed
  End Type flow_type

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

End Module flow_control
