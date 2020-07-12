Module flow_control
  Use kinds, Only: wi,&
                   wp

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
  !> Simulation type keys: MD
  Integer(Kind=wi), Parameter, Public :: MD = 1
  !> Simulation type keys: EVB
  Integer(Kind=wi), Parameter, Public :: EVB = 2
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
