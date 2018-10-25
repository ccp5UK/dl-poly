!> Module for electrostatic data and routines, common between Ewald and 
!> non-Ewald methods
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J.Madge July 2018
!> Amended - J.S. Wilkins October 2018
Module electrostatic
  Use kinds, Only : wi,wp
  Use mpole, Only : mpole_type
  Implicit None

  Private

  ! Electrostatic potential keys
  !> No electrostatics
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_NULL = 0
  !> Ewald Sum
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_EWALD = 1
  !> Distance dependent dielectric potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_DDDP = 2
  !> Direct real-space Coulomb potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_COULOMB = 3
  !> Force-shifted and damped Coulomb potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_COULOMB_FORCE_SHIFT = 4
  !> Reaction field and damped Coulomb potential
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_COULOMB_REACTION_FIELD = 5
  !> Direct space Poisson solver
  Integer( Kind = wi ), Parameter, Public :: ELECTROSTATIC_POISSON = 6

  !> Type containing electrostatic potential data
  Type, Public :: electrostatic_type
    Private

    !> No electrostatics switch
    Logical, Public :: no_elec

    !> Electrostatic potential key
    Integer( Kind = wi ), Public :: key = ELECTROSTATIC_NULL

    Logical,           Public :: initialised
    Logical,           Public :: multipolar = .false.
    Type( mpole_type), Public :: mpoles
    Integer,           Public :: num_mpoles = 0

    !> How many mpole derivatives
    Integer, Dimension(0:0), Public :: nmpole_derivs = [1]
    !> My mpole derivatives
    Integer, Dimension(3,1,0:0), Public :: mpole_derivs = reshape([0,0,0],[3,1,1])

    !> Damped or not?
    Logical,           Public :: damp
    !> Damping distance
    Real( Kind = wp ), Public :: damping
    
    !> Relative dielectric constant
    Real( Kind = wp ), Public :: eps

    Logical, Public :: lecx
    
    Integer, Public :: nstfce
    
    Real( Kind = wp ), Public :: force_shift = 0.0_wp, energy_shift = 0.0_wp
    
    Real( Kind = wp ), Dimension(0:2), Public :: reaction_field

  End Type electrostatic_type

End Module electrostatic
