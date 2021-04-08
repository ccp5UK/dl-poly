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
  Use numerics, Only : interp_table, calc_erfc_n, calc_erfc_deriv_n
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

    !> Electrostatic potential key
    Integer( Kind = wi ), Public :: key = ELECTROSTATIC_NULL

    !> No electrostatics switch
    Logical, Public :: no_elec = .false.

    Logical,           Public :: initialised = .false.
    Logical,           Public :: multipolar = .false.
    Type( mpole_type), Public :: mpoles
    Integer,           Public :: num_mpoles = 0

    !> How many mpole derivatives
    Integer, Dimension(0:0), Public :: nmpole_derivs = [1]
    !> My mpole derivatives
    Integer, Dimension(3,1,0:0), Public :: mpole_derivs = reshape([0,0,0],[3,1,1])

    !> Damped or not?
    Logical,           Public :: damp = .false.
    !> Damping distance
    Real( Kind = wp ), Public :: damping = 0.0_wp

    !> Relative dielectric constant
    Real( Kind = wp ), Public :: eps = 1.0_wp

    Logical, Public :: lecx = .false.

    Integer, Public :: nstfce = 1

    Real( Kind = wp ), Public :: force_shift = 0.0_wp, energy_shift = 0.0_wp

    Real( Kind = wp ), Dimension(0:2), Public :: reaction_field = 0.0_wp

    type ( interp_table ), Public :: erfc, erfc_deriv

  Contains
    Procedure, Public :: init_erf_tables
    Procedure, Public :: erfcgen
    ! If an exact erfc is desired
    Procedure, NoPass, Public :: calc_erfc => calc_erfc_n
    Procedure, NoPass, Public :: calc_erfc_deriv => calc_erfc_deriv_n
  End Type electrostatic_type

contains
  Subroutine init_erf_tables(electro, nsamples)
    Class( electrostatic_type ), Intent( InOut ) :: electro
    Integer,                     Intent( In    ) :: nsamples

    electro%erfc%nsamples = nsamples
    electro%erfc_deriv%nsamples = nsamples

  End Subroutine init_erf_tables

  Subroutine erfcgen(electro,rcut,alpha)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 routine for generating interpolation tables for erfc/r and its
    !! derivative - for use with Ewald sum
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester december 1994
    !! amended   - i.t.todorov february 2016
    !! amended   - j.s.wilkins september 2019
    !!-----------------------------------------------------------------------
    Implicit None

    Class( electrostatic_type ), Intent ( InOut ) :: electro

    Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha

    if (electro%erfc%initialised .and. electro%erfc_deriv%initialised) return

    call electro%erfc%init(rcut, erfc_ar_over_r)
    call electro%erfc_deriv%init(rcut, erfc_ar_over_r_deriv)

  contains
    Function erfc_ar_over_r(rrr)
      Real( Kind = wp ) :: rrr
      Real( Kind = wp ) :: erfc_ar_over_r

      erfc_ar_over_r = calc_erfc_n(alpha*rrr) / rrr
    end Function erfc_ar_over_r

    Function erfc_ar_over_r_deriv(rrr)
      Real( Kind = wp ) :: rrr, rsq
      Real( Kind = wp ) :: erfc_ar_over_r_deriv

      rsq = rrr ** 2
      erfc_ar_over_r_deriv = (erfc_ar_over_r(rrr) + alpha*calc_erfc_deriv_n(alpha*rrr))/rsq
    end Function erfc_ar_over_r_deriv


  End Subroutine erfcgen

End Module electrostatic
