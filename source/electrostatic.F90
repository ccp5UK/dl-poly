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
  Use numerics, Only : interp_table, calc_erfc, calc_erfc_deriv
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

    type ( interp_table ), Public :: erfc, erfc_deriv

  Contains
    Procedure, Public :: init_erf_tables
    Procedure, Public :: erfcgen
    ! If an exact erfc is desired
    Procedure, NoPass, Public :: calc_erfc => calc_erfc
    Procedure, NoPass, Public :: calc_erfc_deriv => calc_erfc_deriv
  End Type electrostatic_type

contains
  Subroutine init_erf_tables(electro, nsamples)
    Class( electrostatic_type ), Intent( InOut ) :: electro
    Integer,                     Intent( In    ) :: nsamples

    electro%erfc%nsamples = nsamples
    electro%erfc_deriv%nsamples = nsamples

  End Subroutine init_erf_tables

  Subroutine erfcgen(electro,rcut,alpha) !,erc,fer)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 routine for generating interpolation tables for erfc and its
    !! derivative - for use with Ewald sum
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester december 1994
    !! amended   - i.t.todorov february 2016
    !! amended   - j.s.wilkins october  2018
    !!-----------------------------------------------------------------------
    Implicit None

    Class( electrostatic_type ), Intent ( InOut ) :: electro
    ! Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
    ! Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
    ! Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
    ! Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
    ! Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
    ! Real( Kind = wp ), Parameter :: pp =  0.3275911_wp

    Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha
    ! Type( interp_table ),                     Intent( InOut ) :: erc, fer

    Integer           :: i
    Real( Kind = wp ) :: rrr, rsq
    ! Real( Kind = wp ) :: exp1,rrr,rsq,tt
    Integer           :: fail

    if (electro%erfc%initialised .and. electro%erfc_deriv%initialised) return

    call electro%erfc%init(rcut, erfc_ar_over_r)
    call electro%erfc_deriv%init(rcut, erfc_ar_over_r_deriv)

    ! allocate(electro%erfc%table(0:electro%erfc%nsamples), stat=fail)
    ! if (fail > 0) call error_alloc('electro%erfc%table','erfcgen')
    ! allocate(electro%erfc_deriv%table(0:electro%erfc_deriv%nsamples), stat=fail)
    ! if (fail > 0) call error_alloc('electro%erfc_deriv%table','erfcgen')



    ! look-up tables for real space part of ewald sum

    ! electro%erfc%spacing = rcut/Real(electro%erfc%nsamples-4,wp)
    ! electro%erfc_deriv%spacing = electro%erfc%spacing

    ! electro%erfc%recip_spacing = 1.0_wp / electro%erfc%spacing
    ! electro%erfc_deriv%recip_spacing = 1.0_wp / electro%erfc_deriv%spacing

    ! Do i=1,electro%erfc%nsamples
    !   rrr = Real(i,wp)*electro%erfc%spacing
    !   rsq = rrr*rrr

    !   electro%erfc%table(i) =
    !   electro%erfc_deriv%table(i) = (electro%erfc%table(i) + alpha*calc_erfc_deriv(alpha*rrr))/rsq

    !   ! tt = 1.0_wp/(1.0_wp + pp*alpha*rrr)

    !   ! exp1 = Exp(-(alpha*rrr)**2)

    !   ! electro%erfc%table(i) = tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rrr
    !   ! electro%erfc_deriv%table(i) = (electro%erfc%table(i) + 2.0_wp*(alpha/sqrpi)*exp1)/rsq
    ! End Do

    ! ! extrapolation for grid point 0 at distances close to 0

    ! electro%erfc%table(0) = Huge(1.0_wp)
    ! electro%erfc_deriv%table(0) = Huge(1.0_wp+2.0_wp*(alpha/sqrpi))

    ! electro%erfc%end_sample = electro%erfc%table(electro%erfc%nsamples-4)
    ! electro%erfc_deriv%end_sample = electro%erfc_deriv%table(electro%erfc_deriv%nsamples-4)

    ! electro%erfc%initialised = .true.
    ! electro%erfc_deriv%initialised = .true.

  contains
    Function erfc_ar_over_r(rrr)
      Real( Kind = wp ) :: rrr
      Real( Kind = wp ) :: erfc_ar_over_r

      erfc_ar_over_r = calc_erfc(alpha*rrr) / rrr
    end Function erfc_ar_over_r

    Function erfc_ar_over_r_deriv(rrr)
      Real( Kind = wp ) :: rrr
      Real( Kind = wp ) :: erfc_ar_over_r_deriv

      erfc_ar_over_r_deriv = (erfc_ar_over_r(rrr) + alpha*calc_erfc_deriv(alpha*rrr))/rsq
    end Function erfc_ar_over_r_deriv


  End Subroutine erfcgen




End Module electrostatic
