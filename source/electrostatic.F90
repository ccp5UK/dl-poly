!> Module for electrostatic data and routines, common between Ewald and
!> non-Ewald methods
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J.Madge July 2018
!> Amended - J.S. Wilkins October 2018
!> Contrib - B.T Speake July 2024 - charge smearing
Module electrostatic
  Use charge_smearing, Only : linear_smearing_pot, linear_smearing_force, &
                              slater_exp_smearing_pot, slater_exp_smearing_force, &
                              BETA_ORIGINAL, slater_apprx_smearing_pot, &
                              slater_apprx_smearing_force
  Use kinds, Only : wi,wp
  Use mpole, Only : mpole_type
  Use numerics, Only : interp_table, calc_erfc_n, calc_erfc_deriv_n
  Implicit None

  Private

  ! Electrostatic potential keys
  !> No electrostatics
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_NULL = 0
  !> Ewald Sum
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_SPME = 1
  !> Distance dependent dielectric potential
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_DDDP = 2
  !> Direct real-space Coulomb potential
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_COULOMB = 3
  !> Force-shifted and damped Coulomb potential
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_COULOMB_FORCE_SHIFT = 4
  !> Reaction field and damped Coulomb potential
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_COULOMB_REACTION_FIELD = 5
  !> Direct space Poisson solver
  Integer(Kind=wi), Parameter, Public :: ELECTROSTATIC_POISSON = 6

  ! Charge smearing keys
  Integer(Kind=wi), Parameter, Public :: SMEARING_NULL = 0
  Integer(Kind=wi), Parameter, Public :: SMEARING_LINEAR = 1
  Integer(Kind=wi), Parameter, Public :: SMEARING_SLATER_TRUNCATED = 2
  Integer(Kind=wi), Parameter, Public :: SMEARING_SLATER_EXP = 3
  Integer(Kind=wi), Parameter, Public :: SMEARING_GAUSSIAN = 4
  Integer(Kind=wi), Parameter, Public :: SMEARING_GAUSSIAN_EQUAL = 5

  !> Type containing electrostatic potential data
  Type, Public :: electrostatic_type
    Private

    !> Electrostatic potential key
    Integer(Kind=wi), Public              :: key = ELECTROSTATIC_NULL
    !> Charge smearing key
    Integer(Kind=wi), Public              :: smear = SMEARING_NULL
    !> Charge smearing length
    Real(Kind=wp), Public                 :: r_smear = 0.0_wp
    !> Slater charge smearing beta relation
    Integer(Kind=wi), Public              :: b_smear = BETA_ORIGINAL
    !> No electrostatics switch
    Logical, Public                       :: no_elec = .false.
    !> specifies if the correction terms (force shift/reaction field) have been initialised
    Logical, Public                       :: initialised = .false.
    Logical, Public                       :: multipolar = .false.
    Type(mpole_type), Public              :: mpoles
    Integer, Public                       :: num_mpoles = 0
    !> How many mpole derivatives
    Integer, Dimension(0:0), Public       :: nmpole_derivs = [1]
    !> My mpole derivatives
    Integer, Dimension(3,1,0:0), Public   :: mpole_derivs = reshape([0,0,0],[3,1,1])
    !> Damped or not?
    Logical, Public                       :: damp = .false.
    !> Damping distance
    Real(Kind=wp), Public                 :: damping = 0.0_wp
    !> Relative dielectric constant
    Real(Kind=wp), Public                 :: eps = 1.0_wp
    !> Bjerrum length
    Real(Kind=wp), Public                 :: len_bjer = 0.0_wp
    Logical, Public                       :: lecx = .false.
    Integer, Public                       :: nstfce = 1
    Real(Kind=wp), Public                 :: force_shift = 0.0_wp, energy_shift = 0.0_wp
    Real(Kind=wp), Dimension(0:2), Public :: reaction_field = 0.0_wp
    Type(interp_table), Public             :: erfc_over_r, erfc_gamma

  Contains
    Procedure, Public                     :: init_erf_tables
    Procedure, Public                     :: erfcgen
  End Type electrostatic_type

contains
  Subroutine init_erf_tables(electro, nsamples)
    Class(electrostatic_type), Intent(InOut) :: electro
    Integer,                   Intent(In   ) :: nsamples

    electro%erfc_over_r%nsamples = nsamples
    electro%erfc_gamma%nsamples = nsamples

  End Subroutine init_erf_tables

  Subroutine erfcgen(electro,rcut,alpha)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 routine for generating interpolation tables for erfc/r and
    !! -(d/dr erfc/r) / r. For use with Ewald sum.
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester december 1994
    !! amended   - i.t.todorov february 2016
    !! amended   - j.s.wilkins september 2019
    !! ammended  - h.l.devereux january 2025
    !! contrib   - b.t.speake July 2024 - charge smearing
    !!-----------------------------------------------------------------------
    Implicit None

    Class(electrostatic_type),   Intent (InOut) :: electro
    Real(Kind=wp),               Intent(In   )  :: rcut,alpha

    if (electro%erfc_over_r%initialised .and. electro%erfc_gamma%initialised) return

    call electro%erfc_over_r%init(rcut, erfc_ar_over_r)
    call electro%erfc_gamma%init(rcut, erfc_gamma)

  contains
    Function erfc_ar_over_r(rrr)
      !!-----------------------------------------------------------------------
      !!
      !! dl_poly_4 Function to calculate erfc / r.
      !!
      !! copyright - daresbury laboratory
      !! author    - t.forester december 1994
      !! amended   - i.t.todorov february 2016
      !! amended   - j.s.wilkins september 2019
      !! ammended  - h.l.devereux january 2025
      !! contrib   - b.t.speake July 2024 - charge smearing
      !!-----------------------------------------------------------------------
      Real(Kind=wp) :: rrr
      Real(Kind=wp) :: erfc_ar_over_r

      Select Case (electro%smear)
      Case (SMEARING_LINEAR)
        If (rrr < (2 * electro%r_smear)) Then
          erfc_ar_over_r = calc_erfc_n(alpha*rrr)  - linear_smearing_pot(rrr, electro%r_smear)
        Else
          erfc_ar_over_r = calc_erfc_n(alpha*rrr)
        End If

      Case (SMEARING_SLATER_EXP)
        erfc_ar_over_r = calc_erfc_n(alpha*rrr) - slater_exp_smearing_pot(rrr, electro%r_smear, electro%b_smear)

      Case (SMEARING_SLATER_TRUNCATED)
        erfc_ar_over_r = calc_erfc_n(alpha*rrr) - slater_apprx_smearing_pot(rrr, electro%r_smear, electro%b_smear)

      Case (SMEARING_GAUSSIAN)
        If (electro%r_smear == (1.0_wp / (2.0_wp*alpha))) Then
          erfc_ar_over_r = 0.0_wp
        Else
          erfc_ar_over_r = calc_erfc_n(alpha*rrr) - calc_erfc_n(rrr / (2.0_wp*electro%r_smear))
        End If

      Case (SMEARING_GAUSSIAN_EQUAL)
        erfc_ar_over_r = 0.0_wp
        Return

      Case Default
        erfc_ar_over_r = calc_erfc_n(alpha*rrr)
      End Select

      erfc_ar_over_r = erfc_ar_over_r / rrr
    End Function erfc_ar_over_r

    Function erfc_gamma(rrr)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 function to calculate -(d/dr erfc/r) / r.
    !!
    !! copyright - daresbury laboratory
    !! author    - t.forester december 1994
    !! amended   - i.t.todorov february 2016
    !! amended   - j.s.wilkins september 2019
    !! ammended  - h.l.devereux january 2025
    !! contrib   - b.t.speake July 2024 - charge smearing
    !!-----------------------------------------------------------------------
      Real(Kind=wp) :: rrr, inv_rsq, inv_r
      Real(Kind=wp) :: erfc_gamma

      inv_r = 1 / rrr

      Select Case (electro%smear)
      Case (SMEARING_LINEAR)
        If (rrr < (2.0_wp * electro%r_smear)) Then
          erfc_gamma = ((calc_erfc_n(alpha*rrr) * inv_r) - alpha*calc_erfc_deriv_n(alpha*rrr)) - &
                                    linear_smearing_force(rrr, electro%r_smear) * inv_r
        Else
          erfc_gamma = ((calc_erfc_n(alpha*rrr) * inv_r) - alpha*calc_erfc_deriv_n(alpha*rrr))
        End If

      Case (SMEARING_SLATER_EXP)
        erfc_gamma = (calc_erfc_n(alpha*rrr) - slater_exp_smearing_force(rrr, electro%r_smear, electro%b_smear)) * inv_r
        erfc_gamma = erfc_gamma  - alpha*calc_erfc_deriv_n(alpha*rrr)

      Case (SMEARING_SLATER_TRUNCATED)
        erfc_gamma = (calc_erfc_n(alpha*rrr) - slater_apprx_smearing_force(rrr, electro%r_smear, electro%b_smear)) * inv_r
        erfc_gamma = erfc_gamma  - alpha*calc_erfc_deriv_n(alpha*rrr)

      Case (SMEARING_GAUSSIAN)
        If (electro%r_smear == (1.0_wp / (2.0_wp*alpha))) Then
          erfc_gamma = 0.0_wp
        Else
          erfc_gamma = ((calc_erfc_n(alpha*rrr) * inv_r) - alpha*calc_erfc_deriv_n(alpha*rrr))
          erfc_gamma = erfc_gamma - (calc_erfc_n(rrr / (2.0_wp*electro%r_smear)) * inv_r)
          erfc_gamma = erfc_gamma + (calc_erfc_deriv_n(rrr / (2.0_wp*electro%r_smear)) / (2 * electro%r_smear))
        End If

      Case (SMEARING_GAUSSIAN_EQUAL)
        erfc_gamma = 0.0_wp
        Return

      Case Default
        erfc_gamma = ((calc_erfc_n(alpha*rrr) * inv_r) - alpha*calc_erfc_deriv_n(alpha*rrr))
      End Select

      inv_rsq = inv_r * inv_r
      erfc_gamma = erfc_gamma * inv_rsq
    End Function erfc_gamma


  End Subroutine erfcgen

End Module electrostatic
