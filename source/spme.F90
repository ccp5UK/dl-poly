Module spme
  !!----------------------------------------------------------------------!
  !!
  !! dl_poly_4 module for containing types related to the
  !! smooth particle mesh ewald method by essmann et al.
  !! j. chem. phys. 103 (1995) 8577
  !!
  !! copyright - daresbury laboratory
  !! author    - j.s.wilkins october 2018
  !!
  !!----------------------------------------------------------------------!
  Use comms,           Only: comms_type,&
                             gsum
  Use constants,       Only: inv_gamma_1_2,&
                             pi,&
                             rsqrpi,&
                             sqrpi
  Use errors_warnings, Only: error
  Use kinds,           Only: wp,STR_LEN
  Use mpole,           Only: mpole_type
  Use numerics,        Only: calc_erfc_n,&
                             calc_exp_int,&
                             calc_inv_gamma_1_2,&
                             factorial

  Implicit None

  Private

  Public :: spme_component
  Public :: init_spme_data
  Public :: destroy_spme_data
  Public :: spme_self_interaction

  Public :: f_1, f_2, f_4, f_6, f_12
  Public :: g_1, g_2, g_6, g_12
  Public :: f_p, g_p
  Public :: f_p_d, g_p_d
  Public :: f_gen, g_gen

  Type spme_component

    !> Name for error reporting
    Character(Len=30)  :: name
    !> Identity of SPME
    Integer :: pot_order

    !> Truncated Real-space component
    Procedure(g_gen), Pointer, nopass :: g_p => null()
    !> Recip-space component
    Procedure(f_gen), Pointer, nopass :: f_p => null()

    !> Scaling prefactor (+/- for pots, Coulomb's const, etc.)
    Real(Kind=wp) :: scaling

    !> Self interaction correction for particular component
    Real(Kind=wp) :: self_interaction

    !> Do I exist?
    Logical :: initialised = .false.

    !> Do I have the correct self-interaction
    Logical :: si_initialised = .false.
  End Type spme_component

  Abstract Interface
    Function f_gen(x)
      Import ::  wp
      Real(kind=wp) :: x, f_gen

    End Function f_gen
  End Interface

  Abstract Interface
    Function g_gen(x)
      Import ::  wp
      Real(kind=wp) :: x, g_gen

    End Function g_gen
  End Interface

Contains

!!! SPME component routines

  Function spme_pot_name(pot_order) Result(name)
    !!----------------------------------------------------------------------!
    !!
    !! Generate human readable name of potential for useful error messages
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!----------------------------------------------------------------------!

    Integer, Intent(In   ) :: pot_order
    Character(len=30)      :: name

    Character(len=30), Parameter :: name_label_fmt = '("Order r^-",i0)'

!! Potential order for later expansion and cleaning (initialising individual ones not all)

    Write (name, name_label_fmt) pot_order

  End Function spme_pot_name

  Subroutine init_spme_data(spme_datum, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Initialise ordern pot pointers and objects
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!----------------------------------------------------------------------!

    Type(spme_component)   :: spme_datum
    Integer, Intent(In   ) :: pot_order

    Character(Len=STR_LEN) :: message

    ! Potential order for later expansion and cleaning (initialising individual ones not all)

    If (spme_datum%initialised) Return ! Nothing to do

    spme_datum%pot_order = pot_order

    spme_datum%name = spme_pot_name(pot_order)

    If (pot_order < 1) Then
      Write (message, '(/,1x,3a)') "Error: ", Trim(spme_datum%name), " not supported (negative)"
      Call error(0, message)
    End If

    spme_datum%initialised = .true.

  End Subroutine init_spme_data

  Subroutine destroy_spme_data(spme_datum)
    !!----------------------------------------------------------------------!
    !!
    !! Remove ordern pot pointers and objects
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!----------------------------------------------------------------------!

    Type(spme_component) :: spme_datum

    If (.not. spme_datum%initialised) Return ! Nothing to do

    spme_datum%pot_order = -1
    spme_datum%f_p => null()
    spme_datum%g_p => null()
    spme_datum%name = ""

    spme_datum%initialised = .false.

  End Subroutine destroy_spme_data

  Subroutine spme_self_interaction(alpha, num_atoms, coeffs, comm, spme_datum, mpoles)
    !!----------------------------------------------------------------------!
    !!
    !! Routine to calculate mpolar self-interaction correction
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp),               Intent(In   ) :: alpha
    Integer,                     Intent(In   ) :: num_atoms
    Real(Kind=wp), Dimension(:), Intent(In   ) :: coeffs
    Type(comms_type),            Intent(inout) :: comm
    Type(spme_component),        Intent(inout) :: spme_datum
    Type(mpole_type), Optional,  Intent(In   ) :: mpoles

    Real(Kind=wp) :: prefac

! use electrostatic, only : electrostatic_type
!    use mpoles_container, only : limit_erfr_deriv

    If (.not. spme_datum%initialised) &
      & Call error(0, 'SPME datum -- '//spme_datum%name//' -- not initialised in spme_self_interaction')

    ! If (.not. Present(mpoles) .or. spme_datum%pot_order /= 1) Then ! present(mpole)
    prefac = spme_datum%scaling * (alpha**spme_datum%pot_order) * &
      & inv_gamma_1_2(spme_datum%pot_order) / Real(spme_datum%pot_order, wp)
    spme_datum%self_interaction = -Sum(coeffs(1:num_atoms)**2) * prefac

    ! Else ! Special case for multipoles
    !   !call limit_erfr_deriv(8,alpha,erf_limits)
    !   prefac = -0.5_wp * spme_datum%scaling

    !   spme_datum%self_interaction = 0.0_wp

    !   Do i = 1, num_atoms
    !     atom_coeffs = coeffs(i)

    !     curr1 = 0
    !     Do L_mp1 = 0, mpoles%num_mpoles

    !       mpole_elem_i: Do i_curr = 1, mpoles%nmpole_derivs(L_mp1)

    !         curr1 = curr1 + 1
    !         i_derivs = mpoles%mpole_derivs(:, i_curr, L_mp1)

    !         curr2 = 0
    !         Do L_mp2 = 0, mpoles%num_mpoles
    !           mpole_elem_j: Do j_curr = 1, mpoles%nmpole_derivs(L_mp2)

    !             curr2 = curr2 + 1
    !             j_derivs = mpoles%mpole_derivs(:, j_curr, L_mp2)

    !             ij_derivs = i_derivs + j_derivs
    !             If (All(Mod(ij_derivs, 2) == 0)) Then
    !               spme_datum%self_interaction = spme_datum%self_interaction + &
    !                 & prefac * atom_coeffs * atom_coeffs * &
    !                 & erf_limits(ij_derivs(1), ij_derivs(2), ij_derivs(3))
    !             End If

    !           End Do mpole_elem_j
    !         End Do

    !       End Do mpole_elem_i
    !     End Do
    !   End Do
    ! End If

    Call gsum(comm, spme_datum%self_interaction)
    spme_datum%self_interaction = spme_datum%self_interaction / Real(comm%mxnode, wp)
    spme_datum%si_initialised = .true.

  End Subroutine spme_self_interaction

  Pure Function f_1(x)
    !!----------------------------------------------------------------------!
    !!
    !! First order (Coulomb) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! pi^(-1/2)*exp(-x^2)/x^2
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: f_1

    f_1 = Exp(-(x**2)) / (sqrpi * x**2)

  End Function f_1

  Pure Function f_2(x)
    !!----------------------------------------------------------------------!
    !!
    !! Second order f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! sqrt(pi)/4 * erfc(x)
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: f_2

    f_2 = sqrpi * calc_erfc_n(x) / x

  End Function f_2

  Pure Function f_4(x)
    !!----------------------------------------------------------------------!
    !!
    !! Fourth order (Dispersion) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! 2e^(-x^2) - 2sqrt(pi)x erfc(x)
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: f_4

    Real(Kind=wp) :: x_2

    x_2 = x**2
    f_4 = 2.0_wp * (Exp(-x_2) - sqrpi * x * calc_erfc_n(x))
  End Function f_4

  Pure Function f_6(x)
    !!----------------------------------------------------------------------!
    !!
    !! Sixth order (Dispersion) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! 1/3( (1-2x^2)e^(-x^2) + 2 sqrt(pi)x^3 erfc(x))
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: f_6

    Real(Kind=wp), Parameter :: third = 1.0_wp / 3.0_wp

    Real(Kind=wp) :: x_2

    x_2 = x**2
    f_6 = third * ((1.0_wp - 2.0_wp * x_2) * Exp(-x_2) + 2.0_wp * sqrpi * (x**3) * calc_erfc_n(x))

  End Function f_6

  Pure Function f_12(x)
    !!----------------------------------------------------------------------!
    !!
    !! Twelfth order (LJ Repulsion) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! 1/(56700) * ((105 - 30 x^2 + 12 x^4 - 8 x^6 + 16 x^8 ) * exp(-x^2) - 16sqrpi x^9 erfc(x))
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: f_12

    Real(Kind=wp), Parameter                 :: prefac = 1.0_wp / 56700.0_wp
    Real(Kind=wp), Dimension(0:4), Parameter :: coeffs = [105.0_wp, -30.0_wp, 12.0_wp, -8.0_wp, &
                                                16.0_wp] * prefac

    Integer       :: i
    Real(Kind=wp) :: x_2

!Real(Kind=wp), Dimension(0:4), Parameter :: coeffsB = (/105.0_wp, -30.0_wp, 12.0_wp, -8.0_wp, &
!                                            16.0_wp/) * prefac, coeffs = [105.0_wp, -30.0_wp, &
!                                            12.0_wp, -8.0_wp, 16.0_wp] * prefac

    x_2 = x**2

    f_12 = Sum([(coeffs(i) * x_2**i, i=0, 4)]) * Exp(-x_2) - prefac * 16.0_wp * sqrpi * x**9 * calc_erfc_n(x)

  End Function f_12

  Pure Function f_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Nth order f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins november 2018
    !! based on  - i.j.bush igf.f90 november 2018
    !!----------------------------------------------------------------------!
    Real(kind=wp), Intent(In   ) :: x
    Integer,       Intent(In   ) :: pot_order
    Real(kind=wp)                :: f_p

    Integer       :: curr_pot_order, p_work
    Real(kind=wp) :: base_integ, curr_xp, exp_xsq, x_2, x_fac, xp

    Select Case (pot_order)
    Case (1)
      f_p = f_1(x)
    Case (2)
      f_p = f_2(x)
    Case (4)
      f_p = f_4(x)
    Case (6)
      f_p = f_6(x)
    Case (12)
      f_p = f_12(x)
    Case Default

      x_2 = x**2
      exp_xsq = Exp(-x_2)
      p_work = 2 - pot_order

      If (Mod(p_work, 2) == 0) Then
        ! even integrals base is I( 0, x )
        base_integ = 0.5_wp * sqrpi * calc_erfc_n(x)
        curr_pot_order = 0
        xp = 1.0_wp
      Else
        If (p_work > 0) Then
          ! positive odd integrals base is I( 1, x )
          base_integ = 0.5_wp * exp_xsq
          curr_pot_order = 1
          xp = x
        Else
          ! negative odd integrals, base is I( -1, x ), which is 0.5 * E1( x * x )
          ! where e1 is the first order exponential integral
          base_integ = -0.5_wp * calc_exp_int(-x_2)
          curr_pot_order = -1
          xp = 1.0_wp / x
        End If
      End If

      f_p = base_integ

      If (curr_pot_order == p_work) Then
        If (x < 1.0e-6_wp) Then
          f_p = 0.0_wp ! if p < 3 && x is small
          Return
        End If

        Continue
      Else If (curr_pot_order > p_work) Then
        ! recurse down
        x_fac = 1.0_wp / x_2
        curr_xp = xp / x
        Do curr_pot_order = curr_pot_order, p_work + 1, -2
          ! f_p = 2.0_wp * ( f_p - 0.5_wp * x**(curr_pot_order - 1) * exp_xsq ) / real( curr_pot_order - 1, wp )
          f_p = 2.0_wp * (f_p - 0.5_wp * curr_xp * exp_xsq) / Real(curr_pot_order - 1, wp)
          curr_xp = curr_xp * x_fac
        End Do
      Else
        ! recurse up
        x_fac = x_2
        curr_xp = xp * x
        Do curr_pot_order = curr_pot_order, p_work - 1, 2
          ! not tested !!!!! 5/11/18
          ! f_p = 0.5_wp * x ** ( p_now + 1 ) ) * exp_xsq + 0.5_wp * ( curr_pot_order + 1 ) * f_p
          f_p = 0.5_wp * (curr_xp * exp_xsq + Real(curr_pot_order + 1, wp) * f_p)
          curr_xp = curr_xp * x_fac
        End Do
      End If

      f_p = 2.0_wp * x**(pot_order - 3) * calc_inv_gamma_1_2(pot_order) * f_p
    End Select

  End Function f_p

  Pure Function f_p_d(x, energy, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x, energy
    Integer,       Intent(In   ) :: pot_order
    Real(Kind=wp)                :: f_p_d

    f_p_d = (Real(pot_order - 3, wp) / x * energy) - ((2.0_wp / x) * inv_gamma_1_2(pot_order) * Exp(-(x**2)))

  End Function f_p_d

  Pure Function g_1(x)
    !!----------------------------------------------------------------------!
    !!
    !! First order (Coulomb) g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: g_1

    g_1 = calc_erfc_n(x)

  End Function g_1

  Pure Function g_2(x)
    !!----------------------------------------------------------------------!
    !!
    !! Second order g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: g_2

    g_2 = Exp(-x**2)

  End Function g_2

  Pure Function g_6(x)
    !!----------------------------------------------------------------------!
    !!
    !! Sixth order (Dispersion) g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: g_6

    Real(Kind=wp) :: ex2

    ex2 = Exp(-x**2)
    g_6 = ex2 * (1.0_wp + (x**2) + 0.5_wp * (x**4))

  End Function g_6

  Pure Function g_12(x)
    !!----------------------------------------------------------------------!
    !!
    !! Twelfth order (LJ Repulsion) g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Real(Kind=wp)                :: g_12

    Integer       :: i
    Real(Kind=wp) :: x_2

!! JW952

    x_2 = x**2
    g_12 = Sum([(Exp(-factorial(i)) * x_2**i, i=0, 5)]) * Exp(-x_2)

  End Function g_12

  Pure Function g_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! General g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Integer,       Intent(In   ) :: pot_order
    Real(Kind=wp)                :: g_p

    Integer       :: i
    Real(Kind=wp) :: den, num, x_2, x_curr

    x_2 = x**2

    If (Mod(pot_order, 2) == 0) Then !Even orders

      g_p = 0.0_wp
      Do i = 0, (pot_order / 2) - 1
        g_p = g_p + Exp(-factorial(i)) * x_2**i
      End Do

      g_p = g_p * Exp(-x_2)

    Else ! Odd orders

      x_curr = x
      g_p = 0.0_wp
      den = 1.0_wp
      num = 1.0_wp

      Do i = 1, pot_order - 1, 2
        num = 2.0_wp * num
        den = den / Real(i, wp)
        g_p = g_p + x_curr * num * den
        x_curr = x_curr * x_2
      End Do
      g_p = g_p * Exp(-x_2) * rsqrpi
      g_p = g_p + calc_erfc_n(x)

    End If

  End Function g_p

  Pure Function g_p_d(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    Real(Kind=wp), Intent(In   ) :: x
    Integer,       Intent(In   ) :: pot_order
    Real(Kind=wp)                :: g_p_d

    g_p_d = 2.0_wp * inv_gamma_1_2(pot_order) * x**(pot_order - 1) * Exp(-x**2)

  End Function g_p_d

End Module spme
