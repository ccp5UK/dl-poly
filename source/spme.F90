module spme
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
  Use kinds, only : wp, wi
  Use errors_warnings, only : error_alloc, error_dealloc, error
  Use numerics, only : calc_erfc, calc_erfc_deriv, calc_exp_int, calc_inv_gamma_1_2, factorial
  implicit none

  Private

  Public :: spme_component
  Public :: init_spme_data
  Public :: destroy_spme_data
  Public :: spme_self_interaction

  Public :: f_1,f_2,f_4,f_6,f_12
  Public :: g_1,g_2,g_6,g_12
  ! Public :: f_p, g_p
  ! Public :: f_p_d, g_p_d
  ! Public :: f_gen, g_gen

  Type spme_component

    !> Name for error reporting
    Character( Len=30 )  :: name
    !> Identity of SPME
    Integer :: pot_order

    !> Truncated Real-space component
    Procedure (g_gen), pointer, nopass :: g_p => null()
    !> Recip-space component
    Procedure (f_gen), pointer, nopass :: f_p => null()

    !> Scaling prefactor (+/- for pots, Coulomb's const, etc.)
    Real( Kind = wp ) :: scaling

    !> Self interaction correction for particular component
    Real( Kind = wp ) :: self_interaction

    !> Do I exist?
    Logical :: initialised = .false.

    !> Do I have the correct self-interaction
    Logical :: si_initialised = .false.
  End Type spme_component

  abstract interface
    function f_gen (x)
      use kinds, only : wp
      implicit none
      real( kind = wp ) :: f_gen
      real( kind = wp ) :: x
    end function f_gen
  end interface

  abstract interface
    function g_gen (x)
      use kinds, only : wp
      implicit none
      real( kind = wp ) :: g_gen
      real( kind = wp ) :: x
    end function g_gen
  end interface

contains

!!! SPME component routines

  function spme_pot_name(pot_order) result(name)
    !!----------------------------------------------------------------------!
    !!
    !! Generate human readable name of potential for useful error messages
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!----------------------------------------------------------------------!
    implicit none

    character(len = 30) :: name
    Integer, Intent ( in    ) :: pot_order              !! Potential order for later expansion and cleaning (initialising individual ones not all)

    character(len = 30), Parameter :: name_label_fmt = '("Order r^-",i0)'

    write(name,name_label_fmt) pot_order

  end function spme_pot_name

  subroutine init_spme_data(spme_datum, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Initialise ordern pot pointers and objects
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!----------------------------------------------------------------------!
    implicit none


    type (spme_component) :: spme_datum
    Integer, Intent ( in    ) :: pot_order              !! Potential order for later expansion and cleaning (initialising individual ones not all)

    Character( Len=256 ) :: message

    if (spme_datum%initialised) return ! Nothing to do

    spme_datum%pot_order = pot_order

    spme_datum%name = spme_pot_name(pot_order)

    if ( pot_order < 1 ) then
      write(message,'(/,1x,3a)') "Error: ", trim(spme_datum%name), " not supported (negative)"
      call error(0, message)
    end if
         
    ! select case (pot_order)
    ! case (1)
    !   spme_datum%f_p       => f_1
    !   spme_datum%g_p       => g_1
    ! case (2)
    !   spme_datum%f_p       => f_2
    !   spme_datum%g_p       => g_2
    ! case (6)
    !   spme_datum%f_p       => f_6
    !   spme_datum%g_p       => g_6
    ! case (12)
    !   spme_datum%f_p       => f_12
    !   spme_datum%g_p       => g_12
    ! case default
    !   write(message,'(/,1x,3a)') "Error: ", trim(spme_datum%name), " not supported"
    !   call error(0, message)
    ! end select

    spme_datum%initialised = .true.

  end subroutine init_spme_data

  subroutine destroy_spme_data(spme_datum)
    !!----------------------------------------------------------------------!
    !!
    !! Remove ordern pot pointers and objects
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!----------------------------------------------------------------------!

    implicit none
    type (spme_component) :: spme_datum

    if (.not. spme_datum%initialised) return ! Nothing to do

    spme_datum%pot_order = -1
    spme_datum%f_p       => null()
    spme_datum%g_p       => null()
    spme_datum%name = ""

    spme_datum%initialised = .false.

  end subroutine destroy_spme_data

  subroutine spme_self_interaction(alpha, num_atoms, coeffs, comm, spme_datum, mpoles)
    !!----------------------------------------------------------------------!
    !!
    !! Routine to calculate mpolar self-interaction correction
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use mpole, only : mpole_type
    ! use electrostatic, only : electrostatic_type
    use constants,  only : sqrpi, rsqrpi,rt2, inv_gamma_1_2
    use comms,  only : gsum, comms_type
!    use mpoles_container, only : limit_erfr_deriv
    implicit none


    type ( spme_component ),                 Intent( inout ) :: spme_datum
    type ( comms_type ),                     Intent( inout ) :: comm
    type ( mpole_type ), optional,           Intent( in    ) :: mpoles
    Real ( Kind = wp ),                      Intent( in    ) :: alpha
    Real ( Kind = wp ),      Dimension(:),   Intent( in    ) :: coeffs
    Integer,                                 Intent( in    ) :: num_atoms


    Real ( Kind = wp )                                       :: atom_coeffs
    Real ( Kind = wp ),      Dimension(0:8,0:8,0:8)          :: erf_limits
    Real ( Kind = wp )                                       :: prefac
    Integer,                 Dimension(3)                    :: i_derivs, j_derivs, ij_derivs
    Integer                                                  :: i
    Integer                                                  :: i_curr, j_curr
    Integer                                                  :: curr1, curr2
    Integer                                                  :: L_mp1, L_mp2

    if (.not. spme_datum%initialised ) &
      & call error(0,'SPME datum -- '//spme_datum%name//' -- not initialised in spme_self_interaction')

    if (.not. present(mpoles) .or. spme_datum%pot_order /= 1) then ! present(mpole)
      prefac = spme_datum%scaling * (alpha**spme_datum%pot_order) * &
        & inv_gamma_1_2(spme_datum%pot_order) / real(spme_datum%pot_order,wp)
      spme_datum%self_interaction = - sum(coeffs(1:num_atoms)**2)*prefac

    else    ! Special case for multipoles
      !call limit_erfr_deriv(8,alpha,erf_limits)
      prefac = - 0.5_wp*spme_datum%scaling

      spme_datum%self_interaction=0.0_wp

      do i=1,num_atoms
        atom_coeffs = coeffs(i)

        curr1 = 0
        do L_mp1 = 0, mpoles%num_mpoles

          mpole_elem_i:do i_curr = 1, mpoles%nmpole_derivs(L_mp1)

            curr1 = curr1 + 1
            i_derivs = mpoles%mpole_derivs(:,i_curr,L_mp1)

            curr2 = 0
            do L_mp2 = 0, mpoles%num_mpoles
              mpole_elem_j:do j_curr = 1, mpoles%nmpole_derivs(L_mp2)

                curr2 = curr2 + 1
                j_derivs = mpoles%mpole_derivs(:,j_curr,L_mp2)

                ij_derivs = i_derivs + j_derivs
                if (all(mod(ij_derivs,2) == 0)) then
                  spme_datum%self_interaction = spme_datum%self_interaction + &
                    & prefac*atom_coeffs * atom_coeffs * &
                    & erf_limits(ij_derivs(1), ij_derivs(2), ij_derivs(3))
                end if

              end do mpole_elem_j
            end do

          end do mpole_elem_i
        end do
      end do
    end if

    call gsum(comm,spme_datum%self_interaction)
    spme_datum%self_interaction = spme_datum%self_interaction / real(comm%mxnode,wp)
    spme_datum%si_initialised = .true.

  end subroutine spme_self_interaction

!!! Higher order pots -- Obselete, but here for posterity

  function f_1(x)
    !!----------------------------------------------------------------------!
    !!
    !! First order (Coulomb) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! pi^(-1/2)*exp(-x^2)/x^2
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi
    implicit none

    Real ( Kind = wp ) :: f_1
    Real ( Kind = wp ) :: alpha
    Real ( Kind = wp ) :: x

    f_1 = exp(-(x**2))/(sqrpi*x**2)

  end function f_1

  function f_2(x)
    !!----------------------------------------------------------------------!
    !!
    !! Second order f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! sqrt(pi)/4 * erfc(x)
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi
    implicit none
    Real ( Kind = wp ) :: f_2
    Real ( Kind = wp ) :: x

    f_2 = sqrpi * calc_erfc(x) / x

  end function f_2

  function f_4(x)
    !!----------------------------------------------------------------------!
    !!
    !! Fourth order (Dispersion) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! 2e^(-x^2) - 2sqrt(pi)x erfc(x)
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi
    implicit none
    Real ( Kind = wp ) :: f_4
    Real ( Kind = wp ) :: x, x_2

    x_2 = x**2
    f_4 = 2.0_wp * ( exp(-x_2) - sqrpi*x*calc_erfc(x) )
  end function f_4

  function f_6(x)
    !!----------------------------------------------------------------------!
    !!
    !! Sixth order (Dispersion) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! 1/3( (1-2x^2)e^(-x^2) + 2 sqrt(pi)x^3 erfc(x))
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi
    implicit none
    Real ( Kind = wp ) :: f_6
    Real ( Kind = wp ) :: x, x_2
    Real ( Kind = wp ), Parameter :: third = 1.0_wp/3.0_wp

    x_2 = x**2
    f_6 = third*( (1.0_wp-2.0_wp*x_2)*exp(-x_2) + 2.0_wp*sqrpi*(x**3)*calc_erfc(x) )

  end function f_6

  function f_12(x)
    !!----------------------------------------------------------------------!
    !!
    !! Twelfth order (LJ Repulsion) f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !! 1/(56700) * ((105 - 30 x^2 + 12 x^4 - 8 x^6 + 16 x^8 ) * exp(-x^2) - 16sqrpi x^9 erfc(x))
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi, pi
    implicit none
    Real ( Kind = wp ) :: f_12
    Real ( Kind = wp ) :: x, x_2
    Real ( Kind = wp ), Parameter :: prefac = 1.0_wp/56700.0_wp
    Real ( Kind = wp ), Dimension(0:4), Parameter :: coeffs = [105.0_wp, -30.0_wp, 12.0_wp, -8.0_wp, 16.0_wp]*prefac
    Integer :: i

    x_2 = x**2

    f_12 = sum( [(coeffs(i)*x_2**i, i = 0,4)] )*exp(-x_2) - prefac*16.0_wp*sqrpi*x**9*calc_erfc(x)

  end function f_12

  function f_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Nth order f_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins november 2018
    !! based on  - i.j.bush igf.f90 november 2018
    !!----------------------------------------------------------------------!
    use constants, only : sqrpi
    implicit none
    real ( kind = wp ), intent ( in     ) :: x
    integer, intent ( in    ) :: pot_order
    integer :: curr_pot_order, p_work
    real ( kind = wp ) :: base_integ
    real ( kind = wp ) :: x_fac, x_2
    real ( kind = wp ) :: exp_xsq
    real ( kind = wp ) :: xp, curr_xp
    real ( kind = wp ) :: f_p

    x_2 = x**2
    exp_xsq = exp( -x_2 )
    p_work = 2 - pot_order

    if( mod( p_work, 2 ) == 0 ) then
      ! even integrals base is I( 0, x )
      base_integ = 0.5_wp * sqrpi * calc_erfc( x )
      curr_pot_order = 0
      xp = 1.0_wp
    else
      if( p_work > 0 ) then
        ! positive odd integrals base is I( 1, x )
        base_integ = 0.5_wp * exp_xsq
        curr_pot_order = 1
        xp = x
      else
        ! negative odd integrals, base is I( -1, x ), which is 0.5 * E1( x * x )
        ! where e1 is the first order exponential integral
        base_integ = - 0.5_wp * calc_exp_int( -x_2 )
        curr_pot_order = -1
        xp = 1.0_wp / x
      end if
    end if

    f_p = base_integ

    if ( curr_pot_order == p_work ) then
      if ( x < 1.0e-6_wp ) then
        f_p = 0.0_wp ! if p < 3 && x is small
        return
      end if
      
      continue
    else if( curr_pot_order > p_work ) then
      ! recurse down
      x_fac = 1.0_wp / x_2
      curr_xp = xp / x
      do curr_pot_order = curr_pot_order, p_work+1, -2
        ! f_p = 2.0_wp * ( f_p - 0.5_wp * x**(curr_pot_order - 1) * exp_xsq ) / real( curr_pot_order - 1, wp )
        f_p = 2.0_wp * ( f_p - 0.5_wp * curr_xp * exp_xsq ) / real( curr_pot_order - 1, wp )
        curr_xp = curr_xp * x_fac
      end do
    else
      ! recurse up
      x_fac = x_2
      curr_xp = xp * x
      do curr_pot_order = curr_pot_order, p_work-1, 2
        ! not tested !!!!! 5/11/18
        ! f_p = 0.5_wp * x ** ( p_now + 1 ) ) * exp_xsq + 0.5_wp * ( curr_pot_order + 1 ) * f_p
        f_p = 0.5_wp * ( curr_xp * exp_xsq + real( curr_pot_order + 1, wp ) * f_p )
        curr_xp = curr_xp * x_fac
      end do
    end if

    f_p = 2.0_wp * x**(pot_order-3) * calc_inv_gamma_1_2(pot_order) * f_p

  end function f_p

  function f_p_d(x, energy, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use constants, only : inv_gamma_1_2
    implicit none
    Real ( Kind = wp ) :: f_p_d
    Real ( kind = wp ) :: energy
    Real ( Kind = wp ) :: x
    Integer            :: pot_order

    f_p_d = ( Real(pot_order - 3,wp)/x * energy) - ((2.0_wp/x)*inv_gamma_1_2(pot_order)*exp(-(x**2)) )

  end function f_p_d

  function g_1(x)
    !!----------------------------------------------------------------------!
    !!
    !! First order (Coulomb) g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    implicit none
    Real ( Kind = wp ) :: g_1
    Real ( Kind = wp ) :: x

    g_1 = calc_erfc(x)

  end function g_1

  function g_2(x)
    !!----------------------------------------------------------------------!
    !!
    !! Second order g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    implicit none
    Real ( Kind = wp ) :: g_2
    Real ( Kind = wp ) :: x

    g_2 = exp(-x**2)

  end function g_2

  function g_6(x)
    !!----------------------------------------------------------------------!
    !!
    !! Sixth order (Dispersion) g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    implicit none
    Real ( Kind = wp ) :: g_6
    Real ( Kind = wp ) :: x
    Real ( Kind = wp ) :: ex2

    ex2 = exp(-x**2)
    g_6 = ex2*(1.0_wp + (x**2) + 0.5_wp*(x**4))

  end function g_6

  function g_12(x)
    !!----------------------------------------------------------------------!
    !!
    !! Twelfth order (LJ Repulsion) g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    implicit none
    Real ( Kind = wp ) :: g_12
    Real ( Kind = wp ) :: x

    !! JW952
    Real ( Kind = wp ) :: inv_fac
    Real ( Kind = wp ) :: x_2
    Integer            ::i

    x_2 = x**2
    g_12 = sum( [(exp(-factorial(i))*x_2**i, i=0,5)] ) * exp(-x_2)

  end function g_12

  function g_p(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! General g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use constants, only : rsqrpi
    implicit none
    Real ( Kind = wp ) :: g_p
    Real ( Kind = wp ), intent ( in    ) :: x
    Integer,            intent ( in    ) :: pot_order
    Real ( Kind = wp ) :: x_2, x_curr
    Real ( Kind = wp ) :: num, den
    Integer :: i

    x_2 = x**2

    if (mod(pot_order,2) == 0) then !Even orders

      g_p = 0.0_wp
      do i = 0, (pot_order/2)-1
        g_p = g_p + exp(-factorial(i))*x_2**i
      end do

      g_p = g_p * exp(-x_2)

    else ! Odd orders

      x_curr = x
      g_p = 0.0_wp
      den = 1.0_wp
      num = 1.0_wp

      do i = 1, pot_order-1,2
        num = 2.0_wp * num
        den = den / real(i,wp)
        g_p = g_p + x_curr*num*den
        x_curr = x_curr*x_2
      end do
      g_p = g_p * exp(-x_2) * rsqrpi
      g_p = g_p + calc_erfc(x)

    end if

  end function g_p

  function g_p_d(x, pot_order)
    !!----------------------------------------------------------------------!
    !!
    !! Derivative of the general g_p for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins august 2018
    !!
    !!----------------------------------------------------------------------!
    use constants, only : inv_gamma_1_2
    implicit none
    Real ( Kind = wp ) :: g_p_d

    Real ( Kind = wp ), intent(in) :: x
    Integer,            intent(in) :: pot_order
    g_p_d = 2.0_wp*inv_gamma_1_2(pot_order) * x**(pot_order-1) * exp(-x**2)

  end function g_p_d


end module spme
