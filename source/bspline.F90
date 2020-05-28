Module bspline
  !!----------------------------------------------------------------------!
  !!
  !! dl_poly_4 module for containing types and functions
  !! relating to the creation of bsplines
  !!
  !! copyright - daresbury laboratory
  !! author    - j.s.wilkins october 2018
  !!
  !!----------------------------------------------------------------------!
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc
  Use kinds,           Only: wi,&
                             wp
  Use kspace,          Only: kspace_type

  Implicit None
  Private

  Public :: bspline_coeffs_gen, bspline_splines_gen

  !! JW952
  ! Attach to Ewald type?
  Type, Public :: bspline_type
    !> Number of required derivatives
    Integer(Kind=wi), Public :: num_deriv
    !> SPME FFT B-spline order
    Integer(Kind=wi), Public :: num_splines
    !> SPME FFT B-spline order when padding radius > 0
    Integer(Kind=wi), Public :: num_spline_pad
    !> And another one
    Integer(Kind=wi), Public :: num_spline_padded
    !> B-spline coefficients
    Complex(Kind=wp), Dimension(:, :), Allocatable, Public :: coefficients
    !> Precalculated bb*
    Real(Kind=wp), Dimension(:, :), Allocatable, Public :: norm2
    !> Spline derivative
    Real(Kind=wp), Dimension(:, :, :, :), Allocatable, Public :: derivs

    !> Is this bspline initialised correctly
    Logical :: derivs_initialised = .false., coeffs_initialised = .false.
  End Type bspline_type

  Real(Kind=wp), Dimension(:, :), Save, Allocatable                                                 :: ncombk !! Combinations
  Real(Kind=wp), Dimension(:), Save, Allocatable                                                   :: real_no, inv_no !! Real variants to avoid type-casting

Contains
!!! Bspline Routines

  Subroutine bspline_coeffs_gen(kspace, bspline)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine to calculate B-spline coefficients for
    !! Euler exponential splines
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith july 1998
    !!
    !!-----------------------------------------------------------------------

    Type(kspace_type),  Intent(In   ) :: kspace
    Type(bspline_type), Intent(InOut) :: bspline

    Character(Len=256)                          :: message
    Complex(Kind=wp)                            :: temp_spline
    Complex(Kind=wp), Allocatable, Dimension(:) :: ww1, ww2, ww3
    Integer                                     :: fail, i, j, k, n
    Real(Kind=wp), Allocatable, Dimension(:)    :: cspline

    ! Nothing to do
    If (bspline%coeffs_initialised) Return

    ! Setup constants for bspline_splines_gen

    If (Allocated(real_no) .and. Allocated(inv_no)) Then
      Deallocate (real_no, inv_no, Stat=fail)
      If (fail /= 0) Call error_alloc('real_no and inv_no', 'bspline_splines_gen')
    End If
    Allocate (real_no(1:bspline%num_splines), inv_no(1:bspline%num_splines), Stat=fail)
    If (fail /= 0) Call error_alloc('real_no and inv_no', 'bspline_splines_gen')

    Do i = 1, bspline%num_splines
      real_no(i) = Real(i, wp)
      inv_no(i) = 1.0_wp / real_no(i)
    End Do

    If (Allocated(ncombk)) Then
      Deallocate (ncombk, stat=fail)
      If (fail /= 0) Call error_dealloc('ncombk', 'bspline_coeffs_gen')
    End If
    Allocate (ncombk(0:bspline%num_splines, bspline%num_splines), stat=fail)
    If (fail /= 0) Call error_alloc('ncombk', 'bspline_coeffs_gen')

    ncombk = 0.0_wp
    Do n = 1, bspline%num_splines ! If we change spline number, need to recompute coeffs anyway
      Do k = 0, n
        ncombk(k, n) = Product([(real_no(i), i=n - k + 1, n)]) * Product([(inv_no(i), i=1, Max(1, k))])
      End Do
    End Do

    ! Perform validity checks in routine which is only called once!
    If (bspline%num_splines < 2 .or. bspline%num_deriv < 0) Then
      Call error(0, 'Error: illegal spline order in bspline_coeffs')
    End If

    If (bspline%num_deriv > bspline%num_splines - 2) Then ! Can't return this many
      Write (message, '(a,/,2(a,i0))') 'Error: b-spline not high enough order for derivatives.', 'Num Derivatives: ', &
           & bspline%num_deriv, 'Bspline order: ', bspline%num_splines
      Call error(0, message)
    End If

    Allocate (ww1(1:kspace%k_vec_dim(1)), ww2(1:kspace%k_vec_dim(2)), ww3(1:kspace%k_vec_dim(3)), stat=fail)
    If (fail /= 0) Call error_alloc('ww arrays', 'bspline_coeffs_gen')

    Allocate (bspline%coefficients(3, kspace%k_vec_max), stat=fail)
    If (fail /= 0) Call error_alloc('bspline coefficients', 'bspline_coeffs_gen')
    bspline%coefficients = 0.0_wp
    Allocate (bspline%norm2(3, kspace%k_vec_max), stat=fail)
    If (fail /= 0) Call error_alloc('bspline norms', 'bspline_coeffs_gen')

    ! initialise the complex exponential arrays

    Call spl_cexp(kspace%k_vec_dim(1), kspace%k_vec_dim(2), kspace%k_vec_dim(3), ww1, ww2, ww3)

    ! allocate the helper array

    Allocate (cspline(1:bspline%num_splines), stat=fail)
    If (fail /= 0) Call error_alloc('cspline array', 'bspline_coeffs_gen')

    ! calculate B-splines at knots
    cspline = 0.0_wp
    cspline(2) = 1.0_wp

    Do k = 3, bspline%num_splines
      Do j = k, 2, -1
        cspline(j) = (real_no(j - 1) * cspline(j) + real_no(k - j + 1) * cspline(j - 1)) * inv_no(k - 1)
      End Do
    End Do

    ! calculate B-spline coefficients

    Do i = 0, kspace%k_vec_dim(1) - 1
      temp_spline = (0.0_wp, 0.0_wp)

      Do k = 0, bspline%num_splines - 2
        temp_spline = temp_spline + cspline(k + 2) * ww1(Mod(i * k, kspace%k_vec_dim(1)) + 1)
      End Do

      bspline%coefficients(1, i + 1) = ww1(Mod(i * (bspline%num_splines - 1), kspace%k_vec_dim(1)) + 1) / temp_spline
    End Do

    Do i = 0, kspace%k_vec_dim(2) - 1
      temp_spline = (0.0_wp, 0.0_wp)

      Do k = 0, bspline%num_splines - 2
        temp_spline = temp_spline + cspline(k + 2) * ww2(Mod(i * k, kspace%k_vec_dim(2)) + 1)
      End Do

      bspline%coefficients(2, i + 1) = ww2(Mod(i * (bspline%num_splines - 1), kspace%k_vec_dim(2)) + 1) / temp_spline
    End Do

    Do i = 0, kspace%k_vec_dim(3) - 1
      temp_spline = (0.0_wp, 0.0_wp)

      Do k = 0, bspline%num_splines - 2
        temp_spline = temp_spline + cspline(k + 2) * ww3(Mod(i * k, kspace%k_vec_dim(3)) + 1)
      End Do

      bspline%coefficients(3, i + 1) = ww3(Mod(i * (bspline%num_splines - 1), kspace%k_vec_dim(3)) + 1) / temp_spline
    End Do

    ! Calculate magnitude of bspline coefficients

    bspline%norm2 = Real(bspline%coefficients * Conjg(bspline%coefficients), wp)

    Deallocate (cspline, stat=fail)
    If (fail /= 0) Call error_dealloc('cspline arrays', 'bspline_coeffs_gen')
    Deallocate (ww1, ww2, ww3, stat=fail)
    If (fail /= 0) Call error_dealloc('ww arrays', 'bspline_coeffs_gen')

    bspline%derivs_initialised = .false.
    ! end if

    ! nsplines_old = max(nsplines_old,bspline%num_splines)

  End Subroutine bspline_coeffs_gen

  Pure Subroutine bspline_splines_gen(num_atoms, recip_coords, bspline)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine to calculate B-splines for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith july 1998
    !! amended   - i.t.todorov april 2015
    !! reordered - j.s.wilkins november 2018
    !!-----------------------------------------------------------------------
    Integer,                                Intent(In   ) :: num_atoms
    Real(Kind=wp), Dimension(3, num_atoms), Intent(In   ) :: recip_coords
    Type(bspline_type),                     Intent(InOut) :: bspline
    Integer                     :: i, j, k, l, n, s
    Real(Kind=wp)               :: jm1_r, k_r, km1_rr
    Real(Kind=wp), Dimension(3) :: current_bspline_centre, current_bspline_point
    ! Size one larger to avoid overrun (optimising loop)
    Real(Kind=wp), Dimension(3, 1:bspline%num_splines+1) :: current_zero_deriv, current_first_deriv

!! Number of atoms to fill
!! Coordinates of charge centres in reciprocal cell
!! Bspline to be created
!! Origin of current bspline being calculated
!! Current point of the bspline being calculated
!! Sign alternation
!! Real variants to avoid type-casting
!! Loop counters

    ! Perform validity checks in routine which is only called once -- bspline_coeffs_gen !
    s = bspline%num_splines

    ! construct B-splines
    ! Zero initial array -- Avoid hassle
    ! bspline%derivs = 0.0_wp

    ! Reversed order of array for memory access efficiency
    Do i = 1, num_atoms

      current_zero_deriv = 0.0_wp
      ! initializing 2nd order B-spline
      ! for u where (0<u<1) and (1<u<2)

      current_zero_deriv(:, s) = recip_coords(:, i) - Aint(recip_coords(:, i), wp)
      current_zero_deriv(:, s - 1) = 1.0_wp - current_zero_deriv(:, s)

      ! Now on to calculate order k B-spline values at k
      ! points where (0<u<k)

      current_bspline_centre(:) = current_zero_deriv(:, s)

      Do k = s - 2, bspline%num_deriv + 1, -1 ! 3,bspline%num_splines-bspline%num_deriv ! Order of B-spline

        k_r = real_no(s - k + 1)
        km1_rr = inv_no(s - k)
        Do j = k, s - 1 ! Compute order k B-spline at points {k,k-1,...,1}

          jm1_r = real_no(s - j)
          current_bspline_point(:) = current_bspline_centre(:) + jm1_r

          current_zero_deriv(:, j) = (current_bspline_point(:) * current_zero_deriv(:, j) + &
               (k_r - current_bspline_point(:)) * current_zero_deriv(:, j + 1)) * km1_rr
          ! bspline%derivs(:, 0, j, i) = (current_bspline_point(:) * bspline%derivs(:, 0, j, i) &
          !   & + (k_r - current_bspline_point(:)) * bspline%derivs(:, 0, j + 1, i)) * km1_rr

        End Do

        current_zero_deriv(:, s) = current_zero_deriv(:, s) * current_bspline_centre(:) * km1_rr
      End Do

      ! Now compute B-splines for order bspline%num_splines at k points where
      ! (0<u<bspline%num_splines)

      Do l = bspline%num_deriv, 1, -1

        k_r = real_no(s - l + 1)
        km1_rr = inv_no(s - l)
        current_first_deriv(:,:) = 0.0_wp

        Do j = 1, s - 1 !bspline%num_splines,2,-1

          Derivatives of B-splines with order nospl at k-1 points
          sgn = 1.0_wp
          do n = 0,Min(l,s-j)
            current_first_deriv(:, j) = current_first_deriv(:, j) + sgn*ncombk(n, l) * current_zero_deriv(:, j + n)
            sgn = -sgn
          end do
          ! Do n = 0, Min(l, s - j), 2
          !   current_first_deriv(:, j) = current_first_deriv(:, j) + ncombk(n, l) * current_zero_deriv(:, j + n)
          !   current_first_deriv(:, j) = current_first_deriv(:, j) - ncombk(n+1, l) * current_zero_deriv(:, j + n+1)
          ! End Do

          ! Generate current point at a lag behind derivs
          jm1_r = real_no(s - j)
          current_bspline_point(:) = current_bspline_centre(:) + jm1_r
          current_zero_deriv(:, j) = (current_bspline_point(:) * current_zero_deriv(:, j) &
            & + (k_r - current_bspline_point(:)) * current_zero_deriv(:, j + 1)) * km1_rr

        End Do

        bspline%derivs(:, l, 1:s-1, i) = current_first_deriv(:, 1:s-1)
        bspline%derivs(:, l, s, i) = current_zero_deriv(:, s)

        current_zero_deriv(:, s) = current_zero_deriv(:, s) * current_bspline_centre(:) * km1_rr

      End Do

      bspline%derivs(:, 0, :, i) = current_zero_deriv(:, :)

    End Do

  End Subroutine bspline_splines_gen

  Pure Subroutine spl_cexp(ndiv1, ndiv2, ndiv3, ww1, ww2, ww3)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 routine to create complex exponential arrays for b-splines
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith october 1998
    !! amended   - i.t.todorov october 2006
    !!
    !!-----------------------------------------------------------------------
    Use constants, Only: twopi
    Integer,          Intent(In   ) :: ndiv1, ndiv2, ndiv3
    Complex(Kind=wp), Intent(  Out) :: ww1(1:ndiv1), ww2(1:ndiv2), ww3(1:ndiv3)

    Integer       :: i
    Real(Kind=wp) :: arg

    ! initialise complex exponential factors

    ww1(1) = (1.0_wp, 0.0_wp)

    Do i = 1, ndiv1 / 2
      arg = (twopi / Real(ndiv1, wp)) * Real(i, wp)
      ww1(i + 1) = Cmplx(Cos(arg), Sin(arg), Kind=wp)
      ww1(ndiv1 + 1 - i) = Conjg(ww1(i + 1))

    End Do
    ww2(1) = (1.0_wp, 0.0_wp)

    Do i = 1, ndiv2 / 2
      arg = (twopi / Real(ndiv2, wp)) * Real(i, wp)

      ww2(i + 1) = Cmplx(Cos(arg), Sin(arg), Kind=wp)
      ww2(ndiv2 + 1 - i) = Conjg(ww2(i + 1))
    End Do

    ww3(1) = (1.0_wp, 0.0_wp)

    Do i = 1, ndiv3 / 2
      arg = (twopi / Real(ndiv3, wp)) * Real(i, wp)

      ww3(i + 1) = Cmplx(Cos(arg), Sin(arg), Kind=wp)
      ww3(ndiv3 + 1 - i) = Conjg(ww3(i + 1))
    End Do

  End Subroutine spl_cexp

End Module bspline
