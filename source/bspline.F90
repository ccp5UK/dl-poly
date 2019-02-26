module bspline
  !!----------------------------------------------------------------------!
  !!
  !! dl_poly_4 module for containing types and functions
  !! relating to the creation of bsplines
  !!
  !! copyright - daresbury laboratory
  !! author    - j.s.wilkins october 2018
  !!
  !!----------------------------------------------------------------------!
  Use kinds, only : wp, wi
  Use kspace,  only : kspace_type
  Use errors_warnings, only : error,error_alloc,error_dealloc
  Private

  Public :: bspline_coeffs_gen, bspline_splines_gen

  !! JW952
  ! Attach to Ewald type?
  Type, Public :: bspline_type
    !> Number of required derivatives
    Integer( Kind = wi ), Public :: num_deriv
    !> SPME FFT B-spline order
    Integer( Kind = wi ), Public :: num_splines
    !> SPME FFT B-spline order when padding radius > 0
    Integer( Kind = wi ), Public :: num_spline_pad
    !> And another one
    Integer( Kind = wi ), Public :: num_spline_padded
    !> B-spline coefficients
    Complex( Kind = wp ), Dimension( :,: ),     Allocatable, Public :: coefficients
    !> Precalculated bb*
    Real( Kind = wp ),    Dimension( :,: ),     Allocatable, Public :: norm2
    !> Spline derivative
    Real( Kind = wp ),    Dimension( :,:,:,: ), Allocatable, Public :: derivs

    !> Is this bspline initialised correctly
    Logical :: derivs_initialised = .false., coeffs_initialised = .false.
  End Type bspline_type

  Real( Kind = wp ), Dimension(:,:), Save, Allocatable                                                 :: ncombk                 !! Combinations
  Real( Kind = wp ), Dimension(:), Save, Allocatable                                                   :: real_no,inv_no         !! Real variants to avoid type-casting

contains
!!! Bspline Routines

  Subroutine bspline_coeffs_gen(kspace,bspline)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine to calculate B-spline coefficients for
    !! Euler exponential splines
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith july 1998
    !!
    !!-----------------------------------------------------------------------

    Implicit None

    Type( kspace_type ),                                        Intent( In    ) :: kspace
    Type( bspline_type ),                                       Intent( InOut ) :: bspline
    Complex( Kind = wp ),    dimension(:), allocatable                          :: ww1,ww2,ww3
    Complex( Kind = wp )                                                        :: temp_spline
    Real( Kind = wp ),       dimension(:), allocatable                          :: cspline
    Integer                                                                     :: i,j,k,n
!    Integer, save                                                               :: nsplines_old = -1

    Character( Len = 256 )                                                      :: message
    Integer                                                                     :: fail

    ! Nothing to do
    if (bspline%coeffs_initialised) return

    ! Perform validity checks in routine which is only called once!
    if (bspline%num_splines < 2 .or. bspline%num_deriv < 0) then
       Write(message,'(a)') 'Error: illegal spline order in bspline_coeffs'
       Call error(0, message)
    end if
    
    if (bspline%num_deriv > bspline%num_splines - 2) then ! Can't return this many
       Write(message,'(a,/,2(a,i0))') 'Error: b-spline not high enough order for derivatives.','Num Derivatives: ', &
            & bspline%num_deriv,'Bspline order: ',bspline%num_splines
       Call error(0, message)
    End If

    allocate (ww1(1:kspace%k_vec_dim(1)),ww2(1:kspace%k_vec_dim(2)),ww3(1:kspace%k_vec_dim(3)), stat = fail)
    if (fail > 0) call error_alloc('ww arrays','bspline_coeffs_gen')

    allocate (bspline%coefficients(3,kspace%k_vec_max), stat = fail)
    if (fail > 0) call error_alloc('bspline coefficients','bspline_coeffs_gen')
    allocate (bspline%norm2(3,kspace%k_vec_max), stat = fail)
    if (fail > 0) call error_alloc('bspline norms','bspline_coeffs_gen')

    ! initialise the complex exponential arrays

    call spl_cexp(kspace%k_vec_dim(1),kspace%k_vec_dim(2),kspace%k_vec_dim(3),ww1,ww2,ww3)

    ! allocate the helper array

    allocate (cspline(1:bspline%num_splines), stat = fail)
    if (fail > 0) call error_alloc('cspline array','bspline_coeffs_gen')

    ! calculate B-splines at knots
    cspline(1)=0.0_wp
    cspline(2)=1.0_wp

    Do k=3,bspline%num_splines
      cspline(k)=0.0_wp

      Do j=k,2,-1
        cspline(j)=(Real(j-1,wp)*cspline(j)+Real(k-j+1,wp)*cspline(j-1))/Real(k-1,wp)
      End Do
    End Do

    ! calculate B-spline coefficients

    Do i=0,kspace%k_vec_dim(1)-1
      temp_spline=(0.0_wp,0.0_wp)

      Do k=0,bspline%num_splines-2
        temp_spline=temp_spline+cspline(k+2)*ww1(Mod(i*k,kspace%k_vec_dim(1))+1)
      End Do

      bspline%coefficients(1,i+1)=ww1(Mod(i*(bspline%num_splines-1),kspace%k_vec_dim(1))+1)/temp_spline
    End Do

    Do i=0,kspace%k_vec_dim(2)-1
      temp_spline=(0.0_wp,0.0_wp)

      Do k=0,bspline%num_splines-2
        temp_spline=temp_spline+cspline(k+2)*ww2(Mod(i*k,kspace%k_vec_dim(2))+1)
      End Do

      bspline%coefficients(2,i+1)=ww2(Mod(i*(bspline%num_splines-1),kspace%k_vec_dim(2))+1)/temp_spline
    End Do

    Do i=0,kspace%k_vec_dim(3)-1
      temp_spline=(0.0_wp,0.0_wp)

      Do k=0,bspline%num_splines-2
        temp_spline=temp_spline+cspline(k+2)*ww3(Mod(i*k,kspace%k_vec_dim(3))+1)
      End Do

      bspline%coefficients(3,i+1)=ww3(Mod(i*(bspline%num_splines-1),kspace%k_vec_dim(3))+1)/temp_spline
    End Do

    ! Calculate magnitude of bspline coefficients
    bspline%norm2 = Real ( bspline%coefficients*conjg(bspline%coefficients), wp)

    deallocate (cspline,         stat = fail)
    if (fail > 0) call error_dealloc('cspline arrays','bspline_coeffs_gen')
    deallocate (ww1,ww2,ww3, stat = fail)
    if (fail > 0) call error_dealloc('ww arrays','bspline_coeffs_gen')

    ! Setup constants for bspline_splines_gen

    ! if (bspline%num_splines > nsplines_old) then
       if (allocated(real_no) .and. allocated(inv_no)) then
          Deallocate (real_no,inv_no, Stat=fail)
          If (fail > 0) call error_alloc('real_no and inv_no','bspline_splines_gen')
       end if
       Allocate (real_no(1:bspline%num_splines),inv_no(1:bspline%num_splines), Stat=fail)
       If (fail > 0) call error_alloc('real_no and inv_no','bspline_splines_gen')

       Do i=1,bspline%num_splines
          real_no(i) = Real(i,wp)
          inv_no(i)  = 1.0_wp / real_no(i)
       End Do

       if (allocated(ncombk)) then
          Deallocate(ncombk, stat = fail)
          if (fail > 0) call error_dealloc('ncombk','bspline_coeffs_gen')
       end if
           allocate(ncombk(bspline%num_splines,0:bspline%num_splines), stat = fail)
       if (fail > 0) call error_alloc('ncombk','bspline_coeffs_gen')

       ncombk = 0.0_wp
       Do n = 1, bspline%num_splines ! If we change spline number, need to recompute coeffs anyway
         Do k = 0, n
           ncombk(n,k) = product([(real_no(i),i=n-k+1,n)])*product([(inv_no(i),i=1,max(1,k))])
         End Do
       End Do
       bspline%derivs_initialised = .false.
    ! end if

    ! nsplines_old = max(nsplines_old,bspline%num_splines)
    
  End Subroutine bspline_coeffs_gen

  Subroutine bspline_splines_gen(num_atoms,recip_coords,bspline)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine to calculate B-splines for SPME method
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith july 1998
    !! amended   - i.t.todorov april 2015
    !! reordered - j.s.wilkins november 2018
    !!-----------------------------------------------------------------------

    Implicit None

    Integer,                                                                       Intent( In    ) :: num_atoms              !! Number of atoms to fill
    Real( Kind = wp ), Dimension(3,num_atoms),                                     Intent( In    ) :: recip_coords           !! Coordinates of charge centres in reciprocal cell
    Type( bspline_type ),                                                          Intent( InOut ) :: bspline                !! Bspline to be created
    Real( Kind = wp ), Dimension(3)                                                                :: current_bspline_centre !! Origin of current bspline being calculated
    Real( Kind = wp ), Dimension(3)                                                                :: current_bspline_point  !! Current point of the bspline being calculated
    Real( Kind = wp )                                                                              :: sgn                    !! Sign alternation
    Real( Kind = wp )                                                                              :: jm1_r,k_r,km1_rr       !! Real variants to avoid type-casting
    Integer                                                                                        :: i,j,k,l,n              !! Loop counters
    Integer                                                                                        :: s
    
    ! Perform validity checks in routine which is only called once -- bspline_coeffs_gen !

    s = bspline%num_splines
    ! construct B-splines

    !! JW952
    ! Reversed order of array for memory access efficiency
    do i = 1, num_atoms

      ! initializing 2nd order B-spline
      ! for u where (0<u<1) and (1<u<2)

      bspline%derivs(:,0,s,i) = recip_coords(:,i)-Aint(recip_coords(:,i),wp)
      bspline%derivs(:,0,s-1,i) = 1.0_wp - bspline%derivs(:,0,s,i)

      ! Now on to calculate order k B-spline values at k
      ! points where (0<u<k)

      current_bspline_centre(:)=bspline%derivs(:,0,s,i)

      Do k=s-2,bspline%num_deriv+1,-1 ! 3,bspline%num_splines-bspline%num_deriv ! Order of B-spline

        bspline%derivs(:,0,k,i) = 0.0_wp
        k_r   =real_no(s-k+1)
        km1_rr=inv_no(s-k)

        Do j=k,s-1 ! Compute order k B-spline at points {k,k-1,...,1}

          jm1_r=real_no(s-j)
          current_bspline_point(:)=current_bspline_centre(:)+jm1_r

          !! Possibly wrong
          bspline%derivs(:,0,j,i)=(current_bspline_point(:)*bspline%derivs(:,0,j,i) &
            & + (k_r-current_bspline_point(:))*bspline%derivs(:,0,j+1,i))*km1_rr

        End Do

        bspline%derivs(:,0,s,i)=bspline%derivs(:,0,s,i)*current_bspline_centre(:)*km1_rr
      End Do

      ! Now compute B-splines for order bspline%num_splines at k points where
      ! (0<u<bspline%num_splines)

      !! JW952
      ! New loop here!
      Do l = bspline%num_deriv, 1, -1

        k= l
        bspline%derivs(:,0,k,i) = 0.0_wp
        k_r   = real_no(s-k+1)
        km1_rr= inv_no(s-k)

        Do j=1,s-1!bspline%num_splines,2,-1

          bspline%derivs(:,l,j,i) = 0.0_wp
          ! Derivatives of B-splines with order nospl at k-1 points
          sgn = 1.0_wp
          do n = 0,Min(l,s-j)
            bspline%derivs(:,l,j,i)=bspline%derivs(:,l,j,i) + sgn*ncombk(l,n)*bspline%derivs(:,0,j+n,i)
            sgn = -sgn
          end do

          ! Generate current point at a lag behind derivs
          jm1_r=real_no(s-j)
          current_bspline_point(:)=current_bspline_centre(:)+jm1_r
          bspline%derivs(:,0,j,i)=(current_bspline_point(:)*bspline%derivs(:,0,j,i) &
            & + (k_r-current_bspline_point(:))*bspline%derivs(:,0,j+1,i))*km1_rr

        End Do

        bspline%derivs(:,l,s,i)=bspline%derivs(:,0,s,i)

        bspline%derivs(:,0,s,i)=bspline%derivs(:,0,s,i)*current_bspline_centre(:)*km1_rr

      end do

    end do
    
  end Subroutine bspline_splines_gen

  Subroutine spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 routine to create complex exponential arrays for b-splines
    !!
    !! copyright - daresbury laboratory
    !! author    - w.smith october 1998
    !! amended   - i.t.todorov october 2006
    !!
    !!-----------------------------------------------------------------------
    Use constants, only : twopi
    implicit none

    Integer,              Intent( In    ) :: ndiv1,ndiv2,ndiv3
    Complex( Kind = wp ), Intent(   Out ) :: ww1(1:ndiv1),ww2(1:ndiv2),ww3(1:ndiv3)

    Integer           :: i
    Real( Kind = wp ) :: arg

    ! initialise complex exponential factors

    ww1(1)=(1.0_wp,0.0_wp)

    Do i=1,ndiv1/2
      arg=(twopi/Real(ndiv1,wp))*Real(i,wp)
      ww1(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
      ww1(ndiv1+1-i)=Conjg(ww1(i+1))

    End Do
    ww2(1)=(1.0_wp,0.0_wp)

    Do i=1,ndiv2/2
      arg=(twopi/Real(ndiv2,wp))*Real(i,wp)

      ww2(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
      ww2(ndiv2+1-i)=Conjg(ww2(i+1))
    End Do

    ww3(1)=(1.0_wp,0.0_wp)

    Do i=1,ndiv3/2
      arg=(twopi/Real(ndiv3,wp))*Real(i,wp)

      ww3(i+1)=Cmplx(Cos(arg),Sin(arg), Kind = wp)
      ww3(ndiv3+1-i)=Conjg(ww3(i+1))
    End Do

  End Subroutine spl_cexp


end module bspline
