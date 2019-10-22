Module mpole
  Use kinds, Only : wp,wi
  Use constants, Only : zero_plus,nmpldt
  Use site, Only : site_type
  Use parse,      Only : word_2_real, get_line, get_word, lower_case, strip_blanks 
  Use comms, Only : comms_type
  Use numerics, Only : factorial
  Use errors_warnings, Only : error,warning,info
  Use numerics, Only : factorial
  Implicit None

  Private

  !> A derived data type for defining the rotation matrix for multipolar sites
  Type, Public :: rot_mat
    Real( Kind = wp ), Dimension(1:9) :: mtrxa
    Real( Kind = wp ), Dimension(1:3) :: p1,p2
    Integer                           :: flag,mbnd(1:2)
  End Type rot_mat

  ! Polarisation keys
  !> Default polarisation, unscreened & undamped - iAMOEBA like
  Integer( Kind = wi ), Parameter, Public :: POLARISATION_DEFAULT = 0
  !> CHARMM polarisation:
  !>
  !> - $q_\text{shell} = -\lvert q_core \rvert * \sqrt(\alpha * k)
  !> - $k_\text{CHARMM} = 1000 \, \mathrm{kcal} \mathrm{mol}^−1 \mathrm{Å}^−2
  Integer( Kind = wi ), Parameter, Public :: POLARISATION_CHARMM = 1


  Integer,           Dimension(3,35),    Parameter       :: mpole_derivs  = reshape([& !! Mpole derivs in mpole_module order
    & [0,0,0], &
    & [1,0,0], [0,1,0], [0,0,1], &
    & [2,0,0], [1,1,0], [1,0,1], [0,2,0], [0,1,1], [0,0,2], &
    & [3,0,0], [2,1,0], [2,0,1], [1,2,0], [1,1,1], [1,0,2], [0,3,0], [0,2,1], [0,1,2], [0,0,3], &
    & [4,0,0], [3,1,0], [3,0,1], [2,2,0], [2,1,1], [2,0,2], [1,3,0], [1,2,1], [1,1,2], [1,0,3], &
    & [0,4,0], [0,3,1], [0,2,2], [0,1,3], [0,0,4]&
    & ], [3, 35])
  Integer,           Dimension(0:4),     Parameter       :: nmpole_derivs = [1,3,6,10,15]
  
  !> Type containing multipole data
  Type, Public :: mpole_type
    Private

    !> Type of inducible (self-polarisation) scheme
    Integer( Kind = wi ), Public :: key = POLARISATION_DEFAULT

    !> Default thole dumping for CHARMM representation
    Real( Kind = wp ), Public :: thole  = 1.3_wp

    ! Variables for multipolar interactions
    !> Mappings from three indices multipole to a one index multipole
    Integer( Kind = wi ), Allocatable, Public :: map(:,:,:)
    Integer( Kind = wi ), Allocatable, Public :: ltg(:)
    !> Rotation counter flag
    Integer( Kind = wi ), Allocatable, Public :: flg(:)
    !> Bonded connectivity
    Integer( Kind = wi ), Allocatable, Public :: ltp(:,:)
    !> CHARMM core-shell screened electrostatics induction list
    Integer( Kind = wi ), Allocatable, Public :: charmm(:,:)

    !> Local frame
    Real( Kind = wp ), Allocatable, Public :: local_frame(:,:)
    !> Global frame
    Real( Kind = wp ), Allocatable, Public :: global_frame(:,:)

    !> Induced dipole polarisation for sites (inversed if non-zero)
    Real( Kind = wp ), Allocatable, Public :: polarisation_site(:)
    !> Induced dipole polarisation for atoms (inversed if non-zero)
    Real( Kind = wp ), Allocatable, Public :: polarisation_atom(:)

    !> thole dumping coefficient/factor (for self-polarisation) for sites
    Real( Kind = wp ), Allocatable, Public :: dump_site(:)
    !> thole dumping coefficient/factor (for self-polarisation) for atoms
    Real( Kind = wp ), Allocatable, Public :: dump_atom(:)

    !> Rotation matrices
    Type( rot_mat ),   Allocatable, Public :: rotation(:)

    !> Infinitesimal rotations about x
    Real( Kind = wp ), Allocatable, Public :: rotation_x(:,:)
    !> Infinitesimal rotations about y
    Real( Kind = wp ), Allocatable, Public :: rotation_y(:,:)
    !> Infinitesimal rotations about z
    Real( Kind = wp ), Allocatable, Public :: rotation_z(:,:)

    !> Torques due to infinitesimal rotations about x
    Real( Kind = wp ), Allocatable, Public :: torque_x(:)
    !> Torques due to infinitesimal rotations about y
    Real( Kind = wp ), Allocatable, Public :: torque_y(:)
    !> Torques due to infinitesimal rotations about z
    Real( Kind = wp ), Allocatable, Public :: torque_z(:)

    !> n choose k values, used in computing the reciprocal space Ewald sum field
    Real ( Kind = wp ), Allocatable, Public :: n_choose_k(:,:)

    !> How many mpole derivatives
    Integer, Dimension(0:0), Public :: nmpole_derivs = [1]
    !> My mpole derivatives
    Integer, Dimension(3,1,0:0), Public :: mpole_derivs = reshape([0,0,0],[3,1,1])

    Integer, Public :: num_mpoles = 0
    
    !> Maximum number of multipoles
    Integer( Kind = wi ), Public :: max_mpoles
    !> Maximum multipolar order
    Integer( Kind = wi ), Public :: max_order

  Contains
    Private

    Procedure, Public :: init => allocate_mpoles_arrays
    Final :: cleanup
  End Type mpole_type

Contains

  Subroutine allocate_mpoles_arrays(T,max_site,max_exclude,mxatdm,bspline,mxatms)
  Class( mpole_type ) :: T
    Integer( Kind = wi ), Intent( In    ) :: max_site
    Integer( Kind = wi ), Intent( In    ) :: max_exclude
    Integer( Kind = wi ), Intent( In    ) :: mxatdm
    Integer( Kind = wi ), Intent( In    ) :: bspline
    Integer( Kind = wi ), Intent( In    ) :: mxatms

    Integer           :: n,k,om1,numpl,fail(1:9)

    ! No MPOLES file read <= no multipoles directive in FIELD
    If (T%max_mpoles < 1) Return

    om1 = T%max_order + 1
    numpl = (3**om1 - 1)/2

    fail = 0

    Allocate (T%map(0:T%max_order,0:T%max_order,0:T%max_order),T%ltg(1:numpl), stat=fail(1))
    Allocate (T%flg(1:mxatdm),T%ltp(0:max_exclude,1:mxatdm), stat=fail(2))
    If (T%key == 1) Then
      Allocate (T%charmm(0:max_exclude,1:mxatdm), stat=fail(3))
    End If
    Allocate (T%local_frame(1:T%max_mpoles,1:max_site),T%global_frame(1:T%max_mpoles,1:mxatms), stat=fail(4))
    Allocate (T%polarisation_site(1:max_site),T%polarisation_atom(1:mxatms), stat=fail(5))
    Allocate (T%dump_site(1:max_site),T%dump_atom(1:mxatms), stat=fail(6))
    Allocate (T%rotation(1:mxatdm),T%n_choose_k(0:bspline,0:bspline), stat=fail(7))
    Allocate (T%torque_x(1:mxatdm),T%torque_y(1:mxatdm),T%torque_z(1:mxatdm), stat=fail(8))
    Allocate (T%rotation_x(1:T%max_mpoles,1:mxatms), &
      T%rotation_y(1:T%max_mpoles,1:mxatms), &
      T%rotation_z(1:T%max_mpoles,1:mxatms), stat=fail(9))

    If (Any(fail > 0)) Call error(1025)

    T%flg = 0 ; T%ltp = 0
    If (T%key == 1) Then
      T%charmm = 0
    End If

    T%local_frame = 0.0_wp ; T%global_frame = 0.0_wp
    T%polarisation_site = 0.0_wp ; T%polarisation_atom = 0.0_wp
    T%dump_site = 0.0_wp ; T%dump_atom = 0.0_wp

    Do n=1,mxatdm
      T%rotation(n)%mtrxa = 0.0_wp
      T%rotation(n)%p1    = 0.0_wp
      T%rotation(n)%p2    = 0.0_wp
      T%rotation(n)%flag  = 0
      T%rotation(n)%mbnd  = 0
    End Do

    T%torque_x = 0.0_wp ; T%torque_y = 0.0_wp ; T%torque_z = 0.0_wp

    T%rotation_x = 0.0_wp ; T%rotation_y = 0.0_wp ; T%rotation_z = 0.0_wp

    ! Build the multipole map (polymap) and compute the constants T%n_choose_k
    ! Also build the map (T%ltg) that converts between index of a local
    ! multipole to the index of the corresponding global multipole
    T%map(0,0,0)=1

    T%ltg(1)=1

    If (T%max_order >= 1) Then
      T%map(1,0,0)=2 ; T%map(0,1,0)=3 ; T%map(0,0,1)=4

      T%ltg(2)=2 ; T%ltg(3)=3 ; T%ltg(4)=4
    End If

    If (T%max_order >= 2) Then
      T%map(2,0,0)=5 ; T%map(1,1,0)=6 ; T%map(1,0,1)=7
      T%map(0,2,0)=8 ; T%map(0,1,1)=9 ; T%map(0,0,2)=10

      T%ltg(5) =5 ; T%ltg(6) =6 ; T%ltg(7) =7
      T%ltg(8) =6 ; T%ltg(9) =8 ; T%ltg(10)=9
      T%ltg(11)=7 ; T%ltg(12)=9 ; T%ltg(13)=10
    End If

    If (T%max_order >= 3) Then
      T%map(3,0,0)=11 ; T%map(2,1,0)=12 ; T%map(2,0,1)=13
      T%map(1,2,0)=14 ; T%map(1,1,1)=15 ; T%map(1,0,2)=16
      T%map(0,3,0)=17 ; T%map(0,2,1)=18 ; T%map(0,1,2)=19
      T%map(0,0,3)=20

      T%ltg(14)=11 ; T%ltg(15)=12 ; T%ltg(16)=13
      T%ltg(17)=12 ; T%ltg(18)=14 ; T%ltg(19)=15
      T%ltg(20)=13 ; T%ltg(21)=15 ; T%ltg(22)=16
      T%ltg(23)=12 ; T%ltg(24)=14 ; T%ltg(25)=15
      T%ltg(26)=14 ; T%ltg(27)=17 ; T%ltg(28)=18
      T%ltg(29)=15 ; T%ltg(30)=18 ; T%ltg(31)=19
      T%ltg(32)=13 ; T%ltg(33)=15 ; T%ltg(34)=16
      T%ltg(35)=15 ; T%ltg(36)=18 ; T%ltg(37)=19
      T%ltg(38)=16 ; T%ltg(39)=19 ; T%ltg(40)=20
    End If

    If (T%max_order >=4) Then
      T%map(4,0,0)=21 ; T%map(3,1,0)=22 ; T%map(3,0,1)=23
      T%map(2,2,0)=24 ; T%map(2,1,1)=25 ; T%map(2,0,2)=26
      T%map(1,3,0)=27 ; T%map(1,2,1)=28 ; T%map(1,1,2)=29
      T%map(1,0,3)=30 ; T%map(0,4,0)=31 ; T%map(0,3,1)=32
      T%map(0,2,2)=33 ; T%map(0,1,3)=34 ; T%map(0,0,4)=35

      T%ltg(41)=21 ; T%ltg(42)=22 ; T%ltg(43)=23
      T%ltg(44)=22 ; T%ltg(45)=24 ; T%ltg(46)=25
      T%ltg(47)=23 ; T%ltg(48)=25 ; T%ltg(49)=26
      T%ltg(50)=22 ; T%ltg(51)=24 ; T%ltg(52)=25
      T%ltg(53)=24 ; T%ltg(54)=27 ; T%ltg(55)=28
      T%ltg(56)=25 ; T%ltg(57)=28 ; T%ltg(58)=29
      T%ltg(59)=23 ; T%ltg(60)=25 ; T%ltg(61)=26
      T%ltg(62)=25 ; T%ltg(63)=28 ; T%ltg(64)=29
      T%ltg(65)=26 ; T%ltg(66)=29 ; T%ltg(67)=30
      T%ltg(68)=22 ; T%ltg(69)=24 ; T%ltg(70)=25
      T%ltg(71)=24 ; T%ltg(72)=27 ; T%ltg(73)=28
      T%ltg(74)=25 ; T%ltg(75)=28 ; T%ltg(76)=29
      T%ltg(77)=24 ; T%ltg(78)=27 ; T%ltg(79)=28
      T%ltg(80)=27 ; T%ltg(81)=31 ; T%ltg(82)=32
      T%ltg(83)=28 ; T%ltg(84)=32 ; T%ltg(85)=33
      T%ltg(86)=25 ; T%ltg(87)=28 ; T%ltg(88)=29
      T%ltg(89)=28 ; T%ltg(90)=32 ; T%ltg(91)=33
      T%ltg(92)=29 ; T%ltg(93)=33 ; T%ltg(94)=34
      T%ltg(95)=23 ; T%ltg(96)=25 ; T%ltg(97)=26
      T%ltg(98)=25 ; T%ltg(99)=28 ; T%ltg(100)=29

      T%ltg(101)=26 ; T%ltg(102)=29 ; T%ltg(103)=30
      T%ltg(104)=25 ; T%ltg(105)=28 ; T%ltg(106)=29
      T%ltg(107)=28 ; T%ltg(108)=32 ; T%ltg(109)=33
      T%ltg(110)=29 ; T%ltg(111)=33 ; T%ltg(112)=34
      T%ltg(113)=26 ; T%ltg(114)=29 ; T%ltg(115)=30
      T%ltg(116)=29 ; T%ltg(117)=33 ; T%ltg(118)=34
      T%ltg(119)=30 ; T%ltg(120)=34 ; T%ltg(121)=35
    End If

    If (T%max_order < 0 .or. T%max_order >=5) Call error(2071)

    ! compute n choose k
    Do k = 0, bspline
      Do n = 0, bspline
        T%n_choose_k(n,k) = Exp(Factorial(n)-Factorial(n-k)-Factorial(k))
      End Do
    End Do
  End Subroutine allocate_mpoles_arrays

  Subroutine cleanup(T)
    Type( mpole_type ) :: T

    If (Allocated(T%map)) Then
      Deallocate(T%map)
    End If
    If (Allocated(T%ltg)) Then
      Deallocate(T%ltg)
    End If
    If (Allocated(T%flg)) Then
      Deallocate(T%flg)
    End If
    If (Allocated(T%ltp)) Then
      Deallocate(T%ltp)
    End If
    If (Allocated(T%charmm)) Then
      Deallocate(T%charmm)
    End If

    If (Allocated(T%local_frame)) Then
      Deallocate(T%local_frame)
    End If
    If (Allocated(T%global_frame)) Then
      Deallocate(T%global_frame)
    End If

    If (Allocated(T%polarisation_site)) Then
      Deallocate(T%polarisation_site)
    End If
    If (Allocated(T%polarisation_atom)) Then
      Deallocate(T%polarisation_atom)
    End If

    If (Allocated(T%dump_site)) Then
      Deallocate(T%dump_site)
    End If
    If (Allocated(T%dump_atom)) Then
      Deallocate(T%dump_atom)
    End If

    If (Allocated(T%rotation)) Then
      Deallocate(T%rotation)
    End If

    If (Allocated(T%rotation_x)) Then
      Deallocate(T%rotation_x)
    End If
    If (Allocated(T%rotation_y)) Then
      Deallocate(T%rotation_y)
    End If
    If (Allocated(T%rotation_z)) Then
      Deallocate(T%rotation_z)
    End If

    If (Allocated(T%torque_x)) Then
      Deallocate(T%torque_x)
    End If
    If (Allocated(T%torque_y)) Then
      Deallocate(T%torque_y)
    End If
    If (Allocated(T%torque_z)) Then
      Deallocate(T%torque_z)
    End If

    If (Allocated(T%n_choose_k)) Then
      Deallocate(T%n_choose_k)
    End If
  End Subroutine cleanup

End Module mpole
