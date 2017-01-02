Module mpoles_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring configuration variables and arrays for
! multipoles
!
! copyright - daresbury laboratory
! author    - h.a.boateng & i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

! A derived data type for defining the rotation matrix for multipolar sites

  Type rot_mat
     Real( Kind = wp ), Dimension(1:9) :: mtrxa
     Real( Kind = wp ), Dimension(1:3) :: p1,p2
     Integer                           :: flag,mbnd(1:2)
  End Type rot_mat

! Type of inducible (self-polarisation) scheme

  Integer,           Save :: keyind = 0 ! 0 - default :: unscreened & undamped - iAMOEBA like
                                        ! 1 - CHARMM  :: q_shell == -Sign(q_core) * Sqrt(alpha * k) ; k_CHARMM = 1000 kcal*mol^−1*Å^−2

  Real,              Save :: thole  = 1.3_wp  ! default thole dumping for CHARMM representation

! variables for multipolar interactions

  Logical,           Save :: induce=.false. , &
                             gear,aspc,lstsq
  Integer,           Save :: numcof,politer
  Real( Kind = wp ), Save :: convcrit,enepol

  Integer,           Allocatable, Save :: mplmap(:,:,:),mplltg(:)             ! mappings from three indices multipole to a one index multipole
  Integer,           Allocatable, Save :: mplflg(:)                           ! rotation counter flag
  Integer,           Allocatable, Save :: ltpatm(:,:)                         ! bonded connectivity
  Integer,           Allocatable, Save :: lchatm(:,:)                         ! CHARMM core-shell screened electrostatics induction list

  Real( Kind = wp ), Allocatable, Save :: mpllfr(:,:),mplgfr(:,:)             ! local/lab(site) and global(atom) frames
  Real( Kind = wp ), Allocatable, Save :: plrsit(:),plratm(:)                 ! induced dipole polarisation for sites and atoms
  Real( Kind = wp ), Allocatable, Save :: dmpsit(:),dmpatm(:)                 ! sites' and atoms' (thole) dumping coefficient/factor (for self-polarisation)

  Type( rot_mat ),   Allocatable, Save :: mprotm(:)                           ! rotation matrices

  Real( Kind = wp ), Allocatable, Save :: mprotx(:,:),mproty(:,:),mprotz(:,:) ! infinitesimal rotations
  Real( Kind = wp ), Allocatable, Save :: mptrqx(:),mptrqy(:),mptrqz(:)       ! torques due to infinitesimal rotations

  Real( Kind = wp ), Allocatable, Save :: ncombk(:,:)                         ! n combination k values for usage in computing the reciprocal space Ewald sum

  Real( Kind = wp ), Allocatable, Save :: mpfldx(:),mpfldy(:),mpfldz(:)       ! field

  Real( Kind = wp ), Allocatable, Save :: muindx(:),muindy(:),muindz(:)
  Real( Kind = wp ), Allocatable, Save :: indipx(:),indipy(:),indipz(:)

  Real( Kind = wp ), Allocatable, Save :: upidpx(:,:),upidpy(:,:),upidpz(:,:)
  Real( Kind = wp ), Allocatable, Save :: rsdx(:),rsdy(:),rsdz(:)
  Real( Kind = wp ), Allocatable, Save :: polcof(:)


  Public :: allocate_mpoles_arrays

Contains

  Subroutine allocate_mpoles_arrays()

    Use setup_module, Only : mxsite,mxexcl,mxspl,mxompl,mximpl,mxatdm,mxatms

    Implicit None

    Integer           :: n,k,om1,numpl,fail(1:9)
    Real( Kind = wp ) :: Factorial,gearp(1:7),aspcp(1:7)

    If (mximpl < 1) Return ! no MPOLES file read <= no multipoles directive in FIELD

    om1 = mxompl + 1
    numpl = (3**om1 - 1)/2

    fail = 0

    Allocate (mplmap(0:mxompl,0:mxompl,0:mxompl),mplltg(1:numpl),  Stat = fail(1))
    Allocate (mplflg(1:mxatdm),ltpatm(0:mxexcl,1:mxatdm),          Stat = fail(2))
    If (keyind == 1) &
    Allocate (lchatm(0:mxexcl,1:mxatdm),                           Stat = fail(3))
    Allocate (mpllfr(1:mximpl,1:mxsite),mplgfr(1:mximpl,1:mxatms), Stat = fail(4))
    Allocate (plrsit(1:mxsite),plratm(1:mxatms),                   Stat = fail(5))
    Allocate (dmpsit(1:mxsite),dmpatm(1:mxatms),                   Stat = fail(6))
    Allocate (mprotm(1:mxatdm),ncombk(0:mxspl,0:mxspl),            Stat = fail(7))
    Allocate (mptrqx(1:mxatdm),mptrqy(1:mxatdm),mptrqz(1:mxatdm),  Stat = fail(8))
    Allocate (mprotx(1:mximpl,1:mxatms), &
              mproty(1:mximpl,1:mxatms), &
              mprotz(1:mximpl,1:mxatms),                           Stat = fail(9))

    If (Any(fail > 0)) Call error(1025)

    mplflg = 0 ; ltpatm = 0
    If (keyind == 1) lchatm = 0

    mpllfr = 0.0_wp ; mplgfr = 0.0_wp
    plrsit = 0.0_wp ; plratm = 0.0_wp
    dmpsit = 0.0_wp ; dmpatm = 0.0_wp

    Do n=1,mxatdm
       mprotm(n)%mtrxa = 0.0_wp
       mprotm(n)%p1    = 0.0_wp
       mprotm(n)%p2    = 0.0_wp
       mprotm(n)%flag  = 0
       mprotm(n)%mbnd  = 0
    End Do

    mptrqx = 0.0_wp ; mptrqy = 0.0_wp ; mptrqz = 0.0_wp

    mprotx = 0.0_wp ; mproty = 0.0_wp ; mprotz = 0.0_wp

! Build the multipole map (polymap) and compute the constants ncombk
! Also build the map (mplltg) that converts between index of a local
! multipole to the index of the corresponding global multipole

    mplmap(0,0,0)=1

    mplltg(1)=1

    If (mxompl >= 1) Then

       mplmap(1,0,0)=2 ; mplmap(0,1,0)=3 ; mplmap(0,0,1)=4

       mplltg(2)=2 ; mplltg(3)=3 ; mplltg(4)=4

    End If

    If (mxompl >= 2) Then

       mplmap(2,0,0)=5 ; mplmap(1,1,0)=6 ; mplmap(1,0,1)=7
       mplmap(0,2,0)=8 ; mplmap(0,1,1)=9 ; mplmap(0,0,2)=10

       mplltg(5) =5 ; mplltg(6) =6 ; mplltg(7) =7
       mplltg(8) =6 ; mplltg(9) =8 ; mplltg(10)=9
       mplltg(11)=7 ; mplltg(12)=9 ; mplltg(13)=10

    End If

    If (mxompl >= 3) Then

       mplmap(3,0,0)=11 ; mplmap(2,1,0)=12 ; mplmap(2,0,1)=13
       mplmap(1,2,0)=14 ; mplmap(1,1,1)=15 ; mplmap(1,0,2)=16
       mplmap(0,3,0)=17 ; mplmap(0,2,1)=18 ; mplmap(0,1,2)=19
       mplmap(0,0,3)=20

       mplltg(14)=11 ; mplltg(15)=12 ; mplltg(16)=13
       mplltg(17)=12 ; mplltg(18)=14 ; mplltg(19)=15
       mplltg(20)=13 ; mplltg(21)=15 ; mplltg(22)=16
       mplltg(23)=12 ; mplltg(24)=14 ; mplltg(25)=15
       mplltg(26)=14 ; mplltg(27)=17 ; mplltg(28)=18
       mplltg(29)=15 ; mplltg(30)=18 ; mplltg(31)=19
       mplltg(32)=13 ; mplltg(33)=15 ; mplltg(34)=16
       mplltg(35)=15 ; mplltg(36)=18 ; mplltg(37)=19
       mplltg(38)=16 ; mplltg(39)=19 ; mplltg(40)=20

    End If

    If (mxompl >=4) Then

       mplmap(4,0,0)=21 ; mplmap(3,1,0)=22 ; mplmap(3,0,1)=23
       mplmap(2,2,0)=24 ; mplmap(2,1,1)=25 ; mplmap(2,0,2)=26
       mplmap(1,3,0)=27 ; mplmap(1,2,1)=28 ; mplmap(1,1,2)=29
       mplmap(1,0,3)=30 ; mplmap(0,4,0)=31 ; mplmap(0,3,1)=32
       mplmap(0,2,2)=33 ; mplmap(0,1,3)=34 ; mplmap(0,0,4)=35

       mplltg(41)=21 ; mplltg(42)=22 ; mplltg(43)=23
       mplltg(44)=22 ; mplltg(45)=24 ; mplltg(46)=25
       mplltg(47)=23 ; mplltg(48)=25 ; mplltg(49)=26
       mplltg(50)=22 ; mplltg(51)=24 ; mplltg(52)=25
       mplltg(53)=24 ; mplltg(54)=27 ; mplltg(55)=28
       mplltg(56)=25 ; mplltg(57)=28 ; mplltg(58)=29
       mplltg(59)=23 ; mplltg(60)=25 ; mplltg(61)=26
       mplltg(62)=25 ; mplltg(63)=28 ; mplltg(64)=29
       mplltg(65)=26 ; mplltg(66)=29 ; mplltg(67)=30
       mplltg(68)=22 ; mplltg(69)=24 ; mplltg(70)=25
       mplltg(71)=24 ; mplltg(72)=27 ; mplltg(73)=28
       mplltg(74)=25 ; mplltg(75)=28 ; mplltg(76)=29
       mplltg(77)=24 ; mplltg(78)=27 ; mplltg(79)=28
       mplltg(80)=27 ; mplltg(81)=31 ; mplltg(82)=32
       mplltg(83)=28 ; mplltg(84)=32 ; mplltg(85)=33
       mplltg(86)=25 ; mplltg(87)=28 ; mplltg(88)=29
       mplltg(89)=28 ; mplltg(90)=32 ; mplltg(91)=33
       mplltg(92)=29 ; mplltg(93)=33 ; mplltg(94)=34
       mplltg(95)=23 ; mplltg(96)=25 ; mplltg(97)=26
       mplltg(98)=25 ; mplltg(99)=28 ; mplltg(100)=29

       mplltg(101)=26 ; mplltg(102)=29 ; mplltg(103)=30
       mplltg(104)=25 ; mplltg(105)=28 ; mplltg(106)=29
       mplltg(107)=28 ; mplltg(108)=32 ; mplltg(109)=33
       mplltg(110)=29 ; mplltg(111)=33 ; mplltg(112)=34
       mplltg(113)=26 ; mplltg(114)=29 ; mplltg(115)=30
       mplltg(116)=29 ; mplltg(117)=33 ; mplltg(118)=34
       mplltg(119)=30 ; mplltg(120)=34 ; mplltg(121)=35

    End If

    If (mxompl < 0 .or. mxompl >=5) Call error(2071)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute n choose k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Do k = 0, mxspl
       Do n = 0, mxspl
          ncombk(n,k) = Exp(Factorial(n)-Factorial(n-k)-Factorial(k))
       End Do
    End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    If (induce) Then

       Allocate (mpfldx(1:mxatms),mpfldy(1:mxatms),mpfldz(1:mxatms),            Stat = fail(1))
       Allocate (muindx(1:mxatms),muindy(1:mxatms),muindz(1:mxatms),            Stat = fail(2))
       Allocate (indipx(1:mxatms),indipy(1:mxatms),indipz(1:mxatms),            Stat = fail(3))
       Allocate (polcof(numcof),                                                Stat = fail(4))
       Allocate (upidpx(1:numcof,1:mxatms),upidpy(1:numcof,1:mxatms),upidpz(1:numcof,1:mxatms),&
                                                                                Stat = fail(5))
       Allocate (rsdx(1:mxatms),rsdy(1:mxatms),rsdz(1:mxatms),                  Stat = fail(6))

       If (Any(fail > 0)) Call error(1025)

       mpfldx = 0.0_wp; mpfldy = 0.0_wp; mpfldz = 0.0_wp

       muindx = 0.0_wp; muindy = 0.0_wp; muindz = 0.0_wp

       indipx = 0.0_wp; indipy = 0.0_wp; indipz = 0.0_wp

       upidpx = 0.0_wp; upidpy = 0.0_wp; upidpz = 0.0_wp

       rsdx   = 0.0_wp; rsdy   = 0.0_wp; rsdz   = 0.0_wp  ! arrays to store
                                                          ! residuals for
                                                          ! conjugate gradient
                                                          ! minimization of
                                                          ! polarization energy

! Coefficients for polynomial predictor

       polcof = 0.0_wp

! Gear predictor-corrector

       gearp = 0.0_wp

       gearp(1) =   6.0_wp ; gearp(2) = -15.0_wp ; gearp(3) =  20.0_wp
       gearp(4) = -15.0_wp ; gearp(5) =   6.0_wp ; gearp(6) =  -1.0_wp


! Always stable predictor-corrector (aspc)

       aspcp = 0.0_wp

       aspcp(1) =  22.0_wp/ 7.0_wp ; aspcp(2) = -55.0_wp/14.0_wp; aspcp(3) =  55.0_wp/21.0_wp
       aspcp(4) = -22.0_wp/21.0_wp ; aspcp(5) =   5.0_wp/21.0_wp; aspcp(6) =  -1.0_wp/42.0_wp

       If (gear) Then

          Do n = 1, numcof
             polcof(n) = gearp(n)
          End Do

       Else If (aspc) Then

          Do n = 1, numcof
             polcof(n) = aspcp(n)
          End Do

       End If

    End If

  End Subroutine allocate_mpoles_arrays

End Module mpoles_module
