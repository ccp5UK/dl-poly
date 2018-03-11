Module defects1_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global defects variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp

  Implicit None

! REFERENCE1 data

  Logical,                             Save :: l_dfx = .false.

  Integer,                             Save :: nrefs1  = 0 , &
                                               nlrefs1 = 0

  Real( Kind = wp ),                   Save :: celr1(1:9) = 0.0_wp


  Character( Len = 8 ),   Allocatable, Save :: namr1(:)

  Integer,                Allocatable, Save :: lri1(:),lra1(:),indr1(:)

  Real( Kind = wp ),      Allocatable, Save :: xr1(:),yr1(:),zr1(:)

  Public :: allocate_defects1_arrays

Contains

  Subroutine allocate_defects1_arrays()

    Use setup_module, Only : mxatms

    Implicit None

    Integer, Dimension( 1:3 ) :: fail

    fail = 0

    Allocate (namr1(1:mxatms),                               Stat = fail(1))
    Allocate (lri1(1:mxatms),lra1(1:mxatms),indr1(1:mxatms), Stat = fail(2))
    Allocate (xr1(1:mxatms),yr1(1:mxatms),zr1(1:mxatms),     Stat = fail(3))

    If (Any(fail > 0)) Call error(1035)

    namr1 = ' '
    lri1 = 0 ; lra1 = 0 ; indr1 = 0
    xr1 = 0.0_wp ; yr1 = 0.0_wp ; zr1 = 0.0_wp

  End Subroutine allocate_defects1_arrays

End Module defects1_module
