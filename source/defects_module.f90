Module defects_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global defects variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

! REFERENCE data

  Integer,                             Save :: nrefs  = 0 , &
                                               nlrefs = 0

  Real( Kind = wp ),                   Save :: celr(1:9) = 0.0_wp


  Character( Len = 8 ),   Allocatable, Save :: namr(:)

  Integer,                Allocatable, Save :: lri(:),lra(:),indr(:)

  Real( Kind = wp ),      Allocatable, Save :: xr(:),yr(:),zr(:)

  Public :: allocate_defects_arrays

Contains

  Subroutine allocate_defects_arrays()

    Use setup_module, Only : mxatms

    Implicit None

    Integer, Dimension( 1:3 ) :: fail

    fail = 0

    Allocate (namr(1:mxatms),                             Stat = fail(1))
    Allocate (lri(1:mxatms),lra(1:mxatms),indr(1:mxatms), Stat = fail(2))
    Allocate (xr(1:mxatms),yr(1:mxatms),zr(1:mxatms),     Stat = fail(3))

    If (Any(fail > 0)) Call error(1035)

    namr = ' '
    lri = 0 ; lra = 0 ; indr = 0
    xr = 0.0_wp ; yr = 0.0_wp ; zr = 0.0_wp

  End Subroutine allocate_defects_arrays

End Module defects_module
