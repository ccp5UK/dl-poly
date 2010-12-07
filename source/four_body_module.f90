Module four_body_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global four-body potential variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntpfbp = 0


  Logical,           Allocatable, Save :: lfrfbp(:)

  Integer,           Allocatable, Save :: lstfbp(:),ltpfbp(:)

  Real( Kind = wp ), Allocatable, Save :: prmfbp(:,:),rctfbp(:)

  Public :: allocate_four_body_arrays

Contains

  Subroutine allocate_four_body_arrays()

    Use setup_module, Only : mxsite,mxfbp,mxpfbp

    Implicit None

    Integer, Dimension( 1:5 ) :: fail

    fail = 0

    Allocate (lfrfbp(1:Merge(mxsite,0,mxfbp > 0)), Stat = fail(1))
    Allocate (lstfbp(1:mxfbp),                     Stat = fail(2))
    Allocate (ltpfbp(1:mxfbp),                     Stat = fail(3))
    Allocate (prmfbp(1:mxpfbp,1:mxfbp),            Stat = fail(4))
    Allocate (rctfbp(1:mxfbp),                     Stat = fail(5))

    If (Any(fail > 0)) Call error(1024)

    lfrfbp = .false.

    lstfbp = 0
    ltpfbp = 0

    prmfbp = 0.0_wp
    rctfbp = 0.0_wp

  End Subroutine allocate_four_body_arrays

End Module four_body_module
