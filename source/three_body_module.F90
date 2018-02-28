Module three_body_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global three-body potential variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp

  Implicit None

  Integer,                        Save :: ntptbp = 0


  Logical,           Allocatable, Save :: lfrtbp(:)

  Integer,           Allocatable, Save :: lsttbp(:),ltptbp(:)

  Real( Kind = wp ), Allocatable, Save :: prmtbp(:,:),rcttbp(:)

  Public :: allocate_three_body_arrays

Contains

  Subroutine allocate_three_body_arrays()

    Use setup_module, Only : mxsite,mxtbp,mxptbp

    Implicit None

    Integer, Dimension( 1:5 ) :: fail

    fail = 0

    Allocate (lfrtbp(1:Merge(mxsite,0,mxtbp > 0)), Stat = fail(1))
    Allocate (lsttbp(1:mxtbp),                     Stat = fail(2))
    Allocate (ltptbp(1:mxtbp),                     Stat = fail(3))
    Allocate (prmtbp(1:mxptbp,1:mxtbp),            Stat = fail(4))
    Allocate (rcttbp(1:mxtbp),                     Stat = fail(5))

    If (Any(fail > 0)) Call error(1024)

    lfrtbp = .false.

    lsttbp = 0
    ltptbp = 0

    prmtbp = 0.0_wp
    rcttbp = 0.0_wp

  End Subroutine allocate_three_body_arrays

End Module three_body_module
