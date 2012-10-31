Module tersoff_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global tersoff interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntpter = 0, &
                                          potter = 0


  Logical,           Allocatable, Save :: lfrter(:)

  Integer,           Allocatable, Save :: lstter(:),ltpter(:)

  Real( Kind = wp ), Allocatable, Save :: prmter(:,:),prmter2(:,:)
  Real( Kind = wp ), Allocatable, Save :: vmbp(:,:,:),gmbp(:,:,:)

  Public :: allocate_tersoff_arrays

Contains

  Subroutine allocate_tersoff_arrays()

    Use setup_module, Only : mxsite,mxter,mxpter,mxgrid

    Implicit None

    Integer                   :: nprter
    Integer, Dimension( 1:7 ) :: fail

    nprter = (mxter*(mxter+1))/2

    fail = 0

    Allocate (lfrter(1:Merge(mxsite,0,mxter > 0)),    Stat = fail(1))
    Allocate (lstter(1:mxter),                        Stat = fail(2))
    Allocate (ltpter(1:mxter),                        Stat = fail(3))
    Allocate (prmter(1:mxpter,1:mxter),               Stat = fail(4))
    If (potter == 1) Allocate (prmter2(1:nprter,1:2), Stat = fail(5))
    Allocate (vmbp(1:mxgrid,1:nprter,1:3),            Stat = fail(6))
    Allocate (gmbp(1:mxgrid,1:nprter,1:3),            Stat = fail(7))

    If (Any(fail > 0)) Call error(1027)

    lfrter = .false.

    lstter = 0
    ltpter = 0

    prmter  = 0.0_wp
    If (potter == 1) prmter2 = 0.0_wp
    vmbp    = 0.0_wp
    gmbp    = 0.0_wp

  End Subroutine allocate_tersoff_arrays

End Module tersoff_module
