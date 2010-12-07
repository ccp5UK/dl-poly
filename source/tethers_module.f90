Module tethers_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global tether interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntteth = 0


  Integer,           Allocatable, Save :: numteth(:),keytet(:)
  Integer,           Allocatable, Save :: lsttet(:),listtet(:,:),legtet(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmtet(:,:)

  Public :: allocate_tethers_arrays , deallocate_tethers_arrays

Contains

  Subroutine allocate_tethers_arrays()

    Use setup_module, Only : mxtmls,mxtteth,mxteth,mxftet,mxpteth,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numteth(1:mxtmls),           Stat = fail(1))
    Allocate (keytet(1:mxtteth),           Stat = fail(2))
    Allocate (lsttet(1:mxtteth),           Stat = fail(3))
    Allocate (listtet(0:1,1:mxteth),       Stat = fail(4))
    Allocate (legtet(0:mxftet,1:mxatdm),   Stat = fail(5))
    Allocate (prmtet(1:mxpteth,1:mxtteth), Stat = fail(6))

    If (Any(fail > 0)) Call error(1017)

    numteth = 0
    keytet  = 0
    lsttet  = 0
    listtet = 0
    legtet  = 0

    prmtet  = 0.0_wp

  End Subroutine allocate_tethers_arrays

  Subroutine deallocate_tethers_arrays()

    Implicit None

    Integer:: fail

    fail = 0

    Deallocate (numteth,lsttet, Stat = fail)

    If (fail > 0) Call error(1031)

  End Subroutine deallocate_tethers_arrays

End Module tethers_module
