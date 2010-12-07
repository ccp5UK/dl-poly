Module inversions_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global valence invle interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntinv  = 0 , &
                                          ntinv1 = 0


  Integer,           Allocatable, Save :: numinv(:),keyinv(:)
  Integer,           Allocatable, Save :: lstinv(:,:),listinv(:,:),leginv(:,:)

  Real( Kind = wp ), Allocatable, Save :: prminv(:,:)

  Public :: allocate_inversions_arrays , deallocate_inversions_arrays

Contains

  Subroutine allocate_inversions_arrays()

    Use setup_module, Only : mxtmls,mxtinv,mxinv,mxfinv,mxpinv,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numinv(1:mxtmls),          Stat = fail(1))
    Allocate (keyinv(1:mxtinv),          Stat = fail(2))
    Allocate (lstinv(1:4,1:mxtinv),      Stat = fail(3))
    Allocate (listinv(0:4,1:mxinv),      Stat = fail(4))
    Allocate (leginv(0:mxfinv,1:mxatdm), Stat = fail(5))
    Allocate (prminv(1:mxpinv,1:mxtinv), Stat = fail(6))

    If (Any(fail > 0)) Call error(1021)

    numinv  = 0
    keyinv  = 0
    lstinv  = 0
    listinv = 0
    leginv  = 0

    prminv  = 0.0_wp

  End Subroutine allocate_inversions_arrays

  Subroutine deallocate_inversions_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numinv,lstinv, Stat = fail)

    If (fail > 0) Call error(1034)

  End Subroutine deallocate_inversions_arrays

End Module inversions_module
