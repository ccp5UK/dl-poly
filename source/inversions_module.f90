Module inversions_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global valence invle interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: lt_inv=.false. ! no tabulated potentials opted

  Integer,                        Save :: ntinv  = 0 , &
                                          ntinv1 = 0


  Integer,           Allocatable, Save :: numinv(:),keyinv(:)
  Integer,           Allocatable, Save :: lstinv(:,:),listinv(:,:),leginv(:,:)

  Real( Kind = wp ), Allocatable, Save :: prminv(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpinv(:)
  Real( Kind = wp ), Allocatable, Save :: vinv(:,:),ginv(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfinv(:),typinv(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstinv(:,:)

  Public :: allocate_inversions_arrays , deallocate_inversions_arrays , &
            allocate_invr_pot_arrays , allocate_invr_dst_arrays

Contains

  Subroutine allocate_inversions_arrays()

    Use setup_module, Only : mxtmls,mxtinv,mxinv,mxfinv,mxpinv,mxginv,mxatdm

    Implicit None

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numinv(1:mxtmls),          Stat = fail(1))
    Allocate (keyinv(1:mxtinv),          Stat = fail(2))
    Allocate (lstinv(1:4,1:mxtinv),      Stat = fail(3))
    Allocate (listinv(0:4,1:mxinv),      Stat = fail(4))
    Allocate (leginv(0:mxfinv,1:mxatdm), Stat = fail(5))
    Allocate (prminv(1:mxpinv,1:mxtinv), Stat = fail(6))
    If (lt_inv) &
    Allocate (ltpinv(0:mxtinv),          Stat = fail(7))
    If (mxginv > 0) &
    Allocate (ldfinv(0:mxtinv),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1021)

    numinv  = 0
    keyinv  = 0
    lstinv  = 0
    listinv = 0
    leginv  = 0

    prminv  = 0.0_wp

    If (lt_inv) &
    ltpinv  = 0

    If (mxginv > 0) &
    ldfinv  = 0

  End Subroutine allocate_inversions_arrays

  Subroutine deallocate_inversions_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numinv,lstinv, Stat = fail)

    If (fail > 0) Call error(1034)

  End Subroutine deallocate_inversions_arrays

  Subroutine allocate_invr_pot_arrays()

    Use setup_module, Only : mxgrid

    Implicit None

    Integer :: fail(1:2)

    fail = 0

    Allocate (vinv(0:mxgrid,1:ltpinv(0)), Stat = fail(1))
    Allocate (ginv(0:mxgrid,1:ltpinv(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1078)

    vinv = 0.0_wp
    ginv = 0.0_wp

  End Subroutine allocate_invr_pot_arrays

  Subroutine allocate_invr_dst_arrays()

    Use setup_module, Only : mxginv

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (typinv(-11:4,1:ldfinv(0)),dstinv(0:mxginv,1:ldfinv(0)), Stat = fail)

    If (fail > 0) Call error(1079)

    typinv = 0
    dstinv = 0.0_wp

  End Subroutine allocate_invr_dst_arrays

End Module inversions_module
