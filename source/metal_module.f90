Module metal_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global metal interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        save :: ld_met = .false.

  Integer,                        Save :: ntpmet = 0


  Integer,           Allocatable, Save :: lstmet(:),ltpmet(:)

  Real( Kind = wp ), Allocatable, Save :: prmmet(:,:)
  Real( Kind = wp ), Allocatable, Save :: vmet(:,:,:),dmet(:,:,:),fmet(:,:,:)
  Real( Kind = wp ), Allocatable, Save :: elrcm(:),vlrcm(:)

  Public :: allocate_metal_arrays

Contains

  Subroutine allocate_metal_arrays()

    Use setup_module, Only : mxmet,mxpmet,mxgrid,mxatyp

    Implicit None

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (lstmet(1:mxmet),            Stat = fail(1))
    Allocate (ltpmet(1:mxmet),            Stat = fail(2))
    Allocate (prmmet(1:mxpmet,1:mxmet),   Stat = fail(3))
    Allocate (vmet(1:mxgrid,1:mxmet,1:2), Stat = fail(4))
    Allocate (dmet(1:mxgrid,1:mxmet,1:2), Stat = fail(5))
    Allocate (fmet(1:mxgrid,1:mxmet,1:2), Stat = fail(6))
    Allocate (elrcm(0:mxatyp),            Stat = fail(7))
    Allocate (vlrcm(0:mxatyp),            Stat = fail(8))

    If (Any(fail > 0)) Call error(1023)

    lstmet = 0
    ltpmet = 0

    prmmet = 0.0_wp
    vmet   = 0.0_wp
    dmet   = 0.0_wp
    fmet   = 0.0_wp

    elrcm  = 0.0_wp
    vlrcm  = 0.0_wp

  End Subroutine allocate_metal_arrays

End Module metal_module
