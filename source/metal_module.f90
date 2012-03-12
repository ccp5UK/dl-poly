Module metal_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global metal interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        save :: lt_met = .false., &
                                          ld_met = .false.

  Integer,                        Save :: ntpmet = 0


  Integer,           Allocatable, Save :: lstmet(:),ltpmet(:)

  Real( Kind = wp ), Allocatable, Save :: prmmet(:,:)

  Real( Kind = wp ), Allocatable, Save :: elrcm(:),vlrcm(:)

! Possible tabulated calculation arrays

  Real( Kind = wp ), Allocatable, Save :: vmet(:,:,:),dmet(:,:,:),fmet(:,:,:)

  Public :: allocate_metal_arrays, &
            allocate_metal_table_arrays

Contains

  Subroutine allocate_metal_arrays()

    Use setup_module, Only : mxmet,mxpmet,mxatyp

    Implicit None

    Integer, Dimension( 1:5 ) :: fail

    fail = 0

    Allocate (lstmet(1:mxmet),            Stat = fail(1))
    Allocate (ltpmet(1:mxmet),            Stat = fail(2))
    Allocate (prmmet(1:mxpmet,1:mxmet),   Stat = fail(3))
    Allocate (elrcm(0:mxatyp),            Stat = fail(4))
    Allocate (vlrcm(0:mxatyp),            Stat = fail(5))

    If (Any(fail > 0)) Call error(1023)

    lstmet = 0
    ltpmet = 0

    prmmet = 0.0_wp

    elrcm  = 0.0_wp
    vlrcm  = 0.0_wp

  End Subroutine allocate_metal_arrays

  Subroutine allocate_metal_table_arrays()

    Use setup_module, Only : mxmet,mxgrid

    Implicit None

    Integer, Dimension( 1:3 ) :: fail

    fail = 0

    Allocate (vmet(1:mxgrid,1:mxmet,1:2), Stat = fail(1))
    Allocate (dmet(1:mxgrid,1:mxmet,1:2), Stat = fail(2))
    Allocate (fmet(1:mxgrid,1:mxmet,1:2), Stat = fail(3))

    If (Any(fail > 0)) Call error(1069)

    vmet   = 0.0_wp
    dmet   = 0.0_wp
    fmet   = 0.0_wp

  End Subroutine allocate_metal_table_arrays

End Module metal_module
