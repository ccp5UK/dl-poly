Module vdw_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global vdw interaction variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        save :: ld_vdw = .false., &
                                          ls_vdw = .false.

  Integer,                        Save :: ntpvdw = 0


  Integer,           Allocatable, Save :: lstvdw(:),ltpvdw(:)

  Real( Kind = wp ), Allocatable, Save :: prmvdw(:,:)
  Real( Kind = wp ), Allocatable, Save :: vvdw(:,:),gvdw(:,:)

  Public :: allocate_vdw_arrays

Contains

  Subroutine allocate_vdw_arrays()

    Use setup_module, Only : mxvdw,mxpvdw,mxgrid

    Implicit None

    Integer, Dimension( 1:5 ) :: fail

    fail = 0

    Allocate (lstvdw(1:mxvdw),          Stat = fail(1))
    Allocate (ltpvdw(1:mxvdw),          Stat = fail(2))
    Allocate (prmvdw(1:mxpvdw,1:mxvdw), Stat = fail(3))
    Allocate (vvdw(1:mxgrid,1:mxvdw),   Stat = fail(4))
    Allocate (gvdw(1:mxgrid,1:mxvdw),   Stat = fail(5))

    If (Any(fail > 0)) Call error(1022)

    lstvdw = 0
    ltpvdw = 0

    prmvdw = 0.0_wp
    vvdw   = 0.0_wp
    gvdw   = 0.0_wp

  End Subroutine allocate_vdw_arrays

End Module vdw_module
