Module vdw_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global vdw interaction variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: lt_vdw = .false., & ! no tabulated potentials are present
                                          ld_vdw = .false., & ! no direct calculations are opted
                                          ls_vdw = .false.    ! no force-shifting is opted

  Integer,                        Save :: ntpvdw = 0, &       ! number of 2 body interactions
                                          mxtvdw = 0          ! type of mixing


  Integer,           Allocatable, Save :: lstvdw(:),ltpvdw(:)

  Real( Kind = wp ), Allocatable, Save :: prmvdw(:,:),sigeps(:,:)

  Real( Kind = wp ),              Save :: elrc   = 0.0_wp, &
                                          virlrc = 0.0_wp

! Possible tabulated calculation arrays

  Real( Kind = wp ), Allocatable, Save :: vvdw(:,:),gvdw(:,:)

! Possible force-shifting arrays

  Real( Kind = wp ), Allocatable, Save :: afs(:),bfs(:)

  Public :: allocate_vdw_arrays, allocate_vdw_table_arrays, &
            allocate_vdw_direct_fs_arrays

Contains

  Subroutine allocate_vdw_arrays()

    Use setup_module, Only : mxvdw,mxpvdw

    Implicit None

    Integer, Dimension( 1:4 ) :: fail

    fail = 0

    Allocate (lstvdw(1:mxvdw),          Stat = fail(1))
    Allocate (ltpvdw(1:mxvdw),          Stat = fail(2))
    Allocate (prmvdw(1:mxpvdw,1:mxvdw), Stat = fail(3))
    Allocate (sigeps(1:2,1:mxvdw),      Stat = fail(4))

    If (Any(fail > 0)) Call error(1022)

    lstvdw = 0
    ltpvdw = 0

    prmvdw = 0.0_wp
    sigeps = 0.0_wp

  End Subroutine allocate_vdw_arrays

  Subroutine allocate_vdw_table_arrays()

    Use setup_module, Only : mxvdw,mxgvdw

    Implicit None

    Integer, Dimension( 1:2 ) :: fail

    fail = 0

    Allocate (vvdw(0:mxgvdw,1:mxvdw),   Stat = fail(1))
    Allocate (gvdw(0:mxgvdw,1:mxvdw),   Stat = fail(2))

    If (Any(fail > 0)) Call error(1063)

    vvdw = 0.0_wp
    gvdw = 0.0_wp

  End Subroutine allocate_vdw_table_arrays

  Subroutine allocate_vdw_direct_fs_arrays()

    Use setup_module, Only : mxvdw

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (afs(1:mxvdw),bfs(1:mxvdw), Stat = fail)

    If (fail > 0) Call error(1066)

    afs  = 0.0_wp
    bfs  = 0.0_wp

  End Subroutine allocate_vdw_direct_fs_arrays

End Module vdw_module
