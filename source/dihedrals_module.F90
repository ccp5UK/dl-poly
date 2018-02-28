Module dihedrals_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global dihedral interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
! contrib   - a.v.brukhno march 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp

  Implicit None

  Logical,                        Save :: lt_dih=.false. , & ! no tabulated potentials opted
                                          lx_dih=.false.     ! dihedrals with core-shell

  Integer,                        Save :: ntdihd  = 0 , &
                                          ntdihd1 = 0 , &
                                          ncfdih  = 0


  Integer,           Allocatable, Save :: numdih(:),keydih(:)
  Integer,           Allocatable, Save :: lstdih(:,:),listdih(:,:),legdih(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmdih(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpdih(:)
  Real( Kind = wp ), Allocatable, Save :: vdih(:,:),gdih(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfdih(:),typdih(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstdih(:,:)

  Public :: allocate_dihedrals_arrays , deallocate_dihedrals_arrays , &
            allocate_dihd_pot_arrays , allocate_dihd_dst_arrays

Contains

  Subroutine allocate_dihedrals_arrays()

    Use setup_module, Only : mxtmls,mxtdih,mxdihd,mxfdih,mxpdih,mxgdih1,mxatdm

    Implicit None

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numdih(1:mxtmls),          Stat = fail(1))
    Allocate (keydih(1:mxtdih),          Stat = fail(2))
    Allocate (lstdih(1:6,1:mxtdih),      Stat = fail(3))
    Allocate (listdih(0:6,1:mxdihd),     Stat = fail(4))
    Allocate (legdih(0:mxfdih,1:mxatdm), Stat = fail(5))
    Allocate (prmdih(1:mxpdih,1:mxtdih), Stat = fail(6))
    If (lt_dih) &
    Allocate (ltpdih(0:mxtdih),          Stat = fail(7))
    If (mxgdih1 > 0) &
    Allocate (ldfdih(0:mxtdih),          Stat = fail(8))


    If (Any(fail > 0)) Call error(1020)

    numdih  = 0
    keydih  = 0
    lstdih  = 0
    listdih = 0
    legdih  = 0

    prmdih  = 0.0_wp

    If (lt_dih) &
    ltpdih  = 0

    If (mxgdih1 > 0) &
    ldfdih  = 0

  End Subroutine allocate_dihedrals_arrays

  Subroutine deallocate_dihedrals_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numdih,lstdih, Stat = fail)

    If (fail > 0) Call error(1033)

  End Subroutine deallocate_dihedrals_arrays

  Subroutine allocate_dihd_pot_arrays()

    Use setup_module, Only : mxgdih

    Implicit None

    Integer :: fail(1:2)

    fail = 0

    Allocate (vdih(-1:mxgdih,1:ltpdih(0)), Stat = fail(1))
    Allocate (gdih(-1:mxgdih,1:ltpdih(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1076)

    vdih = 0.0_wp
    gdih = 0.0_wp

  End Subroutine allocate_dihd_pot_arrays

  Subroutine allocate_dihd_dst_arrays()

    Use setup_module, Only : mxgdih1

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (typdih(-1:4,1:ldfdih(0)),dstdih(1:mxgdih1,1:ldfdih(0)), Stat = fail)

    If (fail > 0) Call error(1077)

    typdih = 0
    dstdih = 0.0_wp

  End Subroutine allocate_dihd_dst_arrays

End Module dihedrals_module
