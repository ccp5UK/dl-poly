Module dihedrals_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global dihedral interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntdihd  = 0 , &
                                          ntdihd1 = 0

  Logical,                        Save :: lx_dih=.false. ! dihedrals with core-shell

  Integer,           Allocatable, Save :: numdih(:),keydih(:)
  Integer,           Allocatable, Save :: lstdih(:,:),listdih(:,:),legdih(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmdih(:,:)

  Public :: allocate_dihedrals_arrays , deallocate_dihedrals_arrays

Contains

  Subroutine allocate_dihedrals_arrays()

    Use setup_module, Only : mxtmls,mxtdih,mxdihd,mxfdih,mxpdih,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numdih(1:mxtmls),          Stat = fail(1))
    Allocate (keydih(1:mxtdih),          Stat = fail(2))
    Allocate (lstdih(1:6,1:mxtdih),      Stat = fail(3))
    Allocate (listdih(0:6,1:mxdihd),     Stat = fail(4))
    Allocate (legdih(0:mxfdih,1:mxatdm), Stat = fail(5))
    Allocate (prmdih(1:mxpdih,1:mxtdih), Stat = fail(6))

    If (Any(fail > 0)) Call error(1020)

    numdih  = 0
    keydih  = 0
    lstdih  = 0
    listdih = 0
    legdih  = 0

    prmdih  = 0.0_wp

  End Subroutine allocate_dihedrals_arrays

  Subroutine deallocate_dihedrals_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numdih,lstdih, Stat = fail)

    If (fail > 0) Call error(1033)

  End Subroutine deallocate_dihedrals_arrays

End Module dihedrals_module
