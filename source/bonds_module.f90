Module bonds_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global bond interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntbond  = 0 , &
                                          ntbond1 = 0


  Integer,           Allocatable, Save :: numbonds(:),keybnd(:)
  Integer,           Allocatable, Save :: lstbnd(:,:),listbnd(:,:),legbnd(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmbnd(:,:)

  Public :: allocate_bonds_arrays , deallocate_bonds_arrays

Contains

  Subroutine allocate_bonds_arrays()

    Use setup_module, Only : mxtmls,mxtbnd,mxbond,mxfbnd,mxpbnd,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numbonds(1:mxtmls),        Stat = fail(1))
    Allocate (keybnd(1:mxtbnd),          Stat = fail(2))
    Allocate (lstbnd(1:2,1:mxtbnd),      Stat = fail(3))
    Allocate (listbnd(0:2,1:mxbond),     Stat = fail(4))
    Allocate (legbnd(0:mxfbnd,1:mxatdm), Stat = fail(5))
    Allocate (prmbnd(1:mxpbnd,1:mxtbnd), Stat = fail(6))

    If (Any(fail > 0)) Call error(1014)

    numbonds = 0
    keybnd   = 0
    lstbnd   = 0
    listbnd  = 0
    legbnd   = 0

    prmbnd   = 0.0_wp

  End Subroutine allocate_bonds_arrays

  Subroutine deallocate_bonds_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numbonds,lstbnd, Stat = fail)

    If (fail > 0) Call error(1029)

  End Subroutine deallocate_bonds_arrays

End Module bonds_module
