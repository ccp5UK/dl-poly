Module bonds_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global bond interaction variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
! contrib   - a.v.brukhno april 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: lt_bnd = .false. ! no tabulated potentials opted

  Integer,                        Save :: ntbond  = 0 , &
                                          ntbond1 = 0 , &
                                          ncfbnd  = 0

  Real( Kind = wp ),              Save :: rcbnd = 0.0_wp


  Integer,           Allocatable, Save :: numbonds(:),keybnd(:)
  Integer,           Allocatable, Save :: lstbnd(:,:),listbnd(:,:),legbnd(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmbnd(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpbnd(:)
  Real( Kind = wp ), Allocatable, Save :: vbnd(:,:),gbnd(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfbnd(:),typbnd(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstbnd(:,:)

  Public :: allocate_bonds_arrays , deallocate_bonds_arrays , &
            allocate_bond_pot_arrays , allocate_bond_dst_arrays

Contains

  Subroutine allocate_bonds_arrays()

    Use setup_module, Only : mxtmls,mxtbnd,mxbond,mxfbnd,mxpbnd,mxgbnd1,mxatdm

    Implicit None

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numbonds(1:mxtmls),        Stat = fail(1))
    Allocate (keybnd(1:mxtbnd),          Stat = fail(2))
    Allocate (lstbnd(1:2,1:mxtbnd),      Stat = fail(3))
    Allocate (listbnd(0:2,1:mxbond),     Stat = fail(4))
    Allocate (legbnd(0:mxfbnd,1:mxatdm), Stat = fail(5))
    Allocate (prmbnd(1:mxpbnd,1:mxtbnd), Stat = fail(6))
    If (lt_bnd) &
    Allocate (ltpbnd(0:mxtbnd),          Stat = fail(7))
    If (mxgbnd1 > 0) &
    Allocate (ldfbnd(0:mxtbnd),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1014)

    numbonds = 0
    keybnd   = 0
    lstbnd   = 0
    listbnd  = 0
    legbnd   = 0

    prmbnd   = 0.0_wp

    If (lt_bnd) &
    ltpbnd   = 0

    If (mxgbnd1 > 0) &
    ldfbnd   = 0

  End Subroutine allocate_bonds_arrays

  Subroutine deallocate_bonds_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numbonds,lstbnd, Stat = fail)

    If (fail > 0) Call error(1029)

  End Subroutine deallocate_bonds_arrays

  Subroutine allocate_bond_pot_arrays()

    Use setup_module, Only : mxgbnd

    Implicit None

    Integer :: fail(1:2)

    fail = 0

    Allocate (vbnd(-1:mxgbnd,1:ltpbnd(0)), Stat = fail(1))
    Allocate (gbnd(-1:mxgbnd,1:ltpbnd(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1072)

    vbnd = 0.0_wp
    gbnd = 0.0_wp

  End Subroutine allocate_bond_pot_arrays

  Subroutine allocate_bond_dst_arrays()

    Use setup_module, Only : mxgbnd1

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (typbnd(-1:2,1:ldfbnd(0)),dstbnd(1:mxgbnd1,1:ldfbnd(0)), Stat = fail)

    If (fail > 0) Call error(1073)

    typbnd = 0
    dstbnd = 0.0_wp

  End Subroutine allocate_bond_dst_arrays

End Module bonds_module
