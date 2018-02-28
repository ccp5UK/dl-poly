Module angles_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global valence angle interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov april 2014
! contrib   - a.v.brukhno march 2014 (itramolecular TPs & PDFs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp

  Implicit None

  Logical,                        Save :: lt_ang = .false. ! no tabulated potentials opted

  Integer,                        Save :: ntangl  = 0 , &
                                          ntangl1 = 0 , &
                                          ncfang  = 0


  Integer,           Allocatable, Save :: numang(:),keyang(:)
  Integer,           Allocatable, Save :: lstang(:,:),listang(:,:),legang(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmang(:,:)

! Possible tabulated calculation arrays

  Integer,           Allocatable, Save :: ltpang(:)
  Real( Kind = wp ), Allocatable, Save :: vang(:,:),gang(:,:)

! Possible distribution arrays

  Integer,           Allocatable, Save :: ldfang(:),typang(:,:)
  Real( Kind = wp ), Allocatable, Save :: dstang(:,:)

  Public :: allocate_angles_arrays , deallocate_angles_arrays , &
            allocate_angl_pot_arrays , allocate_angl_dst_arrays

Contains

  Subroutine allocate_angles_arrays()

    Use setup_module, Only : mxtmls,mxtang,mxangl,mxfang,mxpang,mxgang1,mxatdm

    Implicit None

    Integer, Dimension( 1:8 ) :: fail

    fail = 0

    Allocate (numang(1:mxtmls),          Stat = fail(1))
    Allocate (keyang(1:mxtang),          Stat = fail(2))
    Allocate (lstang(1:3,1:mxtang),      Stat = fail(3))
    Allocate (listang(0:3,1:mxangl),     Stat = fail(4))
    Allocate (legang(0:mxfang,1:mxatdm), Stat = fail(5))
    Allocate (prmang(1:mxpang,1:mxtang), Stat = fail(6))
    If (lt_ang) &
    Allocate (ltpang(0:mxtang),          Stat = fail(7))
    If (mxgang1 > 0) &
    Allocate (ldfang(0:mxtang),          Stat = fail(8))

    If (Any(fail > 0)) Call error(1013)

    numang  = 0
    keyang  = 0
    lstang  = 0
    listang = 0
    legang  = 0

    prmang  = 0.0_wp

    If (lt_ang) &
    ltpang  = 0

    If (mxgang1 > 0) &
    ldfang  = 0

  End Subroutine allocate_angles_arrays

  Subroutine deallocate_angles_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numang,lstang, Stat = fail)

    If (fail > 0) Call error(1028)

  End Subroutine deallocate_angles_arrays

  Subroutine allocate_angl_pot_arrays()

    Use setup_module, Only : mxgang

    Implicit None

    Integer :: fail(1:2)

    fail = 0

    Allocate (vang(-1:mxgang,1:ltpang(0)), Stat = fail(1))
    Allocate (gang(-1:mxgang,1:ltpang(0)), Stat = fail(2))

    If (Any(fail > 0)) Call error(1074)

    vang = 0.0_wp
    gang = 0.0_wp

  End Subroutine allocate_angl_pot_arrays

  Subroutine allocate_angl_dst_arrays()

    Use setup_module, Only : mxgang1

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (typang(-1:3,1:ldfang(0)),dstang(1:mxgang1,1:ldfang(0)), Stat = fail)

    If (fail > 0) Call error(1075)

    typang = 0
    dstang = 0.0_wp

  End Subroutine allocate_angl_dst_arrays

End Module angles_module
