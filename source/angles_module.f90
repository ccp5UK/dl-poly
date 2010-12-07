Module angles_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global valence angle interaction variables
! and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,                        Save :: ntangl  = 0 , &
                                          ntangl1 = 0


  Integer,           Allocatable, Save :: numang(:),keyang(:)
  Integer,           Allocatable, Save :: lstang(:,:),listang(:,:),legang(:,:)

  Real( Kind = wp ), Allocatable, Save :: prmang(:,:)

  Public :: allocate_angles_arrays , deallocate_angles_arrays

Contains

  Subroutine allocate_angles_arrays()

    Use setup_module, Only : mxtmls,mxtang,mxangl,mxfang,mxpang,mxatdm

    Implicit None

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numang(1:mxtmls),          Stat = fail(1))
    Allocate (keyang(1:mxtang),          Stat = fail(2))
    Allocate (lstang(1:3,1:mxtang),      Stat = fail(3))
    Allocate (listang(0:3,1:mxangl),     Stat = fail(4))
    Allocate (legang(0:mxfang,1:mxatdm), Stat = fail(5))
    Allocate (prmang(1:mxpang,1:mxtang), Stat = fail(6))

    If (Any(fail > 0)) Call error(1013)

    numang  = 0
    keyang  = 0
    lstang  = 0
    listang = 0
    legang  = 0

    prmang  = 0.0_wp

  End Subroutine allocate_angles_arrays

  Subroutine deallocate_angles_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (numang,lstang, Stat = fail)

    If (fail > 0) Call error(1028)

  End Subroutine deallocate_angles_arrays

End Module angles_module
