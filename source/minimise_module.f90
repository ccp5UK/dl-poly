Module minimise_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring minimisation routine arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: l_x = .false.

  Real( Kind = wp ),              Save :: passmin(1:5) = (/ &
                                          0.0_wp         ,  & ! cycles counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average cycles
                                          999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
                                          0.0_wp /)           ! maximum cycles


  Real( Kind = wp ), Allocatable, Save :: oxx(:),oyy(:),ozz(:)

  Public :: allocate_minimise_arrays,deallocate_minimise_arrays

Contains

  Subroutine allocate_minimise_arrays()

    Use setup_module, Only : mxatms

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (oxx(1:mxatms),oyy(1:mxatms),ozz(1:mxatms), Stat = fail)

    If (fail > 0) Call error(1038)

    oxx = 0.0_wp ; oyy = 0.0_wp ; ozz = 0.0_wp

  End Subroutine allocate_minimise_arrays

  Subroutine deallocate_minimise_arrays()

    Implicit None

    Integer :: fail

    fail = 0

    Deallocate (oxx,oyy,ozz, Stat = fail)

    If (fail > 0) Call error(1039)

  End Subroutine deallocate_minimise_arrays

End Module minimise_module
