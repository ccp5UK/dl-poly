Module dpd_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global DPD variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Integer,           Save :: keydpd = 0 ! no DPD

  Real( Kind = wp ), Save :: virdpd      = 0.0_wp , &
                             strdpd(1:9) = 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: gamdpd(:),sigdpd(:)

  Public :: allocate_dpd_arrays

Contains

  Subroutine allocate_dpd_arrays()

    Use setup_module, Only : mxvdw

    Implicit None

    Integer :: fail

    If (keydpd == 0) Return

    fail = 0

    Allocate (gamdpd(0:mxvdw),sigdpd(1:mxvdw), Stat = fail)

    If (fail > 0) Call error(1081)

    gamdpd = 0.0_wp ; sigdpd = 0.0_wp

  End Subroutine allocate_dpd_arrays

End Module dpd_module
