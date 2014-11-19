Module dpd_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global dpd variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,           Save :: l_dpd = .false. ! no dpd

  Real( Kind = wp ), Save :: virdpd      = 0.0_wp , &
                             strdpd(1:9) = 0.0_wp

  Integer,           Allocatable, Save :: lstdpd(:,:)
  Real( Kind = wp ), Allocatable, Save :: gamdpd(:),sigdpd(:)

  Public :: allocate_dpd_arrays

Contains

  Subroutine allocate_dpd_arrays()

    Use setup_module, Only : mxvdw

    Implicit None

    Integer :: fail

    If (.not.l_dpd) Return

    fail = 0

    Allocate (gamdpd(0:mxvdw),sigdpd(1:mxvdw), Stat = fail)

    If (fail > 0) Call error(1081)

    gamdpd = 0.0_wp ; sigdpd = 0.0_wp

  End Subroutine allocate_dpd_arrays

End Module dpd_module
