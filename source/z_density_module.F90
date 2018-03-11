Module z_density_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global z-density variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp

  Implicit None

  Integer,                        Save :: ncfzdn = 0

  Real( Kind = wp ), Allocatable, Save :: zdens(:,:)

  Public :: allocate_z_density_arrays

Contains

  Subroutine allocate_z_density_arrays()

    Use setup_module, Only : mxatyp,mxgrdf

    Implicit None

    Integer :: fail

    fail = 0

    Allocate (zdens(1:mxgrdf,1:mxatyp), Stat = fail)

    If (fail > 0) Call error(1016)

    zdens = 0.0_wp

  End Subroutine allocate_z_density_arrays

End Module z_density_module
