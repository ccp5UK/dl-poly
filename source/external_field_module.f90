Module external_field_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global external field variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

! Only one type of field can be applied on the system (keyfld is a scalar)

  Integer,                        Save :: keyfld = 0


  Real( Kind = wp ), Allocatable, Save :: prmfld(:)

  Public :: allocate_external_field_arrays

Contains

  Subroutine allocate_external_field_arrays()

    Use setup_module, Only : mxpfld

    Implicit None

    Integer, Dimension( 1:1 ) :: fail

    fail = 0

    Allocate (prmfld(mxpfld), Stat = fail(1))

    If (Any(fail > 0)) Call error(1019)

    prmfld = 0.0_wp

  End Subroutine allocate_external_field_arrays

End Module external_field_module
