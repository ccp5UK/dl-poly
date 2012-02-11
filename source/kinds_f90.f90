Module kinds_f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module setting global working precision for Real numbers and
! double precision for Integer numbers @ compile time
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer, Parameter :: sp = Selected_Real_Kind(6,37)    ! single real
  Integer, Parameter :: dp = Selected_Real_Kind(15,307)  ! double real
  Integer, Parameter :: qp = Selected_Real_Kind(33,4931) ! quadrupole real

  Integer, Parameter :: wp = dp                          ! working real

  Integer, Parameter :: ip = Selected_Int_Kind(12)       ! long integer

End Module kinds_f90
