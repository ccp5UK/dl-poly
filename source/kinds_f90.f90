Module kinds_f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module setting global working precision for Real numbers and
! double precision for Integer numbers @ compile time
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Implicit None

  Integer, Parameter :: wp = Selected_Real_Kind(14,300)
  Integer, Parameter :: ip = Selected_Int_Kind(12)

End Module kinds_f90
