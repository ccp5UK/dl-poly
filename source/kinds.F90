Module kinds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module setting global working precision for Real numbers and
! double precision for Integer numbers @ compile time
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Use iso_fortran_env, Only : real32, real64, real128,int64,int32
  Implicit None

  Integer, Parameter :: sp = real32                      
  !! single real
  Integer, Parameter :: dp = real64  
  !! double real
  Integer, Parameter :: qp = real128 
  !! quadrupole real

  Integer, Parameter :: wp = dp
  
  !! working real
  
  Integer, Parameter :: ni = int32       
  !! normal 4byte integer

  Integer, Parameter :: li = int64       
  !! long integer

  Integer, Parameter :: wi = ni
    !! integer precision

End Module kinds
