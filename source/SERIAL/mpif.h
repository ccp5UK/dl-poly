!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mpif.h
!
! dl_poly_4 inclusion for serial compilation of MPI calls
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TYPES

  Integer, Parameter :: MPI_LOGICAL          = 1
  Integer, Parameter :: MPI_CHARACTER        = 1
  Integer, Parameter :: MPI_INTEGER          = 1
  Integer, Parameter :: MPI_REAL             = 1
! MPI_REAL8 is apparently not in the strict MPI2 standard
! It is just an optional data type in the FORTRAN Bindings
  Integer, Parameter :: MPI_DOUBLE_PRECISION = 1
  Integer, Parameter :: MPI_REAL16           = 1

! OPERANDS

  Integer, Parameter :: MPI_LAND = 1
  Integer, Parameter :: MPI_MAX  = 1
  Integer, Parameter :: MPI_MIN  = 1
  Integer, Parameter :: MPI_SUM  = 1

! GLOBAL

  Integer, Parameter :: MPI_UNDEFINED   = 0
  Integer, Parameter :: MPI_COMM_NULL   = 0
  Integer, Parameter :: MPI_COMM_WORLD  = 1
  Integer, Parameter :: MPI_STATUS_SIZE = 1

! MPI_MODES for MPI-I/O

  Integer, Parameter :: MPI_MODE_CREATE          =      1
  Integer, Parameter :: MPI_MODE_EXCL            =      2
  Integer, Parameter :: MPI_MODE_RDONLY          =     10
  Integer, Parameter :: MPI_MODE_WRONLY          =     20
  Integer, Parameter :: MPI_MODE_RDWR            =     30
  Integer, Parameter :: MPI_MODE_DELETE_ON_CLOSE =    100
  Integer, Parameter :: MPI_MODE_SEQUENTIAL      =   1000
  Integer, Parameter :: MPI_MODE_APPEND          =  10000
  Integer, Parameter :: MPI_MODE_UNIQUE_OPEN     = 100000
