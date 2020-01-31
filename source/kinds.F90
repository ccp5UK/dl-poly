!> Module setting the real and integer kinds used in the program
!>
!> Copyright - Daresbury Laboratory
Module kinds

  Use, Intrinsic :: iso_fortran_env, Only: int16,&
                                           int32,&
                                           int64,&
                                           real128,&
                                           real32,&
                                           real64

  Implicit None

  Private

  !> Single real
  Integer, Parameter, Public :: sp = real32
  !> Double real
  Integer, Parameter, Public :: dp = real64
  !> Quadrupole real
  Integer, Parameter, Public :: qp = real128

  !> Working real
  Integer, Parameter, Public :: wp = dp

  !> Normal 4byte integer
  Integer, Parameter, Public :: ni = int32
  !> Long integer
  Integer, Parameter, Public :: li = int64
  !> Short integer
  Integer, Parameter, Public :: si = int16

  !> Working integer
  Integer, Parameter, Public :: wi = ni
End Module kinds
