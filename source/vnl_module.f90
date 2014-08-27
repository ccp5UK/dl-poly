Module vnl_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring VNL conditional update variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
! contrib   - i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: l_vnl = .true., & ! Do update
                                          llvnl = .false.   ! Unconditional VNL

  Real( Kind = wp ),              Save :: skipvnl(1:5) = (/ &
                                          0.0_wp         ,  & ! skips counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average skips
                                          999999999.0_wp ,  & ! minimum skips : ~Huge(1)
                                          0.0_wp /)           ! maximum skips


  Real( Kind = wp ), Allocatable, Save :: xbg(:),ybg(:),zbg(:)

End Module vnl_module
