Module vnl_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring VNL conditional update variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
! contrib   - i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,                        Save :: l_vnl = .true., & ! Do update
                                          llvnl = .false.   ! Unconditional VNL
  Real( Kind = wp ),              Save :: skipvnl(1:5) = 0.0_wp

! skipping accumulator
! skipvnl(1) - cycles counter
! skipvnl(2) - access counter
! skipvnl(3) - average cycles
! skipvnl(4) - minimum cycles
! skipvnl(5) - maximum cycles

  Real( Kind = wp ), Allocatable, Save :: xbg(:),ybg(:),zbg(:)

End Module vnl_module
