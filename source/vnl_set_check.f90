Subroutine vnl_set_check()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine for (re)setting the conditional VNL checkpoint -
! xbg,ybg,zbg at the end of set_halo_particles
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,   Only : idnode
  Use setup_module,   Only : nrite,mxatms
  Use config_module,  Only : nlast,xxx,yyy,zzz
  Use vnl_module,     Only : llvnl,xbg,ybg,zbg

  Implicit None

  Logical, Save :: newjob=.true.

  Integer :: fail

  If (.not.llvnl) Return

  If (newjob) Then ! Init set
     newjob = .false.

     fail = 0
     Allocate (xbg(1:mxatms),ybg(1:mxatms),zbg(1:mxatms), Stat = fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'vnl_set_check allocation failure, node: ', idnode
        Call error(0)
     End If

! CVNL state and skippage accumulators are initialised in vnl_module
!
!    l_vnl = .true.
!    skipvnl(1) - cycles counter
!    skipvnl(2) - access counter
!    skipvnl(3) - average cycles
!    skipvnl(4) - minimum cycles
!    skipvnl(5) - maximum cycles
  End If

! set tracking point

  xbg(1:nlast)=xxx(1:nlast)
  ybg(1:nlast)=yyy(1:nlast)
  zbg(1:nlast)=zzz(1:nlast)

End Subroutine vnl_set_check
