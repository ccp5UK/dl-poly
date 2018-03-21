Module halo

  Use comms,  Only : comms_type,gcheck
  Use deport_data, Only : export_atomic_positions
  Use setup_module,  Only : nrite,mxatms
  Use configuration, Only : natms,nlast,ixyz
  Implicit None
  
  Private
  
  Contains
  
  Subroutine refresh_halo_positions(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to refresh the halo positioning data between
! neighbouring domains/nodes when VNL is skipped
!
! Note: all depends on the ixyz halo array set in set_halo, this assumes
!       that (i) rmet=rcut! as well as (ii) all the error checks in there
!
! copyright - daresbury laboratory
! author    - i.t.todorov & i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Type( comms_type ), Intent(InOut) :: comm
  Logical :: safe
  Integer :: fail,mlast

  Integer, Allocatable :: ixyz0(:)

  fail = 0
  Allocate (ixyz0(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'refresh_halo_ppositions allocation failure, node: ', comm%idnode
     Call error(0)
  End If
  ixyz0(1:nlast) = ixyz(1:nlast)

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call export_atomic_positions(-1,mlast,ixyz0,comm)
  Call export_atomic_positions( 1,mlast,ixyz0,comm)

! exchange atom data in -/+ y directions

  Call export_atomic_positions(-2,mlast,ixyz0,comm)
  Call export_atomic_positions( 2,mlast,ixyz0,comm)

! exchange atom data in -/+ z directions

  Call export_atomic_positions(-3,mlast,ixyz0,comm)
  Call export_atomic_positions( 3,mlast,ixyz0,comm)

  safe=(mlast == nlast)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(138)

  Deallocate (ixyz0, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'referesh_halo_positions deallocation failure, node: ', comm%idnode
     Call error(0)
  End If

End Subroutine refresh_halo_positions

  End Module halo
