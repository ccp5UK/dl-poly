Subroutine metal_ld_set_halo()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of density data between
! neighbouring domains/nodes
!
! Note: all depends on the ixyz halo array set in set_halo, this assumes
!       that (i) rmet=rcut! as well as (ii) all the error checks in there
!
! copyright - daresbury laboratory
! amended   - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module,  Only : nrite,mxatms
  Use config_module, Only : natms,nlast,ixyz

  Implicit None

  Logical :: safe
  Integer :: fail,mlast

  Integer, Allocatable :: ixyz0(:)

  fail = 0
  Allocate (ixyz0(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_set_halo allocation failure, node: ', idnode
     Call error(0)
  End If
  ixyz0(1:nlast) = ixyz(1:nlast)

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call metal_ld_export(-1,mlast,ixyz0)
  Call metal_ld_export( 1,mlast,ixyz0)

! exchange atom data in -/+ y directions

  Call metal_ld_export(-2,mlast,ixyz0)
  Call metal_ld_export( 2,mlast,ixyz0)

! exchange atom data in -/+ z directions

  Call metal_ld_export(-3,mlast,ixyz0)
  Call metal_ld_export( 3,mlast,ixyz0)

! check atom totals after data transfer

  safe=(mlast == nlast)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(96)

  Deallocate (ixyz0, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_set_halo deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine metal_ld_set_halo
