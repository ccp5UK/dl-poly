Subroutine mpoles_rotmat_set_halo()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to set and infinitesimally rotate native multipoles
! rotation matrices and then arrange exchange of halo data between
! neighbouring domains/nodes
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
! contrib   - h.a.boateng february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module,  Only : nrite,mxompl,mxatms
  Use config_module, Only : natms,nlast,ixyz
  Use mpoles_module, Only : mplflg

  Implicit None

  Logical :: safe
  Integer :: fail,i,mlast

  Integer, Allocatable :: ixyz0(:)

  Do i=1,natms
     mplflg(i)=0
     If (mxompl < 3) Then
        Call rotate_mpoles_d(i)
     Else
        Call rotate_mpoles(i)
     End If
  End Do

! Communicate the matrices in the halo

  fail = 0
  Allocate (ixyz0(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'mpoles_rotmat_set_halo allocation failure, node: ', idnode
     Call error(0)
  End If
  ixyz0(1:nlast) = ixyz(1:nlast)

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call mpoles_rotmat_export(-1,mlast,ixyz0)
  Call mpoles_rotmat_export( 1,mlast,ixyz0)

! exchange atom data in -/+ y directions

  Call mpoles_rotmat_export(-2,mlast,ixyz0)
  Call mpoles_rotmat_export( 2,mlast,ixyz0)

! exchange atom data in -/+ z directions

  Call mpoles_rotmat_export(-3,mlast,ixyz0)
  Call mpoles_rotmat_export( 3,mlast,ixyz0)

! check atom totals after data transfer

  safe=(mlast == nlast)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(174)

  Deallocate (ixyz0, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'mpoles_rotmat_set_halo deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine mpoles_rotmat_set_halo
