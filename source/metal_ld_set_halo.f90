Subroutine metal_ld_set_halo(rho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of density data between
! neighbouring domains/nodes
!
! Note: all depends on the ixyz halo array set in set_halo, this
!       assumes that (i) rmet=rcut! & (ii) all the error checks in there
!
! copyright - daresbury laboratory
! amended   - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module,  Only : nrite,mxatms
  Use config_module, Only : natms,nlast,ltg

  Implicit None

  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: rho

  Logical           :: safe
  Integer           :: fail,mlast,i

  Integer,           Dimension( : ), Allocatable :: iwrk

  fail=0
  Allocate (iwrk(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_set_halo allocation failure, node: ', idnode
     Call error(0)
  End If

! Copy domain only ltg in iwrk

  mlast=natms ! No halo
  Do i=1,mlast
     iwrk(i)=ltg(i)
  End Do

! exchange atom data in -/+ x directions

  Call metal_ld_export(-1,mlast,iwrk,rho)
  Call metal_ld_export( 1,mlast,iwrk,rho)

! exchange atom data in -/+ y directions

  Call metal_ld_export(-2,mlast,iwrk,rho)
  Call metal_ld_export( 2,mlast,iwrk,rho)

! exchange atom data in -/+ z directions

  Call metal_ld_export(-3,mlast,iwrk,rho)
  Call metal_ld_export( 3,mlast,iwrk,rho)

! check atom totals after data transfer

  safe=(mlast == nlast)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(96)

! check incoming atomic density assignments

  Do i=natms+1,nlast
     safe=(ltg(i) == iwrk(i))
  End Do
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(98)

  Deallocate (iwrk, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_set_halo deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine metal_ld_set_halo
