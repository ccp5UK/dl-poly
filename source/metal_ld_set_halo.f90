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
! amended   - i.t.todorov june 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,  Only : mxnode,gcheck
  Use config_module, Only : natms,nlast

  Implicit None

  Logical           :: safe
  Integer           :: mlast

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call metal_ld_export(-1,mlast)
  Call metal_ld_export( 1,mlast)

! exchange atom data in -/+ y directions

  Call metal_ld_export(-2,mlast)
  Call metal_ld_export( 2,mlast)

! exchange atom data in -/+ z directions

  Call metal_ld_export(-3,mlast)
  Call metal_ld_export( 3,mlast)

! check atom totals after data transfer

  safe=(mlast == nlast)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(96)

End Subroutine metal_ld_set_halo
