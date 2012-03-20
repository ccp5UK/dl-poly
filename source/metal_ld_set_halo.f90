Subroutine metal_ld_set_halo(rho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of density data between
! neighbouring domains/nodes
!
! Note: all depends on the ixyz halo array set in set_halo, this assumes
!       that (i) rmet=rcut! as well as (ii) all the error checks in there
!
! copyright - daresbury laboratory
! amended   - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : mxnode,gcheck
  Use setup_module,  Only : mxatms
  Use config_module, Only : natms,nlast

  Implicit None

  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: rho

  Logical           :: safe
  Integer           :: mlast

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call metal_ld_export(-1,mlast,rho)
  Call metal_ld_export( 1,mlast,rho)

! exchange atom data in -/+ y directions

  Call metal_ld_export(-2,mlast,rho)
  Call metal_ld_export( 2,mlast,rho)

! exchange atom data in -/+ z directions

  Call metal_ld_export(-3,mlast,rho)
  Call metal_ld_export( 3,mlast,rho)

! check atom totals after data transfer

  safe=(mlast == nlast)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(96)

End Subroutine metal_ld_set_halo
