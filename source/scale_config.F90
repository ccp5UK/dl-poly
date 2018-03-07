Subroutine scale_config(megatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for rescaling the crystallographic information
! from CONFIG to new lattice parameters and saving it in CFGSCL
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use configuration,      Only : cell,natms,xxx,yyy,zzz
  Use development_module, Only : lvcfscl,cels

  Implicit None

  Integer, Intent( In    ) :: megatm

  Character ( Len = 6 ) :: name
  Integer               :: i,nstep
  Real( Kind = wp )     :: rcell(1:9),det,uuu,vvv,www,tstep,time

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Rescale

  Do i=1,natms
     uuu=xxx(i)
     vvv=yyy(i)
     www=zzz(i)

     xxx(i)=rcell(1)*uuu+rcell(4)*vvv+rcell(7)*www
     yyy(i)=rcell(2)*uuu+rcell(5)*vvv+rcell(8)*www
     zzz(i)=rcell(3)*uuu+rcell(6)*vvv+rcell(9)*www

     uuu=xxx(i)
     vvv=yyy(i)
     www=zzz(i)

     xxx(i)=cels(1)*uuu+cels(4)*vvv+cels(7)*www
     yyy(i)=cels(2)*uuu+cels(5)*vvv+cels(8)*www
     zzz(i)=cels(3)*uuu+cels(6)*vvv+cels(9)*www
  End Do

! Write REVCON

  name   = 'CFGSCL' ! file name
  nstep  = 0        ! no steps done
  tstep  = 0.0_wp   ! no step exists
  time   = 0.0_wp   ! time is not relevant

  rcell = cell ; cell = cels
  Call write_config(name,lvcfscl,megatm,nstep,tstep,time)
  cell = rcell

End Subroutine scale_config
