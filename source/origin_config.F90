Subroutine origin_config(megatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for translating the origin of the MD box as
! defined in CONFIG by the (xorg,yorg,zorg) vector and saving it in
! CFGORG
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use config_module,      Only : imcon,cell,natms,xxx,yyy,zzz
  Use development_module, Only : lvcforg,xorg,yorg,zorg

  Implicit None

  Integer, Intent( In    ) :: megatm

  Character ( Len = 6 ) :: name
  Integer               :: i,nstep
  Real( Kind = wp )     :: tstep,time

! Translate

  Do i=1,natms
     xxx(i)=xxx(i)+xorg
     yyy(i)=yyy(i)+yorg
     zzz(i)=zzz(i)+zorg
  End Do

! Restore periodic boundaries

  Call pbcshift(imcon,cell,natms,xxx,yyy,zzz)

! Write REVCON

  name   = 'CFGORG' ! file name
  nstep  = 0        ! no steps done
  tstep  = 0.0_wp   ! no step exists
  time   = 0.0_wp   ! time is not relevant

  Call write_config(name,lvcforg,megatm,nstep,tstep,time)

End Subroutine origin_config
