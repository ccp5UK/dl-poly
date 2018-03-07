Subroutine z_density_collect()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for z-density profile
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module,      Only : mxgrdf
  Use configuration,     Only : cell,natms,ltype,zzz
  Use z_density_module

  Implicit None

  Integer           :: i,k,l
  Real( Kind = wp ) :: zlen,zleno2,rdelr

! accumulator

  ncfzdn=ncfzdn+1

! length of cell in z direction

  zlen=Abs(cell(3))+Abs(cell(6))+Abs(cell(9))

! half of z length

  zleno2=zlen*0.5_wp

! grid interval for density profiles

  rdelr=Real(mxgrdf,wp)/zlen

! set up atom iatm type and accumulate statistic

  Do i=1,natms
     k=ltype(i)

     l=Min(1+Int((zzz(i)+zleno2)*rdelr),mxgrdf)

     zdens(l,k)=zdens(l,k) + 1.0_wp
  End Do

End Subroutine z_density_collect
