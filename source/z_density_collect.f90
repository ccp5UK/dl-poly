Subroutine z_density_collect()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for z-density profile
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,      Only : mxgrdf,zero_plus
  Use site_module,       Only : numtyp
  Use config_module,     Only : cell,natms,ltype,zzz
  Use statistics_module, Only : numzdn,zdens

  Implicit None

  Integer           :: i,k,l
  Real( Kind = wp ) :: zlen,zleno2,rdelr

! accumulator

  numzdn=numzdn+1

! length of cell in z direction

  zlen=Abs(cell(3))+Abs(cell(6))+Abs(cell(9))

! half of z length

  zleno2=zlen*0.5_wp

! grid interval for density profiles

  rdelr=Real(mxgrdf,wp)/zlen

! set up atom iatm type and exclude it if absent crystallographically

  Do i=1,natms
     k=ltype(i)
     If (numtyp(k) > zero_plus) Then

! apply truncation of potential

        l=Int((zzz(i)+zleno2)*rdelr + 1.0_wp)
        l=Max(1,l)
        l=Min(mxgrdf,l)

! accumulate statistic

        zdens(l,k)=zdens(l,k) + 1.0_wp
     End If
  End Do

End Subroutine z_density_collect
