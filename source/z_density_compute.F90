Subroutine z_density_compute()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating z-density profile from
! accumulated data
!
! copyright - daresbury laboratory
! author    - t.forester march 1994
! amended   - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module,  Only : mxgrdf,nrite,nzdndt
  Use site_module,   Only : ntpatm,unqatm
  Use config_module, Only : cfgname,cell,volm
  Use z_density_module

  Implicit None

  Integer           :: j,k
  Real( Kind = wp ) :: zlen,delr,dvolz,factor,rho,rho1,rrr,sum,sum1

  If (idnode == 0) Write(nrite,"(/,/,12x,'Z DENSITY PROFILES',/,/, &
     & 'calculated using ',i8,' configurations')") ncfzdn

! open Z density file and write headers

  If (idnode == 0) Then
     Open(Unit=nzdndt, File='ZDNDAT', Status='replace')
     Write(nzdndt,'(a)') cfgname
     Write(nzdndt,'(2i10)') ntpatm,mxgrdf
  End If

! length of cell in z direction

  zlen=Abs(cell(3))+Abs(cell(6))+Abs(cell(9))

! grid interval for density profiles

  delr=zlen/Real(mxgrdf,wp)

! volume of z-strip

  dvolz=(volm/zlen)*delr

! normalisation factor

  ncfzdn=Max(ncfzdn,1)
  factor=1.0_wp/(Real(ncfzdn,wp)*dvolz)

! for every species

  Do k=1,ntpatm
     If (idnode == 0) Then
        Write(nrite,"(/,'rho(r)  :',a8,/,/,8x,'r',6x,'rho',9x,'n(r)',/)") unqatm(k)
        Write(nzdndt,'(a8)') unqatm(k)
     End If

! global sum of data on all nodes

     If (mxnode > 1) Call gsum(zdens(1:mxgrdf,k))

! running integration of z-density

     sum=0.0_wp

! loop over distances

     Do j=1,mxgrdf
        rrr=(Real(j,wp)-0.5_wp)*delr - zlen*0.5_wp
        rho=zdens(j,k)*factor
        sum=sum + rho*dvolz

! null it if < 1.0e-6_wp

        If (rho < 1.0e-6_wp) Then
           rho1 = 0.0_wp
        Else
           rho1 = rho
        End If

        If (sum < 1.0e-6_wp) Then
           sum1 = 0.0_wp
        Else
           sum1 = sum
        End If

! print out information

        If (idnode == 0) Then
           Write(nzdndt,"(1p,2e14.6)") rrr,rho
           Write(nrite,"(f10.4,1p,2e14.6)") rrr,rho1,sum1
        End If
     End Do
  End Do

  If (idnode == 0) Close(Unit=nzdndt)

End Subroutine z_density_compute
