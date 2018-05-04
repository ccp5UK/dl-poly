Module z_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global z-density variables and arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  Use comms,  Only : comms_type,gsum
  Use setup,  Only : mxgrdf,nrite,nzdndt,mxatyp
  Use site,   Only : ntpatm,unqatm
  Use configuration, Only : cfgname,cell,volm,natms,ltype,zzz
  Use errors_warnings, Only : error,info

  Implicit None

  Integer,                        Save :: ncfzdn = 0

  Real( Kind = wp ), Allocatable, Save :: zdens(:,:)

  Public :: allocate_z_density_arrays

Contains

  Subroutine allocate_z_density_arrays()

    Integer :: fail

    fail = 0

    Allocate (zdens(1:mxgrdf,1:mxatyp), Stat = fail)

    If (fail > 0) Call error(1016)

    zdens = 0.0_wp

  End Subroutine allocate_z_density_arrays

  Subroutine z_density_collect()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for z-density profile
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


Subroutine z_density_compute(comm)

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

  Type( comms_type ), Intent( InOut ) :: comm
  Integer           :: j,k
  Real( Kind = wp ) :: zlen,delr,dvolz,factor,rho,rho1,rrr,sum,sum1
  Character( Len = 256 ) :: messages(2)

  Write(messages(1),'(a)') 'z density profiles:'
  Write(messages(2),'(2x,a,i8,a)') 'calculated using ',ncfzdn,' configurations'
  Call info(messages,2,.true.)

! open Z density file and write headers

  If (comm%idnode == 0) Then
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
     Write(messages(1),'(2x,a,a8)') 'rho(r): ',unqatm(k)
     Write(messages(2),'(9x,a1,6x,a3,9x,a4)') 'r','rho','n(r)'
     Call info(messages,2,.true.)
     If (comm%idnode == 0) Then
        Write(nzdndt,'(a8)') unqatm(k)
     End If

! global sum of data on all nodes

     Call gsum(comm,zdens(1:mxgrdf,k))

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

        Write(messages(1),'(2x,f10.4,1p,2e14.6)') rrr,rho1,sum1
        Call info(messages(1),.true.)
        If (comm%idnode == 0) Then
           Write(nzdndt,"(1p,2e14.6)") rrr,rho
        End If
     End Do
  End Do

  If (comm%idnode == 0) Close(Unit=nzdndt)

End Subroutine z_density_compute

  
End Module z_density
