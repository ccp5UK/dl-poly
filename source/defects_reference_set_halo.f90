Subroutine defects_reference_set_halo &
           (imcon,cut,cwx,cwy,cwz,    &
           dxl,dxr,dyl,dyr,dzl,dzr,   &
           nrefs,nlrefs,namr,lri,lra,indr,xr,yr,zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes for REFERENCE
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module,  Only : nrite,mxatms,zero_plus
  Use config_module, Only : cell
  Use domains_module

  Implicit None

  Integer,              Intent( In    ) :: imcon
  Real( Kind = wp ),    Intent( In    ) :: cut
  Real( Kind = wp ),    Intent(   Out ) :: cwx,cwy,cwz,dxl,dxr,dyl,dyr,dzl,dzr
  Integer,              Intent( In    ) :: nrefs
  Integer,              Intent(   Out ) :: nlrefs
  Character( Len = 8 ), Intent( InOut ) :: namr(1:mxatms)
  Integer,              Intent( InOut ) :: lri(1:mxatms),lra(1:mxatms),indr(1:mxatms)
  Real( Kind = wp ),    Intent( InOut ) :: xr(1:mxatms),yr(1:mxatms),zr(1:mxatms)

  Logical,           Save :: newjob = .true.

  Integer           :: fail,nlx,nly,nlz,i,j,ia,ib
  Real( Kind = wp ) :: celprp(1:10),xdc,ydc,zdc

  Integer, Allocatable :: ixyz(:)

  fail=0
  Allocate (ixyz(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_set_halo allocation failure, node: ', idnode
     Call error(0)
  End If

  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! calculate link cell dimensions per node

  nlx=Int(celprp(7)/(cut*nprx_r))
  nly=Int(celprp(8)/(cut*npry_r))
  nlz=Int(celprp(9)/(cut*nprz_r))

! Get the total number of link-cells in MD cell per direction

  xdc=Real(nlx*nprx,wp)
  ydc=Real(nly*npry,wp)
  zdc=Real(nlz*nprz,wp)

! link-cell widths in reduced space

  cwx=1.0_wp/xdc
  cwy=1.0_wp/ydc
  cwz=1.0_wp/zdc

! Distance from the - edge of this domain

  dxl=Nearest( (-0.5_wp+cwx)+Real(idx,wp)*r_nprx , +1.0_wp)+zero_plus
  dyl=Nearest( (-0.5_wp+cwy)+Real(idy,wp)*r_npry , +1.0_wp)+zero_plus
  dzl=Nearest( (-0.5_wp+cwz)+Real(idz,wp)*r_nprz , +1.0_wp)+zero_plus

! Distance from the + edge of this domain

  dxr=Nearest( (-0.5_wp-cwx)+Real(idx+1,wp)*r_nprx , -1.0_wp)-zero_plus
  dyr=Nearest( (-0.5_wp-cwy)+Real(idy+1,wp)*r_npry , -1.0_wp)-zero_plus
  dzr=Nearest( (-0.5_wp-cwz)+Real(idz+1,wp)*r_nprz , -1.0_wp)-zero_plus

! Convert atomic positions from MD cell centred
! Cartesian coordinates to reduced space ones
! Populate the halo indicator array

  nlrefs=nrefs     ! No halo exists yet
  ixyz(1:nlrefs)=0 ! Initialise halo indicator
  Do i=1,nlrefs
     If (xr(i) <= dxl) ixyz(i)=ixyz(i)+1
     If (xr(i) >= dxr) ixyz(i)=ixyz(i)+2

     If (yr(i) <= dyl) ixyz(i)=ixyz(i)+10
     If (yr(i) >= dyr) ixyz(i)=ixyz(i)+20

     If (zr(i) <= dzl) ixyz(i)=ixyz(i)+100
     If (zr(i) >= dzr) ixyz(i)=ixyz(i)+200
  End Do

! exchange atom data in -/+ x directions

  Call defects_reference_export(-1,ixyz,nlrefs,namr,indr,xr,yr,zr)
  Call defects_reference_export( 1,ixyz,nlrefs,namr,indr,xr,yr,zr)

! exchange atom data in -/+ y directions

  Call defects_reference_export(-2,ixyz,nlrefs,namr,indr,xr,yr,zr)
  Call defects_reference_export( 2,ixyz,nlrefs,namr,indr,xr,yr,zr)

! exchange atom data in -/+ z directions

  Call defects_reference_export(-3,ixyz,nlrefs,namr,indr,xr,yr,zr)
  Call defects_reference_export( 3,ixyz,nlrefs,namr,indr,xr,yr,zr)

  Do i=1,nlrefs
     lri(i)=i
     lra(i)=indr(i)
  End Do
  Call shellsort2(nlrefs,lri,lra)

  Do i=1,nlrefs-1
     j=1
     Do While ((i+j) <= nlrefs)
        If (lra(i) == lra(i+j)) Then
           ia=Min(lri(i),lri(i+j))
           ib=Max(lri(i),lri(i+j))
           lri(i)=ia
           lri(i+j)=ib
           j=j+1
        Else
           Exit
        End If
     End Do
  End Do

! Get domain halo limits in reduced space

! Distance from the - edge of this domain

  dxl=Nearest( (-0.5_wp-cwx)+Real(idx,wp)*r_nprx , -1.0_wp)-zero_plus
  dyl=Nearest( (-0.5_wp-cwy)+Real(idy,wp)*r_npry , -1.0_wp)-zero_plus
  dzl=Nearest( (-0.5_wp-cwz)+Real(idz,wp)*r_nprz , -1.0_wp)-zero_plus

! Distance from the + edge of this domain

  dxr=Nearest( (-0.5_wp+cwx)+Real(idx+1,wp)*r_nprx , +1.0_wp)+zero_plus
  dyr=Nearest( (-0.5_wp+cwy)+Real(idy+1,wp)*r_npry , +1.0_wp)+zero_plus
  dzr=Nearest( (-0.5_wp+cwz)+Real(idz+1,wp)*r_nprz , +1.0_wp)+zero_plus

  Deallocate (ixyz, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'defects_reference_set_halo deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine defects_reference_set_halo
