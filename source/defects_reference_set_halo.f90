Subroutine defects_reference_set_halo  &
           (imcon,cut,cwx,cwy,cwz,dxl, &
           dxr,dyl,dyr,dzl,dzr,        &
           nrefs,nlrefs,namr,lri,lra,indr,xr,yr,zr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes for REFERENCE
!
! copyright - daresbury laboratory
! author    - i.t.todorov may 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode
  Use setup_module,   Only : mxatms
  Use config_module,  Only : cell
  Use domains_module, Only : nprx,npry,nprz

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
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: sidex,sidey,sidez,dispx,dispy,dispz

  Integer           :: i,j,ia,ib
  Real( Kind = wp ) :: celprp(1:10)


  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Get this node's (domain's) coordinates

     idz=idnode/(nprx*npry)
     idy=idnode/nprx-idz*npry
     idx=Mod(idnode,nprx)

! Get the domains' dimensions in reduced space
! (domains are geometrically equvalent)

     sidex=1.0_wp/Real(nprx,wp)
     sidey=1.0_wp/Real(npry,wp)
     sidez=1.0_wp/Real(nprz,wp)

! Calculate the displacements from the origin of the MD cell
! to the origin of this domain in reduced space

! First term (0.5_wp) = move to the bottom left corner of MD cell
! Second term, first term (side) = scale by the number of domains
! in the given direction
! Second term, second term, first term (id) = move to the bottom
! left corner of this domain in the given direction
! Second term, second term, second term (0.5_wp) = move to the
! middle of this domain

     dispx=0.5_wp-sidex*(Real(idx,wp)+0.5_wp)
     dispy=0.5_wp-sidey*(Real(idy,wp)+0.5_wp)
     dispz=0.5_wp-sidez*(Real(idz,wp)+0.5_wp)
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Calculate a link-cell(cut) width in every direction in the
! reduced space of the domain
! First term = the width of the domain in reduced space
! Second term = number of link-cells per domain per direction

  cwx=sidex/Real( Int(sidex*celprp(7)/cut),wp )
  cwy=sidey/Real( Int(sidey*celprp(8)/cut),wp )
  cwz=sidez/Real( Int(sidez*celprp(9)/cut),wp )

! Convert site positions from MD cell centred reduced space coordinates
! to reduced space coordinates of this node

  nlrefs=nrefs
  Do i=1,nlrefs
     xr(i)=xr(i)+dispx
     yr(i)=yr(i)+dispy
     zr(i)=zr(i)+dispz
  End Do

! exchange atom data in -/+ x directions

  Call defects_reference_export &
           (-1,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)
  Call defects_reference_export &
           ( 1,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)

! exchange atom data in -/+ y directions

  Call defects_reference_export &
           (-2,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)
  Call defects_reference_export &
           ( 2,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)

! exchange atom data in -/+ z directions

  Call defects_reference_export &
           (-3,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)
  Call defects_reference_export &
           ( 3,sidex,sidey,sidez,cwx,cwy,cwz,nlrefs,namr,indr,xr,yr,zr)

! restore site coordinates to MD centred reduced space coordinates

  Do i=1,nlrefs
     xr(i)=xr(i)-dispx
     yr(i)=yr(i)-dispy
     zr(i)=zr(i)-dispz
  End Do

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

  dxl=-dispx-sidex*0.5_wp-cwx ; dxr=-dispx+sidex*0.5_wp+cwx
  dyl=-dispy-sidey*0.5_wp-cwy ; dyr=-dispy+sidey*0.5_wp+cwy
  dzl=-dispz-sidez*0.5_wp-cwz ; dzr=-dispz+sidez*0.5_wp+cwz

End Subroutine defects_reference_set_halo
