Subroutine metal_ld_set_halo(imcon,rmet,keyfce,rho)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of density data between
! neighbouring domains/nodes
!
! copyright - daresbury laboratory
! author    - w.smith  april 1999
! amended   - i.t.todorov october 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gcheck
  Use setup_module
  Use domains_module, Only : nprx,npry,nprz
  Use config_module,  Only : cell,natms,nlast,ltg,xxx,yyy,zzz

  Implicit None

  Integer,                                  Intent( In    ) :: imcon,keyfce
  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( InOut ) :: rho

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: cut,sidex,sidey,sidez,dispx,dispy,dispz

  Logical           :: safe
  Integer           :: fail,i,mlast
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9), &
                       cwx,cwy,cwz,ecwx,ecwy,ecwz,  &
                       uuu,vvv,www

  Integer, Dimension( : ), Allocatable :: iwrk

  fail=0
  Allocate (iwrk(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_set_halo allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Define cut

     cut=rmet+1.0e-6_wp

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

! Calculate a link-cell width in every direction in the
! reduced space of the domain
! First term = the width of the domain in reduced space
! Second term = number of link-cells per domain per direction

  cwx=sidex/Real( Int(sidex*celprp(7)/cut),wp )
  cwy=sidey/Real( Int(sidey*celprp(8)/cut),wp )
  cwz=sidez/Real( Int(sidez*celprp(9)/cut),wp )

! "Positive halo" widths in reduced space as needed by SPME for
! b-splines. To be used in halo transport in NEGATIVE DIRECTIONS!!!

  If (keyfce == 2) Then
     ecwx=Real(mxspl,wp)/Real(kmaxa,wp)
     ecwy=Real(mxspl,wp)/Real(kmaxb,wp)
     ecwz=Real(mxspl,wp)/Real(kmaxc,wp)

     ecwx=Max(cwx,ecwx)
     ecwy=Max(cwy,ecwy)
     ecwz=Max(cwz,ecwz)
  Else
     ecwx=cwx
     ecwy=cwy
     ecwz=cwz
  End If

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions from MD cell centred Cartesian coordinates
! to reduced space coordinates of this node

  Do i=1,nlast
     uuu=xxx(i)
     vvv=yyy(i)
     www=zzz(i)

     xxx(i)=rcell(1)*uuu+rcell(4)*vvv+rcell(7)*www+dispx
     yyy(i)=rcell(2)*uuu+rcell(5)*vvv+rcell(8)*www+dispy
     zzz(i)=rcell(3)*uuu+rcell(6)*vvv+rcell(9)*www+dispz

     iwrk(i)=ltg(i)
  End Do

! save the number of atoms here

  mlast=natms

! exchange atom data in -/+ x directions

  Call metal_ld_export(-1,sidex,sidey,sidez,ecwx,ecwy,ecwz,mlast,iwrk,rho)
  Call metal_ld_export( 1,sidex,sidey,sidez,cwx,cwy,cwz,mlast,iwrk,rho)

! exchange atom data in -/+ y directions

  Call metal_ld_export(-2,sidex,sidey,sidez,ecwx,ecwy,ecwz,mlast,iwrk,rho)
  Call metal_ld_export( 2,sidex,sidey,sidez,cwx,cwy,cwz,mlast,iwrk,rho)

! exchange atom data in -/+ z directions

  Call metal_ld_export(-3,sidex,sidey,sidez,ecwx,ecwy,ecwz,mlast,iwrk,rho)
  Call metal_ld_export( 3,sidex,sidey,sidez,cwx,cwy,cwz,mlast,iwrk,rho)

! check atom totals after data transfer

  safe=(mlast == nlast)
  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(96)

! check incoming atomic density assignments

  Do i=natms+1,nlast
     safe=(ltg(i) == iwrk(i))
  End Do

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe) Call error(98)

! restore atomic coordinates to real coordinates

  Do i=1,nlast
     uuu=xxx(i)-dispx
     vvv=yyy(i)-dispy
     www=zzz(i)-dispz

     xxx(i)=cell(1)*uuu+cell(4)*vvv+cell(7)*www
     yyy(i)=cell(2)*uuu+cell(5)*vvv+cell(8)*www
     zzz(i)=cell(3)*uuu+cell(6)*vvv+cell(9)*www
  End Do

  Deallocate (iwrk, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'metal_ld_set_halo deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine metal_ld_set_halo
