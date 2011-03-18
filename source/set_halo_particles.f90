Subroutine set_halo_particles(imcon,rcut,keyfce,lbook)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov march 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode
  Use setup_module
  Use domains_module,      Only : nprx,npry,nprz
  Use site_module
  Use config_module
  Use rigid_bodies_module, Only : m_rgd,rgdxxx,rgdyyy,rgdzzz

  Implicit None

  Logical,           Intent( In    ) :: lbook
  Integer,           Intent( In    ) :: imcon,keyfce
  Real( Kind = wp ), Intent( In    ) :: rcut

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: cut,sidex,sidey,sidez,dispx,dispy,dispz

  Logical           :: oldjob
  Integer           :: fail,nlx,nly,nlz,i,j,ia,ib,nlast_tmp1,nlast_tmp2
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9), &
                       cwx,cwy,cwz,ecwx,ecwy,ecwz,  &
                       uuu,vvv,www

  Real( Kind = wp ), Allocatable :: ott(:)

  Allocate (ott(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'set_halo_particles allocation failure, node: ', idnode
     Call error(0)
  End If

  If (newjob) Then
     newjob = .false.
     oldjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Define cut

     cut=rcut+1.0e-6_wp

! Get this node's (domain's) coordinates

     idz=idnode/(nprx*npry)
     idy=idnode/nprx-idz*npry
     idx=Mod(idnode,nprx)

! Get the domains' dimensions in reduced space
! (domains are geometrically equivalent)

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
  Else
     oldjob = .true.
  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! calculate link cell dimensions per node

  nlx=Int(sidex*celprp(7)/cut)
  nly=Int(sidey*celprp(8)/cut)
  nlz=Int(sidez*celprp(9)/cut)

! Calculate a link-cell width in every direction in the
! reduced space of the domain
! First term = the width of the domain in reduced space
! Second term = number of link-cells per domain per direction

  cwx=sidex/Real(nlx,wp)
  cwy=sidey/Real(nly,wp)
  cwz=sidez/Real(nlz,wp)

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

! Rebound coordinates just in case

  If (.not.oldjob) Then
     Do i=1,natms
        uuu=xxx(i)
        vvv=yyy(i)
        www=zzz(i)

        xxx(i)=rcell(1)*uuu+rcell(4)*vvv+rcell(7)*www+dispx
        yyy(i)=rcell(2)*uuu+rcell(5)*vvv+rcell(8)*www+dispy
        zzz(i)=rcell(3)*uuu+rcell(6)*vvv+rcell(9)*www+dispz

        uuu=xxx(i)-dispx
        vvv=yyy(i)-dispy
        www=zzz(i)-dispz

        xxx(i)=cell(1)*uuu+cell(4)*vvv+cell(7)*www
        yyy(i)=cell(2)*uuu+cell(5)*vvv+cell(8)*www
        zzz(i)=cell(3)*uuu+cell(6)*vvv+cell(9)*www
     End Do
  End If

! Convert atomic positions from MD cell centred Cartesian coordinates
! to reduced space coordinates of this node

  nlast=natms

  Do i=1,nlast
     uuu=xxx(i)
     vvv=yyy(i)
     www=zzz(i)

     xxx(i)=rcell(1)*uuu+rcell(4)*vvv+rcell(7)*www+dispx
     yyy(i)=rcell(2)*uuu+rcell(5)*vvv+rcell(8)*www+dispy
     zzz(i)=rcell(3)*uuu+rcell(6)*vvv+rcell(9)*www+dispz
  End Do

! exchange atom data in -/+ x directions

  If (nprx == 2 .and. nlx == 1) nlast_tmp1 = nlast ! Put tab on nlast
  Call export_atomic_data(-1,sidex,sidey,sidez,ecwx,ecwy,ecwz)
  If (nprx == 2 .and. nlx == 1) nlast_tmp2 = nlast ! Put tab on nlast
  Call export_atomic_data( 1,sidex,sidey,sidez,cwx,cwy,cwz)
  If (nprx == 2 .and. nlx == 1 .and. idx == 1) Then ! Handle exception
     Do i=nlast_tmp1+1,nlast
        ott(i)=xxx(i)
     End Do
     Do i=nlast_tmp1+1,nlast_tmp1+nlast-nlast_tmp2
        xxx(i)=ott(i+nlast_tmp2-nlast_tmp1)
     End Do
     Do i=nlast_tmp1+nlast-nlast_tmp2+1,nlast
        xxx(i)=ott(i+nlast_tmp2-nlast)
     End Do
  End If

! exchange atom data in -/+ y directions

  If (npry == 2 .and. nly == 1) nlast_tmp1 = nlast ! Put tab on nlast
  Call export_atomic_data(-2,sidex,sidey,sidez,ecwx,ecwy,ecwz)
  If (npry == 2 .and. nly == 1) nlast_tmp2 = nlast ! Put tab on nlast
  Call export_atomic_data( 2,sidex,sidey,sidez,cwx,cwy,cwz)
  If (npry == 2 .and. nly == 1 .and. idy == 1) Then ! Handle exception
     Do i=nlast_tmp1+1,nlast
        ott(i)=yyy(i)
     End Do
     Do i=nlast_tmp1+1,nlast_tmp1+nlast-nlast_tmp2
        yyy(i)=ott(i+nlast_tmp2-nlast_tmp1)
     End Do
     Do i=nlast_tmp1+nlast-nlast_tmp2+1,nlast
        yyy(i)=ott(i+nlast_tmp2-nlast)
     End Do
  End If

! exchange atom data in -/+ z directions

  If (nprz == 2 .and. nlz == 1) nlast_tmp1 = nlast ! Put tab on nlast
  Call export_atomic_data(-3,sidex,sidey,sidez,ecwx,ecwy,ecwz)
  If (nprz == 2 .and. nlz == 1) nlast_tmp2 = nlast ! Put tab on nlast
  Call export_atomic_data( 3,sidex,sidey,sidez,cwx,cwy,cwz)
  If (nprz == 2 .and. nlz == 1 .and. idz == 1) Then ! Handle exception
     Do i=nlast_tmp1+1,nlast
        ott(i)=zzz(i)
     End Do
     Do i=nlast_tmp1+1,nlast_tmp1+nlast-nlast_tmp2
        zzz(i)=ott(i+nlast_tmp2-nlast_tmp1)
     End Do
     Do i=nlast_tmp1+nlast-nlast_tmp2+1,nlast
        zzz(i)=ott(i+nlast_tmp2-nlast)
     End Do
  End If

! restore atomic coordinates to real coordinates

  Do i=1,nlast
     uuu=xxx(i)-dispx
     vvv=yyy(i)-dispy
     www=zzz(i)-dispz

     xxx(i)=cell(1)*uuu+cell(4)*vvv+cell(7)*www
     yyy(i)=cell(2)*uuu+cell(5)*vvv+cell(8)*www
     zzz(i)=cell(3)*uuu+cell(6)*vvv+cell(9)*www
  End Do

! assign incoming atom properties

  Do i=natms+1,nlast
     ltype(i)=typsit(lsite(i))
     chge(i)=chgsit(lsite(i))
     weight(i)=wgtsit(lsite(i))
     lfrzn(i)=frzsit(lsite(i))
     lfree(i)=fresit(lsite(i))
  End Do

! Refresh gtl or record global atom indices for local+halo sorting
! and sort multiple entries

  If (gtl_b > 0) Then
     Call get_gtl(lbook)
  Else
     Do i=1,nlast
        lsi(i)=i
        lsa(i)=ltg(i)
     End Do
     Call shellsort2(nlast,lsi,lsa)

     Do i=1,nlast-1
        j=1
        Do While ((i+j) <= nlast)
           If (lsa(i) == lsa(i+j)) Then
              ia=Min(lsi(i),lsi(i+j))
              ib=Max(lsi(i),lsi(i+j))
              lsi(i)=ia
              lsi(i+j)=ib
              j=j+1
           Else
              Exit
           End If
        End Do
     End Do
  End If

  nfree=0
  Do i=1,natms
     If (lfree(i) == 0) Then
        nfree=nfree+1
        lstfre(nfree)=i
     End If
  End Do

! Retag RBs when called again after the very first time
! when it's done in rigid_bodies_setup <- build_book_intra

  If (oldjob .and. m_rgd > 0) Then
     Call rigid_bodies_tags()
     Call rigid_bodies_coms(imcon,xxx,yyy,zzz,rgdxxx,rgdyyy,rgdzzz)
  End If

  Deallocate (ott, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'set_halo_particles deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine set_halo_particles
