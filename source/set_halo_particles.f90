Subroutine set_halo_particles(imcon,rcut,keyfce,lbook)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes
!
! copyright - daresbury laboratory
! amended   - i.t.todorov march 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,        Only : idnode
  Use setup_module
  Use domains_module
  Use site_module
  Use config_module
  Use rigid_bodies_module, Only : m_rgd,rgdxxx,rgdyyy,rgdzzz

  Implicit None

  Logical,           Intent( In    ) :: lbook
  Integer,           Intent( In    ) :: imcon,keyfce
  Real( Kind = wp ), Intent( In    ) :: rcut

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: cut

  Logical           :: oldjob
  Integer           :: fail,nlx,nly,nlz,i,j,ia,ib
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9), &
                       xdc,ydc,zdc,cwx,cwy,cwz,ecwx,ecwy,ecwz

  Real( Kind = wp ), Dimension( : ), Allocatable :: xxt,yyt,zzt

  If (newjob) Then
     newjob = .false.
     oldjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Define cut

     cut=rcut+1.0e-6_wp
  Else
     oldjob = .true.
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

! Larger widths may be needed by SPME for the b-splines -
! used in the halo transport in NEGATIVE DIRECTIONS ONLY!!!

  If (keyfce == 2) Then
     ecwx=Real(kmaxa,wp)/Real(mxspl,wp)
     ecwy=Real(kmaxb,wp)/Real(mxspl,wp)
     ecwz=Real(kmaxc,wp)/Real(mxspl,wp)

     ecwx=Min(cwx,ecwx)
     ecwy=Min(cwy,ecwy)
     ecwz=Min(cwz,ecwz)
  Else
     ecwx=cwx
     ecwy=cwy
     ecwz=cwz
  End If

! Distance from the - edge of this domain

  ecwx=Nearest( (-0.5_wp+ecwx)+Real(idx,wp)*r_nprx , +1.0_wp)+zero_plus
  ecwy=Nearest( (-0.5_wp+ecwy)+Real(idy,wp)*r_npry , +1.0_wp)+zero_plus
  ecwz=Nearest( (-0.5_wp+ecwz)+Real(idz,wp)*r_nprz , +1.0_wp)+zero_plus

! Distance from the + edge of this domain with a possible
! extension strip for the one liked cell per domain scenario

  cwx=Nearest( (-0.5_wp-cwx)+Real(idx+1,wp)*r_nprx , -1.0_wp)-zero_plus-Merge( cwx*1.0e-10_wp , 0.0_wp , nlx == 1 )
  cwy=Nearest( (-0.5_wp-cwy)+Real(idy+1,wp)*r_npry , -1.0_wp)-zero_plus-Merge( cwy*1.0e-10_wp , 0.0_wp , nly == 1 )
  cwz=Nearest( (-0.5_wp-cwz)+Real(idz+1,wp)*r_nprz , -1.0_wp)-zero_plus-Merge( cwz*1.0e-10_wp , 0.0_wp , nlz == 1 )

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions from MD cell centred
! Cartesian coordinates to reduced space ones
! Populate the halo indicator array

  nlast=natms                 ! No halo exists yet
  If (oldjob) ixyz(1:nlast)=0 ! Initialise halo indicator

  fail=0
  Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'set_halo_particles allocation failure, node: ', idnode
     Call error(0)
  End If

  Do i=1,nlast
     xxt(i)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     yyt(i)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     zzt(i)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

     If (xxt(i) <= ecwx) ixyz(i)=ixyz(i)+1
     If (xxt(i) >=  cwx) ixyz(i)=ixyz(i)+2

     If (yyt(i) <= ecwy) ixyz(i)=ixyz(i)+10
     If (yyt(i) >=  cwy) ixyz(i)=ixyz(i)+20

     If (zzt(i) <= ecwz) ixyz(i)=ixyz(i)+100
     If (zzt(i) >=  cwz) ixyz(i)=ixyz(i)+200
  End Do

  Deallocate (xxt,yyt,zzt, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'set_halo_particles deallocation failure, node: ', idnode
     Call error(0)
  End If

! exchange atom data in -/+ x directions

  Call export_atomic_data(-1)
  Call export_atomic_data( 1)

! exchange atom data in -/+ y directions

  Call export_atomic_data(-2)
  Call export_atomic_data( 2)

! exchange atom data in -/+ z directions

  Call export_atomic_data(-3)
  Call export_atomic_data( 3)

! assign incoming atom properties (of the halo only)

  Do i=natms+1,nlast
     ltype(i)=typsit(lsite(i))
     chge(i)=chgsit(lsite(i))
     weight(i)=wgtsit(lsite(i))
     lfrzn(i)=frzsit(lsite(i))
     lfree(i)=fresit(lsite(i))
  End Do

! Record global atom indices for local+halo sorting
! and sort multiple entries

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

End Subroutine set_halo_particles
