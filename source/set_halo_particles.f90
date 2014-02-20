Subroutine set_halo_particles(imcon,rlnk,keyfce)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use domains_module
  Use site_module
  Use config_module

  Implicit None

  Integer,           Intent( In    ) :: imcon,keyfce
  Real( Kind = wp ), Intent( In    ) :: rlnk

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: cut

  Integer           :: nlx,nly,nlz,i,j,ia,ib
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9),x,y,z, &
                       xdc,ydc,zdc,cwx,cwy,cwz,ecwx,ecwy,ecwz

  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Define cut

     cut=rlnk+1.0e-6_wp
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
! extension strip for the one linked cell per domain scenario

  cwx=Nearest( (-0.5_wp-cwx)+Real(idx+1,wp)*r_nprx , -1.0_wp)-zero_plus-Merge( cwx*1.0e-10_wp , 0.0_wp , nlx == 1 )
  cwy=Nearest( (-0.5_wp-cwy)+Real(idy+1,wp)*r_npry , -1.0_wp)-zero_plus-Merge( cwy*1.0e-10_wp , 0.0_wp , nly == 1 )
  cwz=Nearest( (-0.5_wp-cwz)+Real(idz+1,wp)*r_nprz , -1.0_wp)-zero_plus-Merge( cwz*1.0e-10_wp , 0.0_wp , nlz == 1 )

! Get the inverse cell matrix

  Call invert(cell,rcell,det)

! Convert atomic positions from MD cell centred
! Cartesian coordinates to reduced space ones
! Populate the halo indicator array

  nlast=natms     ! No halo exists yet
  ixyz(1:nlast)=0 ! Initialise halo indicator
  Do i=1,nlast
     x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
     y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
     z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

     If (x <= ecwx) ixyz(i)=ixyz(i)+1
     If (x >=  cwx) ixyz(i)=ixyz(i)+2

     If (y <= ecwy) ixyz(i)=ixyz(i)+10
     If (y >=  cwy) ixyz(i)=ixyz(i)+20

     If (z <= ecwz) ixyz(i)=ixyz(i)+100
     If (z >=  cwz) ixyz(i)=ixyz(i)+200
  End Do

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

End Subroutine set_halo_particles
