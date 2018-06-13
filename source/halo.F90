Module halo

  Use comms,  Only : comms_type,gcheck
  Use deport_data, Only : export_atomic_positions, export_atomic_data
  Use setup,  Only : nrite,mxatms,kmaxa,kmaxb,kmaxc,mxspl1
  Use configuration 

  Use domains
  Use site, Only : site_type
  Use mpole

  Use neighbours,       Only : neighbours_type,vnl_set_check
  Use errors_warnings,  Only : error
  
  Implicit None


  Private
  Public :: refresh_halo_positions
  Public :: set_halo_particles
 
  Contains
  
  Subroutine refresh_halo_positions(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to refresh the halo positioning data between
! neighbouring domains/nodes when VNL is skipped
!
! Note: all depends on the ixyz halo array set in set_halo, this assumes
!       that (i) met%rcut=rcut! as well as (ii) all the error checks in there
!
! copyright - daresbury laboratory
! author    - i.t.todorov & i.j.bush february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Type( comms_type ), Intent(InOut) :: comm
  Logical :: safe
  Integer :: fail,mlast

  Integer, Allocatable :: ixyz0(:)
  Character ( Len = 256 )  :: message


  fail = 0
  Allocate (ixyz0(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'refresh_halo_ppositions allocation failure'
     Call error(0,message)
  End If
  ixyz0(1:nlast) = ixyz(1:nlast)

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call export_atomic_positions(-1,mlast,ixyz0,comm)
  Call export_atomic_positions( 1,mlast,ixyz0,comm)

! exchange atom data in -/+ y directions

  Call export_atomic_positions(-2,mlast,ixyz0,comm)
  Call export_atomic_positions( 2,mlast,ixyz0,comm)

! exchange atom data in -/+ z directions

  Call export_atomic_positions(-3,mlast,ixyz0,comm)
  Call export_atomic_positions( 3,mlast,ixyz0,comm)

  safe=(mlast == nlast)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(138)

  Deallocate (ixyz0, Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'referesh_halo_positions deallocation failure'
     Call error(0,message)
  End If

End Subroutine refresh_halo_positions


Subroutine set_halo_particles(keyfce,neigh,site_data,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of halo data between
! neighbouring domains/nodes
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: keyfce
  Type( neighbours_type ), Intent( InOut ) :: neigh
  Type( site_type ), Intent( In    ) :: site_data
  Type ( comms_type ), Intent( InOut  ) :: comm

  Real( Kind = wp ), Save :: cut

  Integer           :: nlx,nly,nlz,i,j,ia,ib
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9),x,y,z, &
                       xdc,ydc,zdc,cwx,cwy,cwz,ecwx,ecwy,ecwz

! Define cut

  cut=neigh%cutoff_extended+1.0e-6_wp

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
     ecwx=Real(mxspl1,wp)/Real(kmaxa,wp)
     ecwy=Real(mxspl1,wp)/Real(kmaxb,wp)
     ecwz=Real(mxspl1,wp)/Real(kmaxc,wp)

! I.e. take the smaller width in reduced space!!!

     ecwx=Max(cwx,ecwx)
     ecwy=Max(cwy,ecwy)
     ecwz=Max(cwz,ecwz)
  Else
     ecwx=cwx
     ecwy=cwy
     ecwz=cwz
  End If

! Distance from the - edge of this domain (larger positive halo)

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

  Call export_atomic_data(-1,comm)
  Call export_atomic_data( 1,comm)

! exchange atom data in -/+ y directions

  Call export_atomic_data(-2,comm)
  Call export_atomic_data( 2,comm)

! exchange atom data in -/+ z directions

  Call export_atomic_data(-3,comm)
  Call export_atomic_data( 3,comm)

! assign incoming atom properties (of the halo only)

  Do i=natms+1,nlast
     ltype(i)=site_data%type_site(lsite(i))
     chge(i)=site_data%charge_site(lsite(i))
     weight(i)=site_data%weight_site(lsite(i))
     lfrzn(i)=site_data%freeze_site(lsite(i))
     lfree(i)=site_data%free_site(lsite(i))
  End Do

! Assign polarisation and dumping factor

  If (mximpl > 0) Then
     Do i=natms+1,nlast
        plratm(i)=plrsit(lsite(i))
        dmpatm(i)=dmpsit(lsite(i))
     End Do
  End If

! Set VNL checkpoint

  Call vnl_set_check(neigh,comm)

! Record global atom indices for local+halo sorting
! and sort multiple entries

  Do i=1,nlast
     lsi(i)=i
     lsa(i)=ltg(i)
  End Do
  Call shellsort2(nlast,lsi,lsa)

! Sort multiple entries

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



  End Module halo
