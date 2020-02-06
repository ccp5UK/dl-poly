Module halo

  Use comms,           Only: comms_type,&
                             gcheck
  Use configuration,   Only: configuration_type
  Use constants,       Only: zero_plus
  Use deport_data,     Only: export_atomic_data,&
                             export_atomic_positions
  Use domains,         Only: domains_type
  Use electrostatic,   Only: ELECTROSTATIC_EWALD
  Use errors_warnings, Only: error
  Use ewald,           Only: ewald_type
  Use kim,             Only: kim_type
  Use kinds,           Only: wp
  Use mpole,           Only: mpole_type
  Use neighbours,      Only: neighbours_type,&
                             vnl_set_check
  Use numerics,        Only: dcell,&
                             invert,&
                             shellsort2
  Use site,            Only: site_type

  Implicit None

  Private

  Public :: refresh_halo_positions
  Public :: set_halo_particles

Contains

  Subroutine refresh_halo_positions(domain, config, kim_data, comm)

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
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(kim_type),           Intent(InOut) :: kim_data
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)   :: message
    Integer              :: fail, mlast
    Integer, Allocatable :: ixyz0(:)
    Logical              :: safe

    fail = 0
    Allocate (ixyz0(1:config%mxatms), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'refresh_halo_ppositions allocation failure'
      Call error(0, message)
    End If
    ixyz0(1:config%nlast) = config%ixyz(1:config%nlast)

    ! No halo, start with domain only particles

    mlast = config%natms

    ! exchange atom data in -/+ x directions

    Call export_atomic_positions(-1, mlast, ixyz0, domain, config, kim_data, comm)
    Call export_atomic_positions(1, mlast, ixyz0, domain, config, kim_data, comm)

    ! exchange atom data in -/+ y directions

    Call export_atomic_positions(-2, mlast, ixyz0, domain, config, kim_data, comm)
    Call export_atomic_positions(2, mlast, ixyz0, domain, config, kim_data, comm)

    ! exchange atom data in -/+ z directions

    Call export_atomic_positions(-3, mlast, ixyz0, domain, config, kim_data, comm)
    Call export_atomic_positions(3, mlast, ixyz0, domain, config, kim_data, comm)

    safe = (mlast == config%nlast)
    Call gcheck(comm, safe)
    If (.not. safe) Call error(138)

    Deallocate (ixyz0, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'referesh_halo_positions deallocation failure'
      Call error(0, message)
    End If
  End Subroutine refresh_halo_positions

  Subroutine set_halo_particles(electro_key, neigh, sites, mpoles, domain, config, ewld, kim_data, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange exchange of halo data between
    ! neighbouring domains/nodes
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                  Intent(In   ) :: electro_key
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(site_type),          Intent(In   ) :: sites
    Type(mpole_type),         Intent(InOut) :: mpoles
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Type(ewald_type),         Intent(In   ) :: ewld
    Type(kim_type),           Intent(InOut) :: kim_data
    Type(comms_type),         Intent(InOut) :: comm

    Integer       :: i, ia, ib, j, nlx, nly, nlz
    Real(Kind=wp) :: celprp(1:10), cut, cwx, cwy, cwz, det, ecwx, ecwy, ecwz, rcell(1:9), x, xdc, &
                     y, ydc, z, zdc

    ! Define cut

    cut = neigh%cutoff_extended + 1.0e-6_wp

    ! Get the dimensional properties of the MD cell

    Call dcell(config%cell, celprp)

    ! calculate link cell dimensions per node

    nlx = Int(celprp(7) / (cut * domain%nx_real))
    nly = Int(celprp(8) / (cut * domain%ny_real))
    nlz = Int(celprp(9) / (cut * domain%nz_real))

    ! Get the total number of link-cells in MD cell per direction

    xdc = Real(nlx * domain%nx, wp)
    ydc = Real(nly * domain%ny, wp)
    zdc = Real(nlz * domain%nz, wp)

    ! link-cell widths in reduced space

    cwx = 1.0_wp / xdc
    cwy = 1.0_wp / ydc
    cwz = 1.0_wp / zdc

    ! Larger widths may be needed by SPME for the b-splines -
    ! used in the halo transport in NEGATIVE DIRECTIONS ONLY!!!

    If (electro_key == ELECTROSTATIC_EWALD) Then
      ecwx=Real(ewld%bspline%num_spline_pad,wp)/Real(ewld%kspace%k_vec_dim(1),wp)
      ecwy=Real(ewld%bspline%num_spline_pad,wp)/Real(ewld%kspace%k_vec_dim(2),wp)
      ecwz=Real(ewld%bspline%num_spline_pad,wp)/Real(ewld%kspace%k_vec_dim(3),wp)

      ! I.e. take the smaller width in reduced space!!!

      ecwx = Max(cwx, ecwx)
      ecwy = Max(cwy, ecwy)
      ecwz = Max(cwz, ecwz)
    Else
      ecwx = cwx
      ecwy = cwy
      ecwz = cwz
    End If

    ! Distance from the - edge of this domain (larger positive halo)

    ecwx = Nearest((-0.5_wp + ecwx) + Real(domain%idx, wp) * domain%nx_recip, +1.0_wp) + zero_plus
    ecwy = Nearest((-0.5_wp + ecwy) + Real(domain%idy, wp) * domain%ny_recip, +1.0_wp) + zero_plus
    ecwz = Nearest((-0.5_wp + ecwz) + Real(domain%idz, wp) * domain%nz_recip, +1.0_wp) + zero_plus

    ! Distance from the + edge of this domain with a possible
    ! extension strip for the one linked cell per domain scenario

    cwx=Nearest( (-0.5_wp-cwx)+Real(domain%idx+1,wp)*domain%nx_recip , -1.0_wp)- &
      & zero_plus-Merge( cwx*1.0e-10_wp , 0.0_wp , nlx == 1 )
    cwy=Nearest( (-0.5_wp-cwy)+Real(domain%idy+1,wp)*domain%ny_recip , -1.0_wp)- &
      & zero_plus-Merge( cwy*1.0e-10_wp , 0.0_wp , nly == 1 )
    cwz=Nearest( (-0.5_wp-cwz)+Real(domain%idz+1,wp)*domain%nz_recip , -1.0_wp)- &
      & zero_plus-Merge( cwz*1.0e-10_wp , 0.0_wp , nlz == 1 )

    ! Get the inverse cell matrix

    Call invert(config%cell, rcell, det)

    ! Convert atomic positions from MD cell centred
    ! Cartesian coordinates to reduced space ones
    ! Populate the halo indicator array

    config%nlast = config%natms ! No halo exists yet
    config%ixyz(1:config%nlast) = 0 ! Initialise halo indicator
    Do i = 1, config%nlast
      x = rcell(1) * config%parts(i)%xxx + rcell(4) * config%parts(i)%yyy + rcell(7) * config%parts(i)%zzz
      y = rcell(2) * config%parts(i)%xxx + rcell(5) * config%parts(i)%yyy + rcell(8) * config%parts(i)%zzz
      z = rcell(3) * config%parts(i)%xxx + rcell(6) * config%parts(i)%yyy + rcell(9) * config%parts(i)%zzz

      If (x <= ecwx) config%ixyz(i) = config%ixyz(i) + 1
      If (x >= cwx) config%ixyz(i) = config%ixyz(i) + 2

      If (y <= ecwy) config%ixyz(i) = config%ixyz(i) + 10
      If (y >= cwy) config%ixyz(i) = config%ixyz(i) + 20

      If (z <= ecwz) config%ixyz(i) = config%ixyz(i) + 100
      If (z >= cwz) config%ixyz(i) = config%ixyz(i) + 200
    End Do

    ! exchange atom data in -/+ x directions

    Call export_atomic_data(-1, domain, config, kim_data, comm)
    Call export_atomic_data(1, domain, config, kim_data, comm)

    ! exchange atom data in -/+ y directions

    Call export_atomic_data(-2, domain, config, kim_data, comm)
    Call export_atomic_data(2, domain, config, kim_data, comm)

    ! exchange atom data in -/+ z directions

    Call export_atomic_data(-3, domain, config, kim_data, comm)
    Call export_atomic_data(3, domain, config, kim_data, comm)

    ! assign incoming atom properties (of the halo only)

    Do i = config%natms + 1, config%nlast
      config%ltype(i) = sites%type_site(config%lsite(i))
      config%parts(i)%chge = sites%charge_site(config%lsite(i))
      config%weight(i) = sites%weight_site(config%lsite(i))
      config%lfrzn(i) = sites%freeze_site(config%lsite(i))
      config%lfree(i) = sites%free_site(config%lsite(i))
    End Do

    ! Assign polarisation and dumping factor

    If (mpoles%max_mpoles > 0) Then
      Do i = config%natms + 1, config%nlast
        mpoles%polarisation_atom(i) = mpoles%polarisation_site(config%lsite(i))
        mpoles%dump_atom(i) = mpoles%dump_site(config%lsite(i))
      End Do
    End If

    ! Set VNL checkpoint

    Call vnl_set_check(neigh, config)

    ! Record global atom indices for local+halo sorting
    ! and sort multiple entries

    Do i = 1, config%nlast
      config%lsi(i) = i
      config%lsa(i) = config%ltg(i)
    End Do
    Call shellsort2(config%nlast, config%lsi, config%lsa)

    ! Sort multiple entries

    Do i = 1, config%nlast - 1
      j = 1
      Do While ((i + j) <= config%nlast)
        If (config%lsa(i) == config%lsa(i + j)) Then
          ia = Min(config%lsi(i), config%lsi(i + j))
          ib = Max(config%lsi(i), config%lsi(i + j))
          config%lsi(i) = ia
          config%lsi(i + j) = ib
          j = j + 1
        Else
          Exit
        End If
      End Do
    End Do

    config%nfree = 0
    Do i = 1, config%natms
      If (config%lfree(i) == 0) Then
        config%nfree = config%nfree + 1
        config%lstfre(config%nfree) = i
      End If
    End Do

  End Subroutine set_halo_particles
End Module halo
