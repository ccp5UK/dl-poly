!> Module containing neighbour list routines and variables
!>
!> Copyright - Daresbury Laboratory
!>
!> Author - J.Madge June 2018
!> Modified - A.B.G Chalk July 2018
!> Modified - i.t.todorov november 2019 (annotation)
!> Modified - A.B. Brukhno march-may 2020 - 'half-halo' VNL
!> Modified - m.a.seaton august 2020 - preprocessing tags for 'half-halo' VNL

Module neighbours

  Use comms,           Only: comms_type,&
                             gcheck,&
                             gmax,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: half_minus,&
                             half_plus, &
                             smalldr
  Use core_shell,      Only: core_shell_type
  Use development,     Only: development_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use kim,             Only: kim_type
  Use kinds,           Only: li,&
                             wi,&
                             wp
  Use mpole,           Only: POLARISATION_CHARMM,&
                             mpole_type
  Use numerics,        Only: dcell,&
                             images,&
                             invert,&
                             match,&
                             equal
  Use statistics,      Only: stats_type
  Use timer,           Only: start_timer,&
                             stop_timer,&
                             timer_type

  Implicit None

  Private

  !> Type containing neighbours data
  Type, Public :: neighbours_type

    Private

    ! Verlet neighbour list data

    !> Update cells flag
    Logical, Public                       :: update = .true.
    !> Unconditional update flag
    Logical, Public                       :: unconditional_update = .false.
    !> Tracking points for Verlet neighbour list
    Real(Kind=wp), Allocatable, Public :: xbg(:), ybg(:), zbg(:)

    !> Largest vdw cutoff, defines Verlet neighbour list radius
    Real(Kind=wp), Public :: cutoff = 0.0_wp
    !> Cutoff squared
    Real(Kind=wp), Public :: cutoff_2 = 0.0_wp
    !> Padding around cutoff
    Real(Kind=wp), Public                 :: padding = 0.0_wp
    !> Actual Verlet neighbour list cutoff (cutoff+padding)
    Real(Kind=wp), Public                 :: cutoff_extended = 0.0_wp
    !> Linked cell list
    Integer(Kind=wi), Allocatable, Public :: list(:, :)
    !> Maximum rank of linked cell list
    Integer(Kind=wi), Public              :: max_list = 0
    !> Maximum number of cells per domain
    Integer(Kind=wi), Public              :: max_cell = 0
    !> Excluded atom list
    Integer(Kind=wi), Allocatable, Public :: list_excl(:, :)
    !> Maximum size of excluded atom list
    Integer(Kind=wi), Public              :: max_exclude = 0
    !> CVNL newstart indicator
    Logical                               :: newstart = .true.
    !> CVNL newsjob indicator
    Logical                               :: newjob = .true.
    !> particle density at which VNL cubcelling stops
    Real(Kind=wp), Public                 :: pdplnc = 50.0_wp

  Contains

    Private

    Procedure, Public :: init_list
    Final             :: cleanup

  End Type neighbours_type

  Public :: vnl_check, vnl_set_check, link_cell_pairs, defects_link_cells

Contains

  !> Initialise the linked cell list
  Subroutine init_list(neigh, mxatdm)

    Class(neighbours_type)          :: neigh
    Integer(Kind=wi), Intent(In   ) :: mxatdm

    Allocate (neigh%list(-3:neigh%max_list, 1:mxatdm))
    Allocate (neigh%list_excl(0:neigh%max_exclude, 1:mxatdm))

    neigh%list(:, :) = 0
    neigh%list_excl(:, :) = 0

  End Subroutine init_list

  !> Perform Verlet neighbour list update check
  !>
  !> Copyright - Daresbury Laboratory
  !>
  !> Author    - I.T.Todorov january 2017
  !> contrib   - i.j.bush february 2014
  !> contrib   - i.t.todorov january 2017
  !> contrib   - i.t.todorov july 2018
  !> contrib   - i.t.todorov november 2019 (magic numbers commenting)

  Subroutine vnl_check(l_str, width, neigh, stat, domain, config, bspline, kim_data, comm)

    Logical,                  Intent(In   ) :: l_str
    Real(Kind=wp),            Intent(InOut) :: width
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(stats_type),         Intent(InOut) :: stat
    Type(domains_type),       Intent(In   ) :: domain
    Type(configuration_type), Intent(InOut) :: config
    Integer(Kind=wi),         Intent(In   ) :: bspline
    Type(kim_type),           Intent(In   ) :: kim_data
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=256)                       :: message
    Integer                                  :: fail, i, ilx, ily, ilz
    Real(Kind=wp)                            :: celprp(1:10), cut, m6, m7, m8, m9, test, tol
    Real(Kind=wp), Allocatable, Dimension(:) :: r, x, y, z

    If (.not. neigh%unconditional_update) Return

    Allocate (x(config%natms), y(config%natms), z(config%natms), r(config%natms), Stat=fail)

    ! Checks

    fail = 0
    If (fail > 0) Then
      Write (message, '(a)') 'vnl_check allocation failure'
      Call error(0, message)
    End If

    Do i = 1, config%natms
      x(i) = config%parts(i)%xxx - neigh%xbg(i)
      y(i) = config%parts(i)%yyy - neigh%ybg(i)
      z(i) = config%parts(i)%zzz - neigh%zbg(i)
    End Do

    Call images(config%imcon, config%cell, config%natms, x, y, z)

    If (config%natms > 0) Then
      r(1:config%natms) = Sqrt(x(1:config%natms)**2 + y(1:config%natms)**2 + z(1:config%natms)**2)
    End If

    ! search for violations: local domain then global

    tol = Merge(Maxval(r(1:config%natms)), 0.0_wp, config%natms > 0)

    If (fail > 0) Then
      Write (message, '(a)') 'vnl_check deallocation failure'
      Call error(0, message)
    End If

    Call gmax(comm, tol)

    ! decide on update if unsafe

    neigh%update = (tol >= half_minus * neigh%padding)

    ! Get the dimensional properties of the MD config%cell

    Call dcell(config%cell, celprp)
    width = Min(celprp(7), celprp(8), celprp(9))

    ! define cut as in link_cell_pairs

    cut = neigh%cutoff_extended + smalldr

    ! calculate link config%cell dimensions per node

    ilx = Int(domain%nx_recip * celprp(7) / cut)
    ily = Int(domain%ny_recip * celprp(8) / cut)
    ilz = Int(domain%nz_recip * celprp(9) / cut)

    m6 = 0.05_wp ! (Angstroms) is the limiting padding distance we bother with for CVNL
    m7 = 0.005_wp ! is the limiting fraction of cutoff (0.5%) we bother with for neigh%padding
    m8 = 0.02_wp ! 2% neigh%padding auto-push (may double)
    m9 = 0.95_wp ! (95%) giving 5% slack

    tol = Min(m6, m7 * neigh%cutoff) ! tolerance
    test = m8 * Merge(1.0_wp, 2.0_wp, bspline > 0) ! 2% (w/ SPME) or 4% (w/o SPME)
    cut = Min(domain%nx_recip * celprp(7), & ! domain size
              domain%ny_recip * celprp(8), &
              domain%nz_recip * celprp(9)) - smalldr ! remove delta r (the smallest distance we care about)

    If (ilx * ily * ilz == 0) Then
      If (cut < neigh%cutoff) Then
        Write (message, '(a)') 'neigh%cutoff <= Min(domain width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
        Call error(307, message, .true.)
      Else ! neigh%padding is defined & in 'no strict' mode
        If (cut < neigh%cutoff_extended) Then
          If (l_str) Then
            Write (message, '(a)') 'neigh%cutoff <= Min(domain width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
            Call error(307, message, .true.)
          Else
            If (cut >= neigh%cutoff) Then ! Re-set neigh%padding with some slack
              neigh%padding = Min(m9 * (cut - neigh%cutoff), test * neigh%cutoff)
              neigh%padding = Real(Int(100.0_wp * neigh%padding, wp)) / 100.0_wp
              If (neigh%padding < tol) neigh%padding = 0.0_wp ! Don't bother
              neigh%cutoff_extended = neigh%cutoff + neigh%padding
              neigh%update = .true.
            End If
          End If
        End If
      End If
    Else ! push the limits when up for update in a 'no strict' regime
      If (neigh%update .and. (.not. l_str)) Then ! Try to re-set neigh%padding with some slack
        If (Int(Real(Min(ilx, ily, ilz), wp) / (1.0_wp + test)) >= 2) Then
          cut = test * neigh%cutoff
        Else
          If (comm%mxnode > 1) Then
            cut = Min(m9 * (Min(domain%nx_recip * celprp(7) / Real(ilx, wp), &
                                domain%ny_recip * celprp(8) / Real(ily, wp), &
                                domain%nz_recip * celprp(9) / Real(ilz, wp)) &
                            - neigh%cutoff - smalldr), test * neigh%cutoff)
          Else ! catch & handle exception
            cut = m9 * (0.5_wp * width - neigh%cutoff - smalldr)
          End If
        End If
        cut = Real(Int(100.0_wp * cut), wp) / 100.0_wp
        If ((.not. (cut < tol)) .and. cut - neigh%padding > 0.005_wp) Then ! Do bother
          Write (message, '(2(a,f0.2),a)') 'cutoff padding reset from ', neigh%padding, &
            ' Angs to ', cut, ' Angs'
          Call info(message, .true.)
          neigh%padding = cut
          neigh%cutoff_extended = neigh%cutoff + neigh%padding
        End If
      End If
    End If

    ! Ensure padding is large enough for KIM model

    If (kim_data%padding_neighbours_required) Then
      If (neigh%padding < kim_data%influence_distance) Then
        neigh%padding = kim_data%influence_distance
        neigh%cutoff_extended = neigh%cutoff + neigh%padding
      End If
    End If

    If (neigh%update) Then ! Deal with skipping statistics
      stat%neighskip(3) = stat%neighskip(2) * stat%neighskip(3)
      stat%neighskip(2) = stat%neighskip(2) + 1.0_wp
      stat%neighskip(3) = stat%neighskip(3) / stat%neighskip(2) + stat%neighskip(1) / stat%neighskip(2)
      If (.not. neigh%newstart) Then ! avoid first compulsory force evaluation
        stat%neighskip(4) = Min(stat%neighskip(1), stat%neighskip(4))
      Else
        neigh%newstart = .false.
      End If
      stat%neighskip(5) = Max(stat%neighskip(1), stat%neighskip(5))

      stat%neighskip(1) = 0.0_wp ! Reset here, checkpoit set by vnl_set_check in set_halo_particles
    Else ! Enjoy telephoning
      stat%neighskip(1) = stat%neighskip(1) + 1.0_wp ! Increment, telephony done for xxx,yyy,zzz in set_halo_positions
    End If

    Deallocate (x, y, z, r, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vnl_check deallocation failure'
      Call error(0, message)
    End If

  End Subroutine vnl_check

  !> (Re)set the conditional Verlet neighbour list checkpoint -
  !> neigh%xbg,neigh%ybg,neigh%zbg at the end of set_halo_particles
  !>
  !> Copyright - Daresbury Laboratory
  !>
  !> Author    - I.T.Todorov january 2017

  Subroutine vnl_set_check(neigh, config)

    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(configuration_type), Intent(InOut) :: config

    Character(Len=256) :: message
    Integer            :: fail, i

    If (.not. neigh%unconditional_update) Return

    If (neigh%newjob) Then ! Init set
      neigh%newjob = .false.

      fail = 0
      Allocate (neigh%xbg(1:config%mxatms), neigh%ybg(1:config%mxatms), neigh%zbg(1:config%mxatms), Stat=fail)
      If (fail > 0) Then
        Write (message, '(a)') 'vnl_set_check allocation failure'
        Call error(0, message)
      End If

      ! CVNL state and skippage accumulators are initialised in vnl_module
      !
      !    neigh%update = .true.
      !    stat%neighskip(1) - cycles counter
      !    stat%neighskip(2) - access counter
      !    stat%neighskip(3) - average cycles
      !    stat%neighskip(4) - minimum cycles
      !    stat%neighskip(5) - maximum cycles
    End If

    ! set tracking point

    Do i = 1, config%nlast
      neigh%xbg(i) = config%parts(i)%xxx
      neigh%ybg(i) = config%parts(i)%yyy
      neigh%zbg(i) = config%parts(i)%zzz
    End Do

  End Subroutine vnl_set_check

  !> Verlet neighbour list based on link-cell method
  !>
  !> Copyright - Daresbury Laboratory
  !>
  !> Author    - I.T.Todorov january 2017
  !> Contrib   - I.J.Bush february 2014
  !> amended   - i.t.todorov november 2019 (nir deallocation fix, and comments)
  !> contrib   - a.v.brukhno march 2020 - 'half-halo' VNL search
  !> amended   - i.t.todorov may 2020 - flickering fix fixed
  !> contrib   - m.a.seaton august 2020 - 'half-halo' VNL search preprocessing tags

  Subroutine link_cell_pairs(rvdw, rmet, lbook, megfrz, cshell, devel, neigh, &
                             mpoles, domain, tmr, config, comm)

    Logical, Intent(In) :: lbook
    Integer, Intent(In) :: megfrz
    Real(Kind=wp), Intent(In) :: rvdw, rmet
    Type(core_shell_type), Intent(InOut) :: cshell
    Type(development_type), Intent(In) :: devel
    Type(neighbours_type), Intent(InOut) :: neigh
    Type(mpole_type), Intent(InOut) :: mpoles
    Type(domains_type), Intent(In) :: domain
    Type(timer_type), Intent(InOut) :: tmr
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type), Intent(InOut) :: comm

    Logical           :: safe, lx0, lx1, ly0, ly1, lz0, lz1

    Integer           :: fail(1:5), l_end, m_end, &
                         icell, ncells, ipass, &
                         kk, ll, ibig, i, ii, j, jj, j_start, &
                         nlx, nly, nlz, nlp, nlp2, nlp3, nlp4, nsbcll, &
                         nlx0s, nly0s, nlz0s, nlx0e, nly0e, nlz0e, &
                         nlx1s, nly1s, nlz1s, nlx1e, nly1e, nlz1e, &
                         ix, iy, iz, ic, ix1, ix2, iy1, iy2, iz1, iz2, &
                         jx, jy, jz, jc

    Real(Kind=wp) :: cut, rcsq, rsq, det, rcell(1:9), celprp(1:10), cnt(0:4), &
                     x, y, z, x1, y1, z1, dispx, dispy, dispz, xdc, ydc, zdc, nlr2

    Logical, Dimension(:), Allocatable :: nir
    Integer, Dimension(:), Allocatable :: nix, niy, niz, &
                                          lct_count, lct_start, &
                                          lct_where, which_cell, &
#ifndef HALF_HALO
                                          cell_dom, cell_bor, at_list
#else /* HALF_HALO */
                                          cell_dom, at_list ! cell_bor
#endif /* HALF_HALO */
    Real(Kind=wp), Dimension(:), Allocatable :: xxt, yyt, zzt

    Character(Len=256) :: message, messages(3)

#ifdef CHRONO
    Call start_timer(tmr, 'Linked Cells')
#endif
    Call dcell(config%cell, celprp)

    ! Get the dimensional properties of the MD config%cell

    det = Min(celprp(7), celprp(8), celprp(9))

    ! halt program if potential cutoff exceeds the minimum half-config%cell width

    If (neigh%cutoff_extended >= det / 2.0_wp) Then
      Call warning(3, neigh%cutoff_extended, det / 2.0_wp, 0.0_wp)
      Call error(95)
    End If

    ! Real space cutoff and squared r.s.c.

    cut = neigh%cutoff_extended + smalldr
    rcsq = neigh%cutoff_extended**2

    ! Calculate the number of link-cells per domain in every direction

    dispx = domain%nx_recip * celprp(7) / cut
    dispy = domain%ny_recip * celprp(8) / cut
    dispz = domain%nz_recip * celprp(9) / cut

    nlx = Int(dispx)
    nly = Int(dispy)
    nlz = Int(dispz)

    ! check for link cell algorithm violations

    If (nlx * nly * nlz == 0) Call error(307)

    ncells = (nlx + 2) * (nly + 2) * (nlz + 2)
    If (ncells > neigh%max_cell) Then
      Call warning(90, Real(ncells, wp), Real(neigh%max_cell, wp), 0.0_wp)
      neigh%max_cell = Nint(1.25_wp * Real(ncells, wp))
      If (ncells > config%mxatms) Call error(69)
    End If

    ! subcelling and new link-cell parameters
    nlp = 1
    nlr2 = Real(config%natms, wp)
#ifdef HALF_HALO
    Call gmax(comm, nlr2) ! Benefit the most loaded domain
#endif /* HALF_HALO */
    det = nlr2 / Real(nlx * nly * nlz, wp)
    Do While (det > neigh%pdplnc)
      nlp = nlp + 1
      rsq = Real(nlp, wp)
      nlx = Int(dispx * rsq)
      nly = Int(dispy * rsq)
      nlz = Int(dispz * rsq)
      det = nlr2 / Real(nlx * nly * nlz, wp)
    End Do
    ncells = (nlx + 2 * nlp) * (nly + 2 * nlp) * (nlz + 2 * nlp)
    nlp2 = (1 + (1 + 2 * nlp)**3) / 2 ! semi-ball size
    nlp3 = nlx * nly * nlz ! domain size
#ifndef HALF_HALO
    nlp4 = nlp3 - & ! border size
           Merge(nlx - 2 * nlp, 0, nlx - 2 * nlp > 0) * &
           Merge(nly - 2 * nlp, 0, nly - 2 * nlp > 0) * &
           Merge(nlz - 2 * nlp, 0, nlz - 2 * nlp > 0)
#endif /* HALF_HALO */

    fail = 0
    Allocate (nix(1:nlp2), niy(1:nlp2), niz(1:nlp2), nir(1:nlp2), Stat=fail(1))
    Allocate (which_cell(1:config%mxatms), at_list(1:config%mxatms), Stat=fail(2))
    Allocate (lct_count(0:ncells), lct_start(0:ncells + 1), lct_where(0:ncells + 1), Stat=fail(3))
#ifndef HALF_HALO
    Allocate (cell_dom(0:nlp3), cell_bor(0:nlp4), Stat=fail(4))
#else /* HALF_HALO */
    Allocate (cell_dom(0:nlp3), Stat=fail(4))
#endif /* HALF_HALO */
    Allocate (xxt(1:config%mxatms), yyt(1:config%mxatms), zzt(1:config%mxatms), Stat=fail(5))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'link_cell_pairs allocation failure'
      Call error(0, message)
    End If
    cell_dom(0) = nlp3 ! save array's limit
#ifndef HALF_HALO
    cell_bor(0) = nlp4 ! save array's limit
#endif /* HALF_HALO */

    ! Create jump-around arrays in link-cell space mapping a
    ! discrete 3D BALL so that ALL two body interactions are
    ! single counted within a domain (but not over-all as
    ! double counting still occurs globally for all shared
    ! inter-domain/semi-hello/cross-domain pairs.

    nix = 0; niy = 0; niz = 0; nir = .false.
    nlp2 = nlp**2
    nlp3 = (nlp - 1)**2
    nsbcll = 0
    Do iz = 0, nlp
      If (iz > 0) Then
        iz1 = (iz - 1)**2
      Else
        iz1 = 0
      End If

      jz = iz**2

      Do iy = -nlp, nlp
        If (iz == 0 .and. iy < 0) Cycle

        nlp4 = Abs(iy)
        If (nlp4 > 0) Then
          iy1 = (nlp4 - 1)**2
        Else
          iy1 = 0
        End If

        ll = iz1 + iy1
        If (ll > nlp2) Cycle

        jy = jz + iy**2

        Do ix = -nlp, nlp
          If (iz == 0 .and. iy == 0 .and. ix < 0) Cycle

          nlp4 = Abs(ix)
          If (nlp4 > 0) Then
            ix1 = (nlp4 - 1)**2
          Else
            ix1 = 0
          End If

          If (ll + ix1 > nlp2) Cycle

          jx = jy + ix**2

          nsbcll = nsbcll + 1

          nix(nsbcll) = ix
          niy(nsbcll) = iy
          niz(nsbcll) = iz
          nir(nsbcll) = (jx < nlp3)
        End Do
      End Do
    End Do
    !  Write(*,*) 'NLP',nlp,nsbcll,nlx,nly,nlz

    ! Get the total number of link-cells in MD config%cell per direction

    xdc = Real(nlx * domain%nx, wp)
    ydc = Real(nly * domain%ny, wp)
    zdc = Real(nlz * domain%nz, wp)

    ! Shifts from global to local link-config%cell space:
    ! (0,0,0) left-most link-config%cell on the domain (halo)
    ! (nlx+2*nlp-1,nly+2*nlp-1,nly+2*nlp-1) right-most
    ! link-cell on the domain (halo)

    jx = nlp - nlx * domain%idx
    jy = nlp - nly * domain%idy
    jz = nlp - nlz * domain%idz

    !***************************************************************
    ! Note(1): Due to numerical inaccuracy it is possible that some
    ! domain particles (1,config%natms) may have link-cell space
    ! coordinates in the halo / at least one coordinate as shown
    ! (nlx+nlp,nly+nlp,nlz+nlp)^(nlp-1,nlp-1,nlp-1) /
    ! as well as halo particles (config%natms+1,config%nlast) may have link-cell
    ! coordinates in the domain / all coordinate as shown
    ! (nlx,nly,nlz)^(nlp,nlp,nlp) / or even outside the
    ! standard one link-cell width, domain-surrounding halo
    ! / at least one coordinate as shown
    ! (>nlx+2*nlp-1,>nly+2*nlp-1,>nlz+2*nlp-1)^(<0,<0,<0) /.
    !
    ! Note(2): In SPME, at high accuracy of the ewald summation
    ! the b-splines may need more positive halo width than the
    ! standard one of one link-cell (in set_halo_particles it is
    ! ensured that such is supplied).  Such large positive halo may
    ! lead to EDGE EFFECTs - link-cells constituting the positive
    ! halo may have larger dimensions than the domain link-cells.
    !***************************************************************

    ! LC limits
    ! halo -start

    nlx0s = 0
    nly0s = 0
    nlz0s = 0

    ! halo -end

    nlx0e = nlp - 1
    nly0e = nlp - 1
    nlz0e = nlp - 1

    ! halo +start

    nlx1s = nlx + nlp
    nly1s = nly + nlp
    nlz1s = nlz + nlp

    ! halo +end

    nlx1e = nlx + 2 * nlp - 1
    nly1e = nly + 2 * nlp - 1
    nlz1e = nlz + 2 * nlp - 1

    ! Form linked list
    ! Initialise cell contents counter

    lct_count = 0

    ! Get the inverse cell matrix

    Call invert(config%cell, rcell, det)

    Do i = 1, config%natms

      ! Convert atomic positions from MD cell centred
      ! Cartesian coordinates to reduced space coordinates

      x = rcell(1) * config%parts(i)%xxx + rcell(4) * config%parts(i)%yyy + rcell(7) * config%parts(i)%zzz
      y = rcell(2) * config%parts(i)%xxx + rcell(5) * config%parts(i)%yyy + rcell(8) * config%parts(i)%zzz
      z = rcell(3) * config%parts(i)%xxx + rcell(6) * config%parts(i)%yyy + rcell(9) * config%parts(i)%zzz

      ! Get config%cell coordinates accordingly

      ix = Int(xdc * (x + 0.5_wp)) + jx
      iy = Int(ydc * (y + 0.5_wp)) + jy
      iz = Int(zdc * (z + 0.5_wp)) + jz

      ! Correction for domain (idnode) only particles (1,config%natms) but due to
      ! some tiny numerical inaccuracy kicked into its halo link-cell space
      ! Put all particles in a bounded link-cell space: lower and upper bounds
      ! as follows nl_coordinate_0e+1 <= i_coordinate <= nl_coordinate_1s-1 !

      ix = Max(Min(ix, nlx1s - 1), nlx0e+1)
      iy = Max(Min(iy, nly1s - 1), nly0e+1)
      iz = Max(Min(iz, nlz1s - 1), nlz0e+1)

      ! Hypercube function transformation (counting starts from one
      ! rather than zero /map_domains/ and 2*nlp more link-cells per
      ! dimension are accounted /coming from the halo/)

      icell = 1 + ix + (nlx + 2 * nlp) * (iy + (nly + 2 * nlp) * iz)

      ! count cell content

      lct_count(icell) = lct_count(icell) + 1

      ! backwards relationship

      which_cell(i) = icell

    End Do

    Do i = config%natms + 1, config%nlast

      ! Convert atomic positions from MD cell centred
      ! Cartesian coordinates to reduced space coordinates

      x = rcell(1) * config%parts(i)%xxx + rcell(4) * config%parts(i)%yyy + rcell(7) * config%parts(i)%zzz
      y = rcell(2) * config%parts(i)%xxx + rcell(5) * config%parts(i)%yyy + rcell(8) * config%parts(i)%zzz
      z = rcell(3) * config%parts(i)%xxx + rcell(6) * config%parts(i)%yyy + rcell(9) * config%parts(i)%zzz

      ! Get cell coordinates accordingly

      If (x > -half_plus) Then
        dispx = xdc * (x + 0.5_wp)
        ix = Int(dispx) + jx
      Else
        dispx = xdc * Abs(x + 0.5_wp)
        ix = -Int(dispx) + jx - 1
      End If
      If (y > -half_plus) Then
        dispy = ydc * (y + 0.5_wp)
        iy = Int(dispy) + jy
      Else
        dispy = ydc * Abs(y + 0.5_wp)
        iy = -Int(dispy) + jy - 1
      End If
      If (z > -half_plus) Then
        dispz = zdc * (z + 0.5_wp)
        iz = Int(dispz) + jz
      Else
        dispz = zdc * Abs(z + 0.5_wp)
        iz = -Int(dispz) + jz - 1
      End If

      ! Exclude any negatively bound residual halo

      If (ix >= nlx0s .and. iy >= nly0s .and. iz >= nlz0s) Then

        ! Correction for halo particles (config%natms+1,config%nlast) of this domain
        ! (idnode) but due to some tiny numerical inaccuracy kicked into
        ! the domain only link-cell space

        lx0 = (ix > nlx0e)
        lx1 = (ix < nlx1s)
        ly0 = (iy > nly0e)
        ly1 = (iy < nly1s)
        lz0 = (iz > nlz0e)
        lz1 = (iz < nlz1s)
        If ((lx0 .and. lx1) .and. &
            (ly0 .and. ly1) .and. &
            (lz0 .and. lz1)) Then

          ! Put on border in the halo

          x     = Abs( dispx - Real( nlx * domain%idx , wp ) )
          x1    = Abs( dispx - Real( nlx * (domain%idx+1) , wp ) )
          dispx = Min( x , x1 )
          y     = Abs( dispy - Real( nly * domain%idy , wp ) )
          y1    = Abs( dispy - Real( nly * (domain%idy+1) , wp ) )
          dispy = Min( y , y1 )
          z     = Abs( dispz - Real( nlz * domain%idz , wp ) )
          z1    = Abs( dispz - Real( nlz * (domain%idz+1) , wp ) )
          dispz = Min( z , z1 )

          If      (dispx <= dispy .and. dispx <= dispz) Then
            If (x < x1) Then
              ix=nlx0e
            Else
              ix=nlx1s
            End If
            If (equal(dispx,dispy)) Then
              If (y < y1) Then
                iy=nly0e
              Else
                iy=nly1s
              End If
            End If
            If (equal(dispx, dispz)) Then
              If (z < z1) Then
                iz=nlz0e
              Else
                iz=nlz1s
              End If
            End If
          Else If (dispy <= dispx .and. dispy <= dispz) Then
            If (y < y1) Then
              iy=nly0e
            Else
              iy=nly1s
            End If
            If (equal(dispy, dispz)) Then
              If (z < z1) Then
                iz=nlz0e
              Else
                iz=nlz1s
              End If
            End If
          Else
            If (z < z1) Then
              iz=nlz0e
            Else
              iz=nlz1s
            End If
          End If

        End If

        ! Check for positively bound residual halo
        lx0 = (ix < nlx0s)
        lx1 = (ix > nlx1e)
        ly0 = (iy < nly0s)
        ly1 = (iy > nly1e)
        lz0 = (iz < nlz0s)
        lz1 = (iz > nlz1e)
        If (.not. &
            (lx0 .or. lx1 .or. &
             ly0 .or. ly1 .or. &
             lz0 .or. lz1)) Then

          ! Hypercube function transformation (counting starts from one
          ! rather than zero /map_domains/ and 2*nlp more link-cells per
          ! dimension are accounted /coming from the halo/)

          icell = 1 + ix + (nlx + 2 * nlp) * (iy + (nly + 2 * nlp) * iz)
        Else

          ! Put possible residual halo in cell=0

          icell = 0

        End If

      Else

        ! Put possible residual halo in cell=0

        icell = 0

      End If

      ! count cell content

      lct_count(icell) = lct_count(icell) + 1

      ! backwards relationship

      which_cell(i) = icell

    End Do

    ! break down local list to list of linked-cell lists

    lct_start(0) = 1
    Do icell = 1, ncells + 1
      lct_start(icell) = lct_start(icell - 1) + lct_count(icell - 1)
    End Do

    ! domain local to linked-lists local mapping

    lct_where = lct_start
    Do i = 1, config%nlast
      !     at_list( lct_where( which_cell( i ) ) ) = i
      ! create a reordered coordinates arrays in the
      ! same manner to speeds up performance later
      j = lct_where(which_cell(i))
      at_list(j) = i

      xxt(j) = config%parts(i)%xxx
      yyt(j) = config%parts(i)%yyy
      zzt(j) = config%parts(i)%zzz

      lct_where(which_cell(i)) = lct_where(which_cell(i)) + 1
    End Do

    nlp3 = 0
    nlp4 = 0
    Do iz = nlz0e+1, nlz1s - 1
      iz1 = iz - nlz0e
      iz2 = iz - nlz1s
      Do iy = nly0e+1, nly1s - 1
        iy1 = iy - nly0e
        iy2 = iy - nly1s
        Do ix = nlx0e+1, nlx1s - 1
          ix1 = ix - nlx0e
          ix2 = ix - nlx1s

          ! index of the primary cell
          ic = 1 + ix + (nlx + 2 * nlp) * (iy + (nly + 2 * nlp) * iz)

          nlp3 = nlp3 + 1
          cell_dom(nlp3) = ic

#ifndef HALF_HALO
          ! loop over the domain's border config%cells only - ipass==2
          If ((ix1 >=  1 .and. ix1 <=  nlp) .or. &
              (ix2 <= -1 .and. ix2 >= -nlp) .or. &
              (iy1 >=  1 .and. iy1 <=  nlp) .or. &
              (iy2 <= -1 .and. iy2 >= -nlp) .or. &
              (iz1 >=  1 .and. iz1 <=  nlp) .or. &
              (iz2 <= -1 .and. iz2 >= -nlp)) Then
            nlp4 = nlp4 + 1
            cell_bor(nlp4) = ic
          End If
#endif /* HALF_HALO */
        End Do
      End Do
    End Do

    ! initialise Verlet neighbourlist (VNL) arrays

    !  neigh%list=0               ! (DEBUG)
#ifndef HALF_HALO
    neigh%list(-2:0, 1:config%natms) = 0 ! (COUNTING DIMENSIONS ONLY)
#else /* HALF_HALO */
    neigh%list(-3:0, 1:config%natms) = 0 ! (COUNTING DIMENSIONS ONLY)
#endif /* HALF_HALO */

    ! initial values of control variables

    ibig = 0
    safe = .true.

    ! primary loop over domain subcells
    ! loop over the domain's cells only

    ipass = 1
    Do icell = 1, cell_dom(0)
      ! index of the primary cell and its coordinates
      ic = cell_dom(icell)
      ix = Mod(ic - 1, nlx + 2 * nlp)
      iz = (ic - 1) / ((nlx + 2 * nlp) * (nly + 2 * nlp))
      iy = (ic - 1) / (nlx + 2 * nlp) - (nly + 2 * nlp) * iz

      ! loop over the primary cell contents

      Do ii = lct_start(ic), lct_start(ic + 1) - 1
        i = at_list(ii) ! get the particle index, by construction [1,config%natms]

        ! secondary loop over neighbouring cells
        ! ipass == 1 - non-negative non-repeatable semi-ball
        Do kk = ipass, nsbcll

          ! be on domain + possible positive halo cells only - ipass==1
          jx = ix + nix(kk)
          jy = iy + niy(kk)
          jz = iz + niz(kk)

          ! index of the secondary cell
          jc = 1 + jx + (nlx + 2 * nlp) * (jy + (nly + 2 * nlp) * jz)

          ! get the secondary list's starting particle index

          If (jc /= ic) Then

            ! if the secondary cell is different from the primary cell
            ! get the head of chain of the secondary cell

            j_start = lct_start(jc)

          Else ! only when ipass==1

            ! if the secondary cell is same as the primary cell
            ! get the next in line from the primary cell running index

            j_start = ii + 1

          End If

          ! bypass on real space cutoff check for safe cells when
          ! atom pairs' distances are guaranteed to be within the cutoff

          If (nir(kk)) Then

            ! check for overflow

            ll = neigh%list(0, i) + lct_start(jc + 1) - j_start
            If (ll <= neigh%max_list) Then

              ! loop over the secondary cell contents

              Do jj = j_start, lct_start(jc + 1) - 1

                ! get the particle index

                j = at_list(jj)

                ! add an entry

                ll = neigh%list(0, i) + 1
                neigh%list(ll, i) = j
                neigh%list(0, i) = ll

              End Do

              ! overflow is to occur

            Else

              ! loop over the secondary cell contents

              Do jj = j_start, lct_start(jc + 1) - 1

                ! get the particle index

                j = at_list(jj)

                ! check for overflow and add an entry

                ll = neigh%list(0, i) + 1
                If (ll <= neigh%max_list) Then
                  neigh%list(ll, i) = j
                Else
                  ibig = Max(ibig, ll)
                  safe = .false.
                End If
                neigh%list(0, i) = ll

                ! end of loop over the secondary cell contents

              End Do

            End If

            ! no bypass on real space cutoff check for safe cells when
            ! atom pairs' distances are not guaranteed to be within the cutoff
            ! distances in real space are needed for checking the cutoff criterion

          Else

            ! loop over the secondary cell contents

            Do jj = j_start, lct_start(jc + 1) - 1
              ! get the particle index

              j = at_list(jj)

              ! check cutoff criterion (all atom pairs MUST BE within the cutoff)
              !                 rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2
              ! use reordered coordinate list in order to get "ordered" rather than "random" access

              rsq = (xxt(jj) - config%parts(i)%xxx)**2 + (yyt(jj) - config%parts(i)%yyy)**2 + (zzt(jj) - config%parts(i)%zzz)**2
              If (rsq <= rcsq) Then

                ! check for overflow and add an entry

                ll = neigh%list(0, i) + 1
                If (ll <= neigh%max_list) Then
                  neigh%list(ll, i) = j
                Else
                  safe = .false.
                  ibig = Max(ibig, ll)
                End If
                neigh%list(0, i) = ll

                ! end of cutoff criterion check

              End If

              ! end of cutoff criterion check

            End Do

            ! end of bypass on real space cutoff check for safe cells

          End If

          ! secondary loop over neighbouring cells

        End Do

        ! end of loop over the primary cell contents

      End Do

      ! end of loops over domain cells

    End Do

#ifndef HALF_HALO
    ! secondary loop over domain's border cells only

    ipass = 2
    Do icell = 1, cell_bor(0)

      ! index of the primary cell and its coordinates

      ic = cell_bor(icell)
      ix = Mod(ic - 1, nlx + 2 * nlp)
      iz = (ic - 1) / ((nlx + 2 * nlp) * (nly + 2 * nlp))
      iy = (ic - 1) / (nlx + 2 * nlp) - (nly + 2 * nlp) * iz

      ! loop over the primary cell contents

      Do ii = lct_start(ic), lct_start(ic + 1) - 1

        i = at_list(ii) ! get the particle index, by construction [1,config%natms]

        ! secondary loop over neighbouring cells,
        ! when ipass==2 exclude self-self (i.e. domain - kk=1)

        Do kk = ipass, nsbcll

          ! describe a negative non-repeatable semi-ball - ipass==2

          jx = ix - nix(kk)
          jy = iy - niy(kk)
          jz = iz - niz(kk)

          ! be on halo cells only - ipass==2

          If ((jx <= nlx0e) .or. (jx >= nlx1s) .or. &
              (jy <= nly0e) .or. (jy >= nly1s) .or. &
              (jz <= nlz0e) .or. (jz >= nlz1s)) Then

            ! index of the secondary cell

            jc = 1 + jx + (nlx + 2 * nlp) * (jy + (nly + 2 * nlp) * jz)

            ! get the secondary list's starting particle index as
            ! ipass==2 the secondary cell is always different from the
            ! primary cell get the head of chain of the secondary cell

            If (jc /= ic) j_start = lct_start(jc)

            ! bypass on real space cutoff check for safe cells when
            ! atom pairs' distances are guaranteed to be within the cutoff

            If (nir(kk)) Then

              ! check for overflow

              ll = neigh%list(0, i) + lct_start(jc + 1) - j_start
              If (ll <= neigh%max_list) Then

                ! loop over the secondary cell contents

                Do jj = j_start, lct_start(jc + 1) - 1

                  ! get the particle index

                  j = at_list(jj)

                  ! add an entry

                  ll = neigh%list(0, i) + 1
                  neigh%list(ll, i) = j
                  neigh%list(0, i) = ll

                End Do

              Else

                ! loop over the secondary cell contents

                Do jj = j_start, lct_start(jc + 1) - 1

                  ! get the particle index

                  j = at_list(jj)

                  ! check for overflow and add an entry

                  ll = neigh%list(0, i) + 1
                  If (ll <= neigh%max_list) Then
                    neigh%list(ll, i) = j
                  Else
                    ibig = Max(ibig, ll)
                    safe = .false.
                  End If
                  neigh%list(0, i) = ll

                  ! end of loop over the secondary cell contents

                End Do

              End If

              ! no bypass on real space cutoff check for safe config%cells when
              ! atom pairs' distances are not guaranteed to be within the cutoff
              ! distances in real space are needed for checking the cutoff criterion

            Else

              ! loop over the secondary cell contents

              Do jj = j_start, lct_start(jc + 1) - 1

                ! get the particle index

                j = at_list(jj)

                ! check cutoff criterion (all atom pairs MUST BE within the cutoff)
                !                    rsq=(xxx(j)-xxx(i))**2+(yyy(j)-yyy(i))**2+(zzz(j)-zzz(i))**2
                ! use reordered coordinate list in order to get "ordered" rather than "random" access

                rsq = (xxt(jj) - config%parts(i)%xxx)**2 + (yyt(jj) - config%parts(i)%yyy)**2 + (zzt(jj) - config%parts(i)%zzz)**2
                If (rsq <= rcsq) Then
                  ! check for overflow and add an entry
                  ll = neigh%list(0, i) + 1
                  If (ll <= neigh%max_list) Then
                    neigh%list(ll, i) = j
                  Else
                    safe = .false.
                    ibig = Max(ibig, ll)
                  End If
                  neigh%list(0, i) = ll

                  ! end of cutoff criterion check

                End If

                ! end of loop over the secondary cell contents

              End Do

              ! end of bypass on real space cutoff check for safe cells

            End If

            ! end of inner if-block on cells and borders

          End If

          ! secondary loop over neighbouring cells

        End Do

        ! end of loop over the primary cell contents

      End Do

      ! end of loops over subconfig%cells only

    End Do
#endif /* HALF_HALO */

    ! terminate job if neighbour list array exceeded

    Call gcheck(comm, safe)
    If (.not. safe) Then
      Call gmax(comm, ibig)
      Call warning(290, Real(ibig, wp), Real(neigh%max_list, wp), 0.0_wp)
      Call error(106)
    End If

    ! Rear down frozen pairs

    If (megfrz > 1) Then
      Do i = 1, config%natms
        l_end = neigh%list(0, i)
        m_end = l_end

        ii = config%lfrzn(i)
        If (ii > 0) Then
          Do kk = l_end, 1, -1
            j = neigh%list(kk, i)
            jj = config%lfrzn(j)
            If (jj > 0) Then
              If (kk < m_end) Then
                neigh%list(kk, i) = neigh%list(m_end, i)
                neigh%list(m_end, i) = j
              End If
              m_end = m_end - 1
            End If
          End Do
        End If

        neigh%list(-2, i) = neigh%list(0, i) ! End of full non-repeatable half VNL (FNRH VNL)
        neigh%list(0, i) = m_end ! End of new neigh%list with no frozen pairs (NFP)
      End Do
    Else
      Do i = 1, config%natms
        neigh%list(-2, i) = neigh%list(0, i) ! End of full non-repeatable half VNL (FNRH VNL)
      End Do
    End If

    ! Rear down excluded pairs on top of the frozen ones

    If (lbook) Then
      Do i = 1, config%natms
        l_end = neigh%list(0, i)
        m_end = l_end

        ii = neigh%list_excl(0, i)
        If (ii > 0) Then
          Do kk = l_end, 1, -1
            j = neigh%list(kk, i)
            jj = config%ltg(j)
            If (match(jj, ii, neigh%list_excl(1:ii, i))) Then
              If (kk < m_end) Then
                neigh%list(kk, i) = neigh%list(m_end, i)
                neigh%list(m_end, i) = j
              End If
              m_end = m_end - 1
            End If
          End Do
        End If

        neigh%list(-1, i) = neigh%list(0, i) ! End of NFP FNRH VNL
        neigh%list(0, i) = m_end ! End of new list with no excluded interactions (NXI)
      End Do

      ! CHARMM core-shell screened electrostatic induction interactions
      ! Push up CHARMM pairs at the top of the bonded part of the neigh%list

      If (mpoles%key == POLARISATION_CHARMM) Then
        Do i = 1, config%natms
          l_end = neigh%list(-1, i) ! search marker to move up
          m_end = neigh%list(0, i) ! CHARMM marker to move down

          ii = mpoles%charmm(0, i)
          If (ii > 0) Then ! find what the local sublist CHARMM marker is
            outside: Do While (l_end > m_end + 1) ! Only when space for swap exists

              ! Check for space at the top

              j = neigh%list(m_end + 1, i)
              jj = config%ltg(j)
              If (match(jj, ii, mpoles%charmm(1:ii, i))) Then
                m_end = m_end + 1 ! move down CHARMM marker
                Cycle outside
              End If

              ! if a swap can be made at the top space m_end+1
              ! check for a qualifier (l_end) at the bottom

              inside: Do While (l_end > m_end + 1) ! Only when space for swap exists
                j = neigh%list(l_end, i)
                jj = config%ltg(j)
                If (match(jj, ii, mpoles%charmm(1:ii, i))) Then
                  ibig = neigh%list(m_end + 1, i)
                  neigh%list(m_end + 1, i) = neigh%list(l_end, i)
                  neigh%list(l_end, i) = ibig

                  l_end = l_end - 1 ! move up search marker
                  m_end = m_end + 1 ! move down CHARMM marker
                  Cycle outside
                End If
                l_end = l_end - 1 ! move up search marker
                Cycle inside
              End Do inside
            End Do outside
          End If
          neigh%list(-3, i) = m_end ! CHARMM end within NFP FNRH VNL, offset wrt NXI end
        End Do
      Else
        Do i = 1, config%natms
          neigh%list(-3, i) = neigh%list(0, i) ! CHARMM end coincides with NXI end
        End Do
      End If
    Else
      Do i = 1, config%natms
        neigh%list(-1, i) = neigh%list(0, i) ! End of NFP FNRH VNL
        neigh%list(-3, i) = neigh%list(0, i) ! CHARMM end coincides with NXI end
      End Do
    End If

    ! check on minimum separation distance between VNL pairs at re/start

    If (devel%l_dis) Then
      ! devel%l_dis=.false. ! at re/start ONLY

      cnt = 0.0_wp
      Do i = 1, config%natms
        ii = config%ltg(i)
        ll = cshell%legshl(0, i)

        ! iz=(which_cell(i)-1)/((nlx + 2*nlp)*(nlx + 2*nlp))
        ! iy=(which_cell(i)-1)/(nlx + 2*nlp) - (nly + 2*nlp)*iz
        ! ix=Mod(which_cell(i)-1,nlx + 2*nlp)

        Do kk = 1, neigh%list(-2, i)
          j = neigh%list(kk, i)
          jj = config%ltg(j)

          ! jz=(which_cell(j)-1)/((nlx + 2*nlp)*(nlx + 2*nlp))
          ! jy=(which_cell(j)-1)/(nlx + 2*nlp) - (nly + 2*nlp)*jz
          ! jx=Mod(which_cell(j)-1,nlx + 2*nlp)

          ! Exclude core-shell units' pairs from the check

          If (ll /= 0) Then ! the primary particle is part of a unit
            If (j <= config%natms) Then ! can check directly if the pair is part of the same unit
              If (cshell%legshl(0, j) /= 0) Then ! the secondary particle is part of a unit
                If (cshell%legshl(1, i) == cshell%legshl(1, j)) Cycle ! both are part of the same unit
              End If
            Else ! cannot check directly
              If (ll > 0) Then ! the primary particle is a core
                If (cshell%listshl(2, cshell%legshl(1, i)) == jj) Cycle ! jj is the shell of that core
              Else ! the primary particle is a shell
                If (cshell%listshl(1, cshell%legshl(1, i)) == jj) Cycle ! jj is the core of that shell
              End If
            End If
          End If

          If (j <= config%natms .or. ii < jj) Then
            cnt(1) = cnt(1) + 1.0_wp ! sum up all pairs (neigh%cutoff_extended=neigh%cutoff+neigh%padding)

            det = Sqrt((config%parts(i)%xxx - config%parts(j)%xxx)**2 + &
                       (config%parts(i)%yyy - config%parts(j)%yyy)**2 + &
                       (config%parts(i)%zzz - config%parts(j)%zzz)**2)

            If (det < devel%r_dis) Then
              safe = .false.
              Write (message, '(a,2(i10,a),f5.3,a)') &
                ' the pair with global indeces: ', ii, '  &', jj, &
                '  violates minimum separation distance (', det, ' Angs)'
              Call warning(message)
              cnt(0) = cnt(0) + 1.0_wp ! sum up violators
            End If

            If (kk <= neigh%list(0, i)) Then
              If (det < neigh%cutoff) cnt(2) = cnt(2) + 1.0_wp ! sum up all pairs (neigh%cutoff, electrostatics)
              If (det < rvdw) cnt(3) = cnt(3) + 1.0_wp ! sum up all pairs (rvdw, vdw)
              If (det <= rmet) cnt(4) = cnt(4) + 1.0_wp ! sum up all pairs (rmet, metal)
            End If
          End If
        End Do
      End Do

      Call gcheck(comm, safe, "enforce")
      Call gsum(comm, cnt)

      If (.not. safe) Then
        Write (message, '(i20,2a,f7.3,a)') &
          Int(cnt(0), li), ' pair(s) of particles in CONFIG ', &
          'violate(s) the minimum separation distance of ', devel%r_dis, ' Angs'
        Call warning(message, .true.)
      End If

      Call info('Pair totals of short range interactions over cutoffs (in Angstroms):', .true.)
      If (Abs(neigh%cutoff_extended - neigh%cutoff) > smalldr) Then
        Write (message, '(2x,a,i20,a,f7.3)') &
          'extended       -  ', Int(cnt(1), li), '  within neigh%cutoff_extended = ', neigh%cutoff_extended
        Call info(message, .true.)
      End If
      Write (messages(1), '(2x,a,i20,a,f7.3)') &
        'electrostatics -  ', Int(cnt(2), li), '  within r = ', neigh%cutoff
      Write (messages(2), '(2x,a,i20,a,f7.3)') &
        'van der Waals  -  ', Int(cnt(3), li), '  within r = ', rvdw
      Write (messages(3), '(2x,a,i20,a,f7.3)') &
        'metal          -  ', Int(cnt(4), li), '  within r = ', rmet
      Call info(messages, 3, .true.)
    End If

    Deallocate (nix, niy, niz, nir, Stat=fail(1))
    Deallocate (which_cell, at_list, Stat=fail(2))
    Deallocate (lct_count, lct_start, lct_where, Stat=fail(3))
#ifndef HALF_HALO
    Deallocate (cell_dom, cell_bor, Stat=fail(4))
#else /* HALF_HALO */
    Deallocate (cell_dom, Stat=fail(4))
#endif /* HALF_HALO */
    Deallocate (xxt, yyt, zzt, Stat=fail(5))
    If (Any(fail > 0)) Then
      Write (message, '(a)') 'link_cell_pairs deallocation failure'
      Call error(0, message)
    End If

#ifdef CHRONO
    Call stop_timer(tmr, 'Linked Cells')
#endif

  End Subroutine link_cell_pairs

  !> Verlet neighbour list based on link-cell method
  !>
  !> Copyright - Daresbury Laboratory
  !>
  !> Author    - I.T.Todorov january 2015
  !> amended   - i.t.todorov may 2020 - flickering fix fixed

  Subroutine defects_link_cells(config, cut, mxcldef, na, nl, xxt, yyt, zzt, nlx, nly, nlz, link, &
                                lct, domain)

    Type(configuration_type), Intent(In   ) :: config
    Real(Kind=wp),            Intent(In   ) :: cut
    Integer,                  Intent(In   ) :: mxcldef, na, nl
    Real(Kind=wp),            Intent(In   ) :: xxt(1:), yyt(1:), zzt(1:)
    Integer,                  Intent(  Out) :: nlx, nly, nlz, link(1:), lct(0:mxcldef)
    Type(domains_type),       Intent(In   ) :: domain

    Integer       :: i, icell, ix, iy, iz, jx, jy, jz, ncells
    Logical       :: lx0, lx1, ly0, ly1, lz0, lz1
    Real(Kind=wp) :: celprp(1:10), dispx, dispy, dispz, xdc, ydc, zdc, x,y,z, x1,y1,z1

    ! Get the dimensional properties of the MD cell

    Call dcell(config%cell, celprp)

    ! Calculate the number of link-cells per domain in every direction

    dispx = celprp(7) / (cut * domain%nx_real)
    dispy = celprp(8) / (cut * domain%ny_real)
    dispz = celprp(9) / (cut * domain%nz_real)

    nlx = Int(dispx)
    nly = Int(dispy)
    nlz = Int(dispz)

    ! check for link cell algorithm violations

    If (nlx * nly * nlz == 0) Call error(307)

    ncells = (nlx + 2) * (nly + 2) * (nlz + 2)
    If (ncells > mxcldef) Then
      Call warning(90, Real(ncells, wp), Real(mxcldef, wp), 0.0_wp)
      Call error(392)
    End If

    ! Get the total number of link-cells in MD cell per direction

    xdc = Real(nlx * domain%nx, wp)
    ydc = Real(nly * domain%ny, wp)
    zdc = Real(nlz * domain%nz, wp)

    ! Shifts from global to local link-cell space:
    ! (0,0,0) left-most link-cell on the domain (halo)
    ! (nlx+1,nly+1,nly+1) right-most
    ! link-cell on the domain (halo)

    jx = 1 - nlx * domain%idx
    jy = 1 - nly * domain%idy
    jz = 1 - nlz * domain%idz

    !***************************************************************
    ! Note(1): Due to numerical inaccuracy it is possible that some
    ! domain particles (1,na) may have link-cell space
    ! coordinates in the halo / at least one coordinate as shown
    ! (nlx+1,nly+1,nlz+1)^(0,0,0) / as well as halo particles
    ! (na+1,nl) may have link-cell coordinates in the domain
    ! / all coordinate as shown (nlx,nly,nlz)^(1,1,1) / or even
    ! outside the standard one link-cell width, domain-surrounding
    ! halo / at least one coordinate as shown
    ! (>nlx+1,>nly+1,>nlz+1)^(<0,<0,<0) /.
    !
    ! Note(2): In SPME, at high accuracy of the ewald summation
    ! the b-splines may need more positive halo width than the
    ! standard one of one link-cell (in set_halo_particles it is
    ! insured that such is supplied).  Such large positive halo may
    ! lead to EDGE EFFECTs - link-cells constituting the positive
    ! halo may have larger dimensions than the domain link-cells.
    !***************************************************************

    ! Form linked list
    ! Initialise link arrays

    Do i = 1, nl
      link(i) = 0
    End Do
    Do icell = 0, ncells
      lct(icell) = 0
    End Do

    Do i = nl, na + 1, -1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

      ! Get cell coordinates accordingly

      If (xxt(i) > -half_plus) Then
        dispx = xdc * (xxt(i) + 0.5_wp)
        ix = Int(dispx) + jx
      Else
        dispx = xdc * Abs(xxt(i) + 0.5_wp)
        ix = -Int(dispx) + jx - 1
      End If
      If (yyt(i) > -half_plus) Then
        dispy = ydc * (yyt(i) + 0.5_wp)
        iy = Int(dispy) + jy
      Else
        dispy = ydc * Abs(yyt(i) + 0.5_wp)
        iy = -Int(dispy) + jy - 1
      End If
      If (zzt(i) > -half_plus) Then
        dispz = zdc * (zzt(i) + 0.5_wp)
        iz = Int(dispz) + jz
      Else
        dispz = zdc * Abs(zzt(i) + 0.5_wp)
        iz = -Int(dispz) + jz - 1
      End If

      ! Correction for halo particles (na+1,nl) of this domain
      ! (idnode) but due to some tiny numerical inaccuracy kicked into
      ! the domain only link-cell space

      lx0 = (ix == 1)
      lx1 = (ix == nlx)
      ly0 = (iy == 1)
      ly1 = (iy == nly)
      lz0 = (iz == 1)
      lz1 = (iz == nlz)
      If ((lx0 .or. lx1) .and. &
          (ly0 .or. ly1) .and. &
          (lz0 .or. lz1)) Then
        x     = Abs( dispx - Real( nlx * domain%idx , wp ) )
        x1    = Abs( dispx - Real( nlx * (domain%idx+1) , wp ) )
        dispx = Min( x , x1 )
        y     = Abs( dispy - Real( nly * domain%idy , wp ) )
        y1    = Abs( dispy - Real( nly * (domain%idy+1) , wp ) )
        dispy = Min( y , y1 )
        z     = Abs( dispz - Real( nlz * domain%idz , wp ) )
        z1    = Abs( dispz - Real( nlz * (domain%idz+1) , wp ) )
        dispz = Min( z , z1 )

        If      (dispx <= dispy .and. dispx <= dispz) Then
          If (x < x1) Then
            ix=0
          Else
            ix=nlx+1
          End If
          If (equal(dispx, dispy)) Then
            If (y < y1) Then
              iy=0
            Else
              iy=nly+1
            End If
          End If
          If (equal(dispx, dispz)) Then
            If (z < z1) Then
              iz=0
            Else
              iz=nlz+1
            End If
          End If
        Else If (dispy <= dispx .and. dispy <= dispz) Then
          If (y < y1) Then
            iy=0
          Else
            iy=nly+1
          End If
          If (equal(dispy,dispz)) Then
            If (z < z1) Then
              iz=0
            Else
              iz=nlz+1
            End If
          End If
        Else
          If (z < z1) Then
            iz=0
          Else
            iz=nlz+1
          End If
        End If
      End If

      ! Check for residual halo

      lx0 = (ix < 0)
      lx1 = (ix > nlx + 1)
      ly0 = (iy < 0)
      ly1 = (iy > nly + 1)
      lz0 = (iz < 0)
      lz1 = (iz > nlz + 1)
      If (.not. &
          (lx0 .or. lx1 .or. &
           ly0 .or. ly1 .or. &
           lz0 .or. lz1)) Then

        ! Hypercube function transformation (counting starts from one
        ! rather than zero /map_domains/ and two more link-cells per
        ! dimension are accounted /coming from the halo/)

        icell = 1 + ix + (nlx + 2) * (iy + (nly + 2) * iz)

      Else

        ! Put possible residual halo in cell=0

        icell = 0

      End If

      ! link points to the next in chain or zero if end of chain occurs
      ! this is the old lct(icell)

      link(i) = lct(icell)

      ! at the end of the do-loop lct will point to the head of chain
      ! for this link-cell (update of lct(icell))

      lct(icell) = i

    End Do
    Do i = na, 1, -1 !!! BACKWARDS ORDER IS ESSENTIAL !!!

      ! Get cell coordinates accordingly

      ix = Int(xdc * (xxt(i) + 0.5_wp)) + jx
      iy = Int(ydc * (yyt(i) + 0.5_wp)) + jy
      iz = Int(zdc * (zzt(i) + 0.5_wp)) + jz

      ! Correction for domain (idnode) only particles (1,na) but due to
      ! some tiny numerical inaccuracy kicked into its halo link-cell space
      ! Put all particles in bounded link-cell space: lower and upper
      ! bounds as 1 <= i_coordinate <= nl_coordinate

      ix = Max(Min(ix, nlx), 1)
      iy = Max(Min(iy, nly), 1)
      iz = Max(Min(iz, nlz), 1)

      ! Hypercube function transformation (counting starts from one
      ! rather than zero /map_domains/ and two more link-cells per
      ! dimension are accounted /coming from the halo/)

      icell = 1 + ix + (nlx + 2) * (iy + (nly + 2) * iz)

      ! link points to the next in chain or zero if end of chain occurs
      ! this is the old lct(icell)

      link(i) = lct(icell)

      ! at the end of the do-loop lct will point to the head of chain
      ! for this link-cell (update of lct(icell))

      lct(icell) = i

    End Do

  End Subroutine defects_link_cells

  !> Dealloate neighbours type

  Subroutine cleanup(neigh)

    Type(neighbours_type) :: neigh

    If (Allocated(neigh%xbg)) Then
      Deallocate (neigh%xbg)
    End If
    If (Allocated(neigh%ybg)) Then
      Deallocate (neigh%ybg)
    End If
    If (Allocated(neigh%zbg)) Then
      Deallocate (neigh%zbg)
    End If

  End Subroutine cleanup

End Module neighbours
