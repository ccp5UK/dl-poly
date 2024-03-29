Module domains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring fundamental domain decomposition variables
  ! and arrays for the entire package
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2016
  ! contrib   - i.j.bush august 2010
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: comms_type,&
                             gsync
  Use errors_warnings, Only: error, error_alloc, error_dealloc
  Use kinds,           Only: wi,&
                             wp
  Use numerics,        Only: factor,&
                             get_nth_prime

  Implicit None

  Private

  !> Size of array for factorization. It allows for the first max_factor-1 prime numbers in the
  !> factorization, so 10 allows for the factors 2, 3, 5, 7, 11, 13, 17, 19, 23
  Integer, Parameter :: max_factor = 10

  Type, Public :: domains_type
    Private

    !> Dimensions of the 3D processor/domain grid
    Integer(Kind=wi), Public :: nx, ny, nz
    !> Real values of the dimensions nx,ny and nz
    Real(Kind=wp), Public    :: nx_real, ny_real, nz_real
    !> Reciprocal values of the dimensions domain%nx,domain%ny and domain%nz
    Real(Kind=wp), Public    :: nx_recip, ny_recip, nz_recip
    !> This domain's coordinates on the grid
    Integer(Kind=wi), Public :: idx, idy, idz
    !> Neighbourhood coordinates domain%map of the grid
    Integer(Kind=wi), Public :: map(1:26)
    !> Unique neighbour domain list
    !>
    !> - 0 if neighbour is unique
    !> - 1 if neighbour is repeated
    Integer(Kind=wi), Public :: map_unique(1:26)
    !> Number of neighbours (0 when serial, 26 otherwise). Used as the dimension
    !> of some shared unit arrays
    Integer(Kind=wi), Public :: neighbours
    Integer(Kind=wi), Public :: mxbfdp, mxbfxp, mxbfsh
  End Type domains_type

  Public :: map_domains, idcube, exchange_grid

Contains

  Subroutine map_domains(imcon, wx, wy, wz, domain, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for obtaining chatersoff%cutoffistics of a parallel
    ! computer and constructing a domain mapping
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & i.j.bush september 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent(In   ) :: imcon
    Real(Kind=wp),      Intent(In   ) :: wx, wy, wz
    Type(domains_type), Intent(InOut) :: domain
    Type(comms_type),   Intent(InOut) :: comm

    Real(Kind=wp), Parameter :: tol = 1.0e-6_wp

    Integer                          :: i, j, jdx, jdy, jdz, limx, limy, limz, nfacsx, nfacsy, nx, &
                                        ny, nz, P, Pyz
    Integer, Dimension(1:max_factor) :: pfacsx, pfacsy
    Real(Kind=wp)                    :: dx, dy, dz, min_S, S

!> Limits on the (x,y,z) sizes of the processor grid
!> processor factorisations of and factorised grid vectors over the x and y
!> sides of the Cartesian rectangular parallelepiped approximating the MD config%cell
!> Target grid product and running value of the (y,z) processor grid product
!> Running values of the (x,y,z) sizes of the processor producing P and Pyz
!> Running values of the (x,y,z) sizes of the grid config%cell volume and surface
!> Minimum surface
!> Search tolerance

    ! DD SEARCH
    If (comm%mxnode == 1) Then
      domain%nx = 1
      domain%ny = 1
      domain%nz = 1
    Else ! mxnode > 1

      ! Impose imcon driven limits on decomposition
      limx = Merge(Huge(1), 2, imcon /= 0)
      limy = Merge(Huge(1), 2, imcon /= 0)
      limz = Merge(Huge(1), 2, imcon /= 0 .and. imcon /= 6)

      ! Minimum surface (the largest possible) and decomposition (none)
      min_S = Huge(1.0_wp)
      domain%nx = -1
      domain%ny = -1
      domain%nz = -1

      P = comm%mxnode
      Call factor(P, pfacsx)
      nfacsx = get_n_factors(pfacsx)
      Do i = 1, nfacsx
        nx = get_nth_factor(pfacsx, i)
        If (nx > limx) Cycle
        dx = wx / Real(nx, wp)

        Pyz = P / nx
        Call factor(Pyz, pfacsy)
        nfacsy = get_n_factors(pfacsy)
        Do j = 1, nfacsy
          ny = get_nth_factor(pfacsy, j)
          If (ny > limy) Cycle
          dy = wy / Real(ny, wp)

          nz = Pyz / ny
          If (nz > limz) Cycle
          dz = wz / Real(nz, wp)

          S = 2.0_wp * (dx * dy + dy * dz + dz * dx)

          If (min_S - S > tol) Then

            ! Qualifier
            min_S = S
            domain%nx = nx
            domain%ny = ny
            domain%nz = nz

          Else If (Abs(min_S - S) < tol) Then

            ! For degenerate cases choose case where have least procs down any dimension (it should help FFT)
            If (Max(nx, ny, nz) < Max(domain%nx, domain%ny, domain%nz)) Then

              min_S = S
              domain%nx = nx
              domain%ny = ny
              domain%nz = nz

            Else If (Max(nx, ny, nz) == Max(domain%nx, domain%ny, domain%nz)) Then

              ! The new case has the same number of procs
              If (nx < domain%nx) Then

                ! Choose first the case which has the least down x
                min_S = S
                domain%nx = nx
                domain%ny = ny
                domain%nz = nz

              Else If (nx == domain%nx .and. ny < domain%ny) Then

                ! If max the same AND x the same choose it if y is less
                min_S = S
                domain%nx = nx
                domain%ny = ny
                domain%nz = nz

              End If
            End If
          End If
        End Do
      End Do

      Call gsync(comm)
      If (domain%nx == -1 .or. domain%ny == -1 .or. domain%nz == -1) Call error(520)

    End If

    ! DD MAPPING
    domain%map = 0
    domain%map_unique = 0

    domain%nx_real = Real(domain%nx, wp); domain%nx_recip = 1.0_wp / domain%nx_real
    domain%ny_real = Real(domain%ny, wp); domain%ny_recip = 1.0_wp / domain%ny_real
    domain%nz_real = Real(domain%nz, wp); domain%nz_recip = 1.0_wp / domain%nz_real

    ! construct domain%map of neighbouring nodes and domains
    domain%idz = comm%idnode / (domain%nx * domain%ny)
    domain%idy = comm%idnode / domain%nx - domain%idz * domain%ny
    domain%idx = Mod(comm%idnode, domain%nx)

    jdz = domain%nz + domain%idz
    jdy = domain%ny + domain%idy
    jdx = domain%nx + domain%idx

    domain%map(1) = idcube(Mod(jdx - 1, domain%nx), domain%idy, domain%idz, domain)
    domain%map(2) = idcube(Mod(domain%idx + 1, domain%nx), domain%idy, domain%idz, domain)
    domain%map(3) = idcube(domain%idx, Mod(jdy - 1, domain%ny), domain%idz, domain)
    domain%map(4) = idcube(domain%idx, Mod(domain%idy + 1, domain%ny), domain%idz, domain)
    domain%map(5) = idcube(domain%idx, domain%idy, Mod(jdz - 1, domain%nz), domain)
    domain%map(6) = idcube(domain%idx, domain%idy, Mod(domain%idz + 1, domain%nz), domain)

    !domain%map(1) points to the node number that is in -x direction to me (idnode)
    !domain%map(2) points to the node number that is in +x direction to me (idnode)
    !domain%map(3) points to the node number that is in -y direction to me (idnode)
    !domain%map(4) points to the node number that is in +y direction to me (idnode)
    !domain%map(5) points to the node number that is in -z direction to me (idnode)
    !domain%map(6) points to the node number that is in +z direction to me (idnode)

    domain%map(7) = idcube(Mod(jdx - 1, domain%nx), Mod(domain%idy + 1, domain%ny), domain%idz, domain)
    domain%map(8) = idcube(Mod(domain%idx + 1, domain%nx), Mod(jdy - 1, domain%ny), domain%idz, domain)
    domain%map(9) = idcube(Mod(jdx - 1, domain%nx), Mod(jdy - 1, domain%ny), domain%idz, domain)
    domain%map(10) = idcube(Mod(domain%idx + 1, domain%nx), Mod(domain%idy + 1, domain%ny), domain%idz, domain)

    domain%map(11) = idcube(Mod(jdx - 1, domain%nx), domain%idy, Mod(domain%idz + 1, domain%nz), domain)
    domain%map(12) = idcube(Mod(domain%idx + 1, domain%nx), domain%idy, Mod(jdz - 1, domain%nz), domain)
    domain%map(13) = idcube(Mod(jdx - 1, domain%nx), domain%idy, Mod(jdz - 1, domain%nz), domain)
    domain%map(14) = idcube(Mod(domain%idx + 1, domain%nx), domain%idy, Mod(domain%idz + 1, domain%nz), domain)

    domain%map(15) = idcube(domain%idx, Mod(jdy - 1, domain%ny), Mod(domain%idz + 1, domain%nz), domain)
    domain%map(16) = idcube(domain%idx, Mod(domain%idy + 1, domain%ny), Mod(jdz - 1, domain%nz), domain)
    domain%map(17) = idcube(domain%idx, Mod(jdy - 1, domain%ny), Mod(jdz - 1, domain%nz), domain)
    domain%map(18) = idcube(domain%idx, Mod(domain%idy + 1, domain%ny), Mod(domain%idz + 1, domain%nz), domain)

    domain%map(19) = idcube(Mod(jdx - 1, domain%nx), Mod(jdy - 1, domain%ny), Mod(jdz - 1, domain%nz), domain)
    domain%map(20) = idcube(Mod(domain%idx + 1, domain%nx), Mod(domain%idy + 1, domain%ny), Mod(domain%idz + 1, domain%nz), domain)
    domain%map(21) = idcube(Mod(jdx - 1, domain%nx), Mod(jdy - 1, domain%ny), Mod(domain%idz + 1, domain%nz), domain)
    domain%map(22) = idcube(Mod(domain%idx + 1, domain%nx), Mod(domain%idy + 1, domain%ny), Mod(jdz - 1, domain%nz), domain)

    domain%map(23) = idcube(Mod(jdx - 1, domain%nx), Mod(domain%idy + 1, domain%ny), Mod(jdz - 1, domain%nz), domain)
    domain%map(24) = idcube(Mod(domain%idx + 1, domain%nx), Mod(jdy - 1, domain%ny), Mod(domain%idz + 1, domain%nz), domain)
    domain%map(25) = idcube(Mod(jdx - 1, domain%nx), Mod(domain%idy + 1, domain%ny), Mod(domain%idz + 1, domain%nz), domain)
    domain%map(26) = idcube(Mod(domain%idx + 1, domain%nx), Mod(jdy - 1, domain%ny), Mod(jdz - 1, domain%nz), domain)

    ! Determine which processors appear more than once
    ! (domain%map_unique(the first unique node)=0 if repeated then =1)
    ! Up to Min(mxnode-1,26) elements of domain%map_unique can be zero
    ! the remainder is images of idnode
    !
    ! NEEDED FOR CATCHING SELF-HALOING
    Do i = 1, 26
      If (comm%idnode == domain%map(i)) domain%map_unique(i) = 1
      Do j = i + 1, 26
        If (domain%map(i) == domain%map(j)) domain%map_unique(j) = 1
      End Do
    End Do

  End Subroutine map_domains

  Function get_n_factors(factors) Result(nfacs)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 function to return the total number of
    ! possible integer factorisations
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush august 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Dimension(:), Intent(In   ) :: factors
    Integer                              :: nfacs

    nfacs = Product(factors(1:Size(factors) - 1) + 1)

    If (factors(Size(factors)) /= 1) nfacs = nfacs * 2
  End Function get_n_factors

  Function get_nth_factor(factors, n)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 function to return the n-th factoriser from the list of
    ! all possible integer factorisers
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush august 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Dimension(:), Intent(In   ) :: factors
    Integer,               Intent(In   ) :: n
    Integer                              :: get_nth_factor

    Integer                             :: dim_prod, i, nfacs, nt
    Integer, Dimension(1:Size(factors)) :: fac_counts

    nfacs = get_n_factors(factors)
    If (n > nfacs) Then
      get_nth_factor = -1
    Else
      If (factors(Size(factors)) /= 1 .and. n > nfacs / 2) Then
        nt = n - nfacs / 2
      Else
        nt = n
      End If
      nt = nt - 1

      dim_prod = Product(factors(1:Size(factors) - 2) + 1)
      Do i = Size(factors) - 1, 2, -1
        fac_counts(i) = nt / dim_prod
        nt = nt - fac_counts(i) * dim_prod
        dim_prod = dim_prod / (factors(i - 1) + 1)
      End Do
      fac_counts(1) = nt

      get_nth_factor = 1
      Do i = 1, Size(factors) - 1
        get_nth_factor = get_nth_factor * (get_nth_prime(i)**fac_counts(i))
      End Do

      If (factors(Size(factors)) /= 1 .and. n > nfacs / 2) Then
        get_nth_factor = get_nth_factor * factors(Size(factors))
      End If
    End If
  End Function get_nth_factor

  Function idcube(i, j, k, domain) Result(id)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 hypercube mapping function
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2006
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent(In   ) :: i, j, k
    Type(domains_type), Intent(In   ) :: domain
    Integer                           :: id

    id = i + domain%nx * (j + domain%ny * k)
  End Function idcube

  Subroutine exchange_grid(qqc_local, local_lb, local_ub, qqc_domain, domain_lb, domain_ub, domain, comm)

    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for exchanging grid data post DaFT
    !!
    !! Get all the data required to calculate the reciprocal space
    !! contribution to the forces on the atoms held by this processor, and
    !! send what is required by the other processors to them.
    !! QQC_LOCAL( IXB:IXT, IYB:IYT, IZB:IZT ) holds all the data I have currently.
    !! QQC_DOMAIN( IXDB:IXDT, IYDB:IYDT, IZDB:IZDT ) will hold all the data I need.
    !! MAP tells me about my neighbouring processors (see the DOMAINS_MODULE).
    !! MXSPL is the order of the spline.
    !! MXSPL1 is the extended spline should the conditional VNL is at work.
    !! MXSPL2 covers for MXSPL1 and if particles are not domain bound.
    !!
    !! copyright - daresbury laboratory
    !! author    - i.j.bush & i.t.todorov june 2014
    !! modified  - j.s.wilkins       september 2018
    !!-----------------------------------------------------------------------
    Use comms, Only: ExchgGrid_tag, comms_type, gsend, gwait, girecv
    Implicit None

    Integer, Dimension(3), Intent(In) :: local_lb, local_ub, domain_lb, domain_ub
    Complex(Kind=wp), Intent(In) :: qqc_local(local_lb(1):local_ub(1),local_lb(2):local_ub(2),local_lb(3):local_ub(3))
    Real(Kind=wp), Intent(Out) :: qqc_domain(domain_lb(1):domain_ub(1),domain_lb(2):domain_ub(2),domain_lb(3):domain_ub(3))
    Integer, Dimension(3) :: offset_lb, offset_ub
    Type(domains_type), Intent(In) :: domain
    Type(comms_type), Intent(InOut) :: comm

    Character(len=256) :: message

    Integer :: me

    ! What's my name?

    me = comm%idnode

    ! Offset -1 => no change
    offset_lb = local_lb - domain_lb - 1
    offset_ub = domain_ub - local_ub - 1

    ! Could strictly be legal?
    If (Any(offset_lb < -1) .or. Any(offset_ub < -1)) Then
      Write (message, '(a)') "Error: reducing grid size in exchanged_grid"
      Call error(0, message)
    End If

    ! Copy over our local data

    qqc_domain(local_lb(1):local_ub(1), local_lb(2):local_ub(2), local_lb(3):local_ub(3)) = Real(qqc_local, wp)

    If (Any(offset_lb >= 0)) Then

      ! Note that because of the way the splines work when particles don't
      ! blur off domains (ewld%bspline1==ewld%bspline), i.e. no conditional VNL updates,
      ! and are bound (ewld%bspline2==ewld%bspline) I only require to receive data from
      ! processors in the negative octant relative to me, and hence only
      ! need to send data to processors in the positive octant.

      ! +X direction face - negative halo

      Call exchange_grid_halo(domain%map(1), domain%map(2), &
                              local_ub(1) - offset_lb(1), local_ub(1), local_lb(2), local_ub(2), local_lb(3), local_ub(3), &
                              domain_lb(1), domain_lb(1) + offset_lb(1), local_lb(2), local_ub(2), local_lb(3), local_ub(3))

      ! +Y direction face (including the +X face extension) - negative halo

      Call exchange_grid_halo(domain%map(3), domain%map(4), &
                              domain_lb(1), local_ub(1), local_ub(2) - offset_lb(2), local_ub(2), local_lb(3), local_ub(3), &
                              domain_lb(1), local_ub(1), domain_lb(2), domain_lb(2) + offset_lb(2), local_lb(3), local_ub(3))

      ! +Z direction face (including the +Y+X faces extensions) - negative halo

      Call exchange_grid_halo(domain%map(5), domain%map(6), &
                              domain_lb(1), local_ub(1), domain_lb(2), local_ub(2), local_ub(3) - offset_lb(3), local_ub(3), &
                              domain_lb(1), local_ub(1), domain_lb(2), local_ub(2), domain_lb(3), domain_lb(3) + offset_lb(3))

    End If

    If (Any(offset_ub >= 0)) Then

      ! -X direction face - positive halo
      Call exchange_grid_halo(domain%map(2), domain%map(1), &
                              local_lb(1), local_lb(1) + offset_ub(1), domain_lb(2), local_ub(2), domain_lb(3), local_ub(3), &
                              domain_ub(1) - offset_ub(1), domain_ub(1), domain_lb(2), local_ub(2), domain_lb(3), local_ub(3))

      ! -Y direction face (including the +&-X faces extensions) - positive halo

      Call exchange_grid_halo(domain%map(4), domain%map(3), &
                              domain_lb(1), domain_ub(1), local_lb(2), local_lb(2) + offset_ub(2), domain_lb(3), local_ub(3), &
                              domain_lb(1), domain_ub(1), domain_ub(2) - offset_ub(2), domain_ub(2), domain_lb(3), local_ub(3))

      ! -Z direction face (including the +&-Y+&-X faces extensions) - positive halo

      Call exchange_grid_halo(domain%map(6), domain%map(5), &
                              domain_lb(1), domain_ub(1), domain_lb(2), domain_ub(2), local_lb(3), local_lb(3) + offset_ub(3), &
                              domain_lb(1), domain_ub(1), domain_lb(2), domain_ub(2), domain_ub(3) - offset_ub(3), domain_ub(3))

    End If

  Contains

    Subroutine exchange_grid_halo(from, to, &
                                  xlb, xlt, ylb, ylt, zlb, zlt, &
                                  xdb, xdt, ydb, ydt, zdb, zdt)

      !!-----------------------------------------------------------------------
      !!
      !! dl_poly_4 subroutine for exchanging grid data post DaFT
      !!
      !! Sends data to a processor      TO from
      !!                                QQC( XLB:XLT, YLB:YLT, ZLB:ZLT )
      !! Receives data from a processor FROM and sticks it in
      !!                                QQC( XDB:XDT, YDB:YDT, ZDB:ZDT )
      !!
      !! Note: Amount of data sent M-U-S-T be the same as that received!!!
      !!
      !! copyright - daresbury laboratory
      !! author    - i.j.bush & i.t.todorov june 2014
      !!
      !!-----------------------------------------------------------------------

    Integer, Intent(In   ) :: from, to, xlb, xlt, ylb, ylt, zlb, zlt, xdb, xdt, ydb, ydt, zdb, zdt

    Integer                                        :: fail(1:2)
    Real(Kind=wp), Allocatable, Dimension(:, :, :) :: recv_buffer, send_buffer

      ! If the processor to receive FROM is actually ME it means there is
      ! only one processor along this axis (so the processor to send TO is
      ! also ME) and so no message passing need be done.  However, there
      ! is a catch.  The domain decomposed spme_forces routine expects data
      ! that would be `seen' through the periodic boundary conditions to be
      ! actually copied from the high positive indices to negative ones and
      ! vice-versa.  The Else clause catches this.

      If (from /= me) Then

        ! Allocate send and receive buffers (of the same size!!!)
        ! so all can be sent and received as one message!!!

        fail = 0
        Allocate (send_buffer(xlb:xlt, ylb:ylt, zlb:zlt), Stat=fail(1))
        Allocate (recv_buffer(xdb:xdt, ydb:ydt, zdb:zdt), Stat=fail(2))
        If (Any(fail > 0)) Call error_alloc('buffers', 'exchange_grid_halo')

        ! Copy the data to be sent

        send_buffer = qqc_domain(xlb:xlt, ylb:ylt, zlb:zlt)

        ! Exchange the data

        Call girecv(comm, recv_buffer(:, :, :), from, ExchgGrid_tag)
        Call gsend(comm, send_buffer(:, :, :), to, ExchgGrid_tag)
        Call gwait(comm)

        ! Copy the received data into the domain halo

        qqc_domain(xdb:xdt, ydb:ydt, zdb:zdt) = recv_buffer

        ! And, as my mum told me, leave things as we found them

        Deallocate (recv_buffer, Stat=fail(1))
        Deallocate (send_buffer, Stat=fail(2))
        If (Any(fail > 0)) Call error_dealloc('buffers', 'exchange_grid_halo')

      Else

        ! Simple on node copy - as sizes are the same as 3D shapes

        qqc_domain(xdb:xdt, ydb:ydt, zdb:zdt) = qqc_domain(xlb:xlt, ylb:ylt, zlb:zlt)

      End If

    End Subroutine exchange_grid_halo

  End Subroutine exchange_grid

End Module domains
