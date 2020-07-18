Module kspace
  !!----------------------------------------------------------------------!
  !!
  !! dl_poly_4 module for containing types and functions
  !! relating to ewald kspace
  !!
  !! copyright - daresbury laboratory
  !! author    - j.s.wilkins october 2018
  !!
  !!----------------------------------------------------------------------!
  Use comms,           Only: comms_type
  Use domains,         Only: domains_type
  Use errors_warnings, Only: error_alloc
  Use kinds,           Only: wi,&
                             wp
  Use parallel_fft,    Only: initialize_fft,&
                             pfft_indices

  Implicit None

  Private
  Public :: setup_kspace

  Type, Public :: kspace_type
    Private

    !> SPME kspace    dimensions
    Integer(Kind=wi), Dimension(3), Public :: k_vec_dim, k_vec_dim_cont
    !> Largest Integer K-Vector Index
    Integer(Kind=wi), Public :: k_vec_max

    !> Real K-Vector indices
    Real(Kind=wp), Dimension(3), Public :: k_vec_dim_real
    !> Real K-Vector indices
    Real(Kind=wp), Dimension(3), Public :: k_vec_dim_real_p_dom
    !> Largest K-Vector index
    Real(Kind=wp), Public :: k_vec_max_real

    !> blocking factors for splines and fft
    Integer, Dimension(3), Public :: block_fac
    !> indexing arrays for x, y & z as used in parallel fft
    Integer, Dimension(:, :), Allocatable, Public :: index
    !> context for parallel fft
    Integer, Public :: context

    !> Domain num cells
    Integer, Dimension(3), Public :: domain_n
    !> Domain location index
    Integer, Dimension(3), Public :: domain_ind

    !> Cells to operate over
    Real(Kind=wp), Dimension(3, 2), Public :: domain_bounds
    !> Cells to operate over
    Integer, Dimension(3, 2), Public :: domain_indices

  End Type kspace_type

Contains

  Subroutine setup_kspace(kspace_in, domain_in, kpoint_grid, comm)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to initialise the kspace type for a particular domain
    !! decomposition
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!----------------------------------------------------------------------!
    Type(kspace_type),     Intent(inout) :: kspace_in
    Type(domains_type),    Intent(In   ) :: domain_in
    Integer, Dimension(3), Intent(In   ) :: kpoint_grid
    Type(comms_type),      Intent(In   ) :: comm

    Integer :: fail, max_block

    kspace_in%k_vec_dim = kpoint_grid

    kspace_in%domain_n = [domain_in%nx, domain_in%ny, domain_in%nz]
    kspace_in%domain_ind = [domain_in%idx, domain_in%idy, domain_in%idz]

    ! Real values of kmax vectors

    kspace_in%k_vec_dim_real = Real(kspace_in%k_vec_dim, wp)
    kspace_in%k_vec_max_real = Maxval(kspace_in%k_vec_dim_real)
    kspace_in%k_vec_max = Maxval(kspace_in%k_vec_dim)

    ! 3d local indices

    kspace_in%domain_indices(:, 1) = kspace_in%domain_ind(:) * (kspace_in%k_vec_dim(:) / kspace_in%domain_n(:)) + 1
    kspace_in%domain_indices(:, 2) = (kspace_in%domain_ind(:) + 1) * (kspace_in%k_vec_dim(:) / kspace_in%domain_n(:))
    ! kspace_in%domain_indices(1, 1) = domain_in%idx * (kspace_in%k_vec_dim(1) / domain_in%nx) + 1
    ! kspace_in%domain_indices(1, 2) = (domain_in%idx + 1) * (kspace_in%k_vec_dim(1) / domain_in%nx)
    ! kspace_in%domain_indices(2, 1) = domain_in%idy * (kspace_in%k_vec_dim(2) / domain_in%ny) + 1
    ! kspace_in%domain_indices(2, 2) = (domain_in%idy + 1) * (kspace_in%k_vec_dim(2) / domain_in%ny)
    ! kspace_in%domain_indices(3, 1) = domain_in%idz * (kspace_in%k_vec_dim(3) / domain_in%nz) + 1
    ! kspace_in%domain_indices(3, 2) = (domain_in%idz + 1) * (kspace_in%k_vec_dim(3) / domain_in%nz)
    kspace_in%domain_bounds(:, 1) = Real(kspace_in%domain_indices(:, 1) - 1, wp)
    kspace_in%domain_bounds(:, 2) = Nearest(Real(kspace_in%domain_indices(:, 2), wp), -1.0_wp)

    ! domain local block limits of kmax space

    kspace_in%block_fac = kspace_in%k_vec_dim / kspace_in%domain_n

    ! set up the indexing arrays for each Dimension (NOT deallocated manually)
    max_block = Maxval(kspace_in%block_fac)
    Allocate (kspace_in%index(3, 1:max_block), Stat=fail)
    If (fail > 0) Call error_alloc('kspace index array', 'setup_kspace')
    kspace_in%index = 0

    ! set up the parallel fft and useful related quantities

    Call initialize_fft(3, kspace_in%k_vec_dim, kspace_in%domain_n, kspace_in%domain_ind, kspace_in%block_fac, &
      & comm%comm, kspace_in%context)

    Call pfft_indices(kspace_in%k_vec_dim(1), kspace_in%block_fac(1), kspace_in%domain_ind(1), &
                      kspace_in%domain_n(1), kspace_in%index(1, 1:kspace_in%block_fac(1)))
    Call pfft_indices(kspace_in%k_vec_dim(2), kspace_in%block_fac(2), kspace_in%domain_ind(2), &
                      kspace_in%domain_n(2), kspace_in%index(2, 1:kspace_in%block_fac(2)))
    Call pfft_indices(kspace_in%k_vec_dim(3), kspace_in%block_fac(3), kspace_in%domain_ind(3), &
                      kspace_in%domain_n(3), kspace_in%index(3, 1:kspace_in%block_fac(3)))

  End Subroutine setup_kspace

End Module kspace
