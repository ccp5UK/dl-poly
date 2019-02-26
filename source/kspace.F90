module kspace
  !!----------------------------------------------------------------------!
  !!
  !! dl_poly_4 module for containing types and functions
  !! relating to ewald kspace
  !!
  !! copyright - daresbury laboratory
  !! author    - j.s.wilkins october 2018
  !!
  !!----------------------------------------------------------------------!
  Use kinds, only : wi, wp
  Use errors_warnings
  implicit none

  Private
  public :: setup_kspace
  
  Type, Public :: kspace_type
     Private

     !> SPME kspace dimensions
     Integer( Kind = wi ), Dimension( 3 ),                  Public :: k_vec_dim, k_vec_dim_cont
     !> Largest Integer K-Vector Index
     Integer( Kind = wi ),                                  Public :: k_vec_max

     !> Real K-Vector indices
     Real( Kind = wp ),    Dimension( 3 ),                  Public :: k_vec_dim_real
     !> Largest K-Vector index
     Real( Kind = wp ),                                     Public :: k_vec_max_real

     !> blocking factors for splines and fft
     Integer,                                               Public :: block_x,block_y,block_z
     !> indexing arrays for x, y & z as used in parallel fft
     Integer,              Dimension( : ),     Allocatable, Public :: index_x,index_y,index_z
     !> context for parallel fft
     Integer,                                               Public :: context

     !> Cells to operate over
     Real( Kind = wp ),    Dimension( 3, 2 ),               Public  :: domain_bounds
     !> Cells to operate over
     Integer          ,    Dimension( 3, 2 ),               Public  :: domain_indices

  End Type kspace_type

contains
  
  subroutine setup_kspace(kspace_in, domain_in, kpoint_grid)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 routine to initialise the kspace type for a particular domain
    !! decomposition
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins february 2019
    !!
    !!----------------------------------------------------------------------!
    use domains, only : domains_type
    type ( kspace_type ),  intent( inout ) :: kspace_in
    type ( domains_type ),    intent( in    ) :: domain_in
    integer, dimension(3), intent( in    ) :: kpoint_grid
    integer :: fail
    
    kspace_in%k_vec_dim = kpoint_grid

    ! Real values of kmax vectors

    kspace_in%k_vec_dim_real = Real(kspace_in%k_vec_dim,wp)
    kspace_in%k_vec_max_real = maxval(kspace_in%k_vec_dim_real)
    kspace_in%k_vec_max = maxval(kspace_in%k_vec_dim)

    ! 3d local indices

    kspace_in%domain_indices(1,1)=domain_in%idx    *(kspace_in%k_vec_dim(1)/domain_in%nx)+1
    kspace_in%domain_indices(1,2)=(domain_in%idx+1)*(kspace_in%k_vec_dim(1)/domain_in%nx)
    kspace_in%domain_indices(2,1)=domain_in%idy    *(kspace_in%k_vec_dim(2)/domain_in%ny)+1
    kspace_in%domain_indices(2,2)=(domain_in%idy+1)*(kspace_in%k_vec_dim(2)/domain_in%ny)
    kspace_in%domain_indices(3,1)=domain_in%idz    *(kspace_in%k_vec_dim(3)/domain_in%nz)+1
    kspace_in%domain_indices(3,2)=(domain_in%idz+1)*(kspace_in%k_vec_dim(3)/domain_in%nz)
    kspace_in%domain_bounds(:,1) = Real   (      kspace_in%domain_indices(:,1)-1,wp)
    kspace_in%domain_bounds(:,2) = nearest( Real(kspace_in%domain_indices(:,2)  ,wp), -1.0_wp)
    
    ! domain local block limits of kmax space

    kspace_in%block_x = kspace_in%k_vec_dim(1) / domain_in%nx
    kspace_in%block_y = kspace_in%k_vec_dim(2) / domain_in%ny
    kspace_in%block_z = kspace_in%k_vec_dim(3) / domain_in%nz

    ! set up the indexing arrays for each Dimension (NOT deallocated manually)

    Allocate ( &
      & kspace_in%index_x( 1:kspace_in%block_x ), &
      & kspace_in%index_y( 1:kspace_in%block_y ), &
      & kspace_in%index_z( 1:kspace_in%block_z ), Stat = fail )
    If (fail > 0) call error_alloc('kspace index array','setup_kspace')

  end subroutine setup_kspace
  
end module kspace
