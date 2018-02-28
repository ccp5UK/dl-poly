Module parallel_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for distributed 3D DFFT commensurate with doman
! decomposition relying on the GPFA 1D DFFT
!
! DaFT - Daresbury advansed Fourier Transforms
!
! copyright - daresbury laboratory
! author    - i.j.bush august 2010
! amended   - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,        Only : wp
  Use setup_module, Only : pi

  Use comms_module ! access to the generalised mpi_wp and
                   ! the intrinsics in mpif.h/mpi-module
  Use gpfa_module, Only : gpfa_set

  Implicit None

  Public :: initialize_fft, summarize_fft, pfft, pfft_indices, pfft_length_ok

  Interface initialize_fft
     Module Procedure initialize_fft
  End Interface

  Interface summarize_fft
     Module Procedure summarize_fft
  End Interface

  Interface pfft
     Module Procedure fft_1d
     Module Procedure fft_3d
  End Interface

  Interface pfft_indices
     Module Procedure generate_indexing
  End Interface

  Private

  ! Size of array for factorization. It allows for the first max_factor-1 prime numbers in the
  ! factorization, so 10 allows for the factors 2, 3, 5, 7, 11, 13, 17, 19, 23

  Integer, Parameter :: max_factor = 10

  ! The number of lowest prime factors that have optimised implementations

  Integer, Parameter :: n_opt_factor = 1


  Type comms_descriptor
     Private
     Integer                              :: communicator
     Integer                              :: n_procs
     Integer                              :: my_proc
     Integer                              :: comms_steps
     Integer                              :: transfer_size
     Integer, Dimension( : ), Allocatable :: exchange
     Integer, Dimension( : ), Allocatable :: trigs_offset
     Logical, Dimension( : ), Allocatable :: first_half
  End Type comms_descriptor


  Type opt_fact_descriptor
     Private
     Integer                                            :: length
     Integer                                            :: Block
     Integer                                            :: sections
     Integer                                            :: local_steps
     Integer                                            :: n_my_sections
     Integer              , Dimension( : ), Allocatable :: my_sec_starts
     Complex( Kind = wp ) , Dimension( : ), Allocatable :: trigs
     Complex( Kind = wp ) , Dimension( : ), Allocatable :: trigs_conjg
     Complex( Kind = wp ) , Dimension( : ), Allocatable :: trigs_short
     Type( comms_descriptor )                           :: communications
  End Type opt_fact_descriptor


  Type other_fact_descriptor
     Private
     Integer                                            :: length
     Integer                                            :: comm
     Integer                                            :: rank
     Integer                                            :: size
     Complex( Kind = wp ) , Dimension( : ), Allocatable :: twiddles
  End Type other_fact_descriptor


  Type dim_descriptor
     Integer                                                :: length

     ! The factorisation of the dimension. The first max_factor - 1 are the
     ! powers of the first max_factor -1 prime numbers, the last number is
     ! the remaining, possibly composite, factor

     Integer,                     Dimension( 1:max_factor ) :: factors

     ! The descriptor for the fully optimised factors, currently only powers of 2

     Type( opt_fact_descriptor ), Dimension( 1:1 )          :: opt_facs

     ! The descriptor for the other factor

     Type( other_fact_descriptor )                          :: other_fac
  End Type dim_descriptor


  Type fft_descriptor
     Private
     Integer                                             :: context
     Integer                                             :: overall_communicator
     Integer                                             :: n_procs
     Integer                                             :: my_proc
     Integer                                             :: dimensionality
     Type( dim_descriptor ), Dimension( : ), Allocatable :: dims
  End Type fft_descriptor


  Integer :: next_context = 0

  Integer, Parameter                                        :: max_ffts = 50
  Type( fft_descriptor ), Dimension( 1:max_ffts )           :: set_up_ffts

  Integer, Parameter                                        :: max_dims  = 3
  Real( Kind = wp ),      Dimension( 1:5, 1:max_dims, 1:2 ) :: fft_times = 0.0_wp
  Integer          ,      Dimension( 1:2 )                  :: n_calls = 0

  Integer, Parameter                                        :: n_strip = 32


  Interface apply_twiddles
     Module Procedure apply_twiddles_3d
  End Interface


Contains


  Subroutine initialize_fft( n_dims, lengths, proc_grid, my_grid_pos, Block, communicator, context )

    ! This initialises things like the comms pattern and trigs
    ! data for the FFT.

    Implicit None

    Integer                       , Intent( In    ) :: n_dims
    Integer, Dimension( 1:n_dims ), Intent( In    ) :: lengths
    Integer, Dimension( 1:n_dims ), Intent( In    ) :: proc_grid
    Integer, Dimension( 1:n_dims ), Intent( In    ) :: my_grid_pos
    Integer, Dimension( 1:n_dims ), Intent( In    ) :: Block
    Integer                       , Intent( In    ) :: communicator
    Integer                       , Intent(   Out ) :: context

    Integer, Dimension( 1:max_factor ) :: factors
    Integer                            :: comms_error
    Integer                            :: error
    Integer                            :: temp
    Integer                            :: i, j

    ! Do some error checking on the input

    Do i = 1, n_dims

       ! Check that each processor will get an equal number
       ! of blocks

       temp = lengths( i ) / ( proc_grid( i ) * Block( i ) )
       temp = temp * proc_grid( i ) * Block( i )
       If ( temp /= lengths( i ) ) Then
          Call fft_error( -3, 'INITIALIZE_FFT', &
               message = 'The blocks are not evenly divided amongst ' // &
               'the processors.' )
       End If

    End Do

    ! Set return context value and update the counter

    next_context = next_context + 1
    context      = next_context

    ! Store the basic info about the FFT and the communications
    ! environment

    set_up_ffts( context )%context        = next_context
    set_up_ffts( context )%dimensionality = n_dims
    Call MPI_COMM_DUP( communicator, set_up_ffts( context )%overall_communicator, comms_error )
    If ( comms_error /= 0 ) Then
       Call fft_error( comms_error, 'INITIALIZE_FFT', &
            called_routine = 'MPI_COMM_DUP' )
    End If

    ! In the overall group find out my name and how
    ! many mates I've got

    Call MPI_COMM_SIZE( set_up_ffts( context )%overall_communicator, set_up_ffts( context )%n_procs, comms_error )
    If ( comms_error /= 0 ) Then
       Call fft_error( comms_error, 'INITIALIZE_FFT', &
            called_routine = 'MPI_COMM_SIZE' )
    End If
    Call MPI_COMM_RANK( set_up_ffts( context )%overall_communicator, set_up_ffts( context )%my_proc, comms_error )
    If ( comms_error /= 0 ) Then
       Call fft_error( comms_error, 'INITIALIZE_FFT', &
            called_routine = 'MPI_COMM_RANK' )
    End If

    ! Check that there are sufficient processors

    If ( Product( proc_grid ) > set_up_ffts( context )%n_procs ) Then
       Call fft_error( -4, 'INITIALIZE_FFT', &
            message = 'There are insufficient processors to ' // &
            'make up the processor cuboid.' )
    End If

    ! Put away space for data describing the FFT along each dimension.

    Allocate ( set_up_ffts( context )%dims( 1:n_dims ), Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'INITIALIZE_FFT', &
            message        = 'Failed to allocate memory for ' // &
            'SET_UP_FFTS( CONTEXT )%DIMS', &
            called_routine = 'ALLOCATE' )
    End If

    ! And set up each dimension in turn

    Do i = 1, n_dims

       set_up_ffts( context )%dims( i )%length = lengths( i )

       ! Factorize the dimension
       Call factor( lengths( i ), factors )

       ! Now set up the optimised factors
       Do j = 1, n_opt_factor
          Call set_opt_factor( i, lengths( i ), proc_grid( i ), my_grid_pos, block( i ), j, factors, &
               communicator, set_up_ffts( context )%dims( i )%opt_facs( j ) )
       End Do
!!$       If ( proc_grid( i ) > 1 ) Then
          Call set_other_factor( i, lengths( i ), proc_grid( i ), my_grid_pos, block( i ), n_opt_factor + 1, factors , &
               communicator, set_up_ffts( context )%dims( i )%other_fac )
!!$       Else
!!$          set_up_ffts( context )%dims( i )%other_fac%length = 1
!!$       End If

    End Do

  Contains

    Subroutine set_trigs( n, trigs )

      ! Simply sets up the trig factors for a forward FFT
      ! of length N

      Implicit None

      Integer                              , Intent( In    ) :: n
      Complex( Kind = wp ), Dimension( 0: ), Intent(   Out ) :: trigs

      Complex( Kind = wp )   :: angle, delta
      Integer                :: i

      delta = Cmplx( 0.0_wp, 2.0_wp * pi / Real( n, Kind = wp ), Kind = wp )
      delta = Exp( delta )

      trigs( 0 ) = Cmplx( 1.0_wp, 0.0_wp, Kind = wp )

      ! Use the recurrence relation to avoid lots of calls to EXP.

      angle = delta
      Do i = 1, n - 1
         trigs( i ) = angle
         angle = angle * delta
      End Do

    End Subroutine set_trigs

    Subroutine set_comms( overall_communicator, grid_pos, this_dim, Block, communications )

      ! Sets up the communications strategy.

      Implicit None

      Integer,                 Intent( In    ) :: overall_communicator
      Integer, Dimension( : ), Intent( In    ) :: grid_pos
      Integer,                 Intent( In    ) :: this_dim
      Integer,                 Intent( In    ) :: Block
      Type( comms_descriptor )                 :: communications

      Integer :: key
      Integer :: split_communicator
      Integer :: last_communicator
      Integer :: comms_error
      Integer :: n_procs
      Integer :: my_proc
      Integer :: comms_steps
      Integer :: bit_filter
      Integer :: trigs_jump
      Integer :: i

      ! First of all generate a communicator so that all
      ! processors that own a bit of a given FFT along this
      ! dimension are in that communicator. The method used
      ! is to succesively chop the cuboid up in directions
      ! perpendicular to the desired dimension until we have
      ! what we want

      Call MPI_COMM_DUP( overall_communicator, last_communicator, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_COMMS', &
              called_routine = 'MPI_COMM_DUP' )
      End If
      split_communicator = last_communicator

      Do i = 1, Size( grid_pos )
         If ( i == this_dim ) Then
            Cycle
         End If
         Call MPI_COMM_RANK( last_communicator, key, comms_error )
         If ( comms_error /= 0 ) Then
            Call fft_error( comms_error, 'SET_COMMS', &
                 called_routine = 'MPI_COMM_RANK' )
         End If
         Call MPI_COMM_SPLIT( last_communicator, grid_pos( i ), key, split_communicator, comms_error )
         If ( comms_error /= 0 ) Then
            Call fft_error( comms_error, 'SET_COMMS', &
                 called_routine = 'MPI_COMM_SPLIT' )
         End If
         Call MPI_COMM_FREE( last_communicator, comms_error )
         If ( comms_error /= 0 ) Then
            Call fft_error( comms_error, 'SET_COMMS', &
                 called_routine = 'MPI_COMM_FREE' )
         End If
         last_communicator = split_communicator
      End Do

      ! SPLIT_COMMUNICATOR now holds the desired communicator.
      ! First find out basic data about it.
      Call MPI_COMM_SIZE( split_communicator, n_procs, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_COMMS', &
              called_routine = 'MPI_COMM_SIZE' )
      End If
      Call MPI_COMM_RANK( split_communicator, my_proc, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_COMMS', &
              called_routine = 'MPI_COMM_RANK' )
      End If

      ! How many of the `butterflys' require communications.

      comms_steps = Int( Log( Real(n_procs,wp) + 0.01_wp ) / Log( 2.0_wp ) )

      ! Store the basic comms environment data about the 1D FFT
      ! down this dimension.

      communications%communicator = split_communicator
      communications%n_procs      = n_procs
      communications%my_proc      = my_proc
      communications%comms_steps  = comms_steps

      ! Now work out the comms strategy. This is a bit complicated,
      ! but basically if we write the processor ranks in binary and
      ! apply bit_filter as below to that it gives you the number
      ! of the processor with which this processor has to exchange data
      ! at the given step.

      bit_filter = Ishft( 1, comms_steps - 1 )
      trigs_jump = Block

      ! Exchange holds the rank of the processor with which to
      ! exchange, first_half flags if this is the first or the
      ! second bit of data at this node of the butterfly ( they
      ! are used slightly differently ). Trigs_offset controls
      ! which trig factors are to be used at this stage.

      Allocate ( communications%exchange( 1:comms_steps ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_COMMS', &
              message        = 'Failed to allocate memory for ' // &
              'COMMUNICATIONS%EXCHANGE', &
              called_routine = 'ALLOCATE' )
      End If
      Allocate ( communications%first_half( 1:comms_steps ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_COMMS', &
              message        = 'Failed to allocate memory for ' // &
              'COMMUNICATIONS%FIRST_HALF', &
              called_routine = 'ALLOCATE' )
      End If
      Allocate ( communications%trigs_offset( 1:comms_steps ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_COMMS', &
              message        = 'Failed to allocate memory for ' // &
              'COMMUNICATIONS%TRIGS_OFFSET', &
              called_routine = 'ALLOCATE' )
      End If

      Do i = 1, comms_steps
         communications%exchange( i ) = Ieor( bit_filter, my_proc )
         communications%first_half( i ) = Ibits( my_proc, comms_steps - i, 1 ) == 0
         communications%trigs_offset( i ) = Iand( my_proc, bit_filter - 1 ) * trigs_jump
         bit_filter = Ishft( bit_filter, -1 )
         trigs_jump = Ishft( trigs_jump,  1 )
      End Do

    End Subroutine set_comms

    Subroutine set_opt_factor( dim, length, n_proc, my_grid_pos, block, n_fac, factors, comm, dim_fac )

      Implicit None

      Integer                    , Intent( In    ) :: dim
      Integer                    , Intent( In    ) :: length
      Integer                    , Intent( In    ) :: n_proc
      Integer, Dimension( : )    , Intent( In    ) :: my_grid_pos
      Integer                    , Intent( In    ) :: block
      Integer                    , Intent( In    ) :: n_fac
      Integer, Dimension( : )    , Intent( In    ) :: factors
      Integer                    , Intent( In    ) :: comm
      Type( opt_fact_descriptor ), Intent(   Out ) :: dim_fac

      Integer, Dimension( : ), Allocatable :: new_pos

      Integer, Dimension( 1:max_factor )   :: proc_factors
      Integer                              :: fac_length
      Integer                              :: n_proc_groups, n_proc_in_group
      Integer                              :: sections, my_sections
      Integer                              :: fac_comm
      Integer                              :: comms_error
      Integer                              :: error
      Integer                              :: i

      ! How long is the bit due to this prime factor?
      fac_length = get_nth_prime( n_fac ) ** factors( n_fac )

      ! Include factors not covered by the processor grid in the first optimised factor
      ! These non powers of 2 will be covered by the serial fft routine
      If ( n_fac == 1 ) Then
         Call factor( n_proc, proc_factors )
         proc_factors = factors - proc_factors
         Do i = 2, max_factor
            fac_length = fac_length * ( get_nth_prime( i ) ** proc_factors( i ) )
         End Do
      End If

      ! And hence how many sections in this block
      sections    = fac_length / block

      ! How many groups will we split the processors into, i.e. how many transforms
      ! of this length will be going on at once
      n_proc_groups = length / fac_length
      n_proc_in_group = n_proc / n_proc_groups

      ! And so how many sections do I own
      my_sections = sections / n_proc_in_group

      ! And so basic info about the optimised factor for this dimension
      dim_fac%length        = fac_length
      dim_fac%block         = block
      dim_fac%sections      = sections
      dim_fac%n_my_sections = my_sections

      ! For this factor how many steps do not require communication?
      dim_fac%local_steps = Nint( Log( Real( my_sections, Kind = wp ) ) / Log( 2.0_wp ) )

      ! Work out where each local block starts
      Allocate ( dim_fac%my_sec_starts( 0:my_sections - 1 ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              message= 'Failed to allocate memory for DIM_FAC%MY_SEC_STARTS', &
              called_routine = 'ALLOCATE' )
      End If
      Do i = 0, my_sections - 1
         dim_fac%my_sec_starts( i ) = i * block
      End Do

      ! Now set up the various sets of trig factors

      ! Trig factors used for steps before the call to the FFT library routine in
      ! the forward transform
      Allocate ( dim_fac%trigs( 0:fac_length - 1 ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              message        = 'Failed to allocate memory for DIM_FAC%TRIGS', &
              called_routine = 'ALLOCATE' )
      End If
      Call set_trigs( fac_length, dim_fac%trigs )

      ! Trig factors used for steps after the call to the FFT library routine in
      ! the backward transform
      Allocate ( dim_fac%trigs_conjg( 0:fac_length - 1 ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              message        = 'Failed to allocate memory for DIM_FAC%TRIGS_CONJG', &
              called_routine = 'ALLOCATE' )
      End If
      dim_fac%trigs_conjg = Conjg( dim_fac%trigs )

      ! Trig factors for the local libray FFT routine
      Allocate ( dim_fac%trigs_short( 0:block - 1 ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              message        = 'Failed to allocate memory for DIM_FAC%TRIGS_SHORT', &
              called_routine = 'ALLOCATE' )
      End If
!!$      Call set_trigs( block, dim_fac%trigs_short )
      Call gpfa_set( dim_fac%trigs_short( 0:block - 1 ), block )

      ! Split the communicator so that each proc group is in a different one
      Call MPI_COMM_SPLIT( comm, ( my_grid_pos( dim ) * block ) / fac_length, my_grid_pos( dim ), fac_comm, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_OPT_FACTOR', &
              called_routine = 'MPI_COMM_SPLIT' )
      End If

      Allocate ( new_pos( 1:Size( my_grid_pos ) ), Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              message        = 'Failed to allocate memory for NEW_POS', &
              called_routine = 'ALLOCATE' )
      End If
      new_pos = my_grid_pos
      Call MPI_COMM_RANK( fac_comm, new_pos( dim ), comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              called_routine = 'MPI_COMM_RANK' )
      End If

      ! And finally set the communication patterns
      Call set_comms( fac_comm, my_grid_pos, dim, block, dim_fac%communications )
      dim_fac%communications%transfer_size = my_sections * block

      Deallocate ( new_pos, Stat = error )
      If ( error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              message        = 'Failed to deallocate memory for NEW_POS', &
              called_routine = 'DEALLOCATE' )
      End If

      Call MPI_COMM_FREE( fac_comm, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( error, 'SET_OPT_FACTOR', &
              called_routine = 'MPI_COMM_FREE' )
      End If

    End Subroutine set_opt_factor

    Subroutine set_other_factor( dim, length, n_proc, my_grid_pos, block, n_fac, factors, comm, dim_fac )

      Implicit None

      Integer                      , Intent( In    ) :: dim
      Integer                      , Intent( In    ) :: length
      Integer                      , Intent( In    ) :: n_proc
      Integer, Dimension( : )      , Intent( In    ) :: my_grid_pos
      Integer                      , Intent( In    ) :: block
      Integer                      , Intent( In    ) :: n_fac
      Integer, Dimension( : )      , Intent( In    ) :: factors
      Integer                      , Intent( In    ) :: comm
      Type( other_fact_descriptor ), Intent(   Out ) :: dim_fac

      Complex( Kind = wp )               :: angle, delta

      Integer, Dimension( 1:max_factor ) :: proc_factors

      Integer                            :: fac_length
      Integer                            :: k2, fac2
      Integer                            :: proc_opt
      Integer                            :: last_comm, split_comm, key
      Integer                            :: comms_error
      Integer                            :: i0
      Integer                            :: i

      Call factor( n_proc, proc_factors )

      ! Determine the length of the transforms for the non-optimised case
      fac_length = factors( Size( factors ) )
      Do i = n_fac, Size( factors ) - 1
!!$         fac_length = fac_length * get_nth_prime( i ) ** factors( i )
         fac_length = fac_length * get_nth_prime( i ) ** proc_factors( i )
      End Do

      dim_fac%length = fac_length

      ! Now the twiddle factors for scaling the transfrom before the optimised transforms
      If ( fac_length /= 1 ) Then
         fac2 = length / fac_length
         Allocate ( dim_fac%twiddles( 0:fac2 - 1 ), Stat = error )
         If ( error /= 0 ) Then
            Call fft_error( error, 'SET_OTHER_FACTOR', &
                 message        = 'Failed to allocate memory for DIM_FAC%TWIDDLES', &
                 called_routine = 'ALLOCATE' )
         End If
         k2 = ( my_grid_pos( dim ) * block ) / fac2
         delta = Cmplx( 0.0_wp, 2.0_wp * pi * Real( k2 , Kind = wp) / Real( length, Kind = wp ), Kind = wp )
         delta = Exp( delta )
         i0 = my_grid_pos( dim ) * block - fac2 * k2
         dim_fac%twiddles( 0 ) = Exp( Cmplx( 0.0_wp, 2.0_wp * pi * Real( k2 * i0 , Kind = wp ) / &
                                                     Real( length, Kind = wp ), Kind = wp ) )
         angle = dim_fac%twiddles( 0 )
         Do i = 1, fac2 - 1
            angle = angle * delta
            dim_fac%twiddles( i ) = angle
         End Do
      Else
         ! Be Careful - allocate to zero size so can pas as argument if need be and simplifies
         ! tidy up code
         Allocate ( dim_fac%twiddles( 0:-1 ), Stat = error )
         If ( error /= 0 ) Then
            Call fft_error( error, 'SET_OTHER_FACTOR', &
                 message        = 'Failed to allocate memory for DIM_FAC%TWIDDLES', &
                 called_routine = 'ALLOCATE' )
         End If
      End If

      ! Now work out a communicator that spans the procs required for the
      ! non-optimised FT (actually a simple DFT). For this need to factor the number of procs.

      ! How many procs in each optimised section?
      proc_opt = 1
      Do i = 1, n_opt_factor
         proc_opt = proc_opt * get_nth_prime( i ) ** proc_factors( i )
      End Do

      ! Now split the overall communicator into ones containing only columns along the dimension
      ! of interest. Do this by chopping up the (hyper-)cuboid along directions orthogonal to the
      ! dimension of interest
      Call MPI_COMM_DUP( comm, last_comm, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
              called_routine = 'MPI_COMM_DUP' )
      End If
      Do i = 1, Size( my_grid_pos )
         If ( i == dim ) Then
            Cycle
         End If
         Call MPI_COMM_RANK( last_comm, key, comms_error )
         If ( comms_error /= 0 ) Then
            Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
                 called_routine = 'MPI_COMM_RANK' )
         End If
         Call MPI_COMM_SPLIT( last_comm, my_grid_pos( i ), key, split_comm, comms_error )
         If ( comms_error /= 0 ) Then
            Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
                 called_routine = 'MPI_COMM_SPLIT' )
         End If
         Call MPI_COMM_FREE( last_comm, comms_error )
         If ( comms_error /= 0 ) Then
            Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
                 called_routine = 'MPI_COMM_FREE' )
         End If
         last_comm = split_comm
      End Do

      ! SPLIT_COMM is now a communicator containing procs along the dimension DIM
      ! Now split that to a comm that just contains procs involved in the unoptimised FTs
      Call MPI_COMM_SPLIT( split_comm, Mod( my_grid_pos( dim ), proc_opt ), my_grid_pos( dim ), dim_fac%comm, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
              called_routine = 'MPI_COMM_SPLIT' )
      End If

      Call MPI_COMM_RANK( dim_fac%comm, dim_fac%rank, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
              called_routine = 'MPI_COMM_RANK' )
      End If
      Call MPI_COMM_SIZE( dim_fac%comm, dim_fac%size, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
              called_routine = 'MPI_COMM_SIZE' )
      End If

      ! And tidy up
      Call MPI_COMM_FREE( split_comm, comms_error )
      If ( comms_error /= 0 ) Then
         Call fft_error( comms_error, 'SET_OTHER_FACTOR', &
              called_routine = 'MPI_COMM_FREE' )
      End If

    End Subroutine set_other_factor

  End Subroutine initialize_fft

  Subroutine summarize_fft( processor, context )

    Implicit None

    Integer, Intent( In    ) :: processor
    Integer, Intent( In    ) :: context

    Integer :: i, j, l

    If ( set_up_ffts( 1 )%my_proc == processor ) Then
       Write( Unit=* , Fmt='('' ---------------------------------------------------------------------- '')' )
       Write( Unit=* , Fmt=* ) 'FFT summary for processor ', processor

       Do l = 1, next_context

          ! Can input a negative context to summarize all the FFTs
          If ( context > 0 .and. l /= context ) Then
             Cycle
          End If

          Write( Unit=* , Fmt=* ) 'FFT with context ', l

          Write( Unit=* , Fmt=* ) 'nodes, dimensionality ', set_up_ffts( l )%n_procs, &
               set_up_ffts( l )%dimensionality

          Do i = 1, set_up_ffts( l )%dimensionality

             Write( Unit=* , Fmt=* )
             Write( Unit=* , Fmt=* ) 'Dim ', i

             Write( Unit=* , Fmt=* ) 'length ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%length, &
                  ' block ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%block, &
                  ' sections ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%sections, &
                  'local steps ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%local_steps, &
                  ' my secs ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%n_my_sections

             Write( Unit=* , Fmt=* ) 'comms'

             Write( Unit=* , Fmt=* ) &
                  'comm ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%communicator, &
                  'procs ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%n_procs, &
                  ' me   ', set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%my_proc, &
                  'comms steps ',set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%comms_steps

             Do j = 1, set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%comms_steps
                Write( Unit=* , Fmt=* ) 'Step ', j, ' Send to ', &
                     set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%exchange( j ), &
                     set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%first_half( j ), &
                     set_up_ffts( l )%dims( i )%opt_facs( 1 )%communications%trigs_offset( j )
             End Do

          End Do

          Write( Unit=* , Fmt=* )

       End Do

       Write( Unit=* , Fmt=* ) ' Total        local        comms         exch        gpfa'
       Do j = 1, 2
          Do i = 1, 3
             Write( Unit=* , Fmt='( 1x, 5( f8.5, 5x ) )' ) fft_times( :, i, j ) / Real(n_calls( 1 ),wp)
          End Do
       End Do
       Write( Unit=* , Fmt='('' ---------------------------------------------------------------------- '')' )
    End If

  End Subroutine summarize_fft

  Subroutine fft_1d( a, work, context, direction )

    Implicit None

    Complex( Kind = wp ), Dimension( 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0: ), Intent( InOut ) :: work
    Integer                              , Intent( In    ) :: context
    Integer                              , Intent( In    ) :: direction

    If ( direction == 1 ) Then
       n_calls( 1 ) = n_calls( 1 ) + 1
       Call forward_1d_fft( a, work, set_up_ffts( context ) )
    Else
       n_calls( 2 ) = n_calls( 2 ) + 1
    End If

  End Subroutine fft_1d

  Subroutine forward_1d_fft( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                 :: desc

    Complex( Kind = wp ), Dimension( : ), Allocatable :: work2

    Complex( Kind = wp )                              :: swap
    Complex( Kind = wp )                              :: trigs_fac

    Integer                                           :: n_in_group
    Integer                                           :: n_groups
    Integer                                           :: group
    Integer                                           :: trigs_ele
    Integer                                           :: start1, start2
    Integer                                           :: i_start, i_end
    Integer                                           :: len_sec
    Integer                                           :: comm
    Integer                                           :: xfer
    Integer                                           :: remote
    Integer                                           :: comms_tag = 10
    Integer                                           :: comms_request
    Integer                                           :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE )           :: comms_status
    Integer                                           :: error
    Integer                                           :: trigs_start, trigs_stride
    Integer                                           :: cut
    Integer                                           :: i, j, k

    Allocate ( work2( 0:desc%dims( 1 )%opt_facs( 1 )%length - 1 ), Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'FORWARD_1D_FFT', &
            message        = 'Failed to allocate memory for WORK2', &
            called_routine = 'ALLOCATE' )
    End If

    len_sec = desc%dims( 1 )%opt_facs( 1 )%block

    n_in_group = Ishft( desc%dims( 1 )%opt_facs( 1 )%n_my_sections, - 1 )
    n_groups   = 1

    trigs_start  = desc%dims( 1 )%opt_facs( 1 )%communications%my_proc * desc%dims( 1 )%opt_facs( 1 )%block
    trigs_stride = desc%dims( 1 )%opt_facs( 1 )%block * ( desc%dims( 1 )%opt_facs( 1 )%communications%n_procs - 1 )

    Do cut = 1, desc%dims( 1 )%opt_facs( 1 )%local_steps
       group = 0
       Do i = 1, n_groups
          trigs_ele = trigs_start
          Do j = 1, n_in_group
             start1 = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( group )
             start2 = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
             Do k = 0, len_sec - 1
                work( k ) = desc%dims( 1 )%opt_facs( 1 )%trigs( trigs_ele )
                trigs_ele = trigs_ele + n_groups
             End Do
             Do k = 0, len_sec - 1
                trigs_fac = work( k )
                swap = a( start1 + k )
                a( start1 + k ) = swap + a( start2 + k )
                a( start2 + k ) = swap - a( start2 + k )
                a( start2 + k ) = trigs_fac * a( start2 + k )
             End Do
             trigs_ele = trigs_ele + trigs_stride
             group = group + 1
          End Do
          group = group + n_in_group
       End Do
       n_groups   = Ishft( n_groups  ,  1 )
       n_in_group = Ishft( n_in_group, -1 )
    End Do

    trigs_stride = n_groups

    comm =     desc%dims( 1 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 1 )%opt_facs( 1 )%communications%transfer_size

    Do cut = 1, desc%dims( 1 )%opt_facs( 1 )%communications%comms_steps
       remote = desc%dims( 1 )%opt_facs( 1 )%communications%exchange( cut )
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag , comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       If ( desc%dims( 1 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          a = work + a
       Else
          Do k = 0, desc%dims( 1 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 1 )%opt_facs( 1 )%communications%trigs_offset( cut )
             i_start = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k )
             i_end   = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
             Do i = i_start, i_end
                work2( i ) = desc%dims( 1 )%opt_facs( 1 )%trigs( trigs_ele )
                trigs_ele = trigs_ele + trigs_stride
             End Do
             Do i = i_start, i_end
                a( i ) = ( work( i ) - a( i ) ) * work2( i )
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, 1 )
    End Do

    Deallocate ( work2, Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'FORWARD_1D_FFT', &
            message        = 'Failed to deallocate memory for WORK2', &
            called_routine = 'DEALLOCATE' )
    End If

    If ( len_sec > 1 ) Then
       Call gpfa_wrap( a, desc%dims( 1 )%opt_facs( 1 )%trigs_short, &
            2, 2 * len_sec, len_sec, &
            1, 1, 1 )
    End If

  End Subroutine forward_1d_fft

  Subroutine fft_3d( a, work, context, direction )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Integer                                      , Intent( In    ) :: context
    Integer                                      , Intent( In    ) :: direction

    If ( direction == 1 ) Then
       n_calls( 1 ) = n_calls( 1 ) + 1
       Call forward_3d_fft_x( a, work, set_up_ffts( context ) )
       Call forward_3d_fft_y( a, work, set_up_ffts( context ) )
       Call forward_3d_fft_z( a, work, set_up_ffts( context ) )
    Else
       n_calls( 2 ) = n_calls( 2 ) + 1
       Call back_3d_fft_z( a, work, set_up_ffts( context ) )
       Call back_3d_fft_y( a, work, set_up_ffts( context ) )
       Call back_3d_fft_x( a, work, set_up_ffts( context ) )
    End If

  End Subroutine fft_3d

  Subroutine forward_3d_fft_x( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                         :: desc

    Complex( Kind = wp )                    :: swap
    Complex( Kind = wp )                    :: trigs_fac

    Integer                                 :: n_in_group
    Integer                                 :: n_groups
    Integer                                 :: group
    Integer                                 :: trigs_ele
    Integer                                 :: start1, start2
    Integer                                 :: i_start, i_end
    Integer                                 :: len_sec
    Integer                                 :: comm
    Integer                                 :: xfer
    Integer                                 :: remote
    Integer                                 :: comms_tag = 10
    Integer                                 :: comms_request
    Integer                                 :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE ) :: comms_status
    Integer                                 :: trigs_start, trigs_stride
    Integer                                 :: n2, n3, size_plane
    Integer                                 :: cut
    Integer                                 :: m, n
    Integer                                 :: m_start, m_finish
    Integer                                 :: i, j, k

    Real( Kind = wp )                       :: start_total, start_exch, start
    Real( Kind = wp )                       :: finish_total, finish_exch, finish

    start_total = fft_time()
    start = start_total

    If ( desc%dims( 1 )%other_fac%length /= 1 ) Then
       Call fft_3d_other( a, work, desc%dims( 1 )%other_fac, 1 )
       Call apply_twiddles( a, desc%dims( 1 )%other_fac, 1, 1 )
    End If

    n2 = Ubound( a, Dim = 2 )
    n3 = Ubound( a, Dim = 3 )

    size_plane = ( n2 + 1 ) * ( n3 + 1 )

    len_sec = desc%dims( 1 )%opt_facs( 1 )%block

    n_groups = 1

    Do n = 0, n3
       Do m_start = 0, n2, n_strip

          m_finish = m_start + n_strip - 1

          n_in_group = Ishft( desc%dims( 1 )%opt_facs( 1 )%n_my_sections, - 1 )
          n_groups   = 1

          trigs_start  = desc%dims( 1 )%opt_facs( 1 )%communications%my_proc * desc%dims( 1 )%opt_facs( 1 )%block
          trigs_stride = desc%dims( 1 )%opt_facs( 1 )%block * ( desc%dims( 1 )%opt_facs( 1 )%communications%n_procs - 1 )

          Do cut = 1, desc%dims( 1 )%opt_facs( 1 )%local_steps
             group = 0
             Do i = 1, n_groups
                trigs_ele = trigs_start
                Do j = 1, n_in_group
                   start1 = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( group )
                   start2 = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
                   Do k = 0, len_sec - 1
                      work( k, 0, 0 ) = desc%dims( 1 )%opt_facs( 1 )%trigs( trigs_ele )
                      trigs_ele = trigs_ele + n_groups
                   End Do
                   Do m = m_start, m_finish
                      Do k = 0, len_sec - 1
                         trigs_fac = work( k, 0, 0 )
                         swap = a( start1 + k, m, n )
                         a( start1 + k, m, n ) = swap + a( start2 + k, m, n )
                         a( start2 + k, m, n ) = swap - a( start2 + k, m, n )
                         a( start2 + k, m, n ) = trigs_fac * a( start2 + k, m, n )
                      End Do
                   End Do
                   trigs_ele = trigs_ele + trigs_stride
                   group = group + 1
                End Do
                group = group + n_in_group
             End Do
             n_groups   = Ishft( n_groups  ,  1 )
             n_in_group = Ishft( n_in_group, -1 )
          End Do

       End Do
    End Do

    finish = fft_time()
    fft_times( 2, 1, 1 ) = fft_times( 2, 1, 1 ) + ( finish - start )
    start = finish

    trigs_stride = n_groups

    comm =     desc%dims( 1 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 1 )%opt_facs( 1 )%communications%transfer_size

    xfer = xfer * size_plane

    Do cut = 1, desc%dims( 1 )%opt_facs( 1 )%communications%comms_steps
       remote = desc%dims( 1 )%opt_facs( 1 )%communications%exchange( cut )
       start_exch = fft_time()
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag, comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       finish_exch = fft_time()
       fft_times( 4, 1, 1 ) = fft_times( 4, 1, 1 ) + ( finish_exch - start_exch )
       If ( desc%dims( 1 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          a = work + a
       Else
          a = work - a
          Do k = 0, desc%dims( 1 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 1 )%opt_facs( 1 )%communications%trigs_offset( cut )
             i_start = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k )
             i_end   = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
             Do i = i_start, i_end
                work( i, 0, 0 ) = desc%dims( 1 )%opt_facs( 1 )%trigs( trigs_ele )
                trigs_ele = trigs_ele + trigs_stride
             End Do
             Do n = 0, n3
                Do m = 0, n2
                   Do i = i_start, i_end
                      a( i, m, n ) = a( i, m, n ) * work( i, 0, 0 )
                   End Do
                End Do
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, 1 )
    End Do

    finish = fft_time()
    fft_times( 3, 1, 1 ) = fft_times( 3, 1, 1 ) + ( finish - start )
    start = finish

    If ( len_sec > 1 ) Then
       Call gpfa_wrap( a, desc%dims( 1 )%opt_facs( 1 )%trigs_short, &
            2, 2 * len_sec, len_sec, &
            desc%dims( 1 )%opt_facs( 1 )%n_my_sections * size_plane, 1, get_start_point( a, 0, 0, 0 ) )
    End If

    finish = fft_time()
    fft_times( 5, 1, 1 ) = fft_times( 5, 1, 1 ) + ( finish - start )
    finish_total = finish
    fft_times( 1, 1, 1 ) = fft_times( 1, 1, 1 ) + ( finish_total - start_total )

  End Subroutine forward_3d_fft_x

  Subroutine forward_3d_fft_y( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                         :: desc

    Complex( Kind = wp )                    :: swap
    Complex( Kind = wp )                    :: trigs_fac

    Integer                                 :: n_in_group
    Integer                                 :: n_groups
    Integer                                 :: group
    Integer                                 :: trigs_ele
    Integer                                 :: start1, start2
    Integer                                 :: len_sec
    Integer                                 :: comm
    Integer                                 :: xfer
    Integer                                 :: remote
    Integer                                 :: comms_tag = 10
    Integer                                 :: comms_request
    Integer                                 :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE ) :: comms_status
    Integer                                 :: trigs_start, trigs_stride
    Integer                                 :: n1, n3, size_plane
    Integer                                 :: cut
    Integer                                 :: m, n
    Integer                                 :: i, j, k

    Real( Kind = wp )                       :: start_total, start_exch, start
    Real( Kind = wp )                       :: finish_total, finish_exch, finish

    start_total = fft_time()
    start = start_total

    If ( desc%dims( 2 )%other_fac%length /= 1 ) Then
       Call fft_3d_other( a, work, desc%dims( 2 )%other_fac, 1 )
       Call apply_twiddles( a, desc%dims( 2 )%other_fac, 1, 2 )
    End If

    n1 = Ubound( a, Dim = 1 )
    n3 = Ubound( a, Dim = 3 )

    size_plane = ( n1 + 1 ) * ( n3 + 1 )

    len_sec = desc%dims( 2 )%opt_facs( 1 )%block

    n_in_group = Ishft( desc%dims( 2 )%opt_facs( 1 )%n_my_sections, - 1 )
    n_groups   = 1

    trigs_start  = desc%dims( 2 )%opt_facs( 1 )%communications%my_proc * desc%dims( 2 )%opt_facs( 1 )%block * n_groups
    trigs_stride = desc%dims( 2 )%opt_facs( 1 )%block * ( desc%dims( 2 )%opt_facs( 1 )%communications%n_procs - 1 )

    Do cut = 1, desc%dims( 2 )%opt_facs( 1 )%local_steps
       group = 0
       Do i = 1, n_groups
          trigs_ele = trigs_start
          Do j = 1, n_in_group
             start1 = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( group )
             start2 = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
             Do k = 0, len_sec - 1
                trigs_fac = desc%dims( 2 )%opt_facs( 1 )%trigs( trigs_ele )
                Do n = 0, n3
                   Do m = 0, n1
                      swap = a( m, start1 + k, n )
                      a( m, start1 + k, n ) = swap + a( m, start2 + k, n )
                      a( m, start2 + k, n ) = swap - a( m, start2 + k, n )
                      a( m, start2 + k, n ) = trigs_fac * a( m, start2 + k, n )
                   End Do
                End Do
                trigs_ele = trigs_ele + n_groups
             End Do
             trigs_ele = trigs_ele + trigs_stride
             group = group + 1
          End Do
          group = group + n_in_group
       End Do
       n_groups   = Ishft( n_groups  ,  1 )
       n_in_group = Ishft( n_in_group, -1 )
    End Do

    finish = fft_time()
    fft_times( 2, 2, 1 ) = fft_times( 2, 2, 1 ) + ( finish - start )
    start = finish

    trigs_stride = n_groups

    comm =     desc%dims( 2 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 2 )%opt_facs( 1 )%communications%transfer_size

    xfer = xfer * size_plane

    Do cut = 1, desc%dims( 2 )%opt_facs( 1 )%communications%comms_steps
       remote = desc%dims( 2 )%opt_facs( 1 )%communications%exchange( cut )
       start_exch = fft_time()
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag, comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       finish_exch = fft_time()
       fft_times( 4, 2, 1 ) = fft_times( 4, 2, 1 ) + ( finish_exch - start_exch )
       If ( desc%dims( 2 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          a = work + a
       Else
          Do k = 0, desc%dims( 2 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 2 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 2 )%opt_facs( 1 )%trigs( trigs_ele )
                Do n = 0, n3
                   Do m = 0, n1
                      a( m, i, n ) = ( - a( m, i, n ) + work( m, i, n ) ) * trigs_fac
                   End Do
                End Do
                trigs_ele = trigs_ele + trigs_stride
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, 1 )
    End Do

    finish = fft_time()
    fft_times( 3, 2, 1 ) = fft_times( 3, 2, 1 ) + ( finish - start )
    start = finish

    If ( len_sec > 1 ) Then
       Do i = 0, n3
          Do j = 0, desc%dims( 2 )%opt_facs( 1 )%n_my_sections - 1
             Call gpfa_wrap( a, desc%dims( 2 )%opt_facs( 1 )%trigs_short, &
                  2 * ( n1 + 1 ), 2, len_sec, &
                  n1 + 1, 1, get_start_point( a, 0, j * len_sec, i ) )
          End Do
       End Do
    End If

    finish = fft_time()
    finish_total = finish
    fft_times( 1, 2, 1 ) = fft_times( 1, 2, 1 ) + ( finish_total - start_total )
    fft_times( 5, 2, 1 ) = fft_times( 5, 2, 1 ) + ( finish - start )

  End Subroutine forward_3d_fft_y

  Subroutine forward_3d_fft_z( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                         :: desc

    Complex( Kind = wp )                    :: swap1, swap2
    Complex( Kind = wp )                    :: trigs_fac

    Integer                                 :: n_in_group
    Integer                                 :: n_groups
    Integer                                 :: group
    Integer                                 :: trigs_ele
    Integer                                 :: start1, start2
    Integer                                 :: len_sec
    Integer                                 :: comm
    Integer                                 :: xfer
    Integer                                 :: remote
    Integer                                 :: comms_tag = 10
    Integer                                 :: comms_request
    Integer                                 :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE ) :: comms_status
    Integer                                 :: trigs_start, trigs_stride
    Integer                                 :: n1, n2, size_plane
    Integer                                 :: cut
    Integer                                 :: m, n
    Integer                                 :: i, j, k

    Real( Kind = wp )                       :: start_total, start_exch, start
    Real( Kind = wp )                       :: finish_total, finish_exch, finish

    start_total = fft_time()
    start = start_total

    If ( desc%dims( 3 )%other_fac%length /= 1 ) Then
       Call fft_3d_other( a, work, desc%dims( 3 )%other_fac, 1 )
       Call apply_twiddles( a, desc%dims( 3 )%other_fac, 1, 3 )
    End If

    n1 = Ubound( a, Dim = 1 )
    n2 = Ubound( a, Dim = 2 )

    size_plane = ( n1 + 1 ) * ( n2 + 1 )

    len_sec = desc%dims( 3 )%opt_facs( 1 )%block

    n_in_group = Ishft( desc%dims( 3 )%opt_facs( 1 )%n_my_sections, - 1 )
    n_groups   = 1

    trigs_start  = desc%dims( 3 )%opt_facs( 1 )%communications%my_proc * desc%dims( 3 )%opt_facs( 1 )%block
    trigs_stride = desc%dims( 3 )%opt_facs( 1 )%block * ( desc%dims( 3 )%opt_facs( 1 )%communications%n_procs - 1 )

    Do cut = 1, desc%dims( 3 )%opt_facs( 1 )%local_steps
       group = 0
       Do i = 1, n_groups
          trigs_ele = trigs_start
          Do j = 1, n_in_group
             start1 = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( group )
             start2 = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
             Do k = 0, len_sec - 1
                trigs_fac = desc%dims( 3 )%opt_facs( 1 )%trigs( trigs_ele )
                Do n = 0, n2
                   Do m = 0, n1
                      swap1 = a( m, n, start1 + k )
                      swap2 = a( m, n, start2 + k )
                      a( m, n, start1 + k ) = swap1 + swap2
                      a( m, n, start2 + k ) = trigs_fac * ( swap1 - swap2 )
                   End Do
                End Do
                trigs_ele = trigs_ele + n_groups
             End Do
             trigs_ele = trigs_ele + trigs_stride
             group = group + 1
          End Do
          group = group + n_in_group
       End Do
       n_groups   = Ishft( n_groups  ,  1 )
       n_in_group = Ishft( n_in_group, -1 )
    End Do

    finish = fft_time()
    fft_times( 2, 3, 1 ) = fft_times( 2, 3, 1 ) + ( finish - start )
    start = finish

    trigs_stride = n_groups

    comm =     desc%dims( 3 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 3 )%opt_facs( 1 )%communications%transfer_size

    xfer = xfer * size_plane

    Do cut = 1, desc%dims( 3 )%opt_facs( 1 )%communications%comms_steps
       remote = desc%dims( 3 )%opt_facs( 1 )%communications%exchange( cut )
       start_exch = fft_time()
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag, comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       finish_exch = fft_time()
       fft_times( 4, 3, 1 ) = fft_times( 4, 3, 1 ) + ( finish_exch - start_exch )
       If ( desc%dims( 3 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          a = work + a
       Else
          a = work - a
          Do k = 0, desc%dims( 3 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 3 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 3 )%opt_facs( 1 )%trigs( trigs_ele )
                a( :, :, i ) = a( :, :, i ) * trigs_fac
                trigs_ele = trigs_ele + trigs_stride
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, 1 )
    End Do

    finish = fft_time()
    fft_times( 3, 3, 1 ) = fft_times( 3, 3, 1 ) + ( finish - start )
    start = finish

    If ( len_sec > 1 ) Then
       Do j = 0, desc%dims( 3 )%opt_facs( 1 )%n_my_sections - 1
          Call gpfa_wrap( a, desc%dims( 3 )%opt_facs( 1 )%trigs_short, &
               2 * size_plane, 2, len_sec, &
               size_plane, 1, get_start_point( a, 0, 0, j * len_sec ) )
       End Do
    End If

    finish = fft_time()
    finish_total = finish
    fft_times( 1, 3, 1 ) = fft_times( 1, 3, 1 ) + ( finish_total - start_total )
    fft_times( 5, 3, 1 ) = fft_times( 5, 3, 1 ) + ( finish - start )

  End Subroutine forward_3d_fft_z

  Subroutine back_3d_fft_x( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                         :: desc

    Complex( Kind = wp )                                      :: swap, tmp
    Complex( Kind = wp ), Dimension( 0:Ubound( a, Dim = 1 ) ) :: work2
    Complex( Kind = wp )                    :: trigs_fac

    Integer                                                   :: n_in_group
    Integer                                                   :: n_groups
    Integer                                                   :: group
    Integer                                                   :: trigs_ele
    Integer                                                   :: start1, start2
    Integer                                                   :: len_sec
    Integer                                                   :: comm
    Integer                                                   :: xfer
    Integer                                                   :: remote
    Integer                                                   :: comms_tag = 10
    Integer                                                   :: comms_request
    Integer                                                   :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE )                   :: comms_status
    Integer                                                   :: trigs_start, trigs_stride
    Integer                                                   :: n2, n3, size_plane
    Integer                                                   :: cut
    Integer                                                   :: m, n
    Integer                                                   :: m_start, m_finish
    Integer                                                   :: i, j, k

    Real( Kind = wp )                                         :: start_total, start_exch, start
    Real( Kind = wp )                                         :: finish_total, finish_exch, finish

    start_total = fft_time()
    start = start_total

    n2 = Ubound( a, Dim = 2 )
    n3 = Ubound( a, Dim = 3 )

    size_plane = ( n2 + 1 ) * ( n3 + 1 )

    len_sec = desc%dims( 1 )%opt_facs( 1 )%block

    If ( len_sec > 1 ) Then
       Call gpfa_wrap( a, desc%dims( 1 )%opt_facs( 1 )%trigs_short, &
            2, 2 * len_sec, len_sec, &
            desc%dims( 1 )%opt_facs( 1 )%n_my_sections * size_plane, -1, get_start_point( a, 0, 0, 0 ) )
    End If

    finish = fft_time()
    fft_times( 5, 1, 2 ) = fft_times( 5, 1, 2 ) + ( finish - start )
    start = finish

    trigs_stride = Ishft( desc%dims( 1 )%opt_facs( 1 )%n_my_sections * &
         desc%dims( 1 )%opt_facs( 1 )%communications%n_procs, -1 )

    comm =     desc%dims( 1 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 1 )%opt_facs( 1 )%communications%transfer_size

    xfer = xfer * size_plane

    Do cut = desc%dims( 1 )%opt_facs( 1 )%communications%comms_steps, 1, -1
       remote = desc%dims( 1 )%opt_facs( 1 )%communications%exchange( cut )
       start_exch = fft_time()
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag, comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       finish_exch = fft_time()
       fft_times( 4, 1, 2 ) = fft_times( 4, 1, 2 ) + ( finish_exch - start_exch )
       If ( desc%dims( 1 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          Do k = 0, desc%dims( 1 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 1 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 1 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                work2( i ) = trigs_fac
                trigs_ele = trigs_ele + trigs_stride
             End Do
             Do n = 0, n3
                Do m = 0, n2
                   Do i = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ), &
                        desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                      trigs_fac = work2( i )
                      a( i, m, n ) = trigs_fac * work( i, m, n ) + a( i, m, n )
                   End Do
                End Do
             End Do
          End Do
       Else
          Do k = 0, desc%dims( 1 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 1 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 1 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                work2( i ) = trigs_fac
                trigs_ele = trigs_ele + trigs_stride
             End Do
             Do n = 0, n3
                Do m = 0, n2
                   Do i = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ), &
                        desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                      trigs_fac = work2( i )
                      a( i, m, n ) = work( i, m, n ) - trigs_fac * a( i, m, n )
                   End Do
                End Do
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, -1 )
    End Do

    finish = fft_time()
    fft_times( 3, 1, 2 ) = fft_times( 3, 1, 2 ) + ( finish - start )
    start = finish

    Do n = 0, n3
       Do m_start = 0, n2, n_strip

          m_finish = m_start + n_strip - 1

          n_in_group = 1
          n_groups   = Ishft( desc%dims( 1 )%opt_facs( 1 )%n_my_sections, - 1 )

          trigs_start  = desc%dims( 1 )%opt_facs( 1 )%communications%my_proc * desc%dims( 1 )%opt_facs( 1 )%block
          trigs_stride = desc%dims( 1 )%opt_facs( 1 )%block * ( desc%dims( 1 )%opt_facs( 1 )%communications%n_procs - 1 )

          Do cut = desc%dims( 1 )%opt_facs( 1 )%local_steps, 1, -1
             group = 0
             Do i = 1, n_groups
                trigs_ele = trigs_start
                Do j = 1, n_in_group
                   start1 = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( group )
                   start2 = desc%dims( 1 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
                   Do k = 0, len_sec - 1
                      work( k, 0, 0 ) = desc%dims( 1 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                      trigs_ele = trigs_ele + n_groups
                   End Do
                   Do m = m_start, m_finish
                      Do k = 0, len_sec - 1
                         trigs_fac = work( k, 0, 0 )
                         swap = a( start1 + k, m, n )
                         tmp  = trigs_fac * a( start2 + k, m, n )
                         a( start1 + k, m, n ) = swap + tmp
                         a( start2 + k, m, n ) = swap - tmp
                      End Do
                   End Do
                   trigs_ele = trigs_ele + trigs_stride
                   group = group + 1
                End Do
                group = group + n_in_group
             End Do
             n_groups   = Ishft( n_groups  , -1 )
             n_in_group = Ishft( n_in_group,  1 )
          End Do

       End Do
    End Do

    If ( desc%dims( 1 )%other_fac%length /= 1 ) Then
       Call apply_twiddles( a, desc%dims( 1 )%other_fac, -1, 1 )
       Call fft_3d_other( a, work, desc%dims( 1 )%other_fac, -1 )
    End If

    finish = fft_time()
    fft_times( 2, 1, 2 ) = fft_times( 2, 1, 2 ) + ( finish - start )
    finish_total = finish
    fft_times( 1, 1, 2 ) = fft_times( 1, 1, 2 ) + ( finish_total - start_total )

  End Subroutine back_3d_fft_x

  Subroutine back_3d_fft_y( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                         :: desc

    Complex( Kind = wp )                    :: swap, tmp
    Complex( Kind = wp )                    :: trigs_fac

    Integer                                 :: n_in_group
    Integer                                 :: n_groups
    Integer                                 :: group
    Integer                                 :: trigs_ele
    Integer                                 :: start1, start2
    Integer                                 :: len_sec
    Integer                                 :: comm
    Integer                                 :: xfer
    Integer                                 :: remote
    Integer                                 :: comms_tag = 10
    Integer                                 :: comms_request
    Integer                                 :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE ) :: comms_status
    Integer                                 :: trigs_start, trigs_stride
    Integer                                 :: n1, n3, size_plane
    Integer                                 :: cut
    Integer                                 :: m, n
    Integer                                 :: i, j, k

    Real( Kind = wp )                       :: start_total, start_exch, start
    Real( Kind = wp )                       :: finish_total, finish_exch, finish

    start_total = fft_time()
    start = start_total

    n1 = Ubound( a, Dim = 1 )
    n3 = Ubound( a, Dim = 3 )

    size_plane = ( n1 + 1 ) * ( n3 + 1 )

    len_sec = desc%dims( 2 )%opt_facs( 1 )%block

    If ( len_sec > 1 ) Then
       Do i = 0, n3
          Do j = 0, desc%dims( 2 )%opt_facs( 1 )%n_my_sections - 1
             Call gpfa_wrap( a, desc%dims( 2 )%opt_facs( 1 )%trigs_short, &
                  2 * ( n1 + 1 ), 2, len_sec, &
                  n1 + 1, -1, get_start_point( a, 0, j * len_sec, i ) )
          End Do
       End Do
    End If

    finish = fft_time()
    fft_times( 5, 2, 2 ) = fft_times( 5, 2, 2 ) + ( finish - start )
    start = finish

    trigs_stride = Ishft( desc%dims( 2 )%opt_facs( 1 )%n_my_sections * &
         desc%dims( 2 )%opt_facs( 1 )%communications%n_procs, -1 )

    comm =     desc%dims( 2 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 2 )%opt_facs( 1 )%communications%transfer_size

    xfer = xfer * size_plane

    Do cut = desc%dims( 2 )%opt_facs( 1 )%communications%comms_steps, 1, -1
       remote = desc%dims( 2 )%opt_facs( 1 )%communications%exchange( cut )
       start_exch = fft_time()
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag, comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       finish_exch = fft_time()
       fft_times( 4, 2, 2 ) = fft_times( 4, 2, 2 ) + ( finish_exch - start_exch )
       If ( desc%dims( 2 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          Do k = 0, desc%dims( 2 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 2 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 2 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                Do n = 0, n3
                   Do m = 0, n1
                      a( m, i, n ) = trigs_fac * work( m, i, n ) + a( m, i, n )
                   End Do
                End Do
                trigs_ele = trigs_ele + trigs_stride
             End Do
          End Do
       Else
          Do k = 0, desc%dims( 2 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 2 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 2 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                Do n = 0, n3
                   Do m = 0, n1
                      a( m, i, n ) = work( m, i, n ) - trigs_fac * a( m, i, n )
                   End Do
                End Do
                trigs_ele = trigs_ele + trigs_stride
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, -1 )
    End Do

    finish = fft_time()
    fft_times( 3, 2, 2 ) = fft_times( 3, 2, 2 ) + ( finish - start )
    start = finish

    n_in_group = 1
    n_groups   = trigs_stride

    trigs_start  = desc%dims( 2 )%opt_facs( 1 )%communications%my_proc * desc%dims( 2 )%opt_facs( 1 )%block
    trigs_stride = desc%dims( 2 )%opt_facs( 1 )%block * ( desc%dims( 2 )%opt_facs( 1 )%communications%n_procs - 1 )

    Do cut = desc%dims( 2 )%opt_facs( 1 )%local_steps, 1, -1
       group = 0
       Do i = 1, n_groups
          trigs_ele = trigs_start
          Do j = 1, n_in_group
             start1 = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( group )
             start2 = desc%dims( 2 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
             Do k = 0, len_sec - 1
                trigs_fac = desc%dims( 2 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                Do n = 0, n3
                   Do m = 0, n1
                      swap = a( m, start1 + k, n )
                      tmp  = trigs_fac * a( m, start2 + k, n )
                      a( m, start1 + k, n ) = swap + tmp
                      a( m, start2 + k, n ) = swap - tmp
                   End Do
                End Do
                trigs_ele = trigs_ele + n_groups
             End Do
             trigs_ele = trigs_ele + trigs_stride
             group = group + 1
          End Do
          group = group + n_in_group
       End Do
       n_groups   = Ishft( n_groups  , -1 )
       n_in_group = Ishft( n_in_group,  1 )
    End Do

    If ( desc%dims( 2 )%other_fac%length /= 1 ) Then
       Call apply_twiddles( a, desc%dims( 2 )%other_fac, -1, 2 )
       Call fft_3d_other( a, work, desc%dims( 2 )%other_fac, -1 )
    End If

    finish = fft_time()
    fft_times( 2, 2, 2 ) = fft_times( 2, 2, 2 ) + ( finish - start )
    finish_total = finish
    fft_times( 1, 2, 2 ) = fft_times( 1, 2, 2 ) + ( finish_total - start_total )

  End Subroutine back_3d_fft_y

  Subroutine back_3d_fft_z( a, work, desc )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( fft_descriptor )                                         :: desc

    Complex( Kind = wp )                    :: swap, tmp
    Complex( Kind = wp )                    :: trigs_fac

    Integer                                 :: n_in_group
    Integer                                 :: n_groups
    Integer                                 :: group
    Integer                                 :: trigs_ele
    Integer                                 :: start1, start2
    Integer                                 :: len_sec
    Integer                                 :: comm
    Integer                                 :: xfer
    Integer                                 :: remote
    Integer                                 :: comms_tag = 10
    Integer                                 :: comms_request
    Integer                                 :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE ) :: comms_status
    Integer                                 :: trigs_start, trigs_stride
    Integer                                 :: n1, n2, size_plane
    Integer                                 :: cut
    Integer                                 :: m, n
    Integer                                 :: i, j, k

    Real( Kind = wp )                       :: start_total, start_exch, start
    Real( Kind = wp )                       :: finish_total, finish_exch, finish

    start_total = fft_time()
    start = start_total

    n1 = Ubound( a, Dim = 1 )
    n2 = Ubound( a, Dim = 2 )

    size_plane = ( n1 + 1 ) * ( n2 + 1 )

    len_sec = desc%dims( 3 )%opt_facs( 1 )%block

    If ( len_sec > 1 ) Then
       Do j = 0, desc%dims( 3 )%opt_facs( 1 )%n_my_sections - 1
          Call gpfa_wrap( a, desc%dims( 3 )%opt_facs( 1 )%trigs_short, &
               2 * size_plane, 2, len_sec, &
               size_plane, -1, get_start_point( a, 0, 0, j * len_sec ) )
       End Do
    End If

    finish = fft_time()
    fft_times( 5, 3, 2 ) = fft_times( 5, 3, 2 ) + ( finish - start )
    start = finish

    trigs_stride = Ishft( desc%dims( 3 )%opt_facs( 1 )%n_my_sections * &
         desc%dims( 3 )%opt_facs( 1 )%communications%n_procs, -1 )

    comm =     desc%dims( 3 )%opt_facs( 1 )%communications%communicator
    xfer = 2 * desc%dims( 3 )%opt_facs( 1 )%communications%transfer_size

    xfer = xfer * size_plane

    Do cut = desc%dims( 3 )%opt_facs( 1 )%communications%comms_steps, 1, -1
       remote = desc%dims( 3 )%opt_facs( 1 )%communications%exchange( cut )
       start_exch = fft_time()
       Call MPI_ISSEND(  a, xfer, wp_mpi, remote, comms_tag, comm, comms_request, comms_error )
       Call MPI_RECV( work, xfer, wp_mpi, remote, comms_tag, comm, comms_status, comms_error )
       Call MPI_WAIT( comms_request, comms_status, comms_error )
       finish_exch = fft_time()
       fft_times( 4, 3, 2 ) = fft_times( 4, 3, 2 ) + ( finish_exch - start_exch )
       If ( desc%dims( 3 )%opt_facs( 1 )%communications%first_half( cut ) ) Then
          Do k = 0, desc%dims( 3 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 3 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 3 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                a( :, :, i ) = a( :, :, i ) + trigs_fac * work( :, :, i )
                trigs_ele = trigs_ele + trigs_stride
             End Do
          End Do
       Else
          Do k = 0, desc%dims( 3 )%opt_facs( 1 )%n_my_sections - 1
             trigs_ele = desc%dims( 3 )%opt_facs( 1 )%communications%trigs_offset( cut )
             Do i = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( k ), &
                  desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( k ) + len_sec - 1
                trigs_fac = desc%dims( 3 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                Do n = 0, n2
                   Do m = 0, n1
                      a( m, n, i ) = work( m, n, i ) - trigs_fac * a( m, n, i )
                   End Do
                End Do
                trigs_ele = trigs_ele + trigs_stride
             End Do
          End Do
       End If
       trigs_stride = Ishft( trigs_stride, -1 )
    End Do

    finish = fft_time()
    fft_times( 3, 3, 2 ) = fft_times( 3, 3, 2 ) + ( finish - start )
    start = finish

    n_in_group = 1
    n_groups   = trigs_stride

    trigs_start  = desc%dims( 3 )%opt_facs( 1 )%communications%my_proc * desc%dims( 3 )%opt_facs( 1 )%block
    trigs_stride = desc%dims( 3 )%opt_facs( 1 )%block * ( desc%dims( 3 )%opt_facs( 1 )%communications%n_procs - 1 )

    Do cut = desc%dims( 3 )%opt_facs( 1 )%local_steps, 1, -1
       group = 0
       Do i = 1, n_groups
          trigs_ele = trigs_start
          Do j = 1, n_in_group
             start1 = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( group )
             start2 = desc%dims( 3 )%opt_facs( 1 )%my_sec_starts( group + n_in_group )
             Do k = 0, len_sec - 1
                trigs_fac = desc%dims( 3 )%opt_facs( 1 )%trigs_conjg( trigs_ele )
                Do n = 0, n2
                   Do m = 0, n1
                      swap = a( m, n, start1 + k )
                      tmp  = trigs_fac * a( m, n, start2 + k )
                      a( m, n, start1 + k ) = swap + tmp
                      a( m, n, start2 + k ) = swap - tmp
                   End Do
                End Do
                trigs_ele = trigs_ele + n_groups
             End Do
             trigs_ele = trigs_ele + trigs_stride
             group = group + 1
          End Do
          group = group + n_in_group
       End Do
       n_groups   = Ishft( n_groups  , -1 )
       n_in_group = Ishft( n_in_group,  1 )
    End Do

    If ( desc%dims( 3 )%other_fac%length /= 1 ) Then
       Call apply_twiddles( a, desc%dims( 3 )%other_fac, -1, 3 )
       Call fft_3d_other( a, work, desc%dims( 3 )%other_fac, -1 )
    End If

    finish = fft_time()
    fft_times( 2, 3, 2 ) = fft_times( 2, 3, 2 ) + ( finish - start )
    finish_total = finish
    fft_times( 1, 3, 2 ) = fft_times( 1, 3, 2 ) + ( finish_total - start_total )

  End Subroutine back_3d_fft_z

  Subroutine fft_3d_other( a, work, desc, sign )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: work
    Type( other_fact_descriptor )                                  :: desc
    Integer                                      , Intent( In    ) :: sign

    Complex( Kind = wp ), Dimension( :, :, : ), Allocatable :: work2

    Complex( Kind = wp )                                    :: t

    Integer                                                 :: nw1, nw2, nw3
    Integer                                                 :: up, down
    Integer                                                 :: which
    Integer                                                 :: rank, rem_rank
    Integer                                                 :: comms_request1, comms_request2
    Integer                                                 :: comms_error
    Integer, Dimension( 1:MPI_STATUS_SIZE )                 :: comms_status
    Integer                                                 :: error
    Integer                                                 :: pulse

    nw1 = Size( a, Dim = 1 )
    nw2 = Size( a, Dim = 2 )
    nw3 = Size( a, Dim = 3 )

    ! Perform the other radix FT by a simple DFT. As the whole design is based upon this
    ! being a small factor should be reasonably efficient unless somebody daft is using DaFT.

    ! Will pass the data around in simple systolic fashion and add to the local data, scaled by the
    ! appropriate phase.

    up   = Modulo( desc%rank - 1, desc%size )
    down = Modulo( desc%rank + 1, desc%size )

    Allocate ( work2( 0:nw1 - 1, 0:nw2 - 1, 0:nw3 - 1 ), Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'FFT_3D_OTHER', &
            message        = 'Failed to allocate memory for WORK2', &
            called_routine = 'ALLOCATE' )
    End If

    which = 0

    work = a

    a = (0.0_wp,0.0_wp)

    rank = desc%rank
    rem_rank = rank

    Do pulse = 0, desc%size - 1

       t = Exp( Cmplx( 0.0_wp, ( 2.0_wp * pi * Real( sign * rank * rem_rank , Kind = wp ) ) / &
                               Real( desc%length , Kind = wp ), Kind = wp ) )

       If ( which == 0 ) Then
          If ( pulse /= desc%size - 1 ) Then
             Call MPI_IRECV(  work2, 2 * Size( work2 ), wp_mpi, down, down, desc%comm, comms_request1, comms_error )
             Call MPI_ISSEND( work , 2 * Size( work  ), wp_mpi,   up, rank, desc%comm, comms_request2, comms_error )
          End If
          a = a + t * work
          If ( pulse /= desc%size - 1 ) Then
             Call MPI_WAIT( comms_request1, comms_status, comms_error )
             Call MPI_WAIT( comms_request2, comms_status, comms_error )
          End If
       Else
          If ( pulse /= desc%size - 1 ) Then
             Call MPI_IRECV(  work , 2 * Size( work  ), wp_mpi, down, down, desc%comm, comms_request1, comms_error )
             Call MPI_ISSEND( work2, 2 * Size( work2 ), wp_mpi,   up, rank, desc%comm, comms_request2, comms_error )
          End If
          a = a + t * work2
          If ( pulse /= desc%size - 1 ) Then
             Call MPI_WAIT( comms_request1, comms_status, comms_error )
             Call MPI_WAIT( comms_request2, comms_status, comms_error )
          End If
       End If

       which = Mod( which + 1, 2 )
       rem_rank = Mod( rem_rank + 1, desc%size )

    End Do

    Deallocate ( work2, Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'FFT_3D_OTHER', &
            message        = 'Failed to deallocate memory for WORK2', &
            called_routine = 'DEALLOCATE' )
    End If
  End Subroutine fft_3d_other

  Subroutine apply_twiddles_3d( a, desc, sign, dim )

    Implicit None

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( InOut ) :: a
    Type( other_fact_descriptor )                                  :: desc
    Integer                                      , Intent( In    ) :: sign
    Integer                                      , Intent( In    ) :: dim

    Complex( Kind = wp ) :: t

    Integer              :: nw1, nw2, nw3
    Integer              :: i, j, k

    nw1 = Size( a, Dim = 1 )
    nw2 = Size( a, Dim = 2 )
    nw3 = Size( a, Dim = 3 )

    Select Case( dim )
    Case( 1 )
       If ( sign == 1 ) Then
          Do k = 0, nw3 - 1
             Do j = 0, nw2 - 1
                Do i = 0, nw1 - 1
                   a( i, j, k ) = a( i, j , k ) * desc%twiddles( i )
                End Do
             End Do
          End Do
       Else
          Do k = 0, nw3 - 1
             Do j = 0, nw2 - 1
                Do i = 0, nw1 - 1
                   a( i, j, k ) = a( i, j , k ) * Conjg( desc%twiddles( i ) )
                End Do
             End Do
          End Do
       End If
    Case( 2 )
       If ( sign == 1 ) Then
          Do k = 0, nw3 - 1
             Do j = 0, nw2 - 1
                t = desc%twiddles( j )
                Do i = 0, nw1 - 1
                   a( i, j, k ) = a( i, j , k ) * t
                End Do
             End Do
          End Do
       Else
          Do k = 0, nw3 - 1
             Do j = 0, nw2 - 1
                t = Conjg( desc%twiddles( j ) )
                Do i = 0, nw1 - 1
                   a( i, j, k ) = a( i, j , k ) * t
                End Do
             End Do
          End Do
       End If
    Case( 3 )
       If ( sign == 1 ) Then
          Do k = 0, nw3 - 1
             t = desc%twiddles( k )
             Do j = 0, nw2 - 1
                Do i = 0, nw1 - 1
                   a( i, j, k ) = a( i, j , k ) * t
                End Do
             End Do
          End Do
       Else
          Do k = 0, nw3 - 1
             t = Conjg( desc%twiddles( k ) )
             Do j = 0, nw2 - 1
                Do i = 0, nw1 - 1
                   a( i, j, k ) = a( i, j , k ) * t
                End Do
             End Do
          End Do
       End If
    Case Default
       Call fft_error( Huge( 1 ), 'APPLY_TWIDDLES_3D', &
            message = 'Internal error - impossible dimension' )
    End Select

  End Subroutine apply_twiddles_3d

  Subroutine generate_indexing( n_tot, block, me_tot, nproc_tot, my_indices )

    Implicit None

    Integer                , Intent( In    ) :: n_tot
    Integer                , Intent( In    ) :: block
    Integer                , Intent( In    ) :: me_tot
    Integer                , Intent( In    ) :: nproc_tot
    Integer, Dimension( : ), Intent(   Out ) :: my_indices

    Integer, Dimension( : ), Allocatable :: all_indices
    Integer, Dimension( : ), Allocatable :: temp_indices

    Integer, Dimension( 1:2 )            :: proc_facs

    Integer                              :: n, me, nproc
    Integer                              :: log_nproc
    Integer                              :: error
    Integer                              :: cuts
    Integer                              :: n_sec, sec_start, sec_length
    Integer                              :: this_index
    Integer                              :: my_start, my_end
    Integer                              :: i, j, k

    n     = n_tot
    me    = me_tot
    nproc = nproc_tot

    Call factor( nproc, proc_facs )

    nproc = 2 ** proc_facs( 1 )

    n = n / proc_facs( 2 )

    me = Mod( me_tot, nproc )

    log_nproc = Nint( Log( Real( nproc, Kind = wp ) ) / Log( 2.0_wp ) )

    Allocate ( all_indices( 1:n_tot ), Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'GENERATE_INDEXING', &
            message        = 'Failed to allocate memory for ALL_INDICES', &
            called_routine = 'ALLOCATE' )
    End If
    Allocate ( temp_indices( 1:n ), Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'GENERATE_INDEXING', &
            message        = 'Failed to allocate memory for TEMP_INDICES', &
            called_routine = 'ALLOCATE' )
    End If

    temp_indices = (/ ( i, i = 1, n ) /)

    If ( proc_facs( 2 ) > 1 ) Then
       k = 1
       Do i = 1, n
          Do j = 0, proc_facs( 2 ) - 1
             all_indices( temp_indices( i ) + j * n ) = k
             k = k + 1
          End Do
       End Do
    Else
       all_indices = temp_indices
    End If

    Do k = 0, proc_facs( 2 ) - 1
       n_sec = 1
       sec_length = n
       Do cuts = 1, log_nproc
          sec_start = 1
          Do i = 1, n_sec
             this_index = sec_start
             Do j = sec_start, sec_start + sec_length / 2 - 1
                temp_indices( j ) = all_indices( this_index + k * n )
                this_index = this_index + 2
             End Do
             this_index = sec_start + 1
             Do j = sec_start + sec_length / 2, sec_start + sec_length - 1
                temp_indices( j ) = all_indices( this_index + k * n )
                this_index = this_index + 2
             End Do
             sec_start = sec_start + sec_length
          End Do
          n_sec = n_sec * 2
          sec_length = sec_length / 2
          all_indices( k * n + 1:( k + 1 ) * n )  = temp_indices
       End Do
    End Do

    my_start = block * me_tot + 1
    my_end   = my_start + block - 1

    my_indices = all_indices( my_start:my_end )

    Deallocate ( temp_indices, Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'GENERATE_INDEXING', &
            message        = 'Failed to deallocate memory for ALL_INDICES', &
            called_routine = 'DEALLOCATE' )
    End If
    Deallocate ( all_indices, Stat = error )
    If ( error /= 0 ) Then
       Call fft_error( error, 'GENERATE_INDEXING', &
            message        = 'Failed to deallocate memory for TEMP_INDICES', &
            called_routine = 'DEALLOCATE' )
    End If
  End Subroutine generate_indexing

  Function get_start_point( a, ix, iy, iz )

    Implicit None

    Integer                                                        :: get_start_point

    Complex( Kind = wp ), Dimension( 0:, 0:, 0: ), Intent( In    ) :: a
    Integer                                      , Intent( In    ) :: ix
    Integer                                      , Intent( In    ) :: iy
    Integer                                      , Intent( In    ) :: iz

    Integer :: l, m, n, lm

    l = Size( a, Dim = 1 )
    m = Size( a, Dim = 2 )
    n = Size( a, Dim = 3 )

    lm = l * m

    get_start_point = ix + iy * l + iz * lm
    get_start_point = get_start_point * 2 + 1

  End Function get_start_point

  Subroutine fft_error( status, calling_routine, message, called_routine )

    Implicit None

    Integer             , Intent( In    )           :: status
    Character( Len = * ), Intent( In    )           :: calling_routine
    Character( Len = * ), Intent( In    ), Optional :: message
    Character( Len = * ), Intent( In    ), Optional :: called_routine

    Integer :: comms_error

    If ( status == 0 ) Then
       Return
    End If
    If ( status > 0 ) Then
       Write( Unit=* , Fmt=* ) 'ERROR detected in routine ' // calling_routine
    Else
       Write( Unit=* , Fmt=* ) 'WARNING from routine ' // calling_routine
    End If
    Write( Unit=* , Fmt=* ) 'The status returned is ', status
    If ( Present( message ) ) Then
       Write( Unit=* , Fmt=* ) 'This corresponds to: '
       Write( Unit=* , Fmt=* ) message
    End If
    If ( Present( called_routine ) ) Then
       Write( Unit=* , Fmt=* ) 'The above status is implementation defined.'
       Write( Unit=* , Fmt=* ) 'The routine/statement flagging the error is ' // &
            called_routine
    End If
    If ( status > 0 ) Then
       Call MPI_ABORT( set_up_ffts( next_context )%overall_communicator, set_up_ffts( next_context )%my_proc, comms_error )
    End If

  End Subroutine fft_error

  Function fft_time()

    Implicit None

    Real( Kind = wp ) :: fft_time

    Integer           :: count, count_rate, count_max

    Call System_clock( count, count_rate, count_max )

    If ( count_rate > 0 ) Then
       fft_time = Real( count , Kind = wp ) / Real( count_rate , Kind = wp )
    Else
       fft_time = 0.0_wp
    End If

  End Function fft_time

  Subroutine factor( n, facs )

    Implicit None

    Integer                , Intent( In    ) :: n
    Integer, Dimension( : ), Intent(   Out ) :: facs

    Integer :: left
    Integer :: p
    Integer :: i

    facs = 0

    left = n
    Do i = 1, Size( facs ) - 1
       p = get_nth_prime( i )

       If ( p <= 0 ) Exit

       Do While ( p * ( left / p ) == left )
          left = left / p

          facs( i ) = facs( i ) + 1
       End Do
    End Do

    facs( Size( facs ) ) = left

  End Subroutine factor

  Function get_nth_prime( n )

    Implicit None

    Integer                  :: get_nth_prime

    Integer, Intent( In    ) :: n

    Integer, Dimension( 1:170 ), Parameter :: primes = (/                             &
           2,      3,      5,      7,     11,     13,     17,     19,     23,     29, &
          31,     37,     41,     43,     47,     53,     59,     61,     67,     71, &
          73,     79,     83,     89,     97,    101,    103,    107,    109,    113, &
         127,    131,    137,    139,    149,    151,    157,    163,    167,    173, &
         179,    181,    191,    193,    197,    199,    211,    223,    227,    229, &
         233,    239,    241,    251,    257,    263,    269,    271,    277,    281, &
         283,    293,    307,    311,    313,    317,    331,    337,    347,    349, &
         353,    359,    367,    373,    379,    383,    389,    397,    401,    409, &
         419,    421,    431,    433,    439,    443,    449,    457,    461,    463, &
         467,    479,    487,    491,    499,    503,    509,    521,    523,    541, &
         547,    557,    563,    569,    571,    577,    587,    593,    599,    601, &
         607,    613,    617,    619,    631,    641,    643,    647,    653,    659, &
         661,    673,    677,    683,    691,    701,    709,    719,    727,    733, &
         739,    743,    751,    757,    761,    769,    773,    787,    797,    809, &
         811,    821,    823,    827,    829,    839,    853,    857,    859,    863, &
         877,    881,    883,    887,    907,    911,    919,    929,    937,    941, &
         947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013 /)

    If ( n <= Size( primes ) ) Then
       get_nth_prime = primes( n )
    Else
       get_nth_prime = -1
    End If

  End Function get_nth_prime

  Function pfft_length_ok( n )

    ! Function to check that the a given length is OK for the serial FFT used
    ! by the module.
    ! THIS VERSION IS FOR GPFA. GPFA supports multiples of 2, 3 and 5.

    Implicit None

    Logical                  :: pfft_length_ok

    Integer, Intent( In    ) :: n

    Integer :: left

    left = n
    Do While ( ( left / 2 ) * 2 == left )
       left = left / 2
    End Do
    Do While ( ( left / 3 ) * 3 == left )
       left = left / 3
    End Do
    Do While ( ( left / 5 ) * 5 == left )
       left = left / 5
    End Do
    pfft_length_ok = left == 1

  End Function pfft_length_ok

End Module parallel_fft
