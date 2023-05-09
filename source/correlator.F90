Module correlators

    !!-----------------------------------------------------------------------
    !!
    !! Module for storing online correlation function calculation data
    !!
    !! copyright - daresbury laboratory
    !! author - h.devereux February 2023 
    !! 
    !!-----------------------------------------------------------------------
  
    Use configuration,    Only: configuration_type
    Use kinds,            Only: wp, wi
    Use comms,            Only: comms_type, mpi_distribution_type, &
                                gatherv_scatterv_index_arrays, &
                                ggatherv, gscatterv, root_id
    Use errors_warnings,  Only: error,&
                                error_alloc,&
                                error_dealloc,&
                                info,&
                                warning
    Implicit None
  
    Private

    ! default parameters
    Integer, Parameter, Public :: DEFAULT_BLOCKS = 1, DEFAULT_POINTS = 100, DEFAULT_WINDOW = 1

    Type, Public :: correlator
      Real(Kind=wp),     Allocatable :: left_accumulator(:,:)
      Real(Kind=wp),     Allocatable :: right_accumulator(:,:)
      Integer(Kind=wi),  Allocatable :: count_accumulated(:)
      Real(Kind=wp),     Allocatable :: left_shift(:,:,:)
      Logical,           Allocatable :: shift_not_null(:,:,:)
      Real(Kind=wp),     Allocatable :: right_shift(:,:,:)
      Integer(Kind=wi),  Allocatable :: shift_index(:)
      Real(Kind=wp),     Allocatable :: correlation(:,:,:,:)
      Integer(Kind=wi),  Allocatable :: count_correlated(:,:)
      Integer(Kind=wi)               :: number_of_blocks, max_block_used, left_dim, &
                                        right_dim, window_size, points_per_block, min_dist, &
                                        buffer_size, count_updated
  
    Contains
      Private
  
      Procedure, Public              :: init => allocate_correlator_arrays
      Procedure, Public              :: update
      Procedure, Private             :: add
      Procedure, Public              :: get_correlation
      Procedure, Public              :: deport_buffer
      Procedure, Public              :: recieve_buffer
      Final                          :: cleanup
  
    End Type

    Type, Public :: correlator_buffer_type
      ! correlators from all processes stored in a packed form
      Real(Kind=wp), Public, Allocatable :: buffer(:)
      ! MPI indexing required to unpack correlators stored in packed form
      Type(mpi_distribution_type), Public :: mpi
    Contains
      Procedure, Public :: initialise => initialise_correlator_buffer_type
      Procedure, Public :: finalise   => finalise_correlator_buffer_type
    End Type correlator_buffer_type

    Type, Public :: indices_buffer_type
      Integer, Public, Allocatable :: buffer(:)
      ! MPI indexing required to unpack correlators stored in packed form
      Type(mpi_distribution_type), Public :: mpi
    Contains
      Procedure, Public :: initialise => initialise_indices_buffer_type
      Procedure, Public :: finalise   => finalise_indices_buffer_type
    End Type indices_buffer_type

  Contains

  Subroutine initialise_correlator_buffer_type(this,comm,size_of)
    Class(correlator_buffer_type), Intent(InOut) :: this
    Type(comms_type),              Intent(In   ) :: comm
    Integer,                       Intent(In   ) :: size_of

    Allocate (this%buffer(size_of))
    Call this%mpi%init(comm)
  End Subroutine initialise_correlator_buffer_type

  Subroutine finalise_correlator_buffer_type(this)
    Class(correlator_buffer_type), Intent(InOut) :: this

    If (Allocated(this%buffer)) Then
      Deallocate(this%buffer)
    End If

    Call this%mpi%finalise()
  End Subroutine finalise_correlator_buffer_type

  Subroutine initialise_indices_buffer_type(this,comm,size_of)
    Class(indices_buffer_type),    Intent(InOut) :: this
    Type(comms_type),              Intent(In   ) :: comm
    Integer,                       Intent(In   ) :: size_of

    Allocate (this%buffer(size_of))
    Call this%mpi%init(comm)
  End Subroutine initialise_indices_buffer_type

  Subroutine finalise_indices_buffer_type(this)
    Class(indices_buffer_type), Intent(InOut) :: this

    If (Allocated(this%buffer)) Then
      Deallocate(this%buffer)
    End If

    Call this%mpi%finalise()
  End Subroutine finalise_indices_buffer_type
  
  Subroutine get_correlation(this, correlation)
    Class(correlator),                                     Intent(InOut)  :: this
    Real(Kind=wp), Dimension(& 
      1:this%points_per_block*this%number_of_blocks, &
      1:this%left_dim, 1:this%right_dim) ,                 Intent(InOut)  :: correlation
    Integer                                                               :: im, i, k
    
    correlation = 0.0

    im = 1
    Do i = 1,this%points_per_block
      If (this%count_correlated(1,i) > 0) Then
        correlation(im,:,:) = correlation(im,:,:) + this%correlation(1,i,:,:) / this%count_correlated(1,i)
        im = im + 1
      End If
    End Do
    Do k = 2,this%max_block_used
      Do i = this%min_dist,this%points_per_block
        If (this%count_correlated(k,i) > 0) Then
          correlation(im,:,:) = correlation(im,:,:) + this%correlation(k,i,:,:) / this%count_correlated(k,i)
          im = im + 1
        EndIf
      End Do
    End Do

  End Subroutine get_correlation

  Recursive Subroutine add(this, data_left, data_right, block_index)
    Class(correlator),   Intent(InOut)  :: this
    Real(Kind=wp),       Intent(In   )  :: data_left(:), data_right(:)
    Integer,             Intent(In   )  :: block_index
    Integer                             :: i, j, s, n, l, r
    Logical                             :: flag
    
    If (block_index > this%number_of_blocks) Then
      Return 
    End If

    s = this%shift_index(block_index) 

    ! indicate blocks at this depth have data
    !   and add to accumulators
    If (block_index > this%max_block_used) Then
      this%max_block_used = block_index
    End If
    ! add new data to the shifts
    !  and accumulate
    this%left_shift(block_index,s,:) = data_left
    this%left_accumulator(block_index,:) = this%left_accumulator(block_index,:) + data_left
    this%shift_not_null(block_index,s,1) = .true.

    this%right_shift(block_index,s,:) = data_right
    this%right_accumulator(block_index,:) = this%right_accumulator(block_index,:) + data_right
    this%shift_not_null(block_index,s,2) = .true.

    this%count_accumulated(block_index) = this%count_accumulated(block_index) + 1

    ! check if we need to move down a block
    If (this%count_accumulated(block_index) == this%window_size) Then
      Call this%add(this%left_accumulator(block_index,:) / this%window_size, &
          this%right_accumulator(block_index,:) / this%window_size, &
          block_index + 1)
      ! the data at this block may be reset
      this%left_accumulator(block_index,:) = 0
      this%right_accumulator(block_index,:) = 0
      this%count_accumulated(block_index) = 0    
    End If

    ! now accumulate the correlation
    i = this%shift_index(block_index)
    If (block_index == 1) Then
      j = i
      Do n = 1, this%points_per_block
        flag = .false.
        If ( this%shift_not_null(block_index,i,1) &
          .and. this%shift_not_null(block_index,j,2) )  Then
          Do l = 1, this%left_dim
            Do r = 1, this%right_dim
              this%correlation(block_index,n,l,r) = this%correlation(block_index,n,l,r) + &
                this%left_shift(block_index,i,l)*this%right_shift(block_index,j,r) 
              flag = .true.
            End Do
          End Do
        End If
        If (flag) Then 
          this%count_correlated(block_index,n) = this%count_correlated(block_index,n) + 1
        End If
        j = j - 1
        If (j <= 0) Then 
          j = this%points_per_block
        End If
      End Do
    Else
      j = i - this%min_dist
      Do n = this%min_dist, this%points_per_block 
        If (j <= 0) Then
          j = this%points_per_block
        End If
        flag = .false.
        If ( this%shift_not_null(block_index,i,1) &
          .and. this%shift_not_null(block_index,j,2) )  Then
          Do l = 1, this%left_dim
            Do r = 1, this%right_dim
              this%correlation(block_index,n,l,r) = this%correlation(block_index,n,l,r) + &
                this%left_shift(block_index,i,l)*this%right_shift(block_index,j,r) 
              flag = .true.
            End Do
          End Do
        End If
        If (flag) Then 
          this%count_correlated(block_index,n) = this%count_correlated(block_index,n) + 1
        End If
        j = j - 1
      End Do
    End If

    this%shift_index(block_index) = this%shift_index(block_index) + 1
    If (this%shift_index(block_index) == this%points_per_block + 1) Then
      this%shift_index(block_index) = 1
    End If

  End Subroutine add

  Subroutine update(this, data_left, data_right)
    Class(correlator), Intent(InOut) :: this
    Real(Kind=wp),     Intent(InOut) :: data_left(:), data_right(:)

    this%count_updated = this%count_updated + 1

    Call add(this,data_left,data_right,1)

  End Subroutine update

  Subroutine allocate_correlator_arrays(this, number_of_blocks, &
                                        points_per_block, window_size, dim_left, dim_right)
    Class(correlator),      Intent(InOut) :: this
    Integer,                Intent(In   ) :: number_of_blocks, points_per_block,&
                                              window_size, dim_left, dim_right
    Integer, Dimension(1:9)               :: fails
    
    this%number_of_blocks = number_of_blocks
    this%left_dim = dim_left
    this%right_dim = dim_right
    this%window_size = window_size
    this%points_per_block = points_per_block
    this%min_dist = points_per_block / window_size
    this%max_block_used = 0

    Allocate(this%left_accumulator(number_of_blocks,dim_left), Stat=fails(1))
    Allocate(this%right_accumulator(number_of_blocks,dim_right), Stat=fails(2))
    Allocate(this%count_accumulated(number_of_blocks), Stat=fails(3))
    Allocate(this%left_shift(number_of_blocks,points_per_block,dim_left), Stat=fails(4))
    Allocate(this%shift_not_null(number_of_blocks,points_per_block,2), Stat=fails(5))
    Allocate(this%right_shift(number_of_blocks,points_per_block,dim_right), Stat=fails(6))
    Allocate(this%shift_index(number_of_blocks), Stat=fails(7))
    Allocate(this%correlation(number_of_blocks,points_per_block,dim_left,dim_right), Stat=fails(8))
    Allocate(this%count_correlated(number_of_blocks,points_per_block), Stat=fails(9))

    this%buffer_size =    number_of_blocks*dim_left                             &
                        + number_of_blocks*dim_right                            &
                        + number_of_blocks                                      &
                        + number_of_blocks*points_per_block*dim_left            &
                        + number_of_blocks*points_per_block*2                   &
                        + number_of_blocks*points_per_block*dim_right           &
                        + number_of_blocks                                      &
                        + number_of_blocks*points_per_block*dim_left*dim_right  &
                        + number_of_blocks*points_per_block+3                    

    this%count_accumulated = 0
    this%left_accumulator = 0
    this%right_accumulator = 0
    this%left_shift = 0
    this%right_shift = 0
    this%shift_not_null = .false.
    this%shift_index = 1
    this%correlation = 0
    this%count_correlated = 0
    this%count_updated = 0

    If (Any(fails > 0)) Call error_alloc("allocate_correlator_arrays","correlator")

  End Subroutine allocate_correlator_arrays

  Subroutine deport_buffer(this, buffer, buffer_index, free)
    Class(correlator),             Intent(InOut) :: this
    Real(Kind=wp), Dimension(:),   Intent(InOut) :: buffer
    Integer,                       Intent(InOut) :: buffer_index
    Integer                                      :: from_index, to_index, elements, i
    Logical, Optional,             Intent(In   ) :: free  

    
    ! pack a buffer to be sent to another process 

    ! left accumulator
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%left_dim
    to_index = from_index+elements - 1
    buffer(from_index:to_index) = Reshape(this%left_accumulator(:,:),(/elements/))
    buffer_index = to_index

    ! right accumulator
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%right_dim 
    to_index = from_index+elements- 1
    buffer(from_index:to_index) = Reshape(this%right_accumulator(:,:),(/elements/))
    buffer_index = to_index

    ! count accumulated
    from_index = buffer_index + 1
    elements = this%number_of_blocks
    to_index = from_index + elements - 1
    buffer(from_index:to_index) = this%count_accumulated(:)
    buffer_index = to_index

    ! left shift
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%points_per_block*this%left_dim
    to_index = from_index + elements - 1
    buffer(from_index:to_index) = Reshape(this%left_shift(:,:,:),(/elements/))
    buffer_index = to_index

    ! right shift
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%points_per_block*this%right_dim 
    to_index = from_index+ elements - 1
    buffer(from_index:to_index) = Reshape(this%right_shift(:,:,:),(/elements/))
    buffer_index = to_index

    ! shift_not_null
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%points_per_block*2
    to_index = from_index + elements - 1
    buffer(from_index:to_index) = 0
    Where(Reshape(this%shift_not_null(:,:,:),(/elements/))) &
      buffer(from_index:to_index) = 1
    buffer_index = to_index

    ! shift index
    from_index = buffer_index + 1
    elements = this%number_of_blocks
    to_index = from_index+ elements - 1
    buffer(from_index:to_index) = this%shift_index(:)
    buffer_index = to_index

    ! correlation
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%points_per_block*this%left_dim*this%right_dim
    to_index = from_index+ elements - 1
    buffer(from_index:to_index) = Reshape(this%correlation(:,:,:,:),(/elements/))
    buffer_index = to_index

    ! count correlated
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%points_per_block 
    to_index = from_index + elements - 1
    buffer(from_index:to_index) = Reshape(this%count_correlated(:,:),(/elements/))
    buffer_index = to_index

    If (Present(free)) Then
      If(free) Then
        this%count_updated = 0
      End If
    End If

  End Subroutine deport_buffer

  Subroutine recieve_buffer(this, buffer, buffer_index)
    Class(correlator),             Intent(InOut)  :: this
    Real(Kind=wp), Dimension(:),   Intent(InOut)  :: buffer
    Integer,                       Intent(InOut)  :: buffer_index
    Integer                                       :: from_index, to_index, s

    ! left accumulator
    from_index = buffer_index + 1
    to_index = from_index + this%number_of_blocks*this%left_dim - 1
    this%left_accumulator(:,:) = Reshape(buffer(from_index:to_index) &
      ,(/this%number_of_blocks,this%left_dim/))
    buffer_index = to_index

    ! right accumulator
    from_index = buffer_index + 1
    to_index = from_index + this%number_of_blocks*this%right_dim - 1
    this%right_accumulator(:,:) = Reshape(buffer(from_index:to_index) &
      ,(/this%number_of_blocks,this%right_dim/))
    buffer_index = to_index

    ! count accumulated
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks - 1
    this%count_accumulated(:) = Floor(buffer(from_index:to_index))
    buffer_index = to_index

    ! left shift
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block*this%left_dim - 1
    this%left_shift(:,:,:) = Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block,this%left_dim/))
    buffer_index = to_index

    ! right shift
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block*this%right_dim - 1
    this%right_shift(:,:,:) = Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block,this%right_dim/))
    buffer_index = to_index

    ! shift_not_null
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block*2 - 1
    this%shift_not_null(:,:,:) = .false.
    Where(Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block,2/))>0) &
      this%shift_not_null(:,:,:) = .true.
    buffer_index = to_index

    ! shift index
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks - 1
    this%shift_index(:) = buffer(from_index:to_index)
    buffer_index = to_index

    ! correlation
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block*this%left_dim*this%right_dim - 1
    this%correlation(:,:,:,:) = Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block,this%left_dim,this%right_dim/))
    buffer_index = to_index

    ! count correlated
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block - 1
    this%count_correlated(:,:) = Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block/))
    buffer_index = to_index

    this%count_updated = 1

  End Subroutine recieve_buffer

  Subroutine cleanup(this)
    Type(correlator), Intent(InOut) :: this
    
    If (Allocated(this%left_accumulator)) Then
      Deallocate (this%left_accumulator)
    End If
    If (Allocated(this%right_accumulator)) Then
      Deallocate (this%right_accumulator)
    End If
    If (Allocated(this%count_accumulated)) Then
      Deallocate (this%count_accumulated)
    End If
    If (Allocated(this%shift_not_null)) Then
      Deallocate (this%shift_not_null)
    End If
    If (Allocated(this%left_shift)) Then
      Deallocate (this%left_shift)
    End If
    If (Allocated(this%right_shift)) Then
      Deallocate (this%right_shift)
    End If
    If (Allocated(this%shift_index)) Then
      Deallocate (this%shift_index)
    End If
    If (Allocated(this%correlation)) Then
      Deallocate (this%correlation)
    End If
    If (Allocated(this%count_correlated)) Then
      Deallocate (this%count_correlated)
    End If

  End Subroutine cleanup

End Module correlators