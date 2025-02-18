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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 Type for calculating Real and Complex correlations, \<AB\>.
    ! Uses the multi-tau algorithm, see add.
    !
    ! Internally both Real and Complex scalars are represented as
    ! Complex. This keeps code duplication minimal in lieu of
    ! templates, but incurs some overhead. This becomes important
    ! for very long lag time, single step accuracy per-atom correlations.
    !
    ! author    - h.l.devereux August 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !> Accumulator for A of \<AB\> for each block
      Complex(Kind=wp), Allocatable :: left_accumulator(:)
      !> Accumulator for B of \<AB\> for each block
      Complex(Kind=wp), Allocatable :: right_accumulator(:)
      !> Count accumulated observable for each block
      Integer(Kind=wi), Allocatable :: count_accumulated(:)
      !> Shift array for A, for each block
      Complex(Kind=wp), Allocatable :: left_shift(:,:)
      Logical,          Allocatable :: shift_not_null(:,:)
      !> Shift array for B, for each block
      Complex(Kind=wp), Allocatable :: right_shift(:,:)
      !> Current shift for each block
      Integer(Kind=wi), Allocatable :: shift_index(:)
      !> Value of correlation at block and point
      Complex(Kind=wp), Allocatable :: correlation(:,:)
      Integer(Kind=wi), Allocatable :: count_correlated(:,:)
      Integer(Kind=wi)              :: number_of_blocks = 0, &
                                       max_block_used = 0,   &
                                       window_size = 0,      &
                                       points_per_block = 0, &
                                       min_dist = 0,         &
                                       buffer_size = 0,      &
                                       count_updated = 0
  
    Contains
      Private
  
      Procedure, Public              :: init => allocate_correlator_arrays
      Procedure, Public              :: update
      Procedure, Private             :: add
      Generic,   Public              :: get_correlation => get_real_correlation, get_complex_correlation
      Procedure, Private             :: get_real_correlation, get_complex_correlation
      Procedure, Private             :: get_correlation_values
      Procedure, Public              :: deport_buffer
      Procedure, Public              :: recieve_buffer
      Final                          :: cleanup
  
    End Type

    Public :: correlator_buffer_size

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
    this%buffer = 0.0_wp
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
    this%buffer = 0
    Call this%mpi%init(comm)
  End Subroutine initialise_indices_buffer_type

  Subroutine finalise_indices_buffer_type(this)
    Class(indices_buffer_type), Intent(InOut) :: this

    If (Allocated(this%buffer)) Then
      Deallocate(this%buffer)
    End If

    Call this%mpi%finalise()
  End Subroutine finalise_indices_buffer_type

  Subroutine get_complex_correlation(this, correlation, timesteps, points_correlated)
    Class(correlator),                                     Intent(InOut)  :: this
    Complex(Kind=wp), Dimension(&
      1:this%points_per_block*this%number_of_blocks),      Intent(InOut)  :: correlation
    Real(Kind=wp), Dimension(&
      1:this%points_per_block*this%number_of_blocks),      Intent(InOut)  :: timesteps
    Integer,                                               Intent(  Out)  :: points_correlated

    Call get_correlation_values(this, correlation, timesteps, points_correlated)
  End Subroutine get_complex_correlation

  Subroutine get_real_correlation(this, correlation, timesteps, points_correlated)
    Class(correlator),                                     Intent(InOut)  :: this
    Real(Kind=wp), Dimension(& 
      1:this%points_per_block*this%number_of_blocks),      Intent(InOut)  :: correlation
    Real(Kind=wp), Dimension(&
      1:this%points_per_block*this%number_of_blocks),      Intent(InOut)  :: timesteps
    Integer,                                               Intent(  Out)  :: points_correlated

    Complex(Kind=wp), Dimension(&
      1:this%points_per_block*this%number_of_blocks) :: complex_correlation

    Call get_correlation_values(this, complex_correlation, timesteps, points_correlated)

    correlation = Real(Real(complex_correlation, Kind=wp), Kind=wp)

  End Subroutine get_real_correlation

  Subroutine get_correlation_values(this, correlation, timesteps, points_correlated)
    Class(correlator),                                     Intent(InOut)  :: this
    Complex(Kind=wp), Dimension(&
      1:this%points_per_block*this%number_of_blocks),      Intent(InOut)  :: correlation
    Real(Kind=wp), Dimension(&
      1:this%points_per_block*this%number_of_blocks),      Intent(InOut)  :: timesteps
    Integer,                                               Intent(  Out)  :: points_correlated
    Integer                                                               :: im, i, k
    
    correlation = 0.0_wp
    timesteps = 0.0_wp

    im = 1
    Do i = 1,this%points_per_block
      If (this%count_correlated(1,i) > 0) Then
        correlation(im) = correlation(im) + this%correlation(1,i) / this%count_correlated(1,i)
        timesteps(im) = i - 1
        im = im + 1
      End If
    End Do
    Do k = 2,this%max_block_used
      Do i = this%min_dist+1,this%points_per_block
        If (this%count_correlated(k,i) > 0) Then
          correlation(im) = correlation(im) + this%correlation(k,i) / this%count_correlated(k,i)
          timesteps(im) = Real(i-1,Kind=wp) * (Real(this%window_size, Kind=wp)**(k-1))
          im = im + 1
        EndIf
      End Do
    End Do

    ! this can be less than blocks * points_per_block 
    !  depending on the averaging parameter
    points_correlated = im

  End Subroutine get_correlation_values

  Recursive Subroutine add(this, data_left, data_right, block_index)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for adding newly observed data to a correlator.
    ! follows the multi-tau algorithm propagating down the block hierachy.
    !
    ! author    - h.l.devereux August 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Class(correlator),   Intent(InOut)  :: this
    Complex(Kind=wp),    Intent(In   )  :: data_left, data_right
    Integer,             Intent(In   )  :: block_index
    Integer                             :: i, j, s, n
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
    this%left_shift(block_index,s) = data_left
    this%left_accumulator(block_index) = this%left_accumulator(block_index) + data_left
    this%shift_not_null(block_index,s) = .true.

    this%right_shift(block_index,s) = data_right
    this%right_accumulator(block_index) = this%right_accumulator(block_index) + data_right
    this%shift_not_null(block_index,s) = .true.

    this%count_accumulated(block_index) = this%count_accumulated(block_index) + 1

    ! check if we need to move down a block
    If (this%count_accumulated(block_index) == this%window_size) Then
      Call this%add(this%left_accumulator(block_index) / this%window_size, &
          this%right_accumulator(block_index) / this%window_size, &
          block_index + 1)
      ! the data at this block may be reset
      this%left_accumulator(block_index) = 0
      this%right_accumulator(block_index) = 0
      this%count_accumulated(block_index) = 0    
    End If

    ! now accumulate the correlation
    i = this%shift_index(block_index)
    If (block_index == 1) Then
      j = i
      Do n = 1, this%points_per_block
        flag = .false.
        If ( this%shift_not_null(block_index,i) &
          .and. this%shift_not_null(block_index,j) )  Then
          this%correlation(block_index,n) = this%correlation(block_index,n) + &
            this%left_shift(block_index,i)*Conjg(this%right_shift(block_index,j))
          flag = .true.
        End If
        If (flag) Then 
          this%count_correlated(block_index,n) = this%count_correlated(block_index,n) + 1
        End If
        j = j - 1
        If (j <= 0) Then 
          j = j + this%points_per_block
        End If
      End Do
    Else
      j = i - this%min_dist+1
      Do n = this%min_dist+1, this%points_per_block 
        If (j <= 0) Then
          j = j + this%points_per_block
        End If
        flag = .false.
        If ( this%shift_not_null(block_index,i) &
          .and. this%shift_not_null(block_index,j) )  Then
          this%correlation(block_index,n) = this%correlation(block_index,n) + &
            this%left_shift(block_index,i)*Conjg(this%right_shift(block_index,j))
          flag = .true.
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_5 subroutine for updating the correlator for newly
    ! observed Real or Complex scalar data. Reals are handled internally
    ! as Complex in leiu of templates.
    !
    ! author    - h.l.devereux August 2024
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Class(correlator), Intent(InOut) :: this
    Complex(Kind=wp),  Intent(In   ) :: data_left, data_right

    this%count_updated = this%count_updated + 1
    Call add(this,data_left,data_right,1)

  End Subroutine update

  !> buffer size for deporting and recieving
  Integer Function correlator_buffer_size(blocks, points)
    Integer, Intent(In   ) :: blocks, points
    correlator_buffer_size = blocks*2  &
      + blocks*2                       &
      + blocks                         &
      + blocks*points*2                &
      + blocks*points                  &
      + blocks*points*2                &
      + blocks                         &
      + blocks*points*2                &
      + blocks*points+4
  End Function correlator_buffer_size

  Subroutine allocate_correlator_arrays(this, number_of_blocks, &
                                        points_per_block, window_size)
    Class(correlator),      Intent(InOut) :: this
    Integer,                Intent(In   ) :: number_of_blocks, points_per_block,&
                                             window_size
    Integer, Dimension(1:9)               :: fails
    
    this%number_of_blocks = number_of_blocks
    this%window_size = window_size
    this%points_per_block = points_per_block
    this%min_dist = points_per_block / window_size
    this%max_block_used = 0

    Allocate(this%left_accumulator(number_of_blocks), Stat=fails(1))
    Allocate(this%right_accumulator(number_of_blocks), Stat=fails(2))
    Allocate(this%count_accumulated(number_of_blocks), Stat=fails(3))
    Allocate(this%left_shift(number_of_blocks,points_per_block), Stat=fails(4))
    Allocate(this%shift_not_null(number_of_blocks,points_per_block), Stat=fails(5))
    Allocate(this%right_shift(number_of_blocks,points_per_block), Stat=fails(6))
    Allocate(this%shift_index(number_of_blocks), Stat=fails(7))
    Allocate(this%correlation(number_of_blocks,points_per_block), Stat=fails(8))
    Allocate(this%count_correlated(number_of_blocks,points_per_block), Stat=fails(9))

    this%buffer_size = correlator_buffer_size(&
      number_of_blocks, points_per_block)                  

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

  !> pack complex data into a real buffer (of 2x length)
  Subroutine pack_complex(dat, buffer, buffer_index)
    Complex(Kind=wp), Intent(In   ) :: dat(:)
    Real(Kind=wp),    Intent(InOut) :: buffer(:)
    Integer,          Intent(InOut) :: buffer_index

    Integer :: i

    Do i = 0, Size(dat)-1
      buffer(buffer_index+i*2+1) = Real(Real(dat(i+1), Kind=wp), Kind=wp)
      buffer(buffer_index+i*2+2) = Real(Aimag(dat(i+1)), Kind=wp)
    End Do
    buffer_index = buffer_index + Size(dat)*2
  End Subroutine pack_complex

  !> upack complex data from a real buffer (of 1/2 length)
  Subroutine unpack_complex(dat, buffer, buffer_index, count)
    Complex(Kind=wp), Allocatable, Intent(InOut) :: dat(:)
    Real(Kind=wp),                 Intent(In   ) :: buffer(:)
    Integer,                       Intent(InOut) :: buffer_index
    Integer,                       Intent(In   ) :: count

    Integer :: i

    If (Allocated(dat)) Deallocate(dat)
    Allocate(dat(1:count))

    Do i = 1, Size(dat)
      dat(i) = Cmplx(buffer(buffer_index+1), buffer(buffer_index+2), Kind=wp)
      buffer_index = buffer_index + 2
    End Do
  End Subroutine

  Subroutine deport_buffer(this, buffer, buffer_index, free)
    Class(correlator),             Intent(InOut) :: this
    Real(Kind=wp), Dimension(:),   Intent(InOut) :: buffer
    Integer,                       Intent(InOut) :: buffer_index
    Logical, Optional,             Intent(In   ) :: free  

    Integer :: from_index, to_index, elements
    
    ! pack a buffer to be sent to another process 

    ! left accumulator
    Call pack_complex(Reshape(this%left_accumulator(:), &
      (/this%number_of_blocks/)), buffer, buffer_index)

    ! right accumulator
    Call pack_complex(Reshape(this%right_accumulator(:), &
      (/this%number_of_blocks/)), buffer, buffer_index)

    !count accumulated
    from_index = buffer_index + 1
    elements = this%number_of_blocks
    to_index = from_index + elements - 1
    buffer(from_index:to_index) = this%count_accumulated(:)
    buffer_index = to_index

    ! left shift
    Call pack_complex(Reshape(this%left_shift(:,:), &
      (/this%number_of_blocks*this%points_per_block/)), &
      buffer, buffer_index)

    ! right shift
    Call pack_complex(Reshape(this%right_shift(:,:), &
      (/this%number_of_blocks*this%points_per_block/)), &
      buffer, buffer_index)

    ! shift_not_null
    from_index = buffer_index + 1
    elements = this%number_of_blocks*this%points_per_block
    to_index = from_index + elements - 1
    buffer(from_index:to_index) = 0
    Where(Reshape(this%shift_not_null(:,:),(/elements/))) &
      buffer(from_index:to_index) = 1
    buffer_index = to_index

    ! shift index
    from_index = buffer_index + 1
    elements = this%number_of_blocks
    to_index = from_index+ elements - 1
    buffer(from_index:to_index) = this%shift_index(:)
    buffer_index = to_index

    ! correlation
    Call pack_complex(Reshape(this%correlation(:,:), &
      (/this%number_of_blocks*this%points_per_block/)), &
      buffer, buffer_index)

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

    Integer                       :: from_index, to_index
    Complex(Kind=wp), Allocatable :: temp(:)

    ! left accumulator
    Call unpack_complex(temp, buffer, buffer_index, this%number_of_blocks)
    this%left_accumulator(:) = Reshape(temp, (/this%number_of_blocks/))

    ! right accumulator
    Call unpack_complex(temp, buffer, buffer_index, this%number_of_blocks)
    this%right_accumulator(:) = Reshape(temp, (/this%number_of_blocks/))

    ! count accumulated
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks - 1
    this%count_accumulated(:) = Floor(buffer(from_index:to_index))
    buffer_index = to_index

    ! left shift
    Call unpack_complex(temp, buffer, buffer_index, &
      this%number_of_blocks*this%points_per_block)
    this%left_shift(:,:) = Reshape(temp, &
      (/this%number_of_blocks,this%points_per_block/))

    ! right shift
    Call unpack_complex(temp, buffer, buffer_index, &
      this%number_of_blocks*this%points_per_block)
    this%right_shift(:,:) = Reshape(temp, &
      (/this%number_of_blocks,this%points_per_block/))

    ! shift_not_null
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block - 1
    this%shift_not_null(:,:) = .false.
    Where(Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block/))>0) &
      this%shift_not_null(:,:) = .true.
    buffer_index = to_index

    ! shift index
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks - 1
    this%shift_index(:) = INT(buffer(from_index:to_index))
    buffer_index = to_index

    ! correlation
    Call unpack_complex(temp, buffer, buffer_index, &
      this%number_of_blocks*this%points_per_block)
    this%correlation(:,:) = Reshape(temp, &
      (/this%number_of_blocks,this%points_per_block/))

    ! count correlated
    from_index = buffer_index + 1
    to_index = from_index+this%number_of_blocks*this%points_per_block - 1
    this%count_correlated(:,:) = INT(Reshape(buffer(from_index:to_index), &
      (/this%number_of_blocks,this%points_per_block/)))
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