Module hash

  !!-----------------------------------------------------------------------
  !!
  !! Module containing hash table routines for reading control input
  !! Hash table is a fixed-size hash table with open addressing
  !!
  !! copyright - daresbury laboratory
  !! author    - j.wilkins march 2020
  !!-----------------------------------------------------------------------
  Use errors_warnings, only : error, error_alloc, error_dealloc
  Implicit None

  Integer, Parameter :: STR_LEN = 256
  Character(Len=*), Parameter :: BAD_VAL = "VAL_NOT_IN_KEYS"

  Type, Public :: control_parameter
     !! Type containing breakdown of control parameter
     Character(Len=STR_LEN) :: key
     Character(Len=STR_LEN) :: val
     Character(Len=STR_LEN) :: unit
  End Type control_parameter

  Type, Public :: parameters_hash_table
     !! Type containing hash table of parameters
     Private
     Type(control_parameter), dimension(:), allocatable :: table_data
     Character(Len=STR_LEN), dimension(:), allocatable :: key_names
     Integer :: used_keys = -1
     Integer :: size = -1
     !> Values in hash table can be overwritten: Default = False
     Logical :: can_overwrite = .false.
     Logical :: allocated = .false.
   Contains

     Private
     Procedure, Public, Pass :: init => allocate_hash_table
     Procedure, Public, Pass :: set => set_hash_value
     Procedure, Public, Pass :: get => get_hash_value
     Procedure, Public, Pass :: hash => hash_value
     Procedure, Public, Pass :: keys => print_keys
     Procedure, Public, Pass :: vals => print_vals
     Procedure, Public, Pass :: keyvals => print_keyvals
     Procedure, Public, Pass(table_to) :: fill => fill_from_table
     Procedure, Public, Pass :: copy => copy_table
     Procedure, Public, Pass :: resize => resize_table
     Procedure, Public, Pass :: expand => expand_table
     Final :: cleanup

  End Type parameters_hash_table

Contains

  Subroutine cleanup(table)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for deallocation of hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Type(parameters_hash_table), Intent( InOut ) :: table
    Integer :: ierr

    if (allocated(table%table_data)) then
       Deallocate(table%table_data, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%table_data", "cleanup hash table")
    end if

    if (allocated(table%key_names)) then
       Deallocate(table%key_names, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%key_names", "cleanup hash table")
    end if

    table%size = -1
    table%used_keys = -1
    table%allocated = .false.

  End Subroutine cleanup

  Subroutine allocate_hash_table(table, size, can_overwrite)
    !!-----------------------------------------------------------------------
    !!
    !! Subroutine for allocation and initialisation of hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------

    Class(parameters_hash_table), Intent( InOut ) :: table

    !> Number of buckets to allocate
    Integer, Intent( In    ) :: size
    Logical, Optional :: can_overwrite
    Integer :: i

    Integer :: ierr

    if (present(can_overwrite)) then
       table%can_overwrite = can_overwrite
    end if

    table%size = size
    table%used_keys = 0
    Allocate(table%table_data(size), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%table_data", "allocate_hash_table")
    Allocate(table%key_names(size), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%key_names", "allocate_hash_table")

    do i = 1, size
       table%table_data(i) = control_parameter(BAD_VAL, "0", "None")
    end do
    table%allocated = .true.

  End Subroutine allocate_hash_table

  Function hash_value(table, input) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Function to hash string using simple sum(ord(input))%max_hash
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: output
    Integer :: i

    if (input == BAD_VAL) call error(0, "Cannot hash value: "//BAD_VAL)
    output = 0
    do i = 1, len_trim(input)
       output = output + ichar(input(i:i))
    end do
    output = mod(output, table%size)

  End Function hash_value

  Function get_hash_value(table, input) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Retrieve stored value from hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( In     ) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: location
    Type(control_parameter) :: output

    if (.not. table%allocated) call error(0, 'Attempting to get from unallocated table')

    location = table%hash(input)

    output = table%table_data(location)
    ! Handle open addressing
    do while (output%key /= input)
       if (output%key == BAD_VAL) then
          exit
       end if
       location = mod(location + 1, table%size)
       output = table%table_data(location)
    end do

  End Function get_hash_value

  Subroutine set_hash_value(table, key, input)
    !!-----------------------------------------------------------------------
    !!
    !! Set table at key to input
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type(control_parameter), Intent( In    ) :: input
    Type(control_parameter) :: output
    Integer :: location

    if (.not. table%allocated) call error(0, 'Attempting to set unallocated table')

    location = table%hash(key)

    output = table%table_data(location)
    if (.not. table%can_overwrite .and. trim(output%key) == trim(key)) &
         call error(0, 'Cannot overwrite key '//trim(key)//' in hash table')

    ! Handle open addressing
    do while (output%key /= BAD_VAL)
       if (.not. table%can_overwrite .and. trim(output%key) == trim(key)) &
            call error(0, 'Cannot overwrite key '//trim(key)//' in hash table')

       location = mod(location + 1, table%size)
       output = table%table_data(location)
    end do

    table%table_data(location) = input
    table%used_keys = table%used_keys + 1
    table%key_names(table%used_keys) = key

  End Subroutine set_hash_value

  Subroutine print_keys(table)
    !!-----------------------------------------------------------------------
    !!
    !! Print all keys in table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( In    ) :: table
    Integer :: i

    do i = 1, table%used_keys
       print('(A)'), table%key_names(i)
    end do

  End Subroutine print_keys

  Subroutine print_vals(table)
    !!-----------------------------------------------------------------------
    !!
    !! Print all values in table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( In    ) :: table
    Type(control_parameter) :: val
    Integer :: i

    do i = 1, table%used_keys
       val = table%get(table%key_names(i))
       print('(2(A,1X))'), val%val, val%unit
    end do

  End Subroutine print_vals

  Subroutine print_keyvals(table)
    !!-----------------------------------------------------------------------
    !!
    !! Print keys and values in table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( In    ) :: table
    Type(control_parameter) :: val
    Integer :: i

    do i = 1, table%used_keys
       val = table%get(table%key_names(i))
       print('(3(A,1X))'), trim(val%key), trim(val%val), trim(val%unit)
    end do

  End Subroutine print_keyvals

  Subroutine fill_from_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Populate table_to with keyvals from table_from
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( In    ) :: table_from
    Class(parameters_hash_table), Intent( InOut ) :: table_to
    Type(control_parameter) :: val
    Integer :: i

    do i = 1, table_from%used_keys
       val = table_from%get(table_from%key_names(i))
       call table_to%set(val%key, val)
    end do

  End Subroutine fill_from_table

  Subroutine resize_table(table, new_size)
    !!-----------------------------------------------------------------------
    !!
    !! Make resize table to new_size or empty table if not provided (destroys table)
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( InOut ) :: table
    Integer, Optional, Intent( In    ) :: new_size
    Integer :: table_size

    if (present(new_size)) then
       table_size = new_size
    else
       table_size = table%size
    end if

    call cleanup(table)
    call table%init(table_size)

  End Subroutine resize_table

  Subroutine copy_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Make table_to a copy of table_from (destroys table_to)
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( In    ) :: table_from
    Class(parameters_hash_table), Intent( InOut ) :: table_to

    call table_to%resize(table_from%size)
    call table_to%fill(table_from)

  End Subroutine copy_table

  Subroutine expand_table(table, new_size)
    !!-----------------------------------------------------------------------
    !!
    !! Non-destructively resize table to new_size (must be larger than old size)
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(parameters_hash_table), Intent( InOut ) :: table
    Class(parameters_hash_table), Pointer :: table_temp
    Integer, Intent( In    ) :: new_size

    if (new_size > table%size) then
       allocate(table_temp)
       call table_temp%copy(table)
       call table%resize(new_size)
       call table%fill(table_temp)
       call cleanup(table_temp)
    else
       call error(0, 'New size for table less than previous size')
    end if

  End Subroutine expand_table


End Module hash
