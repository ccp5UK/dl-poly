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

  Integer, Parameter :: STR_LEN = 256
  Character(Len=*), Parameter :: BAD_VAL = "VAL_NOT_IN_KEYS"

  Type, Public :: control_parameter
     !! Type containing breakdown of control parameter
    Character(Len=STR_LEN) :: name
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

  Contains

    Private
    Procedure, Public, Pass :: init => allocate_hash_table
    Procedure, Public, Pass :: set => set_hash_value
    Procedure, Public, Pass :: get => get_hash_value
    Procedure, Public, Pass :: hash => hash_value
    Procedure, Public, Pass :: keys => print_keys
    Procedure, Public, Pass :: vals => print_vals
    Procedure, Public, Pass :: keyvals => print_keyvals
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

    deallocate(table%table_data, stat=ierr)
    If (ierr /= 0) call error_dealloc("hash%table_data", "cleanup hash table")
    deallocate(table%key_names, stat=ierr)
    If (ierr /= 0) call error_dealloc("hash%key_names", "cleanup hash table")
    table%size = -1
    table%used_keys = -1

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

    location = table%hash(input)

    output = table%table_data(location)
    ! Handle open addressing
    do while (output%name /= input)
      if (output%name == BAD_VAL) then
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

    location = table%hash(key)

    output = table%table_data(location)
    if (.not. table%can_overwrite .and. trim(output%name) == trim(key)) &
       call error(0, 'Cannot overwrite key '//trim(key)//' in hash table')

    ! Handle open addressing
    do while (output%name /= BAD_VAL)
       if (.not. table%can_overwrite .and. trim(output%name) == trim(key)) &
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

    do i = 1, table%used_keys
       val = table%get(table%key_names(i))
       print('(3(A,1X))'), trim(val%name), trim(val%val), trim(val%unit)
    end do

  End Subroutine print_keyvals


End Module hash
