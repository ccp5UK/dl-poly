Module hash

  use errors_warnings, only : error, error_alloc, error_dealloc

  Integer, Parameter :: STR_LEN = 40
  Character(Len=*), Parameter :: BAD_VAL = "VAL_NOT_IN_KEYS"

  Type, Public :: control_parameter

    Character(Len=STR_LEN) :: name
    Character(Len=STR_LEN) :: val
    Character(Len=STR_LEN) :: unit

  End Type control_parameter

  Type, Public :: parameters_hash_table

    Private
    Type(control_parameter), dimension(:), allocatable :: table_data
    Character(Len=STR_LEN), dimension(:) :: keys
    Integer :: used_keys = -1
    Integer :: size = -1

  Contains

    Private
    Procedure, Public :: init => allocate_hash_table
    Procedure, Public :: set => set_hash_value
    Procedure, Public :: get => get_hash_value
    Procedure, Public :: hash => hash_value
    Final :: cleanup

  End Type parameters_hash_table

  Subroutine cleanup(self)
    Type(parameters_hash_table), Intent( InOut ) :: self
    Integer :: ierr

    deallocate(self%table_data, stat=ierr)
    If (ierr /= 0) call error_dealloc("hash%table_data", "cleanup hash table")
    deallocate(self%keys, stat=ierr)
    If (ierr /= 0) call error_dealloc("hash%keys", "cleanup hash table")
    self%size = -1
    self%used_keys = -1

  End Subroutine cleanup

  Subroutine allocate_hash_table(self, size)
    Type(parameters_hash_table), Intent( InOut ) :: self
    Integer, Intent( In    ) :: size
    Integer :: ierr

    self%size = size
    self%used_keys = 0
    Allocate(self%table_data(size), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%table_data", "allocate_hash_table")
    Allocate(self%keys(size), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%keys", "allocate_hash_table")

    do i = 1, size
      self%table_data(i) = control_parameter(BAD_VAL, "0", "None")
    end do

  End Subroutine allocate_hash_table

  Function hash_value(table, input) result(output)
    Type(parameters_hash_table) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: output

    if (input == BAD_VAL) call error("Cannot hash value: "//BAD_VAL, 0)

    do i = 1, len_trim(input)
      output = output + ichar(input(i:i))
    end do
    output = mod(output, table%size)

  End Function hash_value

  Function get_hash_value(table, input) result(output)
    Type(parameters_hash_table), Intent( In     ) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: location
    Type(control_parameter) :: output

    location = table%hash(input)

    output = table%table_data(location)
    ! Handle open addressing
    do while output%name /= input
      if (output%name == BAD_VAL) then
        exit
      end if
      location = mod(location + 1, table%size)
      output = table%table_data(location)
    end do

  End Function get_hash_value

  Subroutine set_hash_value(table, key, value)
    Type(parameters_hash_table), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Type(control_parameter), Intent( In    ) :: value

  End Subroutine set_hash_value


End Module hash
