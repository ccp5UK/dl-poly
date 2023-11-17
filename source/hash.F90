Module hash
  !!-----------------------------------------------------------------------
  !!
  !! Module containing hash table routines for reading control input
  !! Hash table is a fixed-size hash table with open addressing
  !!
  !! copyright - daresbury laboratory
  !! author    - j.wilkins march 2020
  !! contributions - i.j.bush april 2020
  !!-----------------------------------------------------------------------
  Use errors_warnings, Only: error,&
                             error_alloc,&
                             error_dealloc
  Use kinds,           Only: wp

  Implicit None

  Private

  Integer, Parameter :: MAX_KEY = 50
  Character(Len=*), Parameter :: BAD_VAL = "VAL_NOT_IN_KEYS"

  Type, Public :: container
    !! Generic data container
    Private
    Class(*), Pointer, Private :: Data => Null()
  Contains
!     Generic, Public :: Assignment( = ) => set, get
    Procedure, Private :: set => set_container
    Procedure, Pass(C), Private :: get => get_container
  End Type container

  Type, Public :: hash_table
    !! Type containing hash table of parameters
    Private
    Type(container), Dimension(:), Allocatable :: table_data
    Character(Len=MAX_KEY), Dimension(:), Allocatable :: table_keys
    Character(Len=MAX_KEY), Dimension(:), Allocatable :: key_names
    Integer, Public :: collisions = 0
    Integer, Public :: used_keys = -1
    Integer :: size = -1
    !> Values in hash table can be overwritten: Default = True; Set with table%fix
    Logical :: can_overwrite = .true.
    Logical :: allocated = .false.
  Contains

    Private
    Procedure, Public, Pass :: init => allocate_hash_table
    Procedure, Public, Pass :: set => set_hash_value
    Generic, Public  :: get => get_int, get_double, get_complex
    Procedure, Public, Pass :: hash => hash_value
    Procedure, Public, Pass :: print_keys, get_keys
    Procedure, Public, Pass(table_to) :: fill => fill_from_table
    Procedure, Public, Pass :: copy => copy_table
    Procedure, Public, Pass :: resize => resize_table
    Procedure, Public, Pass :: expand => expand_table
    Procedure, Private :: get_int, get_double, get_complex
    Procedure, Public :: get_cont => get_hash_value
    Procedure, Private, Pass :: get_loc => get_loc
    Procedure, Public, Pass :: in => contains_value
    Procedure, Public, Pass :: fix => fix_table
    Procedure, Public, Pass :: destroy => destroy
    Final :: cleanup

  End Type hash_table

  Public :: MAX_KEY
  Public :: get_int, get_double, get_complex

Contains

  Subroutine set_container(C, stuff)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for setting generic data container
    !!
    !! copyright - daresbury laboratory
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Class(container), Intent(InOut) :: C
    Class(*),         Intent(In   ) :: stuff

    If (Associated(C%data)) Deallocate (C%data)
    Allocate (C%data, source=stuff)

  End Subroutine set_container

  Subroutine get_container(stuff, C)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for getting generic data container
    !!
    !! copyright - daresbury laboratory
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Class(*), Pointer, Intent(  Out) :: stuff
    Class(container),  Intent(In   ) :: C

    Allocate (stuff, source=C%data)

  End Subroutine get_container

  Subroutine cleanup(table)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for deallocation of hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Type(hash_table), Intent(InOut) :: table

    Integer :: i, ierr

    if (.not. table%allocated) return

    Do i = 1, table%size
      If (table%table_keys(i) /= BAD_VAL) Then
        Deallocate (table%table_data(i)%data, stat=ierr)
        Nullify (table%table_data(i)%data)
        If (ierr /= 0) Call error_dealloc("hash%table_data element", "cleanup hash table")
      End If
    End Do

    If (Allocated(table%table_data)) Then
      Deallocate (table%table_data, stat=ierr)
      If (ierr /= 0) Call error_dealloc("hash%table_data", "cleanup hash table")
    End If

    If (Allocated(table%key_names)) Then
      Deallocate (table%key_names, stat=ierr)
      If (ierr /= 0) Call error_dealloc("hash%key_names", "cleanup hash table")
    End If

    If (Allocated(table%table_keys)) Then
      Deallocate (table%table_keys, stat=ierr)
      If (ierr /= 0) Call error_dealloc("hash%table_keys", "cleanup hash table")
    End If

    table%size = -1
    table%used_keys = -1
    table%collisions = -1
    table%allocated = .false.

  End Subroutine cleanup

  Subroutine destroy(table)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for deallocation of hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(InOut) :: table

    Integer :: i, ierr

    Do i = 1, table%size
      If (table%table_keys(i) /= BAD_VAL) Then
        Deallocate (table%table_data(i)%data, stat=ierr)
        Nullify (table%table_data(i)%data)
        If (ierr /= 0) Call error_dealloc("hash%table_data element", "cleanup hash table")
      End If
    End Do

    If (Allocated(table%table_data)) Then
      Deallocate (table%table_data, stat=ierr)
      If (ierr /= 0) Call error_dealloc("hash%table_data", "cleanup hash table")
    End If

    If (Allocated(table%key_names)) Then
      Deallocate (table%key_names, stat=ierr)
      If (ierr /= 0) Call error_dealloc("hash%key_names", "cleanup hash table")
    End If

    If (Allocated(table%table_keys)) Then
      Deallocate (table%table_keys, stat=ierr)
      If (ierr /= 0) Call error_dealloc("hash%table_keys", "cleanup hash table")
    End If

    table%size = -1
    table%used_keys = -1
    table%collisions = -1
    table%allocated = .false.

  End Subroutine destroy

  Recursive Subroutine allocate_hash_table(table, nbuckets, can_overwrite)
    !!-----------------------------------------------------------------------
    !!
    !! Subroutine for allocation and initialisation of hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------

    Class(hash_table), Intent(InOut) :: table
    Integer,           Intent(In   ) :: nbuckets
    Logical, Optional                :: can_overwrite

    Integer :: ierr

!> Number of buckets to allocate

    If (table%allocated) Then
      Call table%resize(nbuckets)
      Return
    End If

    If (Present(can_overwrite)) Then
      table%can_overwrite = can_overwrite
    End If

    table%size = nbuckets
    table%used_keys = 0

    Allocate (table%table_data(nbuckets), stat=ierr)
    If (ierr /= 0) Call error_alloc("hash%table_data", "allocate_hash_table")
    Allocate (table%key_names(nbuckets), stat=ierr)
    If (ierr /= 0) Call error_alloc("hash%key_names", "allocate_hash_table")
    Allocate (table%table_keys(nbuckets), stat=ierr)
    If (ierr /= 0) Call error_alloc("hash%table_keys", "allocate_hash_table")

    table%table_keys(:) = BAD_VAL
    table%collisions = 0
    table%allocated = .true.

  End Subroutine allocate_hash_table

  Function hash_value(table, input) Result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Function to hash string using simple sum(ord(input))%max_hash
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(In   ) :: table
    Character(Len=*),  Intent(In   ) :: input
    Integer                          :: output

    Integer :: i

    If (input == BAD_VAL) Call error(0, "Cannot hash value: "//BAD_VAL)
    output = 0
    Do i = 1, Len_trim(input)
      output = output + Ichar(input(i:i))
    End Do
    output = Mod(output, table%size) + 1

  End Function hash_value

  Function get_loc(table, input, must_find) Result(location)
    !!-----------------------------------------------------------------------
    !!
    !! Find location of input or bad_val if not found
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table)                :: table
    Character(Len=*),  Intent(In   ) :: input
    Logical, Optional, Intent(In   ) :: must_find
    Integer                          :: location

    Character(Len=MAX_KEY) :: key

    location = table%hash(input)
    key = table%table_keys(location)
    ! Handle open addressing
    Do While (Trim(key) /= Trim(input))
      If (key == BAD_VAL) Then
        Exit
      End If
      table%collisions = table%collisions + 1
      location = location + 1
      If (location > table%size) location = 1 ! Loop back through safely
      key = table%table_keys(location)
    End Do

    If (Present(must_find)) Then
      If (must_find .and. key == BAD_VAL) Call error(0, 'Key '//input//' not found in table')
    End If

  End Function get_loc

  Function contains_value(table, input) Result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Retrieve stored value from hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(In   ) :: table
    Character(Len=*),  Intent(In   ) :: input
    Logical                          :: output

    Integer :: location

    If (.not. table%allocated) Call error(0, 'Attempting to get from unallocated table')

    location = table%get_loc(input)
    output = table%table_keys(location) == input

  End Function contains_value

  Subroutine get_hash_value(table, input, default, output)
    !!-----------------------------------------------------------------------
    !!
    !! Retrieve stored value from hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table),  Intent(In   ) :: table
    Character(Len=*),   Intent(In   ) :: input
    Class(*), Optional, Intent(In   ) :: default
    Class(*), Pointer,  Intent(  Out) :: output

    Integer :: location

    If (.not. table%allocated) Call error(0, 'Attempting to get from unallocated table')

    location = table%get_loc(input)

    If (table%table_keys(location) == input) Then
      Call table%table_data(location)%get(output)
    Else
      If (Present(default)) Then
        Allocate (output, source=default)
      Else
        Call error(0, 'No data pertaining to key '//Trim(input)//' in table')
      End If
    End If

  End Subroutine get_hash_value

  Subroutine set_hash_value(table, key, input)
    !!-----------------------------------------------------------------------
    !!
    !! Set table at key to input
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(InOut) :: table
    Character(Len=*),  Intent(In   ) :: key
    Class(*),          Intent(In   ) :: input

    Integer :: location

    If (.not. table%allocated) Call error(0, 'Attempting to set unallocated table')

    location = table%get_loc(key)
    If (.not. table%can_overwrite .and. table%table_keys(location) /= BAD_VAL) Then
      Call error(0, 'Cannot overwrite key '//key)
    End If

    If (.not. table%in(key)) Then
      table%used_keys = table%used_keys + 1
      table%key_names(table%used_keys) = key
    End If
    Call table%table_data(location)%set(input)
    table%table_keys(location) = key

  End Subroutine set_hash_value

  Subroutine get_keys(table, keys)
    Class(hash_table),                  Intent(In   ) :: table
    Character(Len=MAX_KEY), Allocatable, Dimension(:) :: keys

    If (Allocated(keys)) Then
      Deallocate (keys)
    End If
    Allocate (keys, source=table%key_names)

  End Subroutine get_keys

  Subroutine print_keys(table)
    !!-----------------------------------------------------------------------
    !!
    !! Print all keys in table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(In   ) :: table

    Integer :: i

    Do i = 1, table%used_keys
      Print('(A)'), table%key_names(i)
    End Do

  End Subroutine print_keys

  Subroutine fill_from_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Populate table_to with keyvals from table_from
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(In   ) :: table_from
    Class(hash_table), Intent(InOut) :: table_to

    Integer :: i, location

    Do i = 1, table_from%used_keys
      location = table_from%get_loc(table_from%key_names(i))
      Call table_to%set(table_from%table_keys(location), table_from%table_data(location))
    End Do

  End Subroutine fill_from_table

  Subroutine resize_table(table, new_size)
    !!-----------------------------------------------------------------------
    !!
    !! Make resize table to new_size or empty table if not provided (destroys table)
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(InOut) :: table
    Integer, Optional, Intent(In   ) :: new_size

    Integer :: table_size

    If (Present(new_size)) Then
      table_size = new_size
    Else
      table_size = table%size
    End If

    Call cleanup(table)
    Call table%init(table_size)

  End Subroutine resize_table

  Subroutine copy_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Make table_to a copy of table_from (destroys table_to)
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(In   ) :: table_from
    Class(hash_table), Intent(InOut) :: table_to

    Call table_to%resize(table_from%size)
    Call table_to%fill(table_from)

  End Subroutine copy_table

  Subroutine expand_table(table, new_size)
    !!-----------------------------------------------------------------------
    !!
    !! Non-destructively resize table to new_size (must be larger than old size)
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent(InOut) :: table
    Integer,           Intent(In   ) :: new_size

    Class(hash_table), Pointer :: table_temp
    Integer                    :: ierr

    If (new_size > table%size) Then
      Allocate (table_temp, stat=ierr)
      If (ierr .ne. 0) Call error_alloc('table_temp', 'expand_table')
      Call table_temp%copy(table)
      Call table%resize(new_size)
      Call table%fill(table_temp)
      Call cleanup(table_temp)
    Else
      Call error(0, 'New size for table less than previous size')
    End If

  End Subroutine expand_table

  Subroutine get_int(table, key, val, default)
    !!-----------------------------------------------------------------------
    !!
    !! get integer type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Class(hash_table), Intent(In   ) :: table
    Character(Len=*),  Intent(In   ) :: key
    Integer,           Intent(  Out) :: val
    Integer, Optional, Intent(In   ) :: default

    Class(*), Pointer :: stuff

    Call table%get_cont(key, default, stuff)

    Select Type (stuff)
    Type is (Integer)
      val = stuff
    Class Default
      Call error(0, 'Trying to get integer from a not integer')
    End Select

    Deallocate (stuff)
    Nullify (stuff)

  End Subroutine get_int

  Subroutine get_double(table, key, val, default)
    !!-----------------------------------------------------------------------
    !!
    !! get double type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Class(hash_table),       Intent(In   ) :: table
    Character(Len=*),        Intent(In   ) :: key
    Real(kind=wp),           Intent(  Out) :: val
    Real(kind=wp), Optional, Intent(In   ) :: default

    Class(*), Pointer :: stuff

    Call table%get_cont(key, default, stuff)

    Select Type (stuff)
    Type is (Real (wp))
      val = stuff
    Class Default
      Call error(0, 'Trying to get real from a not real')
    End Select

    Deallocate (stuff)
    Nullify (stuff)

  End Subroutine get_double

  Subroutine get_complex(table, key, val, default)
    !!-----------------------------------------------------------------------
    !!
    !! get complex type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Class(hash_table),          Intent(In   ) :: table
    Character(Len=*),           Intent(In   ) :: key
    Complex(kind=wp),           Intent(  Out) :: val
    Complex(kind=wp), Optional, Intent(In   ) :: default

    Class(*), Pointer :: stuff

    Call table%get_cont(key, default, stuff)

    Select Type (stuff)
    Type is (Complex (wp))
      val = stuff
    Class Default
      Call error(0, 'Trying to get complex from a not complex')
    End Select

    Deallocate (stuff)
    Nullify (stuff)
  End Subroutine get_complex

  Subroutine fix_table(table)
    Class(hash_table), Intent(InOut) :: table

    table%can_overwrite = .false.
  End Subroutine fix_table

End Module hash
