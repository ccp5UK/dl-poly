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
  Use errors_warnings, only : error_alloc, error_dealloc, error
  Use kinds, only : wp
  Implicit None

  Private

  Integer, Parameter :: STR_LEN = 256
  Integer, Parameter :: MAX_KEY = 50
  Character(Len=*), Parameter :: BAD_VAL = "VAL_NOT_IN_KEYS"

  Type, Public :: container
     !! Generic data container
     Private
     Class( * ), Pointer, Private :: data => Null()
   Contains
!     Generic, Public :: Assignment( = ) => set, get
     Procedure,            Private :: set => set_container
     Procedure, Pass( C ), Private :: get => get_container
  End type container

  Type, Public :: hash_table
     !! Type containing hash table of parameters
     Private
     Type(container), dimension(:), allocatable :: table_data
     Character(Len=MAX_KEY), dimension(:), allocatable :: table_keys
     Character(Len=MAX_KEY), dimension(:), allocatable :: key_names
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
     Generic  , Public  :: get => get_int, get_double, get_complex
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

  Public :: MAX_KEY, STR_LEN
  Public :: get_int, get_double, get_complex

Contains

  Subroutine set_container( C, stuff )
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for setting generic data container
    !!
    !! copyright - daresbury laboratory
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( container ), Intent( InOut ) :: C
    Class( *         ), Intent( In    ) :: stuff

    if (associated(C%data)) Deallocate(C%data)
    Allocate(C%data, source=stuff)

  End Subroutine set_container

  Subroutine get_container( stuff, C )
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for getting generic data container
    !!
    !! copyright - daresbury laboratory
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( container ),              Intent( In    ) :: C
    Class( *         ), Pointer, Intent(   Out ) :: stuff

    Allocate(stuff, source=C%data)

  End Subroutine get_container

  Subroutine cleanup(table)
    !!-----------------------------------------------------------------------
    !!
    !! dl_poly_4 subroutine for deallocation of hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Type(hash_table), Intent( InOut ) :: table
    Integer :: ierr
    Integer :: i

    do i = 1, table%size
      if (table%table_keys(i) /= BAD_VAL) then
        Deallocate(table%table_data(i)%data, stat=ierr)
        Nullify(table%table_data(i)%data)
        If (ierr /= 0) call error_dealloc("hash%table_data element", "cleanup hash table")
      end if
    end do

    if (allocated(table%table_data)) then
       Deallocate(table%table_data, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%table_data", "cleanup hash table")
    end if

    if (allocated(table%key_names)) then
       Deallocate(table%key_names, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%key_names", "cleanup hash table")
    end if

    if (allocated(table%table_keys)) then
       Deallocate(table%table_keys, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%table_keys", "cleanup hash table")
    end if

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
    Class(hash_table), Intent( InOut ) :: table
    Integer :: ierr
    Integer :: i

    do i = 1, table%size
      if (table%table_keys(i) /= BAD_VAL) then
        Deallocate(table%table_data(i)%data, stat=ierr)
        Nullify(table%table_data(i)%data)
        If (ierr /= 0) call error_dealloc("hash%table_data element", "cleanup hash table")
      end if
    end do

    if (allocated(table%table_data)) then
       Deallocate(table%table_data, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%table_data", "cleanup hash table")
    end if

    if (allocated(table%key_names)) then
       Deallocate(table%key_names, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%key_names", "cleanup hash table")
    end if

    if (allocated(table%table_keys)) then
       Deallocate(table%table_keys, stat=ierr)
       If (ierr /= 0) call error_dealloc("hash%table_keys", "cleanup hash table")
    end if

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

    Class(hash_table), Intent( InOut ) :: table

    !> Number of buckets to allocate
    Integer, Intent( In    ) :: nbuckets
    Logical, Optional :: can_overwrite
    Integer :: ierr

    if (table%allocated) then
      call table%resize(nbuckets)
      return
    end if

    if (present(can_overwrite)) then
       table%can_overwrite = can_overwrite
    end if

    table%size = nbuckets
    table%used_keys = 0

    Allocate(table%table_data(nbuckets), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%table_data", "allocate_hash_table")
    Allocate(table%key_names(nbuckets), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%key_names", "allocate_hash_table")
    Allocate(table%table_keys(nbuckets), stat=ierr)
    If (ierr /= 0) call error_alloc("hash%table_keys", "allocate_hash_table")

    table%table_keys(:) = BAD_VAL
    table%collisions = 0
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
    Class(hash_table), Intent( In    ) :: table
    Character(Len=*),  Intent( In    ) :: input
    Integer :: output
    Integer :: i

    if (input == BAD_VAL) call error(0, "Cannot hash value: "//BAD_VAL)
    output = 0
    do i = 1, len_trim(input)
       output = output + ichar(input(i:i))
    end do
    output = mod(output, table%size) + 1

  End Function hash_value

  Function get_loc(table, input, must_find) result(location)
    !!-----------------------------------------------------------------------
    !!
    !! Find location of input or bad_val if not found
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table) :: table
    Character(Len=*), Intent(In) :: input
    Logical, Intent(In), Optional :: must_find
    Integer :: location
    Character(Len=MAX_KEY) :: key


    location = table%hash(input)

    key = table%table_keys(location)
    ! Handle open addressing
    do while (trim(key) /= trim(input))
       if (key == BAD_VAL) then
          exit
       end if
       table%collisions = table%collisions + 1
       location = mod(location + 1, table%size)
       key = table%table_keys(location)
    end do

    if (present(must_find)) then
       if (must_find .and. key == BAD_VAL) call error(0, 'Key '//input//' not found in table')
    end if

  End Function get_loc

  Function contains_value(table, input) result(output)
    !!-----------------------------------------------------------------------
    !!
    !! Retrieve stored value from hash table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In     ) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: location
    Logical :: output

    if (.not. table%allocated) call error(0, 'Attempting to get from unallocated table')

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
    Class(hash_table), Intent( In     ) :: table
    Character(Len=*), Intent( In    ) :: input
    Integer :: location
    Class(*), Intent( In    ), Optional :: default
    Class(*), Intent(   Out ), Pointer :: output

    if (.not. table%allocated) call error(0, 'Attempting to get from unallocated table')

    location = table%get_loc(input)

    if (table%table_keys(location) == input) then
      call table%table_data(location)%get(output)
    else
      if (present(default)) then
        Allocate(output, source = default)
      else
        call error(0, 'No data pertaining to key '//trim(input)//' in table')
      end if
    end if

  End Subroutine get_hash_value

  Subroutine set_hash_value(table, key, input)
    !!-----------------------------------------------------------------------
    !!
    !! Set table at key to input
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( InOut ) :: table
    Character(Len=*), Intent( In    ) :: key
    Class(*), Intent( In    ) :: input
    Integer :: location

    if (.not. table%allocated) call error(0, 'Attempting to set unallocated table')

    location = table%get_loc(key)
    if (.not. table%can_overwrite .and. table%table_keys(location) /= BAD_VAL) then
       call error(0, 'Cannot overwrite key '//key)
    end if

    if (.not. table%in(key)) then
       table%used_keys = table%used_keys + 1
       table%key_names(table%used_keys) = key
    end if
    call table%table_data(location)%set(input)
    table%table_keys(location) = key

  End Subroutine set_hash_value

  Subroutine get_keys(table, keys)
    Class(hash_table), Intent( In    ) :: table
    Character(Len=MAX_KEY), Dimension(:), Allocatable :: keys

    if (allocated(keys)) then
       deallocate(keys)
    end if
    allocate(keys, source=table%key_names)

  End Subroutine get_keys

  Subroutine print_keys(table)
    !!-----------------------------------------------------------------------
    !!
    !! Print all keys in table
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In    ) :: table
    Integer :: i

    do i = 1, table%used_keys
       print('(A)'), table%key_names(i)
    end do

  End Subroutine print_keys

  Subroutine fill_from_table(table_from, table_to)
    !!-----------------------------------------------------------------------
    !!
    !! Populate table_to with keyvals from table_from
    !!
    !! copyright - daresbury laboratory
    !! author    - j.wilkins march 2020
    !!-----------------------------------------------------------------------
    Class(hash_table), Intent( In    ) :: table_from
    Class(hash_table), Intent( InOut ) :: table_to
    Integer :: location
    Integer :: i

    do i = 1, table_from%used_keys
       location = table_from%get_loc(table_from%key_names(i))
       call table_to%set(table_from%table_keys(location), table_from%table_data(location))
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
    Class(hash_table), Intent( InOut ) :: table
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
    Class(hash_table), Intent( In    ) :: table_from
    Class(hash_table), Intent( InOut ) :: table_to

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
    Class(hash_table), Intent( InOut ) :: table
    Class(hash_table), Pointer :: table_temp
    Integer, Intent( In    ) :: new_size
    Integer :: ierr

    if (new_size > table%size) then
       allocate(table_temp, stat=ierr)
       if (ierr.ne.0) call error_alloc('table_temp', 'expand_table')
       call table_temp%copy(table)
       call table%resize(new_size)
       call table%fill(table_temp)
       call cleanup(table_temp)
    else
       call error(0, 'New size for table less than previous size')
    end if

  End Subroutine expand_table

  Subroutine get_int( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get integer type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( hash_table ), Intent( In    ) :: table
    Character(Len=*), Intent( In    ) :: key
    Integer               , Intent(   Out ) :: val
    Integer, Intent( In    ), Optional :: default
    Class( * ), Pointer :: stuff

    call table%get_cont(key, default, stuff)

    Select Type( stuff )
    Type is ( Integer )
      val = stuff
    Class Default
       Call error(0, 'Trying to get integer from a not integer')
    End Select

    deallocate(stuff)
    nullify(stuff)

  End Subroutine get_int

  Subroutine get_double( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get double type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( hash_table ), Intent( In    ) :: table
    Character(Len=*), Intent( In    ) :: key
    Real(kind=wp), Intent( In    ), Optional :: default
    Real(kind=wp)               , Intent(   Out ) :: val
    Class( * ), Pointer :: stuff

    call table%get_cont(key, default, stuff)

    Select Type( stuff )
    Type is ( Real( wp ) )
      val = stuff
    Class Default
       Call error(0, 'Trying to get real from a not real')
    End Select

    deallocate(stuff)
    nullify(stuff)

  End Subroutine get_double

  Subroutine get_complex( table, key, val, default )
    !!-----------------------------------------------------------------------
    !!
    !! get complex type from hash table
    !!
    !! author    - i.j.bush april 2020
    !!-----------------------------------------------------------------------

    Implicit None

    Class( hash_table ), Intent( In    ) :: table
    Character(Len=*), Intent( In    ) :: key
    Complex(kind=wp)               , Intent(   Out ) :: val
    Complex(kind=wp), Intent( In    ), Optional :: default
    Class( * ), Pointer :: stuff

    call table%get_cont(key, default, stuff)

    Select Type( stuff )
    Type is ( Complex( wp ) )
      val = stuff
    Class Default
       Call error(0, 'Trying to get complex from a not complex')
    End Select

    deallocate(stuff)
    nullify(stuff)
  End Subroutine get_complex


  Subroutine fix_table(table)
    Class( hash_table ), Intent( InOut ) :: table

    table%can_overwrite = .false.
  end Subroutine fix_table


end Module hash
