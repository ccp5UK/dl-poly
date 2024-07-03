Module test_hash

  Use asserts,         Only: assert
  Use hash,            Only: hash_table

  Implicit None

  Integer, Parameter :: N = 10

  Contains

    Subroutine run_hash_tests(passed_all)
      Logical, Intent(InOut) :: passed_all

      Type(hash_table) :: table
      Integer          :: i, value
      Character(Len=2) :: key

      Call table%init(N)
      Call table%set("1", 1)
      Call table%get("1", value)
      Call assert(value, 1, "hash table get error", passed_accum = passed_all)

      Call table%get("2", value, 2)
      Call assert(value, 2, "hash table get error", passed_accum = passed_all)

      Do i = 2, N
        Write (key, '(i2)') i
        Call table%set(Trim(key), i)
        Call table%get(Trim(key), value, i)
        Call assert(value, i, "hash table get error", passed_accum = passed_all)
      End Do

      ! check getting nonexistant key on full table
      Call table%get("11", value, 11)
      Call assert(value, 11, "hash table get error", passed_accum = passed_all)
    End Subroutine run_hash_tests

End Module test_hash