Module test_parse

    Use asserts,         Only: assert
    Use kinds,           Only: li
    Use parse,           Only: word_2_integer

    Implicit None

  Contains

    Subroutine run_parse_tests(passed_all)

        Logical, Intent(InOut)            :: passed_all

        Integer           :: i
        Integer(Kind=li)  :: bigint
        Character(Len=19) :: iword

        bigint = 1
        Do i = 0, 10
            bigint = ishft(bigint, i)
            Write (iword, '(i19)') bigint
            Call assert(bigint == word_2_integer(iword), "word_2_integer fail", passed_accum = passed_all)
        End Do

    End Subroutine run_parse_tests

End Module test_parse