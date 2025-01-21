Module test_numerics

    Use asserts,         Only: assert
    Use kinds,           Only: wp
    Use numerics,        Only: sarurnd,&
                               seed_type

    Implicit None

    ! Test vectors obtained from the SARU reference implementation https://doi.org/10.1016/j.cpc.2012.12.003
    Integer,          Parameter, Dimension(1:4) :: s1 = (/51399, 32399, 16674, 6896/), &
                                                   s2 = (/58554, 7739, 65494, 2073/), &
                                                   s3 = (/370, 46222, 27581, 35477/)
    ! Results of using s1's, s2's, s3's in Saru::Saru(...).d(0, 1), https://doi.org/10.1016/j.cpc.2012.12.003
    Real(Kind=wp),    Parameter, Dimension(1:4) :: saruref = (/0.79457369960855373_wp, &
                                                               0.92893400615812716_wp, &
                                                               0.91490455055465536_wp, &
                                                               0.037934117011169652_wp/)

  Contains

    Subroutine run_numerics_tests(passed_all)

        Logical, Intent(InOut)            :: passed_all

        Type(seed_type) :: s
        Integer         :: i

        ! Don't use seed offset
        s%defined = .false.
        Do i = 1, 4
            Call assert(sarurnd(s, s1(i),s2(i),s3(i)), saruref(i), &
              "SARU failed", passed_accum=passed_all, tolerance=1e-6_wp)
        End Do

    End Subroutine run_numerics_tests

End Module test_numerics