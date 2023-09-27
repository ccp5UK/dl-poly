Module test_integrators

    Use asserts,         Only: assert
    Use kinds,           Only: wp
    Use integrators,     Only: trapezium_rule, simpsons_rule, integrator, integrand, integrand_d
  
    Implicit None
  
    Integer,       Parameter :: N = 1024
    Real(Kind=wp), Parameter :: dt = 0.01

  Contains
  
    Subroutine run_integrators_tests(passed_all)
      
        Logical, Intent(InOut)            :: passed_all

        Real(Kind=wp), Allocatable        :: t(:), f(:), tnu(:), fnu(:)
        Real(Kind=wp)                     :: actual, expected, dtau
        Integer                           :: i, k
        Class(integrator), Allocatable    :: inter
        Procedure (integrand_d), Pointer  :: func => NULL()

        Allocate(t(1:N))
        Allocate(f(1:N))

        Do i = 1, N
            t(i) = (Real(i)-1.0_wp) * dt
            f(i) = Cos(t(i))
        End Do

        Allocate(tnu(1:N))
        Allocate(fnu(1:N))

        dtau = 0.0_wp
        tnu(1) = 0.0_wp
        fnu(1) = Cos(tnu(1))
        Do i = 2, N
            ! non uniform pattern
            k = Mod(i*237, 1039)
            k = Mod(k*23 , 1049)
            dtau = (Real(Mod(k,100))/100.0)*(0.02-0.005)+0.005
            tnu(i) = tnu(i-1) + dtau
            fnu(i) = Cos(tnu(i))
        End Do

        expected = Sin(t(Size(t)))

        Allocate(trapezium_rule::inter)

        actual = inter%integrate_uniform(f, dt)

        Call assert(actual, expected,"trapezium rule (uniform) fail", passed_accum = passed_all, tolerance = 1e-3_wp)

        actual = inter%integrate(f, t)

        Call assert(actual, expected, "trapezium rule (non-uniform) on uniform data fail", &
            passed_accum = passed_all, tolerance = 1e-3_wp)

        expected = Sin(tnu(Size(tnu)))

        actual = inter%integrate(fnu, tnu)

        Call assert(actual, expected, "trapezium rule (non-uniform) on non-uniform data fail", &
            passed_accum = passed_all, tolerance = 1e-3_wp)

        expected = DSin(t(Size(t)))

        func => DCos

        actual = inter%integrate_func_d(func, t)

        Call assert(actual, expected, "trapezium rule (func) fail", passed_accum = passed_all, tolerance = 1e-3_wp)

        Deallocate(inter)

        Allocate(simpsons_rule::inter)

        expected = Sin(t(Size(t)))

        actual = inter%integrate_uniform(f, dt)

        Call assert(actual, expected, "simpsons rule (uniform) fail", passed_accum = passed_all, tolerance = 1e-3_wp)
   
        expected = Sin(t(Size(t)-1))

        actual = inter%integrate_uniform(f(1:Size(f)-1), dt)

        Call assert(actual, expected, "simpsons rule (uniform) even intervals fail", passed_accum = passed_all, tolerance = 1e-3_wp)
   
        expected = Sin(t(Size(t)))

        actual = inter%integrate(f, t)

        Call assert(actual, expected, "simpsons rule (non-uniform) on uniform data fail", &
            passed_accum = passed_all, tolerance = 1e-3_wp)

        expected = Sin(tnu(Size(tnu)))

        actual = inter%integrate(fnu, tnu)

        Call assert(actual, expected, "simpsons rule (non-uniform) on non-uniform data fail", &
            passed_accum = passed_all, tolerance = 1e-3_wp)

        Call assert( &
            Abs(inter%integrate_uniform(fnu, dt)-expected) > &
            Abs(inter%integrate(fnu, tnu)-expected), &
            "simpsons rule (non uniform) better on non-uniform data fail", passed_accum = passed_all)

        expected = DSin(t(Size(t)))

        func => DCos

        actual = inter%integrate_func_d(func, t)

        Call assert(actual, expected, "simpsons rule (func) fail", passed_accum = passed_all, tolerance = 1e-3_wp)

        Deallocate(inter)


    End Subroutine run_integrators_tests

End Module test_integrators