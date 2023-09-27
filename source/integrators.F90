Module integrators

    Use kinds, Only: wp

    Use errors_warnings, Only: error

    Implicit None

    Private 

    Type, Abstract, Public :: integrator 
        Contains
            Procedure(integrate),         Deferred :: integrate
            Procedure(integrate_uniform), Deferred :: integrate_uniform 
            Procedure(integrate_func),    Deferred :: integrate_func
            ! without func and func_d there are hard to understand compile errors
            !  when using a Real(4) function 
            !  func => Cos
            !  val = inter%integrate_func(func, t)
            Procedure(integrate_func_d),  Deferred :: integrate_func_d
    End Type integrator 

    Abstract Interface  
        Function integrand(t) Result(v)
            Real, Intent(In   ) :: t
            Real                :: v
        End Function integrand

        Function integrand_d(t) Result(v)
            Import wp
            Real(Kind=wp), Intent(In   ) :: t
            Real(Kind=wp)                :: v
        End Function integrand_d

        Function integrate(i, f, t) Result(v)
            Import integrator, wp
            Class(integrator),              Intent(In   ) :: i
            Real(Kind=wp), Dimension(:),    Intent(In   ) :: f, t
            Real(Kind=wp)                                 :: v
        End Function integrate 

        Function integrate_uniform(i, f, dt) Result(v)
            Import integrator, wp
            Class(integrator),              Intent(In   ) :: i
            Real(Kind=wp), Dimension(:),    Intent(In   ) :: f
            Real(Kind=wp),                  Intent(In   ) :: dt
            Real(Kind=wp)                                 :: v
        End Function integrate_uniform


        Function integrate_func_d(i, f, t) Result(v)
            Import integrator, wp, integrand_d
            Class(integrator),                      Intent(In   ) :: i
            Procedure (integrand_d), Pointer,       Intent(In   ) :: f
            Real(Kind=wp), Dimension(:),            Intent(In   ) :: t
            Real(Kind=wp)                                         :: v
            
        End Function integrate_func_d

        Function integrate_func(i, f, t) Result(v)
            Import integrator, integrand
            Class(integrator),                    Intent(In   ) :: i
            Procedure (integrand), Pointer,       Intent(In   ) :: f
            Real, Dimension(:),                   Intent(In   ) :: t
            Real                                                :: v
            
        End Function integrate_func
    End Interface

    Type, Extends(integrator), Public :: trapezium_rule
        Contains 
            Procedure :: integrate => trapezium_rule_non_uniform
            Procedure :: integrate_uniform => trapezium_rule_uniform 
            Procedure :: integrate_func => trapezium_rule_func
            Procedure :: integrate_func_d => trapezium_rule_func_d
    End Type trapezium_rule

    Type, Extends(integrator), Public :: simpsons_rule
        Contains
            Procedure :: integrate => simpsons_rule_non_uniform
            Procedure :: integrate_uniform => simpsons_rule_uniform
            Procedure :: integrate_func => simpsons_rule_func
            Procedure :: integrate_func_d => simpsons_rule_func_d
    End Type

    Public :: integrand
    Public :: integrand_d

Contains 

    !!!!!!!!!!!!!!!!!!!! Trapezium Rule !!!!!!!!!!!!!!!!!!!!

    Function trapezium_rule_func_d(i, f, t) Result(v)
        Class(trapezium_rule),              Intent(In   ) :: i
        Procedure (integrand_d), Pointer,   Intent(In   ) :: f
        Real(Kind=wp), Dimension(:),        Intent(In   ) :: t
        Real(Kind=wp)                                     :: v, dt
        Integer                                           :: k

        v = 0.0_wp
        Do k = 2, Size(t)
            dt = t(k) - t(k-1)
            v = v + 0.5_wp * dt * (f(t(k))+f(t(k-1))) 
        End Do

    End Function trapezium_rule_func_d

    Function trapezium_rule_func(i, f, t) Result(v)
        Class(trapezium_rule),            Intent(In   ) :: i
        Procedure (integrand), Pointer,   Intent(In   ) :: f
        Real, Dimension(:),               Intent(In   ) :: t
        Real                                            :: v, dt
        Integer                                         :: k

        v = 0.0
        Do k = 2, Size(t)
            dt = t(k) - t(k-1)
            v = v + 0.5 * dt * (f(t(k))+f(t(k-1))) 
        End Do

    End Function trapezium_rule_func

    Function trapezium_rule_uniform(i, f, dt) Result(v)
        Class(trapezium_rule),          Intent(In   ) :: i
        Real(Kind=wp), Dimension(:),    Intent(In   ) :: f
        Real(Kind=wp),                  Intent(In   ) :: dt
        Real(Kind=wp)                                 :: v
        
        v = 0.0_wp
        If (Size(f) == 2) Then
            v = dt * (f(1) + f(2))*0.5_wp
        Else
            v = dt * (Sum(f(2:Size(f)-1)) + 0.5_wp * (f(1) + f(Size(f))))
        End If
    End Function trapezium_rule_uniform

    Function trapezium_rule_non_uniform(i, f, t) Result(v)
        Class(trapezium_rule),          Intent(In   ) :: i
        Real(Kind=wp), Dimension(:),    Intent(In   ) :: f, t
        Real(Kind=wp)                                 :: v, dt
        Integer                                       :: k

        If (Size(f) /= Size(t)) Then
            Call error(0, "unequal data sizes for integrator")
        End If
        
        v = 0.0_wp
        Do k = 2, Size(f)
            dt = t(k) - t(k-1)
            v = v + 0.5_wp * dt * (f(k)+f(k-1)) 
        End Do
    End Function trapezium_rule_non_uniform

    !!!!!!!!!!!!!!!!!!!! Simpson's Rule !!!!!!!!!!!!!!!!!!!!

    Function simpsons_rule_func_d(i, f, t) Result(v)
        Class(simpsons_rule),               Intent(In   ) :: i
        Procedure (integrand_d), Pointer,   Intent(In   ) :: f
        Real(Kind=wp), Dimension(:),        Intent(In   ) :: t
        Real(Kind=wp)                                     :: v, dt
        Integer                                           :: k, n

        n = Size(t)
        dt = (t(Size(t))-t(1)) / Real(n, Kind=wp)
        v = 0.0_wp
        If (Size(t) == 2) Then
            v = dt * (f(t(1)) + f(t(2)))*0.5_wp
        Else
            n = Size(t)-1
            If (Mod(Size(t),2) == 0) Then
                n = n - 3
            End If
            v = 0.0_wp
            Do k = 1, Int(Floor(Real(n,Kind=wp)*0.5_wp))
                v = v + f(t(2*k-1)) + 4.0_wp * f(t(2*k)) + f(t(2*k+1))
            End Do
            v = 1.0_wp/3.0_wp * v * dt
            If (Mod(Size(t),2 ) == 0) Then
                ! uneven interval correction (3/8th rule)
                n = Size(t)-1
                v = v + (3.0_wp/8.0_wp) * dt * ( f(t(n-3)) + 3.0_wp*f(t(n-2)) + 3.0_wp*f(t(n-1)) + f(t(n)))
            End If 
        End If
    End Function simpsons_rule_func_d

    Function simpsons_rule_func(i, f, t) Result(v)
        Class(simpsons_rule),             Intent(In   ) :: i
        Procedure (integrand), Pointer,   Intent(In   ) :: f
        Real, Dimension(:),               Intent(In   ) :: t
        Real                                            :: v, dt
        Integer                                         :: k, n

        dt = (t(Size(t))-t(1)) / Real(Size(t))
        v = 0.0
        If (Size(t) == 2) Then
            v = dt * (f(t(1)) + f(t(2)))*0.5
        Else
            n = Size(t)-1
            If (Mod(Size(t),2) == 0) Then
                n = n - 3
            End If
            v = 0.0
            Do k = 1, Int(Floor(Real(n))*0.5)
                v = v + f(t(2*k-1)) + 4.0 * f(t(2*k)) + f(t(2*k+1))
            End Do
            v = 1.0/3.0 * v * dt
            If (Mod(Size(t),2 ) == 0) Then
                ! uneven interval correction (3/8th rule)
                n = Size(t)-1
                v = v + (3.0/8.0) * dt * ( f(t(n-3)) + 3.0*f(t(n-2)) + 3.0*f(t(n-1)) + f(t(n)))
            End If 
        End If
    End Function simpsons_rule_func

    Function simpsons_rule_uniform(i, f, dt) Result(v)
        Class(simpsons_rule),            Intent(In   ) :: i
        Real(Kind=wp), Dimension(:),     Intent(In   ) :: f
        Real(Kind=wp),                   Intent(In   ) :: dt
        Real(Kind=wp)                                  :: v
        Integer                                        :: k, n
        
        v = 0.0_wp
        If (Size(f) == 2) Then
            v = dt * (f(1) + f(2))*0.5_wp
        Else
            n = Size(f)-1
            If (Mod(Size(f),2) == 0) Then
                n = n - 3
            End If
            v = 0.0_wp
            Do k = 1, Int(Floor(Real(n)*0.5_wp))
                v = v + f(2*k-1) + 4.0_wp * f(2*k) + f(2*k+1)
            End Do
            v = 1.0_wp/3.0_wp * v * dt

            If (Mod(Size(f),2 ) == 0) Then
                ! uneven interval correction (3/8th rule)
                n = Size(f)-1
                v = v + (3.0_wp/8.0_wp) * dt * ( f(n-3) + 3.0_wp*f(n-2) + 3.0_wp*f(n-1) + f(n))
            End If 
        End If
    End Function simpsons_rule_uniform

    Function simpsons_rule_non_uniform(i, f, t) Result(v)
        Class(simpsons_rule),           Intent(In   ) :: i
        Real(Kind=wp), Dimension(:),    Intent(In   ) :: f, t
        Real(Kind=wp)                                 :: v, dt, h0, h1, a, b, c
        Integer                                       :: k, n

        If (Size(f) /= Size(t)) Then
            Call error(0, "unequal data sizes for integrator")
        End If

        ! Shklov, N. (December 1960). 
        ! "Simpson's Rule for Unequally Spaced Ordinates". 
        ! The American Mathematical Monthly. 67 (10): 1022–1023. 
        ! doi:10.2307/2309244. JSTOR 2309244.
        
        n = Size(f)-1
        v = 0.0_wp
        If (Size(f) == 2) Then
            dt = t(2)-t(1)
            v = dt * (f(1) + f(2))*0.5_wp
        Else 
            Do k = 1, Int(Floor(Real(n)*0.5_wp))
                h0 = t(2*k)-t(2*k-1)
                h1 = t(2*k+1)-t(2*k)
                a = h0 + h1 
                b = h0 / h1 
                c = h0 * h1
                v = v + 1.0_wp/6.0_wp * a * ( (2.0_wp-1.0_wp/b)*f(2*k-1)+ a*a/c*f(2*k) + (2.0_wp-b)*f(2*k+1) )
            End Do

            If (Mod(n,2) == 1) Then 
                h0 = t(Size(t)-1)-t(Size(t)-2)
                h1 = t(Size(t))-t(Size(t)-1)
                v = v + f(Size(f))   * (2.0_wp * h1 ** 2.0_wp + 3.0_wp * h0 * h1) / (6.0_wp * (h0 + h1))
                v = v + f(Size(f)-1) * (h1 ** 2.0_wp + 3.0_wp * h1 * h0)          / (6.0_wp * h0)
                v = v - f(Size(f)-2) * h1 ** 3.0_wp                               / (6.0_wp * h0 * (h0 + h1))
            End If
        End If        
    End Function simpsons_rule_non_uniform

End Module integrators

        
