Module two_body_potentials

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Module for defining an abstract potential
  !  Any potential implements and energy function that returns 
  !  a potential energy structure with two fields (energy/ gamma)
  !  this simplifies (makes branchless) core vdw code to act on a 
  !  potential a apposed to selecting dynamically with a case statement 
  !
  ! copyright - daresbury laboratory
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Use kinds,               Only: wp

    Use errors_warnings,     Only: error

    Use constants,           Only: r4pie0

    Implicit None 

    Private 

    Type, Public :: potential_energy
        Real(Kind=wp) :: energy, gamma
    End Type

    Type, Abstract, Public :: potential 
        Contains
            Procedure(energy),         Deferred :: energy
            Procedure(set_parameters), Deferred :: set_parameters
    End Type potential

    Type, Public :: potential_holder
        Class(potential), Allocatable :: p
        Contains 
            Procedure, Public :: energy => holder_energy
    End Type

    Abstract Interface
        Pure Type(potential_energy) Function energy(p, r)
            Import wp, potential_energy, potential
            Class(potential),       Intent(In   ) :: p
            Real(Kind=wp),          Intent(In   ) :: r
        End Function energy

        Subroutine set_parameters(p, param)
            Import wp, potential
            Class(potential),       Intent(InOut) :: p
            Real(Kind=wp),          Intent(In   ) :: param(:)
        End Subroutine set_parameters 
    End Interface

    Type, Extends(potential), Public :: LJ
        Real(Kind=wp) :: sigma, epsilon
        Contains 
            Procedure :: energy => lj_energy 
            Procedure :: set_parameters => lj_set_parameters
    End Type LJ

    Type, Extends(potential), Public :: lj_coh
        Real(Kind=wp) :: sigma, epsilon, coh
        Contains 
            Procedure :: energy => lj_coh_energy
            Procedure :: set_parameters => lj_coh_set_parameters
    End Type lj_coh

    Type, Extends(potential), Public :: LJ126
        Real(Kind=wp) :: a, b
        Contains 
            Procedure :: energy => LJ126_energy
            Procedure :: set_parameters => LJ126_set_parameters
    End Type LJ126

    Type, Extends(potential), Public :: n_m
        Real(Kind=wp) :: e0, n, m, r0
        Contains 
            Procedure :: energy => n_m_energy
            Procedure :: set_parameters => n_m_set_parameters
    End Type n_m

    Type, Extends(potential), Public :: nm_shift
        Real(Kind=wp) :: e0, n, m, r0, r_trunc
        Contains 
            Procedure :: energy => nm_shift_energy
            Procedure :: set_parameters => nm_shift_set_parameters
    End Type nm_shift

    Type, Extends(potential), Public :: morse
        Real(Kind=wp) :: e0, r0, k
        Contains 
            Procedure :: energy => morse_energy
            Procedure :: set_parameters => morse_set_parameters
    End Type morse

    Type, Extends(potential), Public :: morse12
        Real(Kind=wp) :: e0, r0, kk, c
        Contains 
            Procedure :: energy => morse12_energy
            Procedure :: set_parameters => morse12_set_parameters
    End Type morse12

    Type, Extends(potential), Public :: buckingham 
        Real(Kind=wp) :: A, rho, C 
        Contains
            Procedure :: energy => buckingham_energy
            Procedure :: set_parameters => buckingham_set_parameters
    End Type buckingham

    Type, Extends(potential), Public :: bhm
        Real(Kind=wp) :: a, b, sig, c, d
        Contains
            Procedure :: energy => bhm_energy
            Procedure :: set_parameters => bhm_set_parameters
    End Type bhm

    Type, Extends(potential), Public :: hbond
        Real(Kind=wp) :: a, b
        Contains 
            Procedure :: energy => hbond_energy
            Procedure :: set_parameters => hbond_set_parameters
    End Type hbond

    Type, Extends(potential), Public :: wca 
        Real(Kind=wp) :: eps, sig, d, cut 
        Contains 
            Procedure :: energy => wca_energy
            Procedure :: set_parameters => wca_set_parameters
    End Type wca

    Type, Extends(potential), Public :: dpd 
        Real(Kind=wp) :: a, rc 
        Contains 
            Procedure :: energy => dpd_energy
            Procedure :: set_parameters => dpd_set_parameters
    End Type dpd

    Type, Extends(potential), Public :: ndpd
        Real(Kind=wp) ::a, b, n, rc
        Contains 
            Procedure :: energy => ndpd_energy 
            Procedure :: set_parameters => ndpd_set_parameters
    End Type ndpd

    Type, Extends(potential), Public :: amoeba 
        Real(Kind=wp) :: eps, sig 
        Contains 
            Procedure :: energy => amoeba_energy
            Procedure :: set_parameters => amoeba_set_parameters
    End Type amoeba

    Type, Extends(potential), Public :: rydberg 
        Real(Kind=wp) :: a, b, c
        Contains
            Procedure :: energy => rydberg_energy
            Procedure :: set_parameters => rydberg_set_parameters
    End Type rydberg

    Type, Extends(potential), Public :: zbl
        Real(Kind=wp) :: k, ia 


        ! ZBL constants
        Real(wp), Private, Dimension(4) :: zbl_b = [0.18175_wp, 0.50986_wp, 0.28022_wp, 0.02817_wp], &
                                           zbl_c = [3.1998_wp, 0.94229_wp, 0.40290_wp, 0.20162_wp]
        Real(wp), Private               :: zbl_ab = 0.52917721067_wp
        
        Contains 
            Procedure :: energy => zbl_energy
            Procedure :: set_parameters => zbl_set_parameters
    End Type zbl 

    Type, Extends(potential), Public :: fm 
        Real(Kind=wp) :: rm, ic 
        Contains 
            Procedure :: energy => fm_energy
            Procedure :: set_parameters => fm_set_parameters
    End Type

    Type, Extends(potential), Public :: zbls 
        Type(zbl)    :: zbl
        Type(fm)     :: fm
        Type(morse)  :: morse 

        Contains 
            Procedure :: energy => zbls_energy
            Procedure :: set_parameters => zbls_set_parameters
    End Type zbls 

    Type, Extends(potential), Public :: zblb 
        Type(zbl)        :: zbl
        Type(fm)         :: fm 
        Type(buckingham) :: buckingham
        Contains 
            Procedure :: energy => zblb_energy
            Procedure :: set_parameters => zblb_set_parameters
    End Type zblb

    Type, Extends(potential), Public :: sanderson
        Real(Kind=wp) :: A, L, d
        Contains 
            Procedure :: energy => sanderson_energy
            Procedure :: set_parameters => sanderson_set_parameters
    End Type sanderson

    Type, Extends(potential), Public :: MDF 
        Real(Kind=wp) :: ri, rc 
        Contains 
            Procedure :: energy => MDF_energy
            Procedure :: set_parameters => MDF_set_parameters
    End Type MDF

    Type, Extends(potential), Public :: ljf 
        Real(Kind=wp) :: ea, sig2, rc2 
        Contains 
            Procedure :: energy => ljf_energy
            Procedure :: set_parameters => ljf_set_parameters
    End Type ljf

    Type, Extends(potential), Public :: mlj
        Type(LJ)  :: LJ 
        Type(MDF) :: MDF
        Contains 
            Procedure :: energy => mlj_energy
            Procedure :: set_parameters => mlj_set_parameters
    End Type mlj

    Type, Extends(potential), Public :: mbuck 
        Type(buckingham) :: buckingham 
        Type(MDF)        :: MDF 
        Contains 
            Procedure :: energy => mbuck_energy
            Procedure :: set_parameters => mbuck_set_parameters
    End Type mbuck

    Type, Extends(potential), Public :: mlj126
        Type(LJ126) :: LJ126 
        Type(MDF)   :: MDF 
        Contains 
            Procedure :: energy => mlj126_energy
            Procedure :: set_parameters => mlj126_set_parameters
    End Type mlj126

    Type, Extends(potential), Public :: sw 
        Real(Kind=wp) :: eps, A, B, sig, p, q, aa
        Contains 
            Procedure :: energy => sw_energy
            Procedure :: set_parameters => sw_set_parameters
    End Type sw

Contains

    Pure Type(potential_energy) Function holder_energy(h, r)
        Class(potential_holder), Intent(In   ) :: h
        Real(Kind=wp),           Intent(In   ) :: r
        holder_energy = h%p%energy(r)
    End Function

    Pure Type(potential_energy) Function lj_energy(p, r)
        ! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]
        Class(LJ),              Intent(In   ) :: p
        Real(Kind=wp),          Intent(In   ) :: r

        Real(wp) :: sor6

        sor6 = (p%sigma / r)**6
        lj_energy%energy = 4.0_wp * p%epsilon * sor6 * (sor6 - 1.0_wp)
        lj_energy%gamma = 24.0_wp * p%epsilon * sor6 * (2.0_wp * sor6 - 1.0_wp)
    End Function lj_energy

    Subroutine lj_set_parameters(p, param)
        Class(LJ),               Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for lj potential")
        End If
        p%epsilon = param(1)
        p%sigma = param(2)
    End Subroutine lj_set_parameters

    Pure Type(potential_energy) Function lj_coh_energy(p, r)
        ! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]
        Class(lj_coh),           Intent(In   ) :: p
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: sor6

        sor6 = (p%sigma / r)**6
        lj_coh_energy%energy = 4.0_wp * p%epsilon * sor6 * (sor6 - p%coh)
        lj_coh_energy%gamma = 24.0_wp * p%epsilon * sor6 * (2.0_wp * sor6 - p%coh)
    End Function lj_coh_energy 

    Subroutine lj_coh_set_parameters(p, param)
        Class(lj_coh),           Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 3) Then 
            Call error(0, "too few parameters for lj_coh potential")
        End If
        p%epsilon = param(1)
        p%sigma = param(2)
        p%coh = param(3)
    End Subroutine lj_coh_set_parameters

    Pure Type(potential_energy) Function LJ126_energy(p, r)
        ! 12-6 potential :: u=a/r^12-b/r^6
        Class(LJ126),            Intent(In   ) :: p
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: r_6

        r_6 = (1.0_wp / r)**6
        LJ126_energy%energy = (p%a * r_6 - p%b) * r_6
        LJ126_energy%gamma = 6.0_wp * r_6 * (2.0_wp * p%a * r_6 - p%b)
    End Function LJ126_energy

    Subroutine LJ126_set_parameters(p, param)
        Class(LJ126),            Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for LJ126 potential")
        End If
        p%a = param(1)
        p%b = param(2)
    End Subroutine LJ126_set_parameters

    Pure Type(potential_energy) Function n_m_energy(p, r)
        ! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(r0/r)^m]
        Class(n_m),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: a, b, r_n, r_m
        !e0, n, m, r0
        a = p%r0 / r
        b = 1.0_wp / (p%n - p%m)
        r_n = a**int(p%n)
        r_m = a**int(p%m)
    
        n_m_energy%energy = p%e0 * (p%m * r_n - p%n * r_m) * b
        n_m_energy%gamma = p%e0 * p%m * p%n * (r_n - r_m) * b

    End Function n_m_energy

    Subroutine n_m_set_parameters(p, param)
        Class(n_m),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 4) Then 
            Call error(0, "too few parameters for n_m potential")
        End If
        p%e0 = param(1)
        p%n = param(2)
        p%m = param(3)
        p%r0 = param(4)
    End Subroutine n_m_set_parameters

    Pure Type(potential_energy) Function nm_shift_energy(p, r)
        ! shifted and force corrected n-m potential (w.smith)
        Class(nm_shift),         Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: r_inv, r_inv_2
        Real(wp) :: alpha, beta, t, a, b, c, e1, c_inv
        Real(wp) :: n_int, m_int
        Real(wp) :: n, m

        If (r <= p%r_trunc) Then
            r_inv = r**(-1)
            r_inv_2 = r_inv**2
            n = p%n
            n_int = Nint(n)
            m = p%m
            m_int = Nint(m)
        
            t = n - m
        
            b = 1.0_wp / t
            c = p%r_trunc / p%r0
            c_inv = p%r0 / p%r_trunc
        
            beta = c * ((c**(m_int + 1) - 1.0_wp) / (c**(n_int + 1) - 1.0_wp))**b
            alpha = -t / (m * (beta**n_int) * (1.0_wp + (n * c_inv - n - 1.0_wp) * c_inv**n_int) &
                    - n * (beta**m_int) * (1.0_wp + (m * c_inv - m - 1.0_wp) * c_inv**m_int))
            e1 = p%e0 * alpha
        
            a = p%r0 * r_inv
        
            nm_shift_energy%energy = e1 * (m * (beta**n_int) * (a**n_int - (1.0_wp * c_inv)**n_int) &
                    - n * (beta**m_int) * (a**m_int - (1.0_wp * c_inv)**m_int) &
                    + n * m * ((r / p%r_trunc - 1.0_wp) * ((beta * c_inv)**n_int - (beta * c_inv)**m_int))) * b
            nm_shift_energy%gamma = e1 * m * n * ((beta**n_int) * a**n_int - (beta**m_int) * a**m_int &
                    - r / p%r_trunc * ((beta * c_inv)**n_int - (beta * c_inv)**m_int)) * b
        Else
            nm_shift_energy%energy = 0.0_wp
            nm_shift_energy%gamma = 0.0_wp
        End If
    End Function nm_shift_energy

    Subroutine nm_shift_set_parameters(p, param)
        Class(nm_shift),         Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 5) Then 
            Call error(0, "too few parameters for nm_shift potential")
        End If
        p%e0 = param(1)
        p%n = param(2)
        p%m = param(3)
        p%r0 = param(4)
        p%r_trunc = param(5)
    End Subroutine nm_shift_set_parameters

    Pure Type(potential_energy) Function morse_energy(p, r)
        ! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}
        Class(morse),            Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: t

        t = Exp(-p%k * (r - p%r0))

        morse_energy%energy = p%e0 * ((1.0_wp - t)**2 - 1.0_wp)
        morse_energy%gamma = -2.0_wp * r * p%e0 * p%k * (1.0_wp - t) * t

    End Function morse_energy

    Subroutine morse_set_parameters(p, param)
        Class(morse),            Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 3) Then 
            Call error(0, "too few parameters for morse potential")
        End If
        p%e0 = param(1)
        p%r0 = param(2)
        p%k = param(3)
    End Subroutine morse_set_parameters

    Pure Type(potential_energy) Function morse12_energy(p, r)
        ! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}+c/r^12
        Class(morse12),          Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r


        Real(wp) :: t1, t2

        t1 = Exp(-p%kk * (r - p%r0))
        t2 = p%c * r**(-12)

        morse12_energy%energy = p%e0 * t1 * (t1 - 2.0_wp) + t2
        morse12_energy%gamma = -2.0_wp * r * p%e0 * p%kk * (1.0_wp - t1) * t1 + 12.0_wp * t2

    End Function morse12_energy

    Subroutine morse12_set_parameters(p, param)
        Class(morse12),          Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 4) Then 
            Call error(0, "too few parameters for morse12 potential")
        End If
        p%e0 = param(1)
        p%r0 = param(2)
        p%kk = param(3)
        p%c = param(4)
    End Subroutine morse12_set_parameters

    Pure Type(potential_energy) Function buckingham_energy(p, r)
        ! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6
        Class(buckingham),       Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: b, t1, t2

        b = r / p%rho
        t1 = p%A * Exp(-b)
        t2 = -p%C / r**6
    
        buckingham_energy%energy = t1 + t2
        buckingham_energy%gamma = t1 * b + 6.0_wp * t2

    End Function buckingham_energy

    Subroutine buckingham_set_parameters(p, param)
        Class(buckingham),       Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 3) Then 
            Call error(0, "too few parameters for buckingham potential")
        End If
        p%A = param(1)
        p%rho = param(2)
        p%C = param(3)
    End Subroutine buckingham_set_parameters

    Pure Type(potential_energy) Function bhm_energy(p, r)
        ! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8
        Class(bhm),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: r_inv_2, t1, t2, t3

        r_inv_2 = r**(-2)
    
        t1 = p%a * Exp(p%b * (p%sig - r))
        t2 = -p%c * r_inv_2**3
        t3 = -p%d * r_inv_2**4
    
        bhm_energy%energy = t1 + t2 + t3
        bhm_energy%gamma = (t1 * r * p%b + 6.0_wp * t2 + 8.0_wp * t3)
    End Function bhm_energy

    Subroutine bhm_set_parameters(p, param)
        Class(bhm),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 5) Then 
            Call error(0, "too few parameters for bhm potential")
        End If
        p%a = param(1)
        p%b = param(2)
        p%sig = param(3)
        p%c = param(4)
        p%d = param(5)
    End Subroutine bhm_set_parameters

    Pure Type(potential_energy) Function hbond_energy(p, r)
        ! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10
        Class(hbond),            Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: fac12, fac10
        Real(wp) :: r_inv_2
    
        r_inv_2 = r**(-2)
    
        fac12 = p%a * r_inv_2**6
        fac10 = -p%b * r_inv_2**5
    
        hbond_energy%energy = fac12 + fac10
        hbond_energy%gamma = (12.0_wp * fac12 + 10.0_wp * fac10)

    End Function hbond_energy

    Subroutine hbond_set_parameters(p, param)
        Class(hbond),            Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for hbond potential")
        End If
        p%a = param(1)
        p%b = param(2)
    End Subroutine

    Pure Type(potential_energy) Function wca_energy(p, r)
        ! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
        ! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps
        Class(wca),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: sigma_r_6

        If (r < p%cut .or. Abs(r - p%d) < 1.0e-10_wp) Then
          sigma_r_6 = (p%sig / (r - p%d))**6
    
          wca_energy%energy = 4.0_wp * p%eps * sigma_r_6 * (sigma_r_6 - 1.0_wp) + p%eps
          wca_energy%gamma = 24.0_wp * p%eps * sigma_r_6 * (2.0_wp * sigma_r_6 - 1.0_wp) * r / (r - p%d)
        Else
          wca_energy%energy = 0.0_wp
          wca_energy%gamma = 0.0_wp
        end If
    End Function wca_energy

    Subroutine wca_set_parameters(p, param)
        Class(wca),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 4) Then 
            Call error(0, "too few parameters for wca potential")
        End If
        p%eps = param(1)
        p%sig = param(2)
        p%d = param(3)
        p%cut = param(4)
    End Subroutine wca_set_parameters

    Pure Type(potential_energy) Function dpd_energy(p, r)
        ! DPD potential - Groot-Warren (standard) :: u=(1/2).a.rc.(1-r/rc)^2
        Class(dpd),              Intent(In   ) :: p 
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: t1, t2

        If (r < p%rc) Then
    
          t2 = r / p%rc
          t1 = 0.5_wp * p%a * p%rc * (1.0_wp - t2)
    
          dpd_energy%energy = t1 * (1.0_wp - t2)
          dpd_energy%gamma = 2.0_wp * t1 * t2
        Else
    
          dpd_energy%energy = 0.0_wp
          dpd_energy%gamma = 0.0_wp
        End If
    End Function dpd_energy

    Subroutine dpd_set_parameters(p, param)
        Class(dpd),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for dpd potential")
        End If
        p%a = param(1)
        p%rc = param(2)
    End Subroutine dpd_set_parameters

    Pure Type(potential_energy) Function ndpd_energy(p, r)
        ! nDPD potential :: u = (1/(n+1)).a.b.rc.(1-r/rc)^(n+1)-(1/2).a.rc.(1-r/rc)^2
        Class(ndpd),             Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: t0, t1, t2

        If (r < p%rc) Then
          t2 = r / p%rc
          t1 = p%a * p%rc * (1.0_wp - t2)
          t0 = p%b * (1.0_wp - t2) ** (p%n - 1.0_wp)
    
          ndpd_energy%energy = t1 * (1.0_wp - t2) * (t0 / (p%n+1.0_wp) - 0.5_wp)
          ndpd_energy%gamma = t1 * t2 * (t0 - 1.0_wp)
        Else    
          ndpd_energy%energy = 0.0_wp
          ndpd_energy%gamma = 0.0_wp
        End If

    End Function ndpd_energy

    Subroutine ndpd_set_parameters(p, param)
        Class(ndpd),             Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 4) Then 
            Call error(0, "too few parameters for ndpd potential")
        End If
        p%a = param(1)
        p%b = param(2)
        p%n = param(3)
        p%rc = param(4)
    End Subroutine ndpd_set_parameters

    Pure Type(potential_energy) Function amoeba_energy(p, r)
        ! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]
        Class(amoeba),           Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: r_inv_2
        Real(wp) :: rho, t1, t2, t3, t
    
        r_inv_2 = r**(-2)
        rho = r / p%sig
    
        t1 = 1.0_wp / (0.07_wp + rho)
        t2 = 1.0_wp / (0.12_wp + rho**7)
        t3 = p%eps * (1.07_wp * t1)**7
    
        t = t3 * ((1.12_wp * t2) - 2.0_wp)
    
        amoeba_energy%energy = t
        amoeba_energy%gamma = 7.0_wp * (t1 * t + 1.12_wp * t3 * t2**2 * rho**6) * rho
    End Function amoeba_energy

    Subroutine amoeba_set_parameters(p, param)
        Class(amoeba),           Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for amoeba potential")
        End If
        p%eps = param(1)
        p%sig = param(2)
    End Subroutine amoeba_set_parameters

    Pure Type(potential_energy) Function rydberg_energy(p, r)
        ! Rydberg potential:: u=(a+b*r)Exp(-r/c)
        Class(rydberg),          Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: kk, t1

        kk = r / p%c
        t1 = Exp(-kk)
    
        rydberg_energy%energy = (p%a + p%b * r) * t1
        rydberg_energy%gamma = kk * t1 * (p%a - p%b * p%c + p%b * r)

    End Function rydberg_energy

    Subroutine rydberg_set_parameters(p, param)
        Class(rydberg),          Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 3) Then 
            Call error(0, "too few parameters for rydberg potential")
        End If
        p%a = param(1)
        p%b = param(2)
        p%c = param(3)
    End Subroutine rydberg_set_parameters

    Pure Type(potential_energy) Function zbl_energy(p, r)
        ! ZBL potential:: u=Z1Z2/(40r)_{i=1}^4b_ie^{-c_i*r/a}
        Class(zbl),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Integer  :: i
        Real(wp) :: kk, a
        Real(wp) :: ir, t1, x
    
        ! this is in fact inverse a
        a = (p%k**0.23_wp + p%ia**0.23_wp) / (p%zbl_ab * 0.88534_wp)
        kk = p%k * p%ia * r4pie0
    
        zbl_energy%energy = 0.0_wp
        zbl_energy%gamma = 0.0_wp
        x = r * a
        ir = 1.0_wp / r
        Do i = 1, 4
          t1 = p%zbl_b(i) * Exp(-x * p%zbl_c(i))
          zbl_energy%energy = zbl_energy%energy + t1
          zbl_energy%gamma = zbl_energy%gamma - p%zbl_c(i) * t1
        End Do
        zbl_energy%energy = kk * zbl_energy%energy * ir
        ! -rU/r
        zbl_energy%gamma = zbl_energy%energy - a * kk * zbl_energy%gamma
    End Function zbl_energy

    Subroutine zbl_set_parameters(p, param)
        Class(zbl),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for zbl potential")
        End If
        p%k = param(1)
        p%ia = param(2)
    End Subroutine zbl_set_parameters

    Pure Type(potential_energy) Function fm_energy(p, r)
        Class(fm),               Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: t, c

        c = 1.0_wp / p%ic

        If (r < p%rm) Then
            t = Exp(-(p%rm - r) * c) * 0.5_wp
            fm_energy%energy = 1.0_wp - t
            ! -rf/r
            fm_energy%gamma = r * c * t
        Else
            t = Exp(-(r - p%rm) * c) * 0.5_wp
            fm_energy%energy = t
            ! -rf/r
            fm_energy%gamma = r * c * t
        End If

    End Function fm_energy

    Subroutine fm_set_parameters(p, param)
        Class(fm),               Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for fm potential")
        End If
        p%rm = param(1)
        p%ic = param(2)
    End Subroutine fm_set_parameters

    Pure Type(potential_energy) Function zbls_energy(p, r)
        Class(zbls),             Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r
        Type(potential_energy)                 :: zbl_e, fm_e, morse_e

        Real(wp) :: df, dm, dz, f, m, z

        zbl_e = p%zbl%energy(r)
        fm_e = p%fm%energy(r)
        morse_e = p%morse%energy(r)

        z = zbl_e%energy
        dz = zbl_e%gamma
        f = fm_e%energy
        df = fm_e%gamma
        m = morse_e%energy
        dm = morse_e%gamma

        zbls_energy%energy = f * z + (1.0_wp - f) * m
        zbls_energy%gamma = f * dz + df * z + (1.0_wp - f) * dm - df * m

    End Function zbls_energy

    Subroutine zbls_set_parameters(p, param)
        Class(zbls),             Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 7) Then 
            Call error(0, "too few parameters for zbls potential")
        End If
        p%zbl%k = param(1)
        p%zbl%ia = param(2)
        p%fm%rm = param(3)
        p%fm%ic = param(4)
        p%morse%e0 = param(5)
        p%morse%r0 = param(6)
        p%morse%k = param(7)
    End Subroutine zbls_set_parameters

    Pure Type(potential_energy) Function zblb_energy(p, r)
        Class(zblb),             Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r
        Type(potential_energy)                 :: zbl_e, fm_e, buckingham_e

        Real(wp) :: b, db, df, dz, f, z

        zbl_e = p%zbl%energy(r)
        z = zbl_e%energy
        dz = zbl_e%gamma
        
        fm_e = p%fm%energy(r)
        f = fm_e%energy
        df = fm_e%gamma

        buckingham_e = p%buckingham%energy(r)
        b = buckingham_e%energy
        db = buckingham_e%gamma

        zblb_energy%energy = f * z + (1.0_wp - f) * b
        zblb_energy%gamma = f * dz + df * z + (1.0_wp - f) * db - df * b
    End Function zblb_energy

    Subroutine zblb_set_parameters(p, param)
        Class(zblb),             Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 7) Then 
            Call error(0, "too few parameters for zblb potential")
        End If
        p%zbl%k = param(1)
        p%zbl%ia = param(2)
        p%fm%rm = param(3)
        p%fm%ic = param(4)
        p%buckingham%A = param(5)
        p%buckingham%rho = param(6)
        p%buckingham%C = param(7)
    End Subroutine zblb_set_parameters

    Pure Type(potential_energy) Function sanderson_energy(p, r)
        Class(sanderson),        Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: b, t, L, d

        L = p%L
        d = p%d
        b = ((r - L)/d)**2.0_wp
        t = p%A * Exp(-b)

        sanderson_energy%energy = -t
        sanderson_energy%gamma = -2.0_wp*(r-L)*r*t/(d**2.0_wp)

    End Function sanderson_energy

    Subroutine sanderson_set_parameters(p, param)
        Class(sanderson),        Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 3) Then 
            Call error(0, "too few parameters for sanderson potential")
        End If
        p%A = param(1)
        p%L = param(2)
        p%d = param(3)
    End Subroutine sanderson_set_parameters

    Pure Type(potential_energy) Function MDF_energy(p, r)
        Class(MDF),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: rci

        If (r < p%ri) Then
          MDF_energy%energy = 1.0_wp
          MDF_energy%gamma = 0.0_wp
        Else If (r > p%rc) Then
          MDF_energy%energy = 0.0_wp
          MDF_energy%gamma = 0.0_wp
        Else
          rci = (p%rc - p%ri)**5
          MDF_energy%energy = (p%rc - r)**3 * (10.0_wp * p%ri**2 - 5.0_wp * p%rc * p%ri - &
               15.0_wp * r * p%ri + p%rc**2 + 3.0_wp * r * p%rc + 6 * r**2) / rci
          MDF_energy%gamma = 30.0_wp * r * (r - p%rc)**2 * (r - p%ri)**2 / rci
        End If
    End Function MDF_energy

    Subroutine MDF_set_parameters(p, param)
        Class(MDF),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 2) Then 
            Call error(0, "too few parameters for MDF potential")
        End If
        p%ri = param(1)
        p%rc = param(2)
    End Subroutine MDF_set_parameters

    Pure Type(potential_energy) Function ljf_energy(p, r)
        Class(ljf),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r

        Real(wp) :: r2
        Real(wp) :: ir, rct, st, x
    
        r2 = r*r
    
        If (r2 > p%rc2) Then
          ljf_energy%energy = 0.0_wp
          ljf_energy%gamma = 0.0_wp
        Else
          ir = 1.0_wp/r2
          st = p%sig2*ir
          rct = p%rc2*ir
          x = p%ea * (rct - 1.0_wp)**2
          ljf_energy%energy = x * (st - 1.0_wp)
          ljf_energy%gamma = 4.0_wp * p%ea * rct * &
               (rct - 1.0_wp)*(st - 1.0_wp) + 2.0_wp * x * st
        End If
    End Function ljf_energy

    Subroutine ljf_set_parameters(p, param)
        Class(ljf),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 3) Then 
            Call error(0, "too few parameters for ljf potential")
        End If
        p%ea = param(1)
        p%sig2 = param(2)
        p%rc2 = param(3)
    End Subroutine ljf_set_parameters

    Pure Type(potential_energy) Function mlj_energy(p, r)
        Class(mlj),              Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r
        Type(potential_energy)                 :: e_lj, e_mdf

        Real(wp) :: el, em, vl, vm

        e_lj = p%LJ%energy(r)
        el = e_lj%energy
        vl = e_lj%gamma

        e_mdf = p%MDF%energy(r)
        em = e_mdf%energy
        vm = e_mdf%gamma

        mlj_energy%energy = el * em
        mlj_energy%gamma = vl * em + vm * el
    End Function mlj_energy

    Subroutine mlj_set_parameters(p, param)
        Class(mlj),              Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 4) Then 
            Call error(0, "too few parameters for mlj potential")
        End If
        p%LJ%epsilon = param(1)
        p%LJ%sigma = param(2)
        p%MDF%ri = param(3)
        p%MDF%rc = param(4)
    End Subroutine mlj_set_parameters

    Pure Type(potential_energy) Function mbuck_energy(p, r)
        Class(mbuck),            Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r
        Type(potential_energy)                 :: e_buckingham, e_mdf 

        Real(wp) :: eb, em, vb, vm

        e_buckingham = p%buckingham%energy(r)
        eb = e_buckingham%energy
        vb = e_buckingham%gamma
        
        e_mdf = p%MDF%energy(r)
        em = e_mdf%energy
        vm = e_mdf%gamma

        mbuck_energy%energy = eb * em
        mbuck_energy%gamma = vb * em + vm * eb
        
    End Function mbuck_energy

    Subroutine mbuck_set_parameters(p, param)
        Class(mbuck),            Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 5) Then 
            Call error(0, "too few parameters for mbuck potential")
        End If
        p%buckingham%A = param(1)
        p%buckingham%rho = param(2)
        p%buckingham%C = param(3)
        p%MDF%ri = param(4)
        p%MDF%rc = param(5)
    End Subroutine mbuck_set_parameters

    Pure Type(potential_energy) Function mlj126_energy(p, r)
        Class(mlj126),           Intent(In   ) :: p      
        Real(Kind=wp),           Intent(In   ) :: r
        Type(potential_energy)                 :: e_lj126, e_mdf

        Real(Kind=wp) :: el, vl, em, vm

        e_lj126 = p%LJ126%energy(r)
        el = e_lj126%energy
        vl = e_lj126%gamma
 
        e_mdf = p%MDF%energy(r)
        em = e_mdf%energy
        vm = e_mdf%gamma

        mlj126_energy%energy = el * em
        mlj126_energy%gamma = vl * em + vm * el
    End Function mlj126_energy

    Subroutine mlj126_set_parameters(p, param)
        Class(mlj126),           Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 4) Then 
            Call error(0, "too few parameters for mlj126 potential")
        End If
        p%LJ126%a = param(1)
        p%LJ126%b = param(2)
        p%MDF%ri = param(3)
        p%MDF%rc = param(4)
        
    End Subroutine mlj126_set_parameters

    Pure Type(potential_energy) Function sw_energy(p, r)
        ! Stillinger Webber 2 body part potential $u = A*\varepsilon*[B(\sigma/r)^p - (sigma/r)^q]*exp(\sigma/(r-aa*\sigma))}$
        Class(sw),     Intent(In   ) :: p 
        Real(Kind=wp), Intent(In   ) :: r

        Real(Kind=wp) :: t, e, p_r, exp_e, c

        e = p%sig/(r - p%aa*p%sig)

        If (r < p%aa*p%sig) Then
            p_r = p%sig / r
            c = p%B*(p_r)**p%p
            exp_e = Exp(e)

            t = p%A*p%eps*(c - (p_r)**p%q)*exp_e

            sw_energy%energy = t
            sw_energy%gamma = p%A * p%eps * (p%p*c - &
                p%q*(p_r)**p%q)*exp_e + t*r*e/(r - p%aa*p%sig)
        Else
            sw_energy%energy = 0.0_wp
            sw_energy%gamma = 0.0_wp
        End If
        
    End Function 

    Subroutine sw_set_parameters(p, param)
        Class(sw),               Intent(InOut) :: p
        Real(Kind=wp),           Intent(In   ) :: param(:)

        If (Size(param) < 7) Then 
            Call error(0, "too few parameters for sw potential")
        End If
        p%eps = param(1)
        p%A = param(2)
        p%B = param(3)
        p%sig = param(4)
        p%p = param(5)
        p%q = param(6)
        p%aa = param(7)
    End Subroutine sw_set_parameters

End Module two_body_potentials