Module vdw

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global VdW interaction variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov november 2014
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  ! contrib   - a.m.elena march 2019 ! merge potentials.F90 into this module
  ! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
  ! contrib   - j.s.wilkins october 2020 Major cleanup
  ! contrib   - a.m.elena march 2023 add sw potential
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,                    Only: comms_type,&
                                      gbcast,&
                                      gsum
  Use configuration,            Only: configuration_type
  Use constants,                Only: delr_max,&
                                      engunit,&
                                      prsunt,&
                                      r4pie0,&
                                      twopi,&
                                      zero_plus
  Use errors_warnings,          Only: error,&
                                      error_alloc,&
                                      error_dealloc,&
                                      info,&
                                      warning
  Use kinds,                    Only: wi,&
                                      wp,&
                                      STR_LEN
  Use filename,                 Only: FILE_TABVDW, &
                                      file_type
  Use neighbours,               Only: neighbours_type
  Use numerics,                 Only: nequal
  Use parse,                    Only: get_line,&
                                      get_word,&
                                      word_2_real
  Use two_body_potentials,      Only: potential, &
                                      potential_energy, &
                                      potential_holder, &
                                      LJ, lj_coh, LJ126, n_m, &
                                      nm_shift, morse, morse12, &
                                      buckingham, bhm, hbond, &
                                      wca, dpd, ndpd, amoeba, &
                                      rydberg, zbl, zblb, fm, zbls, &
                                      sanderson, MDF, ljf, mlj, &
                                      mbuck, mlj126, sw
  Use site,                      Only: site_type
  Use statistics,                Only: stats_type, calculate_stress
  Use units,                     Only: to_out_units
  Implicit None

  Private

  ! VdW potential parameters
  !> Number of available potential types
  Integer(Kind=wi), Parameter, Public :: NUM_VDW_POTS = 24

  !> No VdW potential
  Integer(Kind=wi), Parameter, Public :: VDW_NULL = -1
  !> Tabulated potential
  Integer(Kind=wi), Parameter, Public :: VDW_TAB = 0
  !> 12-6 potential (Lennard-Jones A B form): $u=A/r^12-B/r^6$
  Integer(Kind=wi), Parameter, Public :: VDW_12_6 = 1
  !> Lennard-Jones potential: $u=4*\epsilon*[(\sigma/r)^12-(\sigma/r)^6]$
  Integer(Kind=wi), Parameter, Public :: VDW_LENNARD_JONES = 2
  !> n-m potential: $u=[e_0/(n-m)]*[m*(r_0/r)^n-n*(d/r)^c]$
  Integer(Kind=wi), Parameter, Public :: VDW_N_M = 3
  !> Buckingham exp-6 potential: $u=a*\exp(-r/\rho)-c/r^6$
  Integer(Kind=wi), Parameter, Public :: VDW_BUCKINGHAM = 4
  !> Born-Huggins-Meyer exp-6-8 potential: $u=a*\exp(b*(\sigma-r)) - c/r^6 - d/r^8$
  Integer(Kind=wi), Parameter, Public :: VDW_BORN_HUGGINS_MEYER = 5
  !> Hydrogen-bond 12-10 potential: $u=a/r^12-b/r^10$
  Integer(Kind=wi), Parameter, Public :: VDW_HYDROGEN_BOND = 6
  !> n-m potential shift forced
  Integer(Kind=wi), Parameter, Public :: VDW_N_M_SHIFT = 7
  !> Morse potential: $u=e_0*([1-\exp(-k(r-r_0))]^2-1)$
  Integer(Kind=wi), Parameter, Public :: VDW_MORSE = 8
  !> Shifted Weeks-Chandler-Anderson potential
  Integer(Kind=wi), Parameter, Public :: VDW_WCA = 9
  !> DPD potential: $u=(A*r_c)/2*(1-r/r_c)^2$
  Integer(Kind=wi), Parameter, Public :: VDW_DPD = 10
  !> AMOEBA 14-7: $u=eps * [1.07/((\sigma/r)+0.07)]^7 * [(1.12/((\sigma/r)^7+0.12))-2]$
  Integer(Kind=wi), Parameter, Public :: VDW_AMOEBA = 11
  !> Lennard-Jones cohesive potential: $u=4*eps*[(\sigma/r)^12-c*(\sigma/r)^6]$
  Integer(Kind=wi), Parameter, Public :: VDW_LENNARD_JONES_COHESIVE = 12
  !> Morse potential with r^12 repuconfig%lsion (mstw): $u=e_0*{[1-\exp(-k(r-r_0))]^2-1}+c/r^12$
  Integer(Kind=wi), Parameter, Public :: VDW_MORSE_12 = 13
  !> Rydberg potential: $u=(a+b*r) \exp(-r/c)$
  Integer(Kind=wi), Parameter, Public :: VDW_RYDBERG = 14
  !> ZBL potential: $u=Z_1 Z_2/(4 \pi \epsilon_0 r) \sum_{i=1}^4 b_i e^{-c_i*r/a}$
  Integer(Kind=wi), Parameter, Public :: VDW_ZBL = 15
  !> ZBL switched with Morse: $u=f(r)\mathrm{zbl}(r)+(1-f(r))*\mathrm{morse}(r)$
  Integer(Kind=wi), Parameter, Public :: VDW_ZBL_SWITCH_MORSE = 16
  !> ZBL switched with Buckingham: $u=f(r)\mathrm{zbl}(r)+(1-f(r))*\mathrm{buckingham}(r)$
  Integer(Kind=wi), Parameter, Public :: VDW_ZBL_SWITCH_BUCKINGHAM = 17
  Integer(Kind=wi), Parameter, Public :: VDW_LJ_MDF = 18
  Integer(Kind=wi), Parameter, Public :: VDW_BUCKINGHAM_MDF = 19
  Integer(Kind=wi), Parameter, Public :: VDW_126_MDF = 20
  Integer(Kind=wi), Parameter, Public :: VDW_LJF = 21
  !> Sanderson potential $u = -A*exp{-[(r-L)/d ]**2}
  Integer(Kind=wi), Parameter, Public :: VDW_SANDERSON = 22
  !> nDPD potential: $u=(A*b*r_c)(n+1)*(1-r/r_c)^(n+1)-(A*r_c)/2*(1-r/r_c)^2$
  Integer(Kind=wi), Parameter, Public :: VDW_NDPD = 23
  ! Stillinger Webber 2 body part potential $u = A*\varepsilon*[B(\sigma/r)^p - (sigma/r)^q]*exp(\sigma/(r-aa*\sigma))}$
  Integer(Kind=wi), Parameter, Public :: VDW_SW = 24

  ! Mixing rule parameters
  !> Null
  Integer(Kind=wi), Parameter, Public :: MIX_NULL = 0
  !> Lorentz-Berthelot: $e_{ij}=(e_i*e_j)^{1/2} \quad s_{ij}=(s_i+s_j)/2$
  Integer(Kind=wi), Parameter, Public :: MIX_LORENTZ_BERTHELOT = 1
  !> Fender-Halsey: $e_{ij}=(2*e_i*e_j)/(e_i+e_j) \quad s_{ij}=(s_i+s_j)/2$
  Integer(Kind=wi), Parameter, Public :: MIX_FENDER_HALSEY = 2
  !> Hogervorst Good-Hope: $e_{ij}=(e_i*e_j)^{1/2} \quad s_{ij}=(s_i*s_j)^{1/2}$
  Integer(Kind=wi), Parameter, Public :: MIX_HOGERVORST = 3
  !> Halgren HHG: $e_{ij}=(4*e_i*e_j)/(e_i^{1/2}+e_j^{1/2})^2 \quad s_{ij}=(s_i^3+s_j^3)/(s_i^2+s_j^2)$
  Integer(Kind=wi), Parameter, Public :: MIX_HALGREN = 4
  !> WaldmanHagler: $e_{ij}=2*(e_i*e_j)^{1/2}*(s_i*s_j)^3/(s_i^6+s_j^6) \quad s_{ij}=[(s_i^6+s_j^6)/2]^{1/6}$
  Integer(Kind=wi), Parameter, Public :: MIX_WALDMAN_HAGLER = 5
  !> Tang-Toennies: $e_{ij}=[(e_i*s_i^6)*(e_j*s_j^6)] / ([(e_i*s_i^12)^{1/13}+(e_j*s_j^12)^{1/13}]/2)^13$
  !>                $s_{ij}=(1/3) \sum_{L=0}^2 [(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^{1/(6-2L)}$
  Integer(Kind=wi), Parameter, Public :: MIX_TANG_TOENNIES = 6
  !> Functional: $e_{ij}=3*(e_i*e_j)^{1/2} * (s_i*s_j)^3 / \sum_{L=0}^2 [(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(6/(6-2L))$
  !>             $s_ij=(1/3) \sum_{L=0}^2 [(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))$
  Integer(Kind=wi), Parameter, Public :: MIX_FUNCTIONAL = 7
  
  !> Type containing Van der Waals data
  Type, Public :: vdw_type
    Private

    !> No Van der Waals switch
    Logical, Public                       :: no_vdw = .false.
    !> Flag for any tabulated potential
    Logical, Public                       :: l_tab = .false.
    !> Direct calculation flag
    Logical, Public                       :: l_direct = .false.
    !> Force shifting flag
    Logical, Public                       :: l_force_shift = .false.
    !> Number of two body interactoins
    Integer(Kind=wi), Public              :: n_vdw = 0
    !> Mixing type
    Integer(Kind=wi), Public              :: mixing = MIX_NULL
    ! list of potentials
    Type(potential_holder), Allocatable, Public :: potentials(:)
    ! list of vdw interactions 
    Integer(Kind=wi), Allocatable, Public :: list(:)
    ! potential type
    Integer(Kind=wi), Allocatable, Public :: ltp(:)
    !> VdW label for pair
    Character(Len=8), Allocatable, Public :: labpair(:, :)
    !> VdW parameters
    Real(Kind=wp), Allocatable, Public    :: param(:, :)
    !> VdW cut off
    Real(Kind=wp), Public                 :: cutoff = 0.0_wp
    !> Energy long range correction
    Real(Kind=wp), Public                 :: elrc = 0.0_wp
    !> Virial long range correction
    Real(Kind=wp), Public                 :: vlrc = 0.0_wp
    ! Possible tabulated calculation arrays
    !> Tabulated potential
    Real(Kind=wp), Allocatable, Public    :: tab_potential(:, :)
    !> Tabulated force
    Real(Kind=wp), Allocatable, Public    :: tab_force(:, :)
    !> Maximum number of grid points
    Integer(Kind=wi), Public              :: max_grid
    ! Possible force-shifting arrays
    Real(Kind=wp), Allocatable, Public    :: afs(:)
    Real(Kind=wp), Allocatable, Public    :: bfs(:)
    !> Maximum number of VdW interations
    Integer(Kind=wi), Public              :: max_vdw
    !> Maximum number of VdW parameters
    Integer(Kind=wi), Public              :: max_param
    !> furst time job
    Logical, Public                       :: newjob = .true.
    Real(Kind=wp), Public                 :: dlrpot, rdr

  Contains
    Private

    Procedure, Public :: init => allocate_vdw_arrays
    Procedure, Public :: init_table => allocate_vdw_table_arrays
    Procedure, Public :: init_direct => allocate_vdw_direct_fs_arrays
    Procedure, Public :: print => dump_vdws
    Procedure, Public :: update_potential_cutoffs => potential_cutoffs
    Final             :: cleanup
  End Type vdw_type

  Public :: vdw_forces_tab, vdw_forces_direct, vdw_generate, vdw_table_read, vdw_lrc, vdw_direct_fs_generate

Contains

  Subroutine potential_cutoffs(vdws)
    Class(vdw_type),  Intent(InOut) :: vdws
    Integer                         :: ivdw, ityp
    Real(Kind=wp)                   :: param(1:5)
 
    Do ivdw = 1, vdws%max_vdw

      ityp = vdws%ltp(ivdw)

      Select Case (ityp)

      Case(VDW_LJ_MDF)

        param(1:3) = vdws%param(1:3, ivdw)
        param(4) = vdws%cutoff
        Call vdws%potentials(ivdw)%p%set_parameters(param(1:4))

      Case(VDW_BUCKINGHAM_MDF)

        param(1:4) = vdws%param(1:4, ivdw)
        param(5) = vdws%cutoff
        Call vdws%potentials(ivdw)%p%set_parameters(param(1:5))

      Case(VDW_126_MDF)

        param(1:3) = vdws%param(1:3, ivdw)
        param(4) = vdws%cutoff
        Call vdws%potentials(ivdw)%p%set_parameters(param(1:4))
      
      End Select

    End Do

  End Subroutine potential_cutoffs

  Subroutine allocate_vdw_arrays(T)
    Class(vdw_type) :: T

    Integer, Dimension(1:5) :: fail

    fail = 0

    Allocate (T%potentials(1:T%max_vdw), Stat=fail(1))
    Allocate (T%list(1:T%max_vdw), Stat=fail(2))
    Allocate (T%ltp(1:T%max_vdw), Stat=fail(3))
    Allocate (T%param(1:T%max_param, 1:T%max_vdw), Stat=fail(4))
    Allocate (T%labpair(1:2, 1:T%max_vdw), Stat=fail(5))

    If (Any(fail > 0)) Call error_alloc('VdW arrays', 'allocate_vdw_arrays')

    T%list = 0
    T%ltp = 0

    T%param = 0.0_wp

  End Subroutine allocate_vdw_arrays

  Subroutine allocate_vdw_table_arrays(T)
    Class(vdw_type) :: T

    Integer, Dimension(1:2) :: fail

    fail = 0

    Allocate (T%tab_potential(0:T%max_grid, 1:T%max_vdw), Stat=fail(1))
    Allocate (T%tab_force(0:T%max_grid, 1:T%max_vdw), Stat=fail(2))

    If (Any(fail > 0)) Call error_alloc('Table arrays', 'allocate_vdw_table_arrays')

    T%tab_potential = 0.0_wp
    T%tab_force = 0.0_wp
  End Subroutine allocate_vdw_table_arrays

  Subroutine allocate_vdw_direct_fs_arrays(T)
    Class(vdw_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%afs(1:T%max_vdw), T%bfs(1:T%max_vdw), Stat=fail)

    If (fail > 0) Call error_alloc('Force-shift arrays', 'allocate_vdw_direct_fs_arrays')

    T%afs = 0.0_wp
    T%bfs = 0.0_wp
  End Subroutine allocate_vdw_direct_fs_arrays

  Pure Real(wp) Function intRadZBL(kk, a, rw, prec)
    Real(wp), Intent(In   ) :: kk, a, rw, prec
    Type(zbl)               :: pot
    Type(potential_energy)  :: e 
    Integer                 :: i, j, n
    Real(wp)                :: f0, f1, f2, h, ie, is, s, sold, x

    n = 10000
    is = rw
    ie = 2 * is
    h = (ie - is) / Real(n, wp)
    s = 0.0_wp
    sold = Huge(1.0_wp)
    j = 1

    pot%k = kk
    pot%ia = a

    Do While (Abs(s - sold) * h / 3.0_wp > prec)
      sold = s

      Do i = (j - 1) * n, j * n, 2

        x = is + i * h

        e = pot%energy(x)
        f0 = e%energy
        
        e = pot%energy(x+h)
        f1 = e%energy
        
        e = pot%energy(x+2.0_wp*h)
        f2 = e%energy

        s = s + x * x * f0 + &
            4.0_wp * (x + h)**2 * f1 + &
            (x + 2.0_wp * h)**2 * f2
      End Do

      j = j + 1
    End Do
    intRadZBL = s * h / 3.0_wp
  End Function intRadZBL

  Real(wp) Function intdRadZBL(kk, a, rw, prec)
    Real(wp), Intent(In   ) :: kk, a, rw, prec
    Type(zbl)               :: pot
    Type(potential_energy)  :: e 
    Integer                 :: i, j, n
    Real(wp)                :: df0, df1, df2, h, ie, is, s, sold, x

    n = 10000
    is = rw
    ie = 2 * is
    h = (ie - is) / Real(n, wp)
    s = 0.0_wp
    sold = Huge(1.0_wp)
    j = 1
    
    pot%k = kk 
    pot%ia = a

    Do While (Abs(s - sold) * h / 3.0_wp > prec)
      sold = s

      Do i = (j - 1) * n, j * n, 2

        x = is + i * h

        e = pot%energy(x)
        df0 = e%gamma
        
        e = pot%energy(x+h)
        df1 = e%gamma
        
        e = pot%energy(x+2.0_wp*h)
        df2 = e%gamma

        s = s + x * x * df0 + &
            4.0_wp * (x + h)**2 * df1 + &
            (x + 2.0_wp * h)**2 * df2
      End Do

      j = j + 1
    End Do
    intdRadZBL = s * h / 3.0_wp
  End Function intdRadZBL

  Pure Real(wp) Function intRadMDF(pot, a, b, c, ri, rw, prec)
    Character(Len=*), Intent(In   ) :: pot
    Real(wp),         Intent(In   ) :: a, b, c, ri, rw, prec
    Class(potential), Allocatable   :: p
    Type(potential_energy)          :: e
    Type(mlj126)                    :: p_mlj126 
    Type(mbuck)                     :: p_mbuck
    Type(mlj)                       :: p_mlj 

    Integer  :: i, j, n
    Real(wp) :: f0, f1, f2, h, ie, is, s, sold, x

    If (pot == 'm126') Then 
      Allocate(mlj126::p)
      p_mlj126%LJ126%a = a 
      p_mlj126%LJ126%b = b
      p_mlj126%MDF%ri = ri 
      p_mlj126%MDF%rc = rw
      p = p_mlj126
    Else If (pot == 'mbuc') Then 
      Allocate(mbuck::p)
      p_mbuck%buckingham%A = a 
      p_mbuck%buckingham%rho = b 
      p_mbuck%buckingham%C = c
      p_mbuck%MDF%ri = ri 
      p_mbuck%MDF%rc = rw 
      p = p_mbuck
    Else If (pot == 'mlj') Then 
      Allocate(mlj::p)
      p_mlj%LJ%sigma = a 
      p_mlj%LJ%epsilon = b
      p_mlj%MDF%ri = ri
      p_mlj%MDF%rc = rw
      p = p_mlj
    Else  
      intRadMDF = 0.0_wp 
      Return 
    End If 

    n = 10000
    is = rw
    ie = 2 * is
    h = (ie - is) / Real(n, wp)
    s = 0.0_wp
    sold = Huge(1.0_wp)
    j = 1

    Do While (Abs(s - sold) * h / 3.0_wp > prec)
      sold = s

      Do i = (j - 1) * n, j * n, 2

        x = is + i * h

        e = p%energy(x)
        f0 = e%energy

        e = p%energy(x+h)
        f1 = e%energy

        e = p%energy(x+2.0_wp*h)
        f2 = e%energy

        s = s + x * x * f0 + &
            4.0_wp * (x + h)**2 * f1 + &
            (x + 2.0_wp * h)**2 * f2
      End Do

      j = j + 1
    End Do
    intRadMDF = s * h / 3.0_wp
  End Function intRadMDF

  Pure Real(wp) Function intdRadMDF(pot, a, b, c, ri, rw, prec)
    Character(Len=*), Intent(In   ) :: pot
    Real(wp),         Intent(In   ) :: a, b, c, ri, rw, prec
    Class(potential), Allocatable   :: p
    Type(potential_energy)          :: e
    Type(mlj126)                    :: p_mlj126 
    Type(mbuck)                     :: p_mbuck
    Type(mlj)                       :: p_mlj 

    Integer  :: i, j, n
    Real(wp) :: df0, df1, df2, h, ie, is, s, sold, x

    If (pot == 'm126') Then 
      Allocate(mlj126::p)
      p_mlj126%LJ126%a = a 
      p_mlj126%LJ126%b = b
      p_mlj126%MDF%ri = ri 
      p_mlj126%MDF%rc = rw
      p = p_mlj126
    Else If (pot == 'mbuc') Then 
      Allocate(mbuck::p)
      p_mbuck%buckingham%A = a 
      p_mbuck%buckingham%rho = b 
      p_mbuck%buckingham%C = c
      p_mbuck%MDF%ri = ri 
      p_mbuck%MDF%rc = rw 
      p = p_mbuck
    Else If (pot == 'mlj') Then 
      Allocate(mlj::p)
      p_mlj%LJ%sigma = a 
      p_mlj%LJ%epsilon = b
      p_mlj%MDF%ri = ri
      p_mlj%MDF%rc = rw
      p = p_mlj
    Else  
      intdRadMDF = 0.0_wp 
      Return 
    End If 

    n = 10000
    is = rw
    ie = 2 * is
    h = (ie - is) / Real(n, wp)
    s = 0.0_wp
    sold = Huge(1.0_wp)
    j = 1

    Do While (Abs(s - sold) * h / 3.0_wp > prec)
      sold = s

      Do i = (j - 1) * n, j * n, 2

        x = is + i * h
       
        e = p%energy(x)
        df0 = e%gamma

        e = p%energy(x+h)
        df1 = e%gamma

        e = p%energy(x+2.0_wp*h)
        df2 = e%gamma

        s = s + x * x * df0 + &
            4.0_wp * (x + h)**2 * df1 + &
            (x + 2.0_wp * h)**2 * df2
      End Do

      j = j + 1
    End Do
    intdRadMDF = s * h / 3.0_wp
  End Function intdRadMDF

  Function mm3(x, A, B, eps)

    Real(Kind=wp), Intent(In   ) :: x, A, B, eps
    Real(Kind=wp)                :: mm3

    Real(Kind=wp) :: ia, ib

    ia = 1.0_wp / (A + x)
    ib = 1.0_wp / (B + x**7)

    mm3 = eps * ((1.0_wp + A) * ia)**7 * ((1.0_wp + B) * ib - 2.0_wp)

  End Function mm3

  Function intRadMM3(r0, A, B, eps, rw, prec)

    Real(Kind=wp), Intent(In   ) :: r0, A, B, eps, rw, prec
    Real(Kind=wp)                :: intRadMM3

    Integer       :: i, j, n
    Real(Kind=wp) :: h, ie, is, s, sold, x

    n = 10000
    is = rw / r0
    ie = 2 * is
    h = (ie - is) / Real(n, wp)
    s = 0.0_wp
    sold = Huge(1.0_wp)
    j = 1

    Do While (Abs(s - sold) * h / 3.0_wp > prec)
      sold = s

      Do i = (j - 1) * n, j * n, 2
        x = is + i * h
        s = s + x * x * mm3(x, A, B, eps) + &
            4.0_wp * (x + h)**2 * mm3(x + h, A, B, eps) + &
            (x + 2.0_wp * h)**2 * mm3(x + 2 * h, A, B, eps)
      End Do

      j = j + 1
    End Do

    intRadMM3 = s * h * r0**3 / 3.0_wp

  End Function intRadMM3

  Function dmm3(x, r0, A, B, eps)

    Real(kind=wp), Intent(In   ) :: x, r0, A, B, eps
    Real(Kind=wp)                :: dmm3

    Real(Kind=wp) :: ia, ib

    ia = 1.0_wp / (A + x)
    ib = 1.0_wp / (B + x**7)

    dmm3 = -7.0_wp * (mm3(x, A, B, eps) * ia - ((1.0_wp + A) * ia)**7 * ((1.0_wp + B) * ib**2) * x**6) / r0

  End Function dmm3

  Function intRaddMM3(r0, A, B, eps, rw, prec)

    Real(kind=wp), Intent(In   ) :: r0, A, B, eps, rw, prec
    Real(Kind=wp)                :: intRaddMM3

    Integer       :: i, j, n
    Real(Kind=wp) :: h, ie, is, s, sold, x

    n = 10000
    is = rw / r0
    ie = 2 * is
    h = (ie - is) / Real(n, wp)
    s = 0.0_wp
    sold = Huge(1.0_wp)
    j = 1

    Do While (Abs(s - sold) * h / 3.0_wp > prec)
      sold = s

      Do i = (j - 1) * n, j * n, 2
        x = is + i * h

        s = s + x**3 * dmm3(x, r0, A, B, eps) + &
            4.0_wp * (x + h)**3 * dmm3(x + h, r0, A, B, eps) + &
            (x + 2.0_wp * h)**3 * dmm3(x + 2 * h, r0, A, B, eps)
      End Do

      j = j + 1
    End Do

    intRaddMM3 = s * h * r0**4 / 3.0_wp

  End Function intRaddMM3

  Subroutine vdw_lrc(sites, vdws, config, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to evaluate vdw long-range corrections to
    ! pressure and energy in a 3D periodic system
    !
    ! copyright - daresbury laboratory
    ! author    - t.forester may 1993
    ! amended   - i.t.todorov september 2016
    ! contrib   - a.m.elena september 2016 (ljc)
    ! contrib   - a.m.elena september 2017 (rydberg)
    ! contrib   - a.m.elena october 2017 (zbl/zbls)
    ! contrib   - a.m.elena december 2017 (zblb)
    ! contrib   - a.m.elena may 2018 (mdf)
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(site_type),          Intent(In   ) :: sites
    Type(vdw_type),           Intent(InOut) :: vdws
    Type(configuration_type), Intent(InOut) :: config
    Type(comms_type),         Intent(InOut) :: comm

    Character(Len=STR_LEN)                       :: messages(3)
    Integer                                  :: fail, i, ivdw, j, k, keypot
    Real(Kind=wp)                            :: a, b, c, d, denprd, e0, eadd, eps, kk, mr, nr, &
                                                padd, plrc, r, r0, s9, sig, t, z1, z2
    Real(Kind=wp), Allocatable, Dimension(:) :: numfrz
    Real( kind = wp ) :: v
    Character(Len=STR_LEN) :: unit

    fail = 0
    Allocate (numfrz(sites%mxatyp), Stat=fail)
    If (fail > 0) Call error_alloc('numfrz', 'vdw_lrc')

    ! initialise long-range corrections to energy and pressure

    plrc = 0.0_wp
    vdws%elrc = 0.0_wp

    If (.not. vdws%l_force_shift) Then

      ! initialise counter arrays and evaluate number density in system

      numfrz = 0.0_wp
      Do i = 1, config%natms
        k = config%ltype(i)
        If (config%lfrzn(i) /= 0) numfrz(k) = numfrz(k) + 1.0_wp
      End Do
      Call gsum(comm, numfrz(1:sites%ntype_atom))

      ! Evaluate only for 3D periodic systems

      If (config%imcon /= 0 .and. config%imcon /= 6) Then
        ivdw = 0

        Do i = 1, sites%ntype_atom
          Do j = 1, i

            eadd = 0.0_wp
            padd = 0.0_wp

            ivdw = ivdw + 1
            k = vdws%list(ivdw)

            keypot = vdws%ltp(k)

            Select case (keypot)
            Case (VDW_TAB)

              ! tabulated energy and pressure lrc

              eadd = vdws%param(1, k)
              padd = -vdws%param(2, k)

            Case (VDW_12_6)

              ! 12-6 potential :: u=a/r^12-b/r^6

              a = vdws%param(1, k)
              b = vdws%param(2, k)
              r = vdws%cutoff

              eadd = a / (9.0_wp * r**9) - b / (3.0_wp * r**3)
              padd = 12.0_wp * a / (9.0_wp * r**9) - 6.0_wp * b / (3.0_wp * r**3)

            Case (VDW_LENNARD_JONES)

              ! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps = vdws%param(1, k)
              sig = vdws%param(2, k)
              r = vdws%cutoff

              eadd = 4.0_wp * eps * (sig**12 / (9.0_wp * r**9) - sig**6 / (3.0_wp * r**3))
              padd = 8.0_wp * eps * (6.0_wp * sig**12 / (9.0_wp * r**9) - sig**6 / (r**3))

            Case (VDW_N_M)

              ! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(r0/r)^m]

              e0 = vdws%param(1, k)
              nr = vdws%param(2, k)
              mr = vdws%param(3, k)
              r0 = vdws%param(4, k)
              r = vdws%cutoff

              eadd = e0 / (nr - mr) * &
                   (mr * r0**nr / ((nr - 3.0_wp) * r**(nr - 3.0_wp)) - nr * r0**mr / ((mr - 3.0_wp) * r**(mr - 3.0_wp)))
              padd = e0 / (nr - mr) * &
                   nr * mr * (r0**nr / ((nr - 3.0_wp) * r**(nr - 3.0_wp)) - r0**mr / ((mr - 3.0_wp) * r**(mr - 3.0_wp)))

            Case (VDW_BUCKINGHAM)

              ! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              c = vdws%param(3, k)
              r = vdws%cutoff

              eadd = -c / (3.0_wp * r**3)
              padd = -2.0_wp * c / (r**3)

            Case (VDW_BORN_HUGGINS_MEYER)

              ! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              c = vdws%param(4, k)
              d = vdws%param(5, k)
              r = vdws%cutoff

              eadd = -c / (3.0_wp * r**3) - d / (5.0_wp * r**5)
              padd = -2.0_wp * c / (r**3) - 8.0_wp * d / (5.0_wp * r**5)

            Case (VDW_HYDROGEN_BOND)

              ! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a = vdws%param(1, k)
              b = vdws%param(2, k)
              r = vdws%cutoff

              eadd = a / (9.0_wp * r**9) - b / (7.0_wp * r**7)
              padd = 12.0_wp * a / (9.0_wp * r**9) - 10.0_wp * b / (7.0_wp * r**7)

            Case (VDW_MORSE)

              ! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

              e0 = vdws%param(1, k)
              r0 = vdws%param(2, k)
              kk = vdws%param(3, k)
              If (kk > Tiny(kk)) Then
                t = Exp(-kk * (vdws%cutoff - r0))

                eadd = -2.0_wp * e0 * t / (kk * kk * kk) * ((kk * vdws%cutoff + 1)**2 + 1) + &
                     e0 * t * t / (4.0_wp * kk * kk * kk) * ((kk * vdws%cutoff + 1)**2 + kk * kk * vdws%cutoff * vdws%cutoff)
                padd = -2.0_wp * e0 * t / (kk * kk * kk) * (kk**3 * vdws%cutoff**3 + &
                     3 * kk**2 * vdws%cutoff**2 + 6 * kk * vdws%cutoff + 6) + &
                     e0 * t * t / (4.0_wp * kk * kk * kk) * &
                     (4.0_wp * kk**3 * vdws%cutoff**3 + 6 * kk**2 * vdws%cutoff**2 + 6 * kk * vdws%cutoff + 3)
              End If

            Case (VDW_AMOEBA)

              ! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

              eps = vdws%param(1, k)
              sig = vdws%param(2, k)

              a = 0.07_wp
              b = 0.12_wp
              e0 = 1.0e-12_wp

              eadd = intRadMM3(sig, a, b, eps, vdws%cutoff, e0)
              padd = -intRaddMM3(sig, a, b, eps, vdws%cutoff, e0)

            Case (VDW_LENNARD_JONES_COHESIVE)

              ! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

              eps = vdws%param(1, k)
              sig = vdws%param(2, k)
              c = vdws%param(3, k)
              r = vdws%cutoff

              eadd = 4.0_wp * eps * (sig**12 / (9.0_wp * r**9) - c * sig**6 / (3.0_wp * r**3))
              padd = 8.0_wp * eps * (6.0_wp * sig**12 / (9.0_wp * r**9) - c * sig**6 / (r**3))

            Case (VDW_MORSE_12)
              ! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

              e0 = vdws%param(1, k)
              r0 = vdws%param(2, k)
              kk = vdws%param(3, k)
              c = vdws%param(4, k)

              If (kk > Tiny(kk)) Then

                t = Exp(-kk * (vdws%cutoff - r0))
                s9 = c / (9.0_wp * vdws%cutoff**9)

                eadd = -2.0_wp * e0 * t / (kk * kk * kk) * ((kk * vdws%cutoff + 1)**2 + 1) + &
                     e0 * t * t / (4.0_wp * kk * kk * kk) * ((kk * vdws%cutoff + 1)**2 + &
                     kk * kk * vdws%cutoff * vdws%cutoff) + s9
                padd = -2.0_wp * e0 * t / (kk * kk * kk) * (kk**3 * vdws%cutoff**3 + &
                     3 * kk**2 * vdws%cutoff**2 + 6 * kk * vdws%cutoff + 6) + &
                     e0 * t * t / (4.0_wp * kk * kk * kk) * (4.0_wp * kk**3 * vdws%cutoff**3 + &
                     6 * kk**2 * vdws%cutoff**2 + 6 * kk * vdws%cutoff + 3) + 12.0_wp * s9
              End If

            Case (VDW_RYDBERG)

              ! Rydberg potential:: u=(a+b*r)Exp(-r/c)

              a = vdws%param(1, k)
              b = vdws%param(2, k)
              c = vdws%param(3, k)
              t = Exp(-vdws%cutoff / c)

            eadd = (b * c * vdws%cutoff**3 + (3 * b * c**2 + a * c) * vdws%cutoff**2 + &
                 (6 * b * c**3 + 2 * a * c**2) * vdws%cutoff + 6 * b * c**4 + 2 * a * c**3) * t
            padd = (b * vdws%cutoff**4 + (3 * b * c + a) * vdws%cutoff**3 + (9 * b * c**2 + 3 * a * c) * vdws%cutoff**2 + &
                    (18 * b * c**3 + 6 * a * c**2) * vdws%cutoff + 18 * b * c**4 + 6 * a * c**3) * t

            Case (VDW_ZBL)

              ! ZBL potential:: u=Z1Z2/(40r)_{i=1}^4b_ie^{-c_i*r/a}

              z1 = vdws%param(1, k)
              z2 = vdws%param(2, k)

              ! this is in fact inverse a
              a = (z1**0.23_wp + z2**0.23_wp) / (0.52917721067_wp * 0.88534_wp)
              kk = z1 * z2 * r4pie0
              eadd = intRadZBL(kk, a, vdws%cutoff, 1e-12_wp)
              padd = intdRadZBL(kk, a, vdws%cutoff, 1e-12_wp)

            Case (VDW_ZBL_SWITCH_MORSE)

              ! ZBL switched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

              e0 = vdws%param(5, k)
              r0 = vdws%param(6, k)
              kk = vdws%param(7, k)

              If (kk > Tiny(kk)) Then
                t = Exp(-kk * (vdws%cutoff - r0))

                eadd = -2.0_wp * e0 * t / (kk * kk * kk) * ((kk * vdws%cutoff + 1)**2 + 1) + &
                     e0 * t * t / (4.0_wp * kk * kk * kk) * ((kk * vdws%cutoff + 1)**2 + kk * kk * vdws%cutoff * vdws%cutoff)
                padd = -2.0_wp * e0 * t / (kk * kk * kk) * (kk**3 * vdws%cutoff**3 + &
                  3 * kk**2 * vdws%cutoff**2 + 6 * kk * vdws%cutoff + 6) + &
                  e0 * t * t / (4.0_wp * kk * kk * kk) * &
                  (4.0_wp * kk**3 * vdws%cutoff**3 + 6 * kk**2 * vdws%cutoff**2 + 6 * kk * vdws%cutoff + 3)
              End If

             Case (VDW_ZBL_SWITCH_BUCKINGHAM)

              ! ZBL switched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

              A = vdws%param(5, k)
              r0 = vdws%param(6, k)
              c = vdws%param(7, k)

              t = A * Exp(-vdws%cutoff / r0)

              eadd = (vdws%cutoff**2 + 2 * r0 * vdws%cutoff + 2 * r0**2) * t * r0 - c / (3.0_wp * vdws%cutoff**3)
              padd = (vdws%cutoff**3 + 3 * r0 * vdws%cutoff**2 + 6 * r0**2 * vdws%cutoff + 6 * r0**3) * t - &
                2.0_wp * c / (vdws%cutoff**3)

             Case (VDW_LJ_MDF)

              ! LJ tapered with MDF:: u=f(r)LJ(r)

              a = vdws%param(1, k)
              b = vdws%param(2, k)
              c = vdws%param(3, k)
              eadd = intRadMDF("mlj", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)
              padd = intdRadMDF("mlj", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)

             Case (VDW_BUCKINGHAM_MDF)

              ! Buckingham tapered with MDF:: u=f(r)Buck(r)

              a = vdws%param(1, k)
              b = vdws%param(2, k)
              c = vdws%param(3, k)
              r0 = vdws%param(4, k)
              eadd = intRadMDF("mbuc", a, b, c, r0, vdws%cutoff, 1e-12_wp)
              padd = intdRadMDF("mbuc", a, b, c, r0, vdws%cutoff, 1e-12_wp)

             Case (VDW_126_MDF)

              ! LJ tapered with MDF:: u=f(r)LJ12-6(r)

              a = vdws%param(1, k)
              b = vdws%param(2, k)
              c = vdws%param(3, k)
              eadd = intRadMDF("m126", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)
              padd = intdRadMDF("m126", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)

             Case (VDW_LJF)
              eadd = 0.0_wp
              padd = 0.0_wp

             Case (VDW_SANDERSON)
              eadd = 0.0_wp ! implement me
              padd = 0.0_wp
             Case (VDW_SW)
              eadd = 0.0_wp ! implement me
              padd = 0.0_wp
            End Select

            ! Self-interaction accounted once, interaction between different species
            ! MUST be accounted twice!!

            If (i /= j) Then
              eadd = eadd * 2.0_wp
              padd = padd * 2.0_wp
            End If

            denprd = twopi * (sites%num_type(i) * sites%num_type(j) - numfrz(i) * numfrz(j)) / config%volm**2

            vdws%elrc = vdws%elrc + config%volm * denprd * eadd
            plrc = plrc + denprd * padd / 3.0_wp

          End Do
        End Do

      End If
    End If

    Write (messages(1), '(a)') 'long-range correction for:'
    Call to_out_units(vdws%elrc,'internal_e',v,unit)
    Write (messages(2), '(2x,a,e15.6)') 'vdw energy('//Trim(unit)//'): ', v
    Call to_out_units(plrc,'internal_p',v,unit)
    Write (messages(3), '(2x,a,e15.6)') 'vdw pressure('//Trim(unit)//'): ', v
    Call info(messages, 3, .true.)
    ! convert plrc to a viral term

    vdws%vlrc = plrc * (-3.0_wp * config%volm)

    Deallocate (numfrz, Stat=fail)
    if (fail /= 0) Call error_dealloc('numfrz', 'vdw_lrc')

  End Subroutine vdw_lrc

  Subroutine vdw_direct_fs_generate(vdws)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for generating force-shifted constant arrays for
    ! direct vdw evaluation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! contrib   - a.m.elena september 2016 (ljc)
    ! contrib   - a.m.elena september 2017 (ryd)
    ! contrib   - a.m.elena october 2017 (zbl/zbls)
    ! contrib   - a.m.elena december 2017 (zblb)
    ! contrib   - a.m.elena april 2018 (mlj/mbuc)
    ! contrib   - a.m.elena may 2018 (m126)
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(vdw_type), Intent(InOut) :: vdws
    Integer                       :: ivdw, keypot
    Real(Kind=wp)                 :: z, dz
    Type(potential_energy)        :: z_dz

    ! allocate arrays for force-shifted corrections

    Call vdws%init_direct()
    if (.not. vdws%l_force_shift) Return

    ! construct arrays for all types of vdw potential

    Do ivdw = 1, vdws%n_vdw

      z_dz = vdws%potentials(ivdw)%p%energy(vdws%cutoff)
      z = z_dz%energy
      dz = z_dz%gamma

      keypot = vdws%ltp(ivdw)

      Select Case (keypot)

        ! a few special cases

        Case (VDW_N_M_SHIFT)

          ! shifted and force corrected n-m potential (w.smith) ::
          z = 0.0_wp
          dz = 0.0_wp

        Case (VDW_DPD) ! all zeroed in vdw

          ! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2
          z = 0.0_wp
          dz = 0.0_wp

        Case (VDW_NDPD)

          z = 0.0_wp
          dz = 0.0_wp

        Case (VDW_LJF)

          z = 0.0_wp
          dz = 0.0_wp

        Case (VDW_126_MDF)

          z = 0.0_wp
          dz = 0.0_wp

      End select

      vdws%afs(ivdw) = dz / vdws%cutoff
      vdws%bfs(ivdw) = -z - dz

    End Do

  End Subroutine vdw_direct_fs_generate

  Subroutine vdw_table_read(vdws, sites, files, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for reading potential energy and force arrays
    ! from TABLE file (for van der waals forces only)
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith march 1994
    ! amended   - i.t.todorov april 2016
    ! amended   - a.m.elena january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(vdw_type),   Intent(InOut) :: vdws
    Type(site_type),  Intent(In   ) :: sites
    Type(file_type),  Intent(InOut) :: files(:)
    Type(comms_type), Intent(InOut) :: comm

    Character(Len=200)                       :: record
    Character(Len=STR_LEN)                       :: message, messages(4)
    Character(Len=40)                        :: word
    Character(Len=8)                         :: atom1, atom2
    Integer                                  :: fail, i, ivdw, j, jtpatm, katom1, katom2, keyvdw, &
                                                l, ngrid, ntable, ierror
    Logical                                  :: remake, safe
    Real(Kind=wp)                            :: cutpot, delpot, dlrpot, ppp, rdr, rrr, t1, t2, &
                                                vk, vk1, vk2
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer


    If (comm%idnode == 0) Then
      Open (Newunit=files(FILE_TABVDW)%unit_no, File=files(FILE_TABVDW)%filename)
      ntable = files(FILE_TABVDW)%unit_no
    End If

    ! skip header record

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Call ferror(24)

    ! read mesh resolution

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Call ferror(24)

    Call get_word(record, word)
    delpot = word_2_real(word)

    Call get_word(record, word)
    cutpot = word_2_real(word)

    Call get_word(record, word)
    ngrid = Nint(word_2_real(word))

    dlrpot = vdws%cutoff / Real(vdws%max_grid - 4, wp)

    ! check grid spacing

    safe = .false.
    If (Abs(delpot - dlrpot) <= 1.0e-8_wp) Then
      safe = .true.
      delpot = dlrpot
    End If
    If (delpot > delr_max .and. (.not. safe)) Then
      Write (messages(1), '(a,1p,e15.7)') 'expected (maximum) radial increment: ', delr_max
      Write (messages(2), '(a,1p,e15.7)') 'TABLE  file actual radial increment: ', delpot
      Write (messages(3), '(a,i10)') 'expected (minimum) number of grid points: ', vdws%max_grid
      Write (messages(4), '(a,i10)') 'TABLE  file actual number of grid points: ', ngrid
      Call info(messages, 4, .true.)
      Call error(22)
    End If
    safe = .true.

    remake = .false.
    If (Abs(1.0_wp - (delpot / dlrpot)) > 1.0e-8_wp) Then
      remake = .true.
      rdr = 1.0_wp / delpot
      Write (message, '(a,i10)') 'TABLE arrays resized for mxgrid = ', vdws%max_grid - 4
      Call info(message, .true.)
    End If

    ! compare grids dimensions

    If (ngrid < vdws%max_grid - 4) Then
      Call warning(270, Real(ngrid, wp), Real(vdws%max_grid - 4, wp), 0.0_wp)
      Call error(0, 'Transfer buffer too small in vdw_table_read')
    End If

    If (cutpot < vdws%cutoff) Call error(0, 'Cutoff too large for TABLE file')

    fail = 0
    Allocate (buffer(0:ngrid), Stat=fail)
    If (fail > 0) call error_alloc('buffer', 'vdw_table_read')

    ! read potential arrays for all pairs

    Do ivdw = 1, vdws%n_vdw

      ! read potential arrays if potential not already defined

      If (vdws%ltp(ivdw) == VDW_TAB) Then

        ! read pair potential labels and long-range corrections

        Call get_line(safe, ntable, record, comm)
        If (.not. safe) Call ferror(24)

        Call get_word(record, atom1)
        Call get_word(record, atom2)

        Call get_word(record, word)
        vdws%param(1, ivdw) = word_2_real(word) * engunit

        Call get_word(record, word)
        vdws%param(2, ivdw) = word_2_real(word) * engunit

        katom1 = 0
        katom2 = 0

        Do jtpatm = 1, sites%ntype_atom
          If (atom1 == sites%unique_atom(jtpatm)) katom1 = jtpatm
          If (atom2 == sites%unique_atom(jtpatm)) katom2 = jtpatm
        End Do

        If (katom1 == 0 .or. katom2 == 0) Then
          Write (message, '(a,a,a,a,a)') '****', atom1, '***', atom2, '**** entry in TABLE'
          Call error(81, message, .true.)
        End If

        keyvdw = (Max(katom1, katom2) * (Max(katom1, katom2) - 1)) / 2 + Min(katom1, katom2)

        ! Only one vdw potential per pair is allowed
        ! (FIELD AND TABLE potentials overlapping)

        If (vdws%list(keyvdw) /= ivdw) &
             Call error(0, 'Incompatible FIELD and TABLE file potentials')

        ! read in potential arrays

        Do i = 1, (ngrid + 3) / 4
          j = Min(4, ngrid - (i - 1) * 4)
          If (comm%idnode == 0) Then
            Read (Unit=ntable, Fmt=*, iostat=ierror) buffer((i - 1) * 4 + 1:(i - 1) * 4 + j)
            If ( ierror > 0 ) Call ferror(24)
          Else
            buffer((i - 1) * 4 + 1:(i - 1) * 4 + j) = 0.0_wp
          End If
        End Do
        Call gbcast(comm, buffer, 0)
        ! linear extrapolation for grid point 0 at distances close to 0

        vdws%tab_potential(0, ivdw) = 2.0_wp * buffer(1) - buffer(2)

        ! reconstruct arrays using 3pt interpolation

        If (remake) Then
          Do i = 1, vdws%max_grid - 4
            rrr = Real(i, wp) * dlrpot
            l = Int(rrr * rdr)
            ppp = rrr * rdr - Real(l, wp)

            vk = buffer(l)

            ! linear extrapolation for the grid points just beyond the cutoff

            If (l + 2 > ngrid) Then
              If (l + 1 > ngrid) Then
                vk1 = 2.0_wp * buffer(l) - buffer(l - 1)
                vk2 = 2.0_wp * vk1 - buffer(l)
              Else
                vk1 = buffer(l + 1)
                vk2 = 2.0_wp * buffer(l + 1) - buffer(l)
              End If
            Else
              vk1 = buffer(l + 1)
              vk2 = buffer(l + 2)
            End If

            t1 = vk + (vk1 - vk) * ppp
            t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)
            vdws%tab_potential(i, ivdw) = t1 + (t2 - t1) * ppp * 0.5_wp
          End Do
        Else
          Do i = 1, vdws%max_grid - 4
            vdws%tab_potential(i, ivdw) = buffer(i)
          End Do

          ! linear extrapolation for the grid point just beyond the cutoff

          vdws%tab_potential(vdws%max_grid - 3, ivdw) = 2.0_wp * vdws%tab_potential(vdws%max_grid - 4, ivdw) - &
                                                        vdws%tab_potential(vdws%max_grid - 5, ivdw)
        End If

        ! linear extrapolation for the grid point at vdws%max_grid-2

        vdws%tab_potential(vdws%max_grid - 2, ivdw) = 2.0_wp * vdws%tab_potential(vdws%max_grid - 3, ivdw) - &
                                                      vdws%tab_potential(vdws%max_grid - 4, ivdw)

        ! read in force arrays

        Do i = 1, (ngrid + 3) / 4
          j = Min(4, ngrid - (i - 1) * 4)
          If (comm%idnode == 0) Then
            Read (Unit=ntable, Fmt=*, Iostat = ierror) buffer((i - 1) * 4 + 1:(i - 1) * 4 + j)
            If ( ierror > 0 ) Call ferror(24)
          Else
            buffer((i - 1) * 4 + 1:(i - 1) * 4 + j) = 0.0_wp
          End If
        End Do
        Call gbcast(comm, buffer, 0)
        ! linear extrapolation for grid point 0 at distances close to 0

        vdws%tab_force(0, ivdw) = (2.0_wp * buffer(1) - 0.5_wp * buffer(2)) / delpot

        ! reconstruct arrays using 3pt interpolation

        If (remake) Then
          Do i = 1, vdws%max_grid - 4
            rrr = Real(i, wp) * dlrpot
            l = Int(rrr * rdr)
            ppp = rrr * rdr - Real(l, wp)

            vk = buffer(l)

            ! linear extrapolation for the grid points just beyond the cutoff

            If (l + 2 > ngrid) Then
              If (l + 1 > ngrid) Then
                vk1 = 2.0_wp * buffer(l) - buffer(l - 1)
                vk2 = 2.0_wp * vk1 - buffer(l)
              Else
                vk1 = buffer(l + 1)
                vk2 = 2.0_wp * buffer(l + 1) - buffer(l)
              End If
            Else
              vk1 = buffer(l + 1)
              vk2 = buffer(l + 2)
            End If

            t1 = vk + (vk1 - vk) * ppp
            t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

            vdws%tab_force(i, ivdw) = t1 + (t2 - t1) * ppp * 0.5_wp
          End Do
        Else
          Do i = 1, vdws%max_grid - 4
            vdws%tab_force(i, ivdw) = buffer(i)
          End Do

          ! linear extrapolation for the grid point just beyond the cutoff

          vdws%tab_force(vdws%max_grid - 3, ivdw) = 2.0_wp * vdws%tab_force(vdws%max_grid - 4, ivdw) - &
                                                    vdws%tab_force(vdws%max_grid - 5, ivdw)
        End If

        ! linear extrapolation for the grid point at vdws%max_grid-2

        vdws%tab_force(vdws%max_grid - 2, ivdw) = 2.0_wp * vdws%tab_force(vdws%max_grid - 3, ivdw) - &
                                                  vdws%tab_force(vdws%max_grid - 4, ivdw)

        ! We must distinguish that something has been defined

        If (Abs(vdws%tab_potential(0, ivdw)) <= zero_plus) Then
          vdws%tab_potential(0, ivdw) = Sign(Tiny(vdws%tab_potential(0, ivdw)), vdws%tab_potential(0, ivdw))
        End If

      End If

    End Do

    Call info('potential tables read from TABLE file', .true.)
    If (comm%idnode == 0) Call files(FILE_TABVDW)%close ()

    ! convert to internal units

    Do ivdw = 1, vdws%n_vdw
      If (vdws%ltp(ivdw) == VDW_TAB) Then

        Do i = 0, vdws%max_grid
          vdws%tab_potential(i, ivdw) = vdws%tab_potential(i, ivdw) * engunit
          vdws%tab_force(i, ivdw) = vdws%tab_force(i, ivdw) * engunit

        End Do
      End If
    End Do

    If (vdws%l_force_shift) Then
      Do ivdw = 1, vdws%n_vdw
        If (vdws%ltp(ivdw) == VDW_TAB) Then

          vdws%tab_potential(vdws%max_grid - 3, ivdw) = 0.0_wp
          vdws%tab_potential(vdws%max_grid - 2, ivdw) = 0.0_wp
          vdws%tab_force(vdws%max_grid - 3, ivdw) = 0.0_wp
          vdws%tab_force(vdws%max_grid - 2, ivdw) = 0.0_wp
        End If
      End Do
    End If

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Call error_dealloc('buffer', 'vdw_table_read')

    ! end of file error exit

    Contains

      Subroutine ferror(code)
        Integer, Intent(In) :: code

        If (comm%idnode == 0) Call files(FILE_TABVDW)%close ()
        Call error(code)

      End Subroutine ferror

  End Subroutine vdw_table_read

  Subroutine dump_vdws(T,comm)
    Class(vdw_type)  :: T
    Type(comms_type), Intent(In) :: comm

    integer :: i, j,u
    Real(Kind = wp ) :: r,dlrpot
    character(len=100) :: filename

    If (comm%idnode ==  0 ) Then


    dlrpot = T%cutoff / Real(T%max_grid - 4, wp)
    Do i = 1, T%n_vdw
        write(filename,'(a,i0,a)')"vdw_",i,".dat"
        Open(Newunit=u,File=Trim(filename),Action="write",Status='unknown')
        Do j = 1, T%max_grid
          r = Real(j, wp) * dlrpot
          Write(u,*)r, T%tab_potential(j, i)/engunit, T%tab_force(j, i)/engunit
        End Do
        Close(u)
    End Do
  End If

  End Subroutine dump_vdws

  Subroutine vdw_generate(vdws)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for generating potential energy and force arrays
    ! for van der waals forces only
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith may 1992
    ! amended   - i.t.todorov march 2016
    ! contrib   - a.m.elena september 2016 (ljc)
    ! contrib   - a.m.elena september 2017 (rydberg)
    ! contrib   - a.m.elena october 2017 (zbl/zbls)
    ! contrib   - a.m.elena december 2017 (zblb)
    ! contrib   - a.m.elena april 2018 (mlj/mbuc)
    ! contrib   - a.m.elena may 2018 (m126)
    ! contrib   - m.a.seaton november 2023 (ndpd)
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(vdw_type), Intent(InOut) :: vdws

    Integer                       :: i, ivdw, keypot
    Real(Kind=wp)                 :: dlrpot, e0, kk, r, r0, t, t1
    Type(potential_energy)        :: z_dz
    Type(mlj)                     :: p_mlj
    Type(mbuck)                   :: p_mbuck
    Type(mlj126)                  :: p_mlj126
    Real(Kind=wp)                 :: param(1:5)

    ! allocate arrays for tabulating

    Call vdws%init_table()

    ! define grid resolution for potential arrays

    dlrpot = vdws%cutoff / Real(vdws%max_grid - 4, wp)

    ! construct arrays for all types of vdw potential

    Do ivdw = 1, vdws%n_vdw

      keypot = vdws%ltp(ivdw)

      If (keypot /= VDW_TAB) Then

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          z_dz = vdws%potentials(ivdw)%p%energy(r)
          vdws%tab_potential(i, ivdw) = z_dz%energy
          vdws%tab_force(i, ivdw) = z_dz%gamma

        End Do

        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))
      
      End If

      select case (keypot)

      Case (VDW_MORSE)

        e0 = vdws%param(1, ivdw)
        r0 = vdws%param(2, ivdw)
        kk = vdws%param(3, ivdw)
        t1 = Exp(+kk * r0)
        vdws%tab_force(0, ivdw) = -2.0_wp * e0 * kk * (1.0_wp - t1) * t1

      Case (VDW_DPD)

        ! DPD potential - Groot-Warren (standard) :: u=(1/2).a.rc.(1-r/rc)^2

        vdws%tab_force(0, ivdw) = vdws%param(1, ivdw)

      Case (VDW_NDPD)

        vdws%tab_force(0, ivdw) = vdws%param(1, ivdw) * (vdws%param(2, ivdw))

      Case (VDW_RYDBERG)

        ! Rydberg potential:: u=(a+b*r)Exp(-r/c)

        vdws%tab_potential(0, ivdw) = vdws%param(1, ivdw)
        vdws%tab_force(0, ivdw) = 0.0_wp

      Case (VDW_LJ_MDF)

        ! LJ tapered with MDF:: u=f(r)LJ(r)
        param(1:3) = vdws%param(1:3, ivdw)
        param(4) = vdws%cutoff
        Call p_mlj%set_parameters(param(1:4))

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          z_dz = p_mlj%energy(r)

          vdws%tab_potential(i, ivdw) = z_dz%energy
          vdws%tab_force(i, ivdw) = z_dz%gamma

        End Do
        
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Case (VDW_BUCKINGHAM_MDF)

        ! Buckingham tapered with MDF:: u=f(r)Buck(r)

        param(1:4) = vdws%param(1:4, ivdw)
        param(5) = vdws%cutoff
        
        Call p_mbuck%set_parameters(param(1:5))

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          z_dz = p_mbuck%energy(r)

          vdws%tab_potential(i, ivdw) = z_dz%energy
          vdws%tab_force(i, ivdw) = z_dz%gamma

        End Do

        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Case (VDW_126_MDF)

        ! LJ tapered with MDF:: u=f(r)LJ(r)

        param(1:3) = vdws%param(1:3, ivdw)
        param(4) = vdws%cutoff

        Call p_mlj126%set_parameters(param(1:4))

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          z_dz = p_mlj126%energy(r)

          vdws%tab_potential(i, ivdw) = z_dz%energy
          vdws%tab_force(i, ivdw) = z_dz%gamma

        End Do

        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))
        
      End select

      ! no shifting to shifted n-m, DPD and nDPD
      If (vdws%l_force_shift .and. (keypot /= VDW_N_M_SHIFT .and. keypot /= VDW_DPD .and. keypot /= VDW_NDPD)) Then

        Do i = 1, vdws%max_grid - 4
          t = vdws%tab_potential(i, ivdw) + &
              vdws%tab_force(vdws%max_grid - 4, ivdw) * (Real(i, wp) * dlrpot / vdws%cutoff - 1.0_wp) - &
              vdws%tab_potential(vdws%max_grid - 4, ivdw)
          t1 = vdws%tab_potential(i - 1, ivdw) + &
               vdws%tab_force(vdws%max_grid - 4, ivdw) * (Real(i - 1, wp) * dlrpot / vdws%cutoff - 1.0_wp) - &
               vdws%tab_potential(vdws%max_grid - 4, ivdw)

        End Do
      End If

      ! Needed to distinguish that something has been defined

      If (Abs(vdws%tab_potential(0, ivdw)) <= zero_plus) Then
        vdws%tab_potential(0, ivdw) = Sign(Tiny(vdws%tab_potential(0, ivdw)), vdws%tab_potential(0, ivdw))
      End If

    End Do

  End Subroutine vdw_generate

  Subroutine vdw_forces_direct(iatm, xxt, yyt, zzt, rrt, engvdw, virvdw, stats, neigh, vdws, config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating vdw energy and force terms using
    ! verlet neighbour list
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1998
    ! amended   - i.t.todorov march 2016
    ! contrib   - a.m.elena september 2016 (ljc)
    ! contrib   - a.m.elena september 2017 (rydberg)
    ! contrib   - a.m.elena october 2017 (zbl/zbls)
    ! contrib   - a.m.elena december 2017 (zblb)
    ! contrib   - a.m.elena april 2018 (mlj/mbuc)
    ! contrib   - a.m.elena may 2018 (m126)
    ! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    ! contrib   - m.a.seaton november 2023 (ndpd)
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                    Intent(In   ) :: iatm
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(stats_type),                           Intent(InOut) :: stats
    Real(Kind=wp),                              Intent(  Out) :: virvdw, engvdw
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: rrt, zzt, yyt, xxt
    Type(vdw_type),                             Intent(InOut) :: vdws
    Type(configuration_type),                   Intent(InOut) :: config

    Integer                 :: ai, aj, idi, ityp, jatm, k, key, mm
    Real(Kind=wp)           :: eng, gamma
    Real(Kind=wp)           :: fix, fiy, fiz, fx, fy, fz
    Real(Kind=wp)           :: r_rrr, r_rrv, r_rsq, r_rvdw, rrr, rscl, rsq
    Real(Kind=wp)           :: strs1, strs2, strs3, strs5, strs6, strs9
    Real(Kind=wp)           :: stress_temp_comp(9)
    Real(Kind=wp)           :: x_temp(3), f_temp(3)
    Type(potential_energy)  :: eng_gamma
    ! define grid resolution for potential arrays and interpolation spacing

    If (vdws%newjob) Then
      vdws%newjob = .false.

      vdws%dlrpot = vdws%cutoff / Real(vdws%max_grid - 4, wp)
      vdws%rdr = 1.0_wp / vdws%dlrpot
    End If

    ! initialise potential energy and virial

    engvdw = 0.0_wp
    virvdw = 0.0_wp

    ! initialise stress tensor accumulators

    strs1 = 0.0_wp
    strs2 = 0.0_wp
    strs3 = 0.0_wp
    strs5 = 0.0_wp
    strs6 = 0.0_wp
    strs9 = 0.0_wp

    stress_temp_comp = 0.0_wp
    ! global identity and type of iatm

    idi = config%ltg(iatm)
    ai = config%ltype(iatm)

    ! load forces

    fix = config%parts(iatm)%fxx
    fiy = config%parts(iatm)%fyy
    fiz = config%parts(iatm)%fzz

    ! start of primary loop for forces evaluation

    Do mm = 1, neigh%list(0, iatm)

      ! atomic and potential function indices

      jatm = neigh%list(mm, iatm)
      aj = config%ltype(jatm)

      If (ai > aj) Then
        key = ai * (ai - 1) / 2 + aj
      Else
        key = aj * (aj - 1) / 2 + ai
      End If

      k = vdws%list(key)

      ! interatomic distance

      rrr = rrt(mm)

      ! validity and truncation of potential

      ityp = vdws%ltp(k)

      If (ityp /= VDW_NULL .and. rrr < vdws%cutoff) Then

        ! Distance derivatives

        r_rrr = 1.0_wp / rrr
        r_rvdw = 1.0_wp / vdws%cutoff
        rsq = rrr**2
        r_rsq = r_rrr**2 !1.0_wp / rsq
        r_rrv = r_rrr * r_rvdw
        rscl = rrr * r_rvdw

        ! Zero energy and force components

        eng = 0.0_wp
        gamma = 0.0_wp

        eng_gamma = vdws%potentials(k)%energy(rrr)

        eng = eng_gamma%energy + vdws%afs(k) * rrr + vdws%bfs(k)
        gamma = eng_gamma%gamma*r_rsq - vdws%afs(k) * r_rrr

        ! calculate forces
        fx = gamma * xxt(mm)
        fy = gamma * yyt(mm)
        fz = gamma * zzt(mm)

        fix = fix + fx
        fiy = fiy + fy
        fiz = fiz + fz

#ifndef HALF_HALO
        If (jatm > config%natms .and. idi >= config%ltg(jatm) .and. .not. stats%collect_pp) &
             eng = 0.0_wp

        If (jatm <= config%natms) Then
#endif /* HALF_HALO */

          config%parts(jatm)%fxx = config%parts(jatm)%fxx - fx
          config%parts(jatm)%fyy = config%parts(jatm)%fyy - fy
          config%parts(jatm)%fzz = config%parts(jatm)%fzz - fz

#ifndef HALF_HALO
        End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
        If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

          ! add interaction energy

          engvdw = engvdw + eng

          ! add virial

          virvdw = virvdw - gamma * rsq

          ! add stress tensor

          strs1 = strs1 + xxt(mm) * fx
          strs2 = strs2 + xxt(mm) * fy
          strs3 = strs3 + xxt(mm) * fz
          strs5 = strs5 + yyt(mm) * fy
          strs6 = strs6 + yyt(mm) * fz
          strs9 = strs9 + zzt(mm) * fz

#ifndef HALF_HALO
        End If
#endif /* HALF_HALO */

        If (stats%collect_pp) Then
          x_temp = (/ xxt(mm), yyt(mm), zzt(mm) /)
          f_temp = (/ fx, fy, fz /)
          stress_temp_comp = calculate_stress( x_temp, f_temp  )
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + eng * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + eng * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */
        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = fix
    config%parts(iatm)%fyy = fiy
    config%parts(iatm)%fzz = fiz

    ! complete stress tensor

    stats%stress(1) = stats%stress(1) + strs1
    stats%stress(2) = stats%stress(2) + strs2
    stats%stress(3) = stats%stress(3) + strs3
    stats%stress(4) = stats%stress(4) + strs2
    stats%stress(5) = stats%stress(5) + strs5
    stats%stress(6) = stats%stress(6) + strs6
    stats%stress(7) = stats%stress(7) + strs3
    stats%stress(8) = stats%stress(8) + strs6
    stats%stress(9) = stats%stress(9) + strs9

  End Subroutine vdw_forces_direct

  Subroutine vdw_forces_tab(iatm, xxt, yyt, zzt, rrt, engvdw, virvdw, stats, neigh, vdws, config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating vdw energy and force terms using
    ! verlet neighbour list
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1998
    ! amended   - i.t.todorov march 2016
    ! contrib   - a.m.elena september 2016 (ljc)
    ! contrib   - a.m.elena september 2017 (rydberg)
    ! contrib   - a.m.elena october 2017 (zbl/zbls)
    ! contrib   - a.m.elena december 2017 (zblb)
    ! contrib   - a.m.elena april 2018 (mlj/mbuc)
    ! contrib   - a.m.elena may 2018 (m126)
    ! contrib   - a.v.brukhno & m.a.seaton august 2020 - 'half-halo' VNL
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                    Intent(In   ) :: iatm
    Type(neighbours_type),                      Intent(In   ) :: neigh
    Type(stats_type),                           Intent(InOut) :: stats
    Real(Kind=wp),                              Intent(  Out) :: virvdw, engvdw
    Real(Kind=wp), Dimension(1:neigh%max_list), Intent(In   ) :: rrt, zzt, yyt, xxt
    Type(vdw_type),                             Intent(InOut) :: vdws
    Type(configuration_type),                   Intent(InOut) :: config

    Real(Kind=wp), Dimension(9) :: stress_temp_comp
    Real(Kind=wp), Dimension(3) :: x_temp, f_temp
    Integer       :: ai, aj, idi, ityp, jatm, k, key, l, mm
    Real(kind=wp) :: t1, t2, vk, vk1, vk2, gk, gk1, gk2, ppp
    Real(kind=wp) :: fix, fiy, fiz, fx, fy, fz, gamma, eng
    Real(kind=wp) :: rrr, rsq, rscl, r_rvdw, r_rrv, r_rsq, r_rrr
    Real(Kind=wp) :: strs1, strs2, strs3, strs5, strs6, strs9

    ! define grid resolution for potential arrays and interpolation spacing

    If (vdws%newjob) Then
      vdws%newjob = .false.

      vdws%dlrpot = vdws%cutoff / Real(vdws%max_grid - 4, wp)
      vdws%rdr = 1.0_wp / vdws%dlrpot
    End If

    ! initialise potential energy and virial

    engvdw = 0.0_wp
    virvdw = 0.0_wp

    ! initialise stress tensor accumulators

    strs1 = 0.0_wp
    strs2 = 0.0_wp
    strs3 = 0.0_wp
    strs5 = 0.0_wp
    strs6 = 0.0_wp
    strs9 = 0.0_wp

    stress_temp_comp = 0.0_wp
    ! global identity and type of iatm

    idi = config%ltg(iatm)
    ai = config%ltype(iatm)

    ! load forces

    fix = config%parts(iatm)%fxx
    fiy = config%parts(iatm)%fyy
    fiz = config%parts(iatm)%fzz

    ! start of primary loop for forces evaluation

    Do mm = 1, neigh%list(0, iatm)

      ! atomic and potential function indices

      jatm = neigh%list(mm, iatm)
      aj = config%ltype(jatm)

      If (ai > aj) Then
        key = ai * (ai - 1) / 2 + aj
      Else
        key = aj * (aj - 1) / 2 + ai
      End If

      k = vdws%list(key)

      If (Abs(vdws%tab_potential(0, k)) < zero_plus) cycle

      ! interatomic distance

      rrr = rrt(mm)

      ! validity and truncation of potential

      ityp = vdws%ltp(k)
      If (ityp /= VDW_NULL .and. rrr < vdws%cutoff) Then

        ! Distance derivatives

        r_rrr = 1.0_wp / rrr
        r_rvdw = 1.0_wp / vdws%cutoff
        rsq = rrr**2
        r_rsq = r_rrr**2 !1.0_wp / rsq
        r_rrv = r_rrr * r_rvdw
        rscl = rrr * r_rvdw

        ! Zero energy and force components

        eng = 0.0_wp
        gamma = 0.0_wp


        l = Int(rrr * vdws%rdr)
        ppp = rrr * vdws%rdr - Real(l, wp)

        ! calculate forces using 3-point interpolation

        gk = vdws%tab_force(l, k); If (l == 0) gk = gk * rrr
        gk1 = vdws%tab_force(l + 1, k)
        gk2 = vdws%tab_force(l + 2, k)

        t1 = gk + (gk1 - gk) * ppp
        t2 = gk1 + (gk2 - gk1) * (ppp - 1.0_wp)

        gamma = (t1 + (t2 - t1) * ppp * 0.5_wp) * r_rsq
        If (vdws%l_force_shift) gamma = gamma - vdws%tab_force(vdws%max_grid - 4, k) * r_rrv ! force-shifting


        ! calculate forces

        fx = gamma * xxt(mm)
        fy = gamma * yyt(mm)
        fz = gamma * zzt(mm)

        fix = fix + fx
        fiy = fiy + fy
        fiz = fiz + fz

#ifndef HALF_HALO
        If (jatm <= config%natms) Then
#endif /* HALF_HALO */

          config%parts(jatm)%fxx = config%parts(jatm)%fxx - fx
          config%parts(jatm)%fyy = config%parts(jatm)%fyy - fy
          config%parts(jatm)%fzz = config%parts(jatm)%fzz - fz

#ifndef HALF_HALO
        End If
#endif /* HALF_HALO */

#ifndef HALF_HALO
        If (jatm <= config%natms .or. idi < config%ltg(jatm)) Then
#endif /* HALF_HALO */

          ! calculate interaction energy using 3-point interpolation

          vk = vdws%tab_potential(l, k)
          vk1 = vdws%tab_potential(l + 1, k)
          vk2 = vdws%tab_potential(l + 2, k)

          t1 = vk + (vk1 - vk) * ppp
          t2 = vk1 + (vk2 - vk1) * (ppp - 1.0_wp)

          eng = t1 + (t2 - t1) * ppp * 0.5_wp
          ! force-shifting
          If (vdws%l_force_shift) Then
            eng = eng + vdws%tab_force(vdws%max_grid - 4, k) * (rscl - 1.0_wp) - &
                 vdws%tab_potential(vdws%max_grid - 4, k)
          End If

          ! add interaction energy

          engvdw = engvdw + eng

          ! add virial

          virvdw = virvdw - gamma * rsq

          ! add stress tensor

          strs1 = strs1 + xxt(mm) * fx
          strs2 = strs2 + xxt(mm) * fy
          strs3 = strs3 + xxt(mm) * fz
          strs5 = strs5 + yyt(mm) * fy
          strs6 = strs6 + yyt(mm) * fz
          strs9 = strs9 + zzt(mm) * fz

#ifndef HALF_HALO
        End If
#endif /* HALF_HALO */
        If (stats%collect_pp) Then
          x_temp = (/ xxt(mm), yyt(mm), zzt(mm) /)
          f_temp = (/ fx, fy, fz /)
          stress_temp_comp = calculate_stress( x_temp, f_temp  )
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + eng * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
          If (jatm <= config%natms) Then
#endif /* HALF_HALO */
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + eng * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
#ifndef HALF_HALO
          End If
#endif /* HALF_HALO */
        End If
      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = fix
    config%parts(iatm)%fyy = fiy
    config%parts(iatm)%fzz = fiz

    ! complete stress tensor

    stats%stress(1) = stats%stress(1) + strs1
    stats%stress(2) = stats%stress(2) + strs2
    stats%stress(3) = stats%stress(3) + strs3
    stats%stress(4) = stats%stress(4) + strs2
    stats%stress(5) = stats%stress(5) + strs5
    stats%stress(6) = stats%stress(6) + strs6
    stats%stress(7) = stats%stress(7) + strs3
    stats%stress(8) = stats%stress(8) + strs6
    stats%stress(9) = stats%stress(9) + strs9

  End Subroutine vdw_forces_tab

  Subroutine cleanup(T)
    Type(vdw_type) :: T

    If (Allocated(T%list)) Then
      Deallocate (T%list)
    End If
    If (Allocated(T%ltp)) Then
      Deallocate (T%ltp)
    End If

    If (Allocated(T%param)) Then
      Deallocate (T%param)
    End If

    If (Allocated(T%labpair)) Then
      Deallocate (T%labpair)
    End If

    If (Allocated(T%tab_potential)) Then
      Deallocate (T%tab_potential)
    End If
    If (Allocated(T%tab_force)) Then
      Deallocate (T%tab_force)
    End If

    If (Allocated(T%afs)) Then
      Deallocate (T%afs)
    End If
    If (Allocated(T%bfs)) Then
      Deallocate (T%bfs)
    End If
  End Subroutine cleanup
End Module vdw
