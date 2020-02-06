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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms,           Only: comms_type,&
                             gbcast,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: delr_max,&
                             engunit,&
                             ntable,&
                             prsunt,&
                             r4pie0,&
                             twopi,&
                             zero_plus
  Use errors_warnings, Only: error,&
                             info,&
                             warning
  Use kinds,           Only: wi,&
                             wp
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: nequal
  Use parse,           Only: get_line,&
                             get_word,&
                             word_2_real
  Use site,            Only: site_type
  Use statistics,      Only: calculate_stress,&
                             stats_type

  Implicit None

  Private

  ! VdW potential parameters
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
  !> DPD potential
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
  !> ZBL swithched with Morse: $u=f(r)\mathrm{zbl}(r)+(1-f(r))*\mathrm{morse}(r)$
  Integer(Kind=wi), Parameter, Public :: VDW_ZBL_SWITCH_MORSE = 16
  !> ZBL swithched with Buckingham: $u=f(r)\mathrm{zbl}(r)+(1-f(r))*\mathrm{buckingham}(r)$
  Integer(Kind=wi), Parameter, Public :: VDW_ZBL_SWITCH_BUCKINGHAM = 17
  Integer(Kind=wi), Parameter, Public :: VDW_LJ_MDF = 18
  Integer(Kind=wi), Parameter, Public :: VDW_BUCKINGHAM_MDF = 19
  Integer(Kind=wi), Parameter, Public :: VDW_126_MDF = 20

  ! Mixing rule parameters
  !> Null
  Integer(Kind=wi), Parameter, Public :: MIX_NULL = 0
  !> Lorentz-Berthelot: $e_{ij}=(e_i*e_j)^{1/2} \quad s_{ij}=(s_i+s_j)/2$
  Integer(Kind=wi), Parameter, Public :: MIX_LORENTZ_BERTHELOT = 1
  !> Fender-Hasley: $e_{ij}=(2*e_i*e_j)/(e_i+e_j) \quad s_{ij}=(s_i+s_j)/2$
  Integer(Kind=wi), Parameter, Public :: MIX_FENDER_HASLEY = 2
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

  ! ZBL constants
  Real(wp), Parameter, Dimension(4) :: b = [0.18175_wp, 0.50986_wp, 0.28022_wp, 0.02817_wp], &
                                       c = [3.1998_wp, 0.94229_wp, 0.40290_wp, 0.20162_wp]
  Real(wp), Parameter, Public       :: ab = 0.52917721067_wp

  !> Type containing Van der Waals data
  Type, Public :: vdw_type
    Private

    !> No Van der Waals switch
    Logical, Public                       :: no_vdw
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
    Integer(Kind=wi), Allocatable, Public :: list(:)
    Integer(Kind=wi), Allocatable, Public :: ltp(:)
    !> VdW parameters
    Real(Kind=wp), Allocatable, Public    :: param(:, :)
    !> VdW cut off
    Real(Kind=wp), Public                 :: cutoff = 0.0_wp
    Real(Kind=wp), Allocatable, Public    :: sigeps(:, :)
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
    Final             :: cleanup
  End Type vdw_type

  Public :: vdw_forces, vdw_generate, vdw_table_read, vdw_lrc, vdw_direct_fs_generate

Contains

  Subroutine allocate_vdw_arrays(T)
    Class(vdw_type) :: T

    Integer, Dimension(1:4) :: fail

    fail = 0

    Allocate (T%list(1:T%max_vdw), Stat=fail(1))
    Allocate (T%ltp(1:T%max_vdw), Stat=fail(2))
    Allocate (T%param(1:T%max_param, 1:T%max_vdw), Stat=fail(3))
    Allocate (T%sigeps(1:2, 1:T%max_vdw), Stat=fail(4))

    If (Any(fail > 0)) Call error(1022)

    T%list = 0
    T%ltp = 0

    T%param = 0.0_wp
    T%sigeps = 0.0_wp
  End Subroutine allocate_vdw_arrays

  Subroutine allocate_vdw_table_arrays(T)
    Class(vdw_type) :: T

    Integer, Dimension(1:2) :: fail

    fail = 0

    Allocate (T%tab_potential(0:T%max_grid, 1:T%max_vdw), Stat=fail(1))
    Allocate (T%tab_force(0:T%max_grid, 1:T%max_vdw), Stat=fail(2))

    If (Any(fail > 0)) Call error(1063)

    T%tab_potential = 0.0_wp
    T%tab_force = 0.0_wp
  End Subroutine allocate_vdw_table_arrays

  Subroutine allocate_vdw_direct_fs_arrays(T)
    Class(vdw_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%afs(1:T%max_vdw), T%bfs(1:T%max_vdw), Stat=fail)

    If (fail > 0) Call error(1066)

    T%afs = 0.0_wp
    T%bfs = 0.0_wp
  End Subroutine allocate_vdw_direct_fs_arrays

  Pure Subroutine zbl(r, k, ia, phi, dphi)
    Real(wp), Intent(In   ) :: r, k, ia
    Real(wp), Intent(  Out) :: phi, dphi

    Integer  :: i
    Real(wp) :: ir, t1, x

    phi = 0.0_wp
    dphi = 0.0_wp
    x = r * ia
    ir = 1.0_wp / r
    Do i = 1, 4
      t1 = b(i) * Exp(-x * c(i))
      phi = phi + t1
      dphi = dphi - c(i) * t1
    End Do
    phi = k * phi * ir
    ! -rU/r
    dphi = phi - ia * k * dphi

  End Subroutine zbl

  Pure Subroutine fm(r, rm, ic, f, df)
    Real(wp), Intent(In   ) :: r, rm, ic
    Real(wp), Intent(  Out) :: f, df

    Real(wp) :: t

    If (r < rm) Then
      t = Exp(-(rm - r) * ic) * 0.5_wp
      f = 1.0_wp - t
      ! -rf/r
      df = r * ic * t
    Else
      t = Exp(-(r - rm) * ic) * 0.5_wp
      f = t
      ! -rf/r
      df = r * ic * t
    End If
  End Subroutine fm

  Pure Subroutine morse(r, d, k, r0, m, dm)
    Real(wp), Intent(In   ) :: r, d, k, r0
    Real(wp), Intent(  Out) :: m, dm

    Real(wp) :: t

    t = Exp(-k * (r - r0))

    m = d * ((1.0_wp - t)**2 - 1.0_wp)
    dm = -2.0_wp * r * d * k * (1.0_wp - t) * t

  End Subroutine morse

  Pure Subroutine buckingham(r, A, r0, C, b, db)
    Real(wp), Intent(In   ) :: r, A, r0, C
    Real(wp), Intent(  Out) :: b, db

    Real(wp) :: ir0, t1, t2

    ir0 = 1.0_wp / r0
    t1 = A * Exp(-r * ir0)
    t2 = -C * (1.0_wp / r)**6

    b = t1 + t2
    db = r * ir0 * t1 + 6.0_wp * t2
  End Subroutine buckingham

  Pure Subroutine zbls(r, kk, ia, rm, ic, d, k, r0, V, dV)
    Real(wp), Intent(In   ) :: r, kk, ia, rm, ic, d, k, r0
    Real(wp), Intent(  Out) :: V, dV

    Real(wp) :: df, dm, dz, f, m, z

    Call zbl(r, kk, ia, z, dz)
    Call fm(r, rm, ic, f, df)
    Call morse(r, d, k, r0, m, dm)
    V = f * z + (1.0_wp - f) * m
    dV = f * dz + df * z + (1.0_wp - f) * dm - df * m
  End Subroutine zbls

  Pure Subroutine zblb(r, kk, ia, rm, ic, A, r0, C, V, dV)
    Real(wp), Intent(In   ) :: r, kk, ia, rm, ic, A, r0, C
    Real(wp), Intent(  Out) :: V, dV

    Real(wp) :: b, db, df, dz, f, z

    Call zbl(r, kk, ia, z, dz)
    Call fm(r, rm, ic, f, df)
    Call buckingham(r, A, r0, C, b, db)
    V = f * z + (1.0_wp - f) * b
    dV = f * dz + df * z + (1.0_wp - f) * db - df * b
  End Subroutine zblb

  Pure Real(wp) Function intRadZBL(kk, a, rw, prec)
    Real(wp), Intent(In   ) :: kk, a, rw, prec

    Integer  :: i, j, n
    Real(wp) :: df0, f0, f1, f2, h, ie, is, s, sold, x

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
        Call zbl(x, kk, a, f0, df0)
        Call zbl(x + h, kk, a, f1, df0)
        Call zbl(x + 2.0_wp * h, kk, a, f2, df0)
        s = s + x * x * f0 + &
            4.0_wp * (x + h)**2 * f1 + &
            (x + 2.0_wp * h)**2 * f2
      End Do

      j = j + 1
    End Do
    intRadZBL = s * h / 3.0_wp
  End Function intRadZBL

  Pure Real(wp) Function intdRadZBL(kk, a, rw, prec)
    Real(wp), Intent(In   ) :: kk, a, rw, prec

    Integer  :: i, j, n
    Real(wp) :: df0, df1, df2, f0, h, ie, is, s, sold, x

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
        Call zbl(x, kk, a, f0, df0)
        Call zbl(x + h, kk, a, f0, df1)
        Call zbl(x + 2.0_wp * h, kk, a, f0, df2)
        s = s + x * x * df0 + &
            4.0_wp * (x + h)**2 * df1 + &
            (x + 2.0_wp * h)**2 * df2
      End Do

      j = j + 1
    End Do
    intdRadZBL = s * h / 3.0_wp
  End Function intdRadZBL

  Pure Subroutine LJ(r, eps, sig, e, v)
    Real(wp), Intent(In   ) :: r, eps, sig
    Real(wp), Intent(  Out) :: e, v

    Real(wp) :: sor6

    sor6 = (sig / r)**6
    e = 4.0_wp * eps * sor6 * (sor6 - 1.0_wp)
    v = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp)

  End Subroutine LJ

  Pure Subroutine LJ126(r, A, B, e, v)
    Real(wp), Intent(In   ) :: r, A, B
    Real(wp), Intent(  Out) :: e, v

    Real(wp) :: sor6

    sor6 = (1.0_wp / r)**6
    e = (A * sor6 - B) * sor6
    v = 6.0_wp * sor6 * (2.0_wp * a * sor6 - b)

  End Subroutine LJ126

  Pure Subroutine MDF(r, ri, rc, e, v)
    Real(wp), Intent(In   ) :: r, ri, rc
    Real(wp), Intent(  Out) :: e, v

    Real(wp) :: rci

    If (r < ri) Then
      e = 1.0_wp
      v = 0.0_wp
    Else If (r > rc) Then
      e = 0.0_wp
      v = 0.0_wp
    Else
      rci = (rc - ri)**5
      e = (rc - r)**3 * (10 * ri**2 - 5 * rc * ri - 15 * r * ri + rc**2 + 3 * r * rc + 6 * r**2) / rci
      v = 30.0_wp * r * (r - rc)**2 * (r - ri)**2 / rci
    End If

  End Subroutine MDF

  Pure Subroutine mlj(r, eps, sig, ri, rc, e, v)
    Real(wp), Intent(In   ) :: r, eps, sig, ri, rc
    Real(wp), Intent(  Out) :: e, v

    Real(wp) :: el, em, vl, vm

    Call LJ(r, eps, sig, el, vl)
    Call MDF(r, ri, rc, em, vm)
    e = el * em
    v = vl * em + vm * el

  End Subroutine mlj

  Pure Subroutine mbuck(r, A, r0, c, ri, rc, e, v)
    Real(wp), Intent(In   ) :: r, A, r0, c, ri, rc
    Real(wp), Intent(  Out) :: e, v

    Real(wp) :: eb, em, vb, vm

    Call buckingham(r, A, r0, c, eb, vb)
    Call MDF(r, ri, rc, em, vm)
    e = eb * em
    v = vb * em + vm * eb

  End Subroutine mbuck

  Pure Subroutine mlj126(r, A, B, ri, rc, e, v)
    Real(wp), Intent(In   ) :: r, A, B, ri, rc
    Real(wp), Intent(  Out) :: e, v

    Real(wp) :: el, em, vl, vm

    Call LJ126(r, A, B, el, vl)
    Call MDF(r, ri, rc, em, vm)
    e = el * em
    v = vl * em + vm * el
  End Subroutine mlj126

  Pure Real(wp) Function intRadMDF(pot, a, b, c, ri, rw, prec)
    Character(Len=*), Intent(In   ) :: pot
    Real(wp),         Intent(In   ) :: a, b, c, ri, rw, prec

    Integer  :: i, j, n
    Real(wp) :: df0, f0, f1, f2, h, ie, is, s, sold, x

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
        ! this is stupid by remember we do it only once
        If (pot == 'm126') Then
          Call mlj126(x, a, b, ri, rw, f0, df0)
          Call mlj126(x + h, a, b, ri, rw, f1, df0)
          Call mlj126(x + 2.0_wp * h, a, b, ri, rw, f2, df0)
        Else If (pot == 'mbuc') Then
          Call mbuck(x, a, b, c, ri, rw, f0, df0)
          Call mbuck(x + h, a, b, c, ri, rw, f1, df0)
          Call mbuck(x + 2.0_wp * h, a, b, c, ri, rw, f2, df0)
        Else If (pot == 'mlj') Then
          Call mlj(x, a, b, ri, rw, f0, df0)
          Call mlj(x + h, a, b, ri, rw, f1, df0)
          Call mlj(x + 2.0_wp * h, a, b, ri, rw, f2, df0)
        End If

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

    Integer  :: i, j, n
    Real(wp) :: df0, df1, df2, f0, h, ie, is, s, sold, x

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
        If (pot == 'm126') Then
          Call mlj126(x, a, b, ri, rw, f0, df0)
          Call mlj126(x + h, a, b, ri, rw, f0, df1)
          Call mlj126(x + 2.0_wp * h, a, b, ri, rw, f0, df2)
        Else If (pot == 'mbuc') Then
          Call mbuck(x, a, b, c, ri, rw, f0, df0)
          Call mbuck(x + h, a, b, c, ri, rw, f0, df1)
          Call mbuck(x + 2.0_wp * h, a, b, c, ri, rw, f0, df2)
        Else If (pot == 'mlj') Then
          Call mlj(x, a, b, ri, rw, f0, df0)
          Call mlj(x + h, a, b, ri, rw, f0, df1)
          Call mlj(x + 2.0_wp * h, a, b, ri, rw, f0, df2)
        End If
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

    Character(Len=256)                       :: message, messages(3)
    Integer                                  :: fail, i, ivdw, j, k, keypot
    Real(Kind=wp)                            :: a, b, c, d, denprd, e0, eadd, eps, kk, mr, nr, &
                                                padd, plrc, r, r0, s9, sig, t, z1, z2
    Real(Kind=wp), Allocatable, Dimension(:) :: numfrz

    fail = 0
    Allocate (numfrz(sites%mxatyp), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vdw_lrc allocation failure'
      Call error(0, message)
    End If

    ! initialise long-range corrections to energy and pressure

    plrc = 0.0_wp
    vdws%elrc = 0.0_wp

    If (vdws%l_force_shift) Go To 10 ! force-shifting

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
          If (keypot == VDW_TAB) Then

            ! tabulated energy and pressure lrc

            eadd = vdws%param(1, k)
            padd = -vdws%param(2, k)

          Else If (keypot == VDW_12_6) Then

            ! 12-6 potential :: u=a/r^12-b/r^6

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            r = vdws%cutoff

            eadd = a / (9.0_wp * r**9) - b / (3.0_wp * r**3)
            padd = 12.0_wp * a / (9.0_wp * r**9) - 6.0_wp * b / (3.0_wp * r**3)

          Else If (keypot == VDW_LENNARD_JONES) Then

            ! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)
            r = vdws%cutoff

            eadd = 4.0_wp * eps * (sig**12 / (9.0_wp * r**9) - sig**6 / (3.0_wp * r**3))
            padd = 8.0_wp * eps * (6.0_wp * sig**12 / (9.0_wp * r**9) - sig**6 / (r**3))

          Else If (keypot == VDW_N_M) Then

            ! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(r0/r)^m]

            e0 = vdws%param(1, k)
            nr = vdws%param(2, k)
            mr = vdws%param(3, k)
            r0 = vdws%param(4, k)
            r = vdws%cutoff

       eadd = e0 / (nr - mr) * (mr * r0**nr / ((nr - 3.0_wp) * r**(nr - 3.0_wp)) - nr * r0**mr / ((mr - 3.0_wp) * r**(mr - 3.0_wp)))
       padd = e0 / (nr - mr) * nr * mr * (r0**nr / ((nr - 3.0_wp) * r**(nr - 3.0_wp)) - r0**mr / ((mr - 3.0_wp) * r**(mr - 3.0_wp)))

          Else If (keypot == VDW_BUCKINGHAM) Then

            ! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

            c = vdws%param(3, k)
            r = vdws%cutoff

            eadd = -c / (3.0_wp * r**3)
            padd = -2.0_wp * c / (r**3)

          Else If (keypot == VDW_BORN_HUGGINS_MEYER) Then

            ! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

            c = vdws%param(4, k)
            d = vdws%param(5, k)
            r = vdws%cutoff

            eadd = -c / (3.0_wp * r**3) - d / (5.0_wp * r**5)
            padd = -2.0_wp * c / (r**3) - 8.0_wp * d / (5.0_wp * r**5)

          Else If (keypot == VDW_HYDROGEN_BOND) Then

            ! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            r = vdws%cutoff

            eadd = a / (9.0_wp * r**9) - b / (7.0_wp * r**7)
            padd = 12.0_wp * a / (9.0_wp * r**9) - 10.0_wp * b / (7.0_wp * r**7)

          Else If (keypot == VDW_MORSE) Then

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

          Else If (keypot == VDW_AMOEBA) Then

            ! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)

            a = 0.07_wp
            b = 0.12_wp
            e0 = 1.0e-12_wp

            eadd = intRadMM3(sig, a, b, eps, vdws%cutoff, e0)
            padd = -intRaddMM3(sig, a, b, eps, vdws%cutoff, e0)

          Else If (keypot == VDW_LENNARD_JONES_COHESIVE) Then

            ! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)
            c = vdws%param(3, k)
            r = vdws%cutoff

            eadd = 4.0_wp * eps * (sig**12 / (9.0_wp * r**9) - c * sig**6 / (3.0_wp * r**3))
            padd = 8.0_wp * eps * (6.0_wp * sig**12 / (9.0_wp * r**9) - c * sig**6 / (r**3))

          Else If (keypot == VDW_MORSE_12) Then
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

          Else If (keypot == VDW_RYDBERG) Then

            ! Rydberg potential:: u=(a+b*r)Exp(-r/c)

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            c = vdws%param(3, k)
            t = Exp(-vdws%cutoff / c)

            eadd = (b * c * vdws%cutoff**3 + (3 * b * c**2 + a * c) * vdws%cutoff**2 + (6 * b * c**3 + 2 * a * c**2) * vdws%cutoff &
                    + 6 * b * c**4 + 2 * a * c**3) * t
            padd = (b * vdws%cutoff**4 + (3 * b * c + a) * vdws%cutoff**3 + (9 * b * c**2 + 3 * a * c) * vdws%cutoff**2 + &
                    (18 * b * c**3 + 6 * a * c**2) * vdws%cutoff + 18 * b * c**4 + 6 * a * c**3) * t

          Else If (keypot == VDW_ZBL) Then

            ! ZBL potential:: u=Z1Z2/(40r)_{i=1}^4b_ie^{-c_i*r/a}

            z1 = vdws%param(1, k)
            z2 = vdws%param(2, k)

            ! this is in fact inverse a
            a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
            kk = z1 * z2 * r4pie0
            eadd = intRadZBL(kk, a, vdws%cutoff, 1e-12_wp)
            padd = intdRadZBL(kk, a, vdws%cutoff, 1e-12_wp)

          Else If (keypot == VDW_ZBL_SWITCH_MORSE) Then

            ! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

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

          Else If (keypot == VDW_ZBL_SWITCH_BUCKINGHAM) Then

            ! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

            A = vdws%param(5, k)
            r0 = vdws%param(6, k)
            c = vdws%param(7, k)

            t = A * Exp(-vdws%cutoff / r0)

            eadd = (vdws%cutoff**2 + 2 * r0 * vdws%cutoff + 2 * r0**2) * t * r0 - c / (3.0_wp * vdws%cutoff**3)
         padd = (vdws%cutoff**3 + 3 * r0 * vdws%cutoff**2 + 6 * r0**2 * vdws%cutoff + 6 * r0**3) * t - 2.0_wp * c / (vdws%cutoff**3)

          Else If (keypot == VDW_LJ_MDF) Then

            ! LJ tappered with MDF:: u=f(r)LJ(r)

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            c = vdws%param(3, k)
            eadd = intRadMDF("mlj", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)
            padd = intdRadMDF("mlj", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)

          Else If (keypot == VDW_BUCKINGHAM_MDF) Then

            ! Buckingham tappered with MDF:: u=f(r)Buck(r)

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            c = vdws%param(3, k)
            r0 = vdws%param(3, k)
            eadd = intRadMDF("mbuc", a, b, c, r0, vdws%cutoff, 1e-12_wp)
            padd = intdRadMDF("mbuc", a, b, c, r0, vdws%cutoff, 1e-12_wp)

          Else If (keypot == VDW_126_MDF) Then

            ! LJ tappered with MDF:: u=f(r)LJ12-6(r)

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            c = vdws%param(3, k)
            eadd = intRadMDF("m126", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)
            padd = intdRadMDF("m126", a, b, 0.0_wp, c, vdws%cutoff, 1e-12_wp)

          End If

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

    10 Continue

    Write (messages(1), '(a)') 'long-range correction for:'
    Write (messages(2), '(2x,a,e15.6)') 'vdw energy ', vdws%elrc / engunit
    Write (messages(3), '(2x,a,e15.6)') 'vdw pressure ', plrc * prsunt
    Call info(messages, 3, .true.)

    ! convert plrc to a viral term

    vdws%vlrc = plrc * (-3.0_wp * config%volm)

    Deallocate (numfrz, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vdw_lrc deallocation failure'
      Call error(0, message)
    End If

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

    Integer       :: ivdw, keypot
    Real(Kind=wp) :: a, b, c, d, dz, e0, eps, ic, k, kk, mr, nr, r0, r0rm, r0rn, r_6, rho, ri, rm, &
                     sig, sor6, t1, t2, t3, z, z1, z2

    ! allocate arrays for force-shifted corrections

    Call vdws%init_direct()

    ! construct arrays for all types of vdw potential

    Do ivdw = 1, vdws%n_vdw

      keypot = vdws%ltp(ivdw)
      If (keypot == VDW_12_6) Then

        ! 12-6 potential :: u=a/r^12-b/r^6

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)

        r_6 = vdws%cutoff**(-6)

        vdws%afs(ivdw) = 6.0_wp * r_6 * (2.0_wp * a * r_6 - b)
        vdws%bfs(ivdw) = -r_6 * (a * r_6 - b) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_LENNARD_JONES) Then

        ! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)

        sor6 = (sig / vdws%cutoff)**6

        vdws%afs(ivdw) = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp)
        vdws%bfs(ivdw) = -4.0_wp * eps * sor6 * (sor6 - 1.0_wp) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_N_M) Then

        ! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0 = vdws%param(1, ivdw)
        nr = vdws%param(2, ivdw)
        mr = vdws%param(3, ivdw)
        r0 = vdws%param(4, ivdw)

        a = r0 / vdws%cutoff
        b = 1.0_wp / Real(nr - mr, wp)
        r0rn = a**nr
        r0rm = a**mr

        vdws%afs(ivdw) = e0 * mr * nr * (r0rn - r0rm) * b
        vdws%bfs(ivdw) = -e0 * (mr * r0rn - nr * r0rm) * b - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_BUCKINGHAM) Then

        ! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a = vdws%param(1, ivdw)
        rho = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)

        If (Abs(rho) <= zero_plus) Then
          If (Abs(a) <= zero_plus) Then
            rho = 1.0_wp
          Else
            Call error(467)
          End If
        End If

        b = vdws%cutoff / rho
        t1 = a * Exp(-b)
        t2 = -c / vdws%cutoff**6

        vdws%afs(ivdw) = (t1 * b + 6.0_wp * t2)
        vdws%bfs(ivdw) = -(t1 + t2) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_BORN_HUGGINS_MEYER) Then

        ! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)
        sig = vdws%param(3, ivdw)
        c = vdws%param(4, ivdw)
        d = vdws%param(5, ivdw)

        t1 = a * Exp(b * (sig - vdws%cutoff))
        t2 = -c / vdws%cutoff**6
        t3 = -d / vdws%cutoff**8

        vdws%afs(ivdw) = (t1 * vdws%cutoff * b + 6.0_wp * t2 + 8.0_wp * t3)
        vdws%bfs(ivdw) = -(t1 + t2 + t3) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_HYDROGEN_BOND) Then

        ! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)

        t1 = a / vdws%cutoff**12
        t2 = -b / vdws%cutoff**10

        vdws%afs(ivdw) = (12.0_wp * t1 + 10.0_wp * t2)
        vdws%bfs(ivdw) = -(t1 + t2) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_N_M_SHIFT) Then

        ! shifted and force corrected n-m potential (w.smith) ::

      Else If (keypot == VDW_MORSE) Then

        ! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0 = vdws%param(1, ivdw)
        r0 = vdws%param(2, ivdw)
        kk = vdws%param(3, ivdw)

        t1 = Exp(-kk * (vdws%cutoff - r0))

        vdws%afs(ivdw) = -2.0_wp * e0 * kk * t1 * (1.0_wp - t1) * vdws%cutoff
        vdws%bfs(ivdw) = -e0 * t1 * (t1 - 2.0_wp) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_WCA) Then

        ! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
        ! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)
        d = vdws%param(3, ivdw)

        sor6 = (sig / (vdws%cutoff - d))**6

        vdws%afs(ivdw) = (24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp) / (vdws%cutoff - d)) * vdws%cutoff
        vdws%bfs(ivdw) = -(4.0_wp * eps * sor6 * (sor6 - 1.0_wp) + eps) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_DPD) Then ! all zeroed in vdw

        ! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

        !       vdws%afs(ivdw) = 0.0_wp !initialised in vdw
        !       vdws%bfs(ivdw) = 0.0_wp !initialised in vdw

      Else If (keypot == VDW_AMOEBA) Then

        ! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)

        rho = sig / vdws%cutoff
        t1 = 1.0_wp / (0.07_wp + rho)
        t2 = 1.0_wp / (0.12_wp + rho**7)
        t3 = eps * (1.07_wp / t1**7)

        vdws%afs(ivdw) = -7.0_wp * t3 * rho * (((1.12_wp / t2) - 2.0_wp) / t1 + (1.12_wp / t2**2) * rho**6)
        vdws%bfs(ivdw) = -t3 * ((1.12_wp / t2) - 2.0_wp) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_LENNARD_JONES_COHESIVE) Then

        ! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)

        sor6 = (sig / vdws%cutoff)**6

        vdws%afs(ivdw) = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - c)
        vdws%bfs(ivdw) = -4.0_wp * eps * sor6 * (sor6 - c) - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_MORSE_12) Then

        ! Morse potential with twelve term:: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0 = vdws%param(1, ivdw)
        r0 = vdws%param(2, ivdw)
        kk = vdws%param(3, ivdw)
        c = vdws%param(4, ivdw)

        t1 = Exp(-kk * (vdws%cutoff - r0))
        sor6 = c / vdws%cutoff**12

        vdws%afs(ivdw) = -2.0_wp * e0 * kk * t1 * (1.0_wp - t1) * vdws%cutoff + 12.0_wp * sor6
        vdws%bfs(ivdw) = -e0 * t1 * (t1 - 2.0_wp) + sor6 - vdws%afs(ivdw)
        vdws%afs(ivdw) = vdws%afs(ivdw) / vdws%cutoff

      Else If (keypot == VDW_RYDBERG) Then

        ! Morse potential with twelve term:: u=(a+b*r)Exp(-r/c)

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)

        kk = 1.0_wp / c
        t1 = Exp(-vdws%cutoff * kk)
        vdws%afs(ivdw) = (a + b * vdws%cutoff) * kk * t1 - b * t1
        vdws%bfs(ivdw) = -(a * c + a * vdws%cutoff + b * vdws%cutoff * vdws%cutoff) * kk * t1

      Else If (keypot == VDW_ZBL) Then

        ! ZBL potential:: u=Z1Z2/(40r)_{i=1}^4b_ie^{-c_i*r/a}

        z1 = vdws%param(1, ivdw)
        z2 = vdws%param(2, ivdw)

        a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
        kk = z1 * z2 * r4pie0

        Call zbl(vdws%cutoff, kk, a, z, dz)
        vdws%afs(ivdw) = dz / vdws%cutoff
        vdws%bfs(ivdw) = -z - dz

      Else If (keypot == VDW_ZBL_SWITCH_MORSE) Then

        ! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = vdws%param(1, ivdw)
        z2 = vdws%param(2, ivdw)
        rm = vdws%param(3, ivdw)
        ic = 1.0_wp / vdws%param(4, ivdw)
        e0 = vdws%param(5, ivdw)
        r0 = vdws%param(6, ivdw)
        k = vdws%param(7, ivdw)

        a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
        kk = z1 * z2 * r4pie0
        Call zbls(vdws%cutoff, kk, a, rm, ic, e0, k, r0, z, dz)
        vdws%afs(ivdw) = dz / vdws%cutoff
        vdws%bfs(ivdw) = -z - dz

      Else If (keypot == VDW_ZBL_SWITCH_BUCKINGHAM) Then

        ! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = vdws%param(1, ivdw)
        z2 = vdws%param(2, ivdw)
        rm = vdws%param(3, ivdw)
        ic = 1.0_wp / vdws%param(4, ivdw)
        e0 = vdws%param(5, ivdw)
        r0 = vdws%param(6, ivdw)
        k = vdws%param(7, ivdw)

        a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
        kk = z1 * z2 * r4pie0
        Call zblb(vdws%cutoff, kk, a, rm, ic, e0, r0, k, z, dz)
        vdws%afs(ivdw) = dz / vdws%cutoff
        vdws%bfs(ivdw) = -z - dz

      Else If (keypot == VDW_LJ_MDF) Then

        ! LJ tappered with MDF:: u=f(r)LJ(r)

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)
        ri = vdws%param(3, ivdw)

        Call mlj(vdws%cutoff, eps, sig, ri, vdws%cutoff, z, dz)
        vdws%afs(ivdw) = dz / vdws%cutoff
        vdws%bfs(ivdw) = -z - dz

      Else If (keypot == VDW_BUCKINGHAM_MDF) Then

        ! Buckingham tappered with MDF:: u=f(r)Buck(r)
        a = vdws%param(1, ivdw)
        rho = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)
        ri = vdws%param(4, ivdw)
        Call mbuck(vdws%cutoff, A, rho, c, ri, vdws%cutoff, z, dz)

        vdws%afs(ivdw) = dz / vdws%cutoff
        vdws%bfs(ivdw) = -z - dz

      Else If (keypot == VDW_126_MDF) Then

        ! LJ tappered with MDF:: u=f(r)LJ(r)
        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)
        ri = vdws%param(3, ivdw)
        Call mlj126(vdws%cutoff, A, B, ri, vdws%cutoff, z, dz)
        vdws%afs(ivdw) = dz / vdws%cutoff
        vdws%bfs(ivdw) = -z - dz

      Else

        Call error(150)

      End If

    End Do

  End Subroutine vdw_direct_fs_generate

  Subroutine vdw_table_read(vdws, sites, comm)

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
    Type(comms_type), Intent(InOut) :: comm

    Character(Len=200)                       :: record
    Character(Len=256)                       :: message, messages(4)
    Character(Len=40)                        :: word
    Character(Len=8)                         :: atom1, atom2
    Integer                                  :: fail, i, ivdw, j, jtpatm, katom1, katom2, keyvdw, &
                                                l, ngrid
    Logical                                  :: remake, safe
    Real(Kind=wp)                            :: cutpot, delpot, dlrpot, ppp, rdr, rrr, t, t1, t2, &
                                                vk, vk1, vk2
    Real(Kind=wp), Allocatable, Dimension(:) :: buffer

    If (comm%idnode == 0) Open (Unit=ntable, File='TABLE')

    ! skip header record

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Go To 100

    ! read mesh resolution

    Call get_line(safe, ntable, record, comm)
    If (.not. safe) Go To 100

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
      Call error(48)
    End If

    If (cutpot < vdws%cutoff) Call error(504)

    fail = 0
    Allocate (buffer(0:ngrid), Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vdw_table_read allocation failure'
      Call error(0, message)
    End If

    ! read potential arrays for all pairs

    Do ivdw = 1, vdws%n_vdw

      ! read potential arrays if potential not already defined

      If (vdws%ltp(ivdw) == VDW_TAB) Then

        ! read pair potential labels and long-range corrections

        Call get_line(safe, ntable, record, comm)
        If (.not. safe) Go To 100

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
          Write (message, '(a,i0,a,i0,a)') '****', atom1, '***', atom2, '**** entry in TABLE'
          Call error(81, message, .true.)
        End If

        keyvdw = (Max(katom1, katom2) * (Max(katom1, katom2) - 1)) / 2 + Min(katom1, katom2)

        ! Only one vdw potential per pair is allowed
        ! (FIELD AND TABLE potentials overlapping)

        If (vdws%list(keyvdw) /= ivdw) Call error(23)

        ! read in potential arrays

        Do i = 1, (ngrid + 3) / 4
          j = Min(4, ngrid - (i - 1) * 4)
          If (comm%idnode == 0) Then
            Read (Unit=ntable, Fmt=*, End=100) buffer((i - 1) * 4 + 1:(i - 1) * 4 + j)
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
            Read (Unit=ntable, Fmt=*, End=100) buffer((i - 1) * 4 + 1:(i - 1) * 4 + j)
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
    If (comm%idnode == 0) Then
      Close (Unit=ntable)
    End If

    ! convert to internal units

    Do ivdw = 1, vdws%n_vdw
      If (vdws%ltp(ivdw) == VDW_TAB) Then

        ! Sigma-epsilon initialisation

        vdws%sigeps(1, ivdw) = -1.0_wp
        vdws%sigeps(2, ivdw) = 0.0_wp

        Do i = 0, vdws%max_grid
          vdws%tab_potential(i, ivdw) = vdws%tab_potential(i, ivdw) * engunit
          vdws%tab_force(i, ivdw) = vdws%tab_force(i, ivdw) * engunit

          ! Sigma-epsilon search

          If ((.not. vdws%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
            If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
              If (Nint(Sign(1.0_wp, vdws%tab_potential(i - 1, ivdw))) == -Nint(Sign(1.0_wp, vdws%tab_potential(i, ivdw)))) &
                vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
            Else ! find epsilon
              If ((vdws%tab_potential(i - 2, ivdw) >= vdws%tab_potential(i - 1, ivdw) .and. &
                   vdws%tab_potential(i - 1, ivdw) <= vdws%tab_potential(i, ivdw)) .and. &
                  (nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i - 1, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 1, ivdw), vdws%tab_potential(i, ivdw)))) Then
                vdws%sigeps(2, ivdw) = -vdws%tab_potential(i - 1, ivdw)
              End If
            End If
          End If
        End Do
      End If
    End Do

    If (vdws%l_force_shift) Then
      Do ivdw = 1, vdws%n_vdw
        If (vdws%ltp(ivdw) == VDW_TAB) Then

          ! Sigma-epsilon initialisation

          vdws%sigeps(1, ivdw) = -1.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp

          ! Sigma-epsilon search

          Do i = 1, vdws%max_grid - 4
            If (i > 20) Then ! Assumes some safety against numeric black holes!!!
              t = vdws%tab_potential(i, ivdw) + vdws%tab_force(vdws%max_grid - 4, ivdw) * &
                  (Real(i, wp) * dlrpot / vdws%cutoff - 1.0_wp) - vdws%tab_potential(vdws%max_grid - 4, ivdw)
              t1 = vdws%tab_potential(i - 1, ivdw) + vdws%tab_force(vdws%max_grid - 4, ivdw) * &
                   (Real(i - 1, wp) * dlrpot / vdws%cutoff - 1.0_wp) - vdws%tab_potential(vdws%max_grid - 4, ivdw)
              If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
                If (Nint(Sign(1.0_wp, t1)) == -Nint(Sign(1.0_wp, t))) &
                  vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
              Else ! find epsilon
                t2 = vdws%tab_potential(i - 2, ivdw) + vdws%tab_force(vdws%max_grid - 4, ivdw) * &
                     (Real(i - 2, wp) * dlrpot / vdws%cutoff - 1.0_wp) - vdws%tab_potential(vdws%max_grid - 4, ivdw)

                If ((t2 >= t1 .and. t1 <= t) .and. &
                    (nequal(t1, t1) .or. nequal(t2, t) .or. nequal(t1, t))) Then
                  vdws%sigeps(2, ivdw) = -t1
                End If
              End If
            End If
          End Do
          vdws%tab_potential(vdws%max_grid - 3, ivdw) = 0.0_wp; vdws%tab_potential(vdws%max_grid - 2, ivdw) = 0.0_wp
          vdws%tab_force(vdws%max_grid - 3, ivdw) = 0.0_wp; vdws%tab_force(vdws%max_grid - 2, ivdw) = 0.0_wp
        End If
      End Do
    End If

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write (message, '(a)') 'vdw_table_read deallocation failure'
      Call error(0, message)
    End If

    Return

    ! end of file error exit

    100 Continue

    If (comm%idnode == 0) Close (Unit=ntable)
    Call error(24)

  End Subroutine vdw_table_read

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
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(vdw_type), Intent(InOut) :: vdws

    Integer       :: i, ivdw, keypot, m, n
    Real(Kind=wp) :: a, alpha, b, beta, c, d, dlrpot, dphi, e0, eps, k, kk, mr, nr, phi, r, r0, &
                     r0rm, r0rn, r_6, rc, rho, ri, rm, sig, sor6, t, t1, t2, t3, z1, z2

    ! allocate arrays for tabulating

    Call vdws%init_table()

    ! define grid resolution for potential arrays

    dlrpot = vdws%cutoff / Real(vdws%max_grid - 4, wp)

    ! construct arrays for all types of vdw potential

    Do ivdw = 1, vdws%n_vdw

      keypot = vdws%ltp(ivdw)
      If (keypot == VDW_12_6) Then

        ! 12-6 potential :: u=a/r^12-b/r^6

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          r_6 = r**(-6)

          vdws%tab_potential(i, ivdw) = r_6 * (a * r_6 - b)
          vdws%tab_force(i, ivdw) = 6.0_wp * r_6 * (2.0_wp * a * r_6 - b)
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          If (a * b > zero_plus) Then
            vdws%sigeps(1, ivdw) = (a / b)**(1.0_wp / 6.0_wp)
            vdws%sigeps(2, ivdw) = b**2 / (4.0_wp * a)
          End If ! else leave undetermined
        End If

      Else If (keypot == VDW_LENNARD_JONES) Then

        ! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          sor6 = (sig / r)**6

          vdws%tab_potential(i, ivdw) = 4.0_wp * eps * sor6 * (sor6 - 1.0_wp)
          vdws%tab_force(i, ivdw) = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp)
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = sig
          vdws%sigeps(2, ivdw) = eps
        End If

      Else If (keypot == VDW_N_M) Then

        ! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0 = vdws%param(1, ivdw)
        nr = vdws%param(2, ivdw)
        mr = vdws%param(3, ivdw)
        r0 = vdws%param(4, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          a = r0 / r
          b = 1.0_wp / (nr - mr)
          r0rn = a**nr
          r0rm = a**mr

          vdws%tab_potential(i, ivdw) = e0 * (mr * r0rn - nr * r0rm) * b
          vdws%tab_force(i, ivdw) = e0 * mr * nr * (r0rn - r0rm) * b
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = r0 * (nr / mr)**(1.0_wp / (mr - nr))
          vdws%sigeps(2, ivdw) = e0
        End If

      Else If (keypot == VDW_BUCKINGHAM) Then

        ! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a = vdws%param(1, ivdw)
        rho = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)

        If (Abs(rho) <= zero_plus) Then
          If (Abs(a) <= zero_plus) Then
            rho = 1.0_wp
          Else
            Call error(467)
          End If
        End If

        ! Sigma-epsilon initialisation

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = -1.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          b = r / rho
          t1 = a * Exp(-b)
          t2 = -c / r**6

          vdws%tab_potential(i, ivdw) = t1 + t2
          vdws%tab_force(i, ivdw) = t1 * b + 6.0_wp * t2

          ! Sigma-epsilon search

          If ((.not. vdws%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
            If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
              If (Nint(Sign(1.0_wp, vdws%tab_potential(i - 1, ivdw))) == -Nint(Sign(1.0_wp, vdws%tab_potential(i, ivdw)))) &
                vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
            Else ! find epsilon
              If ((vdws%tab_potential(i - 2, ivdw) >= vdws%tab_potential(i - 1, ivdw) .and. &
                   vdws%tab_potential(i - 1, ivdw) <= vdws%tab_potential(i, ivdw)) .and. &
                  (nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i - 1, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 1, ivdw), vdws%tab_potential(i, ivdw)))) Then
                vdws%sigeps(2, ivdw) = -vdws%tab_potential(i - 1, ivdw)
              End If
            End If
          End If
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Else If (keypot == VDW_BORN_HUGGINS_MEYER) Then

        ! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)
        sig = vdws%param(3, ivdw)
        c = vdws%param(4, ivdw)
        d = vdws%param(5, ivdw)

        ! Sigma-epsilon initialisation

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = -1.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          t1 = a * Exp(b * (sig - r))
          t2 = -c / r**6
          t3 = -d / r**8

          vdws%tab_potential(i, ivdw) = t1 + t2 + t3
          vdws%tab_force(i, ivdw) = t1 * r * b + 6.0_wp * t2 + 8.0_wp * t3

          ! Sigma-epsilon search

          If ((.not. vdws%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
            If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
              If (Nint(Sign(1.0_wp, vdws%tab_potential(i - 1, ivdw))) == -Nint(Sign(1.0_wp, vdws%tab_potential(i, ivdw)))) &
                vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
            Else ! find epsilon
              If ((vdws%tab_potential(i - 2, ivdw) >= vdws%tab_potential(i - 1, ivdw) .and. &
                   vdws%tab_potential(i - 1, ivdw) <= vdws%tab_potential(i, ivdw)) .and. &
                  (nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i - 1, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 1, ivdw), vdws%tab_potential(i, ivdw)))) Then
                vdws%sigeps(2, ivdw) = -vdws%tab_potential(i - 1, ivdw)
              End If
            End If
          End If
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Else If (keypot == VDW_HYDROGEN_BOND) Then

        ! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          t1 = a / r**12
          t2 = -b / r**10

          vdws%tab_potential(i, ivdw) = t1 + t2
          vdws%tab_force(i, ivdw) = 12.0_wp * t1 + 10.0_wp * t2
        End Do

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = Sqrt(a / b)
          vdws%sigeps(2, ivdw) = ((b / 6.0_wp)**6) * ((5.0_wp / a)**5)
        End If
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Else If (keypot == VDW_N_M_SHIFT) Then

        ! shifted and force corrected n-m potential (w.smith) ::

        e0 = vdws%param(1, ivdw)
        n = Nint(vdws%param(2, ivdw)); nr = Real(n, wp)
        m = Nint(vdws%param(3, ivdw)); mr = Real(m, wp)
        r0 = vdws%param(4, ivdw)
        rc = vdws%param(5, ivdw); If (rc < 1.0e-6_wp) rc = vdws%cutoff

        If (n <= m) Call error(470)

        ! Sigma-epsilon initialisation

        vdws%sigeps(1, ivdw) = -1.0_wp
        vdws%sigeps(2, ivdw) = 0.0_wp

        t = Real(n - m, wp)

        b = 1.0_wp / t
        c = rc / r0; If (c < 1.0_wp) Call error(468)

        beta = c * ((c**(m + 1) - 1.0_wp) / (c**(n + 1) - 1.0_wp))**b
        alpha = -t / (mr * (beta**n) * (1.0_wp + (nr / c - nr - 1.0_wp) / c**n) &
                      - nr * (beta**m) * (1.0_wp + (mr / c - mr - 1.0_wp) / c**m))
        e0 = e0 * alpha

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot
          If (r <= rc) Then
            a = r0 / r

            vdws%tab_potential(i, ivdw) = e0 * (mr * (beta**n) * (a**n - (1.0_wp / c)**n) &
                                                - nr * (beta**m) * (a**m - (1.0_wp / c)**m) &
                                                + nr * mr * ((r / rc - 1.0_wp) * ((beta / c)**n - (beta / c)**m))) * b
            vdws%tab_force(i, ivdw) = e0 * mr * nr * ((beta**n) * a**n - (beta**m) * a**m &
                                                      - r / rc * ((beta / c)**n - (beta / c)**m)) * b

            ! Sigma-epsilon search

            If (i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
                If (Nint(Sign(1.0_wp, vdws%tab_potential(i - 1, ivdw))) == -Nint(Sign(1.0_wp, vdws%tab_potential(i, ivdw)))) &
                  vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
              Else ! find epsilon
                If ((vdws%tab_potential(i - 2, ivdw) >= vdws%tab_potential(i - 1, ivdw) .and. &
                     vdws%tab_potential(i - 1, ivdw) <= vdws%tab_potential(i, ivdw)) .and. &
                    (nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i - 1, ivdw)) .or. &
                     nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i, ivdw)) .or. &
                     nequal(vdws%tab_potential(i - 1, ivdw), vdws%tab_potential(i, ivdw)))) Then
                  vdws%sigeps(2, ivdw) = -vdws%tab_potential(i - 1, ivdw)
                End If
              End If
            End If
          End If ! The else condition is satisfied by the vdw initialisation
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Else If (keypot == VDW_MORSE) Then

        ! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0 = vdws%param(1, ivdw)
        r0 = vdws%param(2, ivdw)
        kk = vdws%param(3, ivdw)

        Do i = 0, vdws%max_grid
          r = Real(i, wp) * dlrpot

          t1 = Exp(-kk * (r - r0))

          vdws%tab_potential(i, ivdw) = e0 * ((1.0_wp - t1)**2 - 1.0_wp)
          vdws%tab_force(i, ivdw) = -2.0_wp * r * e0 * kk * (1.0_wp - t1) * t1
        End Do
        t1 = Exp(+kk * r0)
        vdws%tab_force(0, ivdw) = -2.0_wp * e0 * kk * (1.0_wp - t1) * t1

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = r0 - Log(2.0_wp) / kk
          vdws%sigeps(2, ivdw) = e0
        End If

      Else If (keypot == VDW_WCA) Then

        ! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
        ! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)
        d = vdws%param(3, ivdw)

        ! Sigma-epsilon initialisation

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = -1.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          If (r < vdws%param(4, ivdw) .or. Abs(r - d) < 1.0e-10_wp) Then ! Else leave them zeros
            sor6 = (sig / (r - d))**6

            vdws%tab_potential(i, ivdw) = 4.0_wp * eps * sor6 * (sor6 - 1.0_wp) + eps
            vdws%tab_force(i, ivdw) = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp) * r / (r - d)

            ! Sigma-epsilon search

            If ((.not. vdws%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
                If (Nint(Sign(1.0_wp, vdws%tab_potential(i - 1, ivdw))) == -Nint(Sign(1.0_wp, vdws%tab_potential(i, ivdw)))) &
                  vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
              Else ! find epsilon
                If ((vdws%tab_potential(i - 2, ivdw) >= vdws%tab_potential(i - 1, ivdw) .and. &
                     vdws%tab_potential(i - 1, ivdw) <= vdws%tab_potential(i, ivdw)) .and. &
                    (nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i - 1, ivdw)) .or. &
                     nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i, ivdw)) .or. &
                     nequal(vdws%tab_potential(i - 1, ivdw), vdws%tab_potential(i, ivdw)))) Then
                  vdws%sigeps(2, ivdw) = -vdws%tab_potential(i - 1, ivdw)
                End If
              End If
            End If
          End If
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Else If (keypot == VDW_DPD) Then

        ! DPD potential - Groot-Warren (standard) :: u=(1/2).a.rc.(1-r/rc)^2

        a = vdws%param(1, ivdw)
        rc = vdws%param(2, ivdw)

        Do i = 0, vdws%max_grid
          r = Real(i, wp) * dlrpot

          If (r < rc) Then
            t1 = 0.5_wp * a * rc
            t2 = 1.0_wp - r / rc

            vdws%tab_potential(i, ivdw) = t1 * t2**2
            vdws%tab_force(i, ivdw) = a * t2 * r
          End If
        End Do
        vdws%tab_force(0, ivdw) = a

        vdws%sigeps(1, ivdw) = rc
        vdws%sigeps(2, ivdw) = a

      Else If (keypot == VDW_AMOEBA) Then

        ! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          rho = r / sig
          t1 = 1.0_wp / (0.07_wp + rho)
          t2 = 1.0_wp / (0.12_wp + rho**7)
          t3 = eps * (1.07_wp * t1)**7

          t = t3 * ((1.12_wp * t2) - 2.0_wp)

          vdws%tab_potential(i, ivdw) = t
          vdws%tab_force(i, ivdw) = 7.0_wp * (t1 * t + 1.12_wp * t3 * t2**2 * rho**6) * rho
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = sig
          vdws%sigeps(2, ivdw) = eps
        End If

      Else If (keypot == VDW_LENNARD_JONES_COHESIVE) Then

        ! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          sor6 = (sig / r)**6

          vdws%tab_potential(i, ivdw) = 4.0_wp * eps * sor6 * (sor6 - c)
          vdws%tab_force(i, ivdw) = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - c)
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = sig
          vdws%sigeps(2, ivdw) = eps
        End If

      Else If (keypot == VDW_MORSE_12) Then

        ! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0 = vdws%param(1, ivdw)
        r0 = vdws%param(2, ivdw)
        kk = vdws%param(3, ivdw)
        c = vdws%param(4, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          t1 = Exp(-kk * (r - r0))
          sor6 = c / r**12

          vdws%tab_potential(i, ivdw) = e0 * ((1.0_wp - t1)**2 - 1.0_wp) + sor6
          vdws%tab_force(i, ivdw) = -2.0_wp * r * e0 * kk * (1.0_wp - t1) * t1 + 12.0_wp * sor6
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then !???
          vdws%sigeps(1, ivdw) = r0 - Log(2.0_wp) / kk
          vdws%sigeps(2, ivdw) = e0
        End If

      Else If (keypot == VDW_RYDBERG) Then

        ! Rydberg potential:: u=(a+b*r)Exp(-r/c)

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          kk = r / c
          t1 = Exp(-kk)

          vdws%tab_potential(i, ivdw) = (a + b * r) * t1
          vdws%tab_force(i, ivdw) = t1 * kk * (a - b * c + b * r)
        End Do
        vdws%tab_potential(0, ivdw) = a
        vdws%tab_force(0, ivdw) = 0

        If (.not. vdws%l_force_shift) Then !???
          vdws%sigeps(1, ivdw) = 1.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

      Else If (keypot == VDW_ZBL) Then

        ! ZBL potential:: u=Z1Z2/(40r)_{i=1}^4b_ie^{-c_i*r/a}

        z1 = vdws%param(1, ivdw)
        z2 = vdws%param(2, ivdw)

        ! this is in fact inverse a
        a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
        kk = z1 * z2 * r4pie0

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          Call zbl(r, kk, a, phi, dphi)

          vdws%tab_potential(i, ivdw) = phi
          vdws%tab_force(i, ivdw) = dphi

        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = 0.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

      Else If (keypot == VDW_ZBL_SWITCH_MORSE) Then

        ! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = vdws%param(1, ivdw)
        z2 = vdws%param(2, ivdw)
        rm = vdws%param(3, ivdw)
        c = 1.0_wp / vdws%param(4, ivdw)
        e0 = vdws%param(5, ivdw)
        r0 = vdws%param(6, ivdw)
        k = vdws%param(7, ivdw)

        a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
        kk = z1 * z2 * r4pie0

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          Call zbls(r, kk, a, rm, c, e0, k, r0, phi, dphi)
          vdws%tab_potential(i, ivdw) = phi
          vdws%tab_force(i, ivdw) = dphi
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = 0.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

      Else If (keypot == VDW_ZBL_SWITCH_BUCKINGHAM) Then

        ! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = vdws%param(1, ivdw)
        z2 = vdws%param(2, ivdw)
        rm = vdws%param(3, ivdw)
        c = 1.0_wp / vdws%param(4, ivdw)
        e0 = vdws%param(5, ivdw)
        r0 = vdws%param(6, ivdw)
        k = vdws%param(7, ivdw)

        a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
        kk = z1 * z2 * r4pie0

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          Call zblb(r, kk, a, rm, c, e0, r0, k, phi, dphi)
          vdws%tab_potential(i, ivdw) = phi
          vdws%tab_force(i, ivdw) = dphi
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = 0.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

      Else If (keypot == VDW_LJ_MDF) Then

        ! LJ tappered with MDF:: u=f(r)LJ(r)

        eps = vdws%param(1, ivdw)
        sig = vdws%param(2, ivdw)
        ri = vdws%param(3, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          Call mlj(r, eps, sig, ri, vdws%cutoff, phi, dphi)
          vdws%tab_potential(i, ivdw) = phi
          vdws%tab_force(i, ivdw) = dphi
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = sig
          vdws%sigeps(2, ivdw) = eps
        End If

      Else If (keypot == VDW_BUCKINGHAM_MDF) Then

        ! Buckingham tappered with MDF:: u=f(r)Buck(r)

        a = vdws%param(1, ivdw)
        rho = vdws%param(2, ivdw)
        c = vdws%param(3, ivdw)
        ri = vdws%param(4, ivdw)

        If (Abs(rho) <= zero_plus) Then
          If (Abs(a) <= zero_plus) Then
            rho = 1.0_wp
          Else
            Call error(467)
          End If
        End If

        ! Sigma-epsilon initialisation

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = -1.0_wp
          vdws%sigeps(2, ivdw) = 0.0_wp
        End If

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          Call mbuck(r, A, rho, c, ri, vdws%cutoff, phi, dphi)

          vdws%tab_potential(i, ivdw) = phi
          vdws%tab_force(i, ivdw) = dphi

          ! Sigma-epsilon search

          If ((.not. vdws%l_force_shift) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
            If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
              If (Nint(Sign(1.0_wp, vdws%tab_potential(i - 1, ivdw))) == -Nint(Sign(1.0_wp, vdws%tab_potential(i, ivdw)))) &
                vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
            Else ! find epsilon
              If ((vdws%tab_potential(i - 2, ivdw) >= vdws%tab_potential(i - 1, ivdw) .and. &
                   vdws%tab_potential(i - 1, ivdw) <= vdws%tab_potential(i, ivdw)) .and. &
                  (nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i - 1, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 2, ivdw), vdws%tab_potential(i, ivdw)) .or. &
                   nequal(vdws%tab_potential(i - 1, ivdw), vdws%tab_potential(i, ivdw)))) Then
                vdws%sigeps(2, ivdw) = -vdws%tab_potential(i - 1, ivdw)
              End If
            End If
          End If
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

      Else If (keypot == VDW_126_MDF) Then

        ! LJ tappered with MDF:: u=f(r)LJ(r)

        a = vdws%param(1, ivdw)
        b = vdws%param(2, ivdw)
        ri = vdws%param(3, ivdw)

        Do i = 1, vdws%max_grid
          r = Real(i, wp) * dlrpot

          Call mlj126(r, A, B, ri, vdws%cutoff, phi, dphi)
          vdws%tab_potential(i, ivdw) = phi
          vdws%tab_force(i, ivdw) = dphi
        End Do
        vdws%tab_potential(0, ivdw) = Huge(vdws%tab_potential(1, ivdw))
        vdws%tab_force(0, ivdw) = Huge(vdws%tab_force(1, ivdw))

        If (.not. vdws%l_force_shift) Then
          vdws%sigeps(1, ivdw) = sig
          vdws%sigeps(2, ivdw) = eps
        End If

      Else

        If (.not. vdws%l_tab) Call error(150)

      End If

      ! no shifting to shifted n-m and DPD
      If (vdws%l_force_shift .and. (keypot /= VDW_N_M_SHIFT .and. keypot /= VDW_DPD)) Then

        vdws%sigeps(1, ivdw) = -1.0_wp
        vdws%sigeps(2, ivdw) = 0.0_wp

        Do i = 1, vdws%max_grid - 4
          t = vdws%tab_potential(i, ivdw) + &
              vdws%tab_force(vdws%max_grid - 4, ivdw) * (Real(i, wp) * dlrpot / vdws%cutoff - 1.0_wp) - &
              vdws%tab_potential(vdws%max_grid - 4, ivdw)
          t1 = vdws%tab_potential(i - 1, ivdw) + &
               vdws%tab_force(vdws%max_grid - 4, ivdw) * (Real(i - 1, wp) * dlrpot / vdws%cutoff - 1.0_wp) - &
               vdws%tab_potential(vdws%max_grid - 4, ivdw)

          ! Sigma-epsilon search
          If (i > 20) Then ! Assumes some safety against numeric black holes!!!
            If (Sign(1.0_wp, vdws%sigeps(1, ivdw)) < 0.0_wp) Then ! find sigma
              If (Nint(Sign(1.0_wp, t1)) == -Nint(Sign(1.0_wp, t))) &
                vdws%sigeps(1, ivdw) = (Real(i, wp) - 0.5_wp) * dlrpot
            Else ! find epsilon
              t2 = vdws%tab_potential(i - 2, ivdw) + &
                   vdws%tab_force(vdws%max_grid - 4, ivdw) * (Real(i - 2, wp) * dlrpot / vdws%cutoff - 1.0_wp) - &
                   vdws%tab_potential(vdws%max_grid - 4, ivdw)
              If ((t2 >= t1 .and. t1 <= t) .and. &
                  (nequal(t2, t1) .or. nequal(t2, t) .or. nequal(t1, t))) Then
                vdws%sigeps(2, ivdw) = -t1
              End If
            End If
          End If
        End Do
      End If

      ! Needed to distinguish that something has been defined

      If (Abs(vdws%tab_potential(0, ivdw)) <= zero_plus) Then
        vdws%tab_potential(0, ivdw) = Sign(Tiny(vdws%tab_potential(0, ivdw)), vdws%tab_potential(0, ivdw))
      End If

    End Do

  End Subroutine vdw_generate

  Subroutine vdw_forces(iatm, xxt, yyt, zzt, rrt, engvdw, virvdw, stats, neigh, vdws, config)

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

    Integer                     :: ai, aj, global_id_i, global_id_j, ityp, jatm, k, key, l, m, mm, &
                                   n
    Real(Kind=wp)               :: a, alpha, b, beta, c, d, e0, eng, eps, fix, fiy, fiz, fx, fy, &
                                   fz, gamma, gk, gk1, gk2, kk, mr, nr, ppp, r0, r0rm, r0rn, r_6, &
                                   r_rrr, r_rrv, r_rsq, r_rvdw, rc, rho, ri, rm, rrr, rscl, rsq, &
                                   sig, sor6, t, t1, t2, t3, vk, vk1, vk2, z1, z2
    Real(Kind=wp), Dimension(9) :: stress_temp, stress_temp_comp

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

    stress_temp = 0.0_wp

    ! global identity and type of iatm

    global_id_i = config%ltg(iatm)
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
      global_id_j = config%ltg(jatm)

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
        r_rsq = 1.0_wp / rsq
        r_rrv = r_rrr * r_rvdw
        rscl = rrr * r_rvdw

        ! Zero energy and force components

        eng = 0.0_wp
        gamma = 0.0_wp

        If (vdws%l_direct) Then ! direct calculation

          If (ityp == VDW_12_6) Then

            ! 12-6 potential :: u=a/r^12-b/r^6

            a = vdws%param(1, k)
            b = vdws%param(2, k)

            r_6 = r_rsq**3

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = r_6 * (a * r_6 - b)
            gamma = 6.0_wp * r_6 * (2.0_wp * a * r_6 - b) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_LENNARD_JONES) Then

            ! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)

            sor6 = (sig**2 * r_rsq)**3

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = 4.0_wp * eps * sor6 * (sor6 - 1.0_wp)
            gamma = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_N_M) Then

            ! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(r0/r)^m]

            e0 = vdws%param(1, k)
            nr = vdws%param(2, k)
            mr = vdws%param(3, k)
            r0 = vdws%param(4, k)

            a = r0 * r_rrr
            b = 1.0_wp / (nr - mr)
            r0rn = a**nr
            r0rm = a**mr

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = e0 * (mr * r0rn - nr * r0rm) * b
            gamma = e0 * mr * nr * (r0rn - r0rm) * b * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_BUCKINGHAM) Then

            ! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

            a = vdws%param(1, k)
            rho = vdws%param(2, k)
            c = vdws%param(3, k)
            ! since this involves any constants shall not be tested here, probably read_field is where any sanity checks should
            ! happen and to the GP for the rest - ame
            If (Abs(rho) <= zero_plus) Then
              If (Abs(a) <= zero_plus) Then
                rho = 1.0_wp
              Else
                Call error(467)
              End If
            End If

            b = rrr / rho
            t1 = a * Exp(-b)
            t2 = -c * r_rsq**3

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1 + t2
            gamma = (t1 * b + 6.0_wp * t2) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_BORN_HUGGINS_MEYER) Then

            ! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            sig = vdws%param(3, k)
            c = vdws%param(4, k)
            d = vdws%param(5, k)

            t1 = a * Exp(b * (sig - rrr))
            t2 = -c * r_rsq**3
            t3 = -d * r_rsq**4

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1 + t2 + t3
            gamma = (t1 * rrr * b + 6.0_wp * t2 + 8.0_wp * t3) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_HYDROGEN_BOND) Then

            ! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

            a = vdws%param(1, k)
            b = vdws%param(2, k)

            t1 = a * r_rsq**6
            t2 = -b * r_rsq**5

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1 + t2
            gamma = (12.0_wp * t1 + 10.0_wp * t2) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_N_M_SHIFT) Then

            ! shifted and force corrected n-m potential (w.smith) ::

            e0 = vdws%param(1, k)
            n = Nint(vdws%param(2, k)); nr = Real(n, wp)
            m = Nint(vdws%param(3, k)); mr = Real(m, wp)
            r0 = vdws%param(4, k)
            rc = vdws%param(5, k); If (rc < 1.0e-6_wp) rc = vdws%cutoff

            If (n <= m) Call error(470)

            t = Real(n - m, wp)

            b = 1.0_wp / t
            c = rc / r0; If (c < 1.0_wp) Call error(468)

            beta = c * ((c**(m + 1) - 1.0_wp) / (c**(n + 1) - 1.0_wp))**b
            alpha = -t / (mr * (beta**n) * (1.0_wp + (nr / c - nr - 1.0_wp) / c**n) &
                          - nr * (beta**m) * (1.0_wp + (mr / c - mr - 1.0_wp) / c**m))
            e0 = e0 * alpha

            If (rrr <= rc) Then
              a = r0 * r_rrr

              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = e0 * (mr * (beta**n) * (a**n - (1.0_wp / c)**n) &
                            - nr * (beta**m) * (a**m - (1.0_wp / c)**m) &
                            + nr * mr * ((rrr / rc - 1.0_wp) * ((beta / c)**n - (beta / c)**m))) * b
              gamma = e0 * Real(m * n, wp) * ((beta**n) * a**n - (beta**m) * a**m &
                                              - rrr / rc * ((beta / c)**n - (beta / c)**m)) * b * r_rsq
            End If

          Else If (ityp == VDW_MORSE) Then

            ! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}

            e0 = vdws%param(1, k)
            r0 = vdws%param(2, k)
            kk = vdws%param(3, k)

            t1 = Exp(-kk * (rrr - r0))

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = e0 * t1 * (t1 - 2.0_wp)
            gamma = -2.0_wp * e0 * kk * t1 * (1.0_wp - t1) * r_rrr

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_WCA) Then

            ! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
            ! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)
            d = vdws%param(3, k)

            If (rrr < vdws%param(4, k) .or. Abs(rrr - d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6 = (sig / (rrr - d))**6

              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = 4.0_wp * eps * sor6 * (sor6 - 1.0_wp) + eps
              gamma = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - 1.0_wp) / (rrr * (rrr - d))

              If (vdws%l_force_shift) Then ! force-shifting
                If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                  eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
                gamma = gamma - vdws%afs(k) * r_rrr
              End If
            End If

          Else If (ityp == VDW_DPD) Then

            ! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

            a = vdws%param(1, k)
            rc = vdws%param(2, k)

            If (rrr < rc) Then ! Else leave them zeros
              t2 = rrr / rc
              t1 = 0.5_wp * a * rrr * (1.0_wp - t2)

              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = t1 * (1.0_wp - t2)
              gamma = t1 * (3.0_wp * t2 - 1.0_wp) * r_rsq
            End If

          Else If (ityp == VDW_AMOEBA) Then

            ! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)

            rho = rrr / sig
            t1 = 1.0_wp / (0.07_wp + rho)
            t2 = 1.0_wp / (0.12_wp + rho**7)
            t3 = eps * (1.07_wp * t1)**7

            t = t3 * ((1.12_wp * t2) - 2.0_wp)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t
            gamma = 7.0_wp * (t1 * t + 1.12_wp * t3 * t2**2 * rho**6) * rho * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_LENNARD_JONES_COHESIVE) Then

            ! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)
            c = vdws%param(3, k)

            sor6 = (sig**2 * r_rsq)**3

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = 4.0_wp * eps * sor6 * (sor6 - c)
            gamma = 24.0_wp * eps * sor6 * (2.0_wp * sor6 - c) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_MORSE_12) Then

            ! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}+c/r^12

            e0 = vdws%param(1, k)
            r0 = vdws%param(2, k)
            kk = vdws%param(3, k)
            c = vdws%param(4, k)

            t1 = Exp(-kk * (rrr - r0))
            sor6 = c * r_rsq**6

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = e0 * t1 * (t1 - 2.0_wp) + sor6
            gamma = -2.0_wp * e0 * kk * t1 * (1.0_wp - t1) * r_rrr - 12.0_wp * sor6 * r_rrr

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_RYDBERG) Then

            ! Rydberg potential:: u=(a+b*r)Exp(-r/c)

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            c = vdws%param(3, k)

            kk = rrr / c
            t1 = Exp(-kk)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = (a + b * rrr) * t1
            gamma = kk * t1 * (a - b * c + b * rrr) * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_ZBL) Then

            ! ZBL potential:: u=Z1Z2/(40r)_{i=1}^4b_ie^{-c_i*r/a}

            z1 = vdws%param(1, k)
            z2 = vdws%param(2, k)

            ! this is in fact inverse a
            a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
            kk = z1 * z2 * r4pie0

            Call zbl(rrr, kk, a, t1, gamma)
            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1
            gamma = gamma * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_ZBL_SWITCH_MORSE) Then

            ! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

            z1 = vdws%param(1, k)
            z2 = vdws%param(2, k)
            rm = vdws%param(3, k)
            c = 1.0_wp / vdws%param(4, k)
            e0 = vdws%param(5, k)
            r0 = vdws%param(6, k)
            t2 = vdws%param(7, k)

            ! this is in fact inverse a
            a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
            kk = z1 * z2 * r4pie0

            Call zbls(rrr, kk, a, rm, c, e0, t2, r0, t1, gamma)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1
            gamma = gamma * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_ZBL_SWITCH_BUCKINGHAM) Then

            ! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

            z1 = vdws%param(1, k)
            z2 = vdws%param(2, k)
            rm = vdws%param(3, k)
            c = 1.0_wp / vdws%param(4, k)
            e0 = vdws%param(5, k)
            r0 = vdws%param(6, k)
            t2 = vdws%param(7, k)

            ! this is in fact inverse a
            a = (z1**0.23_wp + z2**0.23_wp) / (ab * 0.88534_wp)
            kk = z1 * z2 * r4pie0

            Call zblb(rrr, kk, a, rm, c, e0, r0, t2, t1, gamma)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1
            gamma = gamma * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_LJ_MDF) Then

            ! LJ tappered with MDF:: u=f(r)LJ(r)

            eps = vdws%param(1, k)
            sig = vdws%param(2, k)
            ri = vdws%param(3, k)
            !rc=vdws%cutoff

            Call mlj(rrr, eps, sig, ri, vdws%cutoff, t1, gamma)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1
            gamma = gamma * r_rsq

            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_BUCKINGHAM_MDF) Then

            ! Buckingham tappered with MDF:: u=f(r)Buck(r)
            a = vdws%param(1, k)
            rho = vdws%param(2, k)
            c = vdws%param(3, k)
            ri = vdws%param(4, k)

            Call mbuck(rrr, a, rho, c, ri, vdws%cutoff, t1, gamma)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1
            gamma = gamma * r_rsq

            ! by construction is zero outside vdws%cutoff so no shifting
            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (ityp == VDW_126_MDF) Then

            ! LJ tappered with MDF:: u=f(r)LJ12-6(r)

            a = vdws%param(1, k)
            b = vdws%param(2, k)
            ri = vdws%param(3, k)
            !rc=vdws%cutoff

            Call mlj126(rrr, a, b, ri, vdws%cutoff, t1, gamma)

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
              eng = t1
            gamma = gamma * r_rsq

            ! by construction is zero outside vdws%cutoff so no shifting
            If (vdws%l_force_shift) Then ! force-shifting
              If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) &
                eng = eng + vdws%afs(k) * rrr + vdws%bfs(k)
              gamma = gamma - vdws%afs(k) * r_rrr
            End If

          Else If (Abs(vdws%tab_potential(0, k)) > zero_plus) Then ! potential read from TABLE - (ityp == VDW_TAB)

            l = Int(rrr * vdws%rdr)
            ppp = rrr * vdws%rdr - Real(l, wp)

            ! calculate interaction energy using 3-point interpolation

            If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) Then
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
            End If

            ! calculate forces using 3-point interpolation

            gk = vdws%tab_force(l, k); If (l == 0) gk = gk * rrr
            gk1 = vdws%tab_force(l + 1, k)
            gk2 = vdws%tab_force(l + 2, k)

            t1 = gk + (gk1 - gk) * ppp
            t2 = gk1 + (gk2 - gk1) * (ppp - 1.0_wp)

            gamma = (t1 + (t2 - t1) * ppp * 0.5_wp) * r_rsq
            If (vdws%l_force_shift) gamma = gamma - vdws%tab_force(vdws%max_grid - 4, k) * r_rrv ! force-shifting

          End If

        Else If (Abs(vdws%tab_potential(0, k)) > zero_plus) Then ! no direct = fully tabulated calculation

          l = Int(rrr * vdws%rdr)
          ppp = rrr * vdws%rdr - Real(l, wp)

          ! calculate interaction energy using 3-point interpolation

          If (stats%collect_pp .or. jatm <= config%natms .or. global_id_i < global_id_j) Then
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
          End If

          ! calculate forces using 3-point interpolation

          gk = vdws%tab_force(l, k); If (l == 0) gk = gk * rrr
          gk1 = vdws%tab_force(l + 1, k)
          gk2 = vdws%tab_force(l + 2, k)

          t1 = gk + (gk1 - gk) * ppp
          t2 = gk1 + (gk2 - gk1) * (ppp - 1.0_wp)

          gamma = (t1 + (t2 - t1) * ppp * 0.5_wp) * r_rsq
          If (vdws%l_force_shift) gamma = gamma - vdws%tab_force(vdws%max_grid - 4, k) * r_rrv ! force-shifting

        End If

        ! calculate forces

        fx = gamma * xxt(mm)
        fy = gamma * yyt(mm)
        fz = gamma * zzt(mm)

        fix = fix + fx
        fiy = fiy + fy
        fiz = fiz + fz

        If (jatm <= config%natms) Then

          config%parts(jatm)%fxx = config%parts(jatm)%fxx - fx
          config%parts(jatm)%fyy = config%parts(jatm)%fyy - fy
          config%parts(jatm)%fzz = config%parts(jatm)%fzz - fz

        End If

        If (jatm <= config%natms .or. global_id_i < global_id_j) Then

          ! add interaction energy

          engvdw = engvdw + eng

          ! add virial

          virvdw = virvdw - gamma * rsq

          ! add stress tensor

          stress_temp_comp = calculate_stress([xxt(mm), yyt(mm), zzt(mm)], [fx, fy, fz])
          stress_temp = stress_temp + stress_temp_comp

        End If

        If (stats%collect_pp) Then
          stress_temp_comp = calculate_stress([xxt(mm), yyt(mm), zzt(mm)], [fx, fy, fz])
          stats%pp_energy(iatm) = stats%pp_energy(iatm) + eng * 0.5_wp
          stats%pp_stress(:, iatm) = stats%pp_stress(:, iatm) + stress_temp_comp * 0.5_wp
          If (jatm <= config%natms) Then
            stats%pp_energy(jatm) = stats%pp_energy(jatm) + eng * 0.5_wp
            stats%pp_stress(:, jatm) = stats%pp_stress(:, jatm) + stress_temp_comp * 0.5_wp
          End If
        End If

      End If

    End Do

    ! load back forces

    config%parts(iatm)%fxx = fix
    config%parts(iatm)%fyy = fiy
    config%parts(iatm)%fzz = fiz

    ! complete stress tensor

    stats%stress = stats%stress + stress_temp

  End Subroutine vdw_forces

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

    If (Allocated(T%sigeps)) Then
      Deallocate (T%sigeps)
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
