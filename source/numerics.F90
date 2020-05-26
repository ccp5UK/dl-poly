Module numerics

  Use comms,           Only: comms_type
  Use constants,       Only: epsilon_wp,&
                             half_minus,&
                             rt2,&
                             rt3,&
                             sqrpi,&
                             zero_plus
  Use errors_warnings, Only: error
  Use kinds,           Only: li,&
                             wi,&
                             wp
  Use particle,        Only: corePart

  Implicit None
  Private

  !> Random seed type
  Type, Public :: seed_type
    Private

    !> Flag indicating whether the seed has been initialised
    Logical, Public          :: defined = .false.
    !> The seed
    Integer(Kind=wi), Public :: seed(1:3)
    !> state variables for uni random number generator. In long run one wants to move to a better random number generaot
    Logical                  :: newjob = .true.
    Logical                  :: newjob_bm = .true.
    Integer                  :: ir, jr
    Real(Kind=wp)            :: c, cd, cm, u(1:97)

  Contains
    Private
    Procedure, Public :: init => init_seed
  End Type seed_type

  ! Private copies to avoid inheritence loop
  Integer, Parameter :: IMCON_NOPBC = 0
  Integer, Parameter :: IMCON_CUBIC = 1
  Integer, Parameter :: IMCON_ORTHORHOMBIC = 2
  Integer, Parameter :: IMCON_PARALLELOPIPED = 3
  Integer, Parameter :: IMCON_SLAB = 6
  ! REMOVED -- DL_POLY 2 ONLY
  Integer, Parameter :: IMCON_TRUNC_OCTO = 4
  Integer, Parameter :: IMCON_RHOMBIC_DODEC = 5
  Integer, Parameter :: IMCON_HEXAGONAL = 7


  !!!!!!!!!!!!!!!!!!!!!!!! THIS IS NUMERIC_CONTAINER !!!!!!!!!!!!!!!!!!!!!
  !
  ! Function uni - two seeded random number generator
  !
  ! Function sarurnd - three seeded random number generator based on SARU
  !
  ! Subroutine box_mueller_saru - generates gaussian random numbers of unit
  !                               variance (with zero mean and standard
  !                               variation of 1)
  !
  ! Subroutine box_mueller_uni - generates gaussian random numbers of unit
  !                              variance (with zero mean and standard
  !                              variation of 1)
  !
  ! Subroutine gauss_1 - constructs velocity arrays with a gaussian
  !                      distribution of unit variance (zero mean) by
  !                      an approximation of the Central Limit Theorem
  !
  ! Subroutine gauss_2 - constructs velocity arrays with a gaussian
  !                      distribution of unit variance (zero mean) using
  !                      the box-mueller method
  !
  ! Subroutine erfcgen - generates interpolation tables for erfc and erfc/r
  !                      derivative
  !
  ! Function match - determines a match between integer value 'n' and an
  !                  array of integers in ascending order
  !
  ! Subroutine shellsort - sorts an integer array in ascending order
  !
  ! Subroutine shellsort2 - sorts an integer array in ascending order,
  !                         keeping the original ranking of the array
  !
  ! Function local_index - finds the local atom number given the global
  !                        atom number
  !
  ! Subroutine dcell - calculates the dimensional properties of a
  !                    simulation cell
  !
  ! Subroutine invert - calculates the inverse of a 3x3 matrix using
  !                     cofactors
  !
  ! Subroutine images - calculates the minimum image distance of
  !                     atom pairs within a specified MD cell
  !
  ! Subroutine images_s - calculates the minimum image distance of
  !                       a single atom pair within a specified MD cell
  !
  ! Subroutine pbcshift - calculates the minimum image of atoms within
  !                       a specified MD cell in accordance with the DD
  !                       boundary convention
  !
  ! Subroutine pbcshfrc - calculates the minimum image of atoms within
  !                       a specified MD cell in accordance with the DD
  !                       boundary convention and returns reduced coordinates
  !
  ! Subroutine pbcshfrl - calculates the minimum image of atoms within
  !                       a specified MD cell in accordance with the DD
  !                       boundary convention from reduced coordinates and
  !                       returns real coordinates
  !
  ! Subroutine jacobi - diagonalises real symmetric matrices by the
  !                     Jacobi method
  !
  ! Subroutine mat_mul - calculates product of two 3x3 matrices written
  !                     in a DL_POLY format as vectors
  !
  ! Function Factorial - computes the factorial (n!) of an input n
  !
  ! MM3_Module - a module to calculate LTC due to AMOEBA 14-7 buffered vdw
  !              interactions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Public :: uni
  Public :: sarurnd
  Public :: match
  Public :: local_index
  Public :: Factorial
  Public :: box_mueller_saru1
  Public :: box_mueller_saru2
  Public :: box_mueller_saru3
  Public :: box_mueller_saru6
  Public :: box_mueller_uni
  Public :: gauss_1
  Public :: gauss_2
  Public :: erfcgen
  Public :: shellsort
  Public :: shellsort2
  Public :: dcell
  Public :: invert
  Public :: images
  Public :: images_s
  Public :: pbcshift
  Public :: pbcshfrc
  Public :: pbcshfrl
  Public :: jacobi
  Public :: mat_mul
  Public :: factor
  Public :: get_nth_prime

  Interface pbcshfrc
    Module Procedure pbcshfrc_parts
    Module Procedure pbcshfrc_arrays
  End Interface pbcshfrc

  Interface pbcshfrl
    Module Procedure pbcshfrl_parts
    Module Procedure pbcshfrl_arrays
  End Interface pbcshfrl

  Interface pbcshift
    Module Procedure pbcshift_parts
    Module Procedure pbcshift_arrays
  End Interface

  Interface equal
    Module Procedure equal_real_wp
  End Interface equal

  Interface nequal
    Module Procedure nequal_real_wp
  End Interface nequal

  Public :: equal, nequal

Contains

  !> Initialise a seed
  Subroutine init_seed(T, seed)
    Class(seed_type)                :: T
    Integer(Kind=wi), Intent(In   ) :: seed(1:3)

    T%defined = .true.
    T%seed(1:3) = seed(1:3)
  End Subroutine init_seed

  Function uni(seed, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 random number generator based on the universal random number
    ! generator of marsaglia, zaman and tsang.
    !
    ! Ref: stats. and prob. lett. 8 (1990) 35-39.)
    !
    ! Note: It returns in [0,1)
    !
    ! This random number generator originally appeared in "Toward a
    ! Universal Random Number Generator" by George Marsaglia, Arif Zaman and
    ! W.W. Tsang in Florida State University Report: FSU-SCRI-87-50 (1987).
    ! It was later modified by F. James and published in "A Review of
    ! Pseudo-random Number Generators".
    ! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
    ! It passes ALL of the tests for random number generators and has a
    ! period of 2^144, is completely portable (gives bit identical results
    ! on all machines with at least 24-bit mantissas in the floating point
    ! representation).
    ! The algorithm is a combination of a Fibonacci sequence (with lags of
    ! 97 and 33, and operation "subtraction plus one, modulo one") and an
    ! "arithmetic sequence" (using subtraction).
    ! Use IJ = 1802 & KL = 9373 (idnode=0) to test the random number
    ! generator. The subroutine RANMAR should be used to generate 20000
    ! random numbers.  Then display the next six random numbers generated
    ! multiplied by 4096*4096.  If the random number generator is working
    ! properly, the random numbers should be:
    !         6533892.0  14220222.0  7275067.0
    !         6172232.0  8354498.0   10633180.0
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov april 2008
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type),  Intent(InOut) :: seed
    Type(comms_type), Intent(In   ) :: comm
    Real(Kind=wp)                   :: uni

    Integer       :: i, ii, ij, j, jj, k, kl, l, m
    Real(Kind=wp) :: s, t

    ! initialise parameters u,c,cd,cm

    If (seed%newjob .or. seed%defined) Then
      seed%newjob = .false.

      ! If no seeding is specified then default to DL_POLY scheme

      If (seed%defined) Then

        seed%defined = .false.

        ! First random number seed must be between 0 and 31328
        ! Second seed must have a value between 0 and 30081

        ij = Mod(Abs(seed%seed(1) + comm%idnode), 31328)
        i = Mod(ij / 177, 177) + 2;
        j = Mod(ij, 177) + 2;
        kl = Mod(Abs(seed%seed(2) + comm%idnode), 30081)
        k = Mod(kl / 169, 178) + 1
        l = Mod(kl, 169)

      Else

        ! initial values of i,j,k must be in range 1 to 178 (not all 1)
        ! initial value of l must be in range 0 to 168

        i = Mod(comm%idnode, 166) + 12
        j = Mod(comm%idnode, 144) + 34
        k = Mod(comm%idnode, 122) + 56
        l = Mod(comm%idnode, 90) + 78

      End If

      seed%ir = 97
      seed%jr = 33

      Do ii = 1, 97

        s = 0.0_wp
        t = 0.5_wp

        Do jj = 1, 24

          m = Mod(Mod(i * j, 179) * k, 179)
          i = j
          j = k
          k = m
          l = Mod(53 * l + 1, 169)
          If (Mod(l * m, 64) >= 32) s = s + t
          t = 0.5_wp * t

        End Do

        seed%u(ii) = s

      End Do

      seed%c = 362436.0_wp / 16777216.0_wp
      seed%cd = 7654321.0_wp / 16777216.0_wp
      seed%cm = 16777213.0_wp / 16777216.0_wp

    End If

    ! calculate random number

    uni = seed%u(seed%ir) - seed%u(seed%jr)
    If (uni < 0.0_wp) uni = uni + 1.0_wp

    seed%u(seed%ir) = uni

    seed%ir = seed%ir - 1
    If (seed%ir == 0) seed%ir = 97

    seed%jr = seed%jr - 1
    If (seed%jr == 0) seed%jr = 97

    seed%c = seed%c - seed%cd
    If (seed%c < 0.0_wp) seed%c = seed%c + seed%cm

    uni = uni - seed%c
    If (uni < 0.0_wp) uni = uni + 1.0_wp

  End Function uni

  Function sarurnd(seed, seeda, seedb, seedc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine random number generator based on the saru random
    ! number generator of Steve Worley with three integer seeds.
    !
    ! Ref: Comp. Phys. Comms. 184 (2013) 1119-1128
    !
    ! Note: It returns in [0,1)
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type), Intent(InOut) :: seed
    Integer,         Intent(In   ) :: seeda, seedb, seedc
    Real(Kind=wp)                  :: sarurnd

    Integer(Kind=li), Parameter :: two32 = 2_li**32
    Real(Kind=wp), Parameter    :: rtwo32 = (2.0_wp)**(-32), two32r = 2.0_wp**32

    Integer(Kind=li) :: itmp, seed1, seed2, seed3, state, u32, v, wstate
    Real(Kind=wp)    :: statepart1, statepart2, tmp

    ! Apply possible shifting - usually (0,0,0)

    If (.not. seed%defined) Then
      seed1 = Int(seeda, Kind=li)
      seed2 = Int(seedb, Kind=li)
      seed3 = Int(seedc, Kind=li)
    Else
      seed1 = Int(seeda, Kind=li) + Int(seed%seed(1), Kind=li)
      seed2 = Int(seedb, Kind=li) + Int(seed%seed(2), Kind=li)
      seed3 = Int(seedc, Kind=li) + Int(seed%seed(3), Kind=li)
    End If

    ! Wrap up

    seed1 = Mod(seed1, two32); If (seed1 < 0) seed1 = seed1 + two32
    seed2 = Mod(seed2, two32); If (seed2 < 0) seed2 = seed2 + two32
    seed3 = Mod(seed3, two32); If (seed3 < 0) seed3 = seed3 + two32

    ! apply premixing to seeds

    ! seed3 ^= (seed1<<7)^(seed2>>6);

    seed3 = Ieor(seed3, Ieor(Ishft(seed1, 7_li), Ishft(seed2, -6_li)))
    seed3 = Mod(seed3, two32)

    ! seed2 += (seed1>>4)^(seed3>>15);

    itmp = Ieor(Ishft(seed1, -4_li), Ishft(seed3, -15_li))
    seed2 = seed2 + itmp
    seed2 = Mod(seed2, two32)

    ! seed1 ^= (seed2<<9)+(seed3<<8);

    seed1 = Ieor(seed1, Ishft(seed2, 9_li) + Ishft(seed3, 8_li))
    seed1 = Mod(seed1, two32)

    ! seed3 ^= 0xA5366B4D*((seed2>>11) ^ (seed1<<1));

    itmp = Ieor(Ishft(seed2, -11_li), Ishft(seed1, 1_li))
    tmp = 2771807053.0_wp * Real(itmp, Kind=wp)
    itmp = Int(Mod(tmp, two32r), Kind=li)
    seed3 = Ieor(seed3, itmp)
    seed3 = Mod(seed3, two32)

    ! seed2 += 0x72BE1579*((seed1<<4) ^ (seed3>>16));

    itmp = Ieor(Ishft(seed1, 4_li), Ishft(seed3, -16_li))
    tmp = 1925059961.0_wp * Real(itmp, Kind=wp)
    itmp = Int(Mod(tmp, two32r), Kind=li)
    seed2 = seed2 + itmp
    seed2 = Mod(seed2, two32)

    ! seed1 ^= 0X3F38A6ED*((seed3>>5) ^ (((signed int)seed2)>>22));

    itmp = Ieor(Ishft(seed3, -5_li), Ishft(Int(Int(seed2), Kind=li), -22_li))
    tmp = 1060677357.0_wp * Real(itmp, Kind=wp)
    itmp = Int(Mod(tmp, two32r), Kind=li)
    If (itmp < 0) itmp = itmp + two32
    seed1 = Ieor(seed1, itmp)
    seed1 = Mod(seed1, two32)

    ! seed2 += seed1*seed3;

    tmp = Real(seed1, Kind=wp) * Real(seed3, Kind=wp)
    itmp = Int(Mod(tmp, two32r), Kind=li)
    seed2 = seed2 + itmp
    seed2 = Mod(seed2, two32)

    ! seed1 += seed3 ^ (seed2>>2);

    seed1 = seed1 + Ieor(seed3, Ishft(seed2, -2_li))
    seed1 = Mod(seed1, two32)

    ! seed2 ^= ((signed int)seed2)>>17;

    seed2 = Ieor(seed2, Ishft(Int(Int(seed2), Kind=li), -17_li))
    If (seed2 < 0) seed2 = seed2 + two32

    ! convert seeds to state values

    ! state = 0x79dedea3*(seed1^(((signed int)seed1)>>14));

    itmp = Ieor(seed1, Ishft(Int(Int(seed1), Kind=li), -14_li))
    tmp = 2044649123.0_wp * Real(itmp, Kind=wp)
    state = Int(Mod(tmp, two32r), Kind=li)
    If (state < 0) state = state + two32

    ! wstate = (state + seed2) ^ (((signed int)state)>>8);

    wstate = Ieor((state + seed2), Ishft(Int(Int(state), Kind=li), -8_li))
    wstate = Mod(wstate, two32); If (wstate < 0) wstate = wstate + two32

    ! state = state + (wstate*(wstate^0xdddf97f5));

    itmp = Ieor(wstate, 3722418165_li)
    tmp = Real(wstate, Kind=wp) * Real(itmp, Kind=wp)
    itmp = Int(Mod(tmp, two32r), Kind=li)
    state = state + itmp
    state = Mod(state, two32)

    ! wstate = 0xABCB96F7 + (wstate>>1);

    wstate = 2882246391_li + Ishft(wstate, -1_li)
    wstate = Mod(wstate, two32)

    ! advance LCG state by 1
    ! state = 0x4beb5d59*state + 0x2600e1f7; // LCG

    tmp = 1273716057.0_wp * Real(state, Kind=wp)
    itmp = Int(Mod(tmp, two32r), Kind=li)
    state = itmp + 637592055_li
    state = Mod(state, two32)

    ! advance Weyl state by 1, for oWeylOffset=0x8009d14b & oWeylPeriod=0xda879add
    ! wstate = wstate + oWeylOffset + ((((signed int)wstate)>>31)&oWeylPeriod); // OWS

    wstate = wstate + 2148127051_li + Iand(Ishft(Int(Int(wstate), Kind=li), -31_li), 3666320093_li)
    wstate = Mod(wstate, two32)

    ! calculate 32-bit pseudo-random number

    ! v = (state ^ (state>>26))+wstate;

    v = Ieor(state, Ishft(state, -26_li)) + wstate
    v = Mod(v, two32)

    ! u32 = (v^(v>>20))*0x6957f5a7;

    itmp = Ieor(v, Ishft(v, -20_li))
    tmp = Real(itmp, Kind=wp) * 1767372199.0_wp
    u32 = Int(Mod(tmp, two32r), Kind=li)

    ! convert to real (double-precision) number between 0 and 1

    ! statep1 = ((signed int)u32)*TWO_N32+(0.5+0.5*TWO_N32);

    statepart1 = Real(Int(Int(u32), Kind=li), Kind=wp) * rtwo32 + (0.5_wp + 0.5_wp * rtwo32)

    ! statep2 = state*(TWO_N32*TWO_N32);

    statepart2 = Real(state, Kind=wp) * rtwo32**2

    ! sarurand=statep1+statep2

    sarurnd = statepart1 + statepart2

  End Function sarurnd

  Subroutine box_mueller_saru1(seed, i, j, gauss1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine using the box-mueller method for generating a
    ! gaussian random number of unit variance (with zero mean and standard
    ! variation of 1).
    !
    ! dependent on sarurnd
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type), Intent(InOut) :: seed
    Integer,         Intent(In   ) :: i, j
    Real(Kind=wp),   Intent(  Out) :: gauss1

    Integer       :: k
    Real(Kind=wp) :: ran0, ran1, ran2

    ! Initialise counter

    k = 0

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random numbers 1 & 2

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss1 = ran0 * ran1

  End Subroutine box_mueller_saru1

  Subroutine box_mueller_saru2(seed, i, j, n, gauss1, l_str)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine for generating a gaussian random number of unit
    ! variance (with zero mean and standard variation of 1) by using either
    ! the box_mueller method or the simplest CLT approximation
    !
    ! dependent on sarurnd
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type), Intent(InOut) :: seed
    Integer,         Intent(In   ) :: i, j, n
    Real(Kind=wp),   Intent(  Out) :: gauss1
    Logical,         Intent(In   ) :: l_str

    Integer       :: k
    Real(Kind=wp) :: ran0, ran1, ran2

    ! Initialise counter

    k = n

    ! Avoid overflow due to Box_Mueller rejection rate of sampling = (1-p/4) = 0.214

    If (Huge(1) - k <= 50) k = -Huge(1) + 49

    ! CLT approximation

    If (.not. l_str) Then
      ran0 = rt3 * (2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp)
      Return
    End If

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random numbers 1 & 2

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss1 = ran0 * ran1

  End Subroutine box_mueller_saru2

  Subroutine box_mueller_saru3(seed, i, j, gauss1, gauss2, gauss3)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine using the box-mueller method for generating 3
    ! gaussian random numbers of unit variance (with zero mean and standard
    ! variation of 1).
    !
    ! dependent on sarurnd
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type), Intent(InOut) :: seed
    Integer,         Intent(In   ) :: i, j
    Real(Kind=wp),   Intent(  Out) :: gauss1, gauss2, gauss3

    Integer       :: k
    Real(Kind=wp) :: ran0, ran1, ran2

    ! Initialise counter

    k = 0

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random numbers 1 & 2

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss1 = ran0 * ran1
    gauss2 = ran0 * ran2

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random number 3 & 4

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss3 = ran0 * ran1

  End Subroutine box_mueller_saru3

  Subroutine box_mueller_saru6(seed, i, j, gauss)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine using the box-mueller method for generating 6
    ! gaussian random numbers of unit variance (with zero mean and standard
    ! variation of 1).
    !
    ! dependent on sarurnd
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type), Intent(InOut) :: seed
    Integer,         Intent(In   ) :: i, j
    Real(Kind=wp),   Intent(  Out) :: gauss(6)

    Integer       :: k
    Real(Kind=wp) :: ran0, ran1, ran2

    ! Initialise counter

    k = 0

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random numbers 1 & 2

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss(1) = ran0 * ran1
    gauss(2) = ran0 * ran2

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random number 3 & 4

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss(3) = ran0 * ran1
    gauss(4) = ran0 * ran2

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * sarurnd(seed, i, j, k) - 1.0_wp
      ran2 = 2.0_wp * sarurnd(seed, i, j, k + 1) - 1.0_wp
      ran0 = ran1**2 + ran2**2
      k = k + 2
    End Do

    ! calculate gaussian random number 5 & 6
    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss(5) = ran0 * ran1
    gauss(6) = ran0 * ran2

  End Subroutine box_mueller_saru6

  Subroutine box_mueller_uni(seed, gauss1, gauss2, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine using the box-mueller method for generating
    ! gaussian random numbers of unit variance (with zero mean and standard
    ! variation of 1).  Otherwise, an approximation of the Central Limit
    ! Theorem must be used: G = (1/A)*[Sum_i=1,N(Ri) - AN/2]*(12/N)^(1/2),
    ! where A is the number of outcomes from the random throw Ri and N is
    ! the number of tries.
    !
    ! dependent on uni
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith may 2008
    ! amended   - i.t.todorov february 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type),  Intent(InOut) :: seed
    Real(Kind=wp),    Intent(  Out) :: gauss1, gauss2
    Type(comms_type), Intent(In   ) :: comm

    Real(Kind=wp) :: ran0, ran1, ran2

    ! make sure uni is initialised

    If (seed%newjob_bm) Then
      seed%newjob_bm = .false.
      ran0 = uni(seed, comm)
    End If

    ! generate uniform random numbers on [-1, 1)

    ran0 = 1.0_wp
    Do While (ran0 <= zero_plus .or. ran0 >= 1.0_wp)
      ran1 = 2.0_wp * uni(seed, comm) - 1.0_wp
      ran2 = 2.0_wp * uni(seed, comm) - 1.0_wp
      ran0 = ran1**2 + ran2**2
    End Do

    ! calculate gaussian random numbers

    ran0 = Sqrt(-2.0_wp * Log(ran0) / ran0)
    gauss1 = ran0 * ran1
    gauss2 = ran0 * ran2

  End Subroutine box_mueller_uni

  Subroutine gauss_1(seed, natms, vxx, vyy, vzz, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine for constructing velocity arrays with a gaussian
    ! distribution of unit variance (zero mean), based on the method
    ! described by Allen and Tildesley in "Computer Simulation of Liquids",
    ! Clarendon Press 1987, P347.  It is based on an approximation of the
    ! Central Limit Theorem : G = (1/A)*[Sum_i=1,N(Ri) - AN/2]*(12/N)^(1/2),
    ! where A is the number of outcomes from the random throw Ri and N is
    ! the number of tries.
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov july 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type),               Intent(InOut) :: seed
    Integer,                       Intent(In   ) :: natms
    Real(Kind=wp), Dimension(1:*), Intent(  Out) :: vxx, vyy, vzz
    Type(comms_type),              Intent(In   ) :: comm

    Real(Kind=wp), Parameter :: a1 = 3.949846138_wp, a3 = 0.252408784_wp, a5 = 0.076542912_wp, &
                                a7 = 0.008355968_wp, a9 = 0.029899776_wp

    Integer       :: i, j
    Real(Kind=wp) :: rr2, rrr

    Do i = 1, natms
      rrr = 0.0_wp
      Do j = 1, 12
        rrr = rrr + uni(seed, comm)
      End Do
      rrr = (rrr - 6.0_wp) / 4.0_wp
      rr2 = rrr * rrr
      vxx(i) = rrr * (a1 + rr2 * (a3 + rr2 * (a5 + rr2 * (a7 + rr2 * a9))))

      rrr = 0.0_wp
      Do j = 1, 12
        rrr = rrr + uni(seed, comm)
      End Do
      rrr = (rrr - 6.0_wp) / 4.0_wp
      rr2 = rrr * rrr
      vyy(i) = rrr * (a1 + rr2 * (a3 + rr2 * (a5 + rr2 * (a7 + rr2 * a9))))

      rrr = 0.0_wp
      Do j = 1, 12
        rrr = rrr + uni(seed, comm)
      End Do
      rrr = (rrr - 6.0_wp) / 4.0_wp
      rr2 = rrr * rrr
      vzz(i) = rrr * (a1 + rr2 * (a3 + rr2 * (a5 + rr2 * (a7 + rr2 * a9))))
    End Do

  End Subroutine gauss_1

  Subroutine gauss_2(seed, natms, vxx, vyy, vzz, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine for constructing velocity arrays with a gaussian
    ! distribution of unit variance (zero mean), based on the box-mueller
    ! method
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(seed_type),               Intent(InOut) :: seed
    Integer,                       Intent(In   ) :: natms
    Real(Kind=wp), Dimension(1:*), Intent(  Out) :: vxx, vyy, vzz
    Type(comms_type),              Intent(InOut) :: comm

    Integer       :: i, j
    Real(Kind=wp) :: gauss1, gauss2

    Do i = 1, (natms + 1) / 2
      j = natms + 1 - i

      Call box_mueller_uni(seed, gauss1, gauss2, comm)
      vxx(i) = gauss1
      vxx(j) = gauss2

      Call box_mueller_uni(seed, gauss1, gauss2, comm)
      vyy(i) = gauss1
      vyy(j) = gauss2

      Call box_mueller_uni(seed, gauss1, gauss2, comm)
      vzz(i) = gauss1
      vzz(j) = gauss2
    End Do

  End Subroutine gauss_2

  Subroutine erfcgen(rcut, alpha, ewald_exclusion_grid, erc, fer)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine for generating interpolation tables for erfc and its
    ! derivative - for use with Ewald sum
    !
    ! copyright - daresbury laboratory
    ! author    - t.forester december 1994
    ! amended   - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp),                                    Intent(In   ) :: rcut, alpha
    Integer,                                          Intent(In   ) :: ewald_exclusion_grid
    Real(Kind=wp), Dimension(0:ewald_exclusion_grid), Intent(  Out) :: erc, fer

    Real(Kind=wp), Parameter :: a1 = 0.254829592_wp, a2 = -0.284496736_wp, a3 = 1.421413741_wp, &
                                a4 = -1.453152027_wp, a5 = 1.061405429_wp, pp = 0.3275911_wp

    Integer       :: i
    Real(Kind=wp) :: drewd, exp1, rrr, rsq, tt

    ! look-up tables for real space part of ewald sum

    drewd = rcut / Real(ewald_exclusion_grid - 4, wp)

    Do i = 1, ewald_exclusion_grid
      rrr = Real(i, wp) * drewd
      rsq = rrr * rrr

      tt = 1.0_wp / (1.0_wp + pp * alpha * rrr)
      exp1 = Exp(-(alpha * rrr)**2)

      erc(i) = tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1 / rrr
      fer(i) = (erc(i) + 2.0_wp * (alpha / sqrpi) * exp1) / rsq
    End Do

    ! extrapolation for grid point 0 at distances close to 0

    erc(0) = Huge(1.0_wp)
    fer(0) = Huge(1.0_wp + 2.0_wp * (alpha / sqrpi))

  End Subroutine erfcgen

  Function match(n, ind_top, list)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 function to determine a match between a positive integer
    ! 'n' and an array of positive integer 'list(1:ind_top)' sorted in
    ! ascending order
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov october 2006
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent(In   ) :: n, ind_top, list(1:*)
    Logical                :: match

    Character(Len=256) :: message
    Integer            :: ind_now, ind_old

    If (n < 1) Then
      Write (message, '(a)') 'Wrong value of n has entered in Function match, called when defining ltg arrays in link_cells'
      Call error(0, message)
    End If

    match = .false.

    If (ind_top < 1) Return

    ind_old = 1
    ind_now = 1

    Do
      If (n == list(ind_now)) Then
        match = .true.
        Return
      Else If (n > list(ind_now)) Then
        If (ind_old == ind_top) Return
        ind_old = ind_now
        ind_now = (ind_old + ind_top + 1) / 2
      Else If (n < list(ind_now)) Then
        ind_now = (ind_old + ind_now) / 2
        If (ind_now == ind_old) Return
      End If
    End Do

  End Function match

  Subroutine shellsort(n, list)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 shell sort routine.  Sorts an array of integers into
    ! ascending order.
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                 Intent(In   ) :: n
    Integer, Dimension(1:*), Intent(InOut) :: list

    Integer :: i, index, nl, value
    Logical :: go

    ! set up sort

    If (n > 1) Then

      ! number of lists

      nl = n / 2

      ! iterate shell sort until there is a list

      Do While (nl > 0)

        ! for all lists from next-to-ground-level up to their end

        Do i = nl + 1, n

          value = list(i)
          index = i

          ! Antibubble down between levels of the same list

          go = .true.
          Do While (index > nl .and. go)
            go = (list(index - nl) > value)

            If (go) Then
              list(index) = list(index - nl)
              index = index - nl
            End If
          End Do

          ! Last insertion as close to the ground as it gets

          list(index) = value

        End Do

        ! Decrease the number of lists

        nl = nl / 2

      End Do

    End If

  End Subroutine shellsort

  Subroutine shellsort2(n, rank, list)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 shell sort routine.  Sorts an array of integers (list) into
    ! ascending order.  The original rank of array list is kept in rank.
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                 Intent(In   ) :: n
    Integer, Dimension(1:*), Intent(InOut) :: rank, list

    Integer :: i, index, nl, rang, value
    Logical :: go

    ! set up sort

    If (n > 1) Then

      ! number of lists

      nl = n / 2

      ! iterate shell sort until there is a list

      Do While (nl > 0)

        ! for all lists from next-to-ground-level up to their end

        Do i = nl + 1, n

          value = list(i)
          rang = rank(i)
          index = i

          ! Antibubble down between levels of the same list

          go = .true.
          Do While (index > nl .and. go)
            go = (list(index - nl) > value)

            If (go) Then
              list(index) = list(index - nl)
              rank(index) = rank(index - nl)
              index = index - nl
            End If
          End Do

          ! Last insertion as close to the ground as it gets

          list(index) = value
          rank(index) = rang

        End Do

        ! Decrease the number of lists

        nl = nl / 2

      End Do

    End If

  End Subroutine shellsort2

  Function local_index(global_index, search_limit, rank, list)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 function to find the local atom number given the
    ! global atom number.  (For use with DD codes only)
    !
    ! If multiple copies present it returns the lowest local atom number
    ! If no copy is present it returns zero
    !
    ! rank(1,*) - array of local atom indices, ranking list
    ! list(1,*) - array of sorted global atom indices
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                 Intent(In   ) :: global_index, search_limit
    Integer, Dimension(1:*), Intent(In   ) :: rank, list
    Integer                                :: local_index

    Integer :: down, lower_bound, point, upper_bound

    ! Initialise

    local_index = 0

    ! Limits for the search

    lower_bound = 1
    upper_bound = search_limit

    ! Get smart, check for exceptions (whether it's a false pass)
    ! and check for a match on bounds

    If (global_index <= 0) Then

      local_index = 0

      Return

    Else If (global_index == list(lower_bound)) Then

      local_index = rank(lower_bound)

      Return

    Else If (global_index == list(search_limit)) Then

      down = search_limit

      Do While (global_index == list(down) .and. down >= lower_bound)

        local_index = rank(down)
        down = down - 1

      End Do

      Return

    End If

    ! Carry on then

    Do While (upper_bound - lower_bound > 1)

      point = (lower_bound + upper_bound) / 2

      If (global_index < list(point)) Then

        upper_bound = point

      Else If (global_index > list(point)) Then

        lower_bound = point

      Else

        down = point

        Do While (global_index == list(down))

          local_index = rank(down)
          down = down - 1

        End Do

        Return

      End If

    End Do

  End Function local_index

  Subroutine dcell(aaa, bbb)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to calculate the dimensional properties of a
    ! simulation cell specified by the input 3x3 matrix aaa (cell vectors in
    ! rows, the matrix is in the form of one dimensional reading
    ! (row1,row2,row3).
    !
    ! The results are returned in the array bbb, with:
    !
    ! bbb(1 to 3) - lengths of cell vectors: a(x,y,z) , b(x,y,z) , c(x,y,z)
    ! bbb(4 to 6) - cosines of cell angles: gamma(a,b) , beta(a,c) , alpha(b,c)
    ! bbb(7 to 9) - perpendicular cell widths : wx(y,z) , wy(x,z) , wz(x,y)
    ! bbb(10)     - cell volume
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov may 2008
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Dimension(1:9),  Intent(In   ) :: aaa
    Real(Kind=wp), Dimension(1:10), Intent(  Out) :: bbb

    Real(Kind=wp) :: axb1, axb2, axb3, bxc1, bxc2, bxc3, cxa1, cxa2, cxa3, d(1:3), x(1:3), y(1:3), &
                     z(1:3)

    ! calculate lengths of cell vectors

    bbb(1) = Sqrt(aaa(1)**2 + aaa(2)**2 + aaa(3)**2)
    bbb(2) = Sqrt(aaa(4)**2 + aaa(5)**2 + aaa(6)**2)
    bbb(3) = Sqrt(aaa(7)**2 + aaa(8)**2 + aaa(9)**2)

    ! calculate cosines of cell angles

    bbb(4) = (aaa(1) * aaa(4) + aaa(2) * aaa(5) + aaa(3) * aaa(6)) / (bbb(1) * bbb(2))
    bbb(5) = (aaa(1) * aaa(7) + aaa(2) * aaa(8) + aaa(3) * aaa(9)) / (bbb(1) * bbb(3))
    bbb(6) = (aaa(4) * aaa(7) + aaa(5) * aaa(8) + aaa(6) * aaa(9)) / (bbb(2) * bbb(3))

    ! calculate vector products of cell vectors

    axb1 = aaa(2) * aaa(6) - aaa(3) * aaa(5)
    axb2 = aaa(3) * aaa(4) - aaa(1) * aaa(6)
    axb3 = aaa(1) * aaa(5) - aaa(2) * aaa(4)

    bxc1 = aaa(5) * aaa(9) - aaa(6) * aaa(8)
    bxc2 = aaa(6) * aaa(7) - aaa(4) * aaa(9)
    bxc3 = aaa(4) * aaa(8) - aaa(5) * aaa(7)

    cxa1 = aaa(8) * aaa(3) - aaa(9) * aaa(2)
    cxa2 = aaa(9) * aaa(1) - aaa(7) * aaa(3)
    cxa3 = aaa(7) * aaa(2) - aaa(8) * aaa(1)

    ! calculate volume of cell

    bbb(10) = Abs(aaa(1) * bxc1 + aaa(2) * bxc2 + aaa(3) * bxc3)

    ! calculate cell perpendicular widths

    d(1) = bbb(10) / Sqrt(bxc1 * bxc1 + bxc2 * bxc2 + bxc3 * bxc3)
    d(2) = bbb(10) / Sqrt(cxa1 * cxa1 + cxa2 * cxa2 + cxa3 * cxa3)
    d(3) = bbb(10) / Sqrt(axb1 * axb1 + axb2 * axb2 + axb3 * axb3)

    x(1) = Abs(aaa(1)) / bbb(1); y(1) = Abs(aaa(2)) / bbb(1); z(1) = Abs(aaa(3)) / bbb(1)
    x(2) = Abs(aaa(4)) / bbb(2); y(2) = Abs(aaa(5)) / bbb(2); z(2) = Abs(aaa(6)) / bbb(2)
    x(3) = Abs(aaa(7)) / bbb(3); y(3) = Abs(aaa(8)) / bbb(3); z(3) = Abs(aaa(9)) / bbb(3)

    ! distribute widths

    If (x(1) >= x(2) .and. x(1) >= x(3)) Then
      bbb(7) = d(1)
      If (y(2) >= y(3)) Then
        bbb(8) = d(2)
        bbb(9) = d(3)
      Else
        bbb(8) = d(3)
        bbb(9) = d(2)
      End If
    Else If (x(2) >= x(1) .and. x(2) >= x(3)) Then
      bbb(7) = d(2)
      If (y(1) >= y(3)) Then
        bbb(8) = d(1)
        bbb(9) = d(3)
      Else
        bbb(8) = d(3)
        bbb(9) = d(1)
      End If
    Else
      bbb(7) = d(3)
      If (y(1) >= y(2)) Then
        bbb(8) = d(1)
        bbb(9) = d(2)
      Else
        bbb(8) = d(2)
        bbb(9) = d(1)
      End If
    End If

  End Subroutine dcell

  Subroutine invert(a, b, d)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to invert a 3x3 matrix using cofactors
    ! matrices are in the form of one dimensional array reading
    ! (row1,row2,row3)
    !
    ! a - input matrix
    ! b - inverted matrix
    ! d - determinant
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov august 2004
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Dimension(1:9), Intent(In   ) :: a
    Real(Kind=wp), Dimension(1:9), Intent(  Out) :: b
    Real(Kind=wp),                 Intent(  Out) :: d

    Real(Kind=wp) :: r

    ! calculate adjoint matrix

    b(1) = a(5) * a(9) - a(6) * a(8)
    b(2) = a(3) * a(8) - a(2) * a(9)
    b(3) = a(2) * a(6) - a(3) * a(5)
    b(4) = a(6) * a(7) - a(4) * a(9)
    b(5) = a(1) * a(9) - a(3) * a(7)
    b(6) = a(3) * a(4) - a(1) * a(6)
    b(7) = a(4) * a(8) - a(5) * a(7)
    b(8) = a(2) * a(7) - a(1) * a(8)
    b(9) = a(1) * a(5) - a(2) * a(4)

    ! calculate determinant

    d = a(1) * b(1) + a(4) * b(2) + a(7) * b(3)
    r = 0.0_wp
    If (Abs(d) > 0.0_wp) r = 1.0_wp / d

    ! complete inverse matrix

    b(1) = r * b(1)
    b(2) = r * b(2)
    b(3) = r * b(3)
    b(4) = r * b(4)
    b(5) = r * b(5)
    b(6) = r * b(6)
    b(7) = r * b(7)
    b(8) = r * b(8)
    b(9) = r * b(9)

  End Subroutine invert

  Subroutine images(imcon, cell, pairs, xxx, yyy, zzz)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image vector of
    ! atom pairs within a specified MD cell.  The cell matrix is in the form
    ! of one dimensional array reading (row1,row2,row3).
    !
    ! Image conditions
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov march 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                       Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9), Intent(In   ) :: cell
    Integer,                       Intent(In   ) :: pairs
    Real(Kind=wp), Dimension(1:*), Intent(InOut) :: xxx, yyy, zzz

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, det, rcell(1:9), xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      aaa = 1.0_wp / cell(1)

      Do i = 1, pairs
        xxx(i) = xxx(i) - cell(1) * Anint(aaa * xxx(i))
        yyy(i) = yyy(i) - cell(1) * Anint(aaa * yyy(i))
        zzz(i) = zzz(i) - cell(1) * Anint(aaa * zzz(i))
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then

      ! rectangular (slab) boundary conditions

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(5)
      ccc = 1.0_wp / cell(9)

      Do i = 1, pairs
        xxx(i) = xxx(i) - cell(1) * Anint(aaa * xxx(i))
        yyy(i) = yyy(i) - cell(5) * Anint(bbb * yyy(i))
        zzz(i) = zzz(i) - cell(9) * Anint(ccc * zzz(i))
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Call invert(cell, rcell, det)

      Do i = 1, pairs
        xss = rcell(1) * xxx(i) + rcell(4) * yyy(i) + rcell(7) * zzz(i)
        yss = rcell(2) * xxx(i) + rcell(5) * yyy(i) + rcell(8) * zzz(i)
        zss = rcell(3) * xxx(i) + rcell(6) * yyy(i) + rcell(9) * zzz(i)

        xss = xss - Anint(xss)
        yss = yss - Anint(yss)
        zss = zss - Anint(zss)

        xxx(i) = cell(1) * xss + cell(4) * yss + cell(7) * zss
        yyy(i) = cell(2) * xss + cell(5) * yss + cell(8) * zss
        zzz(i) = cell(3) * xss + cell(6) * yss + cell(9) * zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, pairs
        xxx(i) = xxx(i) - cell(1) * Anint(aaa * xxx(i))
        yyy(i) = yyy(i) - cell(1) * Anint(aaa * yyy(i))
        zzz(i) = zzz(i) - cell(1) * Anint(aaa * zzz(i))

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(zzz(i))) >= 0.75_wp * cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(1), zzz(i))
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, pairs
        xxx(i) = xxx(i) - cell(1) * Anint(aaa * xxx(i))
        yyy(i) = yyy(i) - cell(1) * Anint(aaa * yyy(i))
        zzz(i) = zzz(i) - cell(9) * Anint(bbb * zzz(i))

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(rt2 * zzz(i))) >= cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(9), zzz(i))
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions

      det = cell(1) * cell(5) - cell(2) * cell(4)

      If (Abs(det) < 1.0e-6_wp) Call error(120)

      det = 1.0_wp / det

      rcell(1) = det * cell(5)
      rcell(2) = -det * cell(2)
      rcell(4) = -det * cell(4)
      rcell(5) = det * cell(1)

      Do i = 1, pairs
        xss = rcell(1) * xxx(i) + rcell(4) * yyy(i)
        yss = rcell(2) * xxx(i) + rcell(5) * yyy(i)

        xss = xss - Anint(xss)
        yss = yss - Anint(yss)

        xxx(i) = cell(1) * xss + cell(4) * yss
        yyy(i) = cell(2) * xss + cell(5) * yss
      End Do

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, pairs
        yyy(i) = yyy(i) - bbb * Anint(ccc * yyy(i))
        zzz(i) = zzz(i) - cell(9) * Anint(ddd * zzz(i))

        If ((Abs(yyy(i)) + Abs(rt3 * xxx(i))) >= bbb) Then
          xxx(i) = xxx(i) - rt3 * Sign(aaa, xxx(i))
          yyy(i) = yyy(i) - Sign(aaa, yyy(i))
        End If
      End Do

    End If

  End Subroutine images

  Subroutine images_s(imcon, cell, xxx, yyy, zzz)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image vector of an
    ! atom pair within a specified MD cell.  The cell matrix is in the form
    ! of one dimensional array reading (row1,row2,row3).
    !
    ! Image conditions
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith july 1992
    ! amended   - i.t.todorov march 2015
    ! tweaked   - h.a.boateng january 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                       Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9), Intent(In   ) :: cell
    Real(Kind=wp),                 Intent(InOut) :: xxx, yyy, zzz

    Real(Kind=wp) :: aaa, bbb, ccc, ddd, det, rcell(1:9), xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      aaa = 1.0_wp / cell(1)

      xxx = xxx - cell(1) * Anint(aaa * xxx)
      yyy = yyy - cell(1) * Anint(aaa * yyy)
      zzz = zzz - cell(1) * Anint(aaa * zzz)

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular (slab) boundary conditions

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(5)
      ccc = 1.0_wp / cell(9)
      xxx = xxx - cell(1) * Anint(aaa * xxx)
      yyy = yyy - cell(5) * Anint(bbb * yyy)
      zzz = zzz - cell(9) * Anint(ccc * zzz)

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Call invert(cell, rcell, det)

      xss = rcell(1) * xxx + rcell(4) * yyy + rcell(7) * zzz
      yss = rcell(2) * xxx + rcell(5) * yyy + rcell(8) * zzz
      zss = rcell(3) * xxx + rcell(6) * yyy + rcell(9) * zzz

      xss = xss - Anint(xss)
      yss = yss - Anint(yss)
      zss = zss - Anint(zss)

      xxx = cell(1) * xss + cell(4) * yss + cell(7) * zss
      yyy = cell(2) * xss + cell(5) * yss + cell(8) * zss
      zzz = cell(3) * xss + cell(6) * yss + cell(9) * zss

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      xxx = xxx - cell(1) * Anint(aaa * xxx)
      yyy = yyy - cell(1) * Anint(aaa * yyy)
      zzz = zzz - cell(1) * Anint(aaa * zzz)

      If ((Abs(xxx) + Abs(yyy) + Abs(zzz)) >= 0.75_wp * cell(1)) Then
        xxx = xxx - 0.5_wp * Sign(cell(1), xxx)
        yyy = yyy - 0.5_wp * Sign(cell(1), yyy)
        zzz = zzz - 0.5_wp * Sign(cell(1), zzz)
      End If

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      xxx = xxx - cell(1) * Anint(aaa * xxx)
      yyy = yyy - cell(1) * Anint(aaa * yyy)
      zzz = zzz - cell(9) * Anint(bbb * zzz)

      If ((Abs(xxx) + Abs(yyy) + Abs(rt2 * zzz)) >= cell(1)) Then
        xxx = xxx - 0.5_wp * Sign(cell(1), xxx)
        yyy = yyy - 0.5_wp * Sign(cell(1), yyy)
        zzz = zzz - 0.5_wp * Sign(cell(9), zzz)
      End If

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions

      det = cell(1) * cell(5) - cell(2) * cell(4)

      If (Abs(det) < 1.0e-6_wp) Call error(120)

      det = 1.0_wp / det

      rcell(1) = det * cell(5)
      rcell(2) = -det * cell(2)
      rcell(4) = -det * cell(4)
      rcell(5) = det * cell(1)

      xss = rcell(1) * xxx + rcell(4) * yyy
      yss = rcell(2) * xxx + rcell(5) * yyy

      xss = xss - Anint(xss)
      yss = yss - Anint(yss)

      xxx = cell(1) * xss + cell(4) * yss
      yyy = cell(2) * xss + cell(5) * yss

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      yyy = yyy - bbb * Anint(ccc * yyy)
      zzz = zzz - cell(9) * Anint(ddd * zzz)

      If ((Abs(yyy) + Abs(rt3 * xxx)) >= bbb) Then
        xxx = xxx - rt3 * Sign(aaa, xxx)
        yyy = yyy - Sign(aaa, yyy)
      End If

    End If

  End Subroutine images_s

  Subroutine pbcshift_parts(imcon, cell, natms, parts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image of atoms within
    ! a specified MD cell in accordance with the domain decomposition
    ! boundary convention for fractional coordinates: every coordinate must
    ! be intervalled as [-0.5,+0.5)
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                        Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9),  Intent(In   ) :: cell
    Integer,                        Intent(In   ) :: natms
    Type(corePart), Dimension(1:*), Intent(InOut) :: parts

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, det, rcell(1:9), xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = aaa * parts(i)%yyy
        zss = aaa * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(1) * zss
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular boundary conditions

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(5)
      ccc = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = bbb * parts(i)%yyy
        zss = ccc * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(5) * yss
        parts(i)%zzz = cell(9) * zss
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Call invert(cell, rcell, det)

      Do i = 1, natms
        xss = rcell(1) * parts(i)%xxx + rcell(4) * parts(i)%yyy + rcell(7) * parts(i)%zzz
        yss = rcell(2) * parts(i)%xxx + rcell(5) * parts(i)%yyy + rcell(8) * parts(i)%zzz
        zss = rcell(3) * parts(i)%xxx + rcell(6) * parts(i)%yyy + rcell(9) * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss + cell(4) * yss + cell(7) * zss
        parts(i)%yyy = cell(2) * xss + cell(5) * yss + cell(8) * zss
        parts(i)%zzz = cell(3) * xss + cell(6) * yss + cell(9) * zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = aaa * parts(i)%yyy
        zss = aaa * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(1) * zss

        If ((Abs(parts(i)%xxx) + Abs(parts(i)%yyy) + Abs(parts(i)%zzz)) >= 0.75_wp * cell(1)) Then
          parts(i)%xxx = parts(i)%xxx - 0.5_wp * Sign(cell(1), parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - 0.5_wp * Sign(cell(1), parts(i)%yyy)
          parts(i)%zzz = parts(i)%zzz - 0.5_wp * Sign(cell(1), parts(i)%zzz)

          xss = aaa * parts(i)%xxx
          yss = aaa * parts(i)%yyy
          zss = aaa * parts(i)%zzz

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          parts(i)%xxx = cell(1) * xss
          parts(i)%yyy = cell(1) * yss
          parts(i)%zzz = cell(1) * zss
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = aaa * parts(i)%yyy
        zss = bbb * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(9) * zss

        If ((Abs(parts(i)%xxx) + Abs(parts(i)%yyy) + Abs(rt2 * parts(i)%zzz)) >= cell(1)) Then
          parts(i)%xxx = parts(i)%xxx - 0.5_wp * Sign(cell(1), parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - 0.5_wp * Sign(cell(1), parts(i)%yyy)
          parts(i)%zzz = parts(i)%zzz - 0.5_wp * Sign(cell(9), parts(i)%zzz)

          xss = aaa * parts(i)%xxx
          yss = aaa * parts(i)%yyy
          zss = bbb * parts(i)%zzz

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          parts(i)%xxx = cell(1) * xss
          parts(i)%yyy = cell(1) * yss
          parts(i)%zzz = cell(9) * zss
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions (SLAB)

      det = cell(1) * cell(5) - cell(2) * cell(4)

      If (Abs(det) < 1.0e-6_wp) Call error(120)

      det = 1.0_wp / det

      rcell(1) = det * cell(5)
      rcell(2) = -det * cell(2)
      rcell(4) = -det * cell(4)
      rcell(5) = det * cell(1)

      Do i = 1, natms
        xss = rcell(1) * parts(i)%xxx + rcell(4) * parts(i)%yyy
        yss = rcell(2) * parts(i)%xxx + rcell(5) * parts(i)%yyy

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss

        parts(i)%xxx = cell(1) * xss + cell(4) * yss
        parts(i)%yyy = cell(2) * xss + cell(5) * yss
      End Do

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, natms
        zss = ddd * parts(i)%zzz
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss
        parts(i)%zzz = cell(9) * zss

        yss = ccc * parts(i)%yyy
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        parts(i)%yyy = bbb * yss

        If ((Abs(parts(i)%yyy) + Abs(rt3 * parts(i)%xxx)) >= bbb) Then
          parts(i)%xxx = parts(i)%xxx - rt3 * Sign(aaa, parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - Sign(aaa, parts(i)%yyy)

          yss = ccc * parts(i)%yyy
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          parts(i)%yyy = bbb * yss
        End If
      End Do

    End If

  End Subroutine pbcshift_parts

  Subroutine pbcshift_arrays(imcon, cell, natms, xxx, yyy, zzz)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image of atoms within
    ! a specified MD cell in accordance with the domain decomposition
    ! boundary convention for fractional coordinates: every coordinate must
    ! be intervalled as [-0.5,+0.5)
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2015
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                       Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9), Intent(In   ) :: cell
    Integer,                       Intent(In   ) :: natms
    Real(Kind=wp), Dimension(1:*), Intent(InOut) :: xxx, yyy, zzz

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, det, rcell(1:9), xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = aaa * yyy(i)
        zss = aaa * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(1) * zss
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular boundary conditions

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(5)
      ccc = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = bbb * yyy(i)
        zss = ccc * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(5) * yss
        zzz(i) = cell(9) * zss
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Call invert(cell, rcell, det)

      Do i = 1, natms
        xss = rcell(1) * xxx(i) + rcell(4) * yyy(i) + rcell(7) * zzz(i)
        yss = rcell(2) * xxx(i) + rcell(5) * yyy(i) + rcell(8) * zzz(i)
        zss = rcell(3) * xxx(i) + rcell(6) * yyy(i) + rcell(9) * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss + cell(4) * yss + cell(7) * zss
        yyy(i) = cell(2) * xss + cell(5) * yss + cell(8) * zss
        zzz(i) = cell(3) * xss + cell(6) * yss + cell(9) * zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = aaa * yyy(i)
        zss = aaa * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(1) * zss

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(zzz(i))) >= 0.75_wp * cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(1), zzz(i))

          xss = aaa * xxx(i)
          yss = aaa * yyy(i)
          zss = aaa * zzz(i)

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          xxx(i) = cell(1) * xss
          yyy(i) = cell(1) * yss
          zzz(i) = cell(1) * zss
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = aaa * yyy(i)
        zss = bbb * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(9) * zss

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(rt2 * zzz(i))) >= cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(9), zzz(i))

          xss = aaa * xxx(i)
          yss = aaa * yyy(i)
          zss = bbb * zzz(i)

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          xxx(i) = cell(1) * xss
          yyy(i) = cell(1) * yss
          zzz(i) = cell(9) * zss
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions (SLAB)

      det = cell(1) * cell(5) - cell(2) * cell(4)

      If (Abs(det) < 1.0e-6_wp) Call error(120)

      det = 1.0_wp / det

      rcell(1) = det * cell(5)
      rcell(2) = -det * cell(2)
      rcell(4) = -det * cell(4)
      rcell(5) = det * cell(1)

      Do i = 1, natms
        xss = rcell(1) * xxx(i) + rcell(4) * yyy(i)
        yss = rcell(2) * xxx(i) + rcell(5) * yyy(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss

        xxx(i) = cell(1) * xss + cell(4) * yss
        yyy(i) = cell(2) * xss + cell(5) * yss
      End Do

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, natms
        zss = ddd * zzz(i)
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss
        zzz(i) = cell(9) * zss

        yss = ccc * yyy(i)
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        yyy(i) = bbb * yss

        If ((Abs(yyy(i)) + Abs(rt3 * xxx(i))) >= bbb) Then
          xxx(i) = xxx(i) - rt3 * Sign(aaa, xxx(i))
          yyy(i) = yyy(i) - Sign(aaa, yyy(i))

          yss = ccc * yyy(i)
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          yyy(i) = bbb * yss
        End If
      End Do

    End If

  End Subroutine pbcshift_arrays

  Subroutine pbcshfrc_parts(imcon, cell, natms, parts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image of atoms within
    ! a specified MD cell in accordance with the domain decomposition
    ! boundary convention for fractional coordinates: every coordinate must
    ! be intervalled as [-0.5,+0.5)
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                        Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9),  Intent(In   ) :: cell
    Integer,                        Intent(In   ) :: natms
    Type(corePart), Dimension(1:*), Intent(InOut) :: parts

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, det, rcell(1:9), xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = aaa * parts(i)%yyy
        zss = aaa * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = xss
        parts(i)%yyy = yss
        parts(i)%zzz = zss
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular boundary conditions

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(5)
      ccc = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = bbb * parts(i)%yyy
        zss = ccc * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = xss
        parts(i)%yyy = yss
        parts(i)%zzz = zss
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Call invert(cell, rcell, det)

      Do i = 1, natms
        xss = rcell(1) * parts(i)%xxx + rcell(4) * parts(i)%yyy + rcell(7) * parts(i)%zzz
        yss = rcell(2) * parts(i)%xxx + rcell(5) * parts(i)%yyy + rcell(8) * parts(i)%zzz
        zss = rcell(3) * parts(i)%xxx + rcell(6) * parts(i)%yyy + rcell(9) * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = xss
        parts(i)%yyy = yss
        parts(i)%zzz = zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = aaa * parts(i)%yyy
        zss = aaa * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(1) * zss

        If ((Abs(parts(i)%xxx) + Abs(parts(i)%yyy) + Abs(parts(i)%zzz)) >= 0.75_wp * cell(1)) Then
          parts(i)%xxx = parts(i)%xxx - 0.5_wp * Sign(cell(1), parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - 0.5_wp * Sign(cell(1), parts(i)%yyy)
          parts(i)%zzz = parts(i)%zzz - 0.5_wp * Sign(cell(1), parts(i)%zzz)

          xss = aaa * parts(i)%xxx
          yss = aaa * parts(i)%yyy
          zss = aaa * parts(i)%zzz

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          parts(i)%xxx = xss
          parts(i)%yyy = yss
          parts(i)%zzz = zss
        Else
          parts(i)%xxx = xss / cell(1)
          parts(i)%yyy = yss / cell(1)
          parts(i)%zzz = zss / cell(1)
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * parts(i)%xxx
        yss = aaa * parts(i)%yyy
        zss = bbb * parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(9) * zss

        If ((Abs(parts(i)%xxx) + Abs(parts(i)%yyy) + Abs(rt2 * parts(i)%zzz)) >= cell(1)) Then
          parts(i)%xxx = parts(i)%xxx - 0.5_wp * Sign(cell(1), parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - 0.5_wp * Sign(cell(1), parts(i)%yyy)
          parts(i)%zzz = parts(i)%zzz - 0.5_wp * Sign(cell(9), parts(i)%zzz)

          xss = aaa * parts(i)%xxx
          yss = aaa * parts(i)%yyy
          zss = bbb * parts(i)%zzz

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          parts(i)%xxx = xss
          parts(i)%yyy = yss
          parts(i)%zzz = zss
        Else
          parts(i)%xxx = xss / cell(1)
          parts(i)%yyy = yss / cell(1)
          parts(i)%zzz = zss / cell(9)
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions (SLAB)

      det = cell(1) * cell(5) - cell(2) * cell(4)

      If (Abs(det) < 1.0e-6_wp) Call error(120)

      det = 1.0_wp / det

      rcell(1) = det * cell(5)
      rcell(2) = -det * cell(2)
      rcell(4) = -det * cell(4)
      rcell(5) = det * cell(1)

      Do i = 1, natms
        xss = rcell(1) * parts(i)%xxx + rcell(4) * parts(i)%yyy
        yss = rcell(2) * parts(i)%xxx + rcell(5) * parts(i)%yyy

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss

        parts(i)%xxx = xss
        parts(i)%yyy = yss
      End Do ! note zzz remains in real space

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, natms
        zss = ddd * parts(i)%zzz
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss
        parts(i)%zzz = zss

        yss = ccc * parts(i)%yyy
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        parts(i)%yyy = bbb * yss

        If ((Abs(parts(i)%yyy) + Abs(rt3 * parts(i)%xxx)) >= bbb) Then
          parts(i)%xxx = parts(i)%xxx - rt3 * Sign(aaa, parts(i)%xxx)
          parts(i)%xxx = parts(i)%xxx / aaa
          parts(i)%yyy = parts(i)%yyy - Sign(aaa, parts(i)%yyy)

          yss = ccc * parts(i)%yyy
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          parts(i)%yyy = yss
        Else
          parts(i)%xxx = parts(i)%xxx / aaa
          parts(i)%yyy = parts(i)%yyy / bbb
        End If
      End Do

    End If

  End Subroutine pbcshfrc_parts

  Subroutine pbcshfrc_arrays(imcon, cell, natms, xxx, yyy, zzz)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image of atoms within
    ! a specified MD cell in accordance with the domain decomposition
    ! boundary convention for fractional coordinates: every coordinate must
    ! be intervalled as [-0.5,+0.5)
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                       Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9), Intent(In   ) :: cell
    Integer,                       Intent(In   ) :: natms
    Real(Kind=wp), Dimension(1:*), Intent(InOut) :: xxx, yyy, zzz

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, det, rcell(1:9), xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = aaa * yyy(i)
        zss = aaa * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = xss
        yyy(i) = yss
        zzz(i) = zss
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular boundary conditions

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(5)
      ccc = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = bbb * yyy(i)
        zss = ccc * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = xss
        yyy(i) = yss
        zzz(i) = zss
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Call invert(cell, rcell, det)

      Do i = 1, natms
        xss = rcell(1) * xxx(i) + rcell(4) * yyy(i) + rcell(7) * zzz(i)
        yss = rcell(2) * xxx(i) + rcell(5) * yyy(i) + rcell(8) * zzz(i)
        zss = rcell(3) * xxx(i) + rcell(6) * yyy(i) + rcell(9) * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = xss
        yyy(i) = yss
        zzz(i) = zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = aaa * yyy(i)
        zss = aaa * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(1) * zss

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(zzz(i))) >= 0.75_wp * cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(1), zzz(i))

          xss = aaa * xxx(i)
          yss = aaa * yyy(i)
          zss = aaa * zzz(i)

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          xxx(i) = xss
          yyy(i) = yss
          zzz(i) = zss
        Else
          xxx(i) = xss / cell(1)
          yyy(i) = yss / cell(1)
          zzz(i) = zss / cell(1)
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = aaa * xxx(i)
        yss = aaa * yyy(i)
        zss = bbb * zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(9) * zss

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(rt2 * zzz(i))) >= cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(9), zzz(i))

          xss = aaa * xxx(i)
          yss = aaa * yyy(i)
          zss = bbb * zzz(i)

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          xxx(i) = xss
          yyy(i) = yss
          zzz(i) = zss
        Else
          xxx(i) = xss / cell(1)
          yyy(i) = yss / cell(1)
          zzz(i) = zss / cell(9)
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions (SLAB)

      det = cell(1) * cell(5) - cell(2) * cell(4)

      If (Abs(det) < 1.0e-6_wp) Call error(120)

      det = 1.0_wp / det

      rcell(1) = det * cell(5)
      rcell(2) = -det * cell(2)
      rcell(4) = -det * cell(4)
      rcell(5) = det * cell(1)

      Do i = 1, natms
        xss = rcell(1) * xxx(i) + rcell(4) * yyy(i)
        yss = rcell(2) * xxx(i) + rcell(5) * yyy(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss

        xxx(i) = xss
        yyy(i) = yss
      End Do ! note zzz remains in real space

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, natms
        zss = ddd * zzz(i)
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss
        zzz(i) = zss

        yss = ccc * yyy(i)
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        yyy(i) = bbb * yss

        If ((Abs(yyy(i)) + Abs(rt3 * xxx(i))) >= bbb) Then
          xxx(i) = xxx(i) - rt3 * Sign(aaa, xxx(i))
          xxx(i) = xxx(i) / aaa
          yyy(i) = yyy(i) - Sign(aaa, yyy(i))

          yss = ccc * yyy(i)
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          yyy(i) = yss
        Else
          xxx(i) = xxx(i) / aaa
          yyy(i) = yyy(i) / bbb
        End If
      End Do

    End If

  End Subroutine pbcshfrc_arrays

  Subroutine pbcshfrl_parts(imcon, cell, natms, parts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image of atoms within
    ! a specified MD cell in accordance with the domain decomposition
    ! boundary convention for fractional coordinates: every coordinate must
    ! be intervalled as [-0.5,+0.5)
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                        Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9),  Intent(In   ) :: cell
    Integer,                        Intent(In   ) :: natms
    Type(corePart), Dimension(1:*), Intent(InOut) :: parts

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      Do i = 1, natms
        xss = parts(i)%xxx
        yss = parts(i)%yyy
        zss = parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(1) * zss
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular boundary conditions

      Do i = 1, natms
        xss = parts(i)%xxx
        yss = parts(i)%yyy
        zss = parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(5) * yss
        parts(i)%zzz = cell(9) * zss
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Do i = 1, natms
        xss = parts(i)%xxx
        yss = parts(i)%yyy
        zss = parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss + cell(4) * yss + cell(7) * zss
        parts(i)%yyy = cell(2) * xss + cell(5) * yss + cell(8) * zss
        parts(i)%zzz = cell(3) * xss + cell(6) * yss + cell(9) * zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = parts(i)%xxx
        yss = parts(i)%yyy
        zss = parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(1) * zss

        If ((Abs(parts(i)%xxx) + Abs(parts(i)%yyy) + Abs(parts(i)%zzz)) >= 0.75_wp * cell(1)) Then
          parts(i)%xxx = parts(i)%xxx - 0.5_wp * Sign(cell(1), parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - 0.5_wp * Sign(cell(1), parts(i)%yyy)
          parts(i)%zzz = parts(i)%zzz - 0.5_wp * Sign(cell(1), parts(i)%zzz)

          xss = aaa * parts(i)%xxx
          yss = aaa * parts(i)%yyy
          zss = aaa * parts(i)%zzz

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          parts(i)%xxx = cell(1) * xss
          parts(i)%yyy = cell(1) * yss
          parts(i)%zzz = cell(1) * zss
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = parts(i)%xxx
        yss = parts(i)%yyy
        zss = parts(i)%zzz

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        parts(i)%xxx = cell(1) * xss
        parts(i)%yyy = cell(1) * yss
        parts(i)%zzz = cell(9) * zss

        If ((Abs(parts(i)%xxx) + Abs(parts(i)%yyy) + Abs(rt2 * parts(i)%zzz)) >= cell(1)) Then
          parts(i)%xxx = parts(i)%xxx - 0.5_wp * Sign(cell(1), parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - 0.5_wp * Sign(cell(1), parts(i)%yyy)
          parts(i)%zzz = parts(i)%zzz - 0.5_wp * Sign(cell(9), parts(i)%zzz)

          xss = aaa * parts(i)%xxx
          yss = aaa * parts(i)%yyy
          zss = bbb * parts(i)%zzz

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          parts(i)%xxx = cell(1) * xss
          parts(i)%yyy = cell(1) * yss
          parts(i)%zzz = cell(9) * zss
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions (SLAB)

      Do i = 1, natms
        xss = parts(i)%xxx
        yss = parts(i)%yyy

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss

        parts(i)%xxx = cell(1) * xss + cell(4) * yss
        parts(i)%yyy = cell(2) * xss + cell(5) * yss
      End Do ! note zzz remains unchanged

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, natms
        zss = parts(i)%zzz
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss
        parts(i)%zzz = cell(9) * zss

        yss = parts(i)%yyy
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        parts(i)%yyy = bbb * yss
        parts(i)%xxx = aaa * parts(i)%xxx

        If ((Abs(parts(i)%yyy) + Abs(rt3 * parts(i)%xxx)) >= bbb) Then
          parts(i)%xxx = parts(i)%xxx - rt3 * Sign(aaa, parts(i)%xxx)
          parts(i)%yyy = parts(i)%yyy - Sign(aaa, parts(i)%yyy)

          yss = ccc * parts(i)%yyy
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          parts(i)%yyy = bbb * yss
        End If
      End Do

    End If

  End Subroutine pbcshfrl_parts

  Subroutine pbcshfrl_arrays(imcon, cell, natms, xxx, yyy, zzz)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the minimum image of atoms within
    ! a specified MD cell in accordance with the domain decomposition
    ! boundary convention for fractional coordinates: every coordinate must
    ! be intervalled as [-0.5,+0.5)
    !
    ! Note: in all cases the centre of the MD cell is at (0,0,0)
    !
    ! imcon=0 no boundary conditions apply
    ! imcon=1 standard cubic boundaries apply
    ! imcon=2 orthorhombic boundaries apply
    ! imcon=3 parallelepiped boundaries apply
    ! imcon=4 truncated octahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=5 rhombic dodecahedron boundaries apply NOT AVAILABLE in DD !!!
    ! imcon=6 x-y parallelogram boundary conditions : no periodicity in z
    ! imcon=7 hexagonal prism boundaries apply      NOT AVAILABLE in DD !!!
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                       Intent(In   ) :: imcon
    Real(Kind=wp), Dimension(1:9), Intent(In   ) :: cell
    Integer,                       Intent(In   ) :: natms
    Real(Kind=wp), Dimension(1:*), Intent(InOut) :: xxx, yyy, zzz

    Integer       :: i
    Real(Kind=wp) :: aaa, bbb, ccc, ddd, xss, yss, zss

    If (imcon == IMCON_CUBIC) Then

      ! standard cubic boundary conditions

      Do i = 1, natms
        xss = xxx(i)
        yss = yyy(i)
        zss = zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(1) * zss
      End Do

    Else If (imcon == IMCON_ORTHORHOMBIC .or. imcon == IMCON_NOPBC) Then ! no PBC box wrapping exception

      ! rectangular boundary conditions

      Do i = 1, natms
        xss = xxx(i)
        yss = yyy(i)
        zss = zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(5) * yss
        zzz(i) = cell(9) * zss
      End Do

    Else If (imcon == IMCON_PARALLELOPIPED) Then

      ! parallelepiped boundary conditions

      Do i = 1, natms
        xss = xxx(i)
        yss = yyy(i)
        zss = zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss + cell(4) * yss + cell(7) * zss
        yyy(i) = cell(2) * xss + cell(5) * yss + cell(8) * zss
        zzz(i) = cell(3) * xss + cell(6) * yss + cell(9) * zss
      End Do

    Else If (imcon == IMCON_TRUNC_OCTO) Then

      ! truncated octahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(5) - cell(9)) < 1.0e-6_wp)) Call error(130)

      aaa = 1.0_wp / cell(1)

      Do i = 1, natms
        xss = xxx(i)
        yss = yyy(i)
        zss = zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(1) * zss

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(zzz(i))) >= 0.75_wp * cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(1), zzz(i))

          xss = aaa * xxx(i)
          yss = aaa * yyy(i)
          zss = aaa * zzz(i)

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          xxx(i) = cell(1) * xss
          yyy(i) = cell(1) * yss
          zzz(i) = cell(1) * zss
        End If
      End Do

    Else If (imcon == IMCON_RHOMBIC_DODEC) Then

      ! rhombic Dodecahedral boundary conditions

      If (.not. (Abs(cell(1) - cell(5)) < 1.0e-6_wp .and. Abs(cell(9) - cell(1) * rt2) < 1.0e-6_wp)) Call error(140)

      aaa = 1.0_wp / cell(1)
      bbb = 1.0_wp / cell(9)

      Do i = 1, natms
        xss = xxx(i)
        yss = yyy(i)
        zss = zzz(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

        xxx(i) = cell(1) * xss
        yyy(i) = cell(1) * yss
        zzz(i) = cell(9) * zss

        If ((Abs(xxx(i)) + Abs(yyy(i)) + Abs(rt2 * zzz(i))) >= cell(1)) Then
          xxx(i) = xxx(i) - 0.5_wp * Sign(cell(1), xxx(i))
          yyy(i) = yyy(i) - 0.5_wp * Sign(cell(1), yyy(i))
          zzz(i) = zzz(i) - 0.5_wp * Sign(cell(9), zzz(i))

          xss = aaa * xxx(i)
          yss = aaa * yyy(i)
          zss = bbb * zzz(i)

          xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss

          xxx(i) = cell(1) * xss
          yyy(i) = cell(1) * yss
          zzz(i) = cell(9) * zss
        End If
      End Do

    Else If (imcon == IMCON_SLAB) Then

      ! x-y boundary conditions (SLAB)

      Do i = 1, natms
        xss = xxx(i)
        yss = yyy(i)

        xss = xss - Anint(xss); If (xss >= half_minus) xss = -xss
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss

        xxx(i) = cell(1) * xss + cell(4) * yss
        yyy(i) = cell(2) * xss + cell(5) * yss
      End Do ! note zzz remains unchanged

    Else If (imcon == IMCON_HEXAGONAL) Then

      ! hexagonal prism boundary conditions

      If (Abs(cell(1) - rt3 * cell(5)) > 1.0e-6_wp) Call error(135)

      aaa = cell(1) / (rt3 * 2.0_wp)
      bbb = cell(1) / rt3
      ccc = rt3 / cell(1)
      ddd = 1.0_wp / cell(9)

      Do i = 1, natms
        zss = zzz(i)
        zss = zss - Anint(zss); If (zss >= half_minus) zss = -zss
        zzz(i) = cell(9) * zss

        yss = yyy(i)
        yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
        yyy(i) = bbb * yss
        xxx(i) = aaa * xxx(i)

        If ((Abs(yyy(i)) + Abs(rt3 * xxx(i))) >= bbb) Then
          xxx(i) = xxx(i) - rt3 * Sign(aaa, xxx(i))
          yyy(i) = yyy(i) - Sign(aaa, yyy(i))

          yss = ccc * yyy(i)
          yss = yss - Anint(yss); If (yss >= half_minus) yss = -yss
          yyy(i) = bbb * yss
        End If
      End Do

    End If

  End Subroutine pbcshfrl_arrays

  Subroutine jacobi(n, aaa, vvv)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Diagonalisation of real square symmetric matrices by the Jacobi method:
    ! a sequence of Jacobi rotations
    !
    ! Users must ensure the symmetry of the input matrix
    !
    ! input parameters: n   - matrix dimension
    !                   aaa - the matrix to be diagonalised
    !                   vvv - the (diagonalised) eigenvector matrix
    !
    ! Jacobi processes lower triangle only - strictly upper triangle
    !                                        remains unchanged
    !
    ! Variable rho sets absolute tolerance on convergence
    ! Variable test is a moving tolerance that diminishes on each pass
    ! until true convergence test<rho
    ! ZERO matrices are accepted and returned
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov july 2008
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                            Intent(In   ) :: n
    Real(Kind=wp), Dimension(1:n, 1:n), Intent(InOut) :: aaa, vvv

    Real(Kind=wp), Parameter :: rho = 1.0e-20_wp

    Integer       :: i, j, k
    Logical       :: pass
    Real(Kind=wp) :: c, c_c, omg, s, s_c, s_s, scale, test, tmp, v_d_hor, v_d_mid, v_d_off, v_d_ver

! ,l

    !  l=0 ! Iteration counter

    ! Rescale (lower triangle) matrix for optimal accuracy
    ! by the largest by magnitude diagonal element

    scale = 0.0_wp
    Do i = 1, n
      If (Abs(aaa(i, i)) > scale) scale = Abs(aaa(i, i))
    End Do
    If (scale <= zero_plus) Then
      vvv = aaa
      Return ! Accept & Return zero matrices
    Else
      Do i = 1, n
        Do j = 1, i
          aaa(j, i) = aaa(j, i) / scale
        End Do
      End Do
    End If

    ! Set initial value of moving tolerance
    ! Sum of all off-diagonal elements (strictly lower triangle)

    test = 0.0_wp
    Do i = 2, n
      Do j = 1, i - 1
        test = test + aaa(j, i)**2
      End Do
    End Do
    test = Sqrt(2.0_wp * test)

    ! Initialise eigenvectors

    vvv = 0.0_wp
    Do i = 1, n
      vvv(i, i) = 1.0_wp
    End Do

    ! Accept & Return already diagonalised matrices
    ! (as well as zero matrices)

    If (test < rho) Return

    ! Recycle until absolute tolerance satisfied

    Do While (test > rho)
      test = test / Real(n, wp)
      If (test < rho) test = rho

      ! Jacobi diagonalisation

      pass = .true.

      ! Recycle until moving tolerance satisfied

      Do While (pass)
        pass = .false.

        ! Loop around the strictly lower triangle matrix

        Do i = 2, n
          Do j = 1, i - 1
            If (Abs(aaa(j, i)) >= test) Then
              !                 l=l+1
              pass = .true.

              v_d_hor = aaa(i, i)
              v_d_ver = aaa(j, j)
              v_d_off = aaa(j, i)
              v_d_mid = 0.5_wp * (v_d_ver - v_d_hor)
              If (Abs(v_d_mid) < rho) Then
                omg = -1.0_wp
              Else
                omg = -v_d_off / Sqrt(v_d_off**2 + v_d_mid**2)
                If (v_d_mid < 0.0_wp) omg = -omg
              End If
              s = omg / Sqrt(2.0_wp * (1.0_wp + Sqrt(1.0_wp - omg**2)))
              s_s = s * s; c_c = 1.0_wp - s_s; c = Sqrt(c_c); s_c = s * c

              Do k = 1, n
                If (k <= j) Then
                  tmp = aaa(k, j) * c - aaa(k, i) * s
                  aaa(k, i) = aaa(k, j) * s + aaa(k, i) * c
                  aaa(k, j) = tmp
                Else If (k > i) Then
                  tmp = aaa(j, k) * c - aaa(i, k) * s
                  aaa(i, k) = aaa(j, k) * s + aaa(i, k) * c
                  aaa(j, k) = tmp
                Else
                  tmp = aaa(j, k) * c - aaa(k, i) * s
                  aaa(k, i) = aaa(j, k) * s + aaa(k, i) * c
                  aaa(j, k) = tmp
                End If

                tmp = vvv(k, j) * c - vvv(k, i) * s
                vvv(k, i) = vvv(k, j) * s + vvv(k, i) * c
                vvv(k, j) = tmp
              End Do

              aaa(i, i) = v_d_hor * c_c + v_d_ver * s_s + 2.0_wp * v_d_off * s_c
              aaa(j, j) = v_d_hor * s_s + v_d_ver * c_c - 2.0_wp * v_d_off * s_c
              aaa(j, i) = (v_d_ver - v_d_hor) * s_c + v_d_off * (c_c - s_s)
            End If
          End Do
        End Do
      End Do
    End Do

    ! Rescale back the lower triangle matrix

    Do i = 1, n
      Do j = 1, i
        aaa(j, i) = aaa(j, i) * scale
      End Do
    End Do

  End Subroutine jacobi

  Subroutine mat_mul(aaa, bbb, ccc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! matrix multiply routine: A*B=C (note order!)
    !
    ! Note: A, B and C are 3x3 matrices in linear arrays as used in dl_poly
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith april 2009
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(In   ) :: aaa(1:9), bbb(1:9)
    Real(Kind=wp), Intent(  Out) :: ccc(1:9)

    ccc(1) = aaa(1) * bbb(1) + aaa(4) * bbb(2) + aaa(7) * bbb(3)
    ccc(2) = aaa(2) * bbb(1) + aaa(5) * bbb(2) + aaa(8) * bbb(3)
    ccc(3) = aaa(3) * bbb(1) + aaa(6) * bbb(2) + aaa(9) * bbb(3)

    ccc(4) = aaa(1) * bbb(4) + aaa(4) * bbb(5) + aaa(7) * bbb(6)
    ccc(5) = aaa(2) * bbb(4) + aaa(5) * bbb(5) + aaa(8) * bbb(6)
    ccc(6) = aaa(3) * bbb(4) + aaa(6) * bbb(5) + aaa(9) * bbb(6)

    ccc(7) = aaa(1) * bbb(7) + aaa(4) * bbb(8) + aaa(7) * bbb(9)
    ccc(8) = aaa(2) * bbb(7) + aaa(5) * bbb(8) + aaa(8) * bbb(9)
    ccc(9) = aaa(3) * bbb(7) + aaa(6) * bbb(8) + aaa(9) * bbb(9)

  End Subroutine mat_mul

  Function Factorial(n)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function to determine the logarithm of a factorial (n!)
    !
    ! copyright - daresbury laboratory
    ! author    - h.a.boateng april 2014
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent(In   ) :: n
    Real(Kind=wp)          :: Factorial

    Integer :: i

    Factorial = 0.0_wp
    Do i = 2, n
      Factorial = Factorial + Log(Real(i, Kind=wp))
    End Do

  End Function Factorial

  Subroutine factor(n, facs)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 prime factorisability decomposition function
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush august 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,               Intent(In   ) :: n
    Integer, Dimension(:), Intent(  Out) :: facs

    Integer :: i, left, p

    facs = 0

    left = n
    Do i = 1, Size(facs) - 1
      p = get_nth_prime(i)

      If (p <= 0) Exit

      Do While (p * (left / p) == left)
        left = left / p

        facs(i) = facs(i) + 1
      End Do
    End Do

    facs(Size(facs)) = left

  End Subroutine factor

  Function get_nth_prime(n)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 return n-th prime function
    !
    ! copyright - daresbury laboratory
    ! author    - i.j.bush august 2010
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent(In   ) :: n
    Integer                :: get_nth_prime

    Integer, Dimension(1:170), Parameter :: primes = (/2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, &
                                            41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101&
                                            , 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157&
                                            , 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223&
                                            , 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277&
                                            , 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349&
                                            , 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419&
                                            , 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479&
                                            , 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563&
                                            , 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619&
                                            , 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691&
                                            , 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769&
                                            , 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853&
                                            , 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929&
                                            , 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, &
                                            1009, 1013/)

    If (n <= Size(primes)) Then
      get_nth_prime = primes(n)
    Else
      get_nth_prime = -1
    End If

  End Function get_nth_prime

  Pure Function equal_real_wp(a, b) Result(equal)
    Real(Kind=wp), Intent(In   ) :: a, b
    Logical                      :: equal

    equal = Abs(a - b) < epsilon_wp
  End Function equal_real_wp

  Pure Function nequal_real_wp(a, b) Result(nequal)
    Real(Kind=wp), Intent(In   ) :: a, b
    Logical                      :: nequal

    nequal = .not. equal_real_wp(a, b)
  End Function nequal_real_wp

End Module numerics
