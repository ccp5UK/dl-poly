Module constants

  !!------------------------------------------------------------------------!
  !!
  !! dl_poly_4 MAIN insert module setting fundamental parameters for the
  !! entire package @ compile time and specifying execution parameters
  !! set @ execution time
  !!
  !! copyright - daresbury laboratory
  !! author    - i.t.todorov august 2016
  !! contrib   - i.t.todorov november 2016
  !! contrib   - a.m.elena february 2017
  !! refactoring:
  !!           - a.m.elena march-october 2018
  !!           - j.madge march-october 2018
  !!           - a.b.g.chalk march-october 2018
  !!           - i.scivetti march-october 2018
  !!
  !! Note(1): The following internal units apply everywhere
  !!
  !! unit of time      (to)    = 1.000000000 x 10**(-12) seconds   (pico-seconds)
  !! unit of length    (lo)    = 1.000000000 x 10**(-10) metres    (Angstroms)
  !! unit of mass      (mo)    = 1.660540200 x 10**(-27) kilograms (Daltons)
  !! unit of charge    (qo)    = 1.602177330 x 10**(-19) coulombs  (electrons)
  !! unit of energy    (eo)    = 1.660540200 x 10**(-23) joules    (10 J mol^-1)
  !! unit of pressure  (po)    = 1.660540200 x 10**(  7) pascals
  !!                           = 1.638825760 x 10**(  2) atmospheres
  !!
  !! Note(2): All modules, defining (and allocating) module specific
  !!          variables, MUST initialise ALL of THEM to ZERO
  !!
  !!------------------------------------------------------------------------!

  Use kinds, Only: wp

  Implicit None

  !> Version particulars

  Character(Len=8), Parameter :: DLP_VERSION = "5.1.0"
  Character(Len=14), Parameter :: DLP_RELEASE = "March 2021"

  !> FIXED PARAMETERS
  !> standard pi related values

  Real(Kind=wp), Parameter ::    pi = 4.0_wp * Atan(1.0_wp)
  Real(Kind=wp), Parameter :: twopi = 2.0_wp * pi
  Real(Kind=wp), Parameter :: fourpi = 4.0_wp * pi
  Real(Kind=wp), Parameter :: sqrpi = Sqrt(pi)
  Real(Kind=wp), Parameter :: rtwopi = 1.0_wp / twopi
  Real(Kind=wp), Parameter :: rsqrpi = 1.0_wp / (Sqrt(pi))

  !> Gamma function in steps of 1/2
  Real(Kind=wp), Dimension(24), Parameter :: gamma_1_2 = [ &
                                             1.7724538509055161040_wp, & ! 1
                                             1.0000000000000000000_wp, & ! 2
                                             0.8862269254527580519_wp, & ! 3
                                             1.0000000000000000000_wp, & ! 4
                                             1.3293403881791370225_wp, & ! 5
                                             2.0000000000000000000_wp, & ! 6
                                             3.3233509704478425562_wp, & ! 7
                                             6.0000000000000000000_wp, & ! 8
                                             11.631728396567449835_wp, & ! 9
                                             24.000000000000000000_wp, & ! 10
                                             52.342777784553526033_wp, & ! 11
                                             120.00000000000000000_wp, & ! 12
                                             287.88527781504438963_wp, & ! 13
                                             720.00000000000000000_wp, & ! 14
                                             1871.2543057977884473_wp, & ! 15
                                             5040.0000000000000000_wp, & ! 16
                                             14034.407293483413014_wp, & ! 17
                                             40320.000000000000000_wp, & ! 18
                                             119292.46199460900971_wp, & ! 19
                                             362880.00000000000000_wp, & ! 20
                                             1133278.3889487856068_wp, & ! 21
                                             3628800.0000000000000_wp, & ! 22
                                             11899423.083962248638_wp, & ! 23
                                             39916800.000000000000_wp & ! 24
                                             ]

  !> Reciprocal of the gamma function
  Real(Kind=wp), Dimension(24), Parameter :: inv_gamma_1_2 = 1.0_wp / gamma_1_2

  !> standard square roots

  Real(Kind=wp), Parameter :: rt2 = 1.41421356237309515e0_wp !! Sqrt(2.0_wp)
  Real(Kind=wp), Parameter :: rt3 = 1.73205080756887719e0_wp !! Sqrt(3.0_wp)

  !> conversion factor for coulombic terms in internal units, i.e.
  !> { unit(charge)^2 / [4*pi*eps0*unit(length)] } / unit(energy)

  Real(Kind=wp), Parameter :: r4pie0 = 138935.4835e0_wp

  !> boltzmann constant in internal units

  Real(Kind=wp), Parameter :: boltz = 8.31451115e-1_wp

  !> Energy unit for OUTPUT and STATIS - defined in read_field

  Real(Kind=wp), Parameter :: eu_ev = 9648.530821_wp, & !! for eV       - most used
                              eu_kcpm = 418.4_wp, & !! for kcal/mol - often used
                              eu_kjpm = 100.0_wp !! for kJ/mol   - rarely used
  !                                 en_kpb  =     boltz         ! for K/Boltzmann - very rarely used
  Real(Kind=wp), Save      :: engunit = 1.0_wp !! for 10 J/mol - internal units == default

  Real(Kind=wp), Parameter :: VA_to_dl = 1.037837512e-4_wp
  Real(Kind=wp), Parameter :: tesla_to_dl = 1.037837512e4_wp

  !> conversion factor for pressure from internal units to katm

  Real(Kind=wp), Parameter :: prsunt = 1.63882576e-1_wp

  !> conversion factor for surface tension from internal units to dyn/cm

  Real(Kind=wp), Parameter :: tenunt = 1.660540200_wp

  !> TTM conversion factors
  Real(Kind=wp), Parameter :: JKms_to_kBAps = 10.0_wp / (boltz * tenunt) ! convert W m^-1 K^-1 to kB A^-1 ps^-1
  Real(Kind=wp), Parameter :: Jm3K_to_kBA3 = 1.0e-7_wp / (boltz * tenunt) ! convert J m^-3 K^-1 to kB A^-3
  Real(Kind=wp), Parameter :: kB_to_eV = boltz / eu_ev ! convert kB to eV
  Real(Kind=wp), Parameter :: eV_to_kB = eu_ev / boltz ! convert eV to kB
  Real(Kind=wp), Parameter :: mJcm2_to_eVA2 = 1.0e4_wp / (eu_ev * tenunt) ! convert mJ cm^-2 to eV A^-2


  !> Maximum bin sizes for distance and angle grids

  Real(Kind=wp), Parameter :: delr_max = 0.01_wp !! Angstroms
  Real(Kind=wp), Parameter :: delth_max = 0.20_wp !! degrees

  ! I/O CHANNELS :: STERR = 0 , STINP = 5 , STOUT = 6 , STERR+STOUT = *

  !> reference file input channel

  Integer, Parameter :: nrefdt = 14

  !> defect file output channel

  Integer, Parameter :: ndefdt = 24

  !> z-density file channel number

  Integer, Parameter :: nzdndt = 26

  !> displacements file channel number

  Integer, Parameter :: nrsddt = 27

  !> intramolecular PDF file channels numbers

  Integer, Parameter :: npdfdt = 28, &
                        npdgdt = 29

  !> vaf file channel number
  Integer, Parameter :: nvafdt = 30

  !> multipoles file channel number
  Integer, Parameter :: nmpldt = 31

  !> ICOORD channel number
  Integer, Parameter :: nicrdt = 32

  !> CCORD channel number
  Integer, Parameter :: nccrdt = 33

  !> ADF channel number
  Integer, Parameter :: nchadf = 34

  !> maximum number of  bsplines to be used for spme
  Integer, Parameter :: MAX_BSPLINE = 20

  Integer, Private :: i
  Real(Kind=wp), Dimension(1:MAX_BSPLINE), Parameter :: &
          real_no = [(Real(i,wp), i = 1,MAX_BSPLINE)], &
          inv_no = [(1.0_wp / Real(i, wp), i = 1, MAX_BSPLINE)]  !! Real variants to avoid type-casting

  !> +0.0 in working precision
  Real(Kind=wp), Parameter     :: zero_plus = Tiny(1.0_wp)
#ifdef NVIDIA
  Real(Kind=wp) :: half_plus
  Real(Kind=wp) :: half_minus
#else
  !> Nearest number to 0.5, greater than 0.5, in working precision
  Real(Kind=wp), Parameter     :: half_plus = Nearest(0.5_wp, +1.0_wp)
  !> Nearest number to 0.5, less than 0.5, in working precision
  Real(Kind=wp), Parameter     :: half_minus = Nearest(0.5_wp, -1.0_wp)
#endif
  !> Smallest difference between floats of kind wp
  Real(Kind=wp), Parameter     :: epsilon_wp = Epsilon(epsilon_wp)
  !> complex zero
  Real(Kind=wp), Parameter     :: czero = Cmplx(0.0_wp, 0.0_wp, wp)
  !> smallest distance we care about for link cells
  Real(Kind=wp), Parameter     :: smalldr = 1.0e-6_wp

  !> New line char
  Character, Parameter :: lf = new_line('a')

End Module constants
