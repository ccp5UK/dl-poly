Module constants

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 MAIN insert module setting fundamental parameters for the
  ! entire package @ compile time and specifying execution parameters
  ! set @ execution time
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov august 2016
  ! contrib   - i.t.todorov november 2016
  ! contrib   - a.m.elena february 2017
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  ! Note(1): The following internal units apply everywhere
  !
  ! unit of time      (to)    = 1.000000000 x 10**(-12) seconds   (pico-seconds)
  ! unit of length    (lo)    = 1.000000000 x 10**(-10) metres    (Angstroms)
  ! unit of mass      (mo)    = 1.660540200 x 10**(-27) kilograms (Daltons)
  ! unit of charge    (qo)    = 1.602177330 x 10**(-19) coulombs  (electrons)
  ! unit of energy    (eo)    = 1.660540200 x 10**(-23) joules    (10 J mol^-1)
  ! unit of pressure  (po)    = 1.660540200 x 10**(  7) pascals
  !                           = 1.638825760 x 10**(  2) atmospheres
  !
  ! Note(2): All modules, defining (and allocating) module specific
  !          variables, MUST initialise ALL of THEM to ZERO
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp

  Implicit None


  ! Version particulars

  Character( Len =  8), Parameter :: DLP_VERSION = "4.10alph"
  Character( Len = 14), Parameter :: DLP_RELEASE = "  october 2018"

  ! FIXED PARAMETERS
  ! standard pi related values

  Real( Kind = wp ), Parameter ::    pi  =  3.14159265358979312e0_wp ! 2.0_wp*Asin(1.0_wp)
  Real( Kind = wp ), Parameter :: twopi  =  6.28318530717958623e0_wp ! 4.0_wp*Asin(1.0_wp)
  Real( Kind = wp ), Parameter :: fourpi = 12.56637061435917246e0_wp ! 8.0_wp*Asin(1.0_wp)
  Real( Kind = wp ), Parameter :: sqrpi  =  1.77245385090551588e0_wp ! Sqrt(pi)
  Real( Kind = wp ), Parameter :: rtwopi =  0.15915494309189535e0_wp ! 1.0_wp/twopi

  ! standard square roots

  Real( Kind = wp ), Parameter :: rt2    =  1.41421356237309515e0_wp ! Sqrt(2.0_wp)
  Real( Kind = wp ), Parameter :: rt3    =  1.73205080756887719e0_wp ! Sqrt(3.0_wp)

  ! conversion factor for coulombic terms in internal units, i.e.
  ! { unit(charge)^2 / [4*pi*eps0*unit(length)] } / unit(energy)

  Real( Kind = wp ), Parameter :: r4pie0 = 138935.4835e0_wp

  ! boltzmann constant in internal units

  Real( Kind = wp ), Parameter :: boltz = 8.31451115e-1_wp

  ! Energy unit for OUTPUT and STATIS - defined in read_field

  Real( Kind = wp ), Parameter :: eu_ev   = 9648.530821_wp, & ! for eV       - most used
    eu_kcpm =  418.4_wp,      & ! for kcal/mol - often used
    eu_kjpm =  100.0_wp         ! for kJ/mol   - rarely used
  !                                 en_kpb  =     boltz         ! for K/Boltzmann - very rarely used
  Real( Kind = wp ), Save      :: engunit =    1.0_wp         ! for 10 J/mol - internal units == default

 
  Real( Kind = wp ), Parameter :: VA_to_dl = 1.037837512e-4_wp
  Real( Kind = wp ), Parameter :: tesla_to_dl = 1.037837512e4_wp


  ! conversion factor for pressure from internal units to katm

  Real( Kind = wp ), Parameter :: prsunt = 1.63882576e-1_wp

  ! conversion factor for surface tension from internal units to dyn/cm

  Real( Kind = wp ), Parameter :: tenunt = 1.660540200_wp

  ! Maximum bin sizes for distance and angle grids

  Real( Kind = wp ), Parameter :: delr_max  = 0.01_wp ! Angstroms
  Real( Kind = wp ), Parameter :: delth_max = 0.20_wp ! degrees

  ! I/O CHANNELS :: STERR = 0 , STINP = 5 , STOUT = 6 , STERR+STOUT = *

  ! tabulated potential file channel

  Integer, Parameter :: ntable = 13

  ! reference file input channel

  Integer, Parameter :: nrefdt = 14

  ! defect file output channel

  Integer, Parameter :: ndefdt = 24

  ! rdf%rdf file channel number

  Integer, Parameter :: nrdfdt = 25

  ! z-density file channel number

  Integer, Parameter :: nzdndt = 26

  ! displacements file channel number

  Integer, Parameter :: nrsddt = 27

  ! intramolecular PDF file channels numbers

  Integer, Parameter :: npdfdt = 28, &
    npdgdt = 29

  ! vaf file channel number

  Integer, Parameter :: nvafdt = 30

  ! multipoles file channel number

  Integer, Parameter :: nmpldt = 31

  ! ICOORD channel number
  Integer, parameter :: nicrdt = 32

  ! CCORD channel number 
  Integer, parameter :: nccrdt = 33

  !> +0.0 in working precision
  Real( Kind = wp ), Parameter :: zero_plus = Tiny(1.0_wp)
  !> Nearest number to 0.5, greater than 0.5, in working precision
  Real( Kind = wp ), Parameter :: half_plus = Nearest(0.5_wp,+1.0_wp)
  !> Nearest number to 0.5, less than 0.5, in working precision
  Real( Kind = wp ), Parameter :: half_minus = Nearest(0.5_wp,-1.0_wp)
  !> Smallest difference between floats of kind wp
  Real( Kind = wp), Parameter :: epsilon_wp = Epsilon(epsilon_wp)
End Module constants
