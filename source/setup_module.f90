Module setup_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 MAIN insert module setting fundamental parameters for the
! entire package @ compile time and specifing execution parameters
! set @ execution time
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
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

  Use kinds_f90

  Implicit None


! FIXED PARAMETERS
! standard pi values

  Real( Kind = wp ), Parameter ::    pi = 3.1415926535897932e0_wp
  Real( Kind = wp ), Parameter :: sqrpi = 1.7724538509055160e0_wp

! standard square roots

  Real( Kind = wp ), Parameter :: rt2 = 1.4142135662373095e0_wp
  Real( Kind = wp ), Parameter :: rt3 = 1.7320508075688772e0_wp

! conversion factor for coulombic terms in internal units, i.e.
! { unit(charge)^2 / [4*pi*eps0*unit(length)] } / unit(energy)

  Real( Kind = wp ), Parameter :: r4pie0 = 138935.4835e0_wp

! boltzmann constant in internal units

  Real( Kind = wp ), Parameter :: boltz = 8.31451115e-1_wp

! conversion factor for pressure from internal units to katm

  Real( Kind = wp ), Parameter :: prsunt = 1.63882576e-1_wp

! conversion factor for surface tension from internal units to dyn/cm

  Real( Kind = wp ), Parameter :: tenunt = 1.660540200_wp

! I/O CHANNELS :: STERR = 0 , STINP = 5 , STOUT = 6 , STERR+STOUT = *
! main input channel

  Integer, Parameter :: nread  = 5

! configuration file input channel

  Integer, Parameter :: nconf  = 11

! force field input channel

  Integer, Parameter :: nfield = 12

! tabulated potential file channel

  Integer, Parameter :: ntable = 13

! reference file input channel

  Integer, Parameter :: nrefdt = 14

! main output channel

  Integer, Parameter :: nrite  = 6

! statistical data file output channel

  Integer, Parameter :: nstats = 21

! accumulators restart dump file

  Integer, Parameter :: nrest  = 22

! trajectory history file channel

  Integer, Parameter :: nhist  = 23

! defect file output channel

  Integer, Parameter :: ndefdt = 24

! rdf file channel number

  Integer, Parameter :: nrdfdt = 25

! z-density file channel number

  Integer, Parameter :: nzdndt = 26

! displacements file channel number

  Integer, Parameter :: nrsddt = 27

! intramolecular PDF file channels numbers

  Integer, Parameter :: npdfdt = 28, &
                        npdgdt = 29


! Random seeding

  Logical, Save :: lseed     = .false.
  Integer, Save :: seed(1:3) = 0

! GLOBAL PARAMETERS FOR ARRAYS' BOUNDS LIMITS (set_bounds)

  Integer, Save ::                                          &
    mxsite,mxatyp,mxtmls,mxexcl,                            &
    mxspl,kmaxa,kmaxb,kmaxc,kmaxa1,kmaxb1,kmaxc1,           &
    mxtshl,mxshl,mxfshl,mxtcon,mxcons,mxfcon,mxlshp,mxproc, &
    mxtpmf(1:2),mxpmf,mxfpmf,mxtrgd,mxrgd,mxlrgd,mxfrgd,    &
    mxtteth,mxteth,mxftet,mxpteth,                          &
    mxtbnd, mxbond,mxfbnd,mxpbnd,                           &
    mxtang, mxangl,mxfang,mxpang,                           &
    mxtdih, mxdihd,mxfdih,mxpdih,                           &
    mxtinv, mxinv, mxfinv,mxpinv,                           &
    mxgrid,mxrdf,mxgrdf,mxvdw,mxpvdw,                       &
    mxgana,mxgbnd,mxgang,mxgdih,mxginv,                     &
    mxmet,mxmed,mxmds,mxpmet,mxter,mxpter,                  &
    mxtbp,mx2tbp,mxptbp,mxfbp,mx3fbp,mxpfbp,mxpfld,         &
    mxstak,mxnstk,mxlist,mxcell,mxatms,mxatdm,              &
    mxbfdp,mxbfss,mxbfxp,mxbfsh,mxbuff

! zero+ and half+/- :: defined in set_bounds

  Real( Kind = wp ), Save :: zero_plus,half_plus,half_minus

! ENGUNIT:: defined in read_field
! engunit = 9648.530821_wp for eV - most used
! engunit = 100.0_wp for kJ/mol   - rarely used
! engunit = 1.0_wp for 10 J/mol   - internal units == default
! engunit = boltz for K/Boltzmann - rarely used

  Real( Kind = wp ), Save :: engunit = 1.0_wp

End Module setup_module
