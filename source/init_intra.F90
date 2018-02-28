Subroutine init_intra()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for initialising local bookkeeping of intra-like
! interactions: core-shell, bond constraints, PMF constraints, RBs,
! tethered atoms, chemical bonds, valence angles, torsion and improper
! torsion angles, and inversion angles; with exclusions too at the top
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES

  Use kinds, only : wp

! CONFIG MODULE

  Use config_module

! INTERACTION MODULES

  Use core_shell_module

  Use constraints_module
  Use pmf_module

  Use rigid_bodies_module

  Use tethers_module

  Use bonds_module
  Use angles_module
  Use dihedrals_module
  Use inversions_module

  Implicit None

! exclusions locals

  lexatm = 0

! core-shell locals

  ntshl  = 0 ; ntshl1 = 0 ; ntshl2 = 0
  listshl = 0
  legshl  = 0

! constraints locals

  ntcons  = 0 ; ntcons1 = 0
  listcon = 0
  legcon  = 0

! PMFs locals

  ntpmf  = 0
  listpmf = 0
  legpmf  = 0

! RBs locals

  ntrgd  = 0 ; ntrgd1 = 0
  listrgd = 0
  legrgd  = 0

  rgdfrz = 0 ; rgdind = 0 ; indrgd = 0

  rgdwgt = 0.0_wp ; rgdwg1 = 0.0_wp
  rgdx   = 0.0_wp ; rgdy   = 0.0_wp ; rgdz   = 0.0_wp
  rgdrix = 0.0_wp ; rgdriy = 0.0_wp ; rgdriz = 0.0_wp
  rgdaxs = 0.0_wp

  q0 = 0.0_wp ; q1 = 0.0_wp ; q2 = 0.0_wp ; q3 = 0.0_wp

  rgdxxx = 0.0_wp ; rgdyyy = 0.0_wp ; rgdzzz = 0.0_wp
  rgdvxx = 0.0_wp ; rgdvyy = 0.0_wp ; rgdvzz = 0.0_wp
  rgdoxx = 0.0_wp ; rgdoyy = 0.0_wp ; rgdozz = 0.0_wp

! tethers locals

  ntteth = 0
  listtet = 0
  legtet  = 0

! bonds locals

  ntbond  = 0 ; ntbond1 = 0
  listbnd  = 0
  legbnd   = 0

! angles locals

  ntangl  = 0 ; ntangl1 = 0
  listang = 0
  legang  = 0

! dihedrals locals

  ntdihd  = 0 ; ntdihd1 = 0
  listdih = 0
  legdih  = 0

! inversions locals

  ntinv  = 0 ; ntinv1 = 0
  listinv = 0
  leginv  = 0

End Subroutine init_intra
