Module errors_warnings
!> Module controlling errors and warnings
!>
!> Copyright - Daresbury Laboratory
!
!> Author - a.m.elena March 2018
!> refactoring:
!>           - a.m.elena march-october 2018
!>           - j.madge march-october 2018
!>           - a.b.g.chalk march-october 2018
!>           - i.scivetti march-october 2018
!> contrib - a.m.elena October 2018 - use standard integer for units
!> contrib - a.m.elena March 2019 - remove error 145
!> contrib - a.m.elena March 2019 - fix wrong logic in warning

  Use comms,                         Only: abort_comms,&
                                           comms_type
  Use, Intrinsic :: iso_fortran_env, Only: error_unit,&
                                           input_unit,&
                                           output_unit,&
                                           iostat_end,&
                                           iostat_eor
  Use kinds,                         Only: wp

  Implicit None

  Private

  Type(comms_type), Save :: eworld
  Integer, Save :: print_level = 2
  Integer, Save :: ounit

  Public :: warning
  Public :: error
  Public :: info
  Public :: init_error_system
  Public :: error_alloc, error_dealloc
  Public :: error_read
  Public :: set_print_level
  Public :: get_print_level
  Public :: check_print_level

  Interface warning
    Module Procedure warning_special
    Module Procedure warning_general
  End Interface warning

  Interface info
    Module Procedure info_sl
    Module Procedure info_ml
  End Interface info

Contains

  Subroutine init_error_system(nrite, comm)

    Integer,          Intent(In   ) :: nrite
    Type(comms_type), Intent(In   ) :: comm

    eworld%comm = comm%comm
    eworld%idnode = comm%idnode
    ounit = nrite

  End Subroutine init_error_system

  Subroutine warning_special(kode, a, b, c, triggered, level_once, level_always)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for printing warning messages and returning
    ! control back to the main program
    !
    ! triggered - return logical as to whether this has triggered
    ! level_once - print level for this to trigger only once (only works
    !              in conjunction with triggered)
    ! level_always - print level for this to always trigger
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2017
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,       Intent(In   ) :: kode
    Real(Kind=wp), Intent(In   ) :: a, b, c
    Logical,       Intent(InOut), Optional :: triggered
    Integer,       Intent(In   ), Optional :: level_once, level_always
    Integer :: l1, la
    Integer :: ia, ib, ic

    if (present(level_once)) then
      l1 = level_once
    else
      l1 = 0
    end if

    if (present(level_always)) then
      la = level_always
    else
      la = 3
    end if

    if (present(triggered) .and. .not. check_print_level(la)) then
      if (triggered .or. .not. check_print_level(l1)) then
        return
      else
        triggered = .true.
      end if
    end if

    If (eworld%idnode == 0) Then

      Write (ounit, '(/,1x,a,i6)') 'warning issued ', kode

      select case (kode)
      Case (2)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - DD with ', ia, ' idle nodes, (out of ', ib, ') mapped on vacuum !!! ***'

      Case (3)

        Write (ounit, '(/,1x,3(a,f7.3),a,/)') &
          '*** warning - DD cutoff(+padding) is: ', a, ' Angstroms while minimum half-cell width is: ', b, ' Angstroms !!! ***'

      Case (4)

        Write (ounit, '(/,2(1x,a,/))') &
          '*** warning - system with uncharged particles !!! ***', &
          '*** "no elec" or/and "no strict" directives in CONTROL may speed up simulation !!! ***'

      Case (5)

        Write (ounit, '(/,1x,a,f12.4,a,/)') &
          '*** warning - non-zero total system charge: ', a, ' positrons !!! ***'

      Case (6)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/,1x,a,a,/,1x,a)') &
          '*** warning - maximum length of linked cell list: ', ia, &
          ' + 1 is less than maximum length of particle exclusion list: ', &
          ib, ' !!! ***', &
          '*** this may be due to using too short a cutoff in CONTROL ', &
          'and/or a badly defined intramolecular topology in FIELD !!! ***', &
          '*** further potential problems may be expected !!! ***'

      Case (7)

        Write (ounit, '(/,1x,a,f7.3,a,/,1x,a,/)') &
          '*** warning - DD cutoff is ', a, ' Angstroms !!! ***', &
          '*** Fennell damping is not recommended for cutoffs shorter than 10-12 Angstroms !!! ***'

      Case (8)

        Write (ounit, '(/,1x,a,2(f7.3,a),/)') &
          '*** warning - : detected maximum rigid body width: ', a, ' Angstroms while the DD cutoff is ', b, ' Angstroms !!! ***'

      Case (10)

        Write (ounit, '(/,1x,a,/)') '*** warning - no pair forces in use !!! ***'

      Case (20)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - : 1..4 scale factors reset for molecule: ', ia, ' sites: ', ib, ' & ', ic, ' !!! ***'

      Case (22)

        Write (ounit, '(/,1x,a,/)') '*** warning - : 1..4 scale factors reset for dihedrals in the system !!! ***'

      Case (30)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - Electrostatics requested in a non-periodic system !!! ***'

      Case (34)

        Write (ounit, '(/,1x,a,/)') '*** warning - redundant directive prim(ary cutoff) ignored !!! ***'

      Case (35)

        Write (ounit, '(3(/,1x,a),/)') &
          "*** warning - DL_POLY_2/Classic directive 'delr - Verlet shell strip cutoff' defaulted to ***", &
          "*** DL_POLY_4 directive 'neigh%padding - real space cutoff padding option' ***", &
          "*** neigh%padding=Max(neigh%padding,delr/4) !!! ***"

      Case (36)

        Write (ounit, '(2(/,1x,a),/)') &
          "*** warning - DL_POLY_2/Classic directive 'mult(iple timestep)' defaulted to ***", &
          "*** DL_POLY_4 directive 'infrequent k-space SPME evaluation' !!! ***"

      Case (37)

        Write (ounit, '(/,1x,a,/)') '*** warning - redundant directive all pairs ignored !!! ***'

      Case (38)

        Write (ounit, '(/,1x,a,/)') '*** warning - redundant directive no link ignored !!! ***'

      Case (40)

        Write (ounit, '(/,1x,a,f7.3,a,/)') '*** warning - tentative cutoff reset to ', a, ' Angstroms !!! ***'

      Case (50)

        Write (ounit, '(2(/,1x,a),f7.3,a,/)') &
          '*** warning - short-range interactions cutoff reset ***', '*** new cutoff radius (rvdw) ', a, ' Angstroms !!! ***'

      Case (60)

        Write (ounit, '(/,1x,a,f7.3,a,/)') '*** warning - total system charge ', a, ' positrons !!! ***'

      Case (70)

        Write (ounit, '(/,1x,a,f7.3,a,/)') '*** warning - switching length reset to ', a, ' Angstroms !!! ***'

      Case (80)

        Write (ounit, '(/,1x,a,/)') '*** warning - requested thermostat unavailable !!! ***'

      Case (90)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        If (ic == 0) Then
          Write (ounit, '(/,1x,a,2(i0,a),/,1x,a,/)') &
            '*** warning - allocating more link-cells ', ia, ' than initially envisaged ', ib, ' , in linkcell_pairs !!! ***', &
            '*** System volume has expanded beyond what was safely presumed as physically sensible !!! ***'
        Else If (ic == 1) Then
          Write (ounit, '(/,1x,a,2(i0,a),/,1x,a,/)') &
            '*** warning - allocating more link-cells ', ia, ' than initially envisaged ', ib, ' , in three_body_forces !!! ***', &
            '*** System volume has expanded beyond what was safely presumed as physically sensible !!! ***'
        Else If (ic == 2) Then
          Write (ounit, '(/,1x,a,2(i0,a),/,1x,a,/)') &
            '*** warning - allocating more link-cells ', ia, ' than initially envisaged ', ib, ' , in four_body_forces !!! ***', &
            '*** System volume has expanded beyond what was safely presumed as physically sensible !!! ***'
        Else If (ic == 3) Then
          Write (ounit, '(/,1x,a,2(i0,a),/,1x,a,/)') &
            '*** warning - allocating more link-cells ', ia, ' than initially envisaged ', ib, ' , in tersoff_forces !!! ***', &
            '*** System volume has expanded beyond what was safely presumed as physically sensible !!! ***'
        Else
          Write (ounit, '(/,1x,a,/)') &
            '*** unspecified warning encountered !!! ***'
        End If

      Case (100)

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - primary link cell algorithm has a link cell dimension that is < 3 !!! ***', &
          '*** DL_POLY_4 RUNNING IN LOW EFFICIENCY MODE !!! ***'

      Case (110)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(2(/,1x,a),2(i0,a),/)') &
          '*** warning - image convention incompatible with the set NsT ensemble ***', &
          '*** config%imcon reset from ', ia, ' to ', ib, ' !!! ***'

      Case (120)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - unspecified atom-atom interactions set to zero !!! ***'

      Case (130)

        Write (ounit, '(/,1x,a,/)') '*** warning - no ensemble is specified !!! ***'

      Case (135)

        Write (ounit, '(/,1x,a,/,1x,a,/)') &
          '*** warning - incorrect ensemble specified for two-temperature model !!! ***', &
          '*** replacing ensemble with NVT inhomogeneous Langevin !!! ***'

      Case (140)

        Write (ounit, '(/,1x,2a,2(f8.5,a),/,1x,a,/)') &
          '*** warning - control distances for variable timestep: ', &
          'Dmin = ', a, ' and Dmax = ', b, ' (Angstroms) !!! ***', &
          '*** do not comply with safty condition: Dmax > 2.5 Dmin > 0 !!! ***'

      Case (150)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required export buffer size ', ia, ' and actual: ', ib, ' !!! ***'

      Case (160)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required import array size ', ia, ' and actual: ', ib, ' !!! ***'

      Case (190)

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - REVOLD format mishmash detected (restart requested) !!! ***', &
          '*** restart is abandoned and clean start is assumed !!! ***'

      Case (200)

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - CONFIG contains positions only !!! ***', &
          '*** clean start is assumed !!! ***'

      Case (210)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(3(/,1x,a),2(i0,a),/)') &
          '*** warning - system under great constraint !!! ***', &
          '*** total degrees of freedom <=  total number of particles !!! ***', &
          '*** ', ia, ' <= ', ib, ' !!! ***'

      Case (220)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - Ewald sum requested in a non-periodic system !!! ***'

      Case (230)

        ia = Nint(a)

        Write (ounit, '(/,1x,a,i0,a,/,1x,a,/)') &
          '*** warning - PMF unit ', ia, ' config%weight is detected zero !!! ***', &
          '*** member config%weights defaulted to atom type masses (or units) !!! ***'

      Case (240)

        ia = Nint(a)

        Write (ounit, '(/,1x,a,i0,a,/)') &
          '*** warning - total number of nodes running in parallel: ', ia, ' !!! ***'

      Case (250)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required coulombic exclusion array size ', ia, ' and actual (neigh%max_exclude): ', ib, ' !!! ***'

      Case (260)

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - system volume in non-periodic systems is the MD config%cell volume !!! ***', &
          '*** system pressure is calculated with respect to this volume !!! ***'

      Case (270)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required reading buffer size is ', ia, ' and actual ', ib, ' !!! ***'

      Case (280)

        Write (ounit, '(/,1x,a,/,1x,2(a,f9.3),a,/)') &
          '*** warning - pseudo thermostat cannot be applied for this model system since !!! ***', &
          '*** specified thermostat wall thickness ', a, ' > 1/4 minimum MD cell width ', b, ' (Angstroms) !!! ***'

      Case (290)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required link-config%cell neigh%list size is ', ia, ' and actual (neigh%max_list) ', ib, ' !!! ***'

      Case (295)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - PMF unit ', ia, ' and rigid body unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (296)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, &
          ' on molecular species type ', ib, &
          ' with compromised polarisability !!! ***'

      Case (297)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and angle unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (298)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and dihedral unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (299)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and inversion unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (300)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and PMF unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (301)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and constraint unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (302)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - core-shell unit ', ia, ' and rigid body unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (303)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and tether unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (304)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - constraint unit ', ia, ' and rigid body unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Case (305)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a))') &
          '*** warning - rigid body unit ', ia, ' on molecular species type ', ib, &
          ' forced to freeze !!! ***'

      Case (306)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a))') &
          '*** warning - constraint unit ', ia, ' on molecular species type ', ib, &
          ' forced to freeze !!! ***'

      Case (307)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - rigid body unit ', ia, ' on molecular species type ', ib, &
          ' set to have type ', ic, ' is problematic !!! ***'

      Case (308)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - site ', ia, ' of constraint unit ', ib, &
          ' on molecular species type ', ic, ' is problematic !!! ***'

      Case (309)

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - site ', ia, ' , member ', ib, ' of PMF unit ', ic, ' is problematic !!! ***'

      Case (310)

        Write (ounit, '(/,1x,a,/,1x,a,2(f6.2,a),/)') &
          '*** warning - control distance for defect look-up MUST be in the interval [Min(0.3,rcut/3);Min(1.2,rcut/2)] !!! ***', &
          '*** defects distance condition will default from ', a, ' to ', b, ' (Angstroms) !!! ***'

      Case (320)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - REFERENCE not found, CONFIG to be used as a reference for defects detection !!! ***'

      Case (330)

        Write (ounit, '(/,1x,a,/,1x,a,2(f10.6,a),/)') &
          '*** warning - iteration cycles length limit for conjugate gradient minimisers exceded !!! ***', &
          '*** specified convergence tolerance: ', a, ' , needed one for a pass: ', b, ' !!! ***'

      Case (340)

        Write (ounit, '(3(/,1x,a),2(f7.4,a),/,1x,a,f7.4,a,/)') &
          '*** warning - inconsistent binsize for spatial distribution functions !!! ***', &
          '*** 1.0E-5 (Angstroms) <= binsize <= neigh%cutoff/4 (Angstroms) !!! ***', &
          '*** 1.0E-5 (Angstroms) <= ', a, ' <= ', b, ' (Angstroms) !!! ***', &
          '*** binsize defaults to ', c, ' (Angstroms) !!! ***'

      Case (350)

        Write (ounit, '(/,1x,a)') &
          '*** warning - expansion along z axis not allowed for slab geometry, nz defaults to 1 !!! ***'

      Case (360)

        Write (ounit, '(/,1x,a,2(f12.7,a),3(/,1x,a),/)') &
          '*** warning - minimisation tolerance ', a, ' defaults to ', b, ' !!! ***', &
          '*** force   : 1.0    <= tolerance <= 1000.00, default = 50.00 !!! ***', &
          '*** energy  : 0.0    <  tolerance <=    0.01, default = 0.005 !!! ***', &
          '*** distance: 1.0e-6 <= tolerance <=    0.10, default = 0.005 !!! ***'

      Case (370)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),2(/,1x,a),/)') &
          '*** warning - k-space evaluation interval ', ia, ' defaults to ', ib, ' !!! ***', &
          '*** the interval must be a positive integer beteen 1 and 10 !!! ***', &
          '*** values > 10 default to 4, no value or 0 defaults to 1 !!! ***'

      Case (380)

        ia = Nint(a)

        Write (ounit, '(/,1x,a,i0,a,/)') &
          '*** warning - IMPACT applied before the end of the equlibration period (', ia, ' step) !!! ***'

      Case (390)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - core-shell unit mix-up or duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (400)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - constraint or rigid body duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (410)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,i0,a,/)') &
          '*** warning - tethered atom duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (420)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - chemical bond duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (430)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - bond angle duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (440)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - dihedral angle duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (450)

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - inversion angle duplicate specification (', ia, ',', ib, ') !!! ***'

      Case (460)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - further semi-isotropic barostat option search abandoned !!! ***'

      Case (470)

        Write (ounit, '(/,1x,a,/,1x,a,2(f5.2,a),/)') &
          '*** warning - control distance for diplacement qualification MUST be >= 0.25 Angstroms !!! ***', &
          '*** displacements distance condition will default from ', a, ' to ', b, ' Angstroms !!! ***'

      Case (480)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - "metal direct" option disabled as incompatible with EAM potentials (TABEAM) !!! ***'

      Case (490)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - "metal sqrtrho" option disabled as incompatible with FS type potentials (analytic forms) !!! ***'

      Case (500)

        Write (ounit, '(2(/,1x,a,/))') &
          '*** warning - insufficient electronic temperature cells to allow energy redistribution from inactive cells !!! ***', &
          '*** option switched off - electronic temperature conservation no longer guaranteed !!! ***'

      Case (510)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - no laser energy deposition for two-temperature model specified !!! ***'

      Case (515)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - no initial stopping power for two-temperature model specified !!! ***'

      Case (520)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - mismatch in timestep between REVIVE and provided electronic temperature lattice !!! ***'

      Case (530)

        Write (ounit, '(/,1x,a,f6.2,a,/,1x,a,/)') &
          '*** warning - possible time energy deposition discrepancy of at least ', a, '% !!! ***', &
          '*** discrepancy may be due to inactive cells !!! ***'

      Case (535)

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - spatial energy deposition cutoff in xy-directions greater than available ionic temperature cells !!! ***'

      Case (540)

        Write (ounit, '(/,1x,a,f6.2,a,/)') &
          '*** warning - possible spatial energy deposition discrepancy of at least ', a, '% !!! ***'

      Case Default

        Write (ounit, '(/,1x,a,/)') &
          '*** unspecified warning encountered !!! ***'

      End select

    End If

  End Subroutine warning_special

  Subroutine warning_general(message, master_only, triggered, level_once, level_always)
    Character(Len=*),  Intent(In   ) :: message
    Logical, Optional, Intent(In   ) :: master_only
    Logical,       Intent(InOut), Optional :: triggered
    Integer,       Intent(In   ), Optional :: level_once, level_always
    Integer :: l1, la
    Logical :: zeroOnly


    if (present(level_once)) then
      l1 = level_once
    else
      l1 = 0
    end if

    if (present(level_always)) then
      la = level_always
    else
      la = 3
    end if

    if (present(triggered) .and. .not. check_print_level(la)) then
      if (triggered .or. .not. check_print_level(l1)) then
        return
      else
        triggered = .true.
      end if
    end if

    zeroOnly = .false.
    If (Present(master_only)) zeroOnly = master_only

    If (zeroOnly) Then
      If (eworld%idnode == 0) Then
        Write (ounit, '(a)') "*** warning - "//Trim(message)//" !!! ***"
      End If
    Else
      Write (ounit, '(a,1x,i0,a)') "*** warning - "//Trim(message)//", node: ", eworld%idnode, " !!! ***"
    End If

  End Subroutine warning_general

  Subroutine info_ml(message, n, master_only, level)
    Character(Len=*),  Intent(In   ) :: message(:)
    Integer,           Intent(In   ) :: n
    Logical, Optional, Intent(In   ) :: master_only
    Integer, Optional :: level
    Integer :: print_check

    Integer :: i
    Logical :: zeroOnly


    print_check = 1
    if (present(level)) then
       print_check = level
    end if

    if (print_check > print_level) return

    zeroOnly = .false.
    If (Present(master_only)) zeroOnly = master_only

    If (zeroOnly) Then
      If (eworld%idnode == 0) Then
        Do i = 1, n
          Write (ounit, '(a)') Trim(message(i))
        End Do
      End If
    Else
      Do i = 1, n
        Write (ounit, '(a,1x,i0)') Trim(message(i))//", node: ", eworld%idnode
      End Do
    End If
  End Subroutine info_ml

  Subroutine info_sl(message, master_only, level)
    Character(Len=*),  Intent(In   ) :: message
    Logical, Optional, Intent(In   ) :: master_only
    Integer, Optional :: level
    Integer :: print_check
    Logical :: zeroOnly

    print_check = 1
    if (present(level)) then
       print_check = level
    end if

    if (print_check > print_level) return

    zeroOnly = .false.
    If (Present(master_only)) zeroOnly = master_only

    If (zeroOnly) Then
      If (eworld%idnode == 0) Then
        Write (ounit, '(a)') Trim(message)
      End If
    Else
      Write (ounit, '(a,1x,i0)') Trim(message)//", node: ", eworld%idnode
    End If

  End Subroutine info_sl

  Subroutine set_print_level(level)
    Integer, Intent(In   ) :: level

    print_level = level

  End Subroutine set_print_level

  Integer Function get_print_level()
    get_print_level = print_level
  end Function get_print_level

  Function check_print_level(level)
    Integer, Intent(In   ) :: level
    Logical :: check_print_level

    check_print_level = print_level >= level

  end Function check_print_level


  Subroutine error(kode, message, master_only)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for printing error messages and bringing about a
    ! controlled termination of the program
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov september 2016
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                    Intent(In   ) :: kode
    Character(Len=*), Optional, Intent(In   ) :: message
    Logical, Optional,          Intent(In   ) :: master_only

    Logical :: zeroOnly

    zeroOnly = .false.
    If (Present(master_only)) zeroOnly = master_only

    If (zeroOnly) Then
      If (eworld%idnode == 0) Then
        If (Present(message)) Then
          Write (ounit, '(a)') Trim(message)
        End If
      End If
    Else
      If (Present(message)) Then
        Write (ounit, '(a,1x,i0)') Trim(message)//", node: ", eworld%idnode
      End If
    End If

    If (eworld%idnode == 0) Then

      Write (ounit, '(/,1x,a,i5)') 'DL_POLY_4 terminated due to error ', kode

      select case (kode)
      Case (1)

        Write (ounit, '(/,1x,a)') 'error - word_2_real failure'

      Case (2)

        Write (ounit, '(/,1x,a)') 'error - too many atom types in FIELD (scan_field)'

      Case (3)

        Write (ounit, '(/,1x,a)') 'error - unknown directive found in CONTROL file'

      Case (4)

        Write (ounit, '(/,1x,a)') 'error - unknown directive found in FIELD or MPOLES file'

      Case (5)

        Write (ounit, '(/,1x,a)') 'error - unknown energy unit requested'

      Case (6)

        Write (ounit, '(/,1x,a)') 'error - energy unit not specified'

      Case (7)

        Write (ounit, '(/,1x,a)') 'error - selected external field incompatible with selected ensemble (NVE only!!!)'

      Case (8)

        Write (ounit, '(/,1x,a)') 'error - ewald precision MUST be a POSITIVE real number'

      Case (9)

        Write (ounit, '(/,1x,a)') 'error - ewald sum parameters MUST be well defined'

      Case (10)

        Write (ounit, '(/,1x,a)') 'error - too many molecular types specified'

      Case (11)

        Write (ounit, '(/,1x,a)') 'error - duplicate molecule directive in FIELD file'

      Case (12)

        Write (ounit, '(/,1x,a)') 'error - unknown molecule directive in FIELD or MPOLES file'

      Case (13)

        Write (ounit, '(/,1x,a)') 'error - molecular species not yet specified'

      Case (14)

        Write (ounit, '(/,1x,a)') 'error - too many unique atom types specified'

      Case (15)

        Write (ounit, '(/,1x,a)') 'error - duplicate vdw potential specified'

      Case (16)

        Write (ounit, '(/,1x,a)') 'error - strange exit from FIELD file processing'

      Case (17)

        Write (ounit, '(/,1x,a)') 'error - strange exit from CONTROL file processing'

      Case (18)

        Write (ounit, '(/,1x,a)') 'error - duplicate three-body potential specified'

      Case (19)

        Write (ounit, '(/,1x,a)') 'error - duplicate four-body potential specified'

      Case (20)

        Write (ounit, '(/,1x,a)') 'error - too many molecule sites specified'

      Case (21)

        Write (ounit, '(/,1x,a)') 'error - molecule contains more atoms/sites than declared'

      Case (22)

        Write (ounit, '(/,1x,a)') 'error - unsuitable radial increment in TABLE||TABBND||TABANG||TABDIH||TABINV file'

      Case (23)

        Write (ounit, '(/,1x,a)') 'error - incompatible FIELD and TABLE file potentials'

      Case (24)

        Write (ounit, '(/,1x,a)') 'error - end of file encountered in TABLE||TABEAM||TABBND||TABANG||TABDIH||TABINV file'

      Case (25)

        Write (ounit, '(/,1x,a)') 'error - wrong atom type found in CONFIG file'

      Case (26)

        Write (ounit, '(/,1x,a)') 'error - neutral group option now redundant'

      Case (27)

        Write (ounit, '(/,1x,a)') "error - unit's member indexed outside molecule's site range"

      Case (28)

        Write (ounit, '(/,1x,a)') 'error - wrongly indexed atom entries found in CONFIG file'

      Case (30)

        Write (ounit, '(/,1x,a)') 'error - too many chemical bonds specified'

      Case (31)

        Write (ounit, '(/,1x,a)') 'error - too many chemical bonds per domain'

      Case (32)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in core-shell unit'

      Case (33)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in constraint bond unit'

      Case (34)

        Write (ounit, '(/,1x,a)') 'error - length of constraint bond unit >= real space cutoff (neigh%cutoff)'

      Case (35)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in chemical bond unit'

      Case (36)

        Write (ounit, '(/,1x,a)') 'error - only one *bonds* directive per molecule is allowed'

      Case (38)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in metal_ld_export'

      Case (39)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in metal_ld_export'

      Case (40)

        Write (ounit, '(/,1x,a)') 'error - too many bond constraints specified'

      Case (41)

        Write (ounit, '(/,1x,a)') 'error - too many bond constraints per domain'

      Case (42)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to deport_atomic_data'

      Case (43)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer exceeded in deport_atomic_data'

      Case (44)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in deport_atomic_data'

      Case (45)

        Write (ounit, '(/,1x,a)') 'error - too many atoms in CONFIG file or per domain'

      Case (46)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to export_atomic_data/positions'

      Case (47)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to metal_ld_export'

      Case (48)

        Write (ounit, '(/,1x,a)') 'error - transfer buffer too small in *_table_read'

      Case (49)

        Write (ounit, '(/,1x,a)') 'error - frozen shell (core-shell unit) specified'

      Case (50)

        Write (ounit, '(/,1x,a)') 'error - too many bond angles specified'

      Case (51)

        Write (ounit, '(/,1x,a)') 'error - too many bond angles per domain'

      Case (52)

        Write (ounit, '(/,1x,a)') 'error - end of FIELD or MPOLES file encountered'

      Case (53)

        Write (ounit, '(/,1x,a)') 'error - end of CONTROL file encountered'

      Case (54)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in export_atomic_data/positions'

      Case (55)

        Write (ounit, '(/,1x,a)') 'error - end of CONFIG file encountered'

      Case (56)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in export_atomic_data/positions'

      Case (57)

        Write (ounit, '(/,1x,a)') 'error - too many core-shell units specified'

      Case (58)

        Write (ounit, '(/,1x,a)') 'error - number of atoms in system not conserved'

      Case (59)

        Write (ounit, '(/,1x,a)') 'error - too many core-shell units per domain'

      Case (60)

        Write (ounit, '(/,1x,a)') 'error - too many dihedral angles specified'

      Case (61)

        Write (ounit, '(/,1x,a)') 'error - too many dihedral angles per domain'

      Case (62)

        Write (ounit, '(/,1x,a)') 'error - too many tethered atoms specified'

      Case (63)

        Write (ounit, '(/,1x,a)') 'error - too many tethered atoms per domain'

      Case (64)

        Write (ounit, '(/,1x,a)') 'error - incomplete core-shell unit found in build_book_intra'

      Case (65)

        Write (ounit, '(/,1x,a)') 'error - too many excluded pairs specified'

      Case (66)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in bond angle unit'

      Case (67)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in dihedral unit'

      Case (68)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in inversion unit'

      Case (69)

        Write (ounit, '(/,1x,a)') 'error - too many link cells required'

      Case (70)

        Write (ounit, '(/,1x,a)') 'error - constraint_quench failure'

      Case (71)

        Write (ounit, '(/,1x,a)') 'error - too many metal potentials specified'

      Case (72)

        Write (ounit, '(/,1x,a)') 'error - too many tersoff potentials specified'

      Case (73)

        Write (ounit, '(/,1x,a)') 'error - too many inversion potentials specified'

      Case (74)

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in tersoff potential neigh%list'

      Case (76)

        Write (ounit, '(/,1x,a)') 'error - duplicate tersoff potential specified'

      Case (77)

        Write (ounit, '(/,1x,a)') 'error - too many inversion angles per domain'

      Case (79)

        Write (ounit, '(/,1x,a)') 'error - tersoff potential cutoff undefined'

      Case (80)

        Write (ounit, '(/,1x,a)') 'error - too many pair potentials specified'

      Case (81)

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in pair potential neigh%list'

      Case (82)

        Write (ounit, '(/,1x,a)') 'error - calculated pair potential index too large'

      Case (83)

        Write (ounit, '(/,1x,a)') 'error - too many three-body/angles potentials specified'

      Case (84)

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in three-body/angles potential neigh%list'

      Case (85)

        Write (ounit, '(/,1x,a)') 'error - required velocities not in CONFIG file'

      Case (86)

        Write (ounit, '(/,1x,a)') 'error - calculated three-body potential index too large'

      Case (88)

        Write (ounit, '(/,1x,a)') 'error - legend array exceeded in build_book_intra'

      Case (89)

        Write (ounit, '(/,1x,a)') 'error - too many four-body/dihedrals/inversions potentials specified'

      Case (90)

        Write (ounit, '(/,1x,a)') 'error - specified tersoff potentials have different types'

      Case (91)

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in four-body/dihedrals/iversions potential neigh%list'

      Case (92)

        Write (ounit, '(/,1x,a)') 'error - specified metal potentials have different types'

      Case (93)

        Write (ounit, '(/,1x,a)') 'error - PMFs mixing with rigid bodies not allowed'

      Case (95)

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff (or neigh%cutoff+neigh%padding) > minimum of all half-cell widths'

      Case (96)

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in metal_ld_set_halo'

      Case (97)

        Write (ounit, '(/,1x,a)') 'error - constraints mixing with rigid bodies not allowed'

      Case (99)

        Write (ounit, '(/,1x,a)') 'error - cannot have shells as part of a constraint, rigid body or tether'

      Case (100)

        Write (ounit, '(/,1x,a)') 'error - core-shell unit separation > neigh%cutoff (the system cutoff)'

      Case (101)

        Write (ounit, '(/,1x,a)') 'error - calculated four-body potential index too large'

      Case (102)

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff < 2*tersoffs%cutoff (maximum cutoff for tersoff potentials)'

      Case (103)

        Write (ounit, '(/,1x,a)') 'error - parameter mxlshp exceeded in pass_shared_units'

      Case (104)

        Write (ounit, '(/,1x,a)') 'error - arrays listme and lstout exceeded in pass_shared_units'

      Case (105)

        Write (ounit, '(/,1x,a)') 'error - shake algorithm (constraints_shake) failed to converge'

      Case (106)

        Write (ounit, '(/,1x,a)') 'error - neighbour list array too small in link_cell_pairs'

      Case (107)

        Write (ounit, '(/,1x,a)') 'error - too many pairs for rdf%rdf look up specified'

      Case (108)

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in rdf%rdf look up neigh%list'

      Case (109)

        Write (ounit, '(/,1x,a)') 'error - calculated pair rdf%rdf index too large'

      Case (110)

        Write (ounit, '(/,1x,a)') 'error - duplicate rdf%rdf look up pair specified'

      Case (111)

        Write (ounit, '(/,1x,a)') 'error - bond constraint separation > neigh%cutoff (the system cutoff)'

      Case (112)

        Write (ounit, '(/,1x,a)') 'error - only one *constraints* directive per molecule is allowed'

      Case (113)

        Write (ounit, '(/,1x,a)') 'error - intra-molecular bookkeeping arrays exceeded in deport_atomic_data'

      Case (114)

        Write (ounit, '(/,1x,a)') 'error - legend array exceeded in deport_atomic_data'

      Case (115)

        Write (ounit, '(/,1x,a)') 'error - transfer buffer exceeded in update_shared_units'

      Case (116)

        Write (ounit, '(/,1x,a)') 'error - incorrect atom transfer in update_shared_units'

      Case (118)

        Write (ounit, '(/,1x,a)') 'error - construction error in pass_shared_units'

      Case (120)

        Write (ounit, '(/,1x,a)') 'error - invalid determinant in matrix inversion'

      Case (122)

        Write (ounit, '(/,1x,a)') 'error - FIELD file not found'

      Case (124)

        Write (ounit, '(/,1x,a)') 'error - CONFIG file not found'

      Case (126)

        Write (ounit, '(/,1x,a)') 'error - CONTROL file not found'

      Case (128)

        Write (ounit, '(/,1x,a)') 'error - chemical bond unit separation > neigh%cutoff (the system cutoff)'

      Case (130)

        Write (ounit, '(/,1x,a)') 'error - bond angle unit diameter > neigh%cutoff (the system cutoff)'

      Case (132)

        Write (ounit, '(/,1x,a)') 'error - dihedral angle unit diameter > neigh%cutoff (the system cutoff)'

      Case (134)

        Write (ounit, '(/,1x,a)') 'error - inversion angle unit diameter > neigh%cutoff (the system cutoff)'

      Case (138)

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in refresh_halo_positions'

      Case (141)

        Write (ounit, '(/,1x,a)') 'error - duplicate metal potential specified'

      Case (150)

        Write (ounit, '(/,1x,a)') 'error - unknown van der waals potential selected'

      Case (151)

        Write (ounit, '(/,1x,a)') 'error - unknown EAM keyword in TABEAM'

      Case (152)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to dpd_v_export'

      Case (154)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer exceeded in dpd_v_export'

      Case (156)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer exceeds limit in dpd_v_export'

      Case (158)

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in dpd_v_set_halo'

      Case (160)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to statistics_connect_spread'

      Case (163)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in statistics_connect_spread'

      Case (164)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in statistics_connect_spread'

      Case (170)

        Write (ounit, '(/,1x,a)') 'error - too many variables for statistics array'

      Case (172)

        Write (ounit, '(/,1x,a)') 'error - duplicate intra-molecular entries specified in TABBND||TABANG||TABDIH||TABINV'

      Case (174)

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in mpoles_rotmat_set_halo'

      Case (176)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to mpoles_rotmat_export'

      Case (178)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in mpoles_rotmat_export'

      Case (180)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in mpoles_rotmat_export'

      Case (200)

        Write (ounit, '(/,1x,a)') 'error - rdf%rdf||z-density||bond%dst||angle%dst||dihedral%dst||inversion%dst'// &
          'buffer array too small in system_revive'

      Case (210)

        Write (ounit, '(/,1x,a)') 'error - only one *angles* directive per molecule is allowed'

      Case (220)

        Write (ounit, '(/,1x,a)') 'error - only one *dihedrals* directive per molecule is allowed'

      Case (230)

        Write (ounit, '(/,1x,a)') 'error - only one *inversions* directive per molecule is allowed'

      Case (240)

        Write (ounit, '(/,1x,a)') 'error - only one *tethers* directive per molecule is allowed'

      Case (300)

        Write (ounit, '(/,1x,a)') 'error - incorrect boundary condition for link-cell algorithms'

      Case (305)

        Write (ounit, '(/,1x,a)') 'error - too few link cells per dimension for many-body and tersoff forces subroutines'

      Case (307)

        Write (ounit, '(/,1x,a)') 'error - link cell algorithm violation'

      Case (308)

        Write (ounit, '(/,1x,a)') 'error - link cell algorithm in contention with SPME sum precision'

      Case (321)

        Write (ounit, '(/,1x,a)') 'error - LFV quaternion integrator failed'

      Case (340)

        Write (ounit, '(/,1x,a)') 'error - invalid integration option requested'

      Case (350)

        Write (ounit, '(/,1x,a)') 'error - too few degrees of freedom'

      Case (360)

        Write (ounit, '(/,1x,a)') 'error - degrees of freedom distribution problem'

      Case (380)

        Write (ounit, '(/,1x,a)') 'error - simulation temperature not specified or < 1 K'

      Case (381)

        Write (ounit, '(/,1x,a)') 'error - simulation timestep not specified'

      Case (382)

        Write (ounit, '(/,1x,a)') 'error - simulation cutoff not specified'

      Case (387)

        Write (ounit, '(/,1x,a)') 'error - system pressure not specified'

      Case (390)

        Write (ounit, '(/,1x,a)') 'error - npt/nst ensemble requested in non-periodic system'

      Case (402)

        Write (ounit, '(/,1x,a)') 'error - van der waals cutoff not specified'

      Case (410)

        Write (ounit, '(/,1x,a)') 'error - cell not consistent with image convention'

      Case (414)

        Write (ounit, '(/,1x,a)') 'error - conflicting ensemble options in CONTROL file'

      Case (416)

        Write (ounit, '(/,1x,a)') 'error - conflicting force options in CONTROL file'

      Case (430)

        Write (ounit, '(/,1x,a)') 'error - integration routine not available'

      Case (432)

        Write (ounit, '(/,1x,a)') 'error -  undefined tersoff potential'

      Case (433)

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff MUST be specified for the Ewald sum precision'

      Case (436)

        Write (ounit, '(/,1x,a)') 'error - unrecognised ensemble'

      Case (440)

        Write (ounit, '(/,1x,a)') 'error - undefined angular potential'

      Case (442)

        Write (ounit, '(/,1x,a)') 'error - undefined three-body potential'

      Case (443)

        Write (ounit, '(/,1x,a)') 'error - undefined four-body potential'

      Case (444)

        Write (ounit, '(/,1x,a)') 'error - undefined bond potential'

      Case (445)

        Write (ounit, '(/,1x,a)') 'error - r_14 > neigh%cutoff in dihedrals_forces'

      Case (446)

        Write (ounit, '(/,1x,a)') 'error - undefined electrostatic key in dihedrals_forces'

      Case (448)

        Write (ounit, '(/,1x,a)') 'error - undefined dihedral potential'

      Case (449)

        Write (ounit, '(/,1x,a)') 'error - undefined inversion potential'

      Case (450)

        Write (ounit, '(/,1x,a)') 'error - undefined tethering potential'

      Case (451)

        Write (ounit, '(/,1x,a)') 'error - three-body potential cutoff undefined'

      Case (452)

        Write (ounit, '(/,1x,a)') 'error - undefined pair potential'

      Case (453)

        Write (ounit, '(/,1x,a)') 'error - four-body potential cutoff undefined'

      Case (454)

        Write (ounit, '(/,1x,a)') 'error - undefined external field'

      Case (456)

        Write (ounit, '(/,1x,a)') 'error - external field xpis-ton is applied to a layer with at least one frozen particle'

      Case (461)

        Write (ounit, '(/,1x,a)') 'error - undefined metal potential'

      Case (462)

        Write (ounit, '(/,1x,a)') 'error - thermostat friction constant MUST be > 0'

      Case (463)

        Write (ounit, '(/,1x,a)') 'error - barostat friction constant MUST be > 0'

      Case (464)

        Write (ounit, '(/,1x,a)') 'error - thermostat relaxation time constant MUST be > 0'

      Case (466)

        Write (ounit, '(/,1x,a)') 'error - barostat relaxation time constant MUST be > 0'

      Case (467)

        Write (ounit, '(/,1x,a)') 'error - rho MUST not be zero in valid buckingham potential'

      Case (468)

        Write (ounit, '(/,1x,a)') 'error - r0 too large for snm potential with current cutoff'

      Case (470)

        Write (ounit, '(/,1x,a)') 'error - n < m in definition of n-m potential'

      Case (471)

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff < 2*rctbp (maximum cutoff for three-body potentials)'

      Case (472)

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff < 2*fourbody%cutoff (maximum cutoff for four-body potentials)'

      Case (474)

        Write (ounit, '(/,1x,a)') 'error - conjugate gradient mimimiser cycle limit exceeded'

      Case (476)

        Write (ounit, '(/,1x,a)') 'error - shells MUST all HAVE either zero or non-zero masses'

      Case (477)

        Write (ounit, '(/,1x,a)') 'error - only one *shells* directive per molecule is allowed'

      Case (478)

        Write (ounit, '(/,1x,a)') 'error - shake algorithms (constraints & pmf) failed to converge'

      Case (480)

        Write (ounit, '(/,1x,a)') 'error - PMF length > minimum of all half-cell widths'

      Case (484)

        Write (ounit, '(/,1x,a)') 'error - only one potential of mean force permitted'

      Case (486)

        Write (ounit, '(/,1x,a)') 'error - only one of the PMF units is permitted to have frozen atoms'

      Case (488)

        Write (ounit, '(/,1x,a)') 'error - too many PMF constraints per domain'

      Case (490)

        Write (ounit, '(/,1x,a)') 'error - local PMF constraint not found locally'

      Case (492)

        Write (ounit, '(/,1x,a)') 'error - a diameter of a PMF unit > minimum of all half cell widths'

      Case (494)

        Write (ounit, '(/,1x,a)') 'error - overconstrained PMF units'

      Case (497)

        Write (ounit, '(/,1x,a)') 'error - pmf_quench failure'

      Case (498)

        Write (ounit, '(/,1x,a)') 'error - shake algorithm (pmf_shake) failed to converge'

      Case (499)

        Write (ounit, '(/,1x,a)') 'error - rattle algorithm (pmf_rattle) failed to converge'

      Case (500)

        Write (ounit, '(/,1x,a)') 'error - PMF unit of zero length is not permitted'

      Case (501)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in PMF unit'

      Case (502)

        Write (ounit, '(/,1x,a)') 'error - coincidence of units in PMF constraint'

      Case (504)

        Write (ounit, '(/,1x,a)') 'error - cutoff too large for TABLE file'

      Case (505)

        Write (ounit, '(/,1x,a)') 'error - EAM metal densities or pair crossfunctions out of range'

      Case (506)

        Write (ounit, '(/,1x,a)') 'error - EAM or MBPC metal densities out of range'

      Case (507)

        Write (ounit, '(/,1x,a)') 'error - metal density embedding out of range'

      Case (508)

        Write (ounit, '(/,1x,a)') 'error - EAM metal interaction entry in TABEAM unspecified in FIELD'

      Case (509)

        Write (ounit, '(/,1x,a)') 'error - duplicate entry for a pair interaction detected in TABEAM'

      Case (510)

        Write (ounit, '(/,1x,a)') 'error - duplicate entry for a density function detected in TABEAM'

      Case (511)

        Write (ounit, '(/,1x,a)') 'error - duplicate entry for an embedding function detected in TABEAM'

      Case (512)

        Write (ounit, '(/,1x,a)') 'error - non-definable vdw/dpd interactions detected in FIELD'

      Case (513)

        Write (ounit, '(/,1x,a)') 'error - particle assigned to non-existent domain in read_config'

      Case (514)

        Write (ounit, '(/,1x,a)') 'error - allowed image conventions are: 0, 1, 2, 3 and 6'

      Case (515)

        Write (ounit, '(/,1x,a)') 'error - rattle algorithm (constraints_rattle) failed to converge'

      Case (516)

        Write (ounit, '(/,1x,a)') 'error - the number nodes MUST be a power of 2 series number'

      Case (517)

        Write (ounit, '(/,1x,a)') 'error - allowed configuration information levels are: 0, 1 and 2'

      Case (518)

        Write (ounit, '(/,1x,a)') 'error - control distances for variable timestep not intact'

      Case (519)

        Write (ounit, '(/,1x,a)') 'error - REVOLD is incompatible or does not exist'

      Case (520)

        Write (ounit, '(/,1x,a)') 'error - domain decomposition failed'

      Case (530)

        Write (ounit, '(/,2(1x,a,/))') 'error - pseudo thermostat thickness MUST comply with', &
          '2 Angs <= thickness < a quarter of the minimum MD cell width'

      Case (540)

        Write (ounit, '(/,2(1x,a,/))') 'error - pseudo thermostat can ONLY be used in bulk simulations', &
          'i.e. imcon MUST be 1, 2 or 3'

      Case (551)

        Write (ounit, '(/,1x,a)') 'error - REFERENCE not found !!!'

      Case (552)

        Write (ounit, '(/,1x,a)') 'error - REFERENCE MUST contain cell parameters !!!'

      Case (553)

        Write (ounit, '(/,1x,a)') 'error - REFERENCE is inconsistent !!!'

      Case (554)

        Write (ounit, '(/,1x,a)') "error - REFERENCE's format different from CONFIG's !!!"

      Case (555)

        Write (ounit, '(/,1x,a)') 'error - particle assigned to non-existent domain in defects_read_reference'

      Case (556)

        Write (ounit, '(/,1x,a)') 'error - too many atoms in REFERENCE file'

      Case (557)

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to defects_reference_export'

      Case (558)

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer exceeded in defects_reference_export'

      Case (559)

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in defects_reference_export'

      Case (560)

        Write (ounit, '(/,1x,a)') 'error - rdef found to be > half the shortest interatomic distance in REFERENCE'

      Case (570)

        Write (ounit, '(/,1x,a)') 'error - unsupported image convention (0) for system expansion option nfold'

      Case (580)

        Write (ounit, '(/,1x,a)') 'error - replay (HISTORY) option can only be used for structural property recalculation'

      Case (585)

        Write (ounit, '(/,1x,a)') 'error - HISTORY file does not exist'

      Case (590)

        Write (ounit, '(/,1x,a)') 'error - uknown minimisation type, only "force", "energy" and "distance" are recognised'

      Case (600)

        Write (ounit, '(/,1x,a)') 'error - "impact" option specified more than once in CONTROL'

      Case (610)

        Write (ounit, '(/,1x,a)') &
          'error - "impact" applied on particle that is either frozen, or the shell of a core-shell unit or part of a RB'

      Case (615)

        Write (ounit, '(/,1x,a)') 'error - q(core)*q(shell)*k(core-shell) MUST NOT be zero'

      Case (620)

        Write (ounit, '(/,1x,a)') 'error - duplicate or mixed intra-molecular entries specified in FIELD'

      Case (623)

        Write (ounit, '(/,1x,a)') "error - MPOLES's molecular data mismatched with respect to FIELD's data"

      Case (625)

        Write (ounit, '(/,1x,a)') 'error - only one *rigid* directive per molecule is allowed'

      Case (630)

        Write (ounit, '(/,1x,a)') 'error - too many rigid body units specified'

      Case (632)

        Write (ounit, '(/,1x,a)') 'error - rigid body unit MUST have at least 2 sites'

      Case (634)

        Write (ounit, '(/,1x,a)') 'error - rigid body unit MUST have at least one non-massless site'

      Case (636)

        Write (ounit, '(/,1x,a)') 'error - rigid body unit MUST NOT have any frozen site'

      Case (638)

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in a rigid body unit'

      Case (640)

        Write (ounit, '(/,1x,a)') 'error - too many rigid body units per domain'

      Case (642)

        Write (ounit, '(/,1x,a)') 'error - rigid body unit diameter > neigh%cutoff (the system cutoff)'

      Case (644)

        Write (ounit, '(/,1x,a)') 'error - overconstrained rigid body unit'

      Case (646)

        Write (ounit, '(/,1x,a)') 'error - overconstrained constraint unit'

      Case (648)

        Write (ounit, '(/,1x,a)') 'error - quaternion setup failed'

      Case (650)

        Write (ounit, '(/,1x,a)') 'error - failed to find principal axis system'

      Case (655)

        Write (ounit, '(/,1x,a)') 'error - FENE bond breaking failure'

      Case (660)

        Write (ounit, '(/,1x,a)') 'error - bond length > cutoff in TABBND or cutoff for PDF collection'

      Case (670)

        Write (ounit, '(/,1x,a)') 'error - insufficient electronic temperature cells for TTM heat diffusion'

      Case (680)

        Write (ounit, '(/,1x,a)') 'error - rpad too large for calculation of ionic temperatures'

      Case (681)

        Write (ounit, '(/,1x,a)') 'error - electronic specific heat not fully specified'

      Case (682)

        Write (ounit, '(/,1x,a)') 'error - thermal conductivity of metal not specified'

      Case (683)

        Write (ounit, '(/,1x,a)') 'error - thermal diffusivity of non-metal not specified'

      Case (684)

        Write (ounit, '(/,1x,a)') 'error - cannot find or open thermal conductivity table file (Ke.dat)'

      Case (685)

        Write (ounit, '(/,1x,a)') 'error - no data found in thermal conductivity table file (Ke.dat)'

      Case (686)

        Write (ounit, '(/,1x,a)') 'error - cannot find or open volumetric heat capacity table file (Ce.dat)'

      Case (687)

        Write (ounit, '(/,1x,a)') 'error - no data found in volumetric heat capacity table file (Ce.dat)'

      Case (688)

        Write (ounit, '(/,1x,a)') 'error - cannot find or open thermal diffusivity table file (De.dat)'

      Case (689)

        Write (ounit, '(/,1x,a)') 'error - no data found in thermal diffusivity table file (De.dat)'

      Case (690)

        Write (ounit, '(/,1x,a)') 'error - cannot find or open coupling constant table file (g.dat)'

      Case (691)

        Write (ounit, '(/,1x,a)') 'error - no data found in coupling constant table file (g.dat)'

      Case (692)

        Write (ounit, '(/,1x,a)') 'error - end of file encountered in table file (Ke.dat, Ce.dat, De.dat or g.dat)'

      Case (693)

        Write (ounit, '(/,1x,a)') 'error - negative electronic temperature: instability in electronic heat diffusion equation'

      Case (694)

        Write (ounit, '(/,1x,a)') 'error - electronic temperature restart file (DUMP_E) does not exist'

      Case (695)

        Write (ounit, '(/,1x,a)') &
          'error - mismatch in electronic temperature lattice sizes between restart (DUMP_E) and CONTROL files'

      Case (696)

        Write (ounit, '(/,1x,a)') 'error - cannot read electronic temperature restart (DUMP_E) file'

      Case (1000)

        Write (ounit, '(/,1x,a)') 'error - working precision mismatch between FORTRAN90 and MPI implementation'

      Case (1001)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gcheck_vector'

      Case (1002)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gcheck_vector'

      Case (1003)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gisum_vector'

      Case (1004)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gisum_vector'

      Case (1005)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grsum_vector'

      Case (1006)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grsum_vector'

      Case (1007)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gimax_vector'

      Case (1008)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gimax_vector'

      Case (1009)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grmax_vector'

      Case (1010)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grmax_vector'

      Case (1011)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in setup -> get_record'

      Case (1012)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in setup -> get_record'

      Case (1013)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angles_arrays'

      Case (1014)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bonds_arrays'

      Case (1015)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in core_shell -> allocate_core_shell_arrays'

      Case (1016)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in statistics -> allocate_statitics_arrays'

      Case (1017)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in tethers_module -> allocate_tethers_arrays'

      Case (1018)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in constraints -> allocate_constraints_arrays'

      Case (1019)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in external_field_module -> allocate_external_field_arrays'

      Case (1020)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dihedrals -> allocate_dihedrals_arrays'

      Case (1021)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in inversions -> allocate_inversion_arrays'

      Case (1022)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in vdw -> allocate_vdw_arrays'

      Case (1023)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_arrays'

      Case (1024)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in three_body_module -> allocate_three_body_arrays'

      Case (1025)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in configuration -> allocate_config_arrays'

      Case (1026)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in site -> allocate_site_arrays'

      Case (1027)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in tersoff_module -> allocate_tersoff_arrays'

      Case (1030)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in core_shell -> deallocate_core_shell_arrays'

      Case (1031)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in tethers_module -> deallocate_tethers_arrays'

      Case (1032)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in constraints -> deallocate_constraints_arrays'

      Case (1033)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in dihedrals -> deallocate_dihedrals_arrays'

      Case (1034)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in inversions -> deallocate_inversions_arrays'

      Case (1035)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in defects_module -> allocate_defects_arrays'

      Case (1036)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in pmf -> allocate_pmf_arrays'

      Case (1037)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in pmf -> deallocate_pmf_arrays'

      Case (1038)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in minimise_module -> allocate_minimise_arrays'

      Case (1040)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ewald_module -> ewald_allocate_kall_arrays'

      Case (1041)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in langevin_module -> langevin_allocate_arrays'

      Case (1042)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in rigid_bodies -> allocate_rigid_bodies_arrays'

      Case (1043)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in rigid_bodies -> deallocate_rigid_bodies_arrays'

      Case (1044)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gimin_vector'

      Case (1045)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gimin_vector'

      Case (1046)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grmin_vector'

      Case (1047)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grmin_vector'

      Case (1048)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grsum_matrix'

      Case (1049)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grsum_matrix'

      Case (1050)

        Write (ounit, '(/,1x,a)') 'error - sorted I/O base communicator not set'

      Case (1053)

        Write (ounit, '(/,1x,a)') 'error - sorted I/O allocation error'

      Case (1056)

        Write (ounit, '(/,1x,a)') 'error - unkown write option given to sorted I/O'

      Case (1059)

        Write (ounit, '(/,1x,a)') 'error - unknown write level given to sorted I/O'

      Case (1060)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in statistics -> allocate_statitics_connect'

      Case (1061)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in statistics -> deallocate_statitics_connect'

      Case (1063)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in vdw -> allocate_vdw_table_arrays'

      Case (1066)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in vdw -> allocate_vdw_direct_fs_arrays'

      Case (1069)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_table_arrays'

      Case (1070)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ewald_module -> ewald_allocate_kfrz_arrays'

      Case (1072)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bond_pot_arrays'

      Case (1073)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bond_dst_arrays'

      Case (1074)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angl_pot_arrays'

      Case (1075)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angl_dst_arrays'

      Case (1076)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dihedrals -> allocate_dihd_pot_arrays'

      Case (1077)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dihedrals -> allocate_dihd_dst_arrays'

      Case (1078)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in inversion_module -> allocate_invr_pot_arrays'

      Case (1079)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in inversion_module -> allocate_invr_dst_arrays'

      Case (1080)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in greenkubo_module -> allocate_greenkubo_arrays'

      Case (1081)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dpd -> allocate_dpd_arrays'

      Case (1082)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_erf_arrays'

      Case (1083)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_module -> allocate_ttm_arrays'

      Case (1084)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_module -> deallocate_ttm_arrays'

      Case (1085)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_ion_temperature'

      Case (1086)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_ion_temperature'

      Case (1087)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_thermal_diffusion'

      Case (1088)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_thermal_diffusion'

      Case (1089)

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_track_module -> depoinit'

      Case (1090)

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_track_module -> depoevolve'

      Case Default

        Write (ounit, '(/,1x,a)') 'error - see message above'

      End select

    End If

    ! close all I/O channels
    Call close_unit(ounit)

    ! abort comms

    Call abort_comms(eworld, kode)

  End Subroutine error

  Subroutine error_alloc(array, routine)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for printing standard message for allocation errors
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins october 2018
    !!
    !!----------------------------------------------------------------------!

    Character(Len=*) :: array, routine

    Write (ounit, '(/,1x,a)') 'error - allocation failure in '//Trim(routine)//' -> '//Trim(array)

    Call close_unit(ounit)
    Call abort_comms(eworld, 1001)

  End Subroutine error_alloc

  Subroutine error_dealloc(array, routine)

    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for printing standard message for deallocation errors
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins october 2018
    !!
    !!----------------------------------------------------------------------!
    Character(Len=*) :: array, routine

    Write (ounit, '(/,1x,a)') 'error - deallocation failure in '//Trim(routine)//' -> '//Trim(array)

    Call close_unit(ounit)
    Call abort_comms(eworld, 1002)

  End Subroutine error_dealloc

  Subroutine error_read(ierr, routine, break_eor, break_end)
    !!----------------------------------------------------------------------!
    !!
    !! dl_poly_4 subroutine for checking if a read was successful
    !!
    !! copyright - daresbury laboratory
    !! author    - j.s.wilkins october 2018
    !!
    !!----------------------------------------------------------------------!
    Integer, Intent(In) :: ierr
    Character(Len=*) :: routine
    Logical, Optional, Intent(Out) :: break_eor
    Logical, Optional, Intent(Out) :: break_end

    if (present(break_eor)) then
       break_eor = .false.
    end if
    if (present(break_end)) then
       break_end = .false.
    end if

    If (ierr == IOSTAT_END) then
       if (present(break_end)) then
          break_end = .true.
       else
          Call error(0, 'Unexpected end of file in '//trim(routine))
       end if
    Else if (ierr == IOSTAT_EOR) then
       if (present(break_eor)) then
          break_eor = .true.
       else
          Call error(0, 'Unexpected end of record in '//trim(routine))
       end if
    else if (ierr /= 0) then
       Call error(0, 'Unknown error in read in '//trim(routine))
    end If

  end Subroutine error_read

  !> Close all open file units
  Subroutine close_unit(i)
    Integer, Intent(InOut) :: i

    Integer :: ierr
    Logical :: has_name, is_open

    Inquire (i, opened=is_open, named=has_name, iostat=ierr)

    If (is_open .and. has_name .and. ierr == 0 .and. All(i /= [-1, ERROR_UNIT, INPUT_UNIT, OUTPUT_UNIT])) Then
      Close (i)
    End If
  End Subroutine close_unit

End Module errors_warnings
