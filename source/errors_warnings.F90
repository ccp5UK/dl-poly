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
                                           output_unit
  Use kinds,                         Only: wp

  Implicit None

  Private

  Type(comms_type), Save :: eworld
  Integer, Save :: ounit

  Public :: warning
  Public :: error
  Public :: info
  Public :: init_error_system
  Public :: error_alloc, error_dealloc

  Interface warning
    Module Procedure warning_special
    Module Procedure warning_general
  End Interface

  Interface info
    Module Procedure info_sl
    Module Procedure info_ml
  End Interface

Contains

  Subroutine init_error_system(nrite, comm)

    Integer,          Intent(In   ) :: nrite
    Type(comms_type), Intent(In   ) :: comm

    eworld%comm = comm%comm
    eworld%idnode = comm%idnode
    ounit = nrite

  End Subroutine init_error_system

  Subroutine warning_special(kode, a, b, c)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for printing warning messages and returning
    ! control back to the main program
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

    Integer :: ia, ib, ic

    If (eworld%idnode == 0) Then

      Write (ounit, '(/,1x,a,i6)') 'warning issued ', kode

      If (kode == 2) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - DD with ', ia, ' idle nodes, (out of ', ib, ') mapped on vacuum !!! ***'

      Else If (kode == 3) Then

        Write (ounit, '(/,1x,3(a,f7.3),a,/)') &
          '*** warning - DD cutoff(+padding) is: ', a, ' Angstroms while minimum half-cell width is: ', b, ' Angstroms !!! ***'

      Else If (kode == 4) Then

        Write (ounit, '(/,2(1x,a,/))') &
          '*** warning - system with uncharged particles !!! ***', &
          '*** "no elec" or/and "no strict" directives in CONTROL may speed up simulation !!! ***'

      Else If (kode == 5) Then

        Write (ounit, '(/,1x,a,f12.4,a,/)') &
          '*** warning - non-zero total system charge: ', a, ' positrons !!! ***'

      Else If (kode == 6) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/,1x,a,a,/,1x,a)') &
          '*** warning - maximum length of linked cell list: ', ia, &
          ' + 1 is less than maximum length of particle exclusion list: ', &
          ib, ' !!! ***', &
          '*** this may be due to using too short a cutoff in CONTROL ', &
          'and/or a badly defined intramolecular topology in FIELD !!! ***', &
          '*** further potential problems may be expected !!! ***'

      Else If (kode == 7) Then

        Write (ounit, '(/,1x,a,f7.3,a,/,1x,a,/)') &
          '*** warning - DD cutoff is ', a, ' Angstroms !!! ***', &
          '*** Fennell damping is not recommended for cutoffs shorther than 10-12 Angstroms !!! ***'

      Else If (kode == 8) Then

        Write (ounit, '(/,1x,a,2(f7.3,a),/)') &
          '*** warning - : detected maximum rigid body width: ', a, ' Angstroms while the DD cutoff is ', b, ' Angstroms !!! ***'

      Else If (kode == 10) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - no pair forces in use !!! ***'

      Else If (kode == 20) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - : 1..4 scale factors reset for molecule: ', ia, ' sites: ', ib, ' & ', ic, ' !!! ***'

      Else If (kode == 22) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - : 1..4 scale factors reset for dihedrals in the system !!! ***'

      Else If (kode == 30) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - Electrostatics requested in a non-periodic system !!! ***'

      Else If (kode == 34) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - redundant directive prim(ary cutoff) ignored !!! ***'

      Else If (kode == 35) Then

        Write (ounit, '(3(/,1x,a),/)') &
          "*** warning - DL_POLY_2/Classic directive 'delr - Verlet shell strip cutoff' defaulted to ***", &
          "*** DL_POLY_4 directive 'neigh%padding - real space cutoff padding option' ***", &
          "*** neigh%padding=Max(neigh%padding,delr/4) !!! ***"

      Else If (kode == 36) Then

        Write (ounit, '(2(/,1x,a),/)') &
          "*** warning - DL_POLY_2/Classic directive 'mult(iple timestep)' defaulted to ***", &
          "*** DL_POLY_4 directive 'infrequent k-space SPME evaluation' !!! ***"

      Else If (kode == 37) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - redundant directive all pairs ignored !!! ***'

      Else If (kode == 38) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - redundant directive no link ignored !!! ***'

      Else If (kode == 40) Then

        Write (ounit, '(/,1x,a,f7.3,a,/)') '*** warning - tentative cutoff reset to ', a, ' Angstroms !!! ***'

      Else If (kode == 50) Then

        Write (ounit, '(2(/,1x,a),f7.3,a,/)') &
          '*** warning - short-range interactions cutoff reset ***', '*** new cutoff radius (rvdw) ', a, ' Angstroms !!! ***'

      Else If (kode == 60) Then

        Write (ounit, '(/,1x,a,f7.3,a,/)') '*** warning - total system charge ', a, ' positrons !!! ***'

      Else If (kode == 70) Then

        Write (ounit, '(/,1x,a,f7.3,a,/)') '*** warning - switching length reset to ', a, ' Angstroms !!! ***'

      Else If (kode == 80) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - requested thermostat unavailable !!! ***'

      Else If (kode == 90) Then

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

      Else If (kode == 100) Then

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - primary link cell algorithm has a link cell dimension that is < 3 !!! ***', &
          '*** DL_POLY_4 RUNNING IN LOW EFFICIENCY MODE !!! ***'

      Else If (kode == 110) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(2(/,1x,a),2(i0,a),/)') &
          '*** warning - image convention incompatible with the set NsT ensemble ***', &
          '*** config%imcon reset from ', ia, ' to ', ib, ' !!! ***'

      Else If (kode == 120) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - unspecified atom-atom interactions set to zero !!! ***'

      Else If (kode == 130) Then

        Write (ounit, '(/,1x,a,/)') '*** warning - no ensemble is specified !!! ***'

      Else If (kode == 135) Then

        Write (ounit, '(/,1x,a,/,1x,a,/)') &
          '*** warning - incorrect ensemble specified for two-temperature model !!! ***', &
          '*** replacing ensemble with NVT inhomogeneous Langevin !!! ***'

      Else If (kode == 140) Then

        Write (ounit, '(/,1x,2a,2(f8.5,a),/,1x,a,/)') &
          '*** warning - control distances for variable timestep: ', &
          'Dmin = ', a, ' and Dmax = ', b, ' (Angstroms) !!! ***', &
          '*** do not comply with safty condition: Dmax > 2.5 Dmin > 0 !!! ***'

      Else If (kode == 150) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required export buffer size ', ia, ' and actual: ', ib, ' !!! ***'

      Else If (kode == 160) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required import array size ', ia, ' and actual: ', ib, ' !!! ***'

      Else If (kode == 190) Then

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - REVOLD format mishmash detected (restart requested) !!! ***', &
          '*** restart is abandoned and clean start is assumed !!! ***'

      Else If (kode == 200) Then

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - CONFIG contains positions only !!! ***', &
          '*** clean start is assumed !!! ***'

      Else If (kode == 210) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(3(/,1x,a),2(i0,a),/)') &
          '*** warning - system under great constraint !!! ***', &
          '*** total degrees of freedom <=  total number of particles !!! ***', &
          '*** ', ia, ' <= ', ib, ' !!! ***'

      Else If (kode == 220) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - Ewald sum requested in a non-periodic system !!! ***'

      Else If (kode == 230) Then

        ia = Nint(a)

        Write (ounit, '(/,1x,a,i0,a,/,1x,a,/)') &
          '*** warning - PMF unit ', ia, ' config%weight is detected zero !!! ***', &
          '*** member config%weights defaulted to atom type masses (or units) !!! ***'

      Else If (kode == 240) Then

        ia = Nint(a)

        Write (ounit, '(/,1x,a,i0,a,/)') &
          '*** warning - total number of nodes running in parallel: ', ia, ' !!! ***'

      Else If (kode == 250) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required coulombic exclusion array size ', ia, ' and actual (neigh%max_exclude): ', ib, ' !!! ***'

      Else If (kode == 260) Then

        Write (ounit, '(2(/,1x,a),/)') &
          '*** warning - system volume in non-periodic systems is the MD config%cell volume !!! ***', &
          '*** system pressure is calculated with respect to this volume !!! ***'

      Else If (kode == 270) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required reading buffer size is ', ia, ' and actual ', ib, ' !!! ***'

      Else If (kode == 280) Then

        Write (ounit, '(/,1x,a,/,1x,2(a,f9.3),a,/)') &
          '*** warning - pseudo thermostat cannot be applied for this model system since !!! ***', &
          '*** specified thermostat wall thickness ', a, ' > 1/4 minimum MD cell width ', b, ' (Angstroms) !!! ***'

      Else If (kode == 290) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - required link-config%cell neigh%list size is ', ia, ' and actual (neigh%max_list) ', ib, ' !!! ***'

      Else If (kode == 295) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - PMF unit ', ia, ' and rigid body unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 296) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, &
          ' on molecular species type ', ib, &
          ' with compromised polarisability !!! ***'

      Else If (kode == 297) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and angle unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 298) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and dihedral unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 299) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and inversion unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 300) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and PMF unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 301) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and constraint unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 302) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - core-shell unit ', ia, ' and rigid body unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 303) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a),/)') &
          '*** warning - core-shell unit ', ia, ' and tether unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 304) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - constraint unit ', ia, ' and rigid body unit ', ib, &
          ' on molecular species type ', ic, ' in an illegal configuration !!! ***'

      Else If (kode == 305) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a))') &
          '*** warning - rigid body unit ', ia, ' on molecular species type ', ib, &
          ' forced to freeze !!! ***'

      Else If (kode == 306) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a))') &
          '*** warning - constraint unit ', ia, ' on molecular species type ', ib, &
          ' forced to freeze !!! ***'

      Else If (kode == 307) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - rigid body unit ', ia, ' on molecular species type ', ib, &
          ' set to have type ', ic, ' is problematic !!! ***'

      Else If (kode == 308) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - site ', ia, ' of constraint unit ', ib, &
          ' on molecular species type ', ic, ' is problematic !!! ***'

      Else If (kode == 309) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write (ounit, '(/,1x,a,3(i0,a))') &
          '*** warning - site ', ia, ' , member ', ib, ' of PMF unit ', ic, ' is problematic !!! ***'

      Else If (kode == 310) Then

        Write (ounit, '(/,1x,a,/,1x,a,2(f6.2,a),/)') &
          '*** warning - control distance for defect look-up MUST be in the interval [Min(0.3,rcut/3);Min(1.2,rcut/2)] !!! ***', &
          '*** defects distance condition will default from ', a, ' to ', b, ' (Angstroms) !!! ***'

      Else If (kode == 320) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - REFERENCE not found, CONFIG to be used as a reference for defects detection !!! ***'

      Else If (kode == 330) Then

        Write (ounit, '(/,1x,a,/,1x,a,2(f10.6,a),/)') &
          '*** warning - iteration cycles length limit for conjugate gradient minimisers exceded !!! ***', &
          '*** specified convergence tolerance: ', a, ' , needed one for a pass: ', b, ' !!! ***'

      Else If (kode == 340) Then

        Write (ounit, '(3(/,1x,a),2(f7.4,a),/,1x,a,f7.4,a,/)') &
          '*** warning - inconsistent binsize for spatial distribution functions !!! ***', &
          '*** 1.0E-5 (Angstroms) <= binsize <= neigh%cutoff/4 (Angstroms) !!! ***', &
          '*** 1.0E-5 (Angstroms) <= ', a, ' <= ', b, ' (Angstroms) !!! ***', &
          '*** binsize defaults to ', c, ' (Angstroms) !!! ***'

      Else If (kode == 350) Then

        Write (ounit, '(/,1x,a)') &
          '*** warning - expansion along z axis not allowed for slab geometry, nz defaults to 1 !!! ***'

      Else If (kode == 360) Then

        Write (ounit, '(/,1x,a,2(f12.7,a),3(/,1x,a),/)') &
          '*** warning - minimisation tolerance ', a, ' defaults to ', b, ' !!! ***', &
          '*** force   : 1.0    <= tolerance <= 1000.00, default = 50.00 !!! ***', &
          '*** energy  : 0.0    <  tolerance <=    0.01, default = 0.005 !!! ***', &
          '*** distance: 1.0e-6 <= tolerance <=    0.10, default = 0.005 !!! ***'

      Else If (kode == 370) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),2(/,1x,a),/)') &
          '*** warning - k-space evaluation interval ', ia, ' defaults to ', ib, ' !!! ***', &
          '*** the interval must be a positive integer beteen 1 and 10 !!! ***', &
          '*** values > 10 default to 4, no value or 0 defaults to 1 !!! ***'

      Else If (kode == 380) Then

        ia = Nint(a)

        Write (ounit, '(/,1x,a,i0,a,/)') &
          '*** warning - IMPACT applied before the end of the equlibration period (', ia, ' step) !!! ***'

      Else If (kode == 390) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - core-shell unit mix-up or duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 400) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - constraint or rigid body duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 410) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,i0,a,/)') &
          '*** warning - tethered atom duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 420) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - chemical bond duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 430) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - bond angle duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 440) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - dihedral angle duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 450) Then

        ia = Nint(a)
        ib = Nint(b)

        Write (ounit, '(/,1x,a,2(i0,a),/)') &
          '*** warning - inversion angle duplicate specification (', ia, ',', ib, ') !!! ***'

      Else If (kode == 460) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - further semi-isotropic barostat option search abandoned !!! ***'

      Else If (kode == 470) Then

        Write (ounit, '(/,1x,a,/,1x,a,2(f5.2,a),/)') &
          '*** warning - control distance for diplacement qualification MUST be >= 0.25 Angstroms !!! ***', &
          '*** displacements distance condition will default from ', a, ' to ', b, ' Angstroms !!! ***'

      Else If (kode == 480) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - "metal direct" option disabled as incompatible with EAM potentials (TABEAM) !!! ***'

      Else If (kode == 490) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - "metal sqrtrho" option disabled as incompatible with FS type potentials (analytic forms) !!! ***'

      Else If (kode == 500) Then

        Write (ounit, '(2(/,1x,a,/))') &
          '*** warning - insufficient electronic temperature cells to allow energy redistribution from inactive cells !!! ***', &
          '*** option switched off - electronic temperature conservation no longer guaranteed !!! ***'

      Else If (kode == 510) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - no laser energy deposition for two-temperature model specified !!! ***'

      Else If (kode == 515) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - no initial stopping power for two-temperature model specified !!! ***'

      Else If (kode == 520) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - mismatch in timestep between REVIVE and provided electronic temperature lattice !!! ***'

      Else If (kode == 530) Then

        Write (ounit, '(/,1x,a,f6.2,a,/,1x,a,/)') &
          '*** warning - possible time energy deposition discrepancy of at least ', a, '% !!! ***', &
          '*** discrepancy may be due to inactive cells !!! ***'

      Else If (kode == 535) Then

        Write (ounit, '(/,1x,a,/)') &
          '*** warning - spatial energy deposition cutoff in xy-directions greater than available ionic temperature cells !!! ***'

      Else If (kode == 540) Then

        Write (ounit, '(/,1x,a,f6.2,a,/)') &
          '*** warning - possible spatial energy deposition discrepancy of at least ', a, '% !!! ***'

      Else

        Write (ounit, '(/,1x,a,/)') &
          '*** unspecified warning encountered !!! ***'

      End If

    End If

  End Subroutine warning_special

  Subroutine warning_general(message, master_only)
    Character(Len=*),  Intent(In   ) :: message
    Logical, Optional, Intent(In   ) :: master_only

    Logical :: zeroOnly

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

  Subroutine info_ml(message, n, master_only)
    Character(Len=*),  Intent(In   ) :: message(:)
    Integer,           Intent(In   ) :: n
    Logical, Optional, Intent(In   ) :: master_only

    Integer :: i
    Logical :: zeroOnly

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

  Subroutine info_sl(message, master_only)
    Character(Len=*),  Intent(In   ) :: message
    Logical, Optional, Intent(In   ) :: master_only

    Logical :: zeroOnly

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

      If (kode == 1) Then

        Write (ounit, '(/,1x,a)') 'error - word_2_real failure'

      Else If (kode == 2) Then

        Write (ounit, '(/,1x,a)') 'error - too many atom types in FIELD (scan_field)'

      Else If (kode == 3) Then

        Write (ounit, '(/,1x,a)') 'error - unknown directive found in CONTROL file'

      Else If (kode == 4) Then

        Write (ounit, '(/,1x,a)') 'error - unknown directive found in FIELD or MPOLES file'

      Else If (kode == 5) Then

        Write (ounit, '(/,1x,a)') 'error - unknown energy unit requested'

      Else If (kode == 6) Then

        Write (ounit, '(/,1x,a)') 'error - energy unit not specified'

      Else If (kode == 7) Then

        Write (ounit, '(/,1x,a)') 'error - selected external field incompatible with selected ensemble (NVE only!!!)'

      Else If (kode == 8) Then

        Write (ounit, '(/,1x,a)') 'error - ewald precision MUST be a POSITIVE real number'

      Else If (kode == 9) Then

        Write (ounit, '(/,1x,a)') 'error - ewald sum parameters MUST be well defined'

      Else If (kode == 10) Then

        Write (ounit, '(/,1x,a)') 'error - too many molecular types specified'

      Else If (kode == 11) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate molecule directive in FIELD file'

      Else If (kode == 12) Then

        Write (ounit, '(/,1x,a)') 'error - unknown molecule directive in FIELD or MPOLES file'

      Else If (kode == 13) Then

        Write (ounit, '(/,1x,a)') 'error - molecular species not yet specified'

      Else If (kode == 14) Then

        Write (ounit, '(/,1x,a)') 'error - too many unique atom types specified'

      Else If (kode == 15) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate vdw potential specified'

      Else If (kode == 16) Then

        Write (ounit, '(/,1x,a)') 'error - strange exit from FIELD file processing'

      Else If (kode == 17) Then

        Write (ounit, '(/,1x,a)') 'error - strange exit from CONTROL file processing'

      Else If (kode == 18) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate three-body potential specified'

      Else If (kode == 19) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate four-body potential specified'

      Else If (kode == 20) Then

        Write (ounit, '(/,1x,a)') 'error - too many molecule sites specified'

      Else If (kode == 21) Then

        Write (ounit, '(/,1x,a)') 'error - molecule contains more atoms/sites than declared'

      Else If (kode == 22) Then

        Write (ounit, '(/,1x,a)') 'error - unsuitable radial increment in TABLE||TABBND||TABANG||TABDIH||TABINV file'

      Else If (kode == 23) Then

        Write (ounit, '(/,1x,a)') 'error - incompatible FIELD and TABLE file potentials'

      Else If (kode == 24) Then

        Write (ounit, '(/,1x,a)') 'error - end of file encountered in TABLE||TABEAM||TABBND||TABANG||TABDIH||TABINV file'

      Else If (kode == 25) Then

        Write (ounit, '(/,1x,a)') 'error - wrong atom type found in CONFIG file'

      Else If (kode == 26) Then

        Write (ounit, '(/,1x,a)') 'error - neutral group option now redundant'

      Else If (kode == 27) Then

        Write (ounit, '(/,1x,a)') "error - unit's member indexed outside molecule's site range"

      Else If (kode == 28) Then

        Write (ounit, '(/,1x,a)') 'error - wrongly indexed atom entries found in CONFIG file'

      Else If (kode == 30) Then

        Write (ounit, '(/,1x,a)') 'error - too many chemical bonds specified'

      Else If (kode == 31) Then

        Write (ounit, '(/,1x,a)') 'error - too many chemical bonds per domain'

      Else If (kode == 32) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in core-shell unit'

      Else If (kode == 33) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in constraint bond unit'

      Else If (kode == 34) Then

        Write (ounit, '(/,1x,a)') 'error - length of constraint bond unit >= real space cutoff (neigh%cutoff)'

      Else If (kode == 35) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in chemical bond unit'

      Else If (kode == 36) Then

        Write (ounit, '(/,1x,a)') 'error - only one *bonds* directive per molecule is allowed'

      Else If (kode == 38) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in metal_ld_export'

      Else If (kode == 39) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in metal_ld_export'

      Else If (kode == 40) Then

        Write (ounit, '(/,1x,a)') 'error - too many bond constraints specified'

      Else If (kode == 41) Then

        Write (ounit, '(/,1x,a)') 'error - too many bond constraints per domain'

      Else If (kode == 42) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to deport_atomic_data'

      Else If (kode == 43) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer exceeded in deport_atomic_data'

      Else If (kode == 44) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in deport_atomic_data'

      Else If (kode == 45) Then

        Write (ounit, '(/,1x,a)') 'error - too many atoms in CONFIG file or per domain'

      Else If (kode == 46) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to export_atomic_data/positions'

      Else If (kode == 47) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to metal_ld_export'

      Else If (kode == 48) Then

        Write (ounit, '(/,1x,a)') 'error - transfer buffer too small in *_table_read'

      Else If (kode == 49) Then

        Write (ounit, '(/,1x,a)') 'error - frozen shell (core-shell unit) specified'

      Else If (kode == 50) Then

        Write (ounit, '(/,1x,a)') 'error - too many bond angles specified'

      Else If (kode == 51) Then

        Write (ounit, '(/,1x,a)') 'error - too many bond angles per domain'

      Else If (kode == 52) Then

        Write (ounit, '(/,1x,a)') 'error - end of FIELD or MPOLES file encountered'

      Else If (kode == 53) Then

        Write (ounit, '(/,1x,a)') 'error - end of CONTROL file encountered'

      Else If (kode == 54) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in export_atomic_data/positions'

      Else If (kode == 55) Then

        Write (ounit, '(/,1x,a)') 'error - end of CONFIG file encountered'

      Else If (kode == 56) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in export_atomic_data/positions'

      Else If (kode == 57) Then

        Write (ounit, '(/,1x,a)') 'error - too many core-shell units specified'

      Else If (kode == 58) Then

        Write (ounit, '(/,1x,a)') 'error - number of atoms in system not conserved'

      Else If (kode == 59) Then

        Write (ounit, '(/,1x,a)') 'error - too many core-shell units per domain'

      Else If (kode == 60) Then

        Write (ounit, '(/,1x,a)') 'error - too many dihedral angles specified'

      Else If (kode == 61) Then

        Write (ounit, '(/,1x,a)') 'error - too many dihedral angles per domain'

      Else If (kode == 62) Then

        Write (ounit, '(/,1x,a)') 'error - too many tethered atoms specified'

      Else If (kode == 63) Then

        Write (ounit, '(/,1x,a)') 'error - too many tethered atoms per domain'

      Else If (kode == 64) Then

        Write (ounit, '(/,1x,a)') 'error - incomplete core-shell unit found in build_book_intra'

      Else If (kode == 65) Then

        Write (ounit, '(/,1x,a)') 'error - too many excluded pairs specified'

      Else If (kode == 66) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in bond angle unit'

      Else If (kode == 67) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in dihedral unit'

      Else If (kode == 68) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in inversion unit'

      Else If (kode == 69) Then

        Write (ounit, '(/,1x,a)') 'error - too many link cells required'

      Else If (kode == 70) Then

        Write (ounit, '(/,1x,a)') 'error - constraint_quench failure'

      Else If (kode == 71) Then

        Write (ounit, '(/,1x,a)') 'error - too many metal potentials specified'

      Else If (kode == 72) Then

        Write (ounit, '(/,1x,a)') 'error - too many tersoff potentials specified'

      Else If (kode == 73) Then

        Write (ounit, '(/,1x,a)') 'error - too many inversion potentials specified'

      Else If (kode == 74) Then

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in tersoff potential neigh%list'

      Else If (kode == 76) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate tersoff potential specified'

      Else If (kode == 77) Then

        Write (ounit, '(/,1x,a)') 'error - too many inversion angles per domain'

      Else If (kode == 79) Then

        Write (ounit, '(/,1x,a)') 'error - tersoff potential cutoff undefined'

      Else If (kode == 80) Then

        Write (ounit, '(/,1x,a)') 'error - too many pair potentials specified'

      Else If (kode == 81) Then

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in pair potential neigh%list'

      Else If (kode == 82) Then

        Write (ounit, '(/,1x,a)') 'error - calculated pair potential index too large'

      Else If (kode == 83) Then

        Write (ounit, '(/,1x,a)') 'error - too many three-body/angles potentials specified'

      Else If (kode == 84) Then

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in three-body/angles potential neigh%list'

      Else If (kode == 85) Then

        Write (ounit, '(/,1x,a)') 'error - required velocities not in CONFIG file'

      Else If (kode == 86) Then

        Write (ounit, '(/,1x,a)') 'error - calculated three-body potential index too large'

      Else If (kode == 88) Then

        Write (ounit, '(/,1x,a)') 'error - legend array exceeded in build_book_intra'

      Else If (kode == 89) Then

        Write (ounit, '(/,1x,a)') 'error - too many four-body/dihedrals/inversions potentials specified'

      Else If (kode == 90) Then

        Write (ounit, '(/,1x,a)') 'error - specified tersoff potentials have different types'

      Else If (kode == 91) Then

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in four-body/dihedrals/iversions potential neigh%list'

      Else If (kode == 92) Then

        Write (ounit, '(/,1x,a)') 'error - specified metal potentials have different types'

      Else If (kode == 93) Then

        Write (ounit, '(/,1x,a)') 'error - PMFs mixing with rigid bodies not allowed'

      Else If (kode == 95) Then

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff (or neigh%cutoff+neigh%padding) > minimum of all half-cell widths'

      Else If (kode == 96) Then

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in metal_ld_set_halo'

      Else If (kode == 97) Then

        Write (ounit, '(/,1x,a)') 'error - constraints mixing with rigid bodies not allowed'

      Else If (kode == 99) Then

        Write (ounit, '(/,1x,a)') 'error - cannot have shells as part of a constraint, rigid body or tether'

      Else If (kode == 100) Then

        Write (ounit, '(/,1x,a)') 'error - core-shell unit separation > neigh%cutoff (the system cutoff)'

      Else If (kode == 101) Then

        Write (ounit, '(/,1x,a)') 'error - calculated four-body potential index too large'

      Else If (kode == 102) Then

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff < 2*tersoffs%cutoff (maximum cutoff for tersoff potentials)'

      Else If (kode == 103) Then

        Write (ounit, '(/,1x,a)') 'error - parameter mxlshp exceeded in pass_shared_units'

      Else If (kode == 104) Then

        Write (ounit, '(/,1x,a)') 'error - arrays listme and lstout exceeded in pass_shared_units'

      Else If (kode == 105) Then

        Write (ounit, '(/,1x,a)') 'error - shake algorithm (constraints_shake) failed to converge'

      Else If (kode == 106) Then

        Write (ounit, '(/,1x,a)') 'error - neighbour list array too small in link_cell_pairs'

      Else If (kode == 107) Then

        Write (ounit, '(/,1x,a)') 'error - too many pairs for rdf%rdf look up specified'

      Else If (kode == 108) Then

        Write (ounit, '(/,1x,a)') 'error - unidentified atom in rdf%rdf look up neigh%list'

      Else If (kode == 109) Then

        Write (ounit, '(/,1x,a)') 'error - calculated pair rdf%rdf index too large'

      Else If (kode == 110) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate rdf%rdf look up pair specified'

      Else If (kode == 111) Then

        Write (ounit, '(/,1x,a)') 'error - bond constraint separation > neigh%cutoff (the system cutoff)'

      Else If (kode == 112) Then

        Write (ounit, '(/,1x,a)') 'error - only one *constraints* directive per molecule is allowed'

      Else If (kode == 113) Then

        Write (ounit, '(/,1x,a)') 'error - intra-molecular bookkeeping arrays exceeded in deport_atomic_data'

      Else If (kode == 114) Then

        Write (ounit, '(/,1x,a)') 'error - legend array exceeded in deport_atomic_data'

      Else If (kode == 115) Then

        Write (ounit, '(/,1x,a)') 'error - transfer buffer exceeded in update_shared_units'

      Else If (kode == 116) Then

        Write (ounit, '(/,1x,a)') 'error - incorrect atom transfer in update_shared_units'

      Else If (kode == 118) Then

        Write (ounit, '(/,1x,a)') 'error - construction error in pass_shared_units'

      Else If (kode == 120) Then

        Write (ounit, '(/,1x,a)') 'error - invalid determinant in matrix inversion'

      Else If (kode == 122) Then

        Write (ounit, '(/,1x,a)') 'error - FIELD file not found'

      Else If (kode == 124) Then

        Write (ounit, '(/,1x,a)') 'error - CONFIG file not found'

      Else If (kode == 126) Then

        Write (ounit, '(/,1x,a)') 'error - CONTROL file not found'

      Else If (kode == 128) Then

        Write (ounit, '(/,1x,a)') 'error - chemical bond unit separation > neigh%cutoff (the system cutoff)'

      Else If (kode == 130) Then

        Write (ounit, '(/,1x,a)') 'error - bond angle unit diameter > neigh%cutoff (the system cutoff)'

      Else If (kode == 132) Then

        Write (ounit, '(/,1x,a)') 'error - dihedral angle unit diameter > neigh%cutoff (the system cutoff)'

      Else If (kode == 134) Then

        Write (ounit, '(/,1x,a)') 'error - inversion angle unit diameter > neigh%cutoff (the system cutoff)'

      Else If (kode == 138) Then

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in refresh_halo_positions'

      Else If (kode == 141) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate metal potential specified'

      Else If (kode == 150) Then

        Write (ounit, '(/,1x,a)') 'error - unknown van der waals potential selected'

      Else If (kode == 151) Then

        Write (ounit, '(/,1x,a)') 'error - unknown EAM keyword in TABEAM'

      Else If (kode == 152) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to dpd_v_export'

      Else If (kode == 154) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer exceeded in dpd_v_export'

      Else If (kode == 156) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer exceeds limit in dpd_v_export'

      Else If (kode == 158) Then

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in dpd_v_set_halo'

      Else If (kode == 160) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to statistics_connect_spread'

      Else If (kode == 163) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in statistics_connect_spread'

      Else If (kode == 164) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in statistics_connect_spread'

      Else If (kode == 170) Then

        Write (ounit, '(/,1x,a)') 'error - too many variables for statistics array'

      Else If (kode == 172) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate intra-molecular entries specified in TABBND||TABANG||TABDIH||TABINV'

      Else If (kode == 174) Then

        Write (ounit, '(/,1x,a)') 'error - incorrect atom totals assignments in mpoles_rotmat_set_halo'

      Else If (kode == 176) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to mpoles_rotmat_export'

      Else If (kode == 178) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer size exceeded in mpoles_rotmat_export'

      Else If (kode == 180) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in mpoles_rotmat_export'

      Else If (kode == 200) Then

        Write (ounit, '(/,1x,a)') 'error - rdf%rdf||z-density||bond%dst||angle%dst||dihedral%dst||inversion%dst'// &
          'buffer array too small in system_revive'

      Else If (kode == 210) Then

        Write (ounit, '(/,1x,a)') 'error - only one *angles* directive per molecule is allowed'

      Else If (kode == 220) Then

        Write (ounit, '(/,1x,a)') 'error - only one *dihedrals* directive per molecule is allowed'

      Else If (kode == 230) Then

        Write (ounit, '(/,1x,a)') 'error - only one *inversions* directive per molecule is allowed'

      Else If (kode == 240) Then

        Write (ounit, '(/,1x,a)') 'error - only one *tethers* directive per molecule is allowed'

      Else If (kode == 300) Then

        Write (ounit, '(/,1x,a)') 'error - incorrect boundary condition for link-cell algorithms'

      Else If (kode == 305) Then

        Write (ounit, '(/,1x,a)') 'error - too few link cells per dimension for many-body and tersoff forces subroutines'

      Else If (kode == 307) Then

        Write (ounit, '(/,1x,a)') 'error - link cell algorithm violation'

      Else If (kode == 308) Then

        Write (ounit, '(/,1x,a)') 'error - link cell algorithm in contention with SPME sum precision'

      Else If (kode == 321) Then

        Write (ounit, '(/,1x,a)') 'error - LFV quaternion integrator failed'

      Else If (kode == 340) Then

        Write (ounit, '(/,1x,a)') 'error - invalid integration option requested'

      Else If (kode == 350) Then

        Write (ounit, '(/,1x,a)') 'error - too few degrees of freedom'

      Else If (kode == 360) Then

        Write (ounit, '(/,1x,a)') 'error - degrees of freedom distribution problem'

      Else If (kode == 380) Then

        Write (ounit, '(/,1x,a)') 'error - simulation temperature not specified or < 1 K'

      Else If (kode == 381) Then

        Write (ounit, '(/,1x,a)') 'error - simulation timestep not specified'

      Else If (kode == 382) Then

        Write (ounit, '(/,1x,a)') 'error - simulation cutoff not specified'

      Else If (kode == 387) Then

        Write (ounit, '(/,1x,a)') 'error - system pressure not specified'

      Else If (kode == 390) Then

        Write (ounit, '(/,1x,a)') 'error - npt/nst ensemble requested in non-periodic system'

      Else If (kode == 402) Then

        Write (ounit, '(/,1x,a)') 'error - van der waals cutoff not specified'

      Else If (kode == 410) Then

        Write (ounit, '(/,1x,a)') 'error - cell not consistent with image convention'

      Else If (kode == 414) Then

        Write (ounit, '(/,1x,a)') 'error - conflicting ensemble options in CONTROL file'

      Else If (kode == 416) Then

        Write (ounit, '(/,1x,a)') 'error - conflicting force options in CONTROL file'

      Else If (kode == 430) Then

        Write (ounit, '(/,1x,a)') 'error - integration routine not available'

      Else If (kode == 432) Then

        Write (ounit, '(/,1x,a)') 'error -  undefined tersoff potential'

      Else If (kode == 433) Then

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff MUST be specified for the Ewald sum precision'

      Else If (kode == 436) Then

        Write (ounit, '(/,1x,a)') 'error - unrecognised ensemble'

      Else If (kode == 440) Then

        Write (ounit, '(/,1x,a)') 'error - undefined angular potential'

      Else If (kode == 442) Then

        Write (ounit, '(/,1x,a)') 'error - undefined three-body potential'

      Else If (kode == 443) Then

        Write (ounit, '(/,1x,a)') 'error - undefined four-body potential'

      Else If (kode == 444) Then

        Write (ounit, '(/,1x,a)') 'error - undefined bond potential'

      Else If (kode == 445) Then

        Write (ounit, '(/,1x,a)') 'error - r_14 > neigh%cutoff in dihedrals_forces'

      Else If (kode == 446) Then

        Write (ounit, '(/,1x,a)') 'error - undefined electrostatic key in dihedrals_forces'

      Else If (kode == 448) Then

        Write (ounit, '(/,1x,a)') 'error - undefined dihedral potential'

      Else If (kode == 449) Then

        Write (ounit, '(/,1x,a)') 'error - undefined inversion potential'

      Else If (kode == 450) Then

        Write (ounit, '(/,1x,a)') 'error - undefined tethering potential'

      Else If (kode == 451) Then

        Write (ounit, '(/,1x,a)') 'error - three-body potential cutoff undefined'

      Else If (kode == 452) Then

        Write (ounit, '(/,1x,a)') 'error - undefined pair potential'

      Else If (kode == 453) Then

        Write (ounit, '(/,1x,a)') 'error - four-body potential cutoff undefined'

      Else If (kode == 454) Then

        Write (ounit, '(/,1x,a)') 'error - undefined external field'

      Else If (kode == 456) Then

        Write (ounit, '(/,1x,a)') 'error - external field xpis-ton is applied to a layer with at least one frozen particle'

      Else If (kode == 461) Then

        Write (ounit, '(/,1x,a)') 'error - undefined metal potential'

      Else If (kode == 462) Then

        Write (ounit, '(/,1x,a)') 'error - thermostat friction constant MUST be > 0'

      Else If (kode == 463) Then

        Write (ounit, '(/,1x,a)') 'error - barostat friction constant MUST be > 0'

      Else If (kode == 464) Then

        Write (ounit, '(/,1x,a)') 'error - thermostat relaxation time constant MUST be > 0'

      Else If (kode == 466) Then

        Write (ounit, '(/,1x,a)') 'error - barostat relaxation time constant MUST be > 0'

      Else If (kode == 467) Then

        Write (ounit, '(/,1x,a)') 'error - rho MUST not be zero in valid buckingham potential'

      Else If (kode == 468) Then

        Write (ounit, '(/,1x,a)') 'error - r0 too large for snm potential with current cutoff'

      Else If (kode == 470) Then

        Write (ounit, '(/,1x,a)') 'error - n < m in definition of n-m potential'

      Else If (kode == 471) Then

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff < 2*rctbp (maximum cutoff for three-body potentials)'

      Else If (kode == 472) Then

        Write (ounit, '(/,1x,a)') 'error - neigh%cutoff < 2*fourbody%cutoff (maximum cutoff for four-body potentials)'

      Else If (kode == 474) Then

        Write (ounit, '(/,1x,a)') 'error - conjugate gradient mimimiser cycle limit exceeded'

      Else If (kode == 476) Then

        Write (ounit, '(/,1x,a)') 'error - shells MUST all HAVE either zero or non-zero masses'

      Else If (kode == 477) Then

        Write (ounit, '(/,1x,a)') 'error - only one *shells* directive per molecule is allowed'

      Else If (kode == 478) Then

        Write (ounit, '(/,1x,a)') 'error - shake algorithms (constraints & pmf) failed to converge'

      Else If (kode == 480) Then

        Write (ounit, '(/,1x,a)') 'error - PMF length > minimum of all half-cell widths'

      Else If (kode == 484) Then

        Write (ounit, '(/,1x,a)') 'error - only one potential of mean force permitted'

      Else If (kode == 486) Then

        Write (ounit, '(/,1x,a)') 'error - only one of the PMF units is permitted to have frozen atoms'

      Else If (kode == 488) Then

        Write (ounit, '(/,1x,a)') 'error - too many PMF constraints per domain'

      Else If (kode == 490) Then

        Write (ounit, '(/,1x,a)') 'error - local PMF constraint not found locally'

      Else If (kode == 492) Then

        Write (ounit, '(/,1x,a)') 'error - a diameter of a PMF unit > minimum of all half cell widths'

      Else If (kode == 494) Then

        Write (ounit, '(/,1x,a)') 'error - overconstrained PMF units'

      Else If (kode == 497) Then

        Write (ounit, '(/,1x,a)') 'error - pmf_quench failure'

      Else If (kode == 498) Then

        Write (ounit, '(/,1x,a)') 'error - shake algorithm (pmf_shake) failed to converge'

      Else If (kode == 499) Then

        Write (ounit, '(/,1x,a)') 'error - rattle algorithm (pmf_rattle) failed to converge'

      Else If (kode == 500) Then

        Write (ounit, '(/,1x,a)') 'error - PMF unit of zero length is not permitted'

      Else If (kode == 501) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in PMF unit'

      Else If (kode == 502) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of units in PMF constraint'

      Else If (kode == 504) Then

        Write (ounit, '(/,1x,a)') 'error - cutoff too large for TABLE file'

      Else If (kode == 505) Then

        Write (ounit, '(/,1x,a)') 'error - EAM metal densities or pair crossfunctions out of range'

      Else If (kode == 506) Then

        Write (ounit, '(/,1x,a)') 'error - EAM or MBPC metal densities out of range'

      Else If (kode == 507) Then

        Write (ounit, '(/,1x,a)') 'error - metal density embedding out of range'

      Else If (kode == 508) Then

        Write (ounit, '(/,1x,a)') 'error - EAM metal interaction entry in TABEAM unspecified in FIELD'

      Else If (kode == 509) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate entry for a pair interaction detected in TABEAM'

      Else If (kode == 510) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate entry for a density function detected in TABEAM'

      Else If (kode == 511) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate entry for an embedding function detected in TABEAM'

      Else If (kode == 512) Then

        Write (ounit, '(/,1x,a)') 'error - non-definable vdw/dpd interactions detected in FIELD'

      Else If (kode == 513) Then

        Write (ounit, '(/,1x,a)') 'error - particle assigned to non-existent domain in read_config'

      Else If (kode == 514) Then

        Write (ounit, '(/,1x,a)') 'error - allowed image conventions are: 0, 1, 2, 3 and 6'

      Else If (kode == 515) Then

        Write (ounit, '(/,1x,a)') 'error - rattle algorithm (constraints_rattle) failed to converge'

      Else If (kode == 516) Then

        Write (ounit, '(/,1x,a)') 'error - the number nodes MUST be a power of 2 series number'

      Else If (kode == 517) Then

        Write (ounit, '(/,1x,a)') 'error - allowed configuration information levels are: 0, 1 and 2'

      Else If (kode == 518) Then

        Write (ounit, '(/,1x,a)') 'error - control distances for variable timestep not intact'

      Else If (kode == 519) Then

        Write (ounit, '(/,1x,a)') 'error - REVOLD is incompatible or does not exist'

      Else If (kode == 520) Then

        Write (ounit, '(/,1x,a)') 'error - domain decomposition failed'

      Else If (kode == 530) Then

        Write (ounit, '(/,2(1x,a,/))') 'error - pseudo thermostat thickness MUST comply with', &
          '2 Angs <= thickness < a quarter of the minimum MD cell width'

      Else If (kode == 540) Then

        Write (ounit, '(/,2(1x,a,/))') 'error - pseudo thermostat can ONLY be used in bulk simulations', &
          'i.e. imcon MUST be 1, 2 or 3'

      Else If (kode == 551) Then

        Write (ounit, '(/,1x,a)') 'error - REFERENCE not found !!!'

      Else If (kode == 552) Then

        Write (ounit, '(/,1x,a)') 'error - REFERENCE MUST contain cell parameters !!!'

      Else If (kode == 553) Then

        Write (ounit, '(/,1x,a)') 'error - REFERENCE is inconsistent !!!'

      Else If (kode == 554) Then

        Write (ounit, '(/,1x,a)') "error - REFERENCE's format different from CONFIG's !!!"

      Else If (kode == 555) Then

        Write (ounit, '(/,1x,a)') 'error - particle assigned to non-existent domain in defects_read_reference'

      Else If (kode == 556) Then

        Write (ounit, '(/,1x,a)') 'error - too many atoms in REFERENCE file'

      Else If (kode == 557) Then

        Write (ounit, '(/,1x,a)') 'error - undefined direction passed to defects_reference_export'

      Else If (kode == 558) Then

        Write (ounit, '(/,1x,a)') 'error - outgoing transfer buffer exceeded in defects_reference_export'

      Else If (kode == 559) Then

        Write (ounit, '(/,1x,a)') 'error - incoming data transfer size exceeds limit in defects_reference_export'

      Else If (kode == 560) Then

        Write (ounit, '(/,1x,a)') 'error - rdef found to be > half the shortest interatomic distance in REFERENCE'

      Else If (kode == 570) Then

        Write (ounit, '(/,1x,a)') 'error - unsupported image convention (0) for system expansion option nfold'

      Else If (kode == 580) Then

        Write (ounit, '(/,1x,a)') 'error - replay (HISTORY) option can only be used for structural property recalculation'

      Else If (kode == 585) Then

        Write (ounit, '(/,1x,a)') 'error - HISTORY file does not exist'

      Else If (kode == 590) Then

        Write (ounit, '(/,1x,a)') 'error - uknown minimisation type, only "force", "energy" and "distance" are recognised'

      Else If (kode == 600) Then

        Write (ounit, '(/,1x,a)') 'error - "impact" option specified more than once in CONTROL'

      Else If (kode == 610) Then

        Write (ounit, '(/,1x,a)') &
          'error - "impact" applied on particle that is either frozen, or the shell of a core-shell unit or part of a RB'

      Else If (kode == 615) Then

        Write (ounit, '(/,1x,a)') 'error - q(core)*q(shell)*k(core-shell) MUST NOT be zero'

      Else If (kode == 620) Then

        Write (ounit, '(/,1x,a)') 'error - duplicate or mixed intra-molecular entries specified in FIELD'

      Else If (kode == 623) Then

        Write (ounit, '(/,1x,a)') "error - MPOLES's molecular data mismatched with respect to FIELD's data"

      Else If (kode == 625) Then

        Write (ounit, '(/,1x,a)') 'error - only one *rigid* directive per molecule is allowed'

      Else If (kode == 630) Then

        Write (ounit, '(/,1x,a)') 'error - too many rigid body units specified'

      Else If (kode == 632) Then

        Write (ounit, '(/,1x,a)') 'error - rigid body unit MUST have at least 2 sites'

      Else If (kode == 634) Then

        Write (ounit, '(/,1x,a)') 'error - rigid body unit MUST have at least one non-massless site'

      Else If (kode == 636) Then

        Write (ounit, '(/,1x,a)') 'error - rigid body unit MUST NOT have any frozen site'

      Else If (kode == 638) Then

        Write (ounit, '(/,1x,a)') 'error - coincidence of particles in a rigid body unit'

      Else If (kode == 640) Then

        Write (ounit, '(/,1x,a)') 'error - too many rigid body units per domain'

      Else If (kode == 642) Then

        Write (ounit, '(/,1x,a)') 'error - rigid body unit diameter > neigh%cutoff (the system cutoff)'

      Else If (kode == 644) Then

        Write (ounit, '(/,1x,a)') 'error - overconstrained rigid body unit'

      Else If (kode == 646) Then

        Write (ounit, '(/,1x,a)') 'error - overconstrained constraint unit'

      Else If (kode == 648) Then

        Write (ounit, '(/,1x,a)') 'error - quaternion setup failed'

      Else If (kode == 650) Then

        Write (ounit, '(/,1x,a)') 'error - failed to find principal axis system'

      Else If (kode == 655) Then

        Write (ounit, '(/,1x,a)') 'error - FENE bond breaking failure'

      Else If (kode == 660) Then

        Write (ounit, '(/,1x,a)') 'error - bond length > cutoff in TABBND or cutoff for PDF collection'

      Else If (kode == 670) Then

        Write (ounit, '(/,1x,a)') 'error - insufficient electronic temperature cells for TTM heat diffusion'

      Else If (kode == 680) Then

        Write (ounit, '(/,1x,a)') 'error - rpad too large for calculation of ionic temperatures'

      Else If (kode == 681) Then

        Write (ounit, '(/,1x,a)') 'error - electronic specific heat not fully specified'

      Else If (kode == 682) Then

        Write (ounit, '(/,1x,a)') 'error - thermal conductivity of metal not specified'

      Else If (kode == 683) Then

        Write (ounit, '(/,1x,a)') 'error - thermal diffusivity of non-metal not specified'

      Else If (kode == 684) Then

        Write (ounit, '(/,1x,a)') 'error - cannot find or open thermal conductivity table file (Ke.dat)'

      Else If (kode == 685) Then

        Write (ounit, '(/,1x,a)') 'error - no data found in thermal conductivity table file (Ke.dat)'

      Else If (kode == 686) Then

        Write (ounit, '(/,1x,a)') 'error - cannot find or open volumetric heat capacity table file (Ce.dat)'

      Else If (kode == 687) Then

        Write (ounit, '(/,1x,a)') 'error - no data found in volumetric heat capacity table file (Ce.dat)'

      Else If (kode == 688) Then

        Write (ounit, '(/,1x,a)') 'error - cannot find or open thermal diffusivity table file (De.dat)'

      Else If (kode == 689) Then

        Write (ounit, '(/,1x,a)') 'error - no data found in thermal diffusivity table file (De.dat)'

      Else If (kode == 690) Then

        Write (ounit, '(/,1x,a)') 'error - cannot find or open coupling constant table file (g.dat)'

      Else If (kode == 691) Then

        Write (ounit, '(/,1x,a)') 'error - no data found in coupling constant table file (g.dat)'

      Else If (kode == 692) Then

        Write (ounit, '(/,1x,a)') 'error - end of file encountered in table file (Ke.dat, Ce.dat, De.dat or g.dat)'

      Else If (kode == 693) Then

        Write (ounit, '(/,1x,a)') 'error - negative electronic temperature: instability in electronic heat diffusion equation'

      Else If (kode == 694) Then

        Write (ounit, '(/,1x,a)') 'error - electronic temperature restart file (DUMP_E) does not exist'

      Else If (kode == 695) Then

        Write (ounit, '(/,1x,a)') &
          'error - mismatch in electronic temperature lattice sizes between restart (DUMP_E) and CONTROL files'

      Else If (kode == 696) Then

        Write (ounit, '(/,1x,a)') 'error - cannot read electronic temperature restart (DUMP_E) file'

      Else If (kode == 1000) Then

        Write (ounit, '(/,1x,a)') 'error - working precision mismatch between FORTRAN90 and MPI implementation'

      Else If (kode == 1001) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gcheck_vector'

      Else If (kode == 1002) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gcheck_vector'

      Else If (kode == 1003) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gisum_vector'

      Else If (kode == 1004) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gisum_vector'

      Else If (kode == 1005) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grsum_vector'

      Else If (kode == 1006) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grsum_vector'

      Else If (kode == 1007) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gimax_vector'

      Else If (kode == 1008) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gimax_vector'

      Else If (kode == 1009) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grmax_vector'

      Else If (kode == 1010) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grmax_vector'

      Else If (kode == 1011) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in setup -> get_record'

      Else If (kode == 1012) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in setup -> get_record'

      Else If (kode == 1013) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angles_arrays'

      Else If (kode == 1014) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bonds_arrays'

      Else If (kode == 1015) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in core_shell -> allocate_core_shell_arrays'

      Else If (kode == 1016) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in statistics -> allocate_statitics_arrays'

      Else If (kode == 1017) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in tethers_module -> allocate_tethers_arrays'

      Else If (kode == 1018) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in constraints -> allocate_constraints_arrays'

      Else If (kode == 1019) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in external_field_module -> allocate_external_field_arrays'

      Else If (kode == 1020) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dihedrals -> allocate_dihedrals_arrays'

      Else If (kode == 1021) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in inversions -> allocate_inversion_arrays'

      Else If (kode == 1022) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in vdw -> allocate_vdw_arrays'

      Else If (kode == 1023) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_arrays'

      Else If (kode == 1024) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in three_body_module -> allocate_three_body_arrays'

      Else If (kode == 1025) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in configuration -> allocate_config_arrays'

      Else If (kode == 1026) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in site -> allocate_site_arrays'

      Else If (kode == 1027) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in tersoff_module -> allocate_tersoff_arrays'

      Else If (kode == 1030) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in core_shell -> deallocate_core_shell_arrays'

      Else If (kode == 1031) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in tethers_module -> deallocate_tethers_arrays'

      Else If (kode == 1032) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in constraints -> deallocate_constraints_arrays'

      Else If (kode == 1033) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in dihedrals -> deallocate_dihedrals_arrays'

      Else If (kode == 1034) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in inversions -> deallocate_inversions_arrays'

      Else If (kode == 1035) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in defects_module -> allocate_defects_arrays'

      Else If (kode == 1036) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in pmf -> allocate_pmf_arrays'

      Else If (kode == 1037) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in pmf -> deallocate_pmf_arrays'

      Else If (kode == 1038) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in minimise_module -> allocate_minimise_arrays'

      Else If (kode == 1040) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ewald_module -> ewald_allocate_kall_arrays'

      Else If (kode == 1041) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in langevin_module -> langevin_allocate_arrays'

      Else If (kode == 1042) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in rigid_bodies -> allocate_rigid_bodies_arrays'

      Else If (kode == 1043) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in rigid_bodies -> deallocate_rigid_bodies_arrays'

      Else If (kode == 1044) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> gimin_vector'

      Else If (kode == 1045) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> gimin_vector'

      Else If (kode == 1046) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grmin_vector'

      Else If (kode == 1047) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grmin_vector'

      Else If (kode == 1048) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in comms_module -> grsum_matrix'

      Else If (kode == 1049) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in comms_module -> grsum_matrix'

      Else If (kode == 1050) Then

        Write (ounit, '(/,1x,a)') 'error - sorted I/O base communicator not set'

      Else If (kode == 1053) Then

        Write (ounit, '(/,1x,a)') 'error - sorted I/O allocation error'

      Else If (kode == 1056) Then

        Write (ounit, '(/,1x,a)') 'error - unkown write option given to sorted I/O'

      Else If (kode == 1059) Then

        Write (ounit, '(/,1x,a)') 'error - unknown write level given to sorted I/O'

      Else If (kode == 1060) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in statistics -> allocate_statitics_connect'

      Else If (kode == 1061) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in statistics -> deallocate_statitics_connect'

      Else If (kode == 1063) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in vdw -> allocate_vdw_table_arrays'

      Else If (kode == 1066) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in vdw -> allocate_vdw_direct_fs_arrays'

      Else If (kode == 1069) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_table_arrays'

      Else If (kode == 1070) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ewald_module -> ewald_allocate_kfrz_arrays'

      Else If (kode == 1072) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bond_pot_arrays'

      Else If (kode == 1073) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bond_dst_arrays'

      Else If (kode == 1074) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angl_pot_arrays'

      Else If (kode == 1075) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angl_dst_arrays'

      Else If (kode == 1076) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dihedrals -> allocate_dihd_pot_arrays'

      Else If (kode == 1077) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dihedrals -> allocate_dihd_dst_arrays'

      Else If (kode == 1078) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in inversion_module -> allocate_invr_pot_arrays'

      Else If (kode == 1079) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in inversion_module -> allocate_invr_dst_arrays'

      Else If (kode == 1080) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in greenkubo_module -> allocate_greenkubo_arrays'

      Else If (kode == 1081) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in dpd -> allocate_dpd_arrays'

      Else If (kode == 1082) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_erf_arrays'

      Else If (kode == 1083) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_module -> allocate_ttm_arrays'

      Else If (kode == 1084) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_module -> deallocate_ttm_arrays'

      Else If (kode == 1085) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_ion_temperature'

      Else If (kode == 1086) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_ion_temperature'

      Else If (kode == 1087) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_thermal_diffusion'

      Else If (kode == 1088) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_thermal_diffusion'

      Else If (kode == 1089) Then

        Write (ounit, '(/,1x,a)') 'error - allocation failure in ttm_track_module -> depoinit'

      Else If (kode == 1090) Then

        Write (ounit, '(/,1x,a)') 'error - deallocation failure in ttm_track_module -> depoevolve'

      Else

        Write (ounit, '(/,1x,a)') 'error - unnamed error found'

      End If

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
