Subroutine warning(kode,a,b,c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for printing warning messages and returning
! control back to the main program
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : idnode
  Use setup_module

  Implicit None

  Integer,           Intent( In    ) :: kode

  Real( Kind = wp ), Intent( In    ) :: a,b,c

  Integer                            :: ia,ib,ic

  If (idnode == 0) Then

     Write(nrite,'(/,1x,a,i6)') 'warning issued ',kode

     If      (kode ==   1) Then

        ia=Nint(a)

        Write(nrite,'(/,1x,a,i0,a,/)') &
        '*** warning - node ', ia, ' mapped on vacuum (no particles) !!! ***'

     Else If (kode ==   2) Then

        ia=Nint(a)
        ib=Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - DD with ', ia, ' idle nodes, (out of ', ib, ' ) mapped on vacuum !!! ***'

     Else If (kode ==   3) Then

        Write(nrite,'(/,1x,2(a,f12.6),a,/)') &
        '*** warning - system cutoff is, ', a, ' , minimum half-cell width is, ', b, ' !!! ***'

     Else If (kode ==   4) Then

        Write(nrite,'(/,1x,2(a,/))')                            &
        '*** warning - system with uncharged particles !!! ***', &
        '*** "no elec" or/and "no strict" directives in CONTROL may speed up simulation !!! ***'

     Else If (kode ==   5) Then

        Write(nrite,'(/,1x,a,f12.4,a,/)') &
        '*** warning - non-zero total system charge: ', a, ' !!! ***'

     Else If (kode ==   6) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),2(/,a),a,/)')                         &
        '*** warning - maximum length of linked cell list: ', ia,          &
        ' is less thanmaximum length of particle exclusion list: ',        &
        ib, ' !!! ***',                                                    &
        '*** this may be due to using too short a cutoff in CONTROL ',     &
        'and/or a badly defined intramolecular topology in FIELD !!! ***', &
        '*** further potential problems may be expected !!! ***'

     Else If (kode ==   7) Then

        Write(nrite,'(/,1x,a,f12.6,a,/,1x,a,/)')                    &
        '*** warning - system cutoff is ', a, ' Angstroms !!! ***', &
        '*** Fennell damping is not recommended for cutoffs shorther than 12 Angstroms !!! ***'

     Else If (kode ==   8) Then

        Write(nrite,'(/,1x,a,2(f8.3,a),/)') &
        '*** warning - : detected maximum rigid body width: ', a, ' but system cutoff ', b, ' ***'

     Else If (kode ==  10) Then

        Write(nrite,'(/,1x,a,/)') '*** warning - no pair forces in use ***'

     Else If (kode ==  20) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a),/)') &
        '*** warning - : 1..4 scale factors reset for molecule: ', ia, ' sites: ',ib, ' & ', ic, ' ***'

     Else If (kode ==  30) Then

        Write(nrite,'(/,1x,a,/)') &
        '*** warning - Electrostatics requested in a non-periodic system !!! ***'

     Else If (kode ==  34) Then

        Write(nrite,'(/,1x,a,/)') '*** warning - redundant directive prim(ary cutoff) ignored ***'

     Else If (kode ==  35) Then

        Write(nrite,'(/,1x,a,/)') '*** warning - redundant directive delr ignored ***'

     Else If (kode ==  36) Then

        Write(nrite,'(2(/,1x,a),/)')                                                &
        "*** warning - DL_POLY_2 directive 'mult(iple timestep)' defaulted to ***", &
        "*** DL_POLY_4 'infrequent k-space SPME evaluation' directive ***"

     Else If (kode ==  37) Then

        Write(nrite,'(/,1x,a,/)') '*** warning - redundant directive all pairs ignored ***'

     Else If (kode ==  38) Then

        Write(nrite,'(/,1x,a,/)') '*** warning - redundant directive no link ignored ***'

     Else If (kode ==  40) Then

        Write(nrite,'(2(/,1x,a),f12.6,a,/)') &
        '*** warning - radial cutoff reset ***', '*** new potential cutoff radius ', a, ' ***'

     Else If (kode ==  50) Then

        Write(nrite,'(2(/,1x,a),f12.6,a,/)') &
        '*** warning - short-range cutoff reset ***', '*** new cutoff radius (rvdw) ', a, ' ***'

     Else If (kode ==  60) Then

        Write(nrite,'(/,1x,a,f12.6,a,/)') '*** warning - total system charge ', a,' ***'

     Else If (kode ==  70) Then

        Write(nrite,'(/,1x,a,f12.6,a,/)') '*** warning - switching length reset to ', a,' ***'

     Else If (kode ==  80) Then

        Write(nrite,'(/,1x,a,f12.6,a,/)') '*** warning - requested thermostat unavailable ***'

     Else If (kode ==  90) Then

        ia=Nint(a)
        ib=Nint(b)

        Write(nrite,'(2(/,1x,a),2(i0,a),/)')                  &
        '*** warning - cannot activate link cell option ***', &
        '*** more link-cells ', ia,' than allowed ', ib,' !!! ***'

     Else If (kode == 100) Then

        Write(nrite,'(2(/,1x,a),/)')                                                       &
        '*** warning - link cell algorithm has a link cell dimension that is < 4 !!! ***', &
        '*** DL_POLY_4 RUNNING IN LOW EFFICIENCY MODE !!! ***'

     Else If (kode == 110) Then

        ia=Nint(a)
        ib=Nint(b)

        Write(nrite,'(2(/,1x,a),2(i0,a),/)')                                         &
        '*** warning - image convention incompatible with the set NsT ensemble ***', &
        '*** imcon reset from ', ia,' to ', ib,' !!! ***'

     Else If (kode == 120) Then

        Write(nrite,'(/,1x,a,/)') &
        '*** warning - unspecified atom-atom interactions set to zero ***'

     Else If (kode == 130) Then

        Write(nrite,'(/,1x,a,/)') '*** warning - no ensemble is specified ***'

     Else If (kode == 140) Then

        Write(nrite,'(/,1x,2a,2(f8.5,a),/,1x,a,/)')                &
        '*** warning - control distances for variable timestep: ', &
        'Dmin = ', a, ' (Ang) and Dmax = ', b, ' (Ang) ***',       &
        '*** do not comply with safty condition: Dmax > 2.5 Dmin > 0 !!! ***'

     Else If (kode == 150) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required export buffer array size ', ia, ' and actual: ', ib, ' !!! ***'

     Else If (kode == 160) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required export array size ', ia, ' and actual: ', ib, ' !!! ***'

     Else If (kode == 170) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required export density buffer array size ', ia, ' and actual: ', ib, ' !!! ***'

     Else If (kode == 180) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required export density array size ', ia, ' and actual: ', ib, ' !!! ***'

     Else If (kode == 190) Then

        Write(nrite,'(2(/,1x,a),/)')                                             &
        '*** warning - REVOLD format mishmash detected (restart requested) ***', &
        '*** restart is abandoned and clean start is assumed ***'

     Else If (kode == 200) Then

        Write(nrite,'(2(/,1x,a),/)')                        &
        '*** warning - CONFIG contains positions only ***', &
        '*** clean start is assumed !!! ***'

     Else If (kode == 210) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(3(/,1x,a),2(i0,a),/)')                        &
        '*** warning - system under great constraint !!! ***',      &
        '*** degrees of freedom <=  total number of particles ***', &
        '*** ', ia, ' <= ', ib, ' ***'

     Else If (kode == 220) Then

        Write(nrite,'(/,1x,a,/)') &
        '*** warning - Ewald sum requested in a non-periodic system !!! ***'

     Else If (kode == 230) Then

        ia = Nint(a)

        Write(nrite,'(/,1x,a,i0,a,/,1x,a,/)')                             &
        '*** warning - PMF unit', ia, ' weight is detected zero !!! ***', &
        '*** member weights defaulted to atom type masses (or units) ***'

     Else If (kode == 240) Then

        ia = Nint(a)

        Write(nrite,'(/,1x,a,i0,a,/)') &
        '*** warning - total number of nodes running in parallel: ', ia, ' !!! ***'

     Else If (kode == 250) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required coulombic exclusion array size ', ia, ' and actual (mxexcl): ', ib, ' !!! ***'

     Else If (kode == 260) Then

        Write(nrite,'(2(/,1x,a),/)')                                                         &
        '*** warning - system volume in non-periodic systems is the MD cell volume !!! ***', &
        '*** system pressure is calculated with respect to this volume !!! ***'

     Else If (kode == 270) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required buffer size is ', ia, ' and actual ', ib, ' !!! ***'

     Else If (kode == 280) Then

        Write(nrite,'(/,1x,a,/,1x,2(a,f8.5),a,/)')                              &
        '*** warning - pseudo thermostat cannot be applied for system !!! ***', &
        '*** minimum thermostat wall thickness ', a, ' minimum MD cell width ', b, ' !!! ***'

     Else If (kode == 290) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),/)') &
        '*** warning - required link-cell list size is ', ia, ' and actual (mxlist) ', ib, ' !!! ***'

     Else If (kode == 295) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a))')                           &
        '*** warning - PMF unit', ia, ' and rigid body unit', ib, &
        ' on molecular species type ', ic, ' in illegal configuration !!! ***'

     Else If (kode == 301) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a),/)')                                &
        '*** warning - core-shell unit', ia, ' and constraint unit', ib, &
        ' on molecular species type ', ic, ' in illegal configuration !!! ***'

     Else If (kode == 302) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a))')                                  &
        '*** warning - core-shell unit', ia, ' and rigid body unit', ib, &
        ' on molecular species type ', ic, ' in illegal configuration !!! ***'

     Else If (kode == 303) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a),/)')                            &
        '*** warning - core-shell unit', ia, ' and tether unit', ib, &
        ' on molecular species type ', ic, ' in illegal configuration !!! ***'

     Else If (kode == 304) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a))')                                  &
        '*** warning - constraint unit', ia, ' and rigid body unit', ib, &
        ' on molecular species type ', ic, ' in illegal configuration !!! ***'

     Else If (kode == 305) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a))')                                         &
        '*** warning - rigid body unit', ia, ' on molecular species type ', ib, &
        ' forced to freeze !!! ***'

     Else If (kode == 306) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a))')                                         &
        '*** warning - constraint unit', ia, ' on molecular species type ', ib, &
        ' forced to freeze !!! ***'

     Else If (kode == 307) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a))')                                         &
        '*** warning - rigid body unit', ia, ' on molecular species type ', ib, &
        ' set to have type ', ic, ' is problematic !!! ***'

     Else If (kode == 308) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a))')                       &
        '*** warning - site', ia, ' of constraint unit ', ib, &
        ' on molecular species type', ic, ' is problematic !!! ***'

     Else If (kode == 309) Then

        ia = Nint(a)
        ib = Nint(b)
        ic = Nint(c)

        Write(nrite,'(/,1x,a,3(i0,a))') &
        '*** warning - site', ia, ' , member ', ib, ' of PMF unit', ic, ' is problematic !!! ***'

     Else If (kode == 310) Then

        Write(nrite,'(/,1x,a,/,1x,a,2(f6.2,a),/)')                                                                                 &
        '*** warning - control distance for defect look-up MUST be in the interval [Min(0.3,rcut/3);Min(1.2,rcut/2)] Ang !!! ***', &
        '*** defects distance condition will default from ', a, ' to ', b, ' Ang !!! ***'

     Else If (kode == 320) Then

        Write(nrite,'(/,1x,a,/)') &
        '*** warning - REFERENCE not found, CONFIG to be used as a reference for defects detection !!! ***'

     Else If (kode == 330) Then

        Write(nrite,'(/,1x,a,/,1x,a,2(f10.6,a),/)')                                                      &
        '*** warning - iteration cycles length limit for conjugate gradient minimisers exceded !!! ***', &
        '*** specified convergence tolerance: ', a, ' , needed one for a pass: ', b, ' !!! ***'

     Else If (kode == 340) Then

        Write(nrite,'(3(/,1x,a),2(f7.4,a),/,1x,a,f7.4,a,/)')                             &
        '*** warning - inconsistent binsize for spatial distribution functions !!! ***', &
        '*** 1.0E-5 (Angstroms) <= binsize <= rcut/4 (Angstroms) ***',                   &
        '*** 1.0E-5 (Angstroms) <= ', a, ' <= ', b, ' (Angstroms) ***',                  &
        '*** binsize defaults to ', c, ' (Angstroms) ***'

     Else If (kode == 350) Then

        Write(nrite,'(/,1x,a)') &
        '*** warning - expansion along z axis not allowed for slab geometry, nz defaults to 1!!! ***'

     Else If (kode == 360) Then

        Write(nrite,'(/,1x,a,2(f12.7,a),3(/,1x,a),/)')                            &
        '*** warning - minimisation tolerance', a, ' defaults to', b, ' !!! ***', &
        '*** force   : 1.0    <= tolerance <= 1000.00, default = 50.00 ***',      &
        '*** energy  : 0.0    <  tolerance <=    0.01, default = 0.005 ***',      &
        '*** distance: 1.0e-6 <= tolerance <=    0.10, default = 0.005 ***'

     Else If (kode == 370) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2(i0,a),2(/,1x,a),/)')                                      &
        '*** warning - k-space evaluation interval', ia, ' defaults to', ib, ' !!! ***', &
        '*** the interval must be a positive integer beteen 1 and 10 ***',               &
        '*** values > 10 default to 4, no value or 0 defaults to 1 ***'

     Else If (kode == 380) Then

        ia = Nint(a)

        Write(nrite,'(/,1x,a,i0,a,/)') &
        '*** warning - IMPACT applied before the end of the equlibration period (', ia, ' step) !!! ***'

     Else If (kode == 390) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2i0,a,/)') &
        '*** warning - core-shell unit mix-up or duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 400) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,i0,a,/)') &
        '*** warning - constraint or rigid body duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 410) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,i0,a,/)') &
        '*** warning - tethered atom duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 420) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2i0,a,/)') &
        '*** warning - chemical bond duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 430) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2i0,a,/)') &
        '*** warning - bond angle duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 440) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2i0,a,/)') &
        '*** warning - dihedral angle duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 450) Then

        ia = Nint(a)
        ib = Nint(b)

        Write(nrite,'(/,1x,a,2i0,a,/)') &
        '*** warning - inversion angle duplicate specification (', ia, ib, ' ) !!! ***'

     Else If (kode == 460) Then

        Write(nrite,'(/,1x,a,/)') &
        '*** warning - further semi-isotropic barostat option search abandoned !!! ***'

     Else If (kode == 470) Then

        Write(nrite,'(/,1x,a,/,1x,a,2(f6.2,a),/)')                                                  &
        '*** warning - control distance for diplacement qualification MUST be >= 0.25 Ang !!! ***', &
        '*** displacements distance condition will default from ', a, ' to ', b, ' Ang !!! ***'

     Else If (kode == 480) Then

        Write(nrite,'(/,1x,a,/,1x,a,2(f6.2,a),/)') &
        '*** warning - "metal direct" option disabled as incompatible with EAM (TABEAM)!!! ***'

     Else

        Write(nrite,'(/,1x,a,/)') &
        '*** unspecified warning encountered ***'

     End If

  End If

End Subroutine warning
