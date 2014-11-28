Subroutine error(kode)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for printing error messages and bringing about a
! controlled termination of the program
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module, Only : idnode,abort_comms
  Use setup_module, Only : nread,nconf,nfield,ntable,nrefdt,nrite, &
                           nstats,nrest,nhist,ndefdt,nrdfdt,nzdndt,nrsddt

  Implicit None

  Integer, Intent( In    ) :: kode

  If (idnode == 0) Then

     Write(nrite,'(/,1x,a,i5)') 'DL_POLY_4 terminated due to error ', kode

     If      (kode ==    1) Then

         Write(nrite,'(/,1x,a)') 'error - word_2_real failure'

     Else If (kode ==    2) Then

         Write(nrite,'(/,1x,a)') 'error - too many atom types in FIELD (scan_field)'

     Else If (kode ==    3) Then

        Write(nrite,'(/,1x,a)') 'error - unknown directive found in CONTROL file'

     Else If (kode ==    4) Then

        Write(nrite,'(/,1x,a)') 'error - unknown directive found in FIELD file'

     Else If (kode ==    5) Then

        Write(nrite,'(/,1x,a)') 'error - unknown energy unit requested'

     Else If (kode ==    6) Then

        Write(nrite,'(/,1x,a)') 'error - energy unit not specified'

     Else If (kode ==    7) Then

        Write(nrite,'(/,1x,a)') 'error - selected external field incompatible with selected ensemble (NVE only!!!)'

     Else If (kode ==    8) Then

        Write(nrite,'(/,1x,a)') 'error - ewald precision must be a POSITIVE real number'

     Else If (kode ==   10) Then

        Write(nrite,'(/,1x,a)') 'error - too many molecular types specified'

     Else If (kode ==   11) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate molecule directive in FIELD file'

     Else If (kode ==   12) Then

        Write(nrite,'(/,1x,a)') 'error - unknown molecule directive in FIELD file'

     Else If (kode ==   13) Then

        Write(nrite,'(/,1x,a)') 'error - molecular species not yet specified'

     Else If (kode ==   14) Then

        Write(nrite,'(/,1x,a)') 'error - too many unique atom types specified'

     Else If (kode ==   15) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate vdw potential specified'

     Else If (kode ==   16) Then

        Write(nrite,'(/,1x,a)') 'error - strange exit from FIELD file processing'

     Else If (kode ==   17) Then

        Write(nrite,'(/,1x,a)') 'error - strange exit from CONTROL file processing'

     Else If (kode ==   18) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate three-body potential specified'

     Else If (kode ==   19) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate four-body potential specified'

     Else If (kode ==   20) Then

        Write(nrite,'(/,1x,a)') 'error - too many molecule sites specified'

     Else If (kode ==   21) Then

        Write(nrite,'(/,1x,a)') 'error - molecule contains more atoms/sites than declared'

     Else If (kode ==   22) Then

        Write(nrite,'(/,1x,a)') 'error - unsuitable radial increment in TABLE||TABBND||TABANG||TABDIH||TABINV file'

     Else If (kode ==   23) Then

        Write(nrite,'(/,1x,a)') 'error - incompatible FIELD and TABLE file potentials'

     Else If (kode ==   24) Then

        Write(nrite,'(/,1x,a)') 'error - end of file encountered in TABLE||TABEAM||TABBND||TABANG||TABDIH||TABINV file'

     Else If (kode ==   25) Then

        Write(nrite,'(/,1x,a)') 'error - wrong atom type found in CONFIG file'

     Else If (kode ==   26) Then

        Write(nrite,'(/,1x,a)') 'error - neutral group option now redundant'

     Else If (kode ==   27) Then

        Write(nrite,'(/,1x,a)') "error - unit's member indexed outside molecule's site range"

     Else If (kode ==   28) Then

        Write(nrite,'(/,1x,a)') 'error - wrongly indexed atom entries found in CONFIG file'

     Else If (kode ==   30) Then

        Write(nrite,'(/,1x,a)') 'error - too many chemical bonds specified'

     Else If (kode ==   31) Then

        Write(nrite,'(/,1x,a)') 'error - too many chemical bonds per domain'

     Else If (kode ==   32) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in core-shell unit'

     Else If (kode ==   33) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in constraint bond unit'

     Else If (kode ==   34) Then

        Write(nrite,'(/,1x,a)') 'error - length of constraint bond unit >= real space cutoff (rcut)'

     Else If (kode ==   35) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in chemical bond unit'

     Else If (kode ==   36) Then

        Write(nrite,'(/,1x,a)') 'error - only one *bonds* directive per molecule is allowed'

     Else If (kode ==   38) Then

        Write(nrite,'(/,1x,a)') 'error - outgoing transfer buffer size exceeded in metal_ld_export'

     Else If (kode ==   39) Then

        Write(nrite,'(/,1x,a)') 'error - incoming data transfer size exceeds limit in metal_ld_export'

     Else If (kode ==   40) Then

        Write(nrite,'(/,1x,a)') 'error - too many bond constraints specified'

     Else If (kode ==   41) Then

        Write(nrite,'(/,1x,a)') 'error - too many bond constraints per domain'

     Else If (kode ==   42) Then

        Write(nrite,'(/,1x,a)') 'error - undefined direction passed to deport_atomic_data'

     Else If (kode ==   43) Then

        Write(nrite,'(/,1x,a)') 'error - outgoing transfer buffer exceeded in deport_atomic_data'

     Else If (kode ==   44) Then

        Write(nrite,'(/,1x,a)') 'error - incoming data transfer size exceeds limit in deport_atomic_data'

     Else If (kode ==   45) Then

        Write(nrite,'(/,1x,a)') 'error - too many atoms in CONFIG file or per domain'

     Else If (kode ==   46) Then

        Write(nrite,'(/,1x,a)') 'error - undefined direction passed to export_atomic_data'

     Else If (kode ==   47) Then

        Write(nrite,'(/,1x,a)') 'error - undefined direction passed to metal_ld_export'

     Else If (kode ==   48) Then

        Write(nrite,'(/,1x,a)') 'error - transfer buffer too small in *_table_read'

     Else If (kode ==   49) Then

        Write(nrite,'(/,1x,a)') 'error - frozen shell (core-shell unit) specified'

     Else If (kode ==   50) Then

        Write(nrite,'(/,1x,a)') 'error - too many bond angles specified'

     Else If (kode ==   51) Then

        Write(nrite,'(/,1x,a)') 'error - too many bond angles per domain'

     Else If (kode ==   52) Then

        Write(nrite,'(/,1x,a)') 'error - end of FIELD file encountered'

     Else If (kode ==   53) Then

        Write(nrite,'(/,1x,a)') 'error - end of CONTROL file encountered'

     Else If (kode ==   54) Then

        Write(nrite,'(/,1x,a)') 'error - outgoing transfer buffer size exceeded in export_atomic_data'

     Else If (kode ==   55) Then

        Write(nrite,'(/,1x,a)') 'error - end of CONFIG file encountered'

     Else If (kode ==   56) Then

        Write(nrite,'(/,1x,a)') 'error - incoming data transfer size exceeds limit in export_atomic_data'

     Else If (kode ==   57) Then

        Write(nrite,'(/,1x,a)') 'error - too many core-shell units specified'

     Else If (kode ==   58) Then

        Write(nrite,'(/,1x,a)') 'error - number of atoms in system not conserved'

     Else If (kode ==   59) Then

        Write(nrite,'(/,1x,a)') 'error - too many core-shell units per domain'

     Else If (kode ==   60) Then

        Write(nrite,'(/,1x,a)') 'error - too many dihedral angles specified'

     Else If (kode ==   61) Then

        Write(nrite,'(/,1x,a)') 'error - too many dihedral angles per domain'

     Else If (kode ==   62) Then

        Write(nrite,'(/,1x,a)') 'error - too many tethered atoms specified'

     Else If (kode ==   63) Then

        Write(nrite,'(/,1x,a)') 'error - too many tethered atoms per domain'

     Else If (kode ==   64) Then

        Write(nrite,'(/,1x,a)') 'error - incomplete core-shell unit found in build_book_intra'

     Else If (kode ==   65) Then

        Write(nrite,'(/,1x,a)') 'error - too many excluded pairs specified'

     Else If (kode ==   66) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in bond angle unit'

     Else If (kode ==   67) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in dihedral unit'

     Else If (kode ==   68) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in inversion unit'

     Else If (kode ==   70) Then

        Write(nrite,'(/,1x,a)') 'error - constraint_quench failure'

     Else If (kode ==   71) Then

        Write(nrite,'(/,1x,a)') 'error - too many metal potentials specified'

     Else If (kode ==   72) Then

        Write(nrite,'(/,1x,a)') 'error - too many tersoff potentials specified'

     Else If (kode ==   73) Then

        Write(nrite,'(/,1x,a)') 'error - too many inversion potentials specified'

     Else If (kode ==   74) Then

        Write(nrite,'(/,1x,a)') 'error - unidentified atom in tersoff potential list'

     Else If (kode ==   76) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate tersoff potential specified'

     Else If (kode ==   77) Then

        Write(nrite,'(/,1x,a)') 'error - too many inversion angles per domain'

     Else If (kode ==   79) Then

        Write(nrite,'(/,1x,a)') 'error - tersoff potential cutoff undefined'

     Else If (kode ==   80) Then

        Write(nrite,'(/,1x,a)') 'error - too many pair potentials specified'

     Else If (kode ==   81) Then

        Write(nrite,'(/,1x,a)') 'error - unidentified atom in pair potential list'

     Else If (kode ==   82) Then

        Write(nrite,'(/,1x,a)') 'error - calculated pair potential index too large'

     Else If (kode ==   83) Then

        Write(nrite,'(/,1x,a)') 'error - too many three-body/angles potentials specified'

     Else If (kode ==   84) Then

        Write(nrite,'(/,1x,a)') 'error - unidentified atom in three-body/angles potential list'

     Else If (kode ==   85) Then

        Write(nrite,'(/,1x,a)') 'error - required velocities not in CONFIG file'

     Else If (kode ==   86) Then

        Write(nrite,'(/,1x,a)') 'error - calculated three-body potential index too large'

     Else If (kode ==   88) Then

        Write(nrite,'(/,1x,a)') 'error - legend array exceeded in build_book_intra'

     Else If (kode ==   89) Then

        Write(nrite,'(/,1x,a)') 'error - too many four-body/dihedrals/inversions potentials specified'

     Else If (kode ==   90) Then

        Write(nrite,'(/,1x,a)') 'error - specified tersoff potentials have different types'

     Else If (kode ==   91) Then

        Write(nrite,'(/,1x,a)') 'error - unidentified atom in four-body/dihedrals/iversions potential list'

     Else If (kode ==   92) Then

        Write(nrite,'(/,1x,a)') 'error - specified metal potentials have different types'

     Else If (kode ==   93) Then

        Write(nrite,'(/,1x,a)') 'error - PMFs mixing with rigid bodies not allowed'

     Else If (kode ==   95) Then

        Write(nrite,'(/,1x,a)') 'error - rcut (or rcut+rpad) > minimum of all half-cell widths'

     Else If (kode ==   96) Then

        Write(nrite,'(/,1x,a)') 'error - incorrect atom totals assignments in metal_ld_set_halo'

     Else If (kode ==   97) Then

        Write(nrite,'(/,1x,a)') 'error - constraints mixing with rigid bodies not allowed'

     Else If (kode ==   99) Then

        Write(nrite,'(/,1x,a)') 'error - cannot have shells as part of a constraint, rigid body or tether'

     Else If (kode ==  100) Then

        Write(nrite,'(/,1x,a)') 'error - core-shell unit separation > rcut (the system cutoff)'

     Else If (kode ==  101) Then

        Write(nrite,'(/,1x,a)') 'error - calculated four-body potential index too large'

     Else If (kode ==  102) Then

        Write(nrite,'(/,1x,a)') 'error - rcut < 2*rcter (maximum cutoff for tersoff potentials)'

     Else If (kode ==  103) Then

        Write(nrite,'(/,1x,a)') 'error - parameter mxlshp exceeded in pass_shared_units'

     Else If (kode ==  104) Then

        Write(nrite,'(/,1x,a)') 'error - arrays listme and lstout exceeded in pass_shared_units'

     Else If (kode ==  105) Then

        Write(nrite,'(/,1x,a)') 'error - shake algorithm (constraints_shake) failed to converge'

     Else If (kode ==  106) Then

        Write(nrite,'(/,1x,a)') 'error - neighbour list array too small in link_cell_pairs'

     Else If (kode ==  107) Then

        Write(nrite,'(/,1x,a)') 'error - too many pairs for rdf look up specified'

     Else If (kode ==  108) Then

        Write(nrite,'(/,1x,a)') 'error - unidentified atom in rdf look up list'

     Else If (kode ==  109) Then

        Write(nrite,'(/,1x,a)') 'error - calculated pair rdf index too large'

     Else If (kode ==  110) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate rdf look up pair specified'

     Else If (kode ==  111) Then

        Write(nrite,'(/,1x,a)') 'error - bond constraint separation > rcut (the system cutoff)'

     Else If (kode ==  112) Then

        Write(nrite,'(/,1x,a)') 'error - only one *constraints* directive per molecule is allowed'

     Else If (kode ==  113) Then

        Write(nrite,'(/,1x,a)') 'error - intra-molecular bookkeeping arrays exceeded in deport_atomic_data'

     Else If (kode ==  114) Then

        Write(nrite,'(/,1x,a)') 'error - legend array exceeded in deport_atomic_data'

     Else If (kode ==  115) Then

        Write(nrite,'(/,1x,a)') 'error - transfer buffer exceeded in update_shared_units'

     Else If (kode ==  116) Then

        Write(nrite,'(/,1x,a)') 'error - incorrect atom transfer in update_shared_units'

     Else If (kode ==  118) Then

        Write(nrite,'(/,1x,a)') 'error - construction error in pass_shared_units'

     Else If (kode ==  120) Then

        Write(nrite,'(/,1x,a)') 'error - invalid determinant in matrix inversion'

     Else If (kode ==  122) Then

        Write(nrite,'(/,1x,a)') 'error - FIELD file not found'

     Else If (kode ==  124) Then

        Write(nrite,'(/,1x,a)') 'error - CONFIG file not found'

     Else If (kode ==  126) Then

        Write(nrite,'(/,1x,a)') 'error - CONTROL file not found'

     Else If (kode ==  128) Then

        Write(nrite,'(/,1x,a)') 'error - chemical bond unit separation > rcut (the system cutoff)'

     Else If (kode ==  130) Then

        Write(nrite,'(/,1x,a)') 'error - bond angle unit diameter > rcut (the system cutoff)'

     Else If (kode ==  132) Then

        Write(nrite,'(/,1x,a)') 'error - dihedral angle unit diameter > rcut (the system cutoff)'

     Else If (kode ==  134) Then

        Write(nrite,'(/,1x,a)') 'error - inversion angle unit diameter > rcut (the system cutoff)'

     Else If (kode ==  138) Then

        Write(nrite,'(/,1x,a)') 'error - incorrect atom totals assignments in refresh_halo_positions'

     Else If (kode ==  141) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate metal potential specified'

     Else If (kode ==  145) Then

        Write(nrite,'(/,1x,a)') 'error - no two-body like forces specified'

     Else If (kode ==  150) Then

        Write(nrite,'(/,1x,a)') 'error - unknown van der waals potential selected'

     Else If (kode ==  151) Then

        Write(nrite,'(/,1x,a)') 'error - unknown EAM keyword in TABEAM'

     Else If (kode ==  152) Then

        Write(nrite,'(/,1x,a)') 'error - undefined direction passed to dpd_v_export'

     Else If (kode ==  154) Then

        Write(nrite,'(/,1x,a)') 'error - outgoing transfer buffer exceeded in dpd_v_export'

     Else If (kode ==  156) Then

        Write(nrite,'(/,1x,a)') 'error - incoming data transfer exceeds limit in dpd_v_export'

     Else If (kode ==  158) Then

        Write(nrite,'(/,1x,a)') 'error - incorrect atom totals assignments in dpd_v_set_halo'

     Else If (kode ==  160) Then

        Write(nrite,'(/,1x,a)') 'error - undefined direction passed to statistics_connect_spread'

     Else If (kode ==  163) Then

        Write(nrite,'(/,1x,a)') 'error - outgoing transfer buffer size exceeded in statistics_connect_spread'

     Else If (kode ==  164) Then

        Write(nrite,'(/,1x,a)') 'error - incoming data transfer size exceeds limit in statistics_connect_spread'

     Else If (kode ==  170) Then

        Write(nrite,'(/,1x,a)') 'error - too many variables for statistics array'

     Else If (kode ==  172) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate intra-molecular entries specified in TABBND||TABANG||TABDIH||TABINV'

     Else If (kode ==  200) Then

        Write(nrite,'(/,1x,a)') 'error - rdf||z-density||dstbnd||dstang||dstdih||dstinv buffer array too small in system_revive'

     Else If (kode ==  210) Then

        Write(nrite,'(/,1x,a)') 'error - only one *angles* directive per molecule is allowed'

     Else If (kode ==  220) Then

        Write(nrite,'(/,1x,a)') 'error - only one *dihedrals* directive per molecule is allowed'

     Else If (kode ==  230) Then

        Write(nrite,'(/,1x,a)') 'error - only one *inversions* directive per molecule is allowed'

     Else If (kode ==  240) Then

        Write(nrite,'(/,1x,a)') 'error - only one *tethers* directive per molecule is allowed'

     Else If (kode ==  300) Then

        Write(nrite,'(/,1x,a)') 'error - incorrect boundary condition for link-cell algorithms'

     Else If (kode ==  305) Then

        Write(nrite,'(/,1x,a)') 'error - too few link cells per dimension for many-body and tersoff forces subroutines'

     Else If (kode ==  307) Then

        Write(nrite,'(/,1x,a)') 'error - link cell algorithm violation'

     Else If (kode ==  308) Then

        Write(nrite,'(/,1x,a)') 'error - link cell algorithm in contention with SPME sum precision'

     Else If (kode ==  321) Then

        Write(nrite,'(/,1x,a)') 'error - LFV quaternion integrator failed'

     Else If (kode ==  340) Then

        Write(nrite,'(/,1x,a)') 'error - invalid integration option requested'

     Else If (kode ==  350) Then

        Write(nrite,'(/,1x,a)') 'error - too few degrees of freedom'

     Else If (kode ==  360) Then

        Write(nrite,'(/,1x,a)') 'error - degrees of freedom distribution problem'

     Else If (kode ==  380) Then

        Write(nrite,'(/,1x,a)') 'error - simulation temperature not specified or < 1 K'

     Else If (kode ==  381) Then

        Write(nrite,'(/,1x,a)') 'error - simulation timestep not specified'

     Else If (kode ==  382) Then

        Write(nrite,'(/,1x,a)') 'error - simulation cutoff not specified'

     Else If (kode ==  387) Then

        Write(nrite,'(/,1x,a)') 'error - system pressure not specified'

     Else If (kode ==  390) Then

        Write(nrite,'(/,1x,a)') 'error - npt/nst ensemble requested in non-periodic system'

     Else If (kode ==  402) Then

        Write(nrite,'(/,1x,a)') 'error - van der waals cutoff not specified'

     Else If (kode ==  410) Then

        Write(nrite,'(/,1x,a)') 'error - cell not consistent with image convention'

     Else If (kode ==  414) Then

        Write(nrite,'(/,1x,a)') 'error - conflicting ensemble options in CONTROL file'

     Else If (kode ==  416) Then

        Write(nrite,'(/,1x,a)') 'error - conflicting force options in CONTROL file'

     Else If (kode ==  430) Then

        Write(nrite,'(/,1x,a)') 'error - integration routine not available'

     Else If (kode ==  432) Then

        Write(nrite,'(/,1x,a)') 'error -  undefined tersoff potential'

     Else If (kode ==  433) Then

        Write(nrite,'(/,1x,a)') 'error - rcut must be specified for the Ewald sum precision'

     Else If (kode ==  436) Then

        Write(nrite,'(/,1x,a)') 'error - unrecognised ensemble'

     Else If (kode ==  440) Then

        Write(nrite,'(/,1x,a)') 'error - undefined angular potential'

     Else If (kode ==  442) Then

        Write(nrite,'(/,1x,a)') 'error - undefined three-body potential'

     Else If (kode ==  443) Then

        Write(nrite,'(/,1x,a)') 'error - undefined four-body potential'

     Else If (kode ==  444) Then

        Write(nrite,'(/,1x,a)') 'error - undefined bond potential'

     Else If (kode ==  445) Then

        Write(nrite,'(/,1x,a)') 'error - r_14 > rcut in dihedrals_forces'

     Else If (kode ==  446) Then

        Write(nrite,'(/,1x,a)') 'error - undefined electrostatic key in dihedrals_forces'

     Else If (kode ==  448) Then

        Write(nrite,'(/,1x,a)') 'error - undefined dihedral potential'

     Else If (kode ==  449) Then

        Write(nrite,'(/,1x,a)') 'error - undefined inversion potential'

     Else If (kode ==  450) Then

        Write(nrite,'(/,1x,a)') 'error - undefined tethering potential'

     Else If (kode ==  451) Then

        Write(nrite,'(/,1x,a)') 'error - three-body potential cutoff undefined'

     Else If (kode ==  452) Then

        Write(nrite,'(/,1x,a)') 'error - undefined pair potential'

     Else If (kode ==  453) Then

        Write(nrite,'(/,1x,a)') 'error - four-body potential cutoff undefined'

     Else If (kode ==  454) Then

        Write(nrite,'(/,1x,a)') 'error - undefined external field'

     Else If (kode ==  456) Then

        Write(nrite,'(/,1x,a)') 'error - external field xpis-ton is applied to a layer with at least one frozen particle'

     Else If (kode ==  461) Then

        Write(nrite,'(/,1x,a)') 'error - undefined metal potential'

     Else If (kode ==  462) Then

        Write(nrite,'(/,1x,a)') 'error - thermostat friction constant must be > 0'

     Else If (kode ==  463) Then

        Write(nrite,'(/,1x,a)') 'error - barostat friction constant must be > 0'

     Else If (kode ==  464) Then

        Write(nrite,'(/,1x,a)') 'error - thermostat relaxation time constant must be > 0'

     Else If (kode ==  466) Then

        Write(nrite,'(/,1x,a)') 'error - barostat relaxation time constant must be > 0'

     Else If (kode ==  467) Then

        Write(nrite,'(/,1x,a)') 'error - rho must not be zero in valid buckingham potential'

     Else If (kode ==  468) Then

        Write(nrite,'(/,1x,a)') 'error - r0 too large for snm potential with current cutoff'

     Else If (kode ==  470) Then

        Write(nrite,'(/,1x,a)') 'error - n < m in definition of n-m potential'

     Else If (kode ==  471) Then

        Write(nrite,'(/,1x,a)') 'error - rcut < 2*rctbp (maximum cutoff for three-body potentials)'

     Else If (kode ==  472) Then

        Write(nrite,'(/,1x,a)') 'error - rcut < 2*rcfbp (maximum cutoff for four-body potentials)'

     Else If (kode ==  474) Then

        Write(nrite,'(/,1x,a)') 'error - conjugate gradient mimimiser cycle limit exceeded'

     Else If (kode ==  476) Then

        Write(nrite,'(/,1x,a)') 'error - shells MUST all HAVE either zero or non-zero masses'

     Else If (kode ==  477) Then

        Write(nrite,'(/,1x,a)') 'error - only one *shells* directive per molecule is allowed'

     Else If (kode ==  478) Then

        Write(nrite,'(/,1x,a)') 'error - shake algorithms (constraints & pmf) failed to converge'

     Else If (kode ==  480) Then

        Write(nrite,'(/,1x,a)') 'error - PMF length > minimum of all half-cell widths'

     Else If (kode ==  484) Then

        Write(nrite,'(/,1x,a)') 'error - only one potential of mean force permitted'

     Else If (kode ==  486) Then

        Write(nrite,'(/,1x,a)') 'error - only one of the PMF units is permitted to have frozen atoms'

     Else If (kode ==  488) Then

        Write(nrite,'(/,1x,a)') 'error - too many PMF constraints per domain'

     Else If (kode ==  490) Then

        Write(nrite,'(/,1x,a)') 'error - local PMF constraint not found locally'

     Else If (kode ==  492) Then

        Write(nrite,'(/,1x,a)') 'error - a diameter of a PMF unit > minimum of all half cell widths'

     Else If (kode ==  494) Then

        Write(nrite,'(/,1x,a)') 'error - overconstrained PMF units'

     Else If (kode ==  497) Then

        Write(nrite,'(/,1x,a)') 'error - pmf_quench failure'

     Else If (kode ==  498) Then

        Write(nrite,'(/,1x,a)') 'error - shake algorithm (pmf_shake) failed to converge'

     Else If (kode ==  499) Then

        Write(nrite,'(/,1x,a)') 'error - rattle algorithm (pmf_rattle) failed to converge'

     Else If (kode ==  500) Then

        Write(nrite,'(/,1x,a)') 'error - PMF unit of zero length is not permitted'

     Else If (kode ==  501) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in PMF unit'

     Else If (kode ==  502) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of units in PMF constraint'

     Else If (kode ==  504) Then

        Write(nrite,'(/,1x,a)') 'error - cutoff too large for TABLE file'

     Else If (kode ==  505) Then

        Write(nrite,'(/,1x,a)') 'error - EAM metal densities or pair crossfunctions out of range'

     Else If (kode ==  506) Then

        Write(nrite,'(/,1x,a)') 'error - EAM metal densities out of range'

     Else If (kode ==  507) Then

        Write(nrite,'(/,1x,a)') 'error - metal density embedding out of range'

     Else If (kode ==  508) Then

        Write(nrite,'(/,1x,a)') 'error - EAM metal interaction entry in TABEAM unspecified in FIELD'

     Else If (kode ==  509) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate entry for a pair interaction detected in TABEAM'

     Else If (kode ==  510) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate entry for a density function detected in TABEAM'

     Else If (kode ==  511) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate entry for an embedding function detected in TABEAM'

     Else If (kode ==  512) Then

        Write(nrite,'(/,1x,a)') 'error - non-definable vdw/dpd interactions detected in FIELD'

     Else If (kode ==  513) Then

        Write(nrite,'(/,1x,a)') 'error - particle assigned to non-existent domain in read_config'

     Else If (kode ==  514) Then

        Write(nrite,'(/,1x,a)') 'error - allowed image conventions are: 0, 1, 2, 3 and 6'

     Else If (kode ==  515) Then

        Write(nrite,'(/,1x,a)') 'error - rattle algorithm (constraints_rattle) failed to converge'

     Else If (kode ==  516) Then

        Write(nrite,'(/,1x,a)') 'error - the number nodes MUST be a power of 2 series number'

     Else If (kode ==  517) Then

        Write(nrite,'(/,1x,a)') 'error - allowed configuration information levels are: 0, 1 and 2'

     Else If (kode ==  518) Then

        Write(nrite,'(/,1x,a)') 'error - control distances for variable timestep not intact'

     Else If (kode ==  519) Then

        Write(nrite,'(/,1x,a)') 'error - REVOLD is incompatible or does not exist'

     Else If (kode ==  520) Then

        Write(nrite,'(/,1x,a)') 'error - domain decomposition failed'

     Else If (kode ==  530) Then

        Write(nrite,'(/,2(1x,a,/))') 'error - pseudo thermostat thickness MUST comply with', &
             '2 Angs <= thickness < a quarter of the minimum MD cell width'

     Else If (kode ==  540) Then

        Write(nrite,'(/,2(1x,a,/))') 'error - pseudo thermostat can ONLY be used in bulk simulations', &
             'i.e. imcon MUST be 1, 2 or 3'

     Else If (kode ==  551) Then

        Write(nrite,'(/,1x,a)') 'error - REFERENCE not found !!!'

     Else If (kode ==  552) Then

        Write(nrite,'(/,1x,a)') 'error - REFERENCE must contain cell parameters !!!'

     Else If (kode ==  553) Then

        Write(nrite,'(/,1x,a)') 'error - REFERENCE is inconsistent !!!'

     Else If (kode ==  554) Then

        Write(nrite,'(/,1x,a)') "error - REFERENCE's format different from CONFIG's !!!"

     Else If (kode ==  555) Then

        Write(nrite,'(/,1x,a)') 'error - particle assigned to non-existent domain in defects_read_reference'

     Else If (kode ==  556) Then

        Write(nrite,'(/,1x,a)') 'error - too many atoms in REFERENCE file'

     Else If (kode ==  557) Then

        Write(nrite,'(/,1x,a)') 'error - undefined direction passed to defects_reference_export'

     Else If (kode ==  558) Then

        Write(nrite,'(/,1x,a)') 'error - outgoing transfer buffer exceeded in defects_reference_export'

     Else If (kode ==  559) Then

        Write(nrite,'(/,1x,a)') 'error - incoming data transfer size exceeds limit in defects_reference_export'

     Else If (kode ==  560) Then

        Write(nrite,'(/,1x,a)') 'error - rdef found to be > half the shortest interatomic distance in REFERENCE'

     Else If (kode ==  570) Then

        Write(nrite,'(/,1x,a)') 'error - unsupported image convention (0) for system expansion option nfold'

     Else If (kode ==  580) Then

        Write(nrite,'(/,1x,a)') 'error - replay (HISTORY) option can only be used for structural property recalculation'

     Else If (kode ==  585) Then

        Write(nrite,'(/,1x,a)') 'error - HISTORY file does not exist'

     Else If (kode ==  590) Then

        Write(nrite,'(/,1x,a)') 'error - uknown minimisation type, only "force", "energy" and "distance" are recognised'

     Else If (kode ==  600) Then

        Write(nrite,'(/,1x,a)') 'error - "impact" option specified more than once in CONTROL'

     Else If (kode ==  610) Then

        Write(nrite,'(/,1x,2a)') &
        'error - "impact" applied on particle that is either frozen, or the shell of a core-shell unit or part of a RB'

     Else If (kode ==  620) Then

        Write(nrite,'(/,1x,a)') 'error - duplicate or mixed intra-molecular entries specified in FIELD'

     Else If (kode ==  625) Then

        Write(nrite,'(/,1x,a)') 'error - only one *rigid* directive per molecule is allowed'

     Else If (kode ==  630) Then

        Write(nrite,'(/,1x,a)') 'error - too many rigid body units specified'

     Else If (kode ==  632) Then

        Write(nrite,'(/,1x,a)') 'error - rigid body unit MUST have at least 2 sites'

     Else If (kode ==  634) Then

        Write(nrite,'(/,1x,a)') 'error - rigid body unit MUST have at least one non-massless site'

     Else If (kode ==  636) Then

        Write(nrite,'(/,1x,a)') 'error - rigid body unit MUST NOT have any frozen site'

     Else If (kode ==  638) Then

        Write(nrite,'(/,1x,a)') 'error - coincidence of particles in a rigid body unit'

     Else If (kode ==  640) Then

        Write(nrite,'(/,1x,a)') 'error - too many rigid body units per domain'

     Else If (kode ==  642) Then

        Write(nrite,'(/,1x,a)') 'error - rigid body unit diameter > rcut (the system cutoff)'

     Else If (kode ==  644) Then

        Write(nrite,'(/,1x,a)') 'error - overconstrained rigid body unit'

     Else If (kode ==  646) Then

        Write(nrite,'(/,1x,a)') 'error - overconstrained constraint unit'

     Else If (kode ==  648) Then

        Write(nrite,'(/,1x,a)') 'error - quaternion setup failed'

     Else If (kode ==  650) Then

        Write(nrite,'(/,1x,a)') 'error - failed to find principal axis system'

     Else If (kode ==  655) Then

        Write(nrite,'(/,1x,a)') 'error - FENE bond breaking failure'

     Else If (kode ==  660) Then

        Write(nrite,'(/,1x,a)') 'error - TABBND or PDF bond breaking failure'

     Else If (kode == 1000) Then

        Write(nrite,'(/,1x,a)') 'error - working precision mismatch between FORTRAN90 and MPI implementation'

     Else If (kode == 1001) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> gcheck_vector'

     Else If (kode == 1002) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> gcheck_vector'

     Else If (kode == 1003) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> gisum_vector'

     Else If (kode == 1004) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> gisum_vector'

     Else If (kode == 1005) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> grsum_vector'

     Else If (kode == 1006) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> grsum_vector'

     Else If (kode == 1007) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> gimax_vector'

     Else If (kode == 1008) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> gimax_vector'

     Else If (kode == 1009) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> grmax_vector'

     Else If (kode == 1010) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> grmax_vector'

     Else If (kode == 1011) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in parse_module -> get_record'

     Else If (kode == 1012) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in parse_module -> get_record'

     Else If (kode == 1013) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angles_arrays'

     Else If (kode == 1014) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bonds_arrays'

     Else If (kode == 1015) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in core_shell_module -> allocate_core_shell_arrays'

     Else If (kode == 1016) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in statistics_module -> allocate_statitics_arrays'

     Else If (kode == 1017) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in tethers_module -> allocate_tethers_arrays'

     Else If (kode == 1018) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in constraints_module -> allocate_constraints_arrays'

     Else If (kode == 1019) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in external_field_module -> allocate_external_field_arrays'

     Else If (kode == 1020) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in dihedrals_module -> allocate_dihedrals_arrays'

     Else If (kode == 1021) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in inversions_module -> allocate_inversion_arrays'

     Else If (kode == 1022) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in vdw_module -> allocate_vdw_arrays'

     Else If (kode == 1023) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_arrays'

     Else If (kode == 1024) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in three_body_module -> allocate_three_body_arrays'

     Else If (kode == 1025) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in config_module -> allocate_config_arrays'

     Else If (kode == 1026) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in site_module -> allocate_site_arrays'

     Else If (kode == 1027) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in tersoff_module -> allocate_tersoff_arrays'

     Else If (kode == 1028) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in angles_module -> deallocate_angles_arrays'

     Else If (kode == 1029) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in bonds_module -> deallocate_bonds_arrays'

     Else If (kode == 1030) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in core_shell_module -> deallocate_core_shell_arrays'

     Else If (kode == 1031) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in tethers_module -> deallocate_tethers_arrays'

     Else If (kode == 1032) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in constraints_module -> deallocate_constraints_arrays'

     Else If (kode == 1033) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in dihedrals_module -> deallocate_dihedrals_arrays'

     Else If (kode == 1034) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in inversions_module -> deallocate_inversions_arrays'

     Else If (kode == 1035) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in defects_module -> allocate_defects_arrays'

     Else If (kode == 1036) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in pmf_module -> allocate_pmf_arrays'

     Else If (kode == 1037) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in pmf_module -> deallocate_pmf_arrays'

     Else If (kode == 1038) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in minimise_module -> allocate_minimise_arrays'

     Else If (kode == 1039) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in minimise_module -> deallocate_minimise_arrays'

     Else If (kode == 1040) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in ewald_module -> ewald_allocate_kall_arrays'

     Else If (kode == 1041) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in langevin_module -> langevin_allocate_arrays'

     Else If (kode == 1042) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in rigid_bodies_module -> allocate_rigid_bodies_arrays'

     Else If (kode == 1043) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in rigid_bodies_module -> deallocate_rigid_bodies_arrays'

     Else If (kode == 1044) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> gimin_vector'

     Else If (kode == 1045) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> gimin_vector'

     Else If (kode == 1046) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> grmin_vector'

     Else If (kode == 1047) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> grmin_vector'

     Else If (kode == 1048) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in comms_module -> grsum_matrix'

     Else If (kode == 1049) Then

        Write(nrite,'(/,1x,a)') 'error - deallocation failure in comms_module -> grsum_matrix'

     Else If (kode == 1050) Then

        Write(nrite,'(/,1x,a)') 'error - sorted I/O base communicator not set'

     Else If (kode == 1053) Then

        Write(nrite,'(/,1x,a)') 'error - sorted I/O allocation error'

     Else If (kode == 1056) Then

        Write(nrite,'(/,1x,a)') 'error - unkown write option given to sorted I/O'

     Else If (kode == 1059) Then

        Write(nrite,'(/,1x,a)') 'error - unknown write level given to sorted I/O'

     Else If (kode == 1060) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in statistics_module -> allocate_statitics_connect'

     Else If (kode == 1061) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in statistics_module -> deallocate_statitics_connect'

     Else If (kode == 1063) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in vdw_module -> allocate_vdw_table_arrays'

     Else If (kode == 1066) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in vdw_module -> allocate_vdw_direct_fs_arrays'

     Else If (kode == 1069) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in metal_module -> allocate_metal_table_arrays'

     Else If (kode == 1070) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in ewald_module -> ewald_allocate_kfrz_arrays'

     Else If (kode == 1072) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bond_pot_arrays'

     Else If (kode == 1073) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in bonds_module -> allocate_bond_dst_arrays'

     Else If (kode == 1074) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angl_pot_arrays'

     Else If (kode == 1075) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in angles_module -> allocate_angl_dst_arrays'

     Else If (kode == 1076) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in dihedrals_module -> allocate_dihd_pot_arrays'

     Else If (kode == 1077) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in dihedrals_module -> allocate_dihd_dst_arrays'

     Else If (kode == 1078) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in inversion_module -> allocate_invr_pot_arrays'

     Else If (kode == 1079) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in inversion_module -> allocate_invr_dst_arrays'

     Else If (kode == 1080) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in greenkubo_module -> allocate_greenkubo_arrays'

     Else If (kode == 1081) Then

        Write(nrite,'(/,1x,a)') 'error - allocation failure in dpd_module -> allocate_dpd_arrays'

     Else

        Write(nrite,'(/,1x,a)') 'error - unnamed error found'

     End If

! close all I/O channels

     Close(Unit=nread)
     Close(Unit=nconf)
     Close(Unit=nfield)
     Close(Unit=ntable)
     Close(Unit=nrefdt)
     Close(Unit=nrite)
     Close(Unit=nstats)
     Close(Unit=nrest)
     Close(Unit=nhist)
     Close(Unit=ndefdt)
     Close(Unit=nrdfdt)
     Close(Unit=nzdndt)
     Close(Unit=nrsddt)

  End If

! abort comms

  Call abort_comms()

End Subroutine error
