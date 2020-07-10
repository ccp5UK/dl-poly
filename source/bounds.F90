Module bounds
  Use angles,          Only: angles_type
  Use bonds,           Only: bonds_type
  Use comms,           Only: comms_type
  Use configuration,   Only: configuration_type,&
                             read_config,&
                             scan_config,&
                             IMCON_NOPBC,&
                             IMCON_CUBIC,&
                             IMCON_ORTHORHOMBIC,&
                             IMCON_PARALLELOPIPED,&
                             IMCON_SLAB,&
                             IMCON_TRUNC_OCTO,&
                             IMCON_RHOMBIC_DODEC,&
                             IMCON_HEXAGONAL

  Use constants,       Only: delr_max,&
                             delth_max,&
                             pi,&
                             rt2,&
                             rt3,&
                             zero_plus
  Use constraints,     Only: constraints_type
  Use control,         Only: scan_control,&
                             scan_control_pre
  Use core_shell,      Only: core_shell_type
  Use development,     Only: development_type
  Use dihedrals,       Only: dihedrals_type
  Use domains,         Only: domains_type,&
                             map_domains
  Use electrostatic,   Only: ELECTROSTATIC_NULL,&
                             electrostatic_type
  Use errors_warnings, Only: error,&
                             info,&
                             warning,&
                             check_print_level
  Use ewald,           Only: ewald_type
  Use external_field,  Only: external_field_type,&
                             FIELD_NULL,&
                             FIELD_UMBRELLA
  Use ffield,          Only: scan_field
  Use filename,        Only: file_type
  Use flow_control,    Only: flow_type
  Use four_body,       Only: four_body_type
  Use greenkubo,       Only: greenkubo_type
  Use inversions,      Only: inversions_type
  Use io,              Only: io_type
  Use kim,             Only: kim_type
  Use kinds,           Only: wi,&
                             wp
  Use metal,           Only: metal_type
  Use mpole,           Only: POLARISATION_CHARMM,&
                             mpole_type
  Use msd,             Only: msd_type
  Use neighbours,      Only: neighbours_type
  Use numerics,        Only: dcell
  Use parallel_fft,    Only: adjust_kmax
  Use pmf,             Only: pmf_type
  Use poisson,         Only: poisson_type
  Use rdfs,            Only: rdf_type
  Use rigid_bodies,    Only: rigid_bodies_type
  Use site,            Only: site_type
  Use statistics,      Only: stats_type
  Use tersoff,         Only: tersoff_type,&
                             TERS_TERSOFF,&
                             TERS_KIHS
  Use tethers,         Only: tethers_type
  Use thermostat,      Only: thermostat_type
  Use three_body,      Only: threebody_type
  Use ttm,             Only: ttm_type, ttm_setup_bounds
  Use vdw,             Only: vdw_type
  Use z_density,       Only: z_density_type

  Implicit None

  Private

  Public :: set_bounds

Contains

  Subroutine set_bounds(site, ttm, io, cshell, cons, pmf, stats, thermo, green, devel, &
                        msd_data, met, pois, bond, angle, dihedral, inversion, tether, threebody, zdensity, &
                        neigh, vdws, tersoffs, fourbody, rdf, mpoles, ext_field, rigid, electro, domain, &
                        config, ewld, kim_data, files, flow, comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine to determine various limits as array bounds,
    ! grid sizes, paddings, iterations, etc. as specified in setup
    !
    ! copyright - daresbury laboratory
    !
    ! author    - i.t.todorov december 2016
    ! contrib   - i.j.bush february 2014
    ! contrib   - m.a.seaton june 2014 (VAF)
    ! contrib   - m.a.seaton march 2017 (TTM)
    ! contrib   - i.t.todorov march 2018 (rpad reset if 'strict' & rpad undefined)
    ! contrib   - i.t.todorov june 2018 (spme suggestions)
    ! contrib   - i.t.todorov june 2018 (fdens & mxatms fixes)
    ! contrib   - i.t.todorov august 2018 (rpad refinements)
    ! contrib   - i.t.todorov january 2019 (mxatms increase)
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    ! contrib   - a.m.elena february 2019 (cherry pick 4.09.2)
    ! contrib   - a.m.elena february 2019 (get back to 4.08 magic numbers for mxnode=1)
    ! contrib   - i.t.todorov may 2019 (rpad & rcut feedback, mxatms & mxatdm setting)
    ! contrib   - i.t.todorov july 2019 (SPME b-spline corrected mxatms & buffers,l_trm)
    ! contrib   - i.t.todorov november 2019 (commenting on magic numbers with changes to magic & buffers' sizes)
    ! contrib   - i.t.todorov november 2019 (changing fdens choosing to handle imbalance better)
    ! contrib   - i.t.todorov november 2019 (increasing bond%max_bonds, angle%max_angles, amending domain%mxbfxp)
    ! contrib   - a.m.elena january 2020 avoid division by zero when shifted coulomb is used
    ! contrib   - i.t.todorov january 2020 (25% increasing mxatms cut down estimates and bsplines if on it)
    ! contrib   - i.t.todorov march 2020 (mxatms overflow safety)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(site_type),           Intent(InOut) :: site
    Type(ttm_type),            Intent(InOut) :: ttm
    Type(io_type),             Intent(InOut) :: io
    Type(core_shell_type),     Intent(InOut) :: cshell
    Type(constraints_type),    Intent(InOut) :: cons
    Type(pmf_type),            Intent(InOut) :: pmf
    Type(stats_type),          Intent(InOut) :: stats
    Type(thermostat_type),     Intent(InOut) :: thermo
    Type(greenkubo_type),      Intent(InOut) :: green
    Type(development_type),    Intent(InOut) :: devel
    Type(msd_type),            Intent(InOut) :: msd_data
    Type(metal_type),          Intent(InOut) :: met
    Type(poisson_type),        Intent(InOut) :: pois
    Type(bonds_type),          Intent(InOut) :: bond
    Type(angles_type),         Intent(InOut) :: angle
    Type(dihedrals_type),      Intent(InOut) :: dihedral
    Type(inversions_type),     Intent(InOut) :: inversion
    Type(tethers_type),        Intent(InOut) :: tether
    Type(threebody_type),      Intent(InOut) :: threebody
    Type(z_density_type),      Intent(InOut) :: zdensity
    Type(neighbours_type),     Intent(InOut) :: neigh
    Type(vdw_type),            Intent(InOut) :: vdws
    Type(tersoff_type),        Intent(InOut) :: tersoffs
    Type(four_body_type),      Intent(InOut) :: fourbody
    Type(rdf_type),            Intent(InOut) :: rdf
    Type(mpole_type),          Intent(InOut) :: mpoles
    Type(external_field_type), Intent(InOut) :: ext_field
    Type(rigid_bodies_type),   Intent(InOut) :: rigid
    Type(electrostatic_type),  Intent(InOut) :: electro
    Type(domains_type),        Intent(InOut) :: domain
    Type(configuration_type),  Intent(InOut) :: config
    Type(ewald_type),          Intent(InOut) :: ewld
    Type(kim_type),            Intent(InOut) :: kim_data
    Type(file_type),           Intent(InOut) :: files(:)
    Type(flow_type),           Intent(InOut) :: flow
    Type(comms_type),          Intent(InOut) :: comm

    Character(Len=256) :: message
    Integer, Dimension(3) :: link_cell
    Integer            :: megatm, mtangl, mtbond, mtcons, mtdihd, mtinv, mtrgd, &
         mtshl, mtteth
    Integer(Kind=wi)   :: mxgrid
    Real(Kind=wp), Dimension(10) :: cell_properties
    Real(Kind=wp)      :: ats, dens, cut, dens0, fdvar, &
         padding2, test, xhi, yhi, zhi

    ! scan the FIELD file data

    Call scan_field(megatm, site, neigh%max_exclude, mtshl, &
                    mtcons, mtrgd, mtteth, mtbond, mtangl, mtdihd, mtinv, &
                    ext_field, cshell, cons, pmf, met, bond, angle, dihedral, inversion, tether, threebody, &
                    vdws, tersoffs, fourbody, rdf, mpoles, rigid, kim_data, files, electro, comm)

    ! Get imc_r & set config%dvar

    Call scan_control_pre(config%imc_n, config%dvar, files, comm)

    ! scan CONFIG file data

    Call scan_config(config, megatm, config%dvar, config%levcfg, xhi, yhi, zhi, io, domain, files, comm)

    ! halt execution for unsupported image conditions in DD
    ! checks for some inherited from DL_POLY_2 are though kept

    If (config%imcon == IMCON_TRUNC_OCTO .or. &
        config%imcon == IMCON_RHOMBIC_DODEC .or. &
        config%imcon == IMCON_HEXAGONAL) Call error(514)

    ! scan CONTROL file data

    Call scan_control(tersoffs%cutoff, rigid%max_rigid, config%imcon, config%imc_n, config%cell, &
                      xhi, yhi, zhi, config%mxgana, config%l_ind, electro%nstfce, &
                      ttm, cshell, stats, thermo, green, devel, msd_data, met, pois, bond, angle, dihedral, &
                      inversion, zdensity, neigh, vdws, tersoffs, rdf, mpoles, electro, ewld, kim_data, &
                      files, flow, comm)

    ! check integrity of cell vectors: for cubic cell

    If (config%imcon == IMCON_CUBIC) then

      ats = (Abs(config%cell(1)) + Abs(config%cell(5))) / 2.0_wp
      test = 1.0e-10_wp * ats ! 1.0e-10_wp tolerance in primitive cell type specification of dimensions
      If (any(Abs(config%cell(1:9:4) - ats) > test)) Call error(410)
    end If

    ! check for diagonal cell matrix if appropriate: imcon=1,2

    If (config%imcon /= IMCON_NOPBC .and. config%imcon /= IMCON_PARALLELOPIPED .and. config%imcon /= IMCON_SLAB) Then
      If (Any(Abs(config%cell(2:4)) > zero_plus)) Then
        Call error(410)
      End If
      If (Any(Abs(config%cell(6:8)) > zero_plus)) Then
        Call error(410)
      End If
    End If

    ! calculate dimensional properties of simulation cell
    ! (for use in link-cells) and ttm%volume and define min cell config%width

    Call dcell(config%cell, cell_properties)
    config%width = Min(cell_properties(7), cell_properties(8), cell_properties(9))

    config%volm = cell_properties(10)


    ! check value of cutoff and reset if necessary

    If (config%imcon > 0) Then
      If (config%imcon == IMCON_SLAB) config%width = Min(cell_properties(7), cell_properties(8))

      ! halt program if potential cutoff exceeds the minimum half-cell config%width

      If (neigh%cutoff >= config%width / 2.0_wp) Then
        Call warning(3, neigh%cutoff, config%width / 2.0_wp, 0.0_wp)
        If (.not. devel%l_trm) Then
          Call error(95)
        End If
      End If
    End If

    if (threebody%mxtbp > 0 .and. threebody%cutoff < 1.0e-6_wp) threebody%cutoff = 0.5_wp * neigh%cutoff
    If (fourbody%max_four_body > 0 .and. fourbody%cutoff < 1.0e-6_wp) fourbody%cutoff = 0.5_wp * neigh%cutoff

    ! config%dvar function

    fdvar = config%dvar**1.7_wp ! 1.7_wp arbitrary power factor magnifying densvar effect

    ! config%dvar push of dihedral%max_legend and neigh%max_exclude ranges as the usual suspects

    dihedral%max_legend = Nint(fdvar * Real(dihedral%max_legend, wp))
    neigh%max_exclude = Nint(fdvar * Real(neigh%max_exclude, wp))

    !!! INTRA-LIKE POTENTIAL PARAMETERS !!!

    call setup_potential_parameters(fdvar, comm%mxnode, config, domain, site, vdws, met, tersoffs, &
       cshell, mtshl, cons, mtcons, rigid, mtrgd, tether, mtteth, bond, mtbond, &
       angle, mtangl, dihedral, mtdihd, inversion, mtinv, threebody, fourbody, ext_field)

    !!! GRIDDING PARAMETERS !!!

    call setup_grids(config, neigh, vdws, met, tersoffs, bond, angle, dihedral, inversion, ext_field, &
       electro, ewld, zdensity, rdf, mxgrid)

    ! DD PARAMETERS - by hypercube mapping of MD cell onto machine resources
    ! Dependences: MD cell config%widths (explicit) and machine resources (implicit)

    Call map_domains(config%imc_n, cell_properties(7), cell_properties(8), cell_properties(9), domain, comm)

    Call info(' ', .true.)
    Write (message, '(a,3(i6,1x))') 'node/domain decomposition (x,y,z): ', &
         domain%nx, domain%ny, domain%nz
    Call info(message, .true., level=2)

    ! TTM matters
    padding2 = Huge(1.0_wp) ! Default to huge so fails in mins
    If (ttm%l_ttm) Call ttm_setup_bounds(ttm, config, domain, megatm, padding2)

    ! Link-cell and VNL matters

    call setup_vnl(fdvar, neigh, flow, domain, ewld, devel, config, met, vdws, electro, rdf, &
       tersoffs, threebody, fourbody, cell_properties, link_cell, kim_data, comm, padding2)

    ! decide on MXATMS while reading CONFIG and scan particle density

    Call read_config(config, megatm, config%levcfg, config%l_ind, flow%strict, neigh%cutoff, config%dvar, xhi, yhi, &
         zhi, dens0, dens, io, domain, files, comm)

    Call setup_buffers(fdvar, dens, dens0, megatm, link_cell, mxgrid, config, domain, stats, neigh, &
       green, site, cshell, cons, pmf, rdf, rigid, tether, bond, angle, dihedral, inversion, zdensity, ewld, mpoles, &
       electro%no_elec, msd_data%l_msd, comm)

    ! reset (increase) link-cell maximum (neigh%max_cell)
    ! if tersoff or three- or four-body potentials exist

    If (tersoffs%max_ter > 0 .or. threebody%mxtbp > 0 .or. fourbody%max_four_body > 0) Then
      cut = neigh%cutoff + 1.0e-6_wp ! define cut,
      If (tersoffs%max_ter > 0) cut = Min(cut, tersoffs%cutoff + 1.0e-6_wp)
      If (threebody%mxtbp > 0) cut = Min(cut, threebody%cutoff + 1.0e-6_wp)
      If (fourbody%max_four_body > 0) cut = Min(cut, fourbody%cutoff + 1.0e-6_wp)

      link_cell(1) = Int(domain%nx_recip * cell_properties(7) / cut)
      link_cell(2) = Int(domain%ny_recip * cell_properties(8) / cut)
      link_cell(3) = Int(domain%nz_recip * cell_properties(9) / cut)

      Write (message, '(a,3i6)') "link-cell decomposition 2 (x,y,z): ", link_cell
      Call info(message, .true., level=3)

      If (any(link_cell  < 3)) Call error(305)

      If      (config%imcon == IMCON_NOPBC) Then
        neigh%max_cell = Max(neigh%max_cell, Nint((fdvar**4) * Real(Product(link_cell + 5), wp)))
      Else If (config%imcon == IMCON_SLAB) Then
        neigh%max_cell = Max(neigh%max_cell, Nint((fdvar**3) * Real(Product(link_cell + 5), wp)))
      Else
        neigh%max_cell = Max(neigh%max_cell, Nint((fdvar**2) * Real(Product(link_cell + 5), wp)))
      End If
    End If

    If (any(link_cell < 3)) Call warning(100, 0.0_wp, 0.0_wp, 0.0_wp)

    Write (message, '(a,3i6)') "Final link-cell decomposition (x,y,z): ", link_cell
    Call info(message, .true., level=1)


  end Subroutine set_bounds

  Subroutine setup_potential_parameters(densvar, num_nodes, config, domain, site, vdws, met, &
       tersoffs, cshell, mtshl, cons, mtcons, rigid, mtrgd, tether, mtteth, bond, mtbond, &
       angle, mtangl, dihedral, mtdihd, inversion, mtinv, threebody, fourbody, ext_field)
    !-----------------------------------------------------------------------
    !
    ! dl_poly_4 subroutine to set up potentials
    !
    ! copyright - daresbury laboratory
    !
    !-----------------------------------------------------------------------
    Real(kind=wp),             Intent( In    ) :: densvar
    Integer,                   Intent( In    ) :: num_nodes
    Type(configuration_type),  Intent( InOut ) :: config
    Type(domains_type),        Intent( InOut ) :: domain
    Type(site_type),           Intent( In    ) :: site
    Type(vdw_type),           Intent( InOut ) :: vdws
    Type(metal_type),         Intent( InOut ) :: met
    Type(tersoff_type),        Intent( InOut ) :: tersoffs
    Type(core_shell_type),     Intent( InOut ) :: cshell
    Type(constraints_type),    Intent( InOut ) :: cons
    Type(rigid_bodies_type),   Intent( InOut ) :: rigid
    Type(tethers_type),        Intent( InOut ) :: tether
    Type(bonds_type),          Intent( InOut ) :: bond
    Type(angles_type),         Intent( InOut ) :: angle
    Type(dihedrals_type),      Intent( InOut ) :: dihedral
    Type(inversions_type),     Intent( InOut ) :: inversion
    Type(threebody_type),      Intent( InOut ) :: threebody
    Type(four_body_type),      Intent( InOut ) :: fourbody
    Type(external_field_type), Intent( InOut ) :: ext_field
    Integer,                   Intent( In    ) :: mtshl, mtcons, mtrgd, mtteth, mtbond, mtangl, mtdihd, mtinv

    if (num_nodes > 1) then
      ! maximum number of core-shell units per node

      If (cshell%mxshl > 0) Then
        cshell%mxshl = Max(cshell%mxshl, num_nodes * mtshl)
        cshell%mxshl = (3 * (Nint(densvar * Real(cshell%mxshl, wp)) + num_nodes - 1)) / num_nodes
      End If

      ! maximum number of constraints per node

      If (cons%mxcons > 0) Then
        cons%mxcons = Max(cons%mxcons, num_nodes * mtcons)
        cons%mxcons = (3 * (Nint(densvar * Real(cons%mxcons, wp)) + num_nodes - 1)) / num_nodes
      End If

      ! maximum number of RBs per node

      If (rigid%max_rigid > 0) Then
        rigid%max_rigid = Max(rigid%max_rigid, num_nodes * mtrgd)
        rigid%max_rigid = (3 * (Nint(densvar * Real(rigid%max_rigid, wp)) + num_nodes - 1)) / num_nodes
      End If

      ! dimension of shared atoms arrays for core-shell, constraint and RB units
      ! Max=Max#(members-per-unit)*Max#(units-per-domain)/2
      ! and maximum number of neighbouring domains/nodes in 3D DD (3^3 - 1)
      config%mxlshp = Max(cshell%mxshl, cons%mxcons, (rigid%max_list * rigid%max_rigid) / 2)
      domain%neighbours = 26

      ! maximum number of tethered atoms per node and tether potential parameters
      If (tether%mxteth > 0) Then
        tether%mxteth = Max(tether%mxteth, num_nodes * mtteth)
        tether%mxteth = (3 * (Nint(densvar * Real(tether%mxteth, wp)) + num_nodes - 1)) / num_nodes
        tether%mxpteth = 3
      Else
        tether%mxpteth = 0
      End If

      ! maximum number of chemical bonds per node and bond potential parameters

      If (bond%max_bonds > 0) Then
        bond%max_bonds = Max(bond%max_bonds, num_nodes * mtbond)
        bond%max_bonds = (4 * (Nint(densvar * Real(bond%max_bonds, wp)) + num_nodes - 1)) / num_nodes
        bond%max_param = 4
      Else
        bond%max_param = 0
      End If

      ! maximum number of bond angles per node and angular potential parameters

      If (angle%max_angles > 0) Then
        angle%max_angles = Max(angle%max_angles, num_nodes * mtangl)
        angle%max_angles = (4 * (Nint(densvar * Real(angle%max_angles, wp)) + num_nodes - 1)) / num_nodes
        angle%max_param = 6
      Else
        angle%max_param = 0
      End If

      ! maximum number of torsion angles per node and dihedral potential parameters

      If (dihedral%max_angles > 0) Then
        dihedral%max_angles = Max(dihedral%max_angles, num_nodes * mtdihd)
        dihedral%max_angles = (3 * (Nint(densvar * Real(dihedral%max_angles, wp)) + num_nodes - 1)) / num_nodes
        dihedral%max_angles = dihedral%max_angles + (dihedral%max_angles + 4) / 5 ! allow for 25% higher density
        dihedral%max_param = 7
      Else
        dihedral%max_param = 0
      End If

      ! maximum number of inversions per node and inversion potential parameters

      If (inversion%max_angles > 0) Then
        inversion%max_angles = Max(inversion%max_angles, num_nodes * mtinv)
        inversion%max_angles = (3 * (Nint(densvar * Real(inversion%max_angles, wp)) + num_nodes - 1)) / num_nodes
        inversion%max_angles = inversion%max_angles + (inversion%max_angles + 4) / 5 ! allow for 25% higher density
        inversion%max_param = 3
      Else
        inversion%max_param = 0
      End If

    Else ! nothing is to be shared on one node

      config%mxlshp = 0
      domain%neighbours = 0

      ! maximum number of tether potential parameters
      If (tether%mxteth > 0) Then
        tether%mxpteth = 3
      Else
        tether%mxpteth = 0
      End If

      ! maximum number of bond potential parameters

      If (bond%max_bonds > 0) Then
        bond%max_param = 4
      Else
        bond%max_param = 0
      End If

      ! maximum number of angular potential parameters

      If (angle%max_angles > 0) Then
        angle%max_param = 6
      Else
        angle%max_param = 0
      End If

      ! maximum number of dihedral potential parameters

      If (dihedral%max_angles > 0) Then
        dihedral%max_param = 7
      Else
        dihedral%max_param = 0
      End If

      ! maximum number of inversion potential parameters

      If (inversion%max_angles > 0) Then
        inversion%max_param = 3
      Else
        inversion%max_param = 0
      End If

    end if

    ! if many body interactions exist then their cutoffs must be (threebody%cutoff,fourbody%cutoff,tersoffs%cutoff) > 1.0e-6_wp Angstroms
    ! maximum number of three-body potentials and parameters
    ! TERSOFF DEFINED ELSEWHERE FOR "REASONS" (end of set_bounds)

    If (threebody%mxtbp > 0) Then
      threebody%mx2tbp = (site%mxatyp * (site%mxatyp + 1)) / 2
      threebody%mxtbp  = threebody%mx2tbp * site%mxatyp

      threebody%mxptbp = 5
    Else
      threebody%mx2tbp = 0
      threebody%mxtbp  = 0

      threebody%mxptbp = 0
    End If

    ! maximum number of four-body potentials and parameters

    If (fourbody%max_four_body > 0) Then
      fourbody%mx3fbp        = (site%mxatyp * (site%mxatyp + 1) * (site%mxatyp + 2)) / 6
      fourbody%max_four_body = fourbody%mx3fbp * site%mxatyp

      fourbody%max_param     = 3
    Else
      fourbody%mx3fbp        = 0
      fourbody%max_four_body = 0

      fourbody%max_param     = 0
    End If


    ! maximum number of external field parameters

    If (ext_field%key /= FIELD_NULL) Then
      ext_field%max_param = 6
    Else
      ext_field%max_param = 0
    End If


    if (vdws%no_vdw) then
      vdws%max_param = 0
    else
      vdws%max_vdw = vdws%max_vdw + 1
      vdws%max_param = 7
    end if

    ! maximum number of grid points for metal interactions

    If (met%max_metal > 0) Then
      ! maximum number of metal potentials and parameters
      met%max_metal = met%max_metal + 1
      met%max_param = 9
    Else
      met%max_param = 0
    End If

    ! maximum number of grid points for tersoff interaction arrays

    if (tersoffs%max_ter > 0) then
      ! maximum number of tersoff potentials (tersoffs%max_ter = tersoffs%max_ter) and parameters
      If (tersoffs%key_pot == TERS_TERSOFF) Then
        tersoffs%max_param = 11
      Else If (tersoffs%key_pot == TERS_KIHS) Then
        tersoffs%max_param = 16
      End If
    else
      tersoffs%max_param = 0
    end if

  end Subroutine setup_potential_parameters

  Subroutine setup_grids(config, neigh, vdws, met, tersoffs, bond, angle, dihedral, inversion, ext_field, &
       electro, ewld, zdensity, rdf, mxgrid)
    !-----------------------------------------------------------------------
    !
    ! dl_poly_4 subroutine to set up grids
    !
    ! copyright - daresbury laboratory
    !
    !-----------------------------------------------------------------------
    Type(configuration_type),  Intent( InOut ) :: config
    Type(neighbours_type),     Intent( In    ) :: neigh
    Type(vdw_type),            Intent( InOut ) :: vdws
    Type(metal_type),          Intent( InOut ) :: met
    Type(tersoff_type),        Intent( InOut ) :: tersoffs
    Type(bonds_type),          Intent( InOut ) :: bond
    Type(angles_type),         Intent( InOut ) :: angle
    Type(dihedrals_type),      Intent( InOut ) :: dihedral
    Type(inversions_type),     Intent( InOut ) :: inversion
    Type(external_field_type), Intent( InOut ) :: ext_field
    Type(electrostatic_type),  Intent( InOut ) :: electro
    Type(ewald_type),          Intent( InOut ) :: ewld
    Type(z_density_type),      Intent( InOut ) :: zdensity
    Type(rdf_type),            Intent( InOut ) :: rdf
    Integer,                   Intent(   Out ) :: mxgrid

    ! Set grids for opted intramolecular distribution analysis if unset
    ! SO THEY ARE SWITCHES FOR EXISTENCE TOO

    If (config%mxgana > 0) Then
      If (bond%bin_pdf == -1) Then
        If (bond%bin_tab > 0) Then
          bond%bin_pdf = bond%bin_tab - 4
        Else
          bond%bin_pdf = Nint(bond%rcut / delr_max)
        End If
      End If
      If (angle%bin_adf == -1) Then
        If (angle%bin_tab > 0) Then
          angle%bin_adf = angle%bin_tab - 4
        Else
          angle%bin_adf = Nint(180.0_wp / delth_max)
        End If
      End If
      If (dihedral%bin_adf == -1) Then
        If (dihedral%bin_tab > 0) Then
          dihedral%bin_adf = dihedral%bin_tab - 4
        Else
          dihedral%bin_adf = Nint(360.0_wp / delth_max)
        End If
      End If
      If (inversion%bin_adf == -1) Then
        If (inversion%bin_tab > 0) Then
          inversion%bin_adf = inversion%bin_tab - 4
        Else
          dihedral%bin_adf = Nint(180.0_wp / delth_max)
        End If
      End If
      config%mxgana = Max(bond%bin_pdf, angle%bin_adf, dihedral%bin_adf, inversion%bin_adf)
    End If
    config%mxtana = 0 ! initialise for buffer size purposes, set in read_field

    ! maximum number of rdf potentials (rdf%max_rdf = rdf%max_rdf)
    ! rdf%max_grid - maximum dimension of rdf%rdf and z-density arrays

    if (zdensity%l_collect) then
      zdensity%max_grid = Nint(neigh%cutoff / zdensity%bin_width)
    else
      zdensity%max_grid = 0
    end if

    If (rdf%l_collect) Then
      If ((rdf%l_collect .and. rdf%max_rdf == 0) .and. (vdws%max_vdw > 0 .or. met%max_metal > 0)) &
           rdf%max_rdf = Max(vdws%max_vdw, met%max_metal) ! (vdws,met) == rdf scanning
      rdf%max_grid  = Nint(neigh%cutoff / rdf%rbin)
    Else
      rdf%max_grid  = 0 ! RDF and Z-density function MUST NOT get called!!!
    End If

    ! RDFs particulars for USR (umbrella sampling restraints)

    If (ext_field%key == FIELD_UMBRELLA) Then
      rdf%cutoff_usr   = 0.45_wp * config%width ! 10% reduced half-width (0.45_wp)
      rdf%max_grid_usr = Nint(rdf%cutoff_usr / rdf%rbin) ! allows for up to ~75% system ttm%volume shrinkage
      rdf%cutoff_usr   = Real(rdf%max_grid_usr, wp) * rdf%rbin ! round up and beautify for Andrey Brukhno's sake
    Else
      rdf%cutoff_usr   = 0.0_wp
      rdf%max_grid_usr = 0 ! decider on calling USR RDF
    End If

    ! the number 1004 is minimum default TABLE's gridding - 1000 (equidistant) values
    ! with 2 extra ones on each side (for derivatives) totals.
    ! maximum of all maximum numbers of grid points for all grids - used for mxbuff

    ! maximum number of grid points for bonds

    if (bond%bin_tab > 0) bond%bin_tab = Min(bond%bin_tab, Max(1004, Nint(bond%rcut / delr_max) + 4))

    ! maximum number of grid points for angles

    if (angle%bin_tab > 0) angle%bin_tab = Min(angle%bin_tab, Nint(180.0_wp / delth_max) + 4)

    ! maximum number of grid points for dihedrals

    if (dihedral%bin_tab > 0) dihedral%bin_tab = Min(dihedral%bin_tab, Nint(360.0_wp / delth_max) + 4)

    ! maximum number of grid points for inversions

    if (inversion%bin_tab > 0) inversion%bin_tab = Min(inversion%bin_tab, Nint(180.0_wp / delth_max) + 4)

    ! maximum number of grid points for electrostatics

    if (ewld%direct .or. electro%no_elec) then
      electro%erfc%nsamples = -1
      electro%erfc_deriv%nsamples = -1
    else
      electro%erfc%nsamples = Max(1004, Nint(neigh%cutoff / delr_max) + 4)
      electro%erfc_deriv%nsamples = Max(1004, Nint(neigh%cutoff / delr_max) + 4)
    end if

    ! maximum number of grid points for vdw interactions - overwritten

    if (vdws%no_vdw) then
      vdws%max_grid = -1
    else
      vdws%max_grid = Max(1004, Nint(vdws%cutoff / delr_max) + 4)
    end if

    ! maximum number of grid points for metal interactions

    If (met%max_metal > 0) Then
      met%maxgrid = Max(met%maxgrid, 1004, Nint(met%rcut / delr_max) + 4)
    else
      met%maxgrid = -1
    end If

    ! maximum number of grid points for tersoff interaction arrays

    if (tersoffs%max_ter > 0) then
      tersoffs%max_grid = Max(1004, Nint(tersoffs%cutoff / delr_max) + 4)
    else
      tersoffs%max_grid = -1
    end if

    ! maximum of all maximum numbers of grid points for all grids - used for mxbuff

    mxgrid = Max(config%mxgana, vdws%max_grid, met%maxgrid, zdensity%max_grid, &
         rdf%max_grid, rdf%max_grid_usr, bond%bin_tab, angle%bin_tab, dihedral%bin_tab, &
         inversion%bin_tab, electro%erfc%nsamples, vdws%max_grid, met%maxgrid, tersoffs%max_grid)

  end Subroutine setup_grids

!   Subroutine setup_vnl(fdvar, neigh, flow, domain, ewld, devel, config, met, vdws, electro, rdf, &
!        tersoffs, threebody, fourbody, cell_properties, link_cell, kim_data, comm, ttm_padding)
!     Type(neighbours_type),        Intent( InOut ) :: neigh
!     Type(flow_type),              Intent( In    ) :: flow
!     Type(domains_type),           Intent( In    ) :: domain
!     Type(ewald_type),             Intent( InOut ) :: ewld
!     Type(metal_type),             Intent( In    ) :: met
!     Type(vdw_type),               Intent( In    ) :: vdws
!     Type(electrostatic_type),     Intent( In    ) :: electro
!     Type(rdf_type),               Intent( In    ) :: rdf
    ! Type(tersoff_type),           Intent( In    ) :: tersoffs
    ! Type(threebody_type),         Intent( In    ) :: threebody
    ! Type(four_body_type),         Intent( In    ) :: fourbody
!     Type(development_type),       Intent( In    ) :: devel
!     Type(configuration_type),     Intent( In    ) :: config
!     Type(kim_type),               Intent( In    ) :: kim_data
!     Type(comms_type),             Intent( In    ) :: comm
!     Integer,       Dimension(3),  Intent(   Out ) :: link_cell
!     Real(Kind=wp), Dimension(10), Intent( In    ) :: cell_properties
!     Real(Kind=wp),                Intent( In    ) :: ttm_padding
!     Real(Kind=wp),                Intent( In    ) :: fdvar

!     Real(Kind=wp), Parameter :: minimum_pad = 0.05 ! Ang
!     Real(Kind=wp), Parameter :: minimum_pad_cutoff_frac = 0.005 ! 0.5%
!     Real(Kind=wp), Parameter :: padding_trial_diff = 0.02_wp  ! 2% neigh%padding auto-push (may double)
!     Real(Kind=wp), Parameter :: trial_pct = 0.95_wp  ! (95%) giving 5% slack
!     Real(Kind=wp), Parameter :: delta_r = 1.0e-6_wp
!     Real(Kind=wp)            :: negligible_pad
!     Character(Len=256) :: message, messages(3)
!     Logical            :: no_default_padding
!     Integer,       Dimension(3) :: print_decomp
!     Integer :: bspline_node_check
!     Real(kind=wp) :: padding_tmp
!     Real(Kind=wp) :: test, tol
!     Real(kind=wp) :: cut

!     ! LC and VNL matters
!     ! Reset_padding implies padding in control
!     no_default_padding = neigh%padding <= zero_plus .and. flow%reset_padding

! 5   Continue

!     If (neigh%padding > zero_plus) Then

!       ! define cut adding delta r (1.0e-6_wp the smallest distance we care about)

!       cut = neigh%cutoff + delta_r

!       ! Provide advise on decomposition

!       print_decomp = Int(cell_properties(7:9) / cut)

!       Write (message, '(a,i6,a,3(i0,a))') &
!            'pure cutoff driven limit on largest possible decomposition:', Product(print_decomp), &
!            ' nodes/domains (', print_decomp(1), ',', print_decomp(2), ',', print_decomp(3), ')'
!       Call info(message, .true., level=3)

!       print_decomp = Max(1, print_decomp / 2)

!       Write (message, '(a,i6,a,3(i0,a))') &
!            'pure cutoff driven limit on largest balanced decomposition:', Product(print_decomp), &
!            ' nodes/domains (', print_decomp(1), ',', print_decomp(2), ',', print_decomp(3), ')'
!       Call info(message, .true., level=3)

!     End If

! 10  Continue ! possible neigh%cutoff redefinition...

!     ! Define link-cell cutoff (minimum config%width)

!     neigh%cutoff_extended = neigh%cutoff + neigh%padding

!     ! define cut adding delta r (1.0e-6_wp the smallest distance we care about)

!     cut = neigh%cutoff_extended + delta_r

!     ! Provide advise on decomposition

!     print_decomp = Int(cell_properties(7:9) / cut)

!     Write (message, '(a,i6,a,3(i0,a))') &
!          'cutoffs driven limit on largest possible decomposition:', Product(print_decomp), &
!          ' nodes/domains (', print_decomp(1), ',', print_decomp(2), ',', print_decomp(3), ')'
!     Call info(message, .true., level=3)

!     print_decomp = Max(1, print_decomp / 2)

!     Write (message, '(a,i6,a,3(i0,a))') &
!          'cutoffs driven limit on largest balanced decomposition:', Product(print_decomp), &
!          ' nodes/domains (', print_decomp(1), ',', print_decomp(2), ',', print_decomp(3), ')'
!     Call info(message, .true., level=3)

!     ! calculate link cell dimensions per node

!     link_cell = Int([domain%nx_recip, domain%ny_recip, domain%nz_recip] * cell_properties(7:9) / cut)

!     ! print link cell algorithm and check for violations or...

!     Write (message, '(a,3i6)') "link-cell decomposition 1 (x,y,z): ", link_cell
!     Call info(message, .true., level=3)

!     negligible_pad  = Min(minimum_pad, minimum_pad_cutoff_frac * neigh%cutoff)  ! tolerance

!     if (ewld%bspline > 0) then                                                  ! 2% (w/ SPME/PS) or 4% (w/o SPME/PS)
!       test = padding_trial_diff
!     else
!       test = 2.0_wp * padding_trial_diff
!     end if
!     ! remove delta r (1.0e-6_wp the smallest distance we care about)
!     cut  = Min(config%width / 2.0_wp, domain%nx_recip * cell_properties(7), & ! domain size
!          domain%ny_recip * cell_properties(8), &
!          domain%nz_recip * cell_properties(9)) - delta_r


!     If (Product(link_cell) == 0) Then
!       If (devel%l_trm) Then ! we are prepared to exit gracefully(-:

!         neigh%cutoff = cut  ! - neigh%padding (was zeroed in scan_control)
!         Write (message, '(a)') "real space cutoff reset has occurred, early run termination is due"
!         Call warning(message, .true.)
!         Go To 10
!       Else

!         If (cut < neigh%cutoff) Then
!           Write (message, '(1x,2(a,f0.3),a)') 'user specified or autogenerated padding: ', neigh%padding, &
!                ' , DD+LC suggested maximum value: ', 0.0_wp, ' Angstrom'
!           Call info(message, .true.)
!           Write (message, '(1x,2(a,f0.3),a)') 'user specified cutoff: ', neigh%cutoff, &
!                ' , DD+LC suggested maximum cutoff: ', cut, ' Angstrom'
!           Write (message, '(a)') 'neigh%cutoff <= Min(domain config%width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
!           Call warning(message, .true.)
!           Call error(307)

!         Else ! neigh%padding is defined & in 'no strict' mode

!           If (neigh%padding > zero_plus .and. (.not. flow%strict)) Then ! Re-set neigh%padding with some slack
!             neigh%padding = Min(trial_pct * (cut - neigh%cutoff), test * neigh%cutoff)
!             If (ttm_padding > zero_plus) neigh%padding = Min(neigh%padding, ttm_padding)
!             neigh%padding = Real(Int(100.0_wp * neigh%padding), wp) / 100.0_wp
!             If (neigh%padding < negligible_pad) neigh%padding = 0.0_wp ! Don't bother
!             Go To 10

!           Else
!             Write (message, '(1x,2(a,f0.3),a)') 'user specified padding: ', neigh%padding, &
!                  ' , DD+LC suggested maximum value: ', trial_pct * (cut - neigh%cutoff), ' Angstrom'
!             Call info(message, .true.)
!             Write (message, '(a)') 'neigh%cutoff <= Min(domain config%width) < neigh%cutoff_extended = neigh%cutoff + neigh%padding'
!             Call warning(message, .true.)
!             Call error(307)

!           End If
!         End If
!       End If
!     Else ! push/reset the limits in 'no strict' mode and be somewhat hopeful for undefined rpad in 'strict' when neigh%padding=0
!       If (.not. (met%max_metal == 0 .and. electro%no_elec .and. vdws%no_vdw .and. rdf%max_rdf == 0 .and. kim_data%active)) Then ! 2b link-cells are needed
!         padding_tmp = 0.0_wp
!         If (comm%mxnode == 1 .and. Minval(link_cell) < 2) Then
!           padding_tmp = trial_pct * (0.5_wp * config%width - neigh%cutoff - delta_r)
!           ! remove delta r (1.0e-6_wp the smallest distance we care about)
!         End If

!         If (neigh%padding <= zero_plus) Then ! When neigh%padding is undefined give it some value
!           If (Int(Real(Minval(link_cell), wp) / (1.0_wp + test)) >= 2) Then ! good non-exception
!             neigh%padding = test * neigh%cutoff
!             If (padding_tmp > zero_plus) neigh%padding = Min(neigh%padding, padding_tmp)
!             If (ttm_padding > zero_plus) neigh%padding = Min(neigh%padding, ttm_padding)
!             neigh%padding = Real(Int(100.0_wp * neigh%padding), wp) / 100.0_wp
!             If (neigh%padding > negligible_pad) Go To 10
!           Else ! not so good non-exception
!             neigh%padding = Min(trial_pct * (Min(domain%nx_recip * cell_properties(7) / Real(link_cell(1), wp), &
!                  domain%ny_recip * cell_properties(8) / Real(link_cell(2), wp), &
!                  domain%nz_recip * cell_properties(9) / Real(link_cell(3), wp)) &
!                  - neigh%cutoff - delta_r), test * neigh%cutoff)
!             ! remove delta r (1.0e-6_wp the smallest distance we care about)
!             If (padding_tmp > zero_plus) neigh%padding = Min(neigh%padding, padding_tmp)
!             If (ttm_padding > zero_plus) neigh%padding = Min(neigh%padding, ttm_padding)
!             neigh%padding = Real(Int(100.0_wp * neigh%padding), wp) / 100.0_wp ! round up
!           End If
!         End If

!         If (neigh%padding < negligible_pad) neigh%padding = 0.0_wp ! Don't bother
!       Else
!         neigh%padding = 0.0_wp ! Don't bother
!       End If

!       If (no_default_padding) Then ! forget about it
!         neigh%padding = 0.0_wp
!       Else If (.not. flow%strict) Then ! less hopeful for undefined rpad
!         neigh%padding = trial_pct * neigh%padding
!         neigh%padding = Real(Int(100.0_wp * neigh%padding), wp) / 100.0_wp ! round up
!         If (neigh%padding < negligible_pad) neigh%padding = 0.0_wp ! Don't bother
!       End If

!       neigh%cutoff_extended = neigh%cutoff + neigh%padding ! recalculate neigh%cutoff_extended respectively
!     End If

!     ! Ensure padding is large enough for KIM model

!     If (kim_data%padding_neighbours_required) Then
!       If (neigh%padding < kim_data%influence_distance) Then
!         neigh%padding = kim_data%influence_distance
!         neigh%cutoff_extended = neigh%cutoff + neigh%padding
!       End If
!     End If

!     neigh%unconditional_update = (neigh%padding > zero_plus) ! Determine/Detect conditional VNL updating at start

!     If (any(link_cell < 3)) Call warning(100, 0.0_wp, 0.0_wp, 0.0_wp)

!     ! get total link cells per domain (boundary padding included)
!     ! total link-cells per node/domain is ncells = (link_cell(x)+2)*(link_cell(y)+2)*(link_cell(z)+2)
!     ! allow for more (possible b-spline SPME triggered increase in nlast),
!     ! total link-cells per node/domain is ncells = (link_cell(x)+3)*(link_cell(y)+3)*(link_cell(z)+3)
!     ! allow for thermal expansion of unsettled systems
!     ! total link-cells per node/domain is ncells = (link_cell(x)+4)*(link_cell(y)+4)*(link_cell(z)+4)
!     ! magnify the effect of densvar on neigh%max_cell for different conf%imcon scenarios

!     If      (config%imcon == IMCON_NOPBC) Then
!       neigh%max_cell = Nint(fdvar**4 * Real(Product(link_cell + 4), wp))
!     Else If (config%imcon == IMCON_SLAB) Then
!       neigh%max_cell = Nint(fdvar**3 * Real(Product(link_cell + 4), wp))
!     Else
!       neigh%max_cell = Nint(fdvar**2 * Real(Product(link_cell + 4), wp))
!     End If

!     ! SPME electrostatics particularities

!     ! qlx,qly,qlz - SPME fictional link-cell dimensions postulating that:
!     ! domain%nx <= ewld%fft_dim_a/ewld%bspline, domain%ny <= ewld%fft_dim_b/ewld%bspline, domain%nz <= ewld%fft_dim_c/ewld%bspline.
!     ! Otherwise, this node's b-splines in SPME will need extra 'positive halo'
!     ! that is not on the immediate neighbouring nodes in negative
!     ! directions but beyond them (which may mean self-halo in some cases)

!     ewld%fft_dim_a = ewld%fft_dim_a1
!     ewld%fft_dim_b = ewld%fft_dim_b1
!     ewld%fft_dim_c = ewld%fft_dim_c1

!     ! ewld%bspline = 0 is an indicator for no SPME or Poisson Solver electrostatics in CONTROL

!     If (ewld%bspline /= 0) Then

!       ! ensure (ewld%fft_dim_a,ewld%fft_dim_b,ewld%fft_dim_c) consistency between the DD
!       ! processor grid (map_domains is already called) and the grid
!       ! method or comment out adjustments if using ewald_spme_force~

!       Call adjust_kmax(ewld%fft_dim_a, domain%nx)
!       Call adjust_kmax(ewld%fft_dim_b, domain%ny)
!       Call adjust_kmax(ewld%fft_dim_c, domain%nz)

!       ! Calculate and check ql.

!       print_decomp(1) = Min(link_cell(1), ewld%fft_dim_a / (ewld%bspline * domain%nx))
!       print_decomp(2) = Min(link_cell(2), ewld%fft_dim_b / (ewld%bspline * domain%ny))
!       print_decomp(3) = Min(link_cell(3), ewld%fft_dim_c / (ewld%bspline * domain%nz))

!       If (.not. neigh%unconditional_update) Then
!         ewld%bspline1 = ewld%bspline
!       Else
!         ewld%bspline1 = ewld%bspline + Ceiling((neigh%padding * Real(ewld%bspline, wp)) / neigh%cutoff)

!         ! Redefine ql.

!         print_decomp(1) = Min(link_cell(1), ewld%fft_dim_a / (ewld%bspline1 * domain%nx))
!         print_decomp(2) = Min(link_cell(2), ewld%fft_dim_b / (ewld%bspline1 * domain%ny))
!         print_decomp(3) = Min(link_cell(3), ewld%fft_dim_c / (ewld%bspline1 * domain%nz))
!       End If

!       ! Hard luck, giving up after trying once more

!       If (Product(print_decomp) == 0) Then
!         If ((no_default_padding .eqv. flow%reset_padding) .and. neigh%padding > zero_plus) Then ! defaulted padding must be removed
!           neigh%padding = 0.0_wp
!           no_default_padding = .true.
!           Go To 5
!         Else
!           bspline_node_check = Min(ewld%fft_dim_a / domain%nx, &
!                ewld%fft_dim_b / domain%ny, &
!                ewld%fft_dim_c / domain%nz) - ewld%bspline
!           If (.not. flow%strict) Then
!             If (bspline_node_check == 0 .and. neigh%padding > zero_plus) Then
!               neigh%padding = 0.0_wp
!               no_default_padding = .true.
!               Go To 5
!             Else If (bspline_node_check > 0 .and. neigh%padding > zero_plus) Then
!               padding_tmp = Min(neigh%padding, trial_pct * neigh%cutoff * (Real(bspline_node_check, wp) / Real(ewld%bspline, wp)))
!               padding_tmp = Real(Int(100.0_wp * padding_tmp), wp) / 100.0_wp
!               If (padding_tmp > zero_plus) Then
!                 If (padding_tmp < negligible_pad) padding_tmp = 0.0_wp ! Don't bother
!               End If
!               If (padding_tmp < neigh%padding .and. padding_tmp >= zero_plus) Then
!                 neigh%padding = padding_tmp
!                 Go To 10
!               Else
!                 Write (message, '(2(a,f0.3),a)') 'user specified or autogenerated padding: ', neigh%padding, &
!                      ' , SPME suggested maximum value: ', padding_tmp, ' Angstrom'
!                 Call info(message, .true.)
!               End If
!             End If
!           Else
!             padding_tmp = 0.0_wp
!             If (neigh%padding > zero_plus) Then
!               padding_tmp = trial_pct * neigh%cutoff * (Real(bspline_node_check, wp) / Real(ewld%bspline, wp))
!               padding_tmp = Real(Int(100.0_wp * padding_tmp), wp) / 100.0_wp
!               If (padding_tmp > zero_plus) Then
!                 If (padding_tmp < negligible_pad) padding_tmp = 0.0_wp ! Don't bother
!               End If
!             End If
!             Write (message, '(2(a,f0.3),a)') 'user specified or autogenerated padding: ', neigh%padding, &
!                  ' , SPME suggested maximum value: ', padding_tmp, ' Angstrom'
!             Call info(message, .true., level=2)
!           End If

!           test = Min(Real(ewld%fft_dim_a1, wp) / Real(domain%nx, wp), &
!                Real(ewld%fft_dim_b1, wp) / Real(domain%ny, wp), &
!                Real(ewld%fft_dim_c1, wp) / Real(domain%nz, wp)) / Real(ewld%bspline, wp)
!           tol = Min(Real(ewld%fft_dim_a, wp) / Real(domain%nx, wp), &
!                Real(ewld%fft_dim_b, wp) / Real(domain%ny, wp), &
!                Real(ewld%fft_dim_c, wp) / Real(domain%nz, wp)) / Real(ewld%bspline1, wp)

!           Write (messages(1), '(a,4(i0,a))') &
!                'SPME driven limit on largest possible decomposition:', &
!                (ewld%fft_dim_a / ewld%bspline1) * (ewld%fft_dim_b / ewld%bspline1) * (ewld%fft_dim_c / ewld%bspline1),   &
!                ' nodes/domains (',                                                                                       &
!                ewld%fft_dim_a / ewld%bspline1, ',', ewld%fft_dim_b / ewld%bspline1, ',', ewld%fft_dim_c / ewld%bspline1, &
!                ') for currently specified cutoff (with padding) & Ewald precision (sum parameters)'
!           Write (messages(2), '(a,f0.2,a)') &
!                'SPME suggested factor to decrease currently specified cutoff (with padding) by: ', &
!                1.0_wp / test, ' for currently specified Ewald precision & domain decomposition'
!           Write (messages(3), '(a,f0.2,a)') &
!                'SPME suggested factor to increase current Ewald precision by: ', &
!                (1.0_wp - tol) * 100.0_wp, ' for currently specified cutoff (with padding) & domain decomposition'
!           Call info(messages, 3, .true., level=2)
!         End If

!         If (devel%l_trm) Then
!           print_decomp = link_cell
!           Write (message, '(a)') "SPME partitioning aborted, early run termination is due"
!           Call warning(message, .true.)
!         Else
!           Call error(308)
!         End If
!       End If

!     End If

!     ! reset (increase) link-cell maximum (neigh%max_cell)
!     ! if tersoff or three- or four-body potentials exist

    ! If (tersoffs%max_ter > 0 .or. threebody%mxtbp > 0 .or. fourbody%max_four_body > 0) Then
    !   cut = neigh%cutoff + 1.0e-6_wp ! define cut,
    !   If (tersoffs%max_ter > 0) cut = Min(cut, tersoffs%cutoff + 1.0e-6_wp)
    !   If (threebody%mxtbp > 0) cut = Min(cut, threebody%cutoff + 1.0e-6_wp)
    !   If (fourbody%max_four_body > 0) cut = Min(cut, fourbody%cutoff + 1.0e-6_wp)

    !   link_cell(1) = Int(domain%nx_recip * cell_properties(7) / cut)
    !   link_cell(2) = Int(domain%ny_recip * cell_properties(8) / cut)
    !   link_cell(3) = Int(domain%nz_recip * cell_properties(9) / cut)

    !   Write (message, '(a,3i6)') "link-cell decomposition 2 (x,y,z): ", link_cell
    !   Call info(message, .true., level=3)

    !   If (any(link_cell  < 3)) Call error(305)

    !   If      (config%imcon == IMCON_NOPBC) Then
    !     neigh%max_cell = Max(neigh%max_cell, Nint((fdvar**4) * Real(Product(link_cell + 5), wp)))
    !   Else If (config%imcon == IMCON_SLAB) Then
    !     neigh%max_cell = Max(neigh%max_cell, Nint((fdvar**3) * Real(Product(link_cell + 5), wp)))
    !   Else
    !     neigh%max_cell = Max(neigh%max_cell, Nint((fdvar**2) * Real(Product(link_cell + 5), wp)))
    !   End If
    ! End If

    ! Write (message, '(a,3i6)') "Final link-cell decomposition (x,y,z): ", link_cell
    ! Call info(message, .true., level=1)
!   end Subroutine setup_vnl

  Subroutine setup_buffers(fdvar, maximum_local_density, maximum_domain_density, megatm, link_cell, max_grid, &
       config, domain, stats, neigh, green, site, cshell, cons, pmf, rdf, &
       rigid, tether, bond, angle, dihedral, inversion, zdensity, ewld, mpoles, &
       no_elec, msd, comm)
    Real(Kind=wp),            Intent(In   ) :: fdvar
    Real(Kind=wp),            Intent(In   ) :: maximum_local_density, maximum_domain_density
    Integer,                  Intent(In   ) :: megatm
    Integer, Dimension(3),    Intent(In   ) :: link_cell
    Type(configuration_type), Intent(InOut) :: config
    Type(domains_type),       Intent(InOut) :: domain
    Type(stats_type),         Intent(InOut) :: stats
    Type(neighbours_type),    Intent(InOut) :: neigh
    Type(greenkubo_type),     Intent(In   ) :: green
    Type(site_type),          Intent(In   ) :: site
    Type(core_shell_type),    Intent(In   ) :: cshell
    Type(constraints_type),   Intent(In   ) :: cons
    Type(pmf_type),           Intent(In   ) :: pmf
    Type(rdf_type),           Intent(In   ) :: rdf
    Type(tethers_type),       Intent(In   ) :: tether
    Type(bonds_type),         Intent(In   ) :: bond
    Type(angles_type),        Intent(In   ) :: angle
    Type(dihedrals_type),     Intent(In   ) :: dihedral
    Type(inversions_type),    Intent(In   ) :: inversion
    Type(z_density_type),     Intent(In   ) :: zdensity
    Type(rigid_bodies_type),  Intent(In   ) :: rigid
    Type(ewald_type),         Intent(InOut) :: ewld
    Type(mpole_type),         Intent(In   ) :: mpoles
    Type(comms_type),         Intent(In   ) :: comm
    Logical,                  Intent(In   ) :: no_elec
    Logical,                  Intent(In   ) :: msd
    Integer,                  Intent(In   ) :: max_grid

    Real(Kind=wp) :: outer_surf_vol, inner_surf_vol
    Real(Kind=wp) :: vcell
    ! Real(Kind=Wp) :: xhi, yhi, zhi
    Real(Kind=wp), Dimension(3) :: dims
    Real(Kind=wp) :: fdens, dens
    Real(Kind=wp) :: tol, test, tmp
    Real(Kind=wp), Parameter :: bigint_r = Real(Huge(1), wp)

    ! Create f(fdvar,maximum_local_density,maximum_domain_density) function of density push, maximum 'local' density, maximum domains' density

    If ((comm%mxnode == 1 .or. Any(link_cell < 3)) .or. &
        (config%imcon == IMCON_NOPBC .or. config%imcon == IMCON_SLAB .or. config%imc_n == IMCON_SLAB) .or. &
        (maximum_local_density / maximum_domain_density <= 0.5_wp) .or. (fdvar > 10.0_wp)) Then
      fdens = maximum_domain_density                                    ! for all possibly bad cases resort to max density
    Else
      fdens = fdvar * (0.5_wp * maximum_domain_density + 0.5_wp * maximum_local_density) ! mix 50:50 and push
    End If

    ! Get reasonable to set fdens limit - all particles in one link-cell

    tol   = Real(megatm, wp) / (Real(Product(link_cell), wp) * Real(comm%mxnode, wp))
    fdens = Min(fdens, tol)

    ! density variation affects the link-cell arrays' dimension
    ! more than domains(+halo) arrays' dimensions, in case of
    ! events of extreme collapse in atomic systems (aggregation)

    ! neigh%max_list is the maximum length of link-cell neigh%list (maximum_domain_density * 4/3 pi neigh%cutoff_extended^3)
    ! + 75% extra tolerance - i.e f(maximum_local_density,maximum_domain_density)*(7.5/3)*pi*neigh%cutoff_extended^3

    neigh%max_list = Nint(fdens * (7.5_wp / 3.0_wp) * pi * neigh%cutoff_extended**3)
    neigh%max_list = Min(neigh%max_list, megatm - 1) ! neigh%max_exclude

    If (neigh%max_list < neigh%max_exclude - 1) Then
      Call warning(6, Real(neigh%max_list, wp), Real(neigh%max_exclude, wp), 0.0_wp)
      neigh%max_list = neigh%max_exclude - 1
    End If

    ! get link-cell volume

    vcell = config%volm / (Real(Product(link_cell), wp) * Real(comm%mxnode, wp))

    ! get averaged link-cell particle number, boosted by fdens + 25% extra tolerance

    test = Min(1.25_wp * fdens*vcell, bigint_r) ! bound to largest integer

    ! set dimension of working coordinate arrays using geometrical DD/LC reasoning * fdvar + 25% extra tolerance
    ! and bound to largest integer

    tmp = test * Real(Product(link_cell + 3), wp)
    tmp = Min(tmp, bigint_r)

    config%mxatms = Max(1 , Nint(tmp)) ! Overestimate & then bound down
    If (comm%mxnode == 1 .or. config%imcon == IMCON_NOPBC) Then ! link_cell(x) >= 2 && link_cell(y) >= 2 && link_cell(z) >= 2
      ! maximum of 8 fold increase in of surface thickness (domain+halo) to volume (domain only) as per geometric reasoning
      tmp = 1.25_wp * fdvar * 8.0_wp * Real(megatm, wp)
      tmp = Min(tmp, bigint_r)
      config%mxatms = Min(config%mxatms , Nint(tmp))
    Else If (config%imcon == IMCON_SLAB .or. config%imc_n == IMCON_SLAB) Then ! comm%mxnode >= 4 .or. (link_cell(x) >= 2 && link_cell(y) >= 2)
      ! maximum of 7 fold increase in of surface thickness (domain+halo) to volume (domain only) as per geometric reasoning
      tmp = 1.25_wp * fdvar * 7.0_wp * Real(megatm, wp)
      tmp = Min(tmp, bigint_r)
      config%mxatms = Min(config%mxatms , Nint(tmp))
    Else If (Product(link_cell) < 2) Then ! comm%mxnode >= 8
      ! maximum of 7 fold increase in of surface thickness (domain+halo) to volume (domain only) as per geometric reasoning
      tmp = 1.25_wp * fdvar * 7.0_wp * Real(megatm, wp)
      tmp = Min(tmp, bigint_r)
      config%mxatms = Min(config%mxatms , Nint(tmp))
    Else
      ! maximum of 8 fold increase in of surface thickness (domain+halo) to volume (domain only) as per geometric reasoning
      tmp = test * Real(Product(link_cell + 2) , wp)
      tmp = Min(tmp, bigint_r)
      config%mxatms = Min(config%mxatms , Nint(tmp))
      tmp = 1.25_wp * fdvar * 8.0_wp * Real(megatm , wp)
      tmp = Min(tmp, bigint_r)
      config%mxatms = Min(config%mxatms , Nint(tmp))
    End If

    test = Min(test, Real(config%mxatms, wp) / Real(Product(link_cell + 2), wp))

    ! maximum number of particles per domain (no halo)

    config%mxatdm = Max(1, Nint(test * Real(Product(link_cell + 1), wp)))

    ! SPME b-spline corrected mxatms (extra positive halo)
    ! .hi is the extra requirement to the normal one LC width halo required per direction
    ! with the limiting case of 1 LC per domain sending all/itself to all its neighbours .hi=2
    ! accounted for by the default halo (the product xhi*yhi*zhi<=8); beyond that config%mxatms
    ! must be redefined by the config%mxatdm based density

    If (.not. no_elec) Then
      If (ewld%active) Then
        ewld%kspace%k_vec_dim_real_p_dom = Real(ewld%kspace%k_vec_dim, wp) / &
             Real([domain%nx,domain%ny,domain%nz], wp)
        dims = Max(1.0_wp, Real(ewld%bspline%num_spline_pad, wp) / &
             (ewld%kspace%k_vec_dim_real_p_dom / &
             Real(link_cell, wp))) + 1.0_wp
        ! yhi = Max(1.0_wp, Real(ewld%bspline%num_spline_pad(2), wp) / (Real(ewld%fft_dim_b, wp) / Real(domain%ny, wp) / Real(link_cell(2), wp))) + 1.0_wp
        ! zhi = Max(1.0_wp, Real(ewld%bspline%num_spline_pad(3), wp) / (Real(ewld%fft_dim_c, wp) / Real(domain%nz, wp) / Real(link_cell(3), wp))) + 1.0_wp
        If (Product(dims) > 8.0_wp) Then
          vcell = config%mxatdm / Real(Product(link_cell), wp)
          dims = dims + Real(link_cell, wp)
          ! xhi = xhi + Real(link_cell(1), wp)
          ! yhi = yhi + Real(link_cell(2), wp)
          ! zhi = zhi + Real(link_cell(3), wp)
          ! tmp = vcell * (xhi * yhi * zhi)
          tmp = vcell * Product(dims)
          tmp = Min(tmp, bigint_r)
          config%mxatms = Nint(tmp)
        End If
      End If
    End If

    ! limit mxatdm as it cannot exceed megatm

    config%mxatdm = Min(config%mxatdm, megatm)

    ! maximum number of timesteps in stack arrays

    stats%mxstak = Max(100, stats%mxstak)

    ! maximum number of variables in stack arrays
    ! update number if the MSD option is used, 51=1+27+...+9+9+1+2+2
    ! consult statistic_collect for more information

    stats%mxnstk = 51 + site%mxatyp
    if (msd) stats%mxnstk = stats%mxnstk + 2*config%mxatdm

    ! maximum dimensions of transfer buffers)))

    ! deport_atomic_data & export_atomic_data (+ metal_ld_export.f90
    ! defects_reference_export & statistics_connect_deport if used)
    ! are supposed to be the most MPI/buffer consuming routines

    ! get largest domain's outer surface volume in LC units (as fraction of domain's volume)

    outer_surf_vol = Real(Product(link_cell + 2), wp) / Real((Minval(link_cell) + 2) * Product(link_cell), wp)

    ! exporting, total per atom per direction (times 13 up to 35)

    domain%mxbfxp = 2 * Nint(Real(config%mxatms, wp) * outer_surf_vol) ! included induced dipoles

    ! deporting, total per atom per direction
    ! pass the full LC contents or 0.5*neigh%padding/neigh%cutoff_extended fraction of it (doubled for safety)

    if (neigh%padding > 0.0_wp) then
      dens = outer_surf_vol * neigh%padding / neigh%cutoff_extended
    else
      dens = outer_surf_vol
    end if

    if (comm%mxnode == 1) then
      domain%mxbfdp = 0
      config%mxbfss = 0
      domain%mxbfsh = 0
    else
      ! Per atom contributions
      domain%mxbfdp = 18 + 12 + (neigh%max_exclude + 1)
      ! statistics connect deporting, total per atom per direction
      config%mxbfss = 8

      if (neigh%unconditional_update) domain%mxbfdp = domain%mxbfdp + 3
      if (mpoles%max_mpoles > 0) then
        domain%mxbfdp = domain%mxbfdp + neigh%max_exclude + 1
        if (mpoles%key == POLARISATION_CHARMM) domain%mxbfdp = domain%mxbfdp + neigh%max_exclude + 1 ! CHARMM sends twice the data
      end if
      if (msd) then
        domain%mxbfdp = domain%mxbfdp + 2 * (6 + stats%mxstak)
        config%mxbfss = 2 * (6 + stats%mxstak)
      end if

      domain%mxbfdp = config%mxatdm * domain%mxbfdp
      config%mxbfss = config%mxatdm * config%mxbfss

      ! Bond contribs
      domain%mxbfdp = domain%mxbfdp + &
           3 * green%samp + &
           4 * cshell%mxshl + &
           4 * cons%mxcons + &
           (Sum(pmf%mxtpmf(1:2) + 3) * pmf%mxpmf) + &
           (rigid%max_list + 13) * rigid%max_rigid  + &
           3 * tether%mxteth + &
           4 * bond%max_bonds + &
           5 * angle%max_angles + &
           8 * dihedral%max_angles + &
           6 * inversion%max_angles

      ! Scale by dens
      domain%mxbfdp = Nint(Real(domain%mxbfdp, wp) * dens)
      config%mxbfss = Nint(Real(config%mxbfss, wp) * dens)

      ! Double for send/recv
      domain%mxbfdp = 2 * domain%mxbfdp
      config%mxbfss = 2 * config%mxbfss

      ! get largest domain's inner surface in LC units (as fraction of domain's volume)

      inner_surf_vol = Real(Product(Max(link_cell - 1, 1)), wp) / Real(Minval(link_cell), wp) / Real(Product(link_cell), wp)

      ! shared units updating, per direction

      domain%mxbfsh = Max(2 * cshell%mxshl, 2 * cons%mxcons, rigid%max_list * rigid%max_rigid)
      domain%mxbfsh = 2 * Nint(Real(domain%mxbfsh, wp) * inner_surf_vol)

    end if

    ! largest principal buffer

    config%mxbuff = Max(&
         domain%mxbfdp, &
         35 * domain%mxbfxp, &
         4 * domain%mxbfsh, &
         2 * Product(Int(ewld%kspace%k_vec_dim_real_p_dom)) + 10, &
         stats%mxnstk * stats%mxstak, &
         max_grid, &
         rdf%max_grid, &
         zdensity%max_grid, &
         rigid%max_list * Max(rigid%max_rigid, rigid%max_type), &
         rigid%max_type * (4 + 3 * rigid%max_list), &
         10000)

  end Subroutine setup_buffers

  Subroutine setup_vnl(fdvar, neigh, flow, domain, ewld, devel, config, met, vdws, electro, rdf, &
       tersoffs, threebody, fourbody, cell_properties, link_cell, kim_data, comm, ttm_padding)
    Type(neighbours_type),        Intent( InOut ) :: neigh
    Type(flow_type),              Intent( In    ) :: flow
    Type(domains_type),           Intent( In    ) :: domain
    Type(ewald_type),             Intent( InOut ) :: ewld
    Type(metal_type),             Intent( In    ) :: met
    Type(vdw_type),               Intent( In    ) :: vdws
    Type(electrostatic_type),     Intent( In    ) :: electro
    Type(rdf_type),               Intent( In    ) :: rdf
    Type(tersoff_type),           Intent( In    ) :: tersoffs
    Type(threebody_type),         Intent( In    ) :: threebody
    Type(four_body_type),         Intent( In    ) :: fourbody
    Type(development_type),       Intent( In    ) :: devel
    Type(configuration_type),     Intent( In    ) :: config
    Type(kim_type),               Intent( In    ) :: kim_data
    Type(comms_type),             Intent( In    ) :: comm
    Integer,       Dimension(3),  Intent(   Out ) :: link_cell
    Real(Kind=wp), Dimension(10), Intent( In    ) :: cell_properties
    Real(Kind=wp),                Intent( In    ) :: ttm_padding
    Real(Kind=wp),                Intent( In    ) :: fdvar

    Real(Kind=wp), Parameter                      :: minimum_pad = 0.05 ! Ang
    Real(Kind=wp), Parameter                      :: minimum_pad_cutoff_frac = 0.005 ! 0.5%
    Real(Kind=wp), Parameter                      :: padding_trial_diff = 0.02_wp  ! 2% neigh%padding auto-push (may double)
    Real(Kind=wp)                                 :: cutoff_pct
    Real(Kind=wp), Parameter                      :: trial_pct = 0.95_wp  ! (95%) giving 5% slack
    Real(Kind=wp), Parameter                      :: delta_r = 1.0e-6_wp
    Real(Kind=wp)                                 :: negligible_pad

    Logical                                       :: no_default_padding
    Real(Kind=wp)                                 :: ewald_padding, default_padding, serial_padding, cutoff_padding

    Integer,       Dimension(3)                   :: print_decomp

    Integer                                       :: bspline_node_check
    Integer                                       :: tries
    Real(Kind=wp)                                 :: test, tol

    Real(Kind=wp)                                 :: width_limit
    Real(Kind=wp)                                 :: ewald_limit
    Real(kind=wp), Dimension(3)                   :: lc_det
    Real(kind=wp)                                 :: cut

    Character(Len=256)                            :: message, messages(3)

    ! LC and VNL matters
    ! Reset_padding implies padding in control
    no_default_padding = neigh%padding <= zero_plus .and. flow%reset_padding

    ! Zero potential paddings
    ewald_padding = 0.0_wp
    serial_padding = 0.0_wp
    default_padding = 0.0_wp
    cutoff_padding = 0.0_wp

    If (check_print_level(3)) Then

      ! define cut adding delta r (1.0e-6_wp the smallest distance we care about)

      cut = neigh%cutoff + delta_r

      ! Provide advise on decomposition

      print_decomp = Int(cell_properties(7:9) / cut)

      Write (message, '(a,i6,a,3(i0,a))') &
           'pure cutoff driven limit on largest possible decomposition:', Product(print_decomp), &
           ' nodes/domains (', print_decomp(1), ',', print_decomp(2), ',', print_decomp(3), ')'
      Call info(message, .true., level=3)

      print_decomp = Max(1, print_decomp / 2)

      Write (message, '(a,i6,a,3(i0,a))') &
           'pure cutoff driven limit on largest balanced decomposition:', Product(print_decomp), &
           ' nodes/domains (', print_decomp(1), ',', print_decomp(2), ',', print_decomp(3), ')'
      Call info(message, .true., level=3)

    End If

    negligible_pad  = Min(minimum_pad, minimum_pad_cutoff_frac * neigh%cutoff)  ! tolerance
    width_limit = Min(config%width * 0.5_wp, domain%nx_recip * cell_properties(7), & ! domain size
         domain%ny_recip * cell_properties(8), &
         domain%nz_recip * cell_properties(9))

    if (neigh%cutoff > width_limit) then
      Write (message, '(1x,2(a,f0.3),a)') 'User specified cutoff: ', neigh%cutoff, &
           ' , DD+LC maximum cutoff: ', width_limit, ' Angstrom'
      Call warning(message, .true.)

      if (devel%l_trm) then ! we are prepared to exit gracefully (-:
        Call warning('Fatal cutoff ignored due to termination', .true.)
        neigh%cutoff = width_limit
        neigh%padding = 0.0_wp
        link_cell = [1, 1, 1]
        return
      else
        call error(0, 'Domain decomposition cannot support cutoff')
      end if

    end if

    if (ewld%active) then
      ! 2% (w/ SPME/PS)
      cutoff_pct = padding_trial_diff

      ewld%kspace%k_vec_dim = ewld%kspace%k_vec_dim_cont

      ! ensure (ewld%fft_dim_a,ewld%fft_dim_b,ewld%fft_dim_c) consistency between the DD
      ! processor grid (map_domains is already called) and the grid
      ! method or comment out adjustments if using ewald_spme_force~
      Call adjust_kmax(ewld%kspace%k_vec_dim(1), domain%nx)
      Call adjust_kmax(ewld%kspace%k_vec_dim(2), domain%ny)
      Call adjust_kmax(ewld%kspace%k_vec_dim(3), domain%nz)

      bspline_node_check = Minval(Int(ewld%kspace%k_vec_dim/[domain%nx, domain%ny, domain%nz]))
      print*, bspline_node_check
      if (bspline_node_check == ewld%bspline%num_splines) then
        call warning('LC+DD with SPME grid too small to support padding')
        neigh%padding = 0.0_wp
        no_default_padding = .true.

      else if (bspline_node_check < ewld%bspline%num_splines) then

        Write (message, '(1x,2(a,3(i0.1,1x)),a)') 'User SPME grid: ', ewld%kspace%k_vec_dim(1:3), &
             ', DD+LC minimum grid: ', (ewld%bspline%num_splines+1)*[domain%nx, domain%ny, domain%nz], ' '
        Call warning(message, .true.)

        tol = Real(bspline_node_check, wp) / Real(ewld%bspline%num_splines, wp)
        Write (message, '(a,f0.2,a)') &
             'SPME suggested factor to increase current Ewald precision by: ', &
             (1.0_wp - tol) * 100.0_wp, ' for currently specified cutoff (with padding) & domain decomposition'
        Call info(message, .true.)


        if (devel%l_trm) then ! we are prepared to exit gracefully (-:
          Call warning('Fatal SPME grid error ignored due to termination', .true.)
          neigh%cutoff = width_limit
          neigh%padding = 0.0_wp
          link_cell = [1, 1, 1]
          return
        else if (.not. flow%strict) then
          ewald_padding = Round(trial_pct * neigh%cutoff * &
               (Real(bspline_node_check, wp) / Real(ewld%bspline%num_splines, wp)), 2)
        else
          call error(0, 'Domain decomposition cannot support SPME grid')
        end if
      end if

      if (no_default_padding) then
        ewld%bspline%num_spline_pad = ewld%bspline%num_splines
      else if (neigh%padding > zero_plus) then ! If no user-defined padding assume ~2% difference
        ewld%bspline%num_spline_pad = ewld%bspline%num_splines + &
             Ceiling(cutoff_pct * Real(ewld%bspline%num_splines, wp))
      else
        ewld%bspline%num_spline_pad = ewld%bspline%num_splines + &
             Ceiling((neigh%padding * Real(ewld%bspline%num_splines, wp)) / neigh%cutoff)
      end if

      ewald_limit = Minval(Real(ewld%kspace%k_vec_dim, wp) / ewld%bspline%num_spline_pad)

    else
      !4% (w/o SPME/PS)
      cutoff_pct = 2.0_wp * padding_trial_diff
      ewald_limit = Huge(1.0_wp)
    end if

    if (comm%mxnode == 1 .and. neigh%cutoff < config%width*0.5_wp) &
         serial_padding = trial_pct * (0.5_wp * config%width - neigh%cutoff - delta_r)

    if (.not. no_default_padding) then
      if (ewald_padding < negligible_pad) ewald_padding = Huge(1.0_wp)
      if (default_padding < negligible_pad) default_padding = Huge(1.0_wp)
      if (serial_padding < negligible_pad) serial_padding = Huge(1.0_wp)
      if (cutoff_padding < negligible_pad) cutoff_padding = Huge(1.0_wp)
      if (neigh%padding < zero_plus) then
        neigh%padding = Huge(1.0_wp)
        cutoff_padding = cutoff_pct * neigh%cutoff
      end if

      write(messages(1), '(7(1X,A15,1X))') "neigh%padding", "ewald_padding", "ttm_padding", "default_padding", &
           "serial_padding", "cutoff_padding", "width_limit"
      write(messages(2), '(7(1X,F15.8,1X))') neigh%padding, ewald_padding, ttm_padding, default_padding, &
           serial_padding, cutoff_padding, width_limit - neigh%cutoff
      call info(messages, 2, .true., level=3)
      neigh%padding = Min(neigh%padding, ewald_padding, ttm_padding, default_padding, serial_padding, cutoff_padding, &
           width_limit - neigh%cutoff) ! Cap

      neigh%padding = Round(neigh%padding, 2)
      if (neigh%padding < negligible_pad) neigh%padding = 0.0_wp
    end if

    If (kim_data%padding_neighbours_required) neigh%padding = max(neigh%padding, kim_data%influence_distance)

    If (met%max_metal == 0 .and. electro%no_elec .and. vdws%no_vdw .and. rdf%max_rdf == 0 .and. kim_data%active) then
      call warning('No metals, electrostatics, vdws, rdfs, but kim active: padding reset to zero',.true.)
      neigh%padding = 0.0_wp ! 2b link-cells are needed
    end If

    ! Define link-cell cutoff (minimum config%width)

    neigh%cutoff_extended = neigh%cutoff + neigh%padding

    ! define cut adding delta r (1.0e-6_wp the smallest distance we care about)

    cut = neigh%cutoff_extended + delta_r

    ! Determining factor for link-cell

    lc_det = cell_properties(7:9) / cut !Min(ewald_limit, Minval(cell_properties(7:9) / cut))

    ! calculate link cell dimensions per node
    link_cell = Int([domain%nx_recip, domain%ny_recip, domain%nz_recip] * lc_det)

    ! print link cell algorithm and check for violations or...

    Write (message, '(a,3i6)') "Link-cell decomposition 1 (x,y,z): ", link_cell
    Call info(message, .true., level=3)

    bad_lc: if (Product(link_cell) == 0) then
      strict: if (flow%strict) then

        Write (message, '(1x,2(a,f0.3),a)') 'User specified padding: ', neigh%padding, &
             ' , DD+LC suggested maximum value: ', trial_pct * (width_limit - neigh%cutoff), ' Angstrom'
        Call info(message, .true.)
        Call error(0, 'User padding too large')

      else if (.not. no_default_padding) then

        tries = 0

        do while (Product(link_cell) < 0 .and. tries < 3 .and. neigh%padding > zero_plus)
          tries = tries + 1

          neigh%padding = Round(Min(trial_pct * (width_limit - neigh%cutoff - delta_r), &
               test * neigh%cutoff, &
               ttm_padding, &
               ewald_limit), 2)
          If (neigh%padding < negligible_pad) neigh%padding = 0.0_wp ! Don't bother

          ewld%bspline%num_spline_pad = ewld%bspline%num_splines + &
               Ceiling((neigh%padding * Real(ewld%bspline%num_splines, wp)) / neigh%cutoff)
          ewald_limit = Minval(Real(ewld%kspace%k_vec_dim, wp) / ewld%bspline%num_splines)

          ! Define link-cell cutoff (minimum config%width)

          neigh%cutoff_extended = neigh%cutoff + neigh%padding

          ! define cut adding delta r (1.0e-6_wp the smallest distance we care about)

          cut = neigh%cutoff_extended + delta_r

          lc_det = Min(ewald_limit, Minval(cell_properties(7:9) / cut))

          ! calculate link cell dimensions per node
          link_cell = Int([domain%nx_recip, domain%ny_recip, domain%nz_recip] * lc_det)
        end do

        if (Product(link_cell) == 0) call error(0, 'Link cell padding generation failure')
      end if strict
    else

      if (.not. flow%strict) then
        neigh%padding = Round(trial_pct*neigh%padding, 2)
        neigh%cutoff_extended = neigh%cutoff + neigh%padding
      end if


    end if bad_lc

    neigh%unconditional_update = (neigh%padding > zero_plus) ! Determine/Detect conditional VNL updating at start

    ! get total link cells per domain (boundary padding included)
    ! total link-cells per node/domain is ncells = (link_cell(x)+2)*(link_cell(y)+2)*(link_cell(z)+2)
    ! allow for more (possible b-spline SPME triggered increase in nlast),
    ! total link-cells per node/domain is ncells = (link_cell(x)+3)*(link_cell(y)+3)*(link_cell(z)+3)
    ! allow for thermal expansion of unsettled systems
    ! total link-cells per node/domain is ncells = (link_cell(x)+4)*(link_cell(y)+4)*(link_cell(z)+4)
    ! magnify the effect of densvar on neigh%max_cell for different conf%imcon scenarios

    If (config%imcon == IMCON_NOPBC) Then
      neigh%max_cell = Nint(fdvar**4 * Real(Product(link_cell + 4), wp))
    Else If (config%imcon == IMCON_SLAB) Then
      neigh%max_cell = Nint(fdvar**3 * Real(Product(link_cell + 4), wp))
    Else
      neigh%max_cell = Nint(fdvar**2 * Real(Product(link_cell + 4), wp))
    End If

  Contains

    Pure Function round(val, prec)
      ! Round a value to given precision
      Real(kind=wp) :: round
      Real(Kind=wp), Intent(In   ) :: val
      Integer,       Intent(In   ) :: prec
      round = Real(Int(10.0_wp**(prec) * val), wp) / (10.0_wp**(prec))
    end Function round


  end Subroutine setup_vnl


End Module bounds
