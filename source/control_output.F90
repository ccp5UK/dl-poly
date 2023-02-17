module control_output
  !!-----------------------------------------------------------------------
  !!
  !! Module to handle printing of new parameters
  !!
  !! copyright - daresbury laboratory
  !! author - j.wilkins august 2022
  !!-----------------------------------------------------------------------

  Use angles,                                 Only: angles_type
  Use angular_distribution,                   Only: adf_type
  Use bonds,                                  Only: bonds_type
  Use comms,                    Only: comms_type
  Use configuration,            Only: IMCON_NOPBC,&
                                      IMCON_SLAB,&
                                      configuration_type
  Use constants,                Only: pi,&
                                      tenunt,&
                                      zero_plus
  Use constraints,              Only: constraints_type
  Use control_parameters,       Only: write_param
  Use coord,                    Only: coord_type
  Use core_shell,               Only: core_shell_type
  Use defects,                  Only: defects_type
  Use development,              Only: development_type
  Use dihedrals,                Only: dihedrals_type
  Use electrostatic,            Only: ELECTROSTATIC_COULOMB,&
                                      ELECTROSTATIC_COULOMB_FORCE_SHIFT,&
                                      ELECTROSTATIC_COULOMB_REACTION_FIELD,&
                                      ELECTROSTATIC_DDDP,&
                                      ELECTROSTATIC_EWALD,&
                                      ELECTROSTATIC_NULL,&
                                      electrostatic_type
  Use errors_warnings,          Only: check_print_level,&
                                      error,&
                                      error_units,&
                                      info,&
                                      get_print_level,&
                                      warning
  Use ewald,                    Only: ewald_type
  Use filename,                 Only: &
                                      FILE_CONFIG, FILE_FIELD, FILE_HISTORF, FILE_HISTORY, &
                                      FILE_MSD, FILE_OUTPUT, FILE_RDF, FILE_REVCON, FILE_REVIVE, &
                                      FILE_REVOLD, FILE_STATS, FILE_TABANG, FILE_TABBND, &
                                      FILE_TABDIH, FILE_TABEAM, FILE_TABINV, FILE_TABVDW, file_type
  Use flow_control,             Only: DFTB,&
                                      RESTART_KEY_CLEAN,&
                                      RESTART_KEY_NOSCALE,&
                                      RESTART_KEY_OLD,&
                                      RESTART_KEY_SCALE,&
                                      flow_type
  Use four_body,                Only: four_body_type
  Use greenkubo,                Only: greenkubo_type
  Use hash,                     Only: MAX_KEY
  Use impacts,                  Only: impact_type
  Use inversions,               Only: inversions_type
  Use io,                       Only: &
                                      IO_READ_DIRECT, IO_READ_MASTER, IO_READ_MPIIO, &
                                      IO_WRITE_SORTED_DIRECT, IO_WRITE_SORTED_MASTER, &
                                      IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT, &
                                      IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, &
                                      io_get_parameters, io_set_parameters, io_type
  Use kim,                      Only: kim_type
  Use kinds,                    Only: STR_LEN,&
                                      wp
  Use metal,                    Only: metal_type
  Use minimise,                 Only: MIN_DISTANCE,&
                                      MIN_ENERGY,&
                                      MIN_FORCE,&
                                      minimise_type
  Use mpole,                    Only: POLARISATION_CHARMM,&
                                      POLARISATION_DEFAULT,&
                                      mpole_type
  Use msd,                      Only: msd_type
  Use neighbours,               Only: neighbours_type
  Use numerics,                 Only: dcell,&
                                      seed_type
  Use parse,                    Only: get_line,&
                                      get_word,&
                                      lower_case,&
                                      word_2_real
  Use plumed,                   Only: plumed_type
  Use pmf,                      Only: pmf_type
  Use rdfs,                     Only: rdf_type
  Use rsds,                     Only: rsd_type
  Use statistics,               Only: stats_type
  Use tersoff,                  Only: tersoff_type
  Use thermostat,               Only: &
                                      CONSTRAINT_NONE, CONSTRAINT_SEMI_ORTHORHOMBIC, &
                                      CONSTRAINT_SURFACE_AREA, CONSTRAINT_SURFACE_TENSION, &
                                      DPD_FIRST_ORDER, DPD_NULL, DPD_SECOND_ORDER, &
                                      ENS_NPT_BERENDSEN, ENS_NPT_BERENDSEN_ANISO, &
                                      ENS_NPT_LANGEVIN, ENS_NPT_LANGEVIN_ANISO, ENS_NPT_MTK, &
                                      ENS_NPT_MTK_ANISO, ENS_NPT_NOSE_HOOVER, &
                                      ENS_NPT_NOSE_HOOVER_ANISO, ENS_NVE, ENS_NVT_ANDERSON, &
                                      ENS_NVT_BERENDSEN, ENS_NVT_EVANS, ENS_NVT_GENTLE, &
                                      ENS_NVT_LANGEVIN, ENS_NVT_LANGEVIN_INHOMO, ENS_NVT_NOSE_HOOVER,&
                                      PSEUDO_LANGEVIN_DIRECT, PSEUDO_LANGEVIN, PSEUDO_GAUSSIAN,&
                                      PSEUDO_DIRECT, thermostat_type
  Use three_body,               Only: threebody_type
  Use timer,                    Only: timer_type
  Use trajectory,               Only: TRAJ_KEY_COMPRESSED,&
                                      TRAJ_KEY_COORD,&
                                      TRAJ_KEY_COORD_VEL,&
                                      TRAJ_KEY_COORD_VEL_FORCE,&
                                      trajectory_type
  Use ttm,                      Only: &
                                      TTM_BC_DIRICHLET, TTM_BC_DIRICHLET_XY, TTM_BC_NEUMANN, &
                                      TTM_BC_PERIODIC, TTM_BC_ROBIN, TTM_BC_ROBIN_XY, TTM_CE_CONST, &
                                      TTM_CE_CONST_DYN, TTM_CE_LINEAR, TTM_CE_LINEAR_DYN, &
                                      TTM_CE_TABULATED, TTM_CE_TANH, TTM_CE_TANH_DYN, TTM_DE_CONST, &
                                      TTM_DE_METAL, TTM_DE_RECIP, TTM_DE_TABULATED, TTM_EPVAR_HETERO, &
                                      TTM_EPVAR_HOMO, TTM_EPVAR_NULL, TTM_KE_CONST, TTM_KE_DRUDE, &
                                      TTM_KE_INFINITE, TTM_KE_TABULATED, TTM_SDEPO_EXP, &
                                      TTM_SDEPO_FLAT, TTM_SDEPO_GAUSS, TTM_SDEPO_NULL, &
                                      TTM_TDEPO_DELTA, TTM_TDEPO_EXP, TTM_TDEPO_GAUSS, &
                                      TTM_TDEPO_PULSE, ttm_type
  Use units,                    Only: &
                                      atomic_units, convert_units, current_units => out_units, &
                                      hartree_units, internal_units, kb_units, kcal_units, &
                                      kj_units, set_out_units, set_timestep, si_units, units_scheme, &
                                      to_out_units
  Use vdw,                      Only: MIX_FENDER_HALSEY,&
                                      MIX_FUNCTIONAL,&
                                      MIX_HALGREN,&
                                      MIX_HOGERVORST,&
                                      MIX_LORENTZ_BERTHELOT,&
                                      MIX_NULL,&
                                      MIX_TANG_TOENNIES,&
                                      MIX_WALDMAN_HAGLER,&
                                      vdw_type
  Use z_density,                Only: z_density_type

  Implicit None

  Private

  Public :: write_parameters

Contains

  Subroutine write_parameters(devel, tmr, seed, io_data, files, neigh, config, &
                              link_cell, flow, stats, thermo, ttm, mpoles, vdws, electro, cshell, &
                              ewld, met, impa, minim, plume, cons, pmf, bond, angle, dihedral, &
                              inversion , msd_data, rdf, vaf, zdensity, adf, coords, defect, traj, &
                              displacement)
    Type(development_type),      Intent(In   ) :: devel
    Type(timer_type),            Intent(In   ) :: tmr
    Type(seed_type),             Intent(In   ) :: seed
    Type(io_type),               Intent(In   ) :: io_data
    Type(file_type),             Intent(In   ) :: files(:)
    Type(neighbours_type),       Intent(In   ) :: neigh
    Type(configuration_type),    Intent(In   ) :: config
    Integer, Dimension(3),       Intent(In   ) :: link_cell
    Type(flow_type),             Intent(In   ) :: flow
    Type(stats_type),            Intent(In   ) :: stats
    Type(thermostat_type),       Intent(In   ) :: thermo
    Type(ttm_type),              Intent(In   ) :: ttm
    Type(mpole_type),            Intent(In   ) :: mpoles
    Type(vdw_type),              Intent(In   ) :: vdws
    Type(electrostatic_type),    Intent(In   ) :: electro
    Type(core_shell_type),       Intent(In   ) :: cshell
    Type(ewald_type),            Intent(In   ) :: ewld
    Type(metal_type),            Intent(In   ) :: met
    Type(impact_type),           Intent(In   ) :: impa
    Type(minimise_type),         Intent(In   ) :: minim
    Type(plumed_type),           Intent(In   ) :: plume
    Type(constraints_type),      Intent(In   ) :: cons
    Type(pmf_type),              Intent(In   ) :: pmf
    Type(bonds_type),            Intent(In   ) :: bond
    Type(angles_type),           Intent(In   ) :: angle
    Type(dihedrals_type),        Intent(In   ) :: dihedral
    Type(inversions_type),       Intent(In   ) :: inversion
    Type(msd_type),              Intent(In   ) :: msd_data
    Type(rdf_type),              Intent(InOut) :: rdf
    Type(greenkubo_type),        Intent(In   ) :: vaf
    Type(z_density_type),        Intent(InOut) :: zdensity
    Type(adf_type),              Intent(In   ) :: adf
    Type(coord_type),            Intent(In   ) :: coords
    Type(defects_type),          Intent(In   ) :: defect(:)
    Type(trajectory_type),       Intent(In   ) :: traj
    Type(rsd_type),              Intent(In   ) :: displacement

    Character(Len=80) :: banner(6)

    Write (banner(1), '(a)') ''
    Write (banner(2), '(a)') '#'//Repeat('*', 79)
    Write (banner(3), '(a4,a72,a4)') '#** ', 'title:'//Repeat(' ', 66), ' ***'
    Write (banner(4), '(a4,a72,a4)') '#** ', config%sysname, ' ***'
    Write (banner(5), '(a)') '#'//Repeat('*', 79)
    Write (banner(6), '(a)') ''
    Call info(banner, 6, .true.)

    If (check_print_level(1)) Call write_devel(devel, tmr, seed)
    If (check_print_level(1)) Call write_io(io_data, files)
    If (check_print_level(1)) Call write_units()
    If (check_print_level(1)) Call write_system_parameters(flow, config, stats, thermo, impa, minim, plume, cons, pmf)
    If (check_print_level(1)) Call write_ensemble(thermo, ttm)
    If (check_print_level(1)) Call write_forcefield(link_cell, neigh, vdws, electro, ewld, mpoles, cshell, met)
    If (ttm%l_ttm .and. check_print_level(1)) Call write_ttm(thermo, ttm)
    If (check_print_level(1)) Call write_bond_analysis(stats, flow, bond, angle, dihedral, inversion)
    If (check_print_level(1)) &
      Call write_structure_analysis(stats, msd_data, rdf, vaf, zdensity, adf, coords, traj, defect, displacement)
    Call info('', .true.)

  End Subroutine write_parameters

  Subroutine write_io(io_data, files)
    Type(io_type),   Intent(In   ) :: io_data
    Type(file_type), Intent(In   ) :: files(:)

    Character(len=256) :: message

    !    Integer :: io_read, io_write, procs_write, procs_read, batch_write, batch_read

    Call info('', .true.)
    Call info('I/O Parameters: ', .true.)

    Select Case (io_data%method_read)
    Case (IO_READ_MPIIO)
      Call write_param('I/O read method', 'parallel by using MPI-I/O', indent=1)

    Case (IO_READ_DIRECT)
      Call write_param('I/O read method', 'parallel by using direct access', indent=1)

    Case (IO_READ_MASTER)
      Call write_param('I/O read method', 'serial by using a single master process', indent=1)

    End Select

    Call write_param('I/O readers', io_data%n_io_procs_read, indent=1)

    Select Case (io_data%method_write)
    Case (IO_READ_MPIIO, IO_READ_DIRECT)
      Call write_param('I/O read batch size (assumed) (Bytes)', io_data%batch_size_read, indent=1, level=3)
    End Select


    Call write_param('I/O read buffer size (Bytes)', io_data%buffer_size_read, indent=1, level=3)

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_MPIIO)
      Call write_param('I/O write method', 'parallel by using MPI-I/O', indent=1)

    Case (IO_WRITE_SORTED_DIRECT, IO_WRITE_UNSORTED_DIRECT)
      Call write_param('I/O write method', 'parallel by using direct access', indent=1)
      Call warning('  in parallel this I/O write method has portability issues', .true.)

    Case (IO_WRITE_SORTED_MASTER, IO_WRITE_UNSORTED_MASTER)
      Call write_param('I/O write method', 'serial by using a single master process', indent=1)

    End Select

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MASTER, IO_WRITE_SORTED_MPIIO, IO_WRITE_SORTED_DIRECT)
      Call write_param('Data sorting', .true., indent=2, on_level=3)

    Case (IO_WRITE_UNSORTED_MASTER, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_UNSORTED_DIRECT)
      Call write_param('Data sorting', .false., indent=2)

    End Select

    Call write_param('I/O writers', io_data%n_io_procs_write, indent=1)

    Select Case (io_data%method_write)
    Case (IO_WRITE_SORTED_MPIIO, IO_WRITE_UNSORTED_MPIIO, IO_WRITE_SORTED_DIRECT, IO_WRITE_UNSORTED_DIRECT)
      Call write_param('I/O write batch size (assumed) (Bytes)', io_data%batch_size_write, indent=1, level=3)
    End Select

    Call write_param('I/O write buffer size (Bytes)', io_data%buffer_size_write, indent=1, level=3)

    Call write_param('I/O parallel error checking', io_data%global_error_check, off_level=3, indent=1)

    Call info('', .true.)
    Call info('File outputs (if not default): ', .true.)
    If (files(FILE_OUTPUT)%filename /= "OUTPUT" .or. check_print_level(3)) &
         Call write_param('OUTPUT file', files(FILE_OUTPUT)%filename, indent=1)
    If (files(FILE_CONFIG)%filename /= "CONFIG" .or. check_print_level(3)) &
         Call write_param('CONFIG file', files(FILE_CONFIG)%filename, indent=1)
    If (files(FILE_FIELD)%filename /= "FIELD" .or. check_print_level(3)) &
         Call write_param('FIELD file', files(FILE_FIELD)%filename, indent=1)
    If (files(FILE_STATS)%filename /= "STATIS" .or. check_print_level(3)) &
         Call write_param('STATIS file', files(FILE_STATS)%filename, indent=1)
    If (files(FILE_HISTORY)%filename /= "HISTORY" .or. check_print_level(3)) &
         Call write_param('HISTORY file', files(FILE_HISTORY)%filename, indent=1)
    If (files(FILE_HISTORF)%filename /= "HISTORF" .or. check_print_level(3)) &
         Call write_param('HISTORF file', files(FILE_HISTORF)%filename, indent=1)
    If (files(FILE_REVIVE)%filename /= "REVIVE" .or. check_print_level(3)) &
         Call write_param('REVIVE file', files(FILE_REVIVE)%filename, indent=1)
    If (files(FILE_REVCON)%filename /= "REVCON" .or. check_print_level(3)) &
         Call write_param('REVCON file', files(FILE_REVCON)%filename, indent=1)
    If (files(FILE_REVOLD)%filename /= "REVOLD" .or. check_print_level(3)) &
         Call write_param('REVOLD file', files(FILE_REVOLD)%filename, indent=1)
    If (files(FILE_RDF)%filename /= 'RDFDAT' .or. check_print_level(3)) &
         Call write_param('RDF file', files(FILE_RDF)%filename, indent=1)
    If (files(FILE_MSD)%filename /= 'MSDTMP' .or. check_print_level(3)) &
         Call write_param('MSD file', files(FILE_MSD)%filename, indent=1)
    If (files(FILE_TABBND)%filename /= 'TABBND' .or. check_print_level(3)) &
         Call write_param('TABBND file', files(FILE_TABBND)%filename, indent=1)
    If (files(FILE_TABANG)%filename /= 'TABANG' .or. check_print_level(3)) &
         Call write_param('TABANG file', files(FILE_TABANG)%filename, indent=1)
    If (files(FILE_TABDIH)%filename /= 'TABDIH' .or. check_print_level(3)) &
         Call write_param('TABDIH file', files(FILE_TABDIH)%filename, indent=1)
    If (files(FILE_TABINV)%filename /= 'TABINV' .or. check_print_level(3)) &
         Call write_param('TABINV file', files(FILE_TABINV)%filename, indent=1)
    If (files(FILE_TABVDW)%filename /= 'TABVDW' .or. check_print_level(3)) &
         Call write_param('TABVDW file', files(FILE_TABVDW)%filename, indent=1)
    If (files(FILE_TABEAM)%filename /= 'TABEAM' .or. check_print_level(3)) &
         Call write_param('TABEAM file', files(FILE_TABEAM)%filename, indent=1)

  End Subroutine write_io

  Subroutine write_units()
    !!-----------------------------------------------------------------------
    !!
    !! Write current units information
    !!
    !! copyright - daresbury laboratory
    !! author - j.wilkins april 2022
    !!-----------------------------------------------------------------------

    Call info('', .true.)
    Call info('Current Units Scheme:', .true.)
    Call write_param('Length units', current_units%length_unit%abbrev, indent=1)
    Call write_param('Time units', current_units%time_unit%abbrev, indent=1)
    Call write_param('Mass units', current_units%mass_unit%abbrev, indent=1)
    Call write_param('Charge units', current_units%charge_unit%abbrev, indent=1)
    Call write_param('Energy units', current_units%energy_unit%abbrev, indent=1)
    Call write_param('Temperature units', current_units%temp_unit%abbrev, indent=1)
    Call write_param('Current units', current_units%current_unit%abbrev, indent=1)
    Call write_param('Luminosity units', current_units%luminosity_unit%abbrev, indent=1, level=3)
    Call write_param('Angle units', current_units%angle_unit%abbrev, indent=1)
    Call write_param('Pressure units', current_units%pressure_unit%abbrev, indent=1)
    Call write_param('Force units', current_units%force_unit%abbrev, indent=1)
    Call write_param('Velocity units', current_units%velocity_unit%abbrev, indent=1)
    Call write_param('Power units', current_units%power_unit%abbrev, indent=1)
    Call write_param('Surface Tension units', current_units%surf_ten_unit%abbrev, indent=1)
    Call write_param('EMF units', current_units%emf_unit%abbrev, indent=1)

  end Subroutine write_units

  Subroutine write_devel(devel, tmr, seed)
    Type(development_type),      Intent(In) :: devel
    Type(timer_type),            Intent(In) :: tmr
    Type(seed_type),             Intent(In) :: seed

    Integer :: i,j,k,l,ij,kl
    Character(Len=STR_LEN) :: message

    Call info('', .true.)
    Call info('Software Params:', .true.)
    Call write_param('Verbosity', get_print_level(), indent=1)
    Call write_param('Unsafe comms', devel%l_fast, indent=1, off_level=3)
    Call write_param('Output final energy', devel%l_eng, indent=1, off_level=3)
    Call write_param('Write ASCII REVIVE', devel%l_rout, indent=1, off_level=3)
    Call write_param('Read ASCII REVOLD', devel%l_rin, indent=1, off_level=3)
    if (devel%r_dis > 0.0_wp) Call write_param('Initial minimum separation', devel%r_dis, indent=1)

    Call write_param('Max timer depth', tmr%max_depth, indent=1, level=2)
    Call write_param('Per process timings', tmr%proc_detail, indent=1, off_level=3)

    Call write_param('Job time', tmr%job, 's', 's', indent=1)
    Call write_param('File closing time', tmr%clear_screen, 's', 's', indent=1, level=3)

    if (seed%defined) then

      ij = Mod(Abs(seed%seed(1)), 31328)
      i = Mod(ij / 177, 177) + 2;
      j = Mod(ij, 177) + 2;
      kl = Mod(Abs(seed%seed(2)), 30081)
      k = Mod(kl / 169, 178) + 1
      l = Mod(kl, 169)

      Write (message, "(a, '[', 3(i0.1, ',', 1X), i0.1, ']')") '  Uniform random seed (proc0): ', i,j,k,l
      Call info(message, .true., level=3)
    else
      Call info('  Uniform random seed (proc0): [12, 34, 56, 78]', .true., level=3)
    end if
    Write (message, "(a, '[', 2(i0.1, ',', 1X), i0.1, ']')") '  Saru random seed: ', seed%seed
    Call info(message, .true., level=2)

    Call write_param('Testing DFTB library', devel%test_dftb_library, indent=1, off_level=3)

  End Subroutine write_devel

  Subroutine write_bond_analysis(stats, flow, bond, angle, dihedral, inversion)
    Type(stats_type),      Intent(In   ) :: stats
    Type(flow_type),       Intent(In   ) :: flow
    Type(bonds_type),      Intent(In   ) :: bond
    Type(angles_type),     Intent(In   ) :: angle
    Type(dihedrals_type),  Intent(In   ) :: dihedral
    Type(inversions_type), Intent(In   ) :: inversion

    Character(Len=STR_LEN) :: messages(4)

    Call info('', .true.)

    If (.not. Any([flow%analyse_bond, flow%analyse_ang, flow%analyse_dih, flow%analyse_inv])) Then
      Call info('Bonded statistics collection: OFF', .true., level=3)
    Else
      Call info('Bonded statistics collection requested for:', .true.)

      If (flow%analyse_bond) Then
        Call write_param('Bonds', '', indent=1)
        Call write_param('Collect every (steps)', flow%freq_bond, indent=2)
        Call write_param('Num samples (points)', bond%bin_pdf, indent=2)
        Call write_param('Cutoff', bond%rcut, 'internal_l', indent=2)
      End If

      If (flow%analyse_ang) Then
        Call write_param('Angles', '', indent=1)
        Call write_param('Collect every (steps)', flow%freq_angle, indent=2)
        Call write_param('Num samples (points)', angle%bin_adf, indent=2)
      End If

      If (flow%analyse_dih) Then
        Call write_param('Dihedrals', '', indent=1)
        Call write_param('Collect every (steps)', flow%freq_dihedral, indent=2)
        Call write_param('Num samples (points)', dihedral%bin_adf, indent=2)
      End If

      If (flow%analyse_inv) Then
        Call write_param('Inversions', '', indent=1)
        Call write_param('Collect every (steps)', flow%freq_inversion, indent=2)
        Call write_param('Num samples (points)', inversion%bin_adf, indent=2)
      End If

    End If

    Call write_param('Probability distribution analysis printing', stats%lpana, off_level=3)

  End Subroutine write_bond_analysis

  Subroutine write_structure_analysis(stats, msd_data, rdf, vaf, zdensity, adf, coords, traj, defect, displacement)
    Type(stats_type),      Intent(In   ) :: stats
    Type(msd_type),        Intent(In   ) :: msd_data
    Type(rdf_type),        Intent(InOut) :: rdf
    Type(greenkubo_type),  Intent(In   ) :: vaf
    Type(z_density_type),  Intent(InOut) :: zdensity
    Type(adf_type),        Intent(In   ) :: adf
    Type(coord_type),      Intent(In   ) :: coords
    Type(trajectory_type), Intent(In   ) :: traj
    Type(defects_type),    Intent(In   ) :: defect(:)
    Type(rsd_type),        Intent(In   ) :: displacement

    Character(Len=STR_LEN) :: messages(4)

    Call info('', .true.)
    Call info('Structural statistics:', .true.)


    Call write_param('RDF collection', rdf%l_collect, off_level=3, indent=1)
    If (rdf%l_collect) Then
      Call write_param('Collect every (steps)', rdf%freq, indent=2)
      Call write_param('Bin size', rdf%rbin, 'internal_l', indent=2)

      If (rdf%l_print) Then
        Call info('  -- RDF printing requested', .true., level=3)
      Else If (stats%lpana) Then
        Call info('  -- RDF printing triggered due to a PDA printing request', .true.)
      Else
        Call info('  -- No RDF printing requested', .true.)
      End If

      ! If (rdf%max_rdf == 0) Then
      !   Call info('  -- No RDF pairs specified in FIELD', .true., level=3)
      ! Else
      !   Call info('  -- RDF pairs specified in FIELD', .true.)
      ! End If
    End If

    Call write_param('Z-density profile collection', zdensity%l_collect, off_level=3, indent=1)
    If (zdensity%l_collect) Then
      Call write_param('Collect every (steps)', zdensity%frequency, indent=2)
      Call write_param('Bin size', zdensity%bin_width, 'internal_l', indent=2)
      Call write_param('Z-density printing', zdensity%l_print, indent=2, on_level=3)
    End If

    Call write_param('VAF profile collection', vaf%samp > 0, off_level=3, indent=1)
    If (vaf%samp > 0) Then
      Call write_param('Collect every (steps)', vaf%freq, indent=2)
      Call write_param('Bin size', vaf%binsize, indent=2)
      Call write_param('VAF printing', vaf%l_print, indent=2, on_level=3)
      Call write_param('VAF time-averaging', vaf%l_average, indent=2)
    End If

    Call write_param('MSDTMP profile collection', msd_data%l_msd, off_level=3, indent=1)
    If (msd_data%l_msd) Then
      Call write_param('File start', msd_data%start, indent=2)
      Call write_param('File interval', msd_data%freq, indent=2)
    End If

    Call write_param('Trajectory recording', traj%ltraj, off_level=3, indent=1)
    If (traj%ltraj) Then
      Call write_param('File start', traj%start, indent=2)
      Call write_param('File interval', traj%freq, indent=2)
      Select Case (traj%key)
      Case (TRAJ_KEY_COORD)
        Call write_param('Trajectory file detail', 'COORD', indent=2)
      Case (TRAJ_KEY_COORD_VEL)
        Call write_param('Trajectory file detail', 'COORD, VEL', indent=2)
      Case (TRAJ_KEY_COORD_VEL_FORCE)
        Call write_param('Trajectory file detail', 'COORD, VEL, FORCE', indent=2)
      Case (TRAJ_KEY_COMPRESSED)
        Call write_param('Trajectory file detail', 'COMPRESSED', indent=2)
      End Select
    End If

    Call write_param('Defects analysis', defect(1)%ldef, off_level=3, indent=1)
    If (defect(1)%ldef) Then
      Call write_param('File start', defect(1)%nsdef, indent=2)
      Call write_param('File interval', defect(1)%isdef, indent=2)
      Call write_param('Distance condition', defect(1)%rdef, 'internal_l', indent=2)
      Call write_param('DEFECTS1 file option', defect(2)%ldef, off_level=3)
    End If

    Call write_param('Displacements analysis', displacement%lrsd, off_level=3, indent=1)
    If (displacement%lrsd) Then
      Call write_param('File start', displacement%nsrsd, indent=2)
      Call write_param('File interval', displacement%isrsd, indent=2)
      Call write_param('Distance condition', displacement%rrsd, 'internal_l', indent=2)
    End If

    Call write_param('Coordination analysis', coords%coordon, off_level=3, indent=1)
    If (coords%coordon) Then
      Call info('  -- display to be implemented', .true.)
    End If

    Call write_param('Angular distribution analysis', adf%adfon, off_level=3, indent=1)
    If (adf%adfon) Then
      Call info('  -- display to be implemented', .true.)
    End If

  End Subroutine write_structure_analysis

  Subroutine write_system_parameters(flow, config, stats, thermo, impa, minim, plume, cons, pmf)
    Type(flow_type),          Intent(In   ) :: flow
    Type(configuration_type), Intent(In   ) :: config
    Type(stats_type),         Intent(In   ) :: stats
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(impact_type),        Intent(In   ) :: impa
    Type(minimise_type),      Intent(In   ) :: minim
    Type(plumed_type),        Intent(In   ) :: plume
    Type(constraints_type),   Intent(In   ) :: cons
    Type(pmf_type),           Intent(In   ) :: pmf

    Character(Len=STR_LEN) :: message, messages(4)
    Integer                :: itmp

    Call info('', .true.)
    Call info("System parameters: ", .true.)

    If (.not. thermo%lvar) Then
      Call write_param('Fixed simulation timestep', thermo%tstep, 'internal_t', indent=1)
    Else
      Call write_param('Variable simulation timestep', thermo%tstep, 'internal_t', indent=1)
      Call write_param('Minimum distance Dmin', thermo%mndis, 'internal_l', indent=2)
      Call write_param('Maximum distance Dmin', thermo%mxdis, 'internal_l', indent=2)
      If (thermo%mxstp > zero_plus .and. thermo%mxstp < 1.0e10_wp) Then
        Call write_param('Timestep ceiling max step', thermo%mxstp, 'internal_t', indent=2)
      end If
    End If

    If (flow%run_steps > 0) Call write_param('Run Duration (steps)', flow%run_steps, indent=1)

    If (flow%equil_Steps > 0) Then
      Call write_param('Equilibration period (steps)', flow%equil_steps, indent=1)
      Call write_param('Include equilibration in averages', .not. flow%equilibration, indent=2, off_level=3)

      Call write_param('Temperature regaussing', thermo%freq_tgaus > 0, indent=2, off_level=3)
      If (thermo%freq_tgaus > 0) Call write_param('Temperature regaussing interval (steps)', thermo%freq_tgaus, indent=2)

      Call write_param('Temperature scaling', thermo%freq_tscale > 0, indent=2, off_level=3)
      If (thermo%freq_tgaus > 0) Call write_param('Temperature scaling interval (steps)', thermo%freq_tscale, indent=2)

      Call write_param('Force capping', flow%force_cap, indent=2, off_level=3)
      If (flow%force_cap) Call write_param('Force capping limit', config%fmax, 'k_b.temp/ang', indent=2)

      Call write_param('Minimisation', minim%minimise, indent=2, off_level=3)
      If (minim%minimise) Then
        Select Case (minim%key)
        Case (MIN_FORCE)
          Call write_param('Minimisation criterion', "Force", indent=2)
          message = 'internal_f'
        Case (MIN_ENERGY)
          Call write_param('Minimisation criterion', "Energy", indent=2)
          message = 'internal_e'
        Case (MIN_DISTANCE)
          Call write_param('Minimisation criterion', "Distance", indent=2)
          message = 'internal_l'
        End Select

        Call write_param('Minimisation frequency (steps)', minim%freq, indent=2)
        Call write_param('Minimisation tolerance', minim%tolerance, trim(message), indent=2)
        Call write_param('Minimisation CGM step', minim%step_length, 'internal_l', indent=2)
      End If
    End If

    If (stats%mxstak > 0) Call write_param('Rolling averages length (steps)', stats%mxstak, indent=1)

    If (Any([flow%freq_restart > 0, flow%freq_output > 0, stats%intsta > 0,.not. flow%print_topology])) Then
      Call info('  File output info:', .true.)
    End If

    If (flow%freq_restart > 0) Call write_param('Restart dumping interview (steps)', flow%freq_restart, indent=2)
    If (flow%freq_output > 0) Call write_param('Data printing interval (steps)', flow%freq_output, indent=2)
    If (stats%intsta > 0) Call write_param('Statistics file interval (steps)', stats%intsta, indent=2)

    Call write_param('Extended field topology printing', flow%print_topology, on_level=3, indent=2)

    Call write_param('Computing currents', stats%cur%on, off_level=3, indent=2)
    Call write_param('Computing heat flux', flow%heat_flux, off_level=3, indent=2)
    Call write_param('Write per-particle info', flow%write_per_particle, off_level=3, indent=2)

    Select Case (flow%restart_key)
    Case (RESTART_KEY_SCALE)
      Call info('  Scaled restart requested (starting a new simulation)', .true.)
    Case (RESTART_KEY_NOSCALE)
      Call info('  Unscaled restart requested (starting a new simulation)', .true.)
    Case (RESTART_KEY_OLD)
      Call info('  Restart requested (continuing an old simulation)', .true.)
      Call warning('  Timestep from REVOLD overides specification in CONTROL', .true.)
    Case (RESTART_KEY_CLEAN)
      If (config%levcfg /= 0) Then
        Call info('  Clean start requested, discarding CONFIG velocities', .true.)
      End If
    End Select

    Call write_param('Simulation temperature', thermo%temp, 'K', indent=1)

    Call write_param('Pseudo-thermostat attached to MD cell boundary', thermo%l_stochastic_boundaries, indent=1, off_level=3)
    If (thermo%l_stochastic_boundaries) Then

      Select Case (thermo%key_pseudo)
      Case (PSEUDO_LANGEVIN_DIRECT)
        Call write_param('Thermostat control', 'Langevin + direct temperature scaling', indent=2)
      Case (PSEUDO_LANGEVIN)
        Call write_param('Thermostat control', 'Langevin temperature scaling', indent=2)
      Case (PSEUDO_GAUSSIAN)
        Call write_param('Thermostat control', 'Gaussian temperature scaling', indent=2)
      Case (PSEUDO_DIRECT)
        Call write_param('Thermostat control', 'Direct temperature scaling', indent=2)
      End Select

      Call write_param('Thermostat thickness', thermo%width_pseudo, 'internal_l', indent=2)
      Call write_param('Thermostat temperature', thermo%temp_pseudo, 'K', indent=2)
    End If

    If (thermo%press > 0.0_wp) Then
      Call write_param('Simulation pressure', thermo%press, 'internal_p', indent=1)
    Else If (Any(thermo%stress > 0.0_wp)) Then
      Write (messages(1), '(3A)') '  Simulation pressure (', trim(current_units%pressure), '):'
      Write (messages(2), '(2X, 3(1p g12.5e2))') &
           (convert_units(thermo%stress(itmp), 'internal_p', current_units%pressure), itmp=1, 3)
      Write (messages(3), '(2X, 3(1p g12.5e2))') &
           (convert_units(thermo%stress(itmp), 'internal_p', current_units%pressure), itmp=4, 6)
      Write (messages(4), '(2X, 3(1p g12.5e2))') &
           (convert_units(thermo%stress(itmp), 'internal_p', current_units%pressure), itmp=7, 9)
      Call info(messages, 4, .true.)
    End If

    If (cons%mxcons > 0 .or. pmf%mxpmf > 0) Then
      Call write_param('SHAKE/RATTLE iterations', cons%max_iter_shake, indent=1)
      Call write_param('SHAKE/RATTLE tolerance', cons%tolerance, 'internal_l', indent=1)
    End If

    If (config%l_exp) Then
      Write (message, '(a,9x,3i5)') '  System expansion: ', config%nx, config%ny, config%nz
      Call info(message, .true.)
    End If

    If (config%dvar > 1.0_wp) Then
      Call write_param('Permitted density variance', config%dvar - 1.0_wp, '', '%', indent=1)
    End If

    Call write_param('Impact calculation', impa%active, indent=1, off_level=3)
    If (impa%active) Then
      Call write_param('Particle (index)', impa%imd, indent=2)
      Call write_param('Timestep (steps)', impa%tmd, indent=2)
      Call write_param('Energy', impa%emd, 'ke.V', indent=2)

      Write (message, '(a,3(1p g12.5e2))') '  -- v-r(x,y,z): ', impa%vmx, impa%vmy, impa%vmz
      Call info(message, .true.)

    End If

    Call write_param('PLUMED calculation', plume%l_plumed, indent=1, off_level=3)
    If (plume%l_plumed) Then
      Call write_param('PLUMED input', plume%input, indent=2)
      Call write_param('PLUMED log', plume%logfile, indent=2)
      Call write_param('PLUMED precision', plume%prec, indent=2, level=2)
      Call write_param('PLUMED restart', plume%restart, indent=2)
    End If

  End Subroutine write_system_parameters

  Subroutine write_forcefield(link_cell, neigh, vdws, electro, ewld, mpoles, cshell, met)
    Integer, Dimension(3),    Intent(In   ) :: link_cell
    Type(neighbours_type),    Intent(In   ) :: neigh
    Type(vdw_type),           Intent(In   ) :: vdws
    Type(electrostatic_type), Intent(In   ) :: electro
    Type(ewald_type),         Intent(In   ) :: ewld
    Type(mpole_type),         Intent(In   ) :: mpoles
    Type(core_shell_type),    Intent(In   ) :: cshell
    Type(metal_type),         Intent(In   ) :: met

    Character(Len=STR_LEN) :: message, messages(4)

    Call info('', .true.)
    Call info('Link cell: ', .true.)

    ! ---------------- SUBCELLING ----------------------------------------------

    Call write_param('Subcelling threshold density', neigh%pdplnc, indent=1)
    Call write_param('Cutoff padding', neigh%padding, 'internal_l', indent=1)
    Write (message, '(a,3i6)') "  Final link-cell decomposition (x,y,z): ", link_cell
    Call info(message, .true., level=1)

    ! ---------------- CUTOFF --------------------------------------------------

    Call info('', .true.)
    Call info('Forcefield Parameters: ', .true.)

    Call write_param('Real space cutoff', neigh%cutoff, 'internal_l', indent=1)

    ! ---------------- VDW SETUP -----------------------------------------------

    If (vdws%l_direct) Then
      Call write_param('VdWs', 'Direct calculation', indent=1)
    Else If (vdws%no_vdw) Then
      Call write_param('VdWs', 'Disabled', indent=1)
    Else If (ewld%vdw) Then
      Call write_param('VdWs', 'Ewald', indent=1)
    Else
      Call write_param('VdWs', 'Tabulated', indent=1)
    End If

    If (.not. vdws%no_vdw) Then
      Call write_param('VdW cutoff', vdws%cutoff, 'internal_l', indent=2)

      If (vdws%mixing /= MIX_NULL) Then
        Call info('  -- Vdw cross terms mixing opted (for undefined mixed potentials)', .true.)
        Call info('    mixing is limited to potentials of the same type only', .true.)
        Call info('    mixing restricted to LJ-like potentials (12-6,LJ,WCA,DPD,AMOEBA)', .true.)

        Select Case (vdws%mixing)
        Case (MIX_LORENTZ_BERTHELOT)
          Call write_param('Mixing scheme', 'Lorentz-Berthelot - e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i+s_j)/2', indent=2)

        Case (MIX_FENDER_HALSEY)
          Call write_param('Mixing scheme', 'Fender-Halsey - e_ij=2*e_i*e_j/(e_i+e_j) ; s_ij=(s_i+s_j)/2', indent=2)

        Case (MIX_HOGERVORST)
          Call write_param('Mixing scheme', 'Hogervorst (good hope) - ' &
                    //'e_ij=(e_i*e_j)^(1/2) ; s_ij=(s_i*s_j)^(1/2)', indent=2)

        Case (MIX_HALGREN)
          Call write_param('Mixing scheme', 'Halgren HHG - ' &
                    //'e_ij=4*e_i*e_j/[e_i^(1/2)+e_j^(1/2)]^2 ; s_ij=(s_i^3+s_j^3)/(s_i^2+s_j^2)', indent=2)

        Case (MIX_WALDMAN_HAGLER)
          Call write_param('Mixing scheme', 'Waldman-Hagler - ' &
                    //'e_ij=2*(e_i*e_j)^(1/2)*(s_i*s_j)^3/(s_i^6+s_j^6) ;s_ij=[(s_i^6+s_j^6)/2]^(1/6)', indent=2)

        Case (MIX_TANG_TOENNIES)
          Call write_param('Mixing scheme', 'Tang-Toennies - ' &
                    //' e_ij=[(e_i*s_i^6)*(e_j*s_j^6)] / {[(e_i*s_i^12)^(1/13)+(e_j*s_j^12)^(1/13)]/2}^13', indent=2)
          Call info(Repeat(' ', 43)//'s_ij={[(e_i*s_i^6)*(e_j*s_j^6)]^(1/2) / e_ij}^(1/6)', .true.)

        Case (MIX_FUNCTIONAL)
          Call write_param('Mixing scheme', 'Functional - ' &
                    //'e_ij=3 * (e_i*e_j)^(1/2) * (s_i*s_j)^3 / ' &
                    //'SUM_L=0^2{[(s_i^3+s_j^3)^2 / (4*(s_i*s_j)^L)]^(6/(6-2L))}', indent=2)
          Call info(Repeat(' ', 40)//'s_ij=(1/3) * SUM_L=0^2{[(s_i^3+s_j^3)^2/(4*(s_i*s_j)^L)]^(1/(6-2L))}', .true.)
        End Select
      End If

      Call write_param('VdW force-shifting', vdws%l_force_shift, indent=2, off_level=3)
    End If

    ! ---------------- METALS --------------------------------------------------

    Call write_param('Metal direct', met%l_direct, indent=1, off_level=3)
    Call write_param('Metal sqrtrho', met%l_emb, indent=1, off_level=3)

    ! ELECTROSTATICS

    Call info('', .true.)
    Call info('Electrostatic Parameters: ', .true.)

    ! ---------------- POLARISATION --------------------------------------------


    Call write_param('Multipoles', mpoles%max_mpoles > 0, indent=1, off_level=3)
    If (mpoles%max_mpoles > 0) Then
      Select Case (mpoles%key)
      Case (POLARISATION_CHARMM)
        Call write_param('Polarisation method', 'CHARMM', indent=1)
      Case (POLARISATION_DEFAULT)
        Call write_param('Polarisation method', 'Default', indent=1)
      End Select
    End If

    Call write_param('Core-shell', cshell%mxshl > 0, indent=1, off_level=3)
    If (cshell%mxshl > 0) Then
      Call write_param('Relaxed shell model CGM tolerance', cshell%rlx_tol(1), 'internal_f', indent=2)
      If (cshell%rlx_tol(2) > 0.0_wp) Then
        Call write_param('Relaxed shell model CGM step', cshell%rlx_tol(2), 'internal_l', indent=2)
      End If
    End If

    ! ---------------- ELECTROSTATICS ------------------------------------------

    Select Case (electro%key)
    Case (ELECTROSTATIC_NULL)
      Call write_param('Electrostatics', 'Disabled', indent=1)
    Case (ELECTROSTATIC_EWALD)

      Call write_param('Electrostatics', 'Smooth Particle Mesh Ewald', indent=1)

      If (ewld%precision > 0.0_wp) Call write_param('Ewald sum precision', ewld%precision, indent=2)

      Call write_param('Ewald convergence parameter', ewld%alpha, 'Ang^-1', indent=2)
      Write (message, '(a,3i5)') '  -- Ewald kmax1 kmax2 kmax3   (x2): ', ewld%kspace%k_vec_dim_cont
      Call info(message, .true.)

      If (Any(ewld%kspace%k_vec_dim /= ewld%kspace%k_vec_dim_cont)) Then
        Write (message, '(a,3i5)') '  -- DaFT adjusted kmax values (x2): ', ewld%kspace%k_vec_dim
        Call info(message, .true.)
      End If

      Call write_param('B-Spline interpolation order', ewld%bspline%num_splines, indent=2)

    Case (ELECTROSTATIC_DDDP)
      Call write_param('Electrostatics', 'Distance Dependent Dielectric', indent=1)

    Case (ELECTROSTATIC_COULOMB)
      Call write_param('Electrostatics', 'Coulombic Potential', indent=1)

    Case (ELECTROSTATIC_COULOMB_FORCE_SHIFT)
      Call write_param('Electrostatics', 'Force-Shifted Coulombic Potential', indent=1)

    Case (ELECTROSTATIC_COULOMB_REACTION_FIELD)
      Call write_param('Electrostatics', 'Reaction Field', indent=1)

    End Select

    If (electro%key /= ELECTROSTATIC_NULL) Then
      If (Abs(electro%eps - 1.0_wp) > zero_plus) Call write_param('Relative dielectric constant', electro%eps, indent=2)

      ! Fix with electro merge
      Call write_param('Fennell damping', &
           electro%key /= ELECTROSTATIC_EWALD .and. electro%damping > zero_plus, indent=2, off_level=3)
      If (electro%key /= ELECTROSTATIC_EWALD .and. electro%damping > zero_plus) &
        Call write_param('Damping parameter', electro%damping, 'internal_l^-1', indent=2)

      Call write_param('Extended Coulombic eXclusion', electro%lecx, indent=2)

    End If

  End Subroutine write_forcefield

  Subroutine write_ensemble(thermo, ttm)
    Type(thermostat_type), Intent(In   ) :: thermo
    Type(ttm_type),        Intent(In   ) :: ttm

    Character(Len=STR_LEN)               :: message
    Character(Len=STR_LEN), Dimension(4) :: messages

    Call info('', .true.)
    Call info('Thermostat details:', .true.)

    ensembles:Select Case(thermo%ensemble)
    Case (ENS_NVE)
      Select Case (thermo%key_dpd)
      Case (DPD_NULL)
        Call write_param('Ensemble', 'NVE (Microcanonical)', indent=1)
      Case (DPD_FIRST_ORDER)

        Call write_param('Ensemble', 'NVT dpd (Dissipative Particle Dynamics)', indent=1)
        Call write_param('Ensemble type', "Shardlow's first order splitting (S1)", indent=1)

      Case (DPD_SECOND_ORDER)

        Call write_param('Ensemble', 'NVT dpd (Dissipative Particle Dynamics)', indent=1)
        Call write_param('Ensemble type', "Shardlow's first order splitting (S2)", indent=1)

      End Select

      If (allocated(thermo%gamdpd)) then
        if (thermo%gamdpd(0) > zero_plus) &
             Call write_param('Drag coefficient', thermo%gamdpd(0), 'Da/ps', indent=2)
      end If

    Case (ENS_NVT_EVANS)

      Call write_param('Ensemble', 'NVT Evans (Isokinetic)', indent=1)
      Call info('  Gaussian temperature constraints in use', .true.)

    Case (ENS_NVT_LANGEVIN)

      Call write_param('Ensemble', 'NVT Langevin (Stochastic Dynamics)', indent=1)
      Call write_param('Thermostat friction', thermo%chi, 'internal_t^-1', indent=2)

    Case (ENS_NVT_ANDERSON)

      Call write_param('Ensemble', 'NVT Andersen', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Softness', thermo%soft, indent=2)

    Case (ENS_NVT_BERENDSEN)

      Call write_param('Ensemble', 'NVT Berendsen', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call warning('If you plan to use the Berendsen thermostat, '// &
           'please read https://doi.org/10.1021/acs.jctc.8b00446', .true.)

    Case (ENS_NVT_NOSE_HOOVER)

      Call write_param('Ensemble', 'NVT Nose-Hoover', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)

    Case (ENS_NVT_GENTLE)

      Call write_param('Ensemble', 'NVT gentle stochastic thermostat', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Thermostat friction', thermo%gama, 'internal_t^-1', indent=2)

    Case (ENS_NVT_LANGEVIN_INHOMO)

      Call write_param('Ensemble', 'NVT inhomogeneous Langevin (Stochastic Dynamics)', indent=1)
      Call write_param('e-phonon friction', thermo%chi_ep, 'internal_t^-1', indent=2)
      Call write_param('e-stopping friction', thermo%chi_es, 'internal_t^-1', indent=2)
      Call write_param('e-stopping velocity', thermo%vel_es2, 'internal_l.internal_t^-1', indent=2)
      If (ttm%l_ttm) Then
        Call info('  Applying thermostat with Two-Temperature Model, using local cell electronic temperatures', .true.)
        If (ttm%ttmthvel .and. ttm%ttmthvelz) Then
          Call info('  with only z-components of velocities adjusted using local centre-of-mass corrections for each cell', .true.)
        Else If (ttm%ttmthvel) Then
          Call info('  with velocities adjusted using local centre-of-mass corrections for each cell', .true.)
        Else
          Call info('  without adjusting velocities for local centre-of-mass corrections in each cell', .true.)
        End If
  
      End If

    Case (ENS_NPT_LANGEVIN)

      Call write_param('Ensemble', 'NPT isotropic Langevin (Stochastic Dynamics)', indent=1)
      Call write_param('Thermostat friction', thermo%chi, 'internal_t^-1', indent=2)
      Call write_param('Barostat friction', thermo%tai, 'internal_t^-1', indent=2)

    Case (ENS_NPT_BERENDSEN)

      Call write_param('Ensemble', 'NPT isotropic Berendsen', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Barostat relaxation time', thermo%tau_p, 'internal_t', indent=2)

    Case (ENS_NPT_NOSE_HOOVER)

      Call write_param('Ensemble', 'NPT isotropic Nose-Hoover (Melchionna)', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Barostat relaxation time', thermo%tau_p, 'internal_t', indent=2)

    Case (ENS_NPT_MTK)

      Call write_param('Ensemble', 'NPT isotropic Martyna-Tuckerman-Klein', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Barostat relaxation time', thermo%tau_p, 'internal_t', indent=2)

    Case (ENS_NPT_LANGEVIN_ANISO)

      Call write_param('Ensemble', 'NPT anisotropic Langevin (Stochastic Dynamics)', indent=1)
      Call write_param('Thermostat friction', thermo%chi, 'internal_t^-1', indent=2)
      Call write_param('Barostat friction', thermo%tai, 'internal_t^-1', indent=2)

    Case (ENS_NPT_BERENDSEN_ANISO)

      Call write_param('Ensemble', 'NPT anisotropic Berendsen', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Barostat relaxation time', thermo%tau_p, 'internal_t', indent=2)

    Case (ENS_NPT_NOSE_HOOVER_ANISO)

      Call write_param('Ensemble', 'NPT anisotropic Nose-Hoover (Melchionna)', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Barostat relaxation time', thermo%tau_p, 'internal_t', indent=2)

    Case (ENS_NPT_MTK_ANISO)

      Call write_param('Ensemble', 'NPT anisotropic Martyna-Tuckerman-Klein', indent=1)
      Call write_param('Thermostat relaxation time', thermo%tau_t, 'internal_t', indent=2)
      Call write_param('Barostat relaxation time', thermo%tau_p, 'internal_t', indent=2)

    End Select ensembles

    ! Semi isotropic ensembles

    Select Case (thermo%iso)
    Case (CONSTRAINT_NONE)
      Continue
    Case (CONSTRAINT_SURFACE_AREA)
      Call write_param('Semi-isotropic barostat', 'constant normal pressure (Pn)', indent=1)
      Call write_param('       (N-Pn-A-T)      ', 'constant surface area (A)', indent=1)

    Case (CONSTRAINT_SURFACE_TENSION, CONSTRAINT_SEMI_ORTHORHOMBIC)

      If (thermo%tension > 0.0_wp) Then
        Call write_param('Semi-isotropic barostat', 'constant normal pressure (Pn)', indent=1)
        Call write_param('     (N-Pn-gamma-T)    ', 'constant surface tension (gamma)', indent=1)
        Call write_param('Simulation surface tension', thermo%tension, 'internal_f/internal_l', 'dyn/cm', indent=2)
      Else
        Call write_param('Semi-isotropic barostat', 'orthorhombic MD cell constraints', indent=1)
      End If

      If (thermo%iso == CONSTRAINT_SEMI_ORTHORHOMBIC) Then
        Call write_param('Semi-isotropic barostat', 'semi-orthorhombic MD cell constraints', indent=1)
      End If

    End Select

    If (Any(thermo%iso == [CONSTRAINT_SURFACE_AREA, CONSTRAINT_SURFACE_TENSION])) Then
      Call info('  -- Semi-isotropic ensembles are only correct for infinite', .true.)
      Call info('       interfaces placed perpendicularly to the z axis', .true.)
    End If

  End Subroutine write_ensemble

  Subroutine write_ttm(thermo, ttm)
    Type(thermostat_type), Intent(In   ) :: thermo
    Type(ttm_type),        Intent(In   ) :: ttm

    Character(Len=STR_LEN) :: message, messages(3)

    Call info('', .true.)
    Call info('TTM Parameters: ', .true.)

    Write (messages(1), '(a,3(1x,i8.1))') '  Ionic temperature grid size (x,y,z): ', ttm%ntsys(1:3)
    Write (messages(2), '(a,3(1x,1p g12.5e2))') '  Temperature grid size (x,y,z): ', ttm%delx, ttm%dely, ttm%delz
    Write (messages(3), '(a,1p g12.5e2)') '  Average number of atoms/cell: ', ttm%sysrho * ttm%volume
    Call info(messages, 3, .true.)
    Write (message, '(a,3(1x,i8.1))') '  Electronic temperature grid size (x,y,z): ', ttm%eltsys(1:3)
    Call info(message, .true.)

    If (ttm%ismetal) Then
      Call info('  Electronic subsystem: Metal (thermal conductivity required)', .true.)
    Else
      Call info('  Electronic subsystem: Non-Metal (thermal diffusivity required)', .true.)
    End If

    Call write_param('Dynamic atomic density', ttm%ttmdyndens, indent=1, off_level=3)
    If (.not. ttm%ttmdyndens) Then
      Call write_param('Atomic density', ttm%cellrho, 'internal_l^-3', indent=1)
    Else
      Call write_param('Initial atomic density', ttm%cellrho, 'internal_l^-3', indent=1)
    End If

    Select Case (ttm%cetype)
    Case (TTM_CE_CONST) !Constant
      Call write_param('Electronic specific heat capacity', 'Constant', indent=1)
      Call write_param('Electronic S.H.C. (kB/atom)', ttm%Ce0 / ttm%cellrho, indent=2)

    Case (TTM_CE_CONST_DYN) !Constant
      Call write_param('Electronic specific heat capacity', 'Constant', indent=1)
      Call write_param('Electronic S.H.C. (kB/atom)', ttm%Ce0, indent=2)

    Case (TTM_CE_TANH) !tanh
      Call write_param('Electronic specific heat capacity', 'Hyperbolic tangent', indent=1)
      Call write_param('Constant term A (kB/atom): ', ttm%sh_A / ttm%cellrho, indent=2)
      Call write_param('Temperature term B', ttm%sh_B * 1.0e4_wp, 'K^-1', indent=2)

    Case (TTM_CE_TANH_DYN) !tanh
      Call write_param('Electronic specific heat capacity', 'Hyperbolic tangent', indent=1)
      Call write_param('Constant term A (kB/atom): ', ttm%sh_A, indent=2)
      Call write_param('Temperature term B', ttm%sh_B * 1.0e4_wp, 'K^-1', indent=2)

    Case (TTM_CE_LINEAR, TTM_CE_LINEAR_DYN) !linear
      Call write_param('Electronic specific heat capacity', 'Linear', indent=1)
      Call write_param('Max. electronic S.H.C. (kB/atom): ', ttm%Cemax, indent=2)
      Call write_param('Fermi Temperature', ttm%Tfermi, 'K', indent=2)

    Case (TTM_CE_TABULATED) !tabulated
      Call write_param('Electronic specific heat capacity', 'Tabulated', indent=1)

    End Select

    If (ttm%ismetal) Then
      Select Case (ttm%ketype) ! print thermal conductivity type only for metals
      Case (TTM_KE_INFINITE) ! Infinite
        Call write_param('Electronic thermal conductivity', 'Infinity', indent=1)

      Case (TTM_KE_CONST) ! Constant
        Call write_param('Electronic thermal conductivity', 'Constant', indent=1)
        Call write_param('Electronic T.C.', ttm%Ka0, 'k_b/ps/Ang', 'W/m/K', indent=2)

      Case (TTM_KE_DRUDE) ! Drude
        Call write_param('Electronic thermal conductivity', 'Drude model', indent=1)
        Call write_param('T.C. at system temp.', ttm%Ka0, 'k_b/ps/Ang', 'W/m/K', indent=2)

      Case (TTM_KE_TABULATED) ! Tabulated
        Call write_param('Electronic thermal conductivity', 'Tabulated', indent=1)

      End Select

    Else
      Select Case (ttm%detype) ! print thermal diffusivity type only for non-metals
      Case (TTM_DE_METAL) ! Off (for metals: will probably not get here!)
        Call write_param('Electronic thermal diffusivity', 'Off', indent=1, level=3)
      Case (TTM_DE_CONST) ! Constant
        Call write_param('Electronic thermal diffusivity', 'Constant', indent=1)
        Call write_param('Electronic T.D.', ttm%diff0, 'internal_l^2/internal_t', indent=2)

      Case (TTM_DE_RECIP) ! Recip

        Call write_param('Electronic thermal diffusivity', 'Reciprocal', indent=1)
        Call write_param('Datum Electronic T.D.', ttm%diff0/thermo%temp, 'internal_l^2/internal_t/K', indent=2)
        Call write_param('Fermi Temperature', ttm%Tfermi, 'K', indent=2)

      Case (TTM_DE_TABULATED) ! Tabulated
        Call write_param('Electronic thermal diffusivity', 'Tabulated', indent=1)

      End Select
    End If

    Call write_param('Mininum number of atoms for ionic cells', ttm%amin, indent=1)
    Call write_param('Energy redistribution', ttm%redistribute, indent=1, off_level=3)
    Call write_param('Elec. stopping power', ttm%dedx, 'e.V/ang', indent=1)

    If (ttm%fluence < zero_plus) Then
      Select Case (ttm%sdepoType)
      Case (TTM_SDEPO_GAUSS)
        Call write_param('Spatial energy deposition', 'Gaussian', indent=1)
        Call write_param('Sigma of distribution', ttm%sig, 'ang', indent=2)
        Call write_param('Distribution cutoff', ttm%sigmax * ttm%sig, 'ang', indent=2)

      Case (TTM_SDEPO_FLAT)
        Call write_param('Spatial energy deposition', 'Homogeneous', indent=1)
      End Select

    Else
      Select Case (ttm%sdepoType)
      Case (TTM_SDEPO_FLAT)
        Call write_param('Spatial energy deposition', 'Homogeneous Laser', indent=1)
        Call write_param('Absorbed fluence', ttm%fluence, 'e.V/ang^2', indent=2)
        Call write_param('Penetration depth', ttm%pdepth, 'internal_l', indent=2)

      Case (TTM_SDEPO_EXP)
        Call write_param('Spatial energy deposition', 'Z-exponential decaying Laser', indent=1)
        Call write_param('Absorbed fluence at surface', ttm%fluence, 'e.V/ang^2', indent=2)
        Call write_param('Penetration depth', ttm%pdepth, 'internal_l', indent=2)

      End Select

    End If

    Select Case (ttm%tdepotype)
    Case (TTM_TDEPO_GAUSS)
      Call write_param('Temporal energy deposition', 'Gaussian', indent=1)
      Call write_param('Sigma of distribution', ttm%tdepo, 'ps', indent=2)
      Call write_param('Distribution cutoff', 2.0_wp * ttm%tcdepo * ttm%tdepo, 'ps', indent=2)

    Case (TTM_TDEPO_EXP)
      Call write_param('Temporal energy deposition', 'Decaying exponential', indent=1)
      Call write_param('Tau of distribution', ttm%tdepo, 'ps', indent=2)
      Call write_param('Distribution cutoff', ttm%tcdepo * ttm%tdepo, 'ps', indent=2)

    Case (TTM_TDEPO_DELTA)
      Call write_param('Temporal energy deposition', 'Dirac delta', indent=1)

    Case (TTM_TDEPO_PULSE)
      Call write_param('Temporal energy deposition', 'Square pulse', indent=1)
      Call write_param('Pulse duration', ttm%tdepo, 'ps', indent=2)

    End Select

    Select Case (ttm%gvar)
    Case (TTM_EPVAR_HOMO)
      Call write_param('Variable electron-phonon coupling', 'Homogeneous', indent=1)
      Call info('    (overrides value given for ensemble, required tabulated stopping terms in g.dat file)', .true.)
    Case (TTM_EPVAR_HETERO)
      Call write_param('Variable electron-phonon coupling', 'Heterogeneous', indent=1)
      Call info('    (overrides value given for ensemble, required tabulated stopping terms in g.dat file)', .true.)
    End Select

    Select Case (ttm%bcTypeE)
    Case (TTM_BC_PERIODIC)
      Call write_param('Electronic temperature boundary conditions', 'Periodic', indent=1)

    Case (TTM_BC_DIRICHLET)
      Call write_param('Electronic temperature boundary conditions', 'Dirichlet', indent=1)
      Call info('  -- Boundaries set to system temperature', .true.)

    Case (TTM_BC_NEUMANN)
      Call write_param('Electronic temperature boundary conditions', 'Neumann', indent=1)

    Case (TTM_BC_DIRICHLET_XY)
      Call write_param('Electronic temperature boundary conditions', 'Dirichlet (XY), Neumann (Z)', indent=1)
      Call info('  -- XY boundaries set to system temperature', .true.)

    Case (TTM_BC_ROBIN)
      Call write_param('Electronic temperature boundary conditions', 'Robin', indent=1)
      Call write_param('Temperature leakage', ttm%fluxout, '', '%', indent=2)

    Case (TTM_BC_ROBIN_XY)
      Call write_param('Electronic temperature boundary conditions', 'Robin (XY), Neumann (Z)', indent=1)
      Call write_param('Temperature leakage', ttm%fluxout, '', '%', indent=2)

    End Select

    Call write_param('Electron-ion coupling offset', ttm%ttmoffset, 'ps', indent=1)
    Call write_param('One-way electron-phonon coupling', ttm%oneway, indent=1, off_level=3)

    Call write_param('TTM statistics file', ttm%ttmstats > 0, indent=1, off_level=3)
    If (ttm%ttmstats > 0) Call write_param('TTM statistics file interval (steps)', ttm%ttmstats, indent=2)

    Call write_param('TTM trajectory (temperature profile) file', ttm%ttmtraj > 0, indent=1, off_level=3)
    If (ttm%ttmstats > 0) Call write_param('TTM trajectory file interval (steps)', ttm%ttmtraj, indent=2)

  End Subroutine write_ttm

end module control_output
