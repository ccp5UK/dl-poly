Module test_control

  Use asserts,                  Only: assert
  Use comms,                    Only: comms_type,&
                                      gsync
  Use control_parameters,       Only: parameters_hash_table
  Use kinds,                    Only: STR_LEN,&
                                      wp
  Use new_control,              Only: initialise_control,&
                                      parse_control_file
  Use units,                    Only: destroy_units,&
                                      initialise_units,&
                                      set_timestep

  Implicit None

  Real(Kind=wp), Parameter :: dyn_cm_to_N_m = 1000.0_wp, &
                              katm_to_internal_p = 1000.0_wp / 163.882576_wp

Contains
  Subroutine run_control_tests(comm, passed)
    Type(comms_type), Intent(InOut) :: comm
    Logical,          Intent(  Out) :: passed

    Character(Len=STR_LEN)      :: ctmp
    Integer                     :: itmp, test_unit
    Logical                     :: ltmp
    Real(kind=wp)               :: rtmp, ts
    Real(kind=wp), Allocatable  :: vtmp(:)
    Type(parameters_hash_table) :: params

    passed = .true.

    Call gsync(comm)

    Open (newunit=test_unit, status="scratch")

    Write (test_unit, '(a,1X,a)') "title", "JUNK"
    Write (test_unit, '(a,1X,a)') "simulation_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "random_seed", "[ 6.666 6.666 6.666 ]"
    Write (test_unit, '(a,1X,a)') "density_variance", "6.666 %"
    Write (test_unit, '(a,1X,a)') "data_dump_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "subcell_threshold", "6.666"
    Write (test_unit, '(a,1X,a)') "time_run", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "time_equilibration", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "time_job", "6.666 hr"
    Write (test_unit, '(a,1X,a)') "time_close", "6.666 min"
    Write (test_unit, '(a,1X,a)') "stats_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "stack_size", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "record_equilibration", "on"
    Write (test_unit, '(a,1X,a)') "write_per_particle", "on"
    Write (test_unit, '(a,1X,a)') "print_probability_distribution", "on"
    Write (test_unit, '(a,1X,a)') "analyse_all", "on"
    Write (test_unit, '(a,1X,a)') "analyse_angles", "on"
    Write (test_unit, '(a,1X,a)') "analyse_bonds", "on"
    Write (test_unit, '(a,1X,a)') "analyse_dihedrals", "on"
    Write (test_unit, '(a,1X,a)') "analyse_inversions", "on"
    Write (test_unit, '(a,1X,a)') "analyse_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "analyse_frequency_bonds", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "analyse_frequency_angles", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "analyse_frequency_dihedrals", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "analyse_frequency_inversions", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "analyse_max_dist", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "analyse_num_bins", "66666"
    Write (test_unit, '(a,1X,a)') "analyse_num_bins_bonds", "66666"
    Write (test_unit, '(a,1X,a)') "analyse_num_bins_angles", "66666"
    Write (test_unit, '(a,1X,a)') "analyse_num_bins_dihedrals", "66666"
    Write (test_unit, '(a,1X,a)') "analyse_num_bins_inversions", "66666"
    Write (test_unit, '(a,1X,a)') "msd_calculate", "on"
    Write (test_unit, '(a,1X,a)') "msd_print", "on"
    Write (test_unit, '(a,1X,a)') "msd_start", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "msd_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "traj_calculate", "on"
    Write (test_unit, '(a,1X,a)') "traj_key", "JUNK"
    Write (test_unit, '(a,1X,a)') "traj_start", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "traj_interval", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "defects_calculate", "on"
    Write (test_unit, '(a,1X,a)') "defects_start", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "defects_interval", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "defects_distance", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "defects_backup", "on"
    Write (test_unit, '(a,1X,a)') "displacements_calculate", "on"
    Write (test_unit, '(a,1X,a)') "displacements_start", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "displacements_interval", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "displacements_distance", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "coord_calculate", "on"
    Write (test_unit, '(a,1X,a)') "coord_ops", "66666"
    Write (test_unit, '(a,1X,a)') "coord_start", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "coord_interval", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "adf_calculate", "on"
    Write (test_unit, '(a,1X,a)') "adf_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "adf_precision", "6.666"
    Write (test_unit, '(a,1X,a)') "rdf_calculate", "on"
    Write (test_unit, '(a,1X,a)') "rdf_print", "off"
    Write (test_unit, '(a,1X,a)') "rdf_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "rdf_binsize", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "rdf_error_analysis", "JUNK"
    Write (test_unit, '(a,1X,a)') "rdf_error_analysis_blocks", "66666"
    Write (test_unit, '(a,1X,a)') "zden_calculate", "on"
    Write (test_unit, '(a,1X,a)') "zden_print", "off"
    Write (test_unit, '(a,1X,a)') "zden_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "zden_binsize", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "vaf_calculate", "on"
    Write (test_unit, '(a,1X,a)') "vaf_print", "off"
    Write (test_unit, '(a,1X,a)') "vaf_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "vaf_binsize", "66666"
    Write (test_unit, '(a,1X,a)') "vaf_averaging", "off"
    Write (test_unit, '(a,1X,a)') "currents_calculate", "on"
    Write (test_unit, '(a,1X,a)') "print_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "io_units_scheme", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_length", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_time", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_mass", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_charge", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_energy", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_pressure", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_force", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_velocity", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_power", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_surface_tension", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_units_emf", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_read_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_read_readers", "66666"
    Write (test_unit, '(a,1X,a)') "io_read_batch_size", "66666"
    Write (test_unit, '(a,1X,a)') "io_read_buffer_size", "66666"
    Write (test_unit, '(a,1X,a)') "io_read_error_check", "on"
    Write (test_unit, '(a,1X,a)') "io_read_ascii_revold", "on"
    Write (test_unit, '(a,1X,a)') "io_write_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_write_writers", "66666"
    Write (test_unit, '(a,1X,a)') "io_write_batch_size", "66666"
    Write (test_unit, '(a,1X,a)') "io_write_buffer_size", "66666"
    Write (test_unit, '(a,1X,a)') "io_write_sorted", "off"
    Write (test_unit, '(a,1X,a)') "io_write_error_check", "on"
    Write (test_unit, '(a,1X,a)') "io_write_ascii_revive", "on"
    Write (test_unit, '(a,1X,a)') "io_file_output", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_config", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_field", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_statis", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_history", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_historf", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_revive", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_revold", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_revcon", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_rdf", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_msd", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_tabbnd", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_tabang", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_tabdih", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_tabinv", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_tabvdw", "JUNK"
    Write (test_unit, '(a,1X,a)') "io_file_tabeam", "JUNK"
    Write (test_unit, '(a,1X,a)') "output_energy", "on"
    Write (test_unit, '(a,1X,a)') "ignore_config_indices", "on"
    Write (test_unit, '(a,1X,a)') "print_topology_info", "on"
    Write (test_unit, '(a,1X,a)') "print_level", "66666"
    Write (test_unit, '(a,1X,a)') "timer_depth", "66666"
    Write (test_unit, '(a,1X,a)') "timer_per_mpi", "on"
    Write (test_unit, '(a,1X,a)') "timestep", "6.666 internal_t"
    Write (test_unit, '(a,1X,a)') "timestep_variable", "on"
    Write (test_unit, '(a,1X,a)') "timestep_variable_min_dist", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "timestep_variable_max_dist", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "timestep_variable_max_delta", "6.666 internal_t"
    Write (test_unit, '(a,1X,a)') "ensemble", "JUNK"
    Write (test_unit, '(a,1X,a)') "ensemble_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "ensemble_thermostat_coupling", "6.666 ps"
    Write (test_unit, '(a,1X,a)') "ensemble_dpd_order", "JUNK"
    Write (test_unit, '(a,1X,a)') "ensemble_dpd_drag", "6.666 Da/ps"
    Write (test_unit, '(a,1X,a)') "ensemble_thermostat_friction", "6.666 ps^-1"
    Write (test_unit, '(a,1X,a)') "ensemble_thermostat_softness", "6.666"
    Write (test_unit, '(a,1X,a)') "ensemble_barostat_coupling", "6.666 ps"
    Write (test_unit, '(a,1X,a)') "ensemble_barostat_friction", "6.666 ps^-1"
    Write (test_unit, '(a,1X,a)') "ensemble_semi_isotropic", "JUNK"
    Write (test_unit, '(a,1X,a)') "ensemble_semi_orthorhombic", "on"
    Write (test_unit, '(a,1X,a)') "ensemble_tension", "6.666 N/m"
    Write (test_unit, '(a,1X,a)') "pressure_tensor", "[ 6.666 6.666 6.666 6.666 6.666 6.666 ] katm"
    Write (test_unit, '(a,1X,a)') "pressure_hydrostatic", "6.666 katm"
    Write (test_unit, '(a,1X,a)') "pressure_perpendicular", "[ 6.666 6.666 6.666 ] katm"
    Write (test_unit, '(a,1X,a)') "temperature", "6.666 K"
    Write (test_unit, '(a,1X,a)') "pseudo_thermostat_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "pseudo_thermostat_width", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "pseudo_thermostat_temperature", "6.666 K"
    Write (test_unit, '(a,1X,a)') "impact_part_index", "66666"
    Write (test_unit, '(a,1X,a)') "impact_time", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "impact_energy", "6.666 ke.V"
    Write (test_unit, '(a,1X,a)') "impact_direction", "[ 6.666 6.666 6.666 ]"
    Write (test_unit, '(a,1X,a)') "ttm_calculate", "on"
    Write (test_unit, '(a,1X,a)') "ttm_num_ion_cells", "66666"
    Write (test_unit, '(a,1X,a)') "ttm_num_elec_cells", "[ 6.666 6.666 6.666 ]"
    Write (test_unit, '(a,1X,a)') "ttm_metal", "on"
    Write (test_unit, '(a,1X,a)') "ttm_heat_cap_model", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_heat_cap", "6.666 k_B"
    Write (test_unit, '(a,1X,a)') "ttm_temp_term", "6.666 K^-1"
    Write (test_unit, '(a,1X,a)') "ttm_fermi_temp", "6.666 K"
    Write (test_unit, '(a,1X,a)') "ttm_elec_cond_model", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_elec_cond", "6.666 W/m/K"
    Write (test_unit, '(a,1X,a)') "ttm_diff_model", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_diff", "6.666 m^2/s"
    Write (test_unit, '(a,1X,a)') "ttm_dens_model", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_dens", "6.666 ang^-3"
    Write (test_unit, '(a,1X,a)') "ttm_min_atoms", "66666"
    Write (test_unit, '(a,1X,a)') "ttm_stopping_power", "6.666 e.V/nm"
    Write (test_unit, '(a,1X,a)') "ttm_spatial_dist", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_spatial_sigma", "6.666 nm"
    Write (test_unit, '(a,1X,a)') "ttm_spatial_cutoff", "6.666"
    Write (test_unit, '(a,1X,a)') "ttm_fluence", "6.666 mJ/cm^2"
    Write (test_unit, '(a,1X,a)') "ttm_penetration_depth", "6.666 nm"
    Write (test_unit, '(a,1X,a)') "ttm_laser_type", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_temporal_dist", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_temporal_duration", "6.666 ps"
    Write (test_unit, '(a,1X,a)') "ttm_temporal_cutoff", "6.666"
    Write (test_unit, '(a,1X,a)') "ttm_variable_ep", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_boundary_condition", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_boundary_xy", "on"
    Write (test_unit, '(a,1X,a)') "ttm_boundary_heat_flux", "on %"
    Write (test_unit, '(a,1X,a)') "ttm_time_offset", "6.666 ps"
    Write (test_unit, '(a,1X,a)') "ttm_oneway", "on"
    Write (test_unit, '(a,1X,a)') "ttm_stats_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "ttm_traj_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "ttm_com_correction", "JUNK"
    Write (test_unit, '(a,1X,a)') "ttm_redistribute", "on"
    Write (test_unit, '(a,1X,a)') "ttm_e-phonon_friction", "6.666 ps^-1"
    Write (test_unit, '(a,1X,a)') "ttm_e-stopping_friction", "6.666 ps^-1"
    Write (test_unit, '(a,1X,a)') "ttm_e-stopping_velocity", "6.666 ang/ps"
    Write (test_unit, '(a,1X,a)') "rlx_cgm_step", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "rlx_tol", "6.666 internal_f"
    Write (test_unit, '(a,1X,a)') "shake_max_iter", "66666"
    Write (test_unit, '(a,1X,a)') "shake_tolerance", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "dftb", "on"
    Write (test_unit, '(a,1X,a)') "fixed_com", "off"
    Write (test_unit, '(a,1X,a)') "reset_temperature_interval", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "regauss_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "rescale_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "equilibration_force_cap", "6.666 k_b.temp/ang"
    Write (test_unit, '(a,1X,a)') "minimisation_criterion", "JUNK"
    Write (test_unit, '(a,1X,a)') "minimisation_tolerance", "6.666"
    Write (test_unit, '(a,1X,a)') "minimisation_step_length", "6.666 ang"
    Write (test_unit, '(a,1X,a)') "minimisation_frequency", "6.666 steps"
    Write (test_unit, '(a,1X,a)') "initial_minimum_separation", "6.666 internal_l"
    Write (test_unit, '(a,1X,a)') "restart", "JUNK"
    Write (test_unit, '(a,1X,a)') "nfold", "[ 6.666 6.666 6.666 ]"
    Write (test_unit, '(a,1X,a)') "cutoff", "6.666 internal_l"
    Write (test_unit, '(a,1X,a)') "padding", "6.666 internal_l"
    Write (test_unit, '(a,1X,a)') "coul_damping", "6.666 ang^-1"
    Write (test_unit, '(a,1X,a)') "coul_dielectric_constant", "6.666"
    Write (test_unit, '(a,1X,a)') "coul_extended_exclusion", "on"
    Write (test_unit, '(a,1X,a)') "coul_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "coul_precision", "6.666"
    Write (test_unit, '(a,1X,a)') "ewald_precision", "6.666"
    Write (test_unit, '(a,1X,a)') "ewald_alpha", "6.666 ang^-1"
    Write (test_unit, '(a,1X,a)') "ewald_kvec", "[ 6.666 6.666 6.666 ]"
    Write (test_unit, '(a,1X,a)') "ewald_kvec_spacing", "6.666 ang^-1"
    Write (test_unit, '(a,1X,a)') "ewald_nsplines", "66666"
    Write (test_unit, '(a,1X,a)') "polarisation_model", "JUNK"
    Write (test_unit, '(a,1X,a)') "polarisation_thole", "6.666"
    Write (test_unit, '(a,1X,a)') "metal_direct", "on"
    Write (test_unit, '(a,1X,a)') "metal_sqrtrho", "on"
    Write (test_unit, '(a,1X,a)') "vdw_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "vdw_cutoff", "6.666 internal_l"
    Write (test_unit, '(a,1X,a)') "vdw_mix_method", "JUNK"
    Write (test_unit, '(a,1X,a)') "vdw_force_shift", "on"
    Write (test_unit, '(a,1X,a)') "plumed", "on"
    Write (test_unit, '(a,1X,a)') "plumed_input", "JUNK"
    Write (test_unit, '(a,1X,a)') "plumed_log", "JUNK"
    Write (test_unit, '(a,1X,a)') "plumed_precision", "6.666"
    Write (test_unit, '(a,1X,a)') "plumed_restart", "off"
    Write (test_unit, '(a,1X,a)') "strict_checks", "off"
    Write (test_unit, '(a,1X,a)') "unsafe_comms", "on"

    Rewind (test_unit)

    Call initialise_units()
    Call initialise_control(params)
    Call parse_control_file(test_unit, params, comm)

    Close (test_unit)

    Call params%retrieve('timestep', ts)
    Call set_timestep(ts)

    Call params%retrieve("adf_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of adf_calculate failed", passed_accum = passed)

    Call params%retrieve("adf_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of adf_frequency failed", passed_accum = passed)

    Call params%retrieve("adf_precision", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of adf_precision failed", passed_accum = passed)

    Call params%retrieve("analyse_all", ltmp)
    Call assert(ltmp, "Accurate retrieval of analyse_all failed", passed_accum = passed)

    Call params%retrieve("analyse_angles", ltmp)
    Call assert(ltmp, "Accurate retrieval of analyse_angles failed", passed_accum = passed)

    Call params%retrieve("analyse_bonds", ltmp)
    Call assert(ltmp, "Accurate retrieval of analyse_bonds failed", passed_accum = passed)

    Call params%retrieve("analyse_dihedrals", ltmp)
    Call assert(ltmp, "Accurate retrieval of analyse_dihedrals failed", passed_accum = passed)

    Call params%retrieve("analyse_frequency_angles", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_angles failed", passed_accum = passed)

    Call params%retrieve("analyse_frequency_bonds", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_bonds failed", passed_accum = passed)

    Call params%retrieve("analyse_frequency_dihedrals", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_dihedrals failed", passed_accum = passed)

    Call params%retrieve("analyse_frequency_inversions", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_inversions failed", passed_accum = passed)

    Call params%retrieve("analyse_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency failed", passed_accum = passed)

    Call params%retrieve("analyse_inversions", ltmp)
    Call assert(ltmp, "Accurate retrieval of analyse_inversions failed", passed_accum = passed)

    Call params%retrieve("analyse_max_dist", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_max_dist failed", passed_accum = passed)

    Call params%retrieve("analyse_num_bins_angles", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_angles failed", passed_accum = passed)

    Call params%retrieve("analyse_num_bins_bonds", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_bonds failed", passed_accum = passed)

    Call params%retrieve("analyse_num_bins_dihedrals", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_dihedrals failed", passed_accum = passed)

    Call params%retrieve("analyse_num_bins", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins failed", passed_accum = passed)

    Call params%retrieve("analyse_num_bins_inversions", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_inversions failed", passed_accum = passed)

    Call params%retrieve("coord_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of coord_calculate failed", passed_accum = passed)

    Call params%retrieve("coord_interval", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of coord_interval failed", passed_accum = passed)

    Call params%retrieve("coord_ops", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of coord_ops failed", passed_accum = passed)

    Call params%retrieve("coord_start", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of coord_start failed", passed_accum = passed)

    Call params%retrieve("coul_damping", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of coul_damping failed", passed_accum = passed)

    Call params%retrieve("coul_dielectric_constant", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of coul_dielectric_constant failed", passed_accum = passed)

    Call params%retrieve("coul_extended_exclusion", ltmp)
    Call assert(ltmp, "Accurate retrieval of coul_extended_exclusion failed", passed_accum = passed)

    Call params%retrieve("coul_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of coul_method failed", passed_accum = passed)

    Call params%retrieve("coul_precision", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of coul_precision failed", passed_accum = passed)

    Call params%retrieve("currents_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of currents_calculate failed", passed_accum = passed)

    Call params%retrieve("cutoff", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of cutoff failed", passed_accum = passed)

    Call params%retrieve("data_dump_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of data_dump_frequency failed", passed_accum = passed)

    Call params%retrieve("defects_backup", ltmp)
    Call assert(ltmp, "Accurate retrieval of defects_backup failed", passed_accum = passed)

    Call params%retrieve("defects_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of defects_calculate failed", passed_accum = passed)

    Call params%retrieve("defects_distance", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of defects_distance failed", passed_accum = passed)

    Call params%retrieve("defects_interval", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of defects_interval failed", passed_accum = passed)

    Call params%retrieve("defects_start", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of defects_start failed", passed_accum = passed)

    Call params%retrieve("density_variance", rtmp)
    Call assert(rtmp, 0.06666_wp, "Accurate retrieval of density_variance failed", passed_accum = passed)

    Call params%retrieve("dftb", ltmp)
    Call assert(ltmp, "Accurate retrieval of dftb failed", passed_accum = passed)

    Call params%retrieve("displacements_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of displacements_calculate failed", passed_accum = passed)

    Call params%retrieve("displacements_distance", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of displacements_distance failed", passed_accum = passed)

    Call params%retrieve("displacements_interval", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of displacements_interval failed", passed_accum = passed)

    Call params%retrieve("displacements_start", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of displacements_start failed", passed_accum = passed)

    Call params%retrieve("ensemble_barostat_coupling", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_barostat_coupling failed", passed_accum = passed)

    Call params%retrieve("ensemble_barostat_friction", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_barostat_friction failed", passed_accum = passed)

    Call params%retrieve("ensemble", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ensemble failed", passed_accum = passed)

    Call params%retrieve("ensemble_dpd_drag", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_dpd_drag failed", passed_accum = passed)

    Call params%retrieve("ensemble_dpd_order", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ensemble_dpd_order failed", passed_accum = passed)

    Call params%retrieve("ensemble_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ensemble_method failed", passed_accum = passed)

    Call params%retrieve("ensemble_semi_isotropic", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ensemble_semi_isotropic failed", passed_accum = passed)

    Call params%retrieve("ensemble_semi_orthorhombic", ltmp)
    Call assert(ltmp, "Accurate retrieval of ensemble_semi_orthorhombic failed", passed_accum = passed)

    Call params%retrieve("ensemble_tension", rtmp)
    ! convert N/m to dyn/cm
    Call assert(rtmp, 6.666_wp * dyn_cm_to_N_m, "Accurate retrieval of ensemble_tension failed", passed_accum = passed)

    Call params%retrieve("ensemble_thermostat_coupling", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_thermostat_coupling failed", passed_accum = passed)

    Call params%retrieve("ensemble_thermostat_friction", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_thermostat_friction failed", passed_accum = passed)

    Call params%retrieve("ensemble_thermostat_softness", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_thermostat_softness failed", passed_accum = passed)

    Call params%retrieve("equilibration_force_cap", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of equilibration_force_cap failed", passed_accum = passed)

    Call params%retrieve("ewald_alpha", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ewald_alpha failed", passed_accum = passed)

    Call params%retrieve("ewald_kvec_spacing", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ewald_kvec_spacing failed", passed_accum = passed)

    Call params%retrieve("ewald_kvec", vtmp, 3)
    Call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of ewald_kvec failed", passed_accum = passed)

    Call params%retrieve("ewald_nsplines", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of ewald_nsplines failed", passed_accum = passed)

    Call params%retrieve("ewald_precision", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ewald_precision failed", passed_accum = passed)

    Call params%retrieve("fixed_com", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of fixed_com failed", passed_accum = passed)

    Call params%retrieve("ignore_config_indices", ltmp)
    Call assert(ltmp, "Accurate retrieval of ignore_config_indices failed", passed_accum = passed)

    Call params%retrieve("impact_direction", vtmp, 3)
    Call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of impact_direction failed", passed_accum = passed)

    Call params%retrieve("impact_energy", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of impact_energy failed", passed_accum = passed)

    Call params%retrieve("impact_part_index", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of impact_part_index failed", passed_accum = passed)

    Call params%retrieve("impact_time", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of impact_time failed", passed_accum = passed)

    Call params%retrieve("initial_minimum_separation", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of initial_minimum_separation failed", passed_accum = passed)

    Call params%retrieve("io_file_config", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_config failed", passed_accum = passed)

    Call params%retrieve("io_file_field", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_field failed", passed_accum = passed)

    Call params%retrieve("io_file_historf", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_historyf failed", passed_accum = passed)

    Call params%retrieve("io_file_history", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_history failed", passed_accum = passed)

    Call params%retrieve("io_file_msd", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_msd failed", passed_accum = passed)

    Call params%retrieve("io_file_output", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_output failed", passed_accum = passed)

    Call params%retrieve("io_file_rdf", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_rdf failed", passed_accum = passed)

    Call params%retrieve("io_file_revcon", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_revcon failed", passed_accum = passed)

    Call params%retrieve("io_file_revive", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_revive failed", passed_accum = passed)

    Call params%retrieve("io_file_revold", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_revold failed", passed_accum = passed)

    Call params%retrieve("io_file_statis", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_statis failed", passed_accum = passed)

    Call params%retrieve("io_file_tabang", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_tabang failed", passed_accum = passed)

    Call params%retrieve("io_file_tabbnd", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_tabbnd failed", passed_accum = passed)

    Call params%retrieve("io_file_tabdih", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_tabdih failed", passed_accum = passed)

    Call params%retrieve("io_file_tabeam", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_tabeam failed", passed_accum = passed)

    Call params%retrieve("io_file_tabinv", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_tabinv failed", passed_accum = passed)

    Call params%retrieve("io_file_tabvdw", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of io_file_tabvdw failed", passed_accum = passed)

    Call params%retrieve("io_read_ascii_revold", ltmp)
    Call assert(ltmp, "Accurate retrieval of io_read_ascii_revold failed", passed_accum = passed)

    Call params%retrieve("io_read_batch_size", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of io_read_batch_size failed", passed_accum = passed)

    Call params%retrieve("io_read_buffer_size", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of io_read_buffer_size failed", passed_accum = passed)

    Call params%retrieve("io_read_error_check", ltmp)
    Call assert(ltmp, "Accurate retrieval of io_read_error_check failed", passed_accum = passed)

    Call params%retrieve("io_read_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_read_method failed", passed_accum = passed)

    Call params%retrieve("io_read_readers", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of io_read_readers failed", passed_accum = passed)

    Call params%retrieve("io_units_charge", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_read_readers failed", passed_accum = passed)

    Call params%retrieve("io_units_emf", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_emf failed", passed_accum = passed)

    Call params%retrieve("io_units_energy", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_energy failed", passed_accum = passed)

    Call params%retrieve("io_units_force", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_force failed", passed_accum = passed)

    Call params%retrieve("io_units_length", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_length failed", passed_accum = passed)

    Call params%retrieve("io_units_mass", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_mass failed", passed_accum = passed)

    Call params%retrieve("io_units_power", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_power failed", passed_accum = passed)

    Call params%retrieve("io_units_pressure", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_pressure failed", passed_accum = passed)

    Call params%retrieve("io_units_scheme", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_scheme failed", passed_accum = passed)

    Call params%retrieve("io_units_surface_tension", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_surface_tension failed", passed_accum = passed)

    Call params%retrieve("io_units_time", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_time failed", passed_accum = passed)

    Call params%retrieve("io_units_velocity", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_units_velocity failed", passed_accum = passed)

    Call params%retrieve("io_write_ascii_revive", ltmp)
    Call assert(ltmp, "Accurate retrieval of io_write_ascii_revive failed", passed_accum = passed)

    Call params%retrieve("io_write_batch_size", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of io_write_batch_size failed", passed_accum = passed)

    Call params%retrieve("io_write_buffer_size", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of io_write_buffer_size failed", passed_accum = passed)

    Call params%retrieve("io_write_error_check", ltmp)
    Call assert(ltmp, "Accurate retrieval of io_write_error_check failed", passed_accum = passed)

    Call params%retrieve("io_write_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of io_write_method failed", passed_accum = passed)

    Call params%retrieve("io_write_sorted", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of io_write_sorted failed", passed_accum = passed)

    Call params%retrieve("io_write_writers", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of io_write_writers failed", passed_accum = passed)

    Call params%retrieve("metal_direct", ltmp)
    Call assert(ltmp, "Accurate retrieval of metal_direct failed", passed_accum = passed)

    Call params%retrieve("metal_sqrtrho", ltmp)
    Call assert(ltmp, "Accurate retrieval of metal_sqrtrho failed", passed_accum = passed)

    Call params%retrieve("minimisation_criterion", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of minimisation_criterion failed", passed_accum = passed)

    Call params%retrieve("minimisation_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of minimisation_frequency failed", passed_accum = passed)

    Call params%retrieve("minimisation_step_length", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of minimisation_step_length failed", passed_accum = passed)

    Call params%retrieve("minimisation_tolerance", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of minimisation_tolerance failed", passed_accum = passed)

    Call params%retrieve("msd_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of msd_calculate failed", passed_accum = passed)

    Call params%retrieve("msd_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of msd_frequency failed", passed_accum = passed)

    Call params%retrieve("msd_print", ltmp)
    Call assert(ltmp, "Accurate retrieval of msd_print failed", passed_accum = passed)

    Call params%retrieve("msd_start", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of msd_start failed", passed_accum = passed)

    Call params%retrieve("nfold", vtmp, 3)
    Call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of nfold failed", passed_accum = passed)

    Call params%retrieve("output_energy", ltmp)
    Call assert(ltmp, "Accurate retrieval of output_energy failed", passed_accum = passed)

    Call params%retrieve("padding", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of padding failed", passed_accum = passed)

    Call params%retrieve("plumed_input", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of plumed_input failed", passed_accum = passed)

    Call params%retrieve("plumed_log", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of plumed_log failed", passed_accum = passed)

    Call params%retrieve("plumed", ltmp)
    Call assert(ltmp, "Accurate retrieval of plumed failed", passed_accum = passed)

    Call params%retrieve("plumed_precision", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of plumed_precision failed", passed_accum = passed)

    Call params%retrieve("plumed_restart", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of plumed_restart failed", passed_accum = passed)

    Call params%retrieve("polarisation_model", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of polarisation_model failed", passed_accum = passed)

    Call params%retrieve("polarisation_thole", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of polarisation_thole failed", passed_accum = passed)

    Call params%retrieve("pressure_hydrostatic", rtmp)
    Call assert(rtmp, 6.666_wp * katm_to_internal_p, "Accurate retrieval of pressure_hydrostatic failed", passed_accum = passed)

    Call params%retrieve("pressure_tensor", vtmp, 6)
    Call assert(vtmp, [6.666_wp, 6.666_wp, 6.666_wp, 6.666_wp, 6.666_wp, 6.666_wp]*katm_to_internal_p, &
      "Accurate retrieval of pressure_tensor failed", passed_accum = passed)

    Call params%retrieve("pressure_perpendicular", vtmp, 3)
    Call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp]*katm_to_internal_p, &
      "Accurate retrieval of pressure_perpendicular failed", passed_accum = passed)

    Call params%retrieve("print_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of print_frequency failed", passed_accum = passed)

    Call params%retrieve("print_level", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of print_level failed", passed_accum = passed)

    Call params%retrieve("print_probability_distribution", ltmp)
    Call assert(ltmp, "Accurate retrieval of print_probability_distribution failed", passed_accum = passed)

    Call params%retrieve("print_topology_info", ltmp)
    Call assert(ltmp, "Accurate retrieval of print_topology_info failed", passed_accum = passed)

    Call params%retrieve("pseudo_thermostat_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of pseudo_thermostat_method failed", passed_accum = passed)

    Call params%retrieve("pseudo_thermostat_temperature", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of pseudo_thermostat_temperature failed", passed_accum = passed)

    Call params%retrieve("pseudo_thermostat_width", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of pseudo_thermostat_width failed", passed_accum = passed)

    Call params%retrieve("random_seed", vtmp, 3)
    Call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of random_seed failed", passed_accum = passed)

    Call params%retrieve("rdf_binsize", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of rdf_binsize failed", passed_accum = passed)

    Call params%retrieve("rdf_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of rdf_calculate failed", passed_accum = passed)

    Call params%retrieve("rdf_error_analysis_blocks", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of rdf_error_analysis_blocks failed", passed_accum = passed)

    Call params%retrieve("rdf_error_analysis", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of rdf_error_analysis failed", passed_accum = passed)

    Call params%retrieve("rdf_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of rdf_frequency failed", passed_accum = passed)

    Call params%retrieve("rdf_print", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of rdf_print failed", passed_accum = passed)

    Call params%retrieve("record_equilibration", ltmp)
    Call assert(ltmp, "Accurate retrieval of record_equilibration failed", passed_accum = passed)

    Call params%retrieve("regauss_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of regauss_frequency failed", passed_accum = passed)

    Call params%retrieve("rescale_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of rescale_frequency failed", passed_accum = passed)

    Call params%retrieve("reset_temperature_interval", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of reset_temperature_interval failed", passed_accum = passed)

    Call params%retrieve("restart", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of restart failed", passed_accum = passed)

    Call params%retrieve("rlx_cgm_step", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of rlx_cgm_step failed", passed_accum = passed)

    Call params%retrieve("rlx_tol", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of rlx_tol failed", passed_accum = passed)

    Call params%retrieve("shake_max_iter", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of shake_max_iter failed", passed_accum = passed)

    Call params%retrieve("shake_tolerance", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of shake_tolerance failed", passed_accum = passed)

    Call params%retrieve("simulation_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of simulation_method failed", passed_accum = passed)

    Call params%retrieve("stack_size", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of stack_size failed", passed_accum = passed)

    Call params%retrieve("stats_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of stats_frequency failed", passed_accum = passed)

    Call params%retrieve("strict_checks", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of strict_checks failed", passed_accum = passed)

    Call params%retrieve("subcell_threshold", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of subcell_threshold failed", passed_accum = passed)

    Call params%retrieve("temperature", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of temperature failed", passed_accum = passed)

    Call params%retrieve("time_close", rtmp)
    Call assert(rtmp, 6.666_wp * 60_wp, "Accurate retrieval of time_close failed", passed_accum = passed)

    Call params%retrieve("timer_depth", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of timer_depth failed", passed_accum = passed)

    Call params%retrieve("time_equilibration", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of time_equilibration failed", passed_accum = passed)

    Call params%retrieve("time_job", rtmp)
    Call assert(rtmp, 6.666_wp * 60_wp * 60_wp, "Accurate retrieval of time_job failed", passed_accum = passed)

    Call params%retrieve("timer_per_mpi", ltmp)
    Call assert(ltmp, "Accurate retrieval of timer_per_mpi failed", passed_accum = passed)

    Call params%retrieve("time_run", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of time_run failed", passed_accum = passed)

    Call params%retrieve("timestep", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep failed", passed_accum = passed)

    Call params%retrieve("timestep_variable", ltmp)
    Call assert(ltmp, "Accurate retrieval of timestep_variable failed", passed_accum = passed)

    Call params%retrieve("timestep_variable_max_delta", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep_variable_max_delta failed", passed_accum = passed)

    Call params%retrieve("timestep_variable_max_dist", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep_variable_max_dist failed", passed_accum = passed)

    Call params%retrieve("timestep_variable_min_dist", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep_variable_min_dist failed", passed_accum = passed)

    Call params%retrieve("title", ctmp)
    Call assert(ctmp, "JUNK", "Accurate retrieval of title failed", passed_accum = passed)

    Call params%retrieve("traj_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of traj_calculate failed", passed_accum = passed)

    Call params%retrieve("traj_interval", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of traj_interval failed", passed_accum = passed)

    Call params%retrieve("traj_key", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of traj_key failed", passed_accum = passed)

    Call params%retrieve("traj_start", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of traj_start failed", passed_accum = passed)

    Call params%retrieve("ttm_boundary_condition", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of traj_start failed", passed_accum = passed)

    Call params%retrieve("ttm_boundary_heat_flux", ltmp)
    Call assert(ltmp, "Accurate retrieval of ttm_boundary_heat_flux failed", passed_accum = passed)

    Call params%retrieve("ttm_boundary_xy", ltmp)
    Call assert(ltmp, "Accurate retrieval of ttm_boundary_xy failed", passed_accum = passed)

    Call params%retrieve("ttm_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of ttm_calculate failed", passed_accum = passed)

    Call params%retrieve("ttm_com_correction", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_com_correction failed", passed_accum = passed)

    Call params%retrieve("ttm_dens_model", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_dens_model failed", passed_accum = passed)

    Call params%retrieve("ttm_dens", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_dens failed", passed_accum = passed)

    Call params%retrieve("ttm_diff_model", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_diff_model failed", passed_accum = passed)

    Call params%retrieve("ttm_diff", rtmp)
    Call assert(rtmp, 6.666_wp * (1e-12_wp)/(1e-10_wp * 1e-10_wp), "Accurate retrieval of ttm_diff failed", passed_accum = passed)

    Call params%retrieve("ttm_elec_cond_model", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_diff_model failed", passed_accum = passed)

    Call params%retrieve("ttm_elec_cond", rtmp)
    Call assert(rtmp, 48.281360_wp, "Accurate retrieval of ttm_elec_cond failed", passed_accum = passed)

    Call params%retrieve("ttm_e-phonon_friction", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_e-phonon_friction failed", passed_accum = passed)

    Call params%retrieve("ttm_e-stopping_friction", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_e-stopping_friction failed", passed_accum = passed)

    Call params%retrieve("ttm_e-stopping_velocity", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_e-stopping_velocity failed", passed_accum = passed)

    Call params%retrieve("ttm_fermi_temp", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_fermi_temp failed", passed_accum = passed)

    Call params%retrieve("ttm_fluence", rtmp)
    Call assert(rtmp, 4.1605899490355451_wp, "Accurate retrieval of ttm_fluence failed", passed_accum = passed)

    Call params%retrieve("ttm_heat_cap_model", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_heat_cap_model failed", passed_accum = passed)

    Call params%retrieve("ttm_heat_cap", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_heat_cap failed", passed_accum = passed)

    Call params%retrieve("ttm_laser_type", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_laser_type failed", passed_accum = passed)

    Call params%retrieve("ttm_metal", ltmp)
    Call assert(ltmp, "Accurate retrieval of ttm_metal failed", passed_accum = passed)

    Call params%retrieve("ttm_min_atoms", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of ttm_min_atoms failed", passed_accum = passed)

    Call params%retrieve("ttm_num_elec_cells", vtmp, 3)
    Call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of ttm_num_elec_cells failed", passed_accum = passed)

    Call params%retrieve("ttm_num_ion_cells", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of ttm_num_ion_cells failed", passed_accum = passed)

    Call params%retrieve("ttm_oneway", ltmp)
    Call assert(ltmp, "Accurate retrieval of ttm_oneway failed", passed_accum = passed)

    Call params%retrieve("ttm_penetration_depth", rtmp)
    Call assert(rtmp,66.659999999999997_wp, "Accurate retrieval of ttm_penetration_depth failed", passed_accum = passed)

    Call params%retrieve("ttm_redistribute", ltmp)
    Call assert(ltmp, "Accurate retrieval of ttm_redistribute failed", passed_accum = passed)

    Call params%retrieve("ttm_spatial_cutoff", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_spatial_cutoff failed", passed_accum = passed)

    Call params%retrieve("ttm_spatial_dist", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_spatial_dist failed", passed_accum = passed)

    Call params%retrieve("ttm_spatial_sigma", rtmp)
    Call assert(rtmp, 66.659999999999997_wp, "Accurate retrieval of ttm_spatial_sigma failed", passed_accum = passed)

    Call params%retrieve("ttm_stats_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_stats_frequency failed", passed_accum = passed)

    Call params%retrieve("ttm_stopping_power", rtmp)
    Call assert(rtmp, 0.66660000000000008_wp, "Accurate retrieval of ttm_stopping_power failed", passed_accum = passed)

    Call params%retrieve("ttm_temporal_cutoff", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_temporal_cutoff failed", passed_accum = passed)

    Call params%retrieve("ttm_temporal_dist", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_temporal_dist failed", passed_accum = passed)

    Call params%retrieve("ttm_temporal_duration", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_temporal_duration failed", passed_accum = passed)

    Call params%retrieve("ttm_temp_term", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_temp_term failed", passed_accum = passed)

    Call params%retrieve("ttm_time_offset", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_time_offset failed", passed_accum = passed)

    Call params%retrieve("ttm_traj_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_traj_frequency failed", passed_accum = passed)

    Call params%retrieve("ttm_variable_ep", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of ttm_variable_ep failed", passed_accum = passed)

    Call params%retrieve("unsafe_comms", ltmp)
    Call assert(ltmp, "Accurate retrieval of unsafe_comms failed", passed_accum = passed)

    Call params%retrieve("vaf_averaging", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of vaf_averaging failed", passed_accum = passed)

    Call params%retrieve("vaf_binsize", itmp)
    Call assert(itmp, 66666, "Accurate retrieval of vaf_binsize failed", passed_accum = passed)

    Call params%retrieve("vaf_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of vaf_calculate failed", passed_accum = passed)

    Call params%retrieve("vaf_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of vaf_frequency failed", passed_accum = passed)

    Call params%retrieve("vaf_print", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of vaf_print failed", passed_accum = passed)

    Call params%retrieve("vdw_cutoff", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of vdw_cutoff failed", passed_accum = passed)

    Call params%retrieve("vdw_force_shift", ltmp)
    Call assert(ltmp, "Accurate retrieval of vdw_force_shift failed", passed_accum = passed)

    Call params%retrieve("vdw_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of vdw_method failed", passed_accum = passed)

    Call params%retrieve("vdw_mix_method", ctmp)
    Call assert(ctmp, "junk", "Accurate retrieval of vdw_mix_method failed", passed_accum = passed)

    Call params%retrieve("zden_binsize", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of zden_binsize failed", passed_accum = passed)

    Call params%retrieve("zden_calculate", ltmp)
    Call assert(ltmp, "Accurate retrieval of zden_calculate failed", passed_accum = passed)

    Call params%retrieve("zden_frequency", rtmp)
    Call assert(rtmp, 6.666_wp, "Accurate retrieval of zden_frequency failed", passed_accum = passed)

    Call params%retrieve("zden_print", ltmp)
    Call assert(.not. ltmp, "Accurate retrieval of zden_print failed", passed_accum = passed)


    Call destroy_units()

  End Subroutine run_control_tests
End Module test_control
