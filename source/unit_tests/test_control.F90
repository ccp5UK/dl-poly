Module test_control

  Use asserts, only : assert
  use comms, only : init_comms, comms_type, gsync
  use constants, only : wp
  use filename, only : file_type
  use new_control, only : initialise_control, parse_control_file
  Use errors_warnings, only : error
  Use hash, only : STR_LEN
  use control_parameter_module, only : parameters_hash_table, print_set
  use units, only : initialise_units, set_timestep, convert_units, destroy_units
  implicit none

contains
  Subroutine run_control_tests(comm)
    type(comms_type), intent(inout) :: comm
    type( parameters_hash_table ) :: params
    Real(kind=wp) :: ts
    Real(Kind=wp) :: rtmp
    Character(Len=STR_LEN) :: ctmp
    Real(kind=wp), Dimension(6) :: vtmp
    Integer :: itmp
    Logical :: ltmp
    Integer :: test_unit

    Call gsync(comm)

    open(newunit = test_unit, status = "scratch")

    write(test_unit, '(a,1X,a)') "title", "JUNK"
    write(test_unit, '(a,1X,a)') "simulation_method", "JUNK"
    write(test_unit, '(a,1X,a)') "random_seed", "[ 6.666 6.666 6.666 ]"
    write(test_unit, '(a,1X,a)') "density_variance", "6.666 %"
    write(test_unit, '(a,1X,a)') "data_dump_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "subcell_threshold", "6.666"
    write(test_unit, '(a,1X,a)') "time_run", "6.666 steps"
    write(test_unit, '(a,1X,a)') "time_equilibration", "6.666 steps"
    write(test_unit, '(a,1X,a)') "time_job", "6.666 hr"
    write(test_unit, '(a,1X,a)') "time_close", "6.666 min"
    write(test_unit, '(a,1X,a)') "stats_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "stack_size", "6.666 steps"
    write(test_unit, '(a,1X,a)') "record_equilibration", "on"
    write(test_unit, '(a,1X,a)') "print_per_particle_contrib", "on"
    write(test_unit, '(a,1X,a)') "print_probability_distribution", "on"
    write(test_unit, '(a,1X,a)') "analyse_all", "on"
    write(test_unit, '(a,1X,a)') "analyse_angles", "on"
    write(test_unit, '(a,1X,a)') "analyse_bonds", "on"
    write(test_unit, '(a,1X,a)') "analyse_dihedrals", "on"
    write(test_unit, '(a,1X,a)') "analyse_inversions", "on"
    write(test_unit, '(a,1X,a)') "analyse_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "analyse_frequency_bonds", "6.666 steps"
    write(test_unit, '(a,1X,a)') "analyse_frequency_angles", "6.666 steps"
    write(test_unit, '(a,1X,a)') "analyse_frequency_dihedrals", "6.666 steps"
    write(test_unit, '(a,1X,a)') "analyse_frequency_inversions", "6.666 steps"
    write(test_unit, '(a,1X,a)') "analyse_max_dist", "6.666 ang"
    write(test_unit, '(a,1X,a)') "analyse_num_bins", "66666"
    write(test_unit, '(a,1X,a)') "analyse_num_bins_bonds", "66666"
    write(test_unit, '(a,1X,a)') "analyse_num_bins_angles", "66666"
    write(test_unit, '(a,1X,a)') "analyse_num_bins_dihedrals", "66666"
    write(test_unit, '(a,1X,a)') "analyse_num_bins_inversions", "66666"
    write(test_unit, '(a,1X,a)') "msd_calculate", "on"
    write(test_unit, '(a,1X,a)') "msd_print", "on"
    write(test_unit, '(a,1X,a)') "msd_start", "6.666 steps"
    write(test_unit, '(a,1X,a)') "msd_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "traj_calculate", "on"
    write(test_unit, '(a,1X,a)') "traj_key", "JUNK"
    write(test_unit, '(a,1X,a)') "traj_start", "6.666 steps"
    write(test_unit, '(a,1X,a)') "traj_interval", "6.666 steps"
    write(test_unit, '(a,1X,a)') "defects_calculate", "on"
    write(test_unit, '(a,1X,a)') "defects_start", "6.666 steps"
    write(test_unit, '(a,1X,a)') "defects_interval", "6.666 steps"
    write(test_unit, '(a,1X,a)') "defects_distance", "6.666 ang"
    write(test_unit, '(a,1X,a)') "defects_backup", "on"
    write(test_unit, '(a,1X,a)') "displacements_calculate", "on"
    write(test_unit, '(a,1X,a)') "displacements_start", "6.666 steps"
    write(test_unit, '(a,1X,a)') "displacements_interval", "6.666 steps"
    write(test_unit, '(a,1X,a)') "displacements_distance", "6.666 ang"
    write(test_unit, '(a,1X,a)') "coord_calculate", "on"
    write(test_unit, '(a,1X,a)') "coord_ops", "66666"
    write(test_unit, '(a,1X,a)') "coord_start", "6.666 steps"
    write(test_unit, '(a,1X,a)') "coord_interval", "6.666 steps"
    write(test_unit, '(a,1X,a)') "adf_calculate", "on"
    write(test_unit, '(a,1X,a)') "adf_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "adf_precision", "6.666"
    write(test_unit, '(a,1X,a)') "rdf_calculate", "on"
    write(test_unit, '(a,1X,a)') "rdf_print", "off"
    write(test_unit, '(a,1X,a)') "rdf_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "rdf_binsize", "6.666 ang"
    write(test_unit, '(a,1X,a)') "rdf_error_analysis", "JUNK"
    write(test_unit, '(a,1X,a)') "rdf_error_analysis_blocks", "66666"
    write(test_unit, '(a,1X,a)') "zden_calculate", "on"
    write(test_unit, '(a,1X,a)') "zden_print", "off"
    write(test_unit, '(a,1X,a)') "zden_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "zden_binsize", "6.666 ang"
    write(test_unit, '(a,1X,a)') "vaf_calculate", "on"
    write(test_unit, '(a,1X,a)') "vaf_print", "off"
    write(test_unit, '(a,1X,a)') "vaf_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "vaf_binsize", "66666"
    write(test_unit, '(a,1X,a)') "vaf_averaging", "off"
    write(test_unit, '(a,1X,a)') "currents_calculate", "on"
    write(test_unit, '(a,1X,a)') "print_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "io_units_scheme", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_length", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_time", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_mass", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_charge", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_energy", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_pressure", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_force", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_velocity", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_power", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_surface_tension", "JUNK"
    write(test_unit, '(a,1X,a)') "io_units_emf", "JUNK"
    write(test_unit, '(a,1X,a)') "io_read_method", "JUNK"
    write(test_unit, '(a,1X,a)') "io_read_readers", "66666"
    write(test_unit, '(a,1X,a)') "io_read_batch_size", "66666"
    write(test_unit, '(a,1X,a)') "io_read_buffer_size", "66666"
    write(test_unit, '(a,1X,a)') "io_read_error_check", "on"
    write(test_unit, '(a,1X,a)') "io_read_ascii_revold", "on"
    write(test_unit, '(a,1X,a)') "io_write_method", "JUNK"
    write(test_unit, '(a,1X,a)') "io_write_writers", "66666"
    write(test_unit, '(a,1X,a)') "io_write_batch_size", "66666"
    write(test_unit, '(a,1X,a)') "io_write_buffer_size", "66666"
    write(test_unit, '(a,1X,a)') "io_write_sorted", "off"
    write(test_unit, '(a,1X,a)') "io_write_error_check", "on"
    write(test_unit, '(a,1X,a)') "io_write_netcdf_format", "JUNK"
    write(test_unit, '(a,1X,a)') "io_write_ascii_revive", "on"
    write(test_unit, '(a,1X,a)') "io_file_output", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_config", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_field", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_statis", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_history", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_historf", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_revive", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_revold", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_revcon", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_rdf", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_msd", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_tabbnd", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_tabang", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_tabdih", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_tabinv", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_tabvdw", "JUNK"
    write(test_unit, '(a,1X,a)') "io_file_tabeam", "JUNK"
    write(test_unit, '(a,1X,a)') "output_energy", "on"
    write(test_unit, '(a,1X,a)') "ignore_config_indices", "on"
    write(test_unit, '(a,1X,a)') "print_topology_info", "on"
    write(test_unit, '(a,1X,a)') "print_level", "66666"
    write(test_unit, '(a,1X,a)') "timer_depth", "66666"
    write(test_unit, '(a,1X,a)') "timer_per_mpi", "on"
    write(test_unit, '(a,1X,a)') "timestep", "6.666 internal_t"
    write(test_unit, '(a,1X,a)') "timestep_variable", "on"
    write(test_unit, '(a,1X,a)') "timestep_variable_min_dist", "6.666 ang"
    write(test_unit, '(a,1X,a)') "timestep_variable_max_dist", "6.666 ang"
    write(test_unit, '(a,1X,a)') "timestep_variable_max_delta", "6.666 internal_t"
    write(test_unit, '(a,1X,a)') "ensemble", "JUNK"
    write(test_unit, '(a,1X,a)') "ensemble_method", "JUNK"
    write(test_unit, '(a,1X,a)') "ensemble_thermostat_coupling", "6.666 ps"
    write(test_unit, '(a,1X,a)') "ensemble_dpd_order", "JUNK"
    write(test_unit, '(a,1X,a)') "ensemble_dpd_drag", "6.666 Da/ps"
    write(test_unit, '(a,1X,a)') "ensemble_thermostat_friction", "6.666 ps^-1"
    write(test_unit, '(a,1X,a)') "ensemble_thermostat_softness", "6.666"
    write(test_unit, '(a,1X,a)') "ensemble_barostat_coupling", "6.666 ps"
    write(test_unit, '(a,1X,a)') "ensemble_barostat_friction", "6.666 ps^-1"
    write(test_unit, '(a,1X,a)') "ensemble_semi_isotropic", "JUNK"
    write(test_unit, '(a,1X,a)') "ensemble_semi_orthorhombic", "on"
    write(test_unit, '(a,1X,a)') "ensemble_tension", "6.666 N/m"
    write(test_unit, '(a,1X,a)') "pressure_tensor", "[ 6.666 6.666 6.666 6.666 6.666 6.666 ] katm"
    write(test_unit, '(a,1X,a)') "pressure_hydrostatic", "6.666 katm"
    write(test_unit, '(a,1X,a)') "pressure_perpendicular", "[ 6.666 6.666 6.666 ] katm"
    write(test_unit, '(a,1X,a)') "temperature", "6.666 K"
    write(test_unit, '(a,1X,a)') "pseudo_thermostat_method", "JUNK"
    write(test_unit, '(a,1X,a)') "pseudo_thermostat_width", "6.666 ang"
    write(test_unit, '(a,1X,a)') "pseudo_thermostat_temperature", "6.666 K"
    write(test_unit, '(a,1X,a)') "impact_part_index", "66666"
    write(test_unit, '(a,1X,a)') "impact_time", "6.666 internal_t"
    write(test_unit, '(a,1X,a)') "impact_energy", "6.666 ke.V"
    write(test_unit, '(a,1X,a)') "impact_direction", "[ 6.666 6.666 6.666 ]"
    write(test_unit, '(a,1X,a)') "ttm_calculate", "on"
    write(test_unit, '(a,1X,a)') "ttm_num_ion_cells", "66666"
    write(test_unit, '(a,1X,a)') "ttm_num_elec_cells", "[ 6.666 6.666 6.666 ]"
    write(test_unit, '(a,1X,a)') "ttm_metal", "on"
    write(test_unit, '(a,1X,a)') "ttm_heat_cap_model", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_heat_cap", "6.666 internal_e/internal_m/K"
    write(test_unit, '(a,1X,a)') "ttm_temp_term", "6.666 K^-1"
    write(test_unit, '(a,1X,a)') "ttm_fermi_temp", "6.666 K"
    write(test_unit, '(a,1X,a)') "ttm_elec_cond_model", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_elec_cond", "6.666 W/m/K"
    write(test_unit, '(a,1X,a)') "ttm_diff_model", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_diff", "6.666 m^2/s"
    write(test_unit, '(a,1X,a)') "ttm_dens_model", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_dens", "6.666 ang^-3"
    write(test_unit, '(a,1X,a)') "ttm_min_atoms", "66666"
    write(test_unit, '(a,1X,a)') "ttm_stopping_power", "6.666 e.V/nm"
    write(test_unit, '(a,1X,a)') "ttm_spatial_dist", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_spatial_sigma", "6.666 nm"
    write(test_unit, '(a,1X,a)') "ttm_spatial_cutoff", "6.666 nm"
    write(test_unit, '(a,1X,a)') "ttm_fluence", "6.666 mJ/cm^2"
    write(test_unit, '(a,1X,a)') "ttm_penetration_depth", "6.666 nm"
    write(test_unit, '(a,1X,a)') "ttm_laser_type", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_temporal_dist", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_temporal_duration", "6.666 ps"
    write(test_unit, '(a,1X,a)') "ttm_temporal_cutoff", "6.666 ps"
    write(test_unit, '(a,1X,a)') "ttm_variable_ep", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_boundary_condition", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_boundary_xy", "on"
    write(test_unit, '(a,1X,a)') "ttm_boundary_heat_flux", "on %"
    write(test_unit, '(a,1X,a)') "ttm_time_offset", "6.666 ps"
    write(test_unit, '(a,1X,a)') "ttm_oneway", "on"
    write(test_unit, '(a,1X,a)') "ttm_stats_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "ttm_traj_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "ttm_com_correction", "JUNK"
    write(test_unit, '(a,1X,a)') "ttm_redistribute", "on"
    write(test_unit, '(a,1X,a)') "ttm_e-phonon_friction", "6.666 ps^-1"
    write(test_unit, '(a,1X,a)') "ttm_e-stopping_friction", "6.666 ps^-1"
    write(test_unit, '(a,1X,a)') "ttm_e-stopping_velocity", "6.666 ang/ps"
    write(test_unit, '(a,1X,a)') "rlx_cgm_step", "6.666 ang"
    write(test_unit, '(a,1X,a)') "rlx_tol", "6.666 internal_f"
    write(test_unit, '(a,1X,a)') "shake_max_iter", "66666"
    write(test_unit, '(a,1X,a)') "shake_tolerance", "6.666 ang"
    write(test_unit, '(a,1X,a)') "dftb", "on"
    write(test_unit, '(a,1X,a)') "fixed_com", "off"
    write(test_unit, '(a,1X,a)') "reset_temperature_interval", "6.666 steps"
    write(test_unit, '(a,1X,a)') "regauss_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "rescale_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "equilibration_force_cap", "6.666 k_b.temp/ang"
    write(test_unit, '(a,1X,a)') "minimisation_criterion", "JUNK"
    write(test_unit, '(a,1X,a)') "minimisation_tolerance", "6.666"
    write(test_unit, '(a,1X,a)') "minimisation_step_length", "6.666 ang"
    write(test_unit, '(a,1X,a)') "minimisation_frequency", "6.666 steps"
    write(test_unit, '(a,1X,a)') "initial_minimum_separation", "6.666 internal_l"
    write(test_unit, '(a,1X,a)') "restart", "JUNK"
    write(test_unit, '(a,1X,a)') "nfold", "[ 6.666 6.666 6.666 ]"
    write(test_unit, '(a,1X,a)') "cutoff", "6.666 internal_l"
    write(test_unit, '(a,1X,a)') "padding", "6.666 internal_l"
    write(test_unit, '(a,1X,a)') "coul_damping", "6.666 ang^-1"
    write(test_unit, '(a,1X,a)') "coul_dielectric_constant", "6.666"
    write(test_unit, '(a,1X,a)') "coul_extended_exclusion", "on"
    write(test_unit, '(a,1X,a)') "coul_method", "JUNK"
    write(test_unit, '(a,1X,a)') "coul_precision", "6.666"
    write(test_unit, '(a,1X,a)') "ewald_precision", "6.666"
    write(test_unit, '(a,1X,a)') "ewald_alpha", "6.666 ang^-1"
    write(test_unit, '(a,1X,a)') "ewald_kvec", "[ 6.666 6.666 6.666 ]"
    write(test_unit, '(a,1X,a)') "ewald_kvec_spacing", "6.666 ang^-1"
    write(test_unit, '(a,1X,a)') "ewald_nsplines", "66666"
    write(test_unit, '(a,1X,a)') "polarisation_model", "JUNK"
    write(test_unit, '(a,1X,a)') "polarisation_thole", "6.666"
    write(test_unit, '(a,1X,a)') "metal_direct", "on"
    write(test_unit, '(a,1X,a)') "metal_sqrtrho", "on"
    write(test_unit, '(a,1X,a)') "vdw_method", "JUNK"
    write(test_unit, '(a,1X,a)') "vdw_cutoff", "6.666 internal_l"
    write(test_unit, '(a,1X,a)') "vdw_mix_method", "JUNK"
    write(test_unit, '(a,1X,a)') "vdw_force_shift", "on"
    write(test_unit, '(a,1X,a)') "plumed", "on"
    write(test_unit, '(a,1X,a)') "plumed_input", "JUNK"
    write(test_unit, '(a,1X,a)') "plumed_log", "JUNK"
    write(test_unit, '(a,1X,a)') "plumed_precision", "6.666"
    write(test_unit, '(a,1X,a)') "plumed_restart", "off"
    write(test_unit, '(a,1X,a)') "strict_checks", "off"
    write(test_unit, '(a,1X,a)') "unsafe_comms", "on"

    rewind(test_unit)

    call initialise_units()
    call initialise_control(params)
    Call parse_control_file(test_unit, params, comm)

    close(test_unit)

    call params%retrieve('timestep', ts)
    call set_timestep(ts)

    call params%retrieve("adf_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of adf_calculate failed")
    call params%retrieve("adf_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of adf_frequency failed")
    call params%retrieve("adf_precision", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of adf_precision failed")
    call params%retrieve("analyse_all", ltmp)
    call assert(ltmp, "Accurate retrieval of analyse_all failed")
    call params%retrieve("analyse_angles", ltmp)
    call assert(ltmp, "Accurate retrieval of analyse_angles failed")
    call params%retrieve("analyse_bonds", ltmp)
    call assert(ltmp, "Accurate retrieval of analyse_bonds failed")
    call params%retrieve("analyse_dihedrals", ltmp)
    call assert(ltmp, "Accurate retrieval of analyse_dihedrals failed")
    call params%retrieve("analyse_frequency_angles", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_angles failed")
    call params%retrieve("analyse_frequency_bonds", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_bonds failed")
    call params%retrieve("analyse_frequency_dihedrals", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_dihedrals failed")
    call params%retrieve("analyse_frequency_inversions", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency_inversions failed")
    call params%retrieve("analyse_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_frequency failed")
    call params%retrieve("analyse_inversions", ltmp)
    call assert(ltmp, "Accurate retrieval of analyse_inversions failed")
    call params%retrieve("analyse_max_dist", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of analyse_max_dist failed")
    call params%retrieve("analyse_num_bins_angles", itmp)
    call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_angles failed")
    call params%retrieve("analyse_num_bins_bonds", itmp)
    call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_bonds failed")
    call params%retrieve("analyse_num_bins_dihedrals", itmp)
    call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_dihedrals failed")
    call params%retrieve("analyse_num_bins", itmp)
    call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins failed")
    call params%retrieve("analyse_num_bins_inversions", itmp)
    call assert(itmp, 66666, "Accurate retrieval of analyse_num_bins_inversions failed")
    call params%retrieve("coord_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of coord_calculate failed")
    call params%retrieve("coord_interval", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of coord_interval failed")
    call params%retrieve("coord_ops", itmp)
    call assert(itmp, 66666, "Accurate retrieval of coord_ops failed")
    call params%retrieve("coord_start", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of coord_start failed")
    call params%retrieve("coul_damping", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of coul_damping failed")
    call params%retrieve("coul_dielectric_constant", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of coul_dielectric_constant failed")
    call params%retrieve("coul_extended_exclusion", ltmp)
    call assert(ltmp, "Accurate retrieval of coul_extended_exclusion failed")
    call params%retrieve("coul_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("coul_precision", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of coul_precision failed")
    call params%retrieve("currents_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of currents_calculate failed")
    call params%retrieve("cutoff", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of cutoff failed")
    call params%retrieve("data_dump_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of data_dump_frequency failed")
    call params%retrieve("defects_backup", ltmp)
    call assert(ltmp, "Accurate retrieval of defects_backup failed")
    call params%retrieve("defects_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of defects_calculate failed")
    call params%retrieve("defects_distance", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of defects_distance failed")
    call params%retrieve("defects_interval", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of defects_interval failed")
    call params%retrieve("defects_start", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of defects_start failed")
    call params%retrieve("density_variance", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of density_variance failed")
    call params%retrieve("dftb", ltmp)
    call assert(ltmp, "Accurate retrieval of dftb failed")
    call params%retrieve("displacements_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of displacements_calculate failed")
    call params%retrieve("displacements_distance", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of displacements_distance failed")
    call params%retrieve("displacements_interval", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of displacements_interval failed")
    call params%retrieve("displacements_start", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of displacements_start failed")
    call params%retrieve("ensemble_barostat_coupling", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_barostat_coupling failed")
    call params%retrieve("ensemble_barostat_friction", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_barostat_friction failed")
    call params%retrieve("ensemble", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ensemble_dpd_drag", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_dpd_drag failed")
    call params%retrieve("ensemble_dpd_order", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ensemble_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ensemble_semi_isotropic", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ensemble_semi_orthorhombic", ltmp)
    call assert(ltmp, "Accurate retrieval of ensemble_semi_orthorhombic failed")
    call params%retrieve("ensemble_tension", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_tension failed")
    call params%retrieve("ensemble_thermostat_coupling", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_thermostat_coupling failed")
    call params%retrieve("ensemble_thermostat_friction", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_thermostat_friction failed")
    call params%retrieve("ensemble_thermostat_softness", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ensemble_thermostat_softness failed")
    call params%retrieve("equilibration_force_cap", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of equilibration_force_cap failed")
    call params%retrieve("ewald_alpha", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ewald_alpha failed")
    call params%retrieve("ewald_kvec_spacing", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ewald_kvec_spacing failed")
    call params%retrieve("ewald_kvec", vtmp(1:3))
    call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of ewald_kvec failed")
    call params%retrieve("ewald_nsplines", itmp)
    call assert(itmp, 66666, "Accurate retrieval of ewald_nsplines failed")
    call params%retrieve("ewald_precision", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ewald_precision failed")
    call params%retrieve("fixed_com", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of fixed_com failed")
    call params%retrieve("ignore_config_indices", ltmp)
    call assert(ltmp, "Accurate retrieval of ignore_config_indices failed")
    call params%retrieve("impact_direction", vtmp(1:3))
    call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of impact_direction failed")
    call params%retrieve("impact_energy", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of impact_energy failed")
    call params%retrieve("impact_part_index", itmp)
    call assert(itmp, 66666, "Accurate retrieval of impact_part_index failed")
    call params%retrieve("impact_time", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of impact_time failed")
    call params%retrieve("initial_minimum_separation", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of initial_minimum_separation failed")
    call params%retrieve("io_file_config", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_field", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_historf", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_history", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_msd", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_output", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_rdf", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_revcon", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_revive", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_revold", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_statis", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_tabang", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_tabbnd", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_tabdih", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_tabeam", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_tabinv", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_file_tabvdw", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_read_ascii_revold", ltmp)
    call assert(ltmp, "Accurate retrieval of io_read_ascii_revold failed")
    call params%retrieve("io_read_batch_size", itmp)
    call assert(itmp, 66666, "Accurate retrieval of io_read_batch_size failed")
    call params%retrieve("io_read_buffer_size", itmp)
    call assert(itmp, 66666, "Accurate retrieval of io_read_buffer_size failed")
    call params%retrieve("io_read_error_check", ltmp)
    call assert(ltmp, "Accurate retrieval of io_read_error_check failed")
    call params%retrieve("io_read_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_read_readers", itmp)
    call assert(itmp, 66666, "Accurate retrieval of io_read_readers failed")
    call params%retrieve("io_units_charge", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_emf", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_energy", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_force", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_length", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_mass", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_power", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_pressure", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_scheme", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_surface_tension", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_time", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_units_velocity", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_write_ascii_revive", ltmp)
    call assert(ltmp, "Accurate retrieval of io_write_ascii_revive failed")
    call params%retrieve("io_write_batch_size", itmp)
    call assert(itmp, 66666, "Accurate retrieval of io_write_batch_size failed")
    call params%retrieve("io_write_buffer_size", itmp)
    call assert(itmp, 66666, "Accurate retrieval of io_write_buffer_size failed")
    call params%retrieve("io_write_error_check", ltmp)
    call assert(ltmp, "Accurate retrieval of io_write_error_check failed")
    call params%retrieve("io_write_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_write_netcdf_format", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("io_write_sorted", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of io_write_sorted failed")
    call params%retrieve("io_write_writers", itmp)
    call assert(itmp, 66666, "Accurate retrieval of io_write_writers failed")
    call params%retrieve("metal_direct", ltmp)
    call assert(ltmp, "Accurate retrieval of metal_direct failed")
    call params%retrieve("metal_sqrtrho", ltmp)
    call assert(ltmp, "Accurate retrieval of metal_sqrtrho failed")
    call params%retrieve("minimisation_criterion", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("minimisation_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of minimisation_frequency failed")
    call params%retrieve("minimisation_step_length", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of minimisation_step_length failed")
    call params%retrieve("minimisation_tolerance", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of minimisation_tolerance failed")
    call params%retrieve("msd_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of msd_calculate failed")
    call params%retrieve("msd_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of msd_frequency failed")
    call params%retrieve("msd_print", ltmp)
    call assert(ltmp, "Accurate retrieval of msd_print failed")
    call params%retrieve("msd_start", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of msd_start failed")
    call params%retrieve("nfold", vtmp(1:3))
    call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of nfold failed")
    call params%retrieve("output_energy", ltmp)
    call assert(ltmp, "Accurate retrieval of output_energy failed")
    call params%retrieve("padding", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of padding failed")
    call params%retrieve("plumed_input", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("plumed_log", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("plumed", ltmp)
    call assert(ltmp, "Accurate retrieval of plumed failed")
    call params%retrieve("plumed_precision", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of plumed_precision failed")
    call params%retrieve("plumed_restart", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of plumed_restart failed")
    call params%retrieve("polarisation_model", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("polarisation_thole", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of polarisation_thole failed")
    call params%retrieve("pressure_hydrostatic", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of pressure_hydrostatic failed")
    call params%retrieve("pressure_tensor", vtmp)
    call assert(vtmp,[6.666_wp, 6.666_wp, 6.666_wp, 6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of pressure_tensor failed")
    call params%retrieve("pressure_perpendicular", vtmp(1:3))
    call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of pressure_perpendicular failed")
    call params%retrieve("print_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of print_frequency failed")
    call params%retrieve("print_level", itmp)
    call assert(itmp, 66666, "Accurate retrieval of print_level failed")
    call params%retrieve("print_per_particle_contrib", ltmp)
    call assert(ltmp, "Accurate retrieval of print_per_particle_contrib failed")
    call params%retrieve("print_probability_distribution", ltmp)
    call assert(ltmp, "Accurate retrieval of print_probability_distribution failed")
    call params%retrieve("print_topology_info", ltmp)
    call assert(ltmp, "Accurate retrieval of print_topology_info failed")
    call params%retrieve("pseudo_thermostat_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("pseudo_thermostat_temperature", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of pseudo_thermostat_temperature failed")
    call params%retrieve("pseudo_thermostat_width", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of pseudo_thermostat_width failed")
    call params%retrieve("random_seed", vtmp(1:3))
    call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of random_seed failed")
    call params%retrieve("rdf_binsize", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of rdf_binsize failed")
    call params%retrieve("rdf_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of rdf_calculate failed")
    call params%retrieve("rdf_error_analysis_blocks", itmp)
    call assert(itmp, 66666, "Accurate retrieval of rdf_error_analysis_blocks failed")
    call params%retrieve("rdf_error_analysis", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("rdf_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of rdf_frequency failed")
    call params%retrieve("rdf_print", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of rdf_print failed")
    call params%retrieve("record_equilibration", ltmp)
    call assert(ltmp, "Accurate retrieval of record_equilibration failed")
    call params%retrieve("regauss_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of regauss_frequency failed")
    call params%retrieve("rescale_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of rescale_frequency failed")
    call params%retrieve("reset_temperature_interval", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of reset_temperature_interval failed")
    call params%retrieve("restart", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("rlx_cgm_step", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of rlx_cgm_step failed")
    call params%retrieve("rlx_tol", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of rlx_tol failed")
    call params%retrieve("shake_max_iter", itmp)
    call assert(itmp, 66666, "Accurate retrieval of shake_max_iter failed")
    call params%retrieve("shake_tolerance", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of shake_tolerance failed")
    call params%retrieve("simulation_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("stack_size", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of stack_size failed")
    call params%retrieve("stats_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of stats_frequency failed")
    call params%retrieve("strict_checks", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of strict_checks failed")
    call params%retrieve("subcell_threshold", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of subcell_threshold failed")
    call params%retrieve("temperature", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of temperature failed")
    call params%retrieve("time_close", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of time_close failed")
    call params%retrieve("timer_depth", itmp)
    call assert(itmp, 66666, "Accurate retrieval of timer_depth failed")
    call params%retrieve("time_equilibration", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of time_equilibration failed")
    call params%retrieve("time_job", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of time_job failed")
    call params%retrieve("timer_per_mpi", ltmp)
    call assert(ltmp, "Accurate retrieval of timer_per_mpi failed")
    call params%retrieve("time_run", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of time_run failed")
    call params%retrieve("timestep", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep failed")
    call params%retrieve("timestep_variable", ltmp)
    call assert(ltmp, "Accurate retrieval of timestep_variable failed")
    call params%retrieve("timestep_variable_max_delta", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep_variable_max_delta failed")
    call params%retrieve("timestep_variable_max_dist", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep_variable_max_dist failed")
    call params%retrieve("timestep_variable_min_dist", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of timestep_variable_min_dist failed")
    call params%retrieve("title", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("traj_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of traj_calculate failed")
    call params%retrieve("traj_interval", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of traj_interval failed")
    call params%retrieve("traj_key", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("traj_start", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of traj_start failed")
    call params%retrieve("ttm_boundary_condition", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_boundary_heat_flux", ltmp)
    call assert(ltmp, "Accurate retrieval of ttm_boundary_heat_flux failed")
    call params%retrieve("ttm_boundary_xy", ltmp)
    call assert(ltmp, "Accurate retrieval of ttm_boundary_xy failed")
    call params%retrieve("ttm_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of ttm_calculate failed")
    call params%retrieve("ttm_com_correction", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_dens_model", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_dens", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_dens failed")
    call params%retrieve("ttm_diff_model", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_diff", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_diff failed")
    call params%retrieve("ttm_elec_cond_model", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_elec_cond", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_elec_cond failed")
    call params%retrieve("ttm_e-phonon_friction", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of  failed")
    call params%retrieve("ttm_e-stopping_friction", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of  failed")
    call params%retrieve("ttm_e-stopping_velocity", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of  failed")
    call params%retrieve("ttm_fermi_temp", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_fermi_temp failed")
    call params%retrieve("ttm_fluence", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_fluence failed")
    call params%retrieve("ttm_heat_cap_model", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_heat_cap", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_heat_cap failed")
    call params%retrieve("ttm_laser_type", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_metal", ltmp)
    call assert(ltmp, "Accurate retrieval of ttm_metal failed")
    call params%retrieve("ttm_min_atoms", itmp)
    call assert(itmp, 66666, "Accurate retrieval of ttm_min_atoms failed")
    call params%retrieve("ttm_num_elec_cells", vtmp(1:3))
    call assert(vtmp(1:3), [6.666_wp, 6.666_wp, 6.666_wp], "Accurate retrieval of ttm_num_elec_cells failed")
    call params%retrieve("ttm_num_ion_cells", itmp)
    call assert(itmp, 66666, "Accurate retrieval of ttm_num_ion_cells failed")
    call params%retrieve("ttm_oneway", ltmp)
    call assert(ltmp, "Accurate retrieval of ttm_oneway failed")
    call params%retrieve("ttm_penetration_depth", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_penetration_depth failed")
    call params%retrieve("ttm_redistribute", ltmp)
    call assert(ltmp, "Accurate retrieval of ttm_redistribute failed")
    call params%retrieve("ttm_spatial_cutoff", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_spatial_cutoff failed")
    call params%retrieve("ttm_spatial_dist", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_spatial_sigma", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_spatial_sigma failed")
    call params%retrieve("ttm_stats_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_stats_frequency failed")
    call params%retrieve("ttm_stopping_power", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_stopping_power failed")
    call params%retrieve("ttm_temporal_cutoff", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_temporal_cutoff failed")
    call params%retrieve("ttm_temporal_dist", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("ttm_temporal_duration", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_temporal_duration failed")
    call params%retrieve("ttm_temp_term", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_temp_term failed")
    call params%retrieve("ttm_time_offset", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_time_offset failed")
    call params%retrieve("ttm_traj_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of ttm_traj_frequency failed")
    call params%retrieve("ttm_variable_ep", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("unsafe_comms", ltmp)
    call assert(ltmp, "Accurate retrieval of unsafe_comms failed")
    call params%retrieve("vaf_averaging", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of vaf_averaging failed")
    call params%retrieve("vaf_binsize", itmp)
    call assert(itmp, 66666, "Accurate retrieval of vaf_binsize failed")
    call params%retrieve("vaf_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of vaf_calculate failed")
    call params%retrieve("vaf_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of vaf_frequency failed")
    call params%retrieve("vaf_print", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of vaf_print failed")
    call params%retrieve("vdw_cutoff", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of vdw_cutoff failed")
    call params%retrieve("vdw_force_shift", ltmp)
    call assert(ltmp, "Accurate retrieval of vdw_force_shift failed")
    call params%retrieve("vdw_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("vdw_mix_method", ctmp)
    call assert(ctmp, "JUNK", "Accurate retrieval of  failed")
    call params%retrieve("zden_binsize", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of zden_binsize failed")
    call params%retrieve("zden_calculate", ltmp)
    call assert(ltmp, "Accurate retrieval of zden_calculate failed")
    call params%retrieve("zden_frequency", rtmp)
    call assert(rtmp, 6.666_wp, "Accurate retrieval of zden_frequency failed")
    call params%retrieve("zden_print", ltmp)
    call assert(.not. ltmp, "Accurate retrieval of zden_print failed")

    call destroy_units()

  end Subroutine run_control_tests
end Module test_control
