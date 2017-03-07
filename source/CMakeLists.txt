include_directories("${CMAKE_BINARY_DIR}/source")
include_directories(VV)
include_directories(LFV)

set(dlpoly_SRCS
angles_compute.f90
angles_forces.f90
angles_table_read.f90
bonds_compute.f90
bonds_forces.f90
bonds_table_read.f90
build_book_intra.f90
build_chrm_intra.f90
build_excl_intra.f90
build_tplg_intra.f90
check_config.f90
compress_book_intra.f90
constraints_pseudo_bonds.f90
constraints_quench.f90
constraints_tags.f90
core_shell_forces.f90
core_shell_kinetic.f90
core_shell_on_top.f90
core_shell_quench.f90
core_shell_relax.f90
coul_cp_forces.f90
coul_cp_mforces.f90
coul_chrm_forces.f90
coul_dddp_forces.f90
coul_dddp_mforces.f90
coul_fscp_forces.f90
coul_fscp_mforces.f90
coul_rfp_forces.f90
coul_rfp_mforces.f90
d_ene_trq_mpoles.f90
defects1_write.f90
defects_link_cells.f90
defects_reference_export.f90
defects_reference_read.f90
defects_reference_read_parallel.f90
defects_reference_set_halo.f90
defects_reference_write.f90
defects_write.f90
deport_atomic_data.f90
dihedrals_14_check.f90
dihedrals_14_vdw.f90
dihedrals_compute.f90
dihedrals_forces.f90
dihedrals_table_read.f90
dpd_thermostat.f90
dpd_v_export.f90
dpd_v_set_halo.f90
error.f90
ewald_excl_forces.f90
ewald_excl_mforces.f90
ewald_excl_mforces_d.f90
ewald_frzn_forces.f90
ewald_frzn_mforces.f90
ewald_real_forces.f90
ewald_real_mforces.f90
ewald_real_mforces_d.f90
ewald_spme_forces.f90
ewald_spme_mforces.f90
ewald_spme_mforces_d.f90
exchange_grid.f90
export_atomic_data.f90
export_atomic_positions.f90
external_field_apply.f90
external_field_correct.f90
four_body_forces.f90
impact.f90
init_intra.f90
intra_coul.f90
intra_mcoul.f90
inversions_compute.f90
inversions_forces.f90
inversions_table_read.f90
io_module.f90
langevin_forces.f90
link_cell_pairs.f90
metal_forces.f90
metal_generate.f90
metal_generate_erf.f90
metal_ld_collect_eam.f90
metal_ld_collect_fst.f90
metal_ld_compute.f90
metal_ld_export.f90
metal_ld_set_halo.f90
metal_lrc.f90
metal_table_derivatives.f90
metal_table_read.f90
minimise_relax.f90
mpoles_container.f90
mpoles_rotmat_export.f90
mpoles_rotmat_set_halo.f90
msd_write.f90
numeric_container.f90
origin_config.f90
pass_shared_units.f90
pmf_coms.f90
pmf_pseudo_bonds.f90
pmf_quench.f90
pmf_tags.f90
pmf_units_set.f90
pmf_vcoms.f90
quaternions_container.f90
rdf_collect.f90
rdf_compute.f90
rdf_excl_collect.f90
rdf_frzn_collect.f90
read_config.f90
read_config_parallel.f90
read_control.f90
read_field.f90
read_history.f90
read_mpoles.f90
refresh_halo_positions.f90
regauss_temperature.f90
relocate_particles.f90
report_topology.f90
rigid_bodies_coms.f90
rigid_bodies_move.f90
rigid_bodies_quench.f90
rigid_bodies_setup.f90
rigid_bodies_split_torque.f90
rigid_bodies_stress.f90
rigid_bodies_tags.f90
rigid_bodies_widths.f90
rsd_write.f90
scale_config.f90
scale_temperature.f90
scan_config.f90
scan_control.f90
scan_control_io.f90
scan_control_pre.f90
scan_field.f90
set_bounds.f90
set_halo_particles.f90
set_temperature.f90
spme_container.f90
statistics_collect.f90
statistics_connect_frames.f90
statistics_connect_set.f90
statistics_connect_spread.f90
statistics_result.f90
system_expand.f90
system_init.f90
system_revive.f90
tag_legend.f90
tersoff_forces.f90
tersoff_generate.f90
tethers_forces.f90
three_body_forces.f90
ttm_ion_temperature.f90
ttm_system_init.f90
ttm_system_revive.f90
ttm_table_read.f90
ttm_table_scan.f90
ttm_thermal_diffusion.f90
trajectory_write.f90
two_body_forces.f90
update_shared_units.f90
usr_collect.f90
usr_compute.f90
vaf_collect.f90
vaf_compute.f90
vaf_write.f90
vdw_direct_fs_generate.f90
vdw_forces.f90
vdw_generate.f90
vdw_lrc.f90
vdw_table_read.f90
vnl_check.f90
vnl_set_check.f90
warning.f90
write_config.f90
xscale.f90
z_density_collect.f90
z_density_compute.f90
zero_k_optimise.f90
)

set(dlpoly_CORE_SRCS
kinds_f90.f90
comms_module.F90
setup_module.f90
parse_module.f90
development_module.F90
domains_module.f90
site_module.f90
config_module.f90
vnl_module.f90
defects_module.f90
defects1_module.f90
dpd_module.f90
ewald_module.f90
mpoles_module.f90
poisson_module.f90
vdw_module.f90
metal_module.f90
tersoff_module.f90
three_body_module.f90
four_body_module.f90
rdf_module.f90
z_density_module.f90
core_shell_module.f90
constraints_module.f90
pmf_module.f90
rigid_bodies_module.f90
tethers_module.f90
bonds_module.f90
angles_module.f90
dihedrals_module.f90
inversions_module.f90
external_field_module.f90
langevin_module.f90
minimise_module.f90
msd_module.f90
statistics_module.f90
greenkubo_module.f90
kinetic_module.f90
gpfa_module.f90
parallel_fft.f90
ttm_module.f90
ttm_utils.f90
ttm_track_module.f90
)
if(BUILDER)
  set_property(SOURCE development_module.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "BUILDER='${BUILDER}'")
endif()
set_property(SOURCE development_module.F90 APPEND PROPERTY
  COMPILE_DEFINITIONS "HOST='${HOST}'")

set(dlpoly_THERMO_SRCS
nvt_b0_scl.f90
nvt_b1_scl.f90
nvt_e0_scl.f90
nvt_e1_scl.f90
)

set(netcdf_SRCS netcdf_module.F90)
set(kim_SRCS kim_module.F90)
set(plumed_SRCS plumed_module.F90)

if( NOT WITH_MPI)
  set(dlpoly_SERIAL_SRCS ${src_dir}/SERIAL/mpi_module.f90)
endif()

add_library(dlpoly ${dlpoly_SRCS} ${netcdf_SRCS})
if(WITH_NETCDF)
  set_target_properties(dlpoly PROPERTIES
    COMPILE_DEFINITIONS "NETCDF"
    )
endif()

add_library(dlpolycore ${dlpoly_CORE_SRCS} ${dlpoly_SERIAL_SRCS})
if( NOT WITH_MPI)
  set_target_properties(dlpolycore PROPERTIES
     COMPILE_DEFINITIONS "SERIAL"
     INCLUDE_DIR ${src_dir}/SERIAL
     )
endif()
add_library(dlpolythermo  ${dlpoly_THERMO_SRCS})
add_library(dlpolykim  ${kim_SRCS})
if(WITH_KIM)
  set_target_properties(dlpolykim PROPERTIES
    COMPILE_DEFINITIONS "KIM"
    )
endif()
add_library(dlpolyplumed  ${plumed_SRCS})
if(WITH_PLUMED)
  set_target_properties(dlpolyplumed PROPERTIES
    COMPILE_DEFINITIONS "PLUMED"
    )
endif()
## sets the linking
link_directories (${LIBRARY_OUTPUT_PATH})
add_subdirectory(VV)
add_subdirectory(LFV)
target_link_libraries(dlpolykim dlpolycore)
target_link_libraries(dlpolyplumed dlpolycore)
target_link_libraries(dlpolythermo dlpolycore)
target_link_libraries(dlpoly dlpolycore dlpolykim dlpolyplumed)
add_executable(${target_name} dl_poly.f90)
set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${OMP_LINK_FLAGS} ${MIC_LINK_FLAGS}")
target_link_libraries(${target_name} dlpoly dlpolycore dlpolythermo dlpolyvv dlpolylfv dlpolykim dlpolyplumed ${LIBS})
install(TARGETS ${target_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
INSTALL_SUBLIB(dlpoly)
INSTALL_SUBLIB(dlpolycore)
INSTALL_SUBLIB(dlpolykim)
INSTALL_SUBLIB(dlpolythermo)
INSTALL_SUBLIB(dlpolyplumed)
