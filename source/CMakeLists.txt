include_directories("${CMAKE_BINARY_DIR}/source")
include_directories(VV)

set(dlpoly_SRCS
build_book.F90
build_chrm.F90
build_excl.F90
build_tplg.F90
defects_link_cells.F90
defects.F90
deport_atomic.F90
dihedrals.F90
dpd.F90
errors.F90
external_field.F90
exchange_grid.F90
four_body_forces.F90
impact.F90
init_intra.F90
intra_coul.F90
intra_mcoul.F90
inversions_compute.F90
inversions_forces.F90
inversions_table_read.F90
io_module.F90
langevin_forces.F90
link_cell_pairs.F90
metal_forces.F90
metal_generate.F90
metal_generate_erf.F90
metal_ld_collect_eam.F90
metal_ld_collect_fst.F90
metal_ld_compute.F90
metal_ld_export.F90
metal_ld_set_halo.F90
metal_lrc.F90
metal_table_derivatives.F90
metal_table_read.F90
minimise_relax.F90
mpoles_rotmat_export.F90
mpoles_rotmat_set_halo.F90
msd_write.F90
numerics.F90
mm3lrc.F90
origin_config.F90
pass_shared_units.F90
pmf_coms.F90
pmf_pseudo_bonds.F90
pmf_quench.F90
pmf_tags.F90
pmf_units_set.F90
pmf_vcoms.F90
quaternions_container.F90
rdf_collect.F90
rdf_compute.F90
rdf_excl_collect.F90
rdf_frzn_collect.F90
read_config.F90
read_config_parallel.F90
read_control.F90
read_field.F90
read_history.F90
read_mpoles.F90
refresh_halo_positions.F90
regauss_temperature.F90
relocate_particles.F90
report_topology.F90
rigid_bodies_coms.F90
rigid_bodies_move.F90
rigid_bodies_quench.F90
rigid_bodies_setup.F90
rigid_bodies_split_torque.F90
rigid_bodies_stress.F90
rigid_bodies_tags.F90
rigid_bodies_widths.F90
rsd_write.F90
scale_config.F90
scale_temperature.F90
scan_config.F90
scan_control.F90
scan_control_io.F90
scan_control_pre.F90
scan_field.F90
set_bounds.F90
set_halo_particles.F90
set_temperature.F90
spme_container.F90
coul_spole.F90
statistics_collect.F90
statistics_connect_frames.F90
statistics_connect_set.F90
statistics_connect_spread.F90
statistics_result.F90
system_expand.F90
system_init.F90
system_revive.F90
tag_legend.F90
tersoff_forces.F90
tersoff_generate.F90
tethers_forces.F90
three_body_forces.F90
ttm_ion_temperature.F90
ttm_system_init.F90
ttm_system_revive.F90
ttm_table_read.F90
ttm_table_scan.F90
ttm_thermal_diffusion.F90
trajectory_write.F90
two_body_forces.F90
update_shared_units.F90
usr_collect.F90
usr_compute.F90
vaf_collect.F90
vaf_compute.F90
vaf_write.F90
vdw_direct_fs_generate.F90
vdw_forces.F90
vdw_generate.F90
vdw_lrc.F90
vdw_table_read.F90
vnl_check.F90
vnl_set_check.F90
warning.F90
write_config.F90
xscale.F90
z_density_collect.F90
z_density_compute.F90
zbl.F90
zero_k_optimise.F90
)

set(dlpoly_CORE_SRCS
kinds.F90
comms.F90
setup_module.F90
parse_module.F90
development_module.F90
domains_module.F90
site_module.F90
configuration.F90
vnl_module.F90
ewald.F90
ewald_spole.F90
ewald_mpoles.F90
mpoles_container.F90
mpole.F90
coul_mpoles.F90
poisson_module.F90
vdw_module.F90
metal_module.F90
tersoff_module.F90
three_body_module.F90
four_body_module.F90
rdf_module.F90
z_density_module.F90
core_shell.F90
constraints.F90
pmf_module.F90
rigid_bodies_module.F90
tethers_module.F90
bonds.F90
angles.F90
inversions_module.F90
langevin_module.F90
minimise_module.F90
msd_module.F90
statistics_module.F90
greenkubo_module.F90
kinetic_module.F90
gpfa_module.F90
parallel_fft.F90
ttm_module.F90
ttm_utils.F90
ttm_track_module.F90
)
if(BUILDER)
  set_property(SOURCE development_module.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "BUILDER='${BUILDER}'")
endif()
set_property(SOURCE development_module.F90 APPEND PROPERTY
  COMPILE_DEFINITIONS "HOST='${HOST}'")

set(dlpoly_THERMO_SRCS
nvt_b0_scl.F90
nvt_b1_scl.F90
nvt_e0_scl.F90
nvt_e1_scl.F90
)

set(netcdf_SRCS netcdf_module.F90)
set(kim_SRCS kim_module.F90)
set(plumed_SRCS plumed_module.F90)

if( NOT WITH_MPI)
  set(dlpoly_SERIAL_SRCS ${src_dir}/SERIAL/mpi_api.F90)
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
target_link_libraries(dlpolykim dlpolycore)
target_link_libraries(dlpolyplumed dlpolycore)
target_link_libraries(dlpolythermo dlpolycore)
target_link_libraries(dlpoly dlpolycore dlpolykim dlpolyplumed)
add_executable(${target_name} dl_poly.F90)
set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${OMP_LINK_FLAGS} ${MIC_LINK_FLAGS}")
target_link_libraries(${target_name} dlpoly dlpolycore dlpolythermo dlpolyvv dlpolykim dlpolyplumed ${LIBS})
install(TARGETS ${target_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
INSTALL_SUBLIB(dlpoly)
INSTALL_SUBLIB(dlpolycore)
INSTALL_SUBLIB(dlpolykim)
INSTALL_SUBLIB(dlpolythermo)
INSTALL_SUBLIB(dlpolyplumed)
