include_directories("${CMAKE_BINARY_DIR}/source")
include_directories(VV)

set(dlpoly_SRCS
build_book.F90
build_chrm.F90
build_excl.F90
build_tplg.F90
defects_link_cells.F90
defects.F90
deport_data.F90
dihedrals.F90
external_field.F90
four_body.F90
impact.F90
intra_coul.F90
link_cells.F90
numerics.F90
mm3lrc.F90
origin_config.F90
shared_units.F90
pmf.F90
quaternions.F90
rdfs.F90
control.F90
ffield.F90
trajectory.F90
halo.F90
temperature.F90
rsds.F90
bounds.F90
msd.F90
system.F90
tag_legend.F90
tersoff.F90
tethers.F90
three_body.F90
two_body.F90
vdw.F90
errors_warnings.F90
zbl.F90
zero_k_optimise.F90
ttm_ion_temperature.F90
ttm_system_init.F90
ttm_system_revive.F90
ttm_table_read.F90
ttm_table_scan.F90
ttm_thermal_diffusion.F90
coul_spole.F90
minimise.F90
)

set(dlpoly_CORE_SRCS
kinds.F90
configuration.F90
vnl.F90
comms.F90
io.F90
mpole.F90
core_shell.F90
dpd.F90
statistics.F90
kinetics.F90
ewald.F90
mpoles_container.F90
ttm.F90
ttm_utils.F90
ttm_track.F90
angles.F90
bonds.F90
constraints.F90
coul_mpoles.F90
ewald_spole.F90
ewald_mpoles.F90
inversions.F90
metal.F90
poisson.F90
setup.F90
parse.F90
development.F90
domains.F90
site.F90
z_density.F90
rigid_bodies.F90
langevin.F90
greenkubo.F90
gpfa235.F90
parallel_fft.F90
)
if(BUILDER)
  set_property(SOURCE development.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "BUILDER='${BUILDER}'")
endif()
set_property(SOURCE development.F90 APPEND PROPERTY
  COMPILE_DEFINITIONS "HOST='${HOST}'")

set(dlpoly_THERMO_SRCS
nvt_b0_scl.F90
nvt_b1_scl.F90
nvt_e0_scl.F90
nvt_e1_scl.F90
)

set(netcdf_SRCS netcdf.F90)
set(kim_SRCS kim.F90)
set(plumed_SRCS plumed.F90)

if( NOT WITH_MPI)
  set(dlpoly_SERIAL_SRCS ${src_dir}/SERIAL/mpi_api.F90)
endif()

add_library(dlpoly ${dlpoly_SRCS})

add_library(dlpolycore ${dlpoly_CORE_SRCS} ${dlpoly_SERIAL_SRCS} ${netcdf_SRCS})
if( NOT WITH_MPI)
  set_target_properties(dlpolycore PROPERTIES
     COMPILE_DEFINITIONS "SERIAL"
     INCLUDE_DIR ${src_dir}/SERIAL
     )
endif()
if(WITH_NETCDF)
  set_target_properties(dlpolycore PROPERTIES
    COMPILE_DEFINITIONS "NETCDF"
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
target_link_libraries(dlpolythermo dlpolycore)
target_link_libraries(dlpolyplumed dlpolycore)
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
