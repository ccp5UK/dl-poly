include_directories("${CMAKE_BINARY_DIR}/source")

set(dlpoly_SRCS
  build_book.F90
  build_chrm.F90
  build_excl.F90
  build_tplg.F90
  control.F90
  trajectory.F90
  temperature.F90
  bounds.F90
  system.F90
  tersoff.F90
  three_body.F90
  ffield.F90
  tethers.F90
  deport_data.F90
  halo.F90
  defects.F90
  rsds.F90
  two_body.F90
  flow.F90
  filename.F90
  timer.F90
  external_field.F90
  stochastic_boundary.F90
  four_body.F90
  langevin.F90
  coul_spole.F90
  shared_units.F90
  numerics.F90
  impacts.F90
  errors_warnings.F90
  potentials.F90
  rdfs.F90
  kinds.F90
  vdw.F90
  minimise.F90
  dihedrals.F90
  pmf.F90
  particle.F90
  configuration.F90
  neighbours.F90
  comms.F90
  io.F90
  mpole.F90
  core_shell.F90
  rigid_bodies.F90
  dpd.F90
  statistics.F90
  msd.F90
  kinetics.F90
  ewald.F90
  mpoles_container.F90
  ttm.F90
  ttm_utils.F90
  ttm_track.F90
  angles.F90
  bonds.F90
  constraints.F90
  coul_mpole.F90
  bspline.F90
  kspace.F90
  spme.F90
  ewald_spole.F90
#  ewald_mpole.F90
  inversions.F90
  metal.F90
  poisson.F90
  constants.F90
  parse.F90
  development.F90
  domains.F90
  site.F90
  z_density.F90
  greenkubo.F90
  gpfa235.F90
  parallel_fft.F90
  drivers.F90
  thermostat.F90
  analysis.F90
  electrostatic.F90
  kim.F90
  plumed.F90
  netcdf.F90
  nve.F90
  nvt_anderson.F90
  nvt_berendsen.F90
  nvt_ekin.F90
  nvt_gst.F90
  nvt_nose_hoover.F90
  nvt_langevin.F90
  nst_berendsen.F90
  nst_nose_hoover.F90
  nst_langevin.F90
  nst_mtk.F90
  npt_berendsen.F90
  npt_nose_hoover.F90
  npt_langevin.F90
  npt_mtk.F90
  )
if(BUILDER)
  set_property(SOURCE development.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "BUILDER='${BUILDER}'")
endif()
set_property(SOURCE development.F90 APPEND PROPERTY
  COMPILE_DEFINITIONS "HOST='${HOST}'")

if( NOT WITH_MPI)
  set(dlpoly_SERIAL_SRCS ${src_dir}/SERIAL/mpi_api.F90)
endif()

add_library(dlpoly ${dlpoly_SRCS} ${dlpoly_SERIAL_SRCS})
if( NOT WITH_MPI)
  set_target_properties(dlpoly PROPERTIES
    COMPILE_DEFINITIONS "SERIAL"
    INCLUDE_DIR ${src_dir}/SERIAL
    )
endif()

if(WITH_NETCDF)
  set_target_properties(dlpoly PROPERTIES
    COMPILE_DEFINITIONS "NETCDF"
    )
endif()
if(WITH_PLUMED)
  set_target_properties(dlpoly PROPERTIES
    COMPILE_DEFINITIONS "PLUMED"
    )
endif()

if(WITH_KIM)
  set_target_properties(dlpoly PROPERTIES
    COMPILE_DEFINITIONS "KIM"
    )
endif()
## sets the linking
link_directories (${LIBRARY_OUTPUT_PATH})
add_executable(${target_name} dl_poly.F90)
set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${OMP_LINK_FLAGS} ${MIC_LINK_FLAGS}")
target_link_libraries(${target_name} dlpoly ${LIBS})
install(TARGETS ${target_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
INSTALL_SUBLIB(dlpoly)
