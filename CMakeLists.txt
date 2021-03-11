cmake_minimum_required(VERSION 3.0.2)
project(DL_POLY Fortran)

set(DLPOLY_VERSION_MAJOR "5")
set(DLPOLY_VERSION_MINOR "0")
set(DLPOLY_VERSION_PATCH "0")
set(DLPOLY_VERSION "${DLPOLY_VERSION_MAJOR}.${DLPOLY_VERSION_MINOR}.${DLPOLY_VERSION_PATCH}")
set(AUTHOR "Ilian Todorov;Bill Smith")
set(AUTHOR_DETAILS "ilian.todorov@stfc.ac.uk")
set(DESCRIPTION "DL_POLY4 classical MD code")

set(src_dir ${CMAKE_SOURCE_DIR}/source)
set(target_name DLPOLY.Z)

#check compilers play fine
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "5.9.9" )
    message(STATUS "Compiler version acceptable? Yes")
  else()
    message(STATUS "Compiler version acceptable? No")
    message(FATAL_ERROR "You need a GNU supported compiler, minimum version 6.0.0")
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules")
include(DLPOLYBuildOptions)
include(DLPOLYBuildFunctions)
include(GNUInstallDirs)
include(flags)
if(COMMAND cmake_policy)
  #  cmake_policy(SET CMP0004 OLD)
endif(COMMAND cmake_policy)

if(WITH_MPI)
  ### if we want MPI check for it and set the internal kitchen
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH}")
  set(LIBS ${LIBS} ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
  string(REGEX REPLACE "^ " "" LIBS "${LIBS}")
  message(STATUS "Build with MPI support!")
else()
  message(STATUS "Build without MPI support!")
endif()

if (BUILDER)
  message(STATUS "${BUILDER} is the builder")
else()
  message(STATUS "default value for builder will be used")
endif()
#if (NOT HOST)
#  message(STATUS "determining hostname...")
#  cmake_host_system_information(RESULT AH QUERY FQDN)
#  set(HOST "${AH}" CACHE STRING )
#endif()
message(STATUS "building on ${HOST}")
if(WITH_HALF_HALO)
  add_definitions(-DHALF_HALO)
  message(STATUS "Build with HALF_HALO support")
endif()

if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  add_definitions(${OpenMP_Fortran_FLAGS})
  add_definitions(-D__OPENMP)
  set(OMP_LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  message(STATUS "Build with OPENMP support")
else()
  message(STATUS "Build without OPENMP support")
endif()

if(WITH_EVB)
  find_package(LAPACK REQUIRED)
  set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  message(STATUS "Build with EVB support")
else()
	message(STATUS "Build without EVB support")
endif()

if(WITH_NETCDF)
  set(NETCDF_F90 "YES")
  find_package(NetCDF REQUIRED)
  include_directories(${NETCDF_INCLUDES_F90})
  set(LIBS ${LIBS} ${NETCDF_LIBRARIES})
  message(STATUS "Build with netcdf support")
else()
  message(STATUS "Build without netcdf support")
endif()

if(WITH_KIM)
  if(INTERNAL_KIM)
    include(ExternalProject)
    enable_language(CXX)
    ExternalProject_Add(int-kim-api
      GIT_REPOSITORY https://github.com/openkim/kim-api
      GIT_TAG v2.1.3
      GIT_SHALLOW true
      GIT_PROGRESS true
      PREFIX int-kim-api
      BINARY_DIR int-kim-api/kim-api-build
      SOURCE_DIR int-kim-api/kim-api
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DBUILD_SHARED_LIBS=off -DCMAKE_BUILD_TYPE=Release
      )

    set(KIM_INC ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/kim-api/mod)
    set(KIM_LIB ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
    set(LIBS ${LIBS} kim-api)
    message(STATUS "Build with OpenKIM support")
    message(STATUS "Building OpenKIM for you")
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(KIM REQUIRED libkim-api>=2.1.3)
    set(KIM_INC ${KIM_LIBDIR}/kim-api/mod)
    set(KIM_LIB ${KIM_LIBRARY_DIRS})
    set(LIBS ${LIBS} ${KIM_LIBRARIES})
    message(STATUS "Build with OpenKIM support")
    message(STATUS "OpenKIM include dirs ${KIM_LIBRARIES}")
  endif()
  link_directories(${KIM_LIB})
else()
  message(STATUS "Build without OpenKIM support")
endif()

if(WITH_PLUMED)

  if (INTERNAL_PLUMED)
    include(ExternalProject)
    enable_language(CXX)
    ExternalProject_Add(int-plumed
      GIT_REPOSITORY https://github.com/plumed/plumed2.git
      GIT_TAG v2.5.2
      GIT_SHALLOW true
      GIT_PROGRESS true
      BUILD_IN_SOURCE 1
      PREFIX int-plumed-api
      CONFIGURE_COMMAND CXX=${CMAKE_CXX_COMPILER} FC=${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_BINARY_DIR}/int-plumed-api/plumed-api/configure --prefix=${CMAKE_INSTALL_PREFIX}
      SOURCE_DIR int-plumed-api/plumed-api
      )
    set(PLUMED_LIB ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
    set(LIBS ${LIBS} plumed)
    message(STATUS "Build with PLUMED support")
    message(STATUS "Building PLUMED for you")
  else()
    find_package(PLUMED 2.1 REQUIRED)
    message(STATUS "Build with PLUMED support")
    set(LIBS ${LIBS} ${PLUMED_LIBRARIES})
  endif()
  link_directories(${PLUMED_LIB})
else()
  message(STATUS "Build without PLUMED support")
endif()

if(UNIX)
  add_definitions("-DUNIX")
elseif(WIN32)
  add_definitions("-DWIN32")
endif()

if(WITH_EXTRATIME)
  add_definitions("-DCHRONO")
  message(STATUS "Build with extra timing support")
  message(STATUS "Use only for development")
else()
  message(STATUS "Build without extra timing support")
endif()

if(WITH_DFTBP)

  if(INTERNAL_DFTBP)
    include(ExternalProject)
    ExternalProject_Add(int-dftbp
      GIT_REPOSITORY https://Alex_Buccheri@bitbucket.org/Alex_Buccheri/dftb_v18_lib.git
      GIT_SHALLOW true
      GIT_TAG origin/lib
      GIT_PROGRESS true
      PREFIX int-dftbp-api
      SOURCE_DIR int-dftbp-api/dftbp-api
      BINARY_DIR int-dftbp-api/build-api
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCMAKE_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DLIB_CHOICE="GCC"
    )

    set(DFTBP_INC ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
    set(DFTBP_LIB ${CMAKE_INSTALL_PREFIX}/lib)
    link_directories(${DFTBP_LIB})
    set(LIBS ${LIBS} dftb)
    message(STATUS "Build with OpenKIM support")
    message(STATUS "Building OpenKIM for you")
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(DFTB REQUIRED dftbplus>=19.2)
    link_directories(${DFTB_LIBRARY_DIRS})
    set(LIBS ${LIBS} ${DFTB_LIBRARIES})
    include_directories(${DFTB_INCLUDE_DIRS})
    message(STATUS "using DFTB+ library ${DFTB_LIBRARIES}")
    message(STATUS "using DFTB+ include ${DFTB_INCLUDE_DIRS}")
  endif()

  if(${SCALAPACK} STREQUAL "GCC")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(SCALAPACK REQUIRED scalapack>=2.0.2)
    link_directories(${SCALAPACK_LIBRARY_DIRS})
    set(LIBS ${LIBS} ${SCALAPACK_LIBRARIES})
    pkg_check_modules(OPENBLAS REQUIRED openblas>=0.3.0)
    link_directories(${OPENBLAS_LIBRARY_DIRS})
    set(LIBS ${LIBS} ${OPENBLAS_LIBRARIES})
  elseif(${SCALAPACK} STREQUAL "MKL")
    link_directories($ENV{MKLROOT}/lib/intel64)
    set(LIBS ${LIBS} "mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64;mkl_intel_lp64;mkl_intel_thread;mkl_core;iomp5")
  endif()

else()
    message(STATUS "Build without DFTB+ support")
endif()

if(WITH_ASSERT)
  message(STATUS "Build with asserts")
else()
    message(STATUS "Asserts excluded")
endif()

######################

###documentation
if (DOCS_DOXYGEN)
  message(STATUS "Build Doxygen API Documentation")
  set(doc_dir doc)
  find_package(Doxygen REQUIRED)
endif()

if (DOCS_FORD)
  message(STATUS "Build FORD API Documentation")
  set(doc_dir doc)
  find_package(Ford REQUIRED)
endif()
######################
##set the output folder for libs and bins
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

add_subdirectory(source)
configure_file(${CMAKE_SOURCE_DIR}/cmake/modulefile.cmake modulefile)

if(DOCS_DOXYGEN)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/Doxyfile.cmake Doxyfile)
  add_custom_target(doxygen
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile)
  install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} docs)")
  install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html/ DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif(DOCS_DOXYGEN)

if(DOCS_FORD)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/dlpoly.md.cmake dlpoly.md)
  add_custom_target(ford
    ${FORD_EXECUTABLE} ${PROJECT_BINARY_DIR}/dlpoly.md)

endif(DOCS_FORD)

if (WITH_COVERAGE)
  include(CodeCoverage)
  set(BUILD_TESTING ON)
  setup_coverage(coverage coverage)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/CodeAnalysis)
  if (MPIEXEC)
    set(dlp "${MPIEXEC} -n ${MPI_NPROCS} ${CMAKE_BINARY_DIR}/bin/${target_name}")
  else()
    set(dlp "${CMAKE_BINARY_DIR}/bin/${target_name}")
  endif()
  single_coverage(runcoverage coverage "${CMAKE_BINARY_DIR}/bin/${target_name}" CodeAnalysis)
else()
  message(STATUS "Build without code coverage testing")
endif()

if (BUILD_TESTING AND NOT WITH_MPI)
  message(STATUS "testing is supported only with MPI")
endif()

if (WITH_COVERAGE AND NOT WITH_MPI)
  message(STATUS "code coverage is supported only with MPI")
endif()

if(WITH_FORCHECK)
  message(STATUS "Build with forcheck support support: make forcheck")
  set(FORCHECK_EXECUTABLE forchk)
  set(FORCHECK_FLAGS "-I ${CMAKE_SOURCE_DIR}/source,${CMAKE_SOURCE_DIR}/source/VV,${CMAKE_SOURCE_DIR}/source/LFV")
  if(NOT WITH_MPI)
    set(FORCHECK_FLAGS "${FORCHECK_FLAGS},${CMAKE_SOURCE_DIR}/source/SERIAL")
  endif()
  set(FORCHECK_FLAGS "${FORCHECK_FLAGS} -l ${CMAKE_BINARY_DIR}/${target_name}.lst")
  set(FORCHECK_FLAGS "${FORCHECK_FLAGS} -allc -ff -f08 -ancmpl -anref -rigorous -anprg -cpp")
  unset(fors)
  unset(defs)
  foreach(t dlpoly dlpolyvv dlpolythermo ${target_name})
    get_target_property(srcs ${t} SOURCES)
    get_target_property(fld ${t} SOURCE_DIR)
    get_target_property(def ${t} COMPILE_DEFINITIONS)
    if (def)
      list(APPEND defs ${def})
    endif()
    foreach(f ${srcs})
      list(APPEND fors ${fld}/${f})
    endforeach()
  endforeach()
  set(mpi $ENV{FCKDIR}/share/forcheck/MPI_3.flb)
  if(defs)
    set(define "-define ${defs}")
  endif()
  add_custom_target(forcheck
    ${FORCHECK_EXECUTABLE} ${FORCHECK_FLAGS} ${fors} ${mpi}
    )
else()
  message(STATUS "Build without forcheck support!")
endif()

if(BUILD_TESTING)
  include(CTest)
  message(STATUS "Build with testing support support: make test")
  find_package(PythonInterp 3.0 REQUIRED)
  macro (do_test arg)
    if (EXISTS ${CMAKE_SOURCE_DIR}/data/${arg}.tar.xz)
      message(STATUS "Using already downloaded ${arg}.tar.xz")
    else()
      message(STATUS "Trying to download ${arg}.tar.xz from https://www.ccp5.ac.uk/sites/www.ccp5.ac.uk/files/DL_POLY/data/${DLPOLY_VERSION_MAJOR}.${DLPOLY_VERSION_MINOR}")
      file(DOWNLOAD https://www.ccp5.ac.uk/sites/www.ccp5.ac.uk/files/DL_POLY/data/${DLPOLY_VERSION_MAJOR}.${DLPOLY_VERSION_MINOR}/${arg}.tar.xz
	  ${CMAKE_SOURCE_DIR}/data/${arg}.tar.xz STATUS DL_TEST)
      list(GET DL_TEST 0 ERRNO)
      if ("${ERRNO}" EQUAL "0")
	message(STATUS "Download successful for ${arg}.tar.xz")
      else()
	list(GET DL_TEST 1 ERRMSG)
	message(FATAL_ERROR "Download failed for ${arg}.tar.xz  with message: ${ERRMSG}")
      endif()
    endif()
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/${arg})
    file(COPY ${CMAKE_SOURCE_DIR}/data/${arg}.tar.xz
      DESTINATION ${CMAKE_BINARY_DIR}/Testing/${arg}/)
    add_test(NAME ${arg}
      COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/utils/beetest -f all.xml -v -k -s all -t ${arg} "${dlp}" WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
    set_tests_properties(${arg} PROPERTIES PASS_REGULAR_EXPRESSION "Status: PASSED")
  endmacro (do_test)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/all.xml.cmake Testing/all.xml)

  set(dlp "${MPIEXEC} ${MPIEXEC_PREFLAGS} -n ${MPI_NPROCS} ${CMAKE_BINARY_DIR}/bin/${target_name} >& err.log")
  foreach(test RANGE 1 9 1)
    do_test("TEST0${test}")
  endforeach()
  foreach(test RANGE 10 25 1)
    do_test("TEST${test}")
  endforeach()
  foreach(test RANGE 28 69 1)
    do_test("TEST${test}")
  endforeach()
  foreach(test RANGE 70 141 1)
    do_test("TEST${test}")
  endforeach()
  foreach(test RANGE 144 163 1)
    do_test("TEST${test}")
  endforeach()
  foreach(test RANGE 169 171 1)
    do_test("TEST${test}")
  endforeach()
  foreach(test RANGE 173 179 1)
    do_test("TEST${test}")
  endforeach()
  if(WITH_DFTBP)
    do_test("TEST180")
  endif()
  foreach(test RANGE 181 182 1)
    do_test("TEST${test}")
  endforeach()
  if(WITH_EVB)
    foreach(test RANGE 191 194 1)
      do_test("TEST${test}")
    endforeach()
  endif()

else()
  message(STATUS "Build without testing target")
endif()
