project(DL_POLY Fortran)

set(DLPOLY_VERSION_MAJOR "4")                                                                                                                                   
set(DLPOLY_VERSION_MINOR "07")
set(DLPOLY_VERSION_PATCH "01")
set(DLPOLY_VERSION "${DLPOLY_VERSION_MAJOR}.${DLPOLY_VERSION_MINOR}.${DLPOLY_VERSION_PATCH}")
set(AUTHOR "Ilian Todorov;Bill Smith")
set(AUTHOR_DETAILS "ilian.todorov@stfc.ac.uk")
set(DESCRIPTION "DL_POLY4 classical MD code")

set(src_dir ${CMAKE_SOURCE_DIR}/source) 
set(target_name DLPOLY.Z) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules")
include(DLPOLYBuildOptions)
include(DLPOLYBuildFunctions)
include(GNUInstallDirs)

cmake_minimum_required(VERSION 3.0.2)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0004 OLD)
endif(COMMAND cmake_policy)

if(WITH_MPI)
### if we want MPI check for it and set the internal kitchen
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "MPI_Fortran_INCLUDE_PATH	${MPI_Fortran_INCLUDE_PATH}")
  if(NOT WITH_PHI)
    set(LIBS ${LIBS} ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
  endif()
  message(STATUS "Build with MPI support!")
else()
  message(STATUS "Build without MPI support!")
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

if(WITH_PHI)
   add_definitions("-mmic")
   set(MIC_LINK_FLAGS "-mmic")
   message(STATUS "Build for Xeon Phi Native")
else()
   message(STATUS "Build for CPU host")
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
   find_package(KIM 1.7.1 REQUIRED)
   include_directories(${KIM_MODULE_DIRS})
   message(STATUS "Build with KIM support")
   set(LIBS ${LIBS} ${KIM_LIBRARIES})
else()
   message(STATUS "Build without KIM support")
endif()

if(WITH_EXTRATIME)
   add_definitions("-DCHRONO")
   message(STATUS "Build with extra timing support")
   message(STATUS "Use only for development")
else()
   message(STATUS "Build without extra timing support")
endif()
######################

###documentation
if (DOCS)
   message(STATUS "Build Doxygen API Documentation")
  set(doc_dir doc)
  find_package(Doxygen REQUIRED)
endif()
######################
##set the output folder for libs and bins
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

add_subdirectory(source)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/modulefile.cmake modulefile)

if(DOCS)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/Doxyfile.cmake Doxyfile)
  add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile)
  install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} docs)")
  install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html/ DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif(DOCS)

if(WITH_TESTING)
  include(CTest)
  message(STATUS "Build with testing support support: make test")

  macro (do_test arg)
    execute_process(
      COMMAND wget -nc ftp://ftp.dl.ac.uk/ccp5/DL_POLY/DL_POLY_4.0/DATA/${arg}.tar.gz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data/
      )
    execute_process(
      COMMAND mkdir -p ${arg}
      COMMAND cp ${CMAKE_SOURCE_DIR}/data/${arg}.tar.gz ${arg}/
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/
      )
    add_test(NAME ${arg} 
      COMMAND ${CMAKE_SOURCE_DIR}/utils/beetest -f all.xml -v -k -s all -t ${arg} "${dlp}"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
    set_tests_properties (${arg}
      PROPERTIES PASS_REGULAR_EXPRESSION "Status: PASSED")
  endmacro (do_test)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/all.xml.cmake Testing/all.xml)

  set(dlp "mpirun -n 8 ${CMAKE_BINARY_DIR}/bin/${target_name}")
  foreach(test TEST1 TEST3 TEST5 TEST7)
    do_test(${test})
  endforeach()

  set(dlp "mpirun -n 16 ${CMAKE_BINARY_DIR}/bin/${target_name}")
  foreach(test TEST2 TEST4 TEST6 TEST8)
    do_test(${test})
  endforeach()
else()
  message(STATUS "Build without testing target")
endif()
