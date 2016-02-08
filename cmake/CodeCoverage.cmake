set(COVERAGE_FLAGS "-g -O0 --coverage -fprofile-arcs -ftest-coverage")
foreach(tool gcov lcov genhtml)
  find_program(${tool}_bin ${tool})
  if(NOT ${tool}_bin)
    message(FATAL_ERROR "${tool} not found! That was all for today!")
  else()
    message(STATUS "Using ${tool} from ${${tool}_bin}")
  endif()
endforeach()

if (CMAKE_Fortran_COMPILER_ID)
  if((NOT ${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU"))
    message(FATAL_ERROR "Compiler is not GNU Fortran! No magic today!") 
  else()
    set(CMAKE_Fortran_FLAGS
      "${COVERAGE_FLAGS}"
      CACHE STRING "Fortran compiler flags for coverage builds."
      FORCE )
  endif() 
else()
  message(STATUS "No Fortran compiler used")
endif()

foreach(COMP C CXX)
  if (CMAKE_${COMP}_COMPILER_ID)
    if((NOT CMAKE_${COMP}_COMPILER_ID MATCHES "GNU") AND (NOT CMAKE_${COMP}_COMPILER_ID MATCHES "Clang"))
      message(FATAL_ERROR "Compiler is not GNU/Clang! No magic today!") 
    endif()
    set(CMAKE_${COMP}_FLAGS
      "${COVERAGE_FLAGS}"
      CACHE STRING "Fortran compiler flags for coverage builds."
      FORCE )
  else()
    message(STATUS "No ${COMP} compiler used")
  endif()
endforeach()

function(setup_coverage _targetname _outputname)
	add_custom_target(${_targetname}
    ${lcov_bin} --directory Testing --zerocounters
    COMMAND ${CMAKE_MAKE_PROGRAM} test
    COMMAND ${lcov_bin} --directory . --capture --output-file ${_outputname}.info
    COMMAND ${lcov_bin} --remove ${_outputname}.info 'Testing/*' '/usr/*' --output-file ${_outputname}.info.cleaned
    COMMAND ${lcov_bin} -l ${_outputname}.info.cleaned
    #COMMAND ${genhtml_bin} -o ${_outputname} ${_outputname}.info.cleaned
    COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running the tests... then process the coverage. Output in ${CMAKE_BINARY_DIR}/${_outputname}/index.html"
	)
endfunction() 

function(single_coverage _targetname _outputname _binaryname _wdir)
  add_custom_target(${_targetname}
    ${lcov_bin} --directory ${_wdir}  --zerocounters
    COMMAND  mpirun -n 4 ${_binaryname}
    COMMAND ${lcov_bin} --directory ${CMAKE_BINARY_DIR} --capture --output-file ${_outputname}.info
    COMMAND ${lcov_bin} --remove ${_outputname}.info '${_wdir}/*' '/usr/*' --output-file ${_outputname}.info.cleaned
    COMMAND ${lcov_bin} -l ${_outputname}.info.cleaned
    #COMMAND ${genhtml_bin} -o ${_outputname}  ${_outputname}.info.cleaned
    COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/${_wdir}
    COMMENT "Running the binary in ${_wdir}... then process the coverage. Output in ${CMAKE_BINARY_DIR}/${_wdir}/${_outputname}/index.html"
	)
endfunction() 
