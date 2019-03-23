if(CMAKE_BUILD_TYPE)
  if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
      set(CMAKE_Fortran_FLAGS_DEBUG  "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe")
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -qopt-report=5 -init=snan -init=arrays") 
    endif()
  endif()
endif()
