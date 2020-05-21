if(CMAKE_BUILD_TYPE)
  if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
      if (${CMAKE_Fortran_COMPILER_VERSION} STRGREATER "6.9.9" )
        set(CMAKE_Fortran_FLAGS_DEBUG  "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe")
      else()
        set(CMAKE_Fortran_FLAGS_DEBUG  "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42  -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe")
      endif()
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -qopt-report=5 -init=snan -init=arrays -check noarg_temp_created")
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
    endif()
  endif()
  if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
      set(CMAKE_Fortran_FLAGS_RELEASE  "-g -Ofast -ftree-vectorize -funroll-loops -ffast-math ")
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS_RELEASE "-g -Ofast -mtune=native")
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -hfp3")
    endif()
  endif()
endif()
