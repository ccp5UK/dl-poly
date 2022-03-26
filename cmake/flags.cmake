if (NOT FLAGS_SET)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

    if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "7.9.9" )
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe -fsanitize=undefined,leak,address" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
    else()
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42  -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -ftree-vectorize -funroll-loops -ffast-math" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
    set(CMAKE_Fortran_FLAGS_PROFILE  "-Ofast -ftree-vectorize -funroll-loops -ffast-math -pg" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)

  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -init=snan -init=arrays" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
    set(CMAKE_Fortran_FLAGS_PROFILE "-Ofast -pg -qopt-report=5" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)

  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")

    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -hfp3" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
    set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -hfp3 -h profile_generate" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)

  endif()
  set(FLAGS_SET 1 CACHE INTERNAL "Flags are set")
endif()
