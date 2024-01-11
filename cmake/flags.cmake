if (NOT FLAGS_SET)
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

      if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "7.9.9" )
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2018 -pedantic -fbacktrace    -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived  -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe " CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
        set(CMAKE_Fortran_FLAGS_DEBUG2 "-static-libasan -g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2018 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe -fsanitize=undefined,leak,address" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
      else()
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2018 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647 -finit-real=snan -finit-logical=true -finit-character=42  -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
      endif()
      set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -std=f2018 -ftree-vectorize -funroll-loops -ffast-math" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
      set(CMAKE_Fortran_FLAGS_PROFILE  "-Ofast -std=f2018  -ftree-vectorize -funroll-loops -ffast-math -pg" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)

    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

      set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -stand f18 -traceback -C -fp-stack-check -ftrapuv -init=snan -init=arrays" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
      set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -stand f18" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
      set(CMAKE_Fortran_FLAGS_PROFILE "-Ofast -stand f18 -pg -qopt-report=5" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)
 
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
      # ifx compiler
      set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -O0 -std18 -check all -warn all -fp-model:consistent" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
      set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -std18 -fp-model:consistent" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
      set(CMAKE_Fortran_FLAGS_PROFILE "-g -traceback -O0 -std18 -check all -warn all -fp-model:consistent -qopt-report" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)

    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")

      set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -hfp3" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
      set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -hfp3 -h profile_generate" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)

    endif()
  elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS_DEBUG "/Od /fpp" CACHE STRING "Flags used by the Fortran compiler during DEBUG builds." FORCE)
      set(CMAKE_Fortran_FLAGS_RELEASE "/fast /fpp" CACHE STRING "Flags used by the Fortran compiler during RELEASE builds." FORCE)
      set(CMAKE_Fortran_FLAGS_PROFILE "/fast /fpp" CACHE STRING "Flags used by the Fortran compiler during PROFILE builds." FORCE)
    endif()
  elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  else()
    message(INFO "no flags preset")
  endif()

  set(FLAGS_SET 1 CACHE INTERNAL "Flags are set")
endif()
