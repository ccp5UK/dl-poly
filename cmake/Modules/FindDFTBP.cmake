#.rst:
# FindDFTBP
# --------
#
# Find the native DFTBP++ includes and library.
#
# IMPORTED Targets
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ::
#
#   DFTBP_INCLUDE_DIRS   - where to find the .h files, etc.
#   DFTBP_LIBRARIES      - List of libraries when using zlib.
#   DFTBP_FOUND          - True if zlib found.
#
# ::
#
# Hints
# ^^^^^
#
# A user may set ``DFTBP_ROOT`` to a DFTB++ installation root to tell this
# module where to look.

#=============================================================================
# Copyright 2019 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

set(_DFTBP_SEARCHES)
set(_DFTBP_LIBS_SEARCH)
# Search DFTBP_ROOT first if it is set.
if(DEFINED ENV{DFTB_ROOT})
  set(_DFTBP_SEARCH_ROOT PATHS $ENV{DFTB_ROOT} NO_DEFAULT_PATH)
  list(APPEND _DFTBP_SEARCHES _DFTBP_SEARCH_ROOT)
  list(APPEND _DFTBP_LIBS_SEARCH _DFTBP_SEARCH_ROOT)
endif()

if(EXISTS $ENV{INCLUDE})
  string(REPLACE ":" ";" INC_LIST $ENV{INCLUDE})
  foreach(IT ${INC_LIST})
    set(_DFTBP_SEARCH_INCLUDE_${IT} ${IT})
    list(APPEND _DFTBP_SEARCHES _DFTBP_SEARCH_INCLUDE_${IT})
  endforeach()
endif()


if(EXISTS $ENV{LD_LIBRARY_PATH})
  string(REPLACE ":" ";" LIB_LIST $ENV{LD_LIBRARY_PATH})
  foreach(IT ${LIB_LIST})
    set(_DFTBP_SEARCH_LIB_${IT} ${IT})
    list(APPEND _DFTBP_LIBS_SEARCH _DFTBP_SEARCH_LIB_${IT})
  endforeach()
endif()

# Try each search configuration.
foreach(search ${_DFTBP_SEARCHES})
  find_path(WRAP_INCLUDE_DIR NAMES typegeometry.mod  ${${search}} HINTS ${${search}} PATH_SUFFIXES include)
endforeach() 
get_filename_component(DFTBP_INCLUDE_DIR "${WRAP_INCLUDE_DIR}/include" DIRECTORY)

set(DFTBP_NAMES dftb)
foreach(search ${_DFTBP_LIBS_SEARCH})
  find_library(DFTBP_LIBRARY  NAMES ${DFTBP_NAMES}  ${${search}} HINTS ${${search}} PATH_SUFFIXES lib "lib64")
endforeach() 

mark_as_advanced(DFTBP_LIBRARY DFTBP_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set DFTBP_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DFTBP REQUIRED_VARS DFTBP_LIBRARY DFTBP_INCLUDE_DIR)

if(DFTBP_FOUND)
  set(DFTBP_INCLUDE_DIRS ${DFTBP_INCLUDE_DIR})
  set(DFTBP_LIBRARIES ${DFTBP_LIBRARY})
endif()
