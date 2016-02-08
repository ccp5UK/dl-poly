#.rst:
# FindPLUMED
# --------
#
# Find the native PLUMED includes and library.
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
#   PLUMED_INCLUDE_DIRS   - where to find the .h files, etc.
#   PLUMED_LIBRARIES      - List of libraries when using zlib.
#   PLUMED_FOUND          - True if zlib found.
#
# ::
#
# Hints
# ^^^^^
#
# A user may set ``PLUMED_ROOT`` to a PLUMED installation root to tell this
# module where to look.

#=============================================================================
# Copyright 2015 Kitware, Inc.
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

set(_PLUMED_SEARCHES)
set(_PLUMED_LIBS_SEARCH)
# Search PLUMED_ROOT first if it is set.
if(EXISTS "$ENV{PLUMED_ROOT}")
  set(_PLUMED_SEARCH_ROOT PATHS $ENV{PLUMED_ROOT} NO_DEFAULT_PATH)
  list(APPEND _PLUMED_SEARCHES _PLUMED_SEARCH_ROOT)
  list(APPEND _PLUMED_LIBS_SEARCH _PLUMED_SEARCH_ROOT)
endif()

string(REPLACE ":" ";" INC_LIST $ENV{INCLUDE})
foreach(IT ${INC_LIST})
  set(_PLUMED_SEARCH_INCLUDE_${IT} ${IT})
  list(APPEND _PLUMED_SEARCHES _PLUMED_SEARCH_INCLUDE_${IT})
endforeach()

string(REPLACE ":" ";" LIB_LIST $ENV{LD_LIBRARY_PATH})
foreach(IT ${LIB_LIST})
  set(_PLUMED_SEARCH_LIB_${IT} ${IT})
  list(APPEND _PLUMED_LIBS_SEARCH _PLUMED_SEARCH_LIB_${IT})
endforeach()

# Try each search configuration.
foreach(search ${_PLUMED_SEARCHES})
  find_path(WRAP_INCLUDE_DIR NAMES Plumed.h   ${${search}} HINTS ${${search}} PATH_SUFFIXES include "plumed/wrapper" "include/plumed/wrapper")
endforeach() 
get_filename_component(PLUMED_INCLUDE_DIR "${WRAP_INCLUDE_DIR}" DIRECTORY)

set(PLUMED_NAMES plumed)
foreach(search ${_PLUMED_LIBS_SEARCH})
  find_library(PLUMED_LIBRARY  NAMES ${PLUMED_NAMES}  ${${search}} HINTS ${${search}} PATH_SUFFIXES lib)
endforeach() 

set(VERSION_FILE "${PLUMED_INCLUDE_DIR}/config/version.h")
mark_as_advanced(PLUMED_LIBRARY PLUMED_INCLUDE_DIR)
if(PLUMED_INCLUDE_DIR AND EXISTS "${VERSION_FILE}")
    file(STRINGS "${VERSION_FILE}" PLUMED_H REGEX "^#define PLUMED_VERSION_SHORT \"*.*\"$")
    string(REGEX REPLACE "^.*_VERSION_SHORT \"([0-9]+.[0-9]+)\"$" "\\1" PLUMED_VERSION_STRING "${PLUMED_H}")
endif()
# handle the QUIETLY and REQUIRED arguments and set PLUMED_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PLUMED REQUIRED_VARS PLUMED_LIBRARY PLUMED_INCLUDE_DIR 
                                       VERSION_VAR PLUMED_VERSION_STRING)

if(PLUMED_FOUND)
    set(PLUMED_INCLUDE_DIRS ${PLUMED_INCLUDE_DIR})
    set(PLUMED_LIBRARIES ${PLUMED_LIBRARY})
endif()
