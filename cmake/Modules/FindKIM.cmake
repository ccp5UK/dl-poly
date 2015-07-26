#.rst:
# FindKIM
# --------
#
# Find the native KIM includes and library.
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
#   KIM_INCLUDE_DIRS   - where to find the .h files, etc.
#   KIM_MODULE_DIRS   - where to find the .mod files, etc.
#   KIM_LIBRARIES      - List of libraries when using zlib.
#   KIM_FOUND          - True if zlib found.
#
# ::
#
# Hints
# ^^^^^
#
# A user may set ``KIM_ROOT`` to a KIM installation root to tell this
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

set(_KIM_SEARCHES)
set(_KIM_LIBS_SEARCH)
# Search KIM_ROOT first if it is set.
if(EXISTS "$ENV{KIM_ROOT}")
  set(_KIM_SEARCH_ROOT PATHS $ENV{KIM_ROOT} NO_DEFAULT_PATH)
  list(APPEND _KIM_SEARCHES _KIM_SEARCH_ROOT)
  list(APPEND _KIM_LIBS_SEARCH _KIM_SEARCH_ROOT)
endif()

string(REPLACE ":" ";" INC_LIST $ENV{INCLUDE})
foreach(IT ${INC_LIST})
  set(_KIM_SEARCH_INCLUDE_${IT} ${IT})
  list(APPEND _KIM_SEARCHES _KIM_SEARCH_INCLUDE_${IT})
endforeach()

string(REPLACE ":" ";" LIB_LIST $ENV{LD_LIBRARY_PATH})
foreach(IT ${LIB_LIST})
  set(_KIM_SEARCH_LIB_${IT} ${IT})
  list(APPEND _KIM_LIBS_SEARCH _KIM_SEARCH_LIB_${IT})
endforeach()

# Try each search configuration.
foreach(search ${_KIM_SEARCHES})
  find_path(KIM_INCLUDE_DIR NAMES KIM_API_Version.h   ${${search}} HINTS ${${search}} PATH_SUFFIXES include "include/kim-api")
  find_path(KIM_MODULE_DIR NAMES kim_api_f03.mod   ${${search}} HINTS ${${search}} PATH_SUFFIXES include "include/kim-api")
endforeach() 

set(KIM_NAMES kim-api)
foreach(search ${_KIM_LIBS_SEARCH})
  find_library(KIM_LIBRARY  NAMES ${KIM_NAMES}  ${${search}} HINTS ${${search}} PATH_SUFFIXES lib)
endforeach() 

mark_as_advanced(KIM_LIBRARY KIM_INCLUDE_DIR KIM_MODULE_DIR)
if(KIM_INCLUDE_DIR AND EXISTS "${KIM_INCLUDE_DIR}/KIM_API_Version.h")
    file(STRINGS "${KIM_INCLUDE_DIR}/KIM_API_Version.h" KIM_H REGEX "^#define KIM_API_VERSION_MAJOR [^\"]*$")
    string(REGEX REPLACE "^.*_VERSION_MAJOR ([0-9]+).*.*$" "\\1" KIM_VERSION_MAJOR "${KIM_H}")
    file(STRINGS "${KIM_INCLUDE_DIR}/KIM_API_Version.h" KIM_H REGEX "^#define KIM_API_VERSION_MINOR [^\"]*$")
    string(REGEX REPLACE "^.*_VERSION_MINOR *.([0-9]+).*$" "\\1" KIM_VERSION_MINOR "${KIM_H}")
    file(STRINGS "${KIM_INCLUDE_DIR}/KIM_API_Version.h" KIM_H REGEX "^#define KIM_API_VERSION_PATCH [^\"]*$")
    string(REGEX REPLACE "^.*_VERSION_PATCH *.*.([0-9]+)$" "\\1" KIM_VERSION_PATCH "${KIM_H}")
    set(KIM_VERSION_STRING "${KIM_VERSION_MAJOR}.${KIM_VERSION_MINOR}.${KIM_VERSION_PATCH}")
endif()
# handle the QUIETLY and REQUIRED arguments and set KIM_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(KIM REQUIRED_VARS KIM_LIBRARY KIM_INCLUDE_DIR KIM_MODULE_DIR
                                       VERSION_VAR KIM_VERSION_STRING)

if(KIM_FOUND)
    set(KIM_INCLUDE_DIRS ${KIM_INCLUDE_DIR})
    set(KIM_MODULE_DIRS ${KIM_MODULE_DIR})
    set(KIM_LIBRARIES ${KIM_LIBRARY})
endif()
