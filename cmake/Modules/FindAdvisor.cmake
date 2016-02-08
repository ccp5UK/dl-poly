#.rst:
# FindAdvisor
# --------
#
# Find the native Advisor includes and library.
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines :prop_tgt:`IMPORTED` target ``Advisor::Annotate``, if
# Advisor has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ::
#
#   ADVISOR_INCLUDE_DIRS   - where to find advisor-annotate.h, etc.
#   ADVISOR_MODULE_DIRS   - where to find advisor-annotate.mod, etc.
#   ADVISOR_LIBRARIES      - List of libraries when using zlib.
#   ADVISOR_FOUND          - True if zlib found.
#
# ::
#
# Hints
# ^^^^^
#
# A user may set ``ADVISOR_ROOT`` to a advisor installation root to tell this
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

set(_ADVISOR_SEARCHES)

# Search ADVISOR_ROOT first if it is set.
if(ADVISOR_ROOT)
  set(_ADVISOR_SEARCH_ROOT PATHS ${ADVISOR_ROOT} NO_DEFAULT_PATH)
  list(APPEND _ADVISOR_SEARCHES _ADVISOR_SEARCH_ROOT)
endif()

# Search ADVISOR_XE_2016_DIR first if it is set.
if(DEFINED ENV{ADVISOR_XE_2016_DIR})
  set(_ADVISOR_SEARCH_ROOT PATHS $ENV{ADVISOR_XE_2016_DIR} NO_DEFAULT_PATH)
  list(APPEND _ADVISOR_SEARCHES _ADVISOR_SEARCH_ROOT)
endif()

set(ADVISOR_NAMES advisor)
set(DL_NAMES dl)
set(libsuf "64")
set(modsuf "intel64")
# Try each search configuration.
foreach(search ${_ADVISOR_SEARCHES})
  find_path(ADVISOR_INCLUDE_DIR NAMES advisor-annotate.h   ${${search}} PATH_SUFFIXES include)
  find_path(ADVISOR_MODULE_DIR  NAMES advisor_annotate.mod ${${search}} PATH_SUFFIXES include/${modsuf})
  find_library(ADVISOR_LIBRARY  NAMES ${ADVISOR_NAMES}     ${${search}} PATH_SUFFIXES lib${libsuf})
endforeach()
mark_as_advanced(ADVISOR_LIBRARY ADVISOR_INCLUDE_DIR ADVISOR_MODULE_DIR)
if(ADVISOR_INCLUDE_DIR AND EXISTS "${ADVISOR_INCLUDE_DIR}/advisor-annotate.h")
    file(STRINGS "${ADVISOR_INCLUDE_DIR}/advisor-annotate.h" ADVISOR_H REGEX "^#define INTEL_ADVISOR_ANNOTATION_VERSION [^\"]*$")
    string(REGEX REPLACE "^.*_VERSION ([0-9]+).*$" "\\1" ADVISOR_VERSION_MAJOR "${ADVISOR_H}")
    string(REGEX REPLACE "^.*_VERSION [0-9]+\\.([0-9]+).*$" "\\1" ADVISOR_VERSION_MINOR "${ADVISOR_H}")
    set(ADVISOR_VERSION_STRING "${ADVISOR_VERSION_MAJOR}.${ADVISOR_VERSION_MINOR}")

endif()

# handle the QUIETLY and REQUIRED arguments and set ADVISOR_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ADVISOR REQUIRED_VARS ADVISOR_LIBRARY ADVISOR_INCLUDE_DIR
                                       VERSION_VAR ADVISOR_VERSION_STRING)

if(ADVISOR_FOUND)
    set(ADVISOR_INCLUDE_DIRS ${ADVISOR_INCLUDE_DIR})
    set(ADVISOR_MODULE_DIRS ${ADVISOR_MODULE_DIR})
    set(ADVISOR_LIBRARIES ${ADVISOR_LIBRARY})

    if(NOT TARGET Advisot::Annotate)
      add_library(Advisor::Annotate UNKNOWN IMPORTED)
      set_target_properties(Advisor::Annotate PROPERTIES
        IMPORTED_LOCATION "${ADVISOR_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ADVISOR_INCLUDE_DIRS}")
    endif()
endif()
