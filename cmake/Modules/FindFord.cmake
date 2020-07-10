#[=======================================================================[.rst:
FindFord
-----------

Ford is a documentation generation tool (see https://github.com/cmacmackin/ford/wiki).
This module looks for Doxygen and some optional tools it supports. These
tools are enabled as components in the :command:`find_package` command:

.. code-block:: cmake

find_package(Ford REQUIRED)

The following variables are defined by this module:

.. variable:: FORD_FOUND

True if the `ford`` executable was found.

.. variable:: FORD_VERSION

The version reported by ``ford --version``.

.. variable:: FORD_EXECUTABLE

the ford executable.

#]=======================================================================]

# For backwards compatibility support
if(Ford_FIND_QUIETLY)
  set(Ford_FIND_QUIETLY TRUE)
endif()

#
# Find Ford...
#
find_program(
  FORD_EXECUTABLE
  NAMES ford
  PATHS
  DOC "Ford documentation generation tool (https://github.com/cmacmackin/ford/wiki)"
  )
mark_as_advanced(FORD_EXECUTABLE)

if(FORD_EXECUTABLE)
  execute_process(
    COMMAND "${FORD_EXECUTABLE}" --version
    OUTPUT_VARIABLE fver
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _Ford_version_result
    )
  if(_Ford_version_result)
    message(WARNING "Unable to determine doxygen version: ${_Ford_version_result}")
  endif()
  string(REGEX REPLACE "^FORD, version ([0-9]+.*.*)" "\\1"  FORD_VERSION "${fver}")
endif()


# Verify find results
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Ford
  REQUIRED_VARS FORD_EXECUTABLE
  VERSION_VAR FORD_VERSION
  )
if(FORD_FOUND)
  set(FORD_FOUND "YES")
else()
  set(FORD_FOUND "NO")
endif()

