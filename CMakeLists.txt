#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
#            _    ,-,    _
#     ,--, /: :\/': :`\/: :\
#    |`;  ' `,'   `.;    `: |                  _               _
#    |    |     |  '  |     |.                | |             | |
#    | :  | F E | A S | T++ ||            __ _| | ____ _ _ __ | |_ _   _
#    | :. |  :  |  :  |  :  | \          / _` | |/ / _` | '_ \| __| | | |
#     \__/: :.. : :.. | :.. |  )        | (_| |   < (_| | | | | |_| |_| |
#          `---',\___/,\___/ /'          \__,_|_|\_\__,_|_| |_|\__|\__,_|
#               `==._ .. . /'
#                    `-::-'
#===============================================================================

#===============================================================================
# CMake Project
#===============================================================================

cmake_minimum_required(VERSION 2.6)

project(akantu)

enable_language(CXX)

#===============================================================================
# Misc.
#===============================================================================

set(AKANTU_CMAKE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

include(${AKANTU_CMAKE_DIR}/AkantuMacros.cmake)

#===============================================================================
# Version Number
#===============================================================================

# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 0)
set(AKANTU_MINOR_VERSION 1)
set(AKANTU_PATCH_VERSION 0)

find_package(Subversion)
if (EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  if(SUBVERSION_FOUND)
    subversion_wc_info(${CMAKE_SOURCE_DIR} MY)
    set(AKANTU_BUILD_VERSION ${MY_WC_REVISION})
    set(AKANTU_VERSION
      "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}.${AKANTU_PATCH_VERSION}.${AKANTU_BUILD_VERSION}"
      )
    file(WRITE VERSION "${AKANTU_VERSION}\n")
  else(SUBVERSION_FOUND)
    message("SVN control files were found but no subversion executable is present... ")
    set(AKANTU_VERSION 0)
  endif(SUBVERSION_FOUND)
else(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  if(EXISTS ${CMAKE_SOURCE_DIR}/VERSION)
    file(STRINGS VERSION AKANTU_VERSION)
  else(EXISTS ${CMAKE_SOURCE_DIR}/VERSION)
    message("No SVN control file neither VERSION file could be found. How was this release made ?")
  endif(EXISTS ${CMAKE_SOURCE_DIR}/VERSION)
endif(EXISTS ${CMAKE_SOURCE_DIR}/.svn)


# Append the library version information to the library target properties
if(NOT AKANTU_NO_LIBRARY_VERSION)
  set(AKANTU_LIBRARY_PROPERTIES ${AKANTU_LIBRARY_PROPERTIES}
    VERSION "${AKANTU_VERSION}"
    SOVERSION "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}"
    )
endif(NOT AKANTU_NO_LIBRARY_VERSION)

#===============================================================================
# Options
#===============================================================================
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
option(AKANTU_DOCUMENTATION "Build source documentation using Doxygen." OFF)

# features
add_optional_package(BLAS "Use BLAS for arithmetic operations" OFF)
if(BLAS_mkl_core_LIBRARY)
  set(AKANTU_USE_BLAS_MKL)
endif()

add_optional_package(LAPACK "Use LAPACK for arithmetic operations" OFF)

# Debug
if(NOT AKANTU_DEBUG)
  set(AKANTU_NDEBUG ON)
else()
  set(AKANTU_NDEBUG OFF)
endif(NOT AKANTU_DEBUG)

if(CMAKE_CXX_FLAGS STREQUAL "")
  set(CMAKE_CXX_FLAGS "-Wall" CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()

set(CMAKE_CXX_FLAGS_PROFILING "-g -pg" CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING)
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_PROFILING)
mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_PROFILING)

#===============================================================================
# Dependencies
#===============================================================================
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_all_packages(${CMAKE_SOURCE_DIR}/packages)
if(AKANTU_SCOTCH)
  set(AKANTU_PARTITIONER ON)
else()
  set(AKANTU_PARTITIONER OFF)
endif()

if(AKANTU_MUMPS)
  set(AKANTU_SOLVER ON)
else()
  set(AKANTU_SOLVER OFF)
endif()

#add_akantu_definitions()

#===============================================================================

add_subdirectory(src)

#===============================================================================
# Documentation
#===============================================================================
if(AKANTU_DOCUMENTATION)
  find_package(Doxygen REQUIRED)

  if(DOXYGEN_FOUND)
    set(DOXYGEN_WARNINGS NO)
    set(DOXYGEN_QUIET YES)
    if(CMAKE_VERBOSE_MAKEFILE)
      set(DOXYGEN_WARNINGS YES)
      set(DOXYGEN_QUIET NO)
    endif(CMAKE_VERBOSE_MAKEFILE)

    add_subdirectory(doc/doxygen)
  endif(DOXYGEN_FOUND)
endif(AKANTU_DOCUMENTATION)

#===============================================================================
# Tests
#===============================================================================
ENABLE_TESTING()
INCLUDE(CTest)

INCLUDE(${AKANTU_CMAKE_DIR}/AkantuTestAndExamples.cmake)

option(AKANTU_TESTS "Activate tests" OFF)
if(AKANTU_TESTS)
  include_directories(
    ${AKANTU_INCLUDE_DIRS}
    ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR}
    )
  find_package(GMSH REQUIRED)
  add_subdirectory(test)
endif(AKANTU_TESTS)

#===============================================================================
# Config gen for external packages
#===============================================================================
export(TARGETS akantu
  FILE "${CMAKE_BINARY_DIR}/AkantuLibraryDepends.cmake")

export(PACKAGE Akantu)
configure_file(cmake/AkantuBuildTreeSettings.cmake.in  "${CMAKE_BINARY_DIR}/AkantuBuildTreeSettings.cmake" @ONLY)

file(WRITE "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
#===============================================================================
# @file   AkantuConfigInclude.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

")

foreach(_option ${AKANTU_PACKAGE_NAMES_LIST_ALL})
  list(FIND AKANTU_OPTION_LIST ${_option} _index)
  if (_index EQUAL -1)
    if(NOT "${_option}" STREQUAL "CORE")
      if(NOT AKANTU_${_option})
        set(AKANTU_${_option} OFF)
      endif()
      file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_HAS_${_option} ${AKANTU_${_option}})")
    endif()
  endif()
endforeach()

file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake"
"

set(AKANTU_HAS_PARTITIONER  ${AKANTU_PARTITIONER})
set(AKANTU_HAS_SOLVER       ${AKANTU_SOLVER})

set(AKANTU_OPTION_LIST ${AKANTU_OPTION_LIST})
")

foreach(_option ${AKANTU_OPTION_LIST})
  file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_USE_${_option} ${AKANTU_${_option}})
set(AKANTU_${_option}_LIBRARIES ${AKANTU_${_option}_LIBRARIES})
set(AKANTU_${_option}_INCLUDE_DIR ${AKANTU_${_option}_INCLUDE_DIR})
")
endforeach()

# Create the AkantuConfig.cmake and AkantuConfigVersion files
get_filename_component(CONF_REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
configure_file(cmake/AkantuConfig.cmake.in "${CMAKE_BINARY_DIR}/AkantuConfig.cmake" @ONLY)
configure_file(cmake/AkantuConfigVersion.cmake.in "${CMAKE_BINARY_DIR}/AkantuConfigVersion.cmake" @ONLY)
configure_file(cmake/AkantuUse.cmake "${CMAKE_BINARY_DIR}/AkantuUse.cmake" COPYONLY)

# Install the export set for use with the install-tree
install(EXPORT AkantuLibraryDepends DESTINATION  lib/akantu
  COMPONENT dev)

install(FILES ${CMAKE_BINARY_DIR}/AkantuConfig.cmake
  ${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake
  ${CMAKE_BINARY_DIR}/AkantuConfigVersion.cmake
  ${CMAKE_SOURCE_DIR}/cmake/AkantuUse.cmake
  DESTINATION  lib/akantu
  COMPONENT dev)

install(FILES
  ${CMAKE_SOURCE_DIR}/cmake/FindIOHelper.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindMumps.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindScotch.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindGMSH.cmake
  DESTINATION  lib/akantu/cmake
  COMPONENT dev)

#===============================================================================
# Package builder target module of CPack
#===============================================================================
set(CPACK_GENERATOR "DEB;TGZ;TBZ2;STGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "guillaume.anciaux@epfl.ch, nicolas.richart@epfl.ch")
#message(${CPACK_DEBIAN_PACKAGE_DEPENDS})
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64" CACHE STRING "Architecture of debian package generation")
else()
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386" CACHE STRING "Architecture of debian package generation")
endif()
#message("CPACK_DEBIAN_PACKAGE_ARCHITECTURE = ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Akantu library")
set(CPACK_PACKAGE_VENDOR "LSMS")
set(PACKAGE_FILE_NAME "akantu" CACHE STRING "Name of package to be generated")
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-src")
set(CPACK_PACKAGE_VERSION "${AKANTU_VERSION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
list(APPEND CPACK_SOURCE_IGNORE_FILES ${AKANTU_EXCLUDE_SOURCE_FILE})
list(APPEND CPACK_SOURCE_IGNORE_FILES ${AKANTU_UNACTIVATED_PACKAGE_FILE})
list(APPEND CPACK_SOURCE_IGNORE_FILES "/.*build.*/;/CVS/;/\\\\.svn/;/\\\\.bzr/;/\\\\.hg/;/\\\\.git/;\\\\.swp$;\\\\.#;/#;~")
include(CPack)

