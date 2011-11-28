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

project(AKANTU)

enable_language(CXX)

#===============================================================================
# Misc.
#===============================================================================

set(AKANTU_CMAKE_DIR "${AKANTU_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
INCLUDE(${AKANTU_CMAKE_DIR}/AkantuMacros.cmake)

#===============================================================================
# Version Number
#===============================================================================

# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 0)
set(AKANTU_MINOR_VERSION 1)
execute_process(COMMAND ./get_revision.sh WORKING_DIRECTORY ${AKANTU_SOURCE_DIR} OUTPUT_VARIABLE AKANTU_BUILD_VERSION)
string(STRIP ${AKANTU_BUILD_VERSION} AKANTU_BUILD_VERSION)
MESSAGE(${AKANTU_BUILD_VERSION})
set(AKANTU_VERSION
  "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}.${AKANTU_BUILD_VERSION}"
  )
# Append the library version information to the library target properties
if(NOT AKANTU_NO_LIBRARY_VERSION)
  set(AKANTU_LIBRARY_PROPERTIES ${AKANTU_LIBRARY_PROPERTIES}
    VERSION "${AKANTU_VERSION}"
    SOVERSION "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}"

    )
endif(NOT AKANTU_NO_LIBRARY_VERSION)
#===============================================================================
# Package builder target module of CPack
#===============================================================================
SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Guillaume Anciaux") 
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6")
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Akantu library")
SET(CPACK_PACKAGE_VENDOR "LSMS")
SET(CPACK_PACKAGE_VERSION ${AKANTU_VERSION})
INCLUDE(CPack)

#===============================================================================
# Options
#===============================================================================
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
option(AKANTU_DOCUMENTATION "Build source documentation using Doxygen." OFF)

# features
add_optional_package(IOHelper "Add IOHelper support in akantu" OFF)
add_optional_package(QVIEW "Add QView support in akantu" OFF)
add_optional_package(BLAS "Use BLAS for arithmetic operations" OFF)
add_optional_package(EPSN "Use EPSN streering environment" OFF)


# Debug
if(NOT AKANTU_DEBUG)
  add_definitions(-DAKANTU_NDEBUG)
endif(NOT AKANTU_DEBUG)

#===========================================================================
# Dependencies
#===========================================================================
find_package(Boost REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_all_packages(${CMAKE_SOURCE_DIR}/packages)

#===========================================================================
# Library
#===========================================================================

add_subdirectory(src)

#===========================================================================
# Documentation
#===========================================================================
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
  include_directories(${AKANTU_INCLUDE_DIRS})
  find_package(GMSH REQUIRED)
  add_subdirectory(test)
endif(AKANTU_TESTS)

#===============================================================================
# Examples
#===============================================================================
option(AKANTU_EXAMPLES "Activate examples" OFF)
if(AKANTU_EXAMPLES)
  find_package(GMSH REQUIRED)
  add_subdirectory(examples)
endif(AKANTU_EXAMPLES)
