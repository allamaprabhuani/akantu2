#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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
#                      _               _
#                     | |             | |
#                 __ _| | ____ _ _ __ | |_ _   _
#                / _` | |/ / _` | '_ \| __| | | |
#               | (_| |   < (_| | | | | |_| |_| |
#                \__,_|_|\_\__,_|_| |_|\__|\__,_|
#
#===============================================================================

#===============================================================================
# CMake Project
#===============================================================================
cmake_minimum_required(VERSION 2.6)

project(akantu)

enable_language(CXX)

#===============================================================================
# Misc. config for cmake
#===============================================================================
set(AKANTU_CMAKE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

include(${AKANTU_CMAKE_DIR}/AkantuMacros.cmake)
include(${AKANTU_CMAKE_DIR}/AkantuPackageSystem.cmake)

#===============================================================================
# Version Number
#===============================================================================
# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 1)
set(AKANTU_MINOR_VERSION 0)
set(AKANTU_PATCH_VERSION 0)

include(${AKANTU_CMAKE_DIR}/AkantuVersion.cmake)

#===============================================================================
# Options
#===============================================================================
# Debug
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
mark_as_advanced(AKANTU_DEBUG)

add_cxx_flags("-Wall")

if(NOT AKANTU_DEBUG)
  add_cxx_flags("-Wno-return-type")
  set(AKANTU_NDEBUG ON)
else()
  remove_cxx_flags("-Wno-return-type")
endif()

#Profiling
set(CMAKE_CXX_FLAGS_PROFILING "-g -pg" CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING)
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_PROFILING)
mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_PROFILING)

include(${AKANTU_CMAKE_DIR}/AkantuDetermineCCompiler.cmake)

#===============================================================================
# Dependencies
#===============================================================================
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_all_packages(${CMAKE_SOURCE_DIR}/packages)

## meta option \todo better way to do it when multiple package give enable the
## same feature
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

#===============================================================================
# Akantu library
#===============================================================================
add_subdirectory(src)

#===============================================================================
# Documentation
#===============================================================================
#Documentation
option(AKANTU_DOCUMENTATION "Build source documentation using Doxygen." OFF)
if(AKANTU_DOCUMENTATION)
  add_subdirectory(doc)
endif()

#===============================================================================
# Tests
#===============================================================================
option(AKANTU_TESTS "Activate tests" OFF)
if(AKANTU_TESTS)
  add_subdirectory(test)
endif(AKANTU_TESTS)

#===============================================================================
# Install and Packaging
#===============================================================================
include(${AKANTU_CMAKE_DIR}/AkantuInstall.cmake)
include(${AKANTU_CMAKE_DIR}/AkantuCPack.cmake)
