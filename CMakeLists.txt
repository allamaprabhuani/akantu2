#===============================================================================
# @file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
# @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
# @author Marion Estelle Chambart <marion.chambart@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Alejandro Marcos Aragon <alejandro.aragon@epfl.ch>
# @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
# @author David Simon Kammer <david.kammer@epfl.ch>
#
# @date   Mon Jun 14 19:12:20 2010
#
# @brief  main configuration file
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
#-------------------------------------------------------------------------------
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
cmake_minimum_required(VERSION 2.8)
project(akantu)
enable_language(CXX)

#===============================================================================
# Misc. config for cmake
#===============================================================================
set(AKANTU_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

if(NOT AKANTU_TARGETS_EXPORT)
  set(AKANTU_TARGETS_EXPORT AkantuLibraryDepends)
endif()

include(CMakeVersionGenerator)
include(CMakePackagesSystem)
include(CMakeFlagsHandling)
include(AkantuMacros)

cmake_deactivate_debug_message()

#===============================================================================
# Version Number
#===============================================================================
# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 1)
set(AKANTU_MINOR_VERSION 0)
set(AKANTU_PATCH_VERSION 0)

define_project_version()

#===============================================================================
# Options
#===============================================================================
# Debug
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
mark_as_advanced(AKANTU_DEBUG)

add_flags(cxx "-Wall")

if(NOT AKANTU_DEBUG)
  add_flags(cxx "-Wno-return-type")
  set(AKANTU_NDEBUG ON)
else()
  remove_flags(cxx "-Wno-return-type")
endif()

#Profiling
set(CMAKE_CXX_FLAGS_PROFILING "-g -pg" CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING)
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_PROFILING)
mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_PROFILING)

include(AkantuDetermineCCompiler)

#===============================================================================
# Dependencies
#===============================================================================
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_all_packages(${PROJECT_SOURCE_DIR}/packages)

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
# Examples and tests
#===============================================================================
option(AKANTU_EXAMPLES "Activate examples" OFF)
option(AKANTU_TESTS "Activate tests" OFF)

if(AKANTU_EXAMPLES OR AKANTU_TESTS)
  include(AkantuTestAndExamples)
  find_package(GMSH REQUIRED)
endif()

if(AKANTU_EXAMPLES)
  add_subdirectory(examples)
endif()

if(AKANTU_TESTS)
  enable_testing()
  include(CTest)
  set(AKANTU_TESTS_EXCLUDE_FILES)
  add_subdirectory(test)
endif()

#===============================================================================
# Install and Packaging
#===============================================================================
include(AkantuInstall)
include(AkantuCPack)
