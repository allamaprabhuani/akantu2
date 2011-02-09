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
#    |    |     |  '  |     |.		      | |             | |
#    | :  | F E | A S | T++ ||		  __ _| | ____ _ _ __ | |_ _   _
#    | :. |  :  |  :  |  :  | \		 / _` | |/ / _` | '_ \| __| | | |
#     \__/: :.. : :.. | :.. |  )        | (_| |   < (_| | | | | |_| |_| |
#          `---',\___/,\___/ /'		 \__,_|_|\_\__,_|_| |_|\__|\__,_|
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

#===============================================================================
# Version Number
#===============================================================================

# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 0)
set(AKANTU_MINOR_VERSION 1)
set(AKANTU_BUILD_VERSION 0)
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
# Options
#===============================================================================
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
option(AKANTU_DOCUMENTATION "Build source documentation using Doxygen." OFF)

# compilation switch
option(AKANTU_BUILD_CONTACT "Build the contact algorithm" OFF)

# features
option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" OFF)
option(AKANTU_USE_MPI "Add MPI support in akantu" OFF)
option(AKANTU_USE_SCOTCH "Add Scotch support in akantu" OFF)
option(AKANTU_USE_MUMPS "Add Mumps support in akantu" OFF)
option(AKANTU_USE_BLAS "Use BLAS for arithmetic operations" OFF)
option(AKANTU_USE_EPSN "Use EPSN streering environment" OFF)


#===============================================================================
# Options handling
#===============================================================================

#===========================================================================
# Documentation
#===========================================================================
find_package(Boost REQUIRED)
if(BOOST_FOUND)
  include_directories(BOOST_INCLUDE_PATH)
endif(BOOST_FOUND)


# Debug
if(NOT AKANTU_DEBUG)
  add_definitions(-DAKANTU_NDEBUG)
endif(NOT AKANTU_DEBUG)

# IOHelper
if(AKANTU_USE_IOHELPER)
  find_package(IOHelper REQUIRED)
  if(IOHELPER_FOUND)
    add_definitions(-DAKANTU_USE_IOHELPER)
    set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
      ${IOHELPER_INCLUDE_PATH}
      )
    set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
      ${IOHELPER_LIBRARIES}
      )
  endif(IOHELPER_FOUND)
endif(AKANTU_USE_IOHELPER)

# MPI
set(AKANTU_MPI_ON FALSE)
if(AKANTU_USE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    add_definitions(-DAKANTU_USE_MPI)
    set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
      ${MPI_INCLUDE_PATH}
      )
    set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
      ${MPI_LIBRARY}
      ${MPI_EXTRA_LIBRARY}
      )
    set(AKANTU_MPI_ON TRUE)
  endif(MPI_FOUND)
endif(AKANTU_USE_MPI)

# Scotch
set(AKANTU_SCOTCH_ON FALSE)
if(AKANTU_USE_SCOTCH)
  find_package(Scotch REQUIRED)
  if(SCOTCH_FOUND)
    add_definitions(-DAKANTU_USE_SCOTCH)
    set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
      ${SCOTCH_INCLUDE_PATH}
      )
    set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
      ${SCOTCH_LIBRARIES}
      )
    set(AKANTU_SCOTCH_ON TRUE)
  endif(SCOTCH_FOUND)
endif(AKANTU_USE_SCOTCH)

# MUMPS
set(AKANTU_MUMPS_ON FALSE)
if(AKANTU_USE_MUMPS)
  find_package(Mumps REQUIRED)
  if(MUMPS_FOUND)
    add_definitions(-DAKANTU_USE_MUMPS)
    set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
      ${MUMPS_INCLUDE_PATH}
      )
    set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
      ${MUMPS_LIBRARIES}
      )
    set(AKANTU_MUMPS_ON TRUE)
  endif(MUMPS_FOUND)
endif(AKANTU_USE_MUMPS)

# BLAS
set(AKANTU_BLAS_ON FALSE)
if(AKANTU_USE_BLAS)
  enable_language(Fortran)
  find_package(BLAS)
  if(BLAS_FOUND)
    add_definitions(-DAKANTU_USE_CBLAS)
    set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
      ${BLAS_LIBRARIES}
      )
    set(AKANTU_BLAS_ON TRUE)
  endif(BLAS_FOUND)
endif(AKANTU_USE_BLAS)

# EPSN
set(AKANTU_EPSN_ON FALSE)
if(AKANTU_USE_EPSN)
  find_package(EPSN)
  if(EPSN_FOUND)
    add_definitions(-DAKANTU_USE_EPSN)
    set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
      ${EPSN_INCLUDE_DIR}
      )
    set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
      ${EPSN_LIBRARIES}
      )
    set(AKANTU_EPSN_ON TRUE)
  endif(EPSN_FOUND)
endif(AKANTU_USE_EPSN)

#===============================================================================
# Library
#===============================================================================

set(AKANTU_COMMON_SRC
  common/aka_common.cc
  common/aka_error.cc
  common/aka_extern.cc
  common/aka_static_memory.cc
  common/aka_memory.cc
  common/aka_vector.cc
  common/aka_math.cc
  fem/mesh.cc
  fem/fem.cc
  fem/element_class.cc
#  model/integration_scheme/central_difference.cc
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_boundary.cc
  model/solid_mechanics/solid_mechanics_model_material.cc
  model/solid_mechanics/material.cc
  model/solid_mechanics/materials/material_elastic.cc
  mesh_utils/mesh_io.cc
  mesh_utils/mesh_pbc.cc
  mesh_utils/mesh_io/mesh_io_msh.cc
  mesh_utils/mesh_partition.cc
  mesh_utils/mesh_utils.cc
  solver/sparse_matrix.cc
  solver/solver.cc
  synchronizer/ghost_synchronizer.cc
  synchronizer/synchronizer.cc
  synchronizer/communicator.cc
  synchronizer/static_communicator.cc
  )

set(AKANTU_CONTACT_SRC
  model/solid_mechanics/contact.cc
  model/solid_mechanics/contact_search.cc
  model/solid_mechanics/contact_neighbor_structure.cc
  model/solid_mechanics/contact/contact_2d_explicit.cc
  model/solid_mechanics/contact/contact_search_2d_explicit.cc
  model/solid_mechanics/contact/regular_grid_neighbor_structure.cc
  model/solid_mechanics/contact/contact_search_3d_explicit.cc
  model/solid_mechanics/contact/contact_3d_explicit.cc
  model/solid_mechanics/contact/grid_2d_neighbor_structure.cc
  model/solid_mechanics/contact/contact_rigid_no_friction.cc
  )

if(AKANTU_BUILD_CONTACT)
  set(AKANTU_COMMON_SRC
    ${AKANTU_COMMON_SRC}
    ${AKANTU_CONTACT_SRC})
endif(AKANTU_BUILD_CONTACT)

if(AKANTU_USE_MPI AND MPI_FOUND)
  set(AKANTU_COMMON_SRC
    ${AKANTU_COMMON_SRC}
    synchronizer/static_communicator_mpi.cc
    )
endif(AKANTU_USE_MPI AND MPI_FOUND)


if(AKANTU_SCOTCH_ON)
  set(AKANTU_COMMON_SRC
    ${AKANTU_COMMON_SRC}
    mesh_utils/mesh_partition/mesh_partition_scotch.cc
    )
endif(AKANTU_SCOTCH_ON)


if(AKANTU_MUMPS_ON)
  set(AKANTU_COMMON_SRC
    ${AKANTU_COMMON_SRC}
    solver/solver_mumps.cc
    )
endif(AKANTU_MUMPS_ON)


set(AKANTU_INCLUDE_DIRS
  common
  fem/
  mesh_utils/
  mesh_utils/mesh_io/
  mesh_utils/mesh_partition/
  model/
  model/integration_scheme
  model/solid_mechanics
  model/solid_mechanics/materials
  model/solid_mechanics/contact
  synchronizer/
  solver/
  )

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
  )

add_library(akantu ${AKANTU_COMMON_SRC})
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})


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

option(AKANTU_TESTS "Activate tests" OFF)
if(AKANTU_TESTS)
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
