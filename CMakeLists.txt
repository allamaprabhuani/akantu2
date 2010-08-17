#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# <insert lisence here>
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

# Debug
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
if(NOT AKANTU_DEBUG)
  add_definitions(-DAKANTU_NDEBUG)
endif(NOT AKANTU_DEBUG)

# IOHelper
option(AKANTU_USE_IOHELPER "Add IOHelper support in akantu" OFF)
if(AKANTU_USE_IOHELPER)
  include(${AKANTU_CMAKE_DIR}/FindIOHelper.cmake)
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

#===============================================================================
# Library
#===============================================================================

set(AKANTU_COMMON_SRC
  common/aka_common.cc
  common/aka_extern.cc
  common/aka_static_memory.cc
  common/aka_memory.cc
  common/aka_vector.cc
  common/aka_math.cc
  fem/mesh.cc
  fem/fem.cc
  fem/element_class.cc
  model/solid_mechanics_model.cc
  model/material.cc
  model/materials/material_elastic.cc
  model/integration_scheme/central_difference.cc
  mesh_utils/mesh_io.cc
  mesh_utils/mesh_io/mesh_io_msh.cc
  mesh_utils/mesh_partition.cc
  )

set(AKANTU_INCLUDE_DIRS
  common
  fem/
  mesh_utils/
  mesh_utils/mesh_io
  model/
  )

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
  )

add_library(akantu ${AKANTU_COMMON_SRC})
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})



#===============================================================================
# Tests
#===============================================================================

option(BUILD_UNIT_TEST "Build the unit test")
if(BUILD_UNIT_TEST)
  set(AKANTU_UNIT_TESTS
    test/test_static_memory.cc
    test/test_mesh_io_msh.cc
    test/test_vector.cc
    test/test_fem.cc
    test/test_solid_mechanics_model.cc
    test/test_solid_mechanics_model_cube3d.cc
    test/test_solid_mechanics_model_bar_traction2d.cc
    )

  foreach(_test ${AKANTU_UNIT_TESTS})
    get_filename_component(_target "${_test}" NAME_WE)
    add_executable(${_target} ${_test})
    target_link_libraries(${_target} akantu ${AKANTU_EXTERNAL_LIBRARIES})
  endforeach(_test)
endif(BUILD_UNIT_TEST)