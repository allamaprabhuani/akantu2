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
#    |    |     |  '  |     |.	              | |             | |        
#    | :  | F E | A S | T++ ||	          __ _| | ____ _ _ __ | |_ _   _ 
#    | :. |  :  |  :  |  :  | \	         / _` | |/ / _` | '_ \| __| | | |
#     \__/: :.. : :.. | :.. |  )        | (_| |   < (_| | | | | |_| |_| |
#          `---',\___/,\___/ /'	         \__,_|_|\_\__,_|_| |_|\__|\__,_|
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
# Library
#===============================================================================

#add_definitions(-DAKANTU_NDEBUG)

set(AKANTU_COMMON_SRC
  common/extern.cc
  common/static_memory.cc
  common/memory.cc
  fem/mesh.cc
  mesh_io/mesh_io_msh.cc
  )

set(AKANTU_INCLUDE_DIRS
  common
  fem/
  mesh_io/
  )

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  )

add_library(akantu ${AKANTU_COMMON_SRC})
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})



#===============================================================================
# Tests
#===============================================================================

option(BUILD_UNIT_TEST "Build the unit test")
if(BUILD_UNIT_TEST)
  set(AKANTU_UNIT_TESTS
    test/static_memory.cc
    test/vector.cc
    )

  foreach(_test ${AKANTU_UNIT_TESTS})
    get_filename_component(_target "${_test}" NAME_WE)
    add_executable(test_${_target} ${_test})
    target_link_libraries(test_${_target} akantu)
  endforeach(_test)
endif(BUILD_UNIT_TEST)