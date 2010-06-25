#===========================================================================
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
#===========================================================================


#===========================================================================
# CMake Project
#===========================================================================

cmake_minimum_required(VERSION 2.6)

project(AKANTU)

enable_language(CXX)

#===========================================================================
# Misc.
#===========================================================================

set(AKANTU_CMAKE_DIR "${AKANTU_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

#===========================================================================
# Version Number
#===========================================================================

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



#===========================================================================
# Library
#===========================================================================

#add_definitions(-DAKANTU_NDEBUG)

set(AKANTU_COMMON_SRC
  common/extern.cc
  common/static_memory.cc
  common/memory.cc
  fem/mesh.cc
  mesh_io/mesh_io_msh.cc
  )

include_directories(
  common
  fem/
  mesh_io/
  )

add_library(akantu ${AKANTU_COMMON_SRC})
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})



#===========================================================================
# Tests
#===========================================================================

option(BUILD_UNIT_TEST "Build the unit test")
if(BUILD_UNIT_TEST)
  add_executable(test_static_memory test/static_memory.cc)
  target_link_libraries(test_static_memory akantu)
endif(BUILD_UNIT_TEST)