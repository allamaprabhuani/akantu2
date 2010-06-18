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

project(MYFEM)

enable_language(CXX)

#===========================================================================
# Misc.
#===========================================================================

set(MYFEM_CMAKE_DIR "${MYFEM_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

#===========================================================================
# Version Number
#===========================================================================

# MYFEM version number.  An even minor number corresponds to releases.
set(MYFEM_MAJOR_VERSION 0)
set(MYFEM_MINOR_VERSION 1)
set(MYFEM_BUILD_VERSION 0)
set(MYFEM_VERSION
  "${MYFEM_MAJOR_VERSION}.${MYFEM_MINOR_VERSION}.${MYFEM_BUILD_VERSION}"
  )

# Append the library version information to the library target properties
if(NOT MYFEM_NO_LIBRARY_VERSION)
  set(MYFEM_LIBRARY_PROPERTIES ${MYFEM_LIBRARY_PROPERTIES}
    VERSION "${MYFEM_VERSION}"
    SOVERSION "${MYFEM_MAJOR_VERSION}.${MYFEM_MINOR_VERSION}"
    )
endif(NOT MYFEM_NO_LIBRARY_VERSION)



#===========================================================================
# Library
#===========================================================================

#add_definitions(-DMYFEM_NDEBUG)

set(MYFEM_COMMON_SRC
  common/extern.cc
  common/static_memory.cc
  common/memory.cc
  )

include_directories(
  common
  )

add_library(myfem ${MYFEM_COMMON_SRC})
set_target_properties(myfem PROPERTIES ${MYFEM_LIBRARY_PROPERTIES})



#===========================================================================
# Tests
#===========================================================================

option(BUILD_UNIT_TEST "Build the unit test")
if(BUILD_UNIT_TEST)
  add_executable(test_static_memory test/static_memory.cc)
  target_link_libraries(test_static_memory myfem)
endif(BUILD_UNIT_TEST)