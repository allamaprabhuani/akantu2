#===============================================================================
# @file   scotch.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for scotch
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
#===============================================================================

set(AKANTU_SCOTCH_FILES
  mesh_utils/mesh_partition/mesh_partition_scotch.cc
  )

if(AKANTU_SCOTCH_ON OR AKANTU_PTSCOTCH_ON)
  set(AKANTU_PARTITIONER_ON ON)
else()
  set(AKANTU_PARTITIONER_ON OFF)
endif()

set(AKANTU_SCOTCH_TESTS
  test_mesh_partitionate_scotch
  test_mesh_partitionate_scotch_advanced
  )


option(AKANTU_USE_THIRD_PARTY_SCOTCH "Use the third-party Scotch instead of the one from the system" OFF)
mark_as_advanced(AKANTU_USE_THIRD_PARTY_SCOTCH)

if(AKANTU_USE_THIRD_PARTY_SCOTCH)
  set(AKANTU_USE_SCOTCH ON CACHE BOOL "Add Scotch support in akantu" FORCE)

  if (AKANTU_USE_OBSOLETE_GETTIMEOFDAY)
    set (SCOTCH_TIMMING_OPTION -DCOMMON_TIMING_OLD)
  endif()


  if(TARGET Scotch)
    return()
  endif()

  include(ExternalProject)

  find_package(BISON)
  find_package(FLEX)
  find_package(ZLIB)

  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(SCOTCH_ARCHITECTURE -DIDXSIZE64)
  else()
    set(SCOTCH_ARCHITECTURE)
  endif()

  configure_file(${PROJECT_SOURCE_DIR}/third-party/Scotchmake.inc.cmake
    ${PROJECT_BINARY_DIR}/third-party/Scotchmake.inc)

  set(SCOTCH_URL ${PROJECT_SOURCE_DIR}/third-party/scotch_5.1.12b_esmumps.tar.gz)
  if(NOT EXISTS ${SCOTCH_URL})
    set(SCOTCH_URL https://gforge.inria.fr/frs/download.php/28978/scotch_5.1.12b_esmumps.tar.gz)
  endif()

  ExternalProject_Add(Scotch
    PREFIX ${PROJECT_BINARY_DIR}/third-party/build/scotch
    URL ${SCOTCH_URL}
    URL_HASH MD5=e13b49be804755470b159d7052764dc0
    PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/third-party/scotch.patch
    CONFIGURE_COMMAND cmake -E copy ${PROJECT_BINARY_DIR}/third-party/Scotchmake.inc src/Makefile.inc
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -C src
    INSTALL_DIR ${PROJECT_BINARY_DIR}/third-party/lib
    INSTALL_COMMAND prefix=${PROJECT_BINARY_DIR}/third-party make -C src install
    COMMAND cmake -E copy lib/libesmumps.a ${PROJECT_BINARY_DIR}/third-party/lib
    )

  set(SCOTCH_LIBRARY ${PROJECT_BINARY_DIR}/third-party/lib/libscotch.a CACHE FILEPATH "" FORCE)
  set(SCOTCH_LIBRARY_ERR ${PROJECT_BINARY_DIR}/third-party/lib/libscotcherr.a CACHE FILEPATH "" FORCE)
  set(SCOTCH_LIBRARY_ERREXIT ${PROJECT_BINARY_DIR}/third-party/lib/libscotcherrexit.a CACHE FILEPATH "" FORCE)
  set(SCOTCH_LIBRARY_ESMUMPS ${PROJECT_BINARY_DIR}/third-party/lib/libesmumps.a CACHE FILEPATH "" FORCE)
  set(SCOTCH_INCLUDE_DIR ${PROJECT_BINARY_DIR}/third-party/include CACHE PATH "" FORCE)
  #===============================================================================
  mark_as_advanced(SCOTCH_LIBRARY)
  mark_as_advanced(SCOTCH_LIBRARY_ERR)
  mark_as_advanced(SCOTCH_LIBRARY_ERREXIT)
  mark_as_advanced(SCOTCH_LIBRARY_ESMUMPS)
  mark_as_advanced(SCOTCH_INCLUDE_DIR)
  set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR})
  if(SCOTCH_LIBRARY_ESMUMPS)
    set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY_ESMUMPS} ${SCOTCH_LIBRARIES_ALL})
  endif()
  set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for scotch" FORCE)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES ${SCOTCH_LIBRARIES})
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${SCOTCH_INCLUDE_DIR})
  set(AKANTU_SCOTCH_INCLUDE_DIR ${SCOTCH_INCLUDE_DIR})
  set(AKANTU_SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES})

  list(APPEND AKANTU_OPTION_LIST SCOTCH)
  set(SCOTCH_FOUND TRUE CACHE INTERNAL "" FORCE)
  set(AKANTU_SCOTCH ON)
else()
  add_optional_external_package(Scotch "Add Scotch support in akantu" OFF)
  #add_optional_package(PTScotch "Add PTScotch support in akantu" OFF)

  if(SCOTCH_INCLUDE_DIR)
    file(STRINGS ${SCOTCH_INCLUDE_DIR}/scotch.h SCOTCH_INCLUDE_CONTENT)
    string(REGEX MATCH "_cplusplus" _match ${SCOTCH_INCLUDE_CONTENT})
    if(_match)
      set(AKANTU_SCOTCH_NO_EXTERN ON)
      list(APPEND AKANTU_DEFINITIONS AKANTU_SCOTCH_NO_EXTERN)
    else()
      set(AKANTU_SCOTCH_NO_EXTERN OFF)
    endif()
  endif()
endif()
