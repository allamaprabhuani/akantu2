#===============================================================================
# @file   mumps.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for mumps support
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
set(AKANTU_MUMPS_FILES
  solver/solver_mumps.cc
  solver/solver_mumps.hh
  )

option(AKANTU_USE_THIRD_PARTY_MUMPS "Use the third-party Mumps instead of the one from the system" OFF)
mark_as_advanced(AKANTU_USE_THIRD_PARTY_MUMPS)
if(AKANTU_USE_THIRD_PARTY_MUMPS)
  set(AKANTU_USE_MUMPS ON CACHE BOOL "Add Mumps support in akantu" FORCE)
  set(MUMPS_DEPENDS)

  include(ExternalProject)
  if(AKANTU_USE_MPI)
    string(REPLACE ";" " -I" MUMPS_MPI_INCLUDE_PATH "-I${MPI_C_INCLUDE_PATH}")
    string(REPLACE ";" " " MUMPS_MPI_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES}")

    ExternalProject_Add(ScaLAPACK
      PREFIX ${PROJECT_BINARY_DIR}/third-party/build/scalapack
      INSTALL_DIR ${PROJECT_BINARY_DIR}/third-party/lib
      URL http://www.netlib.org/scalapack/scalapack-2.0.2.tgz
#      URL_HASH MD5=2f75e600a2ba155ed9ce974a1c4b536f
      CMAKE_ARGS ../ScaLAPACK
      CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/third-party -DCMAKE_C_FLAGS:STRING=-fPIC -DCMAKE_Fortran_FLAGS:STRING=-fPIC
      )
    set(SCALAPACK_LIBRARY ${PROJECT_BINARY_DIR}/third-party/lib/libscalapack.a)
    list(APPEND MUMPS_DEPENDS ScaLAPACK)
    set(MUMPS_TYPE par)
  else()
    set(MUMPS_TYPE seq)
  endif()


  set(MUMPS_LIBRARIES_ALL)
  if(AKANTU_USE_THIRD_PARTY_SCOTCH)
    if(NOT TARGET Scotch)
      include(${PROJECT_SOURCE_DIR}/packages/scotch.cmake)
    endif()
    list(APPEND MUMPS_DEPENDS Scotch)
    list(APPEND MUMPS_LIBRARIES_ALL ${SCOTCH_LIBRARIES})
  endif()

  if(SCALAPACK_LIBRARY)
    list(APPEND MUMPS_LIBRARIES_ALL ${SCALAPACK_LIBRARY})
  endif()

  if("${MUMPS_TYPE}" STREQUAL "seq")
    set(MUMPS_PREFIX _seq)
    set(_libmumps_seq COMMAND cmake -E copy libseq/libmpiseq${MUMPS_PREFIX}.a ${PROJECT_BINARY_DIR}/third-party/lib)
    set(MUMPS_LIBRARY_MPI ${PROJECT_BINARY_DIR}/third-party/lib/libmpiseq${MUMPS_PREFIX}.a CACHE FILEPATH "" FORCE)
    mark_as_advanced(MUMPS_LIBRARY_MPI)
  else()
    set(MUMPS_LIBRARY_MPI "" CACHE INTERNAL)
  endif()

  configure_file(${PROJECT_SOURCE_DIR}/third-party/MUMPSmake.inc.cmake
    ${PROJECT_BINARY_DIR}/third-party/MUMPSmake.inc)

  ExternalProject_Add(MUMPS
    DEPENDS ${MUMPS_DEPENDS}
    PREFIX ${PROJECT_BINARY_DIR}/third-party/build/mumps
    URL ${PROJECT_SOURCE_DIR}/third-party/MUMPS_4.9.2.tar.gz
#    URL_HASH MD5=1c896cdb61878cf094b779404c6512fd
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND cmake -E create_symlink ${PROJECT_BINARY_DIR}/third-party/MUMPSmake.inc Makefile.inc
    BUILD_COMMAND make d
    INSTALL_DIR ${PROJECT_BINARY_DIR}/third-party/lib
    INSTALL_COMMAND cmake -E copy lib/libdmumps${MUMPS_PREFIX}.a ${PROJECT_BINARY_DIR}/third-party/lib
    COMMAND cmake -E copy lib/libmumps_common${MUMPS_PREFIX}.a ${PROJECT_BINARY_DIR}/third-party/lib
    COMMAND cmake -E copy lib/libpord${MUMPS_PREFIX}.a ${PROJECT_BINARY_DIR}/third-party/lib
    COMMAND cmake -E copy_directory include ${PROJECT_BINARY_DIR}/third-party/include
    ${_libmumps_seq}
    )

  set(MUMPS_LIBRARY_DMUMPS ${PROJECT_BINARY_DIR}/third-party/lib/libdmumps${MUMPS_PREFIX}.a CACHE FILEPATH "" FORCE)
  set(MUMPS_LIBRARY_COMMON ${PROJECT_BINARY_DIR}/third-party/lib/libmumps_common${MUMPS_PREFIX}.a CACHE FILEPATH "" FORCE)
  set(MUMPS_LIBRARY_PORD ${PROJECT_BINARY_DIR}/third-party/lib/libpord${MUMPS_PREFIX}.a CACHE FILEPATH "" FORCE)
  set(MUMPS_INCLUDE_DIR ${PROJECT_BINARY_DIR}/third-party/include CACHE PATH "" FORCE)

  mark_as_advanced(MUMPS_LIBRARY_COMMON)
  mark_as_advanced(MUMPS_LIBRARY_DMUMPS)
  mark_as_advanced(MUMPS_LIBRARY_PORD)
  mark_as_advanced(MUMPS_INCLUDE_DIR)

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    list(APPEND MUMPS_LIBRARIES_ALL -lgfortran)
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    list(APPEND MUMPS_LIBRARIES_ALL -lifcore)
  endif()

  list(APPEND MUMPS_LIBRARIES_ALL ${MPI_Fortran_LIBRARIES} ${MUMPS_LIBRARY_COMMON} ${MUMPS_LIBRARY_DMUMPS} ${MUMPS_LIBRARY_PORD} ${MUMPS_LIBRARY_MPI})
  set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "Libraries for MUMPS" FORCE)

  list(APPEND AKANTU_EXTERNAL_LIBRARIES ${MUMPS_LIBRARIES})
  list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${MUMPS_INCLUDE_DIR})
  set(AKANTU_MUMPS_INCLUDE_DIR ${MUMPS_INCLUDE_DIR})
  set(AKANTU_MUMPS_LIBRARIES ${MUMPS_LIBRARIES})
  list(APPEND AKANTU_OPTION_LIST MUMPS)
  set(MUMPS_FOUND TRUE CACHE INTERNAL "" FORCE)
  set(AKANTU_MUMPS ON)

else()
  if(AKANTU_USE_MPI)
    set(MUMPS_TYPE par)
    set(AKANTU_MUMPS_DEB_DEPEND
      libmumps-dev
      )
  else()
    set(MUMPS_TYPE seq)
    set(AKANTU_MUMPS_DEB_DEPEND
      libmumps-seq-dev
      )
  endif()

  add_optional_external_package(Mumps "Add Mumps support in akantu" OFF)

endif()

set(AKANTU_MUMPS_TESTS
  test_sparse_matrix_profile
  test_sparse_matrix_assemble
  test_solver_mumps
  test_sparse_matrix_product
  )