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

add_optional_package(Scotch "Add Scotch support in akantu" OFF)
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
  )
