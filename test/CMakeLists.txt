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

INCLUDE(${AKANTU_CMAKE_DIR}/AkantuTestAndExamples.cmake)

#===============================================================================
# List of tests
#===============================================================================
add_akantu_test(test_vector "Test akantu vector")
add_akantu_test(test_static_memory "Test static memory")
add_akantu_test(test_fem "Test finite element functionalties")
add_akantu_test(test_mesh_utils "Test mesh utils")
add_akantu_test(test_model "Test model objects")
add_akantu_test(test_solver "Test solver function")

if(AKANTU_CONTACT)
  #add_akantu_test(test_contact_neighbor_structure "Test contact neighbor structure")
  add_akantu_test(test_surface_extraction "Test mesh utils surface extraction")
endif(AKANTU_CONTACT)

if(AKANTU_MPI)
  add_akantu_test(test_synchronizer "Test synchronizers")
endif(AKANTU_MPI)

