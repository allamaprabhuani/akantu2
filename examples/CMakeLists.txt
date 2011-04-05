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

#===============================================================================
# List of examples
#===============================================================================
add_example(2d_compression "Compression example")

if(AKANTU_SCOTCH_ON AND AKANTU_MPI_ON)
  add_example(3d_cube_compression "Cube compression example")
  add_example(scalability_test "Scalability test")
  add_example(compression_strength "Compression strength")
  add_example(Damage_Hole_3D "Damage Hole 3D")
add_example(test_mazars "test_mazars")
add_example(essai_bresilien "essai_bresilien")
add_example(bresilien_elas "bresilien_elas")

endif()

if(AKANTU_EPSN_ON)
  add_example(epsn "Example of steering akantu with EPSN")
endif(AKANTU_EPSN_ON)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/local_examples.cmake")
  include("local_examples.cmake")
endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/local_examples.cmake")
