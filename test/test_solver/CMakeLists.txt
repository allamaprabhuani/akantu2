#===============================================================================
# @file   CMakeLists.txt
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
# @brief configuration for solver tests
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

add_mesh(test_solver_mesh triangle.geo 2 1)

register_test(test_sparse_matrix_profile
  SOURCES test_sparse_matrix_profile.cc
  DEPENDENCIES test_solver_mesh
  )

register_test(test_sparse_matrix_assemble
  SOURCES test_sparse_matrix_assemble.cc
  DEPENDENCIES test_solver_mesh
  )

register_test(test_solver_mumps SOURCES test_solver_mumps.cc)

register_test(test_sparse_matrix_product
  SOURCES test_sparse_matrix_product.cc
  FILES_TO_COPY bar.msh
  )
