#===============================================================================
# @file   parallel.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Oct 16 14:05:02 2012
#
# @brief  meta package description for parallelization
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

add_meta_package(PARALLEL "Add parallel support in Akantu" OFF MPI SCOTCH)

set(AKANTU_PARALLEL_TESTS
  test_solid_mechanics_model_bar_traction2d_parallel
  test_solid_mechanics_model_segment_parallel
  test_solid_mechanics_model_pbc_parallel
  test_synchronizer_communication
  test_dof_synchronizer
  )

set(AKANTU_PARALLEL_MANUAL_FILES
  manual-parallel.tex
  )
