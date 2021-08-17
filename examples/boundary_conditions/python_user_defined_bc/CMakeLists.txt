#===============================================================================
# @file   CMakeLists.txt
#
# @author Mathias Lebihain <mathias.lebihain@enpc.fr>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Feb 18 2020
# @date last modification: Wed Feb 26 2020
#
# @brief  Python examples
#
#
# @section LICENSE
#
# Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


add_mesh(python_user_defined_bc_mesh fine_mesh.geo 2 1)

register_example(python_user_defined_bc
  SOURCES python_user_defined_bc.cc
  DEPENDS python_user_defined_bc_mesh
  PYTHON
  FILES_TO_COPY boundary_condition.py material.dat
  )

target_include_directories(python_user_defined_bc
  PRIVATE ${PROJECT_SOURCE_DIR}/python)

target_link_libraries(python_user_defined_bc
  PRIVATE pybind11::embed)
