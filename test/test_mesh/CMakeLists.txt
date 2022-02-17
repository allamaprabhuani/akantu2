#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Feb 12 2018
# @date last modification:  Fri Nov 02 2018
#
# @brief  CMakeList of the mesh tests
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_mesh(cube_periodic cube_periodic.geo DIM 3 ORDER 1)
add_mesh(square_periodic square_periodic.geo DIM 2 ORDER 1)

register_test(test_mesh_periodic
  SOURCES test_mesh_periodic.cc
  DEPENDS cube_periodic square_periodic
  PACKAGE core
  UNSTABLE)
