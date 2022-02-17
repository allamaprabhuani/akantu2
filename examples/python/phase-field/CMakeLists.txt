#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Aug 03 2021
# @date last modification: Tue Aug 03 2021
#
# @brief  CMake file for phase-field python example
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


add_mesh(mesh_phasefield plate.geo DIM 2)
add_mesh(mesh_phasefield_static plate_static.geo DIM 2)

register_example(phasefield_dynamic
  SCRIPT phasefield-dynamic.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS mesh_phasefield
  )

register_example(phasefield_static
  SCRIPT phasefield-static.py
  PYTHON
  FILES_TO_COPY material_static.dat
  DEPENDS mesh_phasefield_static
  )
