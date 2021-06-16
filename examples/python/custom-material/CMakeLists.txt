#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Jan 22 2016
# @date last modification: Thu Feb 20 2020
#
# @brief  CMakeLists for custom-material examples
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


add_mesh(square square.geo DIM 2)
register_example(bi-material
  SCRIPT bi-material.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS square)

add_mesh(bar bar.geo DIM 2)
register_example(custom-material
  SCRIPT custom-material.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS bar)
