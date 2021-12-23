#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sat Mar 13 2021
# @date last modification: Sun Mar 14 2021
#
# @brief  CMakeLists for structural mechanics examples
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
register_example(structural_mechanics
  SCRIPT structural_mechanics.py
  PYTHON
  )

register_example(structural_mechanics_softening
  SCRIPT structural_mechanics_softening.py
  PYTHON
  )

register_example(structural_mechanics_dynamics
  SCRIPT structural_mechanics_dynamics.py
  PYTHON
  )
