#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Feb 18 2020
#
# @brief  CMakeLists for the cohesive examples
#
#
# @section LICENSE
#
# Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_example(predefined_bc
  "Example of predefined boundary condition" PACKAGE core)
add_example(user_defined_bc
  "Example of user defined boundary condition" PACKAGE core)
add_example(python_user_defined_bc
  "Example of python user defined boundary condition" PACKAGE pybind11)
