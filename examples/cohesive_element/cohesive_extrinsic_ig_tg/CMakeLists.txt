#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Tue Feb 18 2020
#
# @brief  compilation of Cohesive element examples in extrinsic with 2 different bulk
# materials
#
#
# @section LICENSE
#
# Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_mesh(cohesive_extrinsic_ig_tg_mesh triangle.geo DIM 2 ORDER 1 OUTPUT square.msh)
register_example(cohesive_extrinsic_ig_tg
  SOURCES cohesive_extrinsic_ig_tg.cc
  DEPENDS cohesive_extrinsic_ig_tg_mesh
  FILES_TO_COPY material.dat
  )
