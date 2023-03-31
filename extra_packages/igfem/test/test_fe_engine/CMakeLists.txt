#===============================================================================
# Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
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
# @section DESCRIPTION
#
#===============================================================================


register_test(test_interpolate
  SOURCES test_interpolate.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )

register_test(test_integrate
  SOURCES test_integrate.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )

register_test(test_fe_engine_precomputation
  SOURCES test_fe_engine_precomputation.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )

register_test(test_gradient
  SOURCES test_gradient.cc igfem_mesh_generation.cc
  FILES_TO_COPY mesh_info_igfem_triangle_5.txt mesh_info_igfem_triangle_4.txt
  PACKAGE IGFEM
  )
