#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Dec 11 2020
#
# @brief  configuration for contact detection tests
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
# @section DESCRIPTION
#
#===============================================================================


set(_meshes)

add_mesh(detection_seg data/cohesive_1D.geo
  DIM 2 ORDER 1
  OUTPUT _detection_segment_2.msh)
list(APPEND _meshes detection_seg)

add_mesh(detection_tri data/cohesive_strait_2D.geo
  DIM 2 ORDER 1
  OUTPUT _detection_triangle_3.msh)
list(APPEND _meshes detection_tri)

add_mesh(detection_quad data/cohesive_strait_2D_structured.geo
  DIM 2 ORDER 1
  OUTPUT _detection_quadrangle_4.msh)
list(APPEND _meshes detection_quad)

add_mesh(detection_tri_quad data/cohesive_strait_2D_mixte.geo
  DIM 2 ORDER 1
  OUTPUT _detection_triangle_3_quadrangle_4.msh)
list(APPEND _meshes detection_tri_quad)

add_mesh(detection_tet data/cohesive_strait_3D.geo
  DIM 3 ORDER 1
  OUTPUT _detection_tetrahedron_4.msh)
list(APPEND _meshes detection_tet)

register_gtest_sources(
  SOURCES test_detection.cc
  PACKAGE contact_mechanics
  DEPENDS ${_meshes}
  FILES_TO_COPY material_0.dat material_1.dat
  )

