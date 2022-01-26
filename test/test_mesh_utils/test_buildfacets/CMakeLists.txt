#===============================================================================
# @file   CMakeLists.txt
#
# @author Mauro Corrado <mauro.corrado@epfl.ch>
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date creation: Sun Oct 19 2014
# @date last modification:  Mon Sep 28 2015
#
# @brief  configuration for build facets test
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

set(_meshes
  _triangle_3.msh
  _triangle_6.msh
  _quadrangle_4.msh
  _quadrangle_8.msh
  _quadrangle_4_triangle_3.msh
  _quadrangle_8_triangle_6.msh
  _tetrahedron_10.msh
  _hexahedron_8.msh
  _hexahedron_20.msh
  _pentahedron_6.msh
  _pentahedron_15.msh
  _hexahedron_8_pentahedron_6.msh
  _hexahedron_20_pentahedron_15.msh
  )

register_gtest_sources(
  SOURCES test_buildfacets.cc
  PACKAGE core
  )

register_gtest_test(test_buildfacets
  FILES_TO_COPY ${_meshes}
  )
