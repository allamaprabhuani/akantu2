#===============================================================================
# @file   cgal.cmake
#
# @author Lucas Frérot <lucas.frerot@epfl.ch>
#
# @date creation: Thu Feb 19 2015
# @date last modification: Mon Mar 2 2015
#
# @brief  package description for CGAL
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

package_declare(CGAL EXTERNAL
  DESCRIPTION "Add CGAL support in akantu"
  COMPILE_FLAGS "-frounding-math"
  BOOST_COMPONENTS system
  )

package_declare_sources(CGAL
  geometry/mesh_geom_abstract.hh
  geometry/mesh_geom_abstract.cc
  geometry/mesh_geom_factory.hh
  geometry/mesh_geom_factory_tmpl.hh
  geometry/mesh_geom_container.hh
  geometry/mesh_geom_container.cc

  geometry/tree_type_helper.hh
  geometry/geom_helper_functions.hh

  geometry/aabb_primitives/triangle.hh
  geometry/aabb_primitives/tetrahedron.hh

  geometry/aabb_primitives/aabb_primitive.hh
  )
