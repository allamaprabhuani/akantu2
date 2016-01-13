/**
 * @file   aabb_primitive.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Tue Jun 2 2015
 * @date last modification: Tue Jun 2 2015
 *
 * @brief  Macro classe (primitive) for AABB CGAL algos
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#include "aabb_primitive.hh"

__BEGIN_AKANTU__

Triangle_primitive::Point Triangle_primitive::reference_point() const {
  return primitive.vertex(0);
}

Line_arc_primitive::Point Line_arc_primitive::reference_point() const {
  Real x = to_double(primitive.source().x());
  Real y = to_double(primitive.source().y());
  Real z = to_double(primitive.source().z());
  return Spherical::Point_3(x, y, z);
}

__END_AKANTU__
