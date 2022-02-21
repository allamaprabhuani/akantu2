/**
 * @file   geometry_utils_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Oct 06 2019
 * @date last modification: Wed Sep 16 2020
 *
 * @brief  Geometry utils
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "geometry_utils.hh"

#ifndef __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__
#define __AKANTU_GEOMETRY_UTILS_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool GeometryUtils::isBoundaryElement(const Mesh & mesh,
                                             const Element & subelement) {

  const auto & element_to_subelement =
      mesh.getElementToSubelement(subelement.type)(subelement.element);

  // for regular boundary elements when surfaceselector is set to
  // physical surfaces, the mesh contains only 1 element attached to a
  // boundary subelement
  if (element_to_subelement.size() == 1 and
      element_to_subelement[0].kind() == _ek_regular) {
    return true;
  }

  // for cohesive interface elements when surfaceSelector is set
  // either cohesive surface selector or all surface selector, in this
  // case mesh passed is actually mesh_facet and for boundary or
  // cohesive  interface 2 elements are associated to a subelement
  // we want only one regular element attached to the subelement

  UInt nb_elements_regular = 0;
  UInt nb_elements_cohesive = 0;

  for (auto elem : element_to_subelement) {
    if (elem == ElementNull) {
      continue;
    }

    if (elem.kind() == _ek_regular) {
      ++nb_elements_regular;
    }

    if (elem.kind() == _ek_cohesive) {
      ++nb_elements_cohesive;
    }
  }

  auto nb_elements = element_to_subelement.size();
  return nb_elements_regular < nb_elements;
}

/* -------------------------------------------------------------------------- */
inline bool GeometryUtils::isValidProjection(const Vector<Real> & projection,
                                             Real extension_tolerance) {

  UInt nb_xi_inside = 0;

  for (auto xi : projection) {
    if (xi >= -1.0 - extension_tolerance and xi <= 1.0 + extension_tolerance) {
      nb_xi_inside++;
    }
  }

  return nb_xi_inside == projection.size();
}

} // namespace akantu

#endif
