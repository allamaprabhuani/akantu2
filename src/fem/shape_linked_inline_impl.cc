/**
 * @file   shape_linked_inline_impl.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Apr  7 18:15:24 2011
 *
 * @brief ShapeLinked inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
inline const Vector<Real> & ShapeLinked::getShapes(const ElementType & type,
						   const GhostType & ghost_type,
						   UInt id) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(shapes.exists(type, ghost_type),
		      "No shapes of type "
		      << type << " in " << this->id);
  AKANTU_DEBUG_OUT();
  return *(shapes(type, ghost_type)[id]);
}

/* -------------------------------------------------------------------------- */
inline const Vector<Real> & ShapeLinked::getShapesDerivatives(const ElementType & type,
							      const GhostType & ghost_type,
							      UInt id) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(shapes_derivatives.exists(type, ghost_type),
		      "No shapes_derivatives of type "
		      << type << " in " << this->id);
  AKANTU_DEBUG_OUT();
  return *(shapes_derivatives(type, ghost_type)[id]);
}
