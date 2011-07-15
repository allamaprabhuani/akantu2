/**
 * @file   shape_functions_inline_impl.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Apr  1 13:57:05 2011
 *
 * @brief  
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
inline UInt ShapeFunctions::getShapeSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_size = 0;
#define GET_SHAPE_SIZE(type)				\
  shape_size = ElementClass<type>::getShapeSize()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPE_SIZE);
#undef GET_SHAPE_SIZE

  AKANTU_DEBUG_OUT();
  return shape_size;
}

/* -------------------------------------------------------------------------- */
inline UInt ShapeFunctions::getShapeDerivativesSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_derivatives_size = 0;
#define GET_SHAPE_DERIVATIVES_SIZE(type)				\
  shape_derivatives_size = ElementClass<type>::getShapeDerivativesSize()

  AKANTU_BOOST_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_SIZE);
#undef GET_SHAPE_DERIVATIVES_SIZE

  AKANTU_DEBUG_OUT();
  return shape_derivatives_size;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeFunctions::setControlPointsByType(Vector<Real> & points, const GhostType & ghost_type){
  control_points(type, ghost_type) = &points;
}
