/**
 * @file   fem_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 20 23:40:43 2010
 *
 * @brief  Implementation of the inline functions of the FEM Class
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
inline Mesh & FEM::getMesh() const {
  return *mesh;
}


/* -------------------------------------------------------------------------- */
inline Real FEM::getElementInradius(Real * coord, const ElementType & type) {
  AKANTU_DEBUG_IN();

  Real inradius = 0;

#define GET_INRADIUS(type)						\
  inradius = ElementClass<type>::getInradius(coord);			\

  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_INRADIUS);
#undef GET_INRADIUS

  AKANTU_DEBUG_OUT();
  return inradius;
}

// /* -------------------------------------------------------------------------- */
// inline UInt FEM::getNbQuadraturePoints(const ElementType & type) {
//   AKANTU_DEBUG_IN();

//   UInt nb_quad_points = 0;

// #define GET_NB_QUAD(type)
//   nb_quad_points = ElementClass<type>::getNbQuadraturePoints();

//   AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(GET_NB_QUAD);
// #undef GET_NB_QUAD

//   AKANTU_DEBUG_OUT();
//   return nb_quad_points;
// }

