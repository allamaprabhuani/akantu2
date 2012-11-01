/**
 * @file   contact_search_2d_explicit_inline_impl.cc
 *
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 *
 * @date   Fri Nov 19 14:23:18 2010
 *
 * @brief  Inline functions declaration of class ContactSearch2dExplicit
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
// inline void ContactSearch2d::vector2(const Real * x, const Real * y, Real *vec) {
//   vec[0] = y[0]-x[0];
//   vec[1] = y[1]-x[1];
// }

/* -------------------------------------------------------------------------- */
template <class T>
inline Int ContactSearch2dExplicit::getSign(T v) {
  return (v > 0) - (v < 0);
}

/* -------------------------------------------------------------------------- */
// inline Real ContactSearch2dExplicit::F_LINE(UInt node1, UInt node2, UInt node3) {
//   return (pos_val[node3*2]*(pos_val[node1*2+1]-pos_val[node2*2+1])-pos_val[node3*2+1]*(pos_val[node1*2]-pos_val[node2*2])+pos_val[node1*2]*pos_val[node2*2+1]-pos_val[node2*2]*pos_val[node1*2+1])
// }
