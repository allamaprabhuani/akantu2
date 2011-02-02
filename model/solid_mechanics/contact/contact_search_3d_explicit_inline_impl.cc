/**
 * @file   contact_search_3d_explicit_inline_impl.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Nov 25 15:22:31 2010
 *
 * @brief  Implementation of inline functions of the explicit 3d contact search algorithm
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
inline Real ContactSearch3dExplicit::computeSquareDistanceBetweenNodes(const UInt node_1, const UInt node_2) {
  AKANTU_DEBUG_IN();
  
  Real square_distance = 0.0;
  
  /// get the spatial dimension and current position of nodes
  //UInt spatial_dimension = contact.getModel().getSpatialDimension();
  Real * current_position = contact.getModel().getCurrentPosition().values; 

  /// compute the square distance between the nodes
  for(UInt dim = 0; dim < spatial_dimension; ++dim) {
    Real tmp_value = current_position[node_1 * spatial_dimension + dim]
                   - current_position[node_2 * spatial_dimension + dim];
    square_distance += tmp_value * tmp_value;
  }

  AKANTU_DEBUG_OUT();
  return square_distance;
}
