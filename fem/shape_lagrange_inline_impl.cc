/**
 * @file   shape_lagrange_inline_impl.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Feb 10 21:12:54 2011
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
template <ElementType type>
inline void ShapeLagrange::
computeShapeDerivativesOnCPointsByElement(UInt spatial_dimension,
					  Real * node_coords,
					  UInt nb_nodes_per_element,
					  Real * natural_coords,
					  UInt nb_points,
					  Real * shapesd) {
  // compute dnds
  Real dnds[nb_nodes_per_element * spatial_dimension * nb_points];
  ElementClass<type>::computeDNDS(natural_coords, nb_points, dnds);
  // compute dxds
  Real dxds[spatial_dimension * spatial_dimension * nb_points];
  ElementClass<type>::computeDXDS(dnds, nb_points, node_coords,
				  spatial_dimension, dxds);
  // compute shape derivatives
  ElementClass<type>::computeShapeDerivatives(dxds, dnds, nb_points,
					      spatial_dimension, shapesd);
}
