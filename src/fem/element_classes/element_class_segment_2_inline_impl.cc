/**
 * @file   element_class_segment_2_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _segment_2
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
 * @section DESCRIPTION
 *
 * @verbatim
              q
   --x--------|--------x---> x
    -1        0        1
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * w_1(x) &=& 1/2(1 - x) \\
 * w_2(x) &=& 1/2(1 + x)
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * x_{q}  &=& 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_segment_2>::nb_nodes_per_element;
template<> UInt ElementClass<_segment_2>::nb_quadrature_points;
template<> UInt ElementClass<_segment_2>::spatial_dimension;


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_segment_2>::computeShapes(const Real * natural_coords, 
								Real * shapes){

  /// natural coordinate
  Real c = natural_coords[0];
  /// shape functions
  shapes[0] = 0.5*(1-c);
  shapes[1] =0.5*(1+c);
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_segment_2>::computeDNDS(__attribute__ ((unused)) 
							   const Real * natural_coords,
							      Real * dnds){
  /// dN1/de
  dnds[0] = - .5;
  /// dN2/de
  dnds[1] =   .5;
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_segment_2>::computeJacobian(const Real * dxds,
								   const UInt dimension, 
								   Real & jac){

  if (dimension == spatial_dimension){
    jac = dxds[0];
  }
  else {
    jac = Math::norm2(dxds);
  }
}
 
/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_segment_2>::getInradius(const Real * coord) {
  return sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
}

