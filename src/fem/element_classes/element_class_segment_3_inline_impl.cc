/**
 * @file   element_class_segment_3_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Jul 16 09:09:21 2010
 *
 * @brief  Specialization of the element_class class for the type _segment_3
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
         -1         0         1
  	 -----x---------x---------x-----> x
	      1         3         2
 @endverbatim
 *
 * @subsection coords Nodes coordinates
 *
 * @f[
 * \begin{array}{lll}
 *   x_{1}  = -1   &  x_{2} = 1 & x_{3} = 0
 * \end{array}
 * @f]
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{ll}
 *   w_1(x) = \frac{x}{2}(x - 1) & w'_1(x) = x - \frac{1}{2}\\
 *   w_2(x) = \frac{x}{2}(x + 1) & w'_2(x) = x + \frac{1}{2}\\
 *   w_3(x) =  1-x^2 & w'_3(x) = -2x
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f[
 * \begin{array}{ll}
 * x_{q1}  = -1/\sqrt{3} & x_{q2} = 1/\sqrt{3}
 * \end{array}
 * @f]
 */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_segment_3>::nb_nodes_per_element;
template<> UInt ElementClass<_segment_3>::nb_quadrature_points;
template<> UInt ElementClass<_segment_3>::spatial_dimension;



/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_segment_3>::computeShapes(const Real * natural_coords, 
								Real * shapes){
  Real c = natural_coords[0];
  shapes[0] = (c - 1) * c / 2;
  shapes[1] = (c + 1) * c / 2;
  shapes[2] = 1 - c * c;
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_segment_3>::computeDNDS(const Real * natural_coords,
							      Real * dnds){

  Real c = natural_coords[0];
  dnds[0]  = c - .5;
  dnds[1]  = c + .5;
  dnds[2]  = -2 * c;
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_segment_3>::computeJacobian(const Real * dxds,
							       const UInt dimension, 
							       Real & jac){
  if (dimension == spatial_dimension){
    jac = dxds[0];
  } else {
    jac = Math::norm2(dxds);
  }
}

/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_segment_3>::getInradius(const Real * coord) {
  Real dist1 = sqrt((coord[0] - coord[1])*(coord[0] - coord[1]));
  Real dist2 = sqrt((coord[1] - coord[2])*(coord[1] - coord[2]));
  return dist1 < dist2 ? dist1 : dist2;
}
