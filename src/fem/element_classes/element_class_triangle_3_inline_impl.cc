/**
 * @file   element_class_triangle_3_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _triangle_3
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
       \eta
	 ^
	 |
	 x (0,0,1)
	 |`
	 |  `
	 |  q `
	 |  °   `
	 x--------x----->  \xi
    (1,0,0)      (0,1,0)
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta \\
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 1/3 \qquad  \eta_{q0} = 1/3
 * @f}
 */

  // /// shape functions
  // shape[0] = 1./3.; /// N1(q_0)
  // shape[1] = 1./3.; /// N2(q_0)
  // shape[2] = 1./3.; /// N3(q_0)


/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_3>::nb_nodes_per_element;
template<> UInt ElementClass<_triangle_3>::nb_quadrature_points;
template<> UInt ElementClass<_triangle_3>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_triangle_3>::computeShapes(const Real * natural_coords,
							    Real * shapes){

  /// Natural coordinates
  Real c0 = 1 - natural_coords[0] - natural_coords[1]; /// @f$ c0 = 1 - \xi - \eta @f$
  Real c1 = natural_coords[0];                         /// @f$ c1 = \xi @f$
  Real c2 = natural_coords[1];                         /// @f$ c2 = \eta @f$

  shapes[0] = c0; /// N1(q_0)
  shapes[1] = c1; /// N2(q_0)
  shapes[2] = c2; /// N3(q_0)
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_triangle_3>::computeDNDS(__attribute__ ((unused)) const Real * natural_coords,
							       Real * dnds){

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial \xi}  & \frac{\partial N3}{\partial \xi} \\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial \eta} & \frac{\partial N3}{\partial \eta}
   *          \end{array}
   *        \right)
   * @f]
   */

  dnds[0] = -1.; dnds[1] =  1.; dnds[2] =  0.;
  dnds[3] = -1.; dnds[4] =  0.; dnds[5] =  1.;
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_triangle_3>::computeJacobian(const Real * dxds,
								   const UInt dimension,
								   Real & jac){
  if (dimension == spatial_dimension){
    Real det_dxds = Math::det2(dxds);
    jac = det_dxds;
  }
  else {
    Real vprod[dimension];
    Math::vectorProduct3(dxds,dxds+3,vprod);
    jac = Math::norm3(vprod);
  }
}



/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_triangle_3>::getInradius(const Real * coord) {
  return Math::triangle_inradius(coord, coord+2, coord+4);
}

/* -------------------------------------------------------------------------- */

template<> inline bool ElementClass<_triangle_3>::contains(const types::RVector & natural_coords) {
  if (natural_coords[0] < 0.) return false;
  if (natural_coords[0] > 1.) return false;
  if (natural_coords[1] < 0.) return false;
  if (natural_coords[1] > 1.) return false;
  if (natural_coords[0]+natural_coords[1] > 1.) return false;
  return true;
}
/* -------------------------------------------------------------------------- */


