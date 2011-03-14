/**
 * @file   element_class_quadrangle_4.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Oct 27 17:27:44 2010
 *
 * @brief   Specialization of the element_class class for the type _quadrangle_4
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
 * @verbatim
         \eta
 	  ^
 (-1,1)   |   (1,1)
     x---------x
     |    |    |
     |    |    |
   --|---------|----->  \xi
     |    |    |
     |    |    |
     x---------x
 (-1,-1)  |   (1,-1)
 @endverbatim
 *
 * @subsection shapes Shape functions
 * @f[
 * \begin{array}{lll}
 * N1 = (1 - \xi) (1 - \eta) / 4
 *       & \frac{\partial N1}{\partial \xi}  = - (1 - \eta) / 4
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi) / 4 \\
 * N2 = (1 + \xi) (1 - \eta) / 4 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta) / 4
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi) / 4 \\
 * N3 = (1 + \xi) (1 + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi) / 4 \\
 * N4 = (1 - \xi) (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = - (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = (1 - \xi) / 4 \\
 * \end{array}
 * @f]
 *
 * @subsection quad_points Position of quadrature points
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 0 \qquad  \eta_{q0} = 0
 * @f}
 */


/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_quadrangle_4>::nb_nodes_per_element;
template<> UInt ElementClass<_quadrangle_4>::nb_quadrature_points;
template<> UInt ElementClass<_quadrangle_4>::spatial_dimension;

/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_4>::computeShapes(const Real * natural_coords,
								   Real * shapes) {
  /// Natural coordinates
  const Real * c = natural_coords;

  shapes[0] = .25 * (1 - c[0]) * (1 - c[1]); /// N1(q_0)
  shapes[1] = .25 * (1 + c[0]) * (1 - c[1]); /// N2(q_0)
  shapes[2] = .25 * (1 + c[0]) * (1 + c[1]); /// N3(q_0)
  shapes[3] = .25 * (1 - c[0]) * (1 + c[1]); /// N4(q_0)
}
/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_4>::computeDNDS(const Real * natural_coords,
								 Real * dnds) {

  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial \xi}
   *               & \frac{\partial N3}{\partial \xi}  & \frac{\partial N4}{\partial \xi}\\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial \eta}
   *               & \frac{\partial N3}{\partial \eta} & \frac{\partial N4}{\partial \eta}
   *          \end{array}
   *        \right)
   * @f]
   */

  const Real * c = natural_coords;

  dnds[0] = - .25 * (1 - c[1]);
  dnds[1] =   .25 * (1 - c[1]);
  dnds[2] =   .25 * (1 + c[1]);
  dnds[3] = - .25 * (1 + c[1]);

  dnds[4] = - .25 * (1 - c[0]);
  dnds[5] = - .25 * (1 + c[0]);
  dnds[6] =   .25 * (1 + c[0]);
  dnds[7] =   .25 * (1 - c[0]);
}


/* -------------------------------------------------------------------------- */
template <> inline void ElementClass<_quadrangle_4>::computeJacobian(const Real * dxds,
								   const UInt dimension,
								   Real & jac){
  Real weight = 4.;
  if (dimension == spatial_dimension){
    Real det_dxds = Math::det2(dxds);
    jac = det_dxds * weight;
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
}



/* -------------------------------------------------------------------------- */
template<> inline Real ElementClass<_quadrangle_4>::getInradius(const Real * coord) {
  Real a = Math::distance_2d(coord + 0, coord + 2);
  Real b = Math::distance_2d(coord + 2, coord + 4);
  Real c = Math::distance_2d(coord + 4, coord + 6);
  Real d = Math::distance_2d(coord + 6, coord + 0);

  // Real septimetre = (a + b + c + d) / 2.;

  // Real p = Math::distance_2d(coord + 0, coord + 4);
  // Real q = Math::distance_2d(coord + 2, coord + 6);

  // Real area = sqrt(4*(p*p * q*q) - (a*a + b*b + c*c + d*d)*(a*a + c*c - b*b - d*d)) / 4.;
  // Real h = sqrt(area);  // to get a length
  // Real h = area / septimetre;  // formula of inradius for circumscritable quadrelateral
  Real h = std::min(a, std::min(b, std::min(c, d)));

  return h;
}
